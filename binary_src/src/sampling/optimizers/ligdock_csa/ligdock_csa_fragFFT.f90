!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/sampling/optimizers/ligdock_csa/ligdock_csa_fragFFT.f90
!
! Description:
!
!-------------------------------------------------------------------------------
MODULE LIGDOCK_CSA_FRAGFFT
!-------------------------------------------------------------------------------

use globals
use logger
use in_out_structure
use ran
use string
use sort, only: sort1, sort2
use mathfunctions, only: v_norm, cross, bound_ang
use geometry, only: sample_rotation_by_Fibonacci, rotation_matrix, &
                    rotmat_to_axisangle, gen_uniform_random_quat, rotation, &
                    expmap_quat,inv_expmap_quat, pert_rotation, gen_random_axis
use rmsd, only: ls_rmsd
use FFT_utils
use clustering, only: hierarchical_clustering_dcut, hierarchical_clustering, &
                      DBSCAN, DBSCAN_2
!
use energy_vars
use ligdock_energy, only: receptor_gridization, ligand_gridization, &
                          finalize_FFT_grid, ligdock_energy_using_grid
use Xscore, only: get_hydrophobicity, initialize_X_score
use energy_utils,   only: update_R_for_ligand
!
use transrot_operator, only: transrot_crd
use ligand_operator 
use simplex
!
use mkl_dfti

implicit none
save
private

integer, parameter :: max_hotspot = 20  ! maximum number of binding hot-spots per fragment
integer, parameter :: max_lig_frag = 20 ! maximum number of fragments
integer, parameter :: max_lig_frag_sq = max_lig_frag**2
integer, parameter :: max_frag_per_br = 4  ! maximum number of fragments in each branch
integer, parameter :: max_clique = 1000
integer, parameter :: max_mem = max_lig_frag_sq*max_hotspot
real(dp), parameter :: d_cut = 1.0d0
real(dp) :: cl_cut, cl_cut_sq
real(dp), parameter :: FFT_E_CUTOFF = 0.0d0

!-------------------------------------------------------------------------------
type FFT_result_type
!-------------------------------------------------------------------------------
real :: E
real :: R(3)
!-------------------------------------------------------------------------------
end type FFT_result_type
!-------------------------------------------------------------------------------
type frag_hotspot_type
!-------------------------------------------------------------------------------
integer :: n_spot
real :: R(3, max_hotspot)
real :: E(max_hotspot)
! FOR TEST: TODO remove after method development
real :: err(max_hotspot)
!-------------------------------------------------------------------------------
end type frag_hotspot_type
!-------------------------------------------------------------------------------
type hotspot_graph_type
!-------------------------------------------------------------------------------
integer :: n_mem
integer :: i_mem(max_hotspot*max_lig_frag)
integer :: hotspot_type(max_hotspot*max_lig_frag)
real :: R(3, max_hotspot*max_lig_frag)
real, dimension(max_hotspot*max_lig_frag,max_hotspot*max_lig_frag) :: dist
!TEST
real :: err(max_hotspot*max_lig_frag)
!-------------------------------------------------------------------------------
end type hotspot_graph_type
!-------------------------------------------------------------------------------
type assoc_graph_type
!-------------------------------------------------------------------------------
integer :: n_mem
integer :: i_mem(2,max_hotspot*max_lig_frag_sq)
logical :: connect(max_hotspot*max_lig_frag_sq,max_hotspot*max_lig_frag_sq)
!-------------------------------------------------------------------------------
end type assoc_graph_type
!-------------------------------------------------------------------------------
type clique_type
!-------------------------------------------------------------------------------
integer :: n_mem
integer :: i_mem(max_lig_frag)
!-------------------------------------------------------------------------------
end type clique_type
!-------------------------------------------------------------------------------

integer :: FFT_status
type(DFTI_DESCRIPTOR), pointer :: desc
type(hotspot_graph_type) :: rec_graph
integer, allocatable :: fragment(:,:), n_atm_frag(:)
integer :: n_frag
integer :: core_branch_idx_frag(max_lig_frag)

public :: max_lig_frag
public :: hotspot_graph_type
public :: rec_graph
!
public :: generate_fragment_binding_map
public :: generate_lig_graph
public :: do_pharmdock
public :: finalize_pharmdock

CONTAINS
!-------------------------------------------------------------------------------
subroutine generate_fragment_binding_map(protein, ligand)
!-------------------------------------------------------------------------------
! Detect fragment binding hot-spots using FFT-based fragment docking
! To do that:
!    1. Ligand is splitted into upto max_lig_frag rigid fragments
!    2. For each fragment, run FFT docking
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
!
type(docking_grid_param_type) :: FFT_grid_info
type(frag_hotspot_type) :: hotspot_rec(max_lig_frag)
integer :: i_frag, n_rot, gene_dim, i_cut

gene_dim = 6 + (ligand%n_br-1)
print*, 'GENE DIM: ', gene_dim
if (.not. use_Xscore) then ! hydrophobicity is needed to select ligand fragment!
    call initialize_X_score(protein, ligand)
end if

allocate(fragment(ligand%n_atm,max_lig_frag))
allocate(n_atm_frag(max_lig_frag))
!
! To run energy minimization of bound fragments, simplex should be re-initialize
! (TR sampling only).
call finalize_simplex()
call initialize_simplex(6)
!
call log_p('Prepare FFT-based initial bank generation', level=30)
call ligand_fragmentation(ligand, fragment, n_atm_frag, n_frag)
print*, 'NUMBER OF FRAGMENT: ', n_frag
call initialize_FFT_dock(FFT_grid_info)
!
do i_frag = 1, n_frag
    write(log_msg, '(A,I2,A)') ' -Running FFT-based docking of ', i_frag, 'th fragment'
    call log_p(log_msg, level=30)
    call do_FFT_docking(FFT_grid_info, protein, ligand, fragment(:,i_frag), &
                        n_atm_frag(i_frag), hotspot_rec(i_frag), i_frag)
end do
!
! re-initialize simplex for whole ligand
call finalize_simplex()
call initialize_simplex(gene_dim)
call finalize_FFT_dock()
!
! Generate fragment hot-spot graph
call generate_rec_graph(hotspot_rec, rec_graph, n_frag)
!stop

end subroutine generate_fragment_binding_map
!-------------------------------------------------------------------------------
subroutine finalize_pharmdock()
!-------------------------------------------------------------------------------
deallocate(fragment)
deallocate(n_atm_frag)

end subroutine finalize_pharmdock
!-------------------------------------------------------------------------------
subroutine ligand_fragmentation(ligand, fragment, n_atm_frag, n_frag)
!-------------------------------------------------------------------------------
! Get up to max_lig_frag rigid fragment in ligand
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand
integer, intent(out) :: fragment(:,:), n_atm_frag(:), n_frag
!
integer :: n_atm_piece(ligand%n_br), key(ligand%n_br)
integer :: i_br, j_br, br_idx, i_atm
integer :: n_frag_cand
integer :: frag_cand(ligand%n_atm, ligand%n_br, ligand%n_core_br)
integer :: n_atm_cand(ligand%n_br, ligand%n_core_br)
integer :: tmp_frag(ligand%n_atm), n_tmp_cand
integer :: root_idx, core_idx, n_frag_per_br(ligand%n_core_br)
real(dp), allocatable :: hp(:)

core_branch_idx_frag(:) = 0
root_idx = ligand%piece(ligand%atm_in_br(1,1))
n_frag_per_br(:) = 0
!
! Count number of atoms in each rigid piece
n_atm_piece(:) = 0
do i_atm = 1, ligand%n_atm
    i_br = ligand%piece(i_atm)
    n_atm_piece(i_br) = n_atm_piece(i_br) + 1
end do

call sort1(ligand%n_br, n_atm_piece, key(1:ligand%n_br))

n_frag = 0
n_atm_frag(:) = 0
!
n_atm_cand(:,:) = 0
!
do i_br = ligand%n_br, 1, -1
    br_idx = key(i_br)
    if (br_idx == root_idx) then ! including root fragment always
        n_frag = n_frag + 1
        call get_atoms_in_fragment(ligand, br_idx, n_atm_frag(n_frag), &
                                   fragment(:,n_frag))
    else
        ! Get atom index in current fragment
        n_tmp_cand = 0
        tmp_frag(:) = 0
        call get_atoms_in_fragment(ligand, br_idx, n_tmp_cand, &
                                   tmp_frag)
        call find_core_branch(ligand, tmp_frag, n_tmp_cand, core_idx)
        n_frag_per_br(core_idx) = n_frag_per_br(core_idx) + 1
        n_atm_cand(n_frag_per_br(core_idx), core_idx)   = n_tmp_cand
        frag_cand(:, n_frag_per_br(core_idx), core_idx) = tmp_frag(:)
    end if
end do

do i_br = 1, ligand%n_core_br
    n_frag_cand = n_frag_per_br(i_br)  ! number of fragment candidates per core_branch
    allocate(hp(n_frag_cand)) 
    do j_br = 1, n_frag_cand
        call get_hydrophobicity(ligand%lig_no, frag_cand(:,j_br,i_br), &
                                n_atm_cand(j_br,i_br), hp(j_br))
        hp(j_br) = -1.0d0 * abs(hp(j_br))  ! negative -> hydrophobic
                                           ! positive -> hydrophilic
    end do
    call sort2(n_frag_cand, hp(1:n_frag_cand), key(1:n_frag_cand))
    !
    n_tmp_cand = 0
    do j_br = 1, n_frag_cand
        !if (abs(hp(j_br)) < 0.3d0) cycle ! neither hydrophobic nor hydrophilic
        !                                 ! (tends to bind anywhere?)
        n_tmp_cand = n_tmp_cand + 1
        n_frag = n_frag + 1
        n_atm_frag(n_frag) = n_atm_cand(key(j_br), i_br)
        fragment(:, n_frag) = frag_cand(:,key(j_br),i_br)
        core_branch_idx_frag(n_frag) = i_br
        if (n_frag == max_lig_frag) exit
        if (n_tmp_cand == max_frag_per_br) exit
    end do
    if (n_tmp_cand == 0) then 
        n_frag = n_frag + 1
        n_atm_frag(n_frag) = n_atm_cand(key(1), i_br)
        fragment(:, n_frag) = frag_cand(:,key(1),i_br)
        core_branch_idx_frag(n_frag) = i_br
    end if
    deallocate(hp)
    if (n_frag == max_lig_frag) exit
end do

end subroutine ligand_fragmentation
!-------------------------------------------------------------------------------
subroutine get_atoms_in_fragment(ligand, br_idx, n_atm_frag, frag)
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand
integer, intent(in) :: br_idx
integer, intent(inout) :: n_atm_frag
integer, intent(inout) :: frag(:)
!
integer :: i_atm

! Get atom index in current fragment
do i_atm = 1, ligand%n_atm
    if (ligand%piece(i_atm) == br_idx) then
        n_atm_frag = n_atm_frag + 1
        frag(n_atm_frag) = i_atm
    end if
end do
call attach_bridged_atom(ligand, br_idx, n_atm_frag, frag)

end subroutine get_atoms_in_fragment
!-------------------------------------------------------------------------------
subroutine attach_bridged_atom(ligand, br_idx, n_frag_atm, frag)
!-------------------------------------------------------------------------------
! Find & add atoms connected directly to current fragment
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand
integer, intent(in) :: br_idx
integer, intent(inout) :: n_frag_atm, frag(:)
integer :: i_br

do i_br = 2, ligand%n_br
    if (ligand%piece(ligand%bridge(1,i_br)) == br_idx) then
        n_frag_atm = n_frag_atm + 1
        frag(n_frag_atm) = ligand%bridge(2,i_br)
    else if (ligand%piece(ligand%bridge(2,i_br)) == br_idx) then
        n_frag_atm = n_frag_atm + 1
        frag(n_frag_atm) = ligand%bridge(1,i_br)
    end if
end do

end subroutine attach_bridged_atom
!-------------------------------------------------------------------------------
subroutine find_core_branch(ligand, atom_s, n_atm, core_idx)
!-------------------------------------------------------------------------------
! Find core_branch (which means how to connect to the root) index
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand
integer, intent(in) :: atom_s(:), n_atm
integer, intent(out) :: core_idx
!
integer :: i_atm, j_atm, i_core, br_idx, n_match

core_idx = -1
do i_core = 1, ligand%n_core_br
    br_idx = ligand%core_bridge(i_core)
    n_match = 0
    do i_atm = 1, ligand%n_atm_br(br_idx)
        do j_atm = 1, n_atm
            if (atom_s(j_atm) == ligand%atm_in_br(i_atm, br_idx)) then
                n_match = n_match + 1
                exit
            end if
        end do
        if (n_match > 1 .or. n_match == ligand%n_atm_br(br_idx)) exit
    end do
    if (n_match > 1 .or. n_match == ligand%n_atm_br(br_idx)) then
        core_idx = i_core
        exit
    end if
end do

end subroutine find_core_branch
!-------------------------------------------------------------------------------
subroutine initialize_FFT_dock(FFT_grid_info)
!-------------------------------------------------------------------------------
! To define protein-based pharmacophore, it uses FFT docking.
! In this subroutine, grid parameter and plan for FFT are set.
! Also, FFT grid for receptor is calculated.
!-------------------------------------------------------------------------------
type(docking_grid_param_type), intent(out) :: FFT_grid_info

call setup_FFT_grid_info(FFT_grid_info)
call receptor_gridization(FFT_grid_info)
!
FFT_status=DftiCreateDescriptor(desc, DFTI_SINGLE, DFTI_COMPLEX,&
                                3, FFT_grid_info%n_elem)
FFT_status=DftiCommitDescriptor(desc)
call get_r_conv(r_grid, FFT_grid_info%n_elem)

end subroutine initialize_FFT_dock
!-------------------------------------------------------------------------------
subroutine finalize_FFT_dock()
!-------------------------------------------------------------------------------

FFT_status=DftiFreeDescriptor(desc)
call finalize_FFT_grid()

end subroutine finalize_FFT_dock
!-------------------------------------------------------------------------------
subroutine do_FFT_docking(FFT_grid_info, protein, ligand, fragment, n_atm_frag, &
                          hotspot_rec, i_frag, n_rot, n_rot_z, cutoff_cl)
!-------------------------------------------------------------------------------
type(docking_grid_param_type), intent(in) :: FFT_grid_info
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
integer, intent(in) :: fragment(:), n_atm_frag
type(frag_hotspot_type), intent(inout) :: hotspot_rec
integer, intent(in) :: i_frag
integer, intent(in), optional :: n_rot, n_rot_z
real(dp), intent(in), optional :: cutoff_cl
!
integer :: rot_sample, rot_sample_z
!
real(dp), allocatable :: p(:,:)
real(dp) :: R_frag(3, max_lig_atom), R_frag_rotated(3, max_lig_atom)
real(dp) :: U(3,3), gene(7), axis1(3), axis2(3)
integer :: i_rot, cand_idx, i_cand, i_pool, n_pool, min_idx(1)
integer, parameter :: n_pool_per_rot = 3
!
type(FFT_result_type), allocatable :: cand_s(:)
!DEBUG
real(dp) :: R_frag_debug(3, max_lig_atom)
real(dp) :: tmp_R(3), tmp_vec1(3), dist, error, hp

if (.not. present(n_rot)) then
    !rot_sample = 60
    rot_sample = 100
else
    rot_sample = n_rot
end if
!
if (.not. present(n_rot_z)) then
    !rot_sample_z = 12 ! 30 degree
    rot_sample_z = 24 ! 15 degree
else
    rot_sample = n_rot_z 
end if
!
if (.not. present(cutoff_cl)) then
    !cl_cut = 2.0d0
    cl_cut = 1.5d0
else
    cl_cut = cutoff_cl
end if
cl_cut_sq = cl_cut ** 2
!
allocate(cand_s(rot_sample*rot_sample_z*n_pool_per_rot))
allocate(p(4,rot_sample*rot_sample_z))
!
call update_R_frag(fragment, n_atm_frag, ligand, R_frag)
call get_frag_cntr(protein%ligand(ligand%lig_no), &
                   fragment, n_atm_frag, tmp_R, tmp_vec1)
call get_hydrophobicity(ligand%lig_no, fragment, n_atm_frag, hp) ! For TEST
                                                                 ! TODO: remove 
write(*,'(A,I2,I3,7F8.3)'), 'Fragment: ', i_frag, n_atm_frag, tmp_R(:), tmp_vec1(:), hp
! DEBUG
call update_R_frag(fragment, n_atm_frag, ligand, R_frag_debug, .false.)
call write_native_frag(R_frag_debug, ligand, fragment, n_atm_frag, i_frag)
!
call sample_rotation_by_Fibonacci(p, rot_sample_z, rot_sample, tmp_vec1)
call ligand_gridization(FFT_grid_info, ligand, fragment, &
                        n_atm_frag, R_frag_rotated, .true.)
do i_rot = 1, rot_sample*rot_sample_z
    i_cand = (i_rot-1)*n_pool_per_rot
    call rotation_matrix(3, p(:,i_rot), U)
    call rotate_fragment(R_frag, R_frag_rotated, n_atm_frag, U)
    call ligand_gridization(FFT_grid_info, ligand, fragment, &
                            n_atm_frag, R_frag_rotated, .false.)
    call translate_fragment(r_grid, FFT_grid_info%n_elem, &
                            cand_s(i_cand+1 : i_cand+n_pool_per_rot), &
                            n_pool_per_rot)
end do
!
!! DEBUG
!call write_FFT_result(cand_s, rot_sample*n_pool_per_rot, i_frag)
! 
! get top conf among candidate
call update_hotspot_rec(hotspot_rec, cand_s, rot_sample*rot_sample_z*n_pool_per_rot, n_atm_frag, &
                        i_frag) 
!
do i_rot = 1, hotspot_rec%n_spot
    gene(1:3) = hotspot_rec%R(:,i_rot) - tmp_R(:)
    dist = sqrt(dot_product(gene,gene))
    write(*,'(A, I2,4F8.3,F10.3)') 'Site ', i_rot, hotspot_rec%R(:,i_rot), &
                    dist, hotspot_rec%E(i_rot)
    hotspot_rec%err(i_rot) = dist
end do
!
call write_pharm_result(hotspot_rec, i_frag)
deallocate(cand_s)

end subroutine do_FFT_docking
!-------------------------------------------------------------------------------
subroutine write_native_frag(R_frag, ligand, fragment, n_atm, i_frag)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: R_frag(:,:)
type(ligand_type), intent(in) :: ligand
integer, intent(in) :: fragment(:), n_atm, i_frag
!
character(len=3) :: num_frag
character(len=len_fname) :: pdb_file
integer, parameter :: pdb_unit = 75
integer :: i_atm

call num2str_999(i_frag, num_frag)
write(pdb_file, '(A,A,A)') 'ans_', num_frag, '.pdb'
call open_write_pdb(pdb_unit, pdb_file)

do i_atm = 1, n_atm
    write(pdb_unit, 73) 'HETATM', i_atm, &
        ref_lig(ligand%lig_type)%atom_name(fragment(i_atm)), &
        ref_lig(ligand%lig_type)%lig_name, & 
        ' ', i_frag, ' ', R_frag(:,i_atm) 
end do

call close_write_pdb(pdb_unit)

73  format(a6,i5,1x,a4,1x,a3,1x,a1,i4,a1,3x,3f8.3)

end subroutine write_native_frag
!-------------------------------------------------------------------------------
subroutine write_FFT_result(cand_s, n_cand, i_frag)
!-------------------------------------------------------------------------------
type(FFT_result_type), intent(in) :: cand_s(:)
integer, intent(in) :: n_cand, i_frag
character(len=len_fname) :: pdb_file
integer, parameter :: pdb_unit = 74
integer :: i_cand
character(len=3) :: num_frag

call num2str_999(i_frag, num_frag)
write(pdb_file, '(A,A,A)') 'fragment_', num_frag, '.pdb'
call open_write_pdb(pdb_unit, pdb_file)

do i_cand = 1, n_cand
    write(pdb_unit, 72) 'HETATM', i_cand, ' C  ', 'UNK', & 
        ' ', i_cand, ' ', cand_s(i_cand)%R(:), 1.00, cand_s(i_cand)%E 
end do

call close_write_pdb(pdb_unit)

72  format(a6,i5,1x,a4,1x,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2)
!-------------------------------------------------------------------------------
end subroutine write_FFT_result
!-------------------------------------------------------------------------------
subroutine write_pharm_result(hotspot_rec, i_frag)
!-------------------------------------------------------------------------------
type(frag_hotspot_type), intent(in) :: hotspot_rec
integer, intent(in) :: i_frag
character(len=len_fname) :: pdb_file
integer, parameter :: pdb_unit = 74
integer :: i_cand
character(len=3) :: num_frag

call num2str_999(i_frag, num_frag)
write(pdb_file, '(A,A,A)') 'pharm_', num_frag, '.pdb'
call open_write_pdb(pdb_unit, pdb_file)

do i_cand = 1, hotspot_rec%n_spot
    write(pdb_unit, 72) 'HETATM', i_cand, ' C  ', 'UNK', & 
        ' ', i_cand, ' ', hotspot_rec%R(:, i_cand), 1.00, hotspot_rec%E(i_cand)
end do

call close_write_pdb(pdb_unit)

72  format(a6,i5,1x,a4,1x,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2)
!-------------------------------------------------------------------------------
end subroutine write_pharm_result
!-------------------------------------------------------------------------------
subroutine generate_rec_graph(hotspot_rec, rec_graph, n_frag)
!-------------------------------------------------------------------------------
type(frag_hotspot_type), intent(in) :: hotspot_rec(:)
type(hotspot_graph_type), intent(out) :: rec_graph
integer, intent(in) :: n_frag
!
integer :: n_mem, i_frag, j_frag, i_mem, j_mem, i_idx, j_idx
real(dp) :: dr(3), dist

n_mem = 0
do i_frag = 1, n_frag
    do i_mem = 1, hotspot_rec(i_frag)%n_spot
        n_mem = n_mem + 1
        rec_graph%i_mem(n_mem) = i_mem
        rec_graph%hotspot_type(n_mem) = i_frag
        rec_graph%R(:,n_mem) = hotspot_rec(i_frag)%R(:,i_mem)
        !
        rec_graph%err(n_mem) = hotspot_rec(i_frag)%err(i_mem)
    end do
end do
rec_graph%n_mem = n_mem

rec_graph%dist(:,:) = 0.0
do i_mem = 1, n_mem - 1
    i_frag = rec_graph%hotspot_type(i_mem)
    i_idx = rec_graph%i_mem(i_mem)
    do j_mem = i_mem + 1, n_mem
        j_frag = rec_graph%hotspot_type(j_mem)
        j_idx = rec_graph%i_mem(j_mem)
        !
        dr(:) = hotspot_rec(i_frag)%R(:,i_idx) - hotspot_rec(j_frag)%R(:,j_idx)
        dist = sqrt(dot_product(dr,dr))
        !
        rec_graph%dist(i_mem, j_mem) = dist
        rec_graph%dist(j_mem, i_mem) = dist
    end do
end do

end subroutine generate_rec_graph
!-------------------------------------------------------------------------------
subroutine generate_lig_graph(lig_res, lig_graph)
!-------------------------------------------------------------------------------
type(ligand_residue_type), intent(in) :: lig_res
type(hotspot_graph_type), intent(out) :: lig_graph
!
integer :: n_mem, i_frag, j_frag, i_mem, j_mem
real(dp) :: dr(3), dist, R1(3), R2(3), vec1(3)

n_mem = 0
do i_frag = 1, n_frag
    n_mem = n_mem + 1
    lig_graph%i_mem(n_mem) = n_mem
    lig_graph%hotspot_type(n_mem) = i_frag
end do
lig_graph%n_mem = n_mem
lig_graph%dist(:,:) = 0.0

do i_mem = 1, n_mem
    i_frag = lig_graph%hotspot_type(i_mem)
    call get_frag_cntr(lig_res, fragment(:,i_frag), n_atm_frag(i_frag), R1,&
                       vec1)
    lig_graph%R(:,i_mem) = R1(:)
end do
    
do i_mem = 1, n_mem - 1
    R1(:) = lig_graph%R(:,i_mem)
    do j_mem = i_mem + 1, n_mem
        R2(:) = lig_graph%R(:,j_mem)
        !
        dr(:) = R2(:) - R1(:)
        dist = sqrt(dot_product(dr,dr))
        !
        lig_graph%dist(i_mem, j_mem) = dist
        lig_graph%dist(j_mem, i_mem) = dist
    end do
end do
!

end subroutine generate_lig_graph
!-------------------------------------------------------------------------------
subroutine do_pharmdock(protein, ligand, rec_graph, lig_graph, trans, &
                        gene_s, n_gene, E_cutoff, max_gene)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
type(hotspot_graph_type), intent(in) :: rec_graph, lig_graph
real(dp), intent(in) :: trans(:), E_cutoff
real(dp), intent(inout) :: gene_s(:,:)
integer, intent(out) :: n_gene
integer, intent(in) :: max_gene
!
type(assoc_graph_type) :: assoc_graph
type(clique_type) :: clique(max_clique)
integer :: n_clique, n_mem, gene_dim
integer :: i_cl, i_mem, lig_idx, rec_idx, i
integer :: key(max_gene), n_gene0
real(dp) :: tmp_gene(ligand%n_br+6), gene0_s(ligand%n_br+6, max_gene), E_s(max_gene)
real(dp) :: R1(3,max_lig_frag), R2(3,max_lig_frag), ref_R2(3,max_lig_frag), add
real(dp) :: lig_vec(3), err
real(dp) :: cntr1(3), cntr2(3), vec1(3), vec2(3)
real(dp) :: U(3,3), error, theta, axis(3)
logical :: wrong_orientation
!
real(dp) :: min_gene(ligand%n_br+6)
type(molecule_type) :: tmp_prot
type(energy_type) :: ff, min_ff
real(dp) :: g(3,ligand%n_atm)
logical :: matched_core(0:ligand%n_core_br)
integer :: n_matched_core
integer :: i_iter
integer, parameter :: max_iter = 20

call generate_assoc_graph(rec_graph, lig_graph, assoc_graph)
call find_clique(assoc_graph, clique, n_clique)
!print*, 'N_CLIQUE: ', n_clique
!
n_gene = 0
tmp_prot = protein
gene_dim = ligand%n_br + 5
!
do i_cl = 1, n_clique
    !if (clique(i_cl)%n_mem < min(3, clique(1)%n_mem)) exit
    !matched_core(:) = .false.
    do i_mem = 1, clique(i_cl)%n_mem
        rec_idx = assoc_graph%i_mem(1,clique(i_cl)%i_mem(i_mem))
        lig_idx = assoc_graph%i_mem(2,clique(i_cl)%i_mem(i_mem))
        !
        R1(:,i_mem) = rec_graph%R(:,rec_idx)
        !
        R2(:,i_mem) = lig_graph%R(:,lig_idx)
    end do
    !
    !n_matched_core = 0
    !do i_mem = 0, ligand%n_core_br
    !    if (matched_core(i_mem)) then
    !        n_matched_core = n_matched_core + 1
    !    end if
    !end do
    !if (n_matched_core < min(2, ligand%n_core_br)) cycle
    !
    n_mem = clique(i_cl)%n_mem
    if (clique(i_cl)%n_mem == 2) then
        do i = 1, 3
            vec1(i) = 2.0d0*random()-1.0d0
            vec2(i) = 2.0d0*random()-1.0d0
        end do
        call v_norm(vec1)
        call v_norm(vec2)
        do i = 1, 3
            R1(i,3) = R1(i,2) + vec1(i)
            R2(i,3) = R2(i,2) + vec1(i)
        end do
        do i = 1, 3
            R1(i,4) = R1(i,1) + vec2(i)
            R2(i,4) = R2(i,1) + vec2(i)
        end do
        n_mem = n_mem + 2
    end if
    ref_R2 = R2
    call ls_rmsd(n_mem, R2, R1, 1, U, error)
    !
    do i_mem = 1, 3
        cntr1(i_mem) = sum(R1(i_mem,1:n_mem))/dble(n_mem)
        cntr2(i_mem) = sum(ref_R2(i_mem,1:n_mem))/dble(n_mem)
    end do
    n_gene = n_gene + 1
    gene_s(1:3, n_gene) = cntr1(1:3) - cntr2(1:3) + trans(1:3)
    call rotmat_to_axisangle(U, theta, axis)
    call v_norm(axis)
    !gene_s(4, n_gene) = theta
    gene_s(4:6, n_gene) = axis(1:3)*theta
    !
    call construct_ligand(tmp_prot, ligand, gene_s(:,n_gene))
    call update_R_for_ligand(tmp_prot%ligand(ligand%lig_no), ligand%lig_no)
    call ligdock_energy_using_grid(tmp_prot, ligand, min_ff, g, .false.)
    if (min_ff%total > E_cutoff) then
        n_gene = n_gene - 1
        cycle
    end if
    !
    min_gene(1:gene_dim) = gene_s(1:gene_dim, n_gene)
    do i_iter = 1, max_iter
        gene_s(1:gene_dim,n_gene) = min_gene(1:gene_dim)
        call perturb_gene(gene_s(:,n_gene))
        call construct_ligand(tmp_prot, ligand, gene_s(:,n_gene))
        call update_R_for_ligand(tmp_prot%ligand(ligand%lig_no), ligand%lig_no)
        call ligdock_energy_using_grid(tmp_prot, ligand, ff, g, .false.)
        if (ff%total < min_ff%total) then
            min_gene(1:gene_dim) = gene_s(1:gene_dim,n_gene)
            min_ff = ff
        end if
    end do
    !
    if (n_gene == max_gene) exit
end do
!print*, "N_GENE: ", n_gene

!if (n_gene < max_gene) then
if (n_gene == 0) then
    i_mem = 0
    do i_cl = 1, assoc_graph%n_mem
        rec_idx = assoc_graph%i_mem(1,i_cl)
        lig_idx = assoc_graph%i_mem(2,i_cl)
        !if (rec_graph%hotspot_type(rec_idx) /= 1) cycle 
        !
        n_gene = n_gene + 1
        cntr1(:) = rec_graph%R(:,rec_idx)
        cntr2(:) = lig_graph%R(:,lig_idx)
        gene_s(1:3, n_gene) = cntr1(1:3) - cntr2(1:3) + trans(1:3)
        call gen_random_axis(axis)
        !
        theta = two_pi*random() - 1.0d0*pi
        gene_s(4:6, n_gene) = axis(1:3)*bound_ang(theta)
        !
        if (n_gene == max_gene) exit
    end do
end if

end subroutine do_pharmdock
!-------------------------------------------------------------------------------
subroutine perturb_gene(gene)
!-------------------------------------------------------------------------------
real(dp), intent(inout) :: gene(:)
real(dp), parameter :: l_trans = 0.5d0
real(dp), parameter :: l_axis  = 0.2d0
real(dp), parameter :: l_angle = 10.0d0*deg2rad
integer :: i
real(dp) :: p(3), q(4)

do i = 1, 3
    gene(i) = gene(i) + (2.0d0*random()-1.0d0)*l_trans
end do

p = gene(4:6)
call expmap_quat(p,q)
call pert_rotation(l_angle, q)
call inv_expmap_quat(q,p)
gene(4:6) = p

end subroutine perturb_gene
!-------------------------------------------------------------------------------
subroutine setup_FFT_grid_info(FFT_grid_info)
!-------------------------------------------------------------------------------
type(docking_grid_param_type), intent(out) :: FFT_grid_info
integer, dimension(3) :: added_elem
integer :: i, power_2

do i = 1, 3
    power_2 = 2
    do while (dock_grid_info%n_elem(i) > power_2)
        power_2 = power_2 * 2
    end do
    added_elem(i) = power_2 - dock_grid_info%n_elem(i)
end do

! 1~dock_grid_info%n_elem(i) : same to docking grid
! dock_grid_info%n_elem(i)+1 ~ FFT_grid_info%n_elem(i) : Fill max_energy
FFT_grid_info%n_elem(:) = dock_grid_info%n_elem(:) + added_elem(:)
FFT_grid_info%grid_cntr(:) = dock_grid_info%grid_cntr(:)
FFT_grid_info%grid_width = dock_grid_info%grid_width

end subroutine setup_FFT_grid_info
!-------------------------------------------------------------------------------
subroutine get_r_conv(r_conv, n_elem)
!-------------------------------------------------------------------------------
type(docking_FFT_grid_type), intent(inout) :: r_conv
integer, intent(in) :: n_elem(3)
integer :: i_conv

if (use_atdk3 .or. use_atdk4) then
    do i_conv = 1, r_conv%n_atdk
        FFT_status = DftiComputeForward(desc, r_conv%atdk_grid(:,i_conv))
    end do
end if

if (use_drugscore) then
    do i_conv = 1, r_conv%n_drugscore
        FFT_status = DftiComputeForward(desc, r_conv%drugscore_grid(:,i_conv))
    end do
end if
if (use_Xscore) then
    FFT_status = DftiComputeForward(desc, r_conv%HM_grid)
end if

end subroutine get_r_conv
!-------------------------------------------------------------------------------
subroutine update_R_frag(fragment, n_atm, ligand, R_frag, to_center)
!-------------------------------------------------------------------------------
integer, intent(in) :: fragment(:), n_atm
type(ligand_type), intent(in) :: ligand
real(dp), intent(out) :: R_frag(:,:)
logical, intent(in), optional :: to_center
integer :: i_atm, atm_idx
real(dp) :: trans_v(3), cntr(3)
logical :: move_to_center

if (present(to_center)) then
    move_to_center = to_center
else
    move_to_center = .true.
end if

R_frag = 0.0
cntr(:) = 0.0
do i_atm = 1, n_atm
    atm_idx = fragment(i_atm)
    R_frag(:,i_atm) = ref_lig(ligand%lig_type)%R(:,atm_idx)
    cntr(:) = cntr(:) + R_frag(:,i_atm)
end do

if (.not. move_to_center) return

cntr(:) = cntr(:)/n_atm

trans_v(:) = cntr(:) - dock_grid_info%grid_cntr(:)

do i_atm = 1, n_atm
    R_frag(:,i_atm) = R_frag(:,i_atm) - trans_v(:)
end do

end subroutine update_R_frag
!-------------------------------------------------------------------------------
subroutine rotate_fragment(R_frag, R_frag_rotated, n_atm, U)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: R_frag(:,:)
real(dp), intent(out) :: R_frag_rotated(:,:)
integer, intent(in) :: n_atm
real(dp), intent(in) :: U(3,3)
real(dp) :: cntr(3), trans_v(3)
integer :: i_atm

cntr(:) = 0.0
do i_atm = 1, n_atm
    cntr(:) = cntr(:) + R_frag(:,i_atm)
end do
cntr(:) = cntr(:)/dble(n_atm)

trans_v = dock_grid_info%grid_cntr(:) - cntr(:)

R_frag_rotated = R_frag
call transrot_crd(R_frag_rotated(:,1:n_atm), n_atm, trans_v, cntr, U)

end subroutine rotate_fragment
!-------------------------------------------------------------------------------
subroutine translate_fragment(r_conv, n_elem, cand_s, n_pool)
!-------------------------------------------------------------------------------
type(docking_FFT_grid_type), intent(in) :: r_conv
integer, intent(in) :: n_elem(3)
integer, intent(in) :: n_pool
type(FFT_result_type), intent(out) :: cand_s(:)
!
integer :: i_conv, n_elem_3
real :: score(n_elem(1)*n_elem(2)*n_elem(3))

score = 0.0
n_elem_3 = n_elem(1)*n_elem(2)*n_elem(3)

if (use_atdk3 .or. use_atdk4) then
    do i_conv = 1, l_grid%n_atdk
        FFT_status = DftiComputeForward(desc, l_grid%atdk_grid(:,i_conv))

        l_grid%atdk_grid(:,i_conv) &
            = r_conv%atdk_grid(:,i_conv)*conjg(l_grid%atdk_grid(:,i_conv))
        !
        FFT_status = DftiComputeBackward(desc, l_grid%atdk_grid(:,i_conv))
        score = score &
              + real(l_grid%atdk_grid(:,i_conv))/n_elem_3
    end do
end if
if (use_drugscore) then
    do i_conv = 1, l_grid%n_drugscore
        FFT_status = DftiComputeForward(desc, l_grid%drugscore_grid(:,i_conv))

        l_grid%drugscore_grid(:,i_conv) &
            = conjg(r_conv%drugscore_grid(:,i_conv))*l_grid%drugscore_grid(:,i_conv)

        FFT_status = DftiComputeBackward(desc, l_grid%drugscore_grid(:,i_conv))
        score = score &
              + drugscore_w*real(l_grid%drugscore_grid(:,i_conv))/n_elem_3 
    end do
end if
if (use_Xscore) then
    FFT_status = DftiComputeForward(desc, l_grid%HM_grid(:))

    l_grid%HM_grid(:) &
        = conjg(r_conv%HM_grid(:))*l_grid%HM_grid(:)
    
    FFT_status = DftiComputeBackward(desc, l_grid%HM_grid(:))
    score = score &
          + X_HM_w*real(l_grid%HM_grid(:))/n_elem_3 
end if
    
! (1,1,1) = (cntr, cntr, cntr)
! (i,j,k) = (cntr+(i-1), cntr+(j-1), cntr+(k-1))

call get_possible_candidate(score, n_elem, cand_s, n_pool)

end subroutine translate_fragment
!-------------------------------------------------------------------------------
subroutine get_possible_candidate(score, n_elem, cand_s, n_pool)
!-------------------------------------------------------------------------------
real, intent(inout) :: score(:)
integer, intent(in) :: n_elem(3), n_pool
type(FFT_result_type), intent(out) :: cand_s(:)
!
integer :: add_elem(3)
integer :: diff(3), min_idx(3), min_hash(1), i, i_pool

add_elem(:) = n_elem(:) - dock_grid_info%n_elem(:)
!
do i_pool = 1, n_pool
    min_hash = minloc(score)
    !
    cand_s(i_pool)%E = score(min_hash(1))
    !
    min_idx(1) = mod(min_hash(1), n_elem(1))
    min_idx(2) = mod((min_hash(1)-min_idx(1))/n_elem(1), n_elem(2)) + 1
    min_idx(3) = min_hash(1) - n_elem(1)*(min_idx(2)-1) - min_idx(1)
    min_idx(3) = min_idx(3)/(n_elem(1)*n_elem(2)) + 1
    !
    diff(:) = min_idx(:) - 1
    do i = 1, 3
        if (diff(i) > (dock_grid_info%n_elem(i)-1)/2 + add_elem(i)) then
            diff(i) = diff(i) - n_elem(i)
        end if
    end do
    cand_s(i_pool)%R(:) = dock_grid_info%grid_cntr(:) + diff(:)*dock_grid_info%grid_width
    !
    score(min_hash(1)) = max_energy
end do

end subroutine get_possible_candidate
!-------------------------------------------------------------------------------
subroutine update_hotspot_rec(hotspot_rec, cand_s, n_cand, n_atm_frag, i_frag)
!-------------------------------------------------------------------------------
type(frag_hotspot_type), intent(inout) :: hotspot_rec
type(FFT_result_type), intent(in):: cand_s(:)
integer, intent(in) :: n_cand, n_atm_frag
integer, intent(in) :: i_frag ! for DEBUG
!
real(dp) :: E_s(n_cand)
real(dp) :: dist, d(3), error
integer :: i_cand, j_cand, key(n_cand), n_cand_fin, i_cl, n_cl, n_cl_2, n_count
integer :: minPts
!
real(dp), allocatable :: d_mat(:,:), E_tmp(:)
real, allocatable :: R_s(:,:), R_tmp(:,:)
integer, allocatable :: cl_id(:,:), cl_size(:), k_cl(:)
integer, allocatable :: cl_id_2(:,:), cl_size_2(:), k_cl_2(:)
logical :: is_far
real(dp) :: Ecut

if (n_cand == 0) then
    hotspot_rec%n_spot = 0
    return
end if

do i_cand = 1, n_cand
    E_s(i_cand) = cand_s(i_cand)%E
end do

call sort2(n_cand, E_s(:), key(:))
Ecut = min(FFT_E_CUTOFF, E_s(1)*0.75d0)
!
n_cand_fin = n_cand
do i_cand = 1, n_cand
    if (E_s(i_cand) > Ecut) then
        n_cand_fin = i_cand - 1
        exit
    end if
end do
n_cand_fin = min(n_cand_fin, n_cand/10)
if (n_cand_fin == 0) then
    hotspot_rec%n_spot = 0
    return
end if

!n_cand_fin = n_cand_fin * 3 / 4
allocate(R_s(3,n_cand_fin))
do i_cand = 1, n_cand_fin
    R_s(:,i_cand) = cand_s(key(i_cand))%R(:)
end do


allocate(d_mat(n_cand_fin, n_cand_fin))
allocate(cl_id(n_cand_fin, n_cand_fin))
allocate(cl_size(n_cand_fin))
allocate(k_cl(n_cand_fin))
!
!call hierarchical_clustering_dcut(d_mat, n_cand_fin, cl_cut, &
!                             E_s(1:n_cand_fin), cl_id, cl_size, n_cl)
!call DBSCAN(d_mat, n_cand_fin, cl_cut, 6, cl_id, cl_size, n_cl, .true.)
!call DBSCAN(d_mat, n_cand_fin, cl_cut, 200, cl_id, cl_size, n_cl, .true.)
if (n_cand_fin == 1) then
    n_cl = 1
    cl_id(1,1) = 1
    cl_size(1) = 1
else
    d_mat = 0.0d0
    do i_cand = 1, n_cand_fin - 1
        do j_cand = i_cand + 1, n_cand_fin
            d(:) = R_s(:,i_cand) - R_s(:,j_cand)
            dist = sqrt(dot_product(d,d))
            d_mat(i_cand, j_cand) = dist
            d_mat(j_cand, i_cand) = dist
        end do
    end do
    !
    minPts = min(50, n_cand_fin/3)
    minPts = max(1, minPts)
    !
    call DBSCAN_2(d_mat, n_cand_fin, cl_cut, minPts, cl_id, cl_size, n_cl, .true.)
    !call DBSCAN(d_mat, n_cand_fin, cl_cut, minPts, cl_id, cl_size, n_cl, .true.)
    !print*, n_cand_fin, n_cl
end if
!
if (n_cl < max_hotspot/4) then
    do i_cl = 1, n_cl 
        !call write_cluster_result(R_s, E_s, cl_id(:,i_cl), cl_size(i_cl), i_cl, i_frag) ! DEBUG
        i_cand = cl_id(1,i_cl)
        call get_cluster_center(cl_id(1:cl_size(i_cl),i_cl), &
                                cl_size(i_cl), R_s, hotspot_rec%R(:,i_cl))
        !
        hotspot_rec%E(i_cl) = E_s(i_cand)
    end do
    hotspot_rec%n_spot = n_cl
else
    deallocate(d_mat)
    allocate(d_mat(n_cl, n_cl))
    allocate(R_tmp(3,n_cl))
    allocate(E_tmp(n_cl))
    d_mat = 0.0d0
    do i_cand = 1, n_cl
        call get_cluster_center(cl_id(1:cl_size(i_cand),i_cand), &
                                cl_size(i_cand), R_s, R_tmp(:,i_cand))
        E_tmp(i_cand) = E_s(cl_id(1,i_cand))
        do j_cand = 1, i_cand-1
            d(:) = R_tmp(:,i_cand) - R_tmp(:,j_cand)
            dist = sqrt(dot_product(d,d))
            d_mat(i_cand, j_cand) = dist
            d_mat(j_cand, i_cand) = dist
        end do
    end do
    !
    deallocate(cl_id)
    deallocate(cl_size)
    deallocate(k_cl)
    allocate(cl_id(n_cl, n_cl))
    allocate(cl_size(n_cl))
    allocate(k_cl(n_cl))
    !
    print*, 're-clustering from ', n_cl
    n_cand_fin = n_cl
    call hierarchical_clustering_dcut(d_mat, n_cand_fin, cl_cut*2.0d0, &
                                      E_tmp(1:n_cand_fin), cl_id, cl_size, n_cl)
    do i_cl = 1, min(n_cl, max_hotspot)
        !call write_cluster_result(R_s, E_s, cl_id(:,i_cl), cl_size(i_cl), i_cl, i_frag) ! DEBUG
        i_cand = cl_id(1,i_cl)
        call get_cluster_center(cl_id(1:cl_size(i_cl),i_cl), &
                                cl_size(i_cl), R_tmp, hotspot_rec%R(:,i_cl))
        !
        hotspot_rec%E(i_cl) = E_tmp(i_cand)
    end do
    hotspot_rec%n_spot = min(n_cl, max_hotspot)
    print*, 're-clustering ended ', n_cl
    deallocate(R_tmp)
    deallocate(E_tmp)
end if

deallocate(cl_id)
deallocate(cl_size)
deallocate(k_cl)
deallocate(R_s)

end subroutine update_hotspot_rec
!-------------------------------------------------------------------------------
subroutine write_cluster_result(R_s, E_s, id_s, n_mem, i_cl, i_frag)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: R_s(:,:), E_s(:)
integer, intent(in) :: id_s(:), n_mem, i_cl, i_frag
character(len=len_fname) :: pdb_file
integer, parameter :: pdb_unit = 74
integer :: i_cand
character(len=3) :: num_cl, num_frag

call num2str_999(i_cl, num_cl)
call num2str_999(i_frag, num_frag)
write(pdb_file, '(A,A,A,A,A)') 'cluster_', num_frag, '_', num_cl, '.pdb'
call open_write_pdb(pdb_unit, pdb_file)

do i_cand = 1, n_mem
    write(pdb_unit, 72) 'HETATM', i_cand, ' C  ', 'UNK', & 
        ' ', i_cand, ' ', R_s(:,id_s(i_cand)), 1.00, E_s(id_s(i_cand))
end do

call close_write_pdb(pdb_unit)

72  format(a6,i5,1x,a4,1x,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2)

end subroutine write_cluster_result
!-------------------------------------------------------------------------------
subroutine get_cluster_center(id_s, n_mem, R_s, cntr)
!-------------------------------------------------------------------------------
integer, intent(in) :: id_s(:), n_mem
real, intent(in) :: R_s(:,:)
real, intent(out) :: cntr(:)
!
integer :: i_mem

cntr = 0.0d0
do i_mem = 1, n_mem
    cntr(:) = cntr(:) + R_s(:,id_s(i_mem))
end do
cntr(:) = cntr(:)/n_mem

end subroutine get_cluster_center
!-------------------------------------------------------------------------------
subroutine get_frag_cntr(lig_res, fragment, n_atm_frag, cntr, vec1)
!-------------------------------------------------------------------------------
type(ligand_residue_type), intent(in) :: lig_res
integer, intent(in) :: fragment(:), n_atm_frag
real(dp), intent(out) :: cntr(:), vec1(:)
integer :: i_atm
real(dp) :: vec_1(3), vec_2(3)

cntr(:) = 0.0

do i_atm = 1, n_atm_frag
    cntr(:) = cntr(:) + lig_res%R(:,fragment(i_atm))
end do
cntr(:) = cntr(:) / n_atm_frag
vec_1(:) = lig_res%R(:,fragment(2)) - lig_res%R(:,fragment(1))
vec_2(:) = lig_res%R(:,fragment(3)) - lig_res%R(:,fragment(1))

call cross(vec_1, vec_2, vec1)
call v_norm(vec1)

end subroutine get_frag_cntr
!-------------------------------------------------------------------------------
subroutine generate_assoc_graph(graph1, graph2, assoc_graph)
!-------------------------------------------------------------------------------
type(hotspot_graph_type), intent(in) :: graph1, graph2
type(assoc_graph_type), intent(out) :: assoc_graph
integer :: n_mem, i_mem, j_mem, idx(4)
logical :: compatible

n_mem = 0
do i_mem = 1, graph1%n_mem
    do j_mem = 1, graph2%n_mem
        if (graph1%hotspot_type(i_mem) == graph2%hotspot_type(j_mem)) then
            n_mem = n_mem + 1
            assoc_graph%i_mem(1,n_mem) = i_mem
            assoc_graph%i_mem(2,n_mem) = j_mem
        end if
    end do
end do
assoc_graph%n_mem = n_mem

assoc_graph%connect(:,:) = .false.
do i_mem = 1, assoc_graph%n_mem-1
    do j_mem = i_mem + 1, assoc_graph%n_mem
        idx(1) = assoc_graph%i_mem(1,i_mem)
        idx(2) = assoc_graph%i_mem(2,i_mem)
        idx(3) = assoc_graph%i_mem(1,j_mem)
        idx(4) = assoc_graph%i_mem(2,j_mem)
        call check_compatible_edge(idx, graph1, graph2, compatible)
        if (compatible) then
            assoc_graph%connect(i_mem,j_mem) = .true.
            assoc_graph%connect(j_mem,i_mem) = .true.
        end if
    end do
end do 

!! TEST
!write(*,*) 'TEST'
!do i_mem = 1, assoc_graph%n_mem - 1
!    do j_mem = i_mem + 1, assoc_graph%n_mem
!        if (assoc_graph%connect(i_mem, j_mem)) then
!            write(*,'(I4,I4,I4,I4,I4,I4)') i_mem, j_mem, assoc_graph%i_mem(:,i_mem), assoc_graph%i_mem(:,j_mem)
!        end if
!    end do
!end do

end subroutine generate_assoc_graph
!-------------------------------------------------------------------------------
subroutine check_compatible_edge(idx, graph1, graph2, compatible)
!-------------------------------------------------------------------------------
integer, dimension(4), intent(in):: idx
type(hotspot_graph_type), intent(in) :: graph1, graph2
logical, intent(out) :: compatible
real :: dist1, dist2

compatible = .false.
if (idx(1) == idx(3)) return
if (idx(2) == idx(4)) return

dist1 = graph1%dist(idx(1),idx(3))
dist2 = graph2%dist(idx(2),idx(4))
if (abs(dist1-dist2) < d_cut) then
    compatible = .true.
end if

end subroutine check_compatible_edge
!-------------------------------------------------------------------------------
subroutine find_clique(assoc_graph, clique, n_clique)
!-------------------------------------------------------------------------------
type(assoc_graph_type), intent(in) :: assoc_graph
type(clique_type), intent(out) :: clique(:)
integer, intent(out) :: n_clique
integer :: n_exclude(0:max_lig_frag), n_cand(0:max_lig_frag)
integer :: cand(max_mem, 0:max_lig_frag), exclude(max_mem, 0:max_lig_frag)
integer :: curr_clique(max_lig_frag)
integer :: i_mem, depth

! initialize
n_clique = 0
do i_mem = 1, max_clique
    clique(i_mem)%n_mem = 0
end do

n_exclude(:) = 0
exclude(:,:) = -1
curr_clique(:) = -1
!
n_cand(:) = 0
cand(:,:) = -1
!
depth = 0
do i_mem = 1, assoc_graph%n_mem
    cand(i_mem,depth) = i_mem
end do
n_cand(depth) = assoc_graph%n_mem

call bron_kerbosch(assoc_graph, cand, n_cand, exclude, n_exclude, curr_clique, &
                   clique, n_clique, depth)

end subroutine find_clique
!-------------------------------------------------------------------------------
recursive subroutine bron_kerbosch(assoc_graph, cand, n_cand, exclude, n_ex, &
                                   curr_clique, clique, n_clique, depth)
!-------------------------------------------------------------------------------
type(assoc_graph_type), intent(in) :: assoc_graph
integer, intent(inout) :: cand(max_mem,0:max_lig_frag)
integer, intent(inout) :: exclude(max_mem,0:max_lig_frag), curr_clique(:)
integer, intent(inout) :: n_cand(0:max_lig_frag), n_ex(0:max_lig_frag), depth
type(clique_type), intent(inout) :: clique(:)
integer, intent(inout) :: n_clique
integer :: i_cand, n_cand0
logical :: dead_end, excluded

call check_dead_end(assoc_graph, exclude(:,depth), n_ex(depth), &
                    cand(:,depth), n_cand(depth), dead_end)
if (dead_end) return 

depth = depth + 1

do i_cand = 1, n_cand(depth-1)
    curr_clique(depth) = cand(i_cand, depth-1)
    !
    call is_in_exclude(curr_clique(depth), exclude(:,depth-1), n_ex(depth-1), excluded)
    if (excluded) cycle
    !
    ! udpate candidate member (intersection of candidate and neighbor of selected node)
    call update_member(curr_clique(depth), cand, n_cand, assoc_graph, depth)
    ! udpate exclude member (intersection of exclude and neighbor of selected node)
    call update_member(curr_clique(depth), exclude, n_ex, assoc_graph, depth)
    !
    if (n_cand(depth) == 0 .and. n_ex(depth) == 0 .and. depth > 1) then ! Find clique
        if (n_clique == 0) then
            n_clique = n_clique + 1
            clique(n_clique)%n_mem = depth
            clique(n_clique)%i_mem(1:depth) = curr_clique(1:depth)
        else if (depth > clique(max_clique)%n_mem) then
            call update_clique(curr_clique, depth, clique, n_clique)
        end if
    else
        call bron_kerbosch(assoc_graph, cand, n_cand, exclude, n_ex, &
                           curr_clique, clique, n_clique, depth)
    end if
    !
    n_ex(depth-1) = n_ex(depth-1) + 1
    exclude(n_ex(depth-1), depth-1) = curr_clique(depth)
end do

depth = depth - 1

end subroutine bron_kerbosch
!-------------------------------------------------------------------------------
subroutine update_member(mem_idx, cand, n_cand, assoc_graph, depth)
!-------------------------------------------------------------------------------
integer, intent(in) :: mem_idx, depth
integer, intent(inout) :: cand(max_mem,0:max_lig_frag), n_cand(0:max_lig_frag)
type(assoc_graph_type), intent(in) :: assoc_graph
!
integer :: i_mem, mem_jdx

n_cand(depth) = 0
do i_mem = 1, n_cand(depth-1)
    mem_jdx = cand(i_mem, depth-1)
    if (assoc_graph%connect(mem_jdx, mem_idx)) then
        n_cand(depth) = n_cand(depth) + 1
        cand(n_cand(depth), depth) = mem_jdx
    end if
end do

end subroutine update_member
!-------------------------------------------------------------------------------
subroutine update_clique(curr_clique, clique_size, clique, n_clique)
!-------------------------------------------------------------------------------
integer, intent(in) :: curr_clique(:), clique_size
type(clique_type), intent(inout) :: clique(:)
integer, intent(inout) :: n_clique
type(clique_type) :: clique0(max_clique)
integer :: i_cl, j_cl

clique0 = clique
do i_cl = 1, max_clique
    if (clique_size > clique(i_cl)%n_mem) then 
        clique(i_cl)%i_mem(:) = curr_clique(:)
        clique(i_cl)%n_mem = clique_size
        do j_cl = i_cl+1, min(n_clique+1, max_clique)
            clique(j_cl) = clique0(j_cl-1)
        end do
        exit
    end if
end do
if (n_clique < max_clique) n_clique = n_clique + 1

end subroutine update_clique
!-------------------------------------------------------------------------------
subroutine check_dead_end(assoc_graph, exclude, n_ex, cand, n_cand, dead_end)
!-------------------------------------------------------------------------------
type(assoc_graph_type), intent(in) :: assoc_graph
integer, intent(in) :: exclude(:), n_ex, cand(:), n_cand
logical, intent(out)  :: dead_end
integer :: i, j, v, u

do i = 1, n_ex
    v = exclude(i)
    dead_end = .true.
    do j = 1, n_cand
        u = cand(j)
        if (.not. assoc_graph%connect(u,v)) then
            dead_end = .false.
            exit
        end if
    end do
    if (dead_end) return
end do

dead_end = .false.

end subroutine check_dead_end
!-------------------------------------------------------------------------------
subroutine is_in_exclude(mem, exclude, n_ex, excluded)
!-------------------------------------------------------------------------------
integer, intent(in) :: mem
integer, intent(in) :: exclude(:), n_ex
logical, intent(out) :: excluded
integer :: i_ex

excluded = .false.
do i_ex = 1, n_ex
    if (exclude(i_ex) == mem) then
        excluded = .true.
        return
    end if
end do

end subroutine is_in_exclude
!-------------------------------------------------------------------------------
END MODULE LIGDOCK_CSA_FRAGFFT
!-------------------------------------------------------------------------------
