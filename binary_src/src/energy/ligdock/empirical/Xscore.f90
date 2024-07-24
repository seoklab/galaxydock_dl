!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/ligdock/empirical/Xscore.f90
!
! Description:
!
!-------------------------------------------------------------------------------
MODULE XSCORE
!-------------------------------------------------------------------------------

use globals
use logger
use string, only: parse_string
use mathfunctions, only: sigmoidal1, sigmoidal2
!
use energy_vars
use ligdock_E_utils, only: copy_to_grid, calc_E_using_grid
use xlogp_Score_m, only: get_logp_of_ligand

implicit none
save 
private

integer, parameter :: max_res_info = 40
real(dp), parameter :: HMcutoff = 6.0d0
real(dp), parameter :: HMcutoff2 = HMcutoff**2

!-------------------------------------------------------------------------------
type Xscore_res_info_type
!-------------------------------------------------------------------------------
character(len=4) :: res_name
character(len=4), dimension(max_atm) :: atom_name
real(dp), dimension(max_atm) :: logP
integer :: n_atm
!-------------------------------------------------------------------------------
end type Xscore_res_info_type
!-------------------------------------------------------------------------------

real(dp), allocatable :: logP(:)
logical, allocatable :: is_ligand(:)

public :: initialize_X_score
public :: finalize_X_score
public :: construct_HM_grid
public :: calc_HM_score_using_grid
public :: calc_HM_score
public :: FFT_rec_grid_HM
public :: FFT_lig_grid_HM
public :: get_hydrophobicity

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_X_score(protein, ligand)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
integer :: i_lig, ref_lig_no, lig_atm, n_atm

allocate(logP(tn%atom), is_ligand(tn%atom))

logP(:) = 0.0d0
is_ligand(:) = .false.

! get log_P value
call get_logp_of_protein(protein, logP)

lig_atm = tn%nonligatm
do i_lig = 1, protein%n_lig
    ref_lig_no = protein%ligand(i_lig)%lig_type
    n_atm = ref_lig(ref_lig_no)%n_atm

    if (ligand%lig_no == i_lig) then
        is_ligand(lig_atm+1:lig_atm+n_atm) = .true.
    end if

    call get_logp_of_ligand(ref_lig(ref_lig_no), logP(lig_atm+1:lig_atm+n_atm),&
                            .true.)
    lig_atm = lig_atm + n_atm
end do

end subroutine initialize_X_score
!-------------------------------------------------------------------------------
subroutine finalize_X_score()
!-------------------------------------------------------------------------------

deallocate(logP, is_ligand)

end subroutine finalize_X_score
!-------------------------------------------------------------------------------
subroutine get_logp_of_protein(protein, logP)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
real(dp), intent(inout) :: logP(:)
character(len=len_fname) :: file_name
type(Xscore_res_info_type) :: Xscore_res_info(max_res_info)
integer :: num_res_info

file_name = trim(data_dir) // 'X_score_res.prm'
call read_residue_def_file(file_name, Xscore_res_info, num_res_info)

call assign_logP_of_protein(protein, logP, Xscore_res_info, num_res_info)

!-------------------------------------------------------------------------------
end subroutine get_logp_of_protein
!-------------------------------------------------------------------------------
subroutine read_residue_def_file(file_name, Xscore_res_info, num_res_info)
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: file_name
type(Xscore_res_info_type), intent(out) :: Xscore_res_info(:)
integer, intent(out) :: num_res_info
!
integer :: f_unit, ioerror
integer :: num_word
character(len=len_fname) :: line, word(25)
! 
integer :: res_idx, atm_idx
character(len=4) :: res_name

f_unit = 17
open(unit=f_unit, file=trim(file_name), action='read', iostat=ioerror)
if (ioerror /= 0) then
    write(log_msg,'(A,A)') "ERROR: Cannot find data file to get logP value of protein, ",&
                            file_name
    call terminate_with_error(log_msg)
end if
  
res_idx = 0
atm_idx = 0

do
    read(f_unit, '(A120)', iostat=ioerror) line
    if(ioerror /= 0) exit
   
    if (line(1:1) == "!" .or. line(1:1) == "#") cycle

    call parse_string(line, num_word, word)
    if(word(1) == 'RESI') then 
        res_name = word(2)
        res_idx = res_idx + 1
        atm_idx = 0
        Xscore_res_info(res_idx)%res_name = res_name
    end if

    if(word(1) == 'ATOM') then
        atm_idx = atm_idx + 1
        Xscore_res_info(res_idx)%n_atm = atm_idx
        Xscore_res_info(res_idx)%atom_name(atm_idx) = word(2)
        read(word(9), '(f6.3)') Xscore_res_info(res_idx)%logP(atm_idx)
    end if
end do
num_res_info = res_idx

close(f_unit)

end subroutine read_residue_def_file
!-------------------------------------------------------------------------------
subroutine assign_logP_of_protein(protein, logP, Xscore_res_info, num_res_info)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
real(dp), intent(inout) :: logP(:)
type(Xscore_res_info_type), intent(in) :: Xscore_res_info(:)
integer, intent(in) :: num_res_info
integer :: i_atm, res_no, atm_no, ref_res_no
integer :: res_type, atm_type
logical :: is_ter
character(len=4) :: res_name, atom_name

! protein
do i_atm = 1, tn%stdatm
    res_no = i_R(1, i_atm)
    atm_no = i_R(2, i_atm)
    ref_res_no = protein%residue(res_no)%res_type

    is_ter = .false.
    if (protein%residue(res_no)%ter_type == 'N' .or. &
        protein%residue(res_no)%ter_type == 'C') then
        is_ter = .true.
    end if

    if (is_ter) then
        res_name = protein%residue(res_no)%res_name(2:4)
    else
        res_name = protein%residue(res_no)%res_name
    end if
    call find_Xscore_res_idx(res_name, res_type, Xscore_res_info, num_res_info)

    atom_name = ref_res(ref_res_no)%atom_name(atm_no)
    call find_Xscore_atm_idx(atom_name, res_type, atm_type, Xscore_res_info)
    logP(i_atm) = Xscore_res_info(res_type)%logP(atm_type)
end do

! hetmol (small cofactor)
do i_atm = tn%stdatm + 1, tn%nonligatm
    res_no = i_R(1, i_atm) - tn%stdres
    atm_no = i_R(2, i_atm)
    ref_res_no = protein%hetmol(res_no)%res_type

    res_name = protein%hetmol(res_no)%res_name
    
    call find_Xscore_res_idx(res_name, res_type, Xscore_res_info, num_res_info)

    atom_name = ref_res(ref_res_no)%atom_name(atm_no)
    call find_Xscore_atm_idx(atom_name, res_type, atm_type, Xscore_res_info)
    logP(i_atm) = Xscore_res_info(res_type)%logP(atm_type)
end do

end subroutine assign_logP_of_protein
!===============================================================================
! Grid Construction
!===============================================================================
subroutine construct_HM_grid(grid_info, dock_grid, flex_res)
!-------------------------------------------------------------------------------
type(docking_grid_param_type), intent(in) :: grid_info
type(docking_grid_type), intent(inout) :: dock_grid
logical, intent(in) :: flex_res(:)
integer :: n_hash
real(dp) :: tmp_grid_E
real(dp) :: ref_pt(3), grid_pt(3)
integer :: i, x, y, z

n_hash = grid_info%n_elem(1) * grid_info%n_elem(2) * grid_info%n_elem(3) 
! allocate HMscore grid
allocate(dock_grid%HM_grid(4, n_hash))

tmp_grid_E = 0.0d0

do i = 1, 3
    ref_pt(i) = grid_info%grid_cntr(i) &
              - (dble(grid_info%n_elem(i))-1.0d0)/2.0d0 * grid_info%grid_width
end do

do z = 1, grid_info%n_elem(3)
    grid_pt(3) = (dble(z)-1.0d0)*grid_info%grid_width + ref_pt(3)
    do y = 1, grid_info%n_elem(2)
        grid_pt(2) = (dble(y)-1.0d0)*grid_info%grid_width + ref_pt(2)
        do x = 1, grid_info%n_elem(1)
            grid_pt(1) = (dble(x)-1.0d0)*grid_info%grid_width + ref_pt(1)
            call calc_HM_grid(grid_pt(1:3), tmp_grid_E, flex_res)
            call copy_to_grid(tmp_grid_E, grid_info%n_elem(1:3),&
                              dock_grid%HM_grid(:,:), x, y, z, n_hash)
        end do
    end do
end do

end subroutine construct_HM_grid
!-------------------------------------------------------------------------------
subroutine calc_HM_grid(grid_pt, pt_E, flex_res) 
!-------------------------------------------------------------------------------
real(dp), intent(in) :: grid_pt(3)
real(dp), intent(out) :: pt_E
logical, intent(in) :: flex_res(:)
!
integer :: i_atm, res_no
real(dp) :: dr(3), dist_sqr

pt_E = 0.0d0

do i_atm = 1, tn%atom
    if (is_ligand(i_atm)) cycle
    res_no = i_R(1, i_atm)
    if (res_no <= tn%stdres) then
        if (flex_res(res_no)) cycle
    end if
    
    dr(1:3) = R(1:3, i_atm) - grid_pt(1:3)
    dist_sqr = dot_product(dr, dr)

    if (dist_sqr > HMcutoff2) cycle

    pt_E = pt_E + logP(i_atm)
end do

end subroutine calc_HM_grid
!-------------------------------------------------------------------------------
subroutine calc_HM_score_using_grid(ligand, grid_info, HM_grid, HM_E, g, calc_g)
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand
type(docking_grid_param_type), intent(in) :: grid_info
real(dp), intent(in) :: HM_grid(:,:)
real(dp), intent(out) :: HM_E, g(:,:)
logical, intent(in) :: calc_g
integer :: i_atm, atm_idx, i
real(dp) :: tmp_E, field(2)
!
real(dp), parameter :: E_min = -0.55d0
real(dp), parameter :: E_max = -0.50d0

HM_E = 0.0d0
if (calc_g) g = 0.0d0

do i_atm = 1, ligand%n_atm
    atm_idx = ii_L(i_atm, ligand%lig_no)
    if (logP(atm_idx) < 0.0) cycle
    call calc_E_using_grid(R(1:3, atm_idx), grid_info, HM_grid, tmp_E, &
                           g(:,i_atm), calc_g)
    field(:) = sigmoidal1(tmp_E, E_min, E_max, calc_g)
    HM_E = HM_E + field(1)*logP(atm_idx)
    if (calc_g) then
        g(:,i_atm) = logP(atm_idx)*field(2)*g(:,i_atm)
    end if
end do

end subroutine calc_HM_score_using_grid
!-------------------------------------------------------------------------------
subroutine calc_HM_score(f, g, appl_res, calc_g)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: f, g(3, tn%atom)
logical, intent(in) :: appl_res(:), calc_g
!
integer :: i_atm, j_atm, res_no
real(dp) :: dr(3), dist_sqr, logP_sum

f = 0.0d0
do i_atm = tn%nonligatm+1, tn%atom ! for ligand
    if (.not. is_ligand(i_atm)) cycle
    if (logP(i_atm) < 0.0) cycle
    
    logP_sum = 0.0d0
    do j_atm = 1, tn%atom ! for protein/cofactor
        if (is_ligand(j_atm)) cycle

        res_no = i_R(1, j_atm)
        if (.not. appl_res(res_no)) cycle

        dr(1:3) = R(1:3, j_atm) - R(1:3, i_atm)
        dist_sqr = dot_product(dr,dr)
        
        if (dist_sqr > HMcutoff2) cycle
        logP_sum = logP_sum + logP(j_atm)
    end do
    
    if (logP_sum > -0.50) then
        f = f + logP(i_atm)
    end if
end do

end subroutine calc_HM_score
!===============================================================================
! Misc.
!===============================================================================
subroutine find_Xscore_res_idx(res_name, res_type, Xscore_res_info, num_res_info)
!-------------------------------------------------------------------------------
character(len=4), intent(inout) :: res_name
integer, intent(out) :: res_type
type(Xscore_res_info_type), intent(in) :: Xscore_res_info(:)
integer, intent(in) :: num_res_info
integer :: i_res

res_type = 0

if (trim(res_name) == 'HID' .or. trim(res_name) == 'HIE'&
    .or. trim(res_name) == 'HIP') then
    res_name = 'HIS'
else if (trim(res_name) == 'CYM' .or. trim(res_name) == 'CYX') then
    res_name = 'CYS'
else if (trim(res_name) == 'CPR') then
    res_name = 'PRO'
else if (trim(res_name) == 'ZN' .or. trim(res_name) == 'LI') then
    res_name = 'HET'
else if (trim(res_name) == 'NA' .or. trim(res_name) == 'K') then
    res_name = 'HET'
else if (trim(res_name) == 'CA' .or. trim(res_name) == 'MG') then
    res_name = 'HET'
else if (trim(res_name) == 'AL' .or. trim(res_name) == 'MN') then
    res_name = 'HET'
else if (trim(res_name) == 'FE' .or. trim(res_name) == 'NI') then
    res_name = 'HET'
else if (trim(res_name) == 'CD' .or. trim(res_name) == 'CO') then
    res_name = 'HET'
else if (trim(res_name) == 'CU' .or. trim(res_name) == 'HG') then
    res_name = 'HET'
else if (trim(res_name) == 'F' .or. trim(res_name) == 'CL') then
    res_name = 'HET'
else if (trim(res_name) == 'BR' .or. trim(res_name) == 'I') then
    res_name = 'HET'
end if

do i_res = 1, num_res_info
    if (trim(res_name) == trim(Xscore_res_info(i_res)%res_name)) then
        res_type = i_res
        return
    end if
end do

write(log_msg,'(A,A)') 'Failed to get logP value of this residue, ', res_name
call terminate_with_error(log_msg)

end subroutine find_Xscore_res_idx
!-------------------------------------------------------------------------------
subroutine find_Xscore_atm_idx(atom_name, res_type, atm_type, Xscore_res_info)
!-------------------------------------------------------------------------------
character(len=4), intent(inout) :: atom_name
integer, intent(in) :: res_type
integer, intent(out) :: atm_type
type(Xscore_res_info_type), intent(in) :: Xscore_res_info(:)
character(len=4) :: res_name
integer :: i_atm

atm_type = 0
res_name = Xscore_res_info(res_type)%res_name

if (trim(res_name) == 'ACE' .or. trim(res_name) == 'NME') then
   if(trim(atom_name) == 'CH3') then
      atom_name = 'CA'
   else if(trim(atom_name) == 'HH31') then
      atom_name = 'HA1'
   else if(trim(atom_name) == 'HH32') then
      atom_name = 'HA2'
   else if(trim(atom_name) == 'HH33') then
      atom_name = 'HA3'
   end if
! It's for considering different atom name between amber_topo and
! Xscore_res_info.
else if(trim(res_name) /= 'ALA' .and. trim(res_name) /= 'VAL' &
        .and. trim(res_name) /= 'THR') then
   if(trim(atom_name) == 'HA3') then
      atom_name = 'HA1'
   else if(trim(atom_name) == 'HB3') then
      atom_name = 'HB1'
   else if(trim(atom_name) == 'HG3') then
      atom_name = 'HG1'
   else if(trim(atom_name) == 'HD3') then
      atom_name = 'HD1'
   else if(trim(res_name) == 'LYS' .and. trim(atom_name) == 'HE3') then
      atom_name = 'HE1'
   else if(trim(atom_name) == 'HG13') then
      atom_name = 'HG11'
   end if
end if

do i_atm = 1, Xscore_res_info(res_type)%n_atm 
   if (trim(atom_name) == trim(Xscore_res_info(res_type)%atom_name(i_atm))) then
      atm_type = i_atm
      return
   end if
end do

write(log_msg,'(A,A,A,A)') 'Failed to get logP value of this atom, ', atom_name, &
                           ' in this residue, ', res_name
call terminate_with_error(log_msg)

end subroutine find_Xscore_atm_idx
!-------------------------------------------------------------------------------
subroutine FFT_rec_grid_HM(grid_info, dock_grid, FFT_grid_info, r_grid)
!-------------------------------------------------------------------------------
type(docking_grid_param_type), intent(in) :: grid_info, FFT_grid_info
type(docking_grid_type), intent(in) :: dock_grid
type(docking_FFT_grid_type), intent(inout) :: r_grid
!
integer :: x, y, z, hash, FFT_hash

FFT_hash = FFT_grid_info%n_elem(1)*FFT_grid_info%n_elem(2)*FFT_grid_info%n_elem(3)
allocate(r_grid%HM_grid(FFT_hash))

r_grid%HM_grid = cmplx(0.0,0.0)

do z = 1, FFT_grid_info%n_elem(3)
    do y = 1, FFT_grid_info%n_elem(2)
        do x = 1, FFT_grid_info%n_elem(1)
            !
            FFT_hash = FFT_grid_info%n_elem(2)*FFT_grid_info%n_elem(1)*(z-1) &
                     + FFT_grid_info%n_elem(1)*(y-1) + x
            !
            if ((x > grid_info%n_elem(1)) .or. &
                (y > grid_info%n_elem(2)) .or. &
                (z > grid_info%n_elem(3))) then  ! outside of docking grid
                r_grid%HM_grid(FFT_hash) = cmplx(-1.0*max_energy, 0.0)
                cycle
            end if
            !
            hash = grid_info%n_elem(2)*grid_info%n_elem(1)*(z-1) &
                 + grid_info%n_elem(1)*(y-1) + x
            if (dock_grid%HM_grid(1,hash) > -0.5) then
                r_grid%HM_grid(FFT_hash) = cmplx(1.0, 0.0)
            end if
        end do
    end do
end do

end subroutine FFT_rec_grid_HM
!-------------------------------------------------------------------------------
subroutine FFT_lig_grid_HM(FFT_grid_info, grid_info, ligand, &
                           fragment, n_frag_atm, R_frag, l_grid, ordering)
!-------------------------------------------------------------------------------
type(docking_grid_param_type), intent(in) :: FFT_grid_info, grid_info
type(ligand_type), intent(in) :: ligand
integer, intent(in) :: fragment(:), n_frag_atm
real(dp), intent(in) :: R_frag(:,:)
type(docking_FFT_grid_type), intent(inout) :: l_grid
logical, intent(in) :: ordering
integer :: n_type, i_type, j_type, lig_type
integer :: i_atm, atm_no
integer :: x, y, z, X0(3), FFT_hash
real(dp) :: dr(3), tcrd(3)

if (.not. allocated(l_grid%HM_grid)) then
    FFT_hash = FFT_grid_info%n_elem(1)*FFT_grid_info%n_elem(2)*FFT_grid_info%n_elem(3)
    allocate(l_grid%HM_grid(FFT_hash))
end if

if (ordering) return

l_grid%HM_grid = cmplx(0.0,0.0)

! Lp(i,j,k) = logP or 0
do i_atm = 1, n_frag_atm
    atm_no = ii_L(fragment(i_atm), ligand%lig_no)
    if (logP(atm_no) < 0.0) cycle
    !
    dr(:) = R_frag(:,i_atm) - FFT_grid_info%grid_cntr(:)
    tcrd(:) = dr(:)/FFT_grid_info%grid_width
    tcrd(:) = tcrd(:) + (dble(grid_info%n_elem(:)-1)*0.5d0 + 1.0d0)
    X0(:) = nint(tcrd(:)) ! nearest grid point
    !
    x = X0(1)
    y = X0(2)
    z = X0(3)
    !
    FFT_hash = FFT_grid_info%n_elem(2)*FFT_grid_info%n_elem(1)*(z-1) &
             + FFT_grid_info%n_elem(1)*(y-1) + x
    !
    l_grid%HM_grid(FFT_hash) = l_grid%HM_grid(FFT_hash) + cmplx(logP(atm_no),0.0)
end do

end subroutine FFT_lig_grid_HM
!-------------------------------------------------------------------------------
subroutine get_hydrophobicity(lig_no, fragment, n_frag_atm, hp)
!-------------------------------------------------------------------------------
integer, intent(in) :: fragment(:), n_frag_atm, lig_no
real(dp), intent(out) :: hp
!
integer :: i_frag, i_atm, atm_no

hp = 0.0
do i_atm = 1, n_frag_atm
    atm_no = ii_L(fragment(i_atm), lig_no)
    hp = hp + logp(atm_no)
end do

end subroutine get_hydrophobicity
!-------------------------------------------------------------------------------
END MODULE XSCORE
!-------------------------------------------------------------------------------
