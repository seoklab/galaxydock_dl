!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/modeling/knowledge/knowledge_gb.f90
!
! Description: KGB (Knowledge-based potential with GB) energy
!-------------------------------------------------------------------------------
MODULE KNOWLEDGE_GB
!-------------------------------------------------------------------------------
use globals
use logger, only: log_p, terminate_with_error
use string, only: parse_string
use in_out_utils, only: find_atom_idx
use mathfunctions, only: cubic_spline, v_norm
use convert_res_name, only: convert_to_stdres
!
use energy_vars, only: R, i_R, ii_R, i_P, LRoff, max_neigh2, max_energy,&
                       kgb_water_file, use_water_kgb, use_water_kgb_only, &
                       reinitialize, wkgb_w
use facts, only: reff

implicit none
private

!-------------------------------------------------------------------------------
! PARAMETERS
!-------------------------------------------------------------------------------
! parameters for atom_types
integer, parameter :: tot_n_atom_type = 167
integer :: n_atom_type
integer :: p_atom_type
integer :: w_atom_index
!
! parameters for proteins
integer, parameter :: tot_n_rbin = 30
integer, parameter :: tot_n_abin = 6
integer :: tot_n_solv = 3
!
! parameters for water-potential
integer, parameter :: tot_n_rbin_w = 21
integer, parameter :: tot_n_abin_w = 6
integer, parameter :: tot_n_solv_w = 6
!integer, parameter :: tot_n_solv_w = 1

real(dp), parameter :: dx_radial  = 0.5d0
integer, parameter :: sequence_separation = 7

!-------------------------------------------------------------------------------
! DERIVED TYPES
!-------------------------------------------------------------------------------
type kgb_atom_type
!-------------------------------------------------------------------------------
integer :: i_n, i_p             ! indices for nn, np table
character(len=4) :: resName
character(len=4) :: atmName
integer :: nC, iC(2)            ! # of connected atoms, and their index
logical :: is_polar
!-------------------------------------------------------------------------------
end type kgb_atom_type
!-------------------------------------------------------------------------------
type kgb_score_type
!-------------------------------------------------------------------------------
real :: f(tot_n_rbin)
real :: fpp(tot_n_rbin)

end type kgb_score_type
!-------------------------------------------------------------------------------
type wkgb_score_type
!-------------------------------------------------------------------------------
real :: f(tot_n_rbin_w)
real :: fpp(tot_n_rbin_w)

end type wkgb_score_type
!-------------------------------------------------------------------------------
type kgb_index_type
!-------------------------------------------------------------------------------
logical :: disabled
integer :: ig, i_n, i_p
integer :: nC, iC(2)
logical :: is_polar
logical :: is_water

end type kgb_index_type
!-------------------------------------------------------------------------------
type(kgb_score_type), allocatable :: nn_table(:,:,:)    ! i_solv, ib, ia
type(kgb_score_type), allocatable :: np_table(:,:,:,:)  ! i_ang, i_solv, ib, ia
type(kgb_score_type), allocatable :: pp_table(:,:,:,:)  ! i_ang, i_solv, ib, ia
!
type(wkgb_score_type), allocatable :: wn_table(:,:)     ! i_solv, ia
type(wkgb_score_type), allocatable :: wp_table(:,:,:)   ! i_ang, i_solv, ia

type(kgb_index_type), allocatable :: i_G(:)

integer :: map_index(tot_n_atom_type+1)

public :: initialize_kgb_energy
public :: finalize_kgb_energy
!
public :: calc_kgb_energy
public :: calc_wkgb_energy_grid
public :: calc_wkgb_energy_grid_decompose

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_kgb_energy(protein)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(kgb_atom_type) :: kgb_atoms(tot_n_atom_type)
character(len=len_fname) :: infile_kgb_atoms
character(len=len_fname) :: infile_kgb_score_file
character(len=len_fname) :: infile_wkgb_score_file
integer :: i_het
logical,save :: reset_spline = .true.

if (.not. use_water_kgb) then
    do i_het = 1, protein%n_het
        if (trim(protein%hetmol(i_het)%pdb_res_name) == 'HOH' .or. &
            trim(protein%hetmol(i_het)%pdb_res_name) == 'WAT' ) then
            use_water_kgb = .true.
            exit
        end if
    end do
end if
!
infile_kgb_atoms = trim(data_dir) // 'kgb_atomtype.list'
!
if (use_water_kgb) then
    infile_kgb_score_file = trim(data_dir) // 'kgb_p30r18ex.151112.score'
    if (trim(kgb_water_file) == '') then
        infile_wkgb_score_file = trim(data_dir) // 'kgb_p30r18wat.151205.score'
    else
        infile_wkgb_score_file = trim(kgb_water_file)
    end if
else
    infile_kgb_score_file = trim(data_dir) // 'kgb_p30r18im.151014.score'
end if

tot_n_solv = 3

!
if (.not. reinitialize) then
    write(log_msg, "(A)") "  Reading KGB energy atom list"
    call log_p(log_msg, me=me, level=40)
    call read_kgb_atom_list(infile_kgb_atoms, kgb_atoms)
    !
    if (.not. use_water_kgb_only) then
        write(log_msg, "(A)") "  Reading KGB energy table"
        call log_p(log_msg, me=me, level=40)
        allocate(nn_table(tot_n_solv, n_atom_type, n_atom_type))
        allocate(np_table(tot_n_abin, tot_n_solv, n_atom_type, p_atom_type))
        allocate(pp_table(tot_n_abin, tot_n_solv, p_atom_type, p_atom_type))
        call read_kgb_score_file(infile_kgb_score_file)
        reset_spline = .true.
    end if
end if
!
if (use_water_kgb) then
    if (.not. allocated(wn_table)) then
        write(log_msg, "(A)") "  Reading wKGB energy table"
        call log_p(log_msg, me=me, level=40)
        w_atom_index = n_atom_type+1
        allocate(wn_table(tot_n_solv_w, w_atom_index)) ! +1 for water-water
        allocate(wp_table(tot_n_abin_w, tot_n_solv_w, p_atom_type))
        call read_wkgb_score_file(infile_wkgb_score_file)
        reset_spline = .true.
    end if
end if
!
write(log_msg, "(A)") "  Setting up the protein for KGB energy evaluation"
call log_p(log_msg, me=me, level=40)
if (use_water_kgb) then
    allocate(i_G(tn%stdatm+tn%hetatm))
else
    allocate(i_G(tn%stdatm))
end if
call setup_protein_for_kgb(protein, kgb_atoms)
if (reset_spline) call setup_radial_spline()
reset_spline = .false.

end subroutine initialize_kgb_energy
!-------------------------------------------------------------------------------
subroutine finalize_kgb_energy()
!-------------------------------------------------------------------------------
deallocate(i_G)
!
if (.not. reinitialize) then
    if (allocated(nn_table)) then
        deallocate(nn_table)
        deallocate(np_table)
        deallocate(pp_table)
    end if
end if

if (allocated(wn_table)) then
    deallocate(wn_table)
    deallocate(wp_table)
end if

end subroutine finalize_kgb_energy
!-------------------------------------------------------------------------------
subroutine calc_kgb_energy(ff, gg, appl_respair, calc_g)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: ff(2), gg(3,tn%atom,2)
logical, intent(in) :: appl_respair(tn%residue,tn%residue), calc_g

integer :: i_atm, j_atm, ia, n_pair, i_pair(max_neigh2)
real(dp) :: Ri(3), Rj(3), Rij(3), reff_i, reff_j, reff_ij
real(dp) :: pv(3, tn%stdatm)
real(dp) :: dij, dr, dfdr(3)
real(dp) :: fb(2), fb_pp(2), fg(2)
integer :: ir, is, it, is_s(4)
logical :: use_ang(max_neigh2)
!
ff(:) = 0.0
gg(:,:,:) = 0.0d0

! Define polar vectors
do i_atm = 1, tn%stdatm
    if (i_G(i_atm)%disabled .or. (.not. i_G(i_atm)%is_polar)) cycle
    call define_polar_vector(i_G(i_atm), R(:,i_atm), pv(:,i_atm))
end do

do i_atm = 1, tn%stdatm
    if (i_G(i_atm)%disabled) cycle
    Ri(1:3) = R(1:3, i_atm)
    reff_i = reff(i_atm)
    call get_kgb_energy_pair(i_atm, appl_respair, n_pair, i_pair, use_ang)

    do ia = 1, n_pair
        j_atm = i_pair(ia)
        Rj(1:3) = R(1:3, j_atm)
        reff_j = reff(j_atm)
        Rij = (Rj - Ri)
        dij = sqrt(dot_product(Rij, Rij))
        Rij = Rij/dij
        !
        reff_ij = sqrt(reff_i*reff_j)
        is = max(1, min(tot_n_solv, 1+int((reff_ij-1.0)/3.0)))
        ir = int(2.0*dij - 0.5) + 1
        dr = dble(ir) - 2.0*dij + 0.5
             
        ! distance-dependent term
        fb(1:2)    = nn_table(is, i_G(j_atm)%i_n, i_G(i_atm)%i_n)%f(ir:ir+1)
        fb_pp(1:2) = nn_table(is, i_G(j_atm)%i_n, i_G(i_atm)%i_n)%fpp(ir:ir+1)
        !
        fg = cubic_spline(dx_radial, fb, fb_pp, dr, calc_g)
        ff(1) = ff(1) + fg(1)
        if (calc_g) then
            dfdr = fg(2)*Rij
            gg(:,i_atm,1) = gg(:,i_atm,1) - dfdr
            gg(:,j_atm,1) = gg(:,j_atm,1) + dfdr
        end if

        if (.not. use_ang(ia)) cycle

        if (i_G(i_atm)%is_polar) then
            it = max(1, min(tot_n_abin, int((dot_product(pv(:,i_atm), Rij)+1.0)*3.0)+1))
            fb(1:2)    = np_table(it, is, i_G(j_atm)%i_n, i_G(i_atm)%i_p)%f(ir:ir+1)
            fb_pp(1:2) = np_table(it, is, i_G(j_atm)%i_n, i_G(i_atm)%i_p)%fpp(ir:ir+1)
            !
            fg = cubic_spline(dx_radial, fb, fb_pp, dr, calc_g)
            ff(2) = ff(2) + fg(1)
            if (calc_g) then
                dfdr = fg(2)*Rij
                gg(:,i_atm,2) = gg(:,i_atm,2) - dfdr
                gg(:,j_atm,2) = gg(:,j_atm,2) + dfdr
            end if
            !
            if (i_G(j_atm)%is_polar) then
                it = max(1, min(tot_n_abin, int((-dot_product(pv(:,j_atm), Rij)+1.0)*3.0)+1))
                fb(1:2)    = np_table(it, is, i_G(i_atm)%i_n, i_G(j_atm)%i_p)%f(ir:ir+1)
                fb_pp(1:2) = np_table(it, is, i_G(i_atm)%i_n, i_G(j_atm)%i_p)%fpp(ir:ir+1)
                !
                fg = cubic_spline(dx_radial, fb, fb_pp, dr, calc_g)
                ff(2) = ff(2) + fg(1)
                if (calc_g) then
                    dfdr = fg(2)*Rij
                    gg(:,i_atm,2) = gg(:,i_atm,2) - dfdr
                    gg(:,j_atm,2) = gg(:,j_atm,2) + dfdr
                end if
                !
                it = max(1, min(tot_n_abin, int(( dot_product(pv(:,i_atm), pv(:,j_atm))+1.0)*3.0)+1))
                fb(1:2)    = pp_table(it, is, i_G(j_atm)%i_p, i_G(i_atm)%i_p)%f(ir:ir+1)
                fb_pp(1:2) = pp_table(it, is, i_G(j_atm)%i_p, i_G(i_atm)%i_p)%fpp(ir:ir+1)
                !
                fg = cubic_spline(dx_radial, fb, fb_pp, dr, calc_g)
                ff(2) = ff(2) + fg(1)
                if (calc_g) then
                    dfdr = fg(2)*Rij
                    gg(:,i_atm,2) = gg(:,i_atm,2) - dfdr
                    gg(:,j_atm,2) = gg(:,j_atm,2) + dfdr
                end if
            end if
        else if (i_G(j_atm)%is_polar) then
            it = max(1, min(tot_n_abin, int((-dot_product(pv(:,j_atm), Rij)+1.0)*3.0)+1))
            fb(1:2)    = np_table(it, is, i_G(i_atm)%i_n, i_G(j_atm)%i_p)%f(ir:ir+1)
            fb_pp(1:2) = np_table(it, is, i_G(i_atm)%i_n, i_G(j_atm)%i_p)%fpp(ir:ir+1)
            !
            fg = cubic_spline(dx_radial, fb, fb_pp, dr, calc_g)
            ff(2) = ff(2) + fg(1)
            if (calc_g) then
                dfdr = fg(2)*Rij
                gg(:,i_atm,2) = gg(:,i_atm,2) - dfdr
                gg(:,j_atm,2) = gg(:,j_atm,2) + dfdr
            end if
        end if
    end do
end do

if (.not. use_water_kgb) return

do i_atm = 1, tn%stdatm ! Evaluate protein-water pair
    if (i_G(i_atm)%disabled) cycle
    call get_wkgb_energy_pair(i_atm, appl_respair, n_pair, i_pair)
    if (n_pair == 0) cycle
    !
    Ri(1:3) = R(1:3, i_atm)
    reff_i = reff(i_atm)
    is = max(1, min(tot_n_solv_w, 1+int(reff_i-2.0)))

    do ia = 1, n_pair
        j_atm = i_pair(ia)
        Rj(1:3) = R(1:3, j_atm)
        Rij = (Rj - Ri)
        dij = sqrt(dot_product(Rij, Rij))
        Rij = Rij/dij
        !
        ir = int(2.0*dij - 0.5) + 1
        dr = dble(ir) - 2.0*dij + 0.5
        !
        fb(1:2)    = wn_table(is, i_G(i_atm)%i_n)%f(ir:ir+1)
        fb_pp(1:2) = wn_table(is, i_G(i_atm)%i_n)%fpp(ir:ir+1)
        !
        fg = cubic_spline(dx_radial, fb, fb_pp, dr, calc_g)
        ff(1) = ff(1) + fg(1)*wkgb_w
        if (calc_g) then
            dfdr = fg(2)*Rij*wkgb_w
            gg(:,i_atm,1) = gg(:,i_atm,1) - dfdr
            gg(:,j_atm,1) = gg(:,j_atm,1) + dfdr
        end if
        
        if (i_G(i_atm)%is_polar) then
            it = max(1, min(tot_n_abin_w, int((dot_product(pv(:,i_atm), Rij)+1.0)*3.0)+1))
            fb(1:2)    = wp_table(it, is, i_G(i_atm)%i_p)%f(ir:ir+1)
            fb_pp(1:2) = wp_table(it, is, i_G(i_atm)%i_p)%fpp(ir:ir+1)
            !
            fg = cubic_spline(dx_radial, fb, fb_pp, dr, calc_g)
            ff(2) = ff(2) + fg(1)*wkgb_w
            if (calc_g) then
                dfdr = fg(2)*Rij*wkgb_w
                gg(:,i_atm,2) = gg(:,i_atm,2) - dfdr
                gg(:,j_atm,2) = gg(:,j_atm,2) + dfdr
            end if
        end if
    end do
end do

!return
!
!do i_atm = tn%stdatm+1, tn%stdatm+tn%hetatm-1
!    if (i_G(i_atm)%disabled .or. (.not. i_G(i_atm)%is_water)) cycle
!    Ri(1:3) = R(1:3, i_atm)
!    reff_i = reff(i_atm)
!    is = max(1, min(tot_n_solv_w, 1+int(reff_i-2.0)))
!    !
!    do j_atm = i_atm+1, tn%stdatm+tn%hetatm
!        if (i_G(j_atm)%disabled .or. (.not. i_G(j_atm)%is_water)) cycle
!        Rj(1:3) = R(1:3, j_atm)
!        reff_j = reff(j_atm)
!        js = max(1, min(tot_n_solv_w, 1+int(reff_j-2.0)))
!        !
!        Rij = Rj - Ri
!        dij = sqrt(dot_product(Rij, Rij))
!        if (dij > 10.0d0 .or. dij < 0.25d0) cycle
!        Rij = Rij / dij
!        !
!        ir = int(2.0*dij - 0.5) + 1
!        dr = dble(ir) - 2.0*dij + 0.5
!        !
!        if (is == js) then
!            fb(1:2)    = wn_table(is, w_atom_index)%f(ir:ir+1)
!            fb_pp(1:2) = wn_table(is, w_atom_index)%fpp(ir:ir+1)
!            !
!            fg = cubic_spline(dx_radial, fb, fb_pp, dr, calc_g)
!            ff(1) = ff(1) + fg(1)
!            if (calc_g) then
!                dfdr = fg(2)*Rij
!                gg(:,i_atm,1) = gg(:,i_atm,1) - dfdr
!                gg(:,j_atm,1) = gg(:,j_atm,1) + dfdr
!            end if
!        else
!            fb(1:2)    = wn_table(is, w_atom_index)%f(ir:ir+1)
!            fb_pp(1:2) = wn_table(is, w_atom_index)%fpp(ir:ir+1)
!            !
!            fg = cubic_spline(dx_radial, fb, fb_pp, dr, calc_g)
!            ff(1) = ff(1) + fg(1)*0.5
!            if (calc_g) then
!                dfdr = fg(2)*Rij*0.5
!                gg(:,i_atm,1) = gg(:,i_atm,1) - dfdr
!                gg(:,j_atm,1) = gg(:,j_atm,1) + dfdr
!            end if
!            !
!            fb(1:2)    = wn_table(js, w_atom_index)%f(ir:ir+1)
!            fb_pp(1:2) = wn_table(js, w_atom_index)%fpp(ir:ir+1)
!            !
!            fg = cubic_spline(dx_radial, fb, fb_pp, dr, calc_g)
!            ff(1) = ff(1) + fg(1)*0.5
!            if (calc_g) then
!                dfdr = fg(2)*Rij*0.5
!                gg(:,i_atm,1) = gg(:,i_atm,1) - dfdr
!                gg(:,j_atm,1) = gg(:,j_atm,1) + dfdr
!            end if
!        end if
!    end do
!end do

end subroutine calc_kgb_energy
!-------------------------------------------------------------------------------
subroutine get_kgb_energy_pair(i_atm, appl_respair, n_pair, i_pair, use_ang) 
!-------------------------------------------------------------------------------
integer, intent(in) :: i_atm
logical, intent(in) :: appl_respair(tn%residue,tn%residue)
integer, intent(out) :: n_pair
integer, intent(out) :: i_pair(max_neigh2)
logical, intent(out) :: use_ang(max_neigh2)

integer :: ia_start
integer :: i_res, j_res, j_atm, ia

n_pair = 0
i_res = i_R(1,i_atm)

ia_start = i_P(i_atm)%pair_end_index(2)+1

do ia = ia_start, i_P(i_atm)%n_pair
    j_atm = i_P(i_atm)%i_pair(ia)
    j_res = i_R(1,j_atm)
    if (j_res /= i_res .and. (.not. i_G(j_atm)%disabled) .and. &
        appl_respair(i_res,j_res) .and. (.not. i_G(j_atm)%is_water)) then
        n_pair = n_pair + 1
        i_pair(n_pair) = j_atm
        if (abs(j_res-i_res) >= sequence_separation) then
            use_ang(n_pair) = .true.
        else
            use_ang(n_pair) = .false.
        end if
        if (n_pair == max_neigh2) return
    endif
enddo

do ia = 1, i_P(i_atm)%n_Lpair
    j_atm = i_P(i_atm)%i_Lpair(ia)
    j_res = i_R(1,j_atm)
    if (j_res /= i_res .and. (.not. i_G(j_atm)%disabled) .and. &
        appl_respair(i_res,j_res) .and. (.not. i_G(j_atm)%is_water)) then
        n_pair = n_pair + 1
        i_pair(n_pair) = j_atm
        if (abs(j_res-i_res) >= sequence_separation) then
            use_ang(n_pair) = .true.
        else
            use_ang(n_pair) = .false.
        end if
        if (n_pair == max_neigh2) return
    endif
enddo

end subroutine get_kgb_energy_pair
!-------------------------------------------------------------------------------
subroutine get_wkgb_energy_pair(i_atm, appl_respair, n_pair, i_pair)
!-------------------------------------------------------------------------------
integer, intent(in) :: i_atm
logical, intent(in) :: appl_respair(tn%residue,tn%residue)
integer, intent(out) :: n_pair
integer, intent(out) :: i_pair(max_neigh2)

integer :: ia_start
integer :: i_res, j_res, j_atm, ia

n_pair = 0
i_res = i_R(1,i_atm)

ia_start = i_P(i_atm)%pair_end_index(2)+1

do ia = ia_start, i_P(i_atm)%n_pair
    j_atm = i_P(i_atm)%i_pair(ia)
    j_res = i_R(1, j_atm)
    if (j_res /= i_res .and. (.not. i_G(j_atm)%disabled) .and. &
        appl_respair(i_res,j_res) .and. (i_G(j_atm)%is_water)) then
        n_pair = n_pair + 1
        i_pair(n_pair) = j_atm
        if (n_pair == max_neigh2) return
    endif
end do

do j_atm = tn%stdatm+1, tn%stdatm+tn%hetatm ! iterate on water.
    j_res = i_R(1, j_atm)
    !if ((i_G(j_atm)%disabled) .or. (.not. appl_respair(i_res,j_res))) cycle
    if (i_G(j_atm)%disabled) cycle
    !
    do ia = 1, i_P(j_atm)%n_pair
        !if ((i_P(j_atm)%i_pair(ia) == i_atm) .and. (.not. i_G(i_atm)%disabled)) then
        if (i_P(j_atm)%i_pair(ia) == i_atm) then
            n_pair = n_pair + 1
            i_pair(n_pair) = j_atm
        end if
    end do
end do

end subroutine get_wkgb_energy_pair
!-------------------------------------------------------------------------------
subroutine read_kgb_atom_list(infile_kgb_atoms, kgb_atoms)
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: infile_kgb_atoms
type(kgb_atom_type), intent(out) :: kgb_atoms(:)

character(len=len_fname) :: line, word(10)
integer :: f_unit, io, n_word, ia, i_prev

n_atom_type = 0
p_atom_type = 0

f_unit = 41
open(f_unit, file=trim(infile_kgb_atoms), iostat=io)
if (io < 0) then
    call terminate_with_error("Error: No such file for infile_kgb_atoms")
end if

ia = 0
i_prev = 0
do
    read(f_unit, '(A120)', iostat=io) line
    if (io < 0) exit

    call parse_string(line, n_word, word)

    ia = ia + 1
    kgb_atoms(ia)%resName = word(1)(1:4)
    kgb_atoms(ia)%atmName = word(2)(1:4)
    read(word(3),"(I3)") kgb_atoms(ia)%i_n

    if (trim(word(4)) == 'T') then
        kgb_atoms(ia)%is_polar = .true.

        kgb_atoms(ia)%nC = n_word - 4
        read(word(5), "(I3)") kgb_atoms(ia)%iC(1) 
        if (kgb_atoms(ia)%nC == 2) then
            read(word(6), "(I3)") kgb_atoms(ia)%iC(2) 
        else if (trim(kgb_atoms(ia)%atmName) == 'N') then
            kgb_atoms(ia)%nC = kgb_atoms(ia)%nC + 1
            kgb_atoms(ia)%iC(2) = -1
        end if
    else
        kgb_atoms(ia)%is_polar = .false.
        kgb_atoms(ia)%nC = 0
    endif

    if (i_prev /= kgb_atoms(ia)%i_n) then
        i_prev = kgb_atoms(ia)%i_n
        n_atom_type = n_atom_type + 1
        if (kgb_atoms(ia)%is_polar) then
            p_atom_type = p_atom_type + 1
            map_index(n_atom_type) = p_atom_type
            kgb_atoms(ia)%i_p = p_atom_type
        end if
    else if (kgb_atoms(ia)%is_polar) then
        kgb_atoms(ia)%i_p = p_atom_type
    end if
enddo

close(f_unit)

end subroutine read_kgb_atom_list
!-------------------------------------------------------------------------------
subroutine read_kgb_score_file(infile_kgb_score_file)
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: infile_kgb_score_file

character(len=len_fname) :: line
integer :: f_unit, io
integer :: mode, ia, ib, is, ir
real(dp) :: ff, ff_s(tot_n_abin)

f_unit = 42
open(f_unit, file=trim(infile_kgb_score_file), iostat=io)
if (io < 0) then
    call terminate_with_error("Error: No such file for infile_kgb_score_file")
end if

do
    read(f_unit, "(A120)", iostat=io) line
    if (io < 0) exit

    if (line(1:4) == 'TYPE') then
        if (line(6:7) == 'nn') then
            mode = 0
        else if (line(6:7) == 'np') then
            mode = 1
        else if (line(6:7) == 'pp') then
            mode = 2
        end if
    else
        if (mode == 0) then
            read(line(1:22), '(2(I3,1x),I1,1x,I2,1x,F9.3)') ia,ib,is,ir,ff
            nn_table(is,ib,ia)%f(ir) = ff*0.001

        else
            read(line(1:72), '(2(I3,1x),I1,1x,I2,6(1x,F9.3))') ia,ib,is,ir,ff_s
            if (mode == 1) then
                ia = map_index(ia)
                if (ia == 0) cycle
                np_table(1:tot_n_abin,is,ib,ia)%f(ir) = ff_s(1:tot_n_abin)*0.001
            else
                ia = map_index(ia)
                ib = map_index(ib)
                if (ia == 0 .or. ib == 0) cycle
                pp_table(1:tot_n_abin,is,ib,ia)%f(ir) = ff_s(1:tot_n_abin)*0.001
            end if
        end if
    end if 
enddo

close(f_unit)

end subroutine read_kgb_score_file
!-------------------------------------------------------------------------------
subroutine read_wkgb_score_file(infile_wkgb_score_file)
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: infile_wkgb_score_file

character(len=len_fname) :: line
integer :: f_unit, io
integer :: mode, ia, is, ir, i_ang
real(dp) :: ff_s(tot_n_solv_w)

do ia = 1, n_atom_type
    do is = 1, tot_n_solv_w
        wn_table(is, ia)%f(:) = 0.0
        wn_table(is, ia)%fpp(:) = 0.0
    end do
end do

do ia = 1, p_atom_type
    do is = 1, tot_n_solv_w
        wn_table(is, ia)%f(:) = 0.0
        wn_table(is, ia)%fpp(:) = 0.0
        do i_ang = 1, tot_n_abin_w
            wp_table(i_ang, is, ia)%f(:) = 0.0
            wp_table(i_ang, is, ia)%fpp(:) = 0.0
        end do
    end do
end do

f_unit = 44
open(f_unit, file=trim(infile_wkgb_score_file), iostat=io)
if (io < 0) then
    call terminate_with_error("Error: No such file for infile_wkgb_score_file")
end if

do
    read(f_unit, "(A120)", iostat=io) line
    if (io < 0) exit

    if (line(1:4) == 'TYPE') then
        if (line(6:7) == 'nn') then
            mode = 0
        else if (line(6:7) == 'np') then
            mode = 1
        end if
    else
        if (mode == 0) then
            read(line(1: 3), "(I3)") ia
            read(line(5: 6), "(I2)") ir
            read(line(8:63), "(7(1x,f7.3))") ff_s

            wn_table(1:tot_n_solv_w, ia)%f(ir) = ff_s(1:tot_n_solv_w)

        else
            read(line( 1: 3), "(I3)") ia
            read(line( 5: 6), "(I2)") ir
            read(line( 8: 8), "(I1)") i_ang
            read(line(10:65), "(7(1x,f7.3))") ff_s

            ia = map_index(ia)
            if (ia == 0) cycle

            wp_table(i_ang, 1:tot_n_solv_w, ia)%f(ir) = ff_s(1:tot_n_solv_w)
        end if
    end if 
enddo

close(f_unit)

end subroutine read_wkgb_score_file
!-------------------------------------------------------------------------------
subroutine setup_protein_for_kgb(protein, kgb_atoms)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(kgb_atom_type), intent(in) :: kgb_atoms(:)

integer :: ia, ia_s(max_atm), tna
integer :: i_ref, i_res, i_atm, j_atm, ig, ic, i
integer :: i_ref_prev, i_res_prev
character(len=4) :: resName, atmName
character(len=6), parameter :: error_mode = "ignore"

if (use_water_kgb) then
    tna = tn%stdatm + tn%hetatm
else
    tna = tn%stdatm
end if
do i_atm = 1, tna
    i_G(i_atm)%disabled = .false.
    i_G(i_atm)%ig = 0
    i_G(i_atm)%is_polar = .false.
    i_G(i_atm)%is_water = .false.
    i_G(i_atm)%iC(:) = -1
    i_G(i_atm)%nC  = 0
    i_G(i_atm)%i_n = 0
    i_G(i_atm)%i_p = 0
end do

do i_res = 1, protein%n_res
    i_ref = protein%residue(i_res)%res_type
    resName = ref_res(i_ref)%res_name
    call convert_to_stdres(resName)

    do i_atm = 1, ref_res(i_ref)%n_atm
        ia = ii_R(i_atm, i_res)
        ia_s(i_atm) = ia
        atmName = ref_res(i_ref)%atom_name(i_atm)
        if (atmName(1:3) == 'OXT' .or. atmName(1:2) == 'OT') then
            atmName = 'O'
        end if
        call find_kgb_atom_index(kgb_atoms, resName, atmName, ig)
        if (ig == -1) then
            i_G(ia)%disabled = .true.
            cycle
        end if

        i_G(ia)%ig  = ig
        i_G(ia)%nC  = kgb_atoms(ig)%nC
        i_G(ia)%i_n = kgb_atoms(ig)%i_n
        i_G(ia)%i_p = kgb_atoms(ig)%i_p
        i_G(ia)%is_polar = kgb_atoms(ig)%is_polar
    end do

    do i_atm = 1, ref_res(i_ref)%n_atm
        ia = ia_s(i_atm)
        ig = i_G(ia)%ig
        if (i_G(ia)%disabled) cycle

        do i = 1, i_G(ia)%nC
            ic = kgb_atoms(ig)%iC(i)
            if (ic /= -1) then  ! curr residue atom
                do j_atm = 1, ref_res(i_ref)%n_atm
                    if (i_G(ia_s(j_atm))%i_n == ic) then
                        i_G(ia)%iC(i) = ia_s(j_atm)
                        exit
                    end if
                end do
            else                ! -C
                if (protein%residue(i_res)%ter_type /= 'N') then
                    i_res_prev = i_res-1
                    i_ref_prev = protein%residue(i_res_prev)%res_type
                    call find_atom_idx(i_res_prev, i_ref_prev, "C   ", j_atm, error_mode)
                    i_G(ia)%iC(i) = ii_R(j_atm, i_res_prev)
                end if
            end if

            if (i_G(ia)%iC(i) < 0) then
                i_G(ia)%disabled = .true.
                exit
            end if
        end do
    end do
end do

if (.not. use_water_kgb) return

do ia = tn%stdatm+1, tn%stdatm+tn%hetatm
    i_G(ia)%disabled = .true.
end do

do i_res = 1, protein%n_het
    i_ref = protein%hetmol(i_res)%res_type
    resName = ref_res(i_ref)%res_name
    if (trim(resName) /= 'HOH' .and. trim(resName) /= 'WAT') cycle
    !
    do i_atm = 1, ref_res(i_ref)%n_atm
        atmName = ref_res(i_ref)%atom_name(i_atm)
        if (trim(atmName) == 'O') then
            ia = ii_R(i_atm, i_res+protein%n_res)
            i_G(ia)%ig = w_atom_index
            i_G(ia)%nC = kgb_atoms(w_atom_index)%nC
            i_G(ia)%i_n = kgb_atoms(w_atom_index)%i_n
            i_G(ia)%i_p = kgb_atoms(w_atom_index)%i_p
            i_G(ia)%is_polar = .false.
            i_G(ia)%is_water = .true.
            i_G(ia)%disabled = .false.
        end if
    end do
end do

end subroutine setup_protein_for_kgb
!-------------------------------------------------------------------------------
subroutine setup_radial_spline()
!-------------------------------------------------------------------------------
integer :: ia, ja, it, is
real(dp) :: f(tot_n_rbin), fpp(tot_n_rbin), x(tot_n_rbin)

do ia = 1, tot_n_rbin
    x(ia) = dble(ia)*0.5 - 0.25
end do

if (.not. use_water_kgb_only) then
    ! nn_table
    do ia = 1, n_atom_type
        do ja = 1, n_atom_type
            do is = 1, tot_n_solv
                f(1:tot_n_rbin) = nn_table(is, ja, ia)%f(1:tot_n_rbin)
                call get_spline_derivatives(tot_n_rbin, f, x, fpp)
                nn_table(is, ja, ia)%fpp(:) = fpp(1:tot_n_rbin)
            end do
        end do
    end do

    ! np_table
    do ia = 1, p_atom_type
        do ja = 1, n_atom_type
            do is = 1, tot_n_solv
                do it = 1, tot_n_abin
                    f(1:tot_n_rbin) = np_table(it, is, ja, ia)%f(1:tot_n_rbin)
                    call get_spline_derivatives(tot_n_rbin, f, x, fpp)
                    np_table(it, is, ja, ia)%fpp(1:tot_n_rbin) = fpp(1:tot_n_rbin)
                end do
            end do
        end do
    end do

    ! pp_table
    do ia = 1, p_atom_type
        do ja = 1, p_atom_type
            do is = 1, tot_n_solv
                do it = 1, tot_n_abin
                    f(1:tot_n_rbin) = pp_table(it, is, ja, ia)%f(1:tot_n_rbin)
                    call get_spline_derivatives(tot_n_rbin, f, x, fpp)
                    pp_table(it, is, ja, ia)%fpp(1:tot_n_rbin) = fpp(1:tot_n_rbin)
                end do
            end do
        end do
    end do
end if

if (.not. use_water_kgb) return

! wn_table
do ia = 1, w_atom_index    ! +1 for water-water
    do is = 1, tot_n_solv_w
        f(1:tot_n_rbin_w) = wn_table(is, ia)%f(1:tot_n_rbin_w)
        call get_spline_derivatives(tot_n_rbin_w, f(1:tot_n_rbin_w), x, &
                                    fpp(1:tot_n_rbin_w))
        wn_table(is, ia)%fpp(1:tot_n_rbin_w) = fpp(1:tot_n_rbin_w)
    end do
end do

! wp_table
do ia = 1, p_atom_type
    do it = 1, tot_n_abin
        do is = 1, tot_n_solv_w
            f(1:tot_n_rbin_w) = wp_table(it, is, ia)%f(1:tot_n_rbin_w)
            call get_spline_derivatives(tot_n_rbin_w, f(1:tot_n_rbin_w), x, &
                                        fpp(1:tot_n_rbin_w))
            wp_table(it, is, ia)%fpp(1:tot_n_rbin_w) = fpp(1:tot_n_rbin_w)
        end do
    end do
end do

end subroutine setup_radial_spline
!-------------------------------------------------------------------------------
subroutine get_spline_derivatives(n_bin, f, x, fpp)
!-------------------------------------------------------------------------------
integer, intent(in) :: n_bin
real(dp), intent(in) :: f(n_bin), x(n_bin)
real(dp), intent(out) :: fpp(n_bin)

real(dp) :: u(n_bin), p, sig
integer :: i

fpp(1:n_bin) = 0.0d0
u(1) = 0.0d0

do i = 2, n_bin-1
    sig = (x(i)-x(i-1)) / (x(i+1)-x(i-1))
    p = sig*fpp(i-1) + 2.0d0
    fpp(i) = (sig-1.0d0) / p
    u(i) = (6.0d0*((f(i+1) - f(i))/(x(i+1)-x(i)) - (f(i)-f(i-1))/(x(i)-x(i-1))) / &
           (x(i+1)-x(i-1)) - sig*u(i-1)) / p
enddo

do i = n_bin-1, 1, -1
    fpp(i) = fpp(i)*fpp(i+1) + u(i)
enddo

end subroutine get_spline_derivatives
!-------------------------------------------------------------------------------
subroutine find_kgb_atom_index(kgb_atoms, resName, atmName, i_atm)
!-------------------------------------------------------------------------------
type(kgb_atom_type), intent(in) :: kgb_atoms(:)
character(len=4), intent(in) :: resName, atmName
integer, intent(out) :: i_atm
integer :: ia

i_atm = -1
do ia = 1, tot_n_atom_type
    if (resName == kgb_atoms(ia)%resName .and. &
        atmName == kgb_atoms(ia)%atmName) then
        i_atm = ia
        return
    end if
end do

end subroutine find_kgb_atom_index
!-------------------------------------------------------------------------------
subroutine define_polar_vector(i_G, Ri, v)
!-------------------------------------------------------------------------------
type(kgb_index_type), intent(in) :: i_G
real(dp), intent(in)  :: Ri(3)
real(dp), intent(out) :: v(3)
real(dp) :: v1(3), v2(3)

if (i_G%nC == 1) then
    v(1:3) = R(1:3,i_G%iC(1)) - Ri(1:3)
else
    v1(1:3) = R(1:3,i_G%iC(1)) - Ri(1:3)
    v2(1:3) = R(1:3,i_G%iC(2)) - Ri(1:3)
    call v_norm(v1)
    call v_norm(v2)
    v = v1 + v2
end if

call v_norm(v)

end subroutine define_polar_vector
!-------------------------------------------------------------------------------
subroutine calc_wkgb_energy_grid(ff, n_grid, n_pair, i_pair, d_pair, Rij)
!-------------------------------------------------------------------------------
real, intent(out) :: ff(n_grid)
integer, intent(in) :: n_grid, n_pair(n_grid), i_pair(max_neigh2, n_grid)
real, intent(in) :: d_pair(max_neigh2, n_grid), Rij(3, max_neigh2, n_grid)
!
integer :: is_s(tn%atom)
real(dp) :: pv(3, tn%atom)
integer :: i, ia, i_atm, ir, it
real(dp) :: dr, fb(2), fb_pp(2), fg(2)

ff(:) = 0.0

! Define polar vectors
!do i_atm = 1, tn%stdatm
is_s(:) = 0
do i_atm = 1, tn%atom
    if (i_G(i_atm)%disabled) cycle
    is_s(i_atm) = max(1, min(tot_n_solv_w, 1+int(reff(i_atm)-2.0)))
    if (i_G(i_atm)%is_polar) then
        call define_polar_vector(i_G(i_atm), R(:,i_atm), pv(:,i_atm))
    end if
end do

do ia = 1, n_grid
    do i = 1, n_pair(ia)
        i_atm = i_pair(i, ia)
        if (i_G(i_atm)%disabled) cycle
        !
        ir = min(20, int(2.0*d_pair(i,ia) - 0.5) + 1)
        dr = dble(ir) - 2.0*d_pair(i,ia) + 0.5
        !
        !write(*,'(I4,1x,f12.7)') ir, d_pair(i,ia)
        fb(1:2) = wn_table(is_s(i_atm), i_G(i_atm)%i_n)%f(ir:ir+1)
        fb_pp(1:2) = wn_table(is_s(i_atm), i_G(i_atm)%i_n)%fpp(ir:ir+1)
        !
        fg = cubic_spline(dx_radial, fb, fb_pp, dr, .false.)
        ff(ia) = ff(ia) + fg(1)
        !
        if (i_G(i_atm)%is_polar) then
            it = max(1, min(tot_n_abin_w, int((dot_product(pv(:,i_atm), Rij(:,i,ia))+1.0)*3.0)+1))
            fb(1:2)    = wp_table(it, is_s(i_atm), i_G(i_atm)%i_p)%f(ir:ir+1)
            fb_pp(1:2) = wp_table(it, is_s(i_atm), i_G(i_atm)%i_p)%fpp(ir:ir+1)
            !
            fg = cubic_spline(dx_radial, fb, fb_pp, dr, .false.)
            ff(ia) = ff(ia) + fg(1)
        end if
    end do
end do

end subroutine calc_wkgb_energy_grid
!-------------------------------------------------------------------------------
subroutine calc_wkgb_energy_grid_decompose(protein, n_grid, n_pair, i_pair, d_pair, Rij)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
integer, intent(in) :: n_grid, n_pair(n_grid), i_pair(max_neigh2, n_grid)
real, intent(in) :: d_pair(max_neigh2, n_grid), Rij(3, max_neigh2, n_grid)
!
integer :: i, ia, i_atm, is_s(tn%atom), ir, it, i_res, i_ref
real(dp) :: pv(3, tn%stdatm), dr, fb(2), fb_pp(2), fg(2)
real(dp) :: ff(max_neigh2,n_grid)

ff(:,:) = 0.0d0

! Define polar vectors
do i_atm = 1, tn%stdatm
    if (i_G(i_atm)%disabled) cycle
    is_s(i_atm) = max(1, min(tot_n_solv_w, 1+int(reff(i_atm)-2.0)))
    if (i_G(i_atm)%is_polar) then
        call define_polar_vector(i_G(i_atm), R(:,i_atm), pv(:,i_atm))
    end if
end do

do ia = 1, n_grid
    do i = 1, n_pair(ia)
        i_atm = i_pair(i, ia)
        if (i_G(i_atm)%disabled) cycle
        !
        ir = min(20, int(2.0*d_pair(i,ia) - 0.5) + 1)
        dr = dble(ir) - 2.0*d_pair(i,ia) + 0.5
        !
        fb(1:2) = wn_table(is_s(i_atm), i_G(i_atm)%i_n)%f(ir:ir+1)
        fb_pp(1:2) = wn_table(is_s(i_atm), i_G(i_atm)%i_n)%fpp(ir:ir+1)
        !
        fg = cubic_spline(dx_radial, fb, fb_pp, dr, .false.)
        ff(i,ia) = ff(i,ia) + fg(1)
        !
        if (i_G(i_atm)%is_polar) then
            it = max(1, min(tot_n_abin_w, int((dot_product(pv(:,i_atm), Rij(:,i,ia))+1.0)*3.0)+1))
            fb(1:2)    = wp_table(it, is_s(i_atm), i_G(i_atm)%i_p)%f(ir:ir+1)
            fb_pp(1:2) = wp_table(it, is_s(i_atm), i_G(i_atm)%i_p)%fpp(ir:ir+1)
            !
            fg = cubic_spline(dx_radial, fb, fb_pp, dr, .false.)
            ff(i,ia) = ff(i,ia) + fg(1)
        end if
        !
        i_res = i_R(1, i_pair(i,ia))
        i_ref = protein%residue(i_res)%res_type
        write(*,'(I4,1x,f12.7,1x,I4,1x,A4,3(1x,I3),1x,f7.3,2(1x,f12.7))') ia, ff(i,ia), i_res, &
            ref_res(i_ref)%atom_name(i_R(2,i_pair(i,ia))),&
            i_G(i_atm)%i_n, ir, is_s(i_atm), reff(i_atm), wn_table(is_s(i_atm), i_G(i_atm)%i_n)%f(ir:ir+1)
    end do
    !
    write(*,'(I4,1x,f12.7,1x,A)') ia, sum(ff(1:n_pair(ia),ia)), 'sum'
    write(*,'(A)') '#----------------------------------'
end do

end subroutine calc_wkgb_energy_grid_decompose
!-------------------------------------------------------------------------------
END MODULE KNOWLEDGE_GB
!-------------------------------------------------------------------------------
