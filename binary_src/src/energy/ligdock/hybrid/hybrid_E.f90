!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/ligdock/hybrid/hybrid_E.f90
!
! Description: This module contains subroutines related to
!   1. construct energy table, T_prot(i_rot)
!       T_prot(i_rot) = prot_w * atdk_vdw_E(flex_sc-rigid_prot) 
!                     + ROTA_w * ROTA(flex_sc-rigid_prot)
!   2. calc energy using pre-constructed table
!
!-------------------------------------------------------------------------------
MODULE HYBRID_E
!-------------------------------------------------------------------------------

use globals
use logger, only: log_p, terminate_with_error
use convert_res_name, only: convert_to_stdres
use rotamer, only: get_rotamer_index, place_rotamer, &
                   identify_rotamer_state_single
!
use energy_vars
use pair_list, only: pairlist_update
use energy_utils, only: protein_to_R, update_R_for_sc
!
use autodock3, only: calc_atdk3_prot_vdw
use ROTA, only: calc_ROTA_score

implicit none
save
private

integer, parameter :: max_rotamer_state = 81

!-------------------------------------------------------------------------------
type T_prot_type
!-------------------------------------------------------------------------------
real(dp) :: E(max_rotamer_state)
!-------------------------------------------------------------------------------
end type T_prot_type
!-------------------------------------------------------------------------------
type T_prot_twobody_type
!-------------------------------------------------------------------------------
real(dp) :: E(max_rotamer_state,max_rotamer_state)
!-------------------------------------------------------------------------------
end type T_prot_twobody_type
!-------------------------------------------------------------------------------

type(T_prot_type), allocatable :: T_prot(:)
integer, allocatable :: i_T_prot(:)  ! i_T_prot(i_usc) = i_res

type(T_prot_twobody_type), allocatable :: T_prot_twobody(:,:)

public :: construct_prot_table
public :: finalize_prot_table
public :: calc_prot_E_using_table

CONTAINS
!===============================================================================
! Construct flex_sc-rigid_prot interxn energy table
!===============================================================================
subroutine construct_prot_table(protein, is_usc)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
logical, intent(in) :: is_usc(:)

allocate(T_prot(n_usc))
allocate(i_T_prot(n_usc))

allocate(T_prot_twobody(n_usc,n_usc))

call construct_onebody_table(protein, is_usc, T_prot, i_T_prot)
call construct_twobody_table(protein, i_T_prot, T_prot_twobody)

call protein_to_R(protein)

end subroutine construct_prot_table
!-------------------------------------------------------------------------------
subroutine construct_onebody_table(protein, is_usc, T_prot, i_T_prot)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
logical, intent(in) :: is_usc(:)
type(T_prot_type), intent(inout) :: T_prot(:)
integer, intent(inout) :: i_T_prot(:)
!
type(molecule_type) :: tmp_prot
integer :: i_res, i_usc
integer :: i_rot, i_rot_start, i_rot_end
real(dp) :: ROTA_E, vdw_E, g(3,tn%atom)
logical :: appl_pair(protein%n_res,protein%n_res)

tmp_prot = protein
i_usc = 0
do i_res = 1, protein%n_res
    if (.not. is_usc(i_res)) cycle
    call set_appl_pair_for_onebody(i_res, is_usc, appl_pair, protein%n_res)

    i_usc = i_usc + 1
    i_T_prot(i_usc) = i_res

    call get_rotamer_index(protein%residue(i_res), i_rot_start, i_rot_end)

    do i_rot = i_rot_start, i_rot_end
        call place_rotamer(tmp_prot, i_res, i_rot)
        call update_R_for_sc(tmp_prot%residue(i_res), i_res)
        call pairlist_update()
        call calc_ROTA_score(ROTA_E, g, appl_pair, .false.)
        call calc_atdk3_prot_vdw(vdw_E, g, appl_pair, .false.)
        T_prot(i_usc)%E(i_rot-i_rot_start+1) = rota_w*ROTA_E + prot_w*vdw_E
    end do
end do

end subroutine construct_onebody_table
!-------------------------------------------------------------------------------
subroutine construct_twobody_table(protein, i_T_prot, T_prot_twobody)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
integer, intent(in) :: i_T_prot(:)
type(T_prot_twobody_type), intent(inout) :: T_prot_twobody(:,:)
!
type(molecule_type) :: tmp_prot
integer :: i_usc, j_usc, i_res, j_res
integer :: i_rot, i_rot_start, i_rot_end
integer :: j_rot, j_rot_start, j_rot_end
integer :: i_pair, ii_rot, jj_rot
logical :: appl_pair(protein%n_res,protein%n_res)
real(dp) :: ROTA_E, vdw_E, g(3,tn%atom)

tmp_prot = protein
do i_usc = 1, n_usc-1
    i_res = i_T_prot(i_usc)
    call get_rotamer_index(protein%residue(i_res), i_rot_start, i_rot_end)

    do j_usc = i_usc + 1, n_usc
        j_res = i_T_prot(j_usc)
        call get_rotamer_index(protein%residue(j_res), j_rot_start, j_rot_end)

        i_pair = 0
        do i_rot = i_rot_start, i_rot_end
            ii_rot = i_rot - i_rot_start + 1
            call place_rotamer(tmp_prot, i_res, i_rot)
            call update_R_for_sc(tmp_prot%residue(i_res), i_res)
            do j_rot = j_rot_start, j_rot_end
                jj_rot = j_rot - j_rot_start + 1
                i_pair = i_pair + 1
                call place_rotamer(tmp_prot, j_res, j_rot)
                call update_R_for_sc(tmp_prot%residue(j_res), j_res)
                call pairlist_update()
                call set_appl_pair_for_twobody(i_res, j_res, appl_pair)
                call calc_ROTA_score(ROTA_E, g, appl_pair, .false.)
                call calc_atdk3_prot_vdw(vdw_E, g, appl_pair,  .false.)
                T_prot_twobody(j_usc, i_usc)%E(jj_rot, ii_rot) = rota_w*ROTA_E &
                                                               + prot_w*vdw_E
            end do
        end do
    end do
end do

end subroutine construct_twobody_table
!-------------------------------------------------------------------------------
subroutine finalize_prot_table()
!-------------------------------------------------------------------------------
if (allocated(T_prot)) then
    deallocate(T_prot, i_T_prot)
    deallocate(T_prot_twobody)
end if

end subroutine finalize_prot_table
!-------------------------------------------------------------------------------
subroutine set_appl_pair_for_onebody(i_res, is_usc, appl_pair, n_res)
!-------------------------------------------------------------------------------
logical, intent(in) :: is_usc(:)
logical, intent(out) :: appl_pair(:,:)
integer, intent(in) :: i_res, n_res
integer :: j_res

appl_pair(:,:) = .false.

if (is_usc(i_res)) then
    do j_res = 1, n_res
        if (.not. is_usc(j_res)) then
            appl_pair(i_res, j_res) = .true.
            appl_pair(j_res, i_res) = .true.
        end if
    end do
end if

end subroutine set_appl_pair_for_onebody
!-------------------------------------------------------------------------------
subroutine set_appl_pair_for_twobody(i_res, j_res, appl_pair)
!-------------------------------------------------------------------------------
integer, intent(in) :: i_res, j_res
logical, intent(out) :: appl_pair(:,:)

appl_pair(:,:) = .false.
appl_pair(i_res, j_res) = .true.
appl_pair(j_res, i_res) = .true.

end subroutine set_appl_pair_for_twobody
!-------------------------------------------------------------------------------
subroutine calc_prot_E_using_table(rot_list, prot_E)
!-------------------------------------------------------------------------------
! Calculate ROTA score using pre-constructed energy table.
!-------------------------------------------------------------------------------
integer, intent(in) :: rot_list(:)
real(dp), intent(out) :: prot_E
integer :: i_rot, j_rot
integer :: i_usc, j_usc

prot_E = 0.0d0

! calc onebody energy
do i_usc = 1, n_usc
    i_rot = rot_list(i_usc)
    prot_E = prot_E + T_prot(i_usc)%E(i_rot)
end do

! calc twobody energy
do i_usc = 1, n_usc - 1
    i_rot = rot_list(i_usc)
    do j_usc = i_usc + 1, n_usc
        j_rot = rot_list(j_usc)
        prot_E = prot_E + T_prot_twobody(j_usc,i_usc)%E(j_rot, i_rot)
    end do
end do

end subroutine calc_prot_E_using_table
!-------------------------------------------------------------------------------
END MODULE HYBRID_E
!-------------------------------------------------------------------------------
