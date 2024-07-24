!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/modeling/knowledge/ramachandran_score.f90
!
! Description:
!-------------------------------------------------------------------------------
MODULE RAMACHANDRAN_SCORE
!-------------------------------------------------------------------------------
use globals
use logger,        only: log_p, terminate_with_error
use mathfunctions, only: cross, bicubic_interpolation, multiple_binormal, bound_ang
use geometry,      only: calc_tor_and_grad
use ramachandran,  only: use_rama, initialize_ramachandran, finalize_ramachandran, &
                         protein_rama, protein_ss_prob, delta_angle_rad, &
                         fill_inter_neigh
!
use energy_vars,   only: use_rama_score, R

implicit none
save
private

public :: initialize_rama_score
public :: finalize_rama_score
!
public :: calc_rama_score
!
public :: calc_phipsi_using_R
public :: get_phipsi_using_R

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_rama_score(protein)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein

use_rama = .true.
call initialize_ramachandran(protein)

end subroutine initialize_rama_score
!-------------------------------------------------------------------------------
subroutine finalize_rama_score()
!-------------------------------------------------------------------------------
call finalize_ramachandran()

end subroutine finalize_rama_score
!-------------------------------------------------------------------------------
subroutine calc_rama_score(f, g, use_res, calc_g)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: f, g(3,tn%atom)
logical, intent(in) :: use_res(tn%stdres), calc_g

logical :: rama_use_res(tn%stdres)
integer :: i_res, i_atm, i_ref, i_phi, i_psi, atmno
real(dp) :: phipsi(2), dphipsi(3,4,2)
real(dp) :: ff(4,4), x(2), arg(3)

f = 0.0d0
g(:,:) = 0.0d0

rama_use_res(2:tn%stdres-1) = (use_res(1:tn%stdres-2) .or. use_res(2:tn%stdres-1)) &
                               .or. use_res(3:tn%stdres)

do i_res = 2, tn%stdres-1
    if (.not. rama_use_res(i_res)) cycle
    i_ref = res_index(i_res)%ref_res_no
    if (ref_res(i_ref)%ter_type == 'N' .or. &
        ref_res(i_ref)%ter_type == 'C') cycle

    call calc_phipsi_using_R(i_res, phipsi, dphipsi, calc_g)
    i_phi = min(35, int((phipsi(1)+pi)/delta_angle_rad))
    i_psi = min(35, int((phipsi(2)+pi)/delta_angle_rad))
    !
    call fill_inter_neigh(i_res, i_phi, i_psi, ff(:,:))
    x(1) = phipsi(1) + pi - dble(i_phi)*delta_angle_rad
    x(2) = phipsi(2) + pi - dble(i_psi)*delta_angle_rad
    !
    arg(:) = bicubic_interpolation(x(1), x(2), delta_angle_rad, delta_angle_rad, ff)
    !
    f = f + arg(1)
    if (calc_g) then
        do i_atm = 1, 4
            atmno = res_index(i_res)%phi_atmid(i_atm,1)
            g(:,atmno) = g(:,atmno) + arg(2)*dphipsi(:,i_atm,1)
            atmno = res_index(i_res)%psi_atmid(i_atm,1)
            g(:,atmno) = g(:,atmno) + arg(3)*dphipsi(:,i_atm,2)
        end do
    end if
end do

end subroutine calc_rama_score
!-------------------------------------------------------------------------------
subroutine calc_phipsi_using_R(i_res, phipsi, dphipsi, calc_g)
!-------------------------------------------------------------------------------
integer, intent(in) :: i_res
real(dp), intent(out) :: phipsi(2), dphipsi(3,4,2)
logical, intent(in) :: calc_g
real(dp) :: ra(3,4,2)
integer :: atmid(4,2)

atmid(1:4,1) = res_index(i_res)%phi_atmid(:,1)
atmid(1:4,2) = res_index(i_res)%psi_atmid(:,1)

ra(:,:,1) = R(:,atmid(:,1))
ra(:,:,2) = R(:,atmid(:,2))

call calc_tor_and_grad(ra(:,:,1), phipsi(1), dphipsi(:,:,1), calc_g)
call calc_tor_and_grad(ra(:,:,2), phipsi(2), dphipsi(:,:,2), calc_g)

end subroutine calc_phipsi_using_R
!-------------------------------------------------------------------------------
subroutine get_phipsi_using_R(phi, psi)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: phi(tn%stdres), psi(tn%stdres)
real(dp) :: phipsi(2), dphipsi(3,4,2)
integer :: i_res, ref_res_no

phi(:) = 0.0d0
psi(:) = 0.0d0

do i_res = 2, tn%stdres-1
    ref_res_no = res_index(i_res)%ref_res_no
    if (ref_res(ref_res_no)%ter_type == 'N' .or. ref_res(ref_res_no)%ter_type == 'C') cycle
    call calc_phipsi_using_R(i_res, phipsi, dphipsi, .false.)
    phi(i_res) = phipsi(1)
    psi(i_res) = phipsi(2)
end do

end subroutine get_phipsi_using_R
!-------------------------------------------------------------------------------
END MODULE RAMACHANDRAN_SCORE
!-------------------------------------------------------------------------------
