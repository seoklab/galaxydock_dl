!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/modeling/knowledge/rotamer_score.f90
!
! Description:
!
!-------------------------------------------------------------------------------
MODULE ROTAMER_SCORE
!-------------------------------------------------------------------------------

use globals
use logger,        only: log_p
use rotamer,       only: rotamer_index, rotamer_record, is_symmetric_chi, &
                         rotamer_index_type, rotamer_record_type
use geometry,      only: calc_tor_and_grad, calc_tor_angle
use energy_vars,   only: R
use mathfunctions, only: multiple_multinomial, bound_ang

implicit none
save
private

public :: calc_rotamer_score

CONTAINS
!-------------------------------------------------------------------------------
subroutine calc_rotamer_score(f, g, use_res, calc_g)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: f, g(3,tn%atom)
logical, intent(in) :: use_res(tn%residue)
logical, intent(in) :: calc_g
integer :: i_res, i_ref, i_rot_res, i_atm, ia, i_chi, n_chi
integer :: phi_atmid(4), psi_atmid(4), i_phi, i_psi
integer :: i_rot_idx(2), n_rotamer, i_rot, rot_no
real(dp) :: phi, psi, v_phi(3,3), v_psi(3,3)
real(dp) :: Rchi(3,4), chi(4), dcdr(3,4,4)
integer, parameter :: max_rotamer = 81
real(dp) :: prob(max_rotamer), chi0(4,max_rotamer), sigma(4,max_rotamer)
real(dp) :: fv(5), dfdc(4)

f = 0.0d0
g(:,:) = 0.0d0

do i_res = 1, tn%stdres
    n_chi = res_index(i_res)%n_chi
    if ((.not. use_res(i_res)) .or. n_chi == 0) cycle
    i_ref = res_index(i_res)%ref_res_no
    i_rot_res = ref_res(i_ref)%i_rot_res
    !
    ! Calculate Ramachandran angles
    phi_atmid(:) = res_index(i_res)%phi_atmid(:,1)
    psi_atmid(:) = res_index(i_res)%psi_atmid(:,1)
    do i_atm = 1, 3
        v_phi(:,i_atm) = R(:,phi_atmid(i_atm+1)) - R(:,phi_atmid(i_atm))
        v_psi(:,i_atm) = R(:,psi_atmid(i_atm+1)) - R(:,psi_atmid(i_atm))
    end do
    call calc_tor_angle(v_phi(:,1), v_phi(:,2), v_phi(:,3), phi)
    call calc_tor_angle(v_psi(:,1), v_psi(:,2), v_psi(:,3), psi)
    !
    ! Identify rotamer indices
    i_phi = anint(rad2deg*phi/10.0d0)
    i_psi = anint(rad2deg*psi/10.0d0)
    i_rot_idx(1) = rotamer_index(i_rot_res)%bbdep_id(1,i_phi,i_psi)
    i_rot_idx(2) = rotamer_index(i_rot_res)%bbdep_id(2,i_phi,i_psi)
    n_rotamer = i_rot_idx(2) - i_rot_idx(1) + 1
    !
    ! Calculate Chi-angles
    chi(:) = 0.0d0
    dcdr(:,:,:) = 0.0d0
    do i_chi = 1, n_chi
        Rchi(:,:) = 0.0d0
        do i_atm = 1, 4
            Rchi(:,i_atm) = R(:,res_index(i_res)%chi_atmid(i_atm,i_chi,1))
        end do
        call calc_tor_and_grad(Rchi(:,:), chi(i_chi), dcdr(:,:,i_chi), calc_g)
    end do
    !
    ! Get rotamer records
    do i_rot = i_rot_idx(1), i_rot_idx(2)
        rot_no = i_rot - i_rot_idx(1) + 1
        !
        prob(rot_no)    = rotamer_record(i_rot)%prob
        chi0(:,rot_no)  = rotamer_record(i_rot)%chi(:)
        sigma(:,rot_no) = rotamer_record(i_rot)%sigma(:)
        !
        do i_chi = 1, n_chi
            if (is_symmetric_chi(i_ref, i_chi)) then
                if (bound_ang(chi(i_chi)-chi0(i_chi,rot_no)) > pi_half) then
                    chi0(i_chi,rot_no) = chi0(i_chi,rot_no) + pi
                else if (bound_ang(chi(i_chi)-chi0(i_chi,rot_no)) < -pi_half) then
                    chi0(i_chi,rot_no) = chi0(i_chi,rot_no) - pi
                end if
            end if
        end do
    end do
    !
    ! Get function and derivative
    fv = multiple_multinomial(n_rotamer, n_chi, chi, chi0, sigma, prob, &
                              .true., .true., calc_g)
    f = f + fv(1)
    dfdc(1:4) = fv(2:5)
    if (calc_g) then
        do i_chi = 1, n_chi
            do i_atm = 1, 4
                ia = res_index(i_res)%chi_atmid(i_atm,i_chi,1)
                g(:,ia) = g(:,ia) + dfdc(i_chi)*dcdr(:,i_atm,i_chi)
            end do
        end do
    end if
end do

end subroutine calc_rotamer_score
!-------------------------------------------------------------------------------
END MODULE ROTAMER_SCORE
!-------------------------------------------------------------------------------
