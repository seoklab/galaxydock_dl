!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/ligdock/empirical/PLP_tor.f90
!
! Description:
!
!-------------------------------------------------------------------------------
MODULE PLP_TOR
!-------------------------------------------------------------------------------

use globals
use geometry, only: calc_tor_angle, calc_tor_and_grad
!
use energy_vars

implicit none
save
private

integer, allocatable :: lig_tor_type(:,:)
integer, allocatable :: lig_tor_atm(:,:)
integer :: num_lig_tor

public :: initialize_PLP_tor
public :: finalize_PLP_tor
public :: calc_PLP_tor

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_PLP_tor(ligand)
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand

num_lig_tor = ligand%n_br - 1
allocate(lig_tor_type(2,num_lig_tor))
allocate(lig_tor_atm(4,num_lig_tor))
call setup_PLP_tor(ligand, lig_tor_type, lig_tor_atm)

end subroutine initialize_PLP_tor
!-------------------------------------------------------------------------------
subroutine finalize_PLP_tor()
!-------------------------------------------------------------------------------
deallocate(lig_tor_type)
deallocate(lig_tor_atm)

end subroutine finalize_PLP_tor
!-------------------------------------------------------------------------------
subroutine setup_PLP_tor(ligand, lig_tor_type, lig_tor_atm)
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand
integer, intent(out) :: lig_tor_type(:,:), lig_tor_atm(:,:)
integer :: i, i_tor, tor_idx, atm_idx, atm_1, atm_2
character(len=4) :: hybrid_1, hybrid_2

do i_tor = 1, ligand%n_br - 1
    ! get atom index to define this torsion
    tor_idx = ligand%i_rot_tor(i_tor)
    do i = 1, 4
        atm_idx = ref_lig(ligand%lig_type)%dih(i+1, tor_idx)
        lig_tor_atm(i,i_tor) = ii_L(atm_idx, ligand%lig_no)
    end do

    ! assign torsion type
    lig_tor_type(:,i_tor) = 1
    
    atm_1 = ligand%bridge(1, i_tor+1)
    atm_2 = ligand%bridge(2, i_tor+1)

    hybrid_1 = ref_lig(ligand%lig_type)%mol2_type(atm_1)(3:6)
    hybrid_2 = ref_lig(ligand%lig_type)%mol2_type(atm_2)(3:6)

    if (hybrid_1(1:1) == '2') then
        lig_tor_type(1,i_tor) = 2
    else if(hybrid_1(1:1) == '3') then
        lig_tor_type(1,i_tor) = 3
    else if(hybrid_1(1:2) == 'ar') then
        lig_tor_type(1,i_tor) = 2
    else if(hybrid_1(1:2) == 'am') then
        lig_tor_type(1,i_tor) = 2
    else if(hybrid_1(1:2) == 'pl') then
        lig_tor_type(1,i_tor) = 2
    else if(hybrid_1(1:1) == '4') then
        lig_tor_type(1,i_tor) = 3
    else if(hybrid_1(1:2) == 'o2') then
        lig_tor_type(1,i_tor) = 2
    end if

    if(hybrid_2(1:1) == '2') then
        lig_tor_type(2,i_tor) = 2
    else if(hybrid_2(1:1) == '3') then
        lig_tor_type(2,i_tor) = 3
    else if(hybrid_2(1:2) == 'ar') then
        lig_tor_type(2,i_tor) = 2
    else if(hybrid_2(1:2) == 'am') then
        lig_tor_type(2,i_tor) = 2
    else if(hybrid_2(1:2) == 'pl') then
        lig_tor_type(2,i_tor) = 2
    else if(hybrid_2(1:1) == '4') then
        lig_tor_type(2,i_tor) = 3
    else if(hybrid_2(1:2) == 'o2') then
        lig_tor_type(2,i_tor) = 2
    end if
end do

end subroutine setup_PLP_tor
!-------------------------------------------------------------------------------
subroutine calc_PLP_tor(E_tor, g, calc_g)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: E_tor, g(:,:)
logical, intent(in) :: calc_g
integer :: i_tor, i_atm, atm_idx, L_atm_idx
real(dp) :: tor_ang, i_tor_energy, ra(3,4), dadr(3,4), dEda

E_tor = 0.0d0
if (calc_g) g=0.0d0

do i_tor = 1, num_lig_tor
    ! calc torsion angles
    do i_atm = 1, 4
        ra(:,i_atm) = R(:,lig_tor_atm(i_atm, i_tor))
    end do
    call calc_tor_and_grad(ra, tor_ang, dadr, calc_g)

    ! calc PLP_tor score
    if(lig_tor_type(1, i_tor) == 2 .and. lig_tor_type(2, i_tor) == 3) then
        call sp2_sp3(tor_ang, i_tor_energy, dEda, calc_g)
        E_tor = E_tor + i_tor_energy
      
    else if(lig_tor_type(1, i_tor) == 3 .and. lig_tor_type(2, i_tor) == 2) then
        call sp2_sp3(tor_ang, i_tor_energy, dEda, calc_g)
        E_tor = E_tor + i_tor_energy
              
    else if(lig_tor_type(1, i_tor) == 3 .and. lig_tor_type(2, i_tor) == 3) then
        call sp3_sp3(tor_ang, i_tor_energy, dEda, calc_g)
        E_tor = E_tor + i_tor_energy
    end if
    if (calc_g) then
        do i_atm = 1, 4
            atm_idx = lig_tor_atm(i_atm, i_tor)
            L_atm_idx = i_R(2,atm_idx)
            g(:,L_atm_idx) = g(:,L_atm_idx) + dEda*dadr(:,i_atm)
        end do
    end if
end do

!-------------------------------------------------------------------------------
end subroutine calc_PLP_tor
!-------------------------------------------------------------------------------
subroutine sp3_sp3(tor_ang, E_tor, dEda, calc_g)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: tor_ang
real(dp), intent(out) :: E_tor, dEda
logical, intent(in) :: calc_g

E_tor = 1.5*(1.0 + cos(3.0d0*tor_ang - pi))
if (calc_g) then
    dEda = -1.5d0 * 3.0d0 * sin(3.0d0*tor_ang-pi)
end if

!-------------------------------------------------------------------------------
end subroutine sp3_sp3
!-------------------------------------------------------------------------------
subroutine sp2_sp3(tor_ang, E_tor, dEda, calc_g)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: tor_ang
real(dp), intent(out) :: E_tor, dEda
logical, intent(in) :: calc_g

E_tor = 0.75*(1.0 + cos(6.0d0*tor_ang))
if (calc_g) then
    dEda = -0.75d0 * 6.0d0 * sin(6.0d0*tor_ang)
end if

!-------------------------------------------------------------------------------
end subroutine sp2_sp3
!-------------------------------------------------------------------------------
END MODULE PLP_TOR
!-------------------------------------------------------------------------------
