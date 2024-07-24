!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/sampling/optimizers/ligdock_csa/ligdock_csa_dist.f90
!
! Description:
!
!-------------------------------------------------------------------------------
MODULE LIGDOCK_CSA_DIST
!-------------------------------------------------------------------------------

use globals
!
use rmsd, only: calc_rmsd_simple, calc_rmsd_PM, prep_RMSD_PM_calculation
use rotamer, only: place_rotamer
use xlogp_Score_m, only: get_xlogP_atom_types
!
use ligdock_csa_vars

implicit none
save
public

character(len=20), allocatable:: atom_types(:)
integer :: n_type
integer, allocatable :: n_atm_s(:), atm_id_s(:,:)

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_distance(ligand)
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand

allocate(atom_types(ligand%n_atm))
allocate(n_atm_s(ligand%n_atm))
allocate(atm_id_s(ligand%n_atm, ligand%n_atm))

call get_xlogP_atom_types(ref_lig(ligand%lig_type), atom_types)
call prep_RMSD_PM_calculation(atom_types, ligand%n_atm, n_type, n_atm_s, atm_id_s)

end subroutine initialize_distance
!-------------------------------------------------------------------------------
subroutine finalize_distance()
!-------------------------------------------------------------------------------
deallocate(atom_types)
deallocate(n_atm_s)
deallocate(atm_id_s)

end subroutine finalize_distance
!-------------------------------------------------------------------------------
subroutine calc_dist_btw_bank(bank1, bank2, n_lig_atm, dist)
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(in) :: bank1, bank2
integer, intent(in) :: n_lig_atm
real(dp), intent(out) :: dist
!
real(dp) :: prot_dist, sub_dist
integer :: i_usc, n_atm

dist = 0.0

call calc_rmsd_simple(bank1%R(:, 1:n_lig_atm), bank2%R(:, 1:n_lig_atm), &
                      n_lig_atm, dist)
!call calc_rmsd_PM(bank1%R(:,1:n_lig_atm), bank2%R(:,1:n_lig_atm), &
!                  n_type, n_atm_s, atm_id_s, dist)
if (n_usc > 0) then
    prot_dist = 0.0d0
    do i_usc = 1, n_usc
        n_atm = bank1%flex_sc(i_usc)%n_atm
        call calc_rmsd_simple(bank1%flex_sc(i_usc)%R(:, 1:n_atm), &
                              bank2%flex_sc(i_usc)%R(:, 1:n_atm), &
                              n_atm, sub_dist)
        prot_dist = prot_dist + sub_dist
    end do
    prot_dist = prot_dist / dble(n_usc)
    dist = dist + 0.5*prot_dist
!    dist = dist + prot_dist
end if

end subroutine calc_dist_btw_bank
!-------------------------------------------------------------------------------
subroutine update_distance(bank, n_bank, n_lig_atm, Dij)
!-------------------------------------------------------------------------------
type(ligdock_bank_type), intent(in) :: bank(:)
integer, intent(in) :: n_bank, n_lig_atm
real(dp), intent(out) :: Dij(:,:)
integer :: i_bank, j_bank, i_dist
real(dp) :: dist, sum_dist

sum_dist = 0.0
Dij(:,:) = 0.0d0
i_dist = 0

do i_bank = 1, n_bank - 1
    do j_bank = i_bank + 1, n_bank
        call calc_dist_btw_bank(bank(i_bank), bank(j_bank), n_lig_atm, dist)
        Dij(i_bank, j_bank) = dist
        Dij(j_bank, i_bank) = dist
        sum_dist = sum_dist + dist
        i_dist = i_dist + 1
    end do
end do

D_ave = sum_dist/dble(i_dist)

end subroutine update_distance
!-------------------------------------------------------------------------------
subroutine setup_distance(n_bank_add, n_bank, Dij)
!-------------------------------------------------------------------------------
integer, intent(in) :: n_bank_add, n_bank
real(dp), intent(in) :: Dij(:,:)
!
integer :: i_bank, j_bank, i_dist
real(dp) :: nst, D_ave_of_new_conf, sum_dist

i_dist = 0
sum_dist = 0.0d0
nst = dble(n_opt_to_D_min/n_new_conf)
 
do i_bank = n_bank - n_bank_add + 1, n_bank - 1
    do j_bank = i_bank + 1, n_bank
        i_dist = i_dist + 1
        sum_dist = sum_dist + Dij(i_bank, j_bank)
    end do
end do
  
D_ave_of_new_conf = sum_dist/dble(i_dist)
  
D_cut = D_ave_of_new_conf/factor_init_D_cut            ! The initial cutoff
D_min = D_ave_of_new_conf/factor_min_D_cut         ! the minimum cutoff
  
! The cutoff goes to minimum after n_opt_to_D_min/n_new_conf minimizations
xctdif = (D_cut/D_min)**(-1.0d0/nst) 
  
end subroutine setup_distance
!-------------------------------------------------------------------------------
END MODULE LIGDOCK_CSA_DIST
!-------------------------------------------------------------------------------
