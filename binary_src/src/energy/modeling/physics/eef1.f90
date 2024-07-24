!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/modeling/physics/eef1.f90
!
! Description: EEF1 solvation energy
!
! TODO add comment
! Reference
!
!-------------------------------------------------------------------------------
MODULE EEF1

use globals
use energy_vars, only: R, i_R, i_P, appl_res, Roff_sqr, max_neigh, update_reff, &
                       max_energy
use cutoff_methods, only: nonbonded_switch_function
use molecular_mechanics, only: LJ_type, QQ_type

implicit none
save
private

integer, allocatable :: neighno(:,:), n_neigh(:) ! Neighbor atom no. index (w/o excluded)
real(dp), allocatable :: rsqr_neigh(:,:)
! Parameters for EEF1
integer :: n_grp_eef1
integer, allocatable :: eef1_type(:)
real(dp), allocatable :: vol(:) ! Atomic volume for atom type
real(dp) :: dg_free(17), dg_ref(17) ! delta-free energy
real(dp), allocatable :: i_lambda(:)! inverse lambda
real(dp), allocatable :: Gpair(:)

public :: initialize_eef1
public :: calc_eef1
public :: finalize_eef1

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_eef1(molecule)
!-------------------------------------------------------------------------------
! TODO comment
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: molecule

allocate(n_neigh(tn%atom), neighno(max_neigh,tn%atom))
allocate(rsqr_neigh(max_neigh,tn%atom))

n_neigh(:) = 0
neighno(:,:) = 0

call setup_eef1_parm(molecule)


end subroutine initialize_eef1
!-------------------------------------------------------------------------------
subroutine setup_eef1_parm(molecule)
!-------------------------------------------------------------------------------
! TODO comment
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: molecule
character(len=4) :: atm_name, res_name
integer :: i_atm, resno, atmno, grp
integer :: ref_res_no

if (top_type == 'allh') then
    n_grp_eef1 = 55
else if (top_type == 'eef1') then
    n_grp_eef1 = 31
else if (top_type == 'coarse') then
    n_grp_eef1 = 34
end if

allocate(i_lambda(tn%atom))
allocate(vol(17))
allocate(eef1_type(n_grp_eef1))
allocate(Gpair(tn%atom))

!TODO not used?
!Roff_sqr_pair = Roff_sqr

! eef1_type: Map LJ_type into parameter index defined below
! Atom orderes are consistent with parameters_charmm22ph.in
! This is not consistent with original EEF1 because we are using charmm22.
! This might be modified into mapping by Lazaridis for charmm22 parms
if (top_type == 'allh') then
    eef1_type(:) = &
          (/ 0, 1, 1, 1, 1, 0, 0, 0, 0, 6, &
          4, 6, 6, 6, 0, 0, 0, 0, 0, 0, & !<-- 11th location(CT) is most problematic
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
          0, 0, 0, 7, 7, 7, 0, 9,10, 0, &
         14, 0,13, 0,15, 0,16,17, 0, 0, &
          0, 0, 0, 0, 0/)
else if (top_type == 'eef1') then
    eef1_type(:) = &
          (/ 0, 0, 1, 6, 3,   4, 5, 6, 6, 1, &
          0, 3, 4, 4, 1,  12, 7, 8,10, 7, &
          9,10,11, 7, 0,  14,15,13,16, 0/)
else if (top_type == 'coarse') then
    eef1_type(:) = &
          (/ 1,  2,  3,  4,  0,  0,  5,  0,  0,  6, &
             7,  8,  9,  0, 10, 11, 12,  0, 13, 14, &
            15,  0,  0, 16,  0,  0,  0, 17,  0,  0, &
             0,  0,  0,  0 /)
end if

! Following order in original paper Table I
if (top_type /= 'coarse') then
    vol(:) = &
          (/ 14.7d0,  8.3d0, 23.7d0, 22.4d0, 30.0d0, &
             18.4d0,  4.4d0,  4.4d0, 11.2d0, 11.2d0, &
             11.2d0,  0.0d0, 10.8d0, 10.8d0, 10.8d0, &
             14.7d0, 21.4d0 /)

    dg_ref(:) = &
          (/ 0.000d0, -0.890d0, -0.187d0,  0.372d0,   1.089d0, &
             0.057d0, -5.950d0, -3.820d0, -5.450d0, -20.000d0, &
            -10.000d0, -1.000d0, -5.920d0, -5.330d0, -10.000d0, &
            -3.240d0, -2.050d0 /)
     
    dg_free(:) = &
          (/ 0.00d0, -1.40d0, -0.25d0,  0.52d0,   1.50d0, &
             0.08d0, -8.90d0, -4.00d0, -7.80d0, -20.00d0, &
            -10.00d0, -1.55d0, -6.70d0, -5.85d0, -10.00d0, &
            -4.10d0, -2.70d0 /)

    i_lambda(:) = 1.0d0/3.5d0

    do i_atm = 1, tn%atom
        grp = eef1_type(LJ_type(i_atm))
        
        resno = i_R(1,i_atm)
        atmno = i_R(2,i_atm)
        ref_res_no = molecule%residue(resno)%res_type
        
        atm_name = ref_res(ref_res_no)%atom_name(atmno)
        res_name = trim(ref_res(ref_res_no)%res_name)

        if (grp == 0) cycle
     
        !Correction for pseudo-ionic sidechains
        if (molecule%residue(resno)%ter_type == 'N' .and. atm_name == ' N  ') then
            i_lambda(i_atm) = 1.0d0/6.0d0
        elseif (molecule%residue(resno)%ter_type == 'C' .and. (atm_name == ' C  ' .or. &
             atm_name == ' O  ' .or. atm_name == ' OXT')) then
            i_lambda(i_atm) = 1.0d0/6.0d0
        elseif (res_name == 'ARG' .and. (atm_name == ' CD ' .or. atm_name == ' NE ' .or. &
             atm_name == ' CZ ' .or. atm_name == 'NH1' .or. atm_name == 'NH2')) then
            i_lambda(i_atm) = 1.0d0/6.0d0
        elseif (res_name == 'LYS' .and. (atm_name == ' CE ' .or. atm_name == ' NZ ')) then
            i_lambda(i_atm) = 1.0d0/6.0d0
        elseif (res_name == 'ASP' .and. (atm_name == ' CB ' .or. atm_name == ' CG ' .or. &
             atm_name == ' OD1' .or. atm_name == 'OD2')) then
            i_lambda(i_atm) = 1.0d0/6.0d0
        elseif (res_name == 'GLU' .and. (atm_name == ' CG ' .or. atm_name == ' CD ' .or. &
             atm_name == ' OE1' .or. atm_name == ' OE2')) then
            i_lambda(i_atm) = 1.0d0/6.0d0
        end if

    end do

else
    vol(:) = &
          (/ 4.4d0, 25.5d0, 23.7d0, 53.7d0, 48.9d0, &
            23.7d0, 78.4d0, 89.5d0, 37.1d0, 67.2d0, &
            75.5d0,106.1d0, 35.4d0, 58.7d0, 81.1d0, &
            29.0d0, 46.6d0 /)

    dg_ref(:) = &
          (/  -5.950d0, -2.665d0, -0.187d0,  0.451d0, -2.224d0, &
                -0.187d0, -2.221d0, -0.352d0, -1.434d0,  0.372d0, &
               0.518d0,  0.591d0, -8.241d0, -4.907d0, -3.851d0, &
                -2.427d0,  0.151d0 /)

    dg_free(:) = &
          (/  -8.90d0, -2.92d0, -0.25d0,  0.62d0, -2.46d0, &
                -0.25d0, -2.11d0, -0.39d0, -1.79d0,  0.52d0, &
               0.72d0,  0.82d0, -8.21d0, -4.87d0, -3.79d0, &
                -3.24d0,  0.21d0 /)

    i_lambda(:) = 1.0d0/3.5d0
     
    do i_atm = 1, tn%atom
        grp = eef1_type(LJ_type(i_atm))

        if (grp == 0) cycle
        
        !if it is a charged atom (Qd, QaD, QaE)
        if (grp == 13 .or. grp == 14 .or. grp ==15) then
            i_lambda(i_atm) = 1.0d0/6.0d0
        end if
    end do

end if

end subroutine setup_eef1_parm
!-------------------------------------------------------------------------------
subroutine set_neighbor_pair(status)
!-------------------------------------------------------------------------------
! Search pairlist for pair interaction
!-------------------------------------------------------------------------------
logical, intent(out) :: status
integer :: i_atm, i_j, j_atm, i_start
real(dp) :: dij

status = .true.
n_neigh(:) = 0

do i_atm = 1, tn%atom
    i_start = i_P(i_atm)%pair_end_index(2)+1
    do i_j = i_start, i_P(i_atm)%n_pair
        j_atm = i_P(i_atm)%i_pair(i_j)
        if (eef1_type(LJ_type(j_atm)) == 0) cycle
        dij = i_P(i_atm)%d(i_j)
        n_neigh(i_atm) = n_neigh(i_atm) + 1
        neighno(n_neigh(i_atm),i_atm) = j_atm
        rsqr_neigh(n_neigh(i_atm),i_atm) = dij*dij
    end do
end do

end subroutine set_neighbor_pair
!-------------------------------------------------------------------------------
subroutine finalize_eef1()
!-------------------------------------------------------------------------------
! TODO comment
!-------------------------------------------------------------------------------

deallocate(vol)
deallocate(i_lambda)
deallocate(eef1_type)
deallocate(Gpair)
deallocate(n_neigh)
deallocate(neighno)
deallocate(rsqr_neigh)

end subroutine finalize_eef1
!-------------------------------------------------------------------------------
subroutine calc_eef1(f, g, appl_res, calc_g)
!-------------------------------------------------------------------------------
! TODO comment
!-------------------------------------------------------------------------------
real(dp), intent(out) :: f, g(3,tn%atom)
logical, intent(in) :: appl_res(tn%residue), calc_g
logical :: status
logical :: non0qi, non0qj
real(dp) :: Ri(3), Rij(3), r_sqr, i_rabs
real(dp) :: qi, dGdrij
real(dp) :: argv(3), df_s(3)
real(dp) :: f_tmp
real(dp) :: ftr, gtr
integer :: i_atm, j_atm, i_j

f = 0.0
g(:,:) = 0.0

if (update_reff) then
    call set_neighbor_pair(status)
    if (.not. status) then
        f = max_energy
        return
    end if
end if

!Self term
do i_atm = 1, tn%atom
    non0qi = eng_para%non0q(QQ_type(i_atm))
    qi = eng_para%charge(QQ_type(i_atm))
    
    if (eef1_type(LJ_type(i_atm)) == 0) cycle
        
    if (.not. appl_res(i_R(1,i_atm))) cycle
    f = f + dg_ref(eef1_type(LJ_type(i_atm)))
end do

!Pair term
do i_atm = 1, tn%atom
    if (eef1_type(LJ_type(i_atm)) == 0) cycle
    Ri(:) = R(:,i_atm)
    non0qi = eng_para%non0q(QQ_type(i_atm))
    Gpair(i_atm) = 0.0d0
    do i_j = 1, n_neigh(i_atm)
        j_atm = neighno(i_j,i_atm)
        non0qj = eng_para%non0q(QQ_type(j_atm))
        
        Rij(:) = R(:,j_atm) - Ri(:)
        r_sqr = rsqr_neigh(i_j,i_atm)
        if (r_sqr < small_real) cycle

        i_rabs = 1.0d0/sqrt(r_sqr)
        
        call nonbonded_switch_function(r_sqr, ftr, gtr, calc_g)
        
        call calc_eef1_pair(f_tmp, dGdrij, i_rabs, r_sqr, LJ_type(i_atm), LJ_type(j_atm), &
            i_lambda(i_atm), i_lambda(j_atm), calc_g)
        
        Gpair(i_atm) = Gpair(i_atm) + f_tmp*ftr

        if (calc_g) then
            df_s(:) = -f_tmp*gtr*Rij(:)
            argv(:) = ftr*dGdrij*Rij(:)
            g(:,i_atm) = g(:,i_atm) + argv(:) + df_s(:)
            g(:,j_atm) = g(:,j_atm) - argv(:) - df_s(:)
        end if
    end do ! Iter over j_atm 
    f = f + Gpair(i_atm)
end do ! Iter over i_atm

end subroutine calc_eef1
!-------------------------------------------------------------------------------
subroutine calc_eef1_pair(f, dfdrij, i_rabs, r_sqr, i_type, j_type, lambda_i, lambda_j, calc_g)
!-------------------------------------------------------------------------------
! TODO comment
!-------------------------------------------------------------------------------
integer, intent(in) :: i_type, j_type
real(dp), intent(in) :: i_rabs, r_sqr, lambda_i, lambda_j
logical, intent(in) :: calc_g
real(dp), intent(out) :: f, dfdrij
real(dp), parameter :: cexp_eef1 = 1.0d0/(pi*sqrt(pi))
real(dp) :: r_abs, xij, xji, tmp, f1, f2, r_i, r_j, dri, drj
real(dp) :: a, b, c, dfdrij1, dfdrij2, l2, r2
integer :: i, j

r_abs = 1.0d0/i_rabs
i = eef1_type(i_type) ! 17 types in EEF1
j = eef1_type(j_type) ! 17 types in EEF1
r_i = eng_para%LJ_para(2,i_type)
r_j = eng_para%LJ_para(2,j_type)
tmp = -0.5d0*cexp_eef1/r_sqr

dri = r_abs - r_i
drj = r_abs - r_j

if (dri > 0.0d0) then
    xij = (r_abs - r_i)*lambda_i
    f1 = vol(j)*dg_free(i)*tmp*exp(-xij**2)*lambda_i
    if (calc_g) dfdrij1 = 2.0d0*f1*(xij*lambda_i + i_rabs)
else
    l2 = lambda_i*lambda_i
    r2 = r_i*r_i
    c = -0.5d0*cexp_eef1*vol(j)*dg_free(i)*lambda_i
    a = 2.0d0*c/r2 * (l2 + 1.0d0/r_i)
    b = -c*(2.0d0*l2/r_i + 1.0d0/r2)

    f1 = a*r_abs + b
    if (calc_g) dfdrij1 = -a
end if

if (drj > 0.0d0) then
    xji = (r_abs - r_j)*lambda_j
    f2 = vol(i)*dg_free(j)*tmp*exp(-xji**2)*lambda_j
    if (calc_g) dfdrij2 = 2.0d0*f2*(xji*lambda_j + i_rabs)
else
    l2 = lambda_j*lambda_j
    r2 = r_j*r_j
    c = -0.5d0*cexp_eef1*vol(i)*dg_free(j)*lambda_j
    a = 2.0d0*c/r2 * (l2 + 1.0d0/r_j)
    b = -c*(2.0d0*l2/r_j + 1.0d0/r2)

    f2 = a*r_abs + b
    if (calc_g) dfdrij2 = -a
end if

f = f1 + f2
if (calc_g) dfdrij = (dfdrij1 + dfdrij2)*i_rabs

end subroutine calc_eef1_pair
!-------------------------------------------------------------------------------
END MODULE EEF1
