!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/modeling/physics/facts.f90
!
! Description: FACTS (simplified GB/SA) solvation energy
!              It should be called by solvation.f90 to calculate energy
!              This module is for both polar and nonpolar solvation free energy
!              calculation by using FACTS. (other SA terms are calculated in
!              surface_area.f90)
!
! TODO add comment
! Reference
!
!-------------------------------------------------------------------------------
MODULE FACTS
!-------------------------------------------------------------------------------

use globals
use logger, only: log_p, terminate_with_error
use string, only: parse_string
use mathfunctions, only: cubic_spline
use in_out_vars, only: infile_parameter
use energy_vars, only: use_solv, use_SA, use_factsmem, truncate_by_shift, &
    update_reff, die_const_pro, die_const_slv, R, i_R, i_P, &
    Roff_sqr, max_neigh, max_energy
use cutoff_methods, only: nonbonded_shift_function
use molecular_mechanics, only: perm_vac_inv, Egb_fac, LJ_type, QQ_type, qiqj

implicit none
save
private

! General parameters
integer :: n_grp
integer :: n_grp_facts
!
integer, allocatable :: atm_type(:) ! Atom type index (differs by solvation type)
real(dp), allocatable :: vol(:) ! Atomic volume for atom type
integer, allocatable :: neighno(:,:), n_neigh(:) ! Neighbor atom no. index (w/o excluded)
real(dp) :: Roff_sqr_pair, i_Roff_sqr_pair ! distance cut for pair
!
integer, parameter :: max_neigh_facts = max_neigh + 100
integer, allocatable :: neighno_self(:,:), n_neigh_self(:) ! Self interacting pairs
real(dp), allocatable :: Gpair(:)

real(dp), allocatable :: fctA(:,:), fctB(:,:), fctC(:,:), fctD(:,:), fctF(:,:) ! FACTS parameters
real(dp), allocatable :: fct_self_2ci(:), fct_self_rcut(:), fct_self_rcut2(:) ! Distance cut for self pairs
real(dp) :: fct_ikappa, fct_tau ! inverse kappa, tau
! Below are for efficient evaluation
logical, allocatable :: shell_atm(:)
logical, allocatable :: shell_atm_prv(:)
real(dp), allocatable :: reff_prv(:), Gel_prv(:)  ! Effective Born-radius & Free energy elec. contribution
real(dp), allocatable :: reff(:)
real(dp), allocatable :: Ai(:), Bi(:), Ci(:), Di(:) ! A, B, C, D values in FACTS
real(dp), allocatable :: Bnum(:), Bdenom(:), Bv(:,:)! numerator/denominator/vector sum for B value

! Parameters for FACTSMEM
integer :: n_bin(2)
integer, parameter :: max_bin = 100
real(dp), allocatable :: eff_die(:,:), eff_np(:,:)
real(dp), allocatable :: zfac(:), dzfac(:)
real(dp) :: eff_prof_list(max_bin,2), c_pp_list(max_bin,2)
real(dp) :: interval(2)
logical, allocatable :: outside(:)

! Parameters for SA
real(dp), allocatable :: area_init(:) ! Initial area, used for fixed framework

public :: reff
public :: initialize_facts
public :: calc_facts
public :: finalize_facts

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_facts(molecule)
!-------------------------------------------------------------------------------
! TODO comment
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: molecule
! Dummy variables TODO tn%residue check
logical :: appl_res(tn%residue), status
real(dp) :: dumm1(tn%atom), dumm2(tn%atom)

allocate(area_init(tn%atom))
allocate(n_neigh(tn%atom), neighno(max_neigh,tn%atom))
allocate(Gpair(tn%atom))
allocate(shell_atm(tn%atom))
allocate(shell_atm_prv(tn%atom))
if (use_factsmem) then
    allocate(eff_die(2,tn%atom), eff_np(2,tn%atom))
    eff_die(:,:) = 0.0d0
    eff_np(:,:) = 0.0d0
end if
allocate(outside(tn%atom))
!
allocate(n_neigh_self(tn%atom), neighno_self(max_neigh,tn%atom))
allocate(Ai(tn%atom), Bi(tn%atom))
allocate(Ci(tn%atom), Di(tn%atom))
allocate(Bnum(tn%atom), Bdenom(tn%atom))
allocate(Bv(3,tn%atom))
allocate(Gel_prv(tn%atom))
allocate(reff_prv(tn%atom))
allocate(reff(tn%atom))

area_init(:) = 0.0d0
n_neigh(:) = 0
neighno(:,:) = 0
fct_ikappa = 1.0d0 / 12.0d0
fct_tau = Egb_fac
appl_res(:) = .true.
shell_atm(:) = .true.
shell_atm_prv(:) = .true.
!
truncate_by_shift = .true.
Gel_prv(:) = 0.0d0

call setup_facts_parm()
if (use_factsmem) then
    allocate(zfac(tn%atom), dzfac(tn%atom))
    call setup_factsmem_parm()
end if

if (use_SA) then
    call calc_facts_SA(area_init(:), dumm2, appl_res, .false., status)
end if

! If framework is fixed:
! Pre-calculate initial values into XXX_init then overwrite only if changed
call set_neighbor_pair(status)
call eff_radius_FACTS(Gel_prv, dumm1, reff_prv, status)
call calc_gb_init(reff_prv(:))
reff(:) = reff_prv(:)

end subroutine initialize_facts
!-------------------------------------------------------------------------------
subroutine finalize_facts()
!-------------------------------------------------------------------------------
! TODO comment
!-------------------------------------------------------------------------------

!TODO added by curie. check
deallocate(Gpair)
deallocate(outside)
deallocate(Gel_prv, reff_prv, reff)
deallocate(atm_type)
deallocate(vol)
deallocate(fctA, fctB, fctC, fctD, fctF)
deallocate(fct_self_2ci)
deallocate(fct_self_rcut)
deallocate(fct_self_rcut2)
!
deallocate(area_init)
deallocate(n_neigh, neighno)
deallocate(shell_atm)
deallocate(shell_atm_prv)
!
deallocate(n_neigh_self, neighno_self)
deallocate(Ai, Bi, Ci, Di)
deallocate(Bnum, Bdenom, Bv)
if (use_factsmem) then
    deallocate(eff_die, eff_np)
    deallocate(zfac, dzfac)
end if

end subroutine finalize_facts
!-------------------------------------------------------------------------------
subroutine setup_facts_parm()
!-------------------------------------------------------------------------------
! Set initial parameters for FACTS, can be used with CHARMM19/22ph/22
!-------------------------------------------------------------------------------
character(len=len_fname) :: facts_parm_fname, line
integer :: f_unit = 25, openstat
integer :: i, i_grp_facts
real(dp) :: prefac

if (top_type == 'polarh') then
    facts_parm_fname = trim(data_dir)//'facts_parm_charmm19.list'
    n_grp = eng_para%n_atom_cls
    n_grp_facts = 8
    Roff_sqr_pair = Roff_sqr ! Pair energy calculation cut !Charmm19 default : 7.5d0**2
else if (top_type == 'allh_ch22') then
    facts_parm_fname = trim(data_dir)//'facts_parm_charmm22.list'
    n_grp = eng_para%n_atom_cls
    n_grp_facts = 19
    Roff_sqr_pair = Roff_sqr
end if
i_Roff_sqr_pair = 1.0d0/Roff_sqr_pair

allocate(atm_type(n_grp))
allocate(vol(n_grp_facts))
allocate(fctA(0:3,0:n_grp_facts))
allocate(fctB(2,  0:n_grp_facts))
allocate(fctC(0:3,0:n_grp_facts))
allocate(fctD(2,  0:n_grp_facts))
allocate(fctF(2,  0:n_grp_facts))
allocate(fct_self_2ci(0:n_grp_facts))
allocate(fct_self_rcut(0:n_grp_facts))
allocate(fct_self_rcut2(0:n_grp_facts))

! Map atomic types to the corresponding atom types in FACTS
! CHARMM22ph use CHARMM19 non-bonded parameters and will bring
! FACT parameters derived for CHARMM19.
! CHARMM19: 7 atom types
! 1. H (1.0) 2. O/N (1.6) 3. S (1.89) 4. C/CR (2.1) 5. CH3(2.165)
! 6. CH2 (2.235) 7. CH1 (2.365) 8. P (1.70)
if (infile_parameter == trim(data_dir)//'parameters_charmm22ph.in') then
    atm_type(1:46) = &
          (/ 1, 1, 4, 4, 7, 6, 5, 4, 4, 4, &
             4, 7, 6, 6, 4, 2, 2, 2, 2, 2, &
             2, 2, 2, 2, 2, 2, 2, 2, 3, 2, &
             7, 2, 2, 2, 8, 1, 1, 0, 0, 0, &
             0, 4, 1, 2, 2, 0/)
    if (n_grp > 46) atm_type(47:n_grp) = 0
else if (infile_parameter == trim(data_dir)//'parameters_charmm19.in') then
    atm_type(1:35) = &
          (/ 1, 1, 1, 1, 4, 7, 6, 5, 4, 7, &
             7, 4, 2, 2, 2, 2, 2, 2, 2, 2, &
             2, 2, 2, 2, 2, 2, 0, 0, 3, 3, &
             4, 0, 0, 0, 0 /)
    if (n_grp > 36) atm_type(36:n_grp) = 0
! CHARMM22: 19 atom types
! 1. H (0.2245) 2. HS (0.45)  3. HR2 (0.7)  4. HR1 (0.9)  5. HA(1.32)
! 6. HP (1.3582)  7. HR3 (1.4680)  8. O (1.7)  9. OH1 (1.77) 10. CPA (1.8)
! 11. N (1.85)  12. SM (1.975) 13. CA (1.9924)  14. C (2.0) 15. CT3 (2.06)
! 16. CP2 (2.175)  17. (2.235)  18. CP1 (2.275)  19. (2.365)
! only protein atoms are grouped. currently, HA1, HA2, HA3 mapped as 5
else if (infile_parameter == trim(data_dir)//'parameters_charmm22.in') then
    atm_type(1:70) = &
          (/ 1, 1, 5, 1, 6, 5, 4, 3, 7, 2, &
             0, 0, 5, 5, 5, 5, 0,14,13,18, &
            16,15,10,10,10,13,18,16,16,14, &
            14,10,10,10, 0, 0, 0, 0, 0,18, &
             0, 0, 0, 0, 0,14, 0,11,11,11, &
            11,11,11,11,11,11,11,11,11, 8, &
             8, 8, 9, 9, 9, 8, 0, 8,14,12 /)
    if (n_grp > 70) atm_type(71:n_grp) = 0
end if
     
open(f_unit, file = facts_parm_fname, iostat = openstat)

if (openstat > 0) then
    call log_p('Terminate with error: FACTS parameter file not found.')
    call terminate_with_error('Please check {data_directory}.')
end if

200 format(a100)

!read self energy calculation cut
do
    read(f_unit, 200) line
    if (line == 'SELF RCUT') exit
end do
i_grp_facts = 0
do
    read(f_unit, 200) line
    if (line == 'END') then
        exit
    else
        i_grp_facts = i_grp_facts + 1
        read(line, *) fct_self_rcut(i_grp_facts)
    end if
end do
do i = 1, n_grp_facts
    fct_self_rcut2(i) = fct_self_rcut(i)**2
    fct_self_2ci(i) = (1.0d0/fct_self_rcut(i))**2
end do

! read atomic volume
do
    read(f_unit, 200) line
    if (line == 'VDW RADIUS') exit
end do
i_grp_facts = 0
do
    read(f_unit, 200) line
    if (line == 'END') then
        exit
    else
        i_grp_facts = i_grp_facts + 1
        read(line, *) vol(i_grp_facts)
    end if
end do
prefac = (4.0d0/3.0d0)*pi
do i = 1, n_grp_facts
    vol(i) = (vol(i)**3)*prefac
end do

! read A, B parameters
if (use_solv) then
    do
        read(f_unit, 200) line
        if (line == 'FACTS_A') exit
    end do
    i_grp_facts = 0
    do
        read(f_unit, 200) line
        if (line == 'END') then
            exit
        else
            i_grp_facts = i_grp_facts + 1
            read(line, *) fctA(0:3,i_grp_facts)
        end if
    end do
    do
        read(f_unit, 200) line
        if (line == 'FACTS_B') exit
    end do
    i_grp_facts = 0
    do
        read(f_unit, 200) line
        if (line == 'END') then
            exit
        else
            i_grp_facts = i_grp_facts + 1
            read(line, *) fctB(:,i_grp_facts)
        end if
    end do
end if
     
! read C, D, F parameters (F is not used)
if (use_SA) then
    do
        read(f_unit, 200) line
        if (line == 'FACTS_C') exit
    end do
    i_grp_facts = 0
    do
        read(f_unit, 200) line
        if (line == 'END') then
            exit
        else
            i_grp_facts = i_grp_facts + 1
            read(line, *) fctC(0:3,i_grp_facts)
        end if
    end do
    do
        read(f_unit, 200) line
        if (line == 'FACTS_D') exit
    end do
    i_grp_facts = 0
    do
        read(f_unit, 200) line
        if (line == 'END') then
            exit
        else
            i_grp_facts = i_grp_facts + 1
            read(line, *) fctD(:,i_grp_facts)
        end if
    end do
    do
        read(f_unit, 200) line
        if (line == 'FACTS_F') exit
    end do
    i_grp_facts = 0
    do
        read(f_unit, 200) line
        if (line == 'END') then
            exit
        else
            i_grp_facts = i_grp_facts + 1
            read(line, *) fctF(:,i_grp_facts)
        end if
    end do
end if

close(f_unit)
  
end subroutine setup_facts_parm
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
    i_start = 1

    do i_j = i_start, i_P(i_atm)%n_pair
        j_atm = i_P(i_atm)%i_pair(i_j)

        if (atm_type(LJ_type(j_atm)) == 0) cycle

        dij = i_P(i_atm)%d(i_j)
        n_neigh(i_atm) = n_neigh(i_atm) + 1
        neighno(n_neigh(i_atm),i_atm) = j_atm
    end do
end do

end subroutine set_neighbor_pair
!-------------------------------------------------------------------------------
subroutine set_shell(appl_res)
!-------------------------------------------------------------------------------
! Search pairlist on confined region (given by appl_res)
! Appl_res contains residues of interest (e.g. loop to be modeled)
! This subroutine returns shell_atm, the atomlist influenced by appl_res,
! containing appl_res atoms as well as neighboring atoms at framework
!-------------------------------------------------------------------------------
logical, intent(in) :: appl_res(tn%residue)
integer :: i_atm, j_atm, i_j

! Make shell_atm list used for FACTS
shell_atm(:) = .false.

! Make list of shell_atm
! First append appl_res atoms
do i_atm = 1, tn%atom
    if (appl_res(i_R(1,i_atm))) shell_atm(i_atm) = .true.
end do

do i_atm = 1, tn%atom
    if (appl_res(i_R(1,i_atm))) then
        do i_j = 1, i_P(i_atm)%n_pair
            j_atm = i_P(i_atm)%i_pair(i_j)
            shell_atm(j_atm) = .true.
        end do
    else
        do i_j = 1, i_P(i_atm)%n_pair
            j_atm = i_P(i_atm)%i_pair(i_j)
            if (appl_res(i_R(1,j_atm))) then
                shell_atm(i_atm) = .true.
                exit
            end if
        end do
    end if
end do

end subroutine set_shell
!-------------------------------------------------------------------------------
subroutine calc_facts(f, g, appl_res, calc_g)
!-------------------------------------------------------------------------------
! TODO comment
!-------------------------------------------------------------------------------
real(dp), intent(inout) :: f(2), g(3,tn%atom,2)
logical, intent(in) :: appl_res(tn%residue), calc_g
integer, parameter :: reff_update_iter = 1
logical :: non0qi, non0qj, status
integer :: i_atm, j_atm, k_atm, i_j, i_k
real(dp) :: Ri(3), Rij(3), r_sqr, argv(3), arg
real(dp) :: df_s(3)
real(dp) :: ftr, gtr
real(dp) :: q2, qq2, qi, qi2
real(dp) :: f_tmp, f_sa
real(dp) :: Gel(tn%atom), area(tn%atom)
real(dp) :: darea(3,0:max_neigh)
real(dp) :: dGel(3,0:max_neigh), dGdreff(tn%atom)
real(dp) :: dGdrij, dGdreff_tmp(2)
real(dp) :: z_v(3), fct_tau_z, dft_dzi, dft_dz(2) !for factsmem
real(dp) :: dreff_z(2) !for factsmem
real(dp) :: slf, tmp_z !for factsmem

if (use_SA) then
    area(:) = area_init(:)
    if (use_factsmem) then
        call calc_eff_factsmem_parms(eff_prof_list(:,2), c_pp_list(:,2), &
            n_bin(2), interval(2), eff_np(:,:), z_v, calc_g)
    end if
end if

reff(:) = reff_prv(:)
Gel(:) = Gel_prv(:)
dGdreff(:) = 0.0d0
if (use_factsmem) then
    call calc_eff_factsmem_parms(eff_prof_list(:,1), c_pp_list(:,1), &
        n_bin(1), interval(1), eff_die(:,:), z_v, calc_g)
    zfac(:) = 0.0d0
    dzfac(:) = 0.0d0
    call calc_reff_z(calc_g)
end if

!! First set neighbor pairlist proper for solvation term type
call set_shell(appl_res)

!! Preprocess for Eff. Born radius calculation for GB-style
if (update_reff) then
    call set_neighbor_pair(status)
    if (.not. status) then
        f = max_energy
        return
    end if
    call eff_radius_FACTS(Gel, area, reff, status)
    if (.not. status) then
        f = max_energy
        return
    end if
end if
  
!! Energy evaluation: Polar solvation
if (use_solv) then
    ! Self term
    do i_atm = 1, tn%atom
        non0qi = eng_para%non0q(QQ_type(i_atm))
        qi = eng_para%charge(QQ_type(i_atm))

        if (non0qi) then
            if (atm_type(LJ_type(i_atm)) == 0) cycle
            qi2 = qi*qi
            if (use_factsmem) then
                ! Self term calculation for FACTSMEM
                call calc_one_fct_tau_z(i_atm, fct_tau_z, dft_dzi, .false.)
                f(1) = f(1) -0.5d0*fct_tau_z*qi2/(zfac(i_atm)*reff(i_atm))
            else
                f(1) = f(1) + qi2*Gel(i_atm)
            end if
        end if
    end do
    
    ! Pair-interaction terms
    do i_atm = 1, tn%atom
        if (atm_type(LJ_type(i_atm)) == 0) cycle
        ! Bypass if any of surrrounding Reff is unchanged
        !
        ! TODO: This can be problematic with the current pairlist update method
        ! With the previous pairlist update method,
        !  pair information are sorted in order of their atom indices,
        ! however, with the current pairlist update method,
        !  pair information are not sorted, but random-like ordered,
        !  so, this make it impossible or hard to applying this kind of bypass.
        ! The keypoint is neighbor list for each atom changes randomly,
        !  while all pair is maintained.
        !
        !if (facts_on) then
        !   call check_env_changed(i_atm, reff_changed(:), env_changed)
        !   if (.not. env_changed) then
        !      f(1) = f(1) + Gpair(i_atm)
        !      cycle
        !   end if
        !end if

        ! Recalculate Gpair if environment is changed
        Ri = R(:,i_atm)
        non0qi = eng_para%non0q(QQ_type(i_atm))
        Gpair(i_atm) = 0.0d0
        do i_j = 1, n_neigh(i_atm)
            j_atm = neighno(i_j,i_atm)
            non0qj = eng_para%non0q(QQ_type(j_atm))
            
            Rij(:) = R(:,j_atm) - Ri(:)
            r_sqr = dot_product(Rij,Rij)
            if (r_sqr < small_real) cycle
            
            ! Non-bonded interaction truncation
            if (truncate_by_shift) then  ! FACTS
                call nonbonded_shift_function(r_sqr, Roff_sqr_pair, ftr, gtr, calc_g)
            end if

            if (non0qi .and. non0qj) then
                ! use GB formular
                q2 = qiqj(QQ_type(i_atm), QQ_type(j_atm))
                if (use_factsmem) then
                    call calc_fct_tau_z(i_atm, j_atm, fct_tau_z, dft_dz, calc_g)
                    qq2 = q2*fct_tau_z
                    call calc_gb(f_tmp, dGdrij, dGdreff_tmp, qq2, r_sqr, &
                        zfac(i_atm)*reff(i_atm), zfac(j_atm)*reff(j_atm), calc_g)
                    dGdreff(i_atm) = dGdreff(i_atm) + dGdreff_tmp(1)*ftr
                    dGdreff(j_atm) = dGdreff(j_atm) + dGdreff_tmp(2)*ftr
                    if ((appl_res(i_R(1,i_atm)) .or. appl_res(i_R(1,j_atm))) .and. calc_g) then
                        df_s(:) = dGdrij*ftr*Rij(:) - f_tmp*gtr*Rij(:)
                        g(:,i_atm,1) = g(:,i_atm,1) + df_s(:)
                        g(:,j_atm,1) = g(:,j_atm,1) - df_s(:)
                        !z dependency of tau
                        g(:,i_atm,1) = g(:,i_atm,1) + ftr*dft_dz(1)*f_tmp/fct_tau_z*z_v(:)
                        g(:,j_atm,1) = g(:,j_atm,1) + ftr*dft_dz(2)*f_tmp/fct_tau_z*z_v(:)
                        !z dependency of Reff
                        call calc_gradient_reff_z(i_atm, j_atm, reff(i_atm), reff(j_atm), &
                            r_sqr, f_tmp, dreff_z)
                        g(:,i_atm,1) = g(:,i_atm,1) + dreff_z(1)*ftr*z_v(:)
                        g(:,j_atm,1) = g(:,j_atm,1) + dreff_z(2)*ftr*z_v(:)
                    end if
                else
                    qq2 = q2*fct_tau
                    call calc_gb(f_tmp, dGdrij, dGdreff_tmp, qq2, r_sqr, &
                        reff(i_atm), reff(j_atm), calc_g)
                    dGdreff(i_atm) = dGdreff(i_atm) + dGdreff_tmp(1)*ftr
                    dGdreff(j_atm) = dGdreff(j_atm) + dGdreff_tmp(2)*ftr
                    if ((appl_res(i_R(1,i_atm)) .or. appl_res(i_R(1,j_atm))) .and. calc_g) then
                        df_s(:) = dGdrij*ftr*Rij(:) - f_tmp*gtr*Rij(:)
                        g(:,i_atm,1) = g(:,i_atm,1) + df_s(:)
                        g(:,j_atm,1) = g(:,j_atm,1) - df_s(:)
                    end if
                end if
                Gpair(i_atm) = Gpair(i_atm) + f_tmp*ftr
            end if
        end do ! Iter over j_atm (upper-diagonal pairs)
        f(1) = f(1) + Gpair(i_atm)
    end do ! Iter over i_atm
end if

!! FACTS solvation gradient / SA
do i_atm = 1, tn%atom
    if (atm_type(LJ_type(i_atm)) == 0) cycle
    if (shell_atm(i_atm) .or. shell_atm_prv(i_atm)) then ! Refresh only if included in shell_atm
        if (calc_g) then
            call gradient_for_FACTS(i_atm, n_neigh_self(i_atm), &
            dGel(:,0:n_neigh_self(i_atm)), darea(:,0:n_neigh_self(i_atm)))
            !
            if (use_solv) then
                qi = eng_para%charge(QQ_type(i_atm))
                qi2 = qi*qi
                non0qi = eng_para%non0q(QQ_type(i_atm))
                if (non0qi) then
                    if (use_factsmem) then !z dependency of tau and reff
                        call calc_one_fct_tau_z(i_atm, fct_tau_z, dft_dzi, calc_g)
                        slf = -0.5d0*fct_tau_z*qi2 / (zfac(i_atm)*reff(i_atm))
                        tmp_z = (dft_dzi/zfac(i_atm) - &
                            fct_tau_z*dzfac(i_atm)*((1.0d0/zfac(i_atm))**2.0d0)) *qi2/fct_tau
                        g(:,i_atm,1) = g(:,i_atm,1) + (slf/Gel(i_atm))*dGel(:,0)
                        g(:,i_atm,1) = g(:,i_atm,1) + tmp_z*Gel(i_atm)*z_v(:)
                        do i_j = 1, n_neigh_self(i_atm)
                            j_atm = neighno_self(i_j,i_atm)
                            g(:,j_atm,1) = g(:,j_atm,1) + (slf/Gel(i_atm))*dGel(:,i_j)
                        end do
                    else
                        g(:,i_atm,1) = g(:,i_atm,1) + qi2*dGel(:,0)
                        do i_j = 1, n_neigh_self(i_atm)
                            j_atm = neighno_self(i_j,i_atm)
                            g(:,j_atm,1) = g(:,j_atm,1) + qi2*dGel(:,i_j)
                        end do
                    end if
                end if
                
                if (abs(dGdreff(i_atm)) > small_real) then
                    if (use_factsmem) then
                        arg = dGdreff(i_atm)/Gel(i_atm)*zfac(i_atm)*reff(i_atm)
                    else
                        arg = dGdreff(i_atm)/Gel(i_atm)*reff(i_atm)
                    end if

                    do i_k = 1, n_neigh_self(i_atm)
                        k_atm = neighno_self(i_k,i_atm)
                        argv(:) = arg*dGel(:,i_k)
                        g(:,k_atm,1) = g(:,k_atm,1) + argv(:)
                        g(:,i_atm,1) = g(:,i_atm,1) - argv(:)
                    end do
                end if
            end if
            !
            if (use_SA) then
                f_sa = area(i_atm)
                if (use_factsmem) then
                    f(2) = f(2) + eff_np(1,i_atm)*f_sa
                else
                    f(2) = f(2) + f_sa
                end if
                
                if (use_factsmem) then
                    g(:,i_atm,2) = g(:,i_atm,2) + eff_np(1,i_atm)*darea(:,0)
                    g(:,i_atm,2) = g(:,i_atm,2) + eff_np(2,i_atm)*z_v(:)*f_sa
                    do i_j = 1, n_neigh_self(i_atm)
                        j_atm = neighno_self(i_j,i_atm)
                        g(:,j_atm,2) = g(:,j_atm,2) + eff_np(1,i_atm)*darea(:,i_j)
                    end do
                else
                    g(:,i_atm,2) = g(:,i_atm,2) + darea(:,0)
                    do i_j = 1, n_neigh_self(i_atm)
                        j_atm = neighno_self(i_j,i_atm)
                        g(:,j_atm,2) = g(:,j_atm,2) + darea(:,i_j)
                    end do
                end if
            end if
        else if (use_SA) then
            f_sa = area(i_atm)
            if (use_factsmem) then
                f(2) = f(2) + eff_np(1,i_atm)*f_sa
            else
                f(2) = f(2) + f_sa
            end if
        end if
        
        ! Just copy initial value if out of shell_atm
    else
        if (use_SA) then
            f_sa = area_init(i_atm)
            if (use_factsmem) then
                f(2) = f(2) + eff_np(1,i_atm)*f_sa
            else
                f(2) = f(2) + f_sa
            end if
        end if
    end if
end do

! Save current Born radii
Gel_prv(:) = Gel(:)
reff_prv(:) = reff(:)
shell_atm_prv(:) = shell_atm(:)

end subroutine calc_facts
!-------------------------------------------------------------------------------
subroutine calc_gb_init(reff_in)
!-------------------------------------------------------------------------------
! Pre-calculate initial gb pair energy for future use :)
! This for FACTS so far...
!-----------------------------------------------------------------------------
real(dp), intent(in) :: reff_in(tn%atom)
logical :: non0qi, non0qj
integer :: i_type, j_type, i_j, i_atm, j_atm
real(dp) :: q2, qq2, f_tmp, ftr, gtr
real(dp) :: r_sqr, Rij(3), Ri(3)
real(dp) :: dumm1, dumm2(2)
real(dp) :: z_v(3), fct_tau_z, dft_dz(2) !for factsmem

if (use_factsmem) then
    call calc_eff_factsmem_parms(eff_prof_list(:,1), c_pp_list(:,1), &
        n_bin(1), interval(1), eff_die, z_v, .false.)
    call calc_reff_z(.false.)
end if

Gpair(:) = 0.0d0
do i_atm = 1, tn%atom
    non0qi = eng_para%non0q(QQ_type(i_atm))
    i_type = atm_type(LJ_type(i_atm))
    if (i_type == 0) cycle

    do i_j = 1, n_neigh(i_atm)
        j_atm = neighno(i_j,i_atm)
        j_type = atm_type(LJ_type(j_atm))
        non0qj = eng_para%non0q(QQ_type(j_atm))
     
        if (j_type == 0) cycle
        
        Rij(:) = R(:,j_atm) - Ri(:)
        r_sqr = dot_product(Rij, Rij)

        q2 = qiqj(QQ_type(i_atm), QQ_type(j_atm))
        if (use_factsmem) then
            call calc_fct_tau_z(i_atm, j_atm, fct_tau_z, dft_dz, .false.)
            qq2 = q2*fct_tau_z
            call calc_gb(f_tmp, dumm1, dumm2, qq2, r_sqr, &
                zfac(i_atm)*reff_in(i_atm), zfac(j_atm)*reff_in(j_atm), .false.)
        else
            qq2 = q2*fct_tau
            call calc_gb(f_tmp, dumm1, dumm2, qq2, r_sqr, &
                reff_in(i_atm), reff_in(j_atm), .false.)
        end if

        call nonbonded_shift_function(r_sqr, Roff_sqr_pair, ftr, gtr, .false.)
        Gpair(i_atm) = Gpair(i_atm) + f_tmp*ftr
    end do
end do

end subroutine calc_gb_init
!-------------------------------------------------------------------------------
subroutine calc_facts_SA(area, g, appl_res, calc_g, status)
!-------------------------------------------------------------------------------
! Simple evaluation of SA for current conformation R
!-------------------------------------------------------------------------------
logical, intent(in) :: appl_res(tn%residue), calc_g
real(dp), intent(out) :: area(tn%atom), g(3,tn%atom)
logical, intent(out) :: status
real(dp) :: dumm1(tn%atom), dumm2(tn%atom)

call set_shell(appl_res)
call eff_radius_FACTS(dumm1, area, dumm2, status)

end subroutine calc_facts_SA
!-------------------------------------------------------------------------------
subroutine calc_gb(f, dGdrij, dGdreff, qq2, r_sqr, reffi, reffj, calc_g)
!-------------------------------------------------------------------------------
! Read parameters to calculate GB pair equation (self term is not here)
! reffi, reffj: effective Born radii
! r_sqr, r_abs, Rij: distance & vectors between ij
! qq2: charges multiplicated
! f: interaction energy
! dGdrij: derivative of f about distance
! dGdreff(2): derivative of f about effective born radius (each for i,j)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: f, dGdrij, dGdreff(2)
real(dp), intent(in) :: qq2, r_sqr
real(dp), intent(in) :: reffi, reffj
logical, intent(in) :: calc_g
real(dp) :: g1, g2
real(dp) :: tmp1, tmp2, tmp3, reffm
  
f = 0.0d0
dGdreff(:) = 0.0d0
reffm = reffi*reffj

! If one of effective Born radius is too small (fully buried)
if (reffm < small_real) then
    f = -qq2/sqrt(r_sqr)
     
    if (calc_g) then
        dGdrij = f/r_sqr
        dGdreff(:) = 0.0d0
    end if

else
    tmp1 = r_sqr*fct_ikappa
    tmp2 = exp(-tmp1/reffm)
    tmp3 = r_sqr + reffm*tmp2 ! GB eq. denominator
    f = -qq2/sqrt(tmp3)

    if (calc_g) then
        g1 = 0.5d0*f/tmp3
        g2 = 2.0d0 - 2.0d0*fct_ikappa*tmp2

        dGdrij = g1*g2 ! dGij/drij
        dGdreff(1) = g1*tmp2*(reffj + tmp1/reffi) !dGij/dReffi
        dGdreff(2) = g1*tmp2*(reffi + tmp1/reffj) !dGij/dReffj
    end if
end if
  
end subroutine calc_gb
!-------------------------------------------------------------------------------
subroutine eff_radius_FACTS(Gel, area, reff, status)
!-------------------------------------------------------------------------------
! Calculate effective Born radius, using heuristic method by FACTS
!-------------------------------------------------------------------------------
real(dp), intent(inout) :: Gel(tn%atom), reff(tn%atom)
real(dp), intent(inout) :: area(tn%atom)
logical, intent(out) :: status
real(dp) :: arg, argv(2)
real(dp) :: r_sqr, Rij(3), i_rsqr, i_rabs
real(dp) :: rcuti2, rcutj2, fct_tau_2, Ri(3)
real(dp) :: d2_pair(max_neigh_facts), d_pair(max_neigh_facts)
integer :: i, j, i_atm, j_atm, i_j
integer :: pairs(max_neigh_facts), n_pair

fct_tau_2 = 0.5d0*fct_tau

Bv(:,:) = 0.0d0
Bnum(:) = 0.0d0
Bdenom(:) = 1.0d0 ! 1 + Sum*(Vij/rij*theta)
Ai(:) = 0.0d0
Bi(:) = 0.0d0
Ci(:) = 0.0d0
Di(:) = 0.0d0
n_neigh_self(:) = 0
neighno_self(:,:) = 0
status = .true.

! 1. Calculate A, B and their derivatives
do i_atm = 1, tn%atom
    i = atm_type(LJ_type(i_atm))
    if (i == 0) cycle
    rcuti2 = fct_self_rcut2(i)
    Ri(:) = R(:,i_atm)

    call get_self_pair(i_atm, pairs, n_pair, d2_pair, d_pair, status)
    if (.not. status) return

    ! Eval
    do i_j = 1, n_pair
        j_atm = pairs(i_j)
        j = atm_type(LJ_type(j_atm))
        if (j == 0) cycle

        rcutj2 = fct_self_rcut2(j)
        Rij(:) = R(:,j_atm) - Ri(:)
        r_sqr = d2_pair(i_j)

        if (r_sqr < small_real .or. (r_sqr > rcuti2 .and. r_sqr > rcutj2)) cycle

        ! Inverse distance for efficiency
        i_rsqr = 1.0d0/r_sqr
        i_rabs = 1.0d0/d_pair(i_j)

        if (r_sqr <= rcuti2) then
            n_neigh_self(i_atm) = n_neigh_self(i_atm) + 1
            neighno_self(n_neigh_self(i_atm),i_atm) = j_atm

            arg = vol(j)*theta(r_sqr,fct_self_2ci(i)) ! 1: theta, 2: C for dtheta/dr(:) = CR(:)

            Ai(i_atm) = Ai(i_atm) + arg
            Bv(:,i_atm) = Bv(:,i_atm) + arg*i_rsqr*Rij(:)
            Bdenom(i_atm) = Bdenom(i_atm) + arg*i_rabs
        end if

        if (r_sqr <= rcutj2) then
            n_neigh_self(j_atm) = n_neigh_self(j_atm) + 1
            neighno_self(n_neigh_self(j_atm),j_atm) = i_atm

            arg = vol(i)*theta(r_sqr,fct_self_2ci(j)) ! 1: theta, 2: C for dtheta/dr(:) = CR(:)
           
            Ai(j_atm) = Ai(j_atm) + arg
            Bv(:,j_atm) = Bv(:,j_atm) - arg*i_rsqr*Rij(:)
            Bdenom(j_atm) = Bdenom(j_atm) + arg*i_rabs
        end if
    end do
end do

! 2. Process A and B to get effective Born radius, free energy, and surface area
do i_atm = 1, tn%atom
    if ((.not. shell_atm(i_atm)) .and. (.not. shell_atm_prv(i_atm))) cycle

    i = atm_type(LJ_type(i_atm))
    Bnum(i_atm) = sqrt(dot_product(Bv(:,i_atm),Bv(:,i_atm))) ! |Sum(Vj*Theta*xij/rij)|
    Bi(i_atm) = Bnum(i_atm)/Bdenom(i_atm)

    ! Calculate polar solvation using arrays above
    if (use_solv) then
        Ci(i_atm) = Ai(i_atm) + fctB(1,i)*Bi(i_atm) &
             + fctB(2,i)*Ai(i_atm)*Bi(i_atm)
        argv(:) = sigmoidal_exp(Ci(i_atm), fctA(0:3,i), .false.)
        Gel(i_atm) = argv(1)
        reff(i_atm) = -fct_tau_2/Gel(i_atm)
        !reff_changed(i_atm) = .true.
    end if

    ! Calculate SA using arrays above
    if (use_SA) then
        Di(i_atm) = Ai(i_atm) + fctD(1,i)*Bi(i_atm) &
             + fctD(2,i)*Ai(i_atm)*Bi(i_atm)
        argv(:) = sigmoidal_exp(Di(i_atm), fctC(0:3,i), .false.)
        area(i_atm) = argv(1)
    end if
end do

end subroutine eff_radius_FACTS
!-------------------------------------------------------------------------------
subroutine get_self_pair(i_atm, pairs, n_pair, d2_pair, d_pair, status)
!-------------------------------------------------------------------------------
integer, intent(in) :: i_atm
integer, intent(out) :: pairs(max_neigh_facts), n_pair
real(dp), intent(out) :: d2_pair(max_neigh_facts), d_pair(max_neigh_facts)
logical, intent(out) :: status
integer :: i_j, j_atm
real(dp) :: dij

! Get pair
pairs(:) = 0
n_pair = 0
status = .true.

do i_j = 1, i_P(i_atm)%n_pair
    j_atm = i_P(i_atm)%i_pair(i_j)
    if (j_atm > tn%atom) cycle
    if (.not. ((shell_atm(i_atm) .or. shell_atm_prv(i_atm)) .or.&
                (shell_atm(j_atm) .or. shell_atm_prv(j_atm)))) cycle

    n_pair = n_pair + 1
    if (n_pair > max_neigh_facts) then
        status = .false.
        return
    end if
    pairs(n_pair) = j_atm
    dij = i_P(i_atm)%d(i_j)
    d_pair(n_pair) = dij
    d2_pair(n_pair) =  dij*dij
end do

end subroutine get_self_pair
!-------------------------------------------------------------------------------
subroutine gradient_for_FACTS(i_atm, n_n, dGel, darea)
!-------------------------------------------------------------------------------
! Return derivatives about Cartesian coordinate for a single atom
! dGel: derivative of free energy 
! darea: derivative of surface area 
!-------------------------------------------------------------------------------
integer, intent(in) :: i_atm, n_n
real(dp), intent(out) :: darea(3,0:n_n)
real(dp), intent(out) :: dGel(3,0:n_n)
real(dp) :: dBcd, dBcn, Sarg, Garg
real(dp) :: dAidr(3), dBidr(3), dCidr(3), dDidr(3)
real(dp) :: tmp1, tmp2, arg(2), arg12(3)
real(dp) :: argv1(3), argv2(3,3), argv3(3)
real(dp) :: r_sqr, Rij(3), Ri(3), i_rsqr, i_rabs
integer :: i, j, i_j, j_atm

i = atm_type(LJ_type(i_atm))

darea(:,:) = 0.0d0
dGel(:,:) = 0.0d0

dBcd = -Bi(i_atm)/Bdenom(i_atm)
dBcn = -1.0d0/(Bnum(i_atm)*Bdenom(i_atm))

! dG/dCi
arg(:) = sigmoidal_exp(Ci(i_atm), fctA(0:3,i), .true.)
Garg = arg(2)

! dSA/dDi
arg(:) = sigmoidal_exp(Di(i_atm), fctC(0:3,i), .true.)
Sarg = arg(2)
     
Ri(:) = R(:,i_atm)

do i_j = 1, n_neigh_self(i_atm)
    j_atm = neighno_self(i_j,i_atm)
    j = atm_type(LJ_type(j_atm))
     
    ! First, Calculate derivatives of Bi, Ci, Di about Cartesian coordinate
    Rij(:) = R(:,j_atm) - Ri(:)
    r_sqr = dot_product(Rij,Rij)
    i_rsqr = 1.0d0/r_sqr
    i_rabs = sqrt(i_rsqr)
     
    arg(:) = vol(j)*gtheta(r_sqr,fct_self_2ci(i)) ! 1: theta, 2: C for dtheta/dr(:) = CR(:)
    tmp1 = arg(1)*i_rsqr ! theta*Vj*xij/r
    tmp2 = (tmp1 + arg(2))*i_rabs

    dAidr(:) = arg(2)*Rij(:)

    arg12 = (-2.0d0*arg(1)*i_rsqr - arg(2))
    argv1(:) = arg12*Rij(:)
    argv2(:,1) = argv1(:)*Rij(1)
    argv2(:,2) = argv1(:)*Rij(2)
    argv2(:,3) = argv1(:)*Rij(3)
    argv2(1,1) = argv2(1,1) + arg(1)
    argv2(2,2) = argv2(2,2) + arg(1)
    argv2(3,3) = argv2(3,3) + arg(1)
    argv2(:,:) = argv2(:,:)*i_rsqr

    argv3(1) = dot_product(Bv(:,i_atm),argv2(:,1))
    argv3(2) = dot_product(Bv(:,i_atm),argv2(:,2))
    argv3(3) = dot_product(Bv(:,i_atm),argv2(:,3))
    dBidr(:) = dBcn*argv3(:) + dBcd*tmp2*Rij(:)

    ! Then combine into dGel and darea
    if (use_solv) then
        dCidr(:) = dAidr(:) + fctB(1,i)*dBidr(:) + &
             fctB(2,i)*(Ai(i_atm)*dBidr(:) + Bi(i_atm)*dAidr(:))
     
        dGel(:,0) = dGel(:,0) + Garg*dCidr(:)
        dGel(:,i_j) = -Garg*dCidr(:)
    end if

    ! Calculate SA using arrays above
    if (use_SA) then
        dDidr(:) = dAidr(:) + fctD(1,i)*dBidr(:) + &
             fctD(2,i)*(Ai(i_atm)*dBidr(:) + Bi(i_atm)*dAidr(:))
        
        darea(:,0) = darea(:,0) + Sarg*dDidr(:)
        darea(:,i_j) = -Sarg*dDidr(:)
    end if
end do

end subroutine gradient_for_FACTS
!-------------------------------------------------------------------------------
! Below subroutines are used for factsmem calculation
!-------------------------------------------------------------------------------
subroutine read_factsmem_parm(eff_die_prof, eff_np_prof)
!-------------------------------------------------------------------------------
! Read parameters needed for FACTSMEM calculation (effective dielectric constant
! profiles and prefactors for nonpoloar solvation energy calculation)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: eff_die_prof(max_bin)   ! 1:z_bin, 2: eff_die
real(dp), intent(out) :: eff_np_prof(max_bin)    ! 1:z_bin, 2: prefactor for np
character(len=len_fname) :: fname(2), line, word(2)
integer :: f_unit = 27, openstat, ioerror, num_word
integer :: i_file, i_bin
real(dp) :: eff_prof(max_bin,2)      ! 1: eff_die, 2: prefactor for np
real(dp) :: prev_z, curr_z

fname(1) = trim(data_dir)//'factsmem_die.dat'
fname(2) = trim(data_dir)//'factsmem_np.dat'
  
n_bin(:) = (/ 0, 0 /)

eff_prof(:,:) = 0.0d0
eff_die_prof(:) = 0.0d0
eff_np_prof(:) = 0.0d0

do i_file = 1, 2
    open(f_unit, file = trim(fname(i_file)), iostat=openstat)
    if (openstat > 0) then
        call terminate_with_error('Terminate with error: No factsmem data files.')
    end if
     
    i_bin = 0
    curr_z = 0.0d0
    prev_z = 0.0d0
    do
        read(f_unit,"(A120)", iostat = ioerror) line
        if (ioerror < 0) exit
        call parse_string(line, num_word, word)
        if (num_word /= 2) cycle
        i_bin = i_bin+1
        read(word(2),"(F10.4)") eff_prof(i_bin, i_file)
        prev_z = curr_z
        read(word(1),"(F10.4)") curr_z
    end do

    n_bin(i_file) = i_bin
    interval(i_file) = curr_z - prev_z
    close(f_unit)
end do

eff_die_prof(:) = eff_prof(:,1)
eff_np_prof(:) = eff_prof(:,2)

end subroutine read_factsmem_parm
!-------------------------------------------------------------------------------
subroutine setup_factsmem_parm()
!-------------------------------------------------------------------------------
! Based on the atom's z-coordinate, calculates effective dielectric constant and
! prefactor for nonpolar solvation energy calculation using cubic interpolation.
!-------------------------------------------------------------------------------
real(dp) :: eff_die_prof(max_bin), eff_np_prof(max_bin)
integer :: i
real(dp), allocatable :: c_pp(:)

call read_factsmem_parm(eff_die_prof, eff_np_prof)

eff_prof_list(:,1) = eff_die_prof(:)
eff_prof_list(:,2) = eff_np_prof(:)

! setup for cubic spline interpolation
c_pp_list(:,:) = 0.0d0
do i = 1, 2
    allocate(c_pp(n_bin(i)))
    call build_profile_spline(n_bin(i), interval(i), eff_prof_list(:,i), c_pp)
    c_pp_list(1:n_bin(i),i) = c_pp(1:n_bin(i))
    deallocate(c_pp)
end do

end subroutine setup_factsmem_parm
!-------------------------------------------------------------------------------
subroutine build_profile_spline(num_bin, interval, c, c_pp)
!-------------------------------------------------------------------------------
! setup second derivatives for spline interpolation (FACTSMEM profiles)
!-------------------------------------------------------------------------------
integer, intent(in) :: num_bin
real(dp), intent(in) :: interval
real(dp), intent(in) :: c(max_bin)
real(dp), intent(out) :: c_pp(num_bin)

integer :: i, j
real(dp), allocatable :: x(:), u(:)
real(dp) ::sig, p

allocate(x(num_bin))
allocate(u(num_bin))

c_pp(:) = 0.0d0

! x values
do i = 1, num_bin
    x(i) = 0.0d0 + (i-1)*interval
end do

do i = 2, num_bin-1
    sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
    p = sig*c_pp(i-1) + 2.0d0
    c_pp(i) = (sig - 1.0d0) / p
    u(i) = (6.0d0*((c(i+1) - c(i))/(x(i+1)-x(i)) - (c(i)-c(i-1))/(x(i)-x(i-1))) / &
          (x(i+1)-x(i-1)) - sig*u(i-1)) / p
end do
  
c_pp(num_bin) = 0.0d0
do j = num_bin-1, 1, -1
    c_pp(j) = c_pp(j)*c_pp(j+1) + u(j)
end do

deallocate(x,u)

end subroutine build_profile_spline
!-------------------------------------------------------------------------------
subroutine calc_eff_factsmem_parms(eff_prof, c_pp, num_bin, interval, eff_parm, z_v, calc_g)
!-------------------------------------------------------------------------------
! calculate effective factsmem parameters for all atoms by cubic spline
!-------------------------------------------------------------------------------
real(dp), intent(in) :: eff_prof(max_bin), c_pp(max_bin), interval
integer, intent(in) :: num_bin
logical, intent(in) :: calc_g
real(dp), intent(out) :: eff_parm(2,tn%atom) ! 1: function, 2: gradient
real(dp), intent(out) :: z_v(3)

real(dp) :: arg
real(dp) :: mid_R(3), z_R
real(dp) :: c_j(2), c_jpp(2), a   !for cubic spline interpolation
real(dp) :: spline(2)  ! 1: function 2: gradient
integer :: i_atm, i_bin, int_bin

outside(:) = .true.
  
!input coordinates should be z-axis aligned (as located in the membrane)
z_v(:) = (/0.0d0, 0.0d0, 1.0d0/)
mid_R(:) = (/0.0d0, 0.0d0, 0.0d0/)

eff_parm(:,:) = 0.0d0
do i_atm = 1, tn%atom
    arg = dot_product((R(:,i_atm)-mid_R(:)),z_v(:))
    if (arg < 0.0d0) then
        outside(i_atm) = .false.
    end if
    z_R = abs(arg) !without sign. out-in symmetric
    int_bin = int(z_R/interval)
    i_bin = int_bin + 1
    ! outside boundary (don't do cubic spline)
    if (i_bin >= num_bin) then
        i_bin = num_bin-1
        a = 0.0d0
    else
        a = 1.0-(z_R-int_bin*interval)/interval
    end if
    !DEBUG(sometimes it goes off boundary. need to check)
    c_j(1) = eff_prof(i_bin)
    c_j(2) = eff_prof(i_bin+1)
    c_jpp(1) = c_pp(i_bin)
    c_jpp(2) = c_pp(i_bin+1)
    spline(:) = cubic_spline(interval, c_j(:), c_jpp(:), a, calc_g)
    eff_parm(1, i_atm) = spline(1)
    if (calc_g) then
        if (outside(i_atm)) then
            eff_parm(2,i_atm) = spline(2)
        else
            eff_parm(2,i_atm) = -spline(2)
        end if
    end if
end do

end subroutine calc_eff_factsmem_parms
!-------------------------------------------------------------------------------
subroutine calc_one_fct_tau_z(i_atm, fct_tau_z, dft_dz, calc_g)
!-------------------------------------------------------------------------------
! calculate tau_z  tau = 1/die_p - 1/die_s(zi)
!-------------------------------------------------------------------------------
integer, intent(in) :: i_atm
logical, intent(in) :: calc_g
real(dp), intent(out) :: fct_tau_z, dft_dz

fct_tau_z = perm_vac_inv*(1.0d0/die_const_pro - 1.0d0/eff_die(1,i_atm))

if (calc_g) then
    dft_dz = perm_vac_inv*eff_die(2,i_atm) / (eff_die(1,i_atm)**2.0d0)
end if

end subroutine calc_one_fct_tau_z
!-------------------------------------------------------------------------------
subroutine calc_fct_tau_z(i_atm, j_atm, fct_tau_z, dft_dz, calc_g)
!-------------------------------------------------------------------------------
! calculate tau_z  tau = 1/die_p - 1/((die(zi)+die(zj))/2)
!-------------------------------------------------------------------------------
integer, intent(in) :: i_atm, j_atm
logical, intent(in) :: calc_g
real(dp), intent(out) :: fct_tau_z, dft_dz(2)

real(dp) :: eff_die_s(2), eff_die_solv
real(dp) :: dft_dz_s(2)
integer :: atm_s(2), i

atm_s(1) = i_atm
atm_s(2) = j_atm

eff_die_s(:) = 0.0d0
do i = 1, 2
    eff_die_s(i) = eff_die(1,atm_s(i))
end do

dft_dz_s(:) = 0.0d0
if (calc_g) then
    do i = 1, 2
        dft_dz_s(i) = eff_die(2,atm_s(i))
    end do
end if
  
eff_die_solv = (eff_die_s(1)+eff_die_s(2))/2.0d0
fct_tau_z = perm_vac_inv*(1.0d0/die_const_pro - 1.0d0/eff_die_solv)

dft_dz(:) = 0.0d0
if (calc_g) then
    dft_dz(1) = 2.0d0*perm_vac_inv*dft_dz_s(1)/((eff_die_s(1)+eff_die_s(2))**2)
    dft_dz(2) = 2.0d0*perm_vac_inv*dft_dz_s(2)/((eff_die_s(1)+eff_die_s(2))**2)
end if

end subroutine calc_fct_tau_z
!-------------------------------------------------------------------------------
subroutine calc_reff_z(calc_g)
!-------------------------------------------------------------------------------
! calculate effectvie Born radius based on approximated kirkwood's equation.
!-------------------------------------------------------------------------------
logical, intent(in) :: calc_g

real(dp) :: beta_w, beta_z
real(dp) :: num, denom, tmp1, tmp2
integer :: i_atm

do i_atm = 1, tn%atom
    beta_w = die_const_pro/die_const_slv
    beta_z = die_const_pro/eff_die(1,i_atm)
    num = 1.0d0 + 0.5d0*beta_z
    denom = 1.0d0 + 0.5d0*beta_w
    zfac(i_atm) = num/denom

    if (calc_g) then
        tmp1 = -0.5d0*die_const_pro*eff_die(2,i_atm)
        tmp2 = denom*(eff_die(1,i_atm)**2.0d0)
        dzfac(i_atm) = tmp1/tmp2
    end if
end do

end subroutine calc_reff_z
!-------------------------------------------------------------------------------
subroutine calc_gradient_reff_z(i_atm, j_atm, reff_wi, reff_wj, r_sqr, f_tmp, dreff_z)
!-------------------------------------------------------------------------------
! Calculate gradient term for z dependent factor multipled to Reff in FACTSMEM
!-------------------------------------------------------------------------------
integer, intent(in) :: i_atm, j_atm
real(dp), intent(in) :: reff_wi, reff_wj  !reff in water
real(dp), intent(in) :: r_sqr, f_tmp
real(dp), intent(out) :: dreff_z(2)
  
real(dp) :: reffm, tmp1, tmp2, tmp3

dreff_z(:) = 0.0d0

reffm = zfac(i_atm)*reff_wi*zfac(j_atm)*reff_wj
tmp1 = exp(-r_sqr*fct_ikappa/reffm)
tmp2 = r_sqr + reffm*tmp1
tmp3 = -0.5d0*f_tmp/tmp2

dreff_z(1) = (reff_wi*zfac(j_atm)*reff_wj + r_sqr*fct_ikappa/zfac(i_atm)) * &
       tmp1*tmp3*dzfac(i_atm)
dreff_z(2) = (reff_wj*zfac(i_atm)*reff_wi + r_sqr*fct_ikappa/zfac(j_atm)) * &
       tmp1*tmp3*dzfac(j_atm)


end subroutine calc_gradient_reff_z
!-------------------------------------------------------------------------------
! Below are functions called by FACTS
!-------------------------------------------------------------------------------
function sigmoidal_exp(x, c, calc_g)
!-------------------------------------------------------------------------------
! Sigmoidal function used for FACTS
!-------------------------------------------------------------------------------
real(dp) :: sigmoidal_exp(2)
logical, intent(in) :: calc_g
real(dp), intent(in) :: x, c(4)
real(dp) :: arg, arg2

arg = exp(-c(3)*(x-c(4)))
arg2 = c(2)/(1.0d0 + arg)

sigmoidal_exp(1) = c(1) + arg2
if (calc_g) sigmoidal_exp(2) = c(3)*arg*arg2/(1.0d0 + arg)

end function sigmoidal_exp
!-------------------------------------------------------------------------------
function theta(r_sqr, c)
!-------------------------------------------------------------------------------
! Sigmoidal function used for FACTS
! 1: theta function
! 2~4: dtheta/dr(:)
!-------------------------------------------------------------------------------
real(dp) :: theta
real(dp), intent(in) :: r_sqr, c
real(dp) :: x

x = 1.0d0 - r_sqr*c

if (x >= 0.0d0) then
    theta = x*x
else
    theta = 0.0d0
end if

end function theta
!-------------------------------------------------------------------------------
function gtheta(r_sqr, c)
!-------------------------------------------------------------------------------
! Sigmoidal function used for FACTS
! 1: theta function
! 2~4: dtheta/dr(:)
!-------------------------------------------------------------------------------
real(dp) :: gtheta(2)
real(dp), intent(in) :: r_sqr, c
real(dp) :: x

x = 1.0d0 - r_sqr*c

if (x >= 0.0d0) then
    gtheta(1) = x*x
    gtheta(2) = 4.0d0*x*c ! coeff for dtheta/dr(:) = CR(:)
else
    gtheta(:) = 0.0d0
end if

end function gtheta
!-------------------------------------------------------------------------------
END MODULE FACTS
!-------------------------------------------------------------------------------
