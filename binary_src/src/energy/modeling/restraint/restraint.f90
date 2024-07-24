!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/modeling/restraint/restraint.f90
!
! Description: Subroutines related to generating restraint data types and 
!              calculating restraint energy. Subroutines for updating restraints are
!              also included.
!
!-------------------------------------------------------------------------------
MODULE RESTRAINT

use globals
use in_out_input, only: read_pdblist
use in_out_utils, only: find_atom_idx, find_ligand_atom_idx
use in_out_vars, only: infile_mol2_topo, read_het, n_mol2_top
use logger, only: terminate_with_error, log_p
use string, only: parse_longstring
use mathfunctions
use geometry, only: calc_ang_and_grad, calc_tor_and_grad
use sort, only: sort2
!
use energy_vars, only: energy_type, R, i_R, ii_R, ii_L, rsr_file, use_update_rsr, &
                       update_NO_rsr, update_lig_rsr, use_rsr, rsr_pdb, &
                       update_rsr_dcut, update_rsr_sig, update_lig_rsr_sig, &
                       infile_update_rsr_pdblist, max_energy, &
                       docking_grid_param_type, docking_grid_type
use energy_utils, only: protein_to_R
!
use ligdock_E_utils, only: copy_to_grid, calc_E_using_grid

implicit none
save
private

! Parameters for setting up arrays
integer, parameter :: max_rsrcolumn = 300   ! Maximum # of parameters of a restraint
integer, parameter :: max_rsrlinelen = 3000 ! Maximum length of a line in restraint file
integer, parameter :: max_regparm = 100     ! Maximum # of parameters for regular restraints
integer, parameter :: max_longparm = 400    ! Maximum # of parameters for long restraints


integer :: n_rsr, n_rsr0         ! # of restraints
integer :: n_rsrreg, n_rsrreg0   ! # of restraints, for regular restraints
integer :: n_rsrlong             ! # of restraints, for long restraints

! Restraint type and their parameter types
type restraint_type
    integer :: form      ! Functional form
    integer :: modality  ! Modality of the restraint function
    integer :: feature   ! Feature type (distance, )
    integer :: grp       ! Restraint group
    integer :: nparm     ! # parameter of the restraint
    integer :: natm      ! # atoms involved to the restraint
    integer :: atmno(8)  ! Atom number
    logical :: selected  ! To use or not
    logical :: disabled  ! To use or not
    integer :: parm_id(2)! Regular or Long type
    logical :: is_ang    ! Is angle (cyclic) or not
    logical :: is_het    ! Is ligand restraint
    real(dp) :: E        ! Restraint energy value
    real(dp), allocatable :: g(:,:) ! gradient
    logical :: active    ! turn on for MELD
end type restraint_type

!TODO comment
type regular_parm_type
    real(dp) :: parm(max_regparm)    ! restraint parameters
    real(dp) :: spl_ypp(max_regparm) ! pre-calculated cubic_spline 2nd derivative
end type regular_parm_type
!TODO comment
type long_parm_type
    real(dp) :: parm(max_longparm)   ! restraint parameters
    real(dp) :: spl_ypp(max_longparm)! pre-calculated cubic_spline 2nd derivative
end type long_parm_type

! General info of all the restraints is stored in type rsr,
! while parameters are stored separately depending on the restraint types
type(restraint_type), allocatable :: rsr(:)
type(restraint_type), allocatable :: rsr0(:)
type(regular_parm_type), allocatable :: reg_parm(:)
type(regular_parm_type), allocatable :: reg_parm0(:)
type(long_parm_type), allocatable :: long_parm(:)
logical :: rsr_allocated

! For update_rsr
integer :: n_update_rsr_pdb
character(len=len_fname) :: update_rsr_pdblist(max_pdbfile)

public :: initialize_rsr
public :: finalize_rsr
public :: calc_rsr_energy
public :: calc_rsr_energy_by_meld
public :: select_rsr
public :: selected_res_rsr
public :: n_update_rsr_pdb
public :: update_rsr_pdblist
public :: protein_to_topolid
public :: read_update_rsr_pdb
public :: update_rsr
public :: clear_updated_rsr
public :: build_rsr_backup
public :: construct_ligdock_rsr_grid
public :: calc_ligdock_rsr_using_grid

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_rsr(molecule)
!-------------------------------------------------------------------------------
! Initialize restraint function, for MODELING
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: molecule

! For the case using only segment fixing
if (trim(rsr_file) == '') return

!Read rsr variables
call read_rsr_input(rsr_file)

if (n_rsr > 0) then
    ! Match restraint index with program index
    call rsrid_to_topolid(molecule)
    ! Precalculate spline related variable
    call build_spline()
    ! Make a backup for updates
    if (use_update_rsr) then
        call build_rsr_backup()
    end if
end if

n_update_rsr_pdb = 0
if (use_update_rsr) then
    if (trim(infile_update_rsr_pdblist) /= '') then
        call read_pdblist(infile_update_rsr_pdblist, &
                                  update_rsr_pdblist, n_update_rsr_pdb)
        if (n_update_rsr_pdb == 0) then
            call terminate_with_error("ERROR: Empty update_rsr_pdblist.")
        end if
    end if
else
    if (n_rsr == 0) then
        use_rsr = .false.
        call log_p('- Energy: Restraint energy deactivated, n_rsr=0.', me=me, level=10)
    end if
end if

end subroutine initialize_rsr
!-------------------------------------------------------------------------------
subroutine finalize_rsr()
!-------------------------------------------------------------------------------
if (allocated(rsr)) then
    deallocate(rsr)
    deallocate(reg_parm, long_parm)
end if
if (use_update_rsr) then
    deallocate(rsr0)
    deallocate(reg_parm0)
end if

end subroutine finalize_rsr
!-------------------------------------------------------------------------------
subroutine build_rsr_backup()
!-------------------------------------------------------------------------------
n_rsr0 = n_rsr
n_rsrreg0 = n_rsrreg
allocate(rsr0(n_rsr))
allocate(reg_parm0(n_rsrreg))
rsr0(:) = rsr(:)
reg_parm0(:) = reg_parm(:)

end subroutine build_rsr_backup
!-------------------------------------------------------------------------------
subroutine read_rsr_input(file_name)
!-------------------------------------------------------------------------------
! Read the restraint input file
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: file_name
integer :: f_unit, ioerror, num_word, openstat
integer :: i, i_rsr, i_reg, i_long, nparm
integer :: type_temp(6), atom_temp(8)
integer :: atm1, atm2
character(len=len_fname) :: word(max_rsrcolumn)
character(len=max_rsrlinelen) :: line
real(dp) :: parm_temp(max_longparm)

!1. First get max_rsrparm for allocation
f_unit = 39
open(f_unit, file = trim(file_name), iostat=openstat)
if (openstat > 0) then
    call terminate_with_error('Terminate with error: No restraint file found. Check {rsr_file}.')
end if
  
n_rsr = 0
n_rsrlong = 0
n_rsrreg = 0
do
    read(f_unit,"(A3000)", iostat = ioerror) line
    if (ioerror < 0) exit

    call parse_longstring(line, num_word, word, max_rsrlinelen)

    if (line(1:1) /= 'R') cycle
    read(word(7),"(I10)") nparm
    if (nparm > max_regparm) then
        n_rsrlong = n_rsrlong + 1
    else
        n_rsrreg = n_rsrreg + 1
    end if
    n_rsr = n_rsr + 1
end do

if (n_rsr == 0) then
    call log_p('No restraint was found from the file. Please check.', me=me)
    return
else
    write(log_msg,"(I6,A)") n_rsr, ' Restraints found from {rsr_file}.'
    call log_p(log_msg, me=me, level=10)
    write(log_msg,"(I6,A)") n_rsrreg, ' restraints are set to Regular type with maximum 100 parameters.'
    call log_p(log_msg, me=me, level=20)
    write(log_msg,"(I6,A)") n_rsrlong, ' restraints are set to Long type with maximum 400 parameters.'
    call log_p(log_msg, me=me, level=20)
end if

!2. Allocate array sizes with given input file
allocate(rsr(n_rsr))
allocate(reg_parm(n_rsrreg), long_parm(n_rsrlong))
rsr_allocated = .true.
close(f_unit)

!3. Then read input file again and assign parameters
open(f_unit, file = trim(file_name))

i_rsr = 0
i_reg = 0
i_long = 0
do
    read(f_unit,"(A3000)", iostat = ioerror) line
    if (ioerror < 0) exit

    call parse_longstring(line, num_word, word, max_rsrlinelen)

    if (line(1:1) /= 'R') cycle
    if (num_word < 2) exit

    i_rsr = i_rsr + 1
    ! Basic information of restraint
    do i = 1, 6
        read(word(i+1),"(I6)") type_temp(i)
    end do
    nparm = type_temp(6)

    ! Atom index
    do i = 1, type_temp(5)
        read(word(i+8),"(I6)") atom_temp(i)
    end do
    ! Restraint parameter
    do i = 1, nparm
        read(word(i+8+type_temp(5)),"(F10.4)") parm_temp(i)
    end do
     
    rsr(i_rsr)%form = type_temp(1)
    rsr(i_rsr)%modality = type_temp(2)
    rsr(i_rsr)%feature = type_temp(3)
    rsr(i_rsr)%natm = type_temp(5)
    rsr(i_rsr)%nparm = type_temp(6)
    rsr(i_rsr)%atmno(1:type_temp(5)) = atom_temp(1:type_temp(5))
    rsr(i_rsr)%selected = .true.
    rsr(i_rsr)%disabled = .false.
    allocate(rsr(i_rsr)%g(3,rsr(i_rsr)%natm))

    if (type_temp(4) == 9) then        ! CA-CA
        rsr(i_rsr)%grp = 1
    else if (type_temp(4) == 10) then  ! N-O
        rsr(i_rsr)%grp = 2
    else if (type_temp(4) == 23) then  ! MNCH-SDCH
        rsr(i_rsr)%grp = 3
    else if (type_temp(4) == 26) then  ! SDCH-SDCH
        rsr(i_rsr)%grp = 4
    else if (type_temp(4) == 18) then  ! Disulfide Bond
        rsr(i_rsr)%grp = 4
    else if (type_temp(4) == 22) then  ! HETATM
        rsr(i_rsr)%grp = 5
    else
        rsr(i_rsr)%selected = .false.
        rsr(i_rsr)%disabled = .true.
    end if

    ! Check angle features
    if (rsr(i_rsr)%feature == 2 .or. rsr(i_rsr)%feature == 3 .or. &
        rsr(i_rsr)%feature == 4) then
        rsr(i_rsr)%is_ang = .true.
    else
        rsr(i_rsr)%is_ang = .false.
    end if

    ! Assign to long parameter
    if (nparm > max_regparm) then
        i_long = i_long + 1
        rsr(i_rsr)%parm_id(1) = 2
        rsr(i_rsr)%parm_id(2) = i_long
        long_parm(i_long)%parm(1:nparm) = parm_temp(1:nparm)
    ! Assign to regular parameter
    else
        i_reg = i_reg + 1
        rsr(i_rsr)%parm_id(1) = 1
        rsr(i_rsr)%parm_id(2) = i_reg
        reg_parm(i_reg)%parm(1:nparm) = parm_temp(1:nparm)
    end if
     
    !Check disulfide bond
    if (type_temp(4) == 18) then
        atm1 = atom_temp(1)
        atm2 = atom_temp(2)
    end if
end do

close(f_unit)

end subroutine read_rsr_input
!-------------------------------------------------------------------------------
subroutine rsrid_to_topolid(molecule)
!-------------------------------------------------------------------------------
! Map atom index of rsrs into atom index of the program
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: molecule
integer :: max_rsr_atm, max_rsr_atm_prev
integer :: f_unit, ioerror, openstat
integer :: i_atm, i_rsr, res_no, het_no, atm_no, ref_res_no, resatm_no
integer :: lig_no, i_mol2
integer :: res_no_pdb, res_no_prv
logical :: is_atom, is_ligand
character(len=6) :: atom_error_mode
character(len=4) :: atmtype, het_name
character(len=max_rsrlinelen) :: line
type(restraint_type) :: rsr_i
! Map of (atmno,resno) <-> rsr atom index
integer, allocatable :: map_atomid(:,:), map(:,:)

max_rsr_atm = tn%atom
allocate(map_atomid(2,max_rsr_atm))

! Parse pdb file
i_atm = 0
res_no_prv = -9999
res_no = 0
het_no = 0
lig_no = 0
f_unit = 25
atom_error_mode = 'ignore'

open(f_unit, file = trim(rsr_pdb), iostat = openstat)
if (openstat > 0) then
    call terminate_with_error('Terminate with error: No rsr_pdb file found. Check {rsr_pdb}.')
end if

do
    read(f_unit,"(A1000)", iostat = ioerror) line
    if (ioerror < 0) exit
    if (line(1:4) == 'ATOM') then
        is_atom = .true.
        is_ligand = .false.
    else
        if (read_het .and. line(1:6) == 'HETATM') then
            is_atom = .false.
            is_ligand = .false.
            het_name = adjustl(line(18:20))
            do i_mol2 = 1, n_mol2_top
                if (het_name == trim(infile_mol2_topo(2,i_mol2))) then
                    is_ligand = .true.
                end if
            end do
        else
            cycle
        end if
    end if
     
    i_atm = i_atm + 1

    if (line(13:13) == ' ') then
        read(line(14:16),"(A4)") atmtype
    else !If atmtype is written like "1HH1"
        read(line(13:16),"(A4)") atmtype
        if (.not. is_ligand) then !TODO: To match mol2, do not change name
            atmtype = trim(atmtype(2:4))//atmtype(1:1)
        end if
    end if
    read(line(23:26),"(I4)") res_no_pdb
    if (res_no_pdb /= res_no_prv) then
        if (is_atom) then
            res_no = res_no + 1
        else if (is_ligand) then
            lig_no = lig_no + 1
        else
            het_no = het_no + 1
        end if
    end if
    res_no_prv = res_no_pdb

    if (i_atm > max_rsr_atm) then
        ! reallocate map_atomid by 2*(max_rsr_atm+tn%atom)
        max_rsr_atm_prev = max_rsr_atm
        max_rsr_atm = max_rsr_atm + tn%atom
        allocate(map(2,max_rsr_atm))
        map(:, 1:max_rsr_atm_prev) = map_atomid(:, 1:max_rsr_atm_prev)
        !
        call move_alloc(map, map_atomid)
    end if

    if (is_atom) then
        ref_res_no = molecule%residue(res_no)%res_type
        call find_atom_idx(res_no, ref_res_no, atmtype, resatm_no, atom_error_mode)
        map_atomid(1,i_atm) = res_no
        map_atomid(2,i_atm) = resatm_no
    else if (is_ligand) then
        ref_res_no = molecule%ligand(lig_no)%lig_type
        call find_ligand_atom_idx(lig_no, ref_res_no, atmtype, resatm_no, atom_error_mode)
        map_atomid(1,i_atm) = tn%nonlig + lig_no
        map_atomid(2,i_atm) = resatm_no
    else
        ref_res_no = molecule%hetmol(het_no)%res_type
        call find_atom_idx(het_no, ref_res_no, atmtype, resatm_no, atom_error_mode)
        map_atomid(1,i_atm) = tn%stdres + het_no
        map_atomid(2,i_atm) = resatm_no
    end if
end do

close(f_unit)

!Reset atom indices
do i_rsr = 1, n_rsr
    rsr_i = rsr(i_rsr)
    do i_atm = 1, rsr_i%natm
        res_no = map_atomid(1,rsr_i%atmno(i_atm))
        atm_no = map_atomid(2,rsr_i%atmno(i_atm))
        if (res_no <= tn%nonlig) then
            ! Is initial value of ii_R(:,:) set to -100?
            if (ii_R(atm_no,res_no) == -100) then
                rsr(i_rsr)%selected = .false.
                rsr(i_rsr)%disabled = .true.
                exit
            end if
            rsr(i_rsr)%atmno(i_atm) = ii_R(atm_no,res_no)
        else
            rsr(i_rsr)%atmno(i_atm) = ii_L(atm_no,res_no-tn%nonlig)
        end if
    end do
end do

deallocate(map_atomid)

end subroutine rsrid_to_topolid
!-------------------------------------------------------------------------------
subroutine build_spline()
!-------------------------------------------------------------------------------
! Prepare spline array spl_ypp (second derivative)
!-------------------------------------------------------------------------------
integer :: i, k, i_rsr, n
integer :: i_p(2), nparm
real(dp) :: c(max_longparm), x(max_longparm), c_pp(max_longparm)
real(dp) :: p, sig, u(max_longparm), c_p0, c_pn
real(dp) :: parm(max_longparm)

do i_rsr = 1, n_rsr
    if (rsr(i_rsr)%form /= 10) cycle ! 10: cubic spline type

    n = rsr(i_rsr)%modality ! No. bin points
    i_p(:) = rsr(i_rsr)%parm_id(:)
    parm(:) = 0.0d0
    nparm = rsr(i_rsr)%nparm

    if (i_p(1) == 1) then ! Regular
        parm(1:nparm) = reg_parm(i_p(2))%parm(1:nparm)
    else ! Long
        parm(1:nparm) = long_parm(i_p(2))%parm(1:nparm)
    end if

    ! Build second derivative recursively
    do i = 1, n
        x(i) = parm(2) + (i-1)*parm(4)
    end do

    c(1:n) = parm(7:6+n)
    c_p0 = parm(5)  !given first derivative
    c_pn = parm(6)  !given first derivative

    c_pp(:) = 0.0d0
    ! with a specified first derivative c_p0,
    c_pp(1) = -0.5d0
    u(1) = 3.0d0/(x(2)-x(1)) * ((c(2)-c(1))/(x(2)-x(1))-c_p0)

    do i = 2, n-1
        sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
        p = sig*c_pp(i-1) + 2.0d0
        c_pp(i) = (sig - 1.0d0)/p
        u(i) = (6.0d0*((c(i+1) - c(i))/(x(i+1)-x(i)) - (c(i)-c(i-1))/(x(i)-x(i-1))) / &
             (x(i+1)-x(i-1)) - sig*u(i-1)) / p
    end do

    ! with a specified first derivative c_pn,
    c_pp(n) = 0.5d0
    u(n) = 3.0d0/(x(n)-x(n-1)) * (c_pn-(c(n)-c(n-1))/(x(n)-x(n-1)))
    do k = n-1, 1, -1
        c_pp(k) = c_pp(k)*c_pp(k+1) + u(k)
    end do

    ! Save obtained spl_ypp values
    if (i_p(1) == 1) then
        reg_parm(i_p(2))%spl_ypp(1:n) = c_pp(1:n)
    else
        long_parm(i_p(2))%spl_ypp(1:n) = c_pp(1:n)
    end if
end do

end subroutine build_spline
!-------------------------------------------------------------------------------
subroutine find_spline_value(f, nspl, x_min, interval, c, c_pp, fdiff, c_j, c_jpp, a)
!-------------------------------------------------------------------------------
! Find spline variables from current feature
!-------------------------------------------------------------------------------
integer, intent(in) :: nspl
real(dp), intent(in) :: x_min, interval, f, c(nspl), c_pp(nspl)
real(dp), intent(out) :: fdiff, c_j(2), c_jpp(2), a
integer :: i, n_j
real(dp) :: x(nspl)

! Fill bin points
do i = 1, nspl
    x(i) = x_min + (i-1)*interval
end do

! Find bin location, and make sure location is in range
n_j = max(1,int((f - x_min)/interval) + 1)
n_j = min(nspl,n_j)

c_j(:) = c(n_j:n_j+1) ! Function value
c_jpp(:) = c_pp(n_j:n_j+1)       ! Second derivative
fdiff = x(n_j+1) - f             ! Value difference from prior bin
a = fdiff / (x(n_j+1) - x(n_j))

end subroutine find_spline_value
!-------------------------------------------------------------------------------
subroutine get_rsr_feature(feature_val, multi_val, dfdr, calc_g)
!-------------------------------------------------------------------------------
! Extract information from current protein coordinate R at once
! This subroutine provides information to [calc_restraint]
!-------------------------------------------------------------------------------
real(dp), intent(out) :: feature_val(n_rsr), multi_val(8,n_rsr), dfdr(3,8,n_rsr)
logical, intent(in) :: calc_g
type(restraint_type) :: rsr_i
integer :: i, i_rsr, i_p(2)
integer :: atom_id(8)
real(dp) :: Rtemp(3,4), r_dist(3), arg, parm(max_longparm)
  
feature_val(:) = 0.0d0
multi_val(:,:) = 0.0d0
dfdr(:,:,:) = 0.0d0

do i_rsr = 1, n_rsr
    rsr_i = rsr(i_rsr)
    if (.not. rsr_i%selected) cycle ! Evaluate only those selected (controled by [select_rsr])

    ! Check restraint type
    i_p(:) = rsr(i_rsr)%parm_id(:)
    if (i_p(1) == 1) then
        parm(1:rsr_i%nparm) = reg_parm(i_p(2))%parm(1:rsr_i%nparm)
    else
        parm(1:rsr_i%nparm) = long_parm(i_p(2))%parm(1:rsr_i%nparm)
    end if
    atom_id(:) = rsr_i%atmno(:)

    ! Select feature to be extracted
    select case(rsr_i%feature)

    case(1) !atom-atom distance
        do i = 1, 2
            Rtemp(:,i) = R(:,atom_id(i))
        end do

        r_dist(:) = Rtemp(:,2) - Rtemp(:,1)
        arg = sqrt(dot_product(r_dist(:),r_dist(:)))
        feature_val(i_rsr) = arg

        ! calculate dfdr for feature distance
        if (calc_g) then
            dfdr(:,1,i_rsr) = -r_dist(:)/arg
            dfdr(:,2,i_rsr) = -dfdr(:,1,i_rsr)
        end if

    case(2) !angle
        do i = 1, 3
            Rtemp(:,i) = R(:,atom_id(i))
        end do
        call calc_ang_and_grad(Rtemp(:,1:3), feature_val(i_rsr), dfdr(:,1:3,i_rsr), calc_g)
  
    case(3) !dihedral angle
        do i = 1, 4
            Rtemp(:,i) = R(:,atom_id(i))
        end do
        call calc_tor_and_grad(Rtemp(:,:), feature_val(i_rsr), dfdr(:,1:4,i_rsr), calc_g)

    case(4) !pair of dihedral angles
        do i = 1, 4
            Rtemp(:,i) = R(:,atom_id(i))
        end do
        call calc_tor_and_grad(Rtemp(:,:), multi_val(1,i_rsr), dfdr(:,1:4,i_rsr), calc_g)

        do i = 1, 4
            Rtemp(:,i) = R(:,atom_id(4+i))
        end do
        call calc_tor_and_grad(Rtemp(:,:), multi_val(2,i_rsr), dfdr(:,5:8,i_rsr), calc_g)

    case(5) ! flat-bottomed atom-atom distance ___/
        do i = 1, 2
            Rtemp(:,i) = R(:,atom_id(i))
        end do

        r_dist(:) = Rtemp(:,2) - Rtemp(:,1)
        arg = sqrt(dot_product(r_dist(:),r_dist(:)))
        if (arg < parm(1)) then
            feature_val(i_rsr) = parm(1) ! to make flat-bottom in function
        else
            feature_val(i_rsr) = arg
        end if

        ! calculate dfdr for feature distance
        if (calc_g) then
            dfdr(:,1,i_rsr) = -r_dist(:)/arg
            dfdr(:,2,i_rsr) = -dfdr(:,1,i_rsr)
        end if

    case(9) ! Added : coordinate
        r_dist(:) = R(:,rsr_i%atmno(1)) - parm(3:5)
        arg = sqrt(dot_product(r_dist(:),r_dist(:)))
        feature_val(i_rsr) = arg
        
        if (calc_g) then
            if (arg < small_real) then
                dfdr(:,1,i_rsr) = 0.0d0
            else
                dfdr(:,1,i_rsr) = r_dist(:)/arg
            end if
        end if

    case(10) ! Added : z-coordinate
        r_dist(1:2) = 0.0d0
        r_dist(3) = R(3,rsr_i%atmno(1)) - parm(3)
        arg = abs(r_dist(3))
        feature_val(i_rsr) = arg
        
        if (calc_g) then
            if (arg < small_real) then
                dfdr(:,1,i_rsr) = 0.0d0
            else
                dfdr(:,1,i_rsr) = r_dist(:)/arg
            end if
        end if

    case(11) ! Added : block for specific direction
        r_dist(:) = R(:,rsr_i%atmno(1)) - parm(3:5)
        arg = dot_product(r_dist, parm(6:8))
        if (arg < parm(1)) arg = parm(1)
        feature_val(i_rsr) = arg
        
        if (calc_g) then
            if (arg < small_real) then
                dfdr(:,1,i_rsr) = 0.0d0
            else
                dfdr(:,1,i_rsr) = parm(6:8)
            end if
        end if
    end select
end do

end subroutine get_rsr_feature
!-------------------------------------------------------------------------------
subroutine calc_rsr_energy(ff, g, calc_g)
!-------------------------------------------------------------------------------
! Main subroutine for restraints and gradients
! This subroutine first calls [get_rsr_feature] to extract all the features from
! current coordinate R at once, 
! then evaluate restraint function by restraints.
!-------------------------------------------------------------------------------
real(dp), intent(out) :: ff(1:5)
real(dp), intent(out) :: g(3,tn%atom,5)
logical, intent(in) :: calc_g
type(restraint_type) :: rsr_i
integer :: i_rsr, i, j, i_atm
real(dp) :: feature_val(n_rsr), multi_val(8,n_rsr), dfdr(3,8,n_rsr), dcdf
real(dp) :: x_i, xm_i(8)
real(dp) :: e_val, scale, g_add(3)
! Multiple-binormal
integer :: n_ss
real(dp) :: dcdf_m(0:8)
real(dp) :: w_b(6), rho_b(6), sigma_b(2,6), xm0(2,6)
! Spline variable
integer :: nspl
real(dp) :: spl_ypp(max_longparm), argv(2)
real(dp) :: x_max, x_min, y_max, y_min, interval, dx_min, dx_max
real(dp) :: spl_fdiff, spl_y_j(2), spl_y_jpp(2), spl_a
real(dp) :: parm(max_longparm)
! Multiple gaussian
integer :: n_g
real(dp) :: argm(3)

ff(:) = 0.0d0
g(:,:,:) = 0.0d0

! First extract all the features required
call get_rsr_feature(feature_val, multi_val, dfdr, calc_g)

! Then evaluate restraints
do i_rsr = 1, n_rsr
    rsr_i = rsr(i_rsr)
    if (.not. rsr_i%selected) cycle ! Evaluate only those selected (controled by [select_rsr])

    x_i = feature_val(i_rsr)
    xm_i(1:8) = multi_val(1:8,i_rsr)
    parm(:) = 0.0d0

    ! Check restraint parameter type (long or regular)
    if (rsr_i%parm_id(1) == 1) then
        parm(1:rsr_i%nparm) = reg_parm(rsr_i%parm_id(2))%parm(1:rsr_i%nparm)
    else
        parm(1:rsr_i%nparm) = long_parm(rsr_i%parm_id(2))%parm(1:rsr_i%nparm)
    end if

    ! Calculate restraint function depending on its functional form
    select case(rsr_i%form)

    case(3) !Single gaussian
        argv(:) = RT*single_gaussian(x_i, parm(1), parm(2), rsr_i%is_ang, calc_g)
        e_val = argv(1)
        if (calc_g) dcdf = argv(2)

    case(4) ! Multiple gaussian
        n_g = rsr_i%modality
        argv(:) = RT*multiple_gaussian(n_g, x_i, parm(n_g+1:2*n_g), &
             parm(1:n_g), parm(2*n_g+1:3*n_g), rsr_i%is_ang, .true., calc_g)

        e_val = argv(1)
        if (calc_g) dcdf = argv(2)

    case(5) ! Single lorentzian
        argv(:) = lorentzian(x_i, parm(1), parm(2), calc_g)
        e_val = argv(1)
        if (calc_g) dcdf = argv(2)

    case(7) ! Cosine potential
        argv(:) = cosine_potential(x_i, parm(1), parm(2), rsr_i%modality, calc_g)
        e_val = argv(1)
        if (calc_g) dcdf = argv(2)
        
    case(9) ! Multiple binormal
        ! n_ss : Number of Secondary structure type
        n_ss = rsr_i%modality
        w_b(1:n_ss) = parm(1:n_ss)
        rho_b(1:n_ss) = parm(5*n_ss+1:6*n_ss)

        ! Rearrange parameters
        do i = 1, 2*n_ss, 2 ! For phi angle
            j = int((i+1)/2)
            sigma_b(1,j) = parm(3*n_ss+i)
            xm0(1,j) = parm(n_ss+i)
        end do
        do i = 2, 2*n_ss, 2 ! For psi angle
            j = int(i/2)
            sigma_b(2,j) = parm(3*n_ss+i)
            xm0(2,j) = parm(n_ss+i)
        end do

        argm(:) = RT*multiple_binormal(n_ss, xm_i(1:2), xm0(1:2,1:n_ss), w_b(1:n_ss), &
             sigma_b(1:2,1:n_ss), rho_b(1:n_ss), .true., calc_g)
        e_val = argm(1)

        if (calc_g) dcdf_m(1:2) = argm(2:3)

    case(10) ! Cubic spline
        nspl = rsr_i%modality
        scale = parm(1)
        x_min = parm(2)
        x_max = parm(3)
        interval = parm(4)
        dx_min = parm(5)
        dx_max = parm(6)
        y_min = parm(7)
        y_max = parm(6+nspl)
        if (rsr_i%parm_id(1) == 1) then
            spl_ypp(1:rsr_i%nparm) = reg_parm(rsr_i%parm_id(2))%spl_ypp(1:rsr_i%nparm)
        else
            spl_ypp(1:rsr_i%nparm) = long_parm(rsr_i%parm_id(2))%spl_ypp(1:rsr_i%nparm)
        end if

        ! Check the location of feature in spline table
        ! and interpolate linearly if out of bound
        if (x_i < x_min) then
            e_val = scale*(y_min - (x_min-x_i)*dx_min)
        else if (x_i > x_max) then
            e_val = scale*(y_max + (x_i-x_max)*dx_max)
        else
            call find_spline_value(x_i, nspl, x_min, interval, parm(7:6+nspl), &
                spl_ypp(1:nspl), spl_fdiff, spl_y_j, spl_y_jpp, spl_a)
            argv(:) = scale*cubic_spline(interval, spl_y_j, spl_y_jpp, spl_a, calc_g)
            e_val = argv(1)
        end if

        if (calc_g) then
            if (x_i < x_min) then
                dcdf = dx_min
            else if (x_i > x_max) then
                dcdf = dx_max
            else
                dcdf = argv(2)
            end if
        end if

    ! User defined functional forms
    case(11) ! Sigmoidal
        argv(:) = sigmoidal1(x_i, parm(2), parm(2)+parm(3), calc_g)
        e_val = parm(1)*argv(1)
        dcdf = parm(1)*argv(2)

    case(12) ! Cosine contact
        ! Parm : Scale / Center value / sigma
        argv(:) = cosine_well(x_i, parm(2), parm(3), calc_g)
        e_val = parm(1)*argv(1)
        dcdf = parm(1)*argv(2)

    end select

    rsr(i_rsr)%E = e_val
    ! Sum over energy
    ff(rsr_i%grp) = ff(rsr_i%grp) + e_val

    if (.not. calc_g) cycle
    ! Gradient calculation : by multiplying dcdf & dfdr
    if (rsr_i%form == 9) then ! Exception for bi-function
        do i = 1, 4
            i_atm = rsr_i%atmno(i)
            g_add(:) = dfdr(:,i,i_rsr)*dcdf_m(1)
            rsr_i%g(:,i) = g_add(:)
            g(:,i_atm,rsr_i%grp) = g(:,i_atm,rsr_i%grp) + g_add(:)
           
            i_atm = rsr_i%atmno(i+4)
            g_add(:) = dfdr(:,i+4,i_rsr)*dcdf_m(2)
            rsr_i%g(:,i) = rsr_i%g(:,i) + g_add(:)
            g(:,i_atm,rsr_i%grp) = g(:,i_atm,rsr_i%grp) + g_add(:)
        end do
    else
        do i = 1, rsr_i%natm
            if (rsr_i%atmno(i) /= 0) then
                i_atm = rsr_i%atmno(i)
                rsr_i%g(:,i) = dcdf*dfdr(:,i,i_rsr)
                g(:,i_atm,rsr_i%grp) = g(:,i_atm,rsr_i%grp) + dcdf*dfdr(:,i,i_rsr)
            end if
        end do
    end if
    rsr(i_rsr)%g = rsr_i%g
end do

end subroutine calc_rsr_energy
!-------------------------------------------------------------------------------
subroutine calc_rsr_energy_by_meld(ff, g, n_features, n_grp, active_fraction, calc_g)
!-------------------------------------------------------------------------------
! calculate energy only for restraints defined as active by meld procedure
!-------------------------------------------------------------------------------
real(dp), intent(out) :: ff(1:5)
real(dp), intent(out) :: g(3,tn%atom,5)
real(dp), intent(in) :: active_fraction(:,:)
integer, intent(in) :: n_features
integer, intent(in) :: n_grp
logical, intent(in) :: calc_g
type(restraint_type) :: rsr_i
integer :: i_rsr, i, i_atm

call calc_rsr_energy(ff, g, calc_g)
call update_restraints_as_meld(n_features, n_grp, active_fraction(1:n_grp,1:n_features))
ff(:) = 0.0d0
g(:,:,:) = 0.0d0
!
do i_rsr = 1, n_rsr
    rsr_i = rsr(i_rsr)
    if (.not. rsr_i%selected) cycle
    if (.not. rsr_i%active) cycle
    ff(rsr_i%grp) = ff(rsr_i%grp) + rsr_i%E
end do

if (calc_g) then
    do i_rsr = 1, n_rsr
        rsr_i = rsr(i_rsr)
        if (.not. rsr_i%selected) cycle
        if (.not. rsr_i%active) cycle
        !
        if (rsr_i%form == 9) then
            do i = 1, 4
                i_atm = rsr_i%atmno(i)
                g(:,i_atm,rsr_i%grp) = g(:,i_atm,rsr_i%grp) + rsr_i%g(:,i)
            end do
        else
            do i = 1, rsr_i%natm
                i_atm = rsr_i%atmno(i)
                if (i_atm /= 0) then
                    g(:,i_atm,rsr_i%grp) = g(:,i_atm,rsr_i%grp) + rsr_i%g(:,i)
                end if
            end do
        end if
    end do
end if

end subroutine calc_rsr_energy_by_meld
!-------------------------------------------------------------------------------
subroutine build_spatial_rsr(molecule, rsrres, scale, mode)
!-------------------------------------------------------------------------------
! Build additional restraints (from another route than restraint input file)
! on intra-segment residue pairs given as input rsrres
! To consider Ca-Ca contacts only ('ca') or additionally append N-O contacts ('')
! can be adjusted by setting input 'mode'
! This subroutine can be called without calling [initialize_rsr]
! Used for segment containing ULR modeling.
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: molecule
logical, intent(in) :: rsrres(molecule%n_res,molecule%n_res)
real(dp), intent(in) :: scale
character(len=2), intent(in) :: mode
integer :: i_res, j_res
integer :: i_cand, n_cand, candatm_id(2,molecule%n_res*molecule%n_res*2)
real(dp) :: R1(3), R2(3), R12(3)
type(restraint_type), allocatable :: rsr_cp(:)
type(regular_parm_type), allocatable :: reg_parm_cp(:)
  
i_cand = 0
candatm_id(:,:) = 0

! First check candidates
do i_res = 1, molecule%n_res-1
    do j_res = i_res+1, molecule%n_res
        if (.not. rsrres(i_res,j_res)) cycle
        
        ! Ca-Ca
        i_cand = i_cand + 1
        candatm_id(1,i_cand) = res_index(i_res)%Ca_id(1)
        candatm_id(2,i_cand) = res_index(j_res)%Ca_id(1)

        if (mode == '') then
            !N-O
            i_cand = i_cand + 1
            candatm_id(1,i_cand) = res_index(i_res)%bb_id(1,1)
            candatm_id(2,i_cand) = res_index(j_res)%bb_id(4,1)

            !O-N
            i_cand = i_cand + 1
            candatm_id(1,i_cand) = res_index(i_res)%bb_id(4,1)
            candatm_id(2,i_cand) = res_index(j_res)%bb_id(1,1)
        end if
    end do
end do

n_cand = i_cand

! Reallocate arrays if defined previously
allocate(rsr_cp(n_rsr), reg_parm_cp(n_rsrreg))
if (rsr_allocated) then
    rsr_cp(:) = rsr(:)
    reg_parm_cp(:) = reg_parm(:)
    deallocate(rsr, reg_parm)
end if
  
write(log_msg,"(A,I6,A8,I6)") '   Additional restraints appended:', n_cand, ', Total:', n_rsr+n_cand
call log_p(log_msg, me=me, level=10)
  
allocate(rsr(n_rsr+n_cand))
allocate(reg_parm(n_rsrreg+n_cand))
rsr_allocated = .true.

! Unify previously assigned restraints
rsr(1:n_rsr) = rsr_cp(1:n_rsr)
reg_parm(1:n_rsrreg) = reg_parm_cp(1:n_rsrreg)
deallocate(rsr_cp, reg_parm_cp)

! Assign parameters for restraint candidate pairs
! Harmonic restraints are applied uniformly
do i_cand = 1, n_cand
    n_rsr = n_rsr + 1
    rsr(n_rsr)%form = 3
    rsr(n_rsr)%modality = 1
    rsr(n_rsr)%feature = 1
    !rsr(n_rsr)%grp = 9
    rsr(n_rsr)%grp = 1
    rsr(n_rsr)%natm = 2
    rsr(n_rsr)%nparm = 2
    rsr(n_rsr)%atmno(1) = candatm_id(1,i_cand)
    rsr(n_rsr)%atmno(2) = candatm_id(2,i_cand)
    rsr(n_rsr)%selected = .true.
    rsr(n_rsr)%disabled = .false.

    n_rsrreg = n_rsrreg + 1
    rsr(n_rsr)%parm_id(1) = 1
    rsr(n_rsr)%parm_id(2) = n_rsrreg

    R1(:) = R(:,candatm_id(1,i_cand))
    R2(:) = R(:,candatm_id(2,i_cand))
    R12(:) = R1(:) - R2(:)

    reg_parm(n_rsrreg)%parm(1) = sqrt(dot_product(R12(:),R12(:)))
    reg_parm(n_rsrreg)%parm(2) = 1.0d0/sqrt(scale)
end do

if (use_update_rsr) then
    deallocate(rsr0)
    deallocate(reg_parm0)

    call build_rsr_backup()
end if

end subroutine build_spatial_rsr
!-------------------------------------------------------------------------------
subroutine select_rsr(span_range)
!-------------------------------------------------------------------------------
! Restraint selection by VTFM (Variable Target Function Method)
!-------------------------------------------------------------------------------
integer, intent(in) :: span_range
type(restraint_type) :: rsr_i
integer :: i_rsr, i_atm
integer :: resno, resdiff, max_resno, min_resno

! Set local restraints
if (.not. use_rsr) return
if (span_range < 0) then
    rsr(:)%selected = .not. rsr(:)%disabled
    return
end if

rsr(:)%selected = .false.

do i_rsr = 1, n_rsr
    rsr_i = rsr(i_rsr)
    if (rsr_i%disabled) cycle
    max_resno = i_R(1,rsr_i%atmno(1))
    min_resno = max_resno
    do i_atm = 2, rsr_i%natm
        resno = i_R(1,rsr_i%atmno(i_atm))
        max_resno = max(max_resno, resno)
        min_resno = min(min_resno, resno)
    end do
    resdiff = max_resno - min_resno
 
    if (resdiff < span_range) then
        rsr(i_rsr)%selected = .true.
    end if
end do

end subroutine select_rsr
!-------------------------------------------------------------------------------
subroutine selected_res_rsr(Escale, selected_res)
!-------------------------------------------------------------------------------
! Restraint selection by Escale and selected residues
!-------------------------------------------------------------------------------
type(energy_type), intent(in) :: Escale
logical, intent(in) :: selected_res(tn%residue)
type(restraint_type) :: rsr_i
logical :: allowed(5)
integer :: i_rsr, i_res, i_atm

if (.not. use_rsr) return

rsr(:)%selected = .true.
allowed(:) = .false.
if (Escale%modeling(16) > small_real) allowed(1) = .true.
if (Escale%modeling(17) > small_real) allowed(2) = .true.
if (Escale%modeling(18) > small_real) allowed(3) = .true.
if (Escale%modeling(19) > small_real) allowed(4) = .true.
if (Escale%modeling(20) > small_real) allowed(5) = .true.

do i_rsr = 1, n_rsr
    rsr_i = rsr(i_rsr)
    if (rsr_i%disabled .or. (.not. allowed(rsr_i%grp))) then
        rsr(i_rsr)%selected = .false.
    else
        do i_atm = 1, rsr_i%natm
            i_res = i_R(1,rsr_i%atmno(i_atm))
            if (.not. selected_res(i_res)) then
                rsr(i_rsr)%selected = .false.
                exit
            end if
        end do
    end if
end do

end subroutine selected_res_rsr
!-------------------------------------------------------------------------------
! Below: subroutines related to rsr update
!-------------------------------------------------------------------------------
subroutine protein_to_topolid(molecule, update_lig_rsr, R_updt, use_R_updt)
!-------------------------------------------------------------------------------
! From the given molecule type this fills in atom idx and coordinate information
! for rsr update. Also reads heteroatoms of a ligand (with mol2 information).
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: molecule
logical, intent(in) :: update_lig_rsr
real(dp), intent(out) :: R_updt(3, tn%atom)
logical, intent(out) :: use_R_updt(tn%atom)
!
integer :: i_res, i_ref, i_atm, i_lig, ia
character(len=4) :: atomName

call protein_to_R(molecule)
use_R_updt(:) = .false.

do i_res = 1, tn%stdres
    i_ref = molecule%residue(i_res)%res_type
    do i_atm = 1, ref_res(i_ref)%n_atm
        ia = ii_R(i_atm, i_res)
        R_updt(:,ia) = R(:,ia)
        use_R_updt(ia) = .true.
    end do
end do

if (update_lig_rsr .and. molecule%n_lig /= 0) then
    do i_lig = 1, tn%ligand
        i_ref = molecule%ligand(i_lig)%lig_type
        do i_atm = 1, ref_lig(i_ref)%n_atm
            atomName = ref_lig(i_ref)%atom_name(i_atm)
            if (atomName(1:1) == 'H') cycle
            ia = ii_L(i_atm, i_lig)
            R_updt(:,ia) = R(:,ia)
            use_R_updt(ia) = .true.
        end do
    end do
end if

end subroutine protein_to_topolid
!-------------------------------------------------------------------------------
subroutine read_update_rsr_pdb(molecule, update_rsr_pdb, update_lig_rsr, &
    R_updt, use_R_updt)
!-------------------------------------------------------------------------------
! Fill in rsr update information from 'update_rsr_pdb'
! TODO comment
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: molecule
character(len=len_fname), intent(in) :: update_rsr_pdb
logical, intent(in) :: update_lig_rsr
real(dp), intent(out) :: R_updt(3, tn%atom)
logical, intent(out) :: use_R_updt(tn%atom)
!
character(len=len_fname) :: line
character(len=6) :: atom_error_mode
character(len=4) :: atomName
character(len=3) :: hetName
integer :: f_unit = 25, openstat, ioerror
integer :: i_ref, i_res, i_atm, resNo, resNo_prev
integer :: lig_no, i_mol2
integer :: ia

open(f_unit, file = trim(update_rsr_pdb), iostat = openstat)
if (openstat > 0) then
    call terminate_with_error('Terminate with error: No update_rsr_pdb file found. Check {update_rsr_pdblist}.')
end if

write(log_msg,"(A,A)") "Updating the restraints using ", trim(update_rsr_pdb)
call log_p(log_msg, level=10)

i_res = 0
resNo_prev = -9999
use_R_updt(:) = .false.
lig_no = 0
atom_error_mode = 'ignore'

! Reading update_rsr_pdb file
do
    read(f_unit, "(A120)", iostat=ioerror) line
    if (ioerror < 0) exit
    if (line(1:4) /= 'ATOM' .and. line(1:6) /= 'HETATM') cycle
    if (line(13:13) == ' ') then
        read(line(14:17),"(A4)") atomName
    else
        read(line(13:16),"(A4)") atomname
        atomName = trim(atomName(2:4))//atomName(1:1)
    end if
    !
    if (line(1:4) == 'ATOM') then
        read(line(23:26),"(I4)") resNo
        if (resNo /= resNo_prev) then
            resNo_prev = resNo
            i_res = i_res + 1
        end if
        i_ref = molecule%residue(i_res)%res_type
        call find_atom_idx(i_res, i_ref, atomName, i_atm, atom_error_mode)
        if (i_atm == -100) cycle
        ia = ii_R(i_atm, i_res)
        read(line(31:54),"(3f8.3)") R_updt(:,ia)
        use_R_updt(ia) = .true.

    else if (update_lig_rsr .and. (line(1:6) == 'HETATM')) then
        hetName = line(18:20)
        do i_mol2 = 1, n_mol2_top
            if (hetName == trim(infile_mol2_topo(2,i_mol2))) then
                read(line(23:26),"(I4)") resNo
                if (resNo /= resNo_prev) then
                    resNo_prev = resNo
                    lig_no = lig_no + 1
                end if
            end if
        end do
        if (atomName(1:1) == 'H') cycle
        i_ref = molecule%ligand(lig_no)%lig_type
        call find_ligand_atom_idx(lig_no, i_ref, atomName, i_atm, atom_error_mode)
        ia = ii_L(i_atm, lig_no)
        read(line(31:54),"(3f8.3)") R_updt(:,ia)
        use_R_updt(ia) = .true.
    end if
end do

close (f_unit)

end subroutine read_update_rsr_pdb
!-------------------------------------------------------------------------------
subroutine update_pair_rsr_info(n_res, R_updt, use_R_updt, update_NO_rsr, &
    n_add, rij_s, use_rij)
!-------------------------------------------------------------------------------
! TODO comment
!-------------------------------------------------------------------------------
integer, intent(in) :: n_res
real(dp), intent(in) :: R_updt(3,tn%atom)
logical, intent(in) :: use_R_updt(tn%atom)
logical, intent(in) :: update_NO_rsr
integer, intent(out) :: n_add(2)
real(dp), intent(out) :: rij_s(3,tn%stdres,tn%stdres)
logical, intent(out) :: use_rij(3,tn%stdres,tn%stdres)
integer :: i_res, j_res, i_atm, j_atm
real(dp) :: dr(3), rij

n_add(1) = 0
n_add(2) = 0
use_rij(:,:,:) = .false.

do i_res = 1, n_res-2
    do j_res = i_res+2, n_res
        ! CA-CA
        i_atm = res_index(i_res)%Ca_id(1)
        j_atm = res_index(j_res)%Ca_id(1)
        if (use_R_updt(i_atm) .and. use_R_updt(j_atm)) then
            dr(:) = R_updt(:,i_atm) - R_updt(:,j_atm)
            rij = sqrt(dot_product(dr(:),dr(:)))
            if (rij < update_rsr_dcut) then
                n_add(1) = n_add(1) + 1
                rij_s(1,i_res,j_res) = rij
                use_rij(1,i_res,j_res) = .true.
            end if
        end if
    end do
end do

if (update_NO_rsr) then
    do i_res = 1, n_res-2
        do j_res = i_res+2, n_res
            !N-O
            i_atm = res_index(i_res)%bb_id(1,1)
            j_atm = res_index(j_res)%bb_id(4,1)
            if (use_R_updt(i_atm) .and. use_R_updt(j_atm)) then
                dr(:) = R_updt(:,i_atm) - R_updt(:,j_atm)
                rij = sqrt(dot_product(dr(:),dr(:)))
                if (rij < update_rsr_dcut) then
                    n_add(2) = n_add(2) + 1
                    rij_s(2,i_res,j_res) = rij
                    use_rij(2,i_res,j_res) = .true.
                end if
            end if
            !O-N
            i_atm = res_index(i_res)%bb_id(4,1)
            j_atm = res_index(j_res)%bb_id(1,1)
            if (use_R_updt(i_atm) .and. use_R_updt(j_atm)) then
                dr(:) = R_updt(:,i_atm) - R_updt(:,j_atm)
                rij = sqrt(dot_product(dr(:),dr(:)))
                if (rij < update_rsr_dcut) then
                    n_add(2) = n_add(2) + 1
                    rij_s(3,i_res,j_res) = rij
                    use_rij(3,i_res,j_res) = .true.
                end if
            end if
        end do
    end do
end if

end subroutine update_pair_rsr_info
!-------------------------------------------------------------------------------
subroutine update_rsr(molecule, R_updt, use_R_updt, update_mode, &
    update_NO_rsr, update_lig_rsr, append_mode)
!-------------------------------------------------------------------------------
! TODO comment
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: molecule
integer, intent(in) :: update_mode   ! 0: distance only, 1: cartesian only, 2: both
real(dp), intent(in) :: R_updt(3,tn%atom)
logical, intent(in) :: use_R_updt(tn%atom)
logical, intent(in) :: update_NO_rsr
logical, intent(in) :: update_lig_rsr
logical, optional, intent(in) :: append_mode
!
integer :: n_res
integer :: i_res, j_res, i_lig, i_ref, i_atm, k
integer :: ia, i_rsr, n_new, n_add(4), n_del, n_ligatm
logical :: use_rij(3,tn%stdres,tn%stdres)
real(dp) :: rij_s(3,tn%stdres,tn%stdres)

n_res = molecule%n_res

n_add(:) = 0
if (update_mode == 0 .or. update_mode == 2) then
    use_rij(:,:,:) = .false.
    call update_pair_rsr_info(n_res, R_updt, use_R_updt, update_NO_rsr, &
        n_add(1:2), rij_s, use_rij)
end if
if (update_mode == 1 .or. update_mode == 2) then
    !number of updated CA atoms
    n_add(3) = n_res
end if
if (update_lig_rsr .and. molecule%n_lig /= 0 .and. &
    (update_mode == 1 .or. update_mode == 2)) then
    !number of updated ligand atoms
    n_ligatm = 0
    do i_lig = 1, molecule%n_lig
        i_ref = molecule%ligand(i_lig)%lig_type
        do i_atm = 1, ref_lig(i_ref)%n_atm
            ia = ii_L(i_atm, i_lig)
            if (.not. use_R_updt(ia)) cycle
            n_ligatm = n_ligatm + 1
        end do
    end do
    n_add(4) = n_ligatm
end if

n_new = sum(n_add)

n_rsr = n_rsr0
n_rsrreg = n_rsrreg0
deallocate(rsr)
deallocate(reg_parm)

allocate(rsr(n_rsr+n_new))
allocate(reg_parm(n_rsrreg+n_new))
rsr(1:n_rsr) = rsr0(1:n_rsr)
reg_parm(1:n_rsrreg) = reg_parm0(1:n_rsrreg)

n_del = 0
if (present(append_mode) .and. append_mode) then
    continue
else
    do i_rsr = 1, n_rsr
        if (rsr(i_rsr)%grp == 1) then
            n_del = n_del+1
            rsr(i_rsr)%selected = .false.
            rsr(i_rsr)%disabled = .true.
        else if (update_NO_rsr .and. rsr(i_rsr)%grp == 2) then
            n_del = n_del+1
            rsr(i_rsr)%selected = .false.
            rsr(i_rsr)%disabled = .true.
        end if
    end do
end if

write(log_msg,"(A,I6)") '  - Restraints are updated:', n_new
call log_p(log_msg, level=10)
write(log_msg,"(A,I6)") '  - Restraints are deactivated:', n_del
call log_p(log_msg, level=10)
write(log_msg,"(A,I6)") '  - Total number of active restraints:', n_rsr-n_del+n_new
call log_p(log_msg, level=10)

if (update_mode == 0 .or. update_mode == 2) then
    ! Distance restraint update. k=1: CA-CA 2: N-O 3:O-N
    i_rsr = n_rsr
    do i_res = 1, n_res-2
        do j_res = i_res+2, n_res
            do k = 1, 3
                if (use_rij(k,i_res,j_res)) then
                    i_rsr = i_rsr + 1
                    rsr(i_rsr)%form = 3
                    rsr(i_rsr)%modality = 1
                    rsr(i_rsr)%feature = 1
                    if (k == 1) then
                        rsr(i_rsr)%grp = 1
                    else
                        rsr(i_rsr)%grp = 2
                    end if
                    rsr(i_rsr)%natm = 2
                    allocate(rsr(i_rsr)%g(3,rsr(i_rsr)%natm))
                    rsr(i_rsr)%nparm = 2
                    if (k == 1) then
                        rsr(i_rsr)%atmno(1) = res_index(i_res)%Ca_id(1)
                        rsr(i_rsr)%atmno(2) = res_index(j_res)%Ca_id(1)
                    else if (k == 2) then
                        rsr(i_rsr)%atmno(1) = res_index(i_res)%bb_id(1,1)
                        rsr(i_rsr)%atmno(2) = res_index(j_res)%bb_id(4,1)
                    else if (k == 3) then
                        rsr(i_rsr)%atmno(1) = res_index(i_res)%bb_id(4,1)
                        rsr(i_rsr)%atmno(2) = res_index(j_res)%bb_id(1,1)
                    end if
                    rsr(i_rsr)%selected = .true.
                    rsr(i_rsr)%disabled = .false.
                 
                    n_rsrreg = n_rsrreg + 1
                    rsr(i_rsr)%parm_id(1) = 1
                    rsr(i_rsr)%parm_id(2) = n_rsrreg

                    reg_parm(n_rsrreg)%parm(1) = rij_s(k,i_res,j_res)
                    reg_parm(n_rsrreg)%parm(2) = update_rsr_sig
                end if
            end do
        end do
    end do
end if

if (update_mode == 1 .or. update_mode == 2) then
    ! Cartesian restraint for CA
    do i_res = 1, n_res
        ia = res_index(i_res)%Ca_id(1)
        if (.not. use_R_updt(ia)) cycle
        i_rsr = i_rsr + 1
        rsr(i_rsr)%form = 3
        rsr(i_rsr)%modality = 1
        rsr(i_rsr)%feature = 9
        rsr(i_rsr)%grp = 1
        rsr(i_rsr)%natm = 1
        allocate(rsr(i_rsr)%g(3,rsr(i_rsr)%natm))
        rsr(i_rsr)%nparm = 5
        rsr(i_rsr)%atmno(1) = ia
        rsr(i_rsr)%selected = .true.
        rsr(i_rsr)%disabled = .false.
        
        n_rsrreg = n_rsrreg + 1
        rsr(i_rsr)%parm_id(1) = 1
        rsr(i_rsr)%parm_id(2) = n_rsrreg

        reg_parm(n_rsrreg)%parm(1) = 0.0d0
        reg_parm(n_rsrreg)%parm(2) = update_rsr_sig
        reg_parm(n_rsrreg)%parm(3:5) = R_updt(1:3, ia)
    end do

    if (update_lig_rsr) then
        !Cartesian restraint for ligand atoms
        do i_lig = 1, molecule%n_lig
            i_ref = molecule%ligand(i_lig)%lig_type
            do i_atm = 1, ref_lig(i_ref)%n_atm
                ia = ii_L(i_atm, i_lig)
                if (.not. use_R_updt(ia)) cycle
                i_rsr = i_rsr + 1
                rsr(i_rsr)%form = 3
                rsr(i_rsr)%modality = 1
                rsr(i_rsr)%feature = 9
                rsr(i_rsr)%grp = 5
                rsr(i_rsr)%natm = 1
                allocate(rsr(i_rsr)%g(3,rsr(i_rsr)%natm))
                rsr(i_rsr)%nparm = 5
                rsr(i_rsr)%atmno(1) = ia
                rsr(i_rsr)%selected = .true.
                rsr(i_rsr)%disabled = .false.
                
                n_rsrreg = n_rsrreg + 1
                rsr(i_rsr)%parm_id(1) = 1
                rsr(i_rsr)%parm_id(2) = n_rsrreg
                
                reg_parm(n_rsrreg)%parm(1) = 0.0d0
                reg_parm(n_rsrreg)%parm(2) = update_lig_rsr_sig
                reg_parm(n_rsrreg)%parm(3:5) = R_updt(1:3, ia)
            end do
        end do
    end if
end if

n_rsr = n_rsr+n_new

end subroutine update_rsr
!-------------------------------------------------------------------------------
subroutine clear_updated_rsr()
!-------------------------------------------------------------------------------
n_rsr = n_rsr0
n_rsrreg = n_rsrreg0
deallocate(rsr)
deallocate(reg_parm)

allocate(rsr(n_rsr))
allocate(reg_parm(n_rsrreg))
rsr(1:n_rsr) = rsr0(1:n_rsr)
reg_parm(1:n_rsrreg) = reg_parm0(1:n_rsrreg)

end subroutine clear_updated_rsr
!-------------------------------------------------------------------------------
subroutine construct_ligdock_rsr_grid(ligand, grid_info, dock_grid, flex_res)
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand
type(docking_grid_param_type), intent(in) :: grid_info
type(docking_grid_type), intent(inout) :: dock_grid
logical, intent(in) :: flex_res(:)
integer :: n_hash
real(dp), allocatable :: tmp_grid_E(:)
real(dp) :: ref_pt(3), grid_pt(3)
integer :: i_lig, i, x, y, z

n_hash = grid_info%n_elem(1) * grid_info%n_elem(2) * grid_info%n_elem(3) 
! allocate autodock grid
allocate(dock_grid%rsr_grid(4, n_hash, ligand%n_atm))
allocate(tmp_grid_E(ligand%n_atm))

tmp_grid_E(:) = 0.0d0

do i = 1, 3
    ref_pt(i) = grid_info%grid_cntr(i) &
              - (dble(grid_info%n_elem(i))-1.0d0)/2.0d0 * grid_info%grid_width
end do

do z = 1, grid_info%n_elem(3)
    do y = 1, grid_info%n_elem(2)
        do x = 1, grid_info%n_elem(1)
            grid_pt(:) = (dble((/x, y, z/))-1.0d0)*grid_info%grid_width + ref_pt(:)
            call calc_rsr_grid(ligand, grid_pt(1:3), tmp_grid_E, flex_res)
            do i_lig = 1, ligand%n_atm
                call copy_to_grid(tmp_grid_E(i_lig), grid_info%n_elem(1:3),&
                                  dock_grid%rsr_grid(:,:,i_lig), x, y, z, n_hash)
            end do
        end do
    end do
end do

deallocate(tmp_grid_E)

end subroutine construct_ligdock_rsr_grid
!-------------------------------------------------------------------------------
subroutine calc_rsr_grid(ligand, grid_pt, E_pt, flex_res)
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand
real(dp), intent(in) :: grid_pt(3)
real(dp), intent(inout) :: E_pt(:)
logical, intent(in) :: flex_res(:)
integer :: i_atm, i_rsr, i_prot, i_lig
integer :: n_g
real(dp) :: dr(3), x_i
real(dp) :: parm(max_longparm), argv(2)
logical :: is_valid

E_pt(:) = 0.0d0

do i_atm = 1, ligand%n_atm
    i_lig = ii_L(i_atm, ligand%lig_no)
    do i_rsr = 1, n_rsr
        call check_valid_rsr(i_lig, rsr(i_rsr), i_prot, is_valid)
        if (.not. is_valid) cycle
        
        dr(:) = R(:,i_prot) - grid_pt(:)
        x_i = sqrt(dot_product(dr,dr))
        ! Check restraint parameter type (long or regular)
        if (rsr(i_rsr)%parm_id(1) == 1) then
            parm(1:rsr(i_rsr)%nparm) &
                = reg_parm(rsr(i_rsr)%parm_id(2))%parm(1:rsr(i_rsr)%nparm)
        else
            parm(1:rsr(i_rsr)%nparm) &
                = long_parm(rsr(i_rsr)%parm_id(2))%parm(1:rsr(i_rsr)%nparm)
        end if
    
        n_g = rsr(i_rsr)%modality
        argv(:) = RT*multiple_gaussian(n_g, x_i, parm(n_g+1:2*n_g), &
            parm(1:n_g), parm(2*n_g+1:3*n_g), rsr(i_rsr)%is_ang, .true., .false.)
        E_pt(i_atm) = E_pt(i_atm) + argv(1)
    end do
end do

do i_atm = 1, ligand%n_atm
    if (E_pt(i_atm) == 0.0) cycle
    E_pt(i_atm) = min(E_pt(i_atm), 1000.0)
end do

end subroutine calc_rsr_grid
!-------------------------------------------------------------------------------
subroutine check_valid_rsr(atm_idx, rsr, partner, is_valid)
!-------------------------------------------------------------------------------
integer, intent(in) :: atm_idx
type(restraint_type), intent(in) :: rsr
integer, intent(out) :: partner
logical, intent(out) :: is_valid

is_valid = .false.
if (atm_idx == rsr%atmno(1)) then
    is_valid = .true.
    partner = rsr%atmno(2)
else if (atm_idx == rsr%atmno(2)) then
    is_valid = .true.
    partner = rsr%atmno(1)
end if

end subroutine check_valid_rsr
!-------------------------------------------------------------------------------
subroutine calc_ligdock_rsr_using_grid(ligand, grid_info, rsr_grid, rsr_E, &
                                       g, calc_g)
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand
type(docking_grid_param_type), intent(in) :: grid_info
real(dp), intent(in) :: rsr_grid(:,:,:)
real(dp), intent(out) :: rsr_E, g(:,:)
logical, intent(in) :: calc_g
integer :: i_atm, atm_idx
real(dp) :: tmp_E

rsr_E = 0.0d0
do i_atm = 1, ligand%n_atm
    atm_idx = ii_L(i_atm, ligand%lig_no)
    call calc_E_using_grid(R(1:3, atm_idx), grid_info, rsr_grid(:,:,i_atm), &
                           tmp_E, g(:,i_atm), calc_g)
    rsr_E = rsr_E + tmp_E
end do

end subroutine calc_ligdock_rsr_using_grid
!-------------------------------------------------------------------------------
subroutine update_restraints_as_meld(n_features, n_grp, active_fraction)
!-------------------------------------------------------------------------------
! For restraint energies of current R, sort by rsr energy and select given number
! of lowest restraints 
! Reference: Alberto Perez, Justin L. MacCallum, and Ken A. Dill, Accelerating 
! molecular simulations of proteins using Bayesian inference on weak information,
! PNAS, 2015. 
!-------------------------------------------------------------------------------
integer, intent(in) :: n_features
integer, intent(in) :: n_grp
real(dp), intent(in) :: active_fraction(:,:)
type(restraint_type) :: rsr_i
real(dp) :: values(n_rsr) 
integer :: key_rsr(n_rsr), key(n_rsr)
integer :: i, j, i_rsr, n_curr_rsr
integer :: n_sel, i_sel_rsr

key(:) = 0
do i = 1, n_features
    do j = 1, n_grp
        if (active_fraction(j,i) == 0.0d0) cycle
        values(:) = max_energy
        key_rsr(:) = 0
        n_curr_rsr = 0
        do i_rsr = 1, n_rsr
            rsr_i = rsr(i_rsr)
            if (rsr_i%feature == i .and. rsr_i%grp == j) then
                n_curr_rsr = n_curr_rsr + 1
                values(n_curr_rsr) =  rsr_i%E
                key_rsr(n_curr_rsr) = i_rsr
            end if
        end do
        !
        if (n_curr_rsr == 0) cycle
        call sort2(n_curr_rsr, values(1:n_curr_rsr), key(1:n_curr_rsr))
        n_sel = ceiling(n_curr_rsr*active_fraction(j,i))
        do i_rsr = 1, n_curr_rsr
            i_sel_rsr = key_rsr(key(i_rsr))
            if (i_rsr <= n_sel) then
                rsr(i_sel_rsr)%active = .true.
            else
                rsr(i_sel_rsr)%active = .false.
            end if
        end do
    end do
end do

end subroutine update_restraints_as_meld
!-------------------------------------------------------------------------------
END MODULE RESTRAINT
