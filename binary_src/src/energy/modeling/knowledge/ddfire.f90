!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/modeling/knowledge/ddfire.f90
!
! Description: dDFIRE statistical potential
!  This module is for evaluaing dDFIRE score, component of GALAXY global energy.
!  First initialize module by calling [setup_ddfire].
!  The original dDFIRE score having discrete values on certain bin can be 
!   calculated by calling subroutine [ddfire_score].
!  To use dDFIRE as an analytic function, call subroutine [ddfire_static_score]. 
!  This is cubic spline interpolated version for original dDFIRE score.
!
! Reference : Yuedong Yang et al. Specific interactions for ab initio 
!             folding of protein terminal regions with secondary structures, 
!             Proteins 2008, 72, 793-803.
!
! TODO: Have to re-write this module with GOAP.f90 style
!-------------------------------------------------------------------------------
MODULE dDFIRE
!-------------------------------------------------------------------------------

use globals
use logger
!
use in_out_utils,     only: find_atom_idx
use string,           only: parse_longstring
use mathfunctions,    only: cubic_spline
use convert_res_name, only: convert_to_stdres
!
use energy_vars,      only: R, i_R, ii_R, i_P, LRoff, max_neigh2, &
                            ddfire_add_scale, dfire_score_file, &
                            dfire_atype_file, max_energy

implicit none
save
private

integer :: trunc_bin             ! no. of distance bins used for enumeration
integer :: nn, np                ! no. of non-dipole / dipole atom types

type score_pair
    real :: val1
    real :: val2
end type score_pair

type(score_pair), allocatable :: nnscore_pair(:,:,:), npscore_pair(:,:,:,:), ppscore_pair(:,:,:,:)
! Specific pairlist of DFIRE
integer, allocatable :: dfire_atmid(:)    ! mapping atom no. into dDFIRE atom type
logical, allocatable :: is_using_atm(:)   ! whether atom is of interest in dDFIRE

! For dipole enumeration
real(dp), allocatable :: polarv(:,:)      ! current polar vectors
integer, allocatable :: n_connected(:)    ! no. connected to atom of interest
integer, allocatable :: con_atmno(:,:)    ! connected atom index
integer, allocatable :: n_conid(:)        ! no. connected atoms (to form dipole)
logical :: ispolar(192)                   ! whether atom is polar atom type
integer :: atom2polar(192)                ! index for polar atoms
character(len=4), allocatable :: connected_list(:,:)  ! atom name for connected

integer :: sequence_separation

public :: initialize_ddfire
public :: finalize_ddfire
!
public :: calc_ddfire_score
!
public :: pairwise_ddfire_static_score
public :: ddfire_score

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_ddfire(protein)
!-------------------------------------------------------------------------------
! Initial setup of dDFIRE for modeling
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
integer :: i_atm, ref_res_no, ref_res_prv, i, atm_in_res(2)
integer :: resno, atmno, id
character(len=6) :: atomorder(192), polarorder(66), error_mode
character(len=4) :: atm_name, res_name
character(len=1) :: ter_type

error_mode = 'ignore'

! Set truncation site where bin > trunc_bin is truncated
trunc_bin = int(2.0d0*LRoff)

allocate(n_connected(tn%atom), con_atmno(2,tn%atom))
allocate(dfire_atmid(tn%atom), polarv(3,tn%atom))
allocate(is_using_atm(tn%atom))

! Read tables and then set arrays
call read_atomlist(dfire_atype_file, atomorder, polarorder)
allocate(nnscore_pair(30,nn,nn), npscore_pair(30,6,nn,np), ppscore_pair(30,6,np,np))
! nn = 192, np = 69 

allocate(n_conid(np), connected_list(2,np))

call read_ddfire_table(dfire_score_file)
call setup_ddfire_misc(protein, atomorder)

! Setup for dipole
if (trim(top_type) /= 'coarse') then
    if (trim(top_type) /= 'cbeta') then
        ! Define connected atoms for dipole enumerations
        call assign_connected()
    end if
    
    n_connected(:) = 0
    con_atmno(:,:) = 0

    ! Match actual atomic numbers connected to each other
    do i_atm = 1, tn%stdatm
        if (dfire_atmid(i_atm) == 0) cycle
        id = atom2polar(dfire_atmid(i_atm))
        if (id == 0) cycle
        
        resno = i_R(1,i_atm)
        atmno = i_R(2,i_atm)
        if (resno > 1) ref_res_prv = ref_res_no
        ref_res_no = protein%residue(resno)%res_type
        res_name = ref_res(ref_res_no)%res_name
        atm_name = ref_res(ref_res_no)%atom_name(atmno)
        ter_type = protein%residue(resno)%ter_type
        
        call find_ddfire_connected(resno, n_connected(i_atm), atm_in_res, &
                                   ref_res_no, id, atm_name, res_name)
        
        do i = 1, n_connected(i_atm)
            con_atmno(i,i_atm) = ii_R(atm_in_res(i),resno)
        end do
        ! Take special care for -C
        if (atm_name == 'N   ' .and. ter_type /= 'N') then
            n_connected(i_atm) = 2
            call find_atom_idx(resno, ref_res_prv, 'C   ', atm_in_res(2), error_mode)
            con_atmno(2,i_atm) = ii_R(atm_in_res(2),resno-1)
        end if
    end do
end if

call setup_static_score()

if (top_type == 'cbeta') then
    sequence_separation = 7
else
    sequence_separation = 0
end if

end subroutine initialize_ddfire
!-------------------------------------------------------------------------------
subroutine finalize_ddfire()
!-------------------------------------------------------------------------------
deallocate(n_connected, dfire_atmid, con_atmno, polarv)
deallocate(nnscore_pair, npscore_pair, ppscore_pair)
end subroutine finalize_ddfire
!-------------------------------------------------------------------------------
subroutine setup_static_score()
!-------------------------------------------------------------------------------
! Initial call for dDFIRE score and precalculation of the interpolation by cubic spline
!-------------------------------------------------------------------------------
real(dp) :: c(trunc_bin+1), c_pp(0:trunc_bin+1), u(trunc_bin+1), x(trunc_bin+1), p, un
integer :: i_type, j_type, angid, i

nnscore_pair(:,:,:)%val2 = 0.0d0
if (np > 0) then
    npscore_pair(:,:,:,:)%val2 = 0.0d0
    ppscore_pair(:,:,:,:)%val2 = 0.0d0
end if

! Definition of binning by PHB
do i = 1, trunc_bin + 1
    x(i) = 0.5d0*i - 0.25d0 ! Table value -> center of the bin
end do

! Distance term
do i_type = 1, nn
    do j_type = 1, nn
        trunc_bin = int(2.0d0*LRoff)
        if ((i_type >= 168 .and. i_type <= 170) .or. &
            (i_type >= 176 .and. i_type <= 178) .or. &
            (i_type >= 185 .and. i_type <= 187)) & 
            trunc_bin = 19 ! For Phos residue(10.0 cut)
        ! Fill spline arrays recursively
        c(trunc_bin+1) = 0.0d0
        c_pp(0) = 0.0d0
        c_pp(1) = -0.5d0
        c(1:trunc_bin) = nnscore_pair(1:trunc_bin,j_type,i_type)%val1
        u(1) = 12.0d0*(c(2) - c(1))

        do i = 2, trunc_bin
            p = 0.5d0*c_pp(i-1) + 2.0d0
            c_pp(i) = -0.5d0/p
            u(i) = (12.0d0*(c(i+1) + c(i-1) - 2.0d0*c(i)) - 0.5d0*u(i-1)) / p
        end do
        un = -12.0d0*(c(trunc_bin+1) - c(trunc_bin))
        c_pp(trunc_bin+1) = (un - 0.5d0*u(trunc_bin)) / (0.5d0*c_pp(trunc_bin) + 1.0d0)
        
        do i = trunc_bin, 1, -1
            c_pp(i) = c_pp(i)*c_pp(i+1) + u(i)
        end do

        ! Save second derivative
        nnscore_pair(1:trunc_bin+1,j_type,i_type)%val2 = c_pp(1:trunc_bin+1)
    end do
end do
trunc_bin = int(2.0d0*LRoff)

! Dipole - nondipole interaction
do i_type = 1, np
    do j_type = 1, nn
        do angid = 1, 6

            ! Fill spline arrays recursively
            c(trunc_bin+1) = 0.0d0
            c_pp(0) = 0.0d0
            c_pp(1) = -0.5d0
            c(1:trunc_bin) = npscore_pair(1:trunc_bin,angid,j_type,i_type)%val1
            u(1) = 12.0d0*(c(2) - c(1))
           
            do i = 2,trunc_bin
                p = 0.5d0*c_pp(i-1) + 2.0d0
                c_pp(i) = -0.5d0/p
                u(i) = (12.0d0*(c(i+1) + c(i-1) - 2.0d0*c(i)) - 0.5d0*u(i-1)) / p
            end do
            un = -12.0d0*(c(trunc_bin+1) - c(trunc_bin))
            c_pp(trunc_bin+1) = (un - 0.5d0*u(trunc_bin)) / &
                (0.5d0*c_pp(trunc_bin) + 1.0d0)
           
            do i = trunc_bin, 1, -1
                c_pp(i) = c_pp(i)*c_pp(i+1) + u(i)
            end do

            ! Save second derivative
            npscore_pair(1:trunc_bin+1,angid,j_type,i_type)%val2 = c_pp(1:trunc_bin+1)

        end do
    end do
end do

! Dipole-dipole interaction
do i_type = 1, np
    do j_type = 1, np
        do angid = 1, 6
            ! Fill spline arrays recursively
            c(trunc_bin+1) = 0.0d0
            c_pp(0) = 0.0d0
            c_pp(1) = -0.5d0
            c(1:trunc_bin) = ppscore_pair(1:trunc_bin,angid,j_type,i_type)%val1
            u(1) = 12.0d0*(c(2)-c(1))
           
            do i = 2, trunc_bin
                p = 0.5d0*c_pp(i-1) + 2.0d0
                c_pp(i) = -0.5d0/p
                u(i) = (12.0d0*(c(i+1) + c(i-1) - 2.0d0*c(i)) - 0.5d0*u(i-1)) / p
            end do
            un = -12.0d0*(c(trunc_bin+1) - c(trunc_bin))
            c_pp(trunc_bin+1) = (un - 0.5d0*u(trunc_bin)) / &
                (0.5d0*c_pp(trunc_bin) + 1.0d0)
           
            do i = trunc_bin, 1, -1
                c_pp(i) = c_pp(i)*c_pp(i+1) + u(i)
            end do

            ! Save second derivative
            ppscore_pair(1:trunc_bin+1,angid,j_type,i_type)%val2 = c_pp(1:trunc_bin+1)

        end do
    end do
end do

end subroutine setup_static_score
!-------------------------------------------------------------------------------
subroutine calc_ddfire_score(ff, g_dyn4, appl_respair, calc_g)
!-------------------------------------------------------------------------------
! To calculate ddfire energy STATIC (for every function call)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: ff(2), g_dyn4(3,tn%atom,2)
logical, intent(in) :: appl_respair(tn%residue,tn%residue), calc_g
integer :: i_pair, res1, res2
integer :: rid, no1, no2, pno1, pno2, angid, atm1, atm2
integer :: atm2_id(max_neigh2), n_pair
real(dp) :: r_pq(3), pvp(3), pvq(3)
real(dp) :: dist, costheta
real(dp) :: fdiff, a, c_j(2), c_jpp(2), arg(2), dcdfv(3)
logical :: status

ff(:) = 0.0d0
g_dyn4(:,:,:) = 0.0d0

! Calculate polar vector only if dipole term is turned on.
if (ddfire_add_scale > small_real) call build_polarv(R)

! Iter over for pre-defined pairlist (at build_dfire_pair)
do atm1 = 1, tn%atom
    if (.not. is_using_atm(atm1)) cycle
    res1 = i_R(1,atm1)

    call get_dfire_pair(atm1, appl_respair, n_pair, atm2_id, status)
    if (.not. status) then
        ff(:) = max_energy
        return
    end if
    do i_pair = 1, n_pair
        atm2 = atm2_id(i_pair)
        res2 = i_R(1,atm2)

        no1 = dfire_atmid(atm1)
        no2 = dfire_atmid(atm2)

        r_pq(:) = R(:,atm2) - R(:,atm1)
        dist = sqrt(dot_product(r_pq,r_pq))

        r_pq(:) = r_pq(:)/dist
        rid = int(2.0d0*dist - 0.5d0) + 1
        ! Difference of the distance value from former bin point
        fdiff = 0.5d0*rid + 0.25d0 - dist

        a = 2.0d0*fdiff
        c_j(:)   = nnscore_pair(rid:rid+1,no2,no1)%val1
        c_jpp(:) = nnscore_pair(rid:rid+1,no2,no1)%val2

        arg(:) = cubic_spline(0.5d0, c_j, c_jpp, a, calc_g)
        ff(1) = ff(1) + arg(1)
        if (calc_g) then
            dcdfv(:) = arg(2)*r_pq(:)
            g_dyn4(:,atm1,1) = g_dyn4(:,atm1,1) - dcdfv(:)
            g_dyn4(:,atm2,1) = g_dyn4(:,atm2,1) + dcdfv(:)
        end if

        ! From here dipolar terms: Cycle if dipole terms are turned off
        if (ddfire_add_scale < small_real) cycle
        if (abs(res1-res2) < sequence_separation) cycle
        
        if (ispolar(no1)) then
            pno1 = atom2polar(no1)
            pvp(:) = polarv(:,atm1)
        end if
        
        if (ispolar(no2)) then
            pno2 = atom2polar(no2)
            pvq(:) = polarv(:,atm2)
        end if

        !Orientation-dependent terms
        !1. Dipole-nondipole
        if (ispolar(no1) .and. (.not. ispolar(no2))) then
            costheta = dot_product(pvp(:), -r_pq(:))
            angid = min(6,abs(int(3*(costheta+1.0d0))) + 1)
            if (angid < 1) cycle
            c_j(:)   = npscore_pair(rid:rid+1,angid,no2,pno1)%val1
            c_jpp(:) = npscore_pair(rid:rid+1,angid,no2,pno1)%val2
            arg(:) = cubic_spline(0.5d0, c_j, c_jpp, a, calc_g)

            ff(2) = ff(2) + arg(1)
            dcdfv = arg(2)*r_pq(:)
            g_dyn4(:,atm1,2) = g_dyn4(:,atm1,2) - dcdfv(:)
            g_dyn4(:,atm2,2) = g_dyn4(:,atm2,2) + dcdfv(:)
           
        !2. Nondipole-dipole
        else if ((.not. ispolar(no1)) .and. ispolar(no2)) then
            costheta = dot_product(pvq(:), r_pq(:))
            angid = min(6,abs(int(3*(costheta+1.0d0))) + 1)
            if (angid < 1) cycle
            c_j(:)   = npscore_pair(rid:rid+1,angid,no1,pno2)%val1
            c_jpp(:) = npscore_pair(rid:rid+1,angid,no1,pno2)%val2
            arg(:) = cubic_spline(0.5d0, c_j, c_jpp, a, calc_g)
            ff(2) = ff(2) + arg(1)
            dcdfv(:) = arg(2)*r_pq(:)
            g_dyn4(:,atm1,2) = g_dyn4(:,atm1,2) - dcdfv(:)
            g_dyn4(:,atm2,2) = g_dyn4(:,atm2,2) + dcdfv(:)
           
        !3. Dipole-dipole
        else if (ispolar(no1) .and. ispolar(no2)) then
            ! as if dipole-nondipole
            costheta = dot_product(pvp(:), -r_pq(:))
            angid = min(6,abs(int(3*(costheta+1.0d0))) + 1)
            if (angid < 1) cycle
            c_j(:)   = npscore_pair(rid:rid+1,angid,no2,pno1)%val1
            c_jpp(:) = npscore_pair(rid:rid+1,angid,no2,pno1)%val2
            arg(:) = cubic_spline(0.5d0, c_j, c_jpp, a, calc_g)

            ff(2) = ff(2) + arg(1)
            dcdfv(:) = arg(2)*r_pq(:)
            g_dyn4(:,atm1,2) = g_dyn4(:,atm1,2) - dcdfv(:)
            g_dyn4(:,atm2,2) = g_dyn4(:,atm2,2) + dcdfv(:)

            ! as if nondipole-dipole
            costheta = dot_product(pvq(:), r_pq(:))
            angid = min(6,abs(int(3*(costheta+1.0d0))) + 1)
            if (angid < 1) cycle
            c_j(:)   = npscore_pair(rid:rid+1,angid,no1,pno2)%val1
            c_jpp(:) = npscore_pair(rid:rid+1,angid,no1,pno2)%val2
            arg(:) = cubic_spline(0.5d0, c_j, c_jpp, a, calc_g)

            ff(2) = ff(2) + arg(1)
            dcdfv(:) = arg(2)*r_pq(:)
            g_dyn4(:,atm1,2) = g_dyn4(:,atm1,2) - dcdfv(:)
            g_dyn4(:,atm2,2) = g_dyn4(:,atm2,2) + dcdfv(:)
           
            ! as if dipole-dipole
            costheta = dot_product(pvp(:), pvq(:))
            angid = min(6,abs(int(3*(costheta+1.0d0))) + 1)
            if (angid < 1) cycle
            c_j(:)   = ppscore_pair(rid:rid+1,angid,pno2,pno1)%val1
            c_jpp(:) = ppscore_pair(rid:rid+1,angid,pno2,pno1)%val2
            arg = cubic_spline(0.5d0, c_j, c_jpp, a, calc_g)
            ff(2) = ff(2) + arg(1)
            dcdfv(:) = arg(2)*r_pq(:)
            g_dyn4(:,atm1,2) = g_dyn4(:,atm1,2) - dcdfv(:)
            g_dyn4(:,atm2,2) = g_dyn4(:,atm2,2) + dcdfv(:)
        end if
    end do
end do

ff(2) = ff(2)
g_dyn4(:,:,2) = g_dyn4(:,:,2)

end subroutine calc_ddfire_score
!-------------------------------------------------------------------------------
subroutine pairwise_ddfire_static_score(ff, appl_respair)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: ff(tn%residue,tn%residue)
logical, intent(in) :: appl_respair(tn%residue,tn%residue)

integer :: i_pair, res1, res2
integer :: rid, no1, no2, pno1, pno2, angid, atm1, atm2
integer :: atm2_id(max_neigh2), n_pair
real(dp) :: r_pq(3), pvp(3), pvq(3)
real(dp) :: dist, costheta
real(dp) :: fdiff, a, c_j(2), c_jpp(2), arg(2)
logical, parameter :: calc_g = .false.
logical :: status

ff(:,:) = 0.0d0

! Calculate polar vector only if dipole term is turned on.
if (ddfire_add_scale > small_real) call build_polarv(R)

! Iter over for pre-defined pairlist (at build_dfire_pair)
do atm1 = 1, tn%atom
    if (.not. is_using_atm(atm1)) cycle
    res1 = i_R(1,atm1)

    call get_dfire_pair(atm1, appl_respair, n_pair, atm2_id, status)
    if (.not. status) then
        ff(:,:) = max_energy
        return
    end if
    do i_pair = 1, n_pair
        atm2 = atm2_id(i_pair)
        res2 = i_R(1,atm2)

        no1 = dfire_atmid(atm1)
        no2 = dfire_atmid(atm2)

        r_pq(:) = R(:,atm2) - R(:,atm1)
        dist = sqrt(dot_product(r_pq,r_pq))

        r_pq(:) = r_pq(:)/dist
        rid = int(2.0d0*dist - 0.5d0) + 1
        ! Difference of the distance value from former bin point
        fdiff = 0.5d0*rid + 0.25d0 - dist

        a = 2.0d0*fdiff
        c_j(:)   = nnscore_pair(rid:rid+1,no2,no1)%val1
        c_jpp(:) = nnscore_pair(rid:rid+1,no2,no1)%val2

        arg(:) = cubic_spline(0.5d0, c_j, c_jpp, a, calc_g)
        ff(res1,res2) = ff(res1,res2) + arg(1)
        ff(res2,res1) = ff(res2,res1) + arg(1)

        ! From here dipolar terms: Cycle if dipole terms are turned off
        if (ddfire_add_scale < small_real) cycle
        
        if (ispolar(no1)) then
            pno1 = atom2polar(no1)
            pvp(:) = polarv(:,atm1)
        end if
        
        if (ispolar(no2)) then
            pno2 = atom2polar(no2)
            pvq(:) = polarv(:,atm2)
        end if

        !Orientation-dependent terms
        !1. Dipole-nondipole
        if (ispolar(no1) .and. (.not. ispolar(no2))) then
            costheta = dot_product(pvp(:), -r_pq(:))
            angid = min(6,abs(int(3*(costheta+1.0d0))) + 1)

            c_j(:)   = npscore_pair(rid:rid+1,angid,no2,pno1)%val1
            c_jpp(:) = npscore_pair(rid:rid+1,angid,no2,pno1)%val2
            arg(:) = cubic_spline(0.5d0, c_j, c_jpp, a, calc_g)

            ff(res1,res2) = ff(res1,res2) + arg(1)
            ff(res2,res1) = ff(res2,res1) + arg(1)
           
        !2. Nondipole-dipole
        else if ((.not. ispolar(no1)) .and. ispolar(no2)) then
            costheta = dot_product(pvq(:), r_pq(:))
            angid = min(6,abs(int(3*(costheta+1.0d0))) + 1)

            c_j(:)   = npscore_pair(rid:rid+1,angid,no1,pno2)%val1
            c_jpp(:) = npscore_pair(rid:rid+1,angid,no1,pno2)%val2
            arg(:) = cubic_spline(0.5d0, c_j, c_jpp, a, calc_g)

            ff(res1,res2) = ff(res1,res2) + arg(1)
            ff(res2,res1) = ff(res2,res1) + arg(1)
           
        !3. Dipole-dipole
        else if (ispolar(no1) .and. ispolar(no2)) then
            ! as if dipole-nondipole
            costheta = dot_product(pvp(:), -r_pq(:))
            angid = min(6,abs(int(3*(costheta+1.0d0))) + 1)

            c_j(:)   = npscore_pair(rid:rid+1,angid,no2,pno1)%val1
            c_jpp(:) = npscore_pair(rid:rid+1,angid,no2,pno1)%val2
            arg(:) = cubic_spline(0.5d0, c_j, c_jpp, a, calc_g)

            ff(res1,res2) = ff(res1,res2) + arg(1)
            ff(res2,res1) = ff(res2,res1) + arg(1)

            ! as if nondipole-dipole
            costheta = dot_product(pvq(:), r_pq(:))
            angid = min(6,abs(int(3*(costheta+1.0d0))) + 1)

            c_j(:)   = npscore_pair(rid:rid+1,angid,no1,pno2)%val1
            c_jpp(:) = npscore_pair(rid:rid+1,angid,no1,pno2)%val2
            arg(:) = cubic_spline(0.5d0, c_j, c_jpp, a, calc_g)

            ff(res1,res2) = ff(res1,res2) + arg(1)
            ff(res2,res1) = ff(res2,res1) + arg(1)
           
            ! as if dipole-dipole
            costheta = dot_product(pvp(:), pvq(:))
            angid = min(6,abs(int(3*(costheta+1.0d0))) + 1)

            c_j(:)   = ppscore_pair(rid:rid+1,angid,pno2,pno1)%val1
            c_jpp(:) = ppscore_pair(rid:rid+1,angid,pno2,pno1)%val2
            arg = cubic_spline(0.5d0, c_j, c_jpp, a, calc_g)

            ff(res1,res2) = ff(res1,res2) + arg(1)
            ff(res2,res1) = ff(res2,res1) + arg(1)
        end if
    end do
end do

end subroutine pairwise_ddfire_static_score
!-------------------------------------------------------------------------------
subroutine ddfire_score(crd, Edfire, mode)
!-------------------------------------------------------------------------------
! For simple dfire point calculation without interpolation
!-------------------------------------------------------------------------------
real(dp), intent(in) :: crd(3,tn%atom)
character(len=3), intent(in) :: mode
real(dp), intent(out) :: Edfire(0:4)
logical :: appl_respair(tn%residue,tn%residue)
real(dp) :: r_pq(3), pvp(3), pvq(3), dist, costheta
integer :: rid, no1, no2, pno1, pno2, i_pair, angid, atm1, atm2, ndfire
integer :: atm2_id(max_neigh2), n_pair
logical :: status

! First refresh dipole vectors for current coordinate
if (ddfire_add_scale > small_real) call build_polarv(R)
ndfire = 0

Edfire(0:4) = 0.0d0
appl_respair(:,:) = .true.

do atm1 = 1, tn%atom
    if (.not. is_using_atm(atm1)) cycle

    call get_dfire_pair(atm1, appl_respair, n_pair, atm2_id, status)
    if (.not. status) then
        Edfire(0:4) = max_energy
        return
    end if

    do i_pair = 1, n_pair
        atm2 = atm2_id(i_pair)

        no1 = dfire_atmid(atm1)
        no2 = dfire_atmid(atm2)
        
        r_pq(:) = crd(:,atm2) - crd(:,atm1)
        dist = sqrt(dot_product(r_pq,r_pq))
        
        r_pq(:) = r_pq(:)/dist

        rid = int(2.0d0*dist) + 1
        
        if (ispolar(no1)) then
            pno1 = atom2polar(no1)
            pvp(:) = polarv(:,atm1)
        end if
        
        if (ispolar(no2)) then
            pno2 = atom2polar(no2)
            pvq(:) = polarv(:,atm2)
        end if
        
        ! Distance-dependent term
        Edfire(1) = Edfire(1) + nnscore_pair(rid,no2,no1)%val1
        ndfire = ndfire + 1
        
        ! Dipole-nondipole
        if (ispolar(no1) .and. (.not. ispolar(no2))) then
            costheta = dot_product(pvp(:), -r_pq(:))
            angid = min(6,int(3*(costheta+1.0)) + 1)
            Edfire(2) = Edfire(2) + npscore_pair(rid,angid,no2,pno1)%val1
        
            ! Nondipole-dipole
        else if ((.not. ispolar(no1)) .and. ispolar(no2)) then
            costheta = dot_product(pvq(:), r_pq(:))
            angid = min(6,int(3*(costheta+1.0)) + 1)
            Edfire(2) = Edfire(2) + npscore_pair(rid,angid,no1,pno2)%val1
           
            ! Dipole-dipole
        else if (ispolar(no1) .and. ispolar(no2)) then
            costheta = dot_product(pvp(:), -r_pq(:))
            angid = min(6,int(3*(costheta+1.0)) + 1)
            Edfire(3) = Edfire(3) + npscore_pair(rid,angid,no2,pno1)%val1
           
            costheta = dot_product(pvq(:), r_pq(:))
            angid = min(6,int(3*(costheta+1.0)) + 1)
            Edfire(3) = Edfire(3) + npscore_pair(rid,angid,no1,pno2)%val1
           
            costheta = dot_product(pvp(:), pvq(:))
            angid = min(6,int(3*(costheta+1.0)) + 1)
            Edfire(4) = Edfire(4) + ppscore_pair(rid,angid,pno2,pno1)%val1
        end if
    end do
end do

Edfire(0) = sum(Edfire(1:4))

end subroutine ddfire_score
!-------------------------------------------------------------------------------
subroutine read_atomlist(atomfile, atomorder, polarorder)
!-------------------------------------------------------------------------------
! Read atomic indices and orders used for dDFIRE table
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: atomfile
character(len=6), intent(out) :: atomorder(192), polarorder(66)
integer :: f_unit, ioerror, num_word
character(len=len_fname) :: word(13)
character(len=1000) :: line
character(len=6) :: atomid
integer :: i, p, type, polarity

f_unit = 41
open(f_unit, file = trim(atomfile))

i = 0
p = 0
atomorder(:) = ''
polarorder(:) = ''
atom2polar(:) = 0
ispolar(:) = .false.

do
    read(f_unit,"(A1000)", iostat = ioerror) line
    if (ioerror < 0) exit

    if (line(1:4) /= 'ATOM') cycle

    call parse_longstring(line, num_word, word, 1000)
     
    atomid(1:3) = word(2)
    atomid(4:6) = word(3)
    read(word(4),"(I3)") type
    read(word(5),"(I3)") polarity
     
    if (type == 0) then
        i = i + 1
        atomorder(i) = atomid
        if (polarity == 1) then ! Polar atoms
            p = p + 1
            polarorder(p) = atomid
            atom2polar(i) = p
            ispolar(i) = .true.
        end if
    end if
end do
  
nn = i
np = p

close(f_unit)

end subroutine read_atomlist
!-------------------------------------------------------------------------------
subroutine read_ddfire_table(tablefile)
!----------------------------- ------------------------------------------------------------
! Read ddfire table and fill the score table array
!----------------------------- ------------------------------------------------------------
character(len=len_fname), intent(in) :: tablefile
integer :: f_unit, ioerror, num_word
character(len=len_fname) :: word(13)
character(len=1000) :: line
character(len=2) :: type
integer :: id1, id2, dist, ang
real(dp) :: score, angscore(6)

f_unit=42
nnscore_pair(:,:,:)%val1 = 0.0d0
npscore_pair(:,:,:,:)%val1 = 0.0d0
ppscore_pair(:,:,:,:)%val1 = 0.0d0

open(f_unit, file = trim(tablefile))

do
    read(f_unit,"(A1000)", iostat = ioerror) line
    if (ioerror < 0) exit

    if (line(1:2) == 'nn') then
        type = 'nn'
        cycle
    elseif (line(1:2) == 'np') then
        type = 'np'
        cycle
    elseif (line(1:2) == 'pp') then
        type = 'pp'
        cycle
    else if (line(1:3) == 'END') then
        exit
    end if

    call parse_longstring(line, num_word, word, 1000)

    if (type == 'nn') then         !nonpolar-nonpolar interaction
        read(word(3),"(I10)") id1
        read(word(6),"(I10)") id2
        read(word(7),"(I10)") dist
        read(word(8),"(F15.5)") score
!        nnscore(dist,id2,id1) = score
        nnscore_pair(dist,id2,id1)%val1 = score


    else
        read(word(3),"(I10)") id1
        read(word(6),"(I10)") id2
        read(word(7),"(I10)") dist
        read(word(8:13),*) angscore(1:6)
        do ang = 1, 6
            if (type == 'np') then     !nonpolar-dipole interaction
!              npscore(dist,ang,id2,id1) = angscore(ang)
                npscore_pair(dist,ang,id2,id1)%val1 = angscore(ang)
            elseif (type == 'pp') then !dipole-dipole interaction
!              ppscore(dist,ang,id2,id1) = angscore(ang)
                ppscore_pair(dist,ang,id2,id1)%val1 = angscore(ang)
            end if
        end do
    end if
end do

close(f_unit)

end subroutine read_ddfire_table
!-------------------------------------------------------------------------------
subroutine assign_connected()
!-------------------------------------------------------------------------------
! Assign connection information of dipole-related atoms
!-------------------------------------------------------------------------------

n_conid(:) = 0
connected_list(:,:) = ''

!CYS SG
n_conid(5) = 1
connected_list(1,5) = 'CB  '
!ASP OD1
n_conid(8) = 1
connected_list(1,8) = 'CG  '
!GLU OE1
n_conid(11) = 1
connected_list(1,11) = 'CD  '
!HIS ND1
n_conid(18) = 2
connected_list(1,18) = 'CG  '
connected_list(2,18) = 'CE1 '
!HIS NE2
n_conid(19) = 2
connected_list(1,19) = 'CE1 '
connected_list(2,19) = 'CD2 '
!LYS NZ
n_conid(24) = 1
connected_list(1,24) = 'CE  '
!ASN OD1
n_conid(31) = 1
connected_list(1,31) = 'CG  '
!ASN ND2
n_conid(32) = 1
connected_list(1,32) = 'CG  '
!GLN OE1
n_conid(36) = 1
connected_list(1,36) = 'CD  '
!GLN NE2
n_conid(37) = 1
connected_list(1,37) = 'CD  '
!ARG NE
n_conid(40) = 1
connected_list(1,40) = 'CZ  '
!ARG NH1
n_conid(41) = 2
connected_list(1,41) = 'CD  '
connected_list(2,41) = 'CZ  '
!SER OG
n_conid(44) = 1
connected_list(1,44) = 'CB  '
!THR OG1
n_conid(47) = 1
connected_list(1,47) = 'CB  '
!TRP NE1
n_conid(52) = 2
connected_list(1,52) = 'CD1 '
connected_list(2,52) = 'CD2 '
!TYR OH
n_conid(55) = 1
connected_list(1,55) = 'CZ  '
!PTR OH
n_conid(58) = 1
connected_list(1,58) = 'CZ  '
!SEP OG
n_conid(61) = 1
connected_list(1,61) = 'CB  '
!TPO OG1
n_conid(64) = 1
connected_list(1,64) = 'CB  '

end subroutine assign_connected
!-------------------------------------------------------------------------------
subroutine find_ddfire_connected(resno, n_con, con_atm, ref_res_no, id, atom_name, res_name)
!-------------------------------------------------------------------------------
! Get connected atom list for polar atoms
! The connection information is used to build polarizaition vector
!-------------------------------------------------------------------------------
character(len=4), intent(in) :: atom_name, res_name
integer, intent(in) :: resno, ref_res_no, id
integer, intent(out) :: n_con, con_atm(2)
integer :: atmno
character(len=6) :: error_mode

error_mode = 'ignore'
con_atm(:) = 0

if (trim(atom_name) == 'N') then
    if (res_name == 'CNME') then
        call find_atom_idx(resno, ref_res_no, 'CH3 ', atmno, error_mode)
    else
        call find_atom_idx(resno, ref_res_no, 'CA  ', atmno, error_mode)
    end if
    con_atm(1) = atmno
    n_con = 1
     
else if (trim(atom_name) == 'O' .or. trim(atom_name) == 'OXT') then
    call find_atom_idx(resno, ref_res_no, 'C   ', atmno, error_mode)
    con_atm(1) = atmno
    n_con = 1
     
else if (n_conid(id) /= 0) then
    call find_atom_idx(resno, ref_res_no, connected_list(1,id), atmno, error_mode)
    con_atm(1) = atmno
    n_con = n_conid(id)
    if (n_conid(id) == 2) then
        call find_atom_idx(resno, ref_res_no, connected_list(2,id), atmno, error_mode)
        con_atm(2) = atmno
    end if
else
    n_con = 0
end if

end subroutine find_ddfire_connected
!-------------------------------------------------------------------------------
subroutine build_polarv(crd)
!-------------------------------------------------------------------------------
! Build dipole vector for given coordinate
!-------------------------------------------------------------------------------
real(dp), intent(in) :: crd(3,tn%atom)
integer :: atmno, id
real(dp) :: r1(3), r2(3), norm_r1(3), norm_r2(3), addv(3), abs_addv
  
do atmno = 1, tn%atom
    if (.not. is_using_atm(atmno)) cycle
    id = atom2polar(dfire_atmid(atmno))

    ! If connected with one atom
    if (n_connected(atmno) == 1) then
        r1(:) = crd(:,atmno) - crd(:,con_atmno(1,atmno))
        polarv(:,atmno) = r1(:)/sqrt(dot_product(r1,r1))

    ! When connected with two atoms
    elseif (n_connected(atmno) == 2) then
        r1(:) = crd(:,atmno) - crd(:,con_atmno(1,atmno))
        r2(:) = crd(:,atmno) - crd(:,con_atmno(2,atmno))
        norm_r1(:) = r1(:)/sqrt(dot_product(r1(:),r1(:)))
        norm_r2(:) = r2(:)/sqrt(dot_product(r2(:),r2(:)))
        addv(:) = norm_r1(:) + norm_r2(:)
        abs_addv = sqrt(dot_product(addv(:),addv(:)))
        if (abs_addv < small_real) then
            polarv(:,atmno) = 0.0d0
        else
            polarv(:,atmno) = addv(:) / abs_addv
        endif
    else
        cycle
    end if

    ! when not normalized
    if (dot_product(polarv(:,atmno),polarv(:,atmno)) > 1.0001d0) then
        write(log_msg,"(A,I5,F8.3,2I4)") 'Error at building polar vector: ', &
              atmno, r1(:), con_atmno(1,atmno), n_connected(atmno)
        call log_p(log_msg, me=me, level=30)
    end if
end do

end subroutine build_polarv
!-------------------------------------------------------------------------------
subroutine setup_ddfire_misc(protein, atomorder)
!-------------------------------------------------------------------------------
! Assign atom index, type and build available pairlist, For MODELING
! Consider residues only for which appl_respair is turned on 
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
character(len=6), intent(in) :: atomorder(192)
character(len=4) :: resatmname(2,tn%residue,max_atm)
integer :: i_atm, i_res, id, ref_res_no
character(len=4) :: atm_name, res_name

!Build resatmname
do i_res = 1, protein%n_res
    ref_res_no = protein%residue(i_res)%res_type
    res_name = ref_res(ref_res_no)%res_name
    if ((.not. res_name == 'PTR') .and. &
        (.not. res_name == 'SEP') .and. &
        (.not. res_name == 'TPO')) &
        call convert_to_stdres(res_name)

    do i_atm = 1, ref_res(ref_res_no)%n_atm
        resatmname(1,i_res,i_atm) = res_name
        resatmname(2,i_res,i_atm) = ref_res(ref_res_no)%atom_name(i_atm)
    end do
end do

!Match identical pairs
dfire_atmid(:) = 0
is_using_atm(:) = .false.
do i_atm = 1, tn%stdatm
    res_name = resatmname(1,i_R(1,i_atm),i_R(2,i_atm))
    atm_name = resatmname(2,i_R(1,i_atm),i_R(2,i_atm))
    if (atm_name(1:1) == 'H') cycle !Ignore hydrogens
     
    ! dDFIRE has 9 atom types which are identical type to the other among 192
    call find_identical_pair(res_name, atm_name)
    call atomtype2id(res_name, atm_name, atomorder, id)
    dfire_atmid(i_atm) = id
    if (id /= 0) then
        is_using_atm(i_atm) = .true.
    end if
end do

end subroutine setup_ddfire_misc
!-------------------------------------------------------------------------------
subroutine find_identical_pair(res_type, atm_type)
!-------------------------------------------------------------------------------
! Return identical atom name with higher preference
!-------------------------------------------------------------------------------
character(len=4), intent(in) :: res_type
character(len=4), intent(inout) :: atm_type
  
if (trim(res_type) == 'ASP' .and. atm_type == 'OD2 ') then
    atm_type = 'OD1 '
else if (trim(res_type) == 'GLU' .and. atm_type == 'OE2 ') then
    atm_type = 'OE1 '
else if (trim(res_type) == 'PHE' .and. atm_type == 'CD2 ') then
    atm_type = 'CD1 '
else if (trim(res_type) == 'PHE' .and. atm_type == 'CE2 ') then
    atm_type = 'CE1 '
else if (trim(res_type) == 'LEU' .and. atm_type == 'CD2 ') then
    atm_type = 'CD1 '
else if (trim(res_type) == 'ARG' .and. atm_type == 'NH2 ') then
    atm_type = 'NH1 '
else if (trim(res_type) == 'VAL' .and. atm_type == 'CG2 ') then
    atm_type = 'CG1 '
else if (trim(res_type) == 'TYR' .and. atm_type == 'CD2 ') then
    atm_type = 'CD1 '
else if (trim(res_type) == 'TYR' .and. atm_type == 'CE2 ') then
    atm_type = 'CE1 '
else if (trim(res_type) == 'PTR' .and. atm_type == 'CD2 ') then
    atm_type = 'CD1 '
else if (trim(res_type) == 'PTR' .and. atm_type == 'CE2 ') then
    atm_type = 'CE1 '
else if (atm_type == 'O2P ') then
    atm_type = 'O1P '
else if (atm_type == 'O3P ') then
    atm_type = 'O1P '
end if
  
end subroutine find_identical_pair
!-------------------------------------------------------------------------------
subroutine get_dfire_pair(atm1, appl_respair, n_pair, atm2_id, status)
!-------------------------------------------------------------------------------
integer, intent(in) :: atm1
logical, intent(in) :: appl_respair(tn%residue,tn%residue)
integer, intent(out) :: n_pair, atm2_id(max_neigh2)
logical, intent(out) :: status
integer :: i_2, res1, atm2

n_pair = 0
res1 = i_R(1,atm1)
atm2_id(1:i_P(atm1)%n_Lpair) = 0
status = .true.

do i_2 = i_P(atm1)%pair_end_index(2)+1, i_P(atm1)%n_pair   ! 1-4 pair, nb
    atm2 = i_P(atm1)%i_pair(i_2)
    if (i_R(1,atm2) /= res1 .and. is_using_atm(atm2) .and. appl_respair(res1, i_R(1,atm2))) then
        n_pair = n_pair + 1
        if (n_pair > max_neigh2) then
            status = .false.
            return
        end if
        atm2_id(n_pair) = atm2
    end if
end do

do i_2 = 1, i_P(atm1)%n_Lpair
    atm2 = i_P(atm1)%i_Lpair(i_2)
    if (i_R(1,atm2) /= res1 .and. is_using_atm(atm2) .and. appl_respair(res1, i_R(1,atm2))) then
        n_pair = n_pair + 1
        if (n_pair > max_neigh2) then
            status = .false.
            return
        end if
        atm2_id(n_pair) = atm2
    end if
end do

end subroutine get_dfire_pair
!-------------------------------------------------------------------------------
subroutine find_atmno(atom_name,resno,atmtype,atmno)
!-------------------------------------------------------------------------------
character(len=4), intent(in) :: atom_name(tn%atom), atmtype
integer, intent(in) :: resno
integer, intent(out) :: atmno
integer :: i_atm

do i_atm = 1, tn%stdatm
    if (i_R(1,i_atm) /= resno) cycle
    if (trim(atom_name(i_atm)) == trim(atmtype)) then
        atmno = i_atm
        exit
    end if
end do

end subroutine find_atmno
!-------------------------------------------------------------------------------
subroutine atomtype2id(restype, atmtype, atomorder, id)
!-------------------------------------------------------------------------------
character(len=4), intent(in) :: restype, atmtype
character(len=6), intent(in) :: atomorder(192)
integer, intent(out) :: id
integer :: i
logical :: selected

selected = .false.

do i = 1, nn
    if (atomorder(i)(1:3) == restype(1:3) .and. atomorder(i)(4:6) == atmtype(1:3)) then
        id = i
        selected = .true.
        exit
    end if
end do

if (.not. selected) id = 0

end subroutine atomtype2id
!-------------------------------------------------------------------------------
END MODULE dDFIRE
!-------------------------------------------------------------------------------
