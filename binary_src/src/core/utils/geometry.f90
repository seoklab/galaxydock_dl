!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/core/utils/geometry.f90
!
! Description:
!  This module contains various geometric operators related to
!   - Coordinate conversion
!   - Angle & dihedral enumerations
!   - Rotation & quaternion related operations
!   - and so on.
!
!-------------------------------------------------------------------------------
MODULE GEOMETRY
!-------------------------------------------------------------------------------
use globals
use logger, only: terminate_with_error
use sort,   only: sort2
use ran,    only: random
use mathfunctions, only: unit_matrix, cross, v_norm, bound_ang, &
                         jacobi, eigsrt, R3_inverse_matrix

implicit none
public

CONTAINS
!===============================================================================
! Coordinate conversion
!===============================================================================
!-------------------------------------------------------------------------------
subroutine cartesian2internal(res1, res2, residue)
!-------------------------------------------------------------------------------
! Update internal coordinate according to current Cartesian coordinate
! from res1 to res2 (when residue(:) is given)
! Input should be <residue_type> 
!
! Remark by Chaok : 
! Somehow, the following declaration did not work.
! I don't understand why!
! Therefore, res1 = 1 can be used in this setting.
!  type(residue_type), intent(inout) :: residue(res1:res2)
!-------------------------------------------------------------------------------
integer, intent(in) :: res1, res2
type(residue_type), intent(inout) :: residue(:)
integer :: i_res, i_ref_res, i_atm, b_idx(3)
real(dp), dimension(3,3) :: rr
real(dp), dimension(3) :: b
real(dp) :: angle
logical :: bnd_read(3)

call calc_bond_vectors_in_molecule(res1, res2, residue)
  
do i_res = res1, res2
    i_ref_res = residue(i_res)%res_type

    do i_atm = 1, ref_res(i_ref_res)%n_atm

        ! calc bond lengths
        if (residue(i_res)%bnd_read(i_atm)) then
            b(:) = residue(i_res)%b(:,i_atm)
            residue(i_res)%b_len(i_atm) = sqrt(dot_product(b,b))
        else
            residue(i_res)%b_len(i_atm) = ref_res(i_ref_res)%b_len0(i_atm)
        end if

        ! calc bond angles
        b_idx(1:2) = ref_res(i_ref_res)%bnd_in_b_ang(1:2,i_atm)
        bnd_read(1:2) = residue(i_res)%bnd_read(b_idx(1:2))
        if (bnd_read(1) .and. bnd_read(2)) then
            rr(:,1:2) = residue(i_res)%b(:,b_idx(1:2))
            call calc_bnd_angle(rr(:,1), rr(:,2), angle)
            residue(i_res)%b_ang(i_atm) = angle
        else
            residue(i_res)%b_ang(i_atm) = ref_res(i_ref_res)%b_ang0(i_atm)
        end if

        ! calc torsion angles
        b_idx(1:3) = ref_res(i_ref_res)%bnd_in_t_ang(1:3,i_atm)
        bnd_read(1:3) = residue(i_res)%bnd_read(b_idx(1:3))
        if (bnd_read(1) .and. bnd_read(2) .and. bnd_read(3)) then
            rr(:,1:3) = residue(i_res)%b(:,b_idx(1:3))
            call calc_tor_angle(rr(:,1), rr(:,2), rr(:,3), angle)
            residue(i_res)%t_ang(i_atm) = angle
        else
            residue(i_res)%t_ang(i_atm) = ref_res(i_ref_res)%t_ang0(i_atm)
        end if
    end do
end do

end subroutine cartesian2internal
!-------------------------------------------------------------------------------
subroutine calc_bond_vectors_in_molecule(res1, res2, residue)
!-------------------------------------------------------------------------------
! Checks whether bond atoms are read or placed, and
! calculates bond vector (residue_type%b(:))
! Used only for cartesian2internal
!-------------------------------------------------------------------------------
integer, intent(in) :: res1, res2
type(residue_type), intent(inout) :: residue(:)
integer :: i_res, i_ref_res, i_bnd, atm_no(2)
integer :: res_no1(-1:1), atm_no1(-1:1), res_no2(-1:1), atm_no2(-1:1)
logical :: atm_read(1:2), bnd_read, first_bnd_read(2)
real(dp) :: first_bnd(3,2)
real(dp) :: r_curr(3,3), r_ref(3,3), dr(3)
real(dp) :: bnd, ang, rr(3,2)

do i_res = res1, res2
    i_ref_res = residue(i_res)%res_type

    ! regular bonds
    do i_bnd = 2, ref_res(i_ref_res)%n_atm
        atm_no(1:2) = ref_res(i_ref_res)%atm_in_bnd(1:2,i_bnd)
        atm_read(1:2) = residue(i_res)%atm_read(atm_no(1:2)) &
             .or. residue(i_res)%atm_placed(atm_no(1:2))
        if (atm_read(1) .and. atm_read(2)) then
            bnd_read = .true.
            residue(i_res)%b(:,i_bnd) = &
                residue(i_res)%R(:,atm_no(2)) & 
                - residue(i_res)%R(:,atm_no(1))
        else
            bnd_read = .false.
        end if
        residue(i_res)%bnd_read(i_bnd) = bnd_read
    end do

    ! links
    if (i_res == 1) then
        ! first residue
        first_bnd_read(1:2) = residue(i_res)%bnd_read(ref_res(i_ref_res)%i_bnd_prev(1:2))

        if (first_bnd_read(1) .and. first_bnd_read(2)) then
            ! take N-CA and CA-C bonds
            first_bnd(:,1) = residue(i_res)%b(:,ref_res(i_ref_res)%i_bnd_prev(1))
            first_bnd(:,2) = residue(i_res)%b(:,ref_res(i_ref_res)%i_bnd_prev(2))
            ! For N terminal residue, assign bond vectors for the previous atoms (dummy?) 
            ! using the other known vectors 
            residue(i_res)%b(:,-1) = first_bnd(:,2) ! -N -> -CA
            residue(i_res)%b(:,0) = first_bnd(:,1)  !-CA ->  -C
            residue(i_res)%b(:,1) = first_bnd(:,2)  ! -C ->   N
        else
            ! arbitrary
            residue(i_res)%b(:,-1) = (/ 1.0d0, 0.0d0, 0.0d0 /)
            residue(i_res)%b(:,0) = (/ cos(60.0d0*deg2rad), sin(60.0d0*deg2rad), 0.0d0 /)
            residue(i_res)%b(:,1) = (/ 1.0d0, 0.0d0, 0.0d0 /)
        end if

        do i_bnd = -1, 1
            residue(i_res)%bnd_read(i_bnd) = .true.
        end do

        ! Calc R(0) of the residue (-C) using R of N and vecotr -C -> N
        residue(i_res)%R(:,0) = residue(i_res)%R(:,1) - residue(i_res)%b(:,1)

        ! Calc quaternion that locates the body in space.
        ! It is quat(0) of the first residue.
        ! r_ref is on the xy plane, and r_curr is the real coordinate.

        ! r_curr
        r_curr(:,3) = residue(i_res)%R(:,0)
        r_curr(:,2) = r_curr(:,3) - residue(i_res)%b(:,0)
        r_curr(:,1) = r_curr(:,2) - residue(i_res)%b(:,-1)
        
        ! r_ref(2)
        r_ref(:,2) = 0.0d0

        ! r_ref(3)
        dr(:) = residue(i_res)%b(:,0)
        bnd = sqrt(dot_product(dr(:), dr(:)))
        r_ref(:,3) = (/ bnd, 0.0d0, 0.0d0 /)

        ! r_ref(1)
        dr(:) = residue(i_res)%b(:,-1)
        bnd = sqrt(dot_product(dr(:), dr(:)))
        rr(:,1) = residue(i_res)%b(:,0)
        rr(:,2) = residue(i_res)%b(:,-1)
        call calc_bnd_angle(rr(:,1), rr(:,2), ang)
        r_ref(:,1) = bnd * (/ cos(ang), sin(ang), 0.0d0 /)
        
        ! calc quat
        call calc_body_quat(r_ref, r_curr, residue(i_res)%quat(:,0))
        
    else
        ! previous residue exists
        res_no1(-1:1) = residue(i_res)%link_res_no(1:3)
        atm_no1(-1:1) = residue(i_res)%link_atm_no(1:3)
        res_no2(-1:0) = res_no1(0:1)
        atm_no2(-1:0) = atm_no1(0:1)
        res_no2(1) = i_res
        atm_no2(1) = 1 ! atom no. 1 is the one linked to the previous residue
        do i_bnd = -1, 1
            atm_read(1) = residue(res_no1(i_bnd))%atm_read(atm_no1(i_bnd)) &
                .or. residue(res_no1(i_bnd))%atm_placed(atm_no1(i_bnd))
            atm_read(2) = residue(res_no2(i_bnd))%atm_read(atm_no2(i_bnd)) &
                .or. residue(res_no2(i_bnd))%atm_placed(atm_no2(i_bnd))
            if (atm_read(1) .and. atm_read(2)) then
                bnd_read = .true.
                residue(i_res)%b(:,i_bnd) = &
                   residue(res_no2(i_bnd))%R(:,atm_no2(i_bnd)) &
                   - residue(res_no1(i_bnd))%R(:,atm_no1(i_bnd))
            else
                bnd_read = .false.
            end if
            residue(i_res)%bnd_read(i_bnd) = bnd_read
        end do
    end if
end do

end subroutine calc_bond_vectors_in_molecule
!-------------------------------------------------------------------------------
subroutine calc_body_quat(r_ref, r_curr, quat)
!-------------------------------------------------------------------------------
! Calculate quaterion between cartesian system of r_curr(:,:)
! about r_ref(:,:)
! Used only for calc_bond_vectors_in_XXX
!-------------------------------------------------------------------------------
real(dp), intent(in) :: r_ref(3,3), r_curr(3,3)
real(dp), intent(out) :: quat(4)
real(dp) :: r(3,3), r0(3,3)
real(dp) :: ang, rr(3,2), p(3), p0(3)
real(dp) :: axis(3), axis_size
integer :: axis_type, vec_type
real(dp) :: quat1(4), quat2(4), U(3,3), s(3)

! translate so that r(3) is at the origin
r(:,1) = r_curr(:,1) - r_curr(:,3)
r(:,2) = r_curr(:,2) - r_curr(:,3)

! r0's
r0(:,1) = r_ref(:,1) - r_ref(:,3)
r0(:,2) = r_ref(:,2) - r_ref(:,3)

! quaternion that overlays [r_ref(3), r_ref(2)] on [r_curr(3), r_curr(2)]
! It can be problematic when vector from N to CA is (+-1.0, 0.0, 0.0). => solved by mink
call cross(r0(:,2), r(:,2), axis)
axis_size = dot_product(axis,axis)
if (axis_size > small_real) then
    axis = axis/sqrt(axis_size)
    rr(:,1) = r0(:,2)
    rr(:,2) = -r(:,2)
    call calc_bnd_angle(rr(:,1), rr(:,2), ang)
else
    axis = (/1.0d0, 0.0d0, 0.0d0/)
    ! ang is 0 or 180degrees
    axis_size = dot_product(r0(:,2), r(:,2))
    if (axis_size < small_real) then
        ang = pi
    else
        ang = 0.0d0
    end if
end if
call quaternion(axis, ang, quat1)

! rotation by quat1
vec_type = 3
call rotation_matrix(vec_type, quat1, U)
r0(:,1) = matmul(U, r0(:,1))

! quaternion that overlays [r_ref(1), r_ref(2)] on [r_curr(1), r_curr(2)]
s(:) = r(:,2) * dot_product(r(:,1),r(:,2))/dot_product(r(:,2),r(:,2))
p(:) = r(:,1) - s(:)
p0(:) = r0(:,1) - s(:)

call cross(p0, p, axis)
axis_size = dot_product(axis,axis)
if (axis_size > small_real) then
    axis = axis/sqrt(axis_size)
    rr(:,1) = p0(:)
    rr(:,2) = -p(:)
    call calc_bnd_angle(rr(:,1), rr(:,2), ang)
else
    axis = (/1.0d0, 0.0d0, 0.0d0/)
    ! ang is 0 or 180degrees
    axis_size = dot_product(p0, p)
    if (axis_size < small_real) then
        ang = pi
    else
        ang = 0.0d0
    end if
end if
call quaternion(axis, ang, quat2)

! product of the two quat
axis_type = 3
call q_product(axis_type, quat2, quat1, quat)

end subroutine calc_body_quat
!-------------------------------------------------------------------------------
subroutine internal2cartesian(res1, res2, residue, frameres)
!-------------------------------------------------------------------------------
! Update internal coordinate according to current Cartesian coordinate
! from res1 to res2 (when residue(:) is given)
! Input should be <residue_type> 
!
! Note that those bonds that close rings are not considered here,
! but must be considered in energy calculation
!
! Note that those bonds that close rings are not considered here,
! but must be considered in energy calculation
!
! Following algorithm is based on the fact that multiple rotation represented
! by quaternions can be expressed like this:
!       q' = q2*q1 (* means Hamilton product. please refer wikipedia)
!         Rotating system using q' is equivalent to rotating system using q1
!         followed by q2.
!
! It starts from fully extended (all backbone atoms are located along the
! x-axis) conformation and determines location of target atoms like following:
!   1. The target atom is first rotated by bond angle about z-axis.
!        q(b_ang) = (sin(b_ang/2), 0, 0, cos(b_ang/2))
!   2. Then, target atom is rotated by torsion angle about x-axis.
!        q(t_ang) = (cos(t_ang/2), sin(t_ang/2), 0, 0)
!   3. Finally, target atom is rotated using quaternion of linked atom (i_atm_prev).
! When quaterion of linked atom is q_prev, the quaternion equivalent to
! rotation operations described above is following:
!        q_curr = q_prev * q(t_ang) * q(b_ang) (calculated using do-loop, k)
!
! reference:
! C. Seok, E. A. Coutsias, Efficiency of Rotational Operators for Geometric Manipulation of Chain Molecules, Bull. Korean. Chem.
! Soc., 2007
!-------------------------------------------------------------------------------
integer, intent(in) :: res1, res2
type(residue_type), intent(inout) :: residue(:)
integer, optional :: frameres
integer :: i_res, i_ref_res, i_atm, k
integer :: i_res_prev, i_atm_prev
real(dp), dimension(4) :: q
real(dp), dimension(3) :: U, b
real(dp) :: angle, quarter_angle, tan_w, tan_sqr, tan1, cosine, sine
real(dp) :: q1, q2, q3, q4, p1, p2, p4
real(dp) :: b0(3,3), b1(3,3), rtmp(3), U2(3,3), c0(3), c1(3), dot

if (present(frameres)) then
    !b0(1,:) = residue(frameres)%b(:,1)
    b0(:,1) = residue(frameres)%b(:,1)
    call v_norm(b0(:,1))
     
    b0(:,2) = residue(frameres)%b(:,2)
    dot = dot_product(b0(:,2),b0(:,1))
    b0(:,2) = b0(:,2) - dot*b0(:,1)
    call v_norm(b0(:,2))

    call cross(b0(:,1), b0(:,2), rtmp(:))
    b0(:,3) = rtmp(:)

    c0(:) = residue(frameres)%R(:,0)
end if

do i_res = res1, res2
    i_ref_res = residue(i_res)%res_type

    ! copy R and quat of the last atm in the linked residue
    if (i_res > res1) then
        i_atm_prev = residue(i_res)%link_atm_no(3)
        i_res_prev = i_res - 1
        residue(i_res)%R(:,0) = residue(i_res_prev)%R(:,i_atm_prev)
        residue(i_res)%quat(:,0) = residue(i_res_prev)%quat(:,i_atm_prev)
    end if

    do i_atm = 1, ref_res(i_ref_res)%n_atm
        
        do k = 1, 2
            if (k == 1) then ! t_ang rotation
                angle = residue(i_res)%t_ang(i_atm)
            else ! b_ang rotation
                angle = residue(i_res)%b_ang(i_atm)
            end if

            quarter_angle = 0.25d0*angle
            tan_w = tan(quarter_angle)
            tan_sqr = tan_w * tan_w
            tan1 = 1.0d0 + tan_sqr
            cosine = (1.0d0 - tan_sqr)/tan1
            sine = 2.0d0*tan_w/tan1

            if (k == 1) then
                i_atm_prev = ref_res(i_ref_res)%atm_in_bnd(1,i_atm)
                q(1:4) = residue(i_res)%quat(:,i_atm_prev)
            end if
           
            q1 = q(1)
            q2 = q(2)
            q3 = q(3)
            q4 = q(4)

            if (k == 1) then
                p1 = cosine
                p2 = sine

                q(1) = q1*p1 - q2*p2
                q(2) = q1*p2 + p1*q2
                q(3) = p1*q3 + q4*p2
                q(4) = p1*q4 - q3*p2
            else
                p1 = sine
                p4 = cosine

                q(1) = q1*p1 - q4*p4
                q(2) = p1*q2 + q3*p4
                q(3) = p1*q3 - q2*p4
                q(4) = p1*q4 + q1*p4

                residue(i_res)%quat(:,i_atm) = q(:)

                q1 = q(1)
                q2 = q(2)
                q3 = q(3)
                q4 = q(4)

                ! 1/2*U
                U(1) = q1*q1 - 0.5d0 + q2*q2
                U(2) = q2*q3 + q1*q4
                U(3) = q2*q4 - q1*q3

                b(:) = 2.0d0*residue(i_res)%b_len(i_atm)*U(:)
                !
                residue(i_res)%b(:,i_atm) = b(:)
                residue(i_res)%R(:,i_atm) = &
                   residue(i_res)%R(:,i_atm_prev) + b(:)
            end if
        end do
    end do
end do

if (present(frameres)) then
    call calc_bond_vectors_in_molecule(frameres-1, frameres, residue)
    residue(frameres)%R(:,0) = residue(frameres-1)%R(:,residue(frameres)%link_atm_no(3))
    !
    b1(:,1) = residue(frameres)%b(:,1)
    call v_norm(b1(:,1))
     
    b1(:,2) = residue(frameres)%b(:,2)
    dot = dot_product(b1(:,2),b1(:,1))
    b1(:,2) = b1(:,2) - dot*b1(:,1)
    call v_norm(b1(:,2))

    call cross(b1(:,1), b1(:,2), rtmp(:))
    b1(:,3) = rtmp(:)

    c1(:) = residue(frameres)%R(:,0)
     
    call R3_inverse_matrix(b1(:,:))
    U2(:,:) = matmul(b0(:,:), b1(:,:))

    do i_res = res1, res2
        i_ref_res = residue(i_res)%res_type

        do i_atm = 1, ref_res(i_ref_res)%n_atm
            rtmp(:) = residue(i_res)%b(:,i_atm)
            call rotation(rtmp(:), U2(:,:))
            residue(i_res)%b(:,i_atm) = rtmp(:)

            rtmp(:) = residue(i_res)%R(:,i_atm) - c1(:)
            call rotation(rtmp(:), U2(:,:))
            residue(i_res)%R(:,i_atm) = rtmp(:) + c0(:)
        end do
    end do
end if

end subroutine internal2cartesian
!-------------------------------------------------------------------------------
subroutine cartesian2internal_het(hetmol)
!-------------------------------------------------------------------------------
! Cartesian to internal, for hetero-molecule
!-------------------------------------------------------------------------------
type(residue_type), intent(inout) :: hetmol
integer :: i_ref_res, i_atm, b_idx(3)
real(dp), dimension(3,3) :: rr
real(dp), dimension(3) :: b
real(dp) :: angle
logical :: bnd_read(3)

call calc_bond_vectors_in_hetmol(hetmol)
  
i_ref_res = hetmol%res_type

do i_atm = 1, ref_res(i_ref_res)%n_atm

    ! calc bond lengths
    if (hetmol%bnd_read(i_atm)) then
        b(:) = hetmol%b(:,i_atm)
        hetmol%b_len(i_atm) = sqrt(dot_product(b,b))
    else
        hetmol%b_len(i_atm) = ref_res(i_ref_res)%b_len0(i_atm)
    end if

    ! calc bond angles
    b_idx(1:2) = ref_res(i_ref_res)%bnd_in_b_ang(1:2,i_atm)
    bnd_read(1:2) = hetmol%bnd_read(b_idx(1:2))
    if (bnd_read(1) .and. bnd_read(2)) then
        rr(:,1:2) = hetmol%b(:,b_idx(1:2))
        call calc_bnd_angle(rr(:,1), rr(:,2), angle)
        hetmol%b_ang(i_atm) = angle
    else
        hetmol%b_ang(i_atm) = ref_res(i_ref_res)%b_ang0(i_atm)
    end if

    ! calc torsion angles
    b_idx(1:3) = ref_res(i_ref_res)%bnd_in_t_ang(1:3,i_atm)
    bnd_read(1:3) = hetmol%bnd_read(b_idx(1:3))
    if (bnd_read(1) .and. bnd_read(2) .and. bnd_read(3)) then
        rr(:,1:3) = hetmol%b(:,b_idx(1:3))
        call calc_tor_angle(rr(:,1), rr(:,2), rr(:,3), angle)
        hetmol%t_ang(i_atm) = angle
    else
        hetmol%t_ang(i_atm) = ref_res(i_ref_res)%t_ang0(i_atm)
    end if
end do

end subroutine cartesian2internal_het
!-------------------------------------------------------------------------------
subroutine calc_bond_vectors_in_hetmol(hetmol)
!-------------------------------------------------------------------------------
! hetmol version of calc_bond_vectors_in_molecule
!-------------------------------------------------------------------------------
type(residue_type), intent(inout) :: hetmol
integer :: i_ref_res, i_bnd, atm_no(2)
logical :: atm_read(1:2), bnd_read
real(dp) :: r_curr(3,3), r_ref(3,3), dr(3)
real(dp) :: bnd, ang, rr(3,2)

i_ref_res = hetmol%res_type

! regular bonds
do i_bnd = 2, ref_res(i_ref_res)%n_atm
    atm_no(1:2) = ref_res(i_ref_res)%atm_in_bnd(1:2,i_bnd)
    atm_read(1:2) = hetmol%atm_read(atm_no(1:2)) .or. hetmol%atm_placed(atm_no(1:2))
    if (atm_read(1) .and. atm_read(2)) then
        bnd_read = .true.
        hetmol%b(:,i_bnd) = hetmol%R(:,atm_no(2)) - hetmol%R(:,atm_no(1))
    else
        bnd_read = .false.
    end if
    hetmol%bnd_read(i_bnd) = bnd_read
end do
  
! arbitrary
hetmol%b(:,-1) = (/ 1.0d0, 0.0d0, 0.0d0 /)
hetmol%b(:,0) = (/ cos(60.0d0*deg2rad), sin(60.0d0*deg2rad), 0.0d0 /)
hetmol%b(:,1) = (/ 1.0d0, 0.0d0, 0.0d0 /)

do i_bnd = -1, 1
    hetmol%bnd_read(i_bnd) = .true.
end do

! Calc R(0) of the residue
hetmol%R(:,0) = hetmol%R(:,1) - hetmol%b(:,1)

! Calc quaternion that locates the body in space.
! It is quat(0) of the first residue.
! r_ref is on the xy plane, and r_curr is the real coordinate.

! r_curr
r_curr(:,3) = hetmol%R(:,0)
r_curr(:,2) = r_curr(:,3) - hetmol%b(:,0)
r_curr(:,1) = r_curr(:,2) - hetmol%b(:,-1)

! r_ref(2)
r_ref(:,2) = 0.0d0

! r_ref(3)
dr(:) = hetmol%b(:,0)
bnd = sqrt(dot_product(dr(:), dr(:)))
r_ref(:,3) = (/ bnd, 0.0d0, 0.0d0 /)

! r_ref(1)
dr(:) = hetmol%b(:,-1)
bnd = sqrt(dot_product(dr(:), dr(:)))
rr(:,1) = hetmol%b(:,0)
rr(:,2) = hetmol%b(:,-1)
call calc_bnd_angle(rr(:,1), rr(:,2), ang)
r_ref(:,1) = bnd * (/ cos(ang), sin(ang), 0.0d0 /)

! calc quat
call calc_body_quat(r_ref, r_curr, hetmol%quat(:,0))

end subroutine calc_bond_vectors_in_hetmol
!-------------------------------------------------------------------------------
subroutine internal2cartesian_het(hetmol)
!-------------------------------------------------------------------------------
! Internal to cartesian, for hetero molecule
! Note that those bonds that close rings are not considered here,
! but must be considered in energy calculation
!-------------------------------------------------------------------------------
type(residue_type), intent(inout) :: hetmol
integer :: i_ref_res, i_atm, k
integer :: i_atm_prev
real(dp) :: q(4), U(3), b(3)
real(dp) :: angle, quarter_angle, tan_w, tan_sqr, tan1, cosine, sine
real(dp) :: q1, q2, q3, q4, p1, p2, p4

i_ref_res = hetmol%res_type
  
do i_atm = 1, ref_res(i_ref_res)%n_atm
        
    do k = 1, 2
        if (k == 1) then ! t_ang rotation
            angle = hetmol%t_ang(i_atm)
        else ! b_ang rotation
            angle = hetmol%b_ang(i_atm)
        end if

        quarter_angle  = 0.25*angle
        tan_w = tan(quarter_angle)
        tan_sqr = tan_w * tan_w
        tan1 = 1.0d0 + tan_sqr
        cosine = (1.0d0 - tan_sqr)/tan1
        sine = 2.0d0*tan_w/tan1

        if (k == 1) then
            i_atm_prev = ref_res(i_ref_res)%atm_in_bnd(1,i_atm)
            q(1:4) = hetmol%quat(:,i_atm_prev)
        end if
        
        q1 = q(1)
        q2 = q(2)
        q3 = q(3)
        q4 = q(4)

        if (k == 1) then
            p1 = cosine
            p2 = sine
           
            q(1) = q1*p1 - q2*p2
            q(2) = q1*p2 + p1*q2
            q(3) = p1*q3 + q4*p2
            q(4) = p1*q4 - q3*p2
        else
            p1 = sine
            p4 = cosine
           
            q(1) = q1*p1 - q4*p4
            q(2) = p1*q2 + q3*p4
            q(3) = p1*q3 - q2*p4
            q(4) = p1*q4 + q1*p4
           
            hetmol%quat(:,i_atm) = q(:)
           
            q1 = q(1)
            q2 = q(2)
            q3 = q(3)
            q4 = q(4)
           
            ! 1/2*U
            U(1) = q1*q1 - 0.5d0 + q2*q2
            U(2) = q2*q3 + q1*q4
            U(3) = q2*q4 - q1*q3
           
            b(:) = 2.0d0*hetmol%b_len(i_atm)*U(:)
            hetmol%b(:,i_atm) = b(:)
            hetmol%R(:,i_atm) = &
                hetmol%R(:,i_atm_prev) + b(:)
        end if
    end do

end do

end subroutine internal2cartesian_het
!-------------------------------------------------------------------------------
subroutine internal2cartesian_sc(residue)
!-------------------------------------------------------------------------------
! Update internal coordinate according to current Cartesian coordinate
! for single residue sidechain atoms
! Input should be <residue_type> 
!-------------------------------------------------------------------------------
type(residue_type), intent(inout) :: residue
integer :: i_ref_res, i_atm, k, i_atm_prev
real(dp), dimension(4) :: q
real(dp), dimension(3) :: U, b
real(dp) :: angle, quarter_angle, tan_w, tan_sqr, tan1, cosine, sine
real(dp) :: q1, q2, q3, q4, p1, p2, p4

i_ref_res = residue%res_type

do i_atm = 1, ref_res(i_ref_res)%n_atm
    if (.not. ref_res(i_ref_res)%is_sc_atom(i_atm)) cycle
     
    do k = 1, 2
        if (k == 1) then ! t_ang rotation
            angle = residue%t_ang(i_atm)
        else ! b_ang rotation
            angle = residue%b_ang(i_atm)
        end if
        
        quarter_angle  = 0.25d0*angle
        tan_w = tan(quarter_angle)
        tan_sqr = tan_w * tan_w
        tan1 = 1.0d0 + tan_sqr
        cosine = (1.0d0 - tan_sqr)/tan1
        sine = 2.0d0*tan_w/tan1
        
        if (k == 1) then
            i_atm_prev = ref_res(i_ref_res)%atm_in_bnd(1,i_atm)
            q(1:4) = residue%quat(:,i_atm_prev)
        end if
        
        q1 = q(1)
        q2 = q(2)
        q3 = q(3)
        q4 = q(4)
        
        if (k == 1) then
            p1 = cosine
            p2 = sine
           
            q(1) = q1*p1 - q2*p2
            q(2) = q1*p2 + p1*q2
            q(3) = p1*q3 + q4*p2
            q(4) = p1*q4 - q3*p2
        else
            p1 = sine
            p4 = cosine
           
            q(1) = q1*p1 - q4*p4
            q(2) = p1*q2 + q3*p4
            q(3) = p1*q3 - q2*p4
            q(4) = p1*q4 + q1*p4
           
            residue%quat(:,i_atm) = q(:)
           
            q1 = q(1)
            q2 = q(2)
            q3 = q(3)
            q4 = q(4)
           
            ! 1/2*U
            U(1) = q1*q1 - 0.5d0 + q2*q2
            U(2) = q2*q3 + q1*q4
            U(3) = q2*q4 - q1*q3
           
            b(:) = 2.0d0*residue%b_len(i_atm)*U(:)
            residue%b(:,i_atm) = b(:)
            residue%R(:,i_atm) = residue%R(:,i_atm_prev) + b(:)
        end if
    end do
end do

end subroutine internal2cartesian_sc
!-------------------------------------------------------------------------------
subroutine ligand_to_topology(ligand)
!-------------------------------------------------------------------------------
! Defining bond lengths, bond angles, and torsion angles of input ligand.
!-------------------------------------------------------------------------------
type(ref_lig_type), intent(inout) :: ligand
integer :: ia, i, j, i_atm, j_atm, i_bnd, i_ang, i_tor
integer :: n_bnd, n_ang(2)
integer :: bnd_s(2,max_lig_atom*2)
integer :: ang_s(3,max_lig_atom*2, 2)
real(dp) :: R(3, 4), dR(3, 3), b_len, b_ang, d_ang

! Calculating Bond Lengths
do i_bnd = 1, ligand%n_bnd
    R(:,1) = ligand%R(:, ligand%bnd(2,i_bnd))
    R(:,2) = ligand%R(:, ligand%bnd(3,i_bnd))
    dR(:,1) = R(:,2) - R(:,1)
    b_len = sqrt(dot_product(dR(:,1), dR(:,1)))
    ligand%b_len0(i_bnd) = b_len
end do

! Defining Angles
i_ang = 0
do ia = 1, ligand%n_atm
    !i_atm = ligand%atom_no(ia)
    n_bnd = 0
    do i_bnd = 1, ligand%n_bnd
        if (ligand%bnd(2,i_bnd) == ia) then
            n_bnd = n_bnd + 1
            bnd_s(:,n_bnd) = (/ia, ligand%bnd(3,i_bnd)/)
        else if (ligand%bnd(3,i_bnd) == ia) then
            n_bnd = n_bnd + 1
            bnd_s(:,n_bnd) = (/ia, ligand%bnd(2,i_bnd)/)
        end if
    end do

    do i = 1, n_bnd-1
        do j = i+1, n_bnd
            i_ang = i_ang + 1
            ligand%ang(1, i_ang)   = i_ang
            ligand%ang(2:4, i_ang) = (/bnd_s(2,i), ia, bnd_s(2,j)/)
            !
            R(:,1) = ligand%R(:, ligand%ang(2,i_ang))
            R(:,2) = ligand%R(:, ligand%ang(3,i_ang))
            R(:,3) = ligand%R(:, ligand%ang(4,i_ang))
            dR(:,1) = R(:,2) - R(:,1)
            dR(:,2) = R(:,3) - R(:,2)
            call calc_bnd_angle(dR(:,1), dR(:,2), b_ang)
            ligand%b_ang0(i_ang) = b_ang
        end do
    end do
end do
ligand%n_ang = i_ang

! Defining Torsions
i_tor = 0
do ia = 1, ligand%n_bnd
    i_atm = ligand%bnd(2, ia)
    j_atm = ligand%bnd(3, ia)
    n_ang = 0
    do i_ang = 1, ligand%n_ang
        if (ligand%ang(2,i_ang) == i_atm .and. ligand%ang(3,i_ang) == j_atm) then
            n_ang(1) = n_ang(1) + 1
            ang_s(:,n_ang(1),1) = (/i_atm, j_atm, ligand%ang(4,i_ang)/)
        else if (ligand%ang(3,i_ang) == i_atm .and. ligand%ang(2,i_ang) == j_atm) then
            n_ang(2) = n_ang(2) + 1
            ang_s(:,n_ang(2),2) = (/i_atm, j_atm, ligand%ang(4,i_ang)/)
        else if (ligand%ang(3,i_ang) == i_atm .and. ligand%ang(4,i_ang) == j_atm) then
            n_ang(2) = n_ang(2) + 1
            ang_s(:,n_ang(2),2) = (/i_atm, j_atm, ligand%ang(2,i_ang)/)
        else if (ligand%ang(4,i_ang) == i_atm .and. ligand%ang(3,i_ang) == j_atm) then
            n_ang(1) = n_ang(1) + 1
            ang_s(:,n_ang(1),1) = (/i_atm, j_atm, ligand%ang(2,i_ang)/)
        end if
    end do

    do i = 1, n_ang(1)
        do j = 1, n_ang(2)
            i_tor = i_tor + 1
            ligand%dih(1, i_tor) = i_tor
            ligand%dih(2:5, i_tor) = (/ang_s(3,j,2), i_atm, j_atm, ang_s(3,i,1)/)
            !
            R(:,1) = ligand%R(:, ligand%dih(2,i_tor))
            R(:,2) = ligand%R(:, ligand%dih(3,i_tor))
            R(:,3) = ligand%R(:, ligand%dih(4,i_tor))
            R(:,4) = ligand%R(:, ligand%dih(5,i_tor))
            dR(:,1) = R(:,2) - R(:,1)
            dR(:,2) = R(:,3) - R(:,2)
            dR(:,3) = R(:,4) - R(:,3)
            call calc_tor_angle(dR(:,1), dR(:,2), dR(:,3), d_ang)
            ligand%d_ang0(i_tor) = d_ang
        end do
    end do
end do
ligand%n_dih = i_tor

end subroutine ligand_to_topology
!===============================================================================
! Angle and Dihedrals
!===============================================================================
!-------------------------------------------------------------------------------
subroutine calc_bnd_angle(crd1, crd2, angle)
!-------------------------------------------------------------------------------
! Simple calculation of bond angle:
! Input should be atom-connecting vectors
! crd1=Rab, crd2=Rbc
!-------------------------------------------------------------------------------
real(dp), intent(in) :: crd1(3), crd2(3)
real(dp), intent(out) :: angle
real(dp) :: p(3), q(3)
real(dp) :: arg

p(:) = -crd1(:)
q(:) = crd2(:)
arg = dot_product(p,q)/sqrt(dot_product(p,p)*dot_product(q,q))
arg = sign(min(abs(arg),1.0d0),arg) ! make sure abs(arg)<=1
angle = acos(arg)

end subroutine calc_bnd_angle
!-------------------------------------------------------------------------------
subroutine calc_tor_angle(crd1, crd2, crd3, angle)
!-------------------------------------------------------------------------------
! Simple calculation of torsion angle:
! Input should be atom-connecting vectors
! dihedral or improper torsion angle formed by abcd, i.e., 
! crd1(:) = Rab, crd2(:) = Rbc, crd3(:) = Rcd
!-------------------------------------------------------------------------------
real(dp), intent(in) :: crd1(3), crd2(3), crd3(3)
real(dp), intent(out) :: angle
real(dp) :: p(3), q(3), s(3)
real(dp) :: arg

call cross(crd1(:), crd2(:), p)
call cross(crd2(:), crd3(:), q)
call cross(crd3(:), crd1(:), s)
arg = dot_product(p,q)/sqrt(dot_product(p,p)*dot_product(q,q))
arg = sign(min(abs(arg),1.0d0),arg) ! make sure abs(arg)<=1
angle = sign(acos(arg), dot_product(s,crd2(:)))

end subroutine calc_tor_angle
!-------------------------------------------------------------------------------
subroutine calc_ang_and_grad(crd, angle, dadr, calc_g)
!-------------------------------------------------------------------------------
! Calculation of bond angles and derivatives about Cartesian coordinate
! Input should be atomic coordinates
! Returns angle connecting crd(:,1)==crd(:,2)==crd(:,3), and gradient (da/dr)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: crd(3,3)
real(dp), intent(out) :: angle, dadr(3,3)
logical, intent(in) :: calc_g
real(dp) :: p(3), q(3), norm_r1, norm_r2
real(dp) :: dadr_1(3), dadr_3(3), arg, tmp

dadr(:,:) = 0.0d0

p(:) = crd(:,1) - crd(:,2)
q(:) = crd(:,3) - crd(:,2)

norm_r1 = sqrt(dot_product(p(:),p(:)))
norm_r2 = sqrt(dot_product(q(:),q(:)))
  
arg = dot_product(p(:),q(:))/(norm_r1*norm_r2)
arg = sign(min(abs(arg),1.0d0),arg) ! make sure abs(arg)<=1
angle = acos(arg)

! calculate dfdr for feature angle
if (calc_g) then
    tmp = sqrt(1.0d0 - arg**2)
    dadr_1 = (p(:)*arg/norm_r1 - q(:)/norm_r2) / (norm_r1*tmp)
    dadr_3 = (q(:)*arg/norm_r2 - p(:)/norm_r1) / (norm_r2*tmp)
    dadr(:,1) = dadr_1(:)
    dadr(:,3) = dadr_3(:)
    dadr(:,2) = -dadr_1(:) - dadr_3(:)
end if

end subroutine calc_ang_and_grad
!-------------------------------------------------------------------------------
subroutine calc_tor_and_grad(crd, angle, dadr, calc_g)
!-------------------------------------------------------------------------------
! Calculation of dihedral angles and derivatives about Cartesian coordinate
! Input should be atomic coordinates
! Returns dihedral angle connecting crd(:,1)==crd(:,2)==crd(:,3)==crd(:,4), 
! and gradient (da/dr)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: crd(3,4)
real(dp), intent(out) :: angle, dadr(3,4)
logical, intent(in) :: calc_g
real(dp) :: r_tor(3,3), p(3), q(3), s(3)
real(dp) :: dadr_1(3), dadr_4(3), arg, proj1, proj2
real(dp) :: norm_r2, norm_r, norm_p2, norm_q2
  
dadr(:,:) = 0.0d0
  
r_tor(:,1) = crd(:,2) - crd(:,1)
r_tor(:,2) = crd(:,3) - crd(:,2)
r_tor(:,3) = crd(:,4) - crd(:,3)

call cross(r_tor(:,1), r_tor(:,2), p)
call cross(r_tor(:,2), r_tor(:,3), q)
call cross(r_tor(:,3), r_tor(:,1), s)
     
norm_r2 = dot_product(r_tor(:,2),r_tor(:,2))
norm_r = sqrt(norm_r2)
norm_p2 = dot_product(p(:),p(:))
norm_q2 = dot_product(q(:),q(:))
  
arg = dot_product(p(:),q(:))/(sqrt(norm_p2)*sqrt(norm_q2))
arg = sign(min(abs(arg),1.0d0),arg) ! make sure abs(arg)<=1
angle = sign(acos(arg), dot_product(s(:),r_tor(:,2)))

if (calc_g) then
    proj1 = dot_product(-r_tor(:,1),r_tor(:,2))/norm_r2
    proj2 = dot_product(-r_tor(:,3),r_tor(:,2))/norm_r2
     
    dadr_1(:) = -(norm_r/norm_p2)*p(:)
    dadr_4(:) = (norm_r/norm_q2)*q(:)
    dadr(:,1) = dadr_1(:)
    dadr(:,4) = dadr_4(:)
    dadr(:,2) = (proj1 - 1.0d0)*dadr_1(:) - proj2*dadr_4(:)
    dadr(:,3) = (proj2 - 1.0d0)*dadr_4(:) - proj1*dadr_1(:)
end if

end subroutine calc_tor_and_grad
!-------------------------------------------------------------------------------
subroutine calc_tor_deviation(n_ang, tor1, tor2, d_tor)
!-------------------------------------------------------------------------------
! Calculate the squared difference of n_ang successive torsion angles
!-------------------------------------------------------------------------------
integer, intent(in) :: n_ang
real(dp), intent(in) :: tor1(n_ang), tor2(n_ang)
real(dp), intent(out) :: d_tor
real(dp) :: dw, t(2), dt
integer :: i, j

dw = 0.0d0
do i = 1, n_ang
    do j = 1, 2
        if (j == 1) then
            t(j) = tor1(i)
        else
            t(j) = tor2(i)
        end if
        t(j) = bound_ang(t(j))
    end do
    dt = t(1) - t(2)
    dw = dw + bound_ang(dt)**2
end do
d_tor = dw

end subroutine calc_tor_deviation
!===============================================================================
! Rotation & Quaternion
!===============================================================================
!-------------------------------------------------------------------------------
subroutine quaternion(axis, angle, q)
!-------------------------------------------------------------------------------
! calculate quaternion, given rotation axis and angle.
! Original:
!    q(1) = cos(theta/2)
!    q(2:4) = axis(1:3) * sin(theta/2)
! Following algorithm is more complicated but faster version to reduce the
! number of calling trigonometric function.
!-------------------------------------------------------------------------------
real(dp), intent(in) :: axis(3), angle
real(dp), intent(out) :: q(4)
real(dp) :: tan_w, tan_sqr, tan1, cosine, sine

! cos(theta/2) = {1-tan_sqr(theta/4)}/{(1+tan_sqr(theta/4)}
tan_w = tan(0.25d0*angle)
tan_sqr = tan_w * tan_w
tan1 = 1.0d0 + tan_sqr
cosine = (1.0d0 - tan_sqr)/tan1
! sine(theta/2) = 2 * tan(theta/4)/(1+tan_sqr(theta/4))
sine = 2.0d0*tan_w/tan1
q(1) = cosine
q(2:4) = axis(1:3)*sine

end subroutine quaternion
!-------------------------------------------------------------------------------
subroutine slerp(q1, q2, t, q_out)
!-------------------------------------------------------------------------------
!spherical linear interpolation in the context of quaternion interpolation
!-------------------------------------------------------------------------------
real(dp), intent(in) :: q1(4), q2(4), t
real(dp), intent(out) :: q_out(4)
real(dp) :: omega, q1_normed(4), q2_normed(4), so

q1_normed = q1
call v_norm(q1_normed)
q2_normed = q2
call v_norm(q2_normed)

omega = acos(dot_product(q1_normed, q2_normed))
so = sin(omega)

q_out = sin((1.0d0-t)*omega)/so*q1 + sin(t*omega)/so*q2

end subroutine slerp
!-------------------------------------------------------------------------------
subroutine lerp(q1, q2, t, q_out)
!-------------------------------------------------------------------------------
!linear interpolation in the context of quaternion interpolation
!-------------------------------------------------------------------------------
real(dp), intent(in) :: q1(4), q2(4), t
real(dp), intent(out) :: q_out(4)

q_out = (1-t)*q1 + t*q2
call v_norm(q_out)

end subroutine lerp
!-------------------------------------------------------------------------------
subroutine average_quaternion(q_s, n_quat, q_out)
!-------------------------------------------------------------------------------
!Get an average (mean) from more than two quaternions (with two, slerp would be used).
!Note: this only works if all the quaternions are relatively close together.
!-------------------------------------------------------------------------------
real(dp), intent(in) :: q_s(:,:)
integer, intent(in) :: n_quat
real(dp), intent(out) :: q_out(4)
real(dp) :: q_1(4), q_i(4), temp, addamount
integer :: i_quat

q_out(1:4) = 0.0d0
q_1 = q_s(1:4,1)
addamount = 1.0d0/n_quat

do i_quat=1, n_quat
    q_i(1:4) = q_s(1:4,i_quat)

    temp = dot_product(q_1,q_i)
    if (temp < 0) q_i = -q_i

    q_out(1:4) = q_out(1:4) + addamount*q_i(1:4)
end do

call v_norm(q_out)
!-------------------------------------------------------------------------------
end subroutine average_quaternion
!-------------------------------------------------------------------------------
subroutine gen_random_axis(axis)
!-------------------------------------------------------------------------------
! Generates a random 3D unit vector with a uniform spherical distribution
!-------------------------------------------------------------------------------
real(dp), intent(out) :: axis(3)
real(dp) :: phi, costheta, theta

phi = two_pi*random()
costheta = 2.0d0*random()-1.0d0
theta = acos(costheta)
axis(1:3) = (/sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)/)

end subroutine gen_random_axis
!-------------------------------------------------------------------------------
subroutine gen_quaternion(max_pert_w, q)
!-------------------------------------------------------------------------------
! Randomly generate quaternion within given input angle.
!-------------------------------------------------------------------------------
real(dp), intent(in) :: max_pert_w
real(dp), intent(out) :: q(4)
integer :: i
real(dp) :: w, axis(3)
  
w = max_pert_w*(1.0d0*random()-0.5d0)
!do i = 1, 3
!    axis(i) = 2.0d0*random()-1.0d0
!end do
call gen_random_axis(axis)

!call v_norm(axis)
call quaternion(axis(:), w, q)

end subroutine gen_quaternion
!-------------------------------------------------------------------------------
subroutine gen_quaternion_btw_vectors(v1, v2, q)
!-------------------------------------------------------------------------------
! NOTE: v1, v2 should be normalized before call this subroutine!
!-------------------------------------------------------------------------------
real(dp), intent(in) :: v1(3), v2(3)
real(dp), intent(out) :: q(4)
!
real(dp) :: cos_theta, half_cos, half_sin, axis(3)
integer :: i

cos_theta = dot_product(v1, v2)
half_cos = sqrt(0.5d0*(1.0d0+cos_theta))
half_sin = sqrt(0.5d0*(1.0d0-cos_theta))
!
call cross(v1, v2, axis)
call v_norm(axis)
!
q(1) = half_cos
do i = 1, 3
    q(i+1) = half_sin*axis(i)
end do

end subroutine gen_quaternion_btw_vectors
!-------------------------------------------------------------------------------
subroutine pert_rotation(max_pert_w, q)
!-------------------------------------------------------------------------------
! Randomly perturb rotation within given input angle.
! Returns a quaternion that represents the perturbed rotation
!-------------------------------------------------------------------------------
real(dp), intent(in) :: max_pert_w
real(dp), intent(inout) :: q(4)
real(dp) :: q_pert(4), q0(4)

q0 = q
call gen_quaternion(max_pert_w, q_pert)
call q_product(3, q_pert, q0, q)

end subroutine pert_rotation
!-------------------------------------------------------------------------------
subroutine gen_uniform_random_quat(q)
!-------------------------------------------------------------------------------
! Randomly generate quaternion (uniformly spread in space.)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: q(4)
real(dp) :: s, s1, s2, theta1, theta2

s = random()
s1 = sqrt(1-s)
s2 = sqrt(s)
theta1 = 2*pi*random()
theta2 = 2*pi*random()
q(1) = cos(theta2)*s2
q(2) = sin(theta1)*s1
q(3) = cos(theta1)*s1
q(4) = sin(theta2)*s2

end subroutine gen_uniform_random_quat
!-------------------------------------------------------------------------------
subroutine sample_rotation_by_Fibonacci(q, n_sampl_z, n_sampl, init_v, randomize)
!-------------------------------------------------------------------------------
! Generate quaternions using spherical Fibonacci point sets
! P(phi, theta) = (cos(phi)sin(theta), sin(phi)sin(theta), cos(theta))
! SF_i = P(phi_i, arccos(z_i)) when phi_i = 2pi*frac(i/golden_ratio), 
!                                   z_i = 1 - (2i+1)/n
!  cf.) golden_ratio = (sqrt(5)+1)/2
!-------------------------------------------------------------------------------
real(dp), intent(out) :: q(:,:)
integer, intent(in) :: n_sampl_z
integer, intent(in) :: n_sampl
real(dp), intent(in) :: init_v(3)
logical, intent(in), optional :: randomize
!
real(dp), parameter :: golden_ratio = (sqrt(5.0d0)+1.0d0)/2.0d0
real(dp) :: r_ref(3), r_curr(3), i, phi, z, z_sqr, frac, sin_theta, rnd
real(dp) :: p(4), p_z(4), d_theta, theta
integer :: i_sampl, i_sampl_z, idx

rnd = random()
if (present(randomize) .and. (.not. randomize)) then
    rnd = 0.0d0
end if
d_theta = 360.0*deg2rad/n_sampl_z

r_ref(:) = init_v(:)
call v_norm(r_ref)
idx = 0
do i_sampl = 1, n_sampl
    i = dble(i_sampl-1)
    !
    frac = i/golden_ratio - floor(i/golden_ratio)
    phi = 2*pi*(frac+rnd)
    !
    z = 1.0d0 - (2.0d0*i + 1.0d0)/dble(n_sampl)
    z_sqr = z**2
    sin_theta = sqrt(1-z_sqr)
    !
    r_curr(1) = cos(phi) * sin_theta
    r_curr(2) = sin(phi) * sin_theta
    r_curr(3) = z
    !
    call v_norm(r_curr)
    call gen_quaternion_btw_vectors(r_ref, r_curr, p)
    do i_sampl_z = 1, n_sampl_z
        idx = idx+1
        theta = (i_sampl_z-1)*d_theta
        p_z = (/cos(theta/2), 0.0d0, 0.0d0, sin(theta/2)/)
        call q_product(2,p,p_z,q(:,idx))
    end do
end do

end subroutine sample_rotation_by_Fibonacci
!-------------------------------------------------------------------------------
subroutine q_product(axis_type, q, p, s)
!-------------------------------------------------------------------------------
! calculates product of two quaternions p and q.
! axis_type = 1 (p(3:4)=0), 2 (p(2:3)=0), 3 
!-------------------------------------------------------------------------------
integer, intent(in) :: axis_type
real(dp), intent(in) :: p(4), q(4)
real(dp), intent(out) :: s(4)
real(dp) :: p_i(3), q_i(3), qcp(3)
real(dp) :: p1, p2, p4, q1, q2, q3, q4

if (axis_type == 1) then ! p(3:4)=0
    p1 = p(1)
    p2 = p(2)
    q1 = q(1)
    q2 = q(2)
    q3 = q(3)
    q4 = q(4)
    s(1) = q1*p1 - q2*p2
    s(2) = q1*p2 + p1*q2
    s(3) = p1*q3 + q4*p2
    s(4) = p1*q4 - q3*p2
else if (axis_type == 2) then ! p(2:3)=0
    p1 = p(1)
    p4 = p(4)
    q1 = q(1)
    q2 = q(2)
    q3 = q(3)
    q4 = q(4)
    s(1) = q1*p1 - q4*p4
    s(2) = p1*q2 + q3*p4
    s(3) = p1*q3 - q2*p4
    s(4) = p1*q4 + q1*p4
else if (axis_type == 3) then
    q_i = q(2:4)
    p_i = p(2:4)
    s(1) = q(1)*p(1) - dot_product(q_i,p_i)
    call cross(q_i, p_i, qcp)
    s(2:4) = q(1)*p_i + p(1)*q_i + qcp
end if

end subroutine q_product
!-------------------------------------------------------------------------------
subroutine rotation_matrix(vec_type, q, U)
!-------------------------------------------------------------------------------
! constructs rotation matrix U from quaternion q.
! vec_type = 1 (the vector rotated by U is on the x axis), 3 (3-d)
! vec_type = 1 is unused. vec_type = 3 is only used.
!-------------------------------------------------------------------------------
integer, intent(in) :: vec_type
real(dp), intent(in) :: q(4)
real(dp), intent(out) :: U(3,3)
real(dp) :: q0, q1, q2, q3, b0, b1, b2, b3
real(dp) :: q00, q01, q02, q03, q11, q12, q13, q22, q23, q33

q0 = q(1)
q1 = q(2)
q2 = q(3)
q3 = q(4)
b0 = 2.0d0*q0
b1 = 2.0d0*q1
q00 = b0*q0 - 1.0d0
q02 = b0*q2
q03 = b0*q3
q11 = b1*q1
q12 = b1*q2
q13 = b1*q3
if (vec_type == 1) then ! need only U(:,1)
    U(1,1) = q00 + q11
    U(2,1) = q12 + q03
    U(3,1) = q13 - q02
else if (vec_type == 3) then
    b2 = 2.0d0*q2
    b3 = 2.0d0*q3
    q01 = b0*q1
    q22 = b2*q2
    q23 = b2*q3
    q33 = b3*q3
    U(1,1) = q00 + q11
    U(1,2) = q12 - q03
    U(1,3) = q13 + q02
    U(2,1) = q12 + q03
    U(2,2) = q00 + q22
    U(2,3) = q23 - q01
    U(3,1) = q13 - q02
    U(3,2) = q23 + q01
    U(3,3) = q00 + q33
end if

end subroutine rotation_matrix
!-------------------------------------------------------------------------------
subroutine rotation(crd, U)
!-------------------------------------------------------------------------------
! Rotate crd about rotation matrix U
! x'(3) = U(3,3) * x(3)
!-------------------------------------------------------------------------------
real(dp), intent(inout) :: crd(3)
real(dp), intent(in) :: U(3,3)
real(dp) :: s(3)

s(1:3) = crd(1:3)
  
crd(1) = dot_product(U(1, 1:3), s(1:3))
crd(2) = dot_product(U(2, 1:3), s(1:3))
crd(3) = dot_product(U(3, 1:3), s(1:3))

end subroutine rotation
!-------------------------------------------------------------------------------
subroutine calc_g_rotmatrix(q, T)
!-------------------------------------------------------------------------------
! Calculate derivative of rotation matrix about quaternion
!-------------------------------------------------------------------------------
real(dp), intent(in) :: q(4)
real(dp), intent(out) :: T(3,3,4)
  
T(:,1,1) = (/ q(1), -q(4),  q(3) /)
T(:,2,1) = (/ q(4),  q(1), -q(2) /)
T(:,3,1) = (/-q(3),  q(2),  q(1) /)

T(:,1,2) = (/ q(2),  q(3),  q(4) /)
T(:,2,2) = (/ q(3), -q(2), -q(1) /)
T(:,3,2) = (/ q(4),  q(1), -q(2) /)

T(:,1,3) = (/-q(3),  q(2),  q(1) /)
T(:,2,3) = (/ q(2),  q(3),  q(4) /)
T(:,3,3) = (/-q(1),  q(4), -q(3) /)

T(:,1,4) = (/-q(4), -q(1),  q(2) /)
T(:,2,4) = (/ q(1), -q(4),  q(3) /)
T(:,3,4) = (/ q(2),  q(3),  q(4) /)

T(:,:,:) = 2.0d0*T(:,:,:)

end subroutine calc_g_rotmatrix
!-------------------------------------------------------------------------------
subroutine inv_expmap_quat(q, p)
!-------------------------------------------------------------------------------
! Do inverse-exponential mapping : q(4) -> p(3)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: q(4)
real(dp), intent(out) :: p(3)
real(dp) :: p_abs, sinw, arg

arg = 1.0d0 - q(1)
p_abs = 2.0d0*acos(q(1))

if (arg < 1.0d-12/8.0d0) then
    p(1:3) = q(2:4)/(0.5d0+(p_abs**2)/48.0d0)
    !p(1:3) = q(2:4)*p_abs/(q(1)+(q(1)**3)/6.0d0)
else
    sinw = sqrt(1.0d0 - q(1)**2)
    p(1:3) = q(2:4)*p_abs/sinw
end if

if (p_abs > pi) then
    p(:) = (1.0d0 - 2.0*pi/p_abs)*p(:)
end if
  
end subroutine inv_expmap_quat
!-------------------------------------------------------------------------------
subroutine expmap_quat(p, q)
!-------------------------------------------------------------------------------
! Do exponential mapping : p(3) -> q(4)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: p(3)
real(dp), intent(out) :: q(4)
real(dp) :: w, p_abs, p_sqr

p_sqr = dot_product(p(:),p(:))
p_abs = sqrt(p_sqr)
w = cos(0.5d0*p_abs)
! Use Taylor expansion when abs approaches to numerical precision
if (dot_product(p(:),p(:)) < 1.0d-6) then
    q(1) = w
    q(2:4) = (0.5d0 + p_sqr/48.0d0)*p(:)
else
    q(1) = w
    q(2:4) = sqrt(1.0d0-w**2)*p(:)/p_abs
end if

end subroutine expmap_quat
!-------------------------------------------------------------------------------
subroutine euler_to_rotmat(p, U)
!-------------------------------------------------------------------------------
! Do from euler-angle to rotation-matrix : p(3) -> U(3,3)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: p(3)
real(dp), intent(out) :: U(3,3)
real(dp) :: c1, c2, c3, s1, s2, s3
  
c1 = cos(p(1))
c2 = cos(p(2))
c3 = cos(p(3))
s1 = sin(p(1))
s2 = sin(p(2))
s3 = sin(p(3))
U(1,1:3) = (/c1*c3-c2*s1*s3, -c1*s3-c2*c3*s1, s1*s2/)
U(2,1:3) = (/c3*s1+c1*c2*s3, c1*c2*c3-s1*s3, -c1*s2/)
U(3,1:3) = (/ s2*s3,  c3*s2,  c2 /)

end subroutine euler_to_rotmat
!-------------------------------------------------------------------------------
subroutine rotmat_to_euler(U, p)
!-------------------------------------------------------------------------------
! Do from rotation-matrix to euler-angle : U(3,3) -> p(3)
! p(1),p(3) range is -pi ~ +pi, p(2) range is 0 ~ +pi
!-------------------------------------------------------------------------------
real(dp), intent(in) :: U(3,3)
real(dp), intent(out) :: p(3)
real(dp) :: beta
  
p(1:3) = (/ atan2(U(1,3),-U(2,3)), acos(U(3,3)), atan2(U(3,1),U(3,2)) /)
beta = p(2)
if (beta < small_real) then ! Gimbal lock
    p(1:3) = (/ atan2(U(2,1),U(2,2)), 0.0d0, 0.0d0 /)
else if (abs(pi-beta) < small_real) then
    p(1:3) = (/ atan2(U(2,1),U(1,1)), pi, 0.0d0 /)
end if
  
end subroutine rotmat_to_euler
!-------------------------------------------------------------------------------
subroutine euler_to_quaternion(p, angle, axis)
!-------------------------------------------------------------------------------
! Do euler-quaternion : p(3) -> rotation axis with rotation angle
! Euler angle notation is given as Tait-Bryan angles
!-------------------------------------------------------------------------------
real(dp), dimension(3), intent(in) :: p
real(dp), dimension(3), intent(out) :: axis
real(dp), intent(out) :: angle
real(dp) :: phi, theta, psi
real(dp) :: q1, q2, q3, q0
real(dp) :: sin_angle

phi = p(1)
theta = p(2)
psi = p(3)

q0 = cos(phi/2.0d0) * cos(theta/2.0d0) * cos(psi/2.0d0)
q0 = q0 + sin(phi/2.0d0) * sin(theta/2.0d0) * sin(psi/2.0d0)

q1 = sin(phi/2.0d0) * cos(theta/2.0d0) * cos(psi/2.0d0)
q1 = q1 - cos(phi/2.0d0) * sin(theta/2.0d0) * sin(psi/2.0d0)

q2 = cos(phi/2.0d0) * sin(theta/2.0d0) * cos(psi/2.0d0)
q2 = q2 + sin(phi/2.0d0) * cos(theta/2.0d0) * sin(psi/2.0d0)

q3 = cos(phi/2.0d0) * cos(theta/2.0d0) * sin(psi/2.0d0)
q3 = q3 - sin(phi/2.0d0) * sin(theta/2.0d0) * cos(psi/2.0d0)

angle = 2.0d0 * acos(q0)
sin_angle = sin(angle/2.0d0)
axis(1) = q1/sin_angle
axis(2) = q2/sin_angle
axis(3) = q3/sin_angle

end subroutine euler_to_quaternion
!-------------------------------------------------------------------------------
subroutine rotmat_to_axisangle(U, theta, axis)
!-------------------------------------------------------------------------------
! Calculates rotation angle and axis from the rotational matrix
!-------------------------------------------------------------------------------
real(dp), intent(in) :: U(3,3)
real(dp), intent(out):: theta, axis(3)
real(dp) :: check_num
  
check_num = (U(1,1) + U(2,2) + U(3,3) -1.0d0)/2.0d0
if (check_num > 1.0d0) then
    check_num = 1.0d0
else if (check_num < -1.0d0) then
    check_num = -1.0d0
end if
theta = acos( check_num )
if ( (theta**2) < small_real ) then
    axis(:) = 1.0d0
else if ( ((theta-pi)**2) < small_real ) then
    axis(:) = -1.0d0
else
    axis(1) = (U(3,2)-U(2,3)) / (2.0d0*sin(theta))
    axis(2) = (U(1,3)-U(3,1)) / (2.0d0*sin(theta))
    axis(3) = (U(2,1)-U(1,2)) / (2.0d0*sin(theta))
end if
  
end subroutine rotmat_to_axisangle
!-------------------------------------------------------------------------------
subroutine reparam_expmap(p)
!-------------------------------------------------------------------------------
! Reparametrize exponential mapping vector p(3) 
! if the vector gets out of boundary
!-------------------------------------------------------------------------------
real(dp), intent(inout) :: p(3)
real(dp) :: p_sqr
  
p_sqr = dot_product(p(:),p(:))
if (p_sqr > pi*pi) then
    p(:) = (1.0d0 - 2.0d0*pi/sqrt(p_sqr))*p(:)
end if

end subroutine reparam_expmap
!-------------------------------------------------------------------------------
subroutine quat_to_yaw(q, p)
!-------------------------------------------------------------------------------
! Do quaternion, yaw-pitch-roll : q(4) -> p(3)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: q(4)
real(dp), intent(out) :: p(3)
real(dp) :: w, x, y, z, Nq, s
real(dp) :: r11, r21, r31, r32, r33
  
w = q(1)
x = q(2)
y = q(3)
z = q(4)

Nq = w*w + x*x + y*y + z*z
s = 0.0d0
if (Nq > small_real) s = 2.0d0/Nq

r11 = 1.0d0 - s*(y*y+z*z)
r21 = s*(x*y + z*w)
r31 = s*(x*z - y*w)
r32 = s*(y*z + x*w)
r33 = 1.0d0 - s*(x*x+y*y)

p(1) = atan2( r21, r11 )
p(2) = atan2( -r31, sqrt((r32**2 + r33**2)) )
p(3) = atan2( r32, r33 )

end subroutine quat_to_yaw
!-------------------------------------------------------------------------------
subroutine yaw_to_rotmat(angles, U)
!-------------------------------------------------------------------------------
! Calculates rotation angle and axis from the rotational matrix
!-------------------------------------------------------------------------------
real(dp), intent(in) :: angles(3)
real(dp), intent(out):: U(:,:)
real(dp) :: sina, cosa, sinb, cosb, sinc, cosc

sina = sin(angles(1))
sinb = sin(angles(2))
sinc = sin(angles(3))
cosa = cos(angles(1))
cosb = cos(angles(2))
cosc = cos(angles(3))
  
U(1,1) = cosa*cosb
U(1,2) = cosa*sinb*sinc - sina*cosc
U(1,3) = cosa*sinb*cosc + sina*sinc
U(2,1) = sina*cosb
U(2,2) = sina*sinb*sinc + cosa*cosc
U(2,3) = sina*sinb*cosc - cosa*sinc
U(3,1) = -sinb
U(3,2) = cosb*sinc
U(3,3) = cosb*cosc

end subroutine yaw_to_rotmat
!-------------------------------------------------------------------------------
!===============================================================================
! Etc.
!===============================================================================
!-------------------------------------------------------------------------------
subroutine get_principal_axis(n_atm, atm_mass, crd, center, axis)
!-------------------------------------------------------------------------------
! Calculate principal axis of coordinates by solving eigenvalue problem
!-------------------------------------------------------------------------------
integer, intent(in) :: n_atm
real(dp), intent(in) :: atm_mass(n_atm), crd(3,n_atm)
real(dp), intent(in) :: center(3)
real(dp), intent(out) :: axis(3,3)
real(dp) :: mass_corr(3,3), x_sqr, m
real(dp) :: crd_s(3,n_atm), e_value(3), e_vector(3,3)
integer :: i, j, k, nrot
real(dp) :: r_N(3), r_C(3), r0(3), dr(3)

do k = 1, n_atm
    crd_s(:,k) = crd(:,k) - center(:)
end do

mass_corr(:,:) = 0.0
do k = 1, n_atm
    x_sqr = dot_product(crd_s(:,k), crd_s(:,k))
    m = atm_mass(k)
    do i = 1, 3
        do j = 1, 3
            if (i == j) then
                mass_corr(i,i) = mass_corr(i,i) + m*(x_sqr - crd_s(i,k)**2)
            else
                mass_corr(i,j) = mass_corr(i,j) - m*crd_s(i,k)*crd_s(j,k)
            end if
        end do
    end do
end do

call jacobi(3, mass_corr, e_value, e_vector, nrot)
call eigsrt(3, e_value, e_vector) ! sorting eigenvectors
axis = dble(e_vector)

! make sure about the directions of the axes
! z axis points from N to C terminal
r_N(:) = crd(:,1)
r_C(:) = crd(:,n_atm)
dr(:) = r_C(:) - r_N(:)
if (dot_product(axis(:,3), dr) < 0.0d0) then
    ! flip x and z axes
    axis(:,1) = -axis(:,1)
    axis(:,3) = -axis(:,3)
end if

! x axis points from center of helix to center of molecule
r0(:) = crd(:, int(n_atm/2))
dr(:) = r0(:) - center(:)
if (dot_product(axis(:,1), dr) < 0.0d0) then
    ! flip y and z axes
    axis(:,1) = -axis(:,1)
    axis(:,2) = -axis(:,2)
end if

end subroutine get_principal_axis
!-------------------------------------------------------------------------------
subroutine get_center_CA(crd, center, start_res, end_res)
!-------------------------------------------------------------------------------
! calc center of mass of C-alpha atoms in  input coordinates.
! since input crd data has no information of atom type (input type: R),
! this subroutine needs residue index of input coordinates.
!
! you should need start and end residue number. 
! (res_index_type don't have total residue number data.)
!-------------------------------------------------------------------------------
  real(dp), intent(in):: crd(:,:)
  real(dp), intent(out):: center(3)
  integer, intent(in) :: start_res, end_res
  integer :: i_res, i_atm
  integer :: tot_atm
  center(:) = 0.0d0
  tot_atm = 0

  do i_res = start_res, end_res
     i_atm = res_index(i_res)%Ca_id(1)
     center(:) = center(:) + crd(:,i_atm)
     tot_atm = tot_atm +1
  enddo
  center(:) = center(:) / dble(tot_atm)

end subroutine get_center_CA
!-------------------------------------------------------------------------------
subroutine get_center_of_mass(n_atm, atm_mass, crd, center)
!-------------------------------------------------------------------------------
! calc center of mass of input coordinates.
!-------------------------------------------------------------------------------
integer, intent(in) :: n_atm
real(dp), intent(in) :: atm_mass(n_atm), crd(3,n_atm)
real(dp), intent(out) :: center(3)
real(dp) :: sum_mass, sum_r(3)
integer :: i

sum_mass = 0.0d0
sum_r(:) = 0.0d0

do i = 1, n_atm
    sum_mass = sum_mass + atm_mass(i)
    sum_r(:) = sum_r(:) + atm_mass(i)*crd(:,i)
end do

center(:) = sum_r(:)/sum_mass

end subroutine get_center_of_mass
!-------------------------------------------------------------------------------
subroutine radius_of_gyration(crd, mass, n, rg)
!-------------------------------------------------------------------------------
! calc radius of gyration of input coordinates.
!-------------------------------------------------------------------------------
integer, intent(in) :: n
real(dp), intent(in) :: crd(3,n), mass(n)
real(dp), intent(out) :: rg
integer :: i
real(dp) :: com(3), Rc(3)

rg = 0.0d0
com(:) = 0.0d0

do i = 1, n
    com(:) = com(:) + mass(i)*crd(:,i)
end do
com(:) = com(:)/sum(mass(:))

do i = 1, n
    Rc(:) = crd(:,i) - com(:)
    rg = rg + mass(i)*dot_product(Rc(:),Rc(:))
end do
  
rg = sqrt(rg/sum(mass(:)))

end subroutine radius_of_gyration
!-------------------------------------------------------------------------------
subroutine identify_t_ang_dependence(ref_res)
!-------------------------------------------------------------------------------
! Identifies torsion angles which rotates together with a selected t_ang
!-------------------------------------------------------------------------------
type(ref_res_type), intent(inout) :: ref_res
integer :: atm_in_t_ang(4), i_ang, i_atm, n_ang
real(dp) :: t_ang0

ref_res%n_ang_dep_on_t_ang(:) = 0
ref_res%ang_dep_on_t_ang(:,:) = 0
ref_res%d_ang_by_t_ang(:,:) = 0.0d0
!print*, ref_res%res_name
if (ref_res%res_name == 'DISU' .or. ref_res%res_name == 'HOH') return

do i_atm = 1, ref_res%n_atm
    atm_in_t_ang(:) = ref_res%atm_in_t_ang(:,i_atm)
    !
    n_ang = 0
    t_ang0 = ref_res%t_ang0(i_atm)
    do i_ang = 1, ref_res%n_atm
        if (atm_in_t_ang(1) == ref_res%atm_in_t_ang(1, i_ang) .and. & 
            atm_in_t_ang(2) == ref_res%atm_in_t_ang(2, i_ang) .and. &
            atm_in_t_ang(3) == ref_res%atm_in_t_ang(3, i_ang)) then

            n_ang = n_ang + 1
            !
!            write(*,'(A,I4,2(2x,I4,4(1x,A4)))') 'DEBUG: ',&
!                n_ang, i_atm, ref_res%atom_name(atm_in_t_ang),&
!                i_ang, ref_res%atom_name(ref_res%atm_in_t_ang(:,i_ang))
            !
            ref_res%ang_dep_on_t_ang(n_ang,i_atm) = i_ang
            ref_res%d_ang_by_t_ang(n_ang, i_atm) = ref_res%t_ang0(i_ang)-t_ang0
        end if
    end do
    !
    if (ref_res%atom_name(atm_in_t_ang(4)) == 'N   ') then
        ! to add torsion defined by -N -CA -C -O
        n_ang = n_ang + 1
        ref_res%ang_dep_on_t_ang(n_ang, i_atm) = -1 ! arbitary number for -O
        ref_res%d_ang_by_t_ang(n_ang, i_atm) = pi
    end if
    !
    ref_res%n_ang_dep_on_t_ang(i_atm) = n_ang
end do

end subroutine identify_t_ang_dependence
!-------------------------------------------------------------------------------
subroutine rotate_torsion_angle(residue, residue_prev, i_ang, t_ang, preserve_dep_in)
!-------------------------------------------------------------------------------
! Rotate torsion angles for a given index i_ang and related t_ang_s
!-------------------------------------------------------------------------------
type(residue_type), intent(inout) :: residue, residue_prev
integer, intent(in) :: i_ang
real(dp), intent(in) :: t_ang
logical, intent(in), optional :: preserve_dep_in
!
type(ref_res_type) :: ref, ref_prev
integer :: i_tor, j_ang, k_ang
real(dp) :: t_ang0, d_ang
logical :: preserve_dep

if (present(preserve_dep_in) .and. preserve_dep_in) then
    preserve_dep = .true.
else
    preserve_dep = .false.
end if

ref = ref_res(residue%res_type)
ref_prev = ref_res(residue_prev%res_type)

do i_tor = 1, ref%n_ang_dep_on_t_ang(i_ang)
    j_ang = ref%ang_dep_on_t_ang(i_tor, i_ang)
    if (preserve_dep) then
        t_ang0 = residue%t_ang(i_ang)
        if (j_ang /= -1) then
            d_ang = residue%t_ang(j_ang) - t_ang0
            residue%t_ang(j_ang) = t_ang + d_ang
        else
            do k_ang = 1, ref_prev%n_atm
                if (ref_prev%atom_name(k_ang) == 'O   ') then
                    j_ang = k_ang
                    exit
                end if
            end do
            d_ang = residue_prev%t_ang(j_ang) - t_ang0
            residue_prev%t_ang(j_ang) = t_ang + d_ang
        end if
    else
        d_ang = ref%d_ang_by_t_ang(i_tor, i_ang)
        !
        if (j_ang /= -1) then
            residue%t_ang(j_ang) = t_ang + d_ang
        else
            do k_ang = 1, ref_prev%n_atm
                if (ref_prev%atom_name(k_ang) == 'O   ') then
                    j_ang = k_ang
                    exit
                end if
            end do
            residue_prev%t_ang(j_ang) = t_ang + d_ang
        end if
    end if
end do

end subroutine rotate_torsion_angle
!-------------------------------------------------------------------------------
END MODULE GEOMETRY
!-------------------------------------------------------------------------------
