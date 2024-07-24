!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/core/utils/ran.f90
!
! Description:
!  This module contains random number initializer & generator.
!  We recommend to use this module instead of Fortran built-in function, since
!  the Fortran built-in function is recognized to invoke some problems.
!-------------------------------------------------------------------------------
MODULE RAN
!-------------------------------------------------------------------------------
use globals
use logger, only: log_p, terminate_with_error

implicit none
private

integer :: seed                        ! Random seed
logical :: seed_given                  ! Whether random seed defined

public :: initialize_random
public :: random
public :: ctrper
public :: seed, seed_given

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_random()
!-------------------------------------------------------------------------------
! Initialize random number (with random seed if given, otherwise as a function of time)
!-------------------------------------------------------------------------------
real(dp) :: sec
integer :: ij, kl, days
character(len=8) :: date
character(len=10) :: time
character(len=5) :: zone
integer, dimension(8) :: values(8), days_in_month(12)

! set random number seed
if (.not. seed_given) then
    days_in_month = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)

    call date_and_time(date, time, zone, values)
   
    if (mod(values(1),4) == 0) days_in_month(2) = 29
   
    days = sum(days_in_month(1:(values(2)-1)))+values(3)
    sec = 1.0d-3 * (    values(8) + 1.0d3*(values(7)+ &
                    60*(values(6) +    60*(values(5) + 24*days))))
    seed = int(sec)+me*10
else 
    call log_p("- Starting with a specified random seed.", me=me, level=20)
    seed = seed+me*10
end if

write(log_msg,"(A,I10)") '  Random seed = ', seed
call log_p(log_msg, me=me, level=20)

ij = mod(seed, 31329)
kl = int(30082*ran3(seed))
call rmarin(ij, kl)

end subroutine initialize_random
! ----------------------------------------------------------------------
function random()
! ----------------------------------------------------------------------
! This is the random number generator proposed by George Marsaglia in
! Florida state university report: fsu-scri-87-50
! It was slightly modified by F. James to produce an array of pseudorandom
! numbers.
! Modified by Chaok 
!-------------------------------------------------------------------------------
real(dp) :: random
real(dp) :: u(97), c, cd, cm
integer :: i97, j97
logical :: test
common /raset1/ u, c, cd, cm, i97, j97, test
real(dp) :: uni

uni = u(i97) - u(j97)
if (uni < 0.0d0) uni = uni + 1.0d0

u(i97) = uni
i97 = i97 - 1
if(i97 == 0) i97 = 97

j97 = j97 - 1
if (j97 == 0) j97 = 97

c = c - cd
if (c < 0.0d0) c = c + cm

uni = uni - c
if (uni < 0.0d0) uni = uni + 1.0d0
random = uni

end function random
!-------------------------------------------------------------------------------
function ran3(idum)
! ------------------------------------------------------------------
! This is a portable random number generator, and
! "rmarin" is to initialize the random number generator "random".
! Two seed numbers should be  0 < iij < 31328,  0 < kkj < 30081.
!-------------------------------------------------------------------------------
real(dp) :: ran3
integer, intent(inout) :: idum
integer, parameter :: mbig = 1000000000, mseed = 161803398, mz = 0
real(dp), parameter :: fac = 1.0d0/mbig
integer :: i, iff, ii, inext, inextp, k
integer :: mj, mk, ma(55)

save iff,inext,inextp,ma
data iff /0/

if (idum < 0 .or. iff == 0) then
    iff = 1
    mj = mseed - iabs(idum)
    mj = mod(mj,mbig)
    ma(55) = mj
    mk = 1
    do i = 1, 54
        ii = mod(21*i,55)
        ma(ii) = mk
        mk = mj - mk
        if (mk <= mz) mk = mk + mbig
        mj = ma(ii)
    end do
   
    do k = 1, 4
        do i = 1, 55
            ma(i) = ma(i) - ma(1+mod(i+30,55))
            if (ma(i) < MZ ) ma(i) = ma(i) + mbig
        end do
    end do
   
    inext = 0
    inextp = 31
    idum = 1
end if

inext = inext + 1
if (inext == 56) inext = 1

inextp = inextp + 1
if (inextp == 56) inextp = 1

mj = ma(inext) - ma(inextp)
if (mj <= mz) mj = mj + mbig

ma(inext) = mj
ran3 = mj*fac

end function ran3
! ----------------------------------------------------------------------
subroutine rmarin(ij,kl)
! ----------------------------------------------------------------------
! This is the initialization routine for the random number generator random()
! note: the seed variables can have values between:    0 <= ij <= 31328
!                                                      0 <= kl <= 30081
! the random number sequences created by these two seeds are of sufficient
! length to complete an entire calculation with. for example, if sveral
! different groups are working on different parts of the same calculation,
! each group could be assigned its own ij seed. this would leave each group
! with 30000 choices for the second seed. that is to say, this random
! number generator can create 900 million different subsequences -- with
! each subsequence having a length of approximately 10^30.
!
! use ij = 1802 & kl = 9373 to test the random number generator. the
! subroutine random should be used to generate 20000 random numbers.
! then display the next six random numbers generated multiplied by 4096*4096
! if the random number generator is working properly, the random numbers
! should be:
!           6533892.0  14220222.0  7275067.0
!           6172232.0  8354498.0   10633180.0
! ----------------------------------------------------------------------
integer, intent(in) :: ij, kl
real(dp) :: u(97), c, cd, cm
integer :: i97, j97
logical :: test
common /raset1/ u, c, cd, cm, i97, j97, test
integer :: i, j, k, l, ii, jj, m
real(dp) :: s, t

test = .false.

if (ij < 0 .or. ij > 31328 .or. kl < 0 .or. kl > 30081) then
    call log_p(' the first random number seed must have a value between 0 and 31328')
    call terminate_with_error(' the second seed must have a value between 0 and 30081')
end if

i = mod(ij/177, 177) + 2
j = mod(ij    , 177) + 2
k = mod(kl/169, 178) + 1
l = mod(kl,     169)
do ii = 1, 97
    s = 0.0d0
    t = 0.5d0
    do jj = 1, 24
        m = mod(mod(i*j, 179)*k, 179)
        i = j
        j = k
        k = m
        l = mod(53*l+1, 169)
        if (mod(l*m, 64) >= 32) then
            s = s + t
        end if
        t = 0.5d0*t
    end do
    u(ii) = s
end do

c = 362436.0d0 / 16777216.0d0
cd = 7654321.0d0 / 16777216.0d0
cm = 16777213.0d0 /16777216.0d0
i97 = 97
j97 = 33
test = .true.

end subroutine rmarin
!-------------------------------------------------------------------------------
subroutine ctrper(n, xdont, pcls)
!-------------------------------------------------------------------------------
! Permute array xdont randomly, 
! but leaving elements close to their initial locations 
! where nearbyness is controled by PCLS.
!
! The relative proportion of initial order and random order
! is 1-PCLS / PCLS, thus when PCLS = 0, there is no change in
! the order whereas the new order is fully random when PCLS = 1.
! Michel Olagnon - May 2000.
!-------------------------------------------------------------------------------
integer, intent(in) :: n
integer, intent(inout)  :: xdont(n)
real(dp), intent(in) :: pcls
real(dp) :: pwrk, xindt(n)
integer :: i, jwrkt(n)

do i = 1, n
    xindt(i) = random()
end do

pwrk = min(max(0.0d0, pcls), 1.0d0)
xindt(:) = dble(n)*xindt(:)
xindt(:) = pwrk*xindt + (1.0d0 - pwrk)*(/ (dble(i), i=1,n) /)

call mrgrnk(n, xindt, jwrkt)
xdont = xdont(jwrkt)

end subroutine ctrper
!-------------------------------------------------------------------------------
subroutine mrgrnk(n, xdont, irngt)
!-------------------------------------------------------------------------------
! Merge-sort ranking of an array
! For performance reasons, the first 2 passes are taken
! out of the standard loop, and use dedicated coding.
!-------------------------------------------------------------------------------
integer, intent(in) :: n
real(dp), intent(in)  :: xdont(n)
integer, intent(out) :: irngt(n)
real(dp) :: xvala, xvalb
integer :: jwrkt(n)
integer :: lmtna, lmtnc, irng1, irng2
integer :: iind, iwrkd, iwrk, iwrkf, jinda, iinda, iindb

!  Fill-in the index array, creating ordered couples
do iind = 2, n, 2
    if (xdont(iind-1) <= xdont(iind)) then
        irngt(iind-1) = iind - 1
        irngt(iind) = iind
    else
        irngt(iind-1) = iind
        irngt(iind) = iind - 1
    end if
end do
if (modulo(n, 2) /= 0) then
    irngt(n) = n
end if

! We will now have ordered subsets A - B - A - B - ...
! and merge A and B couples into     C   -   C   - ...
lmtna = 2
lmtnc = 4

! First iteration. The length of the ordered subsets goes from 2 to 4
do
    if (n <= 2) exit
    !  Loop on merges of A and B into C
    do iwrkd = 0, n-1, 4
        if ((iwrkd+4) > n) then
            if ((iwrkd+2) >= n) exit
            !   1 2 3
            if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) exit
            !   1 3 2
            if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
                irng2 = irngt(iwrkd+2)
                irngt(iwrkd+2) = irngt(iwrkd+3)
                irngt(iwrkd+3) = irng2
                !   3 1 2
            else
                irng1 = irngt(iwrkd+1)
                irngt(iwrkd+1) = irngt(iwrkd+3)
                irngt(iwrkd+3) = irngt(iwrkd+2)
                irngt(iwrkd+2) = irng1
            end if
            exit
        end if
        !   1 2 3 4
        if (xdont(irngt(iwrkd+2)) <= xdont(irngt(iwrkd+3))) cycle
        !   1 3 x x
        if (xdont(irngt(iwrkd+1)) <= xdont(irngt(iwrkd+3))) then
            irng2 = irngt(iwrkd+2)
            irngt(iwrkd+2) = irngt(iwrkd+3)
            if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then !   1 3 2 4
                irngt(iwrkd+3) = irng2
            else !   1 3 4 2
                irngt(iwrkd+3) = irngt(iwrkd+4)
                irngt(iwrkd+4) = irng2
            end if
        else !   3 x x x
            irng1 = irngt(iwrkd+1)
            irng2 = irngt(iwrkd+2)
            irngt(iwrkd+1) = irngt(iwrkd+3)
            if (xdont(irng1) <= xdont(irngt(iwrkd+4))) then
                irngt (iwrkd+2) = irng1
                if (xdont(irng2) <= xdont(irngt(iwrkd+4))) then !   3 1 2 4
                    irngt(iwrkd+3) = irng2
                else  !   3 1 4 2
                    irngt(iwrkd+3) = irngt(iwrkd+4)
                    irngt(iwrkd+4) = irng2
                end if
            else
                !   3 4 1 2
                irngt(iwrkd+2) = irngt(iwrkd+4)
                irngt(iwrkd+3) = irng1
                irngt(iwrkd+4) = irng2
            end if
        end if
    end do
    ! The Cs become As and Bs
    lmtna = 4
    exit
end do

! Iteration loop. Each time, the length of the ordered subsets
! is doubled.
do
    if (lmtna >= n) exit
    iwrkf = 0
    lmtnc = 2 * lmtnc
    
    ! Loop on merges of A and B into C
    do
        iwrk = iwrkf
        iwrkd = iwrkf + 1
        jinda = iwrkf + lmtna
        iwrkf = iwrkf + lmtnc
        if (iwrkf >= n) then
            if (jinda >= n) exit
            iwrkf = n
        end if
        iinda = 1
        iindb = jinda + 1
        ! Shortcut for the case when the max of A is smaller
        ! than the min of B. This line may be activated when the
        ! initial set is already close to sorted.
        
        ! IF (xdont(irngt(jinda)) <= xdont(irngt(iindB))) CYCLE
        ! One steps in the C subset, that we build in the final rank array
        ! Make a copy of the rank array for the merge iteration
        jwrkt(1:lmtna) = irngt(iwrkd:jinda)
        xvala = xdont(jwrkt(iinda))
        xvalb = xdont(irngt(iindb))
        
        do
            iwrk = iwrk + 1
            ! We still have unprocessed values in both A and B
            if (xvala > xvalb) then
                irngt(iwrk) = irngt(iindb)
                iindb = iindb + 1
                if (iindb > iwrkf) then
                    ! Only A still with unprocessed values
                    irngt(iwrk+1:iwrkf) = jwrkt(iinda:lmtna)
                    exit
                end if
                xvalb = xdont(irngt(iindb))
            else
                irngt(iwrk) = jwrkt(iinda)
                iinda = iinda + 1
                if (iinda > lmtna) exit ! Only B still with unprocessed values
                xvala = xdont(jwrkt(iinda))
            end if
        end do
    end do
    !  The Cs become As and Bs
    lmtna = 2 * lmtna
end do
  
end subroutine mrgrnk
!-------------------------------------------------------------------------------
END MODULE RAN
