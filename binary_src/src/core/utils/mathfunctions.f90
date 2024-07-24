!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/core/utils/mathfunctions.f90
!
! Description:
!  This module contains various mathematical functions and methods
!   - Basic functions (cross, normalize, ...)
!   - Statistical functions
!   - Interpolations
!   - Smooth connections
!   - Various Gaussian-function types
!   - Polynomials
!   - Angle related functions
!   - and eigenvalue problem solver.
!
!-------------------------------------------------------------------------------
MODULE MATHFUNCTIONS
!-------------------------------------------------------------------------------
use globals, only: dp, pi, two_pi, top_type

implicit none
public

real(dp), parameter :: i_twopi_2 = 1.0d0/sqrt(two_pi)

! Matrix for bicubic interpolation
integer, private, parameter :: wt2(16*16) = &
    (/ (/ 1, 0,-3, 2, 0, 0, 0, 0,-3, 0, 9,-6, 2, 0,-6, 4 /), &
       (/ 0, 0, 0, 0, 0, 0, 0, 0, 3, 0,-9, 6,-2, 0, 6,-4 /), &
       (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-6, 0, 0,-6, 4 /), &
       (/ 0, 0, 3,-2, 0, 0, 0, 0, 0, 0,-9, 6, 0, 0, 6,-4 /), &
       (/ 0, 0, 0, 0, 1, 0,-3, 2,-2, 0, 6,-4, 1, 0,-3, 2 /), &
       (/ 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 3,-2, 1, 0,-3, 2 /), &
       (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 2, 0, 0, 3,-2 /), &
       (/ 0, 0, 0, 0, 0, 0, 3,-2, 0, 0,-6, 4, 0, 0, 3,-2 /), &
       (/ 0, 1,-2, 1, 0, 0, 0, 0, 0,-3, 6,-3, 0, 2,-4, 2 /), &
       (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 3,-6, 3, 0,-2, 4,-2 /), &
       (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 2,-2 /), &
       (/ 0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 3,-3, 0, 0,-2, 2 /), &
       (/ 0, 0, 0, 0, 0, 1,-2, 1, 0,-2, 4,-2, 0, 1,-2, 1 /), &
       (/ 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 2,-1, 0, 1,-2, 1 /), &
       (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0,-1, 1 /), &
       (/ 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 2,-2, 0, 0,-1, 1 /) /)

integer, private :: wt(16, 16) = reshape(wt2, (/16,16/)) ! To compile with
                                                         ! gfortran

CONTAINS
!===============================================================================
! Basic functions
!===============================================================================
!-------------------------------------------------------------------------------
subroutine cross(p, q, s)
!-------------------------------------------------------------------------------
! Crossproduct p(3) and q(3)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: p(3), q(3)
real(dp), intent(out) :: s(3)

s(1) = p(2)*q(3) - p(3)*q(2)
s(2) = p(3)*q(1) - p(1)*q(3)
s(3) = p(1)*q(2) - p(2)*q(1)

end subroutine cross
!-------------------------------------------------------------------------------
subroutine atan2(y, x, z)
!-------------------------------------------------------------------------------
! Arc tangent for x and its perpendicular y
!-------------------------------------------------------------------------------
real(dp), intent(in) :: x, y
real(dp), intent(out) :: z
real(dp) :: w
real(dp) :: sign, phi

if(y >= 0) then
    sign = 1.0d0
else if(y < 0) then
    sign = -1.0d0
end if

w = y / x
w = abs(w)

phi = atan(w)
if (x > 0) then
    z = phi * sign
else if (x == 0) then
    z = pi / 2.0d0 * sign
else if(x < 0) then
    z = (pi - phi) * sign
end if

end subroutine atan2
!-------------------------------------------------------------------------------
subroutine v_norm(vec)
!-------------------------------------------------------------------------------
! Normalize vector
!-------------------------------------------------------------------------------
real(dp), dimension(:), intent(inout) :: vec
real(dp) :: vec_size

vec_size = sqrt(dot_product(vec,vec))
vec = vec/vec_size
  
end subroutine v_norm
!-------------------------------------------------------------------------------
function projected_to_normal_plane(A,B)
!-------------------------------------------------------------------------------
! projecting vector A to the normal plane of vector B
!-------------------------------------------------------------------------------
real(dp) :: projected_to_normal_plane(3)
real(dp), intent(in) :: A(3), B(3)
real(dp) :: v(3)
real(dp) :: A_B, i_d2B

A_B = dot_product(A(:),B(:))
i_d2B = 1.0d0/dot_product(B(:),B(:))
v(:) = A_B*B(:)*i_d2B

projected_to_normal_plane(:) = A(:)-v(:)

end function projected_to_normal_plane
!-------------------------------------------------------------------------------
subroutine unit_matrix(n, U)
!-------------------------------------------------------------------------------
integer, intent(in) :: n
real(dp), intent(out) :: U(n,n)
integer :: i

U(:,:) = 0.0d0
do i = 1, n
    U(i,i) = 1.0d0
end do

end subroutine unit_matrix
!-------------------------------------------------------------------------------
function hamming_distance(n,X,Y,is_ang)
!-------------------------------------------------------------------------------
! calculate delta-angle distance. 
! returns sum norm or one norm of difference vector, X-Y  ( = ||X-Y||_1)
! the name "hamming distance" is not appropriate word.
! ||X-Y||_1 = |X1 - Y1| + |X2 - Y2| +..... +|Xn - Yn|
! n: dimension of vector
! X,Y: vector
! is_ang: if the value is .true. , caculate absolute angle-diff distance < pi.
!         ( 359 degree is same with -1 degree.)
!-------------------------------------------------------------------------------
real(dp) :: hamming_distance
integer, intent(in) :: n
real(dp), intent(in) :: X(n), Y(n)
logical, intent(in) :: is_ang

integer :: i
real(dp) :: dXY(n)

dXY(:) = abs(X(:)-Y(:))
if (is_ang) then
    do i = 1, n
        dXY(i) = bound_ang(dXY(i))!returns value between -pi and pi, 
                                  !by adding 2n * pi (n: integer)
    end do
end if

hamming_distance = sum(dXY)

end function hamming_distance
!===============================================================================
! Statistical functions
!===============================================================================
!-------------------------------------------------------------------------------
subroutine Pearson_correlation(n, dat1, dat2, corr)
!-------------------------------------------------------------------------------
! linear correlation between two variables dat1 and dat2
!-------------------------------------------------------------------------------
integer, intent(in) :: n
real(dp), intent(in) :: dat1(n), dat2(n)
real(dp), intent(out) :: corr
real(dp) :: av1, av2, stdev1, stdev2, arg
integer :: i

av1 = sum(dat1(:))/dble(n)
av2 = sum(dat2(:))/dble(n)

stdev1 = 0.0d0
stdev2 = 0.0d0
do i = 1, n
    stdev1 = stdev1 + (dat1(i) - av1)*(dat1(i) - av1)
    stdev2 = stdev1 + (dat1(i) - av1)*(dat1(i) - av1)
end do

stdev1 = sqrt(stdev1)
stdev2 = sqrt(stdev2)

corr = 0.0d0
do i = 1, n
    corr = corr + (dat1(i) - av1)*(dat2(i) - av2)
end do

arg = stdev1*stdev2
corr = corr/arg

end subroutine Pearson_correlation
!-------------------------------------------------------------------------------
subroutine Pearson_correlation_with_g(n, dat1, dat2, corr, g)
!-------------------------------------------------------------------------------
! calculates pearson correlation (corr) between two dataset, data1 and data2 
! then, calculates gradient of partial derivative of the correlation.
!
! see M. Strickert et al."Derivatives of pearson Correlation for Gradient-based
! Analysis of Biomedical Data", inteligencia Artificial Vol. 12 No 37, 2008
!-------------------------------------------------------------------------------
integer, intent(in) :: n !number of data
real(dp), intent(in) :: dat1(n), dat2(n) !input data
real(dp), intent(out) :: corr, g(n) !corr: pearson correlation
!g(i): partial derivative of corr in the direction of dat2(i)
real(dp) :: av1, av2, stdev1, stdev2, arg
integer :: i
real(dp) :: B, C, D, d1, d2

av1 = sum(dat1(:))/dble(n) ! average of data1
av2 = sum(dat2(:))/dble(n) ! average of data2

B = 0.0d0
C = 0.0d0
D = 0.0d0
!d1(:) = dat1(:) - av1
!d2(:) = dat2(:) - av2

do i = 1, n
    d1 = dat1(i) - av1 
    d2 = dat2(i) - av2
    B = B + d1*d2 ! covariance * # of data 
    C = C + d1*d1 ! variance of data 1 * # of data
    D = D + d2*d2 ! variance of data 2 * # of data
end do

stdev1 = sqrt(C) 
stdev2 = sqrt(D)

arg = stdev1*stdev2
corr = B/arg !pearson correlation

! Gradient
do i = 1, n
    g(i) = (dat1(i) - av1) - B/D*(dat2(i) - av2)
end do

g(:) = g(:)/arg

end subroutine Pearson_correlation_with_g
!===============================================================================
! Interpolations
!===============================================================================
!-------------------------------------------------------------------------------
function cubic_spline(interval, c_j, c_jpp, a, calc_g)
!-------------------------------------------------------------------------------
! Run cubic spline interpolation for given 
! interval: distance between surrounding bins
! c_j(2): function values at surrounding bins
! c_jpp(2): and the bin's second derivatives
! a: relative location of the point
!    (1.0 at 1st bin, 0.0 at 2nd bin)
! Parameters a, b, c_j, c_jpp should be prepared 
!-------------------------------------------------------------------------------
real(dp) :: cubic_spline(2)
real(dp), intent(in) :: interval, c_j(2), c_jpp(2), a
logical, intent(in) :: calc_g
real(dp) :: b, c, d
real(dp) :: a2, b2, i6, g1, g2, g3

b = 1.0d0 - a
i6 = interval/6.0d0
a2 = a*a
b2 = b*b

c = (a2 - 1.0d0)*a
d = (b2 - 1.0d0)*b

cubic_spline(1) = a*c_j(1) + b*c_j(2) + (c*c_jpp(1) + d*c_jpp(2))*interval*i6

if (calc_g) then
    g1 = (c_j(2) - c_j(1)) / interval
    g2 = (3.0d0*a2 - 1.0d0)*c_jpp(1)
    g3 = (3.0d0*b2 - 1.0d0)*c_jpp(2)
    cubic_spline(2) = g1 - (g2 - g3)*i6
end if

end function cubic_spline
!-------------------------------------------------------------------------------
function bicubic_interpolation(x, y, dx, dy, ff)
!-------------------------------------------------------------------------------
! Run bicubic interpolation for given
! x: x location
! y: y location
! dx, dy: length of bin for x and y direction
! ff(1,:) = f(1:4) function values at surrounding bins on 2-D box
! ff(2,:) = fx(1:4) 
! ff(3,:) = fy(1:4)
! ff(4,:) = fxy(1:4)  x, y, xy direction derivatives
! Note: f, fx, fy, fxy should be numbered as
!                4 3
!                1 2
!-------------------------------------------------------------------------------
real(dp) :: bicubic_interpolation(3)
real(dp), intent(in) :: x, y, dx, dy, ff(4,4)
real(dp) :: c(4,4), cl(16), ansy, ansy1, ansy2, xa(16), dxdy, t, u
integer :: i, j, k, l

dxdy = dx*dy
do i = 1, 4
    xa(i) = ff(1,i)
    xa(i+4) = ff(2,i)*dx
    xa(i+8) = ff(3,i)*dy
    xa(i+12) = ff(4,i)*dxdy
end do

cl(:) = 0.0d0
do i = 1, 16
    do k = 1, 16
        cl(i) = cl(i) + dble(wt(i,k))*xa(k)
    end do
end do

l = 0
do i = 1, 4
    do j = 1, 4
        l = l + 1
        c(i,j) = cl(l)
    end do
end do

t = x/dx
u = y/dy

ansy = 0.0d0
ansy1 = 0.0d0
ansy2 = 0.0d0
do i = 4, 1, -1
    ansy =  t*ansy  + ((c(i,4)*u + c(i,3))*u +c(i,2))*u + c(i,1)
    ansy2 = t*ansy2 + (3.0d0*c(i,4)*u + 2.0d0*c(i,3))*u + c(i,2)
    ansy1 = u*ansy1 + (3.0d0*c(4,i)*t + 2.0d0*c(3,i))*t + c(2,i)
end do

bicubic_interpolation(1) = ansy
bicubic_interpolation(2) = ansy1/dx
bicubic_interpolation(3) = ansy2/dy

end function bicubic_interpolation
!-------------------------------------------------------------------------------
subroutine trilinterp(coord, cell_org, grid_value, width, E, g, calc_g)
!-------------------------------------------------------------------------------
! Run Trilinear interpolration                                  4----8
! coord(3) : current coordinate                                /|   /|  
! cell_org : (x0,y0,z0)                                       3-2--7-6   z y
! grid_value : function values at 8 points in grid            |/   |/    |/
! width : width between grid_point                  (x0,y0,z0)1----5      --x
!-------------------------------------------------------------------------------
real(dp), intent(in) :: coord(3), cell_org(3), grid_value(8), width
real(dp), intent(out) :: E
real(dp), intent(out) :: g(3)
logical, intent(in) :: calc_g
real(dp) :: xd, yd, zd
real(dp) :: c00, c01, c10, c11, c0, c1

xd = (coord(1)-cell_org(1))/width
yd = (coord(2)-cell_org(2))/width
zd = (coord(3)-cell_org(3))/width

c00 = grid_value(1)*(1-yd)+grid_value(2)*yd
c01 = grid_value(3)*(1-yd)+grid_value(4)*yd
c10 = grid_value(5)*(1-yd)+grid_value(6)*yd
c11 = grid_value(7)*(1-yd)+grid_value(8)*yd

c0 = c00*(1-zd)+c01*zd
c1 = c10*(1-zd)+c11*zd

E = c0*(1-xd) + c1*xd

if (calc_g) then ! calc gradient using finite diffrence method
    ! x-direction
    g(1) = (c1-c0)/width
    ! z-direction
    c0 = c00*(1-xd) + c10*xd
    c1 = c01*(1-xd) + c11*xd
    g(3) = (c1-c0)/width
    ! y-direction
    c00 = grid_value(1)*(1-xd)+grid_value(5)*xd
    c01 = grid_value(2)*(1-xd)+grid_value(6)*xd
    c10 = grid_value(3)*(1-xd)+grid_value(7)*xd
    c11 = grid_value(4)*(1-xd)+grid_value(8)*xd
    c0 = c00*(1-zd)+c10*zd
    c1 = c01*(1-zd)+c11*zd
    !
    g(2) = (c1-c0)/width
end if

end subroutine trilinterp
!===============================================================================
! Smooth connections
!===============================================================================
!-------------------------------------------------------------------------------
function sigmoidal1(x, xmin, xmax, calc_g)
!-------------------------------------------------------------------------------
! Connect binnary value (0 and 1) 
! x < xmin: 0.0, x > xmax: 1.0, 
! otherwise continuously change using cosine function
!-------------------------------------------------------------------------------
real(dp) :: sigmoidal1(2)
real(dp), intent(in) :: x, xmin, xmax
logical, intent(in) :: calc_g
real(dp) :: arg1, arg2

sigmoidal1(:) = 0.0d0

if (x < xmin) then
    sigmoidal1(1) = 0.0d0
else if (x > xmax) then 
    sigmoidal1(1) = 1.0d0
else
    arg1 = pi/(xmax-xmin)
    arg2 = cos(arg1*(x-xmin))
    sigmoidal1(1) = 0.5d0 - 0.5d0*arg2
    if (calc_g) sigmoidal1(2) = 0.5d0*sqrt(1.0d0 - arg2**2)*pi/(xmax - xmin)
end if

end function sigmoidal1
!-------------------------------------------------------------------------------
function sigmoidal2(x, xmin, xmax, calc_g)
!-------------------------------------------------------------------------------
! Similar to sigmoidal1, but decrease from 1.0 to 0.0
! otherwise continuously change using cosine function
!-------------------------------------------------------------------------------
real(dp) :: sigmoidal2(2)
real(dp), intent(in) :: x, xmin, xmax
logical, intent(in) :: calc_g
real(dp) :: arg1, arg2

sigmoidal2(:) = 0.0d0

if (x < xmin) then
    sigmoidal2(1) = 1.0d0
else if (x > xmax) then 
    sigmoidal2(1) = 0.0d0
else
    arg1 = pi/(xmax-xmin)
    arg2 = cos(arg1*(x-xmin))
    sigmoidal2(1) = 0.5d0 + 0.5d0*arg2
    if (calc_g) sigmoidal2(2) = -0.5d0*sqrt(1.0d0 - arg2**2)*pi/(xmax - xmin)
end if

end function sigmoidal2
!-------------------------------------------------------------------------------
function cosine_well(x, x0, sigma, calc_g)
!-------------------------------------------------------------------------------
! Ranges from -1 to 0
!-------------------------------------------------------------------------------
real(dp) :: cosine_well(2)
real(dp), intent(in) :: x, x0, sigma
logical, intent(in) :: calc_g
real(dp) :: arg

cosine_well(:) = 0.0d0
arg = (x - x0)/sigma

if (abs(arg) >= 1.0d0) then ! out of boundary
    cosine_well(1) = 0.0d0
else
    cosine_well(1) = -0.5d0*(1.0d0 + cos(arg*pi))
    if (calc_g) cosine_well(2) = 0.5d0*pi*sin(arg*pi)
end if

end function cosine_well
!-------------------------------------------------------------------------------
subroutine fade_function(x, min0, fmin, fmax, max0, f, dfdx)
!-------------------------------------------------------------------------------
! Smoothly connect 4 points with function values either 0 or 1
!                             a b__c d
! (min0,fmin,fmax,max0) =>  ___/    \___
! region defined as         1  2  3  4  5
!-------------------------------------------------------------------------------
real(dp), intent(in) :: x, min0, fmin, fmax, max0
real(dp), intent(out) :: f, dfdx
real(dp) :: dfade_min, dfade_max !Slopes

if (fmin == min0 .and. x == fmin) then
    f = 1.0d0
    dfdx = 0.0d0
elseif (fmax == max0 .and. x == fmax) then
    f = 1.0d0
    dfdx = 0.0d0
else
    dfdx = 0.0d0
    dfade_min = 1.0d0/(fmin-min0) 
    dfade_max = 1.0d0/(max0-fmax)
 
    if (x < fmax) then
        if (x >= fmin) then !region 3
            f = 1.0d0
        elseif (x <= min0) then  !region 1
            f = 0.0d0
        else !region 2
            f = (x - min0)*dfade_min
            dfdx = dfade_min
        end if
    else
        if (x >= max0) then !region 5
            f = 0.0d0
        else !region 4
            f = (max0 - x)*dfade_max
            dfdx = -dfade_max
        end if
    end if
end if

end subroutine fade_function
!===============================================================================
! Gaussian-function related terms
!===============================================================================
!-------------------------------------------------------------------------------
function single_gaussian(x, x0, sigma, is_ang, calc_g)
!-------------------------------------------------------------------------------
! Simple Gaussian, log taken (-log(1/(sigma*sqrt(2*pi))) term is ignored)
! Returns function (1) and dfdx (2)
! Angle difference can be smartly calculated if is_ang turned on
!
! x: current position
! x0: mean value
! sigma: sigma
!-------------------------------------------------------------------------------
real(dp) :: single_gaussian(2)
real(dp), intent(in) :: x, x0, sigma
logical, intent(in) :: is_ang, calc_g
real(dp) :: i_sig2, xdiff

i_sig2 = 1.0d0/(sigma*sigma)

if (is_ang) then
    xdiff = bound_ang(x - x0)
else
    xdiff = x - x0
end if

single_gaussian(1) = 0.5d0*xdiff*xdiff*i_sig2

if (calc_g) single_gaussian(2) = xdiff*i_sig2
  
end function single_gaussian
!-------------------------------------------------------------------------------
function multiple_binormal(n, x, x0, w, sigma, rho, get_log, calc_g)
!-------------------------------------------------------------------------------
! Binormal function (Gaussian with two variables) with multiple modes
! Returns function (1) and dfdx (2:3)
! Alternative functional form instead of exponential is used following MODELLER
!
! x(2): current positions
! x0(2,n): mean values for n-th mode
! w(n): weight for n-th mode
! sigma(2,n): sigma for n-th mode
! rho(n): correlation between two position values for n-th mode
!-------------------------------------------------------------------------------
real(dp) :: multiple_binormal(3)
integer, intent(in) :: n
real(dp), intent(in) :: x(2), x0(2,n), w(n)
real(dp), intent(in) :: sigma(2,n), rho(n)
logical, intent(in) :: get_log, calc_g
integer :: i_m
real(dp) :: x1diff, x2diff, isig1, isig2
real(dp) :: sin1, sin2, cos1, cos2, ss1, ss2
real(dp) :: arg, tmp1, tmp2, tmp3
real(dp) :: g1_tmp, g2_tmp

multiple_binormal(:) = 0.0d0

do i_m = 1, n
    x1diff = x(1) - x0(1,i_m)
    x2diff = x(2) - x0(2,i_m)
    isig1 = 1.0d0/sigma(1,i_m)
    isig2 = 1.0d0/sigma(2,i_m)
    
    sin1 = sin(x1diff)
    sin2 = sin(x2diff)
    cos1 = cos(x1diff)
    cos2 = cos(x2diff)
    tmp1 = 1.0d0 - rho(i_m)**2
    ss1 = sin1*isig1
    ss2 = sin2*isig2
    tmp2 = (1.0d0 - cos1)*isig1*isig1 - rho(i_m)*ss1*ss2 + &
           (1.0d0 - cos2)*isig2*isig2
 
    arg = w(i_m)*exp(-tmp2/tmp1)*isig1*isig2 / (two_pi*sqrt(tmp1))
 
    multiple_binormal(1) = multiple_binormal(1) + arg
    if (calc_g) then
        g1_tmp = ss1 - rho(i_m)*cos1*ss2
        g2_tmp = ss2 - rho(i_m)*cos2*ss1
    
        tmp3 = arg/tmp1
        multiple_binormal(2) = multiple_binormal(2) + tmp3*g1_tmp*isig1
        multiple_binormal(3) = multiple_binormal(3) + tmp3*g2_tmp*isig2
    end if
end do

if (get_log) then
    if (calc_g) multiple_binormal(2:3) = multiple_binormal(2:3)/multiple_binormal(1)
    multiple_binormal(1) = -log(multiple_binormal(1))
else
    if (calc_g) multiple_binormal(2:3) = -1.0d0*multiple_binormal(2:3)
end if

end function multiple_binormal
!-------------------------------------------------------------------------------
function multiple_gaussian(n, x, x0, w, sigma, is_ang, get_log, calc_g)
!-------------------------------------------------------------------------------
! Gaussian function with multiple modes
! Returns function (1) and dfdx (2)
! Angle difference can be smartly calculated if is_ang turned on
!
! x: current positions
! x0(n): mean values for n-th mode
! w(n): weight for n-th mode
! sigma(n): sigma for n-th mode
!-------------------------------------------------------------------------------
real(dp) :: multiple_gaussian(2)
integer, intent(in) :: n
real(dp), intent(in) :: x, x0(n), w(n), sigma(n)
logical, intent(in) :: is_ang, get_log, calc_g
real(dp) :: xdiff, arg, i_sig
integer :: i_m

multiple_gaussian(1:2) = 0.0d0

do i_m = 1, n
    if (is_ang) then
        xdiff = bound_ang(x - x0(i_m))
    else
        xdiff = x - x0(i_m)
    end if
    i_sig = 1.0d0/sigma(i_m)
    
    arg = w(i_m)*i_twopi_2*exp(-0.5d0*(i_sig*xdiff)**2)*i_sig
    multiple_gaussian(1) = multiple_gaussian(1) + arg
    if (calc_g) multiple_gaussian(2) = multiple_gaussian(2) - xdiff*arg*i_sig**2
end do

! Only for if functional form is in log style
if (get_log) then
    if (calc_g) multiple_gaussian(2) = -multiple_gaussian(2) / multiple_gaussian(1)
    multiple_gaussian(1) = -log(multiple_gaussian(1))
end if

end function multiple_gaussian
!-------------------------------------------------------------------------------
function multiple_exp(n, x, x0, c_i, c_exp, is_ang, get_log, calc_g)
!-------------------------------------------------------------------------------
! non-standard exponential function with multiple modes
! Returns function (1) and dfdx (2)
! Angle difference can be smartly calculated if is_ang turned on
!
! x: current positions
! x0(n): mean values for n-th mode
! c_i(n): weight for n-th mode
! c_exp: uniform exponent
!-------------------------------------------------------------------------------
real(dp) :: multiple_exp(2)
integer, intent(in) :: n
real(dp), intent(in) :: x, x0(n), c_i(n), c_exp
logical, intent(in) :: is_ang, get_log, calc_g
real(dp) :: xdiff, sum_ci, arg
integer :: i_m

multiple_exp(:) = 0.0d0
sum_ci = 0.0d0

do i_m = 1, n
    if (is_ang) then
        xdiff = bound_ang(x - x0(i_m))
    else
        xdiff = x - x0(i_m)
    end if
 
    arg = c_i(i_m)*exp(c_exp*xdiff**2)
 
    multiple_exp(1) = multiple_exp(1) + arg
    sum_ci = sum_ci + c_i(i_m)
    if (calc_g) multiple_exp(2) = multiple_exp(2) - 2.0d0*c_exp*xdiff*arg
end do

multiple_exp(:) = multiple_exp(:)/sum_ci

! Only for if functional form is in log style
if (get_log) then
    if (calc_g) multiple_exp(2) = multiple_exp(2) / multiple_exp(1)
    multiple_exp(1) = -log(multiple_exp(1))
end if

end function multiple_exp
!-------------------------------------------------------------------------------
function multiple_biexp(n, x, x0, c_i, c_exp, is_ang, get_log, calc_g)
!-------------------------------------------------------------------------------
! non-standard biexponential function (two variables) with multiple modes
! Returns function (1) and dfdx (2:3)
! Angle difference can be smartly calculated if is_ang turned on
!
! x(2): current positions
! x0(2,n): mean values for n-th mode
! c_i(n): weight for n-th mode
! c_exp: uniform exponent
!-------------------------------------------------------------------------------
real(dp) :: multiple_biexp(3)
integer, intent(in) :: n
real(dp), intent(in) :: x(2), x0(2,n)
real(dp), intent(in) :: c_i(n), c_exp
logical, intent(in) :: is_ang, get_log, calc_g
real(dp) :: x1diff, x2diff, c_tmp, arg1, arg2, sum_ci
integer :: i_m

sum_ci = 0.0d0
multiple_biexp(:) = 0.0d0

do i_m = 1, n
    c_tmp = c_i(i_m)
    sum_ci = sum_ci + c_tmp
    if (is_ang) then
        x1diff = bound_ang(x(1) - x0(1,i_m))
        x2diff = bound_ang(x(2) - x0(2,i_m))
    else
        x1diff = x(1) - x0(1,i_m)
        x2diff = x(2) - x0(2,i_m)
    end if
    arg1 = c_tmp*exp(c_exp*(x1diff**2 + x2diff**2))
  
    multiple_biexp(1) = multiple_biexp(1) + arg1
  
    if (calc_g) then
        arg2 = 2.0d0*c_exp*arg1
        multiple_biexp(2) = multiple_biexp(2) - x1diff*arg2
        multiple_biexp(3) = multiple_biexp(3) - x2diff*arg2
    end if
end do

multiple_biexp(:) = multiple_biexp(:)/sum_ci

! Only for if functional form is in log style
if (get_log) then
    if (calc_g) multiple_biexp(2:3) = multiple_biexp(2:3) / multiple_biexp(1)
    multiple_biexp(1) = -log(multiple_biexp(1))
end if

end function multiple_biexp
!-------------------------------------------------------------------------------
function multiple_multinomial(m, n, x, x0, sigma, w, is_ang, get_log, calc_g)
!-------------------------------------------------------------------------------
! Most general form of Gaussian function,
! with n variables and m modes
! Returns function (1) and dfdx (2:1+n)
! Angle difference can be smartly calculated if is_ang turned on
!
! x(n): current positions
! x0(:,m): mean values for m-th mode
! w(m): weight for m-th mode
! sigma(:,m): sigma for m-th mode
!-------------------------------------------------------------------------------
real(dp) :: multiple_multinomial(1+n)
integer,  intent(in) :: m, n
real(dp), intent(in) :: x(:), x0(:,:), sigma(:,:), w(:)
logical,  intent(in) :: is_ang, get_log, calc_g
integer :: i_m, i_n
real(dp) :: dfdx(n), dfdx_i(n), xdiff(n)
real(dp) :: arg, f_sum, sigmaprod, i_sig2, normfac, tmp
logical :: is_close

f_sum = 0.0d0
dfdx(:) = 0.0d0
normfac = i_twopi_2**n

do i_m = 1, m
    !is_close = .false.
    xdiff(:) = 0.0d0
 
    do i_n = 1, n
        arg = x(i_n) - x0(i_n,i_m)
        if (is_ang) then
            xdiff(i_n) = bound_ang(arg)
        else
            xdiff(i_n) = arg
        end if
  
        ! Below is disabled because of minimization fluctuation (PHB 2011/12/21)
        !if (arg < 5.0d0*sigma(i_n,i_m)) is_close = .true.
    end do
 
    ! Consider modes only if any distance among n variables is at least < 3.0sigma
    !if (.not. is_close) cycle
 
    sigmaprod = 1.0d0
    dfdx_i(:) = 0.0d0
    arg = 0.0d0
 
    do i_n = 1, n
        i_sig2 = 1.0d0/(sigma(i_n,i_m)*sigma(i_n,i_m))
        if (xdiff(i_n) > 5.0d0*sigma(i_n,i_m)) then
            tmp = 10.0d0/sigma(i_n,i_m)
            arg = arg + tmp*xdiff(i_n) - 25.0d0
            if (calc_g) dfdx_i(i_n) = 0.5d0*tmp
        else if (xdiff(i_n) < -5.0d0*sigma(i_n,i_m)) then
            tmp = -10.0d0/sigma(i_n,i_m)
            arg = arg + tmp*xdiff(i_n) - 25.0d0
            if (calc_g) dfdx_i(i_n) = 0.5d0*tmp
        else
            arg = arg + xdiff(i_n)*xdiff(i_n)*i_sig2
            if (calc_g) dfdx_i(i_n) = xdiff(i_n)*i_sig2
        end if
        sigmaprod = sigmaprod*sigma(i_n,i_m)
    end do
 
    arg = w(i_m)*normfac/sigmaprod*exp(-0.5d0*arg)
    f_sum = f_sum + arg
    dfdx(:) = dfdx(:) + arg*dfdx_i(:)
end do

multiple_multinomial(1) = f_sum
multiple_multinomial(2:1+n) = dfdx(1:n)

! Only for if functional form is in log style
if (get_log) then
    multiple_multinomial(1) = -log(f_sum)
    multiple_multinomial(2:1+n) = dfdx(:)/f_sum
end if

end function multiple_multinomial
!-------------------------------------------------------------------------------
function lorentzian(x, x0, gamma, calc_g)
!-------------------------------------------------------------------------------
! 1/gamma - simple lorentzian (Cauchy distribution neglecting 1/pi)
! (1/gamma)(1-(gamma^2/(gamma^2+(x-x0)^2))) 
! The maximum value of this function is 1/gamma. The width is also controlled by gamma
!-------------------------------------------------------------------------------
real(dp) :: lorentzian(2)
real(dp), intent(in) :: x, x0, gamma
logical, intent(in) :: calc_g
real(dp) :: xdiff, tmp
real(dp) :: xdiff_2, gamma_2
real(dp) :: denom, nom

xdiff = x - x0
xdiff_2 = xdiff*xdiff
gamma_2 = gamma**2.0d0
tmp = gamma_2 + xdiff_2
lorentzian(1) = 1.0d0/gamma - gamma/tmp
!
if (calc_g) then
    denom = 2.0d0*gamma*xdiff
    nom = gamma_2*gamma_2 + 2.0d0*gamma_2*xdiff_2 + xdiff_2*xdiff_2
    lorentzian(2) = denom/nom
end if

end function lorentzian
!===============================================================================
! Polynomials
!===============================================================================
!-------------------------------------------------------------------------------
subroutine poly3(x, coeff, min_x, max_x, c, dcdx)
!-------------------------------------------------------------------------------
! Return 3rd order polynomial c and derivative dcdx
! Valid only if min_x < x < max_x
!-------------------------------------------------------------------------------
real(dp), intent(in) :: x, min_x, max_x, coeff(3)
real(dp), intent(out) :: c, dcdx
real(dp) :: x_v(3)

if (x >= max_x .or. x <= min_x) then
    c = 0.0d0
    dcdx = 0.0d0
else
    x_v(1:3) = (/ x*x, x, 1.0d0 /)
   
    c = dot_product(x_v(:), coeff(:))
    dcdx = 2.0d0*coeff(1)*x_v(2) + coeff(2)
end if

end subroutine poly3
!-------------------------------------------------------------------------------
subroutine poly5(x, coeff, min_x, max_x, c, dcdx)
!-------------------------------------------------------------------------------
! Return 5th order polynomial c and derivative dcdx
! Valid only if min_x < x < max_x
!-------------------------------------------------------------------------------
real(dp), intent(in) :: x, min_x, max_x, coeff(5)
real(dp), intent(out) :: c, dcdx
real(dp) :: x_v(5), x2

if (x >= max_x .or. x <= min_x) then
    c = 0.0d0
    dcdx = 0.0d0
else
    x2 = x*x
    x_v(1:5) = (/ x2*x2, x*x2, x2, x, 1.0d0 /)
   
    c = dot_product(x_v(:), coeff(:))
    dcdx = 4.0d0*coeff(1)*x_v(2) + 3.0d0*coeff(2)*x_v(3) &
         + 2.0d0*coeff(3)*x_v(4) + coeff(4)
end if

end subroutine poly5
!-------------------------------------------------------------------------------
subroutine poly7(x, coeff, min_x, max_x, c, dcdx)
!-------------------------------------------------------------------------------
! Return 7th order polynomial c and derivative dcdx
! Valid only if min_x < x < max_x
!-------------------------------------------------------------------------------
real(dp), intent(in) :: x, min_x, max_x, coeff(7)
real(dp), intent(out) :: c, dcdx
real(dp) :: x_v(7), x2, x4

if (x >= max_x .or. x <= min_x) then
    c = 0.0d0
    dcdx = 0.0d0
else
    x2 = x*x
    x4 = x2*x2
    x_v(1:7) = (/ x2*x4, x*x4, x4, x*x2, x2, x, 1.0d0 /)
    
    c = dot_product(x_v(:), coeff(:))
    dcdx = 6.0d0*coeff(1)*x_v(2) + 5.0d0*coeff(2)*x_v(3) &
         + 4.0d0*coeff(3)*x_v(4) + 3.0d0*coeff(4)*x_v(5) &
         + 2.0d0*coeff(5)*x_v(6) + coeff(6)
end if

end subroutine poly7
!-------------------------------------------------------------------------------
subroutine poly8(x, coeff, min_x, max_x, c, dcdx)
!-------------------------------------------------------------------------------
! Return 7th order polynomial c and derivative dcdx
! Valid only if min_x < x < max_x
!-------------------------------------------------------------------------------
real(dp), intent(in) :: x, min_x, max_x, coeff(8)
real(dp), intent(out) :: c, dcdx
real(dp) :: x_v(8), x2, x4, x6

if (x >= max_x .or. x <= min_x) then
    c = 0.0d0
    dcdx = 0.0d0
else
    x2 = x*x
    x4 = x2*x2
    x6 = x2*x4
    x_v(1:8) = (/ x*x6, x6, x*x4, x4, x*x2, x2, x, 1.0d0 /)
    
    c = dot_product(x_v(:), coeff(:))
    dcdx = 7.0d0*coeff(1)*x_v(2) + 6.0d0*coeff(2)*x_v(3) &
         + 5.0d0*coeff(3)*x_v(4) + 4.0d0*coeff(4)*x_v(5) &
         + 3.0d0*coeff(5)*x_v(6) + 2.0d0*coeff(6)*x_v(7) &
         + coeff(7)
end if

end subroutine poly8
!===============================================================================
! Angles
!===============================================================================
!-------------------------------------------------------------------------------
function bound_ang(ang)
!-------------------------------------------------------------------------------
! Modify angle value to be absolute value < pi
!-------------------------------------------------------------------------------
real(dp) :: bound_ang
real(dp), intent(in) :: ang

if (ang < -pi) then
    bound_ang = ang + two_pi
else if (ang > pi) then
    bound_ang = ang - two_pi
else
    bound_ang = ang
end if

end function bound_ang
!-------------------------------------------------------------------------------
function bound_ang90(ang)
!-------------------------------------------------------------------------------
! Modify angle value to be absolute value < pi/2
!-------------------------------------------------------------------------------
real(dp) :: bound_ang90
real(dp), intent(in) :: ang

bound_ang90 = bound_ang(ang)

if (bound_ang90 > 0.5d0*pi) then
    bound_ang90 = pi - bound_ang90
else if (bound_ang90 < -0.5d0*pi) then
    bound_ang90 = pi + bound_ang90
end if

end function bound_ang90
!-------------------------------------------------------------------------------
function cosine_potential(x, a, b, n, calc_g)
!-------------------------------------------------------------------------------
! TODO: top_type would be better to be removed.
! Cosine type potential used for dihedral potential
!-------------------------------------------------------------------------------
real(dp) :: cosine_potential(2)
integer, intent(in) :: n
real(dp), intent(in) :: x, a, b
logical, intent(in) :: calc_g
real(dp) :: f

f = dble(n)*x + a
!curie CHECK this
if (top_type == 'allh_ch22') then
    cosine_potential(1) = b - b*cos(f)
else
    cosine_potential(1) = abs(b) - b*cos(f)     
end if
if (calc_g) cosine_potential(2) = b*n*sin(f)

end function cosine_potential
!===============================================================================
! Matrix and eigenvalue problems
!===============================================================================
!-------------------------------------------------------------------------------
subroutine R3_inverse_matrix(U)
!-------------------------------------------------------------------------------
! Calc inverse matrix of (3,3) matrix
!-------------------------------------------------------------------------------
real(dp), intent(inout):: U(:,:)
real(dp) :: a,b,c,d,e,f,g,h,k, det_u

a = U(1,1)
b = U(1,2)
c = U(1,3)
d = U(2,1)
e = U(2,2)
f = U(2,3)
g = U(3,1)
h = U(3,2)
k = U(3,3)

det_u = a*(e*k-f*h)+b*(f*g-k*d)+c*(d*h-e*g)

U(1,1) = (e*k-f*h)
U(1,2) = (c*h-b*k)
U(1,3) = (b*f-c*e)
U(2,1) = (f*g-d*k)
U(2,2) = (a*k-c*g)
U(2,3) = (c*d-a*f)
U(3,1) = (d*h-e*g)
U(3,2) = (b*g-a*h)
U(3,3) = (a*e-b*d)

U(:,:) = (1.0d0/det_u)*U(:,:)

end subroutine R3_inverse_matrix
!-------------------------------------------------------------------------------
subroutine jacobi(n, a, d, v, nrot)
!-------------------------------------------------------------------------------
! Solve eigenvalue problem using Jacobi method
!-------------------------------------------------------------------------------
integer, intent(in) :: n
integer, intent(out) :: nrot
real(dp), intent(out) :: d(n)
real(dp), intent(inout) :: a(n,n)
real(dp), intent(out) :: v(n,n)
integer :: i, j, ip, iq
real(dp) :: c, g, h, s, sm, t, tau, theta, tresh
real(dp) :: b(n), z(n)

call unit_matrix(n, v(:,:))
do j = 1, n
    b(j) = a(j,j)
end do
d(:) = b(:)
z(:) = 0.0d0
nrot = 0

do i = 1, 50 ! Iter over max trial
    sm = 0.0d0
    do j = 1, n-1
        sm = sm + sum(abs(a(j:n,j)))
    end do
    if (sm == 0.0d0) return
   
    tresh = merge(0.2d0*sm/n**2, 0.0d0, i < 4)
   
    do ip = 1, n-1
        do iq = ip+1, n
            g = 100.0d0*abs(a(ip,iq))
            if ((i > 4) .and. (abs(d(ip))+g == abs(d(ip))) &
                        .and. (abs(d(iq))+g == abs(d(iq)))) then
                a(ip,iq) = 0.0d0
          
            else if (abs(a(ip,iq)) > tresh) then
                h = d(iq) - d(ip)
               
                if (abs(h) + g == abs(h)) then
                    t = a(ip,iq) / h
                else
                    theta = 0.5d0*h/a(ip,iq)
                    t = 1.0d0/(abs(theta) + sqrt(1.0d0 + theta**2))
                    if (theta < 0.0) t = -t
                end if
                c = 1.0d0/sqrt(1.0d0 + t**2)
                s = t*c
                tau = s/(1.0d0 + c)
                h = t*a(ip,iq)
                z(ip) = z(ip) - h
                z(iq) = z(iq) + h
                d(ip) = d(ip) - h
                d(iq) = d(iq) + h
                a(ip,iq) = 0.0d0
               
                call jrotate(ip-1, a(1:ip-1,ip), a(1:ip-1,iq), s, tau)
                call jrotate(iq-ip-1, a(ip,ip+1:iq-1), a(ip+1:iq-1,iq), s, tau)
                call jrotate(n-iq, a(ip,iq+1:n), a(iq,iq+1:n), s, tau)
                call jrotate(n, v(:,ip), v(:,iq), s, tau)
               
                nrot = nrot + 1
            end if
        end do
    end do
   
    b(:) = b(:) + z(:)
    d(:) = b(:)
    z(:) = 0.0d0
end do

end subroutine jacobi
!-------------------------------------------------------------------------------
subroutine jrotate(n, a1, a2, s, tau)
!-------------------------------------------------------------------------------
! Used only for Jacobi
!-------------------------------------------------------------------------------
integer, intent(in) :: n
real(dp), intent(inout) :: a1(n), a2(n)
real(dp), intent(in) :: s, tau
real(dp) :: wk1(n)

wk1(:) = a1(:)
a1(:) = a1(:) - s*(a2(:)  + a1(:)*tau)
a2(:) = a2(:) + s*(wk1(:) - a2(:)*tau)

end subroutine jrotate
!-------------------------------------------------------------------------------
subroutine eigsrt(n, d, v)
!-------------------------------------------------------------------------------
! Sort eigenvectors according to eigenvalues
! Used only for get_principal_axis
!-------------------------------------------------------------------------------
integer, intent(in) :: n
real(dp), intent(inout) :: d(n)
real(dp), intent(inout) :: v(3,n)
integer :: i, j
real(dp) :: di, dj, vi(3), vj(3)

do i = 1, n-1
    j = maxloc(d(i:n), dim=1)
    if (j /= i) then
        di = d(i)
        dj = d(j)
        vi(:) = v(:,i)
        vj(:) = v(:,j)
      
        d(i) = dj
        d(j) = di
        v(:,i) = vj
        v(:,j) = vi
    end if
end do
  
end subroutine eigsrt
!-------------------------------------------------------------------------------
function erfinv(y)
!-------------------------------------------------------------------------------
! Don't ask me what this is :)
! Inverse error function. Please google it
!-------------------------------------------------------------------------------
real(dp) :: erfinv
real(dp), intent(in) :: y
real(dp), parameter :: sqrtpi = 1.772453850905516027d0
real(dp) :: erfx, x, z
real(dp) :: a(4), b(4), c(4), d(2)

! Coefficients for approximation to erfinv in central range
a(:) = (/  0.886226899d0, -1.645349621d0,  0.914624893d0, -0.140543331d0 /)
b(:) = (/ -2.118377725d0,  1.442710462d0, -0.329097515d0,  0.012229801d0 /)
! Coefficients for approximation to erfinv near endpoints
c(:) = (/ -1.970840454d0, -1.624906493d0,  3.429567803d0,  1.641345311d0 /)
d(:) = (/  3.543889200d0,  1.637067800d0 /)

! Get an initial estimate for the inverse error function
if (abs(y) <= 0.7d0) then
    z = y * y
    x = y * (((a(4)*z + a(3))*z + a(2))*z + a(1)) &
         / ((((b(4)*z + b(3))*z + b(2))*z + b(1))*z + 1.0d0)
else if (y > 0.7d0 .and. y < 1.0d0) then
    z = sqrt(-log((1.0d0-y)/2.0d0))
    x = (((c(4)*z + c(3))*z + c(2))*z + c(1)) / ((d(2)*z + d(1))*z + 1.0d0)
else if (y < -0.7d0 .and. y > -1.0d0) then
    z = sqrt(-log((1.0d0+y)/2.0d0))
    x = -(((c(4)*z + c(3))*z + c(2))*z + c(1)) / ((d(2)*z + d(1))*z + 1.0d0)
end if

! Use two steps of Newton-Raphson correction to increase accuracy
call erfcore(x, erfx)
x = x - (erfx - y) / (2.0d0/sqrtpi * exp(-x**2))
call erfcore(x, erfx)
x = x - (erfx - y) / (2.0d0/sqrtpi * exp(-x**2))
erfinv = x
  
end function erfinv
!-------------------------------------------------------------------------------
subroutine erfcore(arg, val)
!-------------------------------------------------------------------------------
! Don't ask me what this is :)
! error function core?? correction??
!-------------------------------------------------------------------------------
real(dp), intent(in) :: arg
real(dp), intent(out) :: val
real(dp), parameter :: sqrpi = 5.6418958354775628695d-1
real(dp), parameter :: thresh = 0.46875d0
real(dp), parameter :: xsmall = 1.11d-16
real(dp), parameter :: xbig = 26.543d0
integer :: i
real(dp) :: x, y, ysq, del, xnum, xden
real(dp) :: a(5), b(4), c(9), d(8), p(6), q(5)

a(:) = (/ 3.16112374387056560d0,  1.13864154151050156d2,  3.77485237685302021d2, &
          3.20937758913846947d3,  1.85777706184603153d-1 /)
b(:) = (/ 2.36012909523441209d1,  2.44024637934444173d2,  1.28261652607737228d3, &
          2.84423683343917062d3 /)

c(:) = (/ 5.64188496988670089d-1, 8.88314979438837594d0,  6.61191906371416295d1, &
          2.98635138197400131d2,  8.81952221241769090d2,  1.71204761263407058d3, &
          2.05107837782607147d3,  1.23033935479799725d3,  2.15311535474403846d-8 /)
d(:) = (/ 1.57449261107098347d1,  1.17693950891312499d2,  5.37181101862009858d2, &
          1.62138957456669019d3,  3.29079923573345963d3,  4.36261909014324716d3, &
          3.43936767414372164d3,  1.23033935480374942d3 /)

p(:) = (/ 3.05326634961232344d-1, 3.60344899949804439d-1, 1.25781726111229246d-1, &
          1.60837851487422766d-2, 6.58749161529837803d-4, 1.63153871373020978d-2 /)
q(:) = (/ 2.56852019228982242d0,  1.87295284992346047d0,  5.27905102951428412d-1, &
          6.05183413124413191d-2, 2.33520497626869185d-3 /)

! Store the argument and its absolute value
x = arg
y = abs(x)

! Evaluate error function for |x| less than 0.46875
if (y <= thresh) then
    ysq = 0.0d0
    if (y > xsmall)  ysq = y * y
    xnum = a(5) * ysq
    xden = ysq
    do i = 1, 3
        xnum = (xnum + a(i)) * ysq
        xden = (xden + b(i)) * ysq
    end do
    val = x * (xnum + a(4)) / (xden + b(4))

! Get complementary error function for 0.46875 <= |x| <= 4.0
else if (y <= 4.0d0) then
    xnum = c(9) * y
    xden = y
    do i = 1, 7
        xnum = (xnum + c(i)) * y
        xden = (xden + d(i)) * y
    end do
    val = (xnum + c(8)) / (xden + d(8))
    ysq = aint(16.0d0*y) / 16.0d0
    del = (y-ysq) * (y+ysq)
    val = exp(-ysq*ysq) * exp(-del) * val
    val = 1.0d0 - val
    if (x < 0.0d0)  val = -val
   
! Get complementary error function for |x| greater than 4.0
else
    val = 0.0d0
    if (y < xbig) then
        ysq = 1.0d0 / (y * y)
        xnum = p(6) * ysq
        xden = ysq
        do i = 1, 4
            xnum = (xnum + p(i)) * ysq
            xden = (xden + q(i)) * ysq
        end do
        val = ysq * (xnum + p(5)) / (xden + q(5))
        val = (sqrpi -  val) / y
        ysq = aint(16.0d0*y) / 16.0d0
        del = (y-ysq) * (y+ysq)
        val = exp(-ysq*ysq) * exp(-del) * val
    end if
    val = 1.0d0 - val
    if (x < 0.0d0)  val = -val
end if

end subroutine erfcore
!-------------------------------------------------------------------------------
END MODULE MATHFUNCTIONS
!-------------------------------------------------------------------------------
