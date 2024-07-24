!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/ppdock/interev.f90
!
! Description:
!   Statistical energy about 2-body and 3-body residue interactions.
!-------------------------------------------------------------------------------
!  Reference : Andreani, J. et. al. InterEvScore: a novel coarse-grained
!              interface scoring function using a multi-body statistical 
!              potential coupled to evolution, Bioinformatics, 2013,29,1742-1749
!-------------------------------------------------------------------------------
MODULE INTEREV
!-------------------------------------------------------------------------------
use globals
use logger
use mathfunctions, only: sigmoidal2
use convert_res_name, only: convert_to_stdres
use energy_vars, only: R

implicit none
save
private

character(3), parameter :: res_type2(20) = &
       (/ 'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'ILE', 'LEU', &
       'LYS', 'MET', 'PHE', 'PRO', 'TRP', 'VAL', 'SER', 'THR', 'TYR', 'HIS' /)
real(dp), parameter :: Rcut21 = 7.5d0
real(dp), parameter :: Rcut22 = 8.5d0
real(dp), parameter :: Rcut22_sqr = Rcut22**2
real(dp), parameter :: Rcut31 = 11.5d0
real(dp), parameter :: Rcut32 = 12.5d0
real(dp), parameter :: Rcut32_sqr = Rcut32**2
integer, allocatable :: res_id2(:), res_id3(:)
real(dp) :: body2(20,20), body3(6,6,6)

public :: initialize_interev
public :: finalize_interev
public :: calc_interev_score

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_interev(protein, use_pssm)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
logical, intent(in), optional :: use_pssm
character(4) :: res_name
integer :: i_res
  
allocate(res_id2(tn%stdres))
allocate(res_id3(tn%stdres))

do i_res = 1, protein%n_res
    res_name = protein%residue(i_res)%res_name
    call convert_to_stdres(res_name)
    call find_res_id2(res_name, res_id2(i_res))
    call find_res_id3(res_name, res_id3(i_res))
end do

body2(:, 1) = (/ -0.785d0, 0.266d0,-0.291d0,-0.044d0,-1.163d0, 0.169d0,&
            0.002d0,-0.689d0,-0.904d0,-0.782d0, 0.136d0,-0.766d0,-0.979d0,&
          -0.125d0,-1.010d0,-0.804d0,-0.286d0,-0.284d0,-0.678d0,-0.226d0 /)
body2(:, 2) = (/  0.266d0, 1.194d0, 0.327d0, 0.029d0,-1.046d0, 0.765d0,&
            0.055d0,-0.020d0, 0.062d0, 0.242d0, 1.561d0,-0.060d0,-0.276d0,&
            0.757d0,-0.662d0, 0.016d0, 0.435d0, 0.465d0,-0.166d0, 0.389d0 /)
body2(:, 3) = (/ -0.291d0, 0.327d0, 0.719d0, 0.619d0, 0.204d0, 0.565d0,&
            0.834d0,-0.073d0,-0.130d0,-0.091d0, 0.484d0, 0.075d0,-0.102d0,&
          -0.017d0,-0.510d0,-0.027d0, 0.322d0, 0.432d0,-0.138d0, 0.381d0 /)
body2(:, 4) = (/ -0.044d0, 0.029d0, 0.619d0, 1.202d0, 0.034d0, 0.771d0,&
            1.557d0, 0.203d0, 0.289d0, 0.310d0, 0.461d0, 0.520d0, 0.086d0,&
            0.375d0, 0.368d0, 0.369d0, 0.365d0, 0.403d0, 0.283d0, 0.106d0 /)
body2(:, 5) = (/ -1.163d0,-1.046d0, 0.204d0, 0.034d0,-2.114d0, 0.065d0,&
            0.431d0,-0.691d0,-1.292d0,-0.812d0,-0.554d0,-0.899d0,-0.934d0,&
          -0.687d0,-0.642d0,-0.655d0,-0.551d0,-0.373d0,-0.496d0,-0.259d0 /)
body2(:, 6) = (/  0.169d0, 0.765d0, 0.565d0, 0.771d0, 0.065d0, 1.418d0,&
            1.006d0, 0.208d0, 0.098d0, 0.312d0, 1.411d0, 0.206d0,-0.075d0,&
            0.550d0, 0.231d0, 0.319d0, 0.718d0, 0.763d0, 0.091d0, 0.643d0 /)
body2(:, 7) = (/  0.002d0, 0.055d0, 0.834d0, 1.557d0, 0.431d0, 1.006d0,&
            1.596d0, 0.388d0,-0.078d0, 0.268d0, 0.507d0, 0.379d0, 0.246d0,&
            0.720d0, 0.593d0, 0.255d0, 0.427d0, 0.603d0, 0.246d0, 0.487d0 /)
body2(:, 8) = (/ -0.689d0,-0.020d0,-0.073d0, 0.203d0,-0.691d0, 0.208d0,&
            0.388d0,-0.312d0,-0.617d0,-0.321d0, 0.224d0,-0.280d0,-0.747d0,&
          -0.198d0,-0.757d0,-0.571d0, 0.020d0, 0.178d0,-0.470d0, 0.091d0 /)
body2(:, 9) = (/ -0.904d0, 0.062d0,-0.130d0, 0.289d0,-1.292d0, 0.098d0,&
          -0.078d0,-0.617d0,-1.205d0,-0.947d0, 0.225d0,-0.802d0,-0.822d0,&
          -0.085d0,-0.888d0,-0.941d0,-0.243d0,-0.506d0,-0.697d0,-0.224d0 /)
body2(:,10) = (/ -0.782d0, 0.242d0,-0.091d0, 0.310d0,-0.812d0, 0.312d0,&
            0.268d0,-0.321d0,-0.947d0,-0.930d0, 0.067d0,-0.957d0,-0.624d0,&
            0.073d0,-0.354d0,-0.956d0,-0.025d0,-0.232d0,-0.520d0,-0.143d0 /)
body2(:,11) = (/  0.136d0, 1.561d0, 0.484d0, 0.461d0,-0.554d0, 1.411d0,&
            0.507d0, 0.224d0, 0.225d0, 0.067d0, 1.602d0, 0.032d0, 0.149d0,&
            0.844d0,-0.404d0, 0.025d0, 0.621d0, 0.678d0, 0.217d0, 0.927d0 /)
body2(:,12) = (/ -0.766d0,-0.060d0, 0.075d0, 0.520d0,-0.899d0, 0.206d0,&
            0.379d0,-0.280d0,-0.802d0,-0.957d0, 0.032d0,-0.171d0,-1.025d0,&
          -0.014d0,-1.091d0,-0.469d0,-0.174d0,-0.172d0,-0.758d0, 0.002d0 /)
body2(:,13) = (/ -0.979d0,-0.276d0,-0.102d0, 0.086d0,-0.934d0,-0.075d0,&
            0.246d0,-0.747d0,-0.822d0,-0.624d0, 0.149d0,-1.025d0,-0.776d0,&
          -0.594d0,-1.003d0,-0.789d0,-0.190d0,-0.378d0,-0.514d0,-0.325d0 /)
body2(:,14) = (/ -0.125d0, 0.757d0,-0.017d0, 0.375d0,-0.687d0, 0.550d0,&
            0.720d0,-0.198d0,-0.085d0, 0.073d0, 0.844d0,-0.014d0,-0.594d0,&
            0.143d0,-1.128d0,-0.182d0,-0.078d0, 0.150d0,-0.697d0,-0.151d0 /)
body2(:,15) = (/ -1.010d0,-0.662d0,-0.510d0, 0.368d0,-0.642d0, 0.231d0,&
            0.593d0,-0.757d0,-0.888d0,-0.354d0,-0.404d0,-1.091d0,-1.003d0,&
          -1.128d0,-1.060d0,-0.844d0,-0.389d0,-0.114d0,-0.591d0,-0.511d0 /)
body2(:,16) = (/ -0.804d0, 0.016d0,-0.027d0, 0.369d0,-0.655d0, 0.319d0,&
            0.255d0,-0.571d0,-0.941d0,-0.956d0, 0.025d0,-0.469d0,-0.789d0,&
          -0.182d0,-0.844d0,-0.801d0,-0.265d0,-0.510d0,-0.627d0,-0.186d0 /)
body2(:,17) = (/ -0.286d0, 0.435d0, 0.322d0, 0.365d0,-0.551d0, 0.718d0,&
            0.427d0, 0.020d0,-0.243d0,-0.025d0, 0.621d0,-0.174d0,-0.190d0,&
          -0.078d0,-0.389d0,-0.265d0, 0.144d0, 0.288d0,-0.047d0, 0.250d0 /)
body2(:,18) = (/ -0.284d0, 0.465d0, 0.432d0, 0.403d0,-0.373d0, 0.763d0,&
            0.603d0, 0.178d0,-0.506d0,-0.232d0, 0.678d0,-0.172d0,-0.378d0,&
            0.150d0,-0.114d0,-0.510d0, 0.288d0, 0.276d0,-0.129d0, 0.102d0 /)
body2(:,19) = (/ -0.678d0,-0.166d0,-0.138d0, 0.283d0,-0.496d0, 0.091d0,&
            0.246d0,-0.470d0,-0.697d0,-0.520d0, 0.217d0,-0.758d0,-0.514d0,&
          -0.697d0,-0.591d0,-0.627d0,-0.047d0,-0.129d0,-0.604d0,-0.145d0 /)
body2(:,20) = (/ -0.226d0, 0.389d0, 0.381d0, 0.106d0,-0.259d0, 0.643d0,&
            0.487d0, 0.091d0,-0.224d0,-0.143d0, 0.927d0, 0.002d0,-0.325d0,&
          -0.151d0,-0.511d0,-0.186d0, 0.250d0, 0.102d0,-0.145d0,-0.033d0 /)

body3(:,1,1) = (/-0.174d0,-0.427d0,-0.610d0,-0.739d0, 0.102d0,-0.035d0/)
body3(:,1,2) = (/-0.427d0,-0.316d0,-0.211d0,-0.366d0, 0.051d0, 0.108d0/)
body3(:,1,3) = (/-0.610d0,-0.211d0,-0.357d0,-0.390d0, 0.138d0,-0.142d0/)
body3(:,1,4) = (/-0.739d0,-0.366d0,-0.390d0,-0.416d0, 0.014d0,-0.234d0/)
body3(:,1,5) = (/ 0.102d0, 0.051d0, 0.138d0, 0.014d0, 0.486d0,-0.078d0/)
body3(:,1,6) = (/-0.035d0, 0.108d0,-0.142d0,-0.234d0,-0.078d0, 0.397d0/)
body3(:,2,1) = (/-0.427d0,-0.316d0,-0.211d0,-0.366d0, 0.051d0, 0.108d0/)
body3(:,2,2) = (/-0.316d0,-0.476d0,-0.103d0,-0.199d0, 0.168d0, 0.127d0/)
body3(:,2,3) = (/-0.211d0,-0.103d0,-0.019d0,-0.061d0, 0.394d0, 0.245d0/)
body3(:,2,4) = (/-0.366d0,-0.199d0,-0.061d0,-0.184d0, 0.273d0, 0.074d0/)
body3(:,2,5) = (/ 0.051d0, 0.168d0, 0.394d0, 0.273d0, 0.848d0, 0.140d0/)
body3(:,2,6) = (/ 0.108d0, 0.127d0, 0.245d0, 0.074d0, 0.140d0, 0.866d0/)
body3(:,3,1) = (/-0.610d0,-0.211d0,-0.357d0,-0.390d0, 0.138d0,-0.142d0/)
body3(:,3,2) = (/-0.211d0,-0.103d0,-0.019d0,-0.061d0, 0.394d0, 0.245d0/)
body3(:,3,3) = (/-0.357d0,-0.019d0,-0.046d0,-0.112d0, 0.570d0, 0.311d0/)
body3(:,3,4) = (/-0.390d0,-0.061d0,-0.112d0,-0.171d0, 0.312d0, 0.016d0/)
body3(:,3,5) = (/ 0.138d0, 0.394d0, 0.570d0, 0.312d0, 1.087d0, 0.305d0/)
body3(:,3,6) = (/-0.142d0, 0.245d0, 0.311d0, 0.016d0, 0.305d0, 0.472d0/)
body3(:,4,1) = (/-0.739d0,-0.366d0,-0.390d0,-0.416d0, 0.014d0,-0.234d0/)
body3(:,4,2) = (/-0.366d0,-0.199d0,-0.061d0,-0.184d0, 0.273d0, 0.074d0/)
body3(:,4,3) = (/-0.390d0,-0.061d0,-0.112d0,-0.171d0, 0.312d0, 0.016d0/)
body3(:,4,4) = (/-0.416d0,-0.184d0,-0.171d0,-0.394d0, 0.074d0,-0.123d0/)
body3(:,4,5) = (/ 0.014d0, 0.273d0, 0.312d0, 0.074d0, 0.980d0,-0.017d0/)
body3(:,4,6) = (/-0.234d0, 0.074d0, 0.016d0,-0.123d0,-0.017d0, 0.483d0/)
body3(:,5,1) = (/ 0.102d0, 0.051d0, 0.138d0, 0.014d0, 0.486d0,-0.078d0/)
body3(:,5,2) = (/ 0.051d0, 0.168d0, 0.394d0, 0.273d0, 0.848d0, 0.140d0/)
body3(:,5,3) = (/ 0.138d0, 0.394d0, 0.570d0, 0.312d0, 1.087d0, 0.305d0/)
body3(:,5,4) = (/ 0.014d0, 0.273d0, 0.312d0, 0.074d0, 0.980d0,-0.017d0/)
body3(:,5,5) = (/ 0.486d0, 0.848d0, 1.087d0, 0.980d0, 2.067d0, 0.469d0/)
body3(:,5,6) = (/-0.078d0, 0.140d0, 0.305d0,-0.017d0, 0.469d0, 0.268d0/)
body3(:,6,1) = (/-0.035d0, 0.108d0,-0.142d0,-0.234d0,-0.078d0, 0.397d0/)
body3(:,6,2) = (/ 0.108d0, 0.127d0, 0.245d0, 0.074d0, 0.140d0, 0.866d0/)
body3(:,6,3) = (/-0.142d0, 0.245d0, 0.311d0, 0.016d0, 0.305d0, 0.472d0/)
body3(:,6,4) = (/-0.234d0, 0.074d0, 0.016d0,-0.123d0,-0.017d0, 0.483d0/)
body3(:,6,5) = (/-0.078d0, 0.140d0, 0.305d0,-0.017d0, 0.469d0, 0.268d0/)
body3(:,6,6) = (/ 0.397d0, 0.866d0, 0.472d0, 0.483d0, 0.268d0, 1.282d0/)
  
end subroutine initialize_interev
!-------------------------------------------------------------------------------
subroutine finalize_interev()
!-------------------------------------------------------------------------------
  
deallocate(res_id2)
deallocate(res_id3)

end subroutine finalize_interev
!-------------------------------------------------------------------------------
subroutine calc_interev_score(f, g, appl_pair, calc_g)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: f(2), g(3,tn%atom,2)
logical, intent(in) :: appl_pair(tn%residue, tn%residue), calc_g
  
call calc_interev2(f(1), g(:,:,1), appl_pair, calc_g)
call calc_interev3(f(2), g(:,:,2), appl_pair, calc_g)

end subroutine calc_interev_score
!-------------------------------------------------------------------------------
subroutine calc_interev2(f, g, appl_pair, calc_g)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: f, g(3,tn%atom)
logical, intent(in) :: appl_pair(tn%residue, tn%residue), calc_g
integer :: res1, res2, atm1, atm2
integer :: type1, type2
real(dp) :: R_diff(3), Rsqr, dist
real(dp) :: score, weight(2), dfdr
  
f = 0.0d0
g(:,:) = 0.0d0
do res1 = 1, tn%recres
    atm1 = res_index(res1)%Ca_id(1)
    type1 = res_id2(res1)

    do res2 = tn%recres+1, tn%stdres
        atm2 = res_index(res2)%Ca_id(1)
        type2 = res_id2(res2)
        if (.not. appl_pair(res1,res2)) cycle
        
        R_diff(:) = R(:,atm2) - R(:,atm1)
        Rsqr = dot_product(R_diff(:), R_diff(:))
        if (Rsqr > Rcut22_sqr) cycle
        dist = sqrt(Rsqr)

        weight(:) = sigmoidal2(dist, Rcut21, Rcut22, calc_g)
        score = weight(1)*body2(type1, type2)
        f = f + score
        
        if (calc_g) then
            dfdr = weight(2)*body2(type1, type2)
            g(:,atm1) = g(:,atm1) - R_diff(:) * dfdr / dist
            g(:,atm2) = g(:,atm2) + R_diff(:) * dfdr / dist
        end if

    end do
end do
        
end subroutine calc_interev2
!-------------------------------------------------------------------------------
subroutine triangle_en_grad(f, g, atm1, atm2, atm3, type1, type2, type3, calc_g)
!-------------------------------------------------------------------------------
real(dp), intent(inout) :: f, g(3,tn%atom)
integer, intent(in) :: atm1, atm2, atm3, type1, type2, type3
logical :: calc_g
real(dp) :: R12(3), R13(3), R23(3), Rsqr, dist, inv_dist
real(dp) :: score, weight(2), dfdr

R12(:) = R(:,atm1) - R(:,atm2)
R13(:) = R(:,atm1) - R(:,atm3)
R23(:) = R(:,atm2) - R(:,atm3)
  
Rsqr = dot_product(R12,R12) + dot_product(R13,R13) + dot_product(R23,R23)
if (Rsqr > Rcut32_sqr) return

dist = sqrt(Rsqr)
inv_dist = 1.0d0 / dist
  
weight(:) = sigmoidal2(dist, Rcut31, Rcut32, calc_g)
score = weight(1)*body3(type1, type2, type3)
f = f + score
  
if (calc_g) then
    dfdr = weight(2)*body3(type1, type2, type3)
    g(:,atm1) = g(:,atm1) + dfdr*inv_dist*(R12(:)+R13(:))
    g(:,atm2) = g(:,atm2) + dfdr*inv_dist*(R23(:)-R12(:))
    g(:,atm3) = g(:,atm3) - dfdr*inv_dist*(R13(:)+R23(:))
end if

end subroutine triangle_en_grad
!-------------------------------------------------------------------------------
subroutine calc_interev3(f, g, appl_pair, calc_g)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: f, g(3,tn%atom)
logical, intent(in) :: appl_pair(tn%residue, tn%residue), calc_g
integer :: res1, res2, res3, atm1, atm2, atm3
integer :: type1, type2, type3
real(dp) :: R12(3), Rsqr

f = 0.0d0
g(:,:) = 0.0d0

do res1 = 1, tn%recres
    atm1 = res_index(res1)%Ca_id(1)
    type1 = res_id3(res1)
     
    do res2 = tn%recres+1, tn%stdres-2
        atm2 = res_index(res2)%Ca_id(1)
        type2 = res_id3(res2)

        R12(:) = R(:,atm1) - R(:,atm2)
        Rsqr = dot_product(R12(:),R12(:))
        if (Rsqr > Rcut32_sqr) cycle

        do res3 = res2+2, tn%stdres
            atm3 = res_index(res3)%Ca_id(1)
            type3 = res_id3(res3)

            call triangle_en_grad(f, g, atm1, atm2, atm3, &
            type1, type2, type3, calc_g)
           
        end do
    end do
end do

do res1 = tn%recres+1, tn%stdres
    atm1 = res_index(res1)%Ca_id(1)
    type1 = res_id3(res1)
     
    do res2 = 1, tn%recres-2
        atm2 = res_index(res2)%Ca_id(1)
        type2 = res_id3(res2)

        R12(:) = R(:,atm1) - R(:,atm2)
        Rsqr = dot_product(R12(:),R12(:))
        if (Rsqr > Rcut32_sqr) cycle

        do res3 = res2+2, tn%recres
            atm3 = res_index(res3)%Ca_id(1)
            type3 = res_id3(res3)

            call triangle_en_grad(f, g, atm1, atm2, atm3, &
            type1, type2, type3, calc_g)
           
        end do
    end do
end do

end subroutine calc_interev3
!-------------------------------------------------------------------------------
subroutine find_res_id2(res_name, res_id)
!-------------------------------------------------------------------------------
character(len=4), intent(in) :: res_name
character(3) :: res
integer, intent(out) :: res_id
integer :: i_res

res = trim(res_name)
do i_res = 1, 20
    if (res_type2(i_res) == res) then
        res_id = i_res
        exit
    end if
end do
  
end subroutine find_res_id2
!-------------------------------------------------------------------------------
subroutine find_res_id3(res_name, res_id)
!-------------------------------------------------------------------------------
character(len=4), intent(in) :: res_name
character(3) :: res
integer, intent(out) :: res_id
  
res = trim(res_name)
res_id = 0
if (res == 'TRP' .or. res == 'TYR' .or. res == 'PHE') then
    res_id = 1
else if (res=='MET'.or.res=='ILE'.or.res=='LEU'.or.res=='VAL'.or.res=='CYS') then
    res_id = 2
else if (res=='ALA'.or.res=='GLY'.or.res=='PRO') then
    res_id = 3
else if (res=='HIS'.or.res=='SER'.or.res=='THR'.or.res=='GLN'.or.res=='ASN') then
    res_id = 4
else if (res=='LYS'.or.res=='ARG') then
    res_id = 5
else if (res=='ASP'.or.res=='GLU') then
    res_id = 6
end if
  
end subroutine find_res_id3
!-------------------------------------------------------------------------------
END MODULE INTEREV
