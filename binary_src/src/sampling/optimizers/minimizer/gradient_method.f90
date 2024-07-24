!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/sampling/optimizers/minimizer/simplex.f90
!
! Description:
!   This module contains subroutines for gradient-based minimization algorithm
!   used in ligand docking.
!
!-------------------------------------------------------------------------------
MODULE GRADIENT_MINIMIZATION
!-------------------------------------------------------------------------------

use globals
use ligdock_energy, only: ligdock_energy_using_grid
use energy_vars
use energy_utils, only: update_R_for_ligand
use energy, only: calc_GALAXY_E
use mathfunctions, only: v_norm, bound_ang, R3_inverse_matrix
use geometry, only: quaternion, calc_g_rotmatrix, expmap_quat, rotation_matrix
use ligand_operator, only: construct_ligand

use in_out_ligand, only: write_mol2

implicit none
public

integer, parameter :: max_gene = 6 + max_lig_br - 1 

CONTAINS
!-------------------------------------------------------------------------------
subroutine convert_grad_to_internal(protein, ligand, gene, R_g, int_g)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: R_g(:,:)
real(dp), intent(in) :: gene(:)
type(ligand_type), intent(in) :: ligand
type(molecule_type), intent(inout) :: protein
real(dp), intent(out) :: int_g(size(gene)) !! current: only torsion angle
integer :: i_atm, atm_idx, k ,j, n
integer :: tor_init, tor_end, gene_dim, i_tor, tor_idx, i, i_gene
real(dp) :: g_tmp, R_init(3), R_end(3), arg1, arg2, p_abs, arg3
real(dp) :: p(3), p_sqr, cntr_R(3)
real(dp) :: q(4), Tq(3,3,4), Tp(3,3,3), dqdp(3,4), v(3), Tv_ij(3)
real(dp) :: axis(3), angle
real(dp) :: U_temp(3,3), p_tmp(3), U_temp2(3,3), q_tmp(4), test_R(3,ligand%n_atm), g_debug ! for debug
real(dp) :: U_inv(3,3), tmp_R(3)

gene_dim = size(gene)
!! Translation gradient
int_g = 0
do i_atm=1, ligand%n_atm
    int_g(1:3) = int_g(1:3) + R_g(1:3, i_atm)
end do

p=gene(4:6)
p_sqr = dot_product(p(:),p(:))

call expmap_quat(p,q)
call calc_g_rotmatrix(q(:), Tq)

!Rotation gradient
!First calculate dqdp, gradient of quaternion about p
if (p_sqr < 1.0d-6) then  ! Taylor expansion (because of singularity..)
    arg1 = 0.5d0 + p_sqr/48.0d0
    arg2 = (p_sqr/40.0d0 - 1.0d0)/24.0d0

    dqdp(:,1) = -0.5d0*p(:)*arg1
    dqdp(:,2) = p(1)*p(:)*arg2
    dqdp(:,3) = p(2)*p(:)*arg2
    dqdp(:,4) = p(3)*p(:)*arg2
    dqdp(1,2) = dqdp(1,2) + arg1
    dqdp(2,3) = dqdp(2,3) + arg1
    dqdp(3,4) = dqdp(3,4) + arg1
else
    p_abs = sqrt(p_sqr)
    arg1 = 0.5d0*cos(0.5d0*p_abs)/p_sqr
    arg2 = sin(0.5d0*p_abs)/p_abs
    arg3 = (arg1 - arg2/p_sqr)

    dqdp(:,1) = -0.5d0*arg2*p(:)
    dqdp(:,2) = arg3*p(1)*p(:)
    dqdp(:,3) = arg3*p(2)*p(:)
    dqdp(:,4) = arg3*p(3)*p(:)
    dqdp(1,2) = dqdp(1,2) + arg2
    dqdp(2,3) = dqdp(2,3) + arg2
    dqdp(3,4) = dqdp(3,4) + arg2
end if

! Then calculate T-matrix, gradient of R-matrix about p-vector
Tp(:,:,:) = 0.0d0
do k = 1, 4
    Tp(:,:,1) = Tp(:,:,1) + Tq(:,:,k)*dqdp(1,k)
    Tp(:,:,2) = Tp(:,:,2) + Tq(:,:,k)*dqdp(2,k)
    Tp(:,:,3) = Tp(:,:,3) + Tq(:,:,k)*dqdp(3,k)
end do

! Finally, calculate drdp
call rotation_matrix(3, q, U_temp)
call R3_inverse_matrix(U_temp)
cntr_R= protein%ligand(ligand%lig_no)%R(1:3, ref_lig(ligand%lig_type)%cntr_atm)
do i_atm=1, ligand%n_atm
    if (i_atm==ref_lig(ligand%lig_type)%cntr_atm) cycle
    v(1:3) = protein%ligand(ligand%lig_no)%R(:, i_atm) - cntr_R(:)
    v(1:3) = matmul(U_temp, v(:))
    do j = 1, 3
        Tv_ij(1:3) = 0.0d0
        do n=1,3
            Tv_ij(n) = dot_product(Tp(:,n,j), v(:))
        end do

        int_g(j+3) = int_g(j+3) + dot_product(Tv_ij(:), R_g(:,i_atm))
    end do
end do

do i_tor=1, ligand%n_br-1
    tor_init = ligand%bridge(1, i_tor+1) 
    tor_end = ligand%bridge(2, i_tor+1)

    R_init = protein%ligand(ligand%lig_no)%R(:,tor_init)
    R_end = protein%ligand(ligand%lig_no)%R(:,tor_end)

    tor_idx = ligand%i_rot_tor(i_tor)
    angle = protein%ligand(ligand%lig_no)%t_ang(tor_idx)
    axis = R_end(:) - R_init(:)
    call v_norm(axis(:))

    do i_atm=1, ligand%n_atm_br(i_tor+1)
        atm_idx =ligand%atm_in_br(i_atm, i_tor+1)
        if (atm_idx == tor_end) cycle ! same as do I_ATM=2

        angle = protein%ligand(ligand%lig_no)%t_ang(tor_idx)
        tmp_R = protein%ligand(ligand%lig_no)%R(:,atm_idx) - R_init
    !
        call cartesian_grad_to_rotation_grad(axis, angle, tmp_R, &
                                             R_g(:,atm_idx), g_tmp)
    !
        int_g(i_tor+6) = int_g(i_tor+6) + g_tmp
    end do
end do
    
end subroutine convert_grad_to_internal
!-------------------------------------------------------------------------------
subroutine cartesian_grad_to_rotation_grad(axis, angle, crd, dfdR, dfd0)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: axis(3), angle, crd(3), dfdR(3)
real(dp), intent(out) :: dfd0
real(dp) :: dqd0(4), q(4), Tq(3,3,4), T0(3,3), U(3,3)
real(dp) :: tan_w, tan1, cosine, tan_sqr, v(3)
integer :: k,i,n
real(dp) :: tmp_angle, q_tmp(4), U_tmp(3,3), crd_after(3)

! convert torsional gradient
call quaternion(axis, angle, q)
call calc_g_rotmatrix(q(:), Tq)

tan_w = tan(0.25d0*angle)
tan_sqr = tan_w*tan_w
tan1 = 1.0d0+tan_sqr
cosine = (1.0d0 - tan_sqr)/tan1

!First calculate dqd0, grad of quaternion about the angle
dqd0(1) = -tan_w/tan1
dqd0(2:4) = 0.5d0*cosine*axis(:)

! Then calculate T-matrix, gradient of R-matrix about the angle
T0 = 0.0d0
do k=1,4
    T0(:,:) = T0(:,:) + Tq(:,:,k)*dqd0(k)
end do

call rotation_matrix(3, q, U)
call R3_inverse_matrix(U)

v(:) = matmul(U,crd)

dfd0 = 0.0d0
do i=1,3
    dfd0 = dfd0 + dot_product(T0(:,i), v)*dfdR(i)
end do

end subroutine cartesian_grad_to_rotation_grad
!-------------------------------------------------------------------------------
subroutine opt_ligand_steepest_descent_intcrd(E, status, gene, protein, ligand, nft)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: protein
type(ligand_type), intent(in) :: ligand
logical, intent(out) :: status
real(dp), intent(inout) :: gene(:)
real(dp), intent(out) :: E(0:n_E_component_ligdock)
integer, intent(out) :: nft
type(energy_type) :: ff, ff_tmp
real(dp) :: f_tmp, cond
real(dp) :: alpha, max_iter, eps, norm2, sigma
logical :: calc_g = .true.
real(dp) :: R_g(3,ligand%n_atm), g_tmp(3, ligand%n_atm)
integer :: i_iter, i_atm, i,j 
real(dp) :: f2, f1, n
integer :: line_iter, gene_dim
real(dp) :: int_g(max_gene), gene_before(max_gene), int_g_tmp(max_gene)
type(ligand_residue_type) :: tmp_ligand

gene_dim = 6+ligand%n_br-1

i_iter = 0
max_iter =200 ! max interation
eps = 0.010d0 !stop condition
cond = 1.0d0 ! initially larger than eps
sigma=0.05d0

call construct_ligand(protein, ligand, gene)
call update_R_for_ligand(protein%ligand(ligand%lig_no), ligand%lig_no)
call ligdock_energy_using_grid(protein, ligand, ff, R_g, calc_g)
call convert_grad_to_internal(protein, ligand, gene, R_g, int_g)
f_tmp = ff%total
gene_before = gene

nft = 1
outer : do while (cond >eps .and. i_iter < max_iter)
    ! do line search (backtracking)
    alpha = 1.0d0
    line_iter = 0

    norm2 = dot_product(int_g(:gene_dim), int_g(:gene_dim))

    tmp_ligand = protein%ligand(ligand%lig_no)

    do 
        protein%ligand(ligand%lig_no) = tmp_ligand
        f2 = f_tmp - sigma*alpha*norm2
            
        !for trans/rot
        do i=1, 6
            gene(i) = gene_before(i) - alpha*int_g(i)
        end do

        !for torsion
        do j=7, gene_dim
            gene(i) = bound_ang(gene_before(i)-alpha*int_g(i))
        end do

        call construct_ligand(protein, ligand, gene)
        call update_R_for_ligand(protein%ligand(ligand%lig_no), ligand%lig_no)
        call ligdock_energy_using_grid(protein, ligand, ff_tmp, g_tmp, .false.)
        nft = nft+1
        f1 = ff_tmp%total

        if (f1>f2) then
            alpha=alpha*0.5d0
            line_iter = line_iter + 1
            if (line_iter > max_iter) then ! line search failed
                call update_R_for_ligand(protein%ligand(ligand%lig_no), ligand%lig_no)
                exit outer
            end if
        else
            exit
        end if
    end do

    gene_before=gene

    call ligdock_energy_using_grid(protein, ligand, ff, R_g, calc_g)
    call convert_grad_to_internal(protein, ligand, gene, R_g, int_g)
    nft = nft+1

!    f = ff%total
    i_iter = i_iter + 1
    cond = abs(ff%total-f_tmp)
    f_tmp = ff%total
end do outer

E(:) = ff%ligdock(:)

end subroutine opt_ligand_steepest_descent_intcrd
!-------------------------------------------------------------------------------
subroutine opt_ligand_steepest_descent_intcrd_nag(E, status, gene, protein, ligand)
!-------------------------------------------------------------------------------
type(molecule_type), intent(inout) :: protein
type(ligand_type), intent(in) :: ligand
logical, intent(out) :: status
real(dp), intent(inout) :: gene(:)
real(dp), intent(out) :: E(0:n_E_component_ligdock)
type(energy_type) :: ff, ff_tmp
real(dp) :: f_tmp, cond
real(dp) :: alpha, max_iter, eps, norm2, sigma
logical :: calc_g = .true.
real(dp) :: R_g(3,ligand%n_atm), g_tmp(3, ligand%n_atm)
integer :: i_iter, i_atm, i,j 
real(dp) :: f2, f1, n
integer :: nft, line_iter, gene_dim
real(dp) :: int_g(max_gene), gene_before(max_gene), int_g_tmp(max_gene)
type(ligand_residue_type) :: tmp_ligand
real(dp) :: v(max_gene), v_before(max_gene), mu

gene_dim = 6+ligand%n_br-1

i_iter = 0
max_iter =200 ! max interation
eps = 0.010d0 !stop condition
cond = 1.0d0 ! initially larger than eps
sigma=0.05d0

call construct_ligand(protein, ligand, gene)
call update_R_for_ligand(protein%ligand(ligand%lig_no), ligand%lig_no)
call ligdock_energy_using_grid(protein, ligand, ff, R_g, calc_g)
call convert_grad_to_internal(protein, ligand, gene, R_g, int_g)
f_tmp = ff%total
gene_before = gene

v=0.0d0
v_before=0.0d0
mu=0.9d0

nft = 0
outer : do while (cond >eps .and. i_iter < max_iter)
    ! do line search (backtracking)
    alpha = 1.0d0
    line_iter = 0

    norm2 = dot_product(int_g(:gene_dim), int_g(:gene_dim))

    tmp_ligand = protein%ligand(ligand%lig_no)

    do 
        protein%ligand(ligand%lig_no) = tmp_ligand
        f2 = f_tmp - sigma*alpha*norm2
            
        v=v_before
        !for trans/rot
        do i=1, 6
            v(i) = mu*v(i) - alpha*int_g(i)
            gene(i) = gene_before(i) - mu*v_before(i) + (1.0d0+mu)*v(i)
        end do

        !for torsion
        do j=7, gene_dim
            v(i) = mu*v(i) - alpha*int_g(i)
            gene(i) = bound_ang(gene_before(i) - mu*v_before(i) + (1.0d0+mu)*v(i))
        end do

        call construct_ligand(protein, ligand, gene)
        call update_R_for_ligand(protein%ligand(ligand%lig_no), ligand%lig_no)
        call ligdock_energy_using_grid(protein, ligand, ff_tmp, g_tmp, .false.)
        nft = nft+1
        f1 = ff_tmp%total

        if (f1>f2) then
            alpha=alpha*0.5d0
            line_iter = line_iter + 1
            if (line_iter > max_iter) then ! line search failed
                call update_R_for_ligand(protein%ligand(ligand%lig_no), ligand%lig_no)
                exit outer
            end if
        else
            exit
        end if
    end do

    gene_before=gene

    call ligdock_energy_using_grid(protein, ligand, ff, R_g, calc_g)
    call convert_grad_to_internal(protein, ligand, gene, R_g, int_g)
    nft = nft+1

!    f = ff%total
    i_iter = i_iter + 1
    cond = abs(ff%total-f_tmp)
    f_tmp = ff%total

    v_before=v
end do outer

E(:) = ff%ligdock(:)

end subroutine opt_ligand_steepest_descent_intcrd_nag
!-------------------------------------------------------------------------------
END MODULE GRADIENT_MINIMIZATION
!-------------------------------------------------------------------------------
