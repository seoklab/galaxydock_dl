!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/core/utils/clustering.f90
!
! Description:
!
!-------------------------------------------------------------------------------
MODULE CLUSTERING
!-------------------------------------------------------------------------------
use globals, only: dp
use ran, only: random
use sort, only: sort1, sort2

implicit none
public

CONTAINS
!-------------------------------------------------------------------------------
subroutine kmeans_clustering(d, n, k, id, cl_size)
!-------------------------------------------------------------------------------
! Run K-means clustering
!  
! k: K (# to be clustered)
! d(n,n): distance between members
! id(n,k): output cluster members
! cl_size(k): output cluster size
!-------------------------------------------------------------------------------
integer, intent(in) :: n, k
real(dp), intent(in) :: d(n,n)
integer, intent(out) :: id(n,k), cl_size(k)
integer, parameter :: max_iter = 30
integer :: iter, i_cl, i, i_mem, j_mem
integer :: cl_no, i_no, j_no
integer :: id_prv(n,k)
integer :: key(n)
real(dp) :: cl_d(n), d_min
logical :: used(n), converge, is_seed(n)

id(:,:) = 0
cl_size(:) = 0
! Initialize randomly
used(:) = .false.
do i_cl = 1, k
    do
        i = int(random()*n) + 1
        if (.not. used(i)) exit
    end do
    used(i) = .true.
    cl_size(i_cl) = 1
    id(1,i_cl) = i
    is_seed(i) = .true.
end do

id_prv(:,:) = id(:,:)

do iter = 1, max_iter
    ! Get member
    do i = 1, n
        if (is_seed(i)) cycle
        d_min = 999999.0d0
        do i_cl = 1, k
            if (d_min > d(id(1,i_cl),i)) then
                d_min = d(id(1,i_cl),i)
                cl_no = i_cl
            end if
        end do
        cl_size(cl_no) = cl_size(cl_no) + 1
        id(cl_size(cl_no),cl_no) = i
    end do
    
    ! Check convergence
    converge = .true.
    do i_cl = 1, k
        do i_mem = 1, cl_size(i_cl)
            if (id_prv(i_mem,i_cl) /= id(i_mem,i_cl)) then
                converge = .false.
                exit
            end if
        end do
        if (.not. converge) exit
    end do
    if (converge .or. iter == max_iter) exit
   
    ! Get seed
    id_prv(:,:) = id(:,:)
    is_seed(:) = .false.
    do i_cl = 1, k
        cl_d(:) = 0.0d0
        do i_mem = 1, cl_size(i_cl)
            i_no = id(i_mem,i_cl)
            do j_mem = 1, cl_size(i_cl)
                j_no = id(j_mem,i_cl)
                cl_d(i_mem) = cl_d(i_mem) + d(i_no,j_no)
            end do
        end do
        call sort2(cl_size(i_cl), cl_d(1:cl_size(i_cl)), key(1:cl_size(i_cl)))
        id(1,i_cl) = id_prv(key(1),i_cl)
        id(2:n,i_cl) = 0
        cl_size(i_cl) = 1
        is_seed(id_prv(key(1),i_cl)) = .true.
    end do
end do

end subroutine kmeans_clustering
!-------------------------------------------------------------------------------
subroutine find_nearest_cluster(n, d, dij, cl_id, cl_size, min_d)
!-------------------------------------------------------------------------------
! Find nearest cluster
! d(n,n): distance between members
! cl_id(n,n): cluster members
! dij(2): cluster index (ith-cluster, jth-cluster)
!-------------------------------------------------------------------------------
integer, intent(in) :: n, cl_id(n,n), cl_size(n)
integer, intent(out) :: dij(2)
real(dp), intent(in) :: d(n,n)
real(dp), intent(out), optional :: min_d
integer :: i, j, ik, jk
real(dp) :: dmin, dsum

dij(:) = 0
do i = 1, n-1
    do j = i+1, n
        dsum = 0.0d0
        do ik = 1, cl_size(i)
            do jk = 1, cl_size(j)
                dsum = dsum + d(cl_id(ik,i),cl_id(jk,j))
            end do
        end do
        if (cl_size(i)*cl_size(j) == 0) cycle
        dsum = dsum / (cl_size(i)*cl_size(j))
        if (dij(1) == 0 .and. dij(2) == 0) then
            dij(1:2) = (/i,j/)
            dmin = dsum
        else if (dsum < dmin) then
            dij(1:2) = (/i,j/)
            dmin = dsum
        end if
    end do
end do

if (present(min_d)) then
    min_d = dmin
end if

end subroutine find_nearest_cluster
!-------------------------------------------------------------------------------
subroutine hierarchical_clustering(d, n, n_cl, w, cl_id, cl_size)
!-------------------------------------------------------------------------------
! Run Hierarchical culstering written by Lim Heo.
!  
! d(n,n): distance between members
! w(n): input weight parameters, highest member in a cluster would be the cl
!       center
! cl_id(n,n): output cluster members
! cl_size(k): output cluster size
!
! CAUTION: output array cl_id is not sorted, while cl_id(1,k) is the cluster center
! CAUTION: output array cl_id and cl_size is n dimensional (not n_cl)
!-------------------------------------------------------------------------------
integer, intent(in) :: n, n_cl
real(dp), intent(in) :: d(n,n), w(n)
integer, intent(out) :: cl_id(n,n), cl_size(n)
integer :: i_cl, i, dij(2), id0(n), k_cl(n)
real(dp) :: d_cl(n,n), w_cl(n)

d_cl(:,:) = d(:,:)
cl_id(:,:) = 0
do i_cl = 1, n
    cl_id(1,i_cl) = i_cl
    cl_size(i_cl) = 1
end do

do i_cl = n, n_cl+1, -1
    call find_nearest_cluster(n, d, dij, cl_id, cl_size)
    do i = 1, cl_size(dij(2))
        cl_size(dij(1)) = cl_size(dij(1)) + 1
        cl_id(cl_size(dij(1)),dij(1)) = cl_id(i,dij(2))
    end do
    cl_size(dij(2)) = 0
end do

do i_cl = 1, n
    if (cl_size(i_cl) == 0) cycle
    do i = 1, cl_size(i_cl)
        w_cl(i) = w(cl_id(i,i_cl))
    end do
    id0(1:cl_size(i_cl)) = cl_id(1:cl_size(i_cl),i_cl)
    call sort2(cl_size(i_cl), w_cl(1:cl_size(i_cl)), k_cl(1:cl_size(i_cl)))
    do i = 1, cl_size(i_cl)
        cl_id(i,i_cl) = id0(k_cl(i))
    end do
end do

!-------------------------------------------------------------------------------
end subroutine hierarchical_clustering
!-------------------------------------------------------------------------------
subroutine hierarchical_clustering_dcut(d, n, d_cut, w, cl_id, cl_size, n_cl)
!-------------------------------------------------------------------------------
! Run Hierarchical culstering written by Lim Heo.
!  
! d(n,n): distance between members
! w(n): input weight parameters, highest member in a cluster would be the cl
!       center
! cl_id(n,n): output cluster members
! cl_size(k): output cluster size
!
! CAUTION: output array cl_id is not sorted, while cl_id(1,k) is the cluster center
! CAUTION: output array cl_id and cl_size is n dimensional (not n_cl)
!-------------------------------------------------------------------------------
integer, intent(in) :: n
real(dp), intent(in) :: d(n,n), w(n), d_cut
integer, intent(out) :: cl_id(n,n), cl_size(n), n_cl
integer :: i_cl, i, dij(2), id0(n), k_cl(n)
real(dp) :: d_cl(n,n), w_cl(n), d_min
integer :: cl_id0(n,n), cl_size0(n), cl_idx

d_cl(:,:) = d(:,:)
cl_id(:,:) = 0
do i_cl = 1, n
    cl_id(1,i_cl) = i_cl
    cl_size(i_cl) = 1
end do

do i_cl = n, 1, -1
    call find_nearest_cluster(n, d, dij, cl_id, cl_size, d_min)
    if (d_min > d_cut) exit
    do i = 1, cl_size(dij(2))
        cl_size(dij(1)) = cl_size(dij(1)) + 1
        cl_id(cl_size(dij(1)),dij(1)) = cl_id(i,dij(2))
    end do
    cl_size(dij(2)) = 0
end do

n_cl = 0
do i_cl = 1, n
    if (cl_size(i_cl) == 0) cycle
    n_cl = n_cl + 1
    do i = 1, cl_size(i_cl)
        w_cl(i) = w(cl_id(i,i_cl))
    end do
    id0(1:cl_size(i_cl)) = cl_id(1:cl_size(i_cl),i_cl)
    call sort2(cl_size(i_cl), w_cl(1:cl_size(i_cl)), k_cl(1:cl_size(i_cl)))
    do i = 1, cl_size(i_cl)
        cl_id(i,i_cl) = id0(k_cl(i))
    end do
end do

call sort1(n, cl_size(1:n), k_cl(1:n))
cl_id0 = cl_id
cl_size0 = cl_size
cl_idx = 0
do i_cl = n, 1, -1
    cl_idx = cl_idx + 1
    cl_size(cl_idx) = cl_size0(i_cl)
    cl_id(:,cl_idx) = cl_id0(:, k_cl(i_cl))
end do

end subroutine hierarchical_clustering_dcut
!-------------------------------------------------------------------------------
subroutine largest_clustering(d, n, w, cl_id, cl_size, dcut, n_cl)
!-------------------------------------------------------------------------------
! Run Largest clustering made by PHB.
!  
! d(n,n): distance between members
! w(n): input weight(?) parameters
! cl_id(n,n): output cluster members
! cl_size(k): output cluster size
! dcut: distance cutoff
! n_cl: # of clusters
!-------------------------------------------------------------------------------
integer, intent(in) :: n
real(dp), intent(in) :: d(n,n), w(n), dcut
integer, intent(out) :: cl_id(n,n), cl_size(n), n_cl
logical :: remain(n)
integer :: i, j, i_mem, n_remain, n_max
integer :: id_tmp(n), n_tmp, id_max(n)
real(dp) :: d_tmp(n,n), w_tmp(n)

! Initialize
cl_id(:,:) = 0
cl_size(:) = 0
remain(:) = .true.
n_cl = 0
n_remain = n

do
    ! Resort index
    i = 0
    id_tmp(:) = 0
    w_tmp(:) = 0.0d0
    do i_mem = 1, n
        if (.not. remain(i_mem)) cycle
        i = i + 1
        id_tmp(i) = i_mem
        w_tmp(i) = w(i_mem)
    end do
   
    n_tmp = i
    d_tmp(:,:) = 0.0d0
    do i = 1, n_tmp-1
        do j = i+1, n_tmp
            d_tmp(i,j) = d(id_tmp(i),id_tmp(j))
            d_tmp(j,i) = d(id_tmp(j),id_tmp(i))
        end do
    end do
   
    ! Finish if singles left
    call get_largest_cluster(d_tmp(1:n_tmp,1:n_tmp), n_tmp, w_tmp(1:n_tmp), &
                             dcut, id_max(1:n_tmp), n_max)
   
    n_cl = n_cl + 1
    n_remain = n_remain - n_max
   
    cl_size(n_cl) = n_max
    do i_mem = 1, n_max
        cl_id(i_mem,n_cl) = id_tmp(id_max(i_mem))
        remain(id_tmp(id_max(i_mem))) = .false.
    end do
    
    if (n_remain == 0) exit
end do

end subroutine largest_clustering
!-------------------------------------------------------------------------------
subroutine get_largest_cluster(d, n, w, dcut, id, nmax) 
!-------------------------------------------------------------------------------
! Find the largest cluster using distance and weights
! First, find close member using dcut.
! Second, weights are summed up.
! Finally, member with the highest weight is the largest cluster.
!
! d(n,n): distance matrix between members
! n: # of members
! w(n): weight parameters
! dcut: cutoff distance
! id(n): members of the largest cluster
! nmax: # of the largest cluster members
!-------------------------------------------------------------------------------
integer, intent(in) :: n
real(dp), intent(in) :: d(n,n), dcut, w(n)
integer, intent(out) :: id(n), nmax
integer :: i_mem, j_mem, i_center
integer :: id_close(n,n), n_close(n)
real(dp) :: mscore, score(n)

n_close(:) = 1
id_close(:,:) = 0
id(:) = 0
score(:) = w(:)

do i_mem = 1, n
    id_close(1,i_mem) = i_mem
end do

do i_mem = 1, n-1
    do j_mem = i_mem+1, n
        if (d(i_mem,j_mem) < dcut) then
            n_close(i_mem) = n_close(i_mem) + 1
            n_close(j_mem) = n_close(j_mem) + 1
     
            score(i_mem) = score(i_mem) + w(j_mem)
            score(j_mem) = score(j_mem) + w(i_mem)
     
            id_close(n_close(i_mem),i_mem) = j_mem
            id_close(n_close(j_mem),j_mem) = i_mem
        end if
    end do
end do

i_center = maxloc(score, dim=1)
mscore = score(i_center)
nmax = n_close(i_center)
id(1:nmax) = id_close(1:nmax,i_center)
     
end subroutine get_largest_cluster
!-------------------------------------------------------------------------------
subroutine NMRclust(d, n, w, cl_id, cl_size, n_cl)
!-------------------------------------------------------------------------------
! Run NMRclust written by Woong-Hee Shin. (Modified by Minkyung Baek)
!  
! d(n,n): distance between members
! w(n): input weight parameters (e.g. energy), lowest member in a cluster 
!       would be the cl center
! cl_id(n,n): output cluster members (cl_id(1,i_cl) is i_cl th cluster center)
! cl_size(n): output cluster size
! n_cl : # of clusters determined by NMRclust
!-------------------------------------------------------------------------------
integer, intent(in) :: n
real(dp), intent(in) :: d(n,n), w(n)
integer, intent(out) :: cl_id(n,n), cl_size(n), n_cl
!
real(dp) :: av_sp(n-1), av_sp_norm(n-1), penalty(n-1)
real(dp) :: w_cl(n)
integer :: h_cl_id(n,n,n-1), h_cl_size(n,n-1)
integer :: i, i_cl, i_step, n_step, dij(2)
integer :: min_loc(1), id0(n), k_cl(n)
integer :: cl_id0(n,n), cl_size0(n), cl_idx

n_step = n-1

! initialize
do i_cl = 1, n
    cl_id(1, i_cl) = i_cl
    cl_size(i_cl) = 1
end do

do i_step = 1, n_step
    ! find two nearest clusters & merge them
    call find_nearest_cluster(n, d, dij, cl_id, cl_size)
    do i = 1, cl_size(dij(2))
        cl_size(dij(1)) = cl_size(dij(1)) + 1
        cl_id(cl_size(dij(1)),dij(1)) = cl_id(i,dij(2))
    end do
    cl_size(dij(2)) = 0
    ! calculate average spread in this step.
    call avg_spread(d, cl_id, cl_size, n, av_sp(i_step))
    ! store history
    h_cl_id(:,:,i_step) = cl_id(:,:)
    h_cl_size(:,i_step) = cl_size(:)
end do

call normalize_spread(av_sp, av_sp_norm, n_step, n)
call step_penalty(av_sp_norm, penalty, n_step)

min_loc = minloc(penalty)
n_cl = 0
do i_cl = 1, n
    if (h_cl_size(i_cl, min_loc(1)) == 0) cycle
    n_cl = n_cl + 1
    cl_size(n_cl) = h_cl_size(i_cl, min_loc(1))
    cl_id(:,n_cl) = h_cl_id(:,i_cl,min_loc(1))
end do

do i_cl = 1, n_cl
    do i = 1, cl_size(i_cl)
        w_cl(i) = w(cl_id(i,i_cl))
    end do
    id0(1:cl_size(i_cl)) = cl_id(1:cl_size(i_cl),i_cl)
    call sort2(cl_size(i_cl), w_cl(1:cl_size(i_cl)), k_cl(1:cl_size(i_cl)))
    do i = 1, cl_size(i_cl)
        cl_id(i,i_cl) = id0(k_cl(i))
    end do
end do

call sort1(n_cl, cl_size(1:n_cl), k_cl(1:n_cl))
cl_id0 = cl_id
cl_size0 = cl_size
cl_idx = 0
do i_cl = n_cl, 1, -1
    cl_idx = cl_idx + 1
    cl_size(cl_idx) = cl_size0(i_cl)
    cl_id(:,cl_idx) = cl_id0(:, k_cl(i_cl))
end do

end subroutine NMRclust
!-------------------------------------------------------------------------------
subroutine avg_spread(d, cl_id, cl_size, n, av_sp)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: d(:,:)
integer, intent(in) :: cl_id(:,:), cl_size(:), n
real(dp), intent(out) :: av_sp
integer :: i_cl, i_mem, j_mem
real(dp) :: spr, n_mem, n_cl

av_sp = 0.0d0

n_cl = 0.0d0
do i_cl = 1, n
    if (cl_size(i_cl) < 2) cycle
    spr = 0.0d0
    n_mem = 0.0d0
    do i_mem = 1, cl_size(i_cl)-1
        do j_mem = i_mem + 1, cl_size(i_cl)
            spr = spr + d(cl_id(i_mem,i_cl), cl_id(j_mem, i_cl))
            n_mem = n_mem + 1.0d0
        end do
    end do
    n_cl = n_cl + 1.0d0
    av_sp = av_sp + spr/n_mem
end do

av_sp = av_sp / n_cl
    
end subroutine avg_spread
!-------------------------------------------------------------------------------
subroutine normalize_spread(av_sp, av_sp_norm, n_step, n)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: av_sp(:)
real(dp), intent(out) :: av_sp_norm(:)
integer, intent(in) :: n_step, n
real(dp) :: min_val, max_val, norm
integer :: i_step

min_val = minval(av_sp(1:n_step))
max_val = maxval(av_sp(1:n_step))

do i_step = 1, n_step
    norm = dble(n-2)
    norm = norm / (max_val - min_val)
    norm = norm * (av_sp(i_step) - min_val)
    norm = norm + 1.0d0
    av_sp_norm(i_step) = norm
end do

end subroutine normalize_spread
!-------------------------------------------------------------------------------
subroutine step_penalty(av_sp_norm, penalty, n_step)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: av_sp_norm(:)
real(dp), intent(out) :: penalty(:)
integer, intent(in) :: n_step
integer :: i_step

do i_step = 1, n_step
    penalty(i_step) = av_sp_norm(i_step) + (n_step - i_step + 1)
end do

end subroutine step_penalty
!-------------------------------------------------------------------------------
subroutine find_centroid(crds, ncrd, natm, cntrd_crd, d2c, i_cntr)
!-------------------------------------------------------------------------------
! Find centroid for multiple copy of structures
!
! crds(3,natm,ncrd): input structures (multiple)
! cntrd_crd: Centroid coordinate
! d2c: distance to centroid for multiple inputs
! i_cntr: Structure index for the closest
!-------------------------------------------------------------------------------
integer, intent(in) :: ncrd, natm
real(dp), intent(in) :: crds(3,natm,ncrd)
integer, intent(out) :: i_cntr
real(dp), intent(out) :: cntrd_crd(3,natm), d2c(ncrd)
real(dp) :: dv(3), d2c_sort(ncrd)
integer :: i_crd, i_atm, key(ncrd)

cntrd_crd(:,:) = 0.0d0

do i_crd = 1, ncrd
    cntrd_crd(:,:) = cntrd_crd(:,:) + crds(:,:,i_crd)
end do
cntrd_crd(:,:) = cntrd_crd(:,:) / dble(ncrd)

d2c(:) = 0.0d0
do i_crd = 1, ncrd
    do i_atm = 1, natm
        dv(:) = crds(:,i_atm,i_crd) - cntrd_crd(:,i_atm)
        d2c(i_crd) = d2c(i_crd) + dot_product(dv(:),dv(:))
    end do
    d2c(i_crd) = d2c(i_crd)/dble(natm)
    d2c(i_crd) = sqrt(d2c(i_crd))
end do

d2c_sort(:) = d2c(:)

call sort2(ncrd, d2c_sort, key)

i_cntr = key(1)

end subroutine find_centroid
!-------------------------------------------------------------------------------
subroutine DBSCAN(d, n, eps, minPts, id, cl_size, n_cl, cl_noise)
!-------------------------------------------------------------------------------
! Run DBSCAN clustering (density-based clustering)
! d(n,n): distance between members (input)
! eps: Distance cutoff to define neighborhood
! minPts: the minimum number of points required to form a dense region within
! eps (input). it should be >= 3.
!
! id(n,n): output cluster members
! cl_size(n): output cluster size
! n_cl: number of cluster
! 
! Please refer wikipedia page for DBSCAN
!-------------------------------------------------------------------------------
integer, intent(in) :: n
real(dp), intent(in) :: d(:,:) ! d(n,n)
real(dp), intent(in) :: eps
integer, intent(in) :: minPts
integer, intent(out) :: id(:,:), cl_size(:) ! id(n,n), cl_size(n)
integer, intent(out) :: n_cl
logical, intent(in), optional :: cl_noise
!
integer, parameter :: UNDEFINED = -1
integer, parameter :: NOISE = -2
integer :: i, j, i_neigh, seed
integer :: label(n), neigh(n), n_neigh, n_queue
integer :: tmp_neigh(n), n_tmp_neigh
logical :: is_in
logical :: noise_to_cluster

if (present(cl_noise)) then
    noise_to_cluster = cl_noise
else
    noise_to_cluster = .false.
end if

label(:) = UNDEFINED
n_cl = 0

do i = 1, n
    if (label(i) /= UNDEFINED) cycle
    call find_neighbors_for_DBSCAN(d, n, i, eps, neigh, n_neigh)
    if (n_neigh < minPts) then
        label(i) = NOISE
        cycle
    end if
    !
    n_cl = n_cl + 1
    n_queue = n_neigh
    i_neigh = 0
    do while (n_queue > 0)
        i_neigh = i_neigh + 1
        seed = neigh(i_neigh)
        n_queue = n_queue - 1 ! extract seed from queue
        if (label(seed) == NOISE) label(seed) = n_cl
        if (label(seed) /= UNDEFINED) cycle
        label(seed) = n_cl
        call find_neighbors_for_DBSCAN(d, n, seed, eps, tmp_neigh, n_tmp_neigh)
        if (n_tmp_neigh >= minPts) then
            do j = 1, n_tmp_neigh
                call check_is_in_neighbor(tmp_neigh(j), neigh, n_neigh, is_in)
                if (.not. is_in) then
                    n_neigh = n_neigh + 1
                    n_queue = n_queue + 1
                    neigh(n_neigh) = tmp_neigh(j)
                end if
            end do
        end if
    end do
end do

if (noise_to_cluster) then
    do i = 1, n
        if (label(i) /= NOISE) cycle
        call find_neighbors_for_DBSCAN(d, n, i, eps, neigh, n_neigh)
        is_in = .false.
        do j = 1, n_neigh
            if ((label(neigh(j)) /= UNDEFINED) .and. &
                (label(neigh(j)) /= NOISE)) then 
                ! neighbor j assigned to cluster
                is_in = .true.
                exit
            end if
        end do
        if (is_in) cycle  ! it has neighbor assigned to some cluster. treat as noise
        ! it is isolated noise clsuters. add it to new cluster
        n_cl = n_cl + 1
        do j = 1, n_neigh
            label(neigh(j)) = n_cl
        end do
    end do
end if

cl_size(:) = 0
do i = 1, n
    if (label(i) == NOISE) cycle
    cl_size(label(i)) = cl_size(label(i)) + 1
    id(cl_size(label(i)), label(i)) = i
end do

end subroutine DBSCAN
!-------------------------------------------------------------------------------
subroutine DBSCAN_2(d, n, eps_in, minPts, id, cl_size, n_cl, cl_noise)
!-------------------------------------------------------------------------------
! Run DBSCAN clustering (density-based clustering)
! d(n,n): distance between members (input)
! eps: Distance cutoff to define neighborhood 
! minPts: the minimum number of points required to form a dense region within
! eps (input). it should be >= 3.
!
! id(n,n): output cluster members
! cl_size(n): output cluster size
! n_cl: number of cluster
! 
! Please refer wikipedia page for DBSCAN
!-------------------------------------------------------------------------------
integer, intent(in) :: n
real(dp), intent(in) :: d(:,:) ! d(n,n)
real(dp), intent(in) :: eps_in
integer, intent(in) :: minPts
integer, intent(out) :: id(:,:), cl_size(:) ! id(n,n), cl_size(n)
integer, intent(out) :: n_cl
logical, intent(in), optional :: cl_noise
!
integer, parameter :: UNDEFINED = -1
integer, parameter :: NOISE = -2
integer :: i, j, i_neigh, seed
integer :: label(n), neigh(n), n_neigh, n_queue
integer :: tmp_neigh(n), n_tmp_neigh
logical :: is_in
logical :: noise_to_cluster
integer :: key(n)
real(dp) :: eps, dist(n), eps_s(n), dist_to_line(n)
real(dp) :: denom

if (present(cl_noise)) then
    noise_to_cluster = cl_noise
else
    noise_to_cluster = .false.
end if

! set eps as elbow points of minPts nearest distance plot
do i = 1, n
    dist(:) = d(:,i)
    call sort2(n, dist, key)
    eps_s(i) = dist(minPts)
end do
call sort2(n, eps_s, key)

! distance to the line connecting tww pts (x1, y1), (x2,y2) from point (x0, y0)
! d = abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)/sqrt((y2-y1)**2+(x2-x1)**2)
denom = sqrt((eps_s(n)-eps_s(1))**2 + dble(n-1)**2)
do i = 1, n
    dist_to_line(n) = abs((eps_s(n)-eps_s(1))*i - (n-1)*eps_s(i) &
                          + n*eps_s(1) - 1*eps_s(n))/denom
end do
call sort2(n, dist_to_line, key)
eps = min(eps_s(key(n)), eps_in)
eps = max(eps, eps_in*0.5d0)
print*, "eps is set to ", eps


label(:) = UNDEFINED
n_cl = 0

do i = 1, n
    if (label(i) /= UNDEFINED) cycle
    call find_neighbors_for_DBSCAN(d, n, i, eps, neigh, n_neigh)
    if (n_neigh < minPts) then
        label(i) = NOISE
        cycle
    end if
    !
    n_cl = n_cl + 1
    n_queue = n_neigh
    i_neigh = 0
    do while (n_queue > 0)
        i_neigh = i_neigh + 1
        seed = neigh(i_neigh)
        n_queue = n_queue - 1 ! extract seed from queue
        if (label(seed) == NOISE) label(seed) = n_cl
        if (label(seed) /= UNDEFINED) cycle
        label(seed) = n_cl
        call find_neighbors_for_DBSCAN(d, n, seed, eps, tmp_neigh, n_tmp_neigh)
        if (n_tmp_neigh >= minPts) then
            do j = 1, n_tmp_neigh
                call check_is_in_neighbor(tmp_neigh(j), neigh, n_neigh, is_in)
                if (.not. is_in) then
                    n_neigh = n_neigh + 1
                    n_queue = n_queue + 1
                    neigh(n_neigh) = tmp_neigh(j)
                end if
            end do
        end if
    end do
end do

if (noise_to_cluster) then
    do i = 1, n
        if (label(i) /= NOISE) cycle
        call find_neighbors_for_DBSCAN(d, n, i, eps, neigh, n_neigh)
        is_in = .false.
        do j = 1, n_neigh
            if ((label(neigh(j)) /= UNDEFINED) .and. &
                (label(neigh(j)) /= NOISE)) then 
                ! neighbor j assigned to cluster
                is_in = .true.
                exit
            end if
        end do
        if (is_in) cycle  ! it has neighbor assigned to some cluster. treat as noise
        ! it is isolated noise clsuters. add it to new cluster
        n_cl = n_cl + 1
        do j = 1, n_neigh
            label(neigh(j)) = n_cl
        end do
    end do
end if

cl_size(:) = 0
do i = 1, n
    if (label(i) == NOISE) cycle
    cl_size(label(i)) = cl_size(label(i)) + 1
    id(cl_size(label(i)), label(i)) = i
end do

end subroutine DBSCAN_2
!-------------------------------------------------------------------------------
subroutine find_neighbors_for_DBSCAN(d, n, P, eps, neigh, n_neigh)
!-------------------------------------------------------------------------------
real(dp), intent(in) :: d(:,:), eps
integer, intent(in) :: n, P
integer, intent(inout) :: neigh(:), n_neigh
!
integer :: Q

n_neigh = 0
do Q = 1, n
    if (d(P,Q) <= eps) then
        n_neigh = n_neigh + 1
        neigh(n_neigh) = Q
    end if
end do

end subroutine find_neighbors_for_DBSCAN
!-------------------------------------------------------------------------------
subroutine check_is_in_neighbor(mem, neigh, n_neigh, is_in)
!-------------------------------------------------------------------------------
integer, intent(in) :: mem, neigh(:), n_neigh
logical, intent(out) :: is_in
integer :: i

is_in = .false.
do i = 1, n_neigh
    if (neigh(i) == mem) then
        is_in = .true.
        exit
    end if
end do

end subroutine check_is_in_neighbor
!-------------------------------------------------------------------------------
END MODULE CLUSTERING
!-------------------------------------------------------------------------------
