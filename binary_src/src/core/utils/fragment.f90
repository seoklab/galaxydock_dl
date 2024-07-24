!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!                     and Prof. Julian Lee @ Soongsil Univ.
!
! File: $GALAXY/src/core/utils/fragment.f90
!
! Description:
!  Implementation of Fragment assembly method developed by Julian Lee 
!   Reference : Julian Lee et. al. Prediction of protein teriary structure
!               using PROFESY, a novel method based on fragment assembly and
!               conformational space annealing, Proteins, 56 704-714.
! Note:
!  1. To control the 'stiffness' of fragment connection criteria, modify
!    [tor_cut1], [tor_cut2] values
! 
!  2. The subroutines in this module would not work for fragment libraries
!    with length other than 9, yet. This might be fixed in the future.
!
!-------------------------------------------------------------------------------
MODULE FRAGMENT
!-------------------------------------------------------------------------------
use globals
use logger, only: log_p, terminate_with_error
use ran,    only: random
use string, only: parse_longstring
use mathfunctions, only: bound_ang, multiple_binormal, hamming_distance
use sort, only: sort2

implicit none
save
private

!===============================================================================
! PARAMETERS
!===============================================================================
integer, parameter :: max_frag_lib   =  10 ! maximum number of fragment library file
integer, parameter :: max_frag_neigh = 200 ! max number of fragments at each residue
!
integer :: n_frag_lib     ! no. fragment library read (only ONE is available yet)
integer :: max_frag_len   ! maximum residue length of fragment
logical :: use_fragment   ! Whether to use fragment library
!
!Torsion difference criteria for connection evaluation
real(dp), parameter :: tor_cut1 = 30.0d0*deg2rad ! criterion for fragment connection
real(dp), parameter :: tor_cut2 = 90.0d0*deg2rad ! Min. criterion for fragment-stem connection
real(dp), parameter :: tor_cut3 = two_pi         ! Max. criterion for fragment-stem connection
!
!===============================================================================
! DERIVED TYPES
!===============================================================================
type frag_type
!-------------------------------------------------------------------------------
! A derived type for fragments
!-------------------------------------------------------------------------------
real(dp) :: phi, psi, omg           ! Torsion angles
real(dp) :: score                   ! ramachandran score of fragment
!-------------------------------------------------------------------------------
end type frag_type
!-------------------------------------------------------------------------------
type frag_lib_type
!-------------------------------------------------------------------------------
! A derived type to store fragment library information
!-------------------------------------------------------------------------------
character(len=len_fname) :: frag_lib_file     ! filename
!
integer :: n_len       ! length of fragment residue. 
integer :: n_neigh     ! number of multiple fragments at single fragment window collect
                       ! from different pdb.
integer :: end_res     ! last residue of the library
integer :: pos_shift   ! residue position shift (generally n_len/2)
!-------------------------------------------------------------------------------
end type frag_lib_type
!-------------------------------------------------------------------------------
type(frag_lib_type) :: frag_lib(max_frag_lib)  ! Fragment library as a global variable
type(frag_type), allocatable :: frag(:,:,:,:)  ! Fragment
!
!===============================================================================
! PUBLIC
!===============================================================================
!public :: frag_type
!public :: frag_lib_type
public :: max_frag_lib
public :: max_frag_len
public :: max_frag_neigh

public :: n_frag_lib
public :: frag_lib

public :: use_fragment
public :: initialize_fragment
public :: finalize_fragment
public :: read_frag_lib
public :: build_tor_bank
public :: assemble_fragments
public :: assign_fragment

public :: select_fragment
public :: assemble_fragments_diverse

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_fragment(tnr)
!-------------------------------------------------------------------------------
! Initialize fragment library arrays
!-------------------------------------------------------------------------------
integer, intent(in) :: tnr

allocate(frag(max_frag_len,max_frag_neigh,tnr,n_frag_lib))
call sort_by_frag_len()

end subroutine initialize_fragment
!-------------------------------------------------------------------------------
subroutine finalize_fragment()
!-------------------------------------------------------------------------------
! Deallocate frag(frag_type)
!-------------------------------------------------------------------------------
deallocate(frag)

end subroutine finalize_fragment
!-------------------------------------------------------------------------------
subroutine read_frag_lib(tnr)
!-------------------------------------------------------------------------------
! Initial setup of fragment library,
!  by reading ROSETTA-format fragment library file
!-------------------------------------------------------------------------------
integer, intent(in) :: tnr
integer :: f_unit, ioerror, num_word, openstat
character(len=len_fname) :: word(30)
character(len=500) :: line
type(frag_type) :: tmp
integer :: i_frag, i_res, i_neigh, k  ! Frag_lib, Position, frag, frag_res

! Read additional info
do i_frag = 1, n_frag_lib
    frag_lib(i_frag)%n_neigh = 0
    frag_lib(i_frag)%end_res = 0
   
    f_unit = 33
    open(f_unit, file = trim(frag_lib(i_frag)%frag_lib_file), iostat = openstat)

    write(log_msg,"(A,A)") '- Reading fragment library file: ', trim(frag_lib(i_frag)%frag_lib_file)
    call log_p(log_msg, me=me, level=20)
   
    if (openstat > 0) then
        call terminate_with_error("Terminate with error: No fragment file found. Please check {infile_frag}.")
    end if
   
    i_res = 0
    i_neigh = 0
    k = 0
   
    do
        read(f_unit, "(A500)", iostat = ioerror) line
        if (ioerror < 0) exit
      
        word(:) = ''
        call parse_longstring(line, num_word, word, 500)
        
        if (word(1) == 'position:') then ! Residue no. divider
            i_neigh = 1
            read(word(2),"(I4)") i_res
            if (frag_lib(i_frag)%n_neigh == 0) then
                read(word(4),"(I4)") frag_lib(i_frag)%n_neigh
            end if
      
        else if (num_word == 0) then ! Frag. neigh. divider
            if (k /= 0) i_neigh = i_neigh + 1
            k = 0
      
        else if (num_word >= 8) then
            k = k + 1
          
            read(word(6),*) tmp%phi
            read(word(7),*) tmp%psi
            read(word(8),*) tmp%omg
          
            ! UPDATED
            if (i_res <= 0 .or. i_res > tnr) cycle
          
            tmp%phi = tmp%phi * deg2rad
            tmp%psi = tmp%psi * deg2rad
            tmp%omg = tmp%omg * deg2rad
            if (i_neigh <= max_frag_neigh) then
                frag(k,i_neigh,i_res,i_frag) = tmp
            end if
        end if
    end do
    
    if (i_res == 0) then
        call terminate_with_error('Error: Problem in fragment library. Please check your library file.')
    end if
   
    frag_lib(i_frag)%end_res = i_res
   
    close(f_unit)
end do

end subroutine read_frag_lib
!-------------------------------------------------------------------------------
subroutine build_tor_bank(i_frag, stem, res1, res2, nres, n_conf, bank, assemble_whole)
!-------------------------------------------------------------------------------
! This is fragment assembly method developed by Julian Lee and
! his co-workers. The code was implemented by PHB as current version.
! This subroutine generates 'n_conf' fragment-assembled conformations 
! from res1 to res2 into 'bank'.
!-------------------------------------------------------------------------------
integer, intent(in) :: i_frag
real(dp), intent(in) :: stem(3,frag_lib(i_frag)%pos_shift,2)
integer, intent(in) :: res1, res2, nres, n_conf 
real(dp), intent(inout) :: bank(3,nres,n_conf)
logical, intent(in) :: assemble_whole

integer :: pos_shift, frag_len
integer :: i_conf, j_conf, k, i_k, curr_k, start_k, iw, resno, curr_res, i_res
integer :: trial, conn_failed, n_cand, i_neigh, frag_id(frag_lib(i_frag)%n_neigh), i_sel
integer :: minw, i1
integer ::  continuity(nres), cont_count, i_seg
real(dp) :: bbtor(3, -frag_lib(i_frag)%pos_shift+1:nres+frag_lib(i_frag)%pos_shift)    ! Save bbtor tmp
real(dp) :: tor_frag(2,frag_lib(i_frag)%n_len), tor_curr(2,frag_lib(i_frag)%n_len)
real(dp) :: conn_tor_cut
logical :: frag_used(frag_lib(i_frag)%n_neigh, -frag_lib(i_frag)%pos_shift+1:nres)
logical :: frag_ok, conn_ok, redundant
character(len=len_fname) :: pdbname(nres), prv_pdb
type(frag_type) :: curr_frag

pos_shift = frag_lib(i_frag)%pos_shift
frag_len = frag_lib(i_frag)%n_len
i_conf = 1
redundant = .false. !Check. for i_conf = 1

do while(i_conf <= n_conf)
    ! Initialization
    bbtor(:,:) = 0.0d0
    bbtor(:,-pos_shift+1:0) = stem(:,:,1)
    bbtor(:,nres+1:nres+pos_shift) = stem(:,:,2)
    frag_used(:,:) = .false.
   
    start_k = max(-pos_shift+1,-res1+2)
    k = start_k
    trial = 0
    conn_failed = 0
    conn_tor_cut = tor_cut2
   
    do while (k <= nres .and. k+res1-1 <= frag_lib(i_frag)%end_res .and. trial < 3*nres)
        resno = k + res1 -1 + pos_shift
        n_cand = 0
        do i_neigh = 1, frag_lib(i_frag)%n_neigh
            if (.not. frag_used(i_neigh,k)) then
                n_cand = n_cand + 1
                frag_id(n_cand) = i_neigh
            end if
        end do
      
        i_sel = frag_id(int(n_cand*random())+1)
        frag_used(i_sel,k) = .true.
      
        curr_frag = frag(1,i_sel,resno,i_frag)
        frag_ok = .true.
        do i_k = 1, frag_len
            tor_frag(1,i_k) = frag(i_k,i_sel,resno-pos_shift,i_frag)%phi
            tor_frag(2,i_k) = frag(i_k,i_sel,resno-pos_shift,i_frag)%psi
        end do
      
        if (resno <= pos_shift) then
            if (curr_frag%phi > two_pi) cycle
        else
            if (n_cand <= 0) then
                frag_used(:,k) = .false.
                cycle
            end if
          
            tor_curr(:,:) = 0.0d0
            do i_k = 1, frag_len
                if (k+i_k-1 <= nres+pos_shift) then
                    tor_curr(1:2,i_k) = bbtor(1:2,k+i_k-1)
                end if
            end do
          
            minw = max(-k+1, 1)
            call examine_end(i_frag, resno, minw, tor_curr(:,:), tor_frag(:,:), i1, frag_ok)
          
            conn_ok = .true.
            if (conn_failed > nres) then
                conn_tor_cut = conn_tor_cut + (tor_cut3 - tor_cut2)/(2.0d0*nres)
            end if
          
            if (.not. assemble_whole) then
                if (k == nres .and. res2+pos_shift < frag_lib(i_frag)%end_res) then
                    call examine_connection(i_frag, tor_frag(:,frag_len-pos_shift+1:frag_len), &
                         stem(:,1:frag_lib(i_frag)%pos_shift,2), conn_tor_cut, conn_ok)
                else if (k == start_k .and. res1-pos_shift > 1) then
                    call examine_connection(i_frag, tor_frag(:,1:pos_shift), &
                         stem(:,1:frag_lib(i_frag)%pos_shift,1), conn_tor_cut, conn_ok)
                end if
            end if
          
            do i_k = 1, frag_len
                if (abs(tor_frag(1,i_k)) < 1.0d-7) then
                    frag_ok = .false.
                end if
            end do
          
            if (.not. frag_ok .or. .not. conn_ok) then
                trial = trial + 1
                if (.not. conn_ok) then
                    conn_failed = conn_failed + 1
                end if
                cycle
            end if
        end if
      
        do iw = i1, frag_len
            curr_k = k+iw-1
            if (curr_k < 1) cycle
            if (curr_k > nres) exit
            curr_res = resno+iw-pos_shift
            if (curr_res >= 1 .and. curr_res <= frag_lib(i_frag)%end_res+frag_len-1) then
                curr_frag = frag(iw, i_sel, resno-pos_shift, i_frag)
                bbtor(1,curr_k) = curr_frag%phi
                bbtor(2,curr_k) = curr_frag%psi
                bbtor(3,curr_k) = curr_frag%omg
            end if
        end do
        
        resno = resno + 1
        k = k + 1
    end do
   
    !Check conformation about previously generated ones
    redundant = .false.
    do j_conf = 1, i_conf-1
        redundant = .true.
        do i_res = 1, nres
            if (abs(bbtor(1,i_res)-bank(1,i_res,j_conf)) > 0.001d0 .or.&
                abs(bbtor(2,i_res)-bank(2,i_res,j_conf)) > 0.001d0) then
                redundant = .false.
                exit
            end if
        end do
        if (redundant) exit
    end do
   
    ! Append only if not redundant against previous conformations
    if (.not. redundant) then
        bank(:,1:nres,i_conf) = bbtor(:,1:nres)
        i_conf = i_conf + 1
      
        !Check continuity of assembled fragment
        continuity(:) = 0
        i_seg = 1
        cont_count = 0
        prv_pdb = ''
        do i_res = 1, nres
            if (trim(pdbname(i_res)) /= trim(prv_pdb)) then
                continuity(i_seg) = cont_count
                i_seg = i_seg + 1
                cont_count = 1
            else
                cont_count = cont_count + 1
            end if
            prv_pdb = pdbname(i_res)
        end do
    end if
end do

end subroutine build_tor_bank
!-------------------------------------------------------------------------------
subroutine examine_end(i_frag, i_res, minw, bbtor1, bbtor2, i1, frag_ok)
!-------------------------------------------------------------------------------
! Compare current torsion between fragment and 
! check if the starting sites are similar to each other
!-------------------------------------------------------------------------------
integer, intent(in) :: i_frag
integer, intent(in) :: i_res, minw
real(dp), intent(inout) :: bbtor1(2,frag_lib(i_frag)%n_len), bbtor2(2,frag_lib(i_frag)%n_len)
integer, intent(out) :: i1
logical, intent(out) :: frag_ok
integer :: frag_len, iw
logical :: is_same(max_frag_neigh)
real(dp) :: diff(2,frag_lib(i_frag)%n_len)
  
frag_len = frag_lib(i_frag)%n_len
i1 = 0
is_same(:) = .false.

do iw = minw, frag_len-1
    call same(bbtor1(:,iw), bbtor2(:,iw), is_same(iw), diff(:,iw), tor_cut1)
end do

do iw = minw, frag_len - 1
    if(is_same(iw)) then
        i1 = iw
        exit
    !If not defined yet
    else if (bbtor1(1,iw) < 1.0d-7 .and. bbtor1(1,iw) < 1.0d-7) then
        i1 = iw
        exit
    end if
end do

if(i1 > 0 .or. (i1 == 0 .and. i_res <= frag_lib(i_frag)%pos_shift)) then
    frag_ok = .true.
else
    frag_ok = .false.
end if

end subroutine examine_end
!-------------------------------------------------------------------------------
subroutine examine_connection(i_frag, tor_frag, stemtor, tor_cut, connect_ok)
!-------------------------------------------------------------------------------
! Check connection of fragment about framework
!-------------------------------------------------------------------------------
integer, intent(in) :: i_frag
real(dp), intent(in) :: tor_frag(2,frag_lib(i_frag)%pos_shift)
real(dp), intent(in) :: stemtor(3,frag_lib(i_frag)%pos_shift), tor_cut
logical, intent(out) :: connect_ok
real(dp) :: diff(2), tor1(2), tor2(2)
logical :: is_same
integer :: i

connect_ok = .true.
do i = 1, frag_lib(i_frag)%pos_shift
    tor1 = stemtor(1:2,i)
    tor2 = tor_frag(:,i)
    call same(tor1, tor2, is_same, diff, tor_cut)
    if (.not. is_same) then
        connect_ok = .false.
        return
    end if
end do

end subroutine examine_connection
!-------------------------------------------------------------------------------
subroutine same(ang1, ang2, is_same, diff, cut)
!-------------------------------------------------------------------------------
! Determine if phi/psi in ang1 and ang2 are similar enough using 'cut' given
!-------------------------------------------------------------------------------
real(dp), intent(in) :: cut
real(dp), intent(inout) :: ang1(2), ang2(2)
logical, intent(out) :: is_same
real(dp), intent(out) :: diff(2)

diff(1) = bound_ang(ang1(1)-ang2(1))
diff(2) = bound_ang(ang1(2)-ang2(2))

if (diff(1)+diff(2) <= cut) then
    is_same = .true.
else
    is_same = .false.
end if
  
end subroutine same
!-------------------------------------------------------------------------------
subroutine sort_by_frag_len()
!-------------------------------------------------------------------------------
! sort fragment library by fragment length
!-------------------------------------------------------------------------------
type(frag_lib_type) :: tmp
integer :: i, j

do i = 1, n_frag_lib
    do j = 1, i-1
        if (frag_lib(i)%n_len > frag_lib(j)%n_len) then
            tmp = frag_lib(i)
            frag_lib(i) = frag_lib(j)
            frag_lib(j) = tmp
        end if
    end do
end do

end subroutine sort_by_frag_len
!-------------------------------------------------------------------------------
subroutine select_fragment(fragment, i_frag, i_res, n_neigh, select_mode, residues)
!-------------------------------------------------------------------------------
! Select a fragment for a given position (i_res) among a fragment library (i_frag)
!  with select_mode:
!  - random : select randomly between 1 and n_neigh
!  - score  : select maximum scored fragment
!  - similar: select the most similar fragment from the given residue torsions
!  - define : use the assigned n_neigh fragment only
! And returns fragment torsion angles phi/psi/omg
!-------------------------------------------------------------------------------
!real(dp), intent(out) :: fragment(3, frag_lib(i_frag)%n_len)
real(dp), intent(out) :: fragment(:,:)
integer, intent(in) :: i_frag, i_res, n_neigh
character(len=len_fname), intent(in) :: select_mode
type(residue_type), intent(in), optional :: residues(:)
integer :: i_f, i_neigh, i
real(dp) :: score_s(n_neigh)

! i_f -> selected fragment index
if (trim(select_mode) == 'random') then
    i_f = min(int(n_neigh*random())+1, n_neigh)
else if (trim(select_mode) == 'score') then
    do i_neigh = 1, n_neigh
        score_s(i_neigh) = frag(1,i_neigh,i_res,i_frag)%score
    end do
    i_f = maxloc(score_s, dim=1)
else if (trim(select_mode) == 'similar') then
    i_f = select_fragment_similar(i_frag,i_res,n_neigh, residues) 
else if (trim(select_mode) == 'define') then
    i_f = n_neigh
end if

do i = 1, frag_lib(i_frag)%n_len
    fragment(1,i) = frag(i, i_f, i_res, i_frag)%phi
    fragment(2,i) = frag(i, i_f, i_res, i_frag)%psi
    fragment(3,i) = frag(i, i_f, i_res, i_frag)%omg
    !
    ! How to apply fragment to residue 
    ! resNo = i_res+i-1
    ! i_ref = protein%residue(resNo)%res_type
    ! if (resNo > 1) then
    !    i_phi = ref_res(i_ref)%iphi_curr
    !    protein%residue(resNo)%t_ang(i_phi) = fragment(1,i)
    !    i_omg = ref_res(i_ref)%iomg_curr
    !    protein%residue(resNo)%t_ang(i_omg) = fragment(3,i)
    ! end if
    ! if (resNo < tnr) then
    !    i_psi = ref_res(i_ref)%ipsi_prev
    !    protein%residue(resNo+1)%t_ang(i_psi) = fragment(2,i)
    ! end if
    !
end do

end subroutine select_fragment
!-------------------------------------------------------------------------------
function select_fragment_similar(i_frag,i_res,n_neigh,residues)
!-------------------------------------------------------------------------------
! select fragment by hamming distance from n_neigh number fragments
!-------------------------------------------------------------------------------
integer :: select_fragment_similar
integer, intent(in) :: i_frag, i_res, n_neigh
type(residue_type), intent(in) :: residues(:)
integer :: i_neigh, n_len, resNo, i_ref, i_ref_next, i, k, ia
real(dp) :: X(frag_lib(i_frag)%n_len*2), Y(frag_lib(i_frag)%n_len*2)
real(dp) :: score_s(n_neigh)

n_len = frag_lib(i_frag)%n_len

i = 0
X(:) = 0.0d0
do resNo = i_res, i_res+n_len-1
    i_ref = residues(resNo)%res_type
    i_ref_next = residues(resNo+1)%res_type

    if (resNo > 1) then
        i = i+1
        ia = ref_res(i_ref)%iphi
        X(i) = residues(resNo)%t_ang(ia)
    end if

    if (resNo < frag_lib(i_frag)%end_res+n_len-1) then
        i = i+1
        ia = ref_res(i_ref_next)%ipsi
        X(i) = residues(resNo+1)%t_ang(ia)
    end if
end do

do i_neigh = 1, n_neigh
    Y(:) = 0.0d0
    do k = 1, n_len
        Y(2*k-1) = frag(k,i_neigh,i_res,i_frag)%phi
        Y(2*k  ) = frag(k,i_neigh,i_res,i_frag)%psi
    end do
    score_s(i_neigh) = hamming_distance(n_len*2, X, Y, .true.)
end do

select_fragment_similar = minloc(score_s, dim=1)

end function select_fragment_similar
!-------------------------------------------------------------------------------
subroutine assemble_fragments(i_frag, resNo, n_mer, n_conf, bank, useTop)
!-------------------------------------------------------------------------------
integer, intent(in) :: i_frag
integer, intent(in) :: resNo, n_mer
integer, intent(in) :: n_conf
integer, intent(in), optional :: useTop
real(dp), intent(out) :: bank(3, n_mer, n_conf)
!
integer :: i,j,k, k_mer, i_res, n_neigh, end_res
integer, allocatable :: index_s(:,:,:), n_index(:)
real(dp), allocatable :: torsion_s(:,:,:)
!
integer :: i_conf, n_trial
real(dp) :: conf(3, n_mer)

bank(:,:,:) = 0.0d0
!
k_mer = frag_lib(i_frag)%n_len  ! fragment length
if (present(useTop)) then
    n_neigh = min(useTop, frag_lib(i_frag)%n_neigh)
else
    n_neigh = frag_lib(i_frag)%n_neigh
end if
end_res = frag_lib(i_frag)%end_res
allocate(n_index(n_mer))
allocate(index_s(3, n_neigh*k_mer, n_mer))
allocate(torsion_s(3, n_neigh*k_mer, n_mer))
!
n_index(:) = 0
do i = 1, n_mer
    do j = 1, k_mer
        ! i_res-th fraglib., j-th torsions are corr. to resNo+(i-1)-th residue
        i_res = resNo + (i-1) - (j-1)
        if (i_res < 1 .or. i_res > end_res) cycle
        !
        do k = 1, n_neigh
            n_index(i) = n_index(i) + 1
            index_s(1:3, n_index(i), i) = (/i_res, k, j/)
            torsion_s(1, n_index(i), i) = frag(j, k, i_res, i_frag)%phi
            torsion_s(2, n_index(i), i) = frag(j, k, i_res, i_frag)%psi
            torsion_s(3, n_index(i), i) = frag(j, k, i_res, i_frag)%omg
        end do
    end do
end do
!
i_conf = 0
n_trial = 0
do while(i_conf < n_conf)
    n_trial = n_trial + 1
    !
    call sample_torsion(torsion_s, index_s, n_index, i_frag, n_mer, k_mer, conf)
    if (is_redundant(bank, conf, i_conf, n_mer)) cycle
    !
    i_conf = i_conf + 1
    bank(:,:, i_conf) = conf
end do
!
deallocate(index_s)
deallocate(n_index)
deallocate(torsion_s)

end subroutine assemble_fragments
!-------------------------------------------------------------------------------
subroutine assemble_fragments_diverse(angle, n_res, n_str)
!-------------------------------------------------------------------------------
! Fragment assembly covering vast torsional conformations
! Search from Top fragments
!-------------------------------------------------------------------------------
integer, intent(in) :: n_res, n_str  
real(dp), intent(out) :: angle(3, n_res, n_str)  

real(dp), parameter :: cutoff = 5.0d0 * deg2rad
integer, parameter :: i_frag = 1
integer :: i_res, j_res, k, n, i, l
integer :: n_fraglen, n_neigh
logical :: redundant
real(dp) :: score(frag_lib(i_frag)%n_neigh*frag_lib(i_frag)%n_len)
real(dp) :: val(3, frag_lib(i_frag)%n_neigh*frag_lib(i_frag)%n_len), val_now(3)
integer :: key(frag_lib(i_frag)%n_neigh*frag_lib(i_frag)%n_len)

n_fraglen = frag_lib(i_frag)%n_len
n_neigh = frag_lib(i_frag)%n_neigh

do i_res = 1, n_res
    ! get all angle pairs from frag_lib
    i = 0
    do k = 1, n_fraglen
        j_res = i_res + 1 - k
        if (j_res < 1 .or. j_res > n_res + 1 - n_fraglen) cycle
        !
        do n = 1, n_neigh
            i = i + 1
            val(1, i) = frag(k, n, j_res, i_frag)%phi
            val(2, i) = frag(k, n, j_res, i_frag)%psi
            val(3, i) = frag(k, n, j_res, i_frag)%omg
            score(i) = n + abs(k - (n_fraglen+1)/2) * 0.01d0
        end do
    end do
    !
    ! get angles from Top frags and centered frags
    call sort2(i, score(1:i), key(1:i))
    !
    n = 1
    angle(:, i_res, n) = val(:, key(1))
    do k = 2, i
        val_now = val(:, key(k))
        ! redundancy
        redundant = .false.
        do l = 1, n
            ! TODO: bound_ang
            if (maxval(abs(val_now(1:2) - angle(1:2, i_res, l))) < cutoff) then
                redundant = .true.
                exit
            end if
        end do
        if (.not. redundant) then
            n = n + 1
            angle(:, i_res, n) = val_now
            if (n == n_str) exit
        end if
    end do
    !
    if (n < n_str) then
        ! duplicate
        do k = n+1, n_str
            angle(:, i_res, k) = angle(:, i_res, mod(k-1,n)+1)
        end do
    end if
end do

end subroutine assemble_fragments_diverse
!-------------------------------------------------------------------------------
subroutine assign_fragment(i_frag, resNo, n_mer, n_conf, n_top, bank)
!-------------------------------------------------------------------------------
! Build torsion angle bank with assigned fragment, naively. Use torsion angles
! from top n_conf fragments.
!-------------------------------------------------------------------------------
integer, intent(in) :: i_frag
integer, intent(in) :: resNo, n_mer !n_mer is loop%n_res
integer, intent(in) :: n_conf
integer, intent(in) :: n_top
real(dp), intent(inout) :: bank(:, :, :)
integer :: n_neigh, k_mer, end_res
integer :: insert_mer
integer :: i_res, i_conf, k
character(len=len_fname) :: select_mode
real(dp), allocatable :: torsion(:,:)

k_mer = frag_lib(i_frag)%n_len !fragment length
n_neigh = min(n_top, frag_lib(i_frag)%n_neigh)
end_res = frag_lib(i_frag)%end_res !last residue of the library
select_mode = 'define'
insert_mer = min(n_mer, k_mer)
!
allocate(torsion(3, k_mer))
!
!Initialize torsion angle bank with given loop information
i_res = resNo
if (i_res < 1 .or. i_res > end_res) return
do i_conf = 1, n_conf
    k = random()*n_neigh
    if (k < 1) k = 1
    !assign k-th neighbor fragment among n_top randomly. Using full length fragment.
    call select_fragment(torsion(1:3, 1:k_mer), i_frag, i_res, k, select_mode)
    bank(:, 1:insert_mer, i_conf) = torsion(:, 1:insert_mer)
end do

deallocate(torsion)

end subroutine assign_fragment
!-------------------------------------------------------------------------------
subroutine sample_torsion(torsion_s, index_s, n_index, i_frag, n_mer, k_mer, conf)
!-------------------------------------------------------------------------------
! - torsion_s: gathered torsion angles
!   torsion_s(1:3, l, kk) = (phi,psi,omg) for kk-th position, l-th angles
! - index_s: indices for gathered torsion angles
!   index_s(1:3, l, kk) = (i_res,k,j)
!   - i_res: i_res-th potision among the fragment library
!   - k: k-th frag. of the i_res-th fraglib.
!   - j: j-th position of the i_res-th fraglib.
!-------------------------------------------------------------------------------
real(dp), intent(in) :: torsion_s(:,:,:)
integer, intent(in) :: index_s(:,:,:)
integer, intent(in) :: n_index(:)
integer, intent(in) :: i_frag, n_mer, k_mer
real(dp), intent(out) :: conf(:,:)
!
integer :: kk, ki, i_res, j, k, i
integer :: idx_s(3, n_mer)
real(dp) :: s(maxval(n_index)), si(maxval(n_index)), x0(2)
!
idx_s(:,:) = 0
do kk = 1, n_mer    ! Assemblying n_mer-peptides
    s(:) = 1.0d0    ! equal prob. for each torsion candidates
    do ki = 1, kk-1 ! Iterating over to calc. conditional prob.
        ! Evaluating only if
        !   prev.selected frag. for pos=ki can cover the pos=kk
        if (idx_s(3,ki)+(kk-ki) < k_mer) then
            i_res = idx_s(1, ki)
            k = idx_s(2, ki) 
            j = idx_s(3, ki) + (kk-ki)
            !
            ! (phi,psi) for the pos=kk from the prev.selected frag. for pos=ki
            x0(1) = frag(j, k, i_res, i_frag)%phi
            x0(2) = frag(j, k, i_res, i_frag)%psi
            !
            si(:) = 0.0d0
            do i = 1, n_index(kk)
                si(i) = torsion_prob(torsion_s(1:2, i, kk), x0)
            end do
            si(:) = si(:) / sum(si)
            s(:) = s(:) * si(:)
        end if
    end do
    !
    s(:) = s(:)/sum(s(1:n_index(kk)))
    i = select_based_on_prob(s, n_index(kk))
    conf(:,kk) = torsion_s(:, i, kk)
    idx_s(:,kk) = index_s(:, i, kk)
end do

end subroutine sample_torsion
!-------------------------------------------------------------------------------
function torsion_prob(x, x0)
!-------------------------------------------------------------------------------
! At first, sigma_sq was intended (30deg)**2, but (30deg) showed better results
!-------------------------------------------------------------------------------
real(dp) :: torsion_prob
real(dp), intent(in) :: x(2), x0(2)
real(dp), parameter :: sigma_sq = 30.0d0*deg2rad
real(dp) :: dx

dx = cos(x(1)-x0(1)) + cos(x(2)-x0(2))
torsion_prob = exp(-(2.0d0 - dx)/sigma_sq)

end function torsion_prob
!-------------------------------------------------------------------------------
function select_based_on_prob(s0, n)
!-------------------------------------------------------------------------------
integer :: select_based_on_prob
real(dp), intent(in) :: s0(:)
integer, intent(in) :: n
integer :: i
real(dp) :: s, r

s = 0.0d0
r = random()
do i = 1, n
    s = s + s0(i)
    if (s > r) then
        select_based_on_prob = i
        return
    end if
end do

select_based_on_prob = i

end function select_based_on_prob
!-------------------------------------------------------------------------------
function is_redundant(bank, conf, n_conf, n_mer)
!-------------------------------------------------------------------------------
logical :: is_redundant
real(dp), intent(in) :: bank(:,:,:)
real(dp), intent(in) :: conf(3, n_mer)
integer, intent(in) :: n_conf, n_mer
integer :: i
real(dp) :: dx(3, n_mer)

do i = 1, n_conf
    dx = abs(conf(:,:)-bank(:,:,i))
    if (maxval(dx) < 0.001d0) then
        is_redundant = .true.
        return
    end if
end do

is_redundant = .false.

end function is_redundant
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
end MODULE FRAGMENT
!-------------------------------------------------------------------------------
