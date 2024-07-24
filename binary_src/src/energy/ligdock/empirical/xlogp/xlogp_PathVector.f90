module xlogp_PathVector_m
    use xlogp_AtomIdVector_m, only : AtomIdVector_t, &
        Assignment_AIdVector, Is_equal_atomIdvector_to, &
        AIdVector_EqualPath, AIdVector_GetVectorSize, &
        AIdVector_Append, AIdVector_Clear, &
        AtomIdVector_Create, AtomIdVector_Destroy, &
        AIdVector_Revert, AIdVector_PathShow

    use xlogp_Atom_m
    use xlogp_constants_m

    implicit none
    !private !! Information hiding
    logical, private, parameter :: DEBUG = .true.

    real, parameter :: pathvector_growth_rate = 1.2
    public :: assignment(=)
    public :: PV_assignment

    type PathVector_t
        !private ! component information hiding
        integer :: num_used != 0
        type(AtomIdVector_t), dimension(XSCORE_MAX_PATHVECTOR_SIZE) :: data !=> null()
        logical :: validity != .false.
    end type PathVector_t

    interface assignment (=) ! overload =
        module procedure PV_assignment
    end interface

!    interface PathVector
!        module procedure newPathVector
!    end interface !PathVector

contains

    function newPathVector(capacity) result(return_value)
        integer, optional :: capacity
        type(PathVector_t) :: return_value

        if(present(capacity)) then
            call PV_create(return_value, capacity)
        else
            call PV_create(return_value, 1)
        endif
    end function newPathVector

    subroutine PV_create(vec, capacity)
        type(PathVector_t) :: vec
        integer, optional :: capacity
        integer :: cap

!        nullify(vec%data)

!        print*, '(*) PV_create : capacity = ', capacity
!        print*, '(*) DEBUG = ', DEBUG

        ! Check that the vector does not have any data left
!        if (associated(vec%data) .AND. size(vec%data) > 0 ) then
!            if(DEBUG) print*, "PathVector_create in before "
!            if(DEBUG) print*, "PathVector_create in before, size = ", size(vec%data)
!
!            call PV_destroy(vec)
!
!            if(DEBUG) print*, "create in ... after "
!        end if

        if (present(capacity)) then
            cap = max(1, capacity)
        else
            cap = 10
        endif

!        allocate(vec%data(1:cap))
        if(cap > XSCORE_MAX_PATHVECTOR_SIZE) then
            print*, '(*) OUT OF BOUND in XSCORE_MAX_PATHVECTOR_SIZE'
            print*, '(*) PATHVECTOR capacity was given to be ', cap
            print*, '(*) XSCORE_MAX_PATHVECTOR_SIZE = ', XSCORE_MAX_PATHVECTOR_SIZE
            stop
        endif

        call PV_initNumUsed(vec)

        vec%validity = .true.

    end subroutine PV_create

    subroutine PV_checkValid(vec)
        type(PathVector_t) :: vec
        if(.not. vec%validity) then
            call PV_create(vec)
        endif
    end subroutine PV_checkValid

    subroutine PV_destroy(vec)
        type(PathVector_t) :: vec

        ! Check that the vector does not have any data left
!        if (associated(vec%data) .AND. vec%num_used > 0) then
!            print*, '(*) Before deallocatge vec%data, vec%num_used = ', &
!                vec%num_used
!
!            deallocate(vec%data)
!        endif

        !vec%num_used = 0
        call PV_initNumUsed(vec)
    end subroutine PV_destroy

    subroutine PV_assignment(vleft, vright)
        type(PathVector_t), intent(inout) :: vleft
        type(PathVector_t), intent(in) :: vright
        integer :: numUsedRight

        integer :: vrightSize
        vrightSize = PV_getCapacitySize(vright)

!        if(.not. associated(vright%data)) then
!            return
!        endif

!        if(.not. associated(vleft%data)) then
!            call PV_create(vleft, vrightSize)
!        endif

        if(PV_getCapacitySize(vleft) < PV_getCapacitySize(vright)) then
            call PV_destroy(vleft)
            call PV_create(vleft, vrightSize)
        endif

        numUsedRight = PV_getNumUsed(vright)
        call PV_setNumUsed(vleft, numUsedRight)
        vleft%data(1:numUsedRight) = vright%data(1:numUsedRight)

    end subroutine PV_assignment

    !!----------------------------------------------------------------------------
    ! Return the number of elements in use
    function PV_getVectorSize(vec) result(return_value)
        type(PathVector_t) :: vec
        integer :: return_value

        return_value = PV_getNumUsed(vec)
    end function PV_getVectorSize

    integer function PV_getCapacitySize(vec)
        type(PathVector_t) :: vec

!        PV_getCapacitySize = size(vec%data)
        PV_getCapacitySize = XSCORE_MAX_PATHVECTOR_SIZE
    end function PV_getCapacitySize

    !!----------------------------------------------------------------------------
    !>
    !!----------------------------------------------------------------------------

    function PV_GetPath(vec, n) result(return_value)
        type(PathVector_t), intent(in out) :: vec
        integer, intent(in) :: n
        type(AtomIdVector_t) :: return_value
        type(AtomIdVector_t) :: empty_path

        if( n .lt. 1 .OR. n .gt. PV_getNumUsed(vec)) then
            print*, "(*) ERROR at xlogp_PathVector_at : Out of boundary. n = ", n, &
                " num_used = ", PV_getNumUsed(vec)
            call AtomIdVector_Create(empty_path)
            call AIdVector_Clear(empty_path)
            return_value = empty_path
        else
            !return_value = vec%data(n)
            call Assignment_AIdVector(return_value,vec%data(n))
        end if

    end function PV_GetPath

    subroutine PV_CopyPath(this, n, path)
        type(PathVector_t), intent(in) :: this
        integer, intent(in) :: n
        type(AtomIdVector_t), intent(out) :: path

        if( n .lt. 1 .OR. n .gt. PV_getNumUsed(this)) then
            print*, "(*) ERROR at xlogp_PathVector_at : Out of boundary. n = ", n, &
                " num_used = ", PV_getNumUsed(this)
            call AtomIdVector_Create(path)
        else
            call Assignment_AIdVector(path,this%data(n))
        end if

    end subroutine PV_CopyPath

!    subroutine PV_GetPathPtr(vec, n, return_value_ptr)
!        type(PathVector_t), intent(in out) :: vec
!        integer, intent(in) :: n
!        type(AtomIdVector_t), pointer :: return_value_ptr
!        type(AtomIdVector_t) :: empty_path
!
!        if( n .lt. 1 .OR. n .gt. PV_getNumUsed(vec)) then
!            print*, "(*) ERROR at xlogp_PathVector_at_ptr : Out of boundary. n = ", n, &
!                " num_used = ", PV_getNumUsed(vec)
!            call AtomIdVector_Create(empty_path)
!            call AIdVector_Clear(empty_path)
!            nullify(return_value_ptr) ! => null()
!        else
!            return_value_ptr => vec%data(n)
!        end if
!
!    end subroutine PV_GetPathPtr

    !!----------------------------------------------------------------------------
    !> Insert one or more empty elements and move original data backward.
    !!
    !! Arguments:
    !!  vec PathVector_t in question
    !!  pos Position to insert the empty elements
    !!  emptyNum Number of empty elements
    !!----------------------------------------------------------------------------
    subroutine PV_insert_empty(vec, pos, emptyNum)
        type(PathVector_t) :: vec
        integer, intent(in) :: pos
        integer, intent(in) :: emptyNum

        integer :: i
        type(AtomIdVector_t) :: empty_path

        call PV_checkValid(vec)

        ! Clear one atom
        call AtomIdVector_Create(empty_path)
        call AIdVector_Clear(empty_path)

        if (emptyNum .lt. 1 .OR. &
            pos .lt. 1 .OR. &
            pos .gt. PV_getNumUsed(vec) ) then
            return
        endif

        if (PV_getNumUsed(vec)+emptyNum .ge. size(vec%data) ) then
            call PV_increase_capacity(vec, PV_getNumUsed(vec) + emptyNum)
        endif

        ! From pos postion, insert emptyNum number of empty atoms
        ! E.g) for [a1,a2,a3,a4,a5,,],
        ! num_used=5
        ! insert_empty(3,2) =>
        ! a5 -> 5+2
        ! a4 -> 5+1
        !
        do i = PV_getNumUsed(vec), pos, -1
            vec%data(i+emptyNum) = vec%data(i)
        enddo

        ! Here, we made some dummy atoms as place holder.
        do i = 1, emptyNum
            vec%data(pos+i-1) = empty_path
        enddo

        call PV_setNumUsed(vec, PV_getNumUsed(vec) + emptyNum)

    end subroutine PV_insert_empty

    subroutine PV_delete_elements(vec, pos, numDelete)
        type(PathVector_t) :: vec
        integer, intent(in) :: pos
        integer, intent(in) :: numDelete

        integer :: i

        if (numDelete .lt. 1 .OR. &
            pos .lt. 1 .OR. &
            pos .gt. PV_getNumUsed(vec) ) then
            return
        endif

        do i=pos,PV_getNumUsed(vec)-numDelete
            vec%data(i) = vec%data(i+numDelete)
        enddo

        call PV_setNumUsed(vec, PV_getNumUsed(vec) - numDelete)
    end subroutine PV_delete_elements

    subroutine PV_delete_one(vec, pos)
        type(PathVector_t) :: vec
        integer, intent(in) :: pos

        integer :: i

        if (pos .lt. 1 .OR. &
            pos .gt. PV_getNumUsed(vec) ) then
            return
        endif

        if(pos == PV_getNumUsed(vec)) then
            call AIdVector_Clear(vec%data(pos))
            call PV_setNumUsed(vec, PV_getNumUsed(vec) - 1)
            return
        endif

        do i=pos, vec%num_used
            vec%data(i) = vec%data(i+1)
        enddo

        call PV_setNumUsed(vec, PV_getNumUsed(vec) - 1)
    end subroutine PV_delete_one

    subroutine PV_clear(vec)
        type(PathVector_t) :: vec
        call PV_clearAllPaths(vec)
    end subroutine PV_clear

    subroutine PV_clearAllPaths(vec)
        type(PathVector_t) :: vec

        call PV_delete_elements(vec, 1, PV_getNumUsed(vec))
    end subroutine PV_clearAllPaths

    !!----------------------------------------------------------------------------
    !> Remove only one atom which corresponds to the given atom in the vector.
    !!----------------------------------------------------------------------------
    subroutine PV_remove(vec, path)
        type(PathVector_t) :: vec
        type(AtomIdVector_t), intent(in) :: path
        type(AtomIdVector_t) :: p1
        integer :: i
        integer :: nLen
        integer :: removeCnt

        removeCnt = 0
        nLen = PV_getVectorSize(vec)
        call AtomIdVector_Create(p1)

        do i=1, nLen
            p1 = vec%data(i)
            if(AIdVector_EqualPath(path, p1)) then
                call PV_delete_one(vec, i)
                removeCnt = removeCnt + 1
            endif
        enddo


    !        if(.true.) then
    !            print*, "(*) Path Remove Cnt = ", removeCnt
    !        endif
    end subroutine PV_remove

    !!----------------------------------------------------------------------------
    !> Remove only one atom which corresponds to the given atom in the vector.
    !!----------------------------------------------------------------------------
    subroutine PV_removeRevert(vec, path)
        type(PathVector_t) :: vec
        type(AtomIdVector_t), intent(in) :: path
        !! END of arguments
        type(AtomIdVector_t) :: pathReverted
        type(AtomIdVector_t) :: p1
        integer :: i
        integer :: nLen
!        integer :: nPathLen1, nPathLen2
        integer :: removeCnt

        removeCnt = 0
        nLen = PV_getVectorSize(vec)
        call AtomIdVector_Create(p1)

        call Assignment_AIdVector(pathReverted,path)
        call AIdVector_Revert(pathReverted)
        do i=1, nLen
            p1 = vec%data(i)
            if(AIdVector_EqualPath(pathReverted, p1)) then
                call PV_delete_one(vec, i)
                removeCnt = removeCnt + 1
            endif
        enddo

    !        if(.true.) then
    !            print*, "(*) Path Remove Cnt = ", removeCnt
    !        endif
    end subroutine PV_removeRevert

    !!----------------------------------------------------------------------------
    !> Append a value to the vector
    !!----------------------------------------------------------------------------
    subroutine PV_append(vec, indata)
        type(PathVector_t) :: vec
        type(AtomIdVector_t), intent(in) :: indata
        type(AtomIdVector_t), save :: prevSave
        integer :: k
        integer, save :: num_used_prev
        integer :: numUsed

        call PV_checkValid(vec)

        numUsed = PV_getNumUsed(vec)

        if (PV_getNumUsed(vec)+1 .ge. size(vec%data) ) then
            call PV_increase_capacity(vec, numUsed+10)
        endif

        if(numUsed .GE. 1) then
            call Assignment_AIdVector(prevSave, vec%data(numUsed))
            num_used_prev = numUsed
        endif

        call PV_setNumUsed(vec, numUsed + 1)
        k = PV_getNumUsed(vec)

        call Assignment_AIdVector(vec%data(k), indata)

!        if(xlogp_PathVector_getNumUsed(vec) .GE. 47) then
!            print*, 'prev, num_used = ', num_used_prev
!            call prevSave%showPath()
!            print*, 'curr, num_used = ', k
!            call vec%data(k)%showPath()
!            print*, '=========================================='
!            print*, ''
!        endif

    end subroutine PV_append

    subroutine PV_append_unique(vec, indata)
        type(PathVector_t) :: vec
        type(AtomIdVector_t), intent(in) :: indata
!        type(AtomIdVector_t) :: temp
        type(AtomIdVector_t), save :: prevSave
        type(AtomIdVector_t) :: tempRevert
        integer :: k
!        integer, save :: num_used_prev
!        integer :: numPath
!        integer :: numUsed

        call PV_checkValid(vec)


        if (PV_getNumUsed(vec)+1 .ge. size(vec%data) ) then
            call PV_increase_capacity(vec, PV_getNumUsed(vec)+10)
        endif

        call Assignment_AIdVector(tempRevert, indata)
        call AIdVector_Revert(tempRevert)

        if(PV_hasElement(vec, indata)) then
            !print*, 'Aleardy existed data(indata). Return without Append'
            !call indata%showPath()
            return
        endif

        if(PV_hasElement(vec, tempRevert)) then
            !print*, 'Aleardy existed data(indata revert). Return without Append'
            !call tempRevert%showPath()
            return
        endif

        if(vec%num_used .GE. 1) then
            call Assignment_AIdVector(prevSave, vec%data(PV_getNumUsed(vec)))
        endif

        call PV_setNumUsed(vec, PV_getNumUsed(vec) + 1)

        k = PV_getNumUsed(vec)

        call Assignment_AIdVector(vec%data(k), indata)

    end subroutine PV_append_unique

    !!----------------------------------------------------------------------------
    !> Append all the elements of the AtomVector to the vector
    !!----------------------------------------------------------------------------
    subroutine PV_appendAll(vec, other)
        type(PathVector_t) :: vec
        type(PathVector_t) :: other
        !! END of arguments
        type(AtomIdVector_t), pointer :: a1
        integer :: i, nLenAppend

        call PV_checkValid(vec)

        nLenAppend = PV_getVectorSize(other)

        if (PV_getNumUsed(vec) + nLenAppend.ge. size(vec%data) ) then
            call PV_increase_capacity(vec, vec%num_used+nLenAppend)
        endif

        !vec%num_used = vec%num_used + 1
        do i=1, nLenAppend
            !call PV_GetPathPtr(other, i, a1)
            call PV_CopyPath(other, i, a1)
            call PV_append(vec, a1)
        enddo
    end subroutine PV_appendAll

    !!----------------------------------------------------------------------------
    !> Put a value at a specific element of the vector (it needs not yet exist)
    !! Arguments:
    !!  vec Vector in question
    !!  n   Index of the element
    !!  data Data to be put in the vector
    !!----------------------------------------------------------------------------
    subroutine PV_put(vec, n, indata)
        type(PathVector_t), intent(inout) :: vec
        integer, intent(in) :: n
        type(AtomIdVector_t), intent(in) :: indata

        call PV_checkValid(vec)

        if (n .lt. 1) return

        if (n .gt. size(vec%data)) then
            call PV_increase_capacity(vec, n)
        endif

        call PV_setNumUsed(vec, max(PV_getNumUsed(vec), n))
        vec%data(n) = indata

    end subroutine PV_put

    !!----------------------------------------------------------------------------
    !> Expand the array holding the data
    !!
    !! Arguments:
    !!  vec PathVector_t in question
    !!  capacity Minimum capacity
    !!----------------------------------------------------------------------------
    subroutine PV_increase_capacity(vec, capacity)
        type(PathVector_t), intent(inout) :: vec
        integer, intent(in) :: capacity
        integer :: new_cap
        type(AtomIdVector_t), dimension(:), pointer :: new_data
        type(AtomIdVector_t) empty_path
        integer :: i
        type(AtomId_t) :: atom1

        call PV_checkValid(vec)

        call AtomIdVector_Create(empty_path, 10)
        call xlogp_AtomId_Clear(atom1)
        do i=1, AIdVector_GetVectorSize(empty_path)
            call AIdVector_Append(empty_path, atom1)
        enddo

        if(capacity .LT. XSCORE_MAX_PATHVECTOR_SIZE) then
            return
        endif

        new_cap = max( capacity, nint(pathvector_growth_rate * size(vec%data)))

        if(new_cap .gt. XSCORE_MAX_PATHVECTOR_SIZE) then
            print*, '(*) New capacity required is greater than XSCORE_MAX_PATHVECTOR_SIZE'
            print*, '(*) new capacity size required = ', new_cap
            print*, '(*) XSCORE_MAX_PATHVECTOR_SIZE = ', XSCORE_MAX_PATHVECTOR_SIZE
            stop
        endif

!        if (new_cap .gt. size(vec%data)) then
!            allocate(new_data(1:new_cap))
!            do i=1, new_cap
!                call AtomIdVector_Create(new_data(i))
!            enddo
!
!            !new_data(1:vec%num_used) = vec%data(1:vec%num_used)
!            do i=1, vec%num_used
!                call Assignment_AIdVector(new_data(i), vec%data(i))
!            enddo
!
!            do i=PV_getNumUsed(vec)+1,new_cap
!                call Assignment_AIdVector(new_data(i), empty_path)
!            enddo
!            !new_data(vec%num_used+1:new_cap) = empty_path
!
!            if(associated(vec%data)) then
!                deallocate(vec%data)
!            endif
!            vec%data => new_data
!        endif

        call AtomIdVector_Destroy(empty_path)

    end subroutine PV_increase_capacity

    !!----------------------------------------------------------------------------
    !> show all the atoms
    !!----------------------------------------------------------------------------
    subroutine PV_showAll(vec)
        type(PathVector_t) :: vec
        integer i, nSize
!        type(AtomIdVector_t) :: a1

        nSize = PV_getVectorSize(vec)
        do i=1, nSize
            call AIdVector_PathShow(vec%data(i), i)
        end do

    end subroutine PV_showAll

    subroutine PV_showMember(vec, number)
        type(PathVector_t) :: vec
        integer, intent(in) :: number

        if(number .LE. 0 .OR. number .GT. PV_getNumUsed(vec)) then
            print*, '(*) ERROR : Out of boundary at xlogp_PathVector_showMember()'
            return
        endif

        call AIdVector_PathShow(vec%data(number))

    end subroutine PV_showMember

    function PV_hasElement(vec, path) result(hasFlag)
        type(PathVector_t) :: vec
        integer i, nSize
        type(AtomIdVector_t) :: path
        logical :: hasFlag

        hasFlag = .false.
        nSize = PV_getVectorSize(vec)
        do i=1, nSize
            if(Is_equal_atomIdvector_to(vec%data(i),path)) then
                hasFlag = .true.
                EXIT
            endif
        enddo
    end function PV_hasElement

    function PV_hasElemIncludingRevert(vec, path) result(hasFlag)
        type(PathVector_t), intent(in) :: vec
        type(AtomIdVector_t), intent(in) :: path
        type(AtomIdVector_t) :: pathRevert
        integer i, nSize
        logical :: hasFlag

        hasFlag = .false.
        nSize = PV_getVectorSize(vec)
        do i=1, nSize
            if(Is_equal_atomIdvector_to(vec%data(i),path)) then
                hasFlag = .true.
                EXIT
            endif
        enddo

        call Assignment_AIdVector(pathRevert, path)
        do i=1, nSize
            if(Is_equal_atomIdvector_to(vec%data(i),pathRevert)) then
                hasFlag = .true.
                EXIT
            endif
        enddo

    end function PV_hasElemIncludingRevert

    subroutine PV_increaseNumUsed(vec)
        type(PathVector_t), intent(inout) :: vec
        call PV_setNumUsed(vec, PV_getNumUsed(vec) + 1)
    end subroutine PV_increaseNumUsed

    subroutine PV_initNumUsed(vec)
        type(PathVector_t), intent(inout) :: vec
        vec%num_used = 0
    end subroutine PV_initNumUsed

    subroutine PV_setNumUsed(vec, num)
        type(PathVector_t), intent(inout) :: vec
        integer, intent(in) :: num
        vec%num_used = num
    end subroutine PV_setNumUsed

    integer function PV_getNumUsed(vec)
        type(PathVector_t), intent(in) :: vec
        PV_getNumUsed = vec%num_used
    end function PV_getNumUsed

end module xlogp_PathVector_m
