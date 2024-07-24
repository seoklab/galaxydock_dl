module xlogp_AtomIdVector_m
    use xlogp_Atom_m, only : AtomId_t, is_equal_atomid_to, &
                    xlogp_AtomId_Clear, xlogp_AtomId_Show
    use xlogp_constants_m
    use xlogp_Bond_m
    implicit none
    !private ! Information hiding

    public :: Assignment_AIdVector
    public :: Is_equal_AtomIdvector_to
    public :: assignment (=)
    public :: operator (==)

    logical, private, parameter :: DEBUG = .false.
    logical, private, parameter :: DEBUG2 = .false.
    logical, private, parameter :: DEBUG3 = .false.

    real, private, parameter :: growth_rate = 1.5

    !!----------------------------------------------------------------------------
    ! Type definition
    type AtomIdVector_t
        !private ! Component hiding
        character(len=100) :: name != 'temp vector'
        integer :: pathId != 0
        integer :: num_used != 0
        !
        !type(AtomId_t), dimension(:), pointer :: data => null()
        !
        ! NOTE : I change allocatable dynamic AtomId_t array into fixed array.
        !           And, We assume that the size of this array corresponds to the number
        !           maximum ligand atoms. ==> 200 (XSCORE_MAX_LIGAND_ATOM_NUMBER)
        !
        type(AtomId_t), dimension(XSCORE_MAX_LIGAND_ATOM_NUMBER) :: data
        logical :: validity != .false.
    end type AtomIdVector_t

    interface AtomIdVector
        module procedure Construct_AtomIdVector
    end interface !AtomIdVector

    interface assignment (=) ! overload =
        module procedure Assignment_AIdVector
    end interface

    interface operator (==) ! overload ==
        module procedure Is_equal_AtomIdvector_to
    end interface

contains

    function Construct_AtomIdVector() result(return_value)
        type(AtomIdVector_t) :: return_value
        call AtomIdVector_Create(return_value)
    end function Construct_AtomIdVector

    subroutine AtomIdVector_Create(this, capacity)
        type(AtomIdVector_t) :: this
        integer, optional :: capacity
        integer :: i
!        type(AtomId_t) :: atom_dummy

        ! Check that the vector does not have any data left
!        if (associated(this%data)) then
!            call this%Destroy()
!        end if

!        if (present(capacity)) then
!            cap = max(1, capacity)
!        else
!            cap = 10
!        endif
!
!        allocate(this%data(1:cap))
!        call atom_dummy%Clear()

        do i=1, XSCORE_MAX_LIGAND_ATOM_NUMBER
            call xlogp_AtomId_Clear(this%data(i))
        enddo

        this%num_used = 0

        this%validity = .true.

    end subroutine AtomIdVector_Create

    subroutine AtomIdVector_Destroy(this)
        type(AtomIdVector_t) :: this
!        integer :: ok

        ! Check that the vector does not have any data left
!        if (associated(this%data) .AND. this%validity) then
!            deallocate(this%data, stat=ok)
!            if(ok /= 0) then
!                print*, '(*) ERROR : deallocation in Destroy()'
!            endif
!        endif

        this%num_used = 0
!        this%data => null()
    end subroutine AtomIdVector_Destroy

    subroutine Assignment_AIdVector(vleft, vright)
        type(AtomIdVector_t), intent(inout) :: vleft
        type(AtomIdVector_t), intent(in) :: vright
        integer :: vrightSize
        integer :: i
        !> END of local variabels

        vrightSize = AIdVector_GetCapacitySize(vright)
        !if(DEBUG2) print*, '>>> assignment called, vrightSize = ', vrightSize

!        if(.not. associated(vright%data)) then
!            !print*, '(*) ERROR: right handed value is not associated'
!        else if (vright%num_used == 0) then
!            !print*, '(*) WARNING : right handed object has zero used member'
!        endif

!        if(.not. associated(vleft%data)) then
!            if(DEBUG) print*, '(*) INFO : Left data reference is null(). Creating left-haded vector ...'
!            call vleft%Create(vrightSize)
!        else if(associated(vleft%data) .AND. &
!                & vleft%GetCapacitySize() < vright%GetCapacitySize()) then
!            if(DEBUG) print*, '(*) INFO : Left and right size is different. Left = ', &
!                                & vleft%GetCapacitySize(), &
!                                ', Right = ', vright%GetCapacitySize()
!            call vleft%Destroy()
!            call vleft%Create(vrightSize)
!        endif
        !> copy ritht object values into left object
!        vleft%data(1:vright%num_used) = vright%data(1:vright%num_used)
        do i=1, vright%num_used
            vleft%data(i) = vright%data(i)
        enddo
        !> change num_used value.
        vleft%num_used = vright%num_used
        vleft%pathId = vright%pathId
    end subroutine Assignment_AIdVector

    logical function Is_equal_AtomIdvector_to(a, b) result(t_f)   ! overload ==
        type(AtomIdVector_t), intent(in) :: a, b ! left & right of ==
        integer :: i
        t_f = .true. ! initialize

        if(AIdVector_GetVectorSize(a) /= AIdVector_GetVectorSize(b)) then
            if(DEBUG2) then
                print*, 'avector size = ', AIdVector_GetVectorSize(a), &
                        ' bvector size = ', AIdVector_GetVectorSize(b)
            endif
            t_f = .false.
            return ! same size?
        endif

        !        t_f = all (a%data == b%data)  ! and all values(AtomId_t) match
        do i=1, AIdVector_GetVectorSize(a)
            if(Is_equal_atomid_to(a%data(i),b%data(i))) then
            else
                t_f = .false.
                if(DEBUG3) then
                    print*, i, '-th elements are different'
                    call xlogp_AtomId_show(a%data(i))
                    call xlogp_AtomId_show(b%data(i))
                endif
                EXIT
            endif
        enddo
    end function Is_equal_AtomIdvector_to

    !!----------------------------------------------------------------------------
    ! Return the number of elements in use
    function AIdVector_GetVectorSize(this) result(return_value)
        type(AtomIdVector_t) :: this
        integer :: return_value

        return_value = this%num_used
    end function AIdVector_GetVectorSize

    !!----------------------------------------------------------------------------
    ! Return the total number of elements which was allowed to be inserted
    function AIdVector_GetCapacitySize(this) result(return_value)
        type(AtomIdVector_t) :: this
        integer :: return_value

        return_value = size(this%data)
    end function AIdVector_GetCapacitySize

    !!----------------------------------------------------------------------------
    !>
    !!----------------------------------------------------------------------------

    function AIdVector_Get(this, n) result(return_value)
        type(AtomIdVector_t), intent(inout) :: this
        integer, intent(in) :: n
        type(AtomId_t) :: return_value
        type(AtomId_t) :: empty_atom

        if( n .lt. 1 .OR. n .gt. this%num_used) then
            print*, "(*) ERROR : Out of boundary, n = ", n, &
                ", num_used = ", this%num_used
            call xlogp_AtomId_clear(empty_atom)
            return_value = empty_atom
        else
            return_value = this%data(n)
        end if

    end function AIdVector_Get

!    subroutine GetAtomIdPtr(this, n, result_ptr)
!        type(AtomIdVector_t), intent(in) :: this
!        integer, intent(in) :: n
!        type(AtomId_t), intent(out), pointer :: result_ptr
!        type(AtomId_t) :: empty_atom
!
!        if( n .lt. 1 .OR. n .gt. this%num_used) then
!            print*, "(*) ERROR : Out of boundary, n = ", n, &
!                ", num_used = ", this%num_used
!            call empty_atom%clear()
!            result_ptr => null()
!        else
!            result_ptr => this%data(n)
!        end if
!
!    end subroutine GetAtomIdPtr


    function AIdVector_AtomId_Value(this, n) result(return_value)
        type(AtomIdVector_t), intent(in) :: this
        integer, intent(in) :: n
        integer :: return_value
        type(AtomId_t) :: empty_atom

        if( n .lt. 1 .OR. n .gt. this%num_used) then
            print*, "(*) ERROR : Out of boundary at GetAtomId_IdValue(), n = ", n, &
                ", num_used = ", this%num_used
            call xlogp_AtomId_clear(empty_atom)
            return_value = 0
        else
            return_value = this%data(n)%id
        end if

    end function AIdVector_AtomId_Value

    function AIdVector_GetPathId(this) result(return_value)
        type(AtomIdVector_t), intent(in) :: this
        integer :: return_value
        return_value = this%pathId
    end function AIdVector_GetPathId

    subroutine AIdVector_SetPathId(this, pathId)
        type(AtomIdVector_t), intent(inout) :: this
        integer, intent(in) :: pathId
        this%pathId = pathId
    end subroutine AIdVector_SetPathId
    !!----------------------------------------------------------------------------
    !> Insert one or more empty elements and move original data backward.
    !!
    !! Arguments:
    !!  this AtomIdVector_t in question
    !!  pos Position to insert the empty elements
    !!  emptyNum Number of empty elements
    !!----------------------------------------------------------------------------
    subroutine AIdVector_Insert_empty(this, pos, emptyNum)
        type(AtomIdVector_t) :: this
        integer, intent(in) :: pos
        integer, intent(in) :: emptyNum

        integer :: i
        type(AtomId_t) :: empty_atom

        ! Clear one atom
        call xlogp_AtomId_clear(empty_atom)

        if (emptyNum .lt. 1 .OR. &
            pos .lt. 1 .OR. &
            pos .gt. this%num_used ) then
            return
        endif

        if (this%num_used+emptyNum .ge. size(this%data) ) then
            call AIdVector_Increase_capacity(this, this%num_used + emptyNum)
        endif

        ! From pos postion, insert emptyNum number of empty atoms
        ! E.g) for [a1,a2,a3,a4,a5,,],
        ! num_used=5
        ! insert_empty(3,2) =>
        ! a5 -> 5+2
        ! a4 -> 5+1
        !
        do i = this%num_used, pos, -1
            this%data(i+emptyNum) = this%data(i)
        enddo

        ! Here, we made some dummy atoms as place holder.
        do i = 1, emptyNum
            this%data(pos+i-1) = empty_atom
        enddo

        this%num_used = this%num_used + emptyNum

    end subroutine AIdVector_Insert_empty

    subroutine AIdVector_Delete_one_element(this, pos)
        type(AtomIdVector_t) :: this
        integer, intent(in) :: pos

        integer :: i

        if (pos .lt. 1 .OR. &
            pos .gt. this%num_used ) then
            print*, '(*) Out of boundary ERROR at Delete_one_element(). pos = ', pos
            return
        endif

        if (pos == this%num_used) then
            call xlogp_AtomId_Clear(this%data(pos))
            this%num_used = this%num_used - 1
            return
        endif

        do i=pos,this%num_used
            this%data(i) = this%data(i+1)
        enddo

        this%num_used = this%num_used - 1
    end subroutine AIdVector_Delete_one_element

    subroutine AIdVector_Delete_elements(this, pos, numDelete)
        type(AtomIdVector_t) :: this
        integer, intent(in) :: pos
        integer, intent(in) :: numDelete

        integer :: i

        if (numDelete .lt. 1 .OR. &
            pos .lt. 1 .OR. &
            pos .gt. this%num_used ) then
            return
        endif

        do i=pos,this%num_used-numDelete
            this%data(i) = this%data(i+numDelete)
            !call this%data(i)%clear()
        enddo

        this%num_used = this%num_used - numDelete
    end subroutine AIdVector_Delete_elements

    subroutine AIdVector_Clear(this)
        type(AtomIdVector_t) :: this

        call AIdVector_Delete_elements(this, 1, this%num_used)
        this%num_used = 0
    end subroutine AIdVector_Clear
    !!----------------------------------------------------------------------------
    !> Remove only one atom which corresponds to the given atom in the vector.
    !!----------------------------------------------------------------------------
    subroutine AIdVector_Remove(this, atom)
        type(AtomIdVector_t) :: this
        type(AtomId_t) :: atom
        !! END of arguments
        type(AtomId_t) :: a1
        integer :: i
        integer :: nLen
        nLen = AIdVector_GetVectorSize(this)

        do i=1, nLen
            a1 = this%data(i)
            if(a1%id == atom%id) then
                call AIdVector_Delete_elements(this, i, 1)
                EXIT
            end if
        enddo
    end subroutine AIdVector_Remove

    !!----------------------------------------------------------------------------
    !> Append a value to the vector
    !!----------------------------------------------------------------------------
    subroutine AIdVector_Append(this, atom)
        type(AtomIdVector_t), intent(inout) :: this
        type(AtomId_t), intent(in) :: atom

!        if(.not. this%validity) then
!            !           print*, '>>>> Since AtomIdVector is not valid, we will call create()'
!            call this%Create()
!        endif

!        if (this%num_used .ge. size(this%data) ) then
!            call this%Increase_capacity(this%num_used+1)
!            if(DEBUG) then
!            !                print*, 'current num_used = ', (this%num_used + 1)
!            endif
!        endif

        this%num_used = this%num_used + 1
        if(DEBUG) then
        !            print*, 'Appending new atom whose id is ', data%id
        endif
        this%data(this%num_used) = atom
    end subroutine AIdVector_Append

    !!----------------------------------------------------------------------------
    !> Append all the elements of the AtomIdVector to the vector
    !!----------------------------------------------------------------------------
    subroutine AIdVector_AppendAll(this, other)
        type(AtomIdVector_t) :: this
        type(AtomIdVector_t) :: other
        !! END of arguments
        type(AtomId_t) :: a1
        integer :: i, nLenAppend

        nLenAppend = AIdVector_GetVectorSize(other)

        if (this%num_used + nLenAppend.ge. size(this%data) ) then
            call AIdVector_Increase_capacity(this, this%num_used+nLenAppend)
        endif

        !this%num_used = this%num_used + 1
        do i=1, nLenAppend
            a1 = other%data(i)
            call AIdVector_Append(this, a1)
        enddo
    end subroutine AIdVector_AppendAll

    !!----------------------------------------------------------------------------
    !> Put a value at a specific element of the vector (it needs not yet exist)
    !! Arguments:
    !!  this Vector in question
    !!  n   Index of the element
    !!  data Data to be Put in the vector
    !!----------------------------------------------------------------------------
    subroutine AIdVector_Put(this, n, data)
        type(AtomIdVector_t), intent(inout) :: this
        integer, intent(in) :: n
        type(AtomId_t), intent(in) :: data

        if (n .lt. 1) return

        if (n .gt. size(this%data)) then
            call AIdVector_Increase_capacity(this, n)
        endif

        this%num_used = max(this%num_used, n)
        this%data(n) = data

    end subroutine AIdVector_Put

    !!----------------------------------------------------------------------------
    !> Expand the array holding the data
    !!
    !! Arguments:
    !!  this AtomIdVector_t in question
    !!  capacity Minimum capacity
    !!----------------------------------------------------------------------------
    subroutine AIdVector_Increase_capacity(this, capacity)
        type(AtomIdVector_t), intent(inout) :: this
        integer, intent(in) :: capacity
        !!------------ END of inPut args
        integer :: new_cap
!        type(AtomId_t), dimension(:), pointer :: new_data
        type(AtomId_t) empty_atom
        !!------------ END of local variabels

        call xlogp_AtomId_clear(empty_atom)
        new_cap = max( capacity, nint(growth_rate * size(this%data)))

        if (new_cap .gt. size(this%data)) then
            print*, '(*) AtomIdVector_t used fixed data size = ', XSCORE_MAX_LIGAND_ATOM_NUMBER
            return
!            allocate(new_data(1:new_cap))
!
!            new_data(1:this%num_used) = this%data(1:this%num_used)
!            new_data(this%num_used+1:new_cap) = empty_atom
!
!            if(associated(this%data)) then
!                deallocate(this%data)
!            endif
!            this%data => new_data
        endif

    end subroutine AIdVector_Increase_capacity

    !!----------------------------------------------------------------------------
    !> show all the atoms
    !!----------------------------------------------------------------------------
    subroutine AIdVector_ShowAll(this)
        type(AtomIdVector_t) :: this
        integer i, nSize
!        type(AtomId_t), pointer :: a1

        nSize = AIdVector_GetVectorSize(this)
        do i=1, nSize
            call xlogp_AtomId_show(this%data(i))
        end do

    end subroutine AIdVector_ShowAll

    !!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    ! Path Related subprocedures

    !!----------------------------------------------------------------------------
    !> Constructs a new Path with two Atoms
    !! Arguments:
    !!  p Path object which will be passed automatically.
    !!  atom1 first atom in the new path
    !!  atom2 second atom in the new path
    !!----------------------------------------------------------------------------
    subroutine AIdVector_ConstructPath(p, atom1, atom2)
        type(AtomIdVector_t) :: p
        type(AtomId_t), intent(in) :: atom1
        type(AtomId_t), intent(in) :: atom2

        call AIdVector_Append(p, atom1)
        call AIdVector_Append(p, atom2)
    end subroutine AIdVector_ConstructPath

    subroutine AIdVector_PathShow(p, rowId)
        implicit none
        type(AtomIdVector_t), intent(in) :: p
        integer, intent(in), optional :: rowId
        integer :: i
        integer :: nMember
!        type(AtomId_t) :: a1
!        character(len=100) :: tmp
!        character(len=500) :: longTmp
        integer :: rowIdNum

        if(present(rowId)) then
            rowIdNum = rowId
        else
            rowIdNum = 0
        endif
        nMember = AIdVector_GetVectorSize(p)
        write(*, FMT="(I3, 1X, I3, 1X, A, I3, A)", ADVANCE="NO") rowIdNum, p%pathId, &
           " [ LEN(Path) = ", nMember, " ] : "
        do i=1, nMember
            write(*, FMT="(I3, 2X)", ADVANCE="NO") p%data(i)%id
        enddo
        write(*,*) ''

    end subroutine AIdVector_PathShow

    subroutine AIdVector_PathShowAsType(p)
        implicit none
        type(AtomIdVector_t), intent(in) :: p
        integer :: i
        integer :: nMember
!        type(AtomId_t) :: a1
!        character(len=100) :: tmp
!        character(len=500) :: longTmp

        nMember = AIdVector_GetVectorSize(p)
        write(*, FMT="(A, I3, A)", ADVANCE="NO") ">>> LEN(Path) = ", nMember, ", "
        do i=1, nMember
            write(*, FMT="(A, 2X)", ADVANCE="NO") trim(p%data(i)%type)
        enddo
        write(*,*) ''

    end subroutine AIdVector_PathShowAsType

    subroutine AIdVector_PathShowAsName(p)
        implicit none
        type(AtomIdVector_t), intent(in) :: p
        integer :: i
        integer :: nMember
!        type(AtomId_t) :: a1
!        character(len=100) :: tmp
!        character(len=500) :: longTmp

        nMember = AIdVector_GetVectorSize(p)
        write(*, FMT="(A, I3, A)", ADVANCE="NO") ">>> LEN(Path) = ", nMember, ", "
        do i=1, nMember
            write(*, FMT="(A, 2X)", ADVANCE="NO") trim(p%data(i)%name)
        enddo
        write(*,*) ''

    end subroutine AIdVector_PathShowAsName


    subroutine AIdVector_Revert(p)
        implicit none
        type(AtomIdVector_t), intent(inout) :: p
        integer :: nSize
        integer :: i, j, mid
        type(AtomId_t) :: a1
        type(AtomId_t) :: a2

        nSize = AIdVector_GetVectorSize(p)
        mid = nSize / 2
        do i=1, mid
            a1 = p%data(i)
            j = nSize - i + 1
            a2 = p%data(j)
            call AIdVector_Put(p, i, a2)
            call AIdVector_Put(p, j, a1)
        enddo

    end subroutine AIdVector_Revert

    function AIdVector_GetIntSectSize(p, other) result(nIntersection)
        type(AtomIdVector_t), intent(in) :: p
        type(AtomIdVector_t), intent(in) :: other
        type(AtomId_t) :: a1, a2
        integer :: nSize1, nSize2, i, j
        integer :: nIntersection

        nSize1 = AIdVector_GetVectorSize(p)
        nSize2 = AIdVector_GetVectorSize(other)
        nIntersection = 0

        do i=1, nSize1
            a1 = p%data(i)
            do j=1, nSize2
                a2 = other%data(j)
                if(a1%id == a2%id) nIntersection = nIntersection + 1
            enddo
        enddo
    end function AIdVector_GetIntSectSize

    subroutine AIdVector_Join(this, other, atom, newPath)
        type(AtomIdVector_t), intent(in) :: this
        type(AtomIdVector_t), intent(in) :: other
        type(AtomId_t), intent(in) :: atom
        type(AtomIdVector_t), intent(out) :: newPath
        type(AtomIdVector_t) tempPath
        integer :: aid1, aid2
        integer :: i, numSum
        !! END of local variables

        call Assignment_AIdVector(newPath, this)
        call Assignment_AIdVector(tempPath, other)

!        print*, "Chk 1"
!        call newPath%ShowPath()
!        print*, "Chk 2"
!        call tempPath%ShowPath()

        ! Make sure that newPath has atom%id at the last position
        aid1 = AIdVector_AtomId_Value(newPath, 1)
        aid2 = AIdVector_AtomId_Value(newPath, &
                AIdVector_GetVectorSize(newPath))
        if (aid1 == atom%id) then
            call AIdVector_Revert(newPath)
        else if(aid2 == atom%id) then
            ! No change
        else
            print*, '(*) ERROR at Join() : No overlapping atom id between first path and atom'
            call AIdVector_PathShow(newPath)
            print*, '(*) atom%id = ', atom%id
        endif

!        print*, 'Chk II '
!        call newPath%ShowPath()

        ! Make sure that second path has atom%id at the first position
        aid1 = AIdVector_AtomId_Value(tempPath, 1)
        aid2 = AIdVector_AtomId_Value(tempPath, &
                AIdVector_GetVectorSize(tempPath))
        if (aid1 == atom%id) then
            ! No change
        else if (aid2 == atom%id) then
            call AIdVector_Revert(tempPath)
        else
            print*, '(*) ERROR at Join() : No overlapping atom id between second path and atom'
            call AIdVector_PathShow(tempPath)
            print*, '(*) atom%id = ', atom%id
        endif

!        print*, 'Chk kk '
!        call tempPath%ShowPath()

        ! Remove overlapping atom element which locates at the first.
        call AIdVector_Delete_one_element(tempPath, 1)

!        print*, 'Chk after delete_one'
!        call tempPath%ShowPath()

        !! Join two pathes
        numSum = AIdVector_GetVectorSize(newPath) + AIdVector_GetVectorSize(tempPath)
!        if(newPath%GetCapacitySize() .LT. numSum) then
!            call newPath%Increase_capacity(numSum+1)
!        endif

        do i=1, AIdVector_GetVectorSize(tempPath)
            call AIdVector_Append(newPath, AIdVector_Get(tempPath, i))
        enddo

!        call newPath%ShowPath()
!        print*, 'Chk end'

    end subroutine AIdVector_Join

    function AIdVector_DetermineRing(p) result(ringFlag)
        type(AtomIdVector_t) :: p
        logical :: ringFlag
        integer :: j, nMemberSize
        type(AtomId_t) :: a1, a2

        ringFlag = .false.
        nMemberSize = AIdVector_GetVectorSize(p) ! member size
        a1 = p%data(1) ! get first element
        a2 = p%data(nMemberSize) ! get last element

        if(nMemberSize > 3 .AND. (a1%id == a2%id)) then !! Check if this path is a correct ring.
            if(DEBUG .eqv. .TRUE.) then
                print*, "New ring found"
                call AIdVector_PathShow(p)
            endif
            ringFlag = .true.
        endif
    end function AIdVector_DetermineRing

    function AIdVector_EqualPath(p, other) result(samePath)
        type(AtomIdVector_t) :: p
        type(AtomIdVector_t) :: other
        logical :: samePath
        integer :: i, nLen1, nLen2
        type(AtomId_t) :: a1, a2

        nLen1 = AIdVector_GetVectorSize(p)
        nLen2 = AIdVector_GetVectorSize(other)

        samePath = .true.

        if (nLen1 /= nLen2) then
            samePath = .false.
        endif

        if (nLen1 == nLen2) then
            do i=1, nLen1
                a1 = p%data(i)
                a2 = other%data(i)
                if(a1%id /= a2%id) then
                    samePath = .false.
                    EXIT
                endif
            enddo
        endif

    end function AIdVector_EqualPath

    function AIdVector_IsEdgeElement(p, atomId) result(edgeFlag)
        type(AtomIdVector_t) :: p
        type(AtomId_t) :: atomId
        logical :: edgeFlag
        integer :: nLen

        nLen = AIdVector_GetVectorSize(p)
        edgeFlag = (p%data(1)%id == atomId%id) .OR. &
                    (p%data(nLen)%id == atomId%id)
    end function AIdVector_IsEdgeElement

end module xlogp_AtomIdVector_m
