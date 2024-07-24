module xlogp_AllRingsFinder_m
    use xlogp_constants_m, only : XSCORE_MAX_ATOM_NEIB, XSCORE_MAX_RING_NUMBER
    use xlogp_Atom_m
    use xlogp_Bond_m
    use xlogp_Group_m
    use xlogp_Ring_m
    use xlogp_AtomIdVector_m
    use xlogp_PathVector_m
    implicit none

    integer, parameter :: ringsearch_timeout = 5000
    logical, private, parameter :: DEBUG = .false.

    type AllRingsFinder_t
        logical :: debug != .false.
        integer :: startTime != 0
        integer :: maxPathLen != 50
        integer :: totIntersectPath != 0
        integer :: totAddedPathCnt != 0
        integer :: totInvolvedPathCntPerAtom != 0
        integer, dimension(:), pointer :: markForPathsWithAtomX !=> null() ! Convenience array for marking
        type(AtomId_t) :: removeAtom !> To communicate with showDetectProcess procedure.

        type(AtomIdVector_t) :: atomList  !> This will store all atoms given from caller
        type(AtomIdVector_t) :: originalAtomList  !> This will store all atoms given from caller
        type(Bond_t), DIMENSION(XSCORE_MAX_LIGAND_BOND_NUMBER) :: bondList !=> null() !> This will store the pointer of bonds
        integer :: bondListSize ! This was added due to that bondList was changed to fixed array
                                ! from dynamic array

        type(PathVector_t) :: pathGraph
        type(PathVector_t) :: newPaths
        type(PathVector_t) :: potentialRings
        type(PathVector_t) :: allRings
        type(PathVector_t) :: removePaths
!        type(Ring_t), DIMENSION(XSCORE_MAX_RING_NUMBER) :: foundRings
!        integer :: foundRingSize

    end type AllRingsFinder_t

    interface AllRingsFinder
        module procedure newAllRingsFinder
    end interface !AllRingsFinder

contains

    !!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    ! Constructor and Destructor

    function newAllRingsFinder(an_atomList, num_atom, bondList, num_bond) result(return_value)
        type(Atom_t), DIMENSION(:), target, intent(in) :: an_atomList
        integer, intent(in) :: num_atom
        type(Bond_t), DIMENSION(:), intent(in) :: bondList
        integer, intent(in) :: num_bond
        type(AllRingsFinder_t) :: return_value

        call xlogp_arf_initialize(return_value, an_atomList, num_atom, bondList, num_bond)
    end function newAllRingsFinder

    subroutine xlogp_arf_initialize(this, an_atomList, num_atom, bondList, num_bond)
        type(AllRingsFinder_t) :: this
        type(Atom_t), DIMENSION(:), target, intent(in) :: an_atomList
        integer, intent(in) :: num_atom
        type(Bond_t), DIMENSION(:), intent(in) :: bondList
        integer, intent(in) :: num_bond

        !> END of arguments
        integer :: i, nBonds
        type(AtomId_t) :: atomId1
        integer :: alloc_stat
        !> END of local variables

        call AtomIdVector_Create(this%atomList) !! => atomList
        call AtomIdVector_Create(this%originalAtomList) !! => atomList

        do i=1, num_atom
            !> We need AtomId_t from Atom_t to make a path graph
            atomId1 = makeAtomId(an_atomList(i))
            call AIdVector_Append(this%atomList, atomId1)
            call AIdVector_Append(this%originalAtomList, atomId1)
        enddo

        nBonds = num_bond ! size(bondList)
        this%bondListSize = nBonds
        if(nBonds > XSCORE_MAX_LIGAND_BOND_NUMBER) then
            print*, '(*) ERROR bond list size out of boundary : max value = ', &
                XSCORE_MAX_LIGAND_BOND_NUMBER, ', given bond list = ', nBonds
            stop
        endif

!        if(.not. associated(this%bondList)) then
!            allocate(this%bondList(nBonds))
!        endif

        do i=1, nBonds
            this%bondList(i) = bondList(i)
        enddo

        call PV_create(this%pathGraph)
        call PV_create(this%newPaths)
        call PV_create(this%potentialRings)
        call PV_create(this%removePaths)
        call PV_create(this%allRings)

!        print*, 'nBonds = ', nBonds

!        call PV_increase_capacity(this%pathGraph, nBonds)
!        call PV_increase_capacity(this%newPaths, nBonds)
!        call PV_increase_capacity(this%potentialRings, nBonds)
!        call PV_increase_capacity(this%removePaths, nBonds)

!        print*, 'nBonds = ', nBonds

        allocate(this%markForPathsWithAtomX(nBonds), stat=alloc_stat)
        if(alloc_stat /= 0) then
            print*, '(*) ALLOC ERROR for markForPathsWithAtomX = ', alloc_stat
        endif
    end subroutine xlogp_arf_initialize

    subroutine xlogp_arf_destroy(this)
        type(AllRingsFinder_t) :: this

        call AtomIdVector_Destroy(this%atomList)
        call AtomIdVector_Destroy(this%originalAtomList)

!        if(associated(this%bondList)) then
!            deallocate(this%bondList)
!        endif

        call PV_destroy(this%pathGraph)
        call PV_destroy(this%newPaths)
        call PV_destroy(this%potentialRings)
        call PV_destroy(this%removePaths)
        call PV_destroy(this%allRings)

        if(associated(this%markForPathsWithAtomX)) then
            deallocate(this%markForPathsWithAtomX)
        endif

    end subroutine xlogp_arf_destroy

    !!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    ! Ring Search Methods

    !!----------------------------------------------------------------------------
    !> This function will search rings in a moleucle
    !! @param maxPathLen
    !! @return allRings (which is a member variable) will gather rings found.
    !!----------------------------------------------------------------------------
    subroutine xlogp_arf_doSearch(this, maxPathLen, debugFlag)
        type(AllRingsFinder_t) :: this
        integer, intent(in) :: maxPathLen
        logical, optional, intent(in) :: debugFlag
!        integer :: i
        type(PathVector_t) :: pathVectorCopy

        if(present(debugFlag)) then
            this%DEBUG = debugFlag
        endif

        this%maxPathLen = maxPathLen

        call PV_create(pathVectorCopy)

        !> First we convert the molecular graph into a path graph using bond list.
        call xlogp_arf_initPathGraphs(this)

        do
            !> Atoms with lower degree will be selected at first.
            this%removeAtom = xlogp_arf_selectAtom(this)

            call xlogp_arf_clearTempPaths(this)

            call arf_prepPathMergeForOneAtom(this, this%removeAtom, &
                & this%pathGraph, this%totIntersectPath)

            if(this%DEBUG) print*, '>>> Before Collapsing '
            if(this%DEBUG) call xlogp_arf_showDetectProcess(this)

            !! Save PathVector for comparison
!            pathVectorCopy = this%pathGraph

            call arf_collapsePathGraphSafe(this, this%pathGraph, maxPathLen)

!            print*, '(##) PathVector Comparison 1 : '
!            call pathVectorCopy%ShowAll()
!            print*, '(##) PathVector Comparison 2 : '
!            call this%pathGraph%ShowAll()

            call xlogp_arf_detectRings(this%potentialRings, this%allRings)
            call removeOneAtomAndPaths(this, this%removeAtom)

            if(this%DEBUG) print*, '>>> After Collapsing '
            if(this%DEBUG) call xlogp_arf_showDetectProcess(this)

            if(PV_getVectorSize(this%pathGraph) <= 0 .OR. &
                AIdVector_GetVectorSize(this%atomList) <= 0) EXIT

            if(this%DEBUG) then
                print*, '!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
            endif
        enddo

    end subroutine xlogp_arf_doSearch

    subroutine xlogp_arf_initPathGraphs(this)
        type(AllRingsFinder_t) :: this
        integer :: i, nBonds, nAtoms
        integer :: atom1, atom2
!        type(Bond_t) :: b1, b2
        type(AtomId_t), allocatable :: a1(:)
        type(AtomId_t), allocatable :: a2(:)
        type(AtomIdVector_t) :: path1

        nBonds = xlogp_arf_getAllBondsNum(this)
        nAtoms = xlogp_arf_getAllAtomsNum(this)

        allocate(a1(nBonds))
        allocate(a2(nBonds))

        call PV_clear(this%pathGraph)

        call PV_increase_capacity(this%pathGraph, nBonds)

        do i=1, nBonds
            atom1 = this%bondList(i)%atom_1
            atom2 = this%bondList(i)%atom_2
            if(atom1 < 1 .OR. atom1 > nAtoms) then
                print*, "(*) ERROR : bond ", i, ", atom1 is abnormal = ", atom1
            endif
            if(atom2 < 1 .OR. atom2 > nAtoms) then
                print*, "(*) ERROR : bond ", i, ", atom2 is abnormal = ", atom2
            endif

            a1(i) = AIdVector_Get(this%atomList, atom1)
            a2(i) = AIdVector_Get(this%atomList, atom2)

            call AIdVector_Clear(path1)
            call AIdVector_SetPathId(path1, i)
            call AIdVector_Append(path1, a1(i))
            call AIdVector_Append(path1, a2(i))

            call PV_append(this%pathGraph, path1)
        enddo

        deallocate(a1)
        deallocate(a2)
    end subroutine xlogp_arf_initPathGraphs

    function xlogp_arf_selectAtom(this) result(minAtom)
        type(AllRingsFinder_t) :: this
        type(AtomId_t) :: minAtom
        integer :: minDegree
        integer :: degree
        type(AtomId_t) :: atom
        integer :: i
        !> END of local variables

        minAtom = AIdVector_Get(this%atomList, 1)
        minDegree = minAtom%num_neib

        do i=2, AIdVector_GetVectorSize(this%atomList)
            atom = AIdVector_Get(this%atomList, i)
            degree = atom%num_neib !! number of connected bonds
            if (degree < minDegree) then
                minAtom = atom
                minDegree = degree
            endif
        enddo
    end function xlogp_arf_selectAtom

    !!----------------------------------------------------------------------------
    !> The conveninece array (markForPathsWithAtomX) should be checked if its
    !! size is sufficient because this%pathGraph may be increased more than
    !! the number of bond list.
    !!----------------------------------------------------------------------------
    subroutine CheckMarkSzForPathsWAtomX(this, nPaths)
        type(AllRingsFinder_t) :: this
        integer, intent(in) :: nPaths
        integer :: alloc_stat

        if(nPaths > size(this%markForPathsWithAtomX)) then
            deallocate(this%markForPathsWithAtomX)
            allocate(this%markForPathsWithAtomX(nPaths+10), stat=alloc_stat)
            if(alloc_stat /= 0) then
                print*, 'ALLOCATION ERROR for this%markForPathsWithAtomX '
                ! @TODO We need some exception handling here.
            endif
        endif
    end subroutine CheckMarkSzForPathsWAtomX

    !!----------------------------------------------------------------------------
    !> Remove an atom from the atom container under certain conditions.
    !! This function will change the member objects of newPaths, potentialRings,
    !!  allRings, and removePaths.
    !!
    !! Actually, pathGraph will be modified during this function.
    !!
    !!----------------------------------------------------------------------------
    subroutine arf_prepPathMergeForOneAtom(this, atom, pathList, totIntersectPath)
        type(AllRingsFinder_t) :: this
        type(AtomId_t), intent(in) :: atom
        type(PathVector_t), intent(inout) :: pathList
        integer, intent(inout) :: totIntersectPath
!        type(AtomIdVector_t), pointer :: path1_ptr, path2_ptr
        type(AtomIdVector_t) :: path1, path2
        type(AtomIdVector_t) :: union !, unionReverted
        integer :: intersectionSize
        integer :: i,j, nPaths
        logical :: firstEdgeElement
        logical :: secondEdgeElement

        nPaths = PV_getVectorSize(pathList)
        totIntersectPath = 0
        this%totInvolvedPathCntPerAtom = 0

        call CheckMarkSzForPathsWAtomX(this, nPaths)

        ! Set flag zero (OFF)
        this%markForPathsWithAtomX = 0
        do i=1, nPaths
            !call PV_GetPathPtr(pathList, i, path1_ptr)
            !path1_ptr => pathList%data(i)

            if(AIdVector_IsEdgeElement(pathList%data(i), atom)) then
                this%markForPathsWithAtomX(i) = 1
                this%totInvolvedPathCntPerAtom = this%totInvolvedPathCntPerAtom + 1
            endif
        enddo

        do i=1, nPaths
            !call PV_GetPathPtr(pathList, i, path1_ptr)
            !path1_ptr => pathList%data(i)

            firstEdgeElement = AIdVector_IsEdgeElement(pathList%data(i), atom)
            if(.not. firstEdgeElement) cycle
            call Assignment_AIdVector(path1, pathList%data(i))
            call PV_append_unique(this%removePaths, path1)
            totIntersectPath = totIntersectPath + 1
            do j=i+1, nPaths
                !call PV_GetPathPtr(pathList, j, path2_ptr)
!                path2_ptr => pathList%data(j)

                secondEdgeElement = AIdVector_IsEdgeElement(pathList%data(j), atom)
                if(.not. secondEdgeElement) cycle
                call Assignment_AIdVector(path2, pathList%data(j))
                totIntersectPath = totIntersectPath + 1
                intersectionSize = AIdVector_GetIntSectSize(path1, path2)
                if (intersectionSize < 3) then
                    !! atom is a link point of two paths.
                    call AIdVector_Join(path1, path2, atom, union)
                    !print*, '(*) UNION PATH ==> '
                    !call union%ShowPath()
                    if(intersectionSize == 1) then
                        call PV_append_unique(this%newPaths, union)
                        !print*, 'IntersectionSize = ', intersectionSize
                    else if(intersectionSize == 2) then
                        ! Two paths meet at two end points => Ring Possibly.
                        call PV_append_unique(this%potentialRings, union)
                        !print*, 'IntersectionSize = ', intersectionSize
                    endif
                    call PV_append_unique(this%removePaths, path2)
                endif
            enddo
        enddo

    end subroutine arf_prepPathMergeForOneAtom

    !!----------------------------------------------------------------------------
    ! Clear removePaths, newPaths, potentialRings
    subroutine xlogp_arf_clearTempPaths(this)
        type(AllRingsFinder_t) :: this
        call PV_clear(this%newPaths)
        call PV_clear(this%removePaths)
        call PV_clear(this%potentialRings)
    end subroutine xlogp_arf_clearTempPaths

    subroutine xlogp_arf_showDetectProcess(this)
        type(AllRingsFinder_t) :: this
        if(.not. this%DEBUG) return
        print*, "(*) totInvolvedPathCntPerAtom = ", this%totInvolvedPathCntPerAtom
        print*, "(*) totIntersectPath = ", this%totIntersectPath
        print*, "(#) Atom to remove : ",  &
            ", id = ", this%removeAtom%Id, ", name = ", this%removeAtom%name, &
            ", type = ", this%removeAtom%type, ", num_neib = ", this%removeAtom%num_neib
        print*, "(#) atomList size = ", AIdVector_GetVectorSize(this%atomList)
        print*, "(#) atomList "
        call AIdVector_PathShow(this%atomList)

        print*, "(@) pathGraph size = ", PV_getVectorSize(this%pathGraph)
        call PV_showAll(this%pathGraph)
        print*, "(@) removePaths size = ", PV_getVectorSize(this%removePaths)
        call PV_showAll(this%removePaths)

        print*, "(@) newPaths size = ", PV_getVectorSize(this%newPaths)
        call PV_showAll(this%newPaths)

        print*, "(#) potentialRings size = ", PV_getVectorSize(this%potentialRings)
        call PV_showAll(this%potentialRings)

        print*, "(#) allRings size = ", PV_getVectorSize(this%allRings)
        call PV_showAll(this%allRings)

    end subroutine xlogp_arf_showDetectProcess

    subroutine xlogp_arf_collapsePathGraph(this, pathList, maxPathLen)
        type(AllRingsFinder_t) :: this
        type(PathVector_t), intent(inout) :: pathList
        integer, intent(in) :: maxPathLen ! Max path length = max ring size detected = max recursion depth
!        type(AtomIdVector_t) :: path1
!        type(AtomIdVector_t), pointer :: path_p
        integer :: i
        integer :: nPathSize
!        integer :: numUsed

        !> Remove previous paths
        do i=1, PV_getVectorSize(this%removePaths)
            !call PV_GetPathPtr(this%removePaths, i, path_p)
!            path_p => this%removePaths%data(i)

            call arf_decreaseAtomConnectivity(this, this%removePaths%data(i))
            call PV_remove(pathList, this%removePaths%data(i))
            call PV_removeRevert(pathList, this%removePaths%data(i))
            !print*, '(*) pathList%remove, remaining num = ', pathList%getVectorSize()
            !call path_p%ShowPath()
        enddo

        !> Add new paths already merged above loop
        do i=1, PV_getVectorSize(this%newPaths)
            !call PV_GetPathPtr(this%newPaths, i, path_p)
!            path_p => this%newPaths%data(i)

            nPathSize = AIdVector_GetVectorSize(this%newPaths%data(i))
            if ( (maxPathLen .LT. 0) .OR. (nPathSize .LE. (maxPathLen+1)) ) then
                call arf_increaseAtomConnectivity(this, this%newPaths%data(i))
                this%totAddedPathCnt = this%totAddedPathCnt + 1
                call AIdVector_SetPathId(this%newPaths%data(i), xlogp_arf_getAllBondsNum(this) + this%totAddedPathCnt)
                !numUsed = pathList%getVectorSize()
                call PV_put(pathList, PV_getVectorSize(pathList)+1, this%newPaths%data(i))
                !call pathList%appendUnique(path_p)
                !print*, '(*) pathList%append '
                !call path_p%ShowPath()
            endif
        enddo
    end subroutine xlogp_arf_collapsePathGraph

    subroutine arf_collapsePathGraphSafe(this, pathList, maxPathLen)
        type(AllRingsFinder_t) :: this
        type(PathVector_t), intent(inout) :: pathList
        integer, intent(in) :: maxPathLen ! Max path length = max ring size detected = max recursion depth
!        type(AtomIdVector_t) :: path1
!        type(AtomIdVector_t), pointer :: path_p
        integer :: i
        integer :: nPathSize
!        integer :: numUsed

        !> Remove previous paths
        do i=1, PV_getVectorSize(this%removePaths)
            !call PV_GetPathPtr(this%removePaths, i, path_p)
!            path_p => this%removePaths%data(i)

            call arf_decreaseAtomConnectivity(this, this%removePaths%data(i))
            call PV_remove(pathList, this%removePaths%data(i))
            call PV_removeRevert(pathList, this%removePaths%data(i))
            !print*, '(*) pathList%remove, remaining num = ', pathList%getVectorSize()
            !call path_p%ShowPath()
        enddo

        !> Add new paths already merged above loop
        do i=1, PV_getVectorSize(this%newPaths)
            !call PV_GetPathPtr(this%newPaths, i, path_p)
!            path_p => this%newPaths%data(i)

            nPathSize = AIdVector_GetVectorSize(this%newPaths%data(i))
            if ( (maxPathLen .LT. 0) .OR. (nPathSize .LE. (maxPathLen+1)) ) then
                call arf_increaseAtomConnectivity(this, this%newPaths%data(i))
                this%totAddedPathCnt = this%totAddedPathCnt + 1
                call AIdVector_SetPathId(this%newPaths%data(i), &
                    xlogp_arf_getAllBondsNum(this) + this%totAddedPathCnt)
                !numUsed = pathList%getVectorSize()
                call PV_put(pathList, PV_getVectorSize(pathList)+1, this%newPaths%data(i))
                !call pathList%appendUnique(path_p)
                !print*, '(*) pathList%append '
                !call path_p%ShowPath()
            endif
        enddo
    end subroutine arf_collapsePathGraphSafe

    subroutine arf_decreaseAtomConnectivity(this, removePath)
        type(AllRingsFinder_t) :: this
        type(AtomIdVector_t), intent(inout) :: removePath
!        type(AtomId_t), pointer :: currAtomId_ptr => null()
        type(AtomId_t) :: currAtomId
        integer :: firstId
        integer :: lastId
        integer :: nVecSize
        integer :: i
        integer :: originalNumNeib
        integer :: newNumNeib

        firstId = AIdVector_AtomId_Value(removePath, 1)
        nVecSize = AIdVector_GetVectorSize(removePath)
        lastId = AIdVector_AtomId_Value(removePath, nVecSize)
        do i=1, AIdVector_GetVectorSize(this%atomList)
            currAtomId = AIdVector_Get(this%atomList, i)
            originalNumNeib = currAtomId%num_neib
            if(firstId == currAtomId%id) then
                currAtomId%num_neib = currAtomId%num_neib - 1
            else if(lastId == currAtomId%id) then
                currAtomId%num_neib = currAtomId%num_neib - 1
            endif
            newNumNeib = currAtomId%num_neib

!            if(this%DEBUG .AND. originalNumNeib /= newNumNeib) then
!                print*, "Decrease Connectivity :: AtomId = ", currAtomId_ptr%id, &
!                    " NumNeib Changed. From ", originalNumNeib, &
!                    " To ", newNumNeib
!            endif
        enddo

    end subroutine arf_decreaseAtomConnectivity

    subroutine arf_increaseAtomConnectivity(this, addPath)
        type(AllRingsFinder_t) :: this
        type(AtomIdVector_t), intent(inout) :: addPath
!        type(AtomId_t), pointer :: currAtomId_ptr => null()
        type(AtomId_t) :: currAtomId
        integer :: firstId
        integer :: lastId
        integer :: nVecSize
        integer :: i
        integer :: originalNumNeib
        integer :: newNumNeib

        firstId = AIdVector_AtomId_Value(addPath, 1)
        nVecSize = AIdVector_GetVectorSize(addPath)
        lastId = AIdVector_AtomId_Value(addPath, nVecSize)
        do i=1, AIdVector_GetVectorSize(this%atomList)
            currAtomId = AIdVector_Get(this%atomList, i)
            originalNumNeib = currAtomId%num_neib
            if(firstId == currAtomId%id) then
                currAtomId%num_neib = currAtomId%num_neib + 1
            else if(lastId == currAtomId%id) then
                currAtomId%num_neib = currAtomId%num_neib + 1
            endif
            newNumNeib = currAtomId%num_neib

!            if(this%DEBUG .AND. originalNumNeib /= newNumNeib) then
!                print*, "Add Connectivity :: AtomId = ", currAtomId_ptr%id, &
!                    " NumNeib Changed. From ", originalNumNeib, &
!                    " To ", newNumNeib
!            endif
        enddo

    end subroutine arf_increaseAtomConnectivity

    !!----------------------------------------------------------------------------
    !> Checks the paths if a ring has been found
    !! @param paths The paths to check for rings
    !! @param ringSet The ringset to add the detected rings to
    !! @param ac The AtomContainer with the original structure
    !!----------------------------------------------------------------------------
    subroutine xlogp_arf_detectRings(potentialRings, rings)
        type(PathVector_t), intent(inout) :: potentialRings
        type(PathVector_t), intent(inout) :: rings
        !> arguments
        integer :: i, nPaths  !, nMemberSize
!        type(AtomIdVector_t), pointer :: path_ptr
        !> local variables

        nPaths = PV_getVectorSize(potentialRings)
        do i=1, nPaths
            !call PV_GetPathPtr(potentialRings, i, path_ptr) ! get each path from pathList
!            path_ptr => potentialRings%data(i)

            if(AIdVector_DetermineRing(potentialRings%data(i))) then !! Check if this path is a correct ring.
                call AIdVector_Delete_elements(potentialRings%data(i), 1, 1)
                if(.not. PV_hasElement(rings, potentialRings%data(i))) then
                    !! Add current path into the ring container
                    call PV_append(rings, potentialRings%data(i))
                endif
            endif
        enddo
    end subroutine xlogp_arf_detectRings

    subroutine removeOneAtomAndPaths(this, atomToRemove)
        type(AllRingsFinder_t) :: this
        type(AtomId_t) :: atomToRemove
!        type(AtomIdVector_t), pointer :: path_p
        integer :: i, nLen
        integer :: nPaths
        integer :: nRelavantPathRemove
!        integer :: alloc_stat

        nRelavantPathRemove = 0
        nPaths = PV_getVectorSize(this%pathGraph)

        call CheckMarkSzForPathsWAtomX(this, nPaths)

        call AIdVector_Remove(this%atomList, atomToRemove)
        this%markForPathsWithAtomX = 0

        do i=1, nPaths
            !call PV_GetPathPtr(this%pathGraph, i, path_p)
!            path_p => this%pathGraph%data(i)

            nLen = AIdVector_GetVectorSize(this%pathGraph%data(i))
            if(atomToRemove%id == AIdVector_AtomId_Value(this%pathGraph%data(i), 1) .OR. &
                atomToRemove%id == AIdVector_AtomId_Value(this%pathGraph%data(i), nLen) ) then
                !print*, '(*) removeOneAtomAndRelavantPaths'
                !call path_p%ShowPath()
                if(AIdVector_GetVectorSize(this%pathGraph%data(i)) == 2) then
                    !call path_p%ShowPath()
                    !call this%pathGraph%remove(path_p)
                    this%markForPathsWithAtomX(i) = 1
                    nRelavantPathRemove = nRelavantPathRemove + 1
                endif
            endif
        enddo

        do i=1, nPaths
            if(this%markForPathsWithAtomX(i) == 1) then
                !call PV_GetPathPtr(this%pathGraph, i, path_p)
!                path_p => this%pathGraph%data(i)
                call PV_remove(this%pathGraph, this%pathGraph%data(i))
            endif
        enddo

!        print*, 'nRelavantPathRemove = ', nRelavantPathRemove
    end subroutine removeOneAtomAndPaths

    subroutine deletePathWithTermAtom(pathGraph, atomToRemove)
!        type(AllRingsFinder_t), intent(inout) :: this
        type(PathVector_t), intent(inout) :: pathGraph
        type(AtomId_t), intent(in) :: atomToRemove
        !! END of arguments
        integer :: i, lenPaths, lenOnePath
        integer :: atIdFirst, atIdLast
!        type(AtomIdVector_t), pointer :: path_p
        type(AtomId_t) :: aid1, aid2
        !! END of local variables

        if(atomToRemove%num_neib > 1) return

        lenPaths = PV_getVectorSize(pathGraph)
        do i=1, lenPaths
            !call PV_GetPathPtr(pathGraph, i, path_p)
!            path_p => pathGraph%data(i)

            lenOnePath = AIdVector_GetVectorSize(pathGraph%data(i))
            aid1 = AIdVector_Get(pathGraph%data(i), 1)
            aid2 = AIdVector_Get(pathGraph%data(i), lenOnePath)
            atIdFirst = aid1%id
            atIdLast = aid2%id

            if(atomToRemove%id == atIdFirst .OR. &
                atomToRemove%id == atIdLast) then

!                if(this%DEBUG) then
!                    print*, "(*) Removing one path as the bonded atom is gone."
!                    print*, "(*) Atom ID to remove = ", atomToRemove%id
!                    call path_p%showPath()
!                    print*, "(*) Remaining pathGraph size = ", pathGraph%getVectorSize()
!                endif

                call PV_delete_elements(pathGraph, i, 1)
            endif
        enddo

    end subroutine deletePathWithTermAtom

    !!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    ! Getter and Setter Methods

    integer function xlogp_arf_getRingSize(this)
        type(AllRingsFinder_t), intent(in) :: this
        xlogp_arf_getRingSize = PV_getVectorSize(this%allRings)
    end function xlogp_arf_getRingSize

    subroutine xlogp_arf_getAllRing(this, nRingSize, return_value)
        type(AllRingsFinder_t), intent(in) :: this
        integer :: nRingSize
        type(Ring_t) :: return_value(nRingSize)
        integer :: alloc_stat
        integer :: i, j, atId1, atId2 !, bnId
        type(AtomIdVector_t) :: path

!        nRingSize = xlogp_arf_getRingSize(this)
!        allocate(return_value(nRingSize), stat=alloc_stat)
!        if(alloc_stat /= 0) then
!            print*, '(*) ERROR : ring array allocation in getAllRing(). ', alloc_stat
!            ! @TODO : Do something to treat this exception.
!        endif

        do i=1, nRingSize
            call PV_CopyPath(this%allRings, i, path)
            !call Assignment_AtomIdVector(path,this%allRings%GetPath(i))
            return_value(i)%num_member = AIdVector_GetVectorSize(path) ! ring size
            do j=1, AIdVector_GetVectorSize(path)
                return_value(i)%atom_id(j) = AIdVector_AtomId_Value(path, j)
            enddo

            ! The first bond id
            atId1 = AIdVector_AtomId_Value(path, &
                AIdVector_GetVectorSize(path))
            atId2 = AIdVector_AtomId_Value(path, 1)
            return_value(i)%bond_id(1) = xlogp_arf_getBondId(this, atId1, atId2)

            ! get from 2nd to the last bond id
            do j=2, AIdVector_GetVectorSize(path)
                atId1 = AIdVector_AtomId_Value(path, j-1)
                atId2 = AIdVector_AtomId_Value(path, j)
                return_value(i)%bond_id(j) = xlogp_arf_getBondId(this, atId1, atId2)
            enddo
        enddo

    end subroutine xlogp_arf_getAllRing

    function xlogp_arf_getBondId(this, atomId1, atomId2) result(return_value)
        type(AllRingsFinder_t) :: this
        integer, intent(in) :: atomId1, atomId2
        integer :: return_value
        integer i, atom_1, atom_2
        return_value = 0

        do i=1, xlogp_arf_getAllBondsNum(this)
            atom_1 = this%bondList(i)%atom_1
            atom_2 = this%bondList(i)%atom_2
            if( (atomId1 == atom_1 .AND. atomId2 == atom_2) .OR. &
                (atomId1 == atom_2 .AND. atomId2 == atom_1) ) then
                return_value = this%bondList(i)%id
                return
            endif
        enddo

    end function xlogp_arf_getBondId

    subroutine xlogp_arf_showAllRings(this)
        type(AllRingsFinder_t) :: this
!        type(AtomIdVector_t), pointer :: path_p
        integer :: i

        print*, '(=) Found Ring Number = [ ', PV_getVectorSize(this%allRings), " ] "
        do i=1, PV_getVectorSize(this%allRings)
            !call PV_GetPathPtr(this%allRings, i, path_p)
!            path_p => this%allRings%data(i)

            call AIdVector_PathShow(this%allRings%data(i))
        enddo
    end subroutine xlogp_arf_showAllRings

    subroutine xlogp_arf_showAllRingsDetail(this)
        type(AllRingsFinder_t) :: this
!        type(AtomIdVector_t), pointer :: path_p
        integer :: i

        if(this%DEBUG) then
            print*, '(*) Total Ring Number Found = ', PV_getVectorSize(this%allRings)
            do i=1, PV_getVectorSize(this%allRings)
                print*, '(#) Ring(',i,') => '
                !call PV_GetPathPtr(this%allRings, i, path_p)
!                path_p => this%allRings%data(i)

                call AIdVector_PathShow(this%allRings%data(i))
                call AIdVector_PathShowAsType(this%allRings%data(i))
                call AIdVector_PathShowAsName(this%allRings%data(i))
            enddo
        endif

    end subroutine xlogp_arf_showAllRingsDetail


    subroutine xlogp_arf_showPathList(this)
        type(AllRingsFinder_t) :: this
!        integer :: i

        if(this%DEBUG) then
            print*, '### END OF REMOVE ATOM ###'
            print*, 'After remove path and append path '
            call PV_showAll(this%pathGraph)
            print*, ' '
        endif

    end subroutine xlogp_arf_showPathList

    function xlogp_arf_getAllAtomsNum(this) result(allAtomNum)
        type(AllRingsFinder_t) :: this
        integer :: allAtomNum

        allAtomNum = AIdVector_GetVectorSize(this%atomList)
    end function xlogp_arf_getAllAtomsNum

    function xlogp_arf_getAllBondsNum(this) result(allBondNum)
        type(AllRingsFinder_t) :: this
        integer :: allBondNum

!        allBondNum = size(this%bondList)
        allBondNum = this%bondListSize
    end function xlogp_arf_getAllBondsNum

    subroutine xlogp_arf_setMaxPathLen(this, maxPathLen_)
        type(AllRingsFinder_t) :: this
        integer, intent(in) :: maxPathLen_

        this%maxPathLen = maxPathLen_
    end subroutine xlogp_arf_setMaxPathLen

end module xlogp_AllRingsFinder_m
