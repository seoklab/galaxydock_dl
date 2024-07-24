module xlogp_AromaticRing_m
    use xlogp_constants_m  !, only : XSCORE_MAX_ATOM_NEIB
    use xlogp_Atom_m
    use xlogp_Bond_m
    use xlogp_Group_m
    use xlogp_Ring_m
    use xlogp_Molecule_m

    use xlogp_AllRingsFinder_m, only : AllRingsFinder_t, xlogp_arf_getAllRing, &
        xlogp_arf_doSearch, xlogp_arf_initialize, newAllRingsFinder, &
        xlogp_arf_getRingSize
    implicit none
    private

    public :: DetectRings, DetectAromaticRings
    public :: Is6MemberAromaticRing
    public :: Is5MemberAromaticRing

    logical, parameter :: mydebug = .false.

contains

    !!----------------------------------------------------------------------------
    !> Find the ring systems in the given molecule.
    !! Define all the atoms in a ring as atom.ring=1; otherwise atom.ring = 0.
    !!----------------------------------------------------------------------------
    subroutine DetectRings(atoms, num_atom, bonds, num_bond, rings, ringSize)
        implicit none
        type(Atom_t), intent(inout) :: atoms(:)
        integer, intent(in) :: num_atom
        type(Bond_t), intent(inout) :: bonds(:)
        integer, intent(in) :: num_bond
        type(Ring_t), intent(inout) :: rings(:)
        integer, intent(inout) :: ringSize

        type(AllRingsFinder_t) :: ringFinderMy
        integer :: i, k
        integer :: atomIdOfRing, bondIdOfRing
        integer :: maxPathLen
        integer :: alloc_stat

!        print*, 'num_atom = ', num_atom
!        print*, 'XSCORE_MAX_LIGAND_ATOM_NUMBER = ', XSCORE_MAX_LIGAND_ATOM_NUMBER
        if(num_atom > XSCORE_MAX_LIGAND_ATOM_NUMBER) then
            print*, 'Out of boundary error in num_atom'
        endif

!        print*, 'num_bond = ', num_bond
!        print*, 'XSCORE_MAX_LIGAND_BOND_NUMBER = ', XSCORE_MAX_LIGAND_BOND_NUMBER
        if(num_bond > XSCORE_MAX_LIGAND_BOND_NUMBER) then
            print*, 'Out of boundary error in num_bond'
        endif

        !> Initialize ring information
        do i=1, num_atom
            atoms(i)%ring = 0
        enddo
        do i=1, num_bond
            bonds(i)%ring = 0
        enddo

        call xlogp_arf_initialize(ringFinderMy, atoms, num_atom, bonds, num_bond)
        maxPathLen = num_atom
        call xlogp_arf_doSearch(ringFinderMy, maxPathLen)
        ringSize = xlogp_arf_getRingSize(ringFinderMy)

!        print*, 'ringSize = ', ringSize
!        print*, 'XSCORE_MAX_RING_NUMBER = ', XSCORE_MAX_RING_NUMBER
        if(ringSize > XSCORE_MAX_RING_NUMBER) then
            print*, 'Out of boundary error in ringSize'
        endif

        call xlogp_arf_getAllRing(ringFinderMy, ringSize, rings)

        do k=1, ringSize
            do i=1, rings(k)%num_member
                atomIdOfRing = rings(k)%atom_id(i)
                atoms(atomIdOfRing)%ring = 1
                bondIdOfRing = rings(k)%bond_id(i)
                bonds(bondIdOfRing)%ring = 1
            enddo
        enddo

    end subroutine DetectRings

    !!----------------------------------------------------------------------------
    !> Detect all the atoms in aromatic rings (6-member and 5-member), and
    !! define them as atom.ring = 2
    !!----------------------------------------------------------------------------
    subroutine DetectAromaticRings(atoms, bonds, rings, ringSize)
        implicit none
        type(Atom_t), intent(inout) :: atoms(:)
        type(Bond_t), intent(inout) :: bonds(:)
        type(Ring_t), intent(inout) :: rings(:)
        integer, intent(inout) :: ringSize

        type(Ring_t) :: ring
        integer :: i

        !> Check all the 6-membered rings first
        do i=1, ringSize
            ring = rings(i)
            if(xlogp_Ring_GetNumMember(ring) /= 6) cycle
            if(Is6MemberAromaticRing(atoms, bonds, ring)) then
                call CalcRingCentroid(atoms, ring)
                call SetRingAsAromatic(atoms, bonds, ring)
            endif
        enddo

        !> Check all the 5-membered rings first
        do i=1, ringSize
            ring = rings(i)
            if(xlogp_Ring_GetNumMember(ring) /= 5) cycle
            if(Is5MemberAromaticRing(atoms, bonds, ring)) then
                call CalcRingCentroid(atoms, ring)
                call SetRingAsAromatic(atoms, bonds, ring)
            endif
        enddo

    end subroutine DetectAromaticRings

    !!----------------------------------------------------------------------------
    !> Set ring and related atom's property as aromatic.
    !!----------------------------------------------------------------------------
    subroutine SetRingAsAromatic(atoms, bonds, ring)
        type(Atom_t), intent(inout) :: atoms(:)
        type(Bond_t), intent(inout) :: bonds(:)
        type(Ring_t), intent(inout) :: ring
        integer :: i, atomid, bondId

        do i=1, xlogp_Ring_GetNumMember(ring)
            atomId = ring%atom_id(i)
            bondId = ring%bond_id(i)
            atoms(atomId)%ring = 2 ! Aromatic ring
            bonds(bondId)%ring = 2 ! Aromatic ring
        enddo
        ring%valid = .true.
        ring%type = 2 ! Aromatic Ring
    end subroutine SetRingAsAromatic

    !!----------------------------------------------------------------------------
    !> calculate the centroid of the given ring
    !!----------------------------------------------------------------------------
    subroutine CalcRingCentroid(atoms, ring)
        type(Atom_t), intent(in) :: atoms(:)
        type(Ring_t), intent(inout) :: ring
        integer :: i, atomId

        ring%centroid = 0.0
        do i=1, xlogp_Ring_GetNumMember(ring)
            atomId = ring%atom_id(i)
            ring%centroid(1) = ring%centroid(1) + atoms(atomId)%coor(1)
            ring%centroid(2) = ring%centroid(2) + atoms(atomId)%coor(2)
            ring%centroid(3) = ring%centroid(3) + atoms(atomId)%coor(3)
        enddo
        ring%centroid(1) = ring%centroid(1) / xlogp_Ring_GetNumMember(ring)
        ring%centroid(2) = ring%centroid(2) / xlogp_Ring_GetNumMember(ring)
        ring%centroid(3) = ring%centroid(3) / xlogp_Ring_GetNumMember(ring)

    end subroutine CalcRingCentroid

    function Is6MemberAromaticRing(atoms, bonds, ring) result(return_value)
        type(Atom_t), intent(in) :: atoms(:)
        type(Bond_t), intent(in) :: bonds(:)
        type(Ring_t), intent(in) :: ring
        logical :: return_value
        integer :: i, iPiElec
        integer :: atomId, bondId, bondIdNext
        logical :: mark, mark1, mark2, markContinuousTwoSingleBonds

        !> Check the number of pi electrons
        mark = .true.
        iPiElec = 0
        do i=1, 6
            atomId = ring%atom_id(i)
            if (atoms(atomId)%type == 'C.ar') then
                iPiElec = iPiElec + 1
            else if(atoms(atomId)%type == 'N.ar') then
                iPiElec = iPiElec + 1
            else if(atoms(atomId)%type == 'C.2') then
                iPiElec = iPiElec + 1
            else if(atoms(atomId)%type == 'N.2') then
                iPiElec = iPiElec + 1
            else
                mark = .false.
                if(mydebug) print*, "(*) No Atoms with pi electrons : ", atoms(atomId)%type
                EXIT
            endif
        enddo

        if( (mark .EQV. .false.) .OR. (iPiElec /= 6) ) then
            if(mydebug) print*, 'mark = ', mark
            if(mydebug) print*, 'iPiElec = ', iPiElec

            return_value = .false.
            return
        endif

        !> Check if there are [two continuous single bonds] on the ring's bond path.
        !! This algorithm is superior to check -=-=-=
        markContinuousTwoSingleBonds = .false.
        CONT_SINGLE_BONDS: do i=1, 5
            bondId = ring%bond_id(i)
            if( (bonds(bondId)%type == '1') .OR. &
                (bonds(bondId)%type == 'am') ) then ! 'am' => amide bond (single bond)
                mark1 = .true.
            else
                mark1 = .false.
            endif

            bondIdNext = ring%bond_id(i+1)
            if( (bonds(bondIdNext)%type == '1') .OR. &
                (bonds(bondIdNext)%type == 'am') ) then ! 'am' => amide bond (single bond)
                mark2 = .true.
            else
                mark2 = .false.
            endif

            if (mark1 .AND. mark2) then
                if(Atom_IsAromaticType(atoms(bonds(bondId)%atom_1))) then
                    cycle CONT_SINGLE_BONDS
                else if(Atom_IsAromaticType(atoms(bonds(bondId)%atom_2))) then
                    cycle CONT_SINGLE_BONDS
                else if(Atom_IsAromaticType(atoms(bonds(bondIdNext)%atom_1))) then
                    cycle CONT_SINGLE_BONDS
                else if(Atom_IsAromaticType(atoms(bonds(bondIdNext)%atom_2))) then
                    cycle CONT_SINGLE_BONDS
                endif

                markContinuousTwoSingleBonds = .true.

                if(mydebug) then
                    print*, '(*) Two continuous single bonds : ', &
                        ' bondid1 = ', bondId, ' bondid2 = ', bondIdNext
                    print*, '(*) i = ', i
                    call Atom_Show(atoms(bonds(bondId)%atom_1))
                    call Atom_Show(atoms(bonds(bondId)%atom_2))
                    call Atom_Show(atoms(bonds(bondIdNext)%atom_1))
                    call Atom_Show(atoms(bonds(bondIdNext)%atom_2))
                endif

                EXIT ! No aromatic ring
            endif
        enddo CONT_SINGLE_BONDS

        if(markContinuousTwoSingleBonds) then
            if(mydebug) print*, 'Two continuous single bonds'
            return_value = .false. ! Not aromaric ring
        else
            return_value = .true.  ! Aromatic Ring
        endif

    !        print*, 'This is the end of Is6MemberAromaticRing'
    end function Is6MemberAromaticRing

    function Is5MemberAromaticRing(atoms, bonds, ring) result(return_value)
        type(Atom_t), intent(in) :: atoms(:)
        type(Bond_t), intent(in) :: bonds(:)
        type(Ring_t), intent(in) :: ring
        logical :: return_value
        integer :: i, iPiElec
        integer :: atomId, bondId, bondIdNext
        logical :: mark, mark1, mark2, markContinuousTwoSingleBonds
        integer :: pi_path(5)
        integer, parameter :: ONE = 1
        integer, parameter :: TWO = 2

        if(mydebug) then
            print*, '(*) Is5MemberAromaticRing was called '
            call xlogp_Ring_Show(ring)
        endif

        !> Check the number of pi electrons
        mark = .true.
        iPiElec = 0
        do i=1, 5
            atomId = ring%atom_id(i)
            if (atoms(atomId)%type == 'C.ar') then
                iPiElec = iPiElec + ONE
                pi_path(i) = ONE
            else if(atoms(atomId)%type == 'N.ar') then
                iPiElec = iPiElec + ONE
                pi_path(i) = ONE
            else if(atoms(atomId)%type == 'C.2') then
                iPiElec = iPiElec + ONE
                pi_path(i) = ONE
            else if(atoms(atomId)%type == 'N.2') then
                iPiElec = iPiElec + ONE
                pi_path(i) = ONE
            else if(atoms(atomId)%type == 'N.pl3') then
                iPiElec = iPiElec + TWO ! Nitrogen trigonal planar (N.pl3) has two pi electrons.
                pi_path(i) = TWO
            else if(atoms(atomId)%type == 'N.pl3.h') then
                iPiElec = iPiElec + TWO ! Nitrogen trigonal planar (N.pl3) has two pi electrons.
                pi_path(i) = TWO
            else if(atoms(atomId)%type == 'N.am') then
                iPiElec = iPiElec + TWO ! Nitrogen amide (N.am) has two pi electrons.
                pi_path(i) = TWO
            else if(atoms(atomId)%type == 'O.3') then
                iPiElec = iPiElec + TWO
                pi_path(i) = TWO
            else if(atoms(atomId)%type == 'S.3') then
                iPiElec = iPiElec + TWO
                pi_path(i) = TWO
            else
                mark = .false.
                if(mydebug) print*, '(*) No aromatic atomtype = ', atoms(atomId)%type
                EXIT
            endif
        enddo

        if(mydebug) print*, 'iPiElec = ', iPiElec

        if(mark .EQV. .false. .OR. iPiElec /= 6) then
            return_value = .false.
            return
        endif

        !> Check if there are two continuous single bonds on the path,
        !! but, it is okay if these two single bonds explain the special atom.
        markContinuousTwoSingleBonds = .false.

        !! Here we iterate from 1 to 4
        CONT_SINGLE_BONDS : do i=1, 4
            bondId = ring%bond_id(i)
            if( (bonds(bondId)%type == '1') .OR. &
                (bonds(bondId)%type == 'am') ) then ! 'am' => amide bond (single bond)
                mark1 = .true.
            else
                mark1 = .false.
            endif

            bondIdNext = ring%bond_id(i+1)
            if( (bonds(bondIdNext)%type == '1') .OR. &
                (bonds(bondIdNext)%type == 'am') ) then ! 'am' => amide bond (single bond)
                mark2 = .true.
            else
                mark2 = .false.
            endif

            if ((.not. mark1) .OR. (.not. mark2)) then
                cycle CONT_SINGLE_BONDS
            else if(pi_path(i) == TWO) then
                cycle CONT_SINGLE_BONDS
            else if(pi_path(i+1) == TWO) then
                cycle CONT_SINGLE_BONDS
            else
                if(Atom_IsAromaticType(atoms(bonds(bondId)%atom_1))) then
                    cycle CONT_SINGLE_BONDS
                else if(Atom_IsAromaticType(atoms(bonds(bondId)%atom_2))) then
                    cycle CONT_SINGLE_BONDS
                else if(Atom_IsAromaticType(atoms(bonds(bondIdNext)%atom_1))) then
                    cycle CONT_SINGLE_BONDS
                else if(Atom_IsAromaticType(atoms(bonds(bondIdNext)%atom_2))) then
                    cycle CONT_SINGLE_BONDS
                endif

                if(mydebug)  then
                    print*, '(*) Two continuous single bonds : ', &
                        ' bondid1 = ', bondId, ' bondid2 = ', bondIdNext
                    print*, '(*) i = ', i
                    call Atom_Show(atoms(bonds(bondId)%atom_1))
                    call Atom_Show(atoms(bonds(bondId)%atom_2))
                    call Atom_Show(atoms(bonds(bondIdNext)%atom_1))
                    call Atom_Show(atoms(bonds(bondIdNext)%atom_2))
                endif

                markContinuousTwoSingleBonds = .true.
                EXIT
            endif
        enddo CONT_SINGLE_BONDS

        if(markContinuousTwoSingleBonds) then
            return_value = .false.
        else
            return_value = .true.
        endif

    end function Is5MemberAromaticRing

end module xlogp_AromaticRing_m
