module xlogp_AtomTypeChecker_m
    use xlogp_constants_m, only : XSCORE_MAX_ATOM_NEIB
    use xlogp_Atom_m
    use xlogp_Bond_m
    use xlogp_Group_m
    implicit none

    private ! Information Hiding

    public :: CheckAllAtomTypes

!SYBYL Atom Type
!C.3 carbon sp3
!C.2 carbon sp2
!C.1 carbon sp1
!
!C.ar carbon aromatic
!C.cat carbocation (C+) used only in a guadinium group
!
!N.3 nitrogen sp3
!N.2 nitrogen sp2
!N.1 nitrogen sp1
!N.ar nitrogen aromatic
!N.am nitrogen amide
!N.pl3  nitrogen trigonal planar
!N.4 nitrogen sp3 positively charged
!
!O.3 oxygen sp3
!O.2 oxygen sp2
!O.co2 oxygen in carboxylate and phosphate groups
!O.spc oxygen in Single Point Charge (SPC) water model
!O.t3p oxygen in Transferable Intermolecular Potential (TIP3P) water model
!S.3 sulfur sp3
!S.2 sulfur sp2
!S.O sulfoxide sulfur (S=O)
!S.O2 sulfone sulfur  (O=S=O)
!
!P.3 phosphorous sp3
!F fluorine

contains

    subroutine CheckCarbon(atomList, numAtom, bondList, mark)
        type(Atom_t), DIMENSION(:), target, intent(inout) :: atomList
        integer, intent(in) :: numAtom
        type(Bond_t), DIMENSION(:), target, intent(inout) :: bondList
        logical, intent(inout) :: mark
        !! end of input arguments
        integer i
        type(Group_t) group

        mark = .true.

        do i=1, numAtom
            if(atomList(i)%type == 'C.2') then
                call xlogp_Group_FindAGroup(group, atomList, bondList, atomList(i)%id)
                !! =C=
                if((group%num_db == 2) .AND. (group%num_nonh == 2) .AND. (group%num_neib == 2)) then
                    atomList(i)%type = 'C.1'
                else if(group%num_db < 1) then
                    mark = .false.
                    !print*, "Waring: carbon atom have a wrong type"
                    !print*, i, ", id = ", atomList(i)%id, ", name = ", atomList(i)%name
                    EXIT
                end if
            else if(atomList(i)%type == 'C.1') then
                call xlogp_Group_FindAGroup(group, atomList, bondList, atomList(i)%id)
                if(group%num_db < 1 .AND. group%num_tb < 1) then
                    mark = .false.
                    !print*, "Waring: carbon atom have a wrong type"
                    !print*, i, ", id = ", atomList(i)%id, ", name = ", atomList(i)%name
                    EXIT
                end if
            end if
        end do

        !if(mark .eqv. .false.) then
        !    print*, "Warning: some carbon atoms have wrong types."
        !end if
    end subroutine CheckCarbon

    subroutine CheckHydrogen(atomList, numAtom, bondList, mark)
        type(Atom_t), DIMENSION(:), target, intent(inout) :: atomList
        integer, intent(in) :: numAtom
        type(Bond_t), DIMENSION(:), target, intent(inout) :: bondList
        logical, intent(inout) :: mark
        !! end of input arguments
        integer i
        integer countH, requireH
        character(len=10) atType

        countH = 0
        requireH = 0
        mark = .true.

        do i=1, numAtom
            atType = TRIM(atomList(i)%type)
            ! Count hydrogen atoms
            if(atType == 'H') then
                if(atomList(atomList(i)%neib(1))%type(1:1) /= 'C') then
                    cycle
                else
                    countH = countH + 1
                end if
            end if

            ! Calculate required hydrogen atoms
            if(atType == 'C.3') then
                requireH = requireH + (4 - atomList(i)%num_nonh)
            else if(atType == 'C.2') then
                requireH = requireH + (3 - atomList(i)%num_nonh)
            else if(atType == 'C.ar') then
                requireH = requireH + (3 - atomList(i)%num_nonh)
            else if(atType == 'C.cat') then !! carbocation
                requireH = requireH + (3 - atomList(i)%num_nonh)
            else if(atType == 'C.1') then
                requireH = requireH + (2 - atomList(i)%num_nonh)
            else
                cycle
            end if
        end do

        if(countH /= requireH) then
       !     print*, 'Warning: hydrogen atoms may not be correctly added.'
       !     print*, 'countH = ', countH, ', requireH = ', requireH
            mark = .false.
        end if

    end subroutine CheckHydrogen

    subroutine CheckOxygen(group)
        type(Group_t), intent(inout) :: group
        integer j

        if(group%center%type(1:1) == 'O') then
            if(group%num_neib == 1 .AND. group%num_nonh == 1) then !! =O
                if(group%center%type == 'O.co2') then   !! carboxylate or phosphate groups
                    group%bond(1)%type = '2'
                else if(group%center%type == 'O.2') then !! sp2 oxygen
                    group%bond(1)%type = '2'
                else
                    group%center%type = 'O.2'
                    group%bond(1)%type = '2'
                end if
            else if(group%num_neib >= 2) then
                group%center%type = 'O.3'
                do j=1, group%num_nonh
                    group%bond(j)%type = '1'
                end do
            end if
        end if
    end subroutine CheckOxygen

    subroutine CheckPhosphorus(group)
        type(Group_t), intent(inout) :: group
        integer j

        if(group%center%type(1:1) == 'P') then
            group%center%type = 'P.3'
        end if
    end subroutine CheckPhosphorus

    subroutine CheckSulfur(group)
        type(Group_t), intent(inout) :: group
        integer j
        logical :: mydebug
        mydebug = .false.

        if(group%center%type(1:1) == 'S') then
           if(mydebug) then
            print*, 'CheckSulfur called'
            print*, 'db_type = ', group%db_type
            print*, 'num_neib = ', group%num_neib
            print*, 'num_nonh = ', group%num_nonh
            print*, 'group%center%type = ', group%center%type
            endif

            if(group%db_type == 3 .AND. group%num_neib == 3) then !! -SO-
                group%center%type = 'S.o'
            else if(group%db_type == 3 .AND. group%num_neib == 4) then !! -SO2
                group%center%type = 'S.o2'
            else if(group%num_nonh >= 4) then !! e.g. -SF5
                group%center%type = 'S.3'
                do j=1, group%num_nonh
                    group%bond(j)%type = '1'
                end do
            else if(group%num_nonh == 1 .AND. group%num_neib == 1) then !! =S
                group%center%type = 'S.2'
                group%bond(1)%type = '2'
            end if
            if(mydebug) then
               print*, 'group%center%type = ', group%center%type
            endif
        end if

    end subroutine CheckSulfur

    subroutine CheckNitrogen(group)
        type(Group_t), intent(inout) :: group
        integer j

        if(group%center%type == 'N.2') then
            if(group%num_db == 2 .AND. &
                group%num_nonh == 2 .AND. &
                group%num_neib == 2) then !! =N=
                group%center%type = 'N.1'
            end if
        else if (group%center%type(1:1) == 'N') then
            if(group%db_type == 3 .AND. group%num_nonh >= 2) then !! -NO, -NO2
                group%center%type = 'N.2'
            else if(group%center%type == 'N.3' .OR. group%center%type == 'N.4') then
                if(group%amide > 0) then !! -NH-SO-, -NH-PO-, or -NH-CO-
                    group%center%type = 'N.am'
                    group%bond(group%amide)%type = 'am'
                else if(group%num_pi == group%num_nonh) then
                    group%center%type = 'N.pl3'
                end if
            else if(group%center%type == 'N.pl3') then
                if(group%amide > 0) then !! -NH-SO-, -NH-PO-, or -NH-CO-
                    group%center%type = 'N.am'
                    group%bond(group%amide)%type = 'am'
                end if
            end if
        end if
    end subroutine CheckNitrogen

    !!----------------------------------------------------------------------------
    !> This function will check N, O, S, P atoms.
    !!----------------------------------------------------------------------------
    subroutine CheckHeteroAtoms(atomList, numAtom, bondList, mark)
        type(Atom_t), DIMENSION(:), target, intent(inout) :: atomList
        integer, intent(in) :: numAtom
        type(Bond_t), DIMENSION(:), target, intent(inout) :: bondList
        logical, intent(inout) :: mark
        !! end of input arguments
        type(Group_t) :: group
        integer i, j, id

        do i=1, numAtom
            call xlogp_Group_FindAGroup(group, atomList, bondList, atomList(i)%id)
            !call group%Show()
            call CheckOxygen(group)
            call CheckPhosphorus(group)
            call CheckSulfur(group)
            call CheckNitrogen(group)

            !! Now assign correct atom type for this atom
            atomList(i)%type = group%center%type

            !! Now assign correct bond type for this atom
            do j=1, group%num_nonh
                id = group%bond(j)%id
                bondList(id)%type = group%bond(j)%type
            end do
        end do
    end subroutine CheckHeteroAtoms


    subroutine CheckAllAtomTypes(atomList, numAtom, bondList, mark)
        type(Atom_t), DIMENSION(:), target, intent(inout) :: atomList
        integer, intent(in) :: numAtom
        type(Bond_t), DIMENSION(:), target, intent(inout) :: bondList
        logical, intent(inout) :: mark

        call CheckCarbon(atomList, numAtom, bondList, mark)

        if(mark .eqv. .false.) then
        !    print*, 'Warning : some carbon atoms have wrong types.'
            return
        end if

        call CheckHydrogen(atomList, numAtom, bondList, mark)
        if(mark .eqv. .false.) then
        !    print*, 'Warning: hydrogen atoms may not be correctly added.'
            return
        end if

        call CheckHeteroAtoms(atomList, numAtom, bondList, mark)

        return
    end subroutine CheckAllAtomTypes


end module xlogp_AtomTypeChecker_m
