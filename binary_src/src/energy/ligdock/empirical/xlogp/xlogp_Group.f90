module xlogp_Group_m
    !!----------------------------------------------------------------------------
    !> Group for X-SCORE
    !!----------------------------------------------------------------------------
    use xlogp_constants_m, only : XSCORE_MAX_ATOM_NEIB
    use xlogp_Atom_m
    use xlogp_Bond_m
    implicit none


    !------------------------------------------------------------
    ! Group Derived Type
    type Group_t

        !sequence ! Let this group type stored as sequential
        logical :: lNeighborData != .false.

        logical :: valid    != .false.
        integer :: num_neib    != 0     ! number of neighboring atoms
        integer :: num_nonh    != 0     !< number of neighboring non-h atoms
        integer :: num_h       != 0     !< number of neighboring hydrogen atoms
        integer :: num_hetero  != 0     !< number of neighboring heteroatoms
        integer :: num_pi      != 0     !< number of neighboring pi atoms

        integer :: num_car     != 0     !< number of neighboring C.ar (Counted only for Find_X_Group)
        integer :: num_nar     != 0     !< number of neighboring N.ar (Counted only for Find_X_Group)

        integer :: num_db      != 0     !< number of double bonds
        integer :: num_tb      != 0     !< number of triple bonds
        integer :: db_type     != 0     !< double bond type, This is used to discriminate
                                    !! which atom is double-bonded to center atom
        integer :: amide       != 0     !< indicator for amino group
        type(Atom_t) :: center      !< center atom of the group
        type(Atom_t), dimension(XSCORE_MAX_ATOM_NEIB) :: neib !< neighboring atoms
        type(Bond_t), dimension(XSCORE_MAX_ATOM_NEIB) :: bond !< neighboring bonds
    end type Group_t

contains

    !------------------------------------------------------------
    ! Clear the group object
    subroutine xlogp_Group_Clear(g)
        implicit none
        type(Group_t), intent(inout) :: g
        integer :: i

        g%lNeighborData = .false.
        g%valid = .false.
        g%num_neib = 0
        g%num_nonh = 0
        g%num_h = 0
        g%num_hetero = 0
        g%num_pi = 0
        g%num_car = 0
        g%num_nar = 0
        g%num_db = 0 ! double bond
        g%db_type = 0
        g%num_tb = 0
        g%amide = 0

        call xlogp_Atom_Clear(g%center)
        do i=1, XSCORE_MAX_ATOM_NEIB
            call xlogp_Atom_Clear(g%neib(i))  ! AtomClear(g%neib(i))
            call xlogp_Bond_Clear(g%bond(i))
        end do

    end subroutine xlogp_Group_Clear

    subroutine xlogp_Group_ShowContents(g)
        implicit none
        type(Group_t), intent(in) :: g
        integer :: i

        print*, "Group center"
        call Atom_Show(g%center) ! AtomShowContents(g%center)

        ! Format Specification
10      FORMAT("Neighbor (", I5,") ", 3X, "Atom ", I5, 3X, A5, 3X, A5, ", Bond ", 3X, A5)

        do i=1, g%num_neib
            write(*, 10) i, g%neib(i)%id, g%neib(i)%name, &
                g%neib(i)%type, g%bond(i)%type
        end do

        ! Format Specification
        ! Iw : integer
        ! Fw.d : w chars and d decimal places
        ! Ew.d : exponent format
        ! Aw : string in w chars
        ! Lw : Print w-1 blanks, then T or F for logical
        ! wX : Print w blanks
        ! Tc : Move ("tab") to absolute position c
        ! TLw : Move left w characters
        ! TRw : Move right w characters
        ! / : Start new line
        !
20      FORMAT("Group Info : ", A, " = ", I3, 3X)
30      FORMAT("Group Info : ", A, " = ", L3, 3X)
        write(*, 30) "valid       ", g%valid
        write(*, 20) "num_neib    ", g%num_neib
        write(*, 20) "num_nonh    ", g%num_nonh
        write(*, 20) "num_h       ", g%num_h
        write(*, 20) "num_hetero  ", g%num_hetero
        write(*, 20) "num_pi      ", g%num_pi
        write(*, 20) "num_car     ", g%num_car
        write(*, 20) "num_nar     ", g%num_nar
        write(*, 20) "num_db      ", g%num_db
        write(*, 20) "num_tb      ", g%num_tb
        write(*, 20) "db_type     ", g%db_type
        write(*, 20) "amide       ", g%amide

    end subroutine xlogp_Group_ShowContents


    subroutine xlogp_Group_SetNeighborInfo(this, atomList, bondList, atomId)
        type(Group_t), intent(inout) :: this
        type(Atom_t), DIMENSION(:), intent(in) :: atomList
        type(Bond_t), DIMENSION(:), intent(in) :: bondList
        integer, intent(in) :: atomId
        integer :: i
        ! END OF DECLARE

        call xlogp_Group_Clear(this)

        ! Specify center atom at first
        this%center = atomList(atomId)


        this%num_neib = this%center%num_neib

        ! Now copy neighboring atoms and bonds into this group builder
        do i=1, this%num_neib
            this%neib(i) = atomList(this%center%neib(i))
            this%bond(i) = bondList(this%center%bond(i))
        end do

        this%lNeighborData = .true.
    end subroutine xlogp_Group_SetNeighborInfo

    subroutine xlogp_Group_CountNumHNonH(this)
        type(Group_t), intent(inout) :: this
        integer :: i

        this%num_h = 0
        this%num_nonh = 0

        DO i=1, this%num_neib
            if(this%neib(i)%type(1:1) == 'H') then
                this%num_h = this%num_h + 1
            else
                this%num_nonh = this%num_nonh + 1
            end if
        end do
    end subroutine xlogp_Group_CountNumHNonH

    !!----------------------------------------------------------------------------
    !> This will count all the neighboring atoms excluding hydrogens
    !!
    !!----------------------------------------------------------------------------
    subroutine xlogp_Group_CountHeteroNum(this)
        type(Group_t), intent(inout) :: this
        integer :: i
        character(len=2) atomType

        this%num_hetero = 0
        DO i=1, this%num_neib

            atomType = TRIM(this%neib(i)%type)

            if (atomType == 'F') then
                this%num_hetero=this%num_hetero+1
            else if(atomType(1:2) == 'Cl') then
                this%num_hetero=this%num_hetero+1
            else if(atomType(1:2) == 'CL') then
                this%num_hetero=this%num_hetero+1
            else if(atomType(1:2) == 'Br') then
                this%num_hetero=this%num_hetero+1
            else if(atomType(1:2) == 'BR') then
                this%num_hetero=this%num_hetero+1
            else if(atomType == 'I') then
                this%num_hetero=this%num_hetero+1
            else if(atomType(1:2) == 'Si') then
                this%num_hetero=this%num_hetero+1
            else if(atomType(1:2) == 'SI') then
                this%num_hetero=this%num_hetero+1
            else if(atomType(1:1) == 'N') then
                this%num_hetero=this%num_hetero+1
            else if(atomType(1:1) == 'O') then
                this%num_hetero=this%num_hetero+1
            else if(atomType(1:1) == 'P') then
                this%num_hetero=this%num_hetero+1
            else if(atomType(1:1) == 'S') then
                this%num_hetero=this%num_hetero+1
            else
                cycle
            end if
        END DO

    end subroutine xlogp_Group_CountHeteroNum

    subroutine xlogp_Group_CountHetNum_XLogp(this)
        type(Group_t), intent(inout) :: this
        integer :: i
        character(len=2) atomType

        this%num_hetero = 0

        DO i=1, this%num_nonh
            if(this%bond(i)%type == '2') then
                cycle
            else if(this%bond(i)%type == '3') then
                cycle
            else if(this%bond(i)%type == 'ar') then
                cycle
            else if(this%neib(i)%type(1:1) == 'N') then
                this%num_hetero = this%num_hetero + 1
            else if(this%neib(i)%type(1:1) == 'O') then
                this%num_hetero = this%num_hetero + 1
            else
                cycle
            end if
        END DO

    end subroutine xlogp_Group_CountHetNum_XLogp

    subroutine Group_CountDbTriBondNumType(this)
        type(Group_t), intent(inout) :: this
        integer :: i
        character(len=10) :: atType
        character(len=3) :: bnType

        this%db_type = 0
        this%num_db = 0
        this%num_tb = 0

        do i=1, this%num_nonh
            if(this%bond(i)%type == '2') then
                this%num_db = this%num_db + 1
                atType = this%neib(i)%type

                ! Determine double bond type to represent
                ! from which atom the double bond was derived.
                if (atType(1:1) == 'C') then
                    this%db_type = 1        ! For atom type 'C', then double bond type is 1
                else if (atType(1:1) == 'N') then
                    this%db_type = 2    ! For 'N' atom, db_type is 2
                else if (atType(1:1) == 'O') then
                    this%db_type = 3
                else if (atType(1:1) == 'S') then
                    this%db_type = 4
                else if (atType(1:1) == 'P') then
                    this%db_type = 5
                else
                    cycle
                end if

            else if(this%bond(i)%type == '1' .AND. &
                this%neib(i)%type == 'O.co2') then

                this%db_type = 3
                this%num_db = this%num_db + 1

            else if(this%bond(i)%type == 'ar' .AND. &
                this%neib(i)%type == 'O.co2') then

                this%db_type = 3
                this%num_db = this%num_db + 1

            else if(this%bond(i)%type == '3') then

                this%num_tb = this%num_tb+1

            else
                cycle
            end if
        end do

    end subroutine Group_CountDbTriBondNumType

    !!----------------------------------------------------------------------------
    !> This will count the number of atoms with pi electrons.
    !!----------------------------------------------------------------------------
    subroutine xlogp_Group_CountPiAtoms(this)
        type(Group_t), intent(inout) :: this
        integer :: i
        character(len=10) :: atType
        character(len=3) :: bnType

        this%num_pi = 0

        do i=1, this%num_nonh
            if(this%bond(i)%type == '2') then
                cycle
            else if(this%bond(i)%type == '3') then
                cycle
            else if(this%neib(i)%type == 'C.ar') then
                this%num_pi = this%num_pi + 1
            else if(this%neib(i)%type == 'C.2') then
                this%num_pi = this%num_pi + 1
            else if(this%neib(i)%type == 'C.1') then
                this%num_pi = this%num_pi + 1
            else if(this%neib(i)%type == 'C.cat') then
                this%num_pi = this%num_pi + 1
            else if(this%neib(i)%type == 'N.2') then
                this%num_pi = this%num_pi + 1
            else if(this%neib(i)%type == 'N.1') then
                this%num_pi = this%num_pi + 1
            else if(this%neib(i)%type == 'N.ar') then
                this%num_pi = this%num_pi + 1
            else
                cycle
            end if
        end do

    end subroutine xlogp_Group_CountPiAtoms


    subroutine xlogp_Group_CountNumNarCar(this)
        type(Group_t), intent(inout) :: this
        integer :: i
        character(len=10) :: atType
        character(len=3) :: bnType

        this%num_nar = 0
        this%num_car = 0

        do i=1, this%num_nonh
            if(this%bond(i)%type /= 'ar') then
                cycle
            else if(this%neib(i)%type == 'N.ar') then
                this%num_nar = this%num_nar + 1
            else if(this%neib(i)%type == 'C.ar') then
                this%num_car = this%num_car + 1
            else if(this%neib(i)%type(1:1) == 'C') then
                this%num_car = this%num_car + 1
            else
                this%num_nar = this%num_nar + 1
            end if
        end do

    end subroutine xlogp_Group_CountNumNarCar

    subroutine xlogp_Group_SetDblType_XLogp(this)
        type(Group_t), intent(inout) :: this
        integer :: i
        character(len=10) :: atType
        character(len=3) :: bnType

        this%db_type = 0

        do i=1, this%num_nonh
            if(this%center%type == 'O.co2' .OR. &
                this%center%type == 'O.2' .OR. &
                this%center%type == 'S.2' ) then

                if(this%neib(i)%type(1:1) == 'C') then
                    this%db_type = 1
                else if(this%neib(i)%type(1:1) == 'N') then
                    this%db_type = 2
                else if(this%neib(i)%type(1:1) == 'O') then
                    this%db_type = 3
                else if(this%neib(i)%type(1:1) == 'S') then
                    this%db_type = 4
                else if(this%neib(i)%type(1:1) == 'P') then
                    this%db_type = 5
                end if

            else if(this%bond(i)%type == '2' .OR. &
                this%neib(i)%type == 'O.co2' .OR. &
                this%neib(i)%type == 'O.2' .OR. &
                this%neib(i)%type == 'S.2' ) then

                if(this%neib(i)%type(1:1) == 'C') then
                    this%db_type = 1
                else if(this%neib(i)%type(1:1) == 'N') then
                    this%db_type = 2
                else if(this%neib(i)%type(1:1) == 'O') then
                    this%db_type = 3
                else if(this%neib(i)%type(1:1) == 'S') then
                    this%db_type = 4
                else if(this%neib(i)%type(1:1) == 'P') then
                    this%db_type = 5
                end if

            else
                cycle
            end if
        end do

    end subroutine xlogp_Group_SetDblType_XLogp

    subroutine xlogp_Group_CheckAmide(this, atomList)
        type(Group_t), intent(inout) :: this
        type(Atom_t), DIMENSION(:), intent(in) :: atomList
        integer :: i, j, id, num
        logical :: mark
        !> END OF DECLARE

        !! Check if the central atom is adjacent to any -SO-, -PO-, or -CO-
        mark = .false.

        OUTLOOP : do i=1, this%num_nonh
            !if(this%neib(i)%type == 'P.3' .OR. &
            !    this%neib(i)%type == 'S.o' .OR. &
            !    this%neib(i)%type == 'S.o2' .OR. &
            !    this%neib(i)%type == 'C.2' ) then
            if(this%neib(i)%type == 'C.2' ) then

                !! Loop over all the neighboring non-hydrogen atoms' neighboring atoms
                !! Here, we search for the Second Shell of center atom.
                num = this%neib(i)%num_nonh
                INLOOP : do j=1, num
                    ! Get a neighbor atom's neighbor id ...
                    id = atomList(this%neib(i)%id)%neib(j)
                    if (id == this%center%id) then
                        cycle INLOOP
                    else if(atomList(id)%type == 'O.2') then
                        mark = .true.
                        exit INLOOP
                    else if(atomList(id)%type == 'O.co2') then
                        mark = .true.
                        exit INLOOP
                    else
                        cycle INLOOP
                    end if
                end do INLOOP !! END OF loop for one 2nd shell

                if(mark) then
                    exit OUTLOOP
                else
                    cycle OUTLOOP
                endif
            else
                cycle OUTLOOP
            end if
        end do OUTLOOP !! END of outer loop

        !! Assign if this group has amide bond or not
        if(mark .eqv. .false.) then
            this%amide = 0
        else
            this%amide = j !! Assign the value of the amide bond
        end if

    end subroutine xlogp_Group_CheckAmide


    subroutine xlogp_Group_SetValidity(this)
        type(Group_t), intent(inout) :: this

        this%valid = .true.
    end subroutine xlogp_Group_SetValidity

    subroutine xlogp_Group_FindAGroup(this, atomList, bondList, atomId)
        type(Group_t), intent(inout) :: this
        type(Atom_t), DIMENSION(:), target, intent(inout) :: atomList
        type(Bond_t), DIMENSION(:), target, intent(inout) :: bondList
        integer, intent(in) :: atomId

        ! Initialize group object
        call xlogp_Group_Clear(this)

        !> Insert Neighboring Atoms and Bonds into the group object
        !! The center atom of group will be also assigned as the atom
        !! corresponding to this%atoms(atomId)
        !!
        call xlogp_Group_SetNeighborInfo(this, atomList, bondList, atomId)

        !> Count number of hydrogen atoms and that of non-hydrogen atoms
        call xlogp_Group_CountNumHNonH(this)

        !> Count number of hetero atoms like as N, O
        call xlogp_Group_CountHeteroNum(this)

        !> Count the number of double and triple bonds
        call Group_CountDbTriBondNumType(this)

        !> Count aromatic (pi with delocalized eletrons) atoms
        call xlogp_Group_CountPiAtoms(this)

        !> Check whether group has neighboring amide bonds.
        call xlogp_Group_CheckAmide(this, atomList)

        !> Now, we achieved all the steps to gather group environment.
        !! So, this group is valid
        call xlogp_Group_SetValidity(this)

    end subroutine xlogp_Group_FindAGroup

    subroutine xlogp_Group_FindXGroup(this, atomList, bondList, atomId)
        type(Group_t), intent(inout) :: this
        type(Atom_t), DIMENSION(:), target, intent(inout) :: atomList
        type(Bond_t), DIMENSION(:), target, intent(inout) :: bondList
        integer, intent(in) :: atomId

        ! Initialize group object
        call xlogp_Group_Clear(this)

        !> Insert Neighboring Atoms and Bonds into the group object
        !! The center atom of group will be also assigned as the atom
        !! corresponding to this%atoms(atomId)
        !!
        call xlogp_Group_SetNeighborInfo(this, atomList, bondList, atomId)

        !> Count number of hydrogen atoms and that of non-hydrogen atoms
        call xlogp_Group_CountNumHNonH(this)

        !> Count number of hetero atoms like as N, O
        call xlogp_Group_CountHetNum_XLogp(this)

        !> Count aromatic (pi with delocalized eletrons) atoms
        call xlogp_Group_CountPiAtoms(this)

        call xlogp_Group_CountNumNarCar(this)

        !> Count the number of double and triple bonds
        call xlogp_Group_SetDblType_XLogp(this)

        !> Now, we achieved all the steps to gather group environment.
        !! So, this group is valid
        call xlogp_Group_SetValidity(this)

    end subroutine xlogp_Group_FindXGroup

end module xlogp_Group_m

