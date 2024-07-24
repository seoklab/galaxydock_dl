module xlogp_LogP_m
    use globals
    use xlogp_Atom_m
    use xlogp_Bond_m
    use xlogp_constants_m
    use xlogp_Group_m
    use xlogp_Molecule_m
    use xlogp_XLogPTyper_m
    use xlogp_XToolTyper_m
    use xlogp_Parameter_m
    use xlogp_Util_m
    use xlogp_Strings_m

    implicit none
    !private ! Information hiding

    REAL, parameter :: LOGP_HYDROPHOBIC_CARBON = 0.211
    REAL, parameter :: LOGP_INTERNAL_HBOND = 0.429
    REAL, parameter :: LOGP_HALOGEN_PAIR = 0.137
    REAL, parameter :: LOGP_NAR_PAIR = 0.485
    REAL, parameter :: LOGP_O3_PAIR = -0.268
    REAL, parameter :: LOGP_ACCEPTOR_PAIR = 0.580
    REAL, parameter :: LOGP_AMINO_ACID = -2.166
    REAL, parameter :: LOGP_SALICYLIC_ACID = 0.554
    REAL, parameter :: LOGP_SULFONIC_ACID = -0.501
    integer, parameter :: M_DP = SELECTED_REAL_KIND(8)

    type LogpFactor_t
        character(len=30) :: symbol != 'Unknown'
        real :: num
        real :: coeff
    end type LogpFactor_t

    type, public :: LogpEvaluate_t
        logical :: DEBUG != .false.
        type(LogpFactor_t) :: logpFactor(10)
        type(Molecule_t), pointer :: mol !=> null()
        type(ForceField_t), pointer :: forceField !=> null()
        integer, DIMENSION(XSCORE_MAX_LIGAND_ATOM_NUMBER) :: recordHbond
        logical :: init != .false.
        real(dp) :: logpValue != 0.0

    end type LogpEvaluate_t


contains

    subroutine Eval_SetDebug(this, debug)
        type(LogpEvaluate_t) :: this
        logical, intent(in) :: debug

        this%DEBUG = debug
    end subroutine Eval_SetDebug

    subroutine Eval_CalcLogP(this, flagLogpForEachAtom)
        type(LogpEvaluate_t) :: this
        !! This flag determines if correction factors are also distributed
        !! onto each atoms or not.
        logical, intent(in) :: flagLogpForEachAtom
        real(M_DP) :: xlogpTmp
        integer :: i
        type(Atom_t), pointer :: atoms(:)
        type(Bond_t), pointer :: bonds(:)
        character(len=20) :: xlogpType


        if(.not. this%init) then
            print*, '(*) You need to initialize Logp with molecule pointer'
            return
        endif

        atoms => this%mol%atoms
        bonds => this%mol%bonds
        do i=1, this%mol%num_atom
            if(.not. atoms(i)%valid) cycle
            call FindXLogPType(atoms, bonds, atoms(i)%id, xlogpType)
            call Param_Get_Atom_LogP(this%forceField, xlogpType, xlogpTmp)
            atoms(i)%xlogptype = xlogpType
            atoms(i)%logp = xlogpTmp

            if(this%DEBUG .AND. atoms(i)%xlogptype(1:1) /= 'H') then
                write(*, '(A, 1x, I3, 1x, A, F8.3, A, A)') &
                    "Atom [",i,"] logp value = ", atoms(i)%logp, &
                    " xlogptype = ", atoms(i)%xlogptype
            endif
        enddo

        do i=1, 10
            this%logpFactor(i)%num = 0
        enddo

        this%logpFactor(1)%num = CountHydrophobicCarbon(this, flagLogpForEachAtom)
        this%logpFactor(2)%num = CountInternalHBond(this, flagLogpForEachAtom)
        this%logpFactor(3)%num = CountHalogen_1_3_Pair(this, flagLogpForEachAtom)
        this%logpFactor(4)%num = CountNar_1_4_Pair(this, flagLogpForEachAtom)
        this%logpFactor(5)%num = CountO3_1_4_Pair(this, flagLogpForEachAtom)
        this%logpFactor(6)%num = CountAcceptor_1_5_Pair(this, flagLogpForEachAtom)
        this%logpFactor(7)%num = 0.0
        this%logpFactor(8)%num = CountAminoAcid(this, flagLogpForEachAtom)
        this%logpFactor(9)%num = CountSalicylicAcid(this, flagLogpForEachAtom)
        this%logpFactor(10)%num = CountSulfonicAcid(this, flagLogpForEachAtom)

    end subroutine Eval_CalcLogP

    subroutine Eval_FindAllXToolType(this)
        type(LogpEvaluate_t) :: this
        !! This flag determines if correction factors are also distributed
        !! onto each atoms or not.
        integer :: i
        type(Atom_t), pointer :: atoms(:)
        type(Bond_t), pointer :: bonds(:)
        character(len=20) :: xType


        if(.not. this%init) then
            print*, '(*) You need to initialize Logp with molecule pointer'
            return
        endif

        atoms => this%mol%atoms
        bonds => this%mol%bonds
        do i=1, this%mol%num_atom
            if(.not. atoms(i)%valid) cycle
            call FindXToolType(atoms, bonds, atoms(i)%id, xType)
            atoms(i)%xtype = xType
        enddo
    end subroutine Eval_FindAllXToolType

    subroutine Eval_FindAllXLogPType(this)
        type(LogpEvaluate_t) :: this
        !! This flag determines if correction factors are also distributed
        !! onto each atoms or not.
        integer :: i
        type(Atom_t), pointer :: atoms(:)
        type(Bond_t), pointer :: bonds(:)
        character(len=20) :: xlogpType


        if(.not. this%init) then
            print*, '(*) You need to initialize Logp with molecule pointer'
            return
        endif

        atoms => this%mol%atoms
        bonds => this%mol%bonds
        do i=1, this%mol%num_atom
            if(.not. atoms(i)%valid) cycle
            call FindXLogPType(atoms, bonds, atoms(i)%id, xlogpType)
            atoms(i)%xlogptype = xlogpType
        enddo
    end subroutine Eval_FindAllXLogPType

    subroutine Eval_InitCoeff(this)
        type(LogpEvaluate_t) :: this
        this%logpFactor(1)%symbol = "Hydrophobic carbon"
        this%logpFactor(2)%symbol = "Internal H-bond"
        this%logpFactor(3)%symbol = "Halogen 1-3 pair"
        this%logpFactor(4)%symbol = "Aromatic nitrogen 1-4 pair"
        this%logpFactor(5)%symbol = "Ortho sp3 oxygen pair"
        this%logpFactor(6)%symbol = "Acceptor 1-5 pair"
        this%logpFactor(7)%symbol = "Paralleled donor pair"
        this%logpFactor(8)%symbol = "Alpha amino acid"
        this%logpFactor(9)%symbol = "Salicylic acid"
        this%logpFactor(10)%symbol = "P-amino sulfnoic acid"

        this%logpFactor(1)%coeff = 0.211
        this%logpFactor(2)%coeff = 0.429
        this%logpFactor(3)%coeff = 0.137
        this%logpFactor(4)%coeff = 0.485
        this%logpFactor(5)%coeff = -0.268
        this%logpFactor(6)%coeff = 0.580
        this%logpFactor(7)%coeff = -0.423
        this%logpFactor(8)%coeff = -2.166
        this%logpFactor(9)%coeff = 0.554
        this%logpFactor(10)%coeff = -0.501
    end subroutine Eval_InitCoeff

    subroutine Eval_InitMolecule(this, molPtr)
        type(LogpEvaluate_t) :: this
        !        type(Molecule_t), intent(inout), pointer :: molPtr
        type(Molecule_t), pointer :: molPtr
        integer :: alloc_stat
        this%mol => molPtr

!        !> Prepare recordHbond
!        if(allocated(this%recordHbond)) then
!            if(size(this%recordHbond) < this%mol%num_atom) then
!                deallocate(this%recordHbond)
!                allocate(this%recordHbond(this%mol%num_atom), stat=alloc_stat)
!                if(alloc_stat /= 0) then
!                    print*, '(*) Allocation ERROR in SetMoleculePointer() = ', alloc_stat
!                endif
!            endif
!        else ! If not allocated
!            allocate(this%recordHbond(this%mol%num_atom), stat=alloc_stat)
!            if(alloc_stat /= 0) then
!                print*, '(*) Allocation ERROR in SetMoleculePointer() = ', alloc_stat
!            endif
!        endif

        call Eval_InitCoeff(this)

        this%init = .true.
    end subroutine Eval_InitMolecule

    subroutine Eval_InitForceField(this, forceFieldPtr)
        type(LogpEvaluate_t) :: this
        !        type(ForceField_t), intent(inout), pointer :: forceFieldPtr
        type(ForceField_t), pointer :: forceFieldPtr
        this%forceField => forceFieldPtr
    end subroutine Eval_InitForceField

    subroutine Eval_Destroy(this)
        type(LogpEvaluate_t) :: this
!        if(allocated(this%recordHbond)) then
!            deallocate(this%recordHbond)
!        endif
    end subroutine Eval_Destroy

    function Connection_1_2_Check(this, id1, id2) result(return_value)
        type(LogpEvaluate_t) :: this
        integer, intent(in) :: id1, id2
        logical :: return_value
        integer :: i
        if (id1 == id2) then
            return_value = .false.
            return
        endif

        !> If one neighbor id of atom2 is equal to id1, then it is 1,2 bond.
        !        associate (atoms => this%mol%atoms)
        do i=1, this%mol%atoms(id2)%num_neib
            if (id1 == this%mol%atoms(id2)%neib(i)) then
                return_value = .true.
                return
            else
                cycle
            endif
        enddo
        !        end associate
        return_value = .false.
    end function Connection_1_2_Check

    function Connection_1_3_Check(this, id1, id2) result(return_value)
        type(LogpEvaluate_t) :: this
        integer, intent(in) :: id1, id2
        logical :: return_value
        integer :: i, j

        if (id1 == id2) then
            return_value = .false.
            return
        endif
        if (Connection_1_2_Check(this, id1, id2)) then
            return_value = .false.
            return
        endif

        !        associate (atoms => this%mol%atoms)
        do i=1, this%mol%atoms(id1)%num_neib
            do j=1, this%mol%atoms(id2)%num_neib
                if (this%mol%atoms(id1)%neib(i) == this%mol%atoms(id2)%neib(j)) then
                    return_value = .true.
                    return
                else
                    cycle
                endif
            enddo
        enddo
        !        end associate
        return_value = .false.
    end function Connection_1_3_Check

    function Connection_1_4_Check(this, id1, id2) result(return_value)
        type(LogpEvaluate_t) :: this
        integer, intent(in) :: id1, id2
        logical :: return_value
        integer :: i, j

        if (id1 == id2) then
            return_value = .false.
            return
        endif
        if (Connection_1_2_Check(this, id1, id2)) then
            return_value = .false.
            return
        endif
        if (Connection_1_3_Check(this, id1, id2)) then
            return_value = .false.
            return
        endif

!        associate (atoms => this%mol%atoms)
            do i=1, this%mol%atoms(id1)%num_neib
                do j=1, this%mol%atoms(id2)%num_neib
                    if (Connection_1_2_Check(this, this%mol%atoms(id1)%neib(i), &
                        this%mol%atoms(id2)%neib(j))) then
                        return_value = .true.
                        return
                    else
                        cycle
                    endif
                enddo
            enddo
!        end associate
        return_value = .false.
    end function Connection_1_4_Check

    function Connection_1_5_Check(this, id1, id2) result(return_value)
        type(LogpEvaluate_t) :: this
        integer, intent(in) :: id1, id2
        logical :: return_value
        integer :: i, j

        if (id1 == id2) then
            return_value = .false.
            return
        endif
        if (Connection_1_2_Check(this, id1, id2)) then
            return_value = .false.
            return
        endif
        if (Connection_1_3_Check(this, id1, id2)) then
            return_value = .false.
            return
        endif
        if (Connection_1_4_Check(this, id1, id2)) then
            return_value = .false.
            return
        endif

!        associate (atoms => this%mol%atoms)
            do i=1, this%mol%atoms(id1)%num_neib
                do j=1, this%mol%atoms(id2)%num_neib
                    if (Connection_1_3_Check(this, this%mol%atoms(id1)%neib(i), &
                        this%mol%atoms(id2)%neib(j))) then
                        return_value = .true.
                        return
                    else
                        cycle
                    endif
                enddo
            enddo
!        end associate
        return_value = .false.
    end function Connection_1_5_Check

    function Connection_1_6_Check(this, id1, id2) result(return_value)
        type(LogpEvaluate_t) :: this
        integer, intent(in) :: id1, id2
        logical :: return_value
        integer :: i, j

        if (id1 == id2) then
            return_value = .false.
            return
        endif
        if (Connection_1_2_Check(this, id1, id2)) then
            return_value = .false.
            return
        endif
        if (Connection_1_3_Check(this, id1, id2)) then
            return_value = .false.
            return
        endif
        if (Connection_1_4_Check(this, id1, id2)) then
            return_value = .false.
            return
        endif
        if (Connection_1_5_Check(this, id1, id2)) then
            return_value = .false.
            return
        endif

!        associate (atoms => this%mol%atoms)
            do i=1, this%mol%atoms(id1)%num_neib
                do j=1, this%mol%atoms(id2)%num_neib
                    if (Connection_1_4_Check(this, this%mol%atoms(id1)%neib(i),&
                        this%mol%atoms(id2)%neib(j))) then
                        return_value = .true.
                        return
                    else
                        cycle
                    endif
                enddo
            enddo
!        end associate
        return_value = .false.
    end function Connection_1_6_Check

    !!----------------------------------------------------------------------------
    !> The author of XLOGP observed that the log P values of compounds with hydrocarbon
    !! chains are often underestimated by atom addition. Such compounds tend to be more
    !! flexible and easier to aggregate in the aqueous phase.
    !!
    !! They define an sp3 or sp2 hybridized carbon atom as a 'hydrophobic carbon atom' if
    !! there is no heteroatom (any atom other than carbon) within the 1-4 relationship.
    !!
    !! Such a carbon atom locates in a hydrophobic micro-environment and
    !! thus its hydrophobicity is reinforced.
    !! The total number of hydrophobic carbon atoms in the given compound is counted.
    !!
    !! This definition of hydrophobic carbons does not include aromatic carbon atoms.
    !!
    !! If the number of hydrophobic carbons is greater than 10, we divide the number by 2.
    !!----------------------------------------------------------------------------
    function CountHydrophobicCarbon(this, flag) result(return_value)
        type(LogpEvaluate_t) :: this
        logical, intent(in) :: flag
        integer :: return_value
        integer :: i
        return_value = 0
        do i=1, this%mol%num_atom
!            associate (atom => this%mol%atoms(i))
                if ((this%mol%atoms(i)%type /= 'C.3') .AND. &
                    (this%mol%atoms(i)%type /= 'C.2') ) then
                    cycle
                else if(.not. HydrophobicNeighborCheck(this, this%mol%atoms(i)%id)) then
                    cycle
                endif
                return_value = return_value + 1
                if(flag) then
                    this%mol%atoms(i)%logp = this%mol%atoms(i)%logp + LOGP_HYDROPHOBIC_CARBON
                else
                    cycle
                endif
!            end associate
        enddo

        if (return_value .GE. 10) then
            return_value = return_value / 2
        endif
    end function CountHydrophobicCarbon

    function HydrophobicNeighborCheck(this, id) result(return_value)
        type(LogpEvaluate_t) :: this
        integer, intent(in) :: id
        logical :: return_value
        integer :: i

        return_value = .true.
        HP: do i=1, this%mol%num_atom
!            associate(atoms => this%mol%atoms)
                if(i==id) then
                    cycle
                else if(this%mol%atoms(i)%type == 'F' .OR. &
                    this%mol%atoms(i)%type == 'Cl' .OR. &
                    this%mol%atoms(i)%type == 'Br' .OR. &
                    this%mol%atoms(i)%type == 'I' .OR. &
                    this%mol%atoms(i)%type == 'Si' .OR. &
                    this%mol%atoms(i)%type(1:1) == 'N' .OR. &
                    this%mol%atoms(i)%type(1:1) == 'O' .OR. &
                    this%mol%atoms(i)%type(1:1) == 'S' .OR. &
                    this%mol%atoms(i)%type(1:1) == 'P') then
                    if (Connection_1_2_Check(this, this%mol%atoms(id)%id, this%mol%atoms(i)%id)) then
                        return_value = .false.
                        EXIT HP
                    else if (Connection_1_3_Check(this, this%mol%atoms(id)%id, this%mol%atoms(i)%id)) then
                        return_value = .false.
                        EXIT HP
                    else if (Connection_1_4_Check(this, this%mol%atoms(id)%id, this%mol%atoms(i)%id)) then
                        return_value = .false.
                        EXIT HP
                    else
                        cycle
                    endif
                else
                    cycle
                endif
!            end associate
        enddo HP

    end function HydrophobicNeighborCheck

    !!----------------------------------------------------------------------------
    !> This will estimate the occurrence of the internal hydrogen bonding.
    !! Defintion to detect an internal hydrogen bond.
    !!
    !! (i) The donor atom could be any nitrogen or oxygen atom with at least
    !!      one hydrogen atom attached, while the acceptor atom could be any
    !!      sp2 oxygen atom or sp3 oxygen atom in a hydroxyl group.
    !!
    !! (ii) Either the donor or the acceptor atom should be linked directly to a ring.
    !!      The ring serves to immobilize the orientations of the donor and the acceptor.
    !!
    !! (iii) If both the donor and the acceptor are linked to a ring, they should be
    !!      1-4 relationship.
    !!
    !! (iv) If only the donor or the acceptor is linked to a ring, they should be of
    !!      1-5 relationship.
    !!
    !!  By using above definitions, only 'reliable' internal hydrogen bonds are taken into
    !!      account.
    !!
    !! The object component recordHbond will be used during this function.
    !!       H
    !!      /
    !!     O
    !!     |   O
    !!    / \//
    !!   |   |
    !!    \ /
    !!----------------------------------------------------------------------------
    function CountInternalHBond(this, flagLogpSum) result(return_value)
        type(LogpEvaluate_t) :: this
        logical, intent(in) :: flagLogpSum
        real :: return_value
        integer :: i, j
        logical :: mark1, mark2
        type(Atom_t), pointer :: atoms(:)
        integer :: k

        this%recordHbond = 0
        return_value = 0.0
        atoms => this%mol%atoms

        DONOR_LOOP: do i=1, this%mol%num_atom
            !> We only consider hydrogen bond donors.
            if(atoms(i)%hb /= 'D' .AND. atoms(i)%hb /= 'DA') then
                cycle
            else if(atoms(i)%ring /= 0) then !! Not allowed in ring
                cycle
            endif

            if(.not. AdjacentRingCheck(this, atoms(i)%id)) then
                mark1 = .false.
            else
                mark1 = .true.
            endif

            ACCEPTOR_LOOP : do j=1, this%mol%num_atom
                !> Now we will find hydrogen bond acceptor only.
                if(i==j) then
                    cycle
                else if(atoms(j)%hb /= 'A' .AND. atoms(j)%hb /= 'DA') then
                    cycle
                else if(atoms(j)%type == 'O.3' .AND. atoms(j)%hb == 'A') then
                    cycle
                else if(atoms(j)%type == 'N.2') then
                    cycle
                else if(atoms(j)%type == 'N.ar') then
                    cycle
                else if(atoms(j)%ring /= 0) then
                    cycle !! It should be not in the ring.
                endif

                if(.not. AdjacentRingCheck(this, atoms(j)%id)) then
                    mark2 = .false.
                else
                    mark2 = .true.
                endif

                if(mark1 .AND. mark2) then
                    if(.not. Connection_1_4_Check(this, atoms(i)%id, atoms(j)%id)) then
                        cycle
                    else
                        this%recordHbond(i) = this%recordHbond(i) + 1
                        this%recordHbond(j) = this%recordHbond(j) + 1
                    endif
                else if(mark1 .AND. .not. mark2) then
                    if(.not. Connection_1_5_Check(this, atoms(i)%id, atoms(j)%id)) then
                        cycle
                    else
                        this%recordHbond(i) = this%recordHbond(i) + 1
                        this%recordHbond(j) = this%recordHbond(j) + 1
                    endif
                else if(.not. mark1 .AND. mark2) then
                    if(.not. Connection_1_5_Check(this, atoms(i)%id, atoms(j)%id)) then
                        cycle
                    else
                        this%recordHbond(i) = this%recordHbond(i) + 1
                        this%recordHbond(j) = this%recordHbond(j) + 1
                    endif
                else
                    cycle
                endif
            enddo ACCEPTOR_LOOP
        enddo DONOR_LOOP

        do i=1, this%mol%num_atom
            if(this%recordHbond(i)==0) cycle
            return_value = return_value + 0.5
            if(flagLogpSum) atoms(i)%logp = atoms(i)%logp + (LOGP_INTERNAL_HBOND / 2.0)
        enddo

    end function CountInternalHBond

    function AdjacentRingCheck(this, atomId) result(return_value)
        type(LogpEvaluate_t) :: this
        integer, intent(in) :: atomId
        logical :: return_value ! Return .true. or .false.
        integer :: i, numNeighbor, neighborAtomId
        type(Atom_t), pointer :: atoms(:)
        atoms => this%mol%atoms

        return_value = .false.
        numNeighbor = atoms(atomId)%num_neib
        do i=1, numNeighbor
            neighborAtomId=atoms(atomId)%neib(i)
            if(atoms(neighborAtomId)%ring /= 0) then
                return_value = .true. ! We see aromartic ring neighborhood.
                EXIT
            endif
        enddo
    end function AdjacentRingCheck

    function CountHalogen_1_3_Pair(this, flag) result(return_value)
        type(LogpEvaluate_t) :: this
        logical, intent(in) :: flag
        real :: return_value
        integer :: i, j
        integer :: num1, num2
        type(Atom_t), pointer :: atoms(:)
        atoms => this%mol%atoms
        num1=0
        num2=0

        do i=1, this%mol%num_atom-1
            if(atoms(i)%type /= 'F') cycle

            do j=i+1, this%mol%num_atom
                if(atoms(j)%type /= 'F') then
                    cycle
                else if(.not. Connection_1_3_Check(this, atoms(i)%id, atoms(j)%id)) then
                    cycle
                endif

                num1 = num1 + 1

                if(flag) then
                    atoms(i)%logp = atoms(i)%logp + ((LOGP_HALOGEN_PAIR)/2.0)
                    atoms(j)%logp = atoms(j)%logp + ((LOGP_HALOGEN_PAIR)/2.0)
                else
                    cycle
                endif
            enddo
        enddo

        do i=1, this%mol%num_atom-1
            if(atoms(i)%type /= 'Cl' .AND. &
                atoms(i)%type /= 'Br' .AND. &
                atoms(i)%type /= 'I') then
                cycle
            endif


            do j=i+1, this%mol%num_atom
                if(atoms(j)%type /= 'Cl' .AND. &
                    atoms(j)%type /= 'Br' .AND. &
                    atoms(j)%type /= 'I') then
                    cycle
                else if(.not. Connection_1_3_Check(this, atoms(i)%id, atoms(j)%id)) then
                    cycle
                endif

                num2 = num2 + 1

                if(flag) then
                    atoms(i)%logp = atoms(i)%logp + ((LOGP_HALOGEN_PAIR)/2.0)
                    atoms(j)%logp = atoms(j)%logp + ((LOGP_HALOGEN_PAIR)/2.0)
                else
                    cycle
                endif
            enddo
        enddo

        return_value = num1+num2

    end function CountHalogen_1_3_Pair

    !!----------------------------------------------------------------------------
    !> This correction factor is triggered when [two aromatic nitrogen atoms] are
    !! in the same aromatic ring and separated by two other atoms.
    !!    N
    !!   / \\
    !! ||   |
    !!   \N//
    !!----------------------------------------------------------------------------
    function CountNar_1_4_Pair(this, flag) result(return_value)
        type(LogpEvaluate_t) :: this
        logical, intent(in) :: flag
        real :: return_value
        integer :: i, j
        integer :: tmp1, tmp2, tmp3, tmp4
        type(Atom_t), pointer :: atoms(:)
        atoms => this%mol%atoms

        return_value = 0.0

        do i=1, this%mol%num_atom
            if(atoms(i)%type /= 'N.ar') cycle
            tmp1 = atoms(i)%neib(1)
            tmp2 = atoms(i)%neib(2)
            do j=i+1, this%mol%num_atom
                if(atoms(j)%type /= 'N.ar') then
                    cycle
                else if(.not. Connection_1_4_Check(this, atoms(i)%id, atoms(j)%id)) then
                    cycle
                else
                    tmp3 = atoms(j)%neib(1)
                    tmp4 = atoms(j)%neib(2)
                    !! Below two 1_2_Check will make sure that N1, N2 is connected by 1_4 in both directions.
                    if(Connection_1_2_Check(this, tmp1, tmp3)) then
                        if(Connection_1_2_Check(this, tmp2, tmp4)) then
                            return_value = return_value + 1
                            if(flag) then
                                atoms(i)%logp = atoms(i)%logp + (LOGP_NAR_PAIR/2.0)
                                atoms(j)%logp = atoms(j)%logp + (LOGP_NAR_PAIR/2.0)
                            endif
                        else
                            cycle
                        endif
                    else if(Connection_1_2_Check(this, tmp1, tmp4)) then
                        if(Connection_1_2_Check(this, tmp2, tmp3)) then
                            return_value = return_value + 1
                            if(flag) then
                                atoms(i)%logp = atoms(i)%logp + (LOGP_NAR_PAIR/2.0)
                                atoms(j)%logp = atoms(j)%logp + (LOGP_NAR_PAIR/2.0)
                            endif
                        else
                            cycle
                        endif
                    else
                        cycle
                    endif
                endif
            enddo
        enddo
    end function CountNar_1_4_Pair

    !!----------------------------------------------------------------------------
    !> Ortho sp3 oxygen pair : This correction factor is used for compounds as next
    !!      OR
    !!   /\\/
    !! ||  |
    !!   \//\
    !!       OR
    !!----------------------------------------------------------------------------
    function CountO3_1_4_Pair(this, flag) result(return_value)
        type(LogpEvaluate_t) :: this
        logical, intent(in) :: flag
        real :: return_value
        integer :: i, j
        integer :: tmp1, tmp2, tmp3, tmp4
        type(Atom_t), pointer :: atoms(:)
        atoms => this%mol%atoms
        return_value = 0.0

        do i=1, this%mol%num_atom
            this%recordHbond(i) = 0
        enddo

        do i=1, this%mol%num_atom-1
            if(atoms(i)%type /= 'O.3') then
                cycle
            else if(atoms(i)%hb == 'DA') then
                cycle
            else if(atoms(i)%ring /= 0) then
                cycle
            else if(.not. AdjacentAromaticCheck(this, atoms(i)%id)) then
                cycle
            endif
            do j=i+1, this%mol%num_atom
                if(atoms(j)%type /= 'O.3') then
                    cycle
                else if(atoms(j)%hb == 'DA') then
                    cycle
                else if(atoms(j)%ring /= 0) then
                    cycle
                else if(.not. AdjacentAromaticCheck(this, atoms(j)%id)) then
                    cycle
                else if(.not. Connection_1_4_Check(this, atoms(i)%id, atoms(j)%id)) then
                    cycle
                else
                    this%recordHbond(i) = this%recordHbond(i) + 1
                    this%recordHbond(j) = this%recordHbond(j) + 1
                endif
            enddo
        enddo

        do i=1, this%mol%num_atom
            if(this%recordHbond(i)==0) cycle
            return_value = return_value + 0.500
            if(flag) then
                atoms(i)%logp = atoms(i)%logp + (LOGP_O3_PAIR/2.0)
            else
                cycle
            endif
        enddo
    end function CountO3_1_4_Pair

    !!----------------------------------------------------------------------------
    !> sp2 oxygen 1-5 pair : It is not counted if such a sub-structure is
    !! part of a ring. LOGP_ACCEPTOR_PAIR (Acceptor 1-5 pair).
    !!    O   O
    !!    ||  ||
    !!   / \ / \
    !!  R   O   R
    !!----------------------------------------------------------------------------
    function CountAcceptor_1_5_Pair(this, flag) result(return_value)
        type(LogpEvaluate_t) :: this
        logical, intent(in) :: flag
        real :: return_value
        integer :: i, j
        integer :: tmp1, tmp2, tmp3, tmp4
        type(Atom_t), pointer :: atoms(:)
        atoms => this%mol%atoms
        return_value = 0.0

        do i=1, this%mol%num_atom-1
            if(atoms(i)%hb /= 'A') then ! Hydrogen bond type must be an acceptor ('A') and atom type be a sp2-oxygen
                cycle
            else if(atoms(i)%type == 'O.3') then
                cycle
            else if(atoms(i)%type == 'N.2') then
                cycle
            else if(atoms(i)%type == 'N.ar') then
                cycle
            endif
            tmp1 = atoms(i)%neib(1)
            if(atoms(tmp1)%type == 'S') then
                cycle
            else if(atoms(tmp1)%type == 'P') then
                cycle
            else if(atoms(tmp1)%ring /= 0) then
                cycle
            endif

            do j=i+1, this%mol%num_atom
                if(atoms(j)%hb /= 'A') then !! Only consider sp2-oxygen acceptor
                    cycle
                else if(atoms(j)%type == 'O.3') then
                    cycle
                else if(atoms(j)%type == 'N.2') then
                    cycle
                else if(atoms(j)%type == 'N.ar') then
                    cycle
                endif
                tmp2 = atoms(j)%neib(1)
                if(atoms(tmp2)%type == 'S') then
                    cycle
                else if(atoms(tmp2)%type == 'P') then
                    cycle
                else if(atoms(tmp2)%ring /= 0) then
                    cycle
                endif

                if(.not. Connection_1_5_Check(this, atoms(i)%id, atoms(j)%id) ) cycle
                return_value = return_value + 1

                if(flag) then
                    atoms(i)%logp = atoms(i)%logp + (LOGP_ACCEPTOR_PAIR/2.0)
                    atoms(j)%logp = atoms(j)%logp + (LOGP_ACCEPTOR_PAIR/2.0)
                else
                    cycle
                endif
            enddo
        enddo
    end function CountAcceptor_1_5_Pair

    !!----------------------------------------------------------------------------
    !> "Salicylic acid and its derivatives are generally more hydrophobic than
    !! it was predicted by atom addition alone. Even after we add the correction
    !! of the internal hydrogen bond, the log P values of such compounds are
    !! still underestimated.
    !!
    !! It seems that the internal hydrogen bond existing in such a compound is much
    !! stronger than an 'average' for salicylic acid and its derivatives.
    !! The indicator could be 0 or 1.
    !!
    !! Reference : Wang,R. et al. (2000) Calculating partition coefficient by
    !! atom-additive method. Perspectives in Drug Discovery and Design, 19, 47-66.
    !!
    !! @TODO : I am not sure if this function is corrrectly implemented or not (HWCHUNG).
    !!----------------------------------------------------------------------------
    function CountSalicylicAcid(this, flag) result(return_value)
        type(LogpEvaluate_t) :: this
        logical, intent(in) :: flag
        real :: return_value
        integer :: i, j
        integer :: tmp
        logical :: mark
        type(Atom_t), pointer :: atoms(:)
        atoms => this%mol%atoms
        return_value = 0.0

        OUTER_LOOP : do i=1, this%mol%num_atom
            if(atoms(i)%type /= 'O.2') cycle
            tmp=atoms(i)%neib(1)
            if(atoms(tmp)%type(1:1) /= 'C') then
                cycle
            else if(atoms(tmp)%ring /= 0) then
                cycle
            else if(.not. AdjacentAromaticCheck(this, atoms(tmp)%id) ) then
                cycle
            endif

            mark = .false.
            do j=1, this%mol%num_atom
                if(i==j) then
                    cycle
                else if(atoms(j)%type /= 'O.3') then
                    cycle
                else if(atoms(j)%ring /= 0) then
                    cycle
                else if(.not. Connection_1_3_Check(this, atoms(i)%id, atoms(j)%id) ) then
                    cycle
                else
                    mark = .true.
                    EXIT
                endif
            enddo

            if(.not. mark) cycle

            mark = .false.

            INNER_LOOP: do j=1, this%mol%num_atom
                if(i==j) then
                    cycle
                else if(atoms(j)%type /= 'O.3') then
                    cycle
                else if(atoms(j)%hb /= 'DA') then
                    cycle
                else if(.not. AdjacentAromaticCheck(this, atoms(j)%id)) then
                    cycle
                else if(.not. Connection_1_5_Check(this, atoms(i)%id, atoms(j)%id)) then
                    cycle
                endif

                return_value = return_value + 1

                if(flag) then
                    atoms(i)%logp = atoms(i)%logp + (LOGP_SALICYLIC_ACID/2.0)
                    atoms(j)%logp = atoms(j)%logp + (LOGP_SALICYLIC_ACID/2.0)
                endif
                EXIT INNER_LOOP !! Why we got exit here?
            enddo INNER_LOOP

            if(return_value /= 0) EXIT OUTER_LOOP
        enddo OUTER_LOOP

    end function CountSalicylicAcid

    function AdjacentAromaticCheck(this,atomId) result(return_value)
        type(LogpEvaluate_t) :: this
        integer :: atomId
        logical :: return_value
        integer :: numNeighbor, neighborId
        type(Atom_t), pointer :: atoms(:)
        integer :: i
        atoms => this%mol%atoms

        return_value = .false.
        numNeighbor = atoms(atomId)%num_neib

        do i=1, numNeighbor
            neighborId=atoms(atomId)%neib(i)
            if(atoms(neighborId)%type == 'C.ar') then
                return_value = .true.
                EXIT
            else
                cycle
            endif
        enddo
    end function AdjacentAromaticCheck

    function AdjacentAromaticId(this,atomId) result(return_value)
        type(LogpEvaluate_t) :: this
        integer :: atomId
        integer :: return_value
        integer :: numNeighbor, neighborId
        type(Atom_t), pointer :: atoms(:)
        integer :: i
        atoms => this%mol%atoms

        return_value = 0
        numNeighbor = atoms(atomId)%num_neib

        do i=1, numNeighbor
            neighborId=atoms(atomId)%neib(i)
            if(atoms(neighborId)%type == 'C.ar') then
                return_value = neighborId
                EXIT
            else
                cycle
            endif
        enddo
    end function AdjacentAromaticId

    !!----------------------------------------------------------------------------
    !> Indicator for alpha-amino acid: Alpha-amino acids do not have free amino
    !! and carboxylic acid groups but rather exist as zwitterions.
    !! The log P value of alpha-amino acid will be largely overestimated
    !! by atom addition alone. Therefore we us a special indicator variable for
    !! alpha-amino acids. 0 or 1.
    !!----------------------------------------------------------------------------
    function CountAminoAcid(this, flag) result(return_value)
        type(LogpEvaluate_t) :: this
        logical, intent(in) :: flag
        real :: return_value
        integer :: i, j
        integer :: tmp
        logical :: mark
        type(Atom_t), pointer :: atoms(:)
        atoms => this%mol%atoms
        return_value = 0.0

        OUT_LOOP : do i=1, this%mol%num_atom
            if(atoms(i)%type /= 'O.2') cycle

            tmp=atoms(i)%neib(1)
            if(atoms(tmp)%type(1:1) /= 'C') then
                cycle
            else if(atoms(tmp)%ring /= 0) then
                cycle
            endif

            mark = .false.

            IN_LOOP: do j=1, this%mol%num_atom
                if(i==j) then
                    cycle
                else if(atoms(j)%type /= 'O.3') then
                    cycle
                else if(atoms(j)%hb /= 'DA') then
                    cycle
                else if(atoms(j)%ring /= 0) then
                    cycle
                else if(.not. Connection_1_3_Check(this, atoms(i)%id, atoms(j)%id)) then
                    cycle
                else
                    mark = .true.
                    EXIT IN_LOOP
                endif
            enddo IN_LOOP

            if(.not. mark) cycle

            do j=1, this%mol%num_atom
                if(atoms(j)%xtype /= 'N.3.h2.pi=0') then
                    cycle
                else if(.not. Connection_1_4_Check(this, atoms(i)%id, atoms(j)%id)) then
                    cycle
                endif

                return_value = return_value + 1

                if(flag) then
                    atoms(i)%logp = atoms(i)%logp + LOGP_AMINO_ACID / 2.0
                    atoms(j)%logp = atoms(j)%logp + LOGP_AMINO_ACID / 2.0
                endif
                EXIT
            enddo

            if(return_value /= 0) EXIT OUT_LOOP

            do j=1, this%mol%num_atom
                if(atoms(j)%type /= 'N.ar') then
                    cycle
                else if(.not. Connection_1_4_Check(this, atoms(i)%id, atoms(j)%id)) then
                    cycle
                endif

                return_value = return_value + 1.0
                if(flag) then
                    atoms(i)%logp = atoms(i)%logp + LOGP_AMINO_ACID / 2.0
                    atoms(j)%logp = atoms(j)%logp + LOGP_AMINO_ACID / 2.0
                endif
                EXIT
            enddo

            if(return_value /= 0) EXIT OUT_LOOP

        enddo OUT_LOOP
    end function CountAminoAcid

    function CountSulfonicAcid(this, flag) result(return_value)
        type(LogpEvaluate_t) :: this
        logical, intent(in) :: flag
        real :: return_value
        integer :: i,j,tmp1,tmp2
        type(Atom_t), pointer :: atoms(:)
        atoms => this%mol%atoms
        return_value = 0

        do i=1, this%mol%num_atom
            if(atoms(i)%type /= 'S.o2') then
                cycle
            else if(atoms(i)%ring /= 0) then
                cycle
            endif

            tmp1 = AdjacentAromaticId(this, atoms(i)%id)
            if(tmp1 == 0) cycle

            do j=1, this%mol%num_atom
                if(atoms(j)%xlogptype /= 'N.3.h2.pi=1') then
                    cycle
                else

                endif

                tmp2 = AdjacentAromaticId(this, atoms(j)%id)
                if(tmp2 == 0) cycle

                if(.not. Connection_1_6_Check(this, atoms(i)%id,atoms(j)%id)) then
                    cycle
                else if(.not. Connection_1_4_Check(this, tmp1, tmp2)) then
                    cycle
                endif
                return_value = return_value + 1
                if(flag) then
                    atoms(i)%logp = atoms(i)%logp + (LOGP_SULFONIC_ACID/2.0)
                    atoms(j)%logp = atoms(j)%logp + (LOGP_SULFONIC_ACID/2.0)
                endif
            enddo
        enddo
    end function CountSulfonicAcid

    function InnerCollisionCheck(this) result(return_value)
        type(LogpEvaluate_t) :: this
        logical :: flag
        logical :: return_value
        integer :: i,j,tmp1,tmp2
        type(Atom_t), pointer :: atoms(:)
        real(dp) :: d, d0
        real(dp) :: dist(3)

        atoms => this%mol%atoms
        return_value = .false.

        do i=1, this%mol%num_atom
            do j=1, this%mol%num_atom
                dist = atoms(i)%coor - atoms(j)%coor
                d = dot_product(dist, dist)
                d0 = atoms(i)%r + atoms(j)%r
                if (d > (d0 - 0.50)) then
                    cycle !! No bump
                else if ( Connection_1_2_Check(this, atoms(i)%id, atoms(j)%id)) then
                    cycle
                else if ( Connection_1_3_Check(this, atoms(i)%id, atoms(j)%id)) then
                    cycle
                else if ( Connection_1_4_Check(this, atoms(i)%id, atoms(j)%id)) then
                    cycle
                else
                    return_value = .true.
                    EXIT
                endif
            enddo
        enddo

    end function InnerCollisionCheck

    subroutine Eval_WriteOutLog(this, fileName)
        implicit none
        type(LogpEvaluate_t) :: this
        character(len=*), intent(in) :: fileName
        integer :: i
        logical :: mark
        integer :: ios
        integer :: unit_no
        character(len=120) :: str1, str2, str3, fmt1, fmt2, fmt3
        integer(4) :: s, h, d, m, y

        call GetUnitNo(unit_no) !! open
        !!----------------------------------------------------------------------------

        open(unit=unit_no, FILE=TRIM(fileName), ACTION='write', iostat=ios)
        if(ios /= 0) then
            print*, '(*) FILE open error at WriteOutLog(). Error code = ', ios
        end if

        write(unit_no, '("#")')
        call GetTimeStamp2(s, h, d, m, y)
        write(str1, '("Thu May ", I2, " ", I2.2, ":", I2.2, ":", I2.2, " ", I4)') d, h, m, s, y

        write(unit_no, '(A, A)') "# LogP calculation for GALAXYDOCK: ", TRIM(str1)
        write(unit_no, '("#")')
        write(unit_no, '(A, A)') "# Name of molecule : ", TRIM(this%mol%name)
        write(unit_no, '(A, A)') "# Molecular formula: ", TRIM(this%mol%formula)
        write(unit_no, '(A, F9.2)') "# Molecular weight: ", this%mol%weight
        write(unit_no, '("#")')

        write(unit_no, '("----------------------------------------------")')
        write(unit_no, '("ID   atom type                  contribution")')
        write(unit_no, '("----------------------------------------------")')

        do i=1, this%mol%num_atom
            if(this%mol%atoms(i)%type == 'H') cycle
            call writenum(this%mol%atoms(i)%id, str1, "I3")
            write(unit_no, '(A5, A10, A13, F10.3)') &
                adjustl(str1), &
                this%mol%atoms(i)%xtype, &
                "", &
                this%mol%atoms(i)%logp
        enddo
        write(unit_no, '("----------------------------------------------")')
        mark = .false.
        do i=1, 10
            if(this%logpFactor(i)%num==0) cycle
            write(unit_no, '(A30,F9.3)') &
                this%logpFactor(i)%symbol, &
                this%logpFactor(i)%num*this%logpFactor(i)%coeff
            mark = .true.
        enddo

        if(mark) then
            write(unit_no, '("----------------------------------------------")')
        endif

        write(unit_no, '(A, F8.2)') "Calculated LogP =", this%logpValue
        write(unit_no, '("")')

        !!----------------------------------------------------------------------------
        close(unit_no) !! close
    end subroutine Eval_WriteOutLog

end module xlogp_LogP_m
