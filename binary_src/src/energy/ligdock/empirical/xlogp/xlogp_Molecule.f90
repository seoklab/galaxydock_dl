module xlogp_Molecule_m
    use globals, only : dp, ref_lig_type, residue_type, molecule_type
    use xlogp_Atom_m
    use xlogp_Bond_m
    use xlogp_Ring_m
    use xlogp_Group_m
    use xlogp_constants_m
    use xlogp_Mol2Reader_m, only : Mol2SimpleMolecule_t
    use xlogp_Util_m
    use xlogp_Strings_m

    !------------------------------------------------------------
    ! Molecule_t Derived Type
    type Molecule_t

        !sequence ! Let this Molecule type stored as sequential
        character(len=256) :: name != 'Unknown molecule'
        character(len=256) :: formula != 'Unknown formula'
        real :: weight != 0.0
        integer :: num_hb_atom != 0
        integer :: num_rotor != 0
        integer :: num_atom != 0
        type(Atom_t), DIMENSION(XSCORE_MAX_LIGAND_ATOM_NUMBER) ::  atoms ! Allocatable atoms container

        integer :: num_bond != 0
        type(Bond_t), DIMENSION(XSCORE_MAX_LIGAND_BOND_NUMBER) :: bonds ! Allocatable bonds container

        integer :: num_subst != 0
        integer :: num_features != 0
        integer :: num_set != 0

        character(len=256) :: mol_type != 'Unknown mol type'
        character(len=256) :: charge_type != 'Unknown charge type'

        integer :: num_ring != 0
        type(Ring_t), DIMENSION(XSCORE_MAX_RING_NUMBER) :: rings ! Allocatable rings container

    end type Molecule_t

    interface Molecule
        module procedure Mol_Construct
    end interface !Molecule

contains

    function Mol_Construct(atom_num, bond_num) result(return_value)
        integer, optional :: atom_num
        integer, optional :: bond_num
        type(Molecule_t) :: return_value
        integer :: atom_num2, bond_num2
        if(present(atom_num)) then
            atom_num2 = atom_num
        else
            atom_num2 = 0
        endif

        if(present(bond_num)) then
            bond_num2 = bond_num
        else
            bond_num2 = 0
        endif

        call Mol_InitWithNumber(return_value, atom_num2, bond_num2)
    end function Mol_Construct

    subroutine Mol_InitWithNoItem(this)
        !implicit none
        !type(Molecule_t), intent(inout) :: this
        type(Molecule_t) this

        print*, 'initWithNoItem was called .... '

!        allocate(this%atoms(0))
!        allocate(this%bonds(0))
!        allocate(this%rings(0))
        this%name = 'Unknown'
        this%formula = 'Unknown'
        this%weight = 0.0
        this%num_hb_atom = 0
        this%num_rotor = 0
        this%num_atom = 0
        this%num_bond = 0
        this%num_subst = 0
        this%num_features = 0
        this%num_set = 0
        this%mol_type = 'Unknown'
        this%charge_type = 'Unknown'
        this%num_ring = 0

    end subroutine Mol_InitWithNoItem

    subroutine Mol_InitWithNumber(this, atom_num, bond_num)
        !implicit none
        !type(Molecule_t), intent(inout) :: this
        type(Molecule_t) this
        integer :: i, j

        integer, intent(in) :: atom_num
        integer, intent(in) :: bond_num

!        allocate(this%atoms(atom_num))
!        allocate(this%bonds(bond_num))
!        allocate(this%rings(0))
        this%num_atom = atom_num
        this%num_bond = bond_num

        this%name = 'Unknown'
        this%formula = 'Unknown'
        this%weight = 0.0
        this%num_hb_atom = 0
        this%num_rotor = 0
        this%num_atom = atom_num
        this%num_bond = bond_num
        this%num_subst = 0
        this%num_features = 0
        this%num_set = 0
        this%mol_type = 'Unknown'
        this%charge_type = 'Unknown'
        this%num_ring = 0

        do i=1, this%num_atom
            call xlogp_Atom_Clear(this%atoms(i))
            this%atoms(i)%id = i
            this%atoms(i)%name = 'Un'
        enddo

        do i=1, this%num_bond
            call xlogp_Bond_Clear(this%bonds(i))
            this%bonds(i)%id = i
        enddo

    end subroutine Mol_InitWithNumber

    subroutine Mol_deallocate(this)
        !implicit none
        !type(Molecule_t), intent(inout) :: mol
        type(Molecule_t) this

!        if (allocated(this%atoms)) deallocate(this%atoms)
!        if (allocated(this%bonds)) deallocate(this%bonds)
!        if (allocated(this%rings)) deallocate(this%rings)
        this%num_atom = 0
        this%num_bond = 0
        this%num_ring = 0

    end subroutine Mol_deallocate

    !------------------------------------------------------------
    ! Clear the Molecule object
    subroutine Mol_Clear(this)
        !type(Molecule_t), intent(inout) :: this
        type(Molecule_t) this

        this%name = 'Unknown'
        this%formula = 'Unknown'
        this%weight = 0.0
        this%num_hb_atom = 0
        this%num_rotor = 0
        this%num_atom = 0
        this%num_bond = 0
        this%num_subst = 0
        this%num_features = 0
        this%num_set = 0
        this%mol_type = 'Unknown'
        this%charge_type = 'Unknown'
        this%num_ring = 0

    end subroutine Mol_Clear

    subroutine Mol_ShowContents(this)
        type(Molecule_t) this

        write (*, '("name = ", A)' ) this%name
        write (*, '("formula = ", A10)' ) this%formula
        write (*, '("weight = ", F8.2)' ) this%weight
        write (*, '("num_hb_atom = ", I5)' ) this%num_hb_atom
        write (*, '("num_rotor = ", I5)' ) this%num_rotor
        write (*, '("num_atom = ", I5)' ) this%num_atom
        write (*, '("num_bond = ", I5)' ) this%num_bond
        write (*, '("num_subst = ", I5)' ) this%num_subst
        write (*, '("num_features = ", I5)' ) this%num_features
        write (*, '("num_set = ", I5)' ) this%num_set
        write (*, '("mol_type = ", A)' ) this%mol_type
        write (*, '("charge_type = ", A)' ) this%charge_type
        write (*, '("num_ring = ", I5)' ) this%num_ring


    end subroutine Mol_ShowContents

    subroutine Mol_ShowAtoms(xd_m1)
        type(Molecule_t) xd_m1

        integer :: i

        call Atom_ShowHeader(xd_m1%atoms(1))

        do i=1, xd_m1%num_atom
            call Atom_Show(xd_m1%atoms(i))
        end do

    end subroutine Mol_ShowAtoms

    subroutine M_ShowAtomNeighbors(xd_m1)
        type(Molecule_t) xd_m1

        integer :: i

        do i=1, xd_m1%num_atom
            call xlogp_Atom_ShowNeighbor(xd_m1%atoms(i))
        end do

    end subroutine M_ShowAtomNeighbors

    subroutine Mol_ShowBonds(xd_m1)
        type(Molecule_t) xd_m1

        integer :: i

        call xlogp_Bond_ShowHeader(xd_m1%bonds(1))

        do i=1, xd_m1%num_bond
            call xlogp_Bond_ShowContents(xd_m1%bonds(i))
        end do

    end subroutine Mol_ShowBonds

    subroutine Mol_SetAtom(this, atomIndex, name, type, hb, ring)
        type(Molecule_t) :: this
        integer, intent(in) :: atomIndex
        character(len=*), intent(in) :: name
        character(len=*), intent(in) :: type
        character(len=*), intent(in) :: hb
        integer, intent(in), optional :: ring

        if(atomIndex < 1 .OR. atomIndex > this%num_atom) then
            print*, '(*) Boundary out of error: atomIndex = ', atomIndex, &
                    ' num_atom = ', this%num_atom
            return
        endif

        this%atoms(atomIndex)%id = atomIndex
        this%atoms(atomIndex)%name = name
        this%atoms(atomIndex)%type = type
        this%atoms(atomIndex)%hb = hb

        if(present(ring)) then
            this%atoms(atomIndex)%ring = ring
        else
            this%atoms(atomIndex)%ring = 0
        endif

        this%atoms(atomIndex)%valid = .true.
    end subroutine Mol_SetAtom

    subroutine Mol_SetBond(this, bondIndex, atomId1, atomId2, bondOrder)
        type(Molecule_t) :: this
        integer, intent(in) :: bondIndex
        integer, intent(in) :: atomId1
        integer, intent(in) :: atomId2
        character(len=*), intent(in) :: bondOrder

        if(bondIndex < 1 .OR. bondIndex > this%num_bond) then
            print*, '(*) Boundary out of error: bondIndex = ', bondIndex, &
                    ' num_bond = ', this%num_bond
            return
        endif

        if(atomId1 < 1 .OR. atomId1 > this%num_atom) then
            print*, '(*) Boundary out of error: atomId1 = ', atomId1, &
                    ' num_atom = ', this%num_atom
            return
        endif

        if(atomId2 < 1 .OR. atomId2 > this%num_atom) then
            print*, '(*) Boundary out of error: atomId2 = ', atomId2, &
                    ' num_atom = ', this%num_atom
            return
        endif

        this%bonds(bondIndex)%atom_1 = atomId1
        this%bonds(bondIndex)%atom_2 = atomId2
        this%bonds(bondIndex)%type = bondOrder ! '1', '2', '3', 'ar'
        !! @TODO : This bond's valid property is not clear, and it
        !!          need to be cleaned up later.
        !!
        if(bondOrder == '1') then
            this%bonds(bondIndex)%valid = 1
        else if(bondOrder == '2') then
            this%bonds(bondIndex)%valid = 2
        else
            this%bonds(bondIndex)%valid = 3
        endif

    end subroutine Mol_SetBond

    subroutine M_mergeAtomsToFragments(this)
        type(Molecule_t) this
        integer i, j, k, mark, tmp1, tmp2, id1, id2
        integer old_part, new_part

        ! 2) Merge the atoms into fragments
        do i=1, this%num_bond
            if (xlogp_Bond_GetValid(this%bonds(i)) .LE. 0) cycle

            tmp1 = this%bonds(i)%atom_1
            tmp2 = this%bonds(i)%atom_2

            if (this%atoms(tmp1)%part >= this%atoms(tmp2)%part) then
                old_part = this%atoms(tmp1)%part
                new_part = this%atoms(tmp2)%part
            else
                old_part = this%atoms(tmp2)%part
                new_part = this%atoms(tmp1)%part
            end if

            do j=1, this%num_atom
                if (this%atoms(j)%part /= old_part) then
                    cycle
                else
                    this%atoms(j)%part = new_part
                end if
            end do
        end do ! End of merge

    end subroutine M_mergeAtomsToFragments

    subroutine M_make_mol_from_ligand(this, ligand)
        type(Molecule_t) this

        type(ref_lig_type), intent(in) :: ligand
        integer :: AlloAtomStatus, AlloBondStatus
        integer :: i

        ! Assign num_atom and num_bond
        this%num_atom = ligand%n_atm
        this%num_bond = ligand%n_bnd
        if(this%num_atom > XSCORE_MAX_LIGAND_ATOM_NUMBER) then
            print*, '(*) ERROR given atom size (', this%num_atom, &
                    ') is greater than possible max atom size', &
                    XSCORE_MAX_LIGAND_ATOM_NUMBER
            stop
        endif

        if(this%num_bond > XSCORE_MAX_LIGAND_BOND_NUMBER) then
            print*, '(*) ERROR given bond size (', this%num_bond, &
                    ') is greater than possible max bond size', &
                    XSCORE_MAX_LIGAND_BOND_NUMBER
            stop
        endif

        do i=1, this%num_atom
            call xlogp_Atom_Clear(this%atoms(i))
            this%atoms(i)%name = ligand%atom_name(i)
            this%atoms(i)%type = ligand%mol2_type(i)
            this%atoms(i)%id = i
            this%atoms(i)%coor(1) = ligand%R(1,i)
            this%atoms(i)%coor(2) = ligand%R(2,i)
            this%atoms(i)%coor(3) = ligand%R(3,i)
            this%atoms(i)%q = ligand%charge(i)
            this%atoms(i)%vdwTypeNumber_ = 0
            this%atoms(i)%transformedVdwTypeNumber_ = 0
            this%atoms(i)%valid = .true.
        end do

        do i=1, this%num_bond
            call xlogp_Bond_Clear(this%bonds(i))
            this%bonds(i)%id = ligand%bnd(1,i)
            this%bonds(i)%atom_1 = ligand%bnd(2,i)
            this%bonds(i)%atom_2 = ligand%bnd(3,i)
            this%bonds(i)%type = ligand%bnd_type(i)
            this%bonds(i)%valid = 1
            call xlogp_Bond_CheckValidity(this%bonds(i))
        end do

    end subroutine M_make_mol_from_ligand

    subroutine M_make_mol_from_SimpleMol(this, simpMol)
        type(Molecule_t) this

        type(Mol2SimpleMolecule_t), intent(in) :: simpMol
        integer :: AlloAtomStatus, AlloBondStatus
        integer :: i

!        if (allocated(this%atoms)) deallocate(this%atoms)
!        if (allocated(this%bonds)) deallocate(this%bonds)
!        if (allocated(this%rings)) deallocate(this%rings)

        ! molecule name
        this%name = simpMol%name

        ! Assign num_atom and num_bond
        this%num_atom = simpMol%num_atom
        this%num_bond = simpMol%num_bond
!        allocate(this%atoms(this%num_atom), stat=AlloAtomStatus)
!        allocate(this%bonds(this%num_bond), stat=AlloBondStatus)
!        if (AlloAtomStatus /= 0 .or. AlloBondStatus /= 0) then
!            print *, 'Allocation Error in xd_molecule'
!        end if
        if(this%num_atom > XSCORE_MAX_LIGAND_ATOM_NUMBER) then
            print*, '(*) ERROR given atom size (', this%num_atom, &
                    ') is greater than possible max atom size', &
                    XSCORE_MAX_LIGAND_ATOM_NUMBER
            stop
        endif

        if(this%num_bond > XSCORE_MAX_LIGAND_BOND_NUMBER) then
            print*, '(*) ERROR given bond size (', this%num_bond, &
                    ') is greater than possible max bond size', &
                    XSCORE_MAX_LIGAND_BOND_NUMBER
            stop
        endif

        do i=1, this%num_atom
            call xlogp_Atom_Clear(this%atoms(i))
            this%atoms(i)%name = simpMol%atoms(i)%name
            this%atoms(i)%type = simpMol%atoms(i)%type
            this%atoms(i)%id = simpMol%atoms(i)%id
            this%atoms(i)%coor(1) = simpMol%atoms(i)%coor(1)
            this%atoms(i)%coor(2) = simpMol%atoms(i)%coor(2)
            this%atoms(i)%coor(3) = simpMol%atoms(i)%coor(3)
            this%atoms(i)%q = simpMol%atoms(i)%q
            this%atoms(i)%vdwTypeNumber_ = 0
            this%atoms(i)%transformedVdwTypeNumber_ = 0
            this%atoms(i)%valid = .true.
        end do

        do i=1, this%num_bond
            call xlogp_Bond_Clear(this%bonds(i))
            this%bonds(i)%id = simpMol%bonds(i)%id
            this%bonds(i)%atom_1 = simpMol%bonds(i)%atom_1
            this%bonds(i)%atom_2 = simpMol%bonds(i)%atom_2
            this%bonds(i)%type = simpMol%bonds(i)%type
            this%bonds(i)%valid = 1
            call xlogp_Bond_CheckValidity(this%bonds(i))
        end do

    end subroutine M_make_mol_from_SimpleMol


    subroutine M_make_mol_from_protein(this, protein)
        type(Molecule_t) this

        type(molecule_type), intent(in) :: protein

    end subroutine M_make_mol_from_protein

    subroutine M_initAtomParts(this)
        type(Molecule_t) this
        integer i

        do i=1, this%num_atom
            this%atoms(i)%part = this%atoms(i)%id
        end do
    end subroutine M_initAtomParts

    !!----------------------------------------------------------------------------
    !> This function requires that Molecule_t has been already incorporated
    !   with its atoms and bonds.
    !
    !!----------------------------------------------------------------------------
    subroutine M_RearrangeFragIdsCont(this)
        type(Molecule_t) this
        integer :: partList(this%num_atom)
        integer :: i, j, numSubst
        logical :: lMark

    !TODO : This function require unittest, but since I could not understand
    ! this procedure, I omit the unittest intentionally.
    !
        numSubst = 1
        OUTER: do i=1, this%num_atom
            lMark = .false.
            INNER: do j=1, numSubst
                if (this%atoms(i)%part /= partList(j)) then
                    cycle
                else
                    lMark = .true.
                    EXIT INNER
                end if
            end do INNER

            if(lMark) then
                this%atoms(i)%part = j + 1 !! An existing fragment
            else !! A new fragment
                partList(numSubst) = this%atoms(i)%part
                numSubst = numSubst + 1
                this%atoms(i)%part = numSubst
            end if
        end do OUTER
    end subroutine M_RearrangeFragIdsCont

    subroutine M_DefinePartsOfBonds(this)
        type(Molecule_t) this
        integer i, atom_1

        do i=1, this%num_bond
            if(this%bonds(i)%valid .GE. 1) then
                this%bonds(i)%part = 1
            else
                atom_1 = this%bonds(i)%atom_1
                this%bonds(i)%part = this%atoms(atom_1)%part
           end if
        end do
    end subroutine M_DefinePartsOfBonds

    !!----------------------------------------------------------------------------
    !> Assign weights to atoms, which is necessary to rank
    !   the neighboring atoms for a central atom.
    !
    !!----------------------------------------------------------------------------

    subroutine M_SetAtomiWeights(this)
        type(Molecule_t), target :: this
        type(Atom_t), DIMENSION(:), pointer :: at
        integer i
        character(len=2) atomType

        at => this%atoms

        do i=1, this%num_atom
            atomType = TRIM(at(i)%type)

            if (atomType == 'F') then
                at(i)%weight = 19.00
            else if(atomType(1:2) == 'Cl') then
                at(i)%weight = 35.45
            else if(atomType(1:2) == 'CL') then
                at(i)%weight = 35.45
            else if(atomType(1:2) == 'Br') then
                at(i)%weight = 79.90
            else if(atomType(1:2) == 'BR') then
                at(i)%weight = 79.90
            else if(atomType == 'I') then
                at(i)%weight = 126.90
            else if(atomType(1:2) == 'Si') then
                at(i)%weight = 28.09
            else if(atomType(1:2) == 'SI') then
                at(i)%weight = 28.09
            else if(atomType(1:1) == 'C') then
                at(i)%weight = 12.01
            else if(atomType(1:1) == 'H') then
                at(i)%weight = 1.00
            else if(atomType(1:1) == 'N') then
                at(i)%weight = 14.01
            else if(atomType(1:1) == 'O') then
                at(i)%weight = 16.00
            else if(atomType(1:1) == 'P') then
                at(i)%weight = 30.97
            else if(atomType(1:1) == 'S') then
                at(i)%weight = 32.06
            else
                at(i)%weight = 0.00
            end if

        end do
    end subroutine M_SetAtomiWeights

    subroutine M_DetectNeighborAtoms(this)
        type(Molecule_t), intent(inout), target :: this
        type(Atom_t), DIMENSION(:), pointer :: at
        type(Bond_t), DIMENSION(:), pointer :: bn
        integer i, id1, id2, tmp1, tmp2

        at => this%atoms
        bn => this%bonds

        ! initialize num_neib and num_nonh values per each atom.
        do i=1, this%num_atom
            at(i)%num_neib = 0
            at(i)%num_nonh = 0
        end do

        do i=1, this%num_bond
            if(bn(i)%valid .LE. 0) cycle

            id1 = bn(i)%atom_1
            id2 = bn(i)%atom_2

            tmp1 = at(id1)%num_neib

            if(tmp1 < XSCORE_MAX_ATOM_NEIB) then
                at(id1)%neib(tmp1+1) = id2
                at(id1)%bond(tmp1+1) = i
                at(id1)%num_neib = at(id1)%num_neib + 1
            end if

            tmp2 = at(id2)%num_neib

            if(tmp2 < XSCORE_MAX_ATOM_NEIB) then
                at(id2)%neib(tmp2+1) = id1
                at(id2)%bond(tmp2+1) = i
                at(id2)%num_neib = at(id2)%num_neib + 1
            end if
        end do

    end subroutine M_DetectNeighborAtoms

    !!----------------------------------------------------------------------------
    !> Arrange the neighboring atoms in a decreasing order
    !! according to atomic weights; Bonds are also re-arranged
    !!----------------------------------------------------------------------------

    subroutine M_ArrangeNeighborAtoms(this)
        type(Molecule_t), intent(inout), target :: this
        type(Atom_t), DIMENSION(:), pointer :: at
        type(Bond_t), DIMENSION(:), pointer :: bn
        integer i, j, k, id1, id2, tmp1, tmp2, mark

        at => this%atoms
        bn => this%bonds

        do i=1, this%num_atom

            if (.not. at(i)%valid ) then
                cycle
            else if(at(i)%num_neib <= 1) then
                cycle
            end if


            do j=1, at(i)%num_neib
                do k=j+1, at(i)%num_neib
                    tmp1 = at(i)%neib(j)
                    tmp2 = at(i)%neib(k)

                    if(at(tmp1)%weight >= at(tmp2)%weight) then
                        cycle
                    else
                        mark = at(i)%neib(j)
                        at(i)%neib(j) = at(i)%neib(k)
                        at(i)%neib(k) = mark

                        mark = at(i)%bond(j)
                        at(i)%bond(j) = at(i)%bond(k)
                        at(i)%bond(k) = mark
                    end if
                end do
            end do
        end do

    end subroutine M_ArrangeNeighborAtoms

    subroutine M_CountNonHydAtomNum(this)
        type(Molecule_t), intent(inout), target :: this
        type(Atom_t), DIMENSION(:), pointer :: at
        integer i, j, k, id1, id2, tmp1, tmp2, mark
        at => this%atoms

        do i=1, this%num_atom
            at(i)%num_nonh = 0 ! First init

            do j=1,at(i)%num_neib
                if(at(at(i)%neib(j))%type(1:1) == 'H') then
                    cycle
                else
                    at(i)%num_nonh = at(i)%num_nonh + 1
                end if
            end do
        end do

    end subroutine M_CountNonHydAtomNum

    subroutine M_GetJointAtom(bond1, bond2, id)
        type(Bond_t), intent(in) :: bond1
        type(Bond_t), intent(in) :: bond2
        integer, intent(out) :: id

        id = bond1%atom_1
        if (id == bond2%atom_1) then
            return
        else if (id == bond2%atom_2) then
            return
        end if

        id = bond1%atom_2

        if(id == bond2%atom_1) then
            return
        else if(id == bond2%atom_2) then
            return
        end if

        id = 0
        return !! Two bonds are not connected

    end subroutine M_GetJointAtom

    subroutine M_DetectNeighborBonds(this)
        type(Molecule_t), intent(inout), target :: this
        type(Atom_t), DIMENSION(:), pointer :: at
        type(Bond_t), DIMENSION(:), pointer :: bn
        integer i, j, k, id1, id2, tmp1, tmp2, mark
        integer jointAtomId

        at => this%atoms
        bn => this%bonds

        do i=1, this%num_bond
            bn(i)%num_neib = 0
        end do

        do i=1, this%num_bond - 1
            do j=i+1, this%num_bond

                if( (bn(i)%valid .LE. 0) .OR. (bn(j)%valid .LE. 0) ) cycle

                call M_GetJointAtom(bn(i), bn(j), jointAtomId)

                if(jointAtomId > 0) then
                    tmp1 = bn(i)%num_neib
                    if (tmp1 < XSCORE_MAX_BOND_NEIB) then
                        bn(i)%neib(tmp1+1) = bn(j)%id
                        bn(i)%num_neib = bn(i)%num_neib + 1
                    end if

                    tmp2 = bn(j)%num_neib
                    if (tmp2 < XSCORE_MAX_BOND_NEIB) then
                        bn(j)%neib(tmp2+1) = bn(i)%id
                        bn(j)%num_neib = bn(j)%num_neib + 1
                    end if
                else
                    cycle
                end if
            end do
        end do
    end subroutine M_DetectNeighborBonds

    subroutine M_DetectConnections(this)
        type(Molecule_t) this
        integer i, j, k, mark, tmp1, tmp2, old_part, new_part, id1, id2

        !!----------------------------------------------------------------------------
        ! Test 1: Detect all the fragments and determines num_subst

        ! 1) Detect how many fragments are in the molecule and initiaize
        !   the molecule as the assembly of isolated atoms.
        !
        call M_initAtomParts(this)

        ! Now all bonds ar valid
        !
        call M_mergeAtomsToFragments(this)
        call M_RearrangeFragIdsCont(this)
        call M_DefinePartsOfBonds(this)

        !!----------------------------------------------------------------------------
        !> 2) Set up connection tables for atoms and bonds
        !!----------------------------------------------------------------------------

        !> Assign atomic weight to atoms, which is necessary to rank
        !! the neighboring atoms for a central atom.
        call M_SetAtomiWeights(this)
        call M_DetectNeighborAtoms(this)
        call M_ArrangeNeighborAtoms(this)
        call M_CountNonHydAtomNum(this)
        call M_DetectNeighborBonds(this)

    end subroutine M_DetectConnections

    subroutine M_FindAtomGrpEnv(this, group, atomId)
        type(Molecule_t), intent(in), target :: this
        integer, intent(in) :: atomId
        type(Group_t), intent(out) :: group
        type(Atom_t), DIMENSION(:), pointer :: at
        type(Bond_t), DIMENSION(:), pointer :: bn

        at => this%atoms
        bn => this%bonds

        call xlogp_Group_FindAGroup(group, at, bn, atomId)

    end subroutine M_FindAtomGrpEnv

    subroutine M_ShowRings(this)
        type(Molecule_t), intent(in) :: this
        integer i

        do i=1, this%num_ring
            call xlogp_Ring_Show(this%rings(i))
        enddo

    end subroutine M_ShowRings


    function M_GetRingSize(this) result(return_value)
        type(Molecule_t), intent(in) :: this
        integer :: return_value

!        if(allocated(this%rings)) then
!            return_value = size(this%rings)
!            !print*, 'return_value1 = ', return_value
!        else
!            return_value = 0
!            !print*, 'return_value2 = ', return_value
!        endif
        return_value = this%num_ring

    end function M_GetRingSize


    function M_GetWeight(this) result(return_value)
        type(Molecule_t), intent(inout) :: this
        real :: return_value
        integer :: i

        return_value = 0.0
        do i=1, this%num_atom
            if(this%atoms(i)%valid .EQV. .false.) then
                cycle
            else
                return_value = return_value + (this%atoms(i)%weight)
            endif
        enddo
    end function M_GetWeight

    function M_GetNumHBAtom(this) result(return_value)
        type(Molecule_t), intent(inout) :: this
        integer :: return_value
        integer :: i

        return_value = 0
        do i=1, this%num_atom
            if(.not. this%atoms(i)%valid) then
                cycle
            else if(this%atoms(i)%hb == 'D') then
                return_value = return_value + 1
            else if(this%atoms(i)%hb == 'DA') then
                return_value = return_value + 1
            else if(this%atoms(i)%hb == 'A') then
                return_value = return_value + 1
            else
                cycle
            endif
        enddo
    end function M_GetNumHBAtom

    function M_GetNumHeavyAtom(this) result(return_value)
        type(Molecule_t), intent(inout) :: this
        integer :: return_value
        integer :: i

        return_value = 0
        do i=1, this%num_atom
            if(.not. this%atoms(i)%valid) then
                cycle
            else if(this%atoms(i)%type == 'H') then
                cycle
            else
                return_value = return_value + 1
            endif
        enddo
    end function M_GetNumHeavyAtom

    subroutine M_Assign_Apparent_Charge(this)
        type(Molecule_t), intent(inout) :: this
        integer i
        this%charge_type = 'FORMAL_CHARGES'
        do i=1, this%num_atom
            if(this%atoms(i)%xtype == 'O.co2') then
                this%atoms(i)%q = -0.500
            else if(this%atoms(i)%xtype == 'N.4') then
                this%atoms(i)%q = 1.0
            else
                this%atoms(i)%q = 0.0
            endif
        enddo

        do i=1, this%num_atom
            if(this%atoms(i)%xtype /= 'N.pl3.h') then
                cycle
            else if(this%atoms(i)%num_nonh /= 1) then
                cycle
            else if(this%atoms(this%atoms(i)%neib(1))%type /= 'C.cat') then
                cycle
            else
                this%atoms(i)%q = 0.5
            endif
        enddo

    end subroutine M_Assign_Apparent_Charge

    !!----------------------------------------------------------------------------
    !> Calculate root position for each atom's hydrogen bonding.
    !!----------------------------------------------------------------------------
    subroutine M_CalculateHBondRoot(this)
        type(Molecule_t), intent(inout) :: this
        integer i, j, tmp, num_nonh
        real(dp) :: tmpx, tmpy, tmpz
        do i=1, this%num_atom
            if(.not. this%atoms(i)%valid) then
                cycle
            else if(this%atoms(i)%hb == 'N') then
                cycle
            else if(this%atoms(i)%hb == 'H') then
                cycle
            else if(this%atoms(i)%hb == 'P') then
                cycle
            endif

            tmpx = 0.0; tmpy = 0.0; tmpz = 0.0
            num_nonh = 0
            do j=1, this%atoms(i)%num_neib
                tmp = this%atoms(i)%neib(j)
                if(this%atoms(tmp)%type(1:1) == 'H') then
                    cycle
                else
                    tmpx = tmpx + this%atoms(tmp)%coor(1)
                    tmpy = tmpy + this%atoms(tmp)%coor(2)
                    tmpz = tmpz + this%atoms(tmp)%coor(3)
                    num_nonh = num_nonh + 1
                endif
            enddo

            if (num_nonh == 0) then
                this%atoms(i)%hb = 'P'
            else
                tmpx = tmpx / num_nonh
                tmpy = tmpy / num_nonh
                tmpz = tmpz / num_nonh
                this%atoms(i)%root(1) = tmpx
                this%atoms(i)%root(2) = tmpy
                this%atoms(i)%root(3) = tmpz
            endif
        enddo

    end subroutine M_CalculateHBondRoot

    subroutine M_GetFormula(this, formula)
        type(Molecule_t), intent(inout) :: this
        character(len=*), intent(out) :: formula
        integer :: i
        integer :: num(12) !! C, H, N, O, P, S, F, Cl, Br, I, Si, Un
        character(3) :: element(12)
        character(10) :: tmp
        character(256) :: tmp_formula, str1, fmt1
        integer :: strLen

        num = 0

        element(1) = "C"
        element(2) = "H"
        element(3) = "N"
        element(4) = "O"
        element(5) = "P"
        element(6) = "S"
        element(7) = "F"
        element(8) = "Cl"
        element(9) = "Br"
        element(10) = "I"
        element(11) = "Si"
        element(12) = "Un"

        do i=1, this%num_atom
            if(.not. this%atoms(i)%valid) then
                num(12) = num(12) + 1
            else if(this%atoms(i)%type == 'C.3') then
                num(1) = num(1) + 1
            else if(this%atoms(i)%type == 'C.2') then
                num(1) = num(1) + 1
            else if(this%atoms(i)%type == 'C.1') then
                num(1) = num(1) + 1
            else if(this%atoms(i)%type == 'C.cat') then
                num(1) = num(1) + 1
            else if(this%atoms(i)%type == 'C.ar') then
                num(1) = num(1) + 1
            else if(this%atoms(i)%type == 'H') then
                num(2) = num(2) + 1
            else if(this%atoms(i)%type == 'H.spc') then
                num(2) = num(2) + 1
            else if(this%atoms(i)%type == 'N.4') then
                num(3) = num(3) + 1
            else if(this%atoms(i)%type == 'N.3') then
                num(3) = num(3) + 1
            else if(this%atoms(i)%type == 'N.2') then
                num(3) = num(3) + 1
            else if(this%atoms(i)%type == 'N.1') then
                num(3) = num(3) + 1
            else if(this%atoms(i)%type == 'N.ar') then
                num(3) = num(3) + 1
            else if(this%atoms(i)%type == 'N.pl3') then
                num(3) = num(3) + 1
            else if(this%atoms(i)%type == 'N.am') then
                num(3) = num(3) + 1
            else if(this%atoms(i)%type == 'O.3') then
                num(4) = num(4) + 1
            else if(this%atoms(i)%type == 'O.2') then
                num(4) = num(4) + 1
            else if(this%atoms(i)%type == 'O.co2') then
                num(4) = num(4) + 1
            else if(this%atoms(i)%type == 'P.3') then
                num(5) = num(5) + 1
            else if(this%atoms(i)%type == 'S.3') then
                num(6) = num(6) + 1
            else if(this%atoms(i)%type == 'S.2') then
                num(6) = num(6) + 1
            else if(this%atoms(i)%type == 'S.o') then
                num(6) = num(6) + 1
            else if(this%atoms(i)%type == 'S.o2') then
                num(6) = num(6) + 1
            else if(this%atoms(i)%type == 'F') then
                num(7) = num(7) + 1
            else if(this%atoms(i)%type == 'Cl') then
                num(8) = num(8) + 1
            else if(this%atoms(i)%type == 'Br') then
                num(9) = num(9) + 1
            else if(this%atoms(i)%type == 'I') then
                num(10) = num(10) + 1
            else if(this%atoms(i)%type == 'Si') then
                num(11) = num(11) + 1
            else
                num(12) = num(12) + 1
            endif
        enddo

        tmp_formula = ""
        do i=1, 12
            if(num(i) == 0) then
                cycle
            else
                tmp_formula = TRIM(tmp_formula)
                strLen = len_trim(tmp_formula)
                write(tmp_formula, '(A, A)') tmp_formula(1:strLen), element(i)
                if(num(i) == 1) then
                cycle
                else
                    !!call writenum(this%mol%atoms(i)%id, str1, "I3")
                    call writenum(num(i), tmp, "I3")
                    strLen = len_trim(tmp_formula)
                    write(tmp_formula, '(A, A)') tmp_formula(1:strLen), TRIM(tmp)
                endif
            endif
        enddo

        formula = trim(tmp_formula)
    end subroutine M_GetFormula

end module xlogp_Molecule_m
