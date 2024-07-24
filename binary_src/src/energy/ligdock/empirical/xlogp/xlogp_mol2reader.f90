module xlogp_Mol2Reader_m
    use xlogp_Atom_m
    use xlogp_Bond_m
    use xlogp_Ring_m
    use xlogp_constants_m

    private ! Hide all

    private :: ParseString ! This method was temporarily copied into this module for orthogonality.
    public :: ReadFromMol2
    public :: Mol2Reader_Clear
    public :: Mol2SimpleMolecule_t

    !------------------------------------------------------------
    ! Mol2Reader Derived Type
    type Mol2SimpleMolecule_t

        ! sequence ! Let this Mol2Reader type stored as sequential
        character(len=256) :: name
        character(len=256) :: formula
        real :: weight
        integer :: num_hb_atom
        integer :: num_rotor
        integer :: num_atom
        type(Atom_t), DIMENSION(XSCORE_MAX_LIGAND_ATOM_NUMBER) :: atoms ! Allocatable atoms container
        integer :: num_bond
        type(Bond_t), DIMENSION(XSCORE_MAX_LIGAND_BOND_NUMBER) :: bonds ! Allocatable bonds container

        integer :: num_subst
        integer :: num_features
        integer :: num_set

        character(len=256) :: mol_type
        character(len=256) :: charge_type

        integer :: num_ring
        type(Ring_t), DIMENSION(XSCORE_MAX_RING_NUMBER) :: rings ! Allocatable rings container

    end type Mol2SimpleMolecule_t


contains
    !------------------------------------------------------------
    ! Clear the Mol2Reader object
    subroutine Mol2Reader_Clear(this)
        type(Mol2SimpleMolecule_t), intent(inout) :: this

        this%name = 'Unknown'
        this%formula = 'No formula'
        this%weight = 0.0
        this%num_hb_atom = 0
        this%num_rotor = 0
        this%num_atom = 0
        this%num_subst = 0
        this%num_features = 0
        this%num_set = 0
        this%mol_type = 'Unknown molecule type'
        this%charge_type = 'Unknown charge type'
        this%num_ring = 0

    end subroutine Mol2Reader_Clear

    subroutine Mol2Reader_Deallocate(this)
        type(Mol2SimpleMolecule_t), intent(inout) :: this

!        if(allocated(this%atoms)) then
!            deallocate(this%atoms)
!        end if
!        if(allocated(this%bonds)) then
!            deallocate(this%bonds)
!        end if
    end subroutine Mol2Reader_Deallocate

    subroutine Mol2Reader_ShowAtoms(this)
        implicit none
        type(Mol2SimpleMolecule_t), intent(in) :: this

        integer :: i

        call Atom_ShowHeader(this%atoms(1))

        do i=1, this%num_atom
            call Atom_Show(this%atoms(i))
        end do

    end subroutine Mol2Reader_ShowAtoms

    subroutine Mol2Reader_ShowBonds(this)
        implicit none
        type(Mol2SimpleMolecule_t), intent(in) :: this

        integer :: i

        do i=1, this%num_bond
            call xlogp_Bond_ShowContents(this%bonds(i))
        end do

    end subroutine Mol2Reader_ShowBonds

    subroutine ReadFromMol2(this, mol2_file, readResult)
        implicit none
        type(Mol2SimpleMolecule_t), intent(inout) :: this
        character(len=120), intent(in) :: mol2_file
        integer :: i_unit, ioerror, n_word, i_line, atm_idx, bnd_idx
        character(len=120) :: word(25), keyword, line
        logical :: molecule_info, atm_info, bnd_info
        logical :: op
        INTEGER, PARAMETER :: MY_UNIT_NUMBER = 31 !> constant unit number
        logical, intent(out) :: readResult
        integer :: AlloAtomStatus, AlloBondStatus
        ! END OF DECLARE

        ! Open file
        i_unit = MY_UNIT_NUMBER
        readResult = .true.

        open(file=trim(mol2_file), unit=i_unit, action='read', status='old', iostat=ioerror)
        if(ioerror /= 0) then
            print*, '(*) ERROR: MOL2 file does not exist! (xd_mol2reader.f90) '
            print*, '(*)        MOL2 file name = ', mol2_file
            readResult = .false.
            return
        end if

        ! Initialize variables
        molecule_info = .false.
        atm_info = .false.
        bnd_info = .false.
        i_line = 0

        ! Read file
        do
            read(i_unit, '(a80)', iostat=ioerror) line
            if(ioerror /= 0) then ! End of file
                exit
            end if

            if(line(1:1) == '#') then ! Skip unnecessary information
                cycle
            end if

            ! Read molecule information
            if(line(1:17) == '@<TRIPOS>MOLECULE') then
                molecule_info = .true.
            end if

            if(molecule_info) then
                i_line = i_line + 1
            end if
            if(molecule_info .and. i_line == 2) then
                call ParseString(line, n_word, word)
                read(word(1), '(A)') this%name
            endif

            if(molecule_info .and. i_line == 3) then
                call ParseString(line, n_word, word)

                read(word(1), '(i3)') this%num_atom
                read(word(2), '(i3)') this%num_bond
                molecule_info = .false.
                i_line = 0

                !> Now we know the number of atoms and bonds
                !! Let me allocate storage for atoms and bonds
                !!
!                if (allocated(this%atoms)) then
!                    deallocate(this%atoms)
!                end if
!                allocate(this%atoms(this%num_atom), stat=AlloAtomStatus)
!
!                if (allocated(this%bonds)) then
!                    deallocate(this%bonds)
!                end if
!                allocate(this%bonds(this%num_bond), stat=AlloBondStatus)

!                if(AlloAtomStatus /= 0 .or. AlloBondStatus /= 0) then
!                    print*, 'Allocation Error in xd_mol2reader'
!                    readResult = .false.
!                    return
!                end if

            end if

            ! Read atom information
            if(line(1:13) == '@<TRIPOS>ATOM') then
                atm_info = .true.
            end if

            if(atm_info .and. line(1:1) /= '@') then

                call ParseString(line, n_word, word)

                i_line = i_line + 1
                read(word(1), '(i3)') this%atoms(i_line)%id         ! atom number
                read(word(2), '(a4)') this%atoms(i_line)%name       ! atom name
                read(word(3), '(f8.3)') this%atoms(i_line)%coor(1)
                read(word(4), '(f8.3)') this%atoms(i_line)%coor(2)
                read(word(5), '(f8.3)') this%atoms(i_line)%coor(3)
                read(word(6), '(a6)') this%atoms(i_line)%type   ! mol2_type
                read(word(9), '(f8.3)') this%atoms(i_line)%q    ! charge
            end if

            ! Read bond information
            if(line(1:13) == '@<TRIPOS>BOND') then
                bnd_info = .true.
                atm_info = .false.
                i_line = 0
            end if

            if(bnd_info .and. line(1:1) /= '@') then

                call ParseString(line, n_word, word)

                i_line = i_line + 1
                read(word(1), '(i6)') this%bonds(i_line)%id
                read(word(2), '(i5)') this%bonds(i_line)%atom_1
                read(word(3), '(i5)') this%bonds(i_line)%atom_2

                ! Modified by HW Chung.
                ! Sometimes bond info has brief description about bond as follows:
                !     32   33   32 2    BACKBONE|DICT
                ! So n_word == 4 -> n_word >=4
                if(n_word .GE. 4) then
                    this%bonds(i_line)%type = word(4)
                end if
            end if

            if(line(1:21) == '@<TRIPOS>SUBSTRUCTURE') then
                bnd_info = .false.
            end if
        end do

        close(i_unit)

        ! Check file consistency
        if(this%atoms(this%num_atom)%id /= this%num_atom) then
            print*, "Error: Inconsistent atom information. Plz check yr mol2 file!"
            readResult = .false.
        end if
        if(this%bonds(this%num_bond)%id /= this%num_bond) then
            print*, "Error: Inconsistent bond information. Plz check yr mol2 file!"
            readResult = .false.
        end if
    end subroutine ReadFromMol2

    subroutine ParseString(long_string, num_word, word)
        !-------------------------------------------------------------------------------
        ! Scan string, get num of "word", and save in array word
        !-------------------------------------------------------------------------------
        character(len=*), intent(in) :: long_string
        integer, intent(out) :: num_word
        character(len=*), intent(out) :: word(:)
        character(len=120) :: string
        integer :: i, n, k
        integer, parameter :: len_fname = 120
        character(len=1) :: prev_char, cur_char

        n = len(long_string)
!        print*, 'long_string = ', long_string

        if (n > len_fname) then
            write(*,"(A,I6)") 'Error. The length of the long string is greater than maximum: ', len_fname
            return
        end if

        num_word = 0
        prev_char = ' '
        string = ' '
        k = 0

        ! go over each character
        do i = 1, n
            cur_char = long_string(i:i)
            if (len(trim(cur_char)) /= 0) then
                if (len(trim(prev_char)) == 0) then
                    ! start of a new word
                    num_word = num_word + 1
                    if (num_word > 1) then
                        ! save previous word and start new
                        word(num_word-1) = string
                        string = ' '
                        k = 0
                    end if
                end if

                if (len(trim(cur_char)) /= 0) then ! add
                    k = k + 1
                    string(k:k) = cur_char
                end if
            end if
            prev_char = cur_char
        end do

        ! save the last word
        if (num_word > 0) then
            word(num_word) = string
        end if
    end subroutine ParseString

end module xlogp_Mol2Reader_m
