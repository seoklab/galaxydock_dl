module xlogp_Parameter_m
    use globals
    use xlogp_Atom_m
    use xlogp_Strings_m
    use xlogp_AtomDefs_m
    use xlogp_XAtomDefs_m
    use xlogp_ResidueDefs_m
    use xlogp_constants_m

    implicit none

    public :: getUnitNumber
    public :: FileExist
    public :: Check_Blank_Line

    integer, private, parameter :: m_dp = 8    !< size of double precision
    logical, public, parameter :: m_lDebug = .false.

    type Ratom_t
        character(len=10) :: name
        character(len=10) :: type
        character(len=10) :: xtype
        real(dp) :: r
        real(dp) :: eps
        real(dp) :: q
        character(len=3) :: hb
        real(dp) :: logp
        real(dp) :: solv
        integer :: ring
        character(len=3) :: pmftype
    end type Ratom_t

    type Rbond_t
        character(len=10) :: atom1
        character(len=10) :: atom2
        character(len=3) :: type
    end type Rbond_t

    type Residue_Def_t
        character(len=10) :: name   !< residue name
        integer :: num_atom         !< number of atoms
        integer :: num_bond         !< number of bonds
        real(m_dp) :: total_charge    !< Total charge of a residue
        character(len=100) :: desc  !< description about residues
        type(Ratom_t) :: atom(50)     !< Residue atoms
        type(Rbond_t) :: bond(50)     !< Residue bonds
    end type Residue_Def_t

    type Atom_Def_t
        character(len=10) :: type
        integer(4) :: vdwTypeNumber_
        real(m_dp) :: weight
        real(m_dp) :: r
        real(m_dp) :: eps
        real(m_dp) :: q
        character(len=3) :: hb
        character(len=50) :: description
    end type Atom_Def_t

    type Xatom_Def
        character(len=20) :: type
        character(len=3) :: hb
        real(m_dp) :: logp
    end type Xatom_Def

    type Bond_Def
        character(len=10) :: atom_1
        character(len=10) :: atom_2
        character(len=3) :: type
        real(m_dp) :: length
    end type Bond_Def

    type Tors_Def
        character(len=10) :: atom_1
        character(len=10) :: atom_2
        character(len=10) :: atom_3
        character(len=10) :: atom_4
        character(len=10) :: type
        real(m_dp) :: V_potential_barrier !< twisting force constant
        integer :: n                    !< periodicity
        integer :: S                    !< sign of torsion angle type
    end type Tors_Def


    type ForceField_t

        integer :: m_num_restype
        integer :: m_num_atomtype
        integer :: m_num_xatomtype
        integer :: m_num_bondtype
        integer :: m_num_tortype

        type(Residue_Def_t), DIMENSION(XSCORE_MAX_RESIDUE_DEF_NUM) :: m_residue
        type(Atom_Def_t), DIMENSION(XSCORE_MAX_ATOM_DEF_NUM) :: m_atom
        type(Xatom_Def), DIMENSION(XSCORE_MAX_XATOM_DEF_NUM) :: m_xatom
        type(Bond_Def), DIMENSION(XSCORE_MAX_BOND_DEF_NUM) :: m_bond
        type(Tors_Def), DIMENSION(XSCORE_MAX_TORS_DEF_NUM) :: m_torsion

        logical :: init != .false.

    end type ForceField_t

contains

    subroutine Param_ShowAllResidueDef(this)
        type(ForceField_t) :: this
        integer :: i
        if(size(this%m_residue) == 0) then
            print*, 'All residue definition size is zero'
            return
        endif
        print*, 'm_residue size = ', this%m_num_restype

        do i=1, this%m_num_restype
            call Param_ShowResidueDef(this%m_residue(i))
        enddo
    end subroutine Param_ShowAllResidueDef

    subroutine Param_ShowAllAtomDef(this)
        type(ForceField_t) :: this
        integer :: i
        if(size(this%m_atom) == 0) then
            print*, 'All atom definition size is zero'
            return
        endif
        print*, 'm_atom size = ', this%m_num_atomtype

        do i=1, this%m_num_atomtype
            call Param_ShowAtomDef(this%m_atom(i), i)
        enddo
    end subroutine Param_ShowAllAtomDef

    subroutine Param_ShowAllXAtomDef(this)
        type(ForceField_t) :: this
        integer :: i
        if(size(this%m_xatom) == 0) then
            print*, 'All xlogp atom definition size is zero'
            return
        endif
        print*, 'm_xatom size = ', this%m_num_xatomtype

        do i=1, this%m_num_xatomtype
            call Param_ShowXAtomDef(this%m_xatom(i), i)
        enddo
    end subroutine Param_ShowAllXAtomDef

    subroutine Param_ClearAtomDef(this)
        type(Atom_Def_t) :: this
        this%type = ''
        this%vdwTypeNumber_ = -1
        this%weight = 0.0
        this%r = 0.0
        this%eps = 0.0
        this%q = 0.0
        this%hb = ''
    end subroutine Param_ClearAtomDef

    subroutine Param_ShowResidueDef(this)
        type(Residue_Def_t) :: this
        integer :: i
100     FORMAT (A4, 2X, A4, 3X, I2, 5X, I2, 3X, F5.2, 1X, A46)
        write(*, 100) "RESI", this%name, this%num_atom, this%num_bond, this%total_charge, this%desc
        do i=1, this%num_atom
            call Param_ShowResidueAtom(this%atom(i))
        enddo
        do i=1, this%num_bond
            call Param_ShowResidueBond(this%bond(i))
        enddo
    end subroutine Param_ShowResidueDef

    subroutine Param_ShowResidueAtom(this)
        type(Ratom_t) :: this
200     FORMAT (A4, 2X, A4, 3X, A7, A10, F6.3, 1X, F6.3, 1X, F6.3, 2X, &
            A2, 2X, F6.3, 2X, F6.3, 2X, I1, 2X, A2, 2X)
        write(*, 200) "ATOM", this%name, this%type, this%xtype, &
            this%r, this%eps, this%q, this%hb, this%logp, &
            this%solv, this%ring, this%pmftype
    end subroutine Param_ShowResidueAtom

    subroutine Param_ShowResidueBond(this)
        type(Rbond_t) :: this
300     FORMAT (A4, 2X, A4, 3X, A4, 3X, A3, A54)
        write(*, 300) "BOND", this%atom1, this%atom2, this%type, " "
    end subroutine Param_ShowResidueBond

    subroutine Param_ShowAtomDef(this, id)
        type(Atom_Def_t) :: this
        integer, intent(in) :: id
        !        character(len=10) :: type
        !        integer(4) :: vdwTypeNumber_
        !        real(m_dp) :: weight
        !        real(m_dp) :: r
        !        real(m_dp) :: eps
        !        real(m_dp) :: q
        !        character(len=3) :: hb
        write(*, '(I3, 2X, A10, 2X, F8.4, 2x, F8.4, 2x, F8.4, 2x, F8.4, 2x, A3, 2x, A35)') &
            id, this%type, this%weight, this%r, this%eps, this%q, this%hb, this%description
    end subroutine Param_ShowAtomDef

    subroutine Param_ShowXAtomDef(this, id)
        type(Xatom_Def) :: this
        integer, intent(in) :: id
        write(*, '(I3, 2X, A20, 2x, F6.3, 2x, A3)') &
            id, this%type, this%logp, this%hb
    end subroutine Param_ShowXAtomDef

    subroutine Param_SetUp(this, residueDefFile, xtoolAtomDefFile, xlogpAtomDefFile)
        type(ForceField_t) :: this
        character(len=*), intent(in) :: residueDefFile
        character(len=*), intent(in) :: xtoolAtomDefFile
        character(len=*), intent(in) :: xlogpAtomDefFile

        !        print*, 'Setup started ...'
        !        print*, 'residueDefFile = ', residueDefFile
        !        print*, 'xtoolAtomDefFile = ', xtoolAtomDefFile
        !        print*, 'xlogpAtomDefFile = ', xlogpAtomDefFile

        call Param_ReadResidueDef(this, residueDefFile)
        call Param_ReadAtomDefinition(this, xtoolAtomDefFile)
        call Param_ReadXAtomDef(this, xlogpAtomDefFile)

    !        print*, 'Setup Ended ...'
    !        call this%ShowAllAtomDef()

    end subroutine Param_SetUp

    subroutine Param_SetUpInMemory(this)
        type(ForceField_t) :: this

        call Param_ReadResDefInMemory(this)
        call Param_ReadAtomDefInMemory(this)
        call Param_ReadXAtomDefInMem(this)

    !        print*, 'Setup Ended ...'
    !        call this%ShowAllAtomDef()

    end subroutine Param_SetUpInMemory

    function FileExist(fileName) result(ex)
        !< Return variable
        logical :: ex

        !< Dummy argument
        character(len=*), intent(in) :: fileName

        !< end of decl

        inquire(FILE=fileName, exist=ex)
    end function FileExist

    subroutine Check_Blank_Line(line, lblank_line)
        !> dummy argument
        character(len=*), intent(inout) :: line
        logical, intent(out) :: lblank_line
        !> local variables
        integer :: length
        character :: ic
        integer :: i, mark
        !> end of decl

        lblank_line = .false.
        line = TRIM(line)
        length = LEN_TRIM(line)
        if(length < 1) then
            lblank_line = .true.
            return
        end if

        !i = ICHAR(unknown_string(1:1)) !> NOTE You can get one character using line(1:1) syntax only.

        mark = 0
        do i=1, length
            if(line(i:i) == CHAR(9)) then
                cycle
            else if(line(i:i) == CHAR(10)) then
                cycle
            else if(line(i:i) == CHAR(13)) then
                exit
            else
                mark = mark + 1
            end if
        end do

        if (mark == 0) then
            lblank_line = .true.
        else
            lblank_line = .false.
        end if
    end subroutine Check_Blank_Line

    function Is_Blank_Line(line) result(lblank_line)
        !> Return value
        logical :: lblank_line

        !> dummy argument
        character(len=*), intent(in) :: line

        !> local variables
        character(len=LEN(line)) :: new_line_str
        integer :: length
        character :: ic
        integer :: i, mark
        !> end of decl

        new_line_str = TRIM(line)
        lblank_line = .false.
        length = LEN_TRIM(new_line_str)
        if(length < 1) then
            lblank_line = .true.
            return
        end if

        !i = ICHAR(unknown_string(1:1)) !> NOTE You can get one character using line(1:1) syntax only.

        mark = 0
        do i=1, length
            if(new_line_str(i:i) == CHAR(9)) then
                cycle
            else if(new_line_str(i:i) == CHAR(10)) then
                cycle
            else if(new_line_str(i:i) == CHAR(13)) then
                exit
            else
                mark = mark + 1
            end if
        end do

        if (mark == 0) then
            lblank_line = .true.
        else
            lblank_line = .false.
        end if
    end function Is_Blank_Line

    function getUnitNumber(unit_num) result(unit_allowed)
        integer, intent(in) :: unit_num
        integer :: unit_allowed
        logical :: lop

        unit_allowed = unit_num
        do
            inquire(unit=unit_allowed, opened=lop)
            if(.not. lop) then
                exit
            else
                unit_allowed = unit_allowed + 1
            end if

            if (unit_allowed > 99) then
                print*, 'ERROR>>> allowable unit numbr exceeds 99'
                exit
            end if
        end do
    end function getUnitNumber

    subroutine Param_ReadResidueDef(this, filename)
        type(ForceField_t) :: this
        !< dummy argument
        character(len=*), intent(in) :: filename
        !< local variable
        integer :: i, num, num_atom, num_bond, resi_count
        character(len=120) :: buf, head
        character(len=120) :: word(20)
        integer :: num_word
        integer :: unit_no, ios, mem_stat
        logical :: lblankLine
        ! END of variable declare

        !!----------------------------------------------------------------------------
        !> File Open
        !!----------------------------------------------------------------------------
        unit_no = getUnitNumber(21)
        open(unit=unit_no, FILE=trim(filename), status='OLD', action='READ', iostat=ios)
        if(ios /= 0) then
            write(*, '(3A, I4)') 'FILE ', filename, ' OPEN Failed with erorr code ', ios
            return  !< TODO Here we need additional treatment for exceptional case.
        end if

        !> Count the number of total residues
        resi_count = 0
        do
            read(unit_no, '(A80)', iostat=ios) buf
            if(ios /= 0) then
                exit
            else if(buf(1:1) == '#') then
                cycle
            else if(Is_Blank_Line(buf)) then
                cycle
            else
                call parse_string(buf, num_word, word)
                if(word(1) == 'RESI') then
                    resi_count = resi_count+1
                end if
            end if
        end do

        REWIND(unit_no)

        this%m_num_restype = resi_count
!        if( .not. allocated(this%m_residue) ) then   !< Check if it is already allocated
!            allocate(this%m_residue(this%m_num_restype), stat=mem_stat)  !< Allocate residue container
!            if(mem_stat /= 0) then
!                print*, 'MEMORY allocation error: We will return right now'
!                return
!            end if
!        else if(size(this%m_residue) /= this%m_num_restype) then  !< Check size of allocated memory
!            deallocate(this%m_residue)
!            allocate(this%m_residue(this%m_num_restype), stat=mem_stat)  !< Allocate residue container
!            if(mem_stat /= 0) then
!                print*, 'MEMORY allocation error: We will return right now'
!                return
!            end if
!        end if

        !> NOTE : Now we will read in the residue templates.
        resi_count = 1
        do
            read(unit_no, '(A120)', iostat=ios) buf
            if(ios /= 0) then
                exit
            else if(buf(1:1) == '#') then
                cycle
            else if(Is_Blank_Line(buf)) then
                cycle
            else
                call parse_string(buf, num_word, word)
                if(word(1) == 'END') then
                    exit
                else if(word(1) == 'RESI') then
                    this%m_residue(resi_count)%name = word(2)
                    read(word(3), '(I2)') this%m_residue(resi_count)%num_atom
                    read(word(4), '(I2)') this%m_residue(resi_count)%num_bond
                    read(word(5), '(F6.3)') this%m_residue(resi_count)%total_charge
                    this%m_residue(resi_count)%desc = TRIM(buf(32:))

                    num_atom = this%m_residue(resi_count)%num_atom

                    !< get information for atoms
                    do i=1, num_atom
                        read(unit_no, '(A120)', iostat=ios) buf
                        if(ios /= 0) exit

                        call parse_string(buf, num_word, word)
                        this%m_residue(resi_count)%atom(i)%name = trim(word(2))
                        this%m_residue(resi_count)%atom(i)%type = trim(word(3))
                        this%m_residue(resi_count)%atom(i)%xtype = trim(word(4))
                        read(word(5), '(F6.3)') this%m_residue(resi_count)%atom(i)%r
                        read(word(6), '(F6.3)') this%m_residue(resi_count)%atom(i)%eps
                        read(word(7), '(F6.3)') this%m_residue(resi_count)%atom(i)%q
                        this%m_residue(resi_count)%atom(i)%hb = trim(word(8))
                        read(word(9), '(F6.3)') this%m_residue(resi_count)%atom(i)%logp
                        read(word(10), '(F6.3)') this%m_residue(resi_count)%atom(i)%solv
                        read(word(11), '(I1)') this%m_residue(resi_count)%atom(i)%ring
                        this%m_residue(resi_count)%atom(i)%pmftype = trim(word(12))
                    end do

                    !< Then, read in the bond information
                    num_bond = this%m_residue(resi_count)%num_bond
                    if(num_bond > 0) then
                        do i=1, num_bond
                            read(unit_no, '(A120)', iostat=ios) buf
                            if(ios /= 0) exit

                            call parse_string(buf, num_word, word)
                            this%m_residue(resi_count)%bond(i)%atom1 = trim(word(2))
                            this%m_residue(resi_count)%bond(i)%atom2 = trim(word(3))
                            this%m_residue(resi_count)%bond(i)%type = trim(word(4)) !< bond type: 1, 2, 3, ar, am
                        end do
                    end if

                    resi_count = resi_count + 1
                end if
            end if
        end do !< END of residue reading

        close(unit_no)
    end subroutine Param_ReadResidueDef

    subroutine Param_ReadResDefInMemory(this)
        type(ForceField_t) :: this
        integer :: i, j, num, num_atom, num_bond, resi_count
        character(len=120) :: buf, head
        character(len=120) :: word(20)
        character(len=30) :: delimiters
        character(len=1) :: tabchar
        integer :: num_word
        logical :: lblankLine
        integer :: mem_stat
        ! END of variable declare

        !> Count the number of total residues
        resi_count = 0
        do i=1, NUM_RESIDUE_DEF_STRS
            buf = residue_def_strs(i)
            if(buf(1:1) == '#') then
                cycle
            else if(Is_Blank_Line(buf)) then
                cycle
            else
                call parse_string(buf, num_word, word)
                if(word(1) == 'RESI') then
                    resi_count = resi_count+1
                end if
            end if
        end do

        this%m_num_restype = resi_count
!        if( .not. allocated(this%m_residue) ) then   !< Check if it is already allocated
!            allocate(this%m_residue(this%m_num_restype), stat=mem_stat)  !< Allocate residue container
!            if(mem_stat /= 0) then
!                print*, 'MEMORY allocation error: We will return right now'
!                return
!            end if
!        else if(size(this%m_residue) /= this%m_num_restype) then  !< Check size of allocated memory
!            deallocate(this%m_residue)
!            allocate(this%m_residue(this%m_num_restype), stat=mem_stat)  !< Allocate residue container
!            if(mem_stat /= 0) then
!                print*, 'MEMORY allocation error: We will return right now'
!                return
!            end if
!        end if

        !> NOTE : Now we will read in the residue templates.
        resi_count = 1
        tabchar = ACHAR(9) !! Ascii code 9 is tab.
        delimiters = ' '//tabchar//' '

        j = 0
        do while (j .LT. NUM_RESIDUE_DEF_STRS)
            j = j + 1
            buf = residue_def_strs(j)
!            print*, 'buf = ', buf
            if(buf(1:1) == '#') then
                cycle
            else if(Is_Blank_Line(buf)) then
                cycle
            else
                call parse(buf, delimiters, word, num_word)
                if(word(1) == 'END') then
                    exit
                else if(word(1) == 'RESI') then
                    this%m_residue(resi_count)%name = word(2)
                    read(word(3), '(I2)') this%m_residue(resi_count)%num_atom
                    read(word(4), '(I2)') this%m_residue(resi_count)%num_bond
                    read(word(5), '(F6.3)') this%m_residue(resi_count)%total_charge
                    this%m_residue(resi_count)%desc = TRIM(buf(32:))

                    num_atom = this%m_residue(resi_count)%num_atom

                    !< get information for atoms
                    do i=1, num_atom
                        !read(unit_no, '(A120)', iostat=ios) buf
                        j = j + 1
                        buf = residue_def_strs(j)
!                        print*, 'buf = ', buf
                        call parse(buf, delimiters, word, num_word)
                        this%m_residue(resi_count)%atom(i)%name = trim(word(2))
                        this%m_residue(resi_count)%atom(i)%type = trim(word(3))
                        this%m_residue(resi_count)%atom(i)%xtype = trim(word(4))
                        read(word(5), '(F6.3)') this%m_residue(resi_count)%atom(i)%r
                        read(word(6), '(F6.3)') this%m_residue(resi_count)%atom(i)%eps
                        read(word(7), '(F6.3)') this%m_residue(resi_count)%atom(i)%q
                        this%m_residue(resi_count)%atom(i)%hb = trim(word(8))
                        read(word(9), '(F6.3)') this%m_residue(resi_count)%atom(i)%logp
                        read(word(10), '(F6.3)') this%m_residue(resi_count)%atom(i)%solv
                        read(word(11), '(I1)') this%m_residue(resi_count)%atom(i)%ring
                        this%m_residue(resi_count)%atom(i)%pmftype = trim(word(12))
                    end do

                    !< Then, read in the bond information
                    num_bond = this%m_residue(resi_count)%num_bond
                    if(num_bond > 0) then
                        do i=1, num_bond
                            !read(unit_no, '(A120)', iostat=ios) buf
                            j = j + 1
                            buf = residue_def_strs(j)

                            call parse_string(buf, num_word, word)
                            this%m_residue(resi_count)%bond(i)%atom1 = trim(word(2))
                            this%m_residue(resi_count)%bond(i)%atom2 = trim(word(3))
                            this%m_residue(resi_count)%bond(i)%type = trim(word(4)) !< bond type: 1, 2, 3, ar, am
                        end do
                    end if

                    resi_count = resi_count + 1
                end if
            end if
        end do !< END of residue reading
    end subroutine Param_ReadResDefInMemory

    subroutine Param_ReadAtomDefinition(this, filename)
        type(ForceField_t) :: this
        !< dummy argument
        character(len=*), intent(in) :: filename
        !< local variable
        integer :: i, num, num_atom
        character(len=120) :: buf, head
        character(len=120) :: word(30)
        character(len=30) :: delimiters
        character(len=1) :: tabchar
        integer :: num_word
        integer :: unit_no, ios, mem_stat
        logical :: lblankLine

        !!----------------------------------------------------------------------------
        !> File Open
        !!----------------------------------------------------------------------------
        unit_no = getUnitNumber(21)
        open(unit=unit_no, FILE=trim(filename), status='OLD', action='READ', iostat=ios)
        if(ios /= 0) then
            write(*, '(3A, I4)') 'FILE ', filename, ' OPEN Failed with erorr code ', ios
            return  !< TODO Here we need additional treatment for exceptional case.
        end if

        !> Count the number of total atoms
        num_atom = 0
        do
            read(unit_no, '(A80)', iostat=ios) buf
            if(ios /= 0) then
                exit
            else if(buf(1:1) == '#') then
                cycle
            else if(Is_Blank_Line(buf)) then
                cycle
            else
                num_atom = num_atom+1
            end if
        end do

        REWIND(unit_no)

        this%m_num_atomtype = num_atom
!        if( .not. allocated(this%m_atom) ) then   !< Check if it is already allocated
!            allocate(this%m_atom(this%m_num_atomtype), stat=mem_stat)  !< Allocate residue container
!            if(mem_stat /= 0) then
!                print*, 'MEMORY allocation error: We will return right now'
!                return
!            end if
!        else if(size(this%m_atom) /= this%m_num_atomtype) then  !< Check size of allocated memory
!            deallocate(this%m_atom)
!            allocate(this%m_atom(this%m_num_atomtype), stat=mem_stat)  !< Allocate residue container
!            if(mem_stat /= 0) then
!                print*, 'MEMORY allocation error: We will return right now'
!                return
!            end if
!        end if

        !> NOTE : Now we will read in the XTool Atom Definitions
        tabchar = ACHAR(9) !! Ascii code 9 is tab.
        delimiters = ' '//tabchar//' '
        num_atom = 1
        do
            read(unit_no, '(A120)', iostat=ios) buf
            if(ios /= 0) then
                exit
            else if(buf(1:1) == '#') then
                cycle
            else if(Is_Blank_Line(buf)) then
                cycle
            else
                !call parse_string(buf, num_word, word)
                word=''
                call parse(buf, delimiters, word, num_word)

                if(ios /= 0) exit
                call Param_ClearAtomDef(this%m_atom(num_atom))
                this%m_atom(num_atom)%type = trim(word(2))
                read(word(3), '(f6.2)' ) this%m_atom(num_atom)%weight
                read(word(4), '(f6.3)' ) this%m_atom(num_atom)%r
                read(word(5), '(f6.3)' ) this%m_atom(num_atom)%eps
                read(word(6), '(f6.3)' ) this%m_atom(num_atom)%q
                this%m_atom(num_atom)%hb = trim(word(7))
                this%m_atom(num_atom)%description = trim(buf(46:))
                call compact(this%m_atom(num_atom)%description)
                call delall(this%m_atom(num_atom)%description,tabchar)
                this%m_atom(num_atom)%vdwTypeNumber_ = -1 ! It is not assigend yet and will be assigned in another step.

                !> Now increase atom number
                num_atom = num_atom + 1

                if(m_lDebug) then
                    print*, 'LINE is => ', buf
                end if
            end if
        end do !< END of atom reading

        close(unit_no)
    end subroutine Param_ReadAtomDefinition

    subroutine Param_ReadAtomDefInMemory(this)
        type(ForceField_t) :: this
        integer :: i, num, num_atom
        character(len=120) :: buf, head
        character(len=120) :: word(30)
        character(len=30) :: delimiters
        character(len=1) :: tabchar
        integer :: num_word
        integer :: mem_stat
        logical :: lblankLine

        !> Count the number of total atoms
        num_atom = 0
        do i=1, NUM_ATOM_DEFS
            buf = atom_defs(i)
            if(buf(1:1) == '#') then
                cycle
            else if(Is_Blank_Line(buf)) then
                cycle
            else
                num_atom = num_atom+1
            end if
        end do

        this%m_num_atomtype = num_atom
!        if( .not. allocated(this%m_atom) ) then   !< Check if it is already allocated
!            allocate(this%m_atom(this%m_num_atomtype), stat=mem_stat)  !< Allocate residue container
!            if(mem_stat /= 0) then
!                print*, 'MEMORY allocation error: We will return right now'
!                return
!            end if
!        else if(size(this%m_atom) /= this%m_num_atomtype) then  !< Check size of allocated memory
!            deallocate(this%m_atom)
!            allocate(this%m_atom(this%m_num_atomtype), stat=mem_stat)  !< Allocate residue container
!            if(mem_stat /= 0) then
!                print*, 'MEMORY allocation error: We will return right now'
!                return
!            end if
!        end if

        !> NOTE : Now we will read in the XTool Atom Definitions
        tabchar = ACHAR(9) !! Ascii code 9 is tab.
        delimiters = ' '//tabchar//' '
        num_atom = 1
        do i=1, NUM_ATOM_DEFS
            buf = atom_defs(i)
            if(buf(1:1) == '#') then
                cycle
            else if(Is_Blank_Line(buf)) then
                cycle
            else
                !call parse_string(buf, num_word, word)
                word=''
                call parse(buf, delimiters, word, num_word)
                !print*, 'buf = ', buf

                call Param_ClearAtomDef(this%m_atom(num_atom))
                this%m_atom(num_atom)%type = trim(word(2))
                read(word(3), '(f6.2)' ) this%m_atom(num_atom)%weight
                read(word(4), '(f6.3)' ) this%m_atom(num_atom)%r
                read(word(5), '(f6.3)' ) this%m_atom(num_atom)%eps
                read(word(6), '(f6.3)' ) this%m_atom(num_atom)%q
                this%m_atom(num_atom)%hb = trim(word(7))
                this%m_atom(num_atom)%description = trim(buf(46:))
                call compact(this%m_atom(num_atom)%description)
                call delall(this%m_atom(num_atom)%description,tabchar)
                this%m_atom(num_atom)%vdwTypeNumber_ = -1 ! It is not assigend yet and will be assigned in another step.

                !> Now increase atom number
                num_atom = num_atom + 1

                if(m_lDebug) then
                    print*, 'LINE is => ', buf
                end if
            end if
        end do !< END of atom reading
    end subroutine Param_ReadAtomDefInMemory

    subroutine Param_ReadXAtomDef(this, filename)
        type(ForceField_t) :: this
        !< dummy argument
        character(len=*), intent(in) :: filename
        !< local variable
        integer :: i, num, num_atom
        character(len=120) :: buf, head
        character(len=120) :: word(20)
        character(len=30) :: delimiters
        character(len=1) :: tabchar

        integer :: num_word
        integer :: unit_no, ios, mem_stat
        logical :: lblankLine

        !!----------------------------------------------------------------------------
        !> File Open
        !!----------------------------------------------------------------------------
        unit_no = getUnitNumber(21)
        open(unit=unit_no, FILE=trim(filename), status='OLD', action='READ', iostat=ios)
        if(ios /= 0) then
            write(*, '(3A, I4)') 'FILE ', filename, ' OPEN Failed with erorr code ', ios
            return  !< TODO Here we need additional treatment for exceptional case.
        end if

        !> Count the number of total atoms
        num_atom = 0
        do
            read(unit_no, '(A80)', iostat=ios) buf
            if(ios /= 0) then
                exit
            else if(buf(1:1) == '#') then
                cycle
            else if(Is_Blank_Line(buf)) then
                cycle
            else
                num_atom = num_atom+1
            end if
        end do

        REWIND(unit_no)

        this%m_num_xatomtype = num_atom
!        if( .not. allocated(this%m_xatom) ) then   !< Check if it is already allocated
!            allocate(this%m_xatom(this%m_num_xatomtype), stat=mem_stat)  !< Allocate residue container
!            if(mem_stat /= 0) then
!                print*, 'MEMORY allocation error: We will return right now'
!                return
!            end if
!        else if(size(this%m_xatom) /= this%m_num_xatomtype) then  !< Check size of allocated memory
!            deallocate(this%m_xatom)
!            allocate(this%m_xatom(this%m_num_xatomtype), stat=mem_stat)  !< Allocate residue container
!            if(mem_stat /= 0) then
!                print*, 'MEMORY allocation error: We will return right now'
!                return
!            end if
!        end if


        !> NOTE : Now we will read in the residue templates.
        tabchar = ACHAR(9) !! Ascii code 9 is tab.
        delimiters = ' '//tabchar//' '
        num_atom = 1
        do
            read(unit_no, '(A120)', iostat=ios) buf
            if(ios /= 0) then
                exit
            else if(buf(1:1) == '#') then
                cycle
            else if(Is_Blank_Line(buf)) then
                cycle
            else
                !call parse_string(buf, num_word, word)
                word=''
                call parse(buf, delimiters, word, num_word)

                if(ios /= 0) exit
                this%m_xatom(num_atom)%type = trim(word(2))
                this%m_xatom(num_atom)%hb = trim(word(3))
                read(word(4), '(f6.3)' ) this%m_xatom(num_atom)%logp

                if(m_lDebug) then
                    write(*, '(A20, 2X, F6.3)') this%m_xatom(num_atom)%type, this%m_xatom(num_atom)%logp
                end if

                !> Now increase atom number
                num_atom = num_atom + 1

            end if
        end do !< END of atom reading

        close(unit_no)

    end subroutine Param_ReadXAtomDef

    subroutine Param_ReadXAtomDefInMem(this)
        type(ForceField_t) :: this
        integer :: i, num, num_atom
        character(len=120) :: buf, head
        character(len=120) :: word(20)
        character(len=30) :: delimiters
        character(len=1) :: tabchar

        integer :: num_word
        integer :: mem_stat
        logical :: lblankLine

        !> Count the number of total atoms
        num_atom = 0
        do i=1, NUM_XATOM_DEFS
            buf = xatom_defs(i)
            if(buf(1:1) == '#') then
                cycle
            else if(Is_Blank_Line(buf)) then
                cycle
            else
                num_atom = num_atom+1
            end if
        end do

        this%m_num_xatomtype = num_atom
!        if( .not. allocated(this%m_xatom) ) then   !< Check if it is already allocated
!            allocate(this%m_xatom(this%m_num_xatomtype), stat=mem_stat)  !< Allocate residue container
!            if(mem_stat /= 0) then
!                print*, 'MEMORY allocation error: We will return right now'
!                return
!            end if
!        else if(size(this%m_xatom) /= this%m_num_xatomtype) then  !< Check size of allocated memory
!            deallocate(this%m_xatom)
!            allocate(this%m_xatom(this%m_num_xatomtype), stat=mem_stat)  !< Allocate residue container
!            if(mem_stat /= 0) then
!                print*, 'MEMORY allocation error: We will return right now'
!                return
!            end if
!        end if


        !> NOTE : Now we will read in the residue templates.
        tabchar = ACHAR(9) !! Ascii code 9 is tab.
        delimiters = ' '//tabchar//' '
        num_atom = 1
        do i=1, NUM_XATOM_DEFS
            buf = xatom_defs(i)
            if(buf(1:1) == '#') then
                cycle
            else if(Is_Blank_Line(buf)) then
                cycle
            else
                !call parse_string(buf, num_word, word)
                word=''
                call parse(buf, delimiters, word, num_word)
!                print*, 'buf = ', buf

                this%m_xatom(num_atom)%type = trim(word(2))
                read(word(3), '(f6.3)' ) this%m_xatom(num_atom)%logp
                this%m_xatom(num_atom)%hb = trim(word(4))

                if(m_lDebug) then
                    write(*, '(A20, 2X, F6.3)') this%m_xatom(num_atom)%type, this%m_xatom(num_atom)%logp
                end if

                !> Now increase atom number
                num_atom = num_atom + 1

            end if
        end do !< END of atom reading

    end subroutine Param_ReadXAtomDefInMem

    subroutine Param_Get_Atom_LogP(this, atom_type, logp)
        type(ForceField_t) :: this
        !> dummy argument which will be replaced as actual argument upon execution.
        character(len=*), intent(in) :: atom_type
        real(m_dp), intent(out) :: logp

        !> local variable
        character(len=20) :: query_atom
        integer :: i
        logical :: mark

        query_atom = trim(ADJUSTL(atom_type))
        if(this%m_num_xatomtype < 0) then
            logp = -9999.0 ! Abnormal case
            print*, 'ERROR: this%m_xatom was not set up yet'
            return
        end if

        mark = .false.
        TYPE_LOOP : do i=1, this%m_num_xatomtype
            if(m_lDebug) then
                write(*, '(A15, 2X, A15)') query_atom, this%m_xatom(i)%type
            endif

            if(query_atom == this%m_xatom(i)%type) then
                logp = this%m_xatom(i)%logp
                mark = .true.
                exit TYPE_LOOP
            end if
        end do TYPE_LOOP

        if(.not. mark) then
            print*, 'Warning : no logP parameter for atom type ', query_atom, ' ...'
            print*, 'Zero value assigned'
            logp = 0.0
        end if

    end subroutine Param_Get_Atom_LogP

    subroutine Param_AssignAll_Atoms(this, atoms, num_atom)
        type(ForceField_t) :: this
        type(Atom_t), intent(inout) :: atoms(:)
        integer, intent(in) :: num_atom
        logical :: assignResult
        integer i
        do i=1, num_atom
            call Param_Assign_OneAtom(this, atoms(i), assignResult)
        enddo
    end subroutine Param_AssignAll_Atoms

    subroutine Param_Assign_OneAtom(this, atm, assignResult)
        type(ForceField_t) :: this
        type(Atom_t), intent(inout) :: atm
        logical, intent(out) :: assignResult
        integer i

        assignResult = .false.

        !> Notice that partial charge is no assigned here
        !! in order not to overwrite the original charge on this atom.
        ATOM_TYPE_LOOP : do i=1, this%m_num_atomtype
            if(atm%xtype /= this%m_atom(i)%type) then
                cycle
            else
                atm%weight = this%m_atom(i)%weight
                atm%R = this%m_atom(i)%r
                atm%r = this%m_atom(i)%r
                atm%eps = this%m_atom(i)%eps
                atm%hb = this%m_atom(i)%hb
                atm%valid = .true.
                atm%vdwTypeNumber_ = this%m_atom(i)%vdwTypeNumber_
                !atm%transformedVdwTypeNumber_ = atm%GetVdwTypeNumberForPotential()
                assignResult = .true.
                EXIT ATOM_TYPE_LOOP
            endif
        enddo ATOM_TYPE_LOOP
        if(.not. assignResult) then
            print*, atm%id, " (*) Assign_Atom_Parameters ERROR : atom%xtype = ", atm%xtype
        else
        !            print*, atm%id, " (#) type = ", atm%type, " xtype = ", atm%xtype, &
        !                    " xlogptype = " , atm%xlogptype
            !print*, "(#) atom%type = ", atm%xtype
        endif
    end subroutine Param_Assign_OneAtom

    subroutine finalize_xd_parameter()
        !< @TODO You must implement deallocation of memory for this module.
        print*, 'DO some deallocation jobs here !!'
    end subroutine finalize_xd_parameter

    !===============================================================================
    ! Data parser
    !===============================================================================
    !-------------------------------------------------------------------------------
    !!----------------------------------------------------------------------------
    !> This was copied temporarily from galaxydock code for convenience in test.
    !!----------------------------------------------------------------------------
    subroutine parse_string(long_string, num_word, word)
        !-------------------------------------------------------------------------------
        ! Scan string, get num of "word", and save in array word
        !-------------------------------------------------------------------------------
        character(len=*), intent(in) :: long_string
        integer, intent(out) :: num_word
        character(len=120), intent(out) :: word(:)
        character(len=120) :: string
        integer :: i, n, k
        character(len=1) :: prev_char, cur_char

        n = len(long_string)
        if (n > 120) then
            write(*,"(A,I6)") 'Error. The length of the long string is greater than maximum: ', 120
            !write(log_msg,"(A,I6)") 'Error. The length of the long string is greater than maximum: ', 120
            !call terminate_with_error(log_msg)
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

    end subroutine parse_string

end module xlogp_Parameter_m
