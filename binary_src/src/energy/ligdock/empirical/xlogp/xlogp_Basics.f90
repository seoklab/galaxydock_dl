module xlogp_constants_m
    !!----------------------------------------------------------------------------
    !> Contants for X-SCORE
    !!----------------------------------------------------------------------------
    implicit none

    logical, parameter :: XSCORE_TRUE = .true.
    logical, parameter :: XSCORE_FALSE = .false.
    real, parameter :: XSCORE_PI = 3.1416
    real, parameter :: XSCORE_TINY = 1.0e-6
    real, parameter :: XSCORE_LARGE = 1.0e+6
    real, parameter :: XSCORE_WATER_R = 1.40
    real, parameter :: XSCORE_POCKET_DEPTH = 4.00
    real, parameter :: XSCORE_LAYER_DEPTH = 3.00
    integer, parameter :: XSCORE_MAX_ATOM_NEIB = 6
    integer, parameter :: XSCORE_MAX_BOND_NEIB = 10
    integer, parameter :: XSCORE_MAX_GRID_NUM = 60000
    integer, parameter :: XSCORE_MAX_ATOMNUM_IN_GRID = 400
    integer, parameter :: XSCORE_MAX_ATOM_IN_RING = 200 ! Is this correct(sufficient) number?
    integer, parameter :: XSCORE_MAX_LIGAND_ATOM_NUMBER = 400
    integer, parameter :: XSCORE_MAX_LIGAND_BOND_NUMBER = 400
    integer, parameter :: XSCORE_MAX_RING_NUMBER = 500
    integer, parameter :: XSCORE_MAX_PATHVECTOR_SIZE = 500

!        type(Residue_Def_t), DIMENSION() :: m_residue(:)
!        type(Atom_Def_t), DIMENSION() :: m_atom(:)
!        type(Xatom_Def), DIMENSION() :: m_xatom(:)
!        type(Bond_Def), DIMENSION() :: m_bond(:)
!        type(Tors_Def), DIMENSION() :: m_torsion(:)
    integer, parameter :: XSCORE_MAX_RESIDUE_DEF_NUM = 50
    integer, parameter :: XSCORE_MAX_ATOM_DEF_NUM = 200
    integer, parameter :: XSCORE_MAX_XATOM_DEF_NUM = 200
    integer, parameter :: XSCORE_MAX_BOND_DEF_NUM = 200
    integer, parameter :: XSCORE_MAX_TORS_DEF_NUM = 200

end module xlogp_constants_m

!!----------------------------------------------------------------------------
!> Atom module for X-SCORE
!!----------------------------------------------------------------------------
module xlogp_Atom_m
    use xlogp_constants_m, only : XSCORE_MAX_ATOM_NEIB
    implicit none

    private ! Information hiding

    public :: Atom_t
    public :: AtomId_t
    public :: radiusForVdwTypeNumber_
    public :: is_equal_atom_to
    public :: is_equal_atomid_to
    public :: operator (==)
    public :: makeAtomId
    public :: Atom_IsAromaticType
    public :: Atom_GetHybType
    public :: Atom_Show
    public :: xlogp_AtomId_Show
    public :: Atom_ShowHeader
    public :: xlogp_Atom_ShowNeighbor
    public :: xlogp_Atom_Clear
    public :: xlogp_AtomId_Clear

    real, parameter, dimension(13) :: radiusForVdwTypeNumber_ = (/ 0.0, &
        2.1, 2.05, 2.0, &
        1.9, 1.8, 1.750, &
        1.650, 1.550, 1.500, &
        1.250, 1.000, 0.000 /)

    integer, parameter :: totalVdwTypeNumber_ = 13

    type Atom_t
        integer :: id !  = 0             ! atom id
        logical :: valid !  = .false.          !< valid indicator
        character(len=10) :: name ! = 'Un'  !< atom name
        character(len=10) :: type != 'Un'  !< basic atom type
        integer :: vdwTypeNumber_  != -1  !<
        integer :: transformedVdwTypeNumber_ != -1 !<

        character(len=10) :: xtype  != ''   !< xtool atom type
        character(len=20) :: xlogptype  != ''   !< xlogp atom type
        character(len=10) :: type2  != 'Un'   !< another atom type when necessary
        character(len=10) :: residue != 'Un'  !< residue name
        character(len=10) :: res_id != '0'   !< residue id: must be strong for PDB
        character(len=1) :: chain != ' '     !< chain label
        real :: coor(3)    != 0.0         !< coordinates
        real :: root(3)    != 0.0         !< HB root's coordinates
        real :: weight     != 0.0         !< atomic weight
        real :: r          != 0.0         !< vdw radius
        real :: eps        != 0.0         !< vdw epsilon value
        real :: q          != 0.0         !< partial atomic charge
        real :: Rx         != 0.0         !< X radius
        real :: logp       != 0.0         !< atomic hydrophobic scale
        real :: solv       != 0.0         !< atomic solvation parameter
        character(len=3) :: hb  != 'N'    !< HB property
        real :: occupancy     != 1.0      !< occupancy probability, for protein atoms
        real :: bfactor       != 0.0      !< B-factor, for protein atoms
        real :: score         != 0.0      !< atomic binding score
        integer :: ring       != 0      !< ring indicator : 1=normal; 2=aromatic
        integer :: origin     != 1      !< atom origin indicator: 1=ligand, 2=protein
        integer :: part       != 1      !< component indicator

        !< For protein atoms: 1=ATOM; 2=HETATM
        integer :: num_neib != 0 !< number of neighboring atoms
        integer :: num_nonh != 0 !< number of non-H neighboring atoms
        integer :: neib(XSCORE_MAX_ATOM_NEIB) != 0 !< ID of neighboring atoms
        integer :: bond(XSCORE_MAX_ATOM_NEIB) != 0 !< ID of neighboring bonds

        integer :: temp != 0 !< for misc uses

    end type Atom_t

    type AtomId_t
        integer :: id != 0
        character(len=10) :: name != 'Un'  !< atom name
        character(len=10) :: type != 'Un'  !< basic atom type
        integer :: num_neib != 0
    end type AtomId_t

    interface operator (==) ! overload ==
        module procedure is_equal_atom_to
    end interface

    interface operator (==) ! overload ==
        module procedure is_equal_atomid_to
    end interface

contains

    !------------------------------------------------------------
    ! Clear the atom object
    subroutine xlogp_AtomId_Clear(this)
        type(AtomId_t) :: this

        this%id = 0
        this%name = 'Un'
        this%type = 'Un'
        this%num_neib = 0
    end subroutine xlogp_AtomId_Clear

    !------------------------------------------------------------
    ! Clear the atom object
    subroutine xlogp_Atom_Clear(this)
        type(Atom_t) :: this

        this%id = 0
        this%valid = .false.
        this%name = 'Un'
        this%type = 'Un'
        this%vdwTypeNumber_ = -1
        this%transformedVdwTypeNumber_ = -1
        this%xtype = ''
        this%xlogptype = ''
        this%type2 = 'Un'
        this%residue = 'Un'
        this%res_id = '0'
        this%chain = ' '
        this%coor = 0.0
        this%root = 0.0
        this%weight = 0.0
        this%r = 0.0
        this%eps = 0.0
        this%q = 0.0
        this%Rx = 0.0
        this%logp = 0.0
        this%solv = 0.0
        this%hb = 'N'
        this%occupancy = 1.0
        this%bfactor = 0.0
        this%score = 0.0
        this%ring = 0
        this%origin = 1
        this%part = 1
        this%num_neib = 0
        this%num_nonh = 0
        this%neib = 0
        this%bond = 0
        this%temp = 0
    end subroutine xlogp_Atom_Clear

    logical function is_equal_atomid_to(a, b) result(t_f) ! overload (==)
        type(AtomId_t), intent(in) :: a, b
        t_f = .false. ! initialize

        if (a%id /= b%id) return ! same id?
        if (a%name /= b%name) return
        if (a%type /= b%type) return
        if (a%num_neib /= b%num_neib) return
        t_f = .true.
    end function is_equal_atomid_to

    type(AtomId_t) function makeAtomId(a) result(newAtomId)
        type(Atom_t), intent(in) :: a

        newAtomId%id = a%id
        newAtomId%name = a%name
        newAtomId%type = a%type
        newAtomId%num_neib = a%num_neib

    end function makeAtomId

    logical function is_equal_atom_to(a, b) result(t_f) ! overload (==)
        type(Atom_t), intent(in) :: a, b
        t_f = .false. ! initialize

        if (a%id /= b%id) return ! same id?
        if (a%name /= b%name) return
        if (a%type /= b%type) return
        if (a%xtype /= b%xtype) return
        !        if (.not. all(a%coor == b%coor)) return
        if (a%num_neib /= b%num_neib) return
        if (a%num_nonh /= b%num_nonh) return
        if (.not. all(a%neib == b%neib)) return
        if (.not. all(a%bond == b%bond)) return
        t_f = .true.
    end function is_equal_atom_to

    subroutine xlogp_AtomId_Show(this)
        type(AtomId_t) :: this

10      FORMAT ('AtomId: ', I5, 2X, A10, 2X, A10, 2X, I5, 2X)
        write (*, 10) this%id, this%name, this%type, this%num_neib

    end subroutine xlogp_AtomId_Show

    subroutine Atom_Show(this)
        type(Atom_t) :: this

10      FORMAT ('Atom: ', I5, 2X, A10, 2X, &
            I1, 2X, F8.4, 2X, F8.4, 2X, F8.4, ' | ', &
            A10, ' | ', A, ' | ', A, ' | ', &
            I5, 2X, I5, 2X, &
            I5, 2X, F8.4, 2X, I1, 2X, F8.4, 2X, A3, 1X)

        write (*, 10) this%id, this%name, &
            this%num_neib, this%coor(1), this%coor(2), this%coor(3), &
            this%type, this%xtype, this%xlogptype, &
            this%vdwTypeNumber_, this%transformedVdwTypeNumber_, &
            this%part, this%weight, this%ring, this%logp, this%hb

    end subroutine Atom_Show

    subroutine Atom_ShowHeader(this)
        implicit none
        type(Atom_t) :: this

        write (*, *) "     ", "id  ", "name ", &
            "num_neib ", "coor(1) ", "coor(2) ", "coor(3) ", &
            "type ", "xtype ", "xlogptype ", &
            "vdwTypeNumber_ ", "transformedVdwTypeNumber_ ", &
            "part ", "weight ", "ring ", "xlogp ", "hb "

    end subroutine Atom_ShowHeader

    integer function Atom_GetDonorType(this)
        type(Atom_t) :: this

        Atom_GetDonorType = 1
    end function Atom_GetDonorType

    integer function Atom_GetAcceptorType(this)
        type(Atom_t) :: this

        Atom_GetAcceptorType = 1
    end function Atom_GetAcceptorType

    integer function GetVdwTypeNumberForPotential(this)
        type(Atom_t) :: this
        integer :: r

        r = -1
        select case (this%vdwTypeNumber_)
            case (1,2,3,35,36,37)
                r = 1 ! vdw Radius 2.1
            case (45)
                r = 2 ! vdw Radius 2.5
            case (7,8,9,38,39,40,41,48)
                r = 3 ! vdw R 2.0
            case (4,5,6,13,44)
                r = 4 ! vdw R 1.9
            case (10,11,12,14,15,16,28)
                r = 5 ! vdw R 1.8
            case (17,18,19,20,21,22,23,24,25,26,27,43,49)
                r = 6 ! vdw R 1.750
            case (29,30,31)
                r = 7 ! vdw R 1.550
            case (32,33,34)
                r = 8 ! vdw R 1.550
            case (42)
                r = 9 ! vdw R 1.500
            case (50)
                r = 19 ! vdw R 1.250
            case (46,47)
                r = 11 ! vdw R 1.0
            case default
                print*, "GetVdwTypeNumberForPotential() Unknown atom !"
        end select
        GetVdwTypeNumberForPotential = r
    end function GetVdwTypeNumberForPotential

    integer function GetVdwTypeNumFromVdwR(this)
        type(Atom_t) :: this

        GetVdwTypeNumFromVdwR = 12
    end function GetVdwTypeNumFromVdwR

    subroutine xlogp_Atom_ShowNeighbor(this)
        type(Atom_t), intent(in) :: this
        integer :: i

10      FORMAT ('Atom: ', I5, 2X, A10, 2X, &
            I5, 2X, &
            I5, 2X, I5, 2X, I5, 2X, &
            I5, 2X, I5, 2X, I5, 2X, F8.4, 2X)

        write (*, 10) this%id, this%name, &
            this%num_neib, &
            (this%neib(i),i=1,6), this%weight

    end subroutine xlogp_Atom_ShowNeighbor

    subroutine xlogp_Atom_SetXType(this, xtype)
        type(Atom_t) :: this
        character(len=*), intent(in) :: xtype
        this%xtype = xtype
    end subroutine xlogp_Atom_SetXType

    function xlogp_Atom_GetXType(this) result(return_value)
        type(Atom_t) :: this
        character(len=20) :: return_value
        return_value = this%xtype
    end function xlogp_Atom_GetXType

    subroutine xlogp_Atom_SetXLogPType(this, xlogptype)
        type(Atom_t) :: this
        character(len=*), intent(in) :: xlogptype
        this%xlogptype = xlogptype
    end subroutine xlogp_Atom_SetXLogPType

    function xlogp_Atom_GetXLogPType(this) result(return_value)
        type(Atom_t) :: this
        character(len=20) :: return_value
        return_value = this%xlogptype
    end function xlogp_Atom_GetXLogPType

    subroutine xlogp_Atom_SetType(this, type)
        type(Atom_t) :: this
        character(len=*), intent(in) :: type
        this%type = type
    end subroutine xlogp_Atom_SetType

    function xlogp_Atom_GetType(this) result(return_value)
        type(Atom_t) :: this
        character(len=10) :: return_value
        return_value = this%type
    end function xlogp_Atom_GetType

    function Atom_GetHybType(this) result(return_value)
        type(Atom_t), intent(in) :: this
        integer :: return_value
        integer, parameter :: SP1=1
        integer, parameter :: SP2=2
        integer, parameter :: SP3=3
        integer, parameter :: SP_NONE=4

        if(this%type == 'C.3') then
            return_value = SP3
        else if(this%type == 'C.2') then
            return_value = SP2
        else if(this%type == 'C.1') then
            return_value = SP1
        else if(this%type == 'C.cat') then ! carbocation
            return_value = SP2
        else if(this%type == 'C.ar') then
            return_value = SP2
        else if(this%type == 'H') then
            return_value = SP_NONE
        else if(this%type == 'H.spc') then
            return_value = SP_NONE
        else if(this%type == 'N.4') then
            return_value = SP3
        else if(this%type == 'N.3') then
            return_value = SP3
        else if(this%type == 'N.2') then
            return_value = SP2
        else if(this%type == 'N.1') then
            return_value = SP1
        else if(this%type == 'N.ar') then
            return_value = SP2
        else if(this%type == 'N.pl3') then ! nitrogen trigonal planar
            return_value = SP2
        else if(this%type == 'N.am') then
            return_value = SP2
        else if(this%type == 'O.3') then
            return_value = SP3
        else if(this%type == 'O.2') then
            return_value = SP2
        else if(this%type == 'O.co2') then
            return_value = SP2
        else if(this%type == 'P.3') then
            return_value = SP3
        else if(this%type == 'S.3') then
            return_value = SP3
        else if(this%type == 'S.2') then
            return_value = SP2
        else if(this%type == 'S.o') then
            return_value = SP3
        else if(this%type == 'S.o2') then
            return_value = SP3
        else if(this%type == 'F') then
            return_value = SP_NONE
        else if(this%type == 'Cl') then
            return_value = SP_NONE
        else if(this%type == 'Br') then
            return_value = SP_NONE
        else if(this%type == 'I') then
            return_value = SP_NONE
        else if(this%type == 'Si') then
            return_value = SP3
        else
            return_value = SP_NONE
        endif
    end function Atom_GetHybType

    function Atom_IsAromaticType(this) result(return_value)
        type(Atom_t) :: this
        logical :: return_value

        if(this%type == 'C.ar' .OR. &
            this%type == 'N.ar' .OR. &
            this%type == 'N.pl3') then
            return_value = .true.
        else
            return_value = .false.
        endif
    end function Atom_IsAromaticType

end module xlogp_Atom_m



!!----------------------------------------------------------------------------
!> Bond module for X-SCORE
!!----------------------------------------------------------------------------
module xlogp_Bond_m
    use xlogp_constants_m, only : XSCORE_MAX_BOND_NEIB
    implicit none

    !------------------------------------------------------------
    ! Bond Derived Type
    type Bond_t

        !sequence ! Let this bond type stored as sequential

        integer :: id != 0! atom id
        integer :: valid != 0 ! valid indicator, 1 (Normal Single Bond), 2 (Rotor Single Bond),
                            ! 3 (Any other valid bond, double, triple, aromatic ...)
        integer :: atom_1 != 0 ! ID of the atom_1
        integer :: atom_2 != 0 ! ID of the atom_2
        character(len=3) :: type != '1' ! basic bond type

        integer :: part != 1 ! component indicator
        integer :: ring != 0 ! ring indicator : 1=normal; 2=aromatic
        real :: length != 0.000 ! bond length
        !! For protein atoms: 1=ATOM; 2=HETATM
        integer :: num_neib != 0 ! number of neighboring bonds
        integer :: neib(XSCORE_MAX_BOND_NEIB) != 0 ! ID of neighboring bonds

    end type Bond_t

contains

    !------------------------------------------------------------
    ! Clear the bond object
    subroutine xlogp_Bond_Clear(b)
        implicit none
        type(Bond_t), intent(inout) :: b

        b%id = 0
        b%valid = 0
        b%part = 1
        b%ring = 0
        b%atom_1 = 0
        b%atom_2 = 0
        b%type = '1'
        b%length = 0.000
        b%num_neib = 0
        b%neib = 0

    end subroutine xlogp_Bond_Clear

    subroutine xlogp_Bond_ShowContents(b)
        implicit none
        type(Bond_t), intent(in) :: b
        integer :: i

10      FORMAT ('Bond: ', I5, 2X, A3, 2X, I1, 2X, I1, 2X, ' between : ', &
            I5, 2X, ' - ', 2X, I5, ' | ', &
            I5, 2X, ' | ', &
            I5, 2X, I5, 2X, I5, 2X, &
            I5, 2X, I5, 2X, I5, 2X)

        write (*, 10) b%id, b%type, b%ring, b%valid, b%atom_1, b%atom_2, &
            b%num_neib, &
            (b%neib(i),i=1,6)

    end subroutine xlogp_Bond_ShowContents

    subroutine xlogp_Bond_ShowHeader(this)
        implicit none
        type(Bond_t) :: this
        integer :: i

        write (*, *) "     ", "ID ", "type  ", "ring ", "valid ", &
            "between atom_1 - ", "atom_2 | ", "num_neib | ", (i,i=1,6)

    end subroutine xlogp_Bond_ShowHeader

    subroutine xlogp_Bond_CheckValidity(this)
        type(Bond_t), intent(inout) :: this

        if(this%atom_1 /= this%atom_2) then
            this%valid = 1
        else
            this%valid = 0
        end if
    end subroutine xlogp_Bond_CheckValidity

    function xlogp_Bond_GetValid(this) result(return_value)
        type(Bond_t), intent(in) :: this
        integer :: return_value

        return_value = this%valid
    end function xlogp_Bond_GetValid



end module xlogp_Bond_m

!!----------------------------------------------------------------------------
!> Dot for X-SCORE
!!----------------------------------------------------------------------------
module xlogp_Dot_m

    !------------------------------------------------------------
    ! Dot Derived Type
    type Dot_t

        !sequence ! Let this Dot type stored as sequential

        integer :: valid ! valid indicator
        character(len=10) :: type ! basic Dot type
        real :: coor(3) ! Coordinates
        real :: unit ! contribution, in either A^2 or A^3
        real :: score ! score on this dot

    end type Dot_t

contains

    !------------------------------------------------------------
    ! Clear the Dot object
    subroutine DotClear(b)
        type(Dot_t), intent(inout) :: b



    end subroutine DotClear

    subroutine DotShowContents(b)
        implicit none
        type(Dot_t), intent(in) :: b

    end subroutine DotShowContents


end module xlogp_Dot_m


!!----------------------------------------------------------------------------
!> Hydrogen Bond for X-SCORE
!!----------------------------------------------------------------------------
module xlogp_HBond_m
    use xlogp_Atom_m
    implicit none

    !------------------------------------------------------------
    ! HBond Derived Type
    type HBond_t

        !sequence ! Let this HBond type stored as sequential

        logical :: valid != .false. ! valid indicator
        integer :: type != 0 !
        ! 1 : donor = latom, acceptor = patom
        ! 2 : donor = patom, acceptor = latom
        ! 3 : donor = metal, acceptor = latom/patom
        ! 4 : donor = latom, acceptor = latom
        ! 5 : donor = patom, acceptor = patom
        ! 0 : invalid, no H-bond
        integer :: d_type != 0 ! donor type: 1=straight, 2=angled
        integer :: a_type != 0 ! acceptor type: 1=straight, 2=angled
        logical :: sb != .false. ! flag for neutral HB or SB (?)
        integer :: latom != 0 ! id of the ligand atom
        integer :: patom != 0 ! id of the protein atom

        type(Atom_t) :: atomH, atomD, atomA
        real :: d != 0.0 ! D-A distance
        real :: a0 != 0.0 ! D-H-A angle
        real :: a1 != 0.0 ! DR-D-A angle
        real :: a2 != 0.0 ! D-A-AR angle

        real :: score != 0.0 ! strength of this HBond


    end type HBond_t

contains

    !------------------------------------------------------------
    ! Clear the HBond object
    subroutine HBondClear(hb)
        implicit none
        type(HBond_t), intent(inout) :: hb
        hb%valid = .false. ; hb%latom = 0; hb%patom = 0
        hb%type = 0; hb%d_type = 0; hb%a_type = 0;
        hb%score = 0
        hb%sb = .false.
        hb%d = 0.0; hb%a0 = 0.0; hb%a1 = 0.0;  hb%a2 = 0.0
        call xlogp_Atom_Clear(hb%atomD) ! AtomClear(hb%atomD)
        call xlogp_Atom_Clear(hb%atomH) ! AtomClear(hb%atomH)
        call xlogp_Atom_Clear(hb%atomA) ! AtomClear(hb%atomA)

    end subroutine HBondClear

    subroutine HBondShowContents(hb)
        implicit none
        type(HBond_t), intent(in) :: hb

10      FORMAT ('HBond type: ', I5, 3X, '----------------------------------------')
        write(*, 10) hb%type
20      FORMAT ('Hydrogen: ', I5, 2X, 5A, 2X, F8.3, 2X, F8.3, 2X, F8.3)
        write(*, 20) hb%atomH%id, hb%atomH%xtype, hb%atomH%coor(1), hb%atomH%coor(2), hb%atomH%coor(3)
30      FORMAT ('Donor: ', I5, 2X, 5A, 2X, F8.3, 2X, F8.3, 2X, F8.3)
        write(*, 30) hb%atomD%id, hb%atomD%xtype, hb%atomD%coor(1), hb%atomD%coor(2), hb%atomD%coor(3)
40      FORMAT ('Acceptor: ', I5, 2X, 5A, 2X, F8.3, 2X, F8.3, 2X, F8.3)
        write(*, 40) hb%atomA%id, hb%atomA%xtype, hb%atomA%coor(1), hb%atomA%coor(2), hb%atomA%coor(3)

        !    30 FORMAT ('Donor: ', I5, 2X, 5A, 2X, F8.3, 2X, F8.3, 2X, F8.3)
        write(*, '("Donor: ", I5, 2X, 5A, 2X, F8.3, 2X, F8.3, 2X, F8.3)') &
            hb%atomD%id, hb%atomD%xtype, hb%atomD%coor(1), hb%atomD%coor(2), hb%atomD%coor(3)

        print '("D-A distance = ", F5.2)', hb%d
        print '("D-H-A distance = ", F5.2)', hb%a0
        print '("DR-D-A distance = ", F5.2)', hb%a1
        print '("D-A-AR distance = ", F5.2)', hb%a2
        print '("Score distance = ", F5.2)', hb%score

    end subroutine HBondShowContents


    integer function Value_HBond(hb)
        implicit none
        type(HBond_t), intent(inout) :: hb

        Value_HBond = 1
    end function Value_HBond

    integer function Value_HBond2(hb)
        implicit none
        type(HBond_t), intent(inout) :: hb

        Value_HBond2 = 1
    end function Value_HBond2

    integer function Value_SBond(hb)
        implicit none
        type(HBond_t), intent(inout) :: hb

        Value_SBond = 1
    end function Value_SBond

    subroutine Determine_DA_Type(hb)
        implicit none
        type(HBond_t), intent(inout) :: hb

    end subroutine Determine_DA_Type


end module xlogp_HBond_m

!!----------------------------------------------------------------------------
!> Residue for X-SCORE
!!----------------------------------------------------------------------------
module xlogp_Residue_m
    use xlogp_constants_m
    implicit none

    real, parameter, dimension(13) :: radiusForVdwTypeNumber_ = (/ 0.0, &
        2.1, 2.05, 2.0, &
        1.9, 1.8, 1.750, &
        1.650, 1.550, 1.500, &
        1.250, 1.000, 0.000 /)

    integer, parameter :: totalVdwTypeNumber_ = 13

    !------------------------------------------------------------
    ! Residue Derived Type
    type Residue_t

        integer :: id
    end type Residue_t

contains

    !------------------------------------------------------------
    ! Clear the Residue object
    subroutine ResidueClear(b)
        type(Residue_t), intent(inout) :: b

    end subroutine ResidueClear

    subroutine ResidueShowContents(b)
        implicit none
        type(Residue_t), intent(in) :: b


    end subroutine ResidueShowContents

end module xlogp_Residue_m



!!----------------------------------------------------------------------------
!> Ring for X-SCORE
!!----------------------------------------------------------------------------
module xlogp_Ring_m
    use xlogp_constants_m, only : XSCORE_MAX_ATOM_IN_RING
    implicit none

    !------------------------------------------------------------
    ! Ring Derived Type
    type Ring_t
        logical :: valid != .false. ! valid indicator
        integer :: type != 0 ! 1 = normal; 2 = aromatic
        integer :: num_member != 0
        integer :: atom_id(XSCORE_MAX_ATOM_IN_RING) != 0 ! atom id in this ring
        integer :: bond_id(XSCORE_MAX_ATOM_IN_RING) != 0 ! bond id in this ring
        real :: centroid(3) != 0.0
    end type Ring_t

contains

    !------------------------------------------------------------
    ! Clear the Ring object
    subroutine xlogp_Ring_Clear(r)
        implicit none
        type(Ring_t), intent(inout) :: r

        r%valid = .false.
        r%type = 0
        r%num_member = 0
        r%atom_id = 0
        r%bond_id = 0

        r%centroid = 0.0

    end subroutine xlogp_Ring_Clear

    !    subroutine CalcCentroid(this)
    !        implicit none
    !        class(Ring_t), intent(in) :: this
    !        integer :: i
    !        integer :: atomId, bondId
    !
    !        do i=1, this%num_member
    !
    !        enddo
    !    end subroutine CalcCentroid

    subroutine xlogp_Ring_Show(r)
        implicit none
        type(Ring_t), intent(in) :: r
        integer :: i

        write(*, FMT="(A, 1X, I3, A)") '[Ring size =', r%num_member, ']'
        if(r%num_member > 0) then
            write(*, FMT="(A, 1X)", ADVANCE="NO") 'atom id   = '
            do i=1, r%num_member
                write(*, FMT="(I3, 2X)", advance="NO") r%atom_id(i)
            enddo
            write(*,*) ''
            write(*, FMT="(A, 1X)", ADVANCE="NO") 'bond id   = '
            do i=1, r%num_member
                write(*, FMT="(I3, 2X)", advance="NO") r%bond_id(i)
            enddo
            write(*,*) ''
        endif

    end subroutine xlogp_Ring_Show

    function xlogp_Ring_GetNumMember(this) result(return_value)
        type(Ring_t), intent(in) :: this
        integer :: return_value
        return_value = this%num_member
    end function xlogp_Ring_GetNumMember


end module xlogp_Ring_m

!!----------------------------------------------------------------------------
!> Torsion for X-SCORE
!!----------------------------------------------------------------------------
module xlogp_Torsion_m
    use xlogp_Atom_m
    implicit none

    !------------------------------------------------------------
    ! Torsion Derived Type
    type Torsion_t

        !sequence ! Let this Torsion type stored as sequential

        type(Atom_t) :: atom_1
        type(Atom_t) :: atom_2
        type(Atom_t) :: atom_3
        type(Atom_t) :: atom_4
        character(len=3) :: type != '' ! basic Torsion type
        integer :: angle != 0 ! torsion angle, in degree
        real :: V_potential_barrier != 0.0 ! potential barrier
        integer :: n != 0 ! periodicity
        integer :: S != 0 ! sign
        real :: e != 0.0 ! torsion energy

    end type Torsion_t

contains

    !------------------------------------------------------------
    ! Clear the Torsion object
    subroutine TorsionClear(b)
        type(Torsion_t), intent(inout) :: b

        b%type = ''
        b%angle = 0
        b%e = 0.0
        b%V_potential_barrier = 0.000
        b%n = 0
        b%S = 0

    end subroutine TorsionClear

    subroutine TorsionShowContents(b)
        implicit none
        type(Torsion_t), intent(in) :: b

        print '("Torsion ", A3, 2X, " between : ")', b%type
        print '(T1, "Atom 1", T1, I5, T1, A5, T1, A5)', &
            b%atom_1%id, b%atom_1%name, b%atom_1%type
        print '(T1, "Atom 2", T1, I5, T1, A5, T1, A5)', &
            b%atom_2%id, b%atom_2%name, b%atom_2%type
        print '(T1, "Atom 3", T1, I5, T1, A5, T1, A5)', &
            b%atom_3%id, b%atom_3%name, b%atom_3%type
        print '(T1, "Atom 4", T1, I5, T1, A5, T1, A5)', &
            b%atom_4%id, b%atom_4%name, b%atom_4%type

    end subroutine TorsionShowContents

end module xlogp_Torsion_m

