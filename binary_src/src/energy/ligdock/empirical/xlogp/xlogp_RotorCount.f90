module xlogp_RotorCount_m
    use globals, only: dp
    use xlogp_constants_m, only : XSCORE_MAX_ATOM_NEIB
    use xlogp_Atom_m
    use xlogp_Bond_m
    use xlogp_Group_m
    use xlogp_Ring_m
    use xlogp_AllRingsFinder_m, only : AllRingsFinder_t
    implicit none
    private

    public :: CountRotor

contains

    !!----------------------------------------------------------------------------
    !> Count the number of rotatable bonds (ROTOR).
    !! If a single bond is normal, bond.valid = 1.
    !! If a single bond is a rotor, bond.valid = 2.
    !!----------------------------------------------------------------------------
    subroutine CountRotor(atoms, num_atom, bonds, num_bond, sumRotor)
        type(Atom_t), intent(inout) :: atoms(:)
        integer, intent(in) :: num_atom
        type(Bond_t), intent(inout) :: bonds(:)
        integer, intent(in) :: num_bond
        real(dp), intent(inout) :: sumRotor
        integer :: i
        integer :: mark, tmp
        integer :: id1, id2
        integer, parameter :: NORMAL_SINGLE_BOND = 1
        integer, parameter :: ROTOR_SINGLE_BOND = 2

        !> Clean atom's scores
        do i=1, num_atom
            atoms(i)%score = 0.0
        enddo

        !> Eliminate all the bonds in ring and all the non-single bonds
        do i=1, num_bond
            if(bonds(i)%ring > 0) then
                cycle
            else if (bonds(i)%type /= '1') then
                cycle
            else
                bonds(i)%valid = ROTOR_SINGLE_BOND  ! If bonds(i)%valid = 2, then the bond is rotor
            endif
        enddo

        !> Eliminate all the R-H, R-OH, R-NH2, R-CH3 bonds
        do i=1, num_bond
            if (bonds(i)%valid /= ROTOR_SINGLE_BOND) cycle
            id1 = bonds(i)%atom_1
            id2 = bonds(i)%atom_2

            !! Filter next cases among ROTOR_SINGLE_BOND
            if(atoms(id1)%num_nonh == 1 .OR. atoms(id2)%num_nonh == 1) then
                bonds(i)%valid = NORMAL_SINGLE_BOND
            else
                cycle
            endif
        enddo

        !> sp2-sp2 rotors
        do i=1, num_bond
            if(bonds(i)%valid /= ROTOR_SINGLE_BOND) cycle
            id1 = bonds(i)%atom_1
            id2 = bonds(i)%atom_2
            mark = 0
            tmp = Atom_GetHybType(atoms(id1))
            if (tmp == 1 .OR. tmp == 2) then
                mark = mark + 1
            endif

            tmp = Atom_GetHybType(atoms(id2))
            if (tmp == 1 .OR. tmp == 2) then
                mark = mark + 1
            endif

            if (mark == 2) then
                bonds(i)%valid = NORMAL_SINGLE_BOND
                cycle
            endif
        enddo

        !> Eliminate terminal rotors, e.g. -PO3, -CF3, -CMe3, -NMe3
        do i=1, num_bond
            if (bonds(i)%valid /= ROTOR_SINGLE_BOND) cycle
            id1 = bonds(i)%atom_1
            id2 = bonds(i)%atom_2

            if(IsTerminalAtom(atoms(id1), id2, atoms)) then
                bonds(i)%valid = NORMAL_SINGLE_BOND
            else if(IsTerminalAtom(atoms(id2), id1, atoms)) then
                bonds(i)%valid = NORMAL_SINGLE_BOND
            else
                cycle
            endif
        enddo

        !> Eliminate abnormal rotors
        do i=1, num_bond
            if(bonds(i)%valid /= ROTOR_SINGLE_BOND) cycle
            id1 = bonds(i)%atom_1
            id2 = bonds(i)%atom_2

            if (atoms(id1)%valid .EQV. .false.) then
                bonds(i)%valid = NORMAL_SINGLE_BOND
            else if(atoms(id2)%valid .EQV. .false.) then
                bonds(i)%valid = NORMAL_SINGLE_BOND
            else
                cycle
            endif
        enddo

        !> Now count the frozen rotors, all the rotors have been labeled as 2
        sumRotor = 0.0
        do i=1, num_bond
            if(bonds(i)%valid /= ROTOR_SINGLE_BOND) then
                cycle
            else
                sumRotor = sumRotor + 1.0
            endif
        enddo

    end subroutine CountRotor

    function IsTerminalAtom(atm, partner, atoms) result(return_value)
        type(Atom_t), intent(in) :: atm
        integer, intent(in) :: partner
        type(Atom_t), intent(in) :: atoms(:)
        logical :: return_value
        integer i, j

        if(Atom_GetHybType(atm) /= 3) then
            return_value = .false.
            return
        endif

        if(atm%num_nonh /= 4) then
            return_value = .false.
            return
        endif

        do i=1, atm%num_nonh
            if(atm%neib(i) == partner) then
                cycle
            else if (atoms(atm%neib(i))%num_nonh /= 1) then
                return_value = .false.
                return
            else
                cycle
            endif
        enddo

        do i=1, atm%num_nonh-1
            if(atm%neib(i) == partner) cycle

            do j=i+1, atm%num_nonh
                if(atm%neib(j) == partner) then
                    cycle
                else if (atoms(atm%neib(j))%xtype == atoms(atm%neib(i))%xtype) then
                    cycle
                else
                    return_value = .false.
                    return
                endif
            enddo
        enddo

        return_value = .true.
    end function IsTerminalAtom

end module xlogp_RotorCount_m
