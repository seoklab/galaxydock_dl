module xlogp_XLogPTyper_m
    use xlogp_constants_m, only : XSCORE_MAX_ATOM_NEIB
    use xlogp_Atom_m
    use xlogp_Bond_m
    use xlogp_Group_m

    implicit none

!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!SYBYL Atom Type
!
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
!
! H hydrogen
! H.spc hydrogen in Single Point Charge (SPC) water model
! H.t3p hydrogen in Transferable intermolecular Potential (TIP3P) water model
!
! LP lone pair
! Du dummy atom
! Du.C dummy carbon
! Any any atom
! Hal halogen
! Het Heteroatom = N, O, S, P
! Hev heavy atom (non hydrogen)
! Li lithium
! Na sodium
! Mg magnesium
! Al aluminum
! Si silicon
! K  potassium
! Ca calcium
! Cr.th chromium (tetrahedral)
! Cr.oh chromium (octahedral)
! Mn manganese
! Fe iron
! Co.oh cobalt (octahedral)
! Cu copper

!! Sybyl Bond Types
! 1 = Single
! 2 = Double
! 3 = Triple
! am = amide
! ar = aromatic
! du = dummy
! un = unknown (cannot be determined from the parameter tables)
! nc = not connected
!!
    type(Group_t) :: g_GroupList(200)

contains

    subroutine ShowAllFoundGroups(num_atom)
        integer :: num_atom
        integer :: i
        do i=1, num_atom
            call xlogp_Group_ShowContents(g_GroupList(i))
        enddo
    end subroutine ShowAllFoundGroups

    subroutine FindXLogPType(atomList, bondList, atomId, aType)
        type(Atom_t), DIMENSION(:), target, intent(inout) :: atomList
        type(Bond_t), DIMENSION(:), target, intent(inout) :: bondList
        integer, intent(in) :: atomId
        character(len=20), intent(out) :: aType
        !! end of input arguments
        type(Group_t) group
        type(Atom_t) atom
        !! end of declaration of local variables

        atom = atomList(atomId)
        call xlogp_Group_FindXGroup(group, atomList, bondList, atomId)
        !! Copy group
        g_GroupList(atomId) = group

        aType = 'Un'

        !> Check if it is a hydrogen type
        !!----------------------------------------------------------------------------
        ! Hydrogen
        if (atom%type == 'H' .OR. atom%type == 'H.spc') then
            if(group%neib(1)%type(1:1) == 'O') then
                aType = 'H.hb'
            else if(group%neib(1)%type(1:1) == 'N') then
                aType = 'H.hb'
            else
                aType = 'H'
            end if
        end if

        !> Check if it is a carbon atom type

        !!----------------------------------------------------------------------------
        ! C.3
        if(atom%type == 'C.3') then
            select case (group%num_nonh)
                case (1)    !! num_nonh == 1
                    if(group%num_hetero == 0) then
                        if(group%num_pi == 0) then
                            aType = 'C.3.h3.pi=0'
                        else
                            aType = 'C.3.h3.pi=1'
                        end if
                    else !! num_hetero > 0
                        aType = 'C.3.h3.x'
                    end if
                case (2) !! num_nonh == 2
                    if(group%num_hetero == 0) then
                        if(group%num_pi == 0) then
                            aType = 'C.3.h2.pi=0'
                        else if(group%num_pi == 1) then
                            aType = 'C.3.h2.pi=1'
                        else
                            aType = 'C.3.h2.pi=2'
                        end if
                    else !! num_hetero > 0
                        if(group%num_pi == 0) then
                            aType = 'C.3.h2.x.pi=0'
                        else if(group%num_pi == 1) then
                            aType = 'C.3.h2.x.pi=1'
                        else
                            aType = 'C.3.h2.x.pi=2'
                        end if
                    end if
                case (3) !! num_nonh == 3
                    if(group%num_hetero == 0) then
                        if(group%num_pi == 0) then
                            aType = 'C.3.h.pi=0'
                        else if(group%num_pi == 1) then
                            aType = 'C.3.h.pi=1'
                        else
                            aType = 'C.3.h.pi>1'
                        end if
                    else !! num_hetero > 0
                        if(group%num_pi == 0) then
                            aType = 'C.3.h.x.pi=0'
                        else if(group%num_pi == 1) then
                            aType = 'C.3.h.x.pi=1'
                        else
                            aType = 'C.3.h.x.pi>1'
                        end if
                    end if

                case (4) !! num_nonh == 4
                    if(group%num_hetero == 0) then
                        if(group%num_pi == 0) then
                            aType = 'C.3.pi=0'
                        else if(group%num_pi == 1) then
                            aType = 'C.3.pi=1'
                        else
                            aType = 'C.3.pi>1'
                        end if
                    else !! num_hetero > 0
                        if(group%num_pi == 0) then
                            aType = 'C.3.x.pi=0'
                        else
                            aType = 'C.3.x.pi>0'
                        end if
                    end if
                case default
                    aType = 'C.3.unknown'
            end select
        end if

        !!----------------------------------------------------------------------------
        ! C.2
        if(atom%type == 'C.2') then
            select case (group%num_nonh)
                case (1) !! num_nonh == 1
                    aType = 'C.2.h2'
                case (2) !! num_nonh == 2
                    if (group%num_hetero == 0) then
                        if(group%num_pi == 0) then
                            aType = 'C.2.h.pi=0'
                        else
                            aType = 'C.2.h.pi=1'
                        end if
                    else !! num_hetero > 0
                        if(group%num_pi == 0) then
                            aType = 'C.2.h.x.pi=0'
                        else
                            aType = 'C.2.h.x.pi=1'
                        end if
                    end if
                case (3)    !! num_nonh == 3
                    if (group%num_hetero == 0) then
                        if(group%num_pi == 0) then
                            aType = 'C.2.pi=0'
                        else
                            aType = 'C.2.pi>0'
                        end if
                    else if(group%num_hetero == 1) then !! num_hetero == 1
                        if(group%num_pi == 0) then
                            aType = 'C.2.x.pi=0'
                        else
                            aType = 'C.2.x.pi>0'
                        end if
                    else
                        if(group%num_pi == 0) then
                            aType = 'C.2.x2.pi=0'
                        else
                            aType = 'C.2.x2.pi>0'
                        end if
                    end if
                case default
                    aType = 'C.2.unknown'
            end select
        end if

        !!----------------------------------------------------------------------------
        ! C.cat
        if(atom%type == 'C.cat') then
            aType = 'C.2.x2.pi>0'
        end if

        !!----------------------------------------------------------------------------
        ! C.ar
        if(atom%type == 'C.ar') then
            select case (group%num_nonh)
                case (2)
                    if(group%num_nar == 0) then
                        aType = 'C.ar.h'
                    else    !! num_nar > 0
                        aType = 'C.ar.h.(X)'
                    end if
                case (3)
                    if(group%num_nar == 0) then
                        if(group%num_hetero == 0) then
                            aType = 'C.ar'
                        else
                            aType = 'C.ar.x'
                        end if
                    else    !! num_nar > 0
                        if(group%num_hetero == 0) then
                            aType = 'C.ar.(X)'
                        else
                            aType = 'C.ar.(X).x'
                        end if
                    end if
                case default
                    aType = 'C.ar.unknown'
            end select
        end if

        !!----------------------------------------------------------------------------
        ! C.1
        if(atom%type == 'C.1') then
            if(group%db_type /= 0) then
                aType = 'C.1.=='
            else if(group%num_nonh == 1) then
                aType = 'C.1.h'
            else if(group%num_nonh == 2) then
                aType = 'C.1'
            else
                aType = 'C.1.unknown'
            end if
        end if


        !!----------------------------------------------------------------------------
        ! N.4 or N.3
        if(atom%type == 'N.4' .OR. &
            atom%type == 'N.3' .OR. &
            atom%type == 'N.pl3') then
            select case (group%num_nonh)
                case (1)
                    if(group%num_hetero == 0) then
                        if(group%num_pi == 0) then
                            aType = 'N.3.h2.pi=0'
                        else
                            aType = 'N.3.h2.pi=1'
                        end if
                    else
                        aType = 'N.3.h2.x'
                    end if
                case (2)
                    if(group%num_hetero == 0) then
                        if(group%num_pi == 0) then
                            aType = 'N.3.h.pi=0'
                        else if(group%num_pi == 1) then
                            aType = 'N.3.h.pi>0'
                        else if(atom%ring == 0) then
                            aType = 'N.3.h.pi>0'
                        else
                            aType = 'N.3.h.ring'
                        end if
                    else    !! num_hetero > 0
                        if(group%num_pi <= 1) then
                            aType = 'N.3.h.x'
                        else if(atom%ring == 0) then
                            aType = 'N.3.h.x'
                        else
                            aType = 'N.3.h.x.ring'
                        end if
                    end if
                case (3) !! num_nonh == 3
                    if(group%num_hetero == 0) then
                        if (group%num_pi == 0) then
                            aType = 'N.3.pi=0'
                        else if(group%num_pi == 1) then
                            aType = 'N.3.pi>0'
                        else if(atom%ring == 0) then
                            aType = 'N.3.pi>0'
                        else
                            aType = 'N.3.ring'
                        end if
                    else
                        if (group%num_pi <= 1) then
                            aType = 'N.3.x'
                        else if (atom%ring == 0) then
                            aType = 'N.3.x'
                        else
                            aType = 'N.3.x.ring'
                        end if
                    end if
                case default
                    aType = 'N.3.unknown'
            end select
        end if

        !!----------------------------------------------------------------------------
        ! N.am
        if(atom%type == 'N.am') then
            if(group%num_nonh==1) then
                aType = 'N.am.h2'
            else if(group%num_nonh==2) then
                if(group%num_hetero==0) then
                    aType = 'N.am.h'
                else
                    aType = 'N.am.h.x'
                 end if
            else if(group%num_nonh==3) then
               if(group%num_hetero==0) then
                  aType = 'N.am'
               else
                  aType = 'N.am.x'
               end if
            else
                aType = 'N.am.unknown'
            end if
        end if

        !!----------------------------------------------------------------------------
        ! N.2
        if(atom%type == 'N.2') then
            select case (group%db_type)
                case (1,4,5) !! N=C, N=S, N=P
                    if(group%num_hetero == 0) then
                        if(group%num_pi == 0) then
                            aType = 'N.2.(=C).pi=0'
                        else
                            aType = 'N.2.(=C).pi=1'
                        end if
                    else
                        if(group%num_pi == 0) then
                            aType = 'N.2.(=C).x.pi=0'
                        else
                            aType = 'N.2.(=C).x.pi=1'
                        end if
                    end if
                case (2)
                    if(group%num_hetero==0) then
                        aType = 'N.2.(=N)'
                    else
                        aType = 'N.2.(=N).x'
                    end if
                case (3)
                    if(group%num_nonh==2) then
                        aType = 'N.2.o'
                    else if(group%num_nonh==3) then
                        aType = 'N.2.o2'
                    else
                        aType = 'N.2.o'
                    end if
                case default
                    aType = 'N.2.unknown'
            end select
        end if

        if(atom%type == 'N.ar') aType = 'N.ar'
        if(atom%type == 'N.1') aType = 'N.1'

        !!----------------------------------------------------------------------------
        ! O.3
        if(atom%type == 'O.3') then
            select case (group%num_nonh)
                case (1)
                    if(group%num_hetero==0) then
                        if(group%num_pi==0) then
                            aType = 'O.3.h.pi=0'
                        else
                            aType = 'O.3.h.pi=1'
                        end if
                    else
                        aType = 'O.3.h.x'
                    end if
                case (2)
                    if(group%num_hetero==0) then
                        if(group%num_pi==0) then
                            aType = 'O.3.pi=0'
                        else
                            aType = 'O.3.pi>0'
                        end if
                    else
                        aType = 'O.3.x'
                    end if
                case default
                    aType = 'O.3.unknown'
            end select
        end if

        if(atom%type == 'O.2') aType = 'O.2'
        if(atom%type == 'O.co2') aType = 'O.co2'

        !!----------------------------------------------------------------------------
        ! Sulfur
        if(atom%type == 'S.3') then
            if(group%num_nonh==1) then
                aType = 'S.3.h'
            else if(group%num_nonh==2) then
                aType = 'S.3'
            else
                aType = 'S.3.unknown'
            end if
        end if

        if(atom%type=='S.2') then
            aType = 'S.2'
        else if(atom%type=='S.o') then
            aType = 'S.o'
        else if(atom%type=='S.o2') then
            aType = 'S.o2'
        endif

        !!----------------------------------------------------------------------------
        ! Phosphorus
        if(atom%type == 'P.3') then
            if(group%db_type == 3) then
                aType = 'P.3.(=O)'
            else if(group%db_type==4) then
                aType = 'P.3.(=S)'
            else
                aType = 'P.3.unknown'
            end if
        end if

        !!----------------------------------------------------------------------------
        ! Fluorine
        if(atom%type == 'F') then
            if(group%num_pi == 0) then
                aType = 'F.pi=0'
            else if(group%num_pi==1) then
                aType = 'F.pi=1'
            else
                aType = 'F.unknown'
            end if
        end if

        !!----------------------------------------------------------------------------
        ! Chlorine
        if(atom%type == 'Cl') then
            if(group%num_pi == 0) then
                aType = 'Cl.pi=0'
            else if(group%num_pi==1) then
                aType = 'Cl.pi=1'
            else
                aType = 'Cl.unknown'
            end if
        end if

        !!----------------------------------------------------------------------------
        ! Bromine
        if(atom%type == 'Br') then
            if(group%num_pi == 0) then
                aType = 'Br.pi=0'
            else if(group%num_pi==1) then
                aType = 'Br.pi=1'
            else
                aType = 'Br.unknown'
            end if
        end if

        !!----------------------------------------------------------------------------
        ! Iodine
        if(atom%type == 'I') then
            if(group%num_pi == 0) then
                aType = 'I.pi=0'
            else if(group%num_pi==1) then
                aType = 'I.pi=1'
            else
                aType = 'I.unknown'
            end if
        end if

        !!----------------------------------------------------------------------------
        ! Silicon
        if(atom%type == 'Si') aType = 'Si'

    end subroutine FindXLogPType

end module xlogp_XLogPTyper_m
