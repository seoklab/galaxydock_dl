module xlogp_XToolTyper_m
    use xlogp_constants_m, only : XSCORE_MAX_ATOM_NEIB
    use xlogp_Atom_m
    use xlogp_Bond_m
    use xlogp_Group_m

    implicit none

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

!################################################################################
!# Atom types and their parameters defined in XTOOL
!# Latest update: 11/13/2002
!################################################################################
!# The 1st column: atom type number
!# The 2nd column: X-Tool atom type
!# The 3rd column: atomic weight
!# The 4th column: X-Tool van der Waals radius
!# The 5th column: X-Tool van der Waals potential
!# The 6th column: X-Tool atomic partial charge
!# The 7th column: hydrogen bonding character
!#              D = H-bond donor
!#              A = H-bond acceptor
!#             DA = H-bond donor/acceptor
!#             DH = polar hydrogen
!#              P = polar atom
!#              H = hydrophobic atom
!#              N = none
!# Note that the ligand O.co2 is set as DA; while the protein O.co2 is set as A
!################################################################################
!1    C.3       12.01  2.100  0.000  0.000  H    carbon sp3 non-polar
!2    C.3.x     12.01  2.100  0.000  0.000  P    carbon sp3 polar
!3    C.3.un    12.01  2.100  0.000  0.000  P    carbon sp3 unknown
!4    C.2       12.01  1.900  0.000  0.000  H    carbon sp2 non-polar
!5    C.2.x     12.01  1.900  0.000  0.000  P    carbon sp2 polar
!6    C.2.un    12.01  1.900  0.000  0.000  P    carbon sp2 unknown
!7    C.ar      12.01  2.000  0.000  0.000  H    carbon aromatic non-polar
!8    C.ar.x    12.01  2.000  0.000  0.000  P    carbon aromatic polar
!9    C.ar.un   12.01  2.000  0.000  0.000  P    carbon aromatic unknown
!10   C.1       12.01  1.800  0.000  0.000  P    carbon sp non-polar
!11   C.1.x     12.01  1.800  0.000  0.000  P    carbon sp polar
!12   C.1.un    12.01  1.800  0.000  0.000  P    carbon sp unknown
!13   C.cat     12.01  1.900  0.000  1.000  P    carboncation (C+)
!14   N.3.h     14.01  1.800  0.000  0.000  D    nitrogen sp3 with H
!15   N.3       14.01  1.800  0.000  0.000  P    nitrogen sp3
!16   N.3.un    14.01  1.800  0.000  0.000  P    nitrogen sp3 unknown
!17   N.pl3.h   14.01  1.750  0.000  0.000  D    nitrogen sp3 planar with H
!18   N.pl3     14.01  1.750  0.000  0.000  P    nitrogen sp3 planar
!19   N.pl3.un  14.01  1.750  0.000  0.000  P    nitrogen sp3 planar unknown
!20   N.2.h     14.01  1.750  0.000  0.000  DA   nitrogen sp2 with H
!21   N.2       14.01  1.750  0.000  0.000  A    nitrogen sp2
!22   N.2.un    14.01  1.750  0.000  0.000  P    nitrogen sp2 unknown
!23   N.ar.h    14.01  1.750  0.000  0.000  D    nitrogen aromatic with H
!24   N.ar      14.01  1.750  0.000  0.000  A    nitrogen aromatic
!25   N.ar.un   14.01  1.750  0.000  0.000  P    nitrogen aromatic unknown
!26   N.1       14.01  1.750  0.000  0.000  A    nitrogen sp
!27   N.1.un    14.01  1.750  0.000  0.000  P    nitrogen sp unknown
!28   N.4       14.01  1.800  0.000  1.000  D    nitrogen sp3 charged (N+)
!29   O.3.h     16.00  1.650  0.000  0.000  DA   oxygen sp3 with H
!30   O.3       16.00  1.650  0.000  0.000  A    oxygen sp3
!31   O.3.un    16.00  1.650  0.000  0.000  P    oxygen sp3 unknown
!32   O.2       16.00  1.550  0.000  0.000  A    oxygen sp2
!33   O.2.un    16.00  1.550  0.000  0.000  P    oxygen sp2 unknown
!34   O.co2     16.00  1.550  0.000 -0.500  DA   oxygen in carboxylate
!35   S.3.h     32.07  2.100  0.000  0.000  H    sulfur sp3 with H
!36   S.3       32.07  2.100  0.000  0.000  H    sulfur sp3
!37   S.3.un    32.07  2.100  0.000  0.000  P    sulfur sp3 unknown
!38   S.2       32.07  2.000  0.000  0.000  P    sulfur sp2
!39   S.2.un    32.07  2.000  0.000  0.000  P    sulfur sp2 unknown
!40   S.o       32.07  2.000  0.000  0.000  P    sulfur in sulfoxide or sulfone
!41   P.3       30.97  2.000  0.000  0.000  P    phosphorous sp3
!42   F         19.00  1.500  0.000  0.000  P    fluorine
!43   Cl        35.45  1.750  0.000  0.000  H    chlorine
!44   Br        79.90  1.900  0.000  0.000  H    bromine
!45   I        126.90  2.050  0.000  0.000  H    iodine
!46   H          1.00  1.000  0.000  0.000  N    hydrogen non-polar
!47   H.hb       1.00  1.000  0.000  0.000  DH   hydrogen polar
!48   Si        28.09  2.000  0.000  0.000  N    silicon sp3
!49   O.w       16.00  1.750  0.000  0.000  DA   water
!50   M+         0.00  1.250  0.000  2.000  M    metal ions
!51   Un         0.00  0.000  0.000  0.000  N    unknown atom
!################################################################################

contains

    subroutine FindXToolType(atomList, bondList, atomId, aType)
        type(Atom_t), DIMENSION(:), target, intent(inout) :: atomList
        type(Bond_t), DIMENSION(:), target, intent(inout) :: bondList
        integer, intent(in) :: atomId
        character(len=20), intent(out) :: aType
        !! end of input arguments
        type(Group_t) group
        type(Atom_t) atom
        !! end of declaration of local variables

        atom = atomList(atomId)
        call xlogp_Group_FindAGroup(group, atomList, bondList, atomId)

        aType = 'Un'

        !> Check if it is a hydrogen type
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
        if(atom%type == 'C.3') then
            if (group%num_hetero == 0) then
                aType = 'C.3'
            else if(group%num_hetero > 0) then
                aType = 'C.3.x'
            else
                aType = 'C.3.un'
            end if
        end if

        if(atom%type == 'C.2' .AND. atom%ring /= 2) then
            if (group%num_hetero == 0) then
                aType = 'C.2'
            else if(group%num_hetero > 0) then
                aType = 'C.2.x'
            else
                aType = 'C.2.un'
            end if
        end if

        if(atom%type == 'C.ar' .OR. ( atom%type == 'C.2' .AND. atom%ring == 2)) then
            if(group%num_hetero == 0) then
                aType = 'C.ar'
            else if(group%num_hetero > 0) then
                aType = 'C.ar.x'
            else
                aType = 'C.ar.un'
            end if
        end if

        if(atom%type == 'C.1') then
            if(group%num_hetero == 0) then
                aType = 'C.1'
            else if(group%num_hetero > 0) then
                aType = 'C.1.x'
            else
                aType = 'C.1.un'
            end if
        end if

        if(atom%type == 'C.cat') then
            aType = 'C.cat'
        end if

        !> Nitrogen
        if(atom%type == 'N.4' .OR. atom%type == 'N.3') then
            if (group%num_nonh <= 2) then
                aType = 'N.4'
            else if(group%num_nonh == 3) then
                aType = 'N.3'
            else
                aType = 'N.3.un'
            end if
        end if

        if( (atom%type == 'N.am' .AND. atom%ring /= 2) .OR. &
            (atom%type == 'N.pl3' .AND. atom%ring /= 2) ) then
            if(group%num_nonh == 1) then
                aType = 'N.pl3.h'
            else if(group%num_nonh == 2) then
                aType = 'N.pl3.h'
            else if(group%num_nonh == 3) then
                aType = 'N.pl3'
            else
                aType = 'N.pl3.un'
            end if
        end if

        if( atom%type == 'N.2' .AND. atom%ring /= 2) then
            if(group%num_nonh == 1) then
                aType = 'N.2.h'
            else if(group%num_nonh == 2) then
                aType = 'N.2'
            else
                aType = 'N.2.un'
            end if
        end if

        if ( atom%type == 'N.ar' .OR. &
            (atom%type == 'N.2' .AND. atom%ring == 2) .OR. &
            (atom%type == 'N.pl3' .AND. atom%ring == 2) .OR. &
            (atom%type == 'N.am' .AND. atom%ring == 2) ) then
            if (group%num_h == 1) then
                aType = 'N.ar.h'
            else if(group%num_h == 0) then
                aType = 'N.ar'
            else
                aType = 'N.ar.un'
            end if
        end if

        if (atom%type == 'N.1') then
            if (group%num_nonh == 1) then
                aType = 'N.1'
            else
                aType = 'N.1.un'
            end if
        end if

        !> Oxygen
        if (atom%type == 'O.3') then
            if (group%num_nonh == 1) then
                aType = 'O.3.h'
            else if(group%num_nonh == 2) then
                aType = 'O.3'
            else
                aType = 'O.3.un'
            end if
        end if

        if(atom%type == 'O.2') then
            aType = 'O.2'
        end if

        if(atom%type == 'O.co2') then
            aType = 'O.co2'
        end if

        !> Sulfur
        if(atom%type == 'S.3') then
            if (group%num_nonh == 1) then
                aType = 'S.3.h'
            else if(group%num_nonh == 2) then
                aType = 'S.3'
            else
                aType = 'S.3.un'
            end if
        end if

        if(atom%type == 'S.2') then
            aType = 'S.2'
        else if(atom%type == 'S.o') then
            aType = 'S.o'
        else if(atom%type == 'S.o2') then
            aType = 'S.o'
        endif

        !> Phosphorus
        if(atom%type == 'P.3') aType = 'P.3'
        if(atom%type == 'F') aType = 'F'
        if(atom%type == 'Cl') aType = 'Cl'
        if(atom%type == 'Br') aType = 'Br'
        if(atom%type == 'I') aType = 'I'
        if(atom%type == 'Si') aType = 'Si'

    end subroutine FindXToolType

end module xlogp_XToolTyper_m
