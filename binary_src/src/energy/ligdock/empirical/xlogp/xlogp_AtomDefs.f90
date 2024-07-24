module xlogp_AtomDefs_m
    implicit none

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
    !#           D = H-bond donor
    !#           A = H-bond acceptor
    !#          DA = H-bond donor/acceptor
    !#          DH = polar hydrogen
    !#           P = polar atom
    !#           H = hydrophobic atom
    !#           N = none
    !# Note that the ligand O.co2 is set as DA; while the protein O.co2 is set as A
    !################################################################################
    integer, parameter :: NUM_ATOM_DEFS = 51
    character(len=97), parameter :: atom_defs(NUM_ATOM_DEFS) = (/ &
        & "  1  C.3          12.0100    2.1000    0.0000    0.0000  H    carbon sp3 non-polar               ", &
        & "  2  C.3.x        12.0100    2.1000    0.0000    0.0000  P    carbon sp3 polar                   ", &
        & "  3  C.3.un       12.0100    2.1000    0.0000    0.0000  P    carbon sp3 unknown                 ", &
        & "  4  C.2          12.0100    1.9000    0.0000    0.0000  H    carbon sp2 non-polar               ", &
        & "  5  C.2.x        12.0100    1.9000    0.0000    0.0000  P    carbon sp2 polar                   ", &
        & "  6  C.2.un       12.0100    1.9000    0.0000    0.0000  P    carbon sp2 unknown                 ", &
        & "  7  C.ar         12.0100    2.0000    0.0000    0.0000  H    carbon aromatic non-polar          ", &
        & "  8  C.ar.x       12.0100    2.0000    0.0000    0.0000  P    carbon aromatic polar              ", &
        & "  9  C.ar.un      12.0100    2.0000    0.0000    0.0000  P    carbon aromatic unknown            ", &
        & " 10  C.1          12.0100    1.8000    0.0000    0.0000  P    carbon sp non-polar                ", &
        & " 11  C.1.x        12.0100    1.8000    0.0000    0.0000  P    carbon sp polar                    ", &
        & " 12  C.1.un       12.0100    1.8000    0.0000    0.0000  P    carbon sp unknown                  ", &
        & " 13  C.cat        12.0100    1.9000    0.0000    1.0000  P    carboncation (C+)                  ", &
        & " 14  N.3.h        14.0100    1.8000    0.0000    0.0000  D    nitrogen sp3 with H                ", &
        & " 15  N.3          14.0100    1.8000    0.0000    0.0000  P    nitrogen sp3                       ", &
        & " 16  N.3.un       14.0100    1.8000    0.0000    0.0000  P    nitrogen sp3 unknown               ", &
        & " 17  N.pl3.h      14.0100    1.7500    0.0000    0.0000  D    nitrogen sp3 planar with H         ", &
        & " 18  N.pl3        14.0100    1.7500    0.0000    0.0000  P    nitrogen sp3 planar                ", &
        & " 19  N.pl3.un     14.0100    1.7500    0.0000    0.0000  P    nitrogen sp3 planar unknown        ", &
        & " 20  N.2.h        14.0100    1.7500    0.0000    0.0000  DA   nitrogen sp2 with H                ", &
        & " 21  N.2          14.0100    1.7500    0.0000    0.0000  A    nitrogen sp2                       ", &
        & " 22  N.2.un       14.0100    1.7500    0.0000    0.0000  P    nitrogen sp2 unknown               ", &
        & " 23  N.ar.h       14.0100    1.7500    0.0000    0.0000  D    nitrogen aromatic with H           ", &
        & " 24  N.ar         14.0100    1.7500    0.0000    0.0000  A    nitrogen aromatic                  ", &
        & " 25  N.ar.un      14.0100    1.7500    0.0000    0.0000  P    nitrogen aromatic unknown          ", &
        & " 26  N.1          14.0100    1.7500    0.0000    0.0000  A    nitrogen sp                        ", &
        & " 27  N.1.un       14.0100    1.7500    0.0000    0.0000  P    nitrogen sp unknown                ", &
        & " 28  N.4          14.0100    1.8000    0.0000    1.0000  D    nitrogen sp3 charged (N+)          ", &
        & " 29  O.3.h        16.0000    1.6500    0.0000    0.0000  DA   oxygen sp3 with H                  ", &
        & " 30  O.3          16.0000    1.6500    0.0000    0.0000  A    oxygen sp3                         ", &
        & " 31  O.3.un       16.0000    1.6500    0.0000    0.0000  P    oxygen sp3 unknown                 ", &
        & " 32  O.2          16.0000    1.5500    0.0000    0.0000  A    oxygen sp2                         ", &
        & " 33  O.2.un       16.0000    1.5500    0.0000    0.0000  P    oxygen sp2 unknown                 ", &
        & " 34  O.co2        16.0000    1.5500    0.0000   -0.5000  DA   oxygen in carboxylate              ", &
        & " 35  S.3.h        32.0700    2.1000    0.0000    0.0000  H    sulfur sp3 with H                  ", &
        & " 36  S.3          32.0700    2.1000    0.0000    0.0000  H    sulfur sp3                         ", &
        & " 37  S.3.un       32.0700    2.1000    0.0000    0.0000  P    sulfur sp3 unknown                 ", &
        & " 38  S.2          32.0700    2.0000    0.0000    0.0000  P    sulfur sp2                         ", &
        & " 39  S.2.un       32.0700    2.0000    0.0000    0.0000  P    sulfur sp2 unknown                 ", &
        & " 40  S.o          32.0700    2.0000    0.0000    0.0000  P    sulfur in sulfoxide or sulfone     ", &
        & " 41  P.3          30.9700    2.0000    0.0000    0.0000  P    phosphorous sp3                    ", &
        & " 42  F            19.0000    1.5000    0.0000    0.0000  P    fluorine                           ", &
        & " 43  Cl           35.4500    1.7500    0.0000    0.0000  H    chlorine                           ", &
        & " 44  Br           79.9000    1.9000    0.0000    0.0000  H    bromine                            ", &
        & " 45  I           126.9000    2.0500    0.0000    0.0000  H    iodine                             ", &
        & " 46  H             1.0000    1.0000    0.0000    0.0000  N    hydrogen non-polar                 ", &
        & " 47  H.hb          1.0000    1.0000    0.0000    0.0000  DH   hydrogen polar                     ", &
        & " 48  Si           28.0900    2.0000    0.0000    0.0000  N    silicon sp3                        ", &
        & " 49  O.w          16.0000    1.7500    0.0000    0.0000  DA   water                              ", &
        & " 50  M+            0.0000    1.2500    0.0000    2.0000  M    metal ions                         ", &
        & " 51  Un            0.0000    0.0000    0.0000    0.0000  N    unknown atom                       " &
        & /)

end module xlogp_AtomDefs_m
