module xlogp_ResidueDefs_m
    implicit none

!The Intel compiler allows up to 511 continuation lines.
!As you know, the standard calls for 39 continuation lines at 132 characters/line
!or 5280 characters.
!
!Our 511 continuation line limit is an extension.
!
!Must be a monster statement you have.

    !################################################################################$
    !# Protein residue templates defined in X-TOOL$
    !# Latest update: 11/20/2003$
    !################################################################################$
    !# The "RESI" section are organized as:$
    !#    1st column: heading$
    !#    2nd column: residue name$
    !#    3rd column: number of atoms$
    !#    4th column: number of bonds$
    !#    5th column: total charge$
    !#    6th column: description$
    !#$
    !# The "ATOM" section are organized as:$
    !#    1st column: heading$
    !#    2nd column: atom name as in PDB file$
    !#    3rd column: basic atom type, Sybyl/Tripos$
    !#    4th column: X-Tool atom type$
    !#    5th column: X-Tool atomic vdw radius$
    !#    6th column: X-Tool atomic vdw potential$
    !#    7th column: X-Tool atomic partial charge$
    !#    8th column: hydrogen bonding character$
    !#    9th column: atomic hydrophobic scale, based on XLOGP2$
    !#    10th column: SAS solvation parameter$
    !#    11th column: ring indicator (0=not ring; 1=normal ring; 2=aromatic ring)$
    !#    12th column: PMF atom type$
    !#$
    !#    Note that O.co2 on the protein side is defined as 'A'; while the O.co2$
    !#    on the ligand side is defined as 'DA'.$
    !#$
    !# The "BOND" section are organized as:$
    !#    1st column: heading$
    !#    2nd column: name of atom 1$
    !#    3rd column: name of atom 2$
    !#    4th column: bond type (1=single; 2=double; 3=triple; ar=aromatic; am=amide)$
    !################################################################################$

    integer, parameter :: NUM_RESIDUE_DEF_STRS = 12
    character(len=77), parameter :: residue_def_strs(NUM_RESIDUE_DEF_STRS) = (/ &
        & "RESI  ACE     6      5    0.00 Acetyl                                        ", &
        & "ATOM  C      C.2    C.2.x      1.900  0.000  0.000  P   -0.030   0.023  0  CP", &
        & "ATOM  O      O.2    O.2        1.550  0.000  0.000  A   -0.399  -0.035  0  OA", &
        & "ATOM  CA     C.3    C.3        2.200  0.000  0.000  H    0.267   0.021  0  CF", &
        & "ATOM  HA1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
        & "ATOM  HA2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
        & "ATOM  HA3    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
        & "BOND  C      O      2                                                        ", &
        & "BOND  C      CA     1                                                        ", &
        & "BOND  CA     HA1    1                                                        ", &
        & "BOND  CA     HA2    1                                                        ", &
        & "BOND  CA     HA3    1                                                        " /)
!        & "RESI  ALA    10      9    0.00 Alanine                                       ", &
!        & "ATOM  N      N.am   N.pl3.h    1.750  0.000  0.000  D   -0.096  -0.052  0  ND", &
!        & "ATOM  H      H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  C      C.2    C.2.x      1.900  0.000  0.000  P   -0.030   0.023  0  CP", &
!        & "ATOM  O      O.2    O.2        1.550  0.000  0.000  A   -0.399  -0.035  0  OA", &
!        & "ATOM  CA     C.3    C.3.x      2.200  0.000  0.000  P   -0.305   0.005  0  CP", &
!        & "ATOM  HA     H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CB     C.3    C.3        2.200  0.000  0.000  H    0.528   0.021  0  CF", &
!        & "ATOM  HB1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HB2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HB3    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "BOND  N      H      1                                                        ", &
!        & "BOND  N      CA     1                                                        ", &
!        & "BOND  C      O      2                                                        ", &
!        & "BOND  C      CA     1                                                        ", &
!        & "BOND  CA     HA     1                                                        ", &
!        & "BOND  CA     CB     1                                                        ", &
!        & "BOND  CB     HB1    1                                                        ", &
!        & "BOND  CB     HB2    1                                                        ", &
!        & "BOND  CB     HB3    1                                                        ", &
!        & "RESI  ARG    24     23    1.00 Arginine                                      ", &
!        & "ATOM  N      N.am   N.pl3.h    1.750  0.000  0.000  D   -0.096  -0.052  0  ND", &
!        & "ATOM  H      H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  C      C.2    C.2.x      1.900  0.000  0.000  P   -0.030   0.023  0  CP", &
!        & "ATOM  O      O.2    O.2        1.550  0.000  0.000  A   -0.399  -0.035  0  OA", &
!        & "ATOM  CA     C.3    C.3.x      2.200  0.000  0.000  P   -0.305   0.005  0  CP", &
!        & "ATOM  HA     H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CB     C.3    C.3        2.200  0.000  0.000  H    0.358   0.021  0  CF", &
!        & "ATOM  HB1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HB2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CG     C.3    C.3        2.200  0.000  0.000  H    0.358   0.021  0  CF", &
!        & "ATOM  HG1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HG2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CD     C.3    C.3.x      2.200  0.000  0.000  P   -0.137   0.005  0  CN", &
!        & "ATOM  HD1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HD2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  NE     N.pl3  N.pl3.h    1.750  0.000  0.000  D   -0.096   0.041  0  NC", &
!        & "ATOM  HE     H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  CZ     C.cat  C.cat      1.800  0.000  0.000  P    0.005   0.023  0  CN", &
!        & "ATOM  NH1    N.pl3  N.pl3.h    1.750  0.000  0.500  D   -0.646   0.041  0  NC", &
!        & "ATOM  HH11   H      H.hb       1.000  0.000  0.000  DH   0.000   0.000  0  HH", &
!        & "ATOM  HH12   H      H.hb       1.000  0.000  0.000  DH   0.000   0.000  0  HH", &
!        & "ATOM  NH2    N.pl3  N.pl3.h    1.750  0.000  0.500  D   -0.646   0.041  0  NC", &
!        & "ATOM  HH21   H      H.hb       1.000  0.000  0.000  DH   0.000   0.000  0  HH", &
!        & "ATOM  HH22   H      H.hb       1.000  0.000  0.000  DH   0.000   0.000  0  HH", &
!        & "BOND  N      H      1                                                        ", &
!        & "BOND  N      CA     1                                                        ", &
!        & "BOND  C      O      2                                                        ", &
!        & "BOND  C      CA     1                                                        ", &
!        & "BOND  CA     HA     1                                                        ", &
!        & "BOND  CA     CB     1                                                        ", &
!        & "BOND  CB     HB1    1                                                        ", &
!        & "BOND  CB     HB2    1                                                        ", &
!        & "BOND  CB     CG     1                                                        ", &
!        & "BOND  CG     HG1    1                                                        ", &
!        & "BOND  CG     HG2    1                                                        ", &
!        & "BOND  CG     CD     1                                                        ", &
!        & "BOND  CD     HD1    1                                                        ", &
!        & "BOND  CD     HD2    1                                                        ", &
!        & "BOND  CD     NE     1                                                        ", &
!        & "BOND  NE     HE     1                                                        ", &
!        & "BOND  NE     CZ     1                                                        ", &
!        & "BOND  CZ     NH1    1                                                        ", &
!        & "BOND  CZ     NH2    1                                                        ", &
!        & "BOND  NH1    HH11   1                                                        ", &
!        & "BOND  NH1    HH12   1                                                        ", &
!        & "BOND  NH2    HH21   1                                                        ", &
!        & "BOND  NH2    HH22   1                                                        ", &
!        & "RESI  ASN    14     13    0.00 Asparagine                                    ", &
!        & "ATOM  N      N.am   N.pl3.h    1.750  0.000  0.000  D   -0.096  -0.052  0  ND", &
!        & "ATOM  H      H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  C      C.2    C.2.x      1.900  0.000  0.000  P   -0.030   0.023  0  CP", &
!        & "ATOM  O      O.2    O.2        1.550  0.000  0.000  A   -0.399  -0.035  0  OA", &
!        & "ATOM  CA     C.3    C.3.x      2.200  0.000  0.000  P   -0.305   0.005  0  CP", &
!        & "ATOM  HA     H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CB     C.3    C.3        2.200  0.000  0.000  H   -0.008   0.021  0  CF", &
!        & "ATOM  HB1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HB2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CG     C.2    C.2.x      1.900  0.000  0.000  P   -0.030   0.023  0  CP", &
!        & "ATOM  OD1    O.2    O.2        1.550  0.000  0.000  A   -0.399  -0.035  0  OA", &
!        & "ATOM  ND2    N.am   N.pl3.h    1.750  0.000  0.000  D   -0.646  -0.052  0  ND", &
!        & "ATOM  HD21   H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  HD22   H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "BOND  N      H      1                                                        ", &
!        & "BOND  N      CA     1                                                        ", &
!        & "BOND  C      O      2                                                        ", &
!        & "BOND  C      CA     1                                                        ", &
!        & "BOND  CA     HA     1                                                        ", &
!        & "BOND  CA     CB     1                                                        ", &
!        & "BOND  CB     HB1    1                                                        ", &
!        & "BOND  CB     HB2    1                                                        ", &
!        & "BOND  CB     CG     1                                                        ", &
!        & "BOND  CG     OD1    2                                                        ", &
!        & "BOND  CG     ND2    am                                                       ", &
!        & "BOND  ND2    HD21   1                                                        ", &
!        & "BOND  ND2    HD22   1                                                        ", &
!        & "RESI  ASP    12     11   -1.00 Aspartic acid                                 ", &
!        & "ATOM  N      N.am   N.pl3.h    1.750  0.000  0.000  D   -0.096  -0.052  0  ND", &
!        & "ATOM  H      H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  C      C.2    C.2.x      1.900  0.000  0.000  P   -0.030   0.023  0  CP", &
!        & "ATOM  O      O.2    O.2        1.550  0.000  0.000  A   -0.399  -0.035  0  OA", &
!        & "ATOM  CA     C.3    C.3.x      2.200  0.000  0.000  P   -0.305   0.005  0  CP", &
!        & "ATOM  HA     H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CB     C.3    C.3        2.200  0.000  0.000  H   -0.008   0.021  0  CF", &
!        & "ATOM  HB1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HB2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CG     C.2    C.2.x      1.900  0.000  0.000  P   -0.030   0.023  0  CO", &
!        & "ATOM  OD1    O.co2  O.co2      1.550  0.000 -0.500  A   -0.880   0.000  0  OC", &
!        & "ATOM  OD2    O.co2  O.co2      1.550  0.000 -0.500  A   -0.880   0.000  0  OC", &
!        & "BOND  N      H      1                                                        ", &
!        & "BOND  N      CA     1                                                        ", &
!        & "BOND  C      O      2                                                        ", &
!        & "BOND  C      CA     1                                                        ", &
!        & "BOND  CA     HA     1                                                        ", &
!        & "BOND  CA     CB     1                                                        ", &
!        & "BOND  CB     HB1    1                                                        ", &
!        & "BOND  CB     HB2    1                                                        ", &
!        & "BOND  CB     CG     1                                                        ", &
!        & "BOND  CG     OD1    2                                                        ", &
!        & "BOND  CG     OD2    2                                                        ", &
!        & "RESI  CYS    11     10    0.00 Cysteine                                      ", &
!        & "ATOM  N      N.am   N.pl3.h    1.750  0.000  0.000  D   -0.096  -0.052  0  ND", &
!        & "ATOM  H      H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  C      C.2    C.2.x      1.900  0.000  0.000  P   -0.030   0.023  0  CP", &
!        & "ATOM  O      O.2    O.2        1.550  0.000  0.000  A   -0.399  -0.035  0  OA", &
!        & "ATOM  CA     C.3    C.3.x      2.200  0.000  0.000  P   -0.305   0.005  0  CP", &
!        & "ATOM  HA     H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CB     C.3    C.3        2.200  0.000  0.000  H    0.358   0.005  0  CP", &
!        & "ATOM  HB1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HB2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  SG     S.3    S.3        2.100  0.000  0.000  H    0.419   0.037  0  SA", &
!        & "ATOM  HG     H      H          1.000  0.000  0.000  N    0.000  -0.106  0  HH", &
!        & "BOND  N      H      1                                                        ", &
!        & "BOND  N      CA     1                                                        ", &
!        & "BOND  C      O      2                                                        ", &
!        & "BOND  C      CA     1                                                        ", &
!        & "BOND  CA     HA     1                                                        ", &
!        & "BOND  CA     CB     1                                                        ", &
!        & "BOND  CB     HB1    1                                                        ", &
!        & "BOND  CB     HB2    1                                                        ", &
!        & "BOND  CB     SG     1                                                        ", &
!        & "BOND  SG     HG     1                                                        ", &
!        & "RESI  GLN    17     16    0.00 Glutamine                                     ", &
!        & "ATOM  N      N.am   N.pl3.h    1.750  0.000  0.000  D   -0.096  -0.052  0  ND", &
!        & "ATOM  H      H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  C      C.2    C.2.x      1.900  0.000  0.000  P   -0.030   0.023  0  CP", &
!        & "ATOM  O      O.2    O.2        1.550  0.000  0.000  A   -0.399  -0.035  0  OA", &
!        & "ATOM  CA     C.3    C.3.x      2.200  0.000  0.000  P   -0.305   0.005  0  CP", &
!        & "ATOM  HA     H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CB     C.3    C.3        2.200  0.000  0.000  H    0.358   0.021  0  CF", &
!        & "ATOM  HB1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HB2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CG     C.3    C.3        2.200  0.000  0.000  H   -0.008   0.021  0  CF", &
!        & "ATOM  HG1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HG2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CD     C.2    C.2.x      1.900  0.000  0.000  P   -0.030   0.023  0  CP", &
!        & "ATOM  OE1    O.2    O.2        1.550  0.000  0.000  A   -0.399  -0.035  0  OA", &
!        & "ATOM  NE2    N.am   N.pl3.h    1.750  0.000  0.000  D   -0.646  -0.052  0  ND", &
!        & "ATOM  HE21   H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  HE22   H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "BOND  N      H      1                                                        ", &
!        & "BOND  N      CA     1                                                        ", &
!        & "BOND  C      O      2                                                        ", &
!        & "BOND  C      CA     1                                                        ", &
!        & "BOND  CA     HA     1                                                        ", &
!        & "BOND  CA     CB     1                                                        ", &
!        & "BOND  CB     HB1    1                                                        ", &
!        & "BOND  CB     HB2    1                                                        ", &
!        & "BOND  CB     CG     1                                                        ", &
!        & "BOND  CG     HG1    1                                                        ", &
!        & "BOND  CG     HG2    1                                                        ", &
!        & "BOND  CG     CD     1                                                        ", &
!        & "BOND  CD     OE1    2                                                        ", &
!        & "BOND  CD     NE2    am                                                       ", &
!        & "BOND  NE2    HE21   1                                                        ", &
!        & "BOND  NE2    HE22   1                                                        ", &
!        & "RESI  GLU    15     14   -1.00 Glutamic acid                                 ", &
!        & "ATOM  N      N.am   N.pl3.h    1.750  0.000  0.000  D   -0.096  -0.052  0  ND", &
!        & "ATOM  H      H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  C      C.2    C.2.x      1.900  0.000  0.000  P   -0.030   0.023  0  CP", &
!        & "ATOM  O      O.2    O.2        1.550  0.000  0.000  A   -0.399  -0.035  0  OA", &
!        & "ATOM  CA     C.3    C.3.x      2.200  0.000  0.000  P   -0.305   0.005  0  CP", &
!        & "ATOM  HA     H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CB     C.3    C.3        2.200  0.000  0.000  H    0.358   0.021  0  CF", &
!        & "ATOM  HB1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HB2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CG     C.3    C.3        2.200  0.000  0.000  H   -0.008   0.021  0  CF", &
!        & "ATOM  HG1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HG2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CD     C.2    C.2.x      1.900  0.000  0.000  P   -0.030   0.023  0  CO", &
!        & "ATOM  OE1    O.co2  O.co2      1.550  0.000 -0.500  A   -0.880   0.000  0  OC", &
!        & "ATOM  OE2    O.co2  O.co2      1.550  0.000 -0.500  A   -0.880   0.000  0  OC", &
!        & "BOND  N      H      1                                                        ", &
!        & "BOND  N      CA     1                                                        ", &
!        & "BOND  C      O      2                                                        ", &
!        & "BOND  C      CA     1                                                        ", &
!        & "BOND  CA     HA     1                                                        ", &
!        & "BOND  CA     CB     1                                                        ", &
!        & "BOND  CB     HB1    1                                                        ", &
!        & "BOND  CB     HB2    1                                                        ", &
!        & "BOND  CB     CG     1                                                        ", &
!        & "BOND  CG     HG1    1                                                        ", &
!        & "BOND  CG     HG2    1                                                        ", &
!        & "BOND  CG     CD     1                                                        ", &
!        & "BOND  CD     OE1    2                                                        ", &
!        & "BOND  CD     OE2    2                                                        ", &
!        & "RESI  GLY     7      6    0.00 Glycine                                       ", &
!        & "ATOM  N      N.am   N.pl3.h    1.750  0.000  0.000  D   -0.096  -0.052  0  ND", &
!        & "ATOM  H      H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  C      C.2    C.2.x      1.900  0.000  0.000  P   -0.030   0.023  0  CP", &
!        & "ATOM  O      O.2    O.2        1.550  0.000  0.000  A   -0.399  -0.035  0  OA", &
!        & "ATOM  CA     C.3    C.3.x      2.200  0.000  0.000  P   -0.303   0.005  0  CP", &
!        & "ATOM  HA1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HA2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "BOND  N      H      1                                                        ", &
!        & "BOND  N      CA     1                                                        ", &
!        & "BOND  C      O      2                                                        ", &
!        & "BOND  C      CA     1                                                        ", &
!        & "BOND  CA     HA1    1                                                        ", &
!        & "BOND  CA     HA2    1                                                        ", &
!        & "RESI  HIS    18     18    0.00 Histidine (all variations)                    ", &
!        & "ATOM  N      N.am   N.pl3.h    1.750  0.000  0.000  D   -0.096  -0.052  0  ND", &
!        & "ATOM  H      H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  C      C.2    C.2.x      1.900  0.000  0.000  P   -0.030   0.023  0  CP", &
!        & "ATOM  O      O.2    O.2        1.550  0.000  0.000  A   -0.399  -0.035  0  OA", &
!        & "ATOM  CA     C.3    C.3.x      2.200  0.000  0.000  P   -0.305   0.005  0  CP", &
!        & "ATOM  HA     H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CB     C.3    C.3        2.200  0.000  0.000  H   -0.008   0.021  0  CF", &
!        & "ATOM  HB1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HB2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CG     C.2    C.ar.x     2.000  0.000  0.000  P   -0.027   0.017  2  cP", &
!        & "ATOM  ND1    N.2    N.ar.h     1.750  0.000  0.500  DA   0.135  -0.078  2  NR", &
!        & "ATOM  HD1    H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  CD2    C.2    C.ar.x     2.000  0.000  0.000  P   -0.310   0.017  2  cP", &
!        & "ATOM  HD2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CE1    C.2    C.ar.x     2.000  0.000  0.000  P   -0.310   0.017  2  cP", &
!        & "ATOM  HE1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  NE2    N.2    N.ar.h     1.750  0.000  0.500  DA   0.135  -0.078  2  NR", &
!        & "ATOM  HE2    H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "BOND  N      H      1                                                        ", &
!        & "BOND  N      CA     1                                                        ", &
!        & "BOND  C      O      2                                                        ", &
!        & "BOND  C      CA     1                                                        ", &
!        & "BOND  CA     HA     1                                                        ", &
!        & "BOND  CA     CB     1                                                        ", &
!        & "BOND  CB     HB1    1                                                        ", &
!        & "BOND  CB     HB2    1                                                        ", &
!        & "BOND  CB     CG     1                                                        ", &
!        & "BOND  CG     ND1    ar                                                       ", &
!        & "BOND  ND1    HD1    1                                                        ", &
!        & "BOND  CG     CD2    ar                                                       ", &
!        & "BOND  CD2    HD2    1                                                        ", &
!        & "BOND  ND1    CE1    ar                                                       ", &
!        & "BOND  CE1    HE1    1                                                        ", &
!        & "BOND  CD2    NE2    ar                                                       ", &
!        & "BOND  NE2    HE2    1                                                        ", &
!        & "BOND  CE1    NE2    ar                                                       ", &
!        & "RESI  ILE    19     18    0.00 Isoleucine                                    ", &
!        & "ATOM  N      N.am   N.pl3.h    1.750  0.000  0.000  D   -0.096  -0.052  0  ND", &
!        & "ATOM  H      H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  C      C.2    C.2.x      1.900  0.000  0.000  P   -0.030   0.023  0  CP", &
!        & "ATOM  O      O.2    O.2        1.550  0.000  0.000  A   -0.399  -0.035  0  OA", &
!        & "ATOM  CA     C.3    C.3.x      2.200  0.000  0.000  P   -0.305   0.005  0  CP", &
!        & "ATOM  HA     H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CB     C.3    C.3        2.200  0.000  0.000  H    0.127   0.021  0  CF", &
!        & "ATOM  HB     H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CG2    C.3    C.3        2.200  0.000  0.000  H    0.528   0.031  0  CF", &
!        & "ATOM  HG21   H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HG22   H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HG23   H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CG1    C.3    C.3        2.200  0.000  0.000  H    0.358   0.031  0  CF", &
!        & "ATOM  HG11   H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HG12   H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CD1    C.3    C.3        2.200  0.000  0.000  H    0.528   0.031  0  CF", &
!        & "ATOM  HD11   H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HD12   H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HD13   H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "BOND  N      H      1                                                        ", &
!        & "BOND  N      CA     1                                                        ", &
!        & "BOND  C      O      2                                                        ", &
!        & "BOND  C      CA     1                                                        ", &
!        & "BOND  CA     HA     1                                                        ", &
!        & "BOND  CA     CB     1                                                        ", &
!        & "BOND  CB     HB     1                                                        ", &
!        & "BOND  CB     CG2    1                                                        ", &
!        & "BOND  CG2    HG21   1                                                        ", &
!        & "BOND  CG2    HG22   1                                                        ", &
!        & "BOND  CG2    HG23   1                                                        ", &
!        & "BOND  CB     CG1    1                                                        ", &
!        & "BOND  CG1    HG11   1                                                        ", &
!        & "BOND  CG1    HG12   1                                                        ", &
!        & "BOND  CG1    CD1    1                                                        ", &
!        & "BOND  CD1    HD11   1                                                        ", &
!        & "BOND  CD1    HD12   1                                                        ", &
!        & "BOND  CD1    HD13   1                                                        ", &
!        & "RESI  LEU    19     18    0.00 Leucine                                       ", &
!        & "ATOM  N      N.am   N.pl3.h    1.750  0.000  0.000  D   -0.096  -0.052  0  ND", &
!        & "ATOM  H      H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  C      C.2    C.2.x      1.900  0.000  0.000  P   -0.030   0.023  0  CP", &
!        & "ATOM  O      O.2    O.2        1.550  0.000  0.000  A   -0.399  -0.035  0  OA", &
!        & "ATOM  CA     C.3    C.3.x      2.200  0.000  0.000  P   -0.305   0.005  0  CP", &
!        & "ATOM  HA     H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CB     C.3    C.3        2.200  0.000  0.000  H    0.358   0.021  0  CF", &
!        & "ATOM  HB1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HB2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CG     C.3    C.3        2.200  0.000  0.000  H    0.127   0.031  0  CF", &
!        & "ATOM  HG     H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CD1    C.3    C.3        2.200  0.000  0.000  H    0.528   0.031  0  CF", &
!        & "ATOM  HD11   H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HD12   H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HD13   H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CD2    C.3    C.3        2.200  0.000  0.000  H    0.528   0.031  0  CF", &
!        & "ATOM  HD21   H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HD22   H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HD23   H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "BOND  N      H      1                                                        ", &
!        & "BOND  N      CA     1                                                        ", &
!        & "BOND  C      O      2                                                        ", &
!        & "BOND  C      CA     1                                                        ", &
!        & "BOND  CA     HA     1                                                        ", &
!        & "BOND  CA     CB     1                                                        ", &
!        & "BOND  CB     HB1    1                                                        ", &
!        & "BOND  CB     HB2    1                                                        ", &
!        & "BOND  CB     CG     1                                                        ", &
!        & "BOND  CG     HG     1                                                        ", &
!        & "BOND  CG     CD1    1                                                        ", &
!        & "BOND  CD1    HD11   1                                                        ", &
!        & "BOND  CD1    HD12   1                                                        ", &
!        & "BOND  CD1    HD13   1                                                        ", &
!        & "BOND  CG     CD2    1                                                        ", &
!        & "BOND  CD2    HD21   1                                                        ", &
!        & "BOND  CD2    HD22   1                                                        ", &
!        & "BOND  CD2    HD23   1                                                        ", &
!        & "RESI  LYS    22     21    1.00 Lysine                                        ", &
!        & "ATOM  N      N.am   N.pl3.h    1.750  0.000  0.000  D   -0.096  -0.052  0  ND", &
!        & "ATOM  H      H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  C      C.2    C.2.x      1.900  0.000  0.000  P   -0.030   0.023  0  CP", &
!        & "ATOM  O      O.2    O.2        1.550  0.000  0.000  A   -0.399  -0.035  0  OA", &
!        & "ATOM  CA     C.3    C.3.x      2.200  0.000  0.000  P   -0.305   0.005  0  CP", &
!        & "ATOM  HA     H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CB     C.3    C.3        2.200  0.000  0.000  H    0.358   0.021  0  CF", &
!        & "ATOM  HB1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HB2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CG     C.3    C.3        2.200  0.000  0.000  H    0.358   0.021  0  CF", &
!        & "ATOM  HG1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HG2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CD     C.3    C.3        2.200  0.000  0.000  H    0.358   0.021  0  CF", &
!        & "ATOM  HD1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HD2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CE     C.3    C.3.x      2.200  0.000  0.000  P   -0.137   0.005  0  CN", &
!        & "ATOM  HE1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HE2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  NZ     N.4    N.4        1.800  0.000  1.000  D   -1.200   0.041  0  NC", &
!        & "ATOM  HZ1    H      H.hb       1.000  0.000  0.000  DH   0.000   0.000  0  HH", &
!        & "ATOM  HZ2    H      H.hb       1.000  0.000  0.000  DH   0.000   0.000  0  HH", &
!        & "ATOM  HZ3    H      H.hb       1.000  0.000  0.000  DH   0.000   0.000  0  HH", &
!        & "BOND  N      H      1                                                        ", &
!        & "BOND  N      CA     1                                                        ", &
!        & "BOND  C      O      2                                                        ", &
!        & "BOND  C      CA     1                                                        ", &
!        & "BOND  CA     HA     1                                                        ", &
!        & "BOND  CA     CB     1                                                        ", &
!        & "BOND  CB     HB1    1                                                        ", &
!        & "BOND  CB     HB2    1                                                        ", &
!        & "BOND  CB     CG     1                                                        ", &
!        & "BOND  CG     HG1    1                                                        ", &
!        & "BOND  CG     HG2    1                                                        ", &
!        & "BOND  CG     CD     1                                                        ", &
!        & "BOND  CD     HD1    1                                                        ", &
!        & "BOND  CD     HD2    1                                                        ", &
!        & "BOND  CD     CE     1                                                        ", &
!        & "BOND  CE     HE1    1                                                        ", &
!        & "BOND  CE     HE2    1                                                        ", &
!        & "BOND  CE     NZ     1                                                        ", &
!        & "BOND  NZ     HZ1    1                                                        ", &
!        & "BOND  NZ     HZ2    1                                                        ", &
!        & "BOND  NZ     HZ3    1                                                        ", &
!        & "RESI  MET    17     16    0.00 Methionine                                    ", &
!        & "ATOM  N      N.am   N.pl3.h    1.750  0.000  0.000  D   -0.096  -0.052  0  ND", &
!        & "ATOM  H      H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  C      C.2    C.2.x      1.900  0.000  0.000  P   -0.030   0.023  0  CP", &
!        & "ATOM  O      O.2    O.2        1.550  0.000  0.000  A   -0.399  -0.035  0  OA", &
!        & "ATOM  CA     C.3    C.3.x      2.200  0.000  0.000  P   -0.305   0.005  0  CP", &
!        & "ATOM  HA     H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CB     C.3    C.3        2.200  0.000  0.000  H    0.358   0.021  0  CF", &
!        & "ATOM  HB1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HB2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CG     C.3    C.3        2.200  0.000  0.000  H    0.358   0.005  0  CP", &
!        & "ATOM  HG1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HG2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  SD     S.3    S.3        2.100  0.000  0.000  H    0.255   0.037  0  SA", &
!        & "ATOM  CE     C.3    C.3        2.200  0.000  0.000  H    0.528   0.005  0  CP", &
!        & "ATOM  HE1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HE2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HE3    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "BOND  N      H      1                                                        ", &
!        & "BOND  N      CA     1                                                        ", &
!        & "BOND  C      O      2                                                        ", &
!        & "BOND  C      CA     1                                                        ", &
!        & "BOND  CA     HA     1                                                        ", &
!        & "BOND  CA     CB     1                                                        ", &
!        & "BOND  CB     HB1    1                                                        ", &
!        & "BOND  CB     HB2    1                                                        ", &
!        & "BOND  CB     CG     1                                                        ", &
!        & "BOND  CG     HG1    1                                                        ", &
!        & "BOND  CG     HG2    1                                                        ", &
!        & "BOND  CG     SD     1                                                        ", &
!        & "BOND  SD     CE     1                                                        ", &
!        & "BOND  CE     HE1    1                                                        ", &
!        & "BOND  CE     HE2    1                                                        ", &
!        & "BOND  CE     HE3    1                                                        ", &
!        & "RESI  NME     6      5    0.00 N-Methyl                                      ", &
!        & "ATOM  N      N.am   N.pl3.h    1.750  0.000  0.000  D   -0.096  -0.052  0  ND", &
!        & "ATOM  H      H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  CA     C.3    C.3.x      2.200  0.000  0.000  P   -0.032   0.005  0  CP", &
!        & "ATOM  HA1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HA2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HA3    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "BOND  N      H      1                                                        ", &
!        & "BOND  N      CA     1                                                        ", &
!        & "BOND  CA     HA1    1                                                        ", &
!        & "BOND  CA     HA2    1                                                        ", &
!        & "BOND  CA     HA3    1                                                        ", &
!        & "RESI  PHE    20     20    0.00 Phenylalanine                                 ", &
!        & "ATOM  N      N.am   N.pl3.h    1.750  0.000  0.000  D   -0.096  -0.052  0  ND", &
!        & "ATOM  H      H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  C      C.2    C.2.x      1.900  0.000  0.000  P   -0.030   0.023  0  CP", &
!        & "ATOM  O      O.2    O.2        1.550  0.000  0.000  A   -0.399  -0.035  0  OA", &
!        & "ATOM  CA     C.3    C.3.x      2.200  0.000  0.000  P   -0.305   0.005  0  CP", &
!        & "ATOM  HA     H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CB     C.3    C.3        2.200  0.000  0.000  H   -0.008   0.021  0  CF", &
!        & "ATOM  HB1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HB2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CG     C.ar   C.ar       2.000  0.000  0.000  H    0.296   0.017  2  cF", &
!        & "ATOM  CD1    C.ar   C.ar       2.000  0.000  0.000  H    0.337   0.017  2  cF", &
!        & "ATOM  HD1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CD2    C.ar   C.ar       2.000  0.000  0.000  H    0.337   0.017  2  cF", &
!        & "ATOM  HD2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CE1    C.ar   C.ar       2.000  0.000  0.000  H    0.337   0.017  2  cF", &
!        & "ATOM  HE1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CE2    C.ar   C.ar       2.000  0.000  0.000  H    0.337   0.017  2  cF", &
!        & "ATOM  HE2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CZ     C.ar   C.ar       2.000  0.000  0.000  H    0.337   0.017  2  cF", &
!        & "ATOM  HZ     H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "BOND  N      H      1                                                        ", &
!        & "BOND  N      CA     1                                                        ", &
!        & "BOND  C      O      2                                                        ", &
!        & "BOND  C      CA     1                                                        ", &
!        & "BOND  CA     HA     1                                                        ", &
!        & "BOND  CA     CB     1                                                        ", &
!        & "BOND  CB     HB1    1                                                        ", &
!        & "BOND  CB     HB2    1                                                        ", &
!        & "BOND  CB     CG     1                                                        ", &
!        & "BOND  CG     CD1    ar                                                       ", &
!        & "BOND  CD1    HD1    1                                                        ", &
!        & "BOND  CG     CD2    ar                                                       ", &
!        & "BOND  CD2    HD2    1                                                        ", &
!        & "BOND  CD1    CE1    ar                                                       ", &
!        & "BOND  CE1    HE1    1                                                        ", &
!        & "BOND  CD2    CE2    ar                                                       ", &
!        & "BOND  CE2    HE2    1                                                        ", &
!        & "BOND  CE1    CZ     ar                                                       ", &
!        & "BOND  CE2    CZ     ar                                                       ", &
!        & "BOND  CZ     HZ     1                                                        ", &
!        & "RESI  PRO    14     15    0.00 Proline                                       ", &
!        & "ATOM  N      N.am   N.pl3      1.750  0.000  0.000  P    0.078  -0.052  1  ND", &
!        & "ATOM  C      C.2    C.2.x      1.900  0.000  0.000  P   -0.030   0.023  0  CP", &
!        & "ATOM  O      O.2    O.2        1.550  0.000  0.000  A   -0.399  -0.035  0  OA", &
!        & "ATOM  CD     C.3    C.3.x      2.200  0.000  0.000  P   -0.137   0.005  1  CP", &
!        & "ATOM  HD1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HD2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CA     C.3    C.3.x      2.200  0.000  0.000  P   -0.305   0.005  1  CP", &
!        & "ATOM  HA     H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CB     C.3    C.3        2.200  0.000  0.000  H    0.358   0.021  1  CF", &
!        & "ATOM  HB1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HB2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CG     C.3    C.3        2.200  0.000  0.000  H    0.358   0.021  1  CF", &
!        & "ATOM  HG1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HG2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "BOND  N      H      1                                                        ", &
!        & "BOND  N      CA     1                                                        ", &
!        & "BOND  C      O      2                                                        ", &
!        & "BOND  C      CA     1                                                        ", &
!        & "BOND  CA     HA     1                                                        ", &
!        & "BOND  CA     CB     1                                                        ", &
!        & "BOND  CB     HB1    1                                                        ", &
!        & "BOND  CB     HB2    1                                                        ", &
!        & "BOND  CB     CG     1                                                        ", &
!        & "BOND  CG     HG1    1                                                        ", &
!        & "BOND  CG     HG2    1                                                        ", &
!        & "BOND  CG     CD     1                                                        ", &
!        & "BOND  CD     HD1    1                                                        ", &
!        & "BOND  CD     HD2    1                                                        ", &
!        & "BOND  N      CD     1                                                        ", &
!        & "RESI  SER    11     10    0.00 Serine                                        ", &
!        & "ATOM  N      N.am   N.pl3.h    1.750  0.000  0.000  D   -0.096  -0.052  0  ND", &
!        & "ATOM  H      H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  C      C.2    C.2.x      1.900  0.000  0.000  P   -0.030   0.023  0  CP", &
!        & "ATOM  O      O.2    O.2        1.550  0.000  0.000  A   -0.399  -0.035  0  OA", &
!        & "ATOM  CA     C.3    C.3.x      2.200  0.000  0.000  P   -0.305   0.005  0  CP", &
!        & "ATOM  HA     H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CB     C.3    C.3.x      2.200  0.000  0.000  P   -0.137   0.005  0  CP", &
!        & "ATOM  HB1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HB2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  OG     O.3    O.3.h      1.650  0.000  0.000  DA  -0.467  -0.014  0  OD", &
!        & "ATOM  HG     H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "BOND  N      H      1                                                        ", &
!        & "BOND  N      CA     1                                                        ", &
!        & "BOND  C      O      2                                                        ", &
!        & "BOND  C      CA     1                                                        ", &
!        & "BOND  CA     HA     1                                                        ", &
!        & "BOND  CA     CB     1                                                        ", &
!        & "BOND  CB     HB1    1                                                        ", &
!        & "BOND  CB     HB2    1                                                        ", &
!        & "BOND  CB     OG     1                                                        ", &
!        & "BOND  OG     HG     1                                                        ", &
!        & "RESI  THR    14     13    0.00 Threonine                                     ", &
!        & "ATOM  N      N.am   N.pl3.h    1.750  0.000  0.000  D   -0.096  -0.052  0  ND", &
!        & "ATOM  H      H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  C      C.2    C.2.x      1.900  0.000  0.000  P   -0.030   0.023  0  CP", &
!        & "ATOM  O      O.2    O.2        1.550  0.000  0.000  A   -0.399  -0.035  0  OA", &
!        & "ATOM  CA     C.3    C.3.x      2.200  0.000  0.000  P   -0.305   0.005  0  CP", &
!        & "ATOM  HA     H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CB     C.3    C.3.x      2.200  0.000  0.000  P   -0.205   0.005  0  CP", &
!        & "ATOM  HB     H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  OG1    O.3    O.3.h      1.650  0.000  0.000  DA  -0.467  -0.014  0  OD", &
!        & "ATOM  HG1    H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  CG2    C.3    C.3        2.200  0.000  0.000  H    0.528   0.021  0  CF", &
!        & "ATOM  HG21   H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HG22   H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HG23   H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "BOND  N      H      1                                                        ", &
!        & "BOND  N      CA     1                                                        ", &
!        & "BOND  C      O      2                                                        ", &
!        & "BOND  C      CA     1                                                        ", &
!        & "BOND  CA     HA     1                                                        ", &
!        & "BOND  CA     CB     1                                                        ", &
!        & "BOND  CB     HB     1                                                        ", &
!        & "BOND  CB     OG1    1                                                        ", &
!        & "BOND  OG1    HG1    1                                                        ", &
!        & "BOND  CB     CG2    1                                                        ", &
!        & "BOND  CG2    HG21   1                                                        ", &
!        & "BOND  CG2    HG22   1                                                        ", &
!        & "BOND  CG2    HG23   1                                                        ", &
!        & "RESI  TRP    24     25    0.00 Tryptophan                                    ", &
!        & "ATOM  N      N.am   N.pl3.h    1.750  0.000  0.000  D   -0.096  -0.052  0  ND", &
!        & "ATOM  H      H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  C      C.2    C.2.x      1.900  0.000  0.000  P   -0.030   0.023  0  CP", &
!        & "ATOM  O      O.2    O.2        1.550  0.000  0.000  A   -0.399  -0.035  0  OA", &
!        & "ATOM  CA     C.3    C.3.x      2.200  0.000  0.000  P   -0.305   0.005  0  CP", &
!        & "ATOM  HA     H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CB     C.3    C.3        2.200  0.000  0.000  H   -0.008   0.021  0  CF", &
!        & "ATOM  HB1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HB2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CG     C.2    C.ar       2.000  0.000  0.000  H    0.013   0.017  2  cF", &
!        & "ATOM  CD2    C.ar   C.ar       2.000  0.000  0.000  H    0.296   0.017  2  cF", &
!        & "ATOM  CE2    C.ar   C.ar.x     2.000  0.000  0.000  P   -0.151   0.017  2  cP", &
!        & "ATOM  CE3    C.ar   C.ar       2.000  0.000  0.000  H    0.337   0.017  2  cF", &
!        & "ATOM  HE3    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CD1    C.2    C.ar.x     2.000  0.000  0.000  P   -0.310   0.017  2  cP", &
!        & "ATOM  HD1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  NE1    N.pl3  N.ar.h     1.750  0.000  0.000  D    0.545  -0.078  2  ND", &
!        & "ATOM  HE1    H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  CZ2    C.ar   C.ar       2.000  0.000  0.000  H    0.337   0.017  2  cF", &
!        & "ATOM  HZ2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CZ3    C.ar   C.ar       2.000  0.000  0.000  H    0.337   0.017  2  cF", &
!        & "ATOM  HZ3    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CH2    C.ar   C.ar       2.000  0.000  0.000  H    0.337   0.017  2  cF", &
!        & "ATOM  HH2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "BOND  N      H      1                                                        ", &
!        & "BOND  N      CA     1                                                        ", &
!        & "BOND  C      O      2                                                        ", &
!        & "BOND  C      CA     1                                                        ", &
!        & "BOND  CA     HA     1                                                        ", &
!        & "BOND  CA     CB     1                                                        ", &
!        & "BOND  CB     HB1    1                                                        ", &
!        & "BOND  CB     HB2    1                                                        ", &
!        & "BOND  CB     CG     1                                                        ", &
!        & "BOND  CG     CD1    ar                                                       ", &
!        & "BOND  CD1    HD1    1                                                        ", &
!        & "BOND  CG     CD2    ar                                                       ", &
!        & "BOND  CD1    NE1    ar                                                       ", &
!        & "BOND  NE1    HE1    1                                                        ", &
!        & "BOND  NE1    CE2    ar                                                       ", &
!        & "BOND  CD2    CE2    ar                                                       ", &
!        & "BOND  CD2    CE3    ar                                                       ", &
!        & "BOND  CE3    HE3    1                                                        ", &
!        & "BOND  CE2    CZ2    ar                                                       ", &
!        & "BOND  CZ2    HZ2    1                                                        ", &
!        & "BOND  CE3    CZ3    ar                                                       ", &
!        & "BOND  CZ3    HZ3    1                                                        ", &
!        & "BOND  CZ2    CH2    ar                                                       ", &
!        & "BOND  CZ3    CH2    ar                                                       ", &
!        & "BOND  CH2    HH2    1                                                        ", &
!        & "RESI  TYR    21     21    0.00 Tyrosine                                      ", &
!        & "ATOM  N      N.am   N.pl3.h    1.750  0.000  0.000  D   -0.096  -0.052  0  ND", &
!        & "ATOM  H      H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  C      C.2    C.2.x      1.900  0.000  0.000  P   -0.030   0.023  0  CP", &
!        & "ATOM  O      O.2    O.2        1.550  0.000  0.000  A   -0.399  -0.035  0  OA", &
!        & "ATOM  CA     C.3    C.3.x      2.200  0.000  0.000  P   -0.305   0.005  0  CP", &
!        & "ATOM  HA     H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CB     C.3    C.3        2.200  0.000  0.000  H   -0.008   0.021  0  CF", &
!        & "ATOM  HB1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HB2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CG     C.ar   C.ar       2.000  0.000  0.000  H    0.296   0.017  2  cF", &
!        & "ATOM  CD1    C.ar   C.ar       2.000  0.000  0.000  H    0.337   0.017  2  cF", &
!        & "ATOM  HD1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CD2    C.ar   C.ar       2.000  0.000  0.000  H    0.337   0.017  2  cF", &
!        & "ATOM  HD2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CE1    C.ar   C.ar       2.000  0.000  0.000  H    0.337   0.017  2  cF", &
!        & "ATOM  HE1    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CE2    C.ar   C.ar       2.000  0.000  0.000  H    0.337   0.017  2  cF", &
!        & "ATOM  HE2    H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CZ     C.ar   C.ar.x     2.000  0.000  0.000  P   -0.151   0.017  2  cP", &
!        & "ATOM  OH     O.3    O.3.h      1.650  0.000  0.000  DA   0.082  -0.014  0  OD", &
!        & "ATOM  HH     H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "BOND  N      H      1                                                        ", &
!        & "BOND  N      CA     1                                                        ", &
!        & "BOND  C      O      2                                                        ", &
!        & "BOND  C      CA     1                                                        ", &
!        & "BOND  CA     HA     1                                                        ", &
!        & "BOND  CA     CB     1                                                        ", &
!        & "BOND  CB     HB1    1                                                        ", &
!        & "BOND  CB     HB2    1                                                        ", &
!        & "BOND  CB     CG     1                                                        ", &
!        & "BOND  CG     CD1    ar                                                       ", &
!        & "BOND  CD1    HD1    1                                                        ", &
!        & "BOND  CG     CD2    ar                                                       ", &
!        & "BOND  CD2    HD2    1                                                        ", &
!        & "BOND  CD1    CE1    ar                                                       ", &
!        & "BOND  CE1    HE1    1                                                        ", &
!        & "BOND  CD2    CE2    ar                                                       ", &
!        & "BOND  CE2    HE2    1                                                        ", &
!        & "BOND  CE1    CZ     ar                                                       ", &
!        & "BOND  CE2    CZ     ar                                                       ", &
!        & "BOND  CZ     OH     1                                                        ", &
!        & "BOND  OH     HH     1                                                        ", &
!        & "RESI  VAL    16     15    0.00 Valine                                        ", &
!        & "ATOM  N      N.am   N.pl3.h    1.750  0.000  0.000  D   -0.096  -0.052  0  ND", &
!        & "ATOM  H      H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  C      C.2    C.2.x      1.900  0.000  0.000  P   -0.030   0.023  0  CP", &
!        & "ATOM  O      O.2    O.2        1.550  0.000  0.000  A   -0.399  -0.035  0  OA", &
!        & "ATOM  CA     C.3    C.3.x      2.200  0.000  0.000  P   -0.305   0.005  0  CP", &
!        & "ATOM  HA     H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CB     C.3    C.3        2.200  0.000  0.000  H    0.127   0.021  0  CF", &
!        & "ATOM  HB     H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CG1    C.3    C.3        2.200  0.000  0.000  H    0.528   0.031  0  CF", &
!        & "ATOM  HG11   H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HG12   H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HG13   H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  CG2    C.3    C.3        2.200  0.000  0.000  H    0.528   0.031  0  CF", &
!        & "ATOM  HG21   H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HG22   H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "ATOM  HG23   H      H          1.000  0.000  0.000  N    0.000   0.000  0  HH", &
!        & "BOND  N      H      1                                                        ", &
!        & "BOND  N      CA     1                                                        ", &
!        & "BOND  C      O      2                                                        ", &
!        & "BOND  C      CA     1                                                        ", &
!        & "BOND  CA     HA     1                                                        ", &
!        & "BOND  CA     CB     1                                                        ", &
!        & "BOND  CB     HB     1                                                        ", &
!        & "BOND  CB     CG1    1                                                        ", &
!        & "BOND  CG1    HG11   1                                                        ", &
!        & "BOND  CG1    HG12   1                                                        ", &
!        & "BOND  CG1    HG13   1                                                        ", &
!        & "BOND  CB     CG2    1                                                        ", &
!        & "BOND  CG2    HG21   1                                                        ", &
!        & "BOND  CG2    HG22   1                                                        ", &
!        & "BOND  CG2    HG23   1                                                        ", &
!        & "RESI  TER     5      0    0.00 N- and C-terminal atoms                       ", &
!        & "ATOM  OXT    O.co2  O.co2      1.550  0.000 -0.500  A   -0.399   0.000  0  OC", &
!        & "ATOM  HOCA   H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  H1     H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  H2     H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  H3     H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "RESI  HOH     3      2    0.00 Water                                         ", &
!        & "ATOM  O      O.w    O.w        1.750  0.000  0.000  DA  -1.000  -0.014  0  OW", &
!        & "ATOM  H1     H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "ATOM  H2     H      H.hb       1.000  0.000  0.000  DH   0.000  -0.115  0  HH", &
!        & "BOND  O      H1     1                                                        ", &
!        & "BOND  O      H2     1                                                        ", &
!        & "RESI  SO4     5      4   -2.00 SO4--                                         ", &
!        & "ATOM  S      S.o2   S.o        2.000  0.000  0.000  P    0.000   0.000  0  UN", &
!        & "ATOM  O1     O.co2  O.co2      1.550  0.000 -0.500  A   -0.399   0.000  0  OC", &
!        & "ATOM  O2     O.co2  O.co2      1.550  0.000 -0.500  A   -0.399   0.000  0  OC", &
!        & "ATOM  O3     O.co2  O.co2      1.550  0.000 -0.500  A   -0.399   0.000  0  OC", &
!        & "ATOM  O4     O.co2  O.co2      1.550  0.000 -0.500  A   -0.399   0.000  0  OC", &
!        & "BOND  S      O1     2                                                        ", &
!        & "BOND  S      O2     2                                                        ", &
!        & "BOND  S      O3     2                                                        ", &
!        & "BOND  S      O4     2                                                        ", &
!        & "RESI  PO4     5      4   -2.00 PO4--                                         ", &
!        & "ATOM  P      P.3    P.3        2.000  0.000  0.000  P    0.000   0.000  0  UN", &
!        & "ATOM  O1     O.co2  O.co2      1.550  0.000 -0.500  A   -0.399   0.000  0  OC", &
!        & "ATOM  O2     O.co2  O.co2      1.550  0.000 -0.500  A   -0.399   0.000  0  OC", &
!        & "ATOM  O3     O.co2  O.co2      1.550  0.000 -0.500  A   -0.399   0.000  0  OC", &
!        & "ATOM  O4     O.co2  O.co2      1.550  0.000 -0.500  A   -0.399   0.000  0  OC", &
!        & "BOND  P      O1     2                                                        ", &
!        & "BOND  P      O2     2                                                        ", &
!        & "BOND  P      O3     2                                                        ", &
!        & "BOND  P      O4     2                                                        ", &
!        & "RESI  HET    19      0    0.00 ions                                          ", &
!        & "ATOM  LI     Li     M+         1.250  0.000  1.000  M   -2.000   0.000  0  UN", &
!        & "ATOM  NA     Na     M+         1.250  0.000  1.000  M   -2.000   0.000  0  UN", &
!        & "ATOM  K      K      M+         1.250  0.000  1.000  M   -2.000   0.000  0  UN", &
!        & "ATOM  CA     Ca     M+         1.250  0.000  2.000  M   -2.000   0.000  0  UN", &
!        & "ATOM  MG     Mg     M+         1.250  0.000  2.000  M   -2.000   0.000  0  UN", &
!        & "ATOM  AL     Al     M+         1.250  0.000  3.000  M   -2.000   0.000  0  UN", &
!        & "ATOM  MN     Mn     M+         1.250  0.000  2.000  M   -2.000   0.000  0  UN", &
!        & "ATOM  FE     Fe     M+         1.250  0.000  2.000  M   -2.000   0.000  0  UN", &
!        & "ATOM  NI     Ni     M+         1.250  0.000  2.000  M   -2.000   0.000  0  UN", &
!        & "ATOM  CD     Cd     M+         1.250  0.000  2.000  M   -2.000   0.000  0  UN", &
!        & "ATOM  CO     Co     M+         1.250  0.000  2.000  M   -2.000   0.000  0  UN", &
!        & "ATOM  CU     Cu     M+         1.250  0.000  2.000  M   -2.000   0.000  0  UN", &
!        & "ATOM  ZN     Zn     M+         1.250  0.000  2.000  M   -2.000   0.000  0  UN", &
!        & "ATOM  HG     Hg     M+         1.250  0.000  2.000  M   -2.000   0.000  0  UN", &
!        & "ATOM  U      U      M+         1.250  0.000  3.000  M   -2.000   0.000  0  UN", &
!        & "ATOM  F      F      F-         1.500  0.000 -1.000  P   -1.000   0.000  0  UN", &
!        & "ATOM  CL     Cl     Cl-        1.750  0.000 -1.000  P   -1.000   0.000  0  UN", &
!        & "ATOM  BR     Br     Br-        1.900  0.000 -1.000  P   -1.000   0.000  0  UN", &
!        & "ATOM  I      I      I-         2.050  0.000 -1.000  P   -1.000   0.000  0  UN" /)


end module xlogp_ResidueDefs_m
