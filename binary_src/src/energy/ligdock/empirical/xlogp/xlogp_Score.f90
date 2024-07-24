!-------------------------------------------------------------------------------
module xlogp_Score_m
!-------------------------------------------------------------------------------

use globals
!
use xlogp_Molecule_m
use xlogp_Mol2Reader_m
use xlogp_AromaticRing_m
use xlogp_Parameter_m
use xlogp_Logp_m
use xlogp_AtomTypeChecker_m
use xlogp_RotorCount_m

implicit none

type(Molecule_t), target :: mol
type(Molecule_t), pointer :: molPtr !=> null()
type(LogpEvaluate_t) :: logpEval
type(ForceField_t), target :: forceField
type(ForceField_t), pointer :: forceFieldPtr !=> null()

character(len=256) :: g_ParameterPathForLogp
logical :: g_ParameterSetup_flag = .false.

CONTAINS
!-------------------------------------------------------------------------------
subroutine get_logp_of_ligand(ligand, logpValues, xlogpCorrection)
!-------------------------------------------------------------------------------
!> This function will calculate the logp values of atoms for a given ligand
! , and it will save the logp values into argument logpValues.
! @param ligand
! @param logpValues array of real types to store logp values.
! @param xlogpCorrection. This flag determines if correction factors are distributed
!      and written onto each atom.
!-------------------------------------------------------------------------------
type(ref_lig_type), intent(inout) :: ligand
real(dp), intent(out) :: logpValues(ligand%n_atm)
logical, intent(in) :: xlogpCorrection
logical :: Success_Fail
logical :: DebugLogpCalc
logical :: paramSetupMark
integer :: i

!> ForceField (Ligand Parameter for xlogp calc) Init
call SetupParameterForLogp(data_dir, paramSetupMark)

forceFieldPtr => forceField
g_ParameterSetup_flag = .true.
DebugLogpCalc = .false.

call CalcLogpOfLigand(ligand, xlogpCorrection, Success_Fail, DebugLogpCalc)
do i=1, ligand%n_atm
    logpValues(i) = mol%atoms(i)%logp
end do

end subroutine get_logp_of_ligand
!-------------------------------------------------------------------------------
subroutine get_rotor_of_ligand(ligand, rotorSum)
!-------------------------------------------------------------------------------
type(ref_lig_type), intent(inout) :: ligand
real(dp), intent(out) :: rotorSum
logical :: Success_Fail
logical :: paramSetupMark

!> ForceField (Ligand Parameter for xlogp calc) Init
call SetupParameterForLogp(data_dir, paramSetupMark)

forceFieldPtr => forceField
g_ParameterSetup_flag = .true.

call Mol_Clear(mol)
call M_make_mol_from_ligand(mol, ligand)
molPtr => mol

Success_Fail = .true.

!> Make connection information
call M_DetectConnections(mol)
call CheckAllAtomTypes(mol%atoms, mol%num_atom, &
    mol%bonds, Success_Fail)

!> Detect rings and aromatic rings
call DetectRings(mol%atoms, mol%num_atom, &
    mol%bonds, mol%num_bond, mol%rings, mol%num_ring)
call DetectAromaticRings(mol%atoms, &
    mol%bonds, mol%rings, mol%num_ring)
mol%num_ring = M_GetRingSize(mol)

!> LogpEvaluate_t init
call Eval_InitMolecule(logpEval, molPtr)
call Eval_InitForceField(logpEval, forceFieldPtr)
!call logpEval%InitCoeff()

!> Classify atom type : xtool type, xlogp type
call Eval_FindAllXToolType(logpEval)
call Eval_FindAllXLogPType(logpEval)

call CheckAllAtomTypes(mol%atoms, mol%num_atom, &
    mol%bonds, Success_Fail)

!! Find atom parameters type, xlogptype, hbond type, vdwType
call Param_AssignAll_Atoms(forceFieldPtr, mol%atoms, mol%num_atom)

!> Calculate Rotor
call CountRotor(mol%atoms, mol%num_atom, &
    mol%bonds, mol%num_bond, rotorSum)

end subroutine get_rotor_of_ligand
!-------------------------------------------------------------------------------
subroutine get_xlogP_atom_types(ligand, atom_types)
!-------------------------------------------------------------------------------
type(ref_lig_type), intent(inout) :: ligand
character(len=20), intent(out) :: atom_types(:)
logical :: Success_Fail
logical :: paramSetupMark
integer :: i

!> ForceField (Ligand Parameter for xlogp calc) Init
call SetupParameterForLogp(data_dir, paramSetupMark)

forceFieldPtr => forceField
g_ParameterSetup_flag = .true.

call Mol_Clear(mol)
call M_make_mol_from_ligand(mol, ligand)
molPtr => mol

Success_Fail = .true.

!> Make connection information
call M_DetectConnections(mol)
call CheckAllAtomTypes(mol%atoms, mol%num_atom, &
    mol%bonds, Success_Fail)

!> Detect rings and aromatic rings
call DetectRings(mol%atoms, mol%num_atom, &
    mol%bonds, mol%num_bond, mol%rings, mol%num_ring)
call DetectAromaticRings(mol%atoms, &
    mol%bonds, mol%rings, mol%num_ring)
mol%num_ring = M_GetRingSize(mol)

!> LogpEvaluate_t init
call Eval_InitMolecule(logpEval, molPtr)
call Eval_InitForceField(logpEval, forceFieldPtr)
!call logpEval%InitCoeff()

!> Classify atom type : xtool type, xlogp type
call Eval_FindAllXToolType(logpEval)
call Eval_FindAllXLogPType(logpEval)

call CheckAllAtomTypes(mol%atoms, mol%num_atom, &
    mol%bonds, Success_Fail)

!! Find atom parameters type, xlogptype, hbond type, vdwType
call Param_AssignAll_Atoms(forceFieldPtr, mol%atoms, mol%num_atom)

do i = 1, ligand%n_atm
    atom_types(i) = mol%atoms(i)%xlogptype
end do

end subroutine get_xlogp_atom_types
!-------------------------------------------------------------------------------
subroutine CalcLogpOfLigand(ligand, xlogpCorrection, Success_Fail, DebugLogpCalc)
!-------------------------------------------------------------------------------
!> This will calculate logp values of the given ligand.
!! @param ligand
!! @param xlogpCorrection. This flag determines if correction factors are distributed
!!      and written onto each atom.
!! @param Success_Fail will say if the overall subroutine has succeeded or not.
!! @param DebugLogpCalc will make verbose output
!-------------------------------------------------------------------------------
type(ref_lig_type), intent(inout) :: ligand
logical, intent(in) :: xlogpCorrection
logical, intent(out) :: Success_Fail
logical, intent(in), optional :: DebugLogpCalc
character(len=256) :: formula

Success_Fail = .false.
if(.not. g_ParameterSetup_flag) then
    print*, '(*) You need to set up parameter !!'
    return
else if(.not. associated(forceFieldPtr)) then
    print*, '(*) You need to set up parameter !!'
    return
endif

call Mol_Clear(mol)
call M_make_mol_from_ligand(mol, ligand)
molPtr => mol

Success_Fail = .true.

!> Make connection information
call M_DetectConnections(mol)
call CheckAllAtomTypes(mol%atoms, mol%num_atom, &
    mol%bonds, Success_Fail)

!> Detect rings and aromatic rings
call DetectRings(mol%atoms, mol%num_atom, &
    mol%bonds, mol%num_bond, mol%rings, mol%num_ring)
call DetectAromaticRings(mol%atoms, &
    mol%bonds, mol%rings, mol%num_ring)

!> LogpEvaluate_t init
call Eval_InitMolecule(logpEval, molPtr)
call Eval_InitForceField(logpEval, forceFieldPtr)
call Eval_InitCoeff(logpEval)

!> Classify atom type : xtool type, xlogp type
call Eval_FindAllXToolType(logpEval)
call Eval_FindAllXLogPType(logpEval)

call CheckAllAtomTypes(mol%atoms, mol%num_atom, &
    mol%bonds, Success_Fail)

!! Find atom parameters type, xlogptype, hbond type, vdwType
call Param_AssignAll_Atoms(forceFieldPtr, mol%atoms, mol%num_atom)

mol%weight = M_GetWeight(mol)
call M_GetFormula(mol, formula)
mol%formula = TRIM(formula)

!> Calculate Logp
if(present(DebugLogpCalc)) then
    call Eval_SetDebug(logpEval, DebugLogpCalc)
else
    call Eval_SetDebug(logpEval, .false.)
endif
call Eval_CalcLogP(logpEval, xlogpCorrection)

end subroutine CalcLogpOfLigand
!-------------------------------------------------------------------------------
subroutine SetupParameterForLogp(paramPath, mark)
!-------------------------------------------------------------------------------
character(len=*), intent(in) :: paramPath
logical, intent(out) :: mark
character(len=120) :: fileNameResidueParam
character(len=120) :: fileNameAtomParam
character(len=120) :: fileNameXlogpParam
character(len=120), parameter :: xscore_residue_param_input = &
    'RESIDUE_DEF_XTOOL'
character(len=120), parameter :: xscore_atom_param_input = &
    "ATOM_DEF_XTOOL"
character(len=120), parameter :: xscore_xlogp_param_input = &
    "ATOM_DEF_XLOGP"
integer :: nLen

g_ParameterPathForLogp = TRIM(paramPath)

g_ParameterSetup_flag = .false.
mark = g_ParameterSetup_flag
nLen = len_trim(g_ParameterPathForLogp)

if(nLen < 1) then
    print*, '(*) Parameter path is not correct: ', paramPath
    return
endif

if(g_ParameterPathForLogp(nLen:nLen) /= '/') then
    g_ParameterPathForLogp(nLen+1:nLen+1) = '/'
    nLen = nLen + 1
endif

write(fileNameResidueParam, '(A,A)') g_ParameterPathForLogp(1:nLen), TRIM(xscore_residue_param_input)
write(fileNameAtomParam, '(A,A)') g_ParameterPathForLogp(1:nLen), TRIM(xscore_atom_param_input)
write(fileNameXlogpParam, '(A,A)') g_ParameterPathForLogp(1:nLen), TRIM(xscore_xlogp_param_input)

!> Check if param files exist or not
inquire(FILE=fileNameResidueParam, exist=mark)
if(.not. mark) then
    print*, '(#) fileNameResidueParam does not exist: ', fileNameResidueParam
    return
endif
inquire(FILE=fileNameAtomParam, exist=mark)
if(.not. mark) then
    print*, '(#) fileNameAtomParam does not exist: ', fileNameAtomParam
    return
endif
inquire(FILE=fileNameXlogpParam, exist=mark)
if(.not. mark) then
    print*, '(#) fileNameXlogpParam does not exist: ', fileNameXlogpParam
    return
endif

!> ForceField (Ligand Parameter for xlogp calc) Init
call Param_SetUp(forceField, fileNameResidueParam, &
    fileNameAtomParam, &
    fileNameXlogpParam)
forceFieldPtr => forceField

g_ParameterSetup_flag = .true.
mark = g_ParameterSetup_flag

end subroutine SetupParameterForLogp
!-------------------------------------------------------------------------------
end module xlogp_Score_m
