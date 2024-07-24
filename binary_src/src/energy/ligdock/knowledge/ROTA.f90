!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/ligdock/knowledge/ROTA.f90
!
! Description: ROTA energy term for protein-protein interaction in
!              flexible docking
!
!-------------------------------------------------------------------------------
MODULE ROTA
!-------------------------------------------------------------------------------

use globals
use logger, only: log_p, terminate_with_error
use convert_res_name, only: convert_to_stdres
use rotamer, only: get_rotamer_index, place_rotamer, &
                   identify_rotamer_state_single
!
use energy_vars
use energy_utils, only: protein_to_R, update_R_for_sc

implicit none
save
private

!===============================================================================
! PARAMETERS
!===============================================================================
integer, parameter :: max_std_aa = 20
integer, parameter :: max_ROTA_bin = 40
integer, parameter :: ALA = 1
integer, parameter :: ARG = 2
integer, parameter :: ASN = 3
integer, parameter :: ASP = 4
integer, parameter :: CYS = 5
integer, parameter :: GLN = 6
integer, parameter :: GLU = 7
integer, parameter :: GLY = 8
integer, parameter :: HIS = 9
integer, parameter :: ILE = 10
integer, parameter :: LEU = 11
integer, parameter :: LYS = 12
integer, parameter :: MET = 13
integer, parameter :: PHE = 14
integer, parameter :: PRO = 15
integer, parameter :: SER = 16
integer, parameter :: THR = 17
integer, parameter :: TRP = 18
integer, parameter :: TYR = 19
integer, parameter :: VAL = 20

!===============================================================================
! DERIVED TYPES
!===============================================================================
type rota_partner_type
!-------------------------------------------------------------------------------
character(len=4) :: atom_name
real(dp) :: bin_score(max_ROTA_bin)
!-------------------------------------------------------------------------------
end type rota_partner_type
!-------------------------------------------------------------------------------
type rota_atom_record_type
!-------------------------------------------------------------------------------
character(len=4) :: atom_name
character(len=3) :: p_res_name(max_std_aa)
type(rota_partner_type) :: partner(max_atm, max_std_aa)
!-------------------------------------------------------------------------------
end type rota_atom_record_type
!-------------------------------------------------------------------------------
type rota_record_type
!-------------------------------------------------------------------------------
character(len=3) :: res_name
type(rota_atom_record_type) :: atom(max_atm)
!-------------------------------------------------------------------------------
end type rota_record_type
!-------------------------------------------------------------------------------

type(rota_record_type) :: rota_record(max_std_aa)
integer, allocatable :: rota_idx(:,:) ! rota_idx(1,i_atm) = rota_res_idx
                                      ! rota_idx(2,i_atm) = rota_atm_idx

public :: initialize_rota
public :: finalize_rota
public :: calc_ROTA_score

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_rota(protein)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein

! 1. read ROTA data file
call initialize_rota_record(rota_record)
call read_rota_data(rota_record)
call log_p('Reading ROTA score data file.... Done!')

! 2. assign index to access rota_record_type
allocate(rota_idx(2, tn%stdatm))
call assign_ROTA_index(protein, rota_idx)

end subroutine initialize_rota
!-------------------------------------------------------------------------------
subroutine finalize_rota()
!-------------------------------------------------------------------------------

deallocate(rota_idx)

end subroutine finalize_rota
!-------------------------------------------------------------------------------
subroutine initialize_rota_record(rota_record)
!-------------------------------------------------------------------------------
type(rota_record_type), intent(inout) :: rota_record(:)
integer :: i_res, i_atm, j_res, j_atm

do i_res = 1, max_std_aa
    rota_record(i_res)%res_name = ''
    do i_atm = 1, max_atm
        rota_record(i_res)%atom(i_atm)%atom_name = ''
        rota_record(i_res)%atom(i_atm)%p_res_name(:) = ''
        do j_res = 1, max_std_aa
            do j_atm = 1, max_atm
                rota_record(i_res)%atom(i_atm)%partner(j_atm,j_res)%atom_name = ''
                rota_record(i_res)%atom(i_atm)%partner(j_atm,j_res)%bin_score(:) = 0.0
            end do
        end do
    end do
end do

end subroutine initialize_rota_record
!-------------------------------------------------------------------------------
subroutine read_rota_data(rota_record)
!-------------------------------------------------------------------------------
type(rota_record_type), intent(inout) :: rota_record(:)
!
character(len=len_fname) :: infile_ROTA
integer :: f_unit, ioerr
character(len=300) :: line
!
character(len=3) :: res_name_prev, res_name, p_res, p_res_prev
character(len=4) :: atm_name_prev, atm_name, p_atm
integer :: i_res, j_res, i_atm, j_atm

infile_ROTA = trim(data_dir) // 'rota.1_2.txt'

f_unit = 17
open(unit=f_unit, file=trim(infile_ROTA), action='read', iostat=ioerr)
if (ioerr /= 0) then
    write(log_msg, '(A,A)') 'ERROR: ROTA data file does not exist!', trim(infile_ROTA)
    call terminate_with_error(log_msg)
end if

res_name_prev = ''
atm_name_prev = ''
p_res_prev = ''

i_res = 0
i_atm = 0
j_res = 0
j_atm = 0

do 
    read(f_unit, '(a300)', iostat=ioerr) line
    if(ioerr /= 0) exit
 
    read(line(1:3), *) res_name
    if(res_name == 'HET' .or. res_name == 'REC') cycle
    if(res_name /= res_name_prev) then
        call get_ROTA_res_idx(res_name, i_res)
        rota_record(i_res)%res_name = res_name
        i_atm = 0
    end if
 
    read(line(5:8), *) atm_name
    if(atm_name(1:1) == '1' .or. atm_name(1:1) == '2' .or. atm_name(1:1) == '3') then
        atm_name = trim(atm_name(2:4)) // atm_name(1:1)
    end if
    if(atm_name /= atm_name_prev) then
        i_atm = i_atm + 1
        rota_record(i_res)%atom(i_atm)%atom_name = atm_name
        j_res = 0
    end if

    read(line(13:15), *) p_res
    if(p_res == 'HET' .or. p_res == 'REC') cycle
    if(p_res /= p_res_prev) then
        call get_ROTA_res_idx(p_res, j_res)
        rota_record(i_res)%atom(i_atm)%p_res_name(j_res) = p_res
        j_atm = 0
    end if
 
    read(line(17:20), *) p_atm
    j_atm = j_atm + 1
    if(p_atm(1:1) == '1' .or. p_atm(1:1) == '2' .or. p_atm(1:1) == '3') then
        p_atm = trim(p_atm(2:4)) // p_atm(1:1)
    end if
    rota_record(i_res)%atom(i_atm)%partner(j_atm, j_res)%atom_name = p_atm
 
    read(line(25:), *) &
        rota_record(i_res)%atom(i_atm)%partner(j_atm, j_res)%bin_score(:)

    res_name_prev = res_name
    atm_name_prev = atm_name
    p_res_prev = p_res
end do

! Debug
!do i_res = 1, max_std_aa
!    write(*,*) 'residue: ',i_res, rota_record(i_res)%res_name
!    do i_atm = 1, max_atm
!        if (rota_record(i_res)%atom(i_atm)%atom_name == '') cycle
!        write(*,*) 'atom: ',i_res, i_atm, rota_record(i_res)%atom(i_atm)%atom_name
!        do j_res = 1, max_std_aa
!            write(*,*) 'partner: ', rota_record(i_res)%atom(i_atm)%p_res_name(j_res)
!            do j_atm = 1, max_atm
!                if (rota_record(i_res)%atom(i_atm)%partner(j_atm,j_res)%atom_name == '') cycle
!                write(*,*) 'p_atm: ', rota_record(i_res)%atom(i_atm)%partner(j_atm,j_res)%atom_name
!                write(*,*) 'score: ', rota_record(i_res)%atom(i_atm)%partner(j_atm,j_res)%bin_score(:)
!            end do
!        end do
!    end do
!end do
!stop

end subroutine read_rota_data
!-------------------------------------------------------------------------------
subroutine assign_ROTA_index(protein, rota_idx)
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
integer, intent(inout) :: rota_idx(:,:)
character(len=4) :: res_name, atom_name
integer :: i_atm, j_atm, res_no, atm_no
integer :: res_idx, atm_idx

rota_idx(:,:) = 0

do res_no = 1, protein%n_res
    do atm_no = 1, protein%residue(res_no)%n_atm
        i_atm = ii_R(atm_no, res_no)

        res_name = protein%residue(res_no)%res_name
        atom_name = ref_res(protein%residue(res_no)%res_type)%atom_name(atm_no)
        call convert_to_stdres(res_name)
        call get_ROTA_res_idx(res_name(1:3), res_idx)
        atm_idx = -1
        do j_atm = 1, max_atm
            if (atom_name == rota_record(res_idx)%atom(j_atm)%atom_name) then
                atm_idx = j_atm
                exit
            end if
        end do
        rota_idx(1, i_atm) = res_idx
        rota_idx(2, i_atm) = atm_idx
    end do
end do

end subroutine assign_ROTA_index
!-------------------------------------------------------------------------------
subroutine get_ROTA_res_idx(res_name, res_idx)
!-------------------------------------------------------------------------------
character(len=3), intent(in) :: res_name
integer, intent(out) :: res_idx

res_idx = -1

if (res_name == 'ALA') then
    res_idx = ALA
else if (res_name == 'ARG') then
    res_idx = ARG
else if (res_name == 'ASN') then
    res_idx = ASN
else if (res_name == 'ASP') then
    res_idx = ASP
else if (res_name == 'CYS') then
    res_idx = CYS
else if (res_name == 'GLN') then
    res_idx = GLN
else if (res_name == 'GLU') then
    res_idx = GLU
else if (res_name == 'GLY') then
    res_idx = GLY
else if (res_name == 'HIS') then
    res_idx = HIS
else if (res_name == 'ILE') then
    res_idx = ILE
else if (res_name == 'LEU') then
    res_idx = LEU
else if (res_name == 'LYS') then
    res_idx = LYS
else if (res_name == 'MET') then
    res_idx = MET
else if (res_name == 'PHE') then
    res_idx = PHE
else if (res_name == 'PRO') then
    res_idx = PRO
else if (res_name == 'SER') then
    res_idx = SER
else if (res_name == 'THR') then
    res_idx = THR
else if (res_name == 'TRP') then
    res_idx = TRP
else if (res_name == 'TYR') then
    res_idx = TYR
else if (res_name == 'VAL') then
    res_idx = VAL
end if

end subroutine get_ROTA_res_idx
!-------------------------------------------------------------------------------
subroutine calc_ROTA_score(f, g, appl_pair, calc_g)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: f, g(:,:)
logical, intent(in) :: appl_pair(:,:), calc_g
!
integer :: i_atm, j_atm, ia, n_pair, i_pair(max_neigh)
integer :: res_1, atm_1, res_2, atm_2
integer :: i_bin, j_bin
real(dp) :: dr(3), dist, tmp_dist, tmp_dist_2, score(2)

f = 0.0d0
g(:,:) = 0.0d0

do i_atm = 1, tn%stdatm
    res_1 = ROTA_idx(1, i_atm)
    atm_1 = ROTA_idx(2, i_atm)
    if (res_1 == -1) cycle
    if (atm_1 == -1) cycle
    call get_ROTA_pair(i_atm, appl_pair, n_pair, i_pair)

    do ia = 1, n_pair
        j_atm = i_pair(ia)
        res_2 = ROTA_idx(1, j_atm)
        atm_2 = ROTA_idx(2, j_atm)

        dr(1:3) = R(1:3, j_atm) - R(1:3, i_atm)
        dist = sqrt(dot_product(dr, dr))
        
        tmp_dist = dist/0.25d0 + 0.5d0

        i_bin = int(tmp_dist)
        
        if (i_bin < 1) i_bin = 1
        if (i_bin >= 40) then
            i_bin = 40
            j_bin = 40
        else
            j_bin = i_bin + 1
        end if

        tmp_dist = tmp_dist - dble(i_bin)
        tmp_dist_2 = 1.0d0 - tmp_dist

        score(1) = rota_record(res_1)%atom(atm_1)%partner(atm_2, res_2)%bin_score(i_bin)
        score(2) = rota_record(res_1)%atom(atm_1)%partner(atm_2, res_2)%bin_score(j_bin)

        f =  f + tmp_dist_2 * score(1) + tmp_dist * score(2)
    end do
end do

end subroutine calc_ROTA_score
!-------------------------------------------------------------------------------
subroutine get_ROTA_pair(i_atm, appl_pair, n_pair, i_pair)
!-------------------------------------------------------------------------------
integer, intent(in) :: i_atm
logical, intent(in) :: appl_pair(:,:)
integer, intent(out) :: n_pair
integer, intent(out) :: i_pair(:)

integer :: ia_start, ia
integer :: i_res, j_res, j_atm, atm_1, atm_2
logical :: is_sc_1, is_sc_2

n_pair = 0
i_res = i_R(1, i_atm)
atm_1 = i_R(2, i_atm)
is_sc_1 = ref_res(res_index(i_res)%ref_res_no)%is_sc_atom(atm_1)

ia_start = i_P(i_atm)%pair_end_index(2)+1

do ia = ia_start, i_P(i_atm)%n_pair
    j_atm = i_P(i_atm)%i_pair(ia)
    j_res = i_R(1, j_atm)
    atm_2 = i_R(2, j_atm)
    if (j_atm > tn%stdatm) cycle
    if (.not. appl_pair(j_res, i_res)) cycle
    is_sc_2 = ref_res(res_index(j_res)%ref_res_no)%is_sc_atom(atm_2)
    if ((.not. is_sc_1) .and. (.not. is_sc_2)) cycle
    if (ROTA_idx(1,j_atm) == -1) cycle
    if (ROTA_idx(2,j_atm) == -1) cycle
    n_pair = n_pair + 1
    i_pair(n_pair) = j_atm
end do

end subroutine get_ROTA_pair
!-------------------------------------------------------------------------------
END MODULE ROTA
!-------------------------------------------------------------------------------
