!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/ligdock/knowledge/drugscore.f90
!
! Description: DrugScore energy term
!
!-------------------------------------------------------------------------------
MODULE DRUGSCORE
!-------------------------------------------------------------------------------

use globals
use logger
use string, only: parse_string
use sort, only: sort1
!
use energy_vars 
!
use ligdock_E_vars
use ligdock_E_utils, only: find_mol2_res_idx, find_mol2_atm_idx, &
                           check_is_new_type, copy_to_grid, calc_E_using_grid

implicit none
save 
private

!===============================================================================
! Parameters
!===============================================================================
! distance cutoff of energy calculation
real(dp), parameter :: DScutoff = 7.0d0
real(dp), parameter :: DScutoff2 = DScutoff**2
! DrugScore types related
integer, parameter :: max_DS_type = 18
integer :: Br, C2, C3, Car, Ccat, Cl, Fluo, Iod, Met
integer :: N3, Nam, Nar, Npl3, O2, O3, Oco2, P3, S3
!
integer, allocatable :: DS_para_idx(:)
logical, allocatable :: is_ligand(:)
integer :: target_lig_no            ! target_lig_no = tn%nonlig+ligand%lig_no

!===============================================================================
! Drived types
!===============================================================================
type DS_record_type
!-------------------------------------------------------------------------------
! Energy table of DrugScore
! If distance index = i, then energy = bin_score(i)
!-------------------------------------------------------------------------------
real(dp) :: bin_score(50)
!-------------------------------------------------------------------------------
end type DS_record_type
!-------------------------------------------------------------------------------
type DS_eng_para_type
!-------------------------------------------------------------------------------
integer :: n_atm_types(3)
integer :: type_idx(max_DS_type, 3)
!-------------------------------------------------------------------------------
end type DS_eng_para_type
type(DS_record_type) :: DS_record(max_DS_type, max_DS_type)
type(DS_eng_para_type) :: DS_para

!===============================================================================
! public
!===============================================================================
public :: initialize_DrugScore
public :: finalize_DrugScore
public :: construct_drugscore_grid
public :: calc_DS_intrxn_E_using_grid
public :: calc_DrugScore
public :: FFT_rec_grid_drugscore
public :: FFT_lig_grid_drugscore

CONTAINS
!===============================================================================
! Initialize/Finalize DrugScore
!===============================================================================
subroutine initialize_DrugScore(protein, ligand, tna, flex_res, &
                                residue_mol2_info, num_mol2_info)
!-------------------------------------------------------------------------------
! setup DrugScore parameters 
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
type(ligand_type), intent(in) :: ligand
integer, intent(in) :: tna
logical, intent(in) :: flex_res(:)
type(residue_mol2_type), intent(in) :: residue_mol2_info(:)
integer, intent(in) :: num_mol2_info
integer :: DS_types(tna)

! read parameter file
call read_DS_param_file(infile_drugscore_prm, DS_record)

! Do atom typing
call assign_DS_atom_types(protein, DS_types, residue_mol2_info, &
                          num_mol2_info)

! store needed parameters to DS_eng_para type
allocate(DS_para_idx(tn%atom))
allocate(is_ligand(tn%atom))
target_lig_no = tn%nonlig + ligand%lig_no
call construct_DS_para(DS_types, DS_para_idx, flex_res)

end subroutine initialize_DrugScore
!-------------------------------------------------------------------------------
subroutine finalize_DrugScore()
!-------------------------------------------------------------------------------
deallocate(DS_para_idx)
deallocate(is_ligand)

end subroutine finalize_DrugScore
!-------------------------------------------------------------------------------
subroutine read_DS_param_file(file_name, DS_record)
!-------------------------------------------------------------------------------
! Read DrugScore data file
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: file_name
type(DS_record_type) :: DS_record(:,:)
integer :: f_unit, ioerror
!
integer :: num_word
character(len=len_fname) :: line, word(25)
!
character(len=len_fname) :: p_atm_type, p_atm_type_prev, l_atm_type, l_atm_type_prev
integer :: i_atm, j_atm, i_bin
real(dp) :: rij, score
character(len=6) :: index_list(max_DS_type)

f_unit = 18
open(unit=f_unit, file=trim(file_name), action='read', iostat=ioerror)
if(ioerror /= 0) then
   write(log_msg, '(A,A)') 'DrugScore_pair data file does not exist! ', &
                             trim(file_name)
   call terminate_with_error(log_msg)
end if
  
l_atm_type_prev = ''
p_atm_type_prev = ''

i_atm = 0
j_atm = 0

do 
   read(f_unit, '(a120)', iostat=ioerror) line
   if(ioerror /= 0) exit
   if (line(1:1) == "!" .or. line(1:1) == "#") cycle
   
   call parse_string(line, num_word, word)
  
   p_atm_type = word(1)
   l_atm_type = word(2)
  
   if (l_atm_type /= l_atm_type_prev) then
      j_atm = j_atm + 1
      if (p_atm_type /= p_atm_type_prev) then
         i_atm = i_atm + 1
         j_atm = 1
         index_list(i_atm) = p_atm_type
      end if
   end if

   read(word(3), '(F4.2)') rij
   read(word(4), '(F8.4)') score

   i_bin = int(10*(rij-1.00) + 1)
   DS_record(i_atm,j_atm)%bin_score(i_bin) = score

   p_atm_type_prev = p_atm_type
   l_atm_type_prev = l_atm_type
end do
close(f_unit)

call set_DS_index(index_list)

end subroutine read_DS_param_file
!-------------------------------------------------------------------------------
subroutine set_DS_index(index_list)
!-------------------------------------------------------------------------------
! set index for atom type in DrugScore
!-------------------------------------------------------------------------------
character(len=6), intent(in) :: index_list(:)
integer :: i_type

do i_type = 1, max_ds_type
    if (trim(index_list(i_type)) == 'Br') then
        Br = i_type
    else if (trim(index_list(i_type)) == 'C.2') then 
        C2 = i_type
    else if (trim(index_list(i_type)) == 'C.3') then 
        C3 = i_type
    else if (trim(index_list(i_type)) == 'C.ar') then 
        Car = i_type
    else if (trim(index_list(i_type)) == 'C.cat') then 
        Ccat = i_type
    else if (trim(index_list(i_type)) == 'Cl') then 
        Cl = i_type
    else if (trim(index_list(i_type)) == 'F') then 
        Fluo = i_type
    else if (trim(index_list(i_type)) == 'I') then 
        Iod = i_type
    else if (trim(index_list(i_type)) == 'Met') then 
        Met = i_type
    else if (trim(index_list(i_type)) == 'N.3') then 
        N3 = i_type
    else if (trim(index_list(i_type)) == 'N.am') then 
        Nam = i_type
    else if (trim(index_list(i_type)) == 'N.ar') then 
        Nar = i_type
    else if (trim(index_list(i_type)) == 'N.pl3') then 
        Npl3 = i_type
    else if (trim(index_list(i_type)) == 'O.2') then 
        O2 = i_type
    else if (trim(index_list(i_type)) == 'O.3') then 
        O3 = i_type
    else if (trim(index_list(i_type)) == 'O.co2') then 
        Oco2 = i_type
    else if (trim(index_list(i_type)) == 'P.3') then 
        P3 = i_type
    else if (trim(index_list(i_type)) == 'S.3') then 
        S3 = i_type
    end if
end do
!-------------------------------------------------------------------------------
end subroutine set_DS_index
!-------------------------------------------------------------------------------
subroutine assign_DS_atom_types(protein, DS_types, &
                                residue_mol2_info, num_mol2_info)
!-------------------------------------------------------------------------------
! Do atom typing of DrugScore based on mol2 type of atom
!-------------------------------------------------------------------------------
type(molecule_type), intent(in) :: protein
integer, intent(inout) :: DS_types(:)
!
type(residue_mol2_type), intent(in) :: residue_mol2_info(:)
integer, intent(in) :: num_mol2_info
!
integer :: i_atm, res_no, atm_no, ref_res_no
integer :: res_type, atm_type, mol2_idx
logical :: is_ter
character(len=4) :: res_name, atom_name
character(len=6) :: mol2_type

DS_types(:) = 0
! for protein
do i_atm = 1, tn%stdatm
    res_no = i_R(1, i_atm)
    atm_no = i_R(2, i_atm)
    ref_res_no = protein%residue(res_no)%res_type
    
    is_ter = .false.
    if (protein%residue(res_no)%ter_type == 'N' .or. &
        protein%residue(res_no)%ter_type == 'C') then
        is_ter = .true.
    end if

    if (is_ter) then
        res_name = protein%residue(res_no)%res_name(2:4)
    else
        res_name = protein%residue(res_no)%res_name
    end if

    call find_mol2_res_idx(res_name, res_type, residue_mol2_info, &
                           num_mol2_info)

    atom_name = ref_res(ref_res_no)%atom_name(atm_no)
    call find_mol2_atm_idx(atom_name, res_type, atm_type, &
                           residue_mol2_info)

    mol2_type = residue_mol2_info(res_type)%mol2_type(atm_type)
    call assign_DS_index(mol2_type, mol2_idx)
    DS_types(i_atm) = mol2_idx
end do

! for hetmol(small cofacor)
do i_atm = tn%stdatm+1, tn%nonligatm
    res_no = i_R(1, i_atm) - protein%n_res
    atm_no = i_R(2, i_atm)
    ref_res_no = protein%hetmol(res_no)%res_type

    res_name = protein%hetmol(res_no)%res_name

    call find_mol2_res_idx(res_name, res_type, residue_mol2_info, &
                           num_mol2_info)

    atom_name = ref_res(ref_res_no)%atom_name(atm_no)
    call find_mol2_atm_idx(atom_name, res_type, atm_type, &
                           residue_mol2_info)

    mol2_type = residue_mol2_info(res_type)%mol2_type(atm_type)
    call assign_DS_index(mol2_type, mol2_idx)
    DS_types(i_atm) = mol2_idx
end do

! for ligand/large cofactor
do i_atm = 1, tn%ligatm
    res_no = i_L(1, i_atm)
    atm_no = i_L(2, i_atm)
    ref_res_no = protein%ligand(res_no)%lig_type

    mol2_type = ref_lig(ref_res_no)%mol2_type(atm_no)
    call assign_DS_index(mol2_type, mol2_idx)
    DS_types(i_atm+tn%nonligatm) = mol2_idx
end do

end subroutine assign_DS_atom_types
!-------------------------------------------------------------------------------
subroutine assign_DS_index(mol2_type, mol2_idx)
!-------------------------------------------------------------------------------
character(len=6), intent(in) :: mol2_type
integer, intent(inout) :: mol2_idx

mol2_idx = -1
if (trim(mol2_type) == 'C.3') then 
    mol2_idx = C3
else if (trim(mol2_type) == 'C.ar') then 
    mol2_idx = Car
else if (trim(mol2_type) == 'C.2') then 
    mol2_idx = C2
else if (trim(mol2_type) == 'C.cat') then 
    mol2_idx = Ccat
else if (trim(mol2_type) == 'N.3') then 
    mol2_idx = N3
else if (trim(mol2_type) == 'N.am') then 
    mol2_idx = Nam
else if (trim(mol2_type) == 'N.ar') then 
    mol2_idx = Nar
else if (trim(mol2_type) == 'N.pl3') then 
    mol2_idx = Npl3
else if (trim(mol2_type) == 'O.2') then 
    mol2_idx = O2
else if (trim(mol2_type) == 'O.3') then 
    mol2_idx = O3
else if (trim(mol2_type) == 'O.co2') then 
    mol2_idx = Oco2
else if (trim(mol2_type) == 'O.w') then
    mol2_idx = O3
else if (trim(mol2_type) == 'Br') then
    mol2_idx = Br
else if (trim(mol2_type) == 'Cl') then 
    mol2_idx = Cl
else if (trim(mol2_type) == 'F') then 
    mol2_idx = Fluo
else if (trim(mol2_type) == 'I') then 
    mol2_idx = Iod
else if (trim(mol2_type) == 'Met') then 
    mol2_idx = Met
else if (trim(mol2_type) == 'P.3') then 
    mol2_idx = P3
else if (trim(mol2_type) == 'S.3') then 
    mol2_idx = S3
else if (trim(mol2_type) == 'N.4') then
    mol2_idx = N3
else if (trim(mol2_type) == 'N.2') then
    mol2_idx = Nar
else if (trim(mol2_type) == 'S.2') then
    mol2_idx = S3
else if (trim(mol2_type) == 'Ca' .or. trim(mol2_type) == 'Mg' &
         .or. trim(mol2_type) == 'Zn' .or. trim(mol2_type) == 'Hg' &
         .or. trim(mol2_type) == 'Mn' .or. trim(mol2_type) == 'Co.oh'&
         .or. trim(mol2_type) == 'Fe' .or. trim(mol2_type) == 'Ni' &
         .or. trim(mol2_type) == 'Cu' .or. trim(mol2_type) == 'Cd') then
    mol2_idx = Met
end if

end subroutine assign_DS_index
!-------------------------------------------------------------------------------
subroutine construct_DS_para(DS_types, DS_para_idx, flex_res)
!-------------------------------------------------------------------------------
! Fill DS_para and DS_para_idx
!-------------------------------------------------------------------------------
integer, intent(in) :: DS_types(:)
integer, intent(out) :: DS_para_idx(:)
logical, intent(in) :: flex_res(:)
integer :: lig_no, res_no, i_atm
integer :: n_prot_type, n_lig_type, type_idx
logical :: new_type

n_prot_type = 0
n_lig_type = 0
is_ligand(:) = .false.

do i_atm = 1, tn%nonligatm ! for protein/small cofactor atoms
    if (DS_types(i_atm) == -1) then
        DS_para_idx(i_atm) = -1
        cycle
    end if
    res_no = i_R(1, i_atm)
    if (res_no <= tn%stdres) then
        if (flex_res(res_no)) then
            call check_is_new_type(DS_para%type_idx(:,2), n_lig_type, &
                                   DS_types(i_atm), new_type, type_idx)
            if (new_type) then
                n_lig_type = n_lig_type + 1
                DS_para%type_idx(n_lig_type, 2) = DS_types(i_atm)
                DS_para_idx(i_atm) = n_lig_type
            else
                DS_para_idx(i_atm) = type_idx
            end if
        else
            call check_is_new_type(DS_para%type_idx(:,1), n_prot_type, &
                                   DS_types(i_atm), new_type, type_idx)
            if (new_type) then
                n_prot_type = n_prot_type + 1
                DS_para%type_idx(n_prot_type, 1) = DS_types(i_atm)
                DS_para_idx(i_atm) = n_prot_type
            else
                DS_para_idx(i_atm) = type_idx
            end if
        end if
    else
        call check_is_new_type(DS_para%type_idx(:,1), n_prot_type, &
                               DS_types(i_atm), new_type, type_idx)
        if (new_type) then
            n_prot_type = n_prot_type + 1
            DS_para%type_idx(n_prot_type, 1) = DS_types(i_atm)
            DS_para_idx(i_atm) = n_prot_type
        else
            DS_para_idx(i_atm) = type_idx
        end if
    end if
end do

do i_atm = tn%nonligatm + 1, tn%atom ! for ligand/large cofactor atoms
    if (DS_types(i_atm) == -1) then
        DS_para_idx(i_atm) = -1
        cycle
    end if
    lig_no = i_L(1, i_atm-tn%nonligatm)
    if (lig_no + tn%nonlig == target_lig_no) then ! for target ligand
        is_ligand(i_atm) = .true.
        call check_is_new_type(DS_para%type_idx(:,2), n_lig_type, &
                               DS_types(i_atm), new_type, type_idx)
        if (new_type) then
            n_lig_type = n_lig_type + 1
            DS_para%type_idx(n_lig_type,2) = DS_types(i_atm)
            DS_para_idx(i_atm) = n_lig_type
        else
            DS_para_idx(i_atm) = type_idx
        end if
    else ! for large cofactor atoms
        call check_is_new_type(DS_para%type_idx(:,1), n_prot_type, &
                               DS_types(i_atm), new_type, type_idx)
        if (new_type) then
            n_prot_type = n_prot_type + 1
            DS_para%type_idx(n_prot_type, 1) = DS_types(i_atm)
            DS_para_idx(i_atm) = n_prot_type
        else 
            DS_para_idx(i_atm) = type_idx
        end if
    end if
end do
DS_para%n_atm_types(1) = n_prot_type
DS_para%n_atm_types(2) = n_lig_type

end subroutine construct_DS_para
!===============================================================================
! Construct DrugScore Grid
!===============================================================================
subroutine construct_drugscore_grid(grid_info, dock_grid, flex_res)
!-------------------------------------------------------------------------------
! Construct drugscore grid (dock_grid%drugscore_grid)
!-------------------------------------------------------------------------------
type(docking_grid_param_type), intent(in) :: grid_info
type(docking_grid_type), intent(inout) :: dock_grid
logical, intent(in) :: flex_res(:)
integer :: n_hash
real(dp), allocatable :: tmp_grid_E(:)
real(dp) :: ref_pt(3), grid_pt(3)
integer :: i_lig, i, x, y, z

n_hash = grid_info%n_elem(1) * grid_info%n_elem(2) * grid_info%n_elem(3) 
! allocate drugscore grid
allocate(dock_grid%drugscore_grid(4, n_hash, DS_para%n_atm_types(2)))
allocate(tmp_grid_E(DS_para%n_atm_types(2)))

tmp_grid_E(:) = 0.0d0

do i = 1, 3
    ref_pt(i) = grid_info%grid_cntr(i) &
              - (dble(grid_info%n_elem(i))-1.0d0)/2.0d0 * grid_info%grid_width
end do

do z = 1, grid_info%n_elem(3)
    grid_pt(3) = (dble(z)-1.0d0)*grid_info%grid_width + ref_pt(3)
    do y = 1, grid_info%n_elem(2)
        grid_pt(2) = (dble(y)-1.0d0)*grid_info%grid_width + ref_pt(2)
        do x = 1, grid_info%n_elem(1)
            grid_pt(1) = (dble(x)-1.0d0)*grid_info%grid_width + ref_pt(1)
            call calc_DS_grid(grid_pt(1:3), tmp_grid_E, flex_res)
            do i_lig = 1, DS_para%n_atm_types(2)
                call copy_to_grid(tmp_grid_E(i_lig), grid_info%n_elem(1:3),&
                                  dock_grid%drugscore_grid(:,:,i_lig), x, y, z,&
                                  n_hash)
            end do
        end do
    end do
end do

deallocate(tmp_grid_E)

end subroutine construct_drugscore_grid
!-------------------------------------------------------------------------------
subroutine calc_DS_grid(grid_pt, E_pt, flex_res) 
!-------------------------------------------------------------------------------
! calculate DrugScore interaction energy on grid point
!-------------------------------------------------------------------------------
real(dp), intent(in) :: grid_pt(3)
real(dp), intent(inout) :: E_pt(:)
logical, intent(in) :: flex_res(:)
integer :: n_atm_types, prot_type, lig_type, i_type
integer :: res_no, atm_no, i_atm
integer :: i_bin
real(dp) :: dr(3), dist
real(dp) :: tmp_E, DS_E

n_atm_types = DS_para%n_atm_types(2)
E_pt(:) = 0.0d0

do res_no = 1, tn%residue
    if (res_no == target_lig_no) cycle
    if (res_no <= tn%stdres) then
        if (flex_res(res_no)) cycle
    end if
    
    do atm_no= 1, res_index(res_no)%n_atm
        if (res_no > tn%nonlig) then
            i_atm = ii_L(atm_no, res_no-tn%nonlig)
        else
            i_atm = ii_R(atm_no, res_no)
        end if
        if (DS_para_idx(i_atm) == -1) cycle

        dr(1:3) = R(1:3, i_atm) - grid_pt(1:3)
        dist = dot_product(dr,dr)

        if (dist > DScutoff2) cycle
        dist = sqrt(dist)
        prot_type = DS_para%type_idx(DS_para_idx(i_atm), 1)
        
        ! calculate grid point energy
        do i_type = 1, n_atm_types
            tmp_E = 0.0d0
            lig_type = DS_para%type_idx(i_type, 2)

            ! bin index of DrugScore
            i_bin = int(10.0*(dist-1.0) + 1.0)

            ! calculate energy at given grid point
            ! if dist> r_max of drugscore, the score decreases linearly.
            if (i_bin >= 50) then
                DS_E = DS_record(prot_type, lig_type)%bin_score(50)
                tmp_E = DS_E * (1.0 - ((dist-5.9) / (DSCutoff-5.9)))
            ! if dist < r_min of drugscore, the score has same value, score(r_min)
            else if (i_bin <= 1) then
                tmp_E = DS_record(prot_type, lig_type)%bin_score(1)
            else
                DS_E = DS_record(prot_type, lig_type)%bin_score(i_bin)
                tmp_E = DS_E + &
                        (DS_record(prot_type, lig_type)%bin_score(i_bin+1) - DS_E) &
                         * (dist - (1.0 + 0.1*(i_bin-1)))/0.1
            end if
            E_pt(i_type) = E_pt(i_type) + tmp_E
        end do
    end do
end do

end subroutine calc_DS_grid
!===============================================================================
! Calculate DrugScore energy
!===============================================================================
subroutine calc_DS_intrxn_E_using_grid(ligand, grid_info, DS_grid, &
                                       drugscore_E, g, calc_g)
!-------------------------------------------------------------------------------
! calculate DrugScore interaction energy between protein and ligand using grid.
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand
type(docking_grid_param_type), intent(in) :: grid_info
real(dp), intent(in) :: DS_grid(:,:,:)
real(dp), intent(out) :: drugscore_E, g(:,:)
logical, intent(in) :: calc_g
integer :: i_atm, atm_idx, lig_type, i
real(dp) :: tmp_E

drugscore_E = 0.0
if (calc_g) g = 0.0d0

do i_atm = 1, ligand%n_atm
    atm_idx = ii_L(i_atm, ligand%lig_no)
    lig_type = DS_para_idx(atm_idx)
    if (lig_type == -1) cycle
    call calc_E_using_grid(R(1:3, atm_idx), grid_info, DS_grid(:,:,lig_type),&
                           tmp_E, g(:,i_atm), calc_g)
    drugscore_E = drugscore_E + tmp_E
end do

end subroutine calc_DS_intrxn_E_using_grid
!-------------------------------------------------------------------------------
subroutine calc_DrugScore(f, g, appl_res, calc_g)
!-------------------------------------------------------------------------------
! calculate DrugScore interaction energy between protein and ligand in
! continuous way.
!-------------------------------------------------------------------------------
real(dp), intent(out) :: f, g(3, tn%atom)
logical, intent(in) :: appl_res(:), calc_g
!
integer :: i_atm, j_atm, ia, atm_no, i_res
integer :: i_bin, prot_type, lig_type
real(dp) :: dist, dr(3), tmp_E

f = 0.0d0
g(:,:) = 0.0d0

do atm_no = 1, res_index(target_lig_no)%n_atm
    i_atm = ii_L(atm_no, target_lig_no-tn%nonlig)
    if (DS_para_idx(i_atm) == -1) cycle
    
    lig_type = DS_para%type_idx(DS_para_idx(i_atm),2)

    do i_res = 1, tn%residue
        if (.not. appl_res(i_res)) cycle
        if (i_res == target_lig_no) cycle
        do ia = 1, res_index(i_res)%n_atm
            j_atm = ii_R(ia,i_res)
            if (DS_para_idx(j_atm) == -1) cycle
            prot_type = DS_para%type_idx(DS_para_idx(j_atm),1)
            dr(1:3) = R(1:3, j_atm) - R(1:3, i_atm)
            dist = dot_product(dr, dr)
            if (dist > DScutoff2) cycle
            dist = sqrt(dist)
            ! calculating DrugScore
            if (dist >= 5.9d0) then
                tmp_E = DS_record(prot_type,lig_type)%bin_score(50)&
                        * (1.0d0-((dist-5.9d0)/(DSCutoff-5.9d0)))
            else if (dist < 1.0d0) then
                tmp_E = DS_record(prot_type, lig_type)%bin_score(1)
            else
                i_bin = int(10.0d0*(dist-1.0d0) + 1.0d0)
                tmp_E = DS_record(prot_type,lig_type)%bin_score(i_bin)&
                        + (DS_record(prot_type,lig_type)%bin_score(i_bin+1)&
                        - DS_record(prot_type,lig_type)%bin_score(i_bin))&
                        * (dist - (1.0d0 + 0.1d0*(i_bin-1)))/0.1d0
            end if
            f = f + tmp_E
        end do
    end do
end do

end subroutine calc_DrugScore
!-------------------------------------------------------------------------------
subroutine FFT_rec_grid_drugscore(grid_info, dock_grid, FFT_grid_info, r_grid)
!-------------------------------------------------------------------------------
type(docking_grid_param_type), intent(in) :: grid_info, FFT_grid_info
type(docking_grid_type), intent(in) :: dock_grid
type(docking_FFT_grid_type), intent(inout) :: r_grid
!
integer :: n_type, i_type, j_type
integer :: x, y, z, hash, FFT_hash
real(dp) :: grid_v

if (mod(DS_para%n_atm_types(2),2) == 0) then
    n_type = DS_para%n_atm_types(2)/2
else
    n_type = DS_para%n_atm_types(2)/2 + 1
end if
r_grid%n_drugscore = n_type

FFT_hash = FFT_grid_info%n_elem(1)*FFT_grid_info%n_elem(2)*FFT_grid_info%n_elem(3)
allocate(r_grid%drugscore_grid(FFT_hash, n_type))
                               
r_grid%drugscore_grid = cmplx(0.0,0.0)

do i_type = 1, DS_para%n_atm_types(2)
    do z = 1, FFT_grid_info%n_elem(3)
        do y = 1, FFT_grid_info%n_elem(2)
            do x = 1, FFT_grid_info%n_elem(1)
                !
                FFT_hash = FFT_grid_info%n_elem(2)*FFT_grid_info%n_elem(1)*(z-1) &
                         + FFT_grid_info%n_elem(1)*(y-1) + x
                !
                if ((x > grid_info%n_elem(1)) .or. &
                    (y > grid_info%n_elem(2)) .or. &
                    (z > grid_info%n_elem(3))) then  ! outside of docking grid
                    if (mod(i_type,2) == 1) then
                        j_type = int(i_type/2) + 1
                        r_grid%drugscore_grid(FFT_hash,j_type) = cmplx(max_energy, max_energy)
                    end if
                    cycle
                end if
                !
                hash = grid_info%n_elem(2)*grid_info%n_elem(1)*(z-1) &
                     + grid_info%n_elem(1)*(y-1) + x
                grid_v = dock_grid%drugscore_grid(1,hash,i_type)
                !
                if (mod(i_type,2) == 1) then
                    j_type = int(i_type/2) + 1
                    r_grid%drugscore_grid(FFT_hash,j_type) = &
                    r_grid%drugscore_grid(FFT_hash,j_type) + cmplx(grid_v, 0.0)
                else
                    j_type = int(i_type/2)
                    r_grid%drugscore_grid(FFT_hash,j_type) = &
                    r_grid%drugscore_grid(FFT_hash,j_type) + cmplx(0.0,grid_v)
                end if
            end do
        end do
    end do
end do

end subroutine FFT_rec_grid_drugscore
!-------------------------------------------------------------------------------
subroutine FFT_lig_grid_drugscore(FFT_grid_info, grid_info, ligand,&
                                  fragment, n_frag_atm, R_frag, l_grid, &
                                  ordering)
!-------------------------------------------------------------------------------
type(docking_grid_param_type), intent(in) :: FFT_grid_info, grid_info
type(ligand_type), intent(in) :: ligand
integer, intent(in) :: fragment(:), n_frag_atm
real(dp), intent(in) :: R_frag(:,:)
type(docking_FFT_grid_type), intent(inout) :: l_grid
logical, intent(in) :: ordering
integer :: n_type, i_type, j_type, lig_type
integer :: i_atm, atm_no
integer :: x, y, z, X0(3), FFT_hash
real(dp) :: dr(3), tcrd(3)
!
integer, allocatable, save :: order_by_type(:), type_s(:)

if (.not. allocated(l_grid%drugscore_grid)) then
    if (mod(DS_para%n_atm_types(2),2) == 0) then
        n_type = int(DS_para%n_atm_types(2)/2)
    else
        n_type = int(DS_para%n_atm_types(2)/2) + 1
    end if
    l_grid%n_drugscore = n_type

    FFT_hash = FFT_grid_info%n_elem(1)*FFT_grid_info%n_elem(2)*FFT_grid_info%n_elem(3)
    allocate(l_grid%drugscore_grid(FFT_hash, n_type))
end if

if (ordering) then
    if (allocated(order_by_type)) then
        deallocate(type_s)
        deallocate(order_by_type)
    end if
    allocate(type_s(n_frag_atm))
    allocate(order_by_type(n_frag_atm))
    do i_atm = 1, n_frag_atm
        atm_no = ii_L(fragment(i_atm), ligand%lig_no)
        type_s(i_atm) = DS_para_idx(atm_no)
    end do
    call sort1(n_frag_atm, type_s, order_by_type)
    return
end if

l_grid%drugscore_grid = cmplx(0.0,0.0)

do atm_no = 1, n_frag_atm
    i_atm = order_by_type(atm_no)
    i_type = type_s(atm_no)
    if (i_type == -1) cycle
    !
    dr(:) = R_frag(:,i_atm) - FFT_grid_info%grid_cntr(:)
    tcrd(:) = dr(:)/FFT_grid_info%grid_width
    tcrd(:) = tcrd(:) + (dble(grid_info%n_elem(:)-1)*0.5d0 + 1.0d0)
    X0(:) = nint(tcrd(:)) ! nearest grid point
    !
    x = X0(1)
    y = X0(2)
    z = X0(3)
    FFT_hash = FFT_grid_info%n_elem(2)*FFT_grid_info%n_elem(1)*(z-1) &
             + FFT_grid_info%n_elem(1)*(y-1) + x
    !
    if (mod(i_type,2) == 1) then
        j_type = int(i_type/2) + 1
        l_grid%drugscore_grid(FFT_hash,j_type) = l_grid%drugscore_grid(FFT_hash,j_type) + & 
                                         cmplx(1.0, 0.0)
    else
        j_type = int(i_type/2)
        l_grid%drugscore_grid(FFT_hash,j_type) = l_grid%drugscore_grid(FFT_hash,j_type) + &
                                         cmplx(0.0, 1.0)
    end if
end do

end subroutine FFT_lig_grid_drugscore
!-------------------------------------------------------------------------------
END MODULE DRUGSCORE
!-------------------------------------------------------------------------------
