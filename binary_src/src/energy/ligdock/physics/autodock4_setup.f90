!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/ligdock/physics/autodock4_setup.f90
!
! Description: Subroutines related setup autodock4 parameters
!
!-------------------------------------------------------------------------------
MODULE AUTODOCK4_SETUP
!-------------------------------------------------------------------------------

use globals
use logger, only: terminate_with_error
use string, only: parse_string
!
use energy_vars
use ligdock_E_utils, only: check_is_new_type
! 
use autodock_vars
use autodock_utils, only: initialize_atdk_para, set_nb_matrix, set_energy_table,&
                          HD_hbond_parameter, NA_hbond_parameter, &
                          OA_SA_hbond_parameter

implicit none
save
private

public :: read_ad4_param_file
public :: construct_atdk4_para
public :: setup_atdk4_E

CONTAINS
!-------------------------------------------------------------------------------
subroutine read_ad4_param_file(file_name, ad4_param)
!-------------------------------------------------------------------------------
! subroutine for reading AutoDock4 parameter file
! Parameters are saved in ad4_param
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: file_name
type(ad4_param_type), intent(out) :: ad4_param
integer :: f_unit, ioerror
!
integer :: num_word
character(len=len_fname) :: line, word(25)
!
integer :: i_atm

f_unit = 19
open(unit=f_unit, file=file_name, status='old', action='read', iostat=ioerror)
if (ioerror /= 0) then
    call terminate_with_error("Error : Can not read autodock energy parameter file")
end if

i_atm = 0
do
    read(f_unit, '(a)', iostat=ioerror) line
    if(ioerror /= 0) exit

    call parse_string(line, num_word, word)
    if(word(1) == 'atom_par') then
        i_atm = i_atm + 1
        ad4_param%atm_type(i_atm) = word(2)
        read(word(3), '(f4.2)') ad4_param%Rii(i_atm)
        read(word(4), '(f5.3)') ad4_param%epsii(i_atm)
        read(word(5), '(f7.4)') ad4_param%vol(i_atm)
        read(word(6), '(f8.5)') ad4_param%solpar(i_atm)
        read(word(7), '(f3.1)') ad4_param%Rij_hb(i_atm)
        read(word(8), '(f3.1)') ad4_param%epsij_hb(i_atm)
        read(word(9), '(i1)') ad4_param%hbond(i_atm)
        read(word(12), '(i1)') ad4_param%bond_index(i_atm)
    end if
end do
close(f_unit)

end subroutine read_ad4_param_file
!-------------------------------------------------------------------------------
subroutine construct_atdk4_para(ad4_param, atdk_types, target_lig_no, &
                                is_Hdon, is_Hacc_N, is_Hacc_O, is_ligand)
!-------------------------------------------------------------------------------
! Construct atdk_para, intern_para which contain parameters for used atoms only.
! In here, atdk_para_idx is filled to access atdk_para.
! atdk_para is used to describe protein-ligand interaction, and intern_para is
! used to describe interaction within ligand.
!-------------------------------------------------------------------------------
type(ad4_param_type), intent(in) :: ad4_param
integer, intent(in) :: atdk_types(:)
integer, intent(in) :: target_lig_no
logical, intent(in) :: is_Hdon(:), is_Hacc_N(:), is_Hacc_O(:)
logical, intent(out) :: is_ligand(:)
integer :: lig_no, res_no, i_atm, atm_no
integer :: n_tot_type, n_lig_type, type_idx
logical :: new_type

call initialize_atdk_para(atdk_para)

n_tot_type = 0
n_lig_type = 0
is_ligand(:) = .false.

do i_atm = 1, tn%nonligatm ! for protein/confactor(hetmol) atoms
    res_no = i_R(1, i_atm)
    call check_is_new_type(atdk_para%type_idx(:), n_tot_type, &
                           atdk_types(i_atm), new_type, type_idx)
    if (new_type) then
        n_tot_type = n_tot_type + 1
        atdk_para%type_idx(n_tot_type) = atdk_types(i_atm)
        atdk_para%type_name(n_tot_type) = ad4_param%atm_type(atdk_types(i_atm))
        atdk_para%is_Hdon(n_tot_type) = is_Hdon(i_atm)
        atdk_para%is_Hacc_N(n_tot_type) = is_Hacc_N(i_atm)
        atdk_para%is_Hacc_O(n_tot_type) = is_Hacc_O(i_atm)
        atdk_para_idx(i_atm) = n_tot_type
    else
        atdk_para_idx(i_atm) = type_idx
    end if
end do

do i_atm = tn%nonligatm + 1, tn%atom ! for ligand/large cofactor atoms
    lig_no = i_L(1, i_atm - tn%nonligatm)
    atm_no = i_L(2, i_atm - tn%nonligatm) 
    call check_is_new_type(atdk_para%type_idx(:), n_tot_type, &
                           atdk_types(i_atm), new_type, type_idx)
    if (new_type) then
        n_tot_type = n_tot_type + 1
        atdk_para%type_idx(n_tot_type) = atdk_types(i_atm)
        atdk_para%type_name(n_tot_type) = ad4_param%atm_type(atdk_types(i_atm))
        atdk_para%is_Hdon(n_tot_type) = is_Hdon(i_atm)
        atdk_para%is_Hacc_N(n_tot_type) = is_Hacc_N(i_atm)
        atdk_para%is_Hacc_O(n_tot_type) = is_Hacc_O(i_atm)
        atdk_para_idx(i_atm) = n_tot_type
        type_idx = n_tot_type
    else
        atdk_para_idx(i_atm) = type_idx
    end if
    if (lig_no == target_lig_no) then ! for target ligand
        is_ligand(i_atm) = .true.
        call check_is_new_type(atdk_para%lig_atom_types(:), n_lig_type, &
                               atdk_para_idx(i_atm), new_type, type_idx)
        if (new_type) then
            n_lig_type = n_lig_type + 1
            atdk_para%lig_atom_types(n_lig_type) = atdk_para_idx(i_atm)
        end if
    end if
end do
atdk_para%n_atm_types(1) = n_tot_type
atdk_para%n_atm_types(2) = n_lig_type

! setup parameter for describing protein-ligand interaction
call setup_atdk4_para(ad4_param, atdk_para)

end subroutine construct_atdk4_para
!-------------------------------------------------------------------------------
subroutine setup_atdk4_para(ad4_param, para)
!-------------------------------------------------------------------------------
! Fill parameters needed in energy calculation of used atom types
!-------------------------------------------------------------------------------
type(ad4_param_type), intent(in) :: ad4_param
type(atdk_eng_para_type), intent(out) :: para
integer :: i_type, j_type
integer :: type_1, type_2
real(dp) :: Rvdw_1, evdw_1, Rhb_1, ehb_1
real(dp) :: Rvdw_2, evdw_2, Rhb_2, ehb_2
integer :: hb_idx_1, hb_idx_2

do i_type = 1, para%n_atm_types(1)
    type_1 = para%type_idx(i_type)

    Rvdw_1 = ad4_param%Rii(type_1)
    evdw_1 = ad4_param%epsii(type_1)
    Rhb_1 = ad4_param%Rij_hb(type_1)
    ehb_1 = ad4_param%epsij_hb(type_1)
    hb_idx_1 = ad4_param%hbond(type_1)
   
    ! save solvation parameter
    para%solpar(i_type) = ad4_param%solpar(type_1)
    para%vol(i_type) = ad4_param%vol(type_1)

    do j_type = 1, para%n_atm_types(1)
        type_2 = para%type_idx(j_type)
        
        Rvdw_2 = ad4_param%Rii(type_2)
        evdw_2 = ad4_param%epsii(type_2)
        Rhb_2 = ad4_param%Rij_hb(type_2)
        ehb_2 = ad4_param%epsij_hb(type_2)
        hb_idx_2 = ad4_param%hbond(type_2)

        ! save vdw parameter
        para%ljrij(j_type, i_type) = (Rvdw_1 + Rvdw_2) / 2.0d0
        para%ljeij(j_type, i_type) = sqrt(evdw_1 * evdw_2)

        ! save hbond parameter
        if (hb_idx_1 > 2 .and. (hb_idx_2 == 1 .or. hb_idx_2 == 2)) then
            ! idx1: H-bond acceptor ; idx2: H-bond donor
            para%hbond_read(j_type, i_type) = .true.
            para%hbrij(j_type, i_type) = Rhb_1
            para%hbeij(j_type, i_type) = ehb_1
        else if ((hb_idx_1 == 1 .or. hb_idx_1 == 2) .and. hb_idx_2 > 2) then
            ! idx1: H-bond donor ; idx2: H-bond acceptor
            para%hbond_read(j_type, i_type) = .true.
            para%hbrij(j_type, i_type) = Rhb_2
            para%hbeij(j_type, i_type) = ehb_2
        else
            ! It's not H-bond pair
            para%hbond_read(j_type, i_type) = .false.
            para%hbrij(j_type, i_type) = 0.0d0
            para%hbeij(j_type, i_type) = 0.0d0
        end if
    end do
end do

do j_type = 1, para%n_atm_types(1)
    type_2 = para%type_idx(j_type)
    para%solpar(j_type) = ad4_param%solpar(type_2)
    para%vol(j_type) = ad4_param%vol(type_2)
end do

end subroutine setup_atdk4_para
!-------------------------------------------------------------------------------
subroutine setup_atdk4_E(ligand, is_Hdon, is_Hacc_N, is_Hacc_O)
!-------------------------------------------------------------------------------
! TODO: comment
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand
logical, intent(in) :: is_Hdon(:), is_Hacc_N(:), is_Hacc_O(:)

call set_nb_matrix(ligand)
call set_parameters_for_intern_E(ligand)
call set_energy_table(atdk_para%n_atm_types(1), atdk_para%n_atm_types(2))

allocate(atdk_hbond(tn%atom))
call initialize_atdk4_Hbond(atdk_hbond, is_Hdon, is_Hacc_N, is_Hacc_O)

end subroutine setup_atdk4_E
!-------------------------------------------------------------------------------
subroutine set_parameters_for_intern_E(ligand)
!-------------------------------------------------------------------------------
type(ligand_type), intent(in) :: ligand
integer :: n_pair, i_atm, j_atm, atm_no_1, atm_no_2, i, j
real(dp) :: q1, q2, SiVj, SjVi

n_pair = 0
do i_atm = 1, ligand%n_atm - 1
    atm_no_1 = ii_L(i_atm, ligand%lig_no)
    q1 = ref_lig(ligand%lig_type)%charge(i_atm)

    do j_atm = i_atm + 1, ligand%n_atm
        if (nb_matrix(j_atm, i_atm) == 0) cycle
        !
        n_pair = n_pair + 1
        !
        atm_no_2 = ii_L(j_atm, ligand%lig_no)
        q2 = ref_lig(ligand%lig_type)%charge(j_atm)
        !
        qiqj(n_pair) = q1*q2
        !
        SiVj = atdk_para%vol(atdk_para_idx(atm_no_2))&
              * (atdk_para%solpar(atdk_para_idx(atm_no_1)) + solpar_q*abs(q1))
        SjVi = atdk_para%vol(atdk_para_idx(atm_no_1)) &
              * (atdk_para%solpar(atdk_para_idx(atm_no_2)) + solpar_q*abs(q2))
        sum_SiVj(n_pair) = SiVj + SjVi
    end do
end do

end subroutine set_parameters_for_intern_E
!-------------------------------------------------------------------------------
subroutine initialize_atdk4_Hbond(atdk_hbond, is_Hdon, is_Hacc_N, is_Hacc_O)
!-------------------------------------------------------------------------------
! TODO: comment
!-------------------------------------------------------------------------------
type(atdk_hbond_type), intent(inout) :: atdk_hbond(:)
logical, intent(in) :: is_Hdon(:), is_Hacc_N(:), is_Hacc_O(:)
integer :: i_atm

! initialize
do i_atm = 1, tn%atom
    atdk_hbond(i_atm)%rexp = 1
    atdk_hbond(i_atm)%rvec(:,:) = 0.0d0
    atdk_hbond(i_atm)%racc = 1.0d0
    atdk_hbond(i_atm)%rdon = 1.0d0
end do

do i_atm = 1, tn%atom ! for protein/small cofactor atom
    if (is_Hdon(i_atm)) then
        call HD_hbond_parameter(atdk_hbond(i_atm), i_atm)
    else if (is_Hacc_N(i_atm)) then
        call NA_hbond_parameter(atdk_hbond(i_atm), i_atm)
    else if (is_Hacc_O(i_atm)) then
        call OA_SA_hbond_parameter(atdk_hbond(i_atm), i_atm)
    end if
end do

end subroutine initialize_atdk4_Hbond
!-------------------------------------------------------------------------------
END MODULE AUTODOCK4_SETUP
!-------------------------------------------------------------------------------
