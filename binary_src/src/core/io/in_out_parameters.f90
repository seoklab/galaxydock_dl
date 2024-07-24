!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/core/io/in_out_parameters.f90
!
! Description:
!   This module contains subroutines for read and setup parameters and
!   topologies
!-------------------------------------------------------------------------------
MODULE IN_OUT_PARAMETERS
!-------------------------------------------------------------------------------
use globals
use logger, only: log_p
use geometry, only: ligand_to_topology, identify_t_ang_dependence
use string, only: parse_string
!
use in_out_vars
use in_out_utils
use in_out_ligand, only: read_mol2, merge_nonpolar_hydrogens
!
use rotamer, only: initialize_rotamer, finalize_rotamer

implicit none
private

public :: initialize_parameters
public :: finalize_parameters

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_parameters(read_rotamer)
!-------------------------------------------------------------------------------
! Read topology/parameters and then Setup reference residue
!-------------------------------------------------------------------------------
logical, intent(in), optional :: read_rotamer
!
call read_biomole_parameters()
call log_p('- Finished with read_parameters.', me=me, level=20)

call read_topology_geo()
call log_p('- Finished with read_topology_geo.', me=me, level=20)

call read_topology_eng()
call log_p('- Finished with read_topology_eng.', me=me, level=20)

if (read_het .and. n_mol2_top /= 0) then
    call read_topology_mol2(ref_lig)
    call log_p('- Finished with read_topology_mol2.', me=me, level=20)
    call read_ligand_parameters(ref_lig)
    call log_p('- Finished with read_parameters for ligands.', me=me, level=20)
    call assign_mol2_ref_res_eng(ref_lig)
end if

call setup_vdw_para()
call setup_ref_res_aux()

if ((present(read_rotamer) .and. (.not. read_rotamer)) .or.&
    top_type == 'coarse' .or. top_type == 'cbeta') return

call initialize_rotamer()


end subroutine initialize_parameters
!-------------------------------------------------------------------------------
subroutine finalize_parameters()
!-------------------------------------------------------------------------------
if (top_type /= 'coarse' .and. top_type /= 'cbeta') then
    call finalize_rotamer()
end if

end subroutine finalize_parameters
!-------------------------------------------------------------------------------
subroutine read_biomole_parameters()
!-------------------------------------------------------------------------------
! Read molecular mechanics energy parameters
!-------------------------------------------------------------------------------
call read_parameters(infile_parameter, eng_para)

end subroutine read_biomole_parameters
!-------------------------------------------------------------------------------
subroutine read_topology_geo()
!-------------------------------------------------------------------------------
! Read topology file and save information in ref_res
!-------------------------------------------------------------------------------
integer :: f_unit = 11, ioerror
character(len=len_fname) :: line
integer :: i_file, i_res, i_atm, i_bnd, j_bnd, bnd_idx
character(len=4) :: res_name, atom_name
character(len=6) :: atom_cls
integer :: i_atm1, i_atm2, i_atm3, i_atm4, n_atm, group
real(dp) :: b_len0, b_ang0, t_ang0, charge
integer :: atm(1:4), atm_in_bnd(2), atm_in_ref_bnd(2)
integer :: openstat
character(len=len_fname) :: error_line
integer :: q_typ_no
character(len=1) :: ter_type
integer :: k_bnd, i_atm_prev
character(len=4) :: atm_name(2), atm_name_ref(2), atom_prev_aa(3)

q_typ_no = 0

! read and save residue topology in ref_res
i_res = 0
do i_file = 1, n_topo_file
    open(f_unit, file = infile_topo(i_file), iostat = openstat)
    if (openstat > 0) then
        write(log_msg,"(A,A)") 'Terminate with error: Topology file not found, ', &
                               trim(infile_topo(i_file))
        call log_p(log_msg)
        call terminate_with_error('Please check {data_directory}.')
    end if

    write(log_msg,"(A,A)") '  Reading topology file: '//trim(infile_topo(i_file))
    call log_p(log_msg,me=me,level=10)

    ter_type = topo_file_ter_type(i_file)

    do
        read(f_unit, "(a100)", iostat = ioerror) line

        if (ioerror < 0) exit
        if (line(1:4) == 'DONE') then
            ref_res(i_res)%res_name = res_name
            ref_res(i_res)%n_atm = n_atm
            ref_res(i_res)%n_gr = n_atm
            ref_res(i_res)%ter_type = ter_type
            if (n_atm > max_atm) then
                write(error_line, *) 'n_atm > max_atm', n_atm, max_atm
                call terminate_with_error(error_line)
            end if

            ! bnd index in bnd angle
            do i_atm = 0, n_atm
                atm(1:3) = ref_res(i_res)%atm_in_b_ang(1:3,i_atm)
                do i_bnd = 1, 2
                    if (i_bnd == 1) then
                        atm_in_bnd(1:2) = atm(1:2)
                    else
                        atm_in_bnd(1:2) = atm(2:3)
                    end if
                    do j_bnd = -1, n_atm
                        atm_in_ref_bnd(1:2) = ref_res(i_res)%atm_in_bnd(1:2,j_bnd)
                        if (atm_in_bnd(1) == atm_in_ref_bnd(1) .and. &
                            atm_in_bnd(2) == atm_in_ref_bnd(2)) then
                            bnd_idx = j_bnd
                            exit
                        end if
                    end do
                    ref_res(i_res)%bnd_in_b_ang(i_bnd,i_atm) = bnd_idx
                end do
            end do

            ! bnd index in tor angle
            do i_atm = 1, n_atm
                atm(1:4) = ref_res(i_res)%atm_in_t_ang(1:4,i_atm)
                do i_bnd = 1, 3
                    if (i_bnd == 1) then
                        atm_in_bnd(1:2) = atm(1:2)
                    else if (i_bnd == 2) then
                        atm_in_bnd(1:2) = atm(2:3)
                    else
                        atm_in_bnd(1:2) = atm(3:4)
                    end if
                    do j_bnd = -1, n_atm
                        atm_in_ref_bnd(1:2) = ref_res(i_res)%atm_in_bnd(1:2,j_bnd)
                        if (atm_in_bnd(1) == atm_in_ref_bnd(1) .and. &
                            atm_in_bnd(2) == atm_in_ref_bnd(2)) then
                            bnd_idx = j_bnd
                            exit
                        end if
                    end do
                    ref_res(i_res)%bnd_in_t_ang(i_bnd,i_atm) = bnd_idx
                end do
            end do

            ! indices for connecting atoms
            if (topo_file_mol_type(i_file) == 'protein') then

                ! N, CA, C, O, atom indices
                if (top_type == 'coarse') then
                    ref_res(i_res)%atom_name(-3:0) = (/ '    ', '-N  ', '-CA ', '-C  ' /)
                else
                    ref_res(i_res)%atom_name(-3:0) = (/ '-O  ', '-N  ', '-CA ', '-C  ' /)
                end if
                ref_res(i_res)%atom_name(ref_res(i_res)%n_atm + 1) = '+N  '
                atom_prev_aa(1:3) = (/ 'N   ', 'CA  ', 'C   ' /)
                if (ref_res(i_res)%res_name == 'NACE') &
                    atom_prev_aa(1:3) = (/ 'O   ', 'CH3 ', 'C   ' /)

                do i_atm = 1, ref_res(i_res)%n_atm
                    do i_atm_prev = 1, 3
                        if (ref_res(i_res)%atom_name(i_atm) == atom_prev_aa(i_atm_prev)) then
                            ref_res(i_res)%i_atm_prev(i_atm_prev) = i_atm
                        end if
                    end do
                    if (ref_res(i_res)%atom_name(i_atm) == 'O') then
                        ref_res(i_res)%i_atm_o = i_atm
                    end if
                end do
     
                ! N-CA, CA-C bond indices
                do k_bnd = 1, 2
                    atm_name(1:2) = atom_prev_aa(k_bnd:k_bnd+1)
           
                    ! go over bonds
                    do i_bnd = 1, ref_res(i_res)%n_atm
                        atm_in_bnd(1:2) = ref_res(i_res)%atm_in_bnd(1:2,i_bnd)
                        atm_name_ref(1:2) = ref_res(i_res)%atom_name(atm_in_bnd(1:2))
                        if (atm_name(1) == atm_name_ref(1) .and. atm_name(2) == atm_name_ref(2)) then
                            bnd_idx = i_bnd
                            exit
                        end if
                    end do
           
                    ref_res(i_res)%i_bnd_prev(k_bnd) = i_bnd
                end do

                ! phi, psi angle indexing
                ! phi angle: torsion angle that places C in the current residue
                ! psi angle: torsion angle that places N in the next residue (+N)
                ! omg angle: torsion angle that places CA in the next residue (+CA)
                ! iphi: phi angle index of the current residue (defines current res C)
                ! ipsi: psi angle index of the previous residue (defines current res N)
                ! iomg: omg angle index of the previous residue (defines current res CA)
                do i_atm = 1, ref_res(i_res)%n_atm
                    if (ref_res(i_res)%atom_name(i_atm) == 'C') then
                        ref_res(i_res)%iphi = i_atm
!                        ref_res(i_res)%iphi_curr = i_atm
                    else if (ref_res(i_res)%atom_name(i_atm) == 'N') then
                        ref_res(i_res)%ipsi = i_atm
!                        ref_res(i_res)%ipsi_prev = i_atm
                    else if (ref_res(i_res)%atom_name(i_atm) == 'CA') then
                        ref_res(i_res)%iomg = i_atm
!                        ref_res(i_res)%iomg_curr = i_atm
                    end if
                end do

            end if

        else if (line(1:2) /= '') then

            ! residue name line
            res_name = line(1:4)
            ! increment residue index
            i_res = i_res + 1

        else
            ! atom line
            read(line,*) i_atm4, atom_name, atom_cls, &
                         i_atm3, i_atm2, i_atm1, b_len0, b_ang0, t_ang0, charge, group 
            i_atm4 = i_atm4 - 3
            i_atm3 = i_atm3 - 3
            i_atm2 = i_atm2 - 3
            i_atm1 = i_atm1 - 3
            n_atm = i_atm4

            if (i_atm4 >= -1) then
                ref_res(i_res)%b_len0(i_atm4) = b_len0
                ref_res(i_res)%atm_in_bnd(1:2, i_atm4) = (/ i_atm3, i_atm4 /)
            end if
            if (i_atm4 >= 0) then
                ref_res(i_res)%b_ang0(i_atm4) = b_ang0*deg2rad
                ref_res(i_res)%atm_in_b_ang(1:3, i_atm4) = (/ i_atm2, i_atm3, i_atm4 /)
            end if
            if (i_atm4 >= 1) then
                ref_res(i_res)%atom_name(i_atm4) = atom_name
                ref_res(i_res)%atm_in_t_ang(1:4, i_atm4) = (/ i_atm1, i_atm2, i_atm3, i_atm4 /)
                ref_res(i_res)%t_ang0(i_atm4) = t_ang0*deg2rad
                ref_res(i_res)%group(i_atm4) = group
                ! charge
                ref_res_eng(i_res)%atom_cls(i_atm4) = atom_cls
                call find_atom_cls(eng_para, atom_cls, ref_res_eng(i_res)%atm_cls(i_atm4))
                q_typ_no = q_typ_no + 1
                eng_para%charge(q_typ_no) = charge
                ref_res_eng(i_res)%qq_idx(i_atm4) = q_typ_no
                if (force_field_type == 'AMBER') then
                    eng_para%charge(q_typ_no) = eng_para%charge(q_typ_no)*electron2kcal
                end if
            end if
        end if
    end do

    close(f_unit)
end do

num_ref_res = i_res
eng_para%n_q_typ = q_typ_no

end subroutine read_topology_geo
!-------------------------------------------------------------------------------
subroutine read_topology_eng()
!-------------------------------------------------------------------------------
! build the reference residue structure for energy calculation. 
!-------------------------------------------------------------------------------
character(len=len_fname) :: line, error_line
type(ref_res_type) :: ref
type(ref_res_eng_type) :: eref
character(len=4) :: string, res_name, atom_name, atom_cls
integer :: f_unit = 19, ref_res_no, atm_no, ang_no, i, idx, imax, ang_type
integer :: indx, i_ang, i_bnd, bnd_no, i_ref_res, i_file
character(len=6), dimension(4) :: atm
character(len=4), dimension(4) :: atom
integer :: ioerror, openstat
character(len=6) :: atom_error_mode = ''
character(len=10) :: topo_mol_type
logical :: status
logical :: pbond_found ! whether +N or -CA is included in angle
character(len=4), dimension(4) :: atom_in_bnd_E_save, atom_in_ang_E_save
integer, dimension(4) :: atm_in_bnd_E_save, atm_in_ang_E_save
integer :: ang_type_save, imax_save
integer :: ref_res_no_PRO, ref_res_no_GLY, ref_res_no_reg

res_name = 'PRO'
call find_res_idx(res_name, ref_res_no_PRO)
res_name = 'GLY'
call find_res_idx(res_name, ref_res_no_GLY)
res_name = 'ALA' ! regular residue
call find_res_idx(res_name, ref_res_no_reg)

do i_file = 1, n_topo_eng_file
    open(f_unit, file = infile_topo_eng(i_file), iostat = openstat)
    write(log_msg,"(A,A)") '  Reading topology file: '//trim(infile_topo_eng(i_file))
    call log_p(log_msg, me=me, level=10)
    if (openstat > 0) then
        call log_p('Terminate with error: Topology_energy file not found.')
        call terminate_with_error('Please check {data_directory}.')
    end if

    topo_mol_type = trim(topo_eng_file_mol_type(i_file))

    ref_res_no = 0
    do
        read(f_unit,"(A100)", iostat = ioerror) line
        if (ioerror < 0) exit

        if (line(1:3) == 'END') then
            ! complete structure of previous residue
            eref%n_ang_E = ang_no
            eref%n_bnd_E = bnd_no
            if (topo_mol_type == 'protein') then
                ! needed for constructing pair list:
                if (top_type == 'coarse') then
                    eref%atom_in_bnd_E(1:2,-1) = (/ '-CA ', '-C  ' /)
                    eref%atom_in_bnd_E(1:2,0) = (/ '-C  ', 'N   ' /)
                else
                    eref%atom_in_bnd_E(1:2,-1)= (/ '-CA ', '-C  ' /)
                    eref%atom_in_bnd_E(1:2,-2)= (/ '-C  ', '-O  ' /)
                    if (res_name == 'NACE') then
                        !CURIE check this
                        if (trim(top_type) /= 'allh' .and. trim(top_type) /= 'allh_ch22') then
                            eref%atom_in_bnd_E(1:2,0) = (/ '-C  ', 'CH3 ' /)
                        else
                            eref%atom_in_bnd_E(1:2,0) = (/ '-C  ', 'HH31' /)
                        end if
                    else
                        eref%atom_in_bnd_E(1:2,0) = (/ '-C  ', 'N   ' /)
                    end if
                end if
            else if (read_het .and. topo_mol_type == 'hetmol') then
            else
                call terminate_with_error('Error. Unknown molecule type '//topo_mol_type)
            end if

            ! For chain connectivity
            if (mol_type /= 'hetmol') then
                if (top_type == 'coarse') then
                    do bnd_no = -1, 0
                        do i = 1, 2
                            call find_atom_idx(0, ref_res_no, eref%atom_in_bnd_E(i, bnd_no), &
                                               eref%atm_in_bnd_E(i, bnd_no), atom_error_mode)
                        end do
                    end do
                else
                    do bnd_no = -2, 0
                        do i = 1, 2
                            call find_atom_idx(0, ref_res_no, eref%atom_in_bnd_E(i, bnd_no), &
                                               eref%atm_in_bnd_E(i, bnd_no), atom_error_mode)
                        end do
                    end do
                end if
            end if
           
            do i_bnd = 1, eref%n_bnd_E
                do i = 1, 2
                    atom(i) = eref%atom_in_bnd_E(i, i_bnd)
                    if (atom(i) == '+N  ') then
                        atom(i) = atom(i)(2:)
                        i_ref_res = ref_res_no_reg ! regular residue
                    else if (atom(i) == '+NP ') then
                        atom(i) = 'N   '
                        i_ref_res = ref_res_no_PRO ! for pro
                    else
                        i_ref_res = ref_res_no
                    end if
                    call find_atom_idx(0, i_ref_res, atom(i), idx, atom_error_mode)
                    atm(i) = ref_res_eng(i_ref_res)%atom_cls(idx)
                end do
                call find_bnd_prm(atm(1:2), eref%bnd_prm(i_bnd), status)
            end do
           
            ! set-up angle parameters
            do i_ang = 1, eref%n_ang_E
                do i = 1, eref%n_atm_in_ang(i_ang)
                    atom(i) = eref%atom_in_ang_E(i, i_ang)
                    if (atom(i) == '+NP ') then
                        atom(i) = 'N   '
                        i_ref_res = ref_res_no_PRO ! for pro
                    else if (atom(i) == '-CAP') then
                        atom(i) = 'CA  '
                        i_ref_res = ref_res_no_PRO ! for pro
                    else if (atom(i) == '-CAG') then
                        atom(i) = 'CA  '
                        i_ref_res = ref_res_no_GLY ! for gly
                    else if (atom(i)(1:1) == '-' .or. atom(i)(1:1) == '+') then
                        atom(i) = atom(i)(2:)
                        i_ref_res = ref_res_no_reg ! regular residue
                    else
                        i_ref_res = ref_res_no
                    end if
                    call find_atom_idx(0, i_ref_res, atom(i), idx, atom_error_mode)
                    atm(i) = ref_res_eng(i_ref_res)%atom_cls(idx)
                end do
                call find_ang_prm(atm(:), eref%ang_E_type(i_ang), eref%ang_prm(i_ang), status)
            end do
            ! save
            ref_res_eng(ref_res_no) = eref
        end if
        
        string = line(1:3)
        if (string == 'RES') then
            ! initialize RES structure
            read(line,*) string, res_name
            call find_res_idx(res_name, ref_res_no)
            ref = ref_res(ref_res_no)
            eref = ref_res_eng(ref_res_no)
            atm_no = 0
            ang_no = 0
            bnd_no = 0

        else if (string == 'ATM') then
            atm_no = atm_no + 1
           
        else if (string == 'BND') then
            ! read bond
            bnd_no = bnd_no + 1
            read(line,*) string, eref%atom_in_bnd_E(1:2, bnd_no)
            do i = 1, 2
                call find_atom_idx(0, ref_res_no, eref%atom_in_bnd_E(i, bnd_no), &
                                   eref%atm_in_bnd_E(i, bnd_no), atom_error_mode)
            end do
            ! to consider two types of +N (std(NH1,NH3), pro(N,NP))
            if (eref%atom_in_bnd_E(2, bnd_no) == '+N  ') then
                atom_in_bnd_E_save(1:2) = eref%atom_in_bnd_E(1:2, bnd_no)
                atm_in_bnd_E_save(1:2) = eref%atm_in_bnd_E(1:2, bnd_no)
                bnd_no = bnd_no + 1
                eref%atom_in_bnd_E(1:2, bnd_no) = atom_in_bnd_E_save(1:2)
                eref%atom_in_bnd_E(2, bnd_no) = '+NP '
                eref%atm_in_bnd_E(1:2, bnd_no) = atm_in_bnd_E_save(1:2)
            end if

        else if (string == 'DIH' .or. string == 'IMP' .or. string == 'ANG') then
            pbond_found = .false.
            ! read angle
            ang_no = ang_no + 1
            if (string == 'DIH') then
                imax = 4
                ang_type = 1
            else if (string == 'IMP') then
                imax = 4
                ang_type = 3
            else if (string == 'ANG') then
                imax = 3
                ang_type = 2
            end if
            eref%ang_E_type(ang_no) = ang_type
            eref%n_atm_in_ang(ang_no) = imax
            read(line,*) string, eref%atom_in_ang_E(1:imax, ang_no)
            do i = 1, imax
                call find_atom_idx(0, ref_res_no, eref%atom_in_ang_E(i, ang_no), &
                                   eref%atm_in_ang_E(i, ang_no), atom_error_mode)
                if (eref%atom_in_ang_E(i, ang_no) == '+N  ' .or. &
                    eref%atom_in_ang_E(i, ang_no) == '-CA ') then
                    pbond_found = .true.
                end if
            end do
            ! to consider two types of +N (std(NH1,NH3), pro(N,NP))
            ! or three types of -CA (std(CT1), pro(CP1), gly(CT2)
            if (pbond_found) then
                !save to copy
                atom_in_ang_E_save(1:imax) = eref%atom_in_ang_E(1:imax, ang_no)
                atm_in_ang_E_save(1:imax) = eref%atm_in_ang_E(1:imax, ang_no)
                ang_type_save = ang_type
                imax_save = imax
                !add for +N of proline
                if ((ang_type == 2 .or. ang_type == 3) .and. &
                    eref%atom_in_ang_E(3, ang_no) == '+N  ') then
                    ang_no = ang_no + 1
                    eref%atom_in_ang_E(1:imax, ang_no) = atom_in_ang_E_save(1:imax)
                    eref%atom_in_ang_E(3, ang_no) = '+NP '
                    eref%atm_in_ang_E(1:imax, ang_no) = atm_in_ang_E_save(1:imax)
                    eref%ang_E_type(ang_no) = ang_type_save
                    eref%n_atm_in_ang(ang_no) = imax_save
                else if (ang_type == 1 .and. eref%atom_in_ang_E(4, ang_no) == '+N  ') then
                    ang_no = ang_no + 1
                    eref%atom_in_ang_E(1:imax, ang_no) = atom_in_ang_E_save(1:imax)
                    eref%atom_in_ang_E(4, ang_no) = '+NP '
                    eref%atm_in_ang_E(1:imax, ang_no) = atm_in_ang_E_save(1:imax)
                    eref%ang_E_type(ang_no) = ang_type_save
                    eref%n_atm_in_ang(ang_no) = imax_save
                else if (ang_type == 1 .and. eref%atom_in_ang_E(1, ang_no) == '-CA ') then
                    !add for -CA of pro
                    ang_no = ang_no + 1
                    eref%atom_in_ang_E(1:imax, ang_no) = atom_in_ang_E_save(1:imax)
                    eref%atom_in_ang_E(1, ang_no) = '-CAP'
                    eref%atm_in_ang_E(1:imax, ang_no) = atm_in_ang_E_save(1:imax)
                    eref%ang_E_type(ang_no) = ang_type_save
                    eref%n_atm_in_ang(ang_no) = imax_save
                    !add for -CA of gly
                    ang_no = ang_no + 1
                    eref%atom_in_ang_E(1:imax, ang_no) = atom_in_ang_E_save(1:imax)
                    eref%atom_in_ang_E(1, ang_no) = '-CAG'
                    eref%atm_in_ang_E(1:imax, ang_no) = atm_in_ang_E_save(1:imax)
                    eref%ang_E_type(ang_no) = ang_type_save
                    eref%n_atm_in_ang(ang_no) = imax_save
                end if
            end if
           
        else if (line(1:4) == 'LINK') then
            ! initialize LINK structure
            read(line,*) string, res_name
            num_ref_res = num_ref_res + 1
            i_ref_res = num_ref_res
            ref_res_no = i_ref_res
            ref_res(ref_res_no)%res_name = res_name
           
            atm_no = 0
            ang_no = 0
            bnd_no = 0
           
            do
                read(f_unit,"(A100)") line
                string = line(1:4)
                if (string == 'END') then
                    eref%n_ang_E = ang_no
                    eref%n_bnd_E = bnd_no
                    ref_res(ref_res_no)%n_atm = atm_no
                 
                    ! set-up angle parameters
                    do i_ang = 1, eref%n_ang_E
                        do i = 1, eref%n_atm_in_ang(i_ang)
                            atom(i) = eref%atom_in_ang_E(i, i_ang)
                            call find_atom_idx(0, ref_res_no, atom(i), idx, atom_error_mode)
                            atm(i) = eref%atom_cls(idx)
                        end do
                        call find_ang_prm(atm(:), eref%ang_E_type(i_ang), eref%ang_prm(i_ang), status)
                    end do
                 
                    ! set-up bond parameters
                    do i_bnd = 1, eref%n_bnd_E
                        do i = 1, 2
                            atom(i) = eref%atom_in_bnd_E(i, i_bnd)
                            call find_atom_idx(0, ref_res_no, atom(i), idx, atom_error_mode)
                            atm(i) = eref%atom_cls(idx)
                        end do
                        call find_bnd_prm(atm(1:2), eref%bnd_prm(i_bnd), status)
                    end do
                 
                    ref_res_eng(ref_res_no) = eref
                    exit
                 
                else if (string == 'ATM') then
                    atm_no = atm_no + 1
                    read(line,*) string, atom_name, atom_cls
                    ref_res(ref_res_no)%atom_name(atm_no) = atom_name
                    ! update number of atoms
                    ref_res(ref_res_no)%n_atm = atm_no
                    eref%atom_cls(atm_no) = atom_cls
                 
                else if (string == 'BND') then
                    bnd_no = bnd_no + 1
                    read(line,*) string, eref%atom_in_bnd_E(1:2, bnd_no)
                    do i = 1, 2
                        call find_atom_idx(0, ref_res_no, eref%atom_in_bnd_E(i, bnd_no), &
                                           eref%atm_in_bnd_E(i, bnd_no), atom_error_mode)
                    end do
                    
                else if (string == 'DIH' .or. string == 'IMP' .or. string == 'ANG') then
                    ang_no = ang_no + 1
                    if (string == 'DIH') then
                        imax = 4
                        ang_type = 1
                    else if (string == 'IMP') then
                        imax = 4
                        ang_type = 3
                    else if (string == 'ANG') then
                        imax = 3
                        ang_type = 2
                    end if
                    eref%ang_E_type(ang_no) = ang_type
                    eref%n_atm_in_ang(ang_no) = imax
                    read(line,*) string, eref%atom_in_ang_E(1:imax, ang_no)
                    do i = 1, imax
                        call find_atom_idx(0, ref_res_no, eref%atom_in_ang_E(i, ang_no), &
                                           eref%atm_in_ang_E(i, ang_no), atom_error_mode)
                    end do
                end if
            end do
        end if
    end do
             
    close(f_unit)
end do

if (num_ref_res > max_ref_res) then
    write(error_line, "(A,I6,I6)") 'num_ref_res > max_ref_res', num_ref_res, max_ref_res
    call terminate_with_error(error_line)
end if

write(log_msg,"(A,I5)") '  Num_ref_res = ', num_ref_res
call log_p(log_msg, me=me, level=30)

end subroutine read_topology_eng
!-------------------------------------------------------------------------------
subroutine read_topology_mol2(ligand)
!-------------------------------------------------------------------------------
! Read ligand mol2 file, and setup ligand topology.
!-------------------------------------------------------------------------------
type(ref_lig_type), intent(out) :: ligand(max_lig)
integer :: i_mol2

do i_mol2 = 1, n_mol2_top
    write(log_msg,"(A,A)") '  Reading topology from mol2: '//trim(infile_mol2_topo(1, i_mol2))
    call log_p(log_msg,me=me,level=20)

    call read_mol2(infile_mol2_topo(1,i_mol2), ligand(i_mol2))
    if (top_type == 'polarh') then
        call log_p('Merging nonpolar hydrogen ligand atoms for top_type polarh.',level=20)
        call merge_nonpolar_hydrogens(ligand(i_mol2))
    end if
    call ligand_to_topology(ligand(i_mol2))
    ligand(i_mol2)%lig_name = infile_mol2_topo(2, i_mol2)(1:3)
end do

end subroutine read_topology_mol2
!-------------------------------------------------------------------------------
subroutine read_ligand_parameters(ref_lig)
!-------------------------------------------------------------------------------
! Read parameter file for ligand, then setup parameters of ref_lig.
!-------------------------------------------------------------------------------
type(ref_lig_type), intent(inout) :: ref_lig(max_lig)
type(eng_para_type) :: lig_para
character(len=6) :: mol2_types(max_lig_atom), mol2_type
integer :: n_mol2_types
integer :: i, i_atm, i_mol2
integer :: atm_no, bnd_no, ang_no
logical :: is_used(4), is_new

! Reading Ligand Parameter File
call read_parameters(infile_parameter_ligand, lig_para)

! Extract mol2_type used in ref_ligs
n_mol2_types = 0
do i_mol2 = 1, n_mol2_top
    do i_atm = 1, ref_lig(i_mol2)%n_atm
        write(mol2_type, '(A,A)') 'L', ref_lig(i_mol2)%mol2_type(i_atm)(1:5)
        ref_lig(i_mol2)%atom_type(i_atm) = mol2_type

        is_new = .true.
        do i = 1, n_mol2_types
            if (mol2_type == mol2_types(i)) then
                is_new = .false.
                exit
            end if
        end do
        if (is_new) then
            n_mol2_types = n_mol2_types + 1
            mol2_types(n_mol2_types) = mol2_type
        end if
    end do
end do

! Extract ATOM parameters used in ref_ligs
atm_no = eng_para%n_atom_cls
do i = 1, lig_para%n_atom_cls
    is_used(1) = mol2_type_is_used(lig_para%atom_cls(i), mol2_types, n_mol2_types)
    if (is_used(1)) then
        atm_no = atm_no + 1
        eng_para%atom_cls(atm_no)   = lig_para%atom_cls(i)
        eng_para%mass(atm_no)       = lig_para%mass(i)
        eng_para%coord_type(atm_no) = lig_para%coord_type(i)
        eng_para%LJ_para(:,atm_no)  = lig_para%LJ_para(:,i)
    end if
end do
eng_para%n_atom_cls = atm_no

bnd_no = eng_para%n_bnd
do i = 1, lig_para%n_bnd
    is_used(1) = mol2_type_is_used(lig_para%atm_in_bnd(1,i),mol2_types,n_mol2_types)
    is_used(2) = mol2_type_is_used(lig_para%atm_in_bnd(2,i),mol2_types,n_mol2_types)
    if (is_used(1) .and. is_used(2)) then
        bnd_no = bnd_no + 1
        eng_para%atm_in_bnd(1:2, bnd_no) = lig_para%atm_in_bnd(1:2, i)
        eng_para%bnd_para(1:2, bnd_no) = lig_para%bnd_para(1:2, i)
    end if
end do
eng_para%n_bnd = bnd_no

ang_no = eng_para%n_ang
do i = 1, lig_para%n_ang
    is_used(1) = mol2_type_is_used(lig_para%atm_in_ang(1,i),mol2_types,n_mol2_types)
    is_used(2) = mol2_type_is_used(lig_para%atm_in_ang(2,i),mol2_types,n_mol2_types)
    is_used(3) = mol2_type_is_used(lig_para%atm_in_ang(3,i),mol2_types,n_mol2_types)
    if (lig_para%ang_type(i) /= 2) then
        is_used(4) = mol2_type_is_used(lig_para%atm_in_ang(4,i),mol2_types,n_mol2_types)
    else
        is_used(4) = .true.
    end if
    if (is_used(1) .and. is_used(2) .and. is_used(3) .and. is_used(4)) then
        ang_no = ang_no + 1
        eng_para%atm_in_ang(:,ang_no) = lig_para%atm_in_ang(:,i)
        eng_para%n_ang_fold(ang_no)   = lig_para%n_ang_fold(i)
        eng_para%ang_type(ang_no)     = lig_para%ang_type(i)
        eng_para%ang_para(:,lig_para%n_ang_fold(i),ang_no) = &
                                        lig_para%ang_para(:,lig_para%n_ang_fold(i),i)
    end if
end do
eng_para%n_ang = ang_no

end subroutine read_ligand_parameters
!-------------------------------------------------------------------------------
function mol2_type_is_used(atm_type, ligand_mol2_types, n_mol2_types)
!-------------------------------------------------------------------------------
! Check atm_type is in ligand or not.
! Variables
!   atm_type: mol2 type in parameter file to check whether this mol2 type is in
!             ligand_mol2_types or not.
!   ligand_mol2_types: mol2 types in ligand 
!   n_mol2_types: No. of mol2 type in ligand
!-------------------------------------------------------------------------------
logical :: mol2_type_is_used
character(len=6), intent(in) :: atm_type, ligand_mol2_types(n_mol2_types)
integer, intent(in) :: n_mol2_types
integer :: i

if (atm_type == 'X' .or. atm_type == 'LX') then
    mol2_type_is_used = .true.
    return
else
    mol2_type_is_used = .false.
end if

do i = 1, n_mol2_types
    if (atm_type == ligand_mol2_types(i)) then
        mol2_type_is_used = .true.
        return
    end if
end do

end function mol2_type_is_used
!-------------------------------------------------------------------------------
subroutine assign_mol2_ref_res_eng(ligand)
!-------------------------------------------------------------------------------
! assign parameters for mol2 type ligand molecule.
!-------------------------------------------------------------------------------
type(ref_lig_type), intent(inout) :: ligand(max_lig)
character(len=6) :: atmName, atom(4)
integer :: i_atm, i_cls, atm_cls, q_typ_no, i_mol2, i
logical :: status

q_typ_no = eng_para%n_q_typ !charge type

do i_mol2 = 1, n_mol2_top !for each ligand (n_mol2_top = number of ligands)
    do i_atm = 1, ligand(i_mol2)%n_atm ! atom index
        atmName = ligand(i_mol2)%atom_type(i_atm) ! ex: CA

        !find atm_cls of atom
        atm_cls = 0
        do i_cls = 1, eng_para%n_atom_cls
            if (eng_para%atom_cls(i_cls) == atmName) then
                atm_cls = i_cls
                exit
            end if
        end do

        ! exeception handling. if atmName is not exist in eng_para%atom_cls,
        ! throw error message
        if (atm_cls == 0) then
            write(log_msg, "(A,2X,A4,A)") 'WARNING: atom type not found in parameter file: ', atmName, ' --> LX'
            call log_p(log_msg, level=0, me=me)

            atmName = 'LX'
            ligand(i_mol2)%atom_type(i_atm) = atmName
            do i_cls = 1, eng_para%n_atom_cls
                if (eng_para%atom_cls(i_cls) == atmName) then
                    atm_cls = i_cls
                    exit
                end if
            end do
        end if
        
        !set atm_cls and charge type
        q_typ_no = q_typ_no + 1
        eng_para%charge(q_typ_no) = ligand(i_mol2)%charge(i_atm)

        ligand(i_mol2)%atm_cls(i_atm) = atm_cls    ! Index for eng_para
        ligand(i_mol2)%q_typ_no(i_atm) = q_typ_no
    end do
    
    !find index for bond parameter in eng_para_structure
    do i = 1, ligand(i_mol2)%n_bnd
        atom(1) = ligand(i_mol2)%atom_type(ligand(i_mol2)%bnd(2,i))
        atom(2) = ligand(i_mol2)%atom_type(ligand(i_mol2)%bnd(3,i))
        call find_bnd_prm(atom(1:2), ligand(i_mol2)%eng_para(i,1), status)
    end do

    !find index for angle parameter in eng_para_structure
    do i = 1, ligand(i_mol2)%n_ang
        atom(1) = ligand(i_mol2)%atom_type(ligand(i_mol2)%ang(2,i))
        atom(2) = ligand(i_mol2)%atom_type(ligand(i_mol2)%ang(3,i))
        atom(3) = ligand(i_mol2)%atom_type(ligand(i_mol2)%ang(4,i))
        call find_ang_prm(atom(1:3), 2, ligand(i_mol2)%eng_para(i,2), status)
    end do

    !find index for dihedral parameter in eng_para_structure
    do i = 1, ligand(i_mol2)%n_dih
        atom(1) = ligand(i_mol2)%atom_type(ligand(i_mol2)%dih(2,i))
        atom(2) = ligand(i_mol2)%atom_type(ligand(i_mol2)%dih(3,i))
        atom(3) = ligand(i_mol2)%atom_type(ligand(i_mol2)%dih(4,i))
        atom(4) = ligand(i_mol2)%atom_type(ligand(i_mol2)%dih(5,i))
        call find_ang_prm(atom(1:4), 1, ligand(i_mol2)%eng_para(i,3), status)
    end do
end do
eng_para%n_q_typ = q_typ_no

end subroutine assign_mol2_ref_res_eng
!-------------------------------------------------------------------------------
subroutine setup_ref_res_aux()
!-------------------------------------------------------------------------------
! Auxiliary things about ref_res are added here...
!-------------------------------------------------------------------------------
integer :: i_ref_res, i_atm
character(len=6) :: atom_name
character(len=4) :: res_name

! is_sc_atom indexing
do i_ref_res = 1, num_ref_res
    do i_atm = 1, ref_res(i_ref_res)%n_atm
        ref_res(i_ref_res)%is_sc_atom(i_atm) = .false.
        atom_name = ref_res(i_ref_res)%atom_name(i_atm)
        if (atom_name == 'N'   .or. atom_name == 'CA' .or. atom_name == 'C' .or. &
            atom_name == 'H'   .or. atom_name == 'CB' .or. atom_name == 'O' .or. &
            atom_name == 'H1'  .or. atom_name == 'H2' .or. atom_name == 'H3' .or. &
            atom_name == 'OXT' .or. atom_name == 'HA') cycle
        res_name = ref_res(i_ref_res)%res_name
        if (res_name == 'NACE' .or. res_name == 'CNME') cycle
        ref_res(i_ref_res)%is_sc_atom(i_atm) = .true.
    end do
    !
    call identify_t_ang_dependence(ref_res(i_ref_res))
end do

!call setup_ref_res_bb_angle()

end subroutine setup_ref_res_aux
!-------------------------------------------------------------------------------
subroutine get_aux_para(eng_para,para_type,i,j)
!-------------------------------------------------------------------------------
! Fill auxiliary parameters in eng_para
!-------------------------------------------------------------------------------
type(eng_para_type), intent(inout) :: eng_para
real(dp) :: para, sigma
integer :: i, ii, j
character(len=2) :: para_type
if (para_type == 'LJ') then
    para = eng_para%LJ_para(j,i)
    if (j == 1 .or. j == 3) then
        eng_para%aux_para(j,i,i) = para
        if (force_field_type == 'AMBER' .and. j==3) then
            eng_para%aux_para(j,i,i) = eng_para%aux_para(j,i,i)*eng_para%V14fac
        end if
    else if (j == 2 .or. j == 4) then
        eng_para%aux_para(j,i,i) = 4.0d0*para*para
    end if

    do ii = 1, eng_para%n_atom_cls
        if (j == 1 .or. j == 3) then ! Epsilon for regular(1) & 1-4(3)
            eng_para%aux_para(j,i,ii) = -sqrt(para*eng_para%LJ_para(j,ii))
            if (force_field_type == 'AMBER' .and. j == 3) then
                eng_para%aux_para(j,i,ii) = eng_para%aux_para(j,i,ii)*eng_para%V14fac
            end if
        else if (j == 2 .or. j == 4) then ! Sigma**2 for regular(2) & 1-4(4)
            sigma = para + eng_para%LJ_para(j,ii)
            eng_para%aux_para(j,i,ii) = sigma*sigma
        end if
        eng_para%aux_para(j,ii,i) = eng_para%aux_para(j,i,ii)
    end do
end if

end subroutine get_aux_para
!-------------------------------------------------------------------------------
subroutine setup_vdw_para()
!-------------------------------------------------------------------------------
! setup van der Waals parameter
!-------------------------------------------------------------------------------
integer :: i, j

! calc sigma_ij and eps_ij for LJ intrxn using the combination rule
! para%aux_para for coarse model may not filled, therefore LJ is not available
if (top_type /= 'coarse') then
    do i = 1, eng_para%n_atom_cls
        do j = 1, 4
            call get_aux_para(eng_para,'LJ',i,j) !get parameters from
                                                 !parameters_XXX.in file
        end do
    end do

    ! calc vdw radii sum of i,j types
    do i = 1, eng_para%n_atom_cls
        do j = 1, eng_para%n_atom_cls
            eng_para%vdwsum(i,j) = eng_para%LJ_para(2,i) + eng_para%LJ_para(2,j)
        end do
    end do
end if

end subroutine setup_vdw_para
!-------------------------------------------------------------------------------
subroutine read_parameters(infile_para, para)
!-------------------------------------------------------------------------------
! Read molecular mechanics energy parameters
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: infile_para
type(eng_para_type), intent(inout) :: para
character(len=len_fname) :: line, word(10)
character(len=10) :: atom
character(len=6) :: string
integer :: f_unit = 23, i, atm_no, atm_no2, ang_no, bnd_no, m
integer :: openstat
real(dp) :: par, b0
real(dp) :: Emin, Rmin, Emin14, Rmin14, k, n, phi0
real(dp) :: tmp1, tmp2
logical :: same_ang
character(len=6) :: atoms(4)
integer :: num_word

open(f_unit, file = infile_para, iostat = openstat)

if (openstat > 0) then
    call log_p('Terminate with error: Infile parameter file not found.')
    call terminate_with_error('Please check {data_directory}.')
end if

200 format(a100)

! read ATOM CLASS
do  
    read(f_unit, 200) line
    if (line == 'ATOM CLASS')exit
end do

atm_no = 0
do 
    read(f_unit, 200) line
    if (line == 'END') exit
    atm_no = atm_no + 1
    read(line, *) para%atom_cls(atm_no), para%mass(atm_no), para%coord_type(atm_no)
end do
para%n_atom_cls = atm_no

! read bond parameters
do  
    read(f_unit, 200) line
    if (line == 'BOND') exit
end do

bnd_no = 0
do 
    read(f_unit, 200) line
    if (line == 'END') exit
    bnd_no = bnd_no + 1
    read(line, *) para%atm_in_bnd(1:2,bnd_no), k, b0
    para%bnd_para(1:2,bnd_no)=(/ k, b0 /)
end do
para%n_bnd = bnd_no

! read bond angle parameters
do  
    read(f_unit, 200) line
    if (line == 'BOND ANGLE') exit
end do

ang_no = 0
do 
    read(f_unit, 200) line
    if (line == 'END') exit
    ang_no = ang_no + 1
    para%ang_type(ang_no) = 2
    read(line, *) para%atm_in_ang(1:3,ang_no), k, phi0
    para%n_ang_fold(ang_no) = 1
    para%ang_para(1:2,1,ang_no)=(/ k, deg2rad*phi0 /)
end do

! read dihedral angle parameters
do  
    read(f_unit, 200) line
    if (line == 'DIHEDRAL ANGLE') exit
end do

do 
    read(f_unit, 200) line
    if (line == 'END') exit
    read(line, *) atoms(1:4), m, k, phi0, n
    same_ang = .false.
    if (ang_no > 0) then
        same_ang = .true.
        do i = 1, 4
            if (atoms(i) /= para%atm_in_ang(i,ang_no)) then
                same_ang = .false.
                exit
            end if
        end do
    end if
    if (same_ang) then
        para%n_ang_fold(ang_no) = para%n_ang_fold(ang_no) + 1
    else
        ang_no = ang_no + 1 
        para%n_ang_fold(ang_no) = 1
        para%ang_type(ang_no) = 1
        para%atm_in_ang(1:4,ang_no) = atoms(1:4) 
    end if
    para%ang_para(1:3,para%n_ang_fold(ang_no),ang_no) &
        = (/ k/dble(m), deg2rad*phi0, n /)
end do

! read improper torsion angle parameters
do  
    read(f_unit, 200) line
    if (line == 'IMPROPER TORSION ANGLE') exit
end do

do 
    read(f_unit, 200) line
    if (line == 'END') exit
    ang_no = ang_no + 1
    para%ang_type(ang_no) = 3
    if (force_field_type == 'CHARMM') then
        read(line, *) para%atm_in_ang(1:4,ang_no), k, phi0
        para%ang_para(1:2,1,ang_no)=(/ k, deg2rad*phi0 /)
    else if (force_field_type == 'AMBER') then
        read(line, *) para%atm_in_ang(1:4,ang_no), k, phi0, n
        para%ang_para(1:3,1,ang_no)=(/ k, deg2rad*phi0, n /)
    end if
    para%n_ang_fold(ang_no) = 1
end do
para%n_ang = ang_no

! read LJ parameters
do  
    read(f_unit, 200) line
    if (line == 'VDW') exit
end do

atm_no = 0
do 
    read(f_unit, 200) line
    if (line == 'END') exit
    if (line(1:6) == 'E14FAC') then
        read(line,*) string, para%E14fac
        if (force_field_type == 'AMBER') para%E14fac = 1.0d0/para%E14fac
    else if (line(1:6)=='V14FAC') then
        read(line,*) string, para%V14fac
        if (force_field_type == 'AMBER') para%V14fac = 1.0d0/para%V14fac
    else
        if (trim(top_type) == 'coarse') then
            read(line,*) atoms(1:2), tmp1, tmp2
            call find_atom_cls(para, atoms(1), atm_no)
            call find_atom_cls(para, atoms(2), atm_no2)
            para%vdwsum(atm_no,atm_no2) = tmp1
            para%vdwsum(atm_no2,atm_no) = tmp1
        else
            read(line,*) atom 
            if (atom /= '!') then
                call find_atom_cls(para, atom, atm_no)
            
                if (force_field_type == 'AMBER') then
                    read(line,*) atom, Emin, Rmin
                    para%LJ_para(1:4,atm_no)=(/ Emin, Rmin, Emin, Rmin /)
                else if (trim(top_type) == 'allh_ch22') then
                    call parse_string(line, num_word, word)
                    if (num_word == 7) then
                        read(line,*) atom, par, Emin, Rmin, par, Emin14, Rmin14
                    else
                        read(line,*) atom, par, Emin, Rmin
                        Emin14 = Emin
                        Rmin14 = Rmin
                    end if
                    para%LJ_para(1:4,atm_no)=(/ Emin, Rmin, Emin14, Rmin14 /)
                else 
                    if (atom(1:1) == 'C' .and. atom(1:2) /= 'C0') then
                        read(line,*) atom, par, Emin, Rmin, par, Emin14, Rmin14
                    else
                        read(line,*) atom, par, Emin, Rmin
                        Emin14 = Emin
                        Rmin14 = Rmin
                    end if
                    para%LJ_para(1:4,atm_no)=(/ Emin, Rmin, Emin14, Rmin14 /)
                end if
            end if
        end if
    end if
end do

close(f_unit)
   
end subroutine read_parameters
!------------------------------------------------------------------------------
END MODULE IN_OUT_PARAMETERS
!-------------------------------------------------------------------------------
