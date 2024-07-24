!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/apps/ligand_docking/ligdock_input.f90
!
! Description: read input file to set parameters related to run ligand docking 
!
!-------------------------------------------------------------------------------
MODULE LIGDOCK_INPUT
!-------------------------------------------------------------------------------
use globals
use logger, only: terminate_with_error
use string, only: parse_string
!
use ligdock_csa_vars
use ligdock_csa_runner, only: set_default_ligdock_csa

implicit none
save
public

CONTAINS
!-------------------------------------------------------------------------------
subroutine read_ligdock_input(input_file, include_init_conf)
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: input_file
integer :: f_unit, ioerr, num_word
character(len=len_fname) :: line, keyword, word(10), string
logical, intent(inout) :: include_init_conf

if (n_usc > 0) then
    call set_default_ligdock_csa(5)
else
    call set_default_ligdock_csa(4)
end if

include_init_conf = .false.
ML_rerank = .true.

f_unit = 20

open(f_unit, file=trim(input_file), iostat=ioerr)
if (ioerr > 0) then
    write(log_msg, '(A,A)') 'ERROR: Cannot open the file: ', input_file
    call terminate_with_error(log_msg)
end if

do
    read(f_unit, '(A120)', iostat=ioerr) line
    if (ioerr < 0) exit

    if (line(1:1) == '!' .or. line(1:1) == '#') cycle

    call parse_string(line, num_word, word)
    keyword = word(1)

    ! Read variables related to first bank generation method
    if (keyword == 'first_bank') then
        if (word(2) == 'rand') then
            lig_conf = 10
        else if (word(2) == 'Voro') then
            lig_conf = 20
        else if (word(2) == 'FFT') then
            lig_conf = 30
        end if
    else if (keyword == 'n_gen_conf') then
        read(word(2),*) n_gen_conf
    else if (keyword == 'n_voro_per_conf') then
        read(word(2),*) n_voro_per_conf
    else if (keyword == 'voro_prefix') then
        voro_prefix = word(2)
    else if (keyword == 'include_initial_conf') then
        if (word(2) == 'yes') then
            include_init_conf = .true.
        end if

    ! Read csa options
    else if (keyword == 'csa_n_csa_cycle') then
        read(line,*) string, n_csa_cycle
    else if (keyword == 'csa_bank') then
        read(line,*) string, max_bank, n_bank_add
    else if (keyword == 'csa_n_opr_x') then
        read(line,*) string, n_opr_1,n_opr_2,n_opr_3,n_opr_4,n_opr_5,n_opr_6,n_opr_7,n_opr_8
    else if (keyword == 'csa_seed') then
        read(line,*) string, n_seed_cycle, n_seed, &
                     max_opt_cycle, cut_pos_seed_size
    else if (keyword == 'csa_is') then
        read(line,*) string, min_tor_select, max_tor_select_ratio
    else if (keyword == 'csa_nran') then
        read(line,*) string, irr
    else if (keyword == 'csa_D_factor') then
        read(line,*) string,factor_init_D_cut, factor_min_D_cut, n_opt_to_D_min
    else if (keyword == 'e0max') then
        read(line,*) string, e0max
    else if (keyword == 'e1max') then
        read(line,*) string, e1max
    else if (keyword == 'max_trial') then
        read(line,*) string, max_trial
    else if (keyword == 'mutation_cut') then
        read(line,*) string, mutation_cut
    else if (keyword == 'rep_conf_seed_cut') then
        read(line,*) string, rep_conf_seed_cut
    else if (keyword == 'print_bank_evol') then 
        if (word(2) == 'yes') then
            print_bank_evol = .true.
        else if(word(2) == 'no') then
            print_bank_evol = .false.
        end if
    else if(keyword == 'print_curr_bank') then
        if(word(2) == 'yes') then
            print_curr_bank = .true.
        end if
    else if(keyword == 'print_bank') then
        if(word(2) == 'yes') then
            print_bank = .true.
        end if

    ! Read output related options
    else if (keyword == 'ligdock_prefix') then
        ligdock_prefix = word(2)

    else if (keyword == 'rerank_mode') then
        if(word(2) == 'no') then
            ML_rerank = .false.
        end if
    end if
end do

close(f_unit)

call check_ligdock_input_consistency()

end subroutine read_ligdock_input
!-------------------------------------------------------------------------------
subroutine check_ligdock_input_consistency()
!-------------------------------------------------------------------------------

! Check consistency between lig_conf and max_trial
if (lig_conf == 20) then
    if (n_gen_conf == 0 .or. n_voro_per_conf == 0) then
        call terminate_with_error('Please check n_gen_conf and n_voro_per_conf')
    end if
    max_trial = n_gen_conf * n_voro_per_conf
end if

n_csa_iter = max_bank/n_bank_add
n_new_conf = n_seed * (n_opr_1 + n_opr_2 + n_opr_3 + n_opr_4 + n_opr_5 + n_opr_6 + n_opr_7 + n_opr_8)

n_gen_str = 0
total_nft = 0

end subroutine check_ligdock_input_consistency
!-------------------------------------------------------------------------------
END MODULE LIGDOCK_INPUT
!-------------------------------------------------------------------------------
