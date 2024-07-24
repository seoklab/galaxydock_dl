!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/energy/energy_input.f90
!
! Description:
!  This file defines reading inputs related to energy calculation.
!-------------------------------------------------------------------------------
MODULE ENERGY_INPUT
!-------------------------------------------------------------------------------
use globals
use logger
use string, only : parse_string
use ramachandran, only: rama_mode
use energy_vars
use energy_weight

implicit none
save
private

public :: read_energy_input

CONTAINS
!-------------------------------------------------------------------------------
subroutine read_energy_input(file_name)
!-------------------------------------------------------------------------------
! Read modeling options from user input file 
! For detailed explanation for various options, refer to documents :)
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: file_name
!
integer, parameter :: f_unit = 20
integer :: ioerror
!
integer :: num_word
character(len=len_fname) :: line, word(11), keyword
!
character(len=len_fname) :: weight_type
logical :: dfire_lib_defined(2)
  
dfire_lib_defined(:) = .false.

! Set default
call set_default_energy_param()

open(f_unit, file = trim(file_name))

do
    read(f_unit, "(A120)", iostat = ioerror) line
    if (ioerror < 0) exit
    if (line(1:1) == '!' .or. line(1:1) == '#') cycle

    call parse_string(line, num_word, word)
    keyword = word(1)
    if (keyword == 'energy_print_level') then
        read(word(2), "(I10)") energy_print_level
        
    ! Pairlist update
    else if (keyword == 'Ron') then
        read(word(2), "(F15.5)") Ron
    else if (keyword == 'Roff') then
        read(word(2), "(F15.5)") Roff
        
    ! Energy general
    else if (keyword == 'weight_type') then
        weight_type = word(2)
        call setup_typical_weight(weight_type)
    else if (keyword == 'report_pairwise') then
        if (word(2) == 'yes') then
            report_pairwise = .true.
        else
            report_pairwise = .false.
        end if
        report_pairwise = .true.
    else if (keyword == 'energy_on_full') then
        if (word(2) == 'yes') then
            energy_on_full = .true.
        else
            energy_on_full = .false.
        end if

    ! Molecular Mechanics
    else if (keyword == 'use_bond') then
        if (word(2) == 'yes') then
            use_bond = .true.
        else if (word(2) == 'no') then
            use_bond = .false.
        end if
    else if (keyword == 'vdw_r_scale') then
        read(word(2),"(F15.5)") vdw_r_scale
    else if (keyword == 'trclash_scale') then
        read(word(2),"(F15.5)") trclash_scale
    else if (keyword == 'scclash_scale') then
        read(word(2),"(F15.5)") scclash_scale
    else if (keyword == 'vdw_weight') then
        read(word(2),"(F15.5)") vdw_w
    else if (keyword == 'elec_weight') then
        read(word(2),"(F15.5)") elec_w
        if (elec_w > small_real) then
            use_elec = .true.
        else
            use_elec = .false.
        end if
    else if (keyword == 'vdw_soft_k') then
        read(word(2),"(F15.5)") vdw_soft_k
    else if (keyword == 'vdw_type') then
        vdw_type = word(2)
        if (trim(vdw_type) == 'LJ') then
            use_LJ = .true.
            use_softsphere = .false.
        else if (trim(vdw_type) == 'soft') then
            use_LJ = .false.
            use_softsphere = .true.
        else if (trim(vdw_type) == 'scwrl4') then
            use_vdw_scwrl4 = .true.
            use_LJ = .false.
            use_softsphere = .false.
        end if
    else if (keyword == 'soften_short') then
        if (word(2) == 'yes') then
            soften_short = .true.
        else if (word(2) == 'no') then
            soften_short = .false.
        end if
    else if (keyword == 'LJ_att_w') then
        read(word(2),"(F15.5)") LJ_att_w
    else if (keyword == 'LJ_rep_w') then
        read(word(2),"(F15.5)") LJ_rep_w
    else if (keyword == 'rdie') then
        if (word(2) == 'yes') then
            rdie = .true.
        else if (word(2) == 'no') then
            rdie = .false.
        end if
    else if (keyword == 'die_const_pro') then
        read(word(2),"(F15.5)") die_const_pro
    else if (keyword == 'die_const_slv') then
        read(word(2),"(F15.5)") die_const_slv
    else if (keyword == 'use_shift_function') then
        if (word(2) == 'no') then
            truncate_by_shift = .false.
        endif
        
    else if (keyword == 'use_rattle') then
        if (word(2) == 'hydrogen') then
            rattle_type = 1
        else if (word(2) == 'all') then
            rattle_type = 2
        else if (word(2) == 'none') then
            rattle_type = -1
        else
            rattle_type = 0
        end if
        
    ! Solvation/SASA
    else if (keyword == 'solv_weight') then
        read(word(2),"(F15.5)") solv_w
        if (solv_w > small_real) then
            use_solv = .true.
        else
            use_solv = .false.
        end if
    else if (keyword == 'solv_type') then
        if (trim(word(2)) == 'facts') then
            solv_type = SOLV_FACTS
        else if (trim(word(2)) == 'factsmem') then
            solv_type = SOLV_FACTS
            use_factsmem = .true.
        else if (trim(word(2)) == 'eef1') then
            solv_type = SOLV_EEF1
        end if
    else if (keyword == 'SA_weight') then
        read(word(2),"(F15.5)") SA_w
        if(SA_w > small_real) then
            use_SA = .true.
        else
            use_SA = .false.
        end if
        if ((use_SA) .and. (tn%recres > 0)) use_SApp = .true.
    else if (keyword == 'SA_type') then
        if (trim(word(2)) == 'facts') then
            SA_type = SOLV_FACTS
        else if (trim(word(2)) == 'factsmem') then
            SA_type = SOLV_FACTS
            use_factsmem = .true.
        else if (trim(word(2)) == 'hasel') then
            SA_type = SOLV_HASEL
        else if (trim(word(2)) == 'haber') then
            SA_type = SOLV_HABER
        end if
    else if (keyword == 'ASP_type') then
        read(word(2),"(I10)") ASP_type
    else if (keyword == 'probe_radius') then
        read(word(2),"(F15.5)") r_probe
        
    ! Rotamer_Score
    else if (keyword == 'rotamer_weight') then
        read(word(2),"(F15.5)") rot_w
        if (rot_w > small_real) then
            use_rotamer_score = .true.
        else
            use_rotamer_score = .false.
        end if
        
    ! Restraint
    else if (keyword == 'rsr_weight') then
        read(word(2),"(F15.5)") rsr_w
        if (rsr_w > small_real) then
            use_rsr = .true.
        else
            use_rsr = .false.
        end if
    else if (keyword == 'rsr_file') then
        rsr_file = word(2)
        use_rsr = .true.
    else if (keyword == 'rsr_pdb') then
        rsr_pdb = word(2)
    else if (keyword == 'update_rsr_pdblist') then
        infile_update_rsr_pdblist = word(2)
        use_update_rsr = .true.
        if (num_word > 2) then
            read(word(3),'(F10.4)') update_rsr_dcut
        end if
        if (num_word > 3) then
            read(word(4),'(F10.4)') update_rsr_sig
        end if
        if (num_word > 4) then
            read(word(5),'(F10.4)') update_lig_rsr_sig
        end if
    else if (keyword == 'update_rsr_self') then
        if (word(2) == 'yes') then
            use_update_rsr = .true.
        end if
    else if (keyword == 'update_NO_rsr') then
        if (word(2) == 'yes') then
            update_NO_rsr = .true.
        end if
    else if (keyword == 'update_lig_rsr') then
        if (word(2) == 'yes') then
            update_lig_rsr = .true.
        end if
    else if (keyword == 'use_meld') then
        if (word(2) == 'yes') then
            use_meld = .true.
        end if
        
    ! Distogram
    else if (keyword == 'distogram_weight') then
        read(word(2),"(F15.5)") distogram_w
        if (distogram_w > small_real) then
            use_distogram = .true.
        else
            use_distogram = .false.
        end if
    else if (keyword == 'distogram_file') then
        distogram_file = word(2)
        use_distogram = .true.
    else if (keyword == 'distogram_pdb') then
        distogram_pdb = word(2)

    ! Backbone torsion angle
    else if (keyword == 'bb_torsion_weight') then
        read(word(2),"(F15.5)") bb_torsion_w
        if (bb_torsion_w > small_real) then
            use_bb_torsion = .true.
        else
            use_bb_torsion = .false.
        end if
    else if (keyword == 'bb_torsion_file') then
        bb_torsion_file = word(2)
        use_bb_torsion = .true.

    ! Machine Learning score
    else if (keyword == 'ml_weight') then
        read(word(2),"(F15.5)") ml_w
        if (ml_w > small_real) then
            use_ml = .true.
        else
            use_ml = .false.
        end if

    ! DFIRE energy related terms
    else if (keyword == 'dfire_score_file') then
        dfire_score_file = trim(data_dir)//word(2)
        dfire_lib_defined(1) = .true.
    else if (keyword == 'dfire_atype_file') then
        dfire_atype_file = trim(data_dir)//word(2)
        dfire_lib_defined(2) = .true.
    else if (keyword == 'dfire_weight') then
        read(word(2),"(F15.5)") dfire_w
        if (dfire_w > small_real) then
            use_ddfire = .true.
        else
            use_ddfire = .false.
        end if
    else if (keyword == 'ddfire_add_scale') then
        read(word(2),"(F15.5)") ddfire_add_scale
    ! Knowledge_GB
    else if (keyword == 'kgb_weight') then
        read(word(2),"(F15.5)") kgb_w
        if (kgb_w > small_real) then
            use_kgb = .true.
        else
            use_kgb = .false.
        end if
    else if (keyword == 'kgb_water_file') then
        kgb_water_file = word(2)
    else if (keyword == 'use_water_kgb') then
        if (word(2) == 'yes') then
            use_water_kgb = .true.
        end if
    else if (keyword == 'wkgb_add_scale') then
        read(word(2),"(F15.5)") wkgb_w
    ! GOAP
    else if (keyword == 'goap_weight') then
        read(word(2),"(F15.5)") goap_w
        if (goap_w > small_real) then
            use_goap = .true.
        else
            use_goap = .false.
        end if
        
    ! Hbond energy related terms
    else if (keyword == 'hbond_weight') then
        read(word(2),"(F15.5)") hbond_w
        if (hbond_w > small_real) then
            use_hbond = .true.
        else
            use_hbond = .false.
        end if
    else if (keyword == 'hbond_SA') then
        if (word(2) == 'yes') then
            use_SAHbond = .true.
        end if
    else if (keyword == 'coop_hbond') then
        if (word(2) == 'yes') then
            use_coHbond = .true.
        end if
        
    ! Ramachandran prefenrece
    else if (keyword == 'rama_weight') then
        read(word(2),"(F15.5)") rama_w
        if (rama_w > small_real) then
            use_rama_score = .true.
        else
            use_rama_score = .false.
        end if
    else if (keyword == 'rama_mode') then
        rama_mode = word(2)

    ! Conservation Score
    else if (keyword == 'conserve_weight') then
        read(word(2),"(F15.5)") conserve_w
        if (conserve_w > small_real) then
            use_conserve = .true.
        end if
    else if (keyword == 'blastfile') then
        blastfile = word(2)
        blastfile_defined = .true.

    ! InterEVScore
    else if (keyword == 'interev_weight') then
        read(word(2),"(F15.5)") interev_w
        if (interev_w > small_real) then
            use_interev = .true.
        end if
        
    ! AutoDock
    else if (keyword == 'infile_atdk_prm') then
        infile_atdk_prm = word(2)
    else if (keyword == 'infile_atdk3_sol_prm') then
        infile_atdk3_sol_prm = word(2)
        
    ! Related docking grid
    else if (keyword == 'grid_box_cntr') then 
        read(word(2:4),*) dock_grid_info%grid_cntr(:)
        use_input_cntr = .true.
    else if (keyword == 'grid_n_elem') then
        read(word(2:4),*) dock_grid_info%n_elem(:)
    else if (keyword == 'grid_width') then
        read(word(2),*) dock_grid_info%grid_width
    else if (keyword == 'soften_dock_E') then
        if (word(2) == 'yes') then
            soften_dock_E = .true.
        end if
    end if
end do

close(f_unit)

call check_energy_input()
call check_energy_input_consistency()

end subroutine read_energy_input
!-------------------------------------------------------------------------------
subroutine check_energy_input()
!-------------------------------------------------------------------------------
use_modeling_E = (use_softsphere .or. use_LJ .or. use_vdw_scwrl4 .or. use_elec .or. use_bond) .or. &
                 (use_solv .or. use_SA) .or. &
                 (use_ddfire .or. use_goap .or. use_kgb) .or. &
                 (use_hbond) .or. (use_rotamer_score) .or. (use_rama_score)  .or. &
                 (use_rsr)

use_ppdock_E   = use_conserve .or. use_interev

use_ligdock_E  = (use_atdk3 .or. use_atdk4) .or. (use_drugscore) .or. (use_Xscore)

end subroutine check_energy_input
!-------------------------------------------------------------------------------
subroutine check_energy_input_consistency()
!-------------------------------------------------------------------------------
if (use_kgb) then
    if (use_ddfire) then
        write(log_msg, '(A)') " ERROR: KGB energy is activated, but dDFIRE is also activated."
        call terminate_with_error(log_msg)
    end if
    !
    if ((.not. use_solv) .or. (solv_type /= SOLV_FACTS)) then
        write(log_msg, '(A)') " ERROR: KGB energy is activated, but FACTS solvation is turned of"
        call terminate_with_error(log_msg)
    end if
end if

end subroutine check_energy_input_consistency
!-------------------------------------------------------------------------------
END MODULE ENERGY_INPUT
!-------------------------------------------------------------------------------
