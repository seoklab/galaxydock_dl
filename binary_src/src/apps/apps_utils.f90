!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/apps/apps_utils.f90
!
! Description:
!
!-------------------------------------------------------------------------------
MODULE APPS_UTILS
!-------------------------------------------------------------------------------
use globals
use logger, only: log_p, terminate_with_error, log_thick_divider, &
                  my_timer, finish_timer
use ran,    only: initialize_random
!
use in_out_parameters, only: initialize_parameters, finalize_parameters
use in_out_vars,  only: infile_pdb, multiple_task, multiple_models, &
                        infile_pdblist
use in_out_input, only: read_global_input, read_pdblist
use in_out_structure, only: n_models_in_pdb, finalize_structure
!
use energy,       only: finalize_energy
use energy_input, only: read_energy_input
!
#ifdef MPI
    use mpi_utils, only: initialize_mpi, finalize_mpi
#endif

implicit none
public

integer :: n_task

type infile_pdb_type
    character(len=len_fname) :: infile_pdb
    integer :: model_no
end type infile_pdb_type
type(infile_pdb_type) :: infile_pdb_s(max_pdbfile)

real(dp), save :: time_start

CONTAINS
!===============================================================================
subroutine initialize_apps(infile_input, use_energy)
!-------------------------------------------------------------------------------
character(len=len_fname), intent(out) :: infile_input
logical, intent(in), optional :: use_energy
!
character(len=len_fname) :: cmd
integer :: n_argc
!
character(len=len_fname) :: fnamelist(max_pdbfile)
integer :: i,j, n_pdb, n_model, model_no(max_pdbfile)

call my_timer(time_start)

! Read arguments given by user
n_argc = iargc()
call getarg(0, cmd)
if (n_argc >= 1) then
    call getarg(1, infile_input)
else 
    write(log_msg,"(A,2X,A)") trim(cmd),"[USER INPUT]"
    call terminate_with_error(log_msg)
end if

#ifdef MPI
    call initialize_mpi()
#endif

! Read the input file
call read_global_input(infile_input)
if ((.not. present(use_energy)) .or. (use_energy)) then
    call read_energy_input(infile_input)
end if

n_task = 0
! Read pdblist if it is multiple_task
if (multiple_task .and. trim(infile_pdblist) /= '') then
    call log_thick_divider(level=20, me=me)
    call log_p('Reading multiple input pdbs.', me=me)
    call read_pdblist(infile_pdblist, fnamelist, n_pdb)
    !
    do i = 1, n_pdb
        if (multiple_models) then
            call n_models_in_pdb(fnamelist(i), n_model, model_no)
            do j = 1, n_model
                n_task = n_task + 1
                infile_pdb_s(n_task)%infile_pdb = fnamelist(i)
                infile_pdb_s(n_task)%model_no = model_no(j)
            end do
        else
            n_task = n_task + 1
            infile_pdb_s(n_task)%infile_pdb = fnamelist(i)
            infile_pdb_s(n_task)%model_no = 1
        end if
    end do
    infile_pdb = infile_pdb_s(1)%infile_pdb
else if (trim(infile_pdb) /= '') then
    if (multiple_models) then
        call n_models_in_pdb(infile_pdb, n_model, model_no)
        if (n_model /= 1) then
            multiple_task = .true.
        end if
        do j = 1, n_model
            n_task = n_task + 1
            infile_pdb_s(n_task)%infile_pdb = infile_pdb
            infile_pdb_s(n_task)%model_no = model_no(j)
        end do
    else
        n_task = n_task + 1
        infile_pdb_s(1)%infile_pdb = infile_pdb
        infile_pdb_s(n_task)%model_no = 1
    end if
end if

call initialize_parameters()
call initialize_random()

end subroutine initialize_apps
!-------------------------------------------------------------------------------
subroutine finalize_apps(use_energy)
!-------------------------------------------------------------------------------
logical, intent(in), optional :: use_energy

if ((.not. present(use_energy)) .or. (use_energy)) then
    call finalize_energy()
end if
call finalize_parameters()

#ifdef MPI
    call finalize_mpi()
#endif

call finalize_structure()

call finish_timer(time_start)

end subroutine finalize_apps
!-------------------------------------------------------------------------------
subroutine distribute_tasks(n_proc, n_task, n_decoy, n_task_proc, task_sch)
!-------------------------------------------------------------------------------
! Setting up the schedule for multiple task with multiple processors
!  - Total number of structures to be generated is n_task*n_decoy
!  - They are distributed as (ie, n_proc=4, n_task=3, n_decoy=3)
!   - proc_0: (1,1), (2,2), (3,3)
!   - proc_1: (1,2), (2,3)
!   - proc_2: (1,3), (3,1)
!   - proc_3: (2,1), (3,2)
!   - The output, task_sch(:,:,:), stores these information as
!   - task_sch(1:2, i_conf, i_proc) = (/k_task, k_decoy/)
!   - use it as
!    do i_task_proc = 1, n_task_proc
!        init_conf = task_sch(1, i_task_proc, me+1)
!        ith_decoy = task_sch(2, i_task_proc, me+1)
!    end do
!-------------------------------------------------------------------------------
integer, intent(in) :: n_proc
integer, intent(in) :: n_task, n_decoy, n_task_proc
integer, intent(inout) :: task_sch(:,:,:)
integer :: i_task, i_proc, k_task, k_decoy

task_sch(:,:,:) = 0

k_task = 1
k_decoy = 0
do i_task = 1, n_task_proc
    do i_proc = 1, n_proc
        k_decoy = k_decoy + 1
        task_sch(1, i_task, i_proc) = k_task
        task_sch(2, i_task, i_proc) = k_decoy
        !
        if (k_decoy == n_decoy) then
            k_decoy = 0
            k_task = k_task + 1
        end if
        if (k_task > n_task) return
    end do
end do

end subroutine distribute_tasks
!-------------------------------------------------------------------------------
END MODULE APPS_UTILS
!-------------------------------------------------------------------------------
