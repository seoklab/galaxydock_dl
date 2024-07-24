!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/apps/mpi_utils.f90
!
! Description:
!
!-------------------------------------------------------------------------------
MODULE MPI_UTILS
!-------------------------------------------------------------------------------
use mpi
!
use globals
use logger, only : log_p, terminate_with_error
use energy_vars, only: energy_type, n_E_component_modeling, &
                       n_E_component_ppdock, n_E_component_ligdock

implicit none

public

!integer, parameter :: idint = 1105, idloop = 1729
integer, parameter :: mpi_tag_signal = 99
integer, parameter :: mpi_tag_R = 100
integer, parameter :: mpi_tag_E = 101
integer, parameter :: mpi_tag_str = 102
integer, parameter :: mpi_tag_n_str = 103

CONTAINS
!-------------------------------------------------------------------------------
subroutine initialize_mpi()
!-------------------------------------------------------------------------------
integer :: ierr
      
! start parallel processing
call mpi_init(ierr)
if (ierr /= 0) then
    call terminate_with_error(' Error. Could not initialize MPI.')
end if

! determine # of nodes and current node
call mpi_comm_rank(mpi_comm_world, me, ierr)
if (ierr /= 0) then
    call terminate_with_error(' Error. Could not determine the ranks of all processes.')
end if

call mpi_comm_size(mpi_comm_world, n_proc, ierr)

write(log_msg, "(A,I4)") "- Initializing Parallel Computing with n_Proc =", n_proc
call log_p(log_msg, me=me, level=10)

if (ierr /= 0) then
    call terminate_with_error(' Error. Could not determine the number of processes.')
end if

end subroutine initialize_mpi
!-------------------------------------------------------------------------------
subroutine finalize_mpi()
!-------------------------------------------------------------------------------
integer :: ierr

call MPI_BARRIER(mpi_comm_world, ierr)
if (ierr /= 0) then
    write(log_msg,"(A,I3)") ' Error in mpi_barrier; IERR : ', ierr
    call terminate_with_error(log_msg)
end if

call mpi_finalize(ierr)
if (ierr /= 0) then
    write(log_msg,"(A,I3)") 'Error in mpi_finalize; IERR : ', ierr
    call terminate_with_error(log_msg)
end if

end subroutine finalize_mpi
!-------------------------------------------------------------------------------
subroutine mpi_send_R(i_proc, R)
!-------------------------------------------------------------------------------
! Send R from me to i_proc
!-------------------------------------------------------------------------------
integer, intent(in) :: i_proc
real(dp), intent(in) :: R(:,:)
!
real(dp), allocatable :: R_send(:)
integer :: tna, ierr
integer :: status(MPI_STATUS_SIZE)

tna = tn%atom
allocate(R_send(tna*3))

R_send = reshape(R, (/tna*3/))

call MPI_SSEND(R_send, tna*3, MPI_DOUBLE_PRECISION, &
              i_proc, mpi_tag_R, MPI_COMM_WORLD, ierr)

write(log_msg, "(A,1x,I4,A,1x,I4,1x,A,I5)") '  send    R from', me, ' to', i_proc, ', ', tna
call log_p(log_msg, me=king, level=40)

deallocate(R_send)

end subroutine mpi_send_R
!-------------------------------------------------------------------------------
subroutine mpi_recv_R(i_proc, R)
!-------------------------------------------------------------------------------
! Receive R from i_proc to me
!-------------------------------------------------------------------------------
integer, intent(in) :: i_proc
real(dp), intent(out) :: R(:,:)
!
real(dp), allocatable :: R_recv(:)
integer :: tna, ierr
integer :: status(MPI_STATUS_SIZE)

tna = tn%atom
allocate(R_recv(tna*3))

call MPI_RECV(R_recv, tna*3, MPI_DOUBLE_PRECISION, &
              i_proc, mpi_tag_R, MPI_COMM_WORLD, status, ierr)

write(log_msg, "(A,1x,I4,A,1x,I4,1x,A,I4)") '  recv  R from', i_proc, ' to', me, ', ', tna
call log_p(log_msg, me=king, level=40)

R = reshape(R_recv, (/3, tna/))

deallocate(R_recv)

end subroutine mpi_recv_R
!-------------------------------------------------------------------------------
subroutine mpi_send_signal(i_proc, signal)
!-------------------------------------------------------------------------------
integer, intent(in) :: i_proc, signal
integer :: ierr

call MPI_SSEND(signal, 1, MPI_INTEGER, i_proc, mpi_tag_signal, MPI_COMM_WORLD, ierr)

write(log_msg, "(A,1x,I4,A,1x,I4)") '  send signal from', me, ' to', i_proc
call log_p(log_msg, me=king, level=40)

end subroutine mpi_send_signal
!-------------------------------------------------------------------------------
subroutine mpi_recv_signal(i_proc, signal)
!-------------------------------------------------------------------------------
integer, intent(in) :: i_proc
integer, intent(out) :: signal
integer :: ierr
integer :: status(MPI_STATUS_SIZE)

call MPI_RECV(signal, 1, MPI_INTEGER, i_proc, mpi_tag_signal, MPI_COMM_WORLD, status, ierr)

write(log_msg, "(A,1x,I4,A,1x,I4)") '  recv signal from', i_proc, ' to', me
call log_p(log_msg, me=king, level=40)

end subroutine mpi_recv_signal
!-------------------------------------------------------------------------------
subroutine mpi_send_energy(i_proc, n_traj, energy)
!-------------------------------------------------------------------------------
integer, intent(in) :: i_proc, n_traj
type(energy_type), intent(in) :: energy(:)
integer :: n_E_component, n_E, ierr, i_traj
real(dp), allocatable :: E_send(:)

n_E_component = n_traj*&
                (n_E_component_modeling + n_E_component_ppdock + n_E_component_ligdock + 4)

allocate(E_send(n_E_component))
!
n_E = 1
do i_traj = 1, n_traj
    E_send(n_E) = energy(i_traj)%total
    n_E = n_E + 1
    E_send(n_E:n_E+n_E_component_modeling) = energy(i_traj)%modeling(0:n_E_component_modeling)
    n_E = n_E+n_E_component_modeling+1
    E_send(n_E:n_E+n_E_component_ppdock)   = energy(i_traj)%ppdock(0:n_E_component_ppdock)
    n_E = n_E+n_E_component_ppdock+1
    E_send(n_E:n_E+n_E_component_ligdock)  = energy(i_traj)%ligdock(0:n_E_component_ligdock)
    n_E = n_E+n_E_component_ligdock+1
end do

call MPI_SSEND(E_send, n_E_component, MPI_DOUBLE_PRECISION, &
               i_proc, mpi_tag_E, MPI_COMM_WORLD, ierr)

write(log_msg, "(A,1x,I4,A,1x,I4)") '  send energy from', me, ' to', i_proc
call log_p(log_msg, me=king, level=40)

deallocate(E_send)

end subroutine mpi_send_energy
!-------------------------------------------------------------------------------
subroutine mpi_recv_energy(i_proc, n_traj, energy)
!-------------------------------------------------------------------------------
integer, intent(in) :: i_proc, n_traj
type(energy_type), intent(out) :: energy(:)
integer :: n_E_component, n_E, ierr, i_traj
integer :: status(MPI_STATUS_SIZE)
real(dp), allocatable :: E_recv(:)

n_E_component = n_traj*&
                (n_E_component_modeling + n_E_component_ppdock + n_E_component_ligdock + 4)

allocate(E_recv(n_E_component))

call MPI_RECV(E_recv, n_E_component, MPI_DOUBLE_PRECISION, &
              i_proc, mpi_tag_E, MPI_COMM_WORLD, status, ierr)

n_E = 1
do i_traj = 1, n_traj
    energy(i_traj)%total = E_recv(n_E)
    n_E = n_E + 1
    energy(i_traj)%modeling(0:n_E_component_modeling) = E_recv(n_E:n_E+n_E_component_modeling)
    n_E = n_E+n_E_component_modeling+1
    energy(i_traj)%ppdock(0:n_E_component_ppdock)     = E_recv(n_E:n_E+n_E_component_ppdock)  
    n_E = n_E+n_E_component_ppdock+1
    energy(i_traj)%ligdock(0:n_E_component_ligdock)   = E_recv(n_E:n_E+n_E_component_ligdock)
    n_E = n_E+n_E_component_ligdock+1
end do

write(log_msg, "(A,1x,I4,A,1x,I4)") '  recv energy from', me, ' to', i_proc
call log_p(log_msg, me=king, level=40)

deallocate(E_recv)

end subroutine mpi_recv_energy
!-------------------------------------------------------------------------------
subroutine mpi_send_string(i_proc, l_string, n_string, string)
!-------------------------------------------------------------------------------
integer, intent(in) :: i_proc, l_string, n_string
character(len=l_string), intent(in) :: string(n_string)
!
integer :: n_char, i,j, ierr
character(len=1), allocatable :: char_send(:)

call MPI_SSEND(n_string, 1, MPI_INTEGER, i_proc, mpi_tag_n_str, &
               MPI_COMM_WORLD, ierr)

n_char = l_string*n_string
allocate(char_send(n_char))
!
do i = 1, n_string
    do j = 1, l_string
        char_send((i-1)*l_string+j) = string(i)(j:j)
    end do
end do
!
call MPI_SSEND(char_send, n_char, MPI_CHARACTER, i_proc, mpi_tag_str, &
               MPI_COMM_WORLD, ierr)
!
deallocate(char_send)

end subroutine mpi_send_string
!-------------------------------------------------------------------------------
subroutine mpi_recv_string(i_proc, l_string, n_string, string)
!-------------------------------------------------------------------------------
integer, intent(in) :: i_proc, l_string
integer, intent(out) :: n_string
character(len=l_string), intent(out) :: string(:)
!
integer :: n_char, i,j, ierr
integer :: status(MPI_STATUS_SIZE)
character(len=1), allocatable :: char_recv(:)

call MPI_RECV(n_string, 1, MPI_INTEGER, i_proc, mpi_tag_n_str, &
              MPI_COMM_WORLD, status, ierr)

n_char = l_string*n_string
allocate(char_recv(n_char))
!
call MPI_RECV(char_recv, n_char, MPI_CHARACTER, i_proc, mpi_tag_str, &
               MPI_COMM_WORLD, status, ierr)
!
do i = 1, n_string
    do j = 1, l_string
        string(i)(j:j) = char_recv((i-1)*l_string+j)
    end do
end do
!
deallocate(char_recv)

end subroutine mpi_recv_string
!-------------------------------------------------------------------------------
END MODULE MPI_UTILS
!-------------------------------------------------------------------------------
