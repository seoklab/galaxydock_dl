!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/core/io/logger.f90
!
! Description:
!  This module controls all the standard GALAXY logs.
!
!-------------------------------------------------------------------------------
MODULE LOGGER
!-------------------------------------------------------------------------------
use globals

implicit none
public

integer :: log_unit                           ! unit of logs

CONTAINS
!-------------------------------------------------------------------------------
subroutine open_logfile(filename, f_unit)
!-------------------------------------------------------------------------------
character(len=len_fname), intent(in) :: filename
integer, intent(in) :: f_unit

if (me == king) then
    open(f_unit, file=filename, status='replace')
end if

end subroutine open_logfile
!-------------------------------------------------------------------------------
subroutine close_logfile(f_unit)
!-------------------------------------------------------------------------------
integer, intent(in) :: f_unit

close(f_unit)

end subroutine close_logfile
!-------------------------------------------------------------------------------
subroutine log_p(msg, me, level) 
!-------------------------------------------------------------------------------
! Report log on the log_unit.
! Processor to be printed out can be controled by optional variable 'me'.
! Log will be printed out only if global_print_level > 'level', the optional variable.
! Unless level specified, log will be printed out always.
!-------------------------------------------------------------------------------
integer, optional :: me, level
character(len=*) :: msg
integer :: temp_me, temp_level
  
if (present(me)) then
    temp_me = me
else 
    temp_me = king
end if

if (present(level)) then
    temp_level = level
else 
    temp_level = -10!print_level_global
end if

if (temp_level > print_level_global) return

if (temp_me == king) then
    write(log_unit,"(A)") trim(msg)
end if

end subroutine log_p
!-------------------------------------------------------------------------------
subroutine log_divider(me, level)
!-------------------------------------------------------------------------------
! Print "-----------------------------------------------------------"
! if current print level is lower than global print level
!-------------------------------------------------------------------------------
integer, optional :: me, level
integer :: temp_me, temp_level

if (present(me)) then
    temp_me = me
else 
    temp_me = king
end if

if (present(level)) then
    temp_level = level
else 
    temp_level = print_level_global
end if

if (temp_level > print_level_global ) return

if (temp_me == king) then
    write(log_unit,"(A)") divider
end if

end subroutine log_divider
!-------------------------------------------------------------------------------
subroutine log_thick_divider(me,level)
!-------------------------------------------------------------------------------
! Print "==========================================================="
! if current print level is lower than global print level
!-------------------------------------------------------------------------------
integer, optional :: me,level
integer :: temp_me, temp_level

if (present(me)) then
    temp_me = me
else 
    temp_me = king
end if

if (present(level)) then
    temp_level = level
else 
    temp_level = print_level_global
end if

if (temp_level > print_level_global ) return

if (temp_me == king) then
    write(log_unit,"(A)") thick_divider
end if

end subroutine log_thick_divider
!-------------------------------------------------------------------------------
subroutine terminate_with_error(error_message)
!-------------------------------------------------------------------------------
character(len=*), intent(in) :: error_message

write(*,"(A)") error_message
stop

end subroutine terminate_with_error
!-------------------------------------------------------------------------------
subroutine my_timer(ttime)
!-------------------------------------------------------------------------------
! Return current time (in seconds)
!-------------------------------------------------------------------------------
real(dp), intent(out) :: ttime
integer(8) :: count, count_rate

call system_clock(count, count_rate)
ttime = real(count, kind = dp) / real(count_rate, kind = dp)

end subroutine my_timer
!-------------------------------------------------------------------------------
subroutine finish_timer(time1)
!-------------------------------------------------------------------------------
! Count current time (time2) and print out time elapsed from time1
!-------------------------------------------------------------------------------
real(dp), intent(in) :: time1
real(dp) :: time2, sec
integer :: day, mn, hour

call my_timer(time2)
call log_divider(me=me)
write(log_msg,"(A,F10.3,A)") 'Total User Time = ', time2 - time1, ' sec'
call log_p(log_msg, me=me)
call convert_time(time2-time1, sec, mn, hour, day)
write(log_msg,"(10X,3(I2,A),F6.3,A)") day, ' dy ', hour, ' hr ', mn, ' min ', sec, ' sec '
call log_p(log_msg, me=me)
call log_divider(me=me)

end subroutine finish_timer
!-------------------------------------------------------------------------------
subroutine convert_time(time, sec, min, hour, day)
!-------------------------------------------------------------------------------
! Convert time (in sec.) into higher-order units
!-------------------------------------------------------------------------------
real(dp), intent(in) :: time
real(dp), intent(out) :: sec
integer, intent(out) :: min,hour,day

sec = time
min = 0; hour = 0; day = 0
if (sec > 60.0d0) then
    min = int(sec/60.0d0)
    sec = sec - min*60
    if (min > 60) then
        hour = int(min/60)
        min = min - hour*60
        if (hour > 24) then
            day = int(hour/24)
            hour = hour - day*24
        end if
    end if
end if

end subroutine convert_time
!-------------------------------------------------------------------------------
subroutine report_time_diff(prefix, time1, time2)
!-------------------------------------------------------------------------------
! Print out on log_unit for time difference (time2 - time1)
!-----------------------------------------------------------------------
character(len=len_fname), intent(in) :: prefix
real(dp), intent(in) :: time1, time2
real(dp) :: time_diff
real(dp) :: sec
integer :: day, mn, hour
character(len=len_fname) :: message

time_diff = time2-time1
call convert_time(time_diff, sec, mn, hour, day)
write(message,"(A,3(I2,A),F4.1,A)") trim(prefix), day, 'D ', hour, 'H ', mn, 'M ', sec, 'S'
call log_p(divider)
call log_p(message)

end subroutine report_time_diff
!-----------------------------------------------------------------------
subroutine report_time_diff2(message, time1, time2)
!-----------------------------------------------------------------------
! Returns message reporting time difference (time2 - time1)
!-----------------------------------------------------------------------
character(len=*), intent(out) :: message
real(dp), intent(in) :: time1, time2
real(dp) :: time_diff
real(dp) :: sec
integer :: day, mn, hour

time_diff = time2-time1
call convert_time(time_diff, sec, mn, hour, day)
write(message, "(3(I2,A),F4.1,A)") day, 'D ', hour, 'H ', mn, 'M ', sec, 'S'

end subroutine report_time_diff2
!-------------------------------------------------------------------------------
END MODULE LOGGER
!-------------------------------------------------------------------------------
