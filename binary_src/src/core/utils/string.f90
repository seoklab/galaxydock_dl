!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/core/utils/string.f90
!
! Description:
!  This module contains string-related subroutines.
!
!-------------------------------------------------------------------------------
MODULE STRING
!-------------------------------------------------------------------------------
use globals
use logger, only: terminate_with_error

implicit none
public

CONTAINS
!-------------------------------------------------------------------------------
subroutine parse_string(long_string, num_word, word)
!-------------------------------------------------------------------------------
! Scan string, get num of "word", and save in array word
!-------------------------------------------------------------------------------
character(len=*), intent(in) :: long_string
integer, intent(out) :: num_word
character(len=len_fname), intent(out) :: word(:)
character(len=len_fname) :: string
integer :: i, n, k
character(len=1) :: prev_char, cur_char

n = len(long_string)
if (n > len_fname) then
    write(log_msg,"(A,I6)") 'Error. The length of the long string is greater than maximum: ', len_fname
    call terminate_with_error(log_msg)
end if

num_word = 0
prev_char = ' '
string = ' '
k = 0

! go over each character
do i = 1, n
    cur_char = long_string(i:i)
    if (len(trim(cur_char)) /= 0) then
        if (len(trim(prev_char)) == 0) then
            ! start of a new word
            num_word = num_word + 1
            if (num_word > 1) then
                ! save previous word and start new
                word(num_word-1) = string
                string = ' '
                k = 0
            end if
        end if
      
        if (len(trim(cur_char)) /= 0) then ! add 
            k = k + 1
            string(k:k) = cur_char
        end if
    end if
    prev_char = cur_char
end do

! save the last word
if (num_word > 0) then
    word(num_word) = string
end if

end subroutine parse_string
!-------------------------------------------------------------------------------
subroutine parse_longstring(long_string, num_word, word, linelen)
!-------------------------------------------------------------------------------
! Same with parse_string, but can treat longer string for which 
! line length is given by linelen
!-------------------------------------------------------------------------------
integer, intent(in) :: linelen
character(len=linelen), intent(in) :: long_string
integer, intent(out) :: num_word
character(len=len_fname), intent(out) :: word(:)
character(len=len_fname) :: string
integer :: i, n, k
character(len=1) :: prev_char, cur_char

n = len(long_string)

num_word = 0
prev_char = ' '
string = ' '
k = 0

! go over each character
do i = 1, n
    cur_char = long_string(i:i)
    if (len(trim(cur_char)) /= 0) then
        if (len(trim(prev_char)) == 0) then
            ! start of a new word
            num_word = num_word + 1
            if (num_word > 1) then
                ! save previous word and start new
                word(num_word-1) = string
                string = ' '
                k = 0
            end if
        end if
        if (len(trim(cur_char)) /= 0) then
            ! add 
            k = k + 1
            string(k:k) = cur_char
        end if
    end if
    prev_char = cur_char
end do

! save the last word
if (num_word > 0) then
    word(num_word) = string
end if

end subroutine parse_longstring
!-------------------------------------------------------------------------------
subroutine strip_string(long_string, n_char, char_range)
!-------------------------------------------------------------------------------
! Same with strip in Python
!-------------------------------------------------------------------------------
character(len=*), intent(in) :: long_string
integer, intent(out) :: n_char, char_range(2)
integer :: i, n

n = len(long_string)

do i = n, 1, -1
    if (len(trim(long_string(i:i))) /= 0) exit
end do
char_range(2) = i

do i = 1, char_range(2)
    if (len(trim(long_string(i:i))) /= 0) exit
end do
char_range(1) = i

n_char = char_range(2)-char_range(1)+1

end subroutine strip_string
!-------------------------------------------------------------------------------
subroutine split_string(long_string, splitter, num_word, word)
!-------------------------------------------------------------------------------
! Same with parse_string, but splits using a character not spaces.
! Same with split in Python
!-------------------------------------------------------------------------------
character(len=*), intent(in) :: long_string
character(len=*), intent(in) :: splitter
integer, intent(out) :: num_word
character(len=len_fname), intent(out) :: word(:)
integer :: i, n, m, k

n = len(long_string)
m = len(splitter)-1
if (n > len_fname) then
    write(log_msg,"(A,I6)") 'Error. The length of the long string is greater than maximum: ', len_fname
    call terminate_with_error(log_msg)
end if

num_word = 1
k = 0

word = ''
do i = 1, n-m
    if ((long_string(i:i+m)) /= splitter(1:(m+1))) then
        k = k + 1
        word(num_word)(k:k) = long_string(i:i)
    else
        k = 0
        num_word = num_word + 1
    end if
end do

end subroutine split_string
!-------------------------------------------------------------------------------
subroutine num2str(num, str)
!-------------------------------------------------------------------------------
! Convert integer into length-4 string
! Unoccupied order of number is filled by '0'
!-------------------------------------------------------------------------------
integer, intent(in) :: num
character(len=4), intent(out) :: str
integer :: i1, i2, i3, i4

if (num > 9999) then
    write(log_msg,"(A,I6)") 'Error in num2str. input number must be less than 10000, but it is', num
    call terminate_with_error(log_msg)
end if

i1 = num/1000
i2 = (num - 1000*i1)/100
i3 = (num - 1000*i1 - 100*i2)/10
i4 = num - 1000*i1 - 100*i2 - i3*10

str = char(i1 + ichar('0'))//char(i2 + ichar('0'))//char(i3 + ichar('0'))//char(i4 + ichar('0'))

end subroutine num2str
!-------------------------------------------------------------------------------
subroutine num2str_999(num, str)
!-------------------------------------------------------------------------------
! Convert integer into length-3 string
! Unoccupied order of number is filled by '0'
!-------------------------------------------------------------------------------
integer, intent(in) :: num
character(len=3), intent(out) :: str
integer :: i1, i2, i3

if (num > 999) then
   write(log_msg,"(A,I6)") 'Error in num2str. input number must be less than 1000, but it is', num
   call terminate_with_error(log_msg)
end if

i1 = num/100
i2 = (num - 100*i1)/10
i3 = num - 100*i1 - 10*i2

str = char(i1 + ichar('0'))//char(i2 + ichar('0'))//char(i3 + ichar('0'))

end subroutine num2str_999
!-------------------------------------------------------------------------------
END MODULE STRING
!-------------------------------------------------------------------------------
