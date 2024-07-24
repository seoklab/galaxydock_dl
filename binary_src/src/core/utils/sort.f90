!-------------------------------------------------------------------------------
! Copyright (C) 2015, Lab of Computational Biology and Biomolecular Engineering
!
! File: $GALAXY/src/core/utils/sort.f90
!
! Description:
!  This module contains data sorting related subroutines.
!
!-------------------------------------------------------------------------------
MODULE SORT
!-------------------------------------------------------------------------------
use globals, only: dp

implicit none
public

CONTAINS
!-------------------------------------------------------------------------------
subroutine sort2(n, list, key)
!-------------------------------------------------------------------------------
!  Takes an input list of reals and sorts it
!  into ascending order using the Heapsort algorithm;
!  it also returns a key into the original ordering
!-------------------------------------------------------------------------------
integer, intent(in) :: n
real(dp), intent(inout) :: list(n)
integer, intent(inout) :: key(n)
integer i, j, k, index, keys
real(dp) :: lists

! initialize index into the original ordering
do i = 1, n
    key(i) = i
end do

! perform the heapsort of the input list
k = n/2 + 1
index = n
do while (n > 1)
    if (k > 1) then
        k = k - 1
        lists = list(k)
        keys = key(k)
    else
        lists = list(index)
        keys = key(index)
        list(index) = list(1)
        key(index) = key(1)
        index = index - 1
        if (index  <= 1) then
            list(1) = lists
            key(1) = keys
            return
        end if
    end if
    i = k
    j = k + k
    do while (j  <= index)
        if (j < index) then
            if (list(j) < list(j+1))  j = j + 1
        end if
        if (lists < list(j)) then
            list(i) = list(j)
            key(i) = key(j)
            i = j
            j = j + j
        else
            j = index + 1
        end if
    end do
    list(i) = lists
    key(i) = keys
end do

end subroutine sort2
!-------------------------------------------------------------------------------
subroutine sort1(n,list,key)
!-------------------------------------------------------------------------------
! Similar to sort1, but integer sorting
!-------------------------------------------------------------------------------
integer, intent(in) :: n
integer, intent(inout) :: list(n)
integer, intent(inout) :: key(n)
integer i, j, k, index, keys
integer :: lists

! initialize index into the original ordering
do i = 1, n
    key(i) = i
end do

! perform the heapsort of the input list
k = n/2 + 1
index = n
do while (n > 1)
    if (k > 1) then
        k = k - 1
        lists = list(k)
        keys = key(k)
    else
        lists = list(index)
        keys = key(index)
        list(index) = list(1)
        key(index) = key(1)
        index = index - 1
        if (index  <= 1) then
            list(1) = lists
            key(1) = keys
            return
        end if
    end if
    i = k
    j = k + k
    do while (j  <= index)
        if (j < index) then
            if (list(j) < list(j+1))  j = j + 1
        end if
        if (lists < list(j)) then
            list(i) = list(j)
            key(i) = key(j)
            i = j
            j = j + j
        else
            j = index + 1
        end if
    end do
    list(i) = lists
    key(i) = keys
end do

end subroutine sort1
!-------------------------------------------------------------------------------
END MODULE SORT
!-------------------------------------------------------------------------------
