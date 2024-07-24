module xlogp_Util_m
    implicit none

contains

    !!----------------------------------------------------------------------------
    !> GetUnitNo will return a free FORTRAN unit number.
    !! This was orginally written by John Burkardt.
    !! @return integer ( kind = 4 ) unitNo
    !! If unitNo = 0, then no free FORTRAN unit could be found, although
    !! all 99 units were checked (except for units 5 and 6).
    !! Otherwise, unitNo is an integer between 1 and 99, representing a
    !! free FORTRAN unit.
    !!
    !! Note that GetUnitNo assumes that units 5 and 6 are special, and
    !!  will never return those values.
    !!
    !!----------------------------------------------------------------------------
    subroutine GetUnitNo(unitNo)
        implicit none
        integer (kind=4), intent(inout) :: unitNo
        integer (kind=4) :: i, ios
        logical :: lopen

        unitNo = 0

        do i=1, 99
            if (i /= 5 .AND. i /= 6 ) then
                inquire(unit=i, opened=lopen, iostat=ios)
                if (ios == 0) then
                    if (.not. lopen) then
                        unitNo = i
                        return
                    endif
                endif
            endif
        enddo

        return
    end subroutine GetUnitNo

    !!----------------------------------------------------------------------------
    !> TimeStamp prints the current YMDHMS date as a time stamp.
    !! E.g) 31 May 2001 9:45:54.872 AM
    !!----------------------------------------------------------------------------
    function GetTimeStamp() result(timeStr)
        implicit none
        character(len=100) :: timeStr
        character (len=8) ampm
        integer (kind=4) :: d, h, m, mm
        integer (kind=4) :: n, s, y
        integer (kind=4) :: values(8)

        character(len=9), parameter, dimension(12) :: month = (/ &
            'January  ', &
            'February ', &
            'March    ', &
            'April    ', &
            'May      ', &
            'June     ', &
            'July     ', &
            'August   ', &
            'September', &
            'October  ', &
            'November ', &
            'December ' /)

        call date_and_time ( values = values )
        y = values(1)
        m = values(2)
        d = values(3)
        h = values(5)
        n = values(6)
        s = values(7)
        mm = values(8)

        if (h < 12 ) then
            ampm = 'AM'
        else if(h==12) then
            if(n==0 .AND. s == 0) then
                ampm = 'Noon'
            else
                ampm = 'PM'
            endif
        else
            h = h -12
            if(h < 12) then
                ampm = 'PM'
            else if(h == 12) then
                if (n==0 .AND. s == 0) then
                    ampm = 'Midnight'
                else
                    ampm = 'AM'
                endif
            endif
        endif

        write(timeStr, '(i2,1x,a,1x,i4,2x,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)') &
            d, trim(month(m)), y, h, ':', s, '.', mm, trim(ampm)
        return
    end function GetTimeStamp

    !!----------------------------------------------------------------------------
    !> TimeStamp prints the current YMDHMS date as a time stamp.
    !! E.g) 31 May 2001 9:45:54.872 AM
    !!----------------------------------------------------------------------------
    subroutine GetTimeStamp2(ss, hh, dd, mm2, yy)
        implicit none
        integer (kind=4) :: ss, hh, dd, mm2, yy
        character(len=100) :: timeStr
        character (len=8) ampm
        integer (kind=4) :: d, h, m, mm
        integer (kind=4) :: n, s, y

        integer (kind=4) :: values(8)

        character(len=9), parameter, dimension(12) :: month = (/ &
            'January  ', &
            'February ', &
            'March    ', &
            'April    ', &
            'May      ', &
            'June     ', &
            'July     ', &
            'August   ', &
            'September', &
            'October  ', &
            'November ', &
            'December ' /)

        call date_and_time ( values = values )
        y = values(1)
        m = values(2)
        d = values(3)
        h = values(5)
        n = values(6)
        s = values(7)
        mm = values(8)

!VALUE(1):   The year
!VALUE(2):   The month
!VALUE(3):   The day of the month
!VAlUE(4):   Time difference with UTC in minutes
!VALUE(5):   The hour of the day
!VALUE(6):   The minutes of the hour
!VALUE(7):   The seconds of the minute
!VALUE(8):   The milliseconds of the second

        yy = values(1)
        dd = values(3)
        hh = values(5)
        mm2 = values(6)
        ss = values(7)

        if (h < 12 ) then
            ampm = 'AM'
        else if(h==12) then
            if(n==0 .AND. s == 0) then
                ampm = 'Noon'
            else
                ampm = 'PM'
            endif
        else
            h = h -12
            if(h < 12) then
                ampm = 'PM'
            else if(h == 12) then
                if (n==0 .AND. s == 0) then
                    ampm = 'Midnight'
                else
                    ampm = 'AM'
                endif
            endif
        endif

        write(timeStr, '(i2,1x,a,1x,i4,2x,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)') &
            d, trim(month(m)), y, h, ':', s, '.', mm, trim(ampm)
        return
    end subroutine GetTimeStamp2

end module xlogp_Util_m

module xlogp_Random_m
! module contains three functions
! ran1 returns a uniform random number between 0-1
! spread returns random number between min - max
! normal returns a normal distribution

    implicit none

    integer, parameter:: b8 = selected_real_kind(14)
    !real(b8), parameter :: pi = 3.141592653589793239_b8
    !integer gene_size

contains
    function ran1()  !returns random number between 0 - 1
        implicit none
        real(b8) ran1,x
        call random_number(x) ! built in fortran 90 random number function
        ran1=x
    end function ran1

    function spread(min,max)  !returns random number between min - max
        implicit none
        real(b8) spread
        real(b8) min,max
        spread=(max - min) * ran1() + min
    end function spread

    function normal(mean,sigma) !returns a normal distribution
        implicit none
        real(b8) normal,tmp
        real(b8) mean,sigma
        integer flag
        real(b8) fac,gsave,rsq,r1,r2
        save flag,gsave
        data flag /0/
        if (flag.eq.0) then
            rsq=2.0_b8
            do while(rsq.ge.1.0_b8.or.rsq.eq.0.0_b8) ! new from for do
                r1=2.0_b8*ran1()-1.0_b8
                r2=2.0_b8*ran1()-1.0_b8
                rsq=r1*r1+r2*r2
            enddo
            fac=sqrt(-2.0_b8*log(rsq)/rsq)
            gsave=r1*fac
            tmp=r2*fac
            flag=1
        else
            tmp=gsave
            flag=0
        endif
        normal=tmp*sigma+mean
        return
    end function normal

end module xlogp_Random_m

module xlogp_Strings_m

    private :: value_dr,value_sr,value_di,value_si
    private :: write_dr,write_sr,write_di,write_si
    private :: writeq_dr,writeq_sr,writeq_di,writeq_si

    ! Real kinds

    integer, parameter :: kr4 = selected_real_kind(6,37)       ! single precision real
    integer, parameter :: kr8 = selected_real_kind(15,307)     ! double precision real

    ! Integer kinds

    integer, parameter :: ki4 = selected_int_kind(9)           ! single precision integer
    integer, parameter :: ki8 = selected_int_kind(18)          ! double precision integer

    !Complex kinds

    integer, parameter :: kc4 = kr4                            ! single precision complex
    integer, parameter :: kc8 = kr8                            ! double precision complex

    !!----------------------------------------------------------------------------
    !> Generic operator for converting a number string to a
    !! number. Calling syntax is 'call value(numstring,number,ios)'
    !! where 'numstring' is a number string and 'number' is a
    !! real number or an integer (single or double precision).
    !!----------------------------------------------------------------------------
    interface value
        module procedure value_dr
        module procedure value_sr
        module procedure value_di
        module procedure value_si
    end interface


    !!----------------------------------------------------------------------------
    !> Generic  interface for writing a number to a string. The
    !! number is left justified in the string. The calling syntax
    !! is 'call writenum(number,string,format)' where 'number' is
    !! a real number or an integer, 'string' is a character string
    !! containing the result, and 'format' is the format desired,
    !! e.g., 'e15.6' or 'i5'.
    !!----------------------------------------------------------------------------
    interface writenum
        module procedure write_dr
        module procedure write_sr
        module procedure write_di
        module procedure write_si
    end interface

    interface writeq  ! Generic interface equating a name to a numerical value. The
                       ! calling syntax is 'call writeq(unit,name,value,format)' where
                       ! unit is the integer output unit number, 'name' is the variable
                       ! name, 'value' is the real or integer value of the variable,
                       ! and 'format' is the format of the value. The result written to
                       ! the output unit has the form <name> = <value>.
        module procedure writeq_dr
        module procedure writeq_sr
        module procedure writeq_di
        module procedure writeq_si
    end interface


!**********************************************************************

contains

    !**********************************************************************

    subroutine parse(str,delims,args,nargs)

        ! Parses the string 'str' into arguments args(1), ..., args(nargs) based on
        ! the delimiters contained in the string 'delims'. Preceding a delimiter in
        ! 'str' by a backslash (\) makes this particular instance not a delimiter.
        ! The integer output variable nargs contains the number of arguments found.

        character(len=*) :: str,delims
        character(len=len_trim(str)) :: strsav
        character(len=*),dimension(:) :: args

        strsav=str
        call compact(str)
        na=size(args)
        do i=1,na
            args(i)=' '
        end do
        nargs=0
        lenstr=len_trim(str)
        if(lenstr==0) return
        k=0

        do
            if(len_trim(str) == 0) exit
            nargs=nargs+1
            call split(str,delims,args(nargs))
            call removebksl(args(nargs))
        end do
        str=strsav

    end subroutine parse

    !**********************************************************************

    subroutine compact(str)

        ! Converts multiple spaces and tabs to single spaces; deletes control characters;
        ! removes initial spaces.

        character(len=*):: str
        character(len=1):: ch
        character(len=len_trim(str)):: outstr

        str=adjustl(str)
        lenstr=len_trim(str)
        outstr=' '
        isp=0
        k=0

        do i=1,lenstr
            ch=str(i:i)
            ich=iachar(ch)

            select case(ich)

                case(9,32)     ! space or tab character
                    if(isp==0) then
                        k=k+1
                        outstr(k:k)=' '
                    end if
                    isp=1

                case(33:)      ! not a space, quote, or control character
                    k=k+1
                    outstr(k:k)=ch
                    isp=0

            end select

        end do

        str=adjustl(outstr)

    end subroutine compact

    !**********************************************************************

    subroutine removesp(str)

        ! Removes spaces, tabs, and control characters in string str

        character(len=*):: str
        character(len=1):: ch
        character(len=len_trim(str))::outstr

        str=adjustl(str)
        lenstr=len_trim(str)
        outstr=' '
        k=0

        do i=1,lenstr
            ch=str(i:i)
            ich=iachar(ch)
            select case(ich)
                case(0:32)  ! space, tab, or control character
                    cycle
                case(33:)
                    k=k+1
                    outstr(k:k)=ch
            end select
        end do

        str=adjustl(outstr)

    end subroutine removesp

    !**********************************************************************

    subroutine value_dr(str,rnum,ios)

        ! Converts number string to a double precision real number

        character(len=*)::str
        real(kr8)::rnum
        integer :: ios

        ilen=len_trim(str)
        ipos=scan(str,'Ee')
        if(.not.is_digit(str(ilen:ilen)) .and. ipos/=0) then
            ios=3
            return
        end if
        read(str,*,iostat=ios) rnum

    end subroutine value_dr

    !**********************************************************************

    subroutine value_sr(str,rnum,ios)

        ! Converts number string to a single precision real number

        character(len=*)::str
        real(kr4) :: rnum
        real(kr8) :: rnumd

        call value_dr(str,rnumd,ios)
        if( abs(rnumd) > huge(rnum) ) then
            ios=15
            return
        end if
        if( abs(rnumd) < tiny(rnum) ) rnum=0.0_kr4
        rnum=rnumd

    end subroutine value_sr

    !**********************************************************************

    subroutine value_di(str,inum,ios)

        ! Converts number string to a double precision integer value

        character(len=*)::str
        integer(ki8) :: inum
        real(kr8) :: rnum

        call value_dr(str,rnum,ios)
        if(abs(rnum)>huge(inum)) then
            ios=15
            return
        end if
        inum=nint(rnum,ki8)

    end subroutine value_di

    !**********************************************************************

    subroutine value_si(str,inum,ios)

        ! Converts number string to a single precision integer value

        character(len=*)::str
        integer(ki4) :: inum
        real(kr8) :: rnum

        call value_dr(str,rnum,ios)
        if(abs(rnum)>huge(inum)) then
            ios=15
            return
        end if
        inum=nint(rnum,ki4)

    end subroutine value_si

    !**********************************************************************

    subroutine shiftstr(str,n)

        ! Shifts characters in in the string 'str' n positions (positive values
        ! denote a right shift and negative values denote a left shift). Characters
        ! that are shifted off the end are lost. Positions opened up by the shift
        ! are replaced by spaces.

        character(len=*):: str

        lenstr=len(str)
        nabs=iabs(n)
        if(nabs>=lenstr) then
            str=repeat(' ',lenstr)
            return
        end if
        if(n<0) str=str(nabs+1:)//repeat(' ',nabs)  ! shift left
        if(n>0) str=repeat(' ',nabs)//str(:lenstr-nabs)  ! shift right
        return

    end subroutine shiftstr

    !**********************************************************************

    subroutine insertstr(str,strins,loc)

        ! Inserts the string 'strins' into the string 'str' at position 'loc'.
        ! Characters in 'str' starting at position 'loc' are shifted right to
        ! make room for the inserted string. Trailing spaces of 'strins' are
        ! removed prior to insertion

        character(len=*):: str,strins
        character(len=len(str))::tempstr

        lenstrins=len_trim(strins)
        tempstr=str(loc:)
        call shiftstr(tempstr,lenstrins)
        tempstr(1:lenstrins)=strins(1:lenstrins)
        str(loc:)=tempstr
        return

    end subroutine insertstr

    !**********************************************************************

    subroutine delsubstr(str,substr)

        ! Deletes first occurrence of substring 'substr' from string 'str' and
        ! shifts characters left to fill hole. Trailing spaces or blanks are
        ! not considered part of 'substr'.

        character(len=*):: str,substr

        lensubstr=len_trim(substr)
        ipos=index(str,substr)
        if(ipos==0) return
        if(ipos == 1) then
            str=str(lensubstr+1:)
        else
            str=str(:ipos-1)//str(ipos+lensubstr:)
        end if
        return

    end subroutine delsubstr

    !**********************************************************************

    subroutine delall(str,substr)

        ! Deletes all occurrences of substring 'substr' from string 'str' and
        ! shifts characters left to fill holes.

        character(len=*):: str,substr

        lensubstr=len_trim(substr)
        do
            ipos=index(str,substr)
            if(ipos == 0) exit
            if(ipos == 1) then
                str=str(lensubstr+1:)
            else
                str=str(:ipos-1)//str(ipos+lensubstr:)
            end if
        end do
        return

    end subroutine delall

    !**********************************************************************

    function uppercase(str) result(ucstr)

        ! convert string to upper case

        character (len=*):: str
        character (len=len_trim(str)):: ucstr

        ilen=len_trim(str)
        ioffset=iachar('A')-iachar('a')
        iquote=0
        ucstr=str
        do i=1,ilen
            iav=iachar(str(i:i))
            if(iquote==0 .and. (iav==34 .or.iav==39)) then
                iquote=1
                iqc=iav
                cycle
            end if
            if(iquote==1 .and. iav==iqc) then
                iquote=0
                cycle
            end if
            if (iquote==1) cycle
            if(iav >= iachar('a') .and. iav <= iachar('z')) then
                ucstr(i:i)=achar(iav+ioffset)
            else
                ucstr(i:i)=str(i:i)
            end if
        end do
        return

    end function uppercase

    !**********************************************************************

    function lowercase(str) result(lcstr)

        ! convert string to lower case

        character (len=*):: str
        character (len=len_trim(str)):: lcstr

        ilen=len_trim(str)
        ioffset=iachar('A')-iachar('a')
        iquote=0
        lcstr=str
        do i=1,ilen
            iav=iachar(str(i:i))
            if(iquote==0 .and. (iav==34 .or.iav==39)) then
                iquote=1
                iqc=iav
                cycle
            end if
            if(iquote==1 .and. iav==iqc) then
                iquote=0
                cycle
            end if
            if (iquote==1) cycle
            if(iav >= iachar('A') .and. iav <= iachar('Z')) then
                lcstr(i:i)=achar(iav-ioffset)
            else
                lcstr(i:i)=str(i:i)
            end if
        end do
        return

    end function lowercase

    !**********************************************************************

    subroutine readline(nunitr,line,ios)

        ! Reads line from unit=nunitr, ignoring blank lines
        ! and deleting comments beginning with an exclamation point(!)

        character (len=*):: line

        do
            read(nunitr,'(a)', iostat=ios) line      ! read input line
            if(ios /= 0) return
            line=adjustl(line)
            ipos=index(line,'!')
            if(ipos == 1) cycle
            if(ipos /= 0) line=line(:ipos-1)
            if(len_trim(line) /= 0) exit
        end do
        return

    end subroutine readline

    !**********************************************************************

    subroutine match(str,ipos,imatch)

        ! Sets imatch to the position in string of the delimiter matching the delimiter
        ! in position ipos. Allowable delimiters are (), [], {}, <>.

        character(len=*) :: str
        character :: delim1,delim2,ch

        lenstr=len_trim(str)
        delim1=str(ipos:ipos)
        select case(delim1)
            case('(')
                idelim2=iachar(delim1)+1
                istart=ipos+1
                iend=lenstr
                inc=1
            case(')')
                idelim2=iachar(delim1)-1
                istart=ipos-1
                iend=1
                inc=-1
            case('[','{','<')
                idelim2=iachar(delim1)+2
                istart=ipos+1
                iend=lenstr
                inc=1
            case(']','}','>')
                idelim2=iachar(delim1)-2
                istart=ipos-1
                iend=1
                inc=-1
            case default
                write(*,*) delim1,' is not a valid delimiter'
                return
        end select
        if(istart < 1 .or. istart > lenstr) then
            write(*,*) delim1,' has no matching delimiter'
            return
        end if
        delim2=achar(idelim2) ! matching delimiter

        isum=1
        do i=istart,iend,inc
            ch=str(i:i)
            if(ch /= delim1 .and. ch /= delim2) cycle
            if(ch == delim1) isum=isum+1
            if(ch == delim2) isum=isum-1
            if(isum == 0) exit
        end do
        if(isum /= 0) then
            write(*,*) delim1,' has no matching delimiter'
            return
        end if
        imatch=i

        return

    end subroutine match

    !**********************************************************************

    subroutine write_dr(rnum,str,fmt)

        ! Writes double precision real number rnum to string str using format fmt

        real(kr8) :: rnum
        character(len=*) :: str,fmt
        character(len=80) :: formt

        formt='('//trim(fmt)//')'
        write(str,formt) rnum
        str=adjustl(str)

    end subroutine write_dr

    !***********************************************************************

    subroutine write_sr(rnum,str,fmt)

        ! Writes single precision real number rnum to string str using format fmt

        real(kr4) :: rnum
        character(len=*) :: str,fmt
        character(len=80) :: formt

        formt='('//trim(fmt)//')'
        write(str,formt) rnum
        str=adjustl(str)

    end subroutine write_sr

    !***********************************************************************

    subroutine write_di(inum,str,fmt)

        ! Writes double precision integer inum to string str using format fmt

        integer(ki8) :: inum
        character(len=*) :: str,fmt
        character(len=80) :: formt

        formt='('//trim(fmt)//')'
        write(str,formt) inum
        str=adjustl(str)

    end subroutine write_di

    !***********************************************************************

    subroutine write_si(inum,str,fmt)

        ! Writes single precision integer inum to string str using format fmt

        integer(ki4) :: inum
        character(len=*) :: str,fmt
        character(len=80) :: formt

        formt='('//trim(fmt)//')'
        write(str,formt) inum
        str=adjustl(str)

    end subroutine write_si

    !***********************************************************************

    subroutine trimzero(str)

        ! Deletes nonsignificant trailing zeroes from number string str. If number
        ! string ends in a decimal point, one trailing zero is added.

        character(len=*) :: str
        character :: ch
        character(len=10) :: exp

        ipos=scan(str,'eE')
        if(ipos>0) then
            exp=str(ipos:)
            str=str(1:ipos-1)
        endif
        lstr=len_trim(str)
        do i=lstr,1,-1
            ch=str(i:i)
            if(ch=='0') cycle
            if(ch=='.') then
                str=str(1:i)//'0'
                if(ipos>0) str=trim(str)//trim(exp)
                exit
            endif
            str=str(1:i)
            exit
        end do
        if(ipos>0) str=trim(str)//trim(exp)

    end subroutine trimzero

    !**********************************************************************

    subroutine writeq_dr(unit,namestr,value,fmt)

        ! Writes a string of the form <name> = value to unit

        real(kr8) :: value
        integer :: unit
        character(len=*) :: namestr,fmt
        character(len=32) :: tempstr

        call writenum(value,tempstr,fmt)
        call trimzero(tempstr)
        write(unit,*) trim(namestr)//' = '//trim(tempstr)

    end subroutine writeq_dr

    !**********************************************************************

    subroutine writeq_sr(unit,namestr,value,fmt)

        ! Writes a string of the form <name> = value to unit

        real(kr4) :: value
        integer :: unit
        character(len=*) :: namestr,fmt
        character(len=32) :: tempstr

        call writenum(value,tempstr,fmt)
        call trimzero(tempstr)
        write(unit,*) trim(namestr)//' = '//trim(tempstr)

    end subroutine writeq_sr

    !**********************************************************************

    subroutine writeq_di(unit,namestr,ivalue,fmt)

        ! Writes a string of the form <name> = ivalue to unit

        integer(ki8) :: ivalue
        integer :: unit
        character(len=*) :: namestr,fmt
        character(len=32) :: tempstr
        call writenum(ivalue,tempstr,fmt)
        call trimzero(tempstr)
        write(unit,*) trim(namestr)//' = '//trim(tempstr)

    end subroutine writeq_di

    !**********************************************************************

    subroutine writeq_si(unit,namestr,ivalue,fmt)

        ! Writes a string of the form <name> = ivalue to unit

        integer(ki4) :: ivalue
        integer :: unit
        character(len=*) :: namestr,fmt
        character(len=32) :: tempstr
        call writenum(ivalue,tempstr,fmt)
        call trimzero(tempstr)
        write(unit,*) trim(namestr)//' = '//trim(tempstr)

    end subroutine writeq_si

    !**********************************************************************

    function is_letter(ch) result(res)

        ! Returns .true. if ch is a letter and .false. otherwise

        character :: ch
        logical :: res

        select case(ch)
            case('A':'Z','a':'z')
                res=.true.
            case default
                res=.false.
        end select
        return

    end function is_letter

    !**********************************************************************

    function is_digit(ch) result(res)

        ! Returns .true. if ch is a digit (0,1,...,9) and .false. otherwise

        character :: ch
        logical :: res

        select case(ch)
            case('0':'9')
                res=.true.
            case default
                res=.false.
        end select
        return

    end function is_digit

    !**********************************************************************

    subroutine split(str,delims,before,sep)

        ! Routine finds the first instance of a character from 'delims' in the
        ! the string 'str'. The characters before the found delimiter are
        ! output in 'before'. The characters after the found delimiter are
        ! output in 'str'. The optional output character 'sep' contains the
        ! found delimiter. A delimiter in 'str' is treated like an ordinary
        ! character if it is preceded by a backslash (\). If the backslash
        ! character is desired in 'str', then precede it with another backslash.

        character(len=*) :: str,delims,before
        character,optional :: sep
        logical :: pres
        character :: ch,cha

        pres=present(sep)
        str=adjustl(str)
        call compact(str)
        lenstr=len_trim(str)
        if(lenstr == 0) return        ! string str is empty
        k=0
        ibsl=0                        ! backslash initially inactive
        before=' '
        do i=1,lenstr
            ch=str(i:i)
            if(ibsl == 1) then          ! backslash active
                k=k+1
                before(k:k)=ch
                ibsl=0
                cycle
            end if
            if(ch == '\') then          ! backslash with backslash inactive
                k=k+1
                before(k:k)=ch
                ibsl=1
                cycle
            end if
            ipos=index(delims,ch)
            if(ipos == 0) then          ! character is not a delimiter
                k=k+1
                before(k:k)=ch
                cycle
            end if
            if(ch /= ' ') then          ! character is a delimiter that is not a space
                str=str(i+1:)
                if(pres) sep=ch
                exit
            end if
            cha=str(i+1:i+1)            ! character is a space delimiter
            iposa=index(delims,cha)
            if(iposa > 0) then          ! next character is a delimiter
                str=str(i+2:)
                if(pres) sep=cha
                exit
            else
                str=str(i+1:)
                if(pres) sep=ch
                exit
            end if
        end do
        if(i >= lenstr) str=''
        str=adjustl(str)              ! remove initial spaces
        return

    end subroutine split

    !**********************************************************************

    subroutine removebksl(str)

        ! Removes backslash (\) characters. Double backslashes (\\) are replaced
        ! by a single backslash.

        character(len=*):: str
        character(len=1):: ch
        character(len=len_trim(str))::outstr

        str=adjustl(str)
        lenstr=len_trim(str)
        outstr=' '
        k=0
        ibsl=0                        ! backslash initially inactive

        do i=1,lenstr
            ch=str(i:i)
            if(ibsl == 1) then          ! backslash active
                k=k+1
                outstr(k:k)=ch
                ibsl=0
                cycle
            end if
            if(ch == '\') then          ! backslash with backslash inactive
                ibsl=1
                cycle
            end if
            k=k+1
            outstr(k:k)=ch              ! non-backslash with backslash inactive
        end do

        str=adjustl(outstr)

    end subroutine removebksl

!**********************************************************************

end module xlogp_Strings_m


module xlogp_Vector_m ! filename : xd_Vector.f90
    ! public, everyting by default, but can specify any
    implicit none

    type Vector
        private
        integer :: size ! vector length
        real, pointer, dimension(:) :: data ! component data
    end type Vector

    ! the use of overloaded operators makes sense.
    ! We overload the addition, subtraction, multiplication, assignment, and
    ! logical equal to operators by defining the correct class members to be used
    ! for different argument types.

    !!----------------------------------------------------------------------------
    ! Overload common operators
    interface operator (+)      ! add others later
        module procedure add_Vector, add_Real_to_Vector
    end interface

    interface operator (-)      ! add unary versions later
        module procedure subtract_Vector, subtract_Real
    end interface

    interface operator (*)  ! overload multiplication (*)
        module procedure dot_Vector, real_mult_Vector, Vector_mult_real
    end interface

    interface assignment (=)    ! overload =
        module procedure equal_Real
    end interface

    interface operator (==)     ! overload ==
        module procedure is_equal_to
    end interface

contains ! functions & operators

    function assign(values) result(name) ! array to vector constructor
        real, intent(in) :: values(:)  ! given rank 1 array
        integer :: length       ! array size
        type(Vector) :: name ! Vector to create

        length = size(values)
        allocate(name%data(length))
        name%size = length
        name%data = values
    end function assign

    function make_Vector(len, values) result(v) ! Optional constructor
        integer, optional, intent(in) :: len ! number of values
        real, optional, intent(in) :: values(:) ! given values
        type(Vector) :: v

        if(present(len)) then   ! create vector data
            v%size = len
            allocate(v%data(len))

            if(present(values)) then; v%data = values  ! vector
            else                    ; v%data = 0.d0     ! null vector
            end if ! values present

        else    ! if no len, then make only one element as default
            v%size = 1
            allocate(v%data(1)) ! default

            if(present(values)) then; v%data(1) = values(1) ! scalar
            else                    ; v%data(1) = 0.0d0     ! null
            end if      ! value present

        end if ! len present
    end function make_Vector

    !!----------------------------------------------------------------------------
    !> This will add r to all the array elements of vector v
    !!----------------------------------------------------------------------------
    function add_Real_to_Vector (v, r) result(new)  ! overload +
        type(Vector), intent(in) :: v
        real, intent(in) :: r
        type(Vector) :: new  ! new = v + r

        if (v%size < 1) then
            print*, "No sizes in add_Real_to_Vector"
            !            stop "No sizes in add_Real_to_Vector"
            ! The STOP statement terminates program execution
            ! before the end of the program unit.
            return
        end if

        allocate(new%data(v%size))
        new%size = v%size
        ! new%data = v%data + r ! as array operation
        new%data(1:v%size) = v%data(1:v%size) + r

    end function add_Real_to_Vector

    function add_Vector (a, b) result(new)      ! vector + vector
        type(Vector), intent(in) :: a, b
        type(Vector) :: new       ! new = a + b
        if (a%size /= b%size) then
            print*, "Sizes differ in add_vector"
        endif

        allocate( new%data(a%size) )
        new%size = a%size
        new%data = a%data + b%data
    end function add_Vector

    function copy_Vector(name) result(new)
        type(Vector), intent(in) :: name
        type(Vector) :: new
        allocate (new%data(name%size) )
        new%size = name%size
        new%data = name%data
    end function copy_Vector

    !!----------------------------------------------------------------------------
    !> delete_Vector is the destructor for this class.
    !! In some sense it is incomplete because it does not delete the size attribute.
    !!----------------------------------------------------------------------------
    subroutine delete_Vector(name)
        ! deallocate allocated items
        type(Vector), intent(inout) :: name
        integer :: ok !! Check deallocate status
        deallocate(name%data, stat=ok)
        if (ok/=0) then
            print*, "Vector not allocated in delete_Vector"
        endif
        name%size = 0
    end subroutine delete_Vector

    function dot_Vector(a, b) result(c) ! overload *
        type(Vector), intent(in) :: a, b
        real :: c
        if (a%size /= b%size) then
            print*, "Sizes differ in dot_Vector"
        endif

        c = dot_product(a%data, b%data)
    end function dot_Vector

    subroutine equal_Real(new, R) ! overload =, real to vector
        type(Vector), intent(inout) :: new
        real, intent(in) :: R
        if ( associated (new%data) ) deallocate (new%data)
        allocate (new%data(1))
        new%size = 1
        new%data = R
    end subroutine equal_Real

    logical function is_equal_to(a, b) result(t_f) ! overload ==
        type(Vector), intent(in) :: a, b  ! left & right of ==
        t_f = .false.  ! initialize
        if (a%size /= b%size) return ! same size ?
        t_f = all (a%data == b%data) ! and all values match
    end function is_equal_to

    function length(name) result(n) !! accessor member
        type(Vector), intent(in) :: name
        integer :: n
        n = name%size
    end function length

    subroutine list(name) !! accessor member, for prettier printing
        type(Vector), intent(in) :: name
        print*, "[", name%data(1:name%size), "]"
    end subroutine list

    function normalize_Vector(name) result(new)
        type(Vector), intent(in) :: name
        type(Vector)            :: new
        real :: total, nil = epsilon(1.0) ! tolerance
        allocate( new%data(name%size))
        new%size = name%size
        total = sqrt( sum(name%data**2)) ! intrinsic functions
        if (total < nil) then
            new%data = 0.d0 ! avoid division by 0
        else
            new%data = new%data/total
        endif
    end function normalize_Vector

    subroutine read_Vector(name) ! read array, assign
        type(Vector), intent(inout) :: name
        integer, parameter :: max = 999
        integer :: length
        read (*,'(i1)', advance='no') length
        if (length <= 0) print*, 'Invalid length in read_Vector'
        if (length >= max) print*, 'Maximum length in read_Vector'
        allocate(name%data(length))
        name%size = length
        read*, name%data(1:length)
    end subroutine read_Vector

    function real_mult_Vector(r, v) result(new) ! overload *
        real, intent(in) :: r
        type(Vector), intent(in) :: v
        type(Vector) :: new  ! new = r * v
        if (v%size < 1) print*, "Zero size in real_mult_vector"
        allocate(new%data(v%size))
        new%size = v%size
        new%data = r * v%data
    end function real_mult_Vector

    function size_Vector(name) result(n) ! accessor member
        type(Vector), intent(in) :: name
        integer :: n
        n = name%size
    end function size_Vector

    function subtract_Real(v, r) result(new) ! vector-real, overload -
        type(Vector), intent(in) :: v
        real, intent(in) :: r
        type(Vector) :: new ! new = v + r

        if (v%size < 1 ) print*, "Zero length in subtract_Real"
        allocate(new%data(v%size))
        new%size = v%size
        new%data = v%data - r
    end function subtract_Real

    function subtract_Vector(a, b) result(new) ! overload -
        type(Vector), intent(in) :: a, b
        type(Vector) :: new
        if (a%size /= b%size) print*, "Sizes differ in subtract_Vector"
        allocate( new%data(a%size))
        new%size = a%size
        new%data = a%data - b%data
    end function subtract_Vector

    function values(name) result(array) ! accessor member
        type(Vector), intent(in) :: name
        real :: array(name%size)
        array = name%data
    end function values

    function Vector_(length, values) result(name) ! constructor
        integer, intent(in) :: length ! array size
        real, target, intent(in) :: values(length) ! given array
        real, pointer           :: pt_to_val(:)     ! pointer to array
        type(Vector) :: name        ! Vector to create
        integer :: get_m            ! allocate flag
        allocate( pt_to_val (length), stat = get_m) ! allocate
        if (get_m /= 0) print*, 'allocate error' ! Check
        pt_to_val = values  ! dereference values
        name = Vector(length, pt_to_val) ! intrinsic constructor
    end function Vector_

    function Vector_max_value(a) result(v) ! accessor member
        type(Vector), intent(in) :: a
        real :: v
        v = maxval( a%data(1:a%size) )
    end function Vector_max_value

    function Vector_min_value(a) result(v) ! accessor member
        type(Vector), intent(in) :: a
        real :: v
        v = minval(a%data(1:a%size))
    end function Vector_min_value

    function Vector_mult_real(v, r) result(new) ! vec*real, overload *
        type(Vector), intent(in) :: v
        real, intent(in) :: r
        type(Vector) :: new ! new = v * r
        if (v%size < 1) print*, "Zero size in Vector_mult_real"
        new = Real_mult_Vector(r, v)
    end function Vector_mult_real

end module xlogp_Vector_m
