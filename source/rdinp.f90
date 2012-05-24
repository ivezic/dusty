

!****************************************************************************!
!                                                                           !
! Replacement of the old RDINP. Simpler operation -- the number is          !
! obtained straight by reading from internal unit instead of the old        !
! method which required parsing and evaluating. Also, introduces            !
! the input variable outUnit, the unit number to which this                 !
! function will write output. In the old version this was hard coded        !
! as unit 12, which was really a bug.                                       !
!                                                    [M.E., March 2012]     !
!                                                                           !
!===========================================================================!
double precision function RDINP(Equal, inUnit, outUnit)
! ==========================================================================!
!  Read lines, up to 255 long, from pre-opened unit inUnit and extract      !
!  all input numbers from them. When EQUAL is set, numeric input data       !
!  must be preceded by an equal sign. All non-numeric data and numbers      !
!  not preceded by = are ignored when EQUAL is on.                          !
!                                                                           !
!  RDINP = next number encountered (after equal sign) and terminated blank. !
!  The blank can be optionally preceded by comma.                           !
!  Input numbers can be written in any FORTRAN allowable format.            !
!  In addition, comma separation is allowed in the input numbers.           !
!                                                                           !
!  All text after % (or !) is ignored, as in TeX (or F90).                  !
!  Lines with * in the first column are echoed to output on pre-opened      !
!  unit outUnit.                                                            !
!                                                                           !
!  The search is conducted between FIRST, which is                          !
!  continuously increased, and LAST.  A new line is read when FIRST         !
!  exceeds LAST, and can be forced by calling with -iUnit.                  !
!                                                                           !
!===========================================================================!
   IMPLICIT None
   Integer :: inUnit, outUnit, ind, First = 1, Last = 0
   Logical, intent(in) :: Equal
   CHARACTER(255) Card, no_comma
   Save Card, First, Last
! -----------------------------------------------------------------------
!
   IF (inUnit.lt.0) Then                         ! force a new line
      First = Last + 1
      inUnit = -inUnit
   END IF

   DO
      If (first > last) then                     ! Time to get a new line
         READ (inUnit, '(A)' , END = 99) Card
         if (len_trim(Card) == 0) cycle          ! ignore empty lines
                                                 ! Echo to output lines that start with * 
         IF (Card(1:1) == '*') &
             WRITE (outUnit,'(A)') TRIM(Card(1:len_trim(card)))
         Card = trim(no_comma(Card))             ! remove commas to allow comma-separated numbers
         first = 1
         last = len_Trim(Card)
         ind = MAX(Index(Card,'%'), Index(Card,'!'))
         if (ind.gt.0) last = ind - 1            ! Everything after % and ! is ignored 
      End If

     !Get past the next '=' when the EQUAL flag is set
      If (Equal) then
        DO WHILE ((Card(first:first) /= '=') .and. (first <= last) )
          first = first + 1
        END DO
      End If
      first = first + 1
      IF (first > last) cycle

     !Find start of next number; necessary when the EQUAL flag is off
      Do While (.not. (Card(first:first).ge.'0'.AND.Card(first:first).le.'9') &
                .and. (first <= last))
          first = first + 1
      End Do
      if (first > last) cycle

     !OK; time to get the number
      READ(card(first - 1:last), *, ERR = 98) RDINP 
      
     !and move past its end
      Do While (Card(first:first) /= ' ')
          first = first + 1
      End Do
      
      return
   End DO

98 WRITE (outUnit,'(3(1x,a,/))')                               &
   ' ****ERROR. RDINP could not read a proper number from', &
   Card(first - 1:last),																			 &
   ' ****Number should be preceded and terminated by spaces'
   RETURN
99 WRITE (outUnit,'(3(1x,a,/))')                                      &
   ' ****TERMINATED. EOF reached by RDINP while looking for input. ', &
   ' *** Last line read:', Card
   RETURN
end function RDINP
! ***********************************************************************


!***********************************************************************
!=======================================================================
function no_comma(str)
!=======================================================================
! Remove all commas from string str. This enables input to RDINP of 
! comma separated numbers
!=======================================================================
   character(len = *), intent(in) :: str
   character(255) temp, no_comma
   integer l, k
   
   ! First create a blank string:      
   do l = 1, len(temp)
      temp(l:l) = ' '  
   end do 

   ! Now fill the blank string with the characters
   ! of str, skipping the commas:
   k = 1
   do l = 1, len_trim(str)
      if (str(l:l) /= ',') then
         temp(k:k) = str(l:l)
         k = k + 1
      end if
   end do

   no_comma = trim(temp)
   return
   
end function no_comma
!***********************************************************************



!***********************************************************************
subroutine rdinps2(equal,iUnit,outUnit,str,LENGTH,UCASE)
!=======================================================================
!     2nd version. Returns also length of meaningful string; if UCASE flag
!     is set, the returned string is in upper case to avoid keyboard entry 
!     problems.
!     Added outUnit for output; original version was hard coded to write
!     to unit 6---a bug.
!=======================================================================
      integer i,iUnit,outUnit,ind,first,last,next,length
      character card*(255),chr,cr
      character*(*) str, ch*1
      logical Equal,bar,UCASE
      save card, first, last
      data first/1/, last/0/
!!-----------------------------------------------------------------------
!     Function statements
      chr(i)    = card(i:i)
      bar(cr)   = cr .eq. ' '
!
      if(iUnit .lt. 0) then
        first = last + 1
        iUnit = -iUnit
      end if
!        start the search for the next string:
  1   continue
      if(first .gt. last) then
        read(iUnit, '(a)', end = 99) card
        first = 1
        if(chr(first) .eq. '*') write(outUnit,'(a)') card
        last = len(card)
        ind = MAX(index(card,'%'), index(card,'!'))
        if(ind .gt. 0) last = ind - 1
      end if
      if(equal) then
        do while (chr(first) .ne. '=')
          first = first + 1
          if(first .gt. last) goto 1
        end do
      end if
      first = first + 1

      do while (bar(chr(first)) .and. first .le. last)
        first = first + 1
      end do
      if(first .gt. last) goto 1

      next=first+1
      do while( .not. bar(chr(next)) )
        next = next + 1
      end do
      str=card(first:next-1)
      length = next - first

      if (UCASE) then
!        convert string to UPPER CASE
         DO i = 1, length
            ch = str(i:i)
            IF ((ch .GE. 'a') .AND. (ch .LE. 'z'))  &
              str(i:i) = Char(IChar(ch)+IChar('A')-IChar('a'))
         END DO
      end if
      return

99    write(outUnit,'(3(1x,a,/))')  &
     ' Terminated. EOF reached by rdinp while looking for input. ', &
     ' last line read:',card
!-----------------------------------------------------------------------
      return
end subroutine rdinps2
!***********************************************************************

!***********************************************************************
INTEGER FUNCTION EMPTY(line)
!=======================================================================
! This function is 1 if string 'line' is empty, or if it contains only
! '%', and 0 otherwise.
!                                                      [Z.I., Nov. 1996]
! =======================================================================
  IMPLICIT none
  INTEGER i, iTeX, l
  CHARACTER ch
  CHARACTER*(*) line
  ! ---------------------------------------------------------------------
  l = LEN(line)
  EMPTY = 1
  iTeX = 0
  DO i = 1, l
     ch = line(i:i)
     IF(EMPTY.EQ.1.AND.ch.EQ.'%') iTeX = 1
     IF (ch.NE.' ') EMPTY = 0
  END DO
  IF (iTeX.EQ.1) EMPTY = 1
  ! ---------------------------------------------------------------------
  RETURN
END FUNCTION EMPTY
!***********************************************************************

!************************************************************************
subroutine clean(strin, strout, length)
  !======================================================================
  ! Find meaningful part of strin without leading and trailing junk
  ! It is returned left-justified in StrOut, right-padded with blanks
  ! The number of meaningful characters is returned in Length.
  ! In case of any problems, StrOut is empty. This sub should be used to
  ! clean every input filename immediately after DUSTY reads it. [ME,'99]
  !======================================================================
  implicit none
  character*(*) strin, strout
  integer i, first, last, length
  !-----------------------------------------------------------------------
  do i = 1, len(strout)
     strout(i:i) = ' '
  end do
  first = 1
  last = len(strin)
  if (first.gt.last) return
  ! Find end of leading junk:
  do while (strin(first:first).le.' ')
     first = first + 1
     if (first.gt.last) return
  end do
  ! Find start of trailing junk:
  do while (strin(last:last).le.' ')
     last = last - 1
     if (last.lt.first) return
  end do
  ! Now trim all junk:
  strout = strin(first:last)
  length = last - first + 1
  !---------------------------------------------------------------------
  return
end subroutine clean
!***********************************************************************

