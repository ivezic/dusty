! ***********************************************************************
! This function is taken from Moshe Elitzur.           [Z.I., Nov. 1995]
! =======================================================================
double precision function RDINP(Equal,iUnit)
  ! =====================================================================
  ! Read lines, up to 232 long, from pre-opened unit iUnit and extract
  ! all input numbers from them. When EQUAL is set, numeric input data
  ! must be preceded by an equal sign.All non-numeric data and numbers
  ! not preceded by = when EQUAL is on are ignored.RDINP = next number
  ! encountered (after equal sign) and terminated by a nun-numeric
  ! symbol (termination with blank is best). Commas and exponential
  ! notation are allowed.  All text after % is ignored, as in TeX.
  ! Lines with * in the first column are echoed to the output device.
  ! The number is comprised of an actual VALUE, decimal part FRAC and
  ! (integer) exponent PWR.  It has a SIGN, and the exponential power
  ! has sign SIGNEX. Logical flag to the decimal part is set in
  ! DECIMAL. The search is conducted between FIRST, which is
  ! continuously increased, and LAST.  A new line is read when FIRST
  ! exceeds LAST, and can be forced by calling with -iUnit.  Actual
  ! extraction of numerical value is done in separate FUNCTION VAL.
  ! ====================================================================
  IMPLICIT None
  Integer iUnit,ind,First, Last
  DOUBLE PRECISION Value,Val,Frac,Pwr,Sign,Signex
  CHARACTER Card*(232),CR,prev,Term,Next
  Logical Equal,digit,minus,sgn,dot,E,decimal
  Save Card,First,Last
  DATA First/1/, Last/0/
  ! --------------------------------------------------------------------
  ! FUNCTION STATEMENTS
  digit(CR) = CR.GE.'0' .AND. CR.LE.'9'
  minus(CR) = CR.EQ.'-'
  sgn(CR)   = CR.EQ.'+' .OR. CR.EQ.'-'
  dot(CR)   = CR.EQ.'.'
  E(CR)     = CR.EQ.'E' .OR. CR.EQ.'e'
  !
  IF (iUnit.lt.0) Then
     First = Last + 1
     iUnit = -iUnit
  END IF
  ! Start the search for the next number:
1 RDINP  = 0.
  VALUE  = 0.
  FRAC   = 0.
  PWR    = 0.
  SIGN   = 1.
  SIGNEX = 1.
  Decimal = .False.
  If (first.gt.last) then
     ! Time to get a new line
     READ (iUnit, '(A)' , END = 99) Card
     first = 1
     last = len(Card)
     ! Find start of trailing junk:
     DO WHILE (Card(last:last).LE.' ')
        last = last - 1
        if (last.lt.first) goto 1
     END DO
     IF (Card(first:first).EQ.'*') WRITE (12,'(A)') Card(1:last)
     ind = Index(Card,'%')
     if (ind.gt.0) last = ind - 1
  End If
  ! Get past the next '=' when the EQUAL flag is set
  If (Equal) then
     DO WHILE (Card(first:first).NE.'=')
        first = first + 1
        IF (first.GT.last) goto 1
     END DO
  End If
  ! OK, start searching for the next digit:
  ! Do While (.not.digit(Card(first:first))) - corrected by ME
  Do While (.not.(Card(first:first).ge.'0'.AND.Card(first:first).le.'9'))
     first = first + 1
     if (first.gt.last) goto 1
  End Do
  ! Check if it is a negative or decimal number
  If (first.gt.1) then
     prev = Card(first-1:first-1)
     ! necessary syntax changes, otherwise the code crushes on ncx[MN]
     ! if (minus(prev)) sign = -1.
     if (prev.eq.'-') sign = -1.
     ! if (dot(prev)) then
     if (prev.eq.'.') then
        decimal = .True.
        if (first.gt.2 .and.(Card(first-2:first-2).eq.'-')) sign = -1.
     end if
  End If
  ! Extract the numerical value
  IF (.not.Decimal) Then
     Value = VAL(Card,first,last,decimal,Term)
     ! Check for a possible decimal part.  Termination with '.E' or
!        '.e' is equivalent to the same ending without '.'
     ! if (first.lt.last.and.dot(Term)) then
     if (first.lt.last.and.Term.eq.'.') then
        first = first + 1
        next = Card(first:first)
        ! if (digit(next)) decimal = .true.
        if (next.GE.'0' .AND. next.LE.'9') decimal = .true.
        ! if (E(next)) Term = 'E'
        if (next.eq.'E' .OR. next.EQ.'e') Term = 'E'
     end if
  END IF
  ! Extract the decimal fraction, when it exists
  IF (Decimal) Frac = VAL(Card,first,last,decimal,Term)
  ! An exponential may exist if any part terminated with 'E' or 'e'
  ! IF (first.lt.last.and.E(term)) then
  IF (first.lt.last.AND.(term.eq.'E'.or.term.eq.'e')) then
     first = first + 1
     next = Card(first:first)
     ! if (first.lt.last.and.sgn(next))then
     if (first.lt.last.AND.(next.EQ.'+'.or.next.EQ.'-')) then
        first = first + 1
        ! if (minus(next)) Signex = -1.
        if (next.eq.'-') Signex = -1.
     end if
     decimal = .False.
     PWR = VAL(Card,first,last,decimal,Term)
  END IF
  ! Finally, put the number together
  RDINP = Sign*(Value + Frac)*10**(Signex*PWR)
  Return
99 WRITE (12,'(3(1x,a,/))')                                           &
        '*** TERMINATED. EOF reached by RDINP while looking for input.', &
        '*** Last line read:',Card
  ! --------------------------------------------------------------------
  RETURN
end function RDINP
! ***********************************************************************

!***********************************************************************
subroutine rdinps2(equal,iUnit,str,LENGTH,UCASE)
!=======================================================================
!     2nd version. Returns also length of meaningful string; if UCASE flag
!     is set, the returned string is in upper case to avoid keyboard entry problems
!=======================================================================
      integer i,iUnit,ind,first,last,next,length
      character card*(232),chr,cr
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
        if(chr(first) .eq. '*') write(12,'(a)') card
        last = len(card)
        ind = index(card,'%')
        if(ind .gt. 0) last = ind-1
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

99    write(16,'(3(1x,a,/))')  &
     ' Terminated. EOF reached by rdinp while looking for input. ', &
     ' last line read:',card
!-----------------------------------------------------------------------
      return
end subroutine rdinps2
!***********************************************************************

!************************************************************************
! This function is taken from Moshe Elitzur            [Z.I., Nov. 1995]
!========================================================================
double precision function Val(Card,first,last,decimal,term)
  !======================================================================
  !  Extract numerical value from CARD, begining at position FIRST up
  !  to the first non-digit encountered.  The terminating character is
  !  returned in TERM, its position in FIRST. Commas are treated as
  !  simple separators.
  !======================================================================
  implicit none
  character card*(*), term, ch
  logical decimal, digit
  integer ival, first, last, first0
  double precision pwr
!-----------------------------------------------------------------------
! function statements
  ival(ch) = ichar(ch) - ichar('0')
  digit(ch) = ch.ge.'0' .and. ch.le.'9'

  val = 0.0d0
  pwr = 1.0d0
  first0 = first
  do 10 first = first0, last
   term = card(first:first)
   if (term.eq.',') goto 10
! if (.not.digit(term)) Return - corrected by ME
   if (.not.(term.ge.'0' .and. term.le.'9')) return
   if (decimal) then
    pwr = pwr*0.1d0
    val = val + pwr*ival(term)
   else
    val = 10.0d0*val + ival(term)
   end if
10 continue
   term = ' '
!-----------------------------------------------------------------------
   return
end function val
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

