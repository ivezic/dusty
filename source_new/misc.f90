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

!***********************************************************************
subroutine attach(root,length,ext,fname)
  !=====================================================================
  ! Attaches extensions to the root cleaned by Clean
  !=====================================================================
  implicit none
  integer i, length
  character*(*) root, ext, fname
  ! ---------------------------------------------------------------------
  do i = 1, len(fname)
     fname(i:i) = ' '
  end do
  fname(:length) = root(:length)
  fname(length + 1:) = ext
  !----------------------------------------------------------------------
  return
end subroutine attach
!************************************************************************

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

!***********************************************************************
subroutine fileMSG(fname,strg)
!=======================================================================
! Prints a message in *.out file in case of error opening the user
! supplied files.
!=======================================================================
  implicit none
  character aux*235, strg*(*), fname*(*)
  integer length, empty
!-----------------------------------------------------------------------
1 read(1,'(a)') aux
  if (empty(aux).eq.1) goto 1
  call clean(aux,fname,length)

  open(10, err=100, file=fname, status='old')
  close(10)
  return

100 write(12,*)' *** Fatal error in dusty! **************************'
  write(12,*)' File with the ',strg
  write(12,'(a2,a)')'  ',fname
  write(12,*)' is missing ?!'
  write(12,*)' ****************************************************'
!  close(12)
!-----------------------------------------------------------------------
  stop
end subroutine fileMSG
!***********************************************************************

!***********************************************************************
subroutine ChkAngle(angle)
!=======================================================================
! Checks if the input angles are in [0,85] degrees interval.  [MN,2005]
!=======================================================================
  implicit none
  integer error
  double precision angle
!-----------------------------------------------------------------------
  if (angle.gt.88.0d0) then
     write(12,*)' ********** Message from input: ***************'
     write(12,*)' * slab illumination angles have to be in the *'
     write(12,*)' * range [0,88] degrees. setting the maximum  *'
     write(12,*)' * illumination angle to 88 degrees.          *'
     write(12,*)' **********************************************'
     angle = 88.0d0
  else
     if (angle.lt.0.0d0) then
        write(12,*)' ********** Message from input: ***************'
        write(12,*)' * slab illumination angles have to be in the *'
        write(12,*)' * range [0,85] degrees. setting the miminum  *'
        write(12,*)' * illumination angle to 0.0 degrees.         *'
        write(12,*)' **********************************************'
     end if
  end if
!-----------------------------------------------------------------------
  return
end subroutine ChkAngle
!***********************************************************************

!***********************************************************************
subroutine MSG(msgno)
!=======================================================================
! This subroutine writes runtime messages to auxiliary file fname.m##
! or to the output file fname.out.             [ZI,Feb'96; MN,Jul'99]
!=======================================================================
  use common
  implicit none

  integer  msgno
!-----------------------------------------------------------------------

  if (msgno.eq.1.and.iX.gt.0) then
   write(18,*)' ************  WARNING  *************'
   write(18,*)' Temperature calculation in FindTemp'
   write(18,*)' achieved the limit of 500 iterations'
  end if
  if (msgno.eq.2.and.iX.gt.0) then
   write(18,*)' ************  WARNING  **************'
   write(18,*)' Energy density iterations in radtransf'
   write(18,*)' achieved the limit of 50000 iterations'
  end if
  if (msgno.eq.3) then
   write(12,*)' **********  INPUT ERROR ***********'
   write(12,*)' * Denstyp is not between 0 and 5! *'
   write(12,*)' * Check input file and try again.  *'
   write(12,*)' ***********************************'
  end if
  if (msgno.eq.4.and.iX.gt.0) then
   write(18,*)' ************  WARNING  *****************'
   write(18,*)' Could not bracket in zbrac (in sub FindTemp).'
   write(18,*)' Something might be wrong in your input.     '
  end if
  if (msgno.eq.5) then
   write(12,*)' ********** INPUT ERROR ***********'
   write(12,*)' * All abundances must be >= 0!  *'
   write(12,*)' * Check input file and try again *'
   write(12,*)' ***********************************'
  end if
  if (msgno.eq.6) then
   write(12,*)' ********** INPUT ERROR ***********'
   write(12,*)' * Wavelengths for the power-law  *'
   write(12,*)' * spectrum must be ascending!    *'
   write(12,*)' * Check input file and try again *'
   write(12,*)' ***********************************'
  end if
  if (msgno.eq.7) then
   write(12,*)' ********** INPUT ERROR ***********'
   write(12,*)' * Relative luminosities must add *'
   write(12,*)' * up to a number >0!!           *'
   write(12,*)' * Check input file and try again. *'
   write(12,*)' ***********************************'
  end if
  if (msgno.eq.8) then
   write(12,*)' ********** A BIG ERROR ***********'
   write(12,*)' *    A black body temperature    *'
   write(12,*)' *        should be > 0 !!        *'
   write(12,*)' * Check input file and try again. *'
   write(12,*)' ***********************************'
  end if
  if (msgno.eq.9) then
   write(12,*)' ********** INPUT ERROR ***********'
   write(12,*)' * The flag for optical properties *'
   write(12,*)' * should be between 1 and 3!!   *'
   write(12,*)' * Check input file and try again. *'
   write(12,*)' ***********************************'
  end if
  if (msgno.eq.10) then
   write(12,*)' ********** INPUT ERROR ***********'
   write(12,*)' * The flag for size distribution  *'
   write(12,*)' * should be between 1 and 3!!    *'
   write(12,*)' * Check input file and try again. *'
   write(12,*)' ***********************************'
  end if
  if (msgno.eq.11) then
   write(12,*)' ********** INPUT ERROR ***********'
   write(12,*)' * The flag for external spectrum *'
   write(12,*)' * should be between 1 and 6 !!  *'
   write(12,*)' * Check input file and try again. *'
   write(12,*)' ***********************************'
  end if
  if (msgno.eq.12) then
   write(12,*)' ********** INPUT ERROR  **********'
   write(12,*)' Only three types of the point spread '
   write(12,*)' function are allowed: 1, 2 or 3 !!  '
   write(12,*)' Check input file and try again.       '
   write(12,*)' ***********************************'
  end if
! msg 14 is not called in this version.
  if (msgno.eq.14.and.iX.gt.0) then
   write(18,*)' ******** message from solve *********'
   write(18,*)' Convergence on en.density is too slow.  '
   write(18,*)' If the accuracy is not reached yet the  '
   write(18,*)' code will increase grid size and try again.'
   write(18,*)' ****************************************'
  end if
  if (msgno.eq.15) then
   write(12,*) ' ***************** WARNING *******************'
   write(12,*) '  No calculation for next model. Parameter npY'
   write(12,*) '  needs to be at least 50. Use of the slab    '
   write(12,*) '  parameters is suggested (see userpar.inc)   '
   write(12,*) ' *********************************************'
  end if
  if (msgno.eq.16) then
   write(12,*)' ****************  WARNING  ******************'
   write(12,*)'  The density profile eta is too steep and the'
   write(12,*)'  code can not handle this. Try decreasing the'
   write(12,*)'  outer radius Yout (see manual, 3.3.3).      '
   write(12,*)' *********************************************'
   if(iX.gt.0) then
    write(18,*)' ****************  WARNING  ******************'
    write(18,*)'  The density profile eta is too steep and the'
    write(18,*)'  code can not handle this. Try decreasing the'
    write(18,*)'  outer radius Yout. (see manual, 3.3.3).     '
    write(18,*)' *********************************************'
   end if
  end if
  if (msgno.eq.17) then
   write(12,*)' *****************  WARNING  ********************'
   write(12,*)'  Eta is too steep and reaches values less than  '
   write(12,*)'  1e-12. Try decreasing the outer radius Yout.   '
   write(12,*)'  (see manual,3.3.3)                             '
   write(12,*)' ************************************************'
   if(iX.gt.0) then
    write(18,*)' *****************  WARNING  ********************'
    write(18,*)'  Eta is too steep and reaches values less than  '
    write(18,*)'  1e-12. Try decreasing the outer radius Yout.   '
    write(18,*)'  (see manual,3.3.3)                             '
    write(18,*)' ************************************************'
   end if
  end if
  if (msgno.eq.18) then
   write(12,*)' ************  WARNING  ************************ '
   write(12,*)'  The dynamical range of eta is more than 1e-12. '
   write(12,*)'  The outer radius Yout must be decreased so that'
   write(12,*)'  eta does not go below 1e-12 (see manual,3.3.3) '
   write(12,*)' *********************************************** '
   if(iX.gt.0) then
    write(18,*)' ************  WARNING  ************************ '
    write(18,*)'  The dynamical range of eta is more than 1e-12. '
    write(18,*)'  The outer radius Yout must be decreased so that'
    write(18,*)'  eta does not go below 1e-12 (see manual,3.3.3)'
    write(18,*)' *********************************************** '
   end if
  end if
  if (msgno.eq.19) then
   write(12,*)' ************ a big error!!************* '
   write(12,*)'  singular matrix in ludcmp when called   '
   write(12,*)'  from analint. stopping the calculation. '
   write(12,*)' **************************************** '
  end if
  if (msgno.eq.20) then
   write(12,*)' ************ a big error!! ************ '
   write(12,*)'  singular matrix in ludcmp when called   '
   write(12,*)'  from invert. stopping the calculation.  '
   write(12,*)' **************************************** '
  end if
  if (msgno.eq.21) then
   write(12,*)' ************ INPUT ERROR ******************'
   write(12,*)' * The type-of-entry value for the central *'
   write(12,*)' * source has to be between 1 and 4 !!     *'
   write(12,*)' * Correct your input and try again.       *'
   write(12,*)' *******************************************'
  end if
  if (msgno.eq.22) then
   write(12,*)' ************ INPUT ERROR ***************'
   write(12,*)' * The typd-of-entry value for external *'
   write(12,*)' * radiation has to be 1 or 2!          *'
   write(12,*)' * Correct your input and try again!    *'
   write(12,*)' ****************************************'
  end if
  if (msgno.eq.23) then
   write(12,*)' ******** message from input: **************'
   write(12,*)' * Left-hand side source is always present *'
   write(12,*)' * in slab geometry. Setting left=1.       *'
   write(12,*)' *******************************************'
  end if
  if (msgno.eq.24) then
   write(12,*)' ******** message from input: ************'
   write(12,*)' * Please choose the version of velocity *'
   write(12,*)' * equation: 1 -linear, 2 -quadratic!    *'
   write(12,*)' * Vorrect your input and try again!     *'
   write(12,*)' *****************************************'
  end if
  if (msgno.eq.25) then
   write(12,*)' ************ INPUT ERROR ***************'
   write(12,*)' * The typd-of-entry value for spectrum *'
   write(12,*)' * from a file can be only 1 or 2!      *'
   write(12,*)' * Correct your input and try again!    *'
   write(12,*)' ****************************************'
  end if
  if (msgno.eq.26) then
   write(12,*)' ********** message from input: *************'
   write(12,*)' * If external spectrum is from a file      *'
   write(12,*)' * type-of-entry is always 2 (normalization *'
   write(12,*)' * factor). Taking factor = 1               *'
   write(12,*)' ********************************************'
  end if
!-----------------------------------------------------------------------

  return
end subroutine MSG
!***********************************************************************

!***********************************************************************
subroutine ReadSpectar(lambdas,Llamstar,Lstar,nLs,is,error)
!=======================================================================
! Reads the source spectrum from a file. This was part of Subroutine Star,
! separated for clarity. is=1 for central s-ce, is=2 for external [MN'01]
! =======================================================================
  use common
  implicit none

  character*235 line
  integer ios1, iLs, nLs, error, nLambdam, Nis, is
! nLambdam is the max number entries for a user supplied stellar spectrum
  parameter (nLambdam = 10000, Nis = 2)
  double precision lambdas(nLambdam), Llamstar(nLambdam),lLs(nLambdam), &
       Ls(nLambdam), Lstar, a, b
! -----------------------------------------------------------------------

! For  startyp.ge.4 stellar spectrum is read from the file 'namestar'
! which is unit=3 (unit=1 is the input file)
! is=1 for the enclosed source, is=2 for the external shell illumination
! possible problems with the file are checked in sub inp_rad
    open(3,file=namestar(is),status='old')
      rewind(3)
      do iLs = 1, 3
        read(3,'(a235)',err=998) line
      end do
      ios1 = 0
      iLs = 0
      do while (ios1.ge.0)
        read(3,*,end=900,err=998,iostat=ios1) a, b
        if(ios1.ge.0) then
          iLs = iLs + 1
          lambdas(iLs) = a
          if (a.le.0.0) goto 998
!         Llamstar is always f_lambda
!         if startyp.eq.4 then file gives lambda*f_lambda
          if (startyp(is).eq.4) Llamstar(iLs) = b/a
!         if startyp.eq.5 then file gives f_lambda
          if (startyp(is).eq.5) Llamstar(iLs) = b
!         if startyp.eq.6 then file gives lnu=lambda**2*f_lambda
          if (startyp(is).eq.6) Llamstar(iLs) = b/(a**2)
        end if
      end do
900   close(3)
      if (iLs.lt.2) goto 998
      nLs = iLs
!     if input wavelengths in descending order turn them around
      if (lambdas(1).gt.lambdas(2)) then
        do iLs = 1, nLs
          lLs(iLs) = lambdas(iLs)
          Ls(iLs) =  Llamstar(iLs)
        end do
        do iLs = 1, nLs
          lambdas(iLs) = lLs(nLs+1-iLs)
          Llamstar(iLs) = Ls(nLs+1-iLs)
        end do
      end if
!     Get the scale, Lstar, of the stellar spectrum
      call Simpson(nLambdam,1,nLs,lambdas,Llamstar,Lstar)
      error = 0
      goto 999
998   write(12,*)' *********** INPUT ERROR *************************'
      write(12,*)' The file with spectral shape of external radiation:'
      write(12,'(a2,a100)')'  ',namestar(1)
      write(12,*)' is missing or not properly formatted?!'
      write(12,*)' ***************************************************'
      error = 3
! -----------------------------------------------------------------------

999 return
end subroutine ReadSpectar
!***********************************************************************

!***********************************************************************
subroutine WriteOut(var1,var2,var3,nG,nameQ,nameNK)
!=======================================================================
! WriteOut prints in fname.out all input parameters,  read before density distribution
! type.
!=======================================================================

  use common
  implicit none

  integer is, iG, nG, i, length
  double precision var1, var2, var3
  character*72 strpow, aux, src, chaux*3
  character*(*) nameQ(npG), nameNK(10)
  logical first
  first = .true.
!-------------------------------------------------------------------------
  is = 1
15 if (SLB) then
   if (is.eq.1) then
    src = 'Left-side source spectrum described by'
   else
    src = 'Right-side source spectrum described by'
   end if
  else
   if (is.eq.1) then
    src = 'Central source spectrum described by'
   else
    src = 'External source spectrum described by'
   end if
  end if
  call Clean(src, aux, length)

  if(Left.eq.0.and.is.eq.1) then
   write(12,*) ' No central source.'
  else
! #1: black body(ies) for startyp=1
   if (startyp(is).eq.1) then
    if (nBB(is).gt.1) then
     call ATTACH(aux, length, ' ', src)
! multiple black bodies
     write(12,'(a2,a37,i2,a13)')'  ', src, nBB(is),' black bodies'
     write(12,'(a27)')' with temperatures (in K):'
     write(12,'(2x,1p,10e10.3)')(Tbb(is,i),i=1,nBB(is))
     write(12,'(a42)')' and relative luminosities, respectively:'
     write(12,'(1p,10e10.1)')(rellum(is,i),i=1,nBB(is))
    else
! for a single black body:
     call ATTACH(aux,length,' a black body',src)
     write(12,'(a2,a)') '  ',src
     write(12,'(a19,1p,e10.3,a2)')' with temperature:',Tbb(is,1),' K'
    end if
   end if

! #2: Engelke-Marengo function for startyp=2
   if (startyp(is).eq.2) then
    call ATTACH(aux, length,' Engelke-Marengo function', src)
    write(12,'(a2,a)') '  ',src
    write(12,'(a13,1p,e10.3,a16)')' with Teff =',Tbb(is,1), ' K and depth of'
    write(12,'(a30,F6.1,a2)')' the SiO absorption feature =', xSiO,' %'
   end if

! #3: power-law(s) for startyp=3
   if (startyp(is).eq.3) then
    if (Nlamtr(is).gt.0) then
     call ATTACH(aux,length,' power law:',src)
     write(12,'(a2,a)') '  ',src
     write(12,*)'    lambda      k'
     do i = 1, Nlamtr(is)
      write(12,'(1x,1p,e10.3)')lamtr(is,i)
      write(12,'(11x,1p,e10.3)')klam(is,i)
     end do
     write(12,'(1x,1p,e10.3)')lamtr(is,Nlamtr(is)+1)
    else
     write(12,*)' Input data for the source spectrum is not good.'
     write(12,*)' Changed to a 10000 K black body'
    end if
   end if

! spectrum from a file for startyp=4,5,6
   if (startyp(is).ge.4.and.startyp(is).le.6) then
    if (is.eq.1) then
     write(12,*)' Stellar spectrum supplied from file:'
    else
     write(12,*)' External spectrum supplied from file:'
    end if
    write(12,'(a2,a100)') '  ',nameStar(is)
    call PrHeader(3,nameStar(is))
   end if
  end if
  write(12,*)' --------------------------------------------'
  if(first) then
! if there is a second source go back to read its parameters
   if(Right.gt.0) then
! repeat printout of source info for the external radiation
! its index is is=2
    is = 2
    first = .false.
    goto 15
   end if
  end if
  is = 1
! -----------------------------------------------------
! Boundary Condition data: typEntry(is) and value
  if(Left.gt.0) then
   if(typEntry(1).eq.5) then
     if (SLB) then
      write(12,'(a45,1p,e9.2,a,i2)')  &
           ' Dust temperature on the slab left boundary:', Tsub(ifidG),' K - Grain:',ifidG
     else
      write(12,'(a41,1p,e9.2,a,i2)')  &
           ' Dust temperature on the inner boundary:', Tsub(ifidG),' K - Grain:',ifidG
     end if
   else if (typEntry(1).eq.1) then
    if (slb) then
     if (startyp(1).gt.3) then
      write(12,'(a33,1p,e9.2,a6)') &
           ' Flux at the slab left boundary:', dilutn*var1,' W/m^2'
     else
      write(12,'(a33,1p,e9.2,a6)') &
           ' Flux at the slab left boundary:', Ji*4*pi,' W/m^2'
     end if
! if input spectrum is from a file its bol.flux is calculated
! and then normalized with the value of 'dilutn' from the input file.
    elseif(sph) then
     if (startyp(1).gt.3) then
      write(12,'(a29,1p,e9.2,a6)') &
           ' Flux at the inner boundary:', dilutn*var1,' W/m^2'
     else
      write(12,'(a29,1p,e9.2,a6)') &
           ' Flux at the inner boundary:', var1,' W/m^2'
     end if
    end if
   else if (typEntry(1).eq.2) then
    write(12,'(a20,1p,e9.2,a17,1p,e9.3,a3)') &
         ' Source luminosity ', var1,' Lo and distance ', var2, ' cm'

   else if (typEntry(1).eq.3.or.typEntry(1).eq.4) then
    write(12,'(a38,1p,e9.2,a6)') &
         ' Mean intensity at the inner surface:', var3,' W/m^2'
  end if
 end if
 if(Right.gt.0) then
  if (sph) then
   if (typEntry(1).eq.1) then
    if (startyp(2).gt.3) then
     write(12,'(a49,1p,e9.2)')  &
          '  Normalization factor of the external radiation:', dilutn
     write(12,'(a54,1p,e9.2,a6)')  &
          '  Mean intensity from the file after renormalization:',var3,' W/m^2'
    else
     write(12,'(a45,1p,e9.2)') &
          '  Dilution factor of the external radiation:', dilutn

    end if
   else
    write(12,'(a44,1p,e9.2,a6)')  &
         '  Mean intensity of the external radiation:', var3,' W/m^2'
   end if
  end if
 end if
 write(12,*)' --------------------------------------------'
! -----------------------------------------------------
! 2) DUST PROPERTIES
!  2.1 Chemical Composition
  if (top.lt.3) then
   write(12,*)' Abundances for supported grains:'
   write(12,*)' Sil-Ow Sil-Oc Sil-DL grf-DL amC-Hn SiC-Pg'
   write(12,'(6f7.3)')(xC(i),i=1,3),xC(4)+xC(5),(xC(i),i=6,7)
   if (top.eq.2) then
    write(12,*)' Abundances for user supplied grains:'
    write(12,'(i6,9i7)')(i,i=1,Nfiles)
    write(12,'(10f7.3)')(xCuser(i),i=1,Nfiles)
    write(12,*)' User supplied n and k from:'
    do i = 1, Nfiles
     write(12,'(a2,i1,a2,a70)')'  ',i,') ',nameNK(i)
    end do
   end if
! user supplied cross-sections:
  else
   do iG = 1, nG
    write(12,*)' Optical properties from file(s):'
    write(12,'(a2,a70)')'  ',nameQ(iG)
    call PrHeader(3,nameQ(iG))
   end do
  end if
  write(12,*)' Sublimation Temperature(s):'
  do iG=1,nG
     write(12,'(a8,i2,a3,f10.3)') '  Grain(',iG,'): ',Tsub(iG) 
  end do
! 2.2 Grain size distribution
  if (top.ne.3) then
   if (szds.eq.3) then
    chaux = 'KMH'
   else
    chaux = 'MRN'
   end if
   write(12,'(a2,a3,a19)')'  ',chaux,'size distribution:'
   call getfs(qsd,1,0,strpow)
   write(12,'(a15,a5)')'      Power q:',strpow
   write(12,'(a15,1p,e9.2,a8)') ' Minimal size:',a1,' microns'
   if (szds.eq.3) then
    write(12,'(a22,1p,e9.2,a8)')' Characteristic size:',a2,' microns'
   else
    write(12,'(a15,1p,e9.2,a8)')' Maximal size:',a2,' microns'
   end if
  end if
  write(12,*)' --------------------------------------------'

!-------------------------------------------------------------------------

  return
end subroutine WriteOut
!***********************************************************************

!***********************************************************************
subroutine getfs(xx,nm,flag,str)
!=======================================================================
! This subroutine writes number xx to a string str according to a format
! f?.nm. Here ? stands for the number of needed places. A blank is
! inserted at the beginning, and for flag.NE.1 another one if number is
! positive. If xx<0 second blank is replaced by '-'. For example, for
! flag=0 and xx = -0.1234D+02, calling this subroutine with nm=1 will
! result in str = ' -12.3', while xx = 0.0123 with nm=3 gives '  0.012'.
! If flag=1 minus will be ignored, for example xx = -0.1234D+02 and nm=1
! will result in str = ' 12.3',                        [Z.I., Nov. 1996]
!=======================================================================
  implicit none

  character ch
  character*(*) str
  integer  flag, nm, db, i, d(20), j, k, dnmp1
  double precision xx, x, rest
  !----------------------------------------------------------------------
  do i = 1, len(str)
     str(i:i) = ' '
  end do
  x = xx
  str(1:1) = ' '
  i = 2
  if (flag.ne.1) then
     if (x.lt.0.0d0) then
        str(i:i) = '-'
     else
        str(i:i) = ' '
     end if
     i = i + 1
  end if
  if (x.lt.0.0d0) x = -x
  ! First check if x will have to be rounded up
  ! Find (nm+1)-th decimal digit
  dnmp1 = int(x*10.0d0**(nm+1)-int(x*10.0d0**nm)*10.0d0)
  if (dnmp1.ge.5) x = x + 1./10.0d0**nm
  if (x.ge.1.0) then
     ! Number of digits before the decimal sign
     db = int(log10(x)) + 1
     ! Copy all these digits to str
     do j = 1, db
        rest = x
        if (j.gt.1) then
           do k = 1, j-1
              rest = rest - d(k)*10.0d0**(db-k)
           end do
        end if
        d(j) = int(rest/10.0d0**(db-j))
        write(ch,'(I1)')d(j)
        str(i:i) = ch
        i = i + 1
     end do
     rest = rest - d(db)
     if (nm.gt.0) then
        str(i:i) = '.'
        i = i + 1
     end if
  else
     str(i:i) = '0'
     i = i + 1
     if (nm.gt.0) then
        str(i:i) = '.'
        i = i + 1
     end if
     rest = x
  end if
  ! Now copy all nm remaining decimal digits to str
  if (nm.gt.0) then
     do j = 1, nm
        d(j) = int(rest*10.0d0**j)
        if (j.gt.1) then
           do k = 1, j-1
              d(j)=d(j)-int(d(k)*10.0d0**(j-k))
           end do
        end if
        write(ch,'(i1)')d(j)
        str(i:i) = ch
        i = i + 1
     end do
  end if
  !--------------------------------------------------------------------
  return
end subroutine getfs
!***********************************************************************

!***********************************************************************
subroutine prHeader(nLines,filenm)
!=======================================================================
! This subroutine prints headers of input data files in fname.out
!=======================================================================
  implicit none

  character filenm*(*), line*80
  integer nLines, i

!----------------------------------------------------------------------

  open(28, file=filenm, status = 'old')
  do i = 1, nLines-1
   read(28,'(a)') line
   write(12,'(a2,a80)') '  ', line
  end do
  rewind(28)
  close(28)
!----------------------------------------------------------------------
  return
end subroutine prheader
!***********************************************************************

! ***********************************************************************
SUBROUTINE LINE(com,typ,unt)
  ! =======================================================================
  ! This subroutine writes a line into file open as unt. For type = 1
  ! the line is '---', and for type = 2 '==='.If com=1 a comment sign # is
  ! added in the beginning (this is when line is used in file headers)
  ! =======================================================================
  IMPLICIT none
  INTEGER com, typ, unt
  ! -----------------------------------------------------------------------
  IF(typ.EQ.1) THEN
     IF(com.eq.1) THEN
        write(unt,'(a50)')'# ------------------------------------------------'
     ELSE
        write(unt,*)'--------------------------------------------------'
     END IF
  ELSE
     IF(com.eq.1) THEN
        write(unt,'(a50)')'# ================================================'
     ELSE
        write(unt,*)'=================================================='
     END IF
  END IF
  ! -----------------------------------------------------------------------
  RETURN
END SUBROUTINE LINE
! ***********************************************************************

!***********************************************************************
subroutine OPPEN(model,rootname,length)
!=======================================================================
! This subroutine prints the results out.              [Z.I., Feb. 1996]
!=======================================================================

  use common
  implicit none
  character ch5*5, rootname*(*), fname*235
  character*72 header1, s3, s4
  integer model, length, i
  !---------------------------------------------------------------------

  ! set up the status indicators
  iError = 0
  iWarning = 0

  if (model.eq.1) iCumm = 0
  ! the following files pertain to all models and are open if model.eq.1
  if (model.eq.1) then
     ! the header to output file *.out is moved to prout [MN,Sep'99]
     ! open file for crosssections
     call attach(rootname,length,'.Xsec',fname)
     open(855,file=fname,status='unknown')
     ! open file for point spread function
     if (iPsf.ne.0.and.psftype.lt.3) then
        ! wavelength dependent psf are also printed out
        call attach(rootname,length,'.psf',fname)
        open(23,file=fname,status='unknown')
        header1 = '    offset'
        write(23,'(a10,20f10.2)')header1,(LambdaOut(i),i=1,nLambdaOut)
     end if
     ! spectra for all models in one file '*.stb' if flag=1
     if(iA.eq.1) then
        call attach(rootname,length,'.stb',fname)
        open(15,file=fname,status='unknown')
     end if
     ! all radial profiles in one file '*.rtb' if flag=1
     if(iB.eq.1) then
        call attach(rootname,length,'.rtb',fname)
        open(16,file=fname,status='unknown')
     end if
     ! all imaging files in  '*.itb' if flag=1
     if (abs(iC).eq.1) then
        call attach(rootname,length,'.itb',fname)
        open(17,file=fname,status='unknown')
     end if
     ! all J-files in  '*.jtb' if flag=1
     if (abs(iJ).eq.1) then
        call attach(rootname,length,'.jtb',fname)
        open(19,file=fname,status='unknown')
     end if
     ! all message files in  '*.mtb' if flag=1
     if(iX.eq.1) then
        call attach(rootname,length,'.mtb',fname)
        open(18,file=fname,status='unknown')
     end if
     ! open the file for error vs. iterations in '*.err'
     if(iInn.eq.1) then
        call attach(rootname,length,'.err',fname)
        open(38,file=fname,status='unknown')
        write(38,'(a12,a235)') 'Input file: ',rootname
     end if
     ! for private rdw option
     if(rdwpr) then
        call attach(rootname,length,'.rdw',fname)
        open(66,file=fname,status='unknown')
     end if
     ! end if for model=1
  end if
!-------------------------------------------------------------
! the following files are open for every model

! (the headers for .s## and .r## files are moved to prout, MN)
! open the spectrum file rootname.s##  (## = model number)
  if(iA.gt.1) then
     write(ch5,'(a2,i3.3)') '.s', model
     call attach(rootname,length,ch5,fname)
     open(15,file=fname,status='unknown')
     if(slb) then
        if(iA.eq.3) then
           write(ch5,'(a2,i3.3)') '.z', model
           call attach(rootname,length,ch5,fname)
           open(25,file=fname,status='unknown')
        end if
     end if
  end if
  ! open the file rootname.r## (y-dependent quantities)
  if(iB.gt.1) then
     write(ch5,'(a2,i3.3)') '.r', model
     call attach(rootname,length,ch5,fname)
     open(16,file=fname,status='unknown')
  end if
  ! open the file rootname.i## (surface brightness)
  if(abs(iC).gt.1) then
     write(ch5,'(a2,i3.3)') '.i', model
     call attach(rootname,length,ch5,fname)
     open(17,file=fname,status='unknown')
  end if
  ! open the file rootname.c## (convolved images) for flag ic<0
  ! if ic=-3 - in a separate file fname.c##
  if(iC.eq.-3.and.iPsf.gt.0) then
     write(ch5,'(a2,i3.3)') '.c', model
     call attach(rootname,length,ch5,fname)
     open(21,file=fname,status='unknown')
  end if
  ! open the file rootname.v## (visibility curves)
  if((abs(iC).eq.3).and.(iV.gt.0)) then
     write(ch5,'(a2,i3.3)') '.v', model
     call attach(rootname,length,ch5,fname)
     open(22,file=fname,status='unknown')
  end if
  ! open the output file rootname.m##
  if (iX.gt.1) then
     write(ch5,'(a2,i3.3)') '.m', model
     call attach(rootname,length,ch5,fname)
     open(18,file=fname,status='unknown')
  end if
  ! open the output file rootname.j##
  if (iJ.gt.1) then
     write(ch5,'(a2,i3.3)') '.j', model
     call attach(rootname,length,ch5,fname)
     open(19,file=fname,status='unknown')
  end if
  !--------------------------------------------------------------------
  return
end subroutine OPPEN
!***********************************************************************

!***********************************************************************
subroutine OccltMSG(us)
!=======================================================================
! Prints a message informing the user about the min Teff required to
! neglect occultation by the central source.
!=======================================================================
  use common
  implicit none
  character*10 tstrg
  integer iL
  double precision qaux(npL), qaux2(npL), res1, res2, Te_min_loc,  &
       us(npL,npY), mx, psitn, Planck, xp
!-----------------------------------------------------------------------
  ! Estimate min Teff required to neglect occultation (eq.(5) in Manual):
  write(18,*) 'Tsub(',ifidG,')=', Tsub(ifidG)
  ! write(18,*) '  lambda(iL)  SigmaA(1,iL) Planck(xP)   xP'
  do iL = 1, nL
     qaux(iL) = sigmaA(ifidG,iL)*us(iL,1)/lambda(iL)
     xP = 14400.0d0/Tsub(ifidG)/lambda(iL)
     qaux2(iL) = sigmaA(1,iL)*Planck(xP)/lambda (iL)
  end do
  call Simpson(npL,1,nL,lambda,qaux,res1)
  call Simpson(npL,1,nL,lambda,qaux2,res2)
  ! approximate psi for opt.thin case:
  psitn = res1/res2
  mx = Tsub(ifidG)*sqrt(sqrt(4.0d0/psitn))
  if(Tsub(ifidG).lt.mx) then
     Te_min_loc = 2.0d0*mx
  else
     Te_min_loc = 2.0d0*Tsub(ifidG)
  end if
  call getfs(Te_min_loc,0,1,tstrg)
  write(12,*) ' ====================================================  '
  write(12,*) ' For compliance with the point-source assumption, the'
  write(12,*) ' following results should only be applied to sources '
  write(12,'(a37,a5,a3)') '  whose effective temperature exceeds ',Tstrg, ' K.'
  write(12,*) ' ===================================================='
!-----------------------------------------------------------------------
  return
end subroutine OccltMSG
!***********************************************************************
