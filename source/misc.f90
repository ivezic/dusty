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



!************************************************************************
!!$
!!$!***********************************************************************
!!$subroutine fileMSG(fname,strg)
!!$!=======================================================================
!!$! Prints a message in *.out file in case of error opening the user
!!$! supplied files.
!!$!=======================================================================
!!$  implicit none
!!$  character aux*235, strg*(*), fname*(*)
!!$  integer length, empty
!!$!-----------------------------------------------------------------------
!!$1 read(1,'(a)') aux
!!$  if (empty(aux).eq.1) goto 1
!!$  call clean(aux,fname,length)
!!$
!!$  open(10, err=100, file=fname, status='old')
!!$  close(10)
!!$  return
!!$
!!$100 write(12,*)' *** Fatal error in dusty! **************************'
!!$  write(12,*)' File with the ',strg
!!$  write(12,'(a2,a)')'  ',fname
!!$  write(12,*)' is missing ?!'
!!$  write(12,*)' ****************************************************'
!!$!  close(12)
!!$!-----------------------------------------------------------------------
!!$  stop
!!$end subroutine fileMSG
!!$!***********************************************************************
!!$
!!$
!!$
!!$
!!$!***********************************************************************
!!$subroutine ChkAngle(angle)
!!$!=======================================================================
!!$! Checks if the input angles are in [0,85] degrees interval.  [MN,2005]
!!$!=======================================================================
!!$  implicit none
!!$  integer error
!!$  double precision angle
!!$!-----------------------------------------------------------------------
!!$  if (angle.gt.88.0d0) then
!!$     write(12,*)' ********** Message from input: ***************'
!!$     write(12,*)' * slab illumination angles have to be in the *'
!!$     write(12,*)' * range [0,88] degrees. setting the maximum  *'
!!$     write(12,*)' * illumination angle to 88 degrees.          *'
!!$     write(12,*)' **********************************************'
!!$     angle = 88.0d0
!!$  else
!!$     if (angle.lt.0.0d0) then
!!$        write(12,*)' ********** Message from input: ***************'
!!$        write(12,*)' * slab illumination angles have to be in the *'
!!$        write(12,*)' * range [0,85] degrees. setting the miminum  *'
!!$        write(12,*)' * illumination angle to 0.0 degrees.         *'
!!$        write(12,*)' **********************************************'
!!$     end if
!!$  end if
!!$!-----------------------------------------------------------------------
!!$  return
!!$end subroutine ChkAngle
!!$!***********************************************************************
!!$
!!$
!!$!***********************************************************************
!!$subroutine ReadSpectar(lambdas,Llamstar,Lstar,nLs,is,error)
!!$!=======================================================================
!!$! Reads the source spectrum from a file. This was part of Subroutine Star,
!!$! separated for clarity. is=1 for central s-ce, is=2 for external [MN'01]
!!$! =======================================================================
!!$  use common
!!$  implicit none
!!$
!!$  character*235 line
!!$  integer ios1, iLs, nLs, error, nLambdam, Nis, is
!!$! nLambdam is the max number entries for a user supplied stellar spectrum
!!$  parameter (nLambdam = 10000, Nis = 2)
!!$  double precision lambdas(nLambdam), Llamstar(nLambdam),lLs(nLambdam), &
!!$       Ls(nLambdam), Lstar, a, b
!!$! -----------------------------------------------------------------------
!!$
!!$! For  startyp.ge.4 stellar spectrum is read from the file 'namestar'
!!$! which is unit=3 (unit=1 is the input file)
!!$! is=1 for the enclosed source, is=2 for the external shell illumination
!!$! possible problems with the file are checked in sub inp_rad
!!$    open(3,file=namestar(is),status='old')
!!$      rewind(3)
!!$      do iLs = 1, 3
!!$        read(3,'(a235)',err=998) line
!!$      end do
!!$      ios1 = 0
!!$      iLs = 0
!!$      do while (ios1.ge.0)
!!$        read(3,*,end=900,err=998,iostat=ios1) a, b
!!$        if(ios1.ge.0) then
!!$          iLs = iLs + 1
!!$          lambdas(iLs) = a
!!$          if (a.le.0.0) goto 998
!!$!         Llamstar is always f_lambda
!!$!         if startyp.eq.4 then file gives lambda*f_lambda
!!$          if (startyp(is).eq.4) Llamstar(iLs) = b/a
!!$!         if startyp.eq.5 then file gives f_lambda
!!$          if (startyp(is).eq.5) Llamstar(iLs) = b
!!$!         if startyp.eq.6 then file gives lnu=lambda**2*f_lambda
!!$          if (startyp(is).eq.6) Llamstar(iLs) = b/(a**2)
!!$        end if
!!$      end do
!!$900   close(3)
!!$      if (iLs.lt.2) goto 998
!!$      nLs = iLs
!!$!     if input wavelengths in descending order turn them around
!!$      if (lambdas(1).gt.lambdas(2)) then
!!$        do iLs = 1, nLs
!!$          lLs(iLs) = lambdas(iLs)
!!$          Ls(iLs) =  Llamstar(iLs)
!!$        end do
!!$        do iLs = 1, nLs
!!$          lambdas(iLs) = lLs(nLs+1-iLs)
!!$          Llamstar(iLs) = Ls(nLs+1-iLs)
!!$        end do
!!$      end if
!!$!     Get the scale, Lstar, of the stellar spectrum
!!$      call Simpson(nLambdam,1,nLs,lambdas,Llamstar,Lstar)
!!$      error = 0
!!$      goto 999
!!$998   write(12,*)' *********** INPUT ERROR *************************'
!!$      write(12,*)' The file with spectral shape of external radiation:'
!!$      write(12,'(a2,a100)')'  ',namestar(1)
!!$      write(12,*)' is missing or not properly formatted?!'
!!$      write(12,*)' ***************************************************'
!!$      error = 3
!!$! -----------------------------------------------------------------------
!!$
!!$999 return
!!$end subroutine ReadSpectar
!!$!***********************************************************************
!!$
!!$!***********************************************************************
!!$subroutine WriteOut(var1,var2,var3,nG,nameQ,nameNK)
!!$!=======================================================================
!!$! WriteOut prints in fname.out all input parameters,  read before density distribution
!!$! type.
!!$!=======================================================================
!!$
!!$  use common
!!$  implicit none
!!$
!!$  integer is, iG, nG, i, length
!!$  double precision var1, var2, var3
!!$  character*72 strpow, aux, src, chaux*3
!!$  character*(*) nameQ(npG), nameNK(10)
!!$  logical first
!!$  first = .true.
!!$!-------------------------------------------------------------------------
!!$  is = 1
!!$15 if (SLB) then
!!$   if (is.eq.1) then
!!$    src = 'Left-side source spectrum described by'
!!$   else
!!$    src = 'Right-side source spectrum described by'
!!$   end if
!!$  else
!!$   if (is.eq.1) then
!!$    src = 'Central source spectrum described by'
!!$   else
!!$    src = 'External source spectrum described by'
!!$   end if
!!$  end if
!!$  call Clean(src, aux, length)
!!$
!!$  if(Left.eq.0.and.is.eq.1) then
!!$   write(12,*) ' No central source.'
!!$  else
!!$! #1: black body(ies) for startyp=1
!!$   if (startyp(is).eq.1) then
!!$    if (nBB(is).gt.1) then
!!$     call ATTACH(aux, length, ' ', src)
!!$! multiple black bodies
!!$     write(12,'(a2,a37,i2,a13)')'  ', src, nBB(is),' black bodies'
!!$     write(12,'(a27)')' with temperatures (in K):'
!!$     write(12,'(2x,1p,10e10.3)')(Tbb(is,i),i=1,nBB(is))
!!$     write(12,'(a42)')' and relative luminosities, respectively:'
!!$     write(12,'(1p,10e10.1)')(rellum(is,i),i=1,nBB(is))
!!$    else
!!$! for a single black body:
!!$     call ATTACH(aux,length,' a black body',src)
!!$     write(12,'(a2,a)') '  ',src
!!$     write(12,'(a19,1p,e10.3,a2)')' with temperature:',Tbb(is,1),' K'
!!$    end if
!!$   end if
!!$
!!$! #2: Engelke-Marengo function for startyp=2
!!$   if (startyp(is).eq.2) then
!!$    call ATTACH(aux, length,' Engelke-Marengo function', src)
!!$    write(12,'(a2,a)') '  ',src
!!$    write(12,'(a13,1p,e10.3,a16)')' with Teff =',Tbb(is,1), ' K and depth of'
!!$    write(12,'(a30,F6.1,a2)')' the SiO absorption feature =', xSiO,' %'
!!$   end if
!!$
!!$! #3: power-law(s) for startyp=3
!!$   if (startyp(is).eq.3) then
!!$    if (Nlamtr(is).gt.0) then
!!$     call ATTACH(aux,length,' power law:',src)
!!$     write(12,'(a2,a)') '  ',src
!!$     write(12,*)'    lambda      k'
!!$     do i = 1, Nlamtr(is)
!!$      write(12,'(1x,1p,e10.3)')lamtr(is,i)
!!$      write(12,'(11x,1p,e10.3)')klam(is,i)
!!$     end do
!!$     write(12,'(1x,1p,e10.3)')lamtr(is,Nlamtr(is)+1)
!!$    else
!!$     write(12,*)' Input data for the source spectrum is not good.'
!!$     write(12,*)' Changed to a 10000 K black body'
!!$    end if
!!$   end if
!!$
!!$! spectrum from a file for startyp=4,5,6
!!$   if (startyp(is).ge.4.and.startyp(is).le.6) then
!!$    if (is.eq.1) then
!!$     write(12,*)' Stellar spectrum supplied from file:'
!!$    else
!!$     write(12,*)' External spectrum supplied from file:'
!!$    end if
!!$    write(12,'(a2,a100)') '  ',nameStar(is)
!!$    call PrHeader(3,nameStar(is))
!!$   end if
!!$  end if
!!$  write(12,*)' --------------------------------------------'
!!$  if(first) then
!!$! if there is a second source go back to read its parameters
!!$   if(Right.gt.0) then
!!$! repeat printout of source info for the external radiation
!!$! its index is is=2
!!$    is = 2
!!$    first = .false.
!!$    goto 15
!!$   end if
!!$  end if
!!$  is = 1
!!$! -----------------------------------------------------
!!$! Boundary Condition data: typEntry(is) and value
!!$  if(Left.gt.0) then
!!$   if(typEntry(1).eq.5) then
!!$     if (SLB) then
!!$      write(12,'(a45,1p,e9.2,a,i2)')  &
!!$           ' Dust temperature on the slab left boundary:', Tsub(ifidG),' K - Grain:',ifidG
!!$     else
!!$      write(12,'(a41,1p,e9.2,a,i2)')  &
!!$           ' Dust temperature on the inner boundary:', Tsub(ifidG),' K - Grain:',ifidG
!!$     end if
!!$   else if (typEntry(1).eq.1) then
!!$    if (slb) then
!!$     if (startyp(1).gt.3) then
!!$      write(12,'(a33,1p,e9.2,a6)') &
!!$           ' Flux at the slab left boundary:', dilutn*var1,' W/m^2'
!!$     else
!!$      write(12,'(a33,1p,e9.2,a6)') &
!!$           ' Flux at the slab left boundary:', Ji*4*pi,' W/m^2'
!!$     end if
!!$! if input spectrum is from a file its bol.flux is calculated
!!$! and then normalized with the value of 'dilutn' from the input file.
!!$    elseif(sph) then
!!$     if (startyp(1).gt.3) then
!!$      write(12,'(a29,1p,e9.2,a6)') &
!!$           ' Flux at the inner boundary:', dilutn*var1,' W/m^2'
!!$     else
!!$      write(12,'(a29,1p,e9.2,a6)') &
!!$           ' Flux at the inner boundary:', var1,' W/m^2'
!!$     end if
!!$    end if
!!$   else if (typEntry(1).eq.2) then
!!$    write(12,'(a20,1p,e9.2,a17,1p,e9.3,a3)') &
!!$         ' Source luminosity ', var1,' Lo and distance ', var2, ' cm'
!!$
!!$   else if (typEntry(1).eq.3.or.typEntry(1).eq.4) then
!!$    write(12,'(a38,1p,e9.2,a6)') &
!!$         ' Mean intensity at the inner surface:', var3,' W/m^2'
!!$  end if
!!$ end if
!!$ if(Right.gt.0) then
!!$  if (sph) then
!!$   if (typEntry(1).eq.1) then
!!$    if (startyp(2).gt.3) then
!!$     write(12,'(a49,1p,e9.2)')  &
!!$          '  Normalization factor of the external radiation:', dilutn
!!$     write(12,'(a54,1p,e9.2,a6)')  &
!!$          '  Mean intensity from the file after renormalization:',var3,' W/m^2'
!!$    else
!!$     write(12,'(a45,1p,e9.2)') &
!!$          '  Dilution factor of the external radiation:', dilutn
!!$
!!$    end if
!!$   else
!!$    write(12,'(a44,1p,e9.2,a6)')  &
!!$         '  Mean intensity of the external radiation:', var3,' W/m^2'
!!$   end if
!!$  end if
!!$ end if
!!$ write(12,*)' --------------------------------------------'
!!$! -----------------------------------------------------
!!$! 2) DUST PROPERTIES
!!$!  2.1 Chemical Composition
!!$  if (top.lt.3) then
!!$   write(12,*)' Abundances for supported grains:'
!!$   write(12,*)' Sil-Ow Sil-Oc Sil-DL grf-DL amC-Hn SiC-Pg'
!!$   write(12,'(6f7.3)')(xC(i),i=1,3),xC(4)+xC(5),(xC(i),i=6,7)
!!$   if (top.eq.2) then
!!$    write(12,*)' Abundances for user supplied grains:'
!!$    write(12,'(i6,9i7)')(i,i=1,Nfiles)
!!$    write(12,'(10f7.3)')(xCuser(i),i=1,Nfiles)
!!$    write(12,*)' User supplied n and k from:'
!!$    do i = 1, Nfiles
!!$     write(12,'(a2,i1,a2,a70)')'  ',i,') ',nameNK(i)
!!$    end do
!!$   end if
!!$! user supplied cross-sections:
!!$  else
!!$   do iG = 1, nG
!!$    write(12,*)' Optical properties from file(s):'
!!$    write(12,'(a2,a70)')'  ',nameQ(iG)
!!$    call PrHeader(3,nameQ(iG))
!!$   end do
!!$  end if
!!$  write(12,*)' Sublimation Temperature(s):'
!!$  do iG=1,nG
!!$     write(12,'(a8,i2,a3,f10.3)') '  Grain(',iG,'): ',Tsub(iG) 
!!$  end do
!!$! 2.2 Grain size distribution
!!$  if (top.ne.3) then
!!$   if (szds.eq.3) then
!!$    chaux = 'KMH'
!!$   else
!!$    chaux = 'MRN'
!!$   end if
!!$   write(12,'(a2,a3,a19)')'  ',chaux,'size distribution:'
!!$   call getfs(qsd,1,0,strpow)
!!$   write(12,'(a15,a5)')'      Power q:',strpow
!!$   write(12,'(a15,1p,e9.2,a8)') ' Minimal size:',a1,' microns'
!!$   if (szds.eq.3) then
!!$    write(12,'(a22,1p,e9.2,a8)')' Characteristic size:',a2,' microns'
!!$   else
!!$    write(12,'(a15,1p,e9.2,a8)')' Maximal size:',a2,' microns'
!!$   end if
!!$  end if
!!$  write(12,*)' --------------------------------------------'
!!$
!!$!-------------------------------------------------------------------------
!!$
!!$  return
!!$end subroutine WriteOut
!!$!***********************************************************************
!!$
!!$
!!$
!!$
!!$
!!$!***********************************************************************
!!$subroutine OccltMSG(us)
!!$!=======================================================================
!!$! Prints a message informing the user about the min Teff required to
!!$! neglect occultation by the central source.
!!$!=======================================================================
!!$  use common
!!$  implicit none
!!$  character*10 tstrg
!!$  integer iL
!!$  double precision qaux(npL), qaux2(npL), res1, res2, Te_min_loc,  &
!!$       us(npL,npY), mx, psitn, Planck, xp
!!$!-----------------------------------------------------------------------
!!$  ! Estimate min Teff required to neglect occultation (eq.(5) in Manual):
!!$  write(18,*) 'Tsub(',ifidG,')=', Tsub(ifidG)
!!$  ! write(18,*) '  lambda(iL)  SigmaA(1,iL) Planck(xP)   xP'
!!$  do iL = 1, nL
!!$     qaux(iL) = sigmaA(ifidG,iL)*us(iL,1)/lambda(iL)
!!$     xP = 14400.0d0/Tsub(ifidG)/lambda(iL)
!!$     qaux2(iL) = sigmaA(1,iL)*Planck(xP)/lambda (iL)
!!$  end do
!!$  call Simpson(npL,1,nL,lambda,qaux,res1)
!!$  call Simpson(npL,1,nL,lambda,qaux2,res2)
!!$  ! approximate psi for opt.thin case:
!!$  psitn = res1/res2
!!$  mx = Tsub(ifidG)*sqrt(sqrt(4.0d0/psitn))
!!$  if(Tsub(ifidG).lt.mx) then
!!$     Te_min_loc = 2.0d0*mx
!!$  else
!!$     Te_min_loc = 2.0d0*Tsub(ifidG)
!!$  end if
!!$  call getfs(Te_min_loc,0,1,tstrg)
!!$  write(12,*) ' ====================================================  '
!!$  write(12,*) ' For compliance with the point-source assumption, the'
!!$  write(12,*) ' following results should only be applied to sources '
!!$  write(12,'(a37,a5,a3)') '  whose effective temperature exceeds ',Tstrg, ' K.'
!!$  write(12,*) ' ===================================================='
!!$!-----------------------------------------------------------------------
!!$  return
!!$end subroutine OccltMSG
!!$!***********************************************************************
