PROGRAM DUSTY
  USE COMMON
  IMPLICIT NONE
  INTEGER :: clock_rate, clock_start, clock_end, io_status, lpath
  CHARACTER(len=3) :: suffix
  CHARACTER(len=4) :: verbosity
  CHARACTER(len=235) :: dustyinpfile, path, apath, stdf(7)

  !-------------------------------------------------------
  ! **************************
  ! *** ABOUT THIS VERSION ***
  ! **************************
  ! version= '4.00' set in common as parameter

  ! get lambda grid
  CALL ReadLambda()
  IF (error.ne.0) THEN 
     PRINT*,'something wrong with lambda grid!'
     STOP
  END IF

  ! timing
  CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) ! Find the rate
  CALL SYSTEM_CLOCK(COUNT=clock_start)

  ! get path of input file, and desired verbosity, as args from command line
  CALL GETARG(1,dustyinpfile)
  CALL GETARG(2,verbosity) ! verbosity if optional, since a default exists
  if (LEN_TRIM(verbosity) == 0) verbosity = "2"  ! and the default verbosity level is "2"
  read(verbosity(1:2),'(I1)') iVerb ! make it an integer
  if (iVerb > 2) iVerb = 2 ! largest implemented value is 2

  ! interpret the input file path, and run the calcs if all os OK
  IF (TRIM(dustyinpfile).eq."") THEN
     PRINT*, "No input file name found on command line."
     PRINT*,'DUSTY STOPPED'
     STOP
  ELSE
     suffix = dustyinpfile(LEN(TRIM(dustyinpfile))-2:) ! the last 3 chars of path

     ! If suffix is 'mas', we have multi-model DUSTY input file. Then
     ! read each line, execute the requested model, and continue with
     ! the next line.
     IF (suffix .eq. 'mas') THEN
        PRINT*, "Found master input file ",TRIM(dustyinpfile), &
             " on on command line."
        OPEN(13,file=trim(dustyinpfile),status='old')
        READ(13,'(a)',iostat=io_status) apath
        DO WHILE (io_status.ge.0)
           CALL clean(apath,path,lpath)
           call rundusty(path,lpath)
           READ(13,'(a)',iostat=io_status) apath
        END DO
        CLOSE(13)

     ! If suffix is 'inp', this is a single-model DUSTY input file.
     ELSE IF (suffix .eq. 'inp') THEN
        PRINT*, "Found normal input file ",TRIM(dustyinpfile), &
             " on on command line."
        apath = dustyinpfile(1:)
        CALL clean(apath,path,lpath)
        call rundusty(path,lpath)

     ! everything else is considered an invalid input file
     ELSE
        PRINT*,'WARNING NOT A PROPER INPUT FILE NAME!!!!'
        PRINT*,'DUSTY STOPPED'
        STOP
     END IF

  END IF

  CALL SYSTEM_CLOCK(COUNT=clock_end) ! stop the timing
  print*,'ellapsed time:',(clock_end-clock_start)/clock_rate,'[s]'

END PROGRAM DUSTY


!***********************************************************************
subroutine rundusty(path,lpath)
!=======================================================================
! This subroutine a single DUSTY process, independently of whether it
! was specified in a single-model .inp file or in a multi-model .mas
! file. All processing of paths and of command-line arguments is done
! in the main code, before this subroutine is invoked.
!                                                      [RN, August 2013]
!=======================================================================
  use common, only: iVerb, error
  implicit none
  integer :: lpath, empty, GridType, Nmodel
  double precision :: tau1, tau2
  double precision, allocatable :: tau(:)
  character(len=235) :: path, nameIn, nameOut, stdf(7)

  interface
     subroutine GetTau(tau1,tau2,GridType,Nmodel,tau)
       integer Nmodel, GridType
       double precision  tau1, tau2, tau(:)
     end subroutine GetTau
     subroutine Kernel(path,lpath,tau,Nmodel)
       integer Nmodel, lpath
       double precision tau(:)
       character(len=235) path
     end subroutine Kernel
  end interface

  call alloc_mem()
  call alloc_mem_nL()

  path = path(1:LEN(TRIM(path))-4)
  lpath = lpath-4

  if (empty(path).ne.1) then
     call attach(path,lpath,'.inp',nameIn)
     call attach(path,lpath,'.out',nameOut)
     call Input(nameIn,nameOut,tau1,tau2,GridType,Nmodel)

     if (iverb.gt.0) then
        print*,'working on input file: ', trim(nameIn)
        if (iVerb.ge.2) print*,'Done with reading input'
     endif

     if (error.ne.3) then 
        if (error.eq.0) then 
           if(allocated(tau)) deallocate(tau)
           allocate(tau(Nmodel))
           tau = 0
           call GetTau(tau1,tau2,GridType,Nmodel,tau)
           if (iVerb.ge.2) print*,'Done with GetTau'
           call Kernel(path,lpath,tau,Nmodel)
        end if
     else
        print*,
     end if

  end if

  call dealloc_mem()

end subroutine rundusty
!***********************************************************************


!***********************************************************************
subroutine ReadLambda()
!=======================================================================
! This subroutine reads and checks that the wavelength grid satisfies
! certain conditions described in the Manual (all wavelengths are given
! in microns). If everything went fine it returns error = 0, and
! fills the wavelength grid in array lambda.
!                                      [ZI,Feb'96; MN,Apr'99; FH,Jan'12]
!=======================================================================
  use common
  implicit none
  integer  iL
  double precision RDINP
  character str*235
  logical Equal
!-----------------------------------------------------------------------
  Equal = .true.
  error = 0
  ! first open the file with lambda grid
  open(4, file='data/lambda_grid.dat', status = 'old')
  call skip_header(4)
  nL = RDINP(Equal,4,6)
  allocate(lambda(nL))
  lambda = 0
  ! initialize lambda array
  do iL = 1, nL
     read(4,*,end=99) lambda(iL)
  end do
99 close(4)
  call sort(lambda,nL)
  do iL = 2, nL
     if (lambda(iL)/lambda(iL-1).gt.1.51d0) then
        write(*,*)' ***************** WARNING!  *******************'
        write(*,*)' the ratio of two consecutive wavelengths in the'
        write(*,*)' grid has to be no bigger than 1.5. you have    '
        write(*,'(2(a4,1p,e8.2))') '    ',lambda(iL)/lambda(iL-1),' at ', lambda(iL)
        write(*,*)' correct this and try again!                    '
        write(*,*)' ***********************************************'
        error = 1
     end if
  end do
  return
end subroutine ReadLambda
!***********************************************************************

!***********************************************************************
subroutine GetTau(tau1,tau2,GridType,Nmodel,tau)
!=======================================================================
! This subroutine generates total optical depth TAUtot.
!                                                      [Z.I., Mar. 1996]
!=======================================================================
  use common
  implicit none
  integer model, Nmodel, iL, GridType
  double precision  q, TAU1, TAU2, tau(:), taumax
!-----------------------------------------------------------------------
  taumax = 0.0d0
  do model = 1, Nmodel
     if(GridType.eq.3) then
        tau(model) = tauIn(model)
     else
        ! Calculate taufid for given model
        if (model.eq.1) then
           tau(model) = tau1
        else
           if (model.eq.Nmodel) then
              tau(model) = tau2
           else
              if(GridType.eq.1) then
                 q =  (tau2 - tau1)/(Nmodel-1.0d0)
                 tau(model) = tau1 + q*(model-1.0d0)
              else
                 q = dexp(dlog(tau2/tau1)/(Nmodel-1.0d0))
                 tau(model) = tau1*q**(model-1.0d0)
              end if
           end if
        end if
     end if
  end do
!-----------------------------------------------------------------------
  return
end subroutine GetTau
!***********************************************************************

subroutine alloc_mem_nL()
use common
implicit none 
integer :: iL, iY, iR, iG, iP
  ! Memory allocation----------------------------------
  allocate(fsL(nL,npY))
  fsL = 0
  allocate(fsR(nL,npY))
  fsR = 0
  allocate(TAUtot(nL))
  TAUtot = 0 
  allocate(TAUslb(nL,npY))
  TAUslb = 0 
  allocate(omega(npG+1,nL))
  omega = 0
  allocate(Utot(nL,npY))
  Utot = 0
  allocate(Utot_old(nL,npY))
  Utot_old = 0
  allocate(Ude(nL,npY))
  Ude = 0
  allocate(Uds(nL,npY))
  Uds = 0
  allocate(fds(nL,npY))
  fds = 0
  allocate(fde(nL,npY))
  fde = 0
  allocate(ftot(nL,npY))
  ftot = 0
  allocate(IstR(nL))
  IstR = 0
  allocate(SLBIntm(npR,nL))
  SLBIntm = 0
  allocate(SLBIntp(npR,nL))
  SLBIntp = 0
  allocate(Intens(nL,npP+2))
  Intens = 0
end subroutine alloc_mem_nL

subroutine alloc_mem()
use common
implicit none 
  integer :: iG, iY,iP
  ! Memory allocation----------------------------------
  allocate(Y(npY))
  Y = 0
  allocate(Yprev(npY))
  Yprev = 0
  allocate(Plast(npY))
  Plast = 0
  allocate(ETAdiscr(npY))
  ETAdiscr = 0
  allocate(fsLbol(npY))
  fsLbol = 0 
  allocate(fsRbol(npY))
  fsRbol = 0
  allocate(fsbol(npY))
  fsbol = 0
  allocate(Jext(npY))
  Jext = 0
  allocate(RPr(npY))
  RPr = 0
  allocate(tauF(npY))
  tauF = 0 
  allocate(eps(npY))
  eps = 0
  allocate(ugas(npY))
  ugas = 0
  allocate(Gamma(npY))
  Gamma = 0
  allocate(fbol(npY))
  fbol = 0
  allocate(ubol(npY))
  ubol = 0
  allocate(fpbol(npY))
  fpbol = 0
  allocate(fmbol(npY))
  fmbol = 0
  allocate(qF(npY))
  qF = 0
  allocate(tauFdyn(npY))
  tauFdyn = 0
  allocate(ETAzp(npP,npY))
  ETAzp = 0 
  allocate(abund(npG+1,npY))
  abund = 0
  allocate(Td_old(npG+1,npY))
  Td_old = 0
  allocate(Td(npG+1,npY))
  Td = 0
  allocate(vrat(npG+1,npY))
  vrat = 0
  allocate(rg(npG+1,npY))
  rg = 0
  allocate(P(npP))
  P = 0
  allocate(iYfirst(npP))
  iYfirst = 0
  allocate(YPequal(npP))
  YPequal = 0
  allocate(IntOut(20,npP+2))
  IntOut = 0
  allocate(bOut(npP+2))
  bout = 0
  allocate(tauZout(npP+2))
  tauZout = 0
  allocate(destroyed(npG,npY))
  destroyed = 0
end subroutine alloc_mem

subroutine dealloc_mem()
use common
implicit none 
  ! Memory allocation----------------------------------
  deallocate(Y)
  deallocate(Yprev)
  deallocate(P)
  deallocate(iYfirst)
  deallocate(YPequal)
  deallocate(Plast)
  deallocate(ETAdiscr)
  deallocate(ETAzp)
  deallocate(abund)
  deallocate(fsL)
  deallocate(fsR)
  deallocate(fsLbol)
  deallocate(fsRbol)
  deallocate(fsbol)
  deallocate(TAUtot)
  deallocate(TAUslb)
  deallocate(omega)
  deallocate(Td)
  deallocate(Td_old)
  deallocate(Utot)
  deallocate(Utot_old)
  deallocate(Ude)
  deallocate(Uds)
  deallocate(fds)
  deallocate(fde)
  deallocate(ftot)
  deallocate(Jext)
  deallocate(RPr)
  deallocate(tauF)
  deallocate(eps)
  deallocate(ugas)
  deallocate(Gamma)
  deallocate(vrat)
  deallocate(SLBIntm)
  deallocate(SLBIntp)
  deallocate(IntOut)
  deallocate(IstR)
  deallocate(bOut)
  deallocate(tauZout)
  deallocate(fbol)
  deallocate(ubol)
  deallocate(qF)
  deallocate(rg)
  deallocate(fpbol)
  deallocate(fmbol)
  deallocate(Intens)
  deallocate(destroyed)
  deallocate(tauFdyn)
  if (allocated(theta)) deallocate(theta)
end subroutine dealloc_mem

subroutine skip_header(iunit)
!Get past Header lines of opened file on unit # iunit
!The end of headers is marked by a line that contains just ">"
  integer iunit,io
  character*128 header
  logical found_ge
  found_ge = .false.
  header = "XXX"
  do while((header(1:1).ne.">").and.(io.ne.-1))
     read(iunit, '(a)',iostat=io) header
     if (header(1:1).eq.">") found_ge = .true.
  end do
  if (.not.found_ge) rewind(iunit)
  return
end subroutine skip_header
