PROGRAM DUSTY
  USE COMMON
  IMPLICIT NONE
  INTEGER :: clock_rate, clock_start, clock_end, io_status, lpath
  INTEGER :: empty, GridType
  INTEGER :: Nmodel
  DOUBLE PRECISION :: RDINP, tau1, tau2
  DOUBLE PRECISION, allocatable :: tau(:)
  CHARACTER(len=4)   :: suffix,verbosity
  CHARACTER(len=235) :: dustyinpfile, path, apath, nameIn, &
       nameOut, stdf(7)
  INTEGER iL
  INTERFACE
     SUBROUTINE GetTau(tau1,tau2,GridType,Nmodel,tau)
       integer Nmodel, GridType
       double precision  TAU1, TAU2, tau(:)
     END SUBROUTINE GETTAU
     SUBROUTINE KERNEL(path,lpath,tau,Nmodel)
       integer  Nmodel, lpath
       double precision tau(:)
       character(len=235) path
     END SUBROUTINE KERNEL
  END INTERFACE
  !-------------------------------------------------------
  ! **************************
  ! *** ABOUT THIS VERSION ***
  ! **************************
  ! version= '4.00' set in common as parameter

  CALL ReadLambda()
  IF (error.ne.0) THEN 
     PRINT*,'something wrong with lambda grid!'
     STOP
  END IF
  CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) ! Find the rate
  CALL SYSTEM_CLOCK(COUNT=clock_start)
  CALL GETARG(1,dustyinpfile)
  ! Get or generate master input file
  IF (TRIM(dustyinpfile).eq."") THEN
     PRINT*, "No input file name found on command line." 
     PRINT*, "Proceeding with default file dusty.master"
     dustyinpfile = "dusty.master"
  ELSE
     suffix = dustyinpfile(LEN(TRIM(dustyinpfile))-6:)
     IF (suffix .eq. '.master') THEN
        PRINT*, "Found master input file ",TRIM(dustyinpfile), &
             " on on command line."
     ELSE IF (suffix .eq. '.inp') THEN
        PRINT*, "Found normal input file ",TRIM(dustyinpfile), &
             " on on command line."
        CALL GETARG(2,verbosity)
        if (TRIM(verbosity).eq."") verbosity = "2"
        !generate temperay master file
        OPEN(unit=100,file="temp.master")
        WRITE(100,*) "% Temporary master input for single input file"
        WRITE(100,*) "% DUSTY version: ",version
        WRITE(100,*) "verbose = ",verbosity
        WRITE(100,*) "% filename:"
        !WRITE(100,*) dustyinpfile(1:LEN(TRIM(dustyinpfile))-4)
        WRITE(100,*) dustyinpfile
        CLOSE(unit=100)
        dustyinpfile = "temp.master"
     ELSE
        PRINT*,'WARNING NOT A PROPER INPUT FILE!!!!'
        PRINT*,'DUSTY STOPPED'
        STOP
     END IF
  END IF
  OPEN(13,file=trim(dustyinpfile),status='old')
  iVerb = RDINP(.true.,13,6)
  READ(13,'(a)',iostat=io_status) apath
  DO WHILE (io_status.ge.0)
     call alloc_mem()
     call alloc_mem_nL()
     CALL clean(apath,path,lpath)
     ! remove .inp from input master file 
     path = path(1:LEN(TRIM(path))-4)
     lpath = lpath-4
     IF (empty(path).ne.1) THEN
        CALL attach(path,lpath,'.inp',nameIn)
        CALL attach(path,lpath,'.out',nameOut)
        print*,'working on file: ',TRIM(nameIn)
        CALL Input(nameIn,nameOut,tau1,tau2,GridType,Nmodel)
        IF (iVerb.gt.0) THEN
           print*,'working on input file: ',TRIM(nameIn)
           IF (iVerb.ge.2) print*,'Done with reading input'
        ENDIF
        IF (error.ne.3) THEN 
           IF (error.eq.0) THEN 
              IF(ALLOCATED(tau)) DEALLOCATE(tau)
              ALLOCATE(tau(Nmodel))
              tau = 0
              CALL GetTau(tau1,tau2,GridType,Nmodel,tau)
              IF (iVerb.ge.2) print*,'Done with GetTau'
              CALL Kernel(path,lpath,tau,Nmodel)
           END IF
        ELSE
           PRINT*,
        END IF
     END IF
     READ(13,'(a)',iostat=io_status) apath
     call dealloc_mem()
  END DO
  CLOSE(13)
  IF (suffix .eq. '.inp') THEN 
     OPEN(unit=100,file="temp.master")
     CLOSE(unit=100,status='delete')
  END IF
  CALL SYSTEM_CLOCK(COUNT=clock_end)
  print*,'ellapsed time:',(clock_end-clock_start)/clock_rate,'[s]'
END PROGRAM DUSTY

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
