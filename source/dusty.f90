PROGRAM DUSTY
  USE COMMON
  IMPLICIT NONE
  INTEGER :: clock_rate, clock_start, clock_end, io_status, lpath
  INTEGER :: empty, GridType
  INTEGER :: Nmodel
  DOUBLE PRECISION :: RDINP, tau1, tau2
  DOUBLE PRECISION, allocatable :: tau(:)
  CHARACTER(len=4)   :: suffix,verbosity
  CHARACTER(len=235) :: dustyinpfile, path, apath, nameIn, nameOut, stdf(7)
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
     PRINT*, "Proceeding with default file dusty.mas"
     dustyinpfile = "dusty.mas"
  ELSE
     suffix = dustyinpfile(LEN(TRIM(dustyinpfile))-3:)
     IF (suffix .eq. '.mas') THEN
        PRINT*, "Found master input file ",TRIM(dustyinpfile), &
             " on on command line."
     ELSE IF (suffix .eq. '.inp') THEN
        PRINT*, "Found normal input file ",TRIM(dustyinpfile), &
             " on on command line."
        CALL GETARG(2,verbosity)
        if (TRIM(verbosity).eq."") verbosity = "1"
        !generate temperay master file
        OPEN(unit=100,file="temp.mas")
        WRITE(100,*) "% Temporary master input for single input file"
        WRITE(100,*) "% DUSTY version: ",version
        WRITE(100,*) "verbose = ",verbosity
        WRITE(100,*) "% filename:"
        WRITE(100,*) dustyinpfile(1:LEN(TRIM(dustyinpfile))-4)
        CLOSE(unit=100)
        dustyinpfile = "temp.mas"
     END IF
  END IF
  OPEN(13,file=trim(dustyinpfile),status='old')
  iVerb = RDINP(.true.,13)
  READ(13,'(a)',iostat=io_status) apath
  DO WHILE (io_status.ge.0)
     CALL clean(apath,path,lpath)
     IF (empty(path).ne.1) THEN
        CALL attach(path,lpath,'.inp',nameIn)
        CALL attach(path,lpath,'.out',nameOut)
        CALL Input(nameIn,nameOut,tau1,tau2,GridType,Nmodel)
        IF (iVerb.gt.0) THEN
           print*,'working on input file: ',TRIM(nameIn)
           IF (iVerb.ge.2) print*,'Done with reading input'
        ENDIF
        IF (error.ne.3) THEN 
           IF (error.eq.0) THEN 
              IF(ALLOCATED(tau)) DEALLOCATE(tau)
              ALLOCATE(tau(Nmodel))
              CALL GetTau(tau1,tau2,GridType,Nmodel,tau)
              IF (iVerb.ge.2) print*,'Done with GetTau'
              CALL Kernel(path,lpath,tau,Nmodel)
!                 CALL Kernel_matrix(nG,path,lpath,tau,Nrec,Nmodel,GridType,error)
!              ELSE
!              END IF
           END IF
        ELSE
           PRINT*,
        END IF
     END IF
     READ(13,'(a)',iostat=io_status) apath
  END DO
  CLOSE(13)
!!$  IF (suffix .eq. '.inp') THEN 
!!$     OPEN(unit=100,file="temp.mas")
!!$     CLOSE(unit=100,status='delete')
!!$  END IF
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
  open(4, file='lambda_grid.dat', status = 'old')
  nL = RDINP(Equal,4)
  allocate(lambda(nL))
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
