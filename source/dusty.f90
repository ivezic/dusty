PROGRAM DUSTY
  USE COMMON
  IMPLICIT NONE
  INTEGER :: clock_rate, clock_start, clock_end, io_status, lpath
  INTEGER :: empty, nG, nameNK(10), GridType, Nrec, error
  INTEGER :: Nmodel,lambdaOK
  PARAMETER (NREC = 1000)
  DOUBLE PRECISION :: RDINP, tau1, tau2, tauIn(Nrec)
  DOUBLE PRECISION, allocatable :: tau(:)
  CHARACTER*4   :: suffix,verbosity
  CHARACTER*12  :: version
  CHARACTER*235 :: dustyinpfile, path, apath, nameIn, nameOut, stdf(7)
  CHARACTER*235 :: nameQ(npG)
  !-------------------------------------------------------
  !**************************
  !*** ABOUT THIS VERSION ***
  !**************************
  version= '4.00'
  CALL ChkLambda(lambdaOK)
  IF (lambdaOK.eq.0) THEN 
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
        WRITE(100,*) "% Temperary master input for single input file"
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
        CALL Input(nameIn,nG,nameOut,nameQ,nameNK,tau1,tau2,&
             tauIn,Nrec,GridType,Nmodel,error,version,stdf)
        IF (iVerb.gt.0) THEN
           print*,'working on input file: ',TRIM(nameIn)
           IF (iVerb.ge.2) print*,'Done with reading input'
        ENDIF
        IF (error.ne.3) THEN 
           call getOptPr(nG,nameQ,nameNK,error,stdf)
           IF (error.eq.0) THEN 
              IF (iVerb.ge.2) print*,'Done with getOptPr'
              IF(ALLOCATED(tau)) DEALLOCATE(tau)
              ALLOCATE(tau(Nmodel))
              CALL GetTau(nG,tau1,tau2,tauIn,Nrec,GridType,Nmodel,tau)
              IF (iVerb.ge.2) print*,'Done with GetTau'
!              IF (SPH) THEN 
                 CALL Kernel_matrix(nG,path,lpath,tauIn,tau,Nrec,Nmodel,GridType,error)
!              ELSE
!                 CALL Kernel(nG,path,lpath,tauIn,tau,Nrec,Nmodel,GridType,error)
!              END IF
           END IF
        ELSE
           PRINT*,
        END IF
     END IF
     READ(13,'(a)',iostat=io_status) apath
  END DO
  CLOSE(13)
  IF (suffix .eq. '.inp') THEN 
     OPEN(unit=100,file="temp.mas")
     CLOSE(unit=100,status='delete')
  END IF
END PROGRAM DUSTY

!***********************************************************************
subroutine ChkLambda(lambdaOK)
!=======================================================================
! This subroutine reads and checks that the wavelength grid satisfies
! certain conditions described in the Manual (all wavelengths are given
! in microns). If everything went fine it returns lambdaOK = 1, and
! fills the wavelength grid in array lambda passed through COMMON grids2.
!                                                  [ZI,Feb'96; MN,Apr'99]
!=======================================================================
  use common
  implicit none
  integer  iL, nLam, lambdaOK
  double precision RDINP
  character str*235
  logical Equal
!-----------------------------------------------------------------------
  Equal = .true.
  ! first open the file with lambda grid
  open(4, file='lambda_grid.dat', status = 'old')
  nLam = RDINP(Equal,4)
  if (nLam.ne.npL) then
     write(*,*)' *************** a big error !! ***************** '
     write(*,*)'  the number of wavelengths in lambda_grid.dat is  '
     write(*,*)'  not equal to the specified npL in userpar.inc    '
     write(*,*)'  make sure the numbers are the same, recompile    '
     write(*,*)'  and try again.                                   '
     write(*,*)' ************************************************* '
     goto 999
  end if
  ! initialize lambda array
  do iL = 1, npL
     read(4,*,end=99) lambda(iL)
  end do
99 close(4)
  call sort(lambda,npL)
  ! check the ends of the lambda grid :
  !  if(lambda(1).gt.0.01d0) then
  !     write(*,*)' *************** WARNING! ********************** '
  !     write(*,*)'  the shortest wavelength in lambda_grid.dat has '
  !     write(*,*)'  to be 0.01 microns. correct this and try again!'
  !     write(*,*)' *********************************************** '
  !     goto 999
  !  end if
  !  if(lambda(npL).lt.36000.0d0) then
  !     write(*,*)' *************** WARNING! ******************* '
  !     write(*,*)'  the longest wavelength in lambda_grid.dat   '
  !     write(*,*)'  has to be 3.6e4 um. correct this and try again!'
  !     write(*,*)' ******************************************** '
  !     goto 999
  !  end if
  ! check the resolution:
  do iL = 2, npL
     if (lambda(iL)/lambda(iL-1).gt.1.51d0) then
        write(*,*)' ***************** WARNING!  *******************'
        write(*,*)' the ratio of two consecutive wavelengths in the'
        write(*,*)' grid has to be no bigger than 1.5. you have    '
        write(*,'(2(a4,1p,e8.2))') '    ',lambda(iL)/lambda(iL-1),' at ', lambda(iL)
        write(*,*)' correct this and try again!                    '
        write(*,*)' ***********************************************'
        goto 999
     end if
  end do
  ! everything is fine
  lambdaOK = 1
  goto 111
  ! something is wrong
999 lambdaOK = 0
!-----------------------------------------------------------------------
111 return
end subroutine ChkLambda
!***********************************************************************
