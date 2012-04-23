! ***********************************************************************
SUBROUTINE solve_matrix(model,taumax,nY,nYprev,itereta,nP,nCav,nIns,initial,deviat,iterfbol,fbolOK)
! =======================================================================
! This subroutine solves the continuum radiative transfer problem for a
! spherically symmetric envelope.                      [Z.I., Nov. 1995]
! =======================================================================
  use common
  use interfaces
  IMPLICIT none
  !--- parameter 
  integer :: model,iterfbol,fbolOK, itereta, nY, nP, nYprev, nCav, nIns
  double precision :: deviat, taumax
  logical initial
  !--- local variables
  integer iPstar, EtaOK, iY
  double precision pstar,TAUlim, delta
  double precision, allocatable :: fs(:,:),u_old(:,:),us(:,:),T4_ext(:)
  double precision, allocatable :: emiss(:,:,:)
!!$
!!$  integer model, error, nG, iterfbol, fbolOK,grid,iY,iL,nY_old,y_incr,imu, &
!!$       iPstar,EtaOK , iP, iZ, nZ, iOut
!!$  double precision pstar,taulim, us(npL,npY),comp_fdiff_bol(npY), &
!!$       comp_fdiff_bol1(npY),calc_fdiff(npY), &
!!$       delta, em(npG,npL,npY), omega(npG+1,npL), iauxl(npL),iauxr(npL), &
!!$       fdsp(npL,npY),fdsm(npL,npY), fdep(npL,npY),fdem(npL,npY),fs(npL,npY), &
!!$       T4_ext(npY), accfbol, fbol_max, fbol_min, aux, &
!!$       Udbol(npY), Usbol(npY), fDebol(npY),fDsbol(npY), maxFerr, deviat
!!$  logical initial
!!$
  ! ----------------------------------------------------------------------
 !!$  IF(iInn.eq.1) THEN
 !!$     write(38,*)'============================================='
 !!$     write(38,'(a7,i5)') ' model= ',model
 !!$     write(38,*)'============================================='
 !!$  END IF
  allocate(fs(nL,npY))
  allocate(us(nL,npY))
  allocate(T4_ext(npY))
  allocate(emiss(nG,nL,npY))
  allocate(u_old(nL,npY))
  IF (iX.NE.0) THEN
     CALL LINE(0,2,18)
     write(18,'(a7,i3,a20)') ' model ',model,'  RUN-TIME MESSAGES '
     CALL LINE(0,1,18)
     write(18,*)'===========  STARTING SOLVE  ==========='
  END IF
  error = 0
  iterfbol = 0
  fbolOK = 0
  itereta = 0
  EtaOK = 0
  !call SetGrids(pstar,iPstar,taumax,nY,nYprev,nP,nCav,nIns,initial,iterfbol+1,itereta+1)
  IF(sph) THEN
     ! Solve for spherical envelope:
     ! temporarily the star is approximated by a point source
     pstar = 0.0
     iPstar = 1
     ! select optical depth for the grid calculation
     ! based on dynamical range
     TAUlim = 0.5D+00*dlog(1.0D+00/dynrange)
     ! if actual maximal tau is smaller use that value
     IF (TAUmax.LT.TAUlim) TAUlim = TAUmax
     ! counter over ETA (for radiatively driven winds only)
     iterETA = 0
     EtaOK = 0
     ! iterations over ETA
     DO WHILE (EtaOK.EQ.0)
        iterETA = iterETA + 1
        IF (iX.NE.0.AND.denstyp.eq.3) THEN !3(RDW)
           write(18,*)'----------------------------------------'
           write(18,*)' ',iterETA,'. iteration over ETA'
        END IF
        IF (iVerb.EQ.2.AND.denstyp.eq.3)  write(*,*) ' ',iterETA,'. iteration over ETA' !3(RDW)
        ! counter for iterations over bolometric flux conservation
        iterFbol = 0
        FbolOK = 0
        DO WHILE (FbolOK.EQ.0)
           iterFbol = iterFbol + 1
           IF (iX.NE.0) THEN
              write(18,*)'  ',iterFbol,'. iteration over Fbol'
           END IF
           IF (iVerb.EQ.2) write(*,*) iterFbol,'. iteration over Fbol'
           ! solve the radiative transfer problem
           Call RADTRANSF_matrix(pstar,iPstar,nY,nYprev,nP,nCav,nIns,TAUmax,FbolOK,initial,deviat,&
                iterFbol,iterEta,model,us,u_old,fs,T4_ext,emiss)
           IF (iVerb.EQ.2) write(*,*) 'Done with RadTransf'
           ! error.EQ.3 : file with stellar spectrum not available
           IF (error.EQ.3) goto 999
           ! error.EQ.5 : Singular matrix
           IF (error.EQ.5) goto 999
           ! error.EQ.6 : Eta exceeds limitations
           IF (error.EQ.6) goto 999
           ! error.EQ.2 : P grid was not produced
           IF (error.EQ.2.AND.iterFbol.EQ.1.AND.iterETA.EQ.1) THEN
              ! if this is the first calculation end this model
              goto 999
           ELSE
                ! if this is a higher iteration use previous solution
                IF (error.EQ.2) THEN
                   IF (iX.NE.0.AND.iterFbol.GT.1) THEN
                      write(18,*)' ======= IMPORTANT WARNING ======== '
                      write(18,*)' In trying to conserve Fbol reached'
                      write(18,*)' the limit for grid sizes. Flux is '
                      write(18,'(a,1p,e9.3)')'  conserved to within ', deviat
                      write(18,*)' Treat all results with caution!'
                   END IF
                   IF (iX.NE.0.AND.iterFbol.EQ.1) THEN
                      write(18,*)' ======== IMPORTANT  WARNING ======== '
                      write(18,*)' In trying to converge on ETA reached'
                      write(18,*)' the limit for grid sizes. Flux is '
                      write(18,'(a,1p,e9.3)')'  conserved to within ', deviat
                      write(18,*)' Treat all results with caution!'
                   END IF
                   error = 0
                   FbolOK = 2
                END IF
           END IF
           ! just in case...
           IF (error.EQ.1) THEN
              IF (iX.NE.0) THEN
                 write(18,*)' *********  FATAL ERROR  *********'
                 write(18,*)' * Something was seriously wrong *'
                 write(18,*)' * Contact Z. Ivezic, M. Elitzur *'
                 write(18,*)' *********************************'
              END IF
              goto 999
           END IF
           ! if Fbol not conserved try again with a finer grid
           IF (FbolOK.EQ.0.AND.iterFbol.LT.20.AND.iX.NE.0) THEN
              write(18,*)'  ******** MESSAGE from SOLVE ********'
              write(18,*)'  Full solution does not conserve Fbol'
              write(18,*)'       Y       TAU/TAUtot        fbol'
              DO iY =1, nY
                 write(18,'(1p,3e13.4)')Y(iY), &
                      ETAzp(1,iY)/ETAzp(1,nY),fbol(iY)
              END DO
              write(18,*)'  Trying again with finer grids'
           END IF
             ! if could not conserve Fbol in 10 trials give it up
             IF (FbolOK.EQ.0.AND.iterFbol.GE.20) THEN
                IF (denstyp.eq.3) THEN !3(RDW)
                   IF (iX.NE.0) THEN
                      write(18,*)' **********  WARNING from SOLVE  **********'
                      write(18,*)' Could not obtain required accuracy in Fbol'
                      write(18,'(a26,1p,e10.3)')'  The achieved accuracy is:',deviat
                      write(18,*)' Will try to converge on the dynamics, but '
                      write(18,*)' treat all results with caution !!         '
                      write(18,*)' If accuracy<=0.01, or TAUmax>1000, this   '
                      write(18,*)' code probably cannot do it. Otherwise,    '
                      write(18,*)' please contact Z. Ivezic or M. Elitzur    '
                      write(18,*)' ******************************************'
                   END IF
                   FbolOK = 1
                ELSE
                   IF (iX.NE.0) THEN
                      write(18,*)' **********  WARNING from SOLVE  **********'
                      write(18,*)' Could not obtain required accuracy in Fbol'
                      write(18,'(a26,1p,e10.3)')'  The achieved accuracy is:',deviat
                      write(18,*)' !!!!  Treat all results with caution  !!!!'
                      write(18,*)' If accuracy<=0.01, or TAUmax>1000, this   '
                      write(18,*)' code probably cannot do it. Otherwise,    '
                      write(18,*)' please contact Z. Ivezic or M. Elitzur    '
                      write(18,*)' ******************************************'
                   END IF
                   FbolOK = 2
                END IF
             END IF
           ! end of loop over flux conservation
        END DO
          ! for winds check if ETA has converged...
          IF ((denstyp.eq.3).AND.FbolOK.NE.2) THEN !3(RDW)
!!$             ! ptr(2) is specified in INPUT and controls converg. crit.
!!$             IF (ptr(2).LT.1.0D-6.AND.iterETA.GT.2)THEN
!!$                EtaOK = 1
!!$             ELSE
!!$                CALL WINDS(nG,EtaOK,ftot)
!!$             END IF
!!$             IF (iterETA.GT.10.AND.EtaOK.EQ.0) THEN
!!$                EtaOK = 2
!!$                iWARNING = iWARNING + 1
!!$                IF (iX.NE.0) THEN
!!$                   write(18,*)' *********  WARNING  *********'
!!$                   write(18,*)' Could not converge on ETA in '
!!$                   write(18,*)' 10 iterations.'
!!$                   write(18,*)' *********************************'
!!$                END IF
!!$             END IF
!!$             ! ...or otherwise finish right away
          ELSE
             EtaOK = 1
          END IF
        ! end of loop over ETA
     END DO
  ELSE
       ! solve for slab case
!!$       CALL SLBsolve(model,nG,error)
!!$       ! error=4 means npY not large enough for oblique illumination grid
!!$       IF (error.eq.4) THEN
!!$          CALL MSG(15)
!!$          PRINT*,'stop MSG(15)'
!!$          goto 999
!!$       END IF
!!$       PRINT*,'Slab case is not implemented for the matrix method!'
!!$       STOP
  END IF
  ! analyze the solution and calculate some auxiliary quantities
  CALL analysis(nY,model,us,T4_ext,delta,deviat)
  IF (iVerb.EQ.2) write(*,*) 'Done with Analysis'
  IF (iX.NE.0) THEN
     write(18,*)' ==== SOLVE successfully completed ====='
     write(18,*)' ======================================='
  END IF
  print*,' ==== SOLVE successfully completed ====='
  ! -----------------------------------------------------------------------
999 deallocate(fs)
  deallocate(us)
  deallocate(T4_ext)
  deallocate(emiss)
  deallocate(u_old)
  RETURN
END SUBROUTINE solve_matrix
!***********************************************************************

!***********************************************************************
SUBROUTINE RADTRANSF_matrix(pstar,iPstar,nY,nYprev,nP,nCav,nIns,TAUmax,&
     FbolOK,initial,deviat,iterFbol,iterEta,model,us,u_old,fs,T4_ext,emiss)
!=======================================================================
! This subroutine solves the continuum radiative transfer problem for a
! spherically symmetric envelope.                      [Z.I., Nov. 1995]
!=======================================================================
  use common
  use interfaces
  IMPLICIT none
  !---parameter
  integer iPstar, FbolOK,iterFbol,iterEta,model,nY,nYprev,nP,nCav,nIns
  double precision :: pstar,TAUmax,deviat
  double precision,allocatable :: u_old(:,:),us(:,:), fs(:,:),T4_ext(:),emiss(:,:,:)
  logical initial
  !---local
  integer iaux,nPok,nYok, itlim, iter, conv, itnum, Fconv,  aconv, &
       iY,iL,Uconv,Uconv1,Uconv2,Uconv3
  double precision,allocatable :: mat0(:,:,:), mat1(:,:,:),&
       mifront(:,:,:), miback(:,:,:), UbolChck(:), Uchck(:,:), fbolold(:)
  double precision ::  BolConv, dmaxF, dmaxU, maxFerr,wgth
  double precision, allocatable :: omat(:,:),emiss_total(:,:),fdebol(:), &
       fdsbol(:)
 
!!$  integer iPstar, nG, FbolOK, error, iterFbol, model, &
!!$       BolConv, Conv, iaux, Fconv, iter,itnum, iY, itlim, uconv
!!$  double precision pstar,taulim,deviat, &
!!$       us(npL,npY),comp_fdiff_bol(npY), &
!!$       comp_fdiff_bol1(npY),calc_fdiff(npY), &
!!$       delta, em(npG,npL,npY), omega(npG+1,npL), iauxl(npL),iauxr(npL), &
!!$       fdsp(npL,npY),fdsm(npL,npY), fdep(npL,npY),fdem(npL,npY),fs(npL,npY), &
!!$       T4_ext(npY), accfbol, fbol_max, fbol_min, aux, &
!!$       Udbol(npY), Usbol(npY), fDebol(npY),fDsbol(npY), maxFerr,  &
!!$       dmaxF, dmaxU, fbolold(npY), mat0(npL,npY,npY), mat1(npL,npY,npY), &
!!$       miback(npL,npP,npY), mifront(npL,npP,npY), Tei, UbolChck(npY), &
!!$       Uchck(npL,npY), Uold(npL,npY), fdbol(npY)
!!$

!!$  allocate(T_old(nG,nY))
  allocate(omat(nL,npY))
  allocate(emiss_total(nL,npY))
  allocate(UbolChck(npY))
  allocate(Uchck(nL,npY))
  allocate(fbolold(npY))
  allocate(fdebol(npY))
  allocate(fdsbol(npY))
  !------------------------------------------------------------------------
  ! generate, or improve, or do not touch the Y and P grids
  IF (iterETA.EQ.1.OR.iterFbol.GT.1) THEN
     IF (iterFbol.EQ.1) THEN
        ! first time generate grids
        CALL SetGrids(pstar,iPstar,TAUmax,nY,nYprev,nP,nCav,nIns,initial,iterfbol,itereta)
        IF (error.NE.0) then 
           print*,' stopping ... something wrong!!!'
           STOP
        END IF
        IF (iVerb.EQ.2) write(*,*) 'Done with SetGrids'
     ELSE
        ! or improve the grid from previous iteration
        !CALL ChkFlux(nY,nYprev,fBol,accFlux,iaux,iterEta)
        call Flux_Consv(nY,nYprev,Ncav,itereta,iterfbol,fbol,fdebol,fdsbol,fbolOK,maxFerr)
        ! added in ver.2.06
        ! IF (maxFerr.GT.0.5D+00) CALL DblYgrid(error) --**--
        IF (error.NE.0) goto 999
        ! generate new impact parameter grid
        ! increase the number of rays through the cavity
        IF (Ncav.LT.80) THEN
           Ncav = 2 * Ncav
           IF (iX.NE.0) write(18,'(a20,i3)')' Ncav increased to:',Ncav
        END IF
        ! increase the number of rays per y-grid interval
        IF (iterFbol.EQ.3.AND.Nins.EQ.2) THEN
           Nins = Nins + 1
           IF (iX.NE.0) write(18,'(a20,i3)')' Nins increased to:',Nins
        END IF
        CALL Pgrid(pstar,iPstar,nY,nP,nCav,nIns)
        ! if P grid is not OK end this model
        IF (error.NE.0) goto 999
        IF (iX.NE.0) THEN
           write(18,'(a23,i3)')' Y grid improved, nY =',nY
           write(18,'(a23,i3)')'                  nP =',nP
           write(18,'(a23,i3)')'                Nins =',Nins
           write(18,'(a23,i3)')'                Ncav =',Ncav
        END IF
     END IF
  ELSE
     IF (iX.NE.0) write(18,*)' Using same Y and P grids'
  END IF
  if (allocated(mat0)) deallocate(mat0)
  if (allocated(mat1)) deallocate(mat1)
  if (allocated(mifront)) deallocate(mifront)
  if (allocated(miback)) deallocate(miback)
  allocate(mat0(nL,nY,nY))
  allocate(mat1(nL,nY,nY))
  allocate(mifront(nL,npP,nY))
  allocate(miback(nL,npP,nY))
  ! generate spline coefficients for ETA
  CALL setupETA(nY,nYprev,itereta)
  ! evaluate ETAzp
  CALL getETAzp(nY,nP)
  ! generate albedo through the envelope
  CALL getOmega(nY)
  do iY=1,nY
     do iL=1,nL
          omat(iL,iY) = SigmaS(nG,iL) / (SigmaA(nG,iL) + SigmaS(nG,iL))        
     end do
  end do
  ! generate stellar moments
  ! CALL Star(pstar,ETAzp,error) --**--
  call Find_Tran(pstar,nY,nP,T4_ext,us,fs)
  IF (iVerb.EQ.2) write(*,*) 'Done with Find_Tran'
  ! issue a message in fname.out about the condition for neglecting
  ! occultation only if T1 is given in input:
  IF(typEntry(1).eq.1.AND.model.eq.1) THEN
     IF(iterFbol.eq.1.AND.iterEta.eq.1.AND.Right.eq.0) CALL OccltMSG(us)
  END IF
  ! generate the first approximation for Td
  ! CALL InitTemp(ETAzp,nG)  --**--
  if (initial.and.iterfbol.eq.1) then
     call Init_Temp(nY,T4_ext,us)
     if(iVerb.eq.2) write(*,*)' Done with initial dust temperature.'
  end if
  !IF (iVerb.EQ.2) write(*,*) 'Done with InitTemp'
  ! find radiative transfer matrices
  IF (iX.NE.0) write(18,*)' Calculating weight matrices'
  IF (iVerb.EQ.2) write(*,*) 'Calculating weight matrices'
  IF (iD.GE.1) THEN
     IF (iVerb.GE.1) write(*,*) 'No disk option in this version'
     ! if disk included:
     ! CALL MatrixD(ETAzp,pstar,iPstar,mat0,mat1,matD,mifront,miback)
  ELSE
     CALL Matrix(pstar,iPstar,mat0,mat1,mifront,miback,nP,nY,nPok,nYok,T4_ext)
  END IF
  Conv = 0
  iter = 0
  ! itlim is an upper limit on number iterations
  itlim = 10000
  IF (iX.NE.0) write(18,*)' Weight matrices OK, calculating Tdust'
  IF (iVerb.EQ.2) write(*,*)' Weight matrices OK, calculating Tdust'
 !!$  IF (iInn.eq.1) THEN
 !!$     write(38,'(a8,i5)') '    nY= ',nY
 !!$     write(38,*) '    iter   maxFerr     dmaxU       dmaxF        T1         Fe1'
 !!$  END IF
  if (initial.and.iterfbol.ne.1.and.typentry(1).eq.5) then
     call find_Text(nY,T4_ext)
  elseif (.not.initial.and.typentry(1).eq.5) then
     call find_Text(nY,T4_ext)
  end if
  do iY = 1, nY
     Jext(iY) = sigma/pi * T4_ext(iY)
!!$     do iL=1,nL
!!$        u_old(iL,iY) = shpR(iL)
!!$        utot(iL,iY) = shpR(iL)
!!$     end do
  end do
  CALL Bolom(Utot,Ubol,nY)
  ! === Iterations over dust temperature =========
  DO WHILE (Conv.EQ.0.AND.iter.LE.itlim)
     iter = iter + 1
     if (initial.and.iterfbol.ne.1.and.typentry(1).eq.5) then
        call find_Text(nY,T4_ext)
     elseif (.not.initial.and.typentry(1).eq.5) then
        call find_Text(nY,T4_ext)
     end if
     ! find emission term
     call Emission(nY,T4_ext,emiss,emiss_total)
     ! solve for Utot
     CALL Invert(nY,mat0,Us,Emiss,U_old,omat)
     !CALL LAMBDA_ITER(nY,mat0,Us,Emiss,U_old,omat)
     call Find_Temp(nY,T4_ext)
     IF(error.NE.0) goto 999
     u_old = utot
     ! find new Td
     ! CALL FindTemp(1,Utot,nG,Td) --**--
     !CALL Bolom(Utot,Ubol,nY)
     ! --------------------------------------
     ! every itnum-th iteration check convergence:
     if (iter.gt.500) then
        itnum = 30
     else iF (iter.GT.80) THEN
        itnum = 10
     ELSE
        itnum = 6
     END IF
     ! first find 'old' flux (i.e. in the previous iteration)
     IF (MOD(iter+1,itnum).EQ.0) THEN
        call Emission(nY,T4_ext,emiss,emiss_total)
        CALL Multiply(1,nY,nY,nL,nL,mat1,Utot,omat,0,fs,fds)
        CALL Multiply(0,nY,nY,nL,nL,mat1,emiss_total,omat,0,fs,fde)
        do iY=1,nY
           fbolold(iY) = 0.
           do iL=1,nL
              fds(iL,iY) = fds(iL,iY)*4*pi
              fde(iL,iY) = fde(iL,iY)*4*pi
              ftot(iL,iY) = fs(iL,iY)+fds(iL,iY) + fde(iL,iY)
              ! begin simpson intgral
              ! weigths
              if (iL.ne.1.and.iL.ne.nL) then
                 wgth = 0.5d0*(lambda(iL+1)-lambda(iL-1))
              else
                 if (iL.eq.1) wgth = 0.5d0*(lambda(1+1)-lambda(1))
                 if (iL.eq.nL) wgth = 0.5d0*(lambda(nL)-lambda(nL-1))
              end if
              ! add contribution to the integral
              fbolold(iY) = fbolold(iY) + ftot(iL,iY)/lambda(iL)*wgth
              ! end simpson integral
           end do
        end do
        !CALL Add(nY,nY,nL,nL,fs,fds,fde,ftot)
        ! find bolometric flux
        !CALL Bolom(ftot,fbolold,nY)
     END IF
     IF (MOD(iter,itnum).EQ.0) THEN
        ! first calculate total flux
        call Emission(nY,T4_ext,emiss,emiss_total)
        CALL Multiply(1,nY,nY,nL,nL,mat1,Utot,omat,0,fs,fds)
        CALL Multiply(0,nY,nY,nL,nL,mat1,emiss_total,omat,0,fs,fde)
        do iY=1,nY
           fbol(iY) = 0.
           do iL=1,nL
              fds(iL,iY) = fds(iL,iY)*4*pi
              fde(iL,iY) = fde(iL,iY)*4*pi
              ftot(iL,iY) = fs(iL,iY)+fds(iL,iY) + fde(iL,iY)
              ! begin simpson intgral
              ! weigths
              if (iL.ne.1.and.iL.ne.nL) then
                 wgth = 0.5d0*(lambda(iL+1)-lambda(iL-1))
              else
                 if (iL.eq.1) wgth = 0.5d0*(lambda(1+1)-lambda(1))
                 if (iL.eq.nL) wgth = 0.5d0*(lambda(nL)-lambda(nL-1))
              end if
              ! add contribution to the integral
              fbol(iY) = fbol(iY) + ftot(iL,iY)/lambda(iL)*wgth
              ! end simpson integral
           end do
        end do
        !CALL Add(nY,nY,nL,nL,fs,fds,fde,ftot)
        ! find bolometric flux
        !CALL Bolom(ftot,fbol,nY)
        ! check convergence of bolometric flux
        ! call Flux_Consv(nY,nYprev,Ncav,itereta,fbol,fbol,fbolOK,maxFerr)
        CALL Converg1(nY,fbolold,fbol,Fconv,dmaxF)
        ! check convergence of energy density
        CALL Converg2(nY,U_old,Utot,Uconv,dmaxU)
        ! find maximal fbol error
        CALL FindErr(nY,fbol,maxFerr)
        !------  printout of errors and convergence with iter.(inner flag): -------
!!$        IF(iInn.EQ.1) THEN
!!$           write(38,'(i7,1p,5e12.4)') iter,maxFerr,dmaxU,dmaxF,Td(1,1),sigma*Tei**4.0D+00
!!$        END IF
        !--------------------------------------------------------------
        IF (abs(maxFerr).LE.accFlux) THEN
           BolConv = 1
        ELSE
           BolConv = 0
        END IF
        !print*,Fconv,Uconv,BolConv,dmaxF,accFlux
        ! total criterion for convergence: Utot must converge, and ftot
        ! must either converge or have the required accuracy
        IF (Uconv*(Fconv+BolConv).GT.0) Conv = 1
     END IF
     ! --------------------------------------
  END DO
  !    === The End of Iterations over Td ===
  IF (iX.NE.0) THEN
     IF (iter.LT.itlim) write(18,*) ' Convergence achieved, number of'
     write(18,'(a34,i4)') ' iterations over energy density: ',iter
     write(18,'(a30,1p,e8.1)') ' Flux conservation OK within:',maxFerr
     IF (iter.GE.itlim) THEN
        CALL MSG(2)
        print*," Reached max number of iterations:",itlim
     END IF
  END IF
  ! calculate the emission term for the converged Td
  call Emission(nY,T4_ext,emiss,emiss_total)
!!$  CALL Emission_matrix(1,0,nG,Us,Em)
  ! calculate flux
  CALL Multiply(1,nY,nY,nL,nL,mat1,Utot,omat,0,fs,fds)
  CALL Multiply(0,nY,nY,nL,nL,mat1,emiss_total,omat,0,fs,fde)
  do iY=1,nY
     fbol(iY) = 0.
     do iL=1,nL
        fds(iL,iY) = fds(iL,iY)*4*pi
        fde(iL,iY) = fde(iL,iY)*4*pi
        ftot(iL,iY) = fs(iL,iY)+fds(iL,iY) + fde(iL,iY)
        ! begin simpson intgral
        ! weigths
        if (iL.ne.1.and.iL.ne.nL) then
           wgth = 0.5d0*(lambda(iL+1)-lambda(iL-1))
        else
           if (iL.eq.1) wgth = 0.5d0*(lambda(1+1)-lambda(1))
           if (iL.eq.nL) wgth = 0.5d0*(lambda(nL)-lambda(nL-1))
        end if
        ! add contribution to the integral
        fbol(iY) = fbol(iY) + ftot(iL,iY)/lambda(iL)*wgth
        ! end simpson integral
     end do
  end do
  call BOLOM(fs,fsbol,nY)
  call BOLOM(fde,fDebol,nY)
  call BOLOM(fds,fDsbol,nY)
  DO iY = 1, nY
     write(7777,*) Y(iY),fDebol(iY),fDsbol(iY),fsbol(iY),Td(1,iY)
  END DO
  close(7777)
  !CALL Add(nY,nY,nL,nL,fs,fds,fde,ftot)
  ! CALL Bolom(ftot,fbol,nY)
  ! check whether, and how well, is bolometric flux conserved
  ! CALL ChkBolom(nY,fbol,accFlux,deviat,FbolOK)
  CALL FindErr(nY,fbol,maxFerr)
  ! added in ver.2.06
  IF (maxFerr.LT.accFlux) FbolOK = 1
  IF (iVerb.EQ.2) write(*,*)' Achieved Flux accuracy:',maxFerr
  !***********************************
  ! calculate additional output quantities
  ! 1) energy densities
  CALL Multiply(1,nY,nY,nL,nL,mat0,Utot,omat,0,Us,Uds)
  !CALL Multiply(0,nY,nY,nL,nL,mat0,emiss_total,omat,0,fs,fde)
  !do iL=1,nL
  !   do iY=1,nY
  !      fde(iL,iY) = fde(iL,iY)*4*pi
  !   end do
  !end do
  CALL Add(nY,nY,nL,nL,Us,Uds,Ude,Uchck)
  CALL Bolom(Utot,Ubol,nY)
  CALL Bolom(Uchck,UbolChck,nY)
  ! 2) scaled radial optical depth, tr
!!$  DO iY = 1, nY
!!$     tr(iY) = ETAzp(1,iY) / ETAzp(1,nY)
!!$  END DO
  ! 3) calculate intensity (at the outer edge) if required
!!$  IF(iC.NE.0) THEN
!!$     IF (iX.NE.0) write(18,*) 'Calculating intensities'
!!$     !CALL FindInt(nG,ETAzp)  --**--
!!$     call sph_int(nG,omega,fs)
!!$     IF (iVerb.EQ.2) write(*,*) 'Done with SPH_INT(FindInt)'
!!$  END IF
!!$  ! if needed convolve intensity with the PSF
!!$  IF (iPSF.NE.0) THEN
!!$     CALL Convolve(IntOut)
!!$     IF (iVerb.EQ.2) write(*,*) 'Done with Convolve'
!!$  END IF
!!$  ! if needed find the visibility function
!!$  IF (iV.NE.0) THEN
!!$     CALL Visibili(IntOut)
!!$     IF (iVerb.EQ.2) write(*,*) 'Done with Visibili'
!!$  END IF
!!$  !============ if the inner flag iInn=1:  =========
!!$  IF(iX.GE.1 .AND. iInn.EQ.1) THEN
!!$     ! if additional output needed in message files when iInn = 1
!!$     CALL Bolom(fs,fsbol,nY)
!!$     CALL Bolom(fs,fsbol,nY)
!!$     CALL ADD2(fds,fde,fdbol,nY)
!!$     write(18,'(a11,1p,E11.3)')'   TAUfid =',TAUfid
!!$     write(18,'(a11,1p,E11.3)')'  MaxFerr =',maxFerr
!!$     write(18,*) '     tr      fbol       fsbol      fdbol       Ubol  '
!!$     !    &'     tr      fbol       fsbol      fdbol     Usbol      Udbol'
!!$     CALL Bolom(Us,Usbol,nY)
!!$     DO iY = 1, nY
!!$        write(18,'(1p,6E11.3)') tr(iY), fbol(iY), fsbol(iY), fdbol(iY), Ubol(iY)
!!$     END DO
!!$  END IF
!!$  !=====================
!!$  !    if needed solve for disk quantities
!!$  !     corr = 0.0
!!$  IF (iD.GE.1) THEN
!!$     IF (iVerb.GE.1) write(*,*) 'No disk option in this version'
!!$  END IF
!!$  !---------------------------------------------------------------------
999 deallocate(mat0)
  deallocate(mat1)
  deallocate(mifront)
  deallocate(miback)
  deallocate(UbolChck)
  deallocate(Uchck)
  deallocate(fbolold)
  deallocate(omat)
  deallocate(emiss_total)
  deallocate(fdebol)
  deallocate(fdsbol)
  RETURN
END SUBROUTINE RADTRANSF_matrix
!***********************************************************************

!***********************************************************************
SUBROUTINE INVERT(nY,mat,Us,Em,Uold,omat)
!=======================================================================
!This subroutine solves the linear system
![Utot] = [Us+Em] + [mat0]*[omega*Utot] by calling LINSYS subroutine.
!       [Z.I., Nov. 1995]
!=======================================================================
  use common
  use interfaces
  IMPLICIT none
  !--- parameter
  integer nY
  double precision,allocatable :: Us(:,:), Uold(:,:), Em(:,:,:),mat(:,:,:),&
       omat(:,:)
  !--- local
  DOUBLE PRECISION  delTAUsc, facc, EtaRat, accFbol
  INTEGER iG,iL,iY, iYaux, Kronecker
  DOUBLE PRECISION,allocatable ::  A(:,:), B(:), X(:)
  !--------------------------------------------------------------------
  allocate(B(nY))
  allocate(A(nY,nY))
  allocate(X(nY))
  error = 0
  ! first copy Utot to Uold
  !  DO iL = 1, nL
  !     DO iY = 1, nY
  !        Uold(iL,iY) = Utot(iL,iY)
  !     END DO
  !  END DO
  ! calculate new energy density
  ! loop over wavelengths
  !$OMP PARALLEL DO PRIVATE(iL,Kronecker,iY,iYaux,iG,B,A,X) 
  DO iL = 1, nL
     ! generate the vector of free coefficients, B, and matrix A
     DO iY = 1, nY
        Uold(iL,iY) = Utot(iL,iY)
        B(iY) = Us(iL,iY)
        DO iYaux = 1, nY
           Kronecker = 0
           IF (iY.EQ.iYaux) Kronecker = 1
           ! loop over grains
           DO iG = 1, nG
              B(iY) = B(iY) + (1.-omat(iL,iYaux))*Em(iG,iL,iYaux)*mat(iL,iY,iYaux)  
              A(iY,iYaux) = Kronecker - omat(iL,iYaux) * mat(iL,iY,iYaux)        
           END DO
        END DO
     END DO
     ! solve the system
     CALL LINSYS(nY,A,B,X)
     IF(error.NE.0) THEN
        CALL MSG(20)
        print*, 'MSG(20)'
        stop
        !RETURN
     END IF
     ! store the result
     DO iY = 1, nY
        IF (X(iY).GE.dynrange) THEN
           Utot(iL,iY) = X(iY)
        ELSE
           Utot(iL,iY) = 0.0D+00
        END IF
     END DO
  END DO
  !$OMP END PARALLEL DO
  !-----------------------------------------------------------------------
  deallocate(B)
  deallocate(A)
  deallocate(X)
  RETURN
END SUBROUTINE INVERT
!***********************************************************************
!!$
!***********************************************************************
SUBROUTINE LAMBDA_ITER(nY,mat,Us,Em,Uold,omat)
!=======================================================================
!This subroutine solves the linear system
![Utot] = [Us+Em] + [mat0]*[omega*Utot] by calling LINSYS subroutine.
!       [Z.I., Nov. 1995]
!=======================================================================
  use common
  use interfaces
  IMPLICIT none
  !--- parameter
  integer nY
  double precision,allocatable :: Us(:,:), Uold(:,:), Em(:,:,:),mat(:,:,:),&
       omat(:,:)
  !--- local
  DOUBLE PRECISION  delTAUsc, facc, EtaRat, accFbol
  INTEGER iG,iL,iY, iYaux, Kronecker
  !--------------------------------------------------------------------
  error = 0
  ! first copy Utot to Uold
  DO iL = 1, nL
     DO iY = 1, nY
        Uold(iL,iY) = Utot(iL,iY)
     END DO
  END DO
  ! calculate new energy density
  ! loop over wavelengths
  !$OMP PARALLEL DO PRIVATE(iL,Kronecker,iY,iYaux,iG) 
  DO iL = 1, nL
     DO iY = 1, nY
        Uold(iL,iY) = Utot(iL,iY)
        Utot(iL,iY) = Us(iL,iY)
        DO iYaux = 1, nY
           Utot(iL,iY) = Utot(iL,iY)+omat(iL,iYaux)*Uold(iL,iYaux)*mat(iL,iY,iYaux)
           do iG=1,nG
              Utot(iL,iY) = Utot(iL,iY)+(1.-omat(iL,iYaux))*Em(iG,iL,iYaux)*mat(iL,iY,iYaux)
           end do
        END DO
        IF (Utot(iL,iY).lt.dynrange) THEN
           Utot(iL,iY) = 0.0D+00
        END IF
     END DO
  END DO
  !$OMP END PARALLEL DO
  !-----------------------------------------------------------------------
  RETURN
END SUBROUTINE LAMBDA_ITER     
!***********************************************************************
!!$
!!$
!***********************************************************************
SUBROUTINE MULTIPLY(type,np1,nr1,np2,nr2,mat,vec1,omat,flag,q1,q2)
!=======================================================================
!This subroutine evaluates the following expression:
![q2] = flag*[q1] + [mat]*[tt*vec1]. Here tt is [omat] for type=1 and
!1-[omat] for type=2. mat is matrix of physical size (np2,np1,np1) and
!real size (nr2,nr1,nr1). omat, vec1, q1 and q2 are matrices of
!physical size (np2,np1) and real size (nr2,nr1).     [Z.I., Nov. 1995]
!=======================================================================
  use common
  use interfaces
  IMPLICIT none
  !---parameter
  integer :: type,np1,nr1,np2,nr2,flag
  double precision, allocatable :: mat(:,:,:), vec1(:,:), omat(:,:), &
       q1(:,:), q2(:,:)
  !---local
  DOUBLE PRECISION delTAUsc, facc, EtaRat, accFbol
  INTEGER i2, i1, idum
  DOUBLE PRECISION aux
  !-----------------------------------------------------------------------
  ! loop over index 2
  !$OMP PARALLEL DO PRIVATE(i1,idum,aux)
  DO i2 = 1, nr2
     ! loop over index 1
     DO i1 = 1, nr1
        q2(i2,i1) = flag * q1(i2,i1)
        ! loop over dummy index (multiplication)
        DO idum = 1, nr1
           IF (type.EQ.1) THEN
              !aux = omat(nG+1,i2)
              aux = omat(i2,idum) 
           ELSE
              !aux = 1.0D+00 - omat(nG+1,i2) 
              aux = 1.0D+00 - omat(i2,idum) 
           END IF
           q2(i2,i1) = q2(i2,i1) + mat(i2,i1,idum)*aux*vec1(i2,idum)
        END DO
        IF (q2(i2,i1).LT.dynrange*dynrange) q2(i2,i1) = 0.0D+00
     END DO
  END DO
  !$OMP END PARALLEL DO
  !---------------------------------------------------------------------
  RETURN
END SUBROUTINE MULTIPLY
!***********************************************************************


!***********************************************************************
SUBROUTINE Converg1(nY,Aold,Anew,Aconv,dmax)
!=======================================================================
!This subroutine checks convergence of an array A(nL,nY) between values
!given in Aold and Anew, when the values are larger than dynrange. If
!the maximum relative difference is smaller than the required accuracy,
!Aconv is assigned 1, otherwise 0.              [Z.I.Jul 96;M.N.Apr.97]
!=======================================================================
  use common
  use interfaces
  IMPLICIT none
  !---parameter
  integer nY,Aconv
  double precision :: dmax
  double precision,allocatable ::  Aold(:), Anew(:)
  !---local
  INTEGER iY
  DOUBLE PRECISION delta
  !-----------------------------------------------------------------------
  Aconv = 1
  dmax = 0.0D+00
  ! loop over radial positions
  DO iY = 1, nY
     ! do it only for elements larger than dynrange
     IF (Anew(iY).GE.dynrange) THEN
        ! find relative difference
        delta = dabs((Anew(iY)-Aold(iY))/Anew(iY))
        IF (delta.GT.dmax) dmax = delta
     END IF
  END DO
  IF (dmax.GT.accFlux/nL*1e-2) Aconv = 0
  !---------------------------------------------------------------------
  RETURN
END SUBROUTINE Converg1
!***********************************************************************


!***********************************************************************
SUBROUTINE Converg2(nY,Aold,Anew,Aconv,dmax)
!=======================================================================
!This subroutine checks convergence of an array A(nL,nY) between values
!given in Aold and Anew, when the values are larger than dynrange. If
!the maximum relative difference is smaller than required accuracy,
!Aconv is assigned 1, otherwise 0.             [Z.I.Jul 96; M.N.Apr.97]
!=======================================================================
  use common
  use interfaces
  IMPLICIT none
  !---parameter
  integer nY,Aconv
  double precision :: dmax
  double precision,allocatable :: Aold(:,:), Anew(:,:)
  !---local
  INTEGER iY, iL
  DOUBLE PRECISION delta
!-----------------------------------------------------------------------
  Aconv = 1
  dmax = 0.0D+00
  ! loop over wavelengths
  DO iL = 1, nL
     ! loop over radial positions
     DO iY = 1, nY
        ! do it only for elements larger than dynrange
        IF (Anew(iL,iY).GE.dynrange) THEN
           ! find relative difference
           delta = dabs((Anew(iL,iY)-Aold(iL,iY))/Anew(iL,iY))
           IF (delta.GT.dmax) dmax = delta
        END IF
     END DO
  END DO
  IF (dmax.GT.accFlux) Aconv = 0
  !-----------------------------------------------------------------------
  RETURN
END SUBROUTINE Converg2
!***********************************************************************


!!$!***********************************************************************
!!$SUBROUTINE ChkBolom(nY,qbol,accur,dev,FbolOK)
!!$!=======================================================================
!!$!This subroutine checks if any element of qbol(i), i=1,nY differs for
!!$!more than accuracy from the median value fmed. If so FbolOK = 0,
!!$!otherwise FbolOK = 1. dev is maximal deviation from fmed. [ZI,'96;MN'00]
!!$!=======================================================================
!!$  use common
!!$  use interfaces
!!$  IMPLICIT none
!!$  !---parameter
!!$  integer nY,FbolOK
!!$  double precision :: accur,dev
!!$  DOUBLE PRECISION,allocatable ::  qBol(:) 
!!$  !---local
!!$  INTEGER iY
!!$  DOUBLE PRECISION fmax,AveDev,RMS
!!$  !-----------------------------------------------------------------------
!!$  FbolOK = 1
!!$  dev = 0.0D+00
!!$  ! loop over iY (radial coordinate)
!!$  if (slb) then
!!$     !CALL SLBmisc(qBol,fmax,fmed,AveDev,RMS,nY)
!!$     PRINT*, 'Matrix method not for slab case'
!!$  END IF
!!$  print*,'---!CALL SLBmisc(qBol,fmax,fmed,AveDev,RMS,nY)---'
!!$  DO iY = 1, nY
!!$     IF (abs(fmed-qBol(iY)).GT.accur) FbolOK = 0
!!$     IF (abs(fmed-qBol(iY)).GT.dev) dev = abs(fmed-qBol(iY))
!!$  END DO
!!$  !-----------------------------------------------------------------------
!!$  RETURN
!!$END SUBROUTINE ChkBolom
!!$!***********************************************************************

!***********************************************************************
SUBROUTINE matrix(pstar,iPstar,m0,m1,mifront,miback,nP,nY,nPok,nYok,T4_ext)
!=======================================================================
!This subroutine evaluates radiative transfer matrix for spherically
!symmetric envelope. Here m is the order of moment (0 for energy dens.,
!1 for flux, 2 for pressure etc.), ETAzp is array of optical depths
!along the line of sight and mat is radiative transfer matrix.
!=======================================================================
  use common
  use interfaces
  IMPLICIT none
  ! --- parameter
  integer nP,nY,nPok,nYok,iPstar
  double precision :: pstar 
  double precision,allocatable :: m0(:,:,:), m1(:,:,:), mifront(:,:,:), &
       miback(:,:,:),T4_ext(:)
  ! --- local
  INTEGER m,iP, flag,iL,iY,iZ,jZ,im,iW
  integer,allocatable :: nZ(:)
  ! double precision :: haux(npP),TAUaux(npL,npP,npY),Tplus(npP,npY,npY), &
  ! Tminus(npP,npY,npY),xN(npP),yN(npP),TAUr(npY),wm(npY), wmT(npY),wp(npY),
  ! alpha(npY,npY),beta(npY,npY),gamma2(npY,npY),delta(npY,npY), &
  ! wgmatp(npY,npY), wgmatm(npY,npY), Yok(npY), Pok(npP)
  DOUBLE PRECISION H, addplus,addminus, &
       result1,result2,resaux, &
       fact, faux
  double precision,allocatable :: haux(:),TAUaux(:,:,:),Tplus(:,:,:), &
       Tminus(:,:,:),xN(:),yN(:),TAUr(:),wm(:),wmT(:),wp(:),&
       alpha(:,:),beta(:,:),gamma2(:,:),delta(:,:),wgmatp(:,:), &
       wgmatm(:,:), Yok(:), Pok(:)
  !---------------------------------------------------------------------
  allocate(nZ(nP))
  allocate(haux(nP))
  allocate(TAUaux(nL,nP,nY))
  allocate(Tplus(nP,nY,nY))
  allocate(Tminus(nP,nY,nY))
  allocate(xN(nP))
  allocate(yN(nP))
  allocate(TAUr(nY))
  allocate(wm(nY))
  allocate(wmT(nY))
  allocate(wp(npY))
  allocate(alpha(nY,nY))
  allocate(beta(nY,nY))
  allocate(gamma2(nY,nY))
  allocate(delta(nY,nY))
  allocate(wgmatp(nY,nY))
  allocate(wgmatm(nY,nY))
  allocate(Yok(nY))
  allocate(Pok(nP))
  error = 0
  ! generate auxiliary arrays haux & nZ
  DO iP = 1, nP
     ! parameter alowing for a line of sight terminating on the star
     ! H(x1,x2) is the step function.
     haux(iP) = H(P(iP),pstar)
     ! if (haux(iP).lt.dynrange) haux(iP)=0.0D0
     ! upper limit for the counter of z position
     nZ(iP) = nY + 1 - iYfirst(iP)
  END DO
  ! Using the local array TAUaux to avoid multiple calculations of the
  ! product
  DO iL = 1, nL
     DO iP = 1, nP
        DO iY = 1, nY
           TAUaux(iL,iP,iY) = ETAzp(iP,iY)*TAUtot(iL)
        END DO
     END DO
  END DO
  ! -- evaluate matrix elements --
  ! loop over wavelengths
  !$OMP PARALLEL DO FIRSTPRIVATE(alpha,beta,gamma2,delta,nZ) &
  !$OMP PRIVATE(iL,iY,iP,iZ,jZ,TAUr,Tplus,Tminus,wgmatp,wgmatm) &
  !$OMP PRIVATE(wmT,fact,wp,wm,addplus,addminus,im,m,xN,yN,faux,resaux,flag,result1,result2)
  DO iL = 1, nL
     ! radial optical depths
     DO iY = 1, nY
        TAUr(iY) = ETAzp(1,iY)*TAUtot(iL)
     END DO
     ! auxiliary arrays for given TAUr
     CALL MYSPLINE(TAUr,nY,alpha,beta,gamma2,delta)
     ! loop over impact parameters
     DO iP = 1, nP-1
        ! set T-s to 0
        DO iY = 1, nY
           DO iW = 1, nY
              Tplus(iP,iY,iW) = 0.0D+00
              Tminus(iP,iY,iW) = 0.0D+00
           END DO
        END DO
        ! generate weights matrices
        CALL WEIGHTS(TAUaux,iP,iL,nZ(iP),alpha,beta,gamma2,delta,wgmatp,wgmatm,nY)
        ! first position on the line of sight
        IF (YPequal(iP).EQ.1) THEN
           iY = iYfirst(iP)
           DO iW = 1, nY
              ! cummulative weights for parts II & III
              wmT(iW) = 0.0D+00
              DO jZ = 1, nZ(iP)-1
                 fact = dexp(-TAUaux(iL,iP,jZ))
                 wmT(iW) = wmT(iW) + fact * wgmatm(jZ,iW)
              END DO
              Tplus(iP,iY,iW) = Tplus(iP,iY,iW) + haux(iP)*wmT(iW)
              Tminus(iP,iY,iW) = Tminus(iP,iY,iW) + wmT(iW)
           END DO
        END IF
        ! loop over positions on the line of sight
        DO iZ = 2, nZ(iP)
           ! increase index for radial position
           iY = iYfirst(iP) + iZ - 1
           ! generate weights for this position
           DO iW = 1, nY
              wp(iW) = 0.0D+00
              wm(iW) = 0.0D+00
              wmT(iW) = 0.0D+00
              ! part I
              DO jZ = 2, iZ
                 fact = dexp(TAUaux(iL,iP,jZ)-TAUaux(iL,iP,iZ))
                 wp(iW) = wp(iW) + fact * wgmatp(jZ,iW)
              END DO
              ! part II & III
              DO jZ = 1, nZ(iP)-1
                 fact = dexp(-(TAUaux(iL,iP,iZ)+TAUaux(iL,iP,jZ)))
                 wmT(iW) = wmT(iW) + fact * wgmatm(jZ,iW)
              END DO
              ! part IV
              IF (iZ.LT.nZ(iP)) THEN
                 DO jZ = iZ, nZ(iP)-1
                    fact = dexp(-(TAUaux(iL,iP,jZ)-TAUaux(iL,iP,iZ)))
                    wm(iW) = wm(iW) + fact * wgmatm(jZ,iW)
                 END DO
              ELSE
                 wm(iW) = 0.0D+00
              END IF
              ! add contribution from this step
              addplus = wp(iW) + haux(iP)*wmT(iW)
              Tplus(iP,iY,iW) = Tplus(iP,iY,iW) + addplus
              addminus = wm(iW)
              Tminus(iP,iY,iW) = Tminus(iP,iY,iW) + addminus
           END DO
           ! end of loop over iZ
        END DO
        ! end of the impact parameter loop, iP
     END DO
     ! add points on the edge
     DO iW = 1, nY
        Tplus(nP,nY,iW) = 0.0D+00
        Tminus(nP,nY,iW) = 0.0D+00
     END DO
     ! ============================
     ! find mat(iL,iY,iW) -> angular (mu) integration
     ! loop over moments (without calculation of rad.pressure)
     DO im = 1, 2
        m = im - 1
        ! loop over radial positions
        DO iY = 1, nY
           ! generate mu arrray
           DO iP = 1, Plast(iY)
              xN(iP) = sqrt(1.0D+00-(P(iP)/Y(iY)*P(iP)/Y(iY)))
           END DO
           ! loop over local (radial) positions
           DO iW = 1, nY
              ! generate intensity array for NORDLUND
              DO iP = 1, Plast(iY)
                 ! 'faux' is a representation of (-1)**m
                 faux = 1.0D+00 - 2.0D+00*MOD(m,2)
                 yN(iP) = Tplus(iP,iY,iW) + faux*Tminus(iP,iY,iW)
                 ! store matrix elements to evaluate intensity (*1/4Pi)
                 IF (im.EQ.1.AND.iY.EQ.nY) THEN
!                    mifront(iL,iP,iW) = 0.0795775D+00 * Tplus(iP,iY,iW)
!                    miback(iL,iP,iW) = 0.0795775D+00 * Tminus(iP,iY,iW)
                    mifront(iL,iP,iW) = Tplus(iP,iY,iW)/(4.0D0*pi)
                    miback(iL,iP,iW) = Tminus(iP,iY,iW)/(4.0D0*pi)
                 END IF
              END DO
              ! angular integration inside cavity
              IF (pstar.GT.0.0D+00) THEN
                 CALL NORDLUND(nY,nP,0,xN,yN,1,iPstar,m,resaux)
                 !IF (error.NE.0) GOTO 999
                 IF (nPcav.GT.iPstar) CALL NORDLUND(nY,nP,0,xN,yN,iPstar+1,nPcav+1,m,result1)
                 !IF (error.NE.0) GOTO 999
                 result1 = result1 + resaux
              ELSE
                 CALL NORDLUND(nY,nP,0,xN,yN,1,nPcav+1,m,result1)
                 !IF (error.NE.0) GOTO 999
              END IF
              ! flag for analytic integration outside cavity
              IF (iY.GT.6) THEN
                 flag = 1
              ELSE
                 flag = 0
              ENDIF
              ! angular integration outside cavity
              IF (iY.GT.1) THEN
                 CALL NORDLUND(nY,nP,flag,xN,yN,nPcav+1,Plast(iY),m,result2)
                 !IF (error.NE.0) GOTO 999
              ELSE
                 result2 = 0.0D+00
              END IF
              ! store current matrix element
              IF (m.EQ.0) m0(iL,iY,iW) = 0.5D+00*Y(iY)*Y(iY)*(result1 + result2)
              IF (m.EQ.1) m1(iL,iY,iW) = 0.5D+00*Y(iY)*Y(iY)*(result1 + result2)
!              IF (m.EQ.0) m0(iL,iY,iW) = 0.5D+00*(result1 + result2)*Y(iY)*Y(iY)/pi ! seems to work for external illumination
!              IF (m.EQ.1) m1(iL,iY,iW) = 0.5D+00*(result1 + result2)*Y(iY)*Y(iY)/pi ! seems to work for external illumination
           END DO
        END DO
     END DO
     ! =============================
     ! end of loop over wavelengths
  END DO
  !$OMP END PARALLEL DO 
  ! save Y and P grids to Yok and Pok, they are needed for analysis
  ! in cases when requirement for finer grids cannot be satisfied and
  ! previous solution is used for output
  nYok = nY
  DO iY = 1, nY
     Yok(iY) = Y(iY)
  END DO
  nPok = nP
  DO iP = 1, nP
     Pok(iP) = P(iP)
  END DO
  !-----------------------------------------------------------------------
999 deallocate(nZ)
  deallocate(haux)
  deallocate(TAUaux)
  deallocate(Tplus)
  deallocate(Tminus)
  deallocate(xN)
  deallocate(yN)
  deallocate(TAUr)
  deallocate(wm)
  deallocate(wmT)
  deallocate(wp)
  deallocate(alpha)
  deallocate(beta)
  deallocate(gamma2)
  deallocate(delta)
  deallocate(wgmatp)
  deallocate(wgmatm)
  deallocate(Yok)
  deallocate(Pok)
  RETURN
END SUBROUTINE matrix
!***********************************************************************

!***********************************************************************
DOUBLE PRECISION FUNCTION H(x1,x2)
!=======================================================================
!This function calculates the step function: H=1 for x1 >= x2 and H=0
!for x1 < x2.           [Z.I., Nov. 1995]
!=======================================================================
  IMPLICIT none
  DOUBLE PRECISION x1, x2
  !-----------------------------------------------------------------------
  IF (x1.GE.x2) THEN
     H = 1.0D+00
  ELSE
     H = 0.0D+00
  END IF
  !-----------------------------------------------------------------------
  RETURN
END FUNCTION H
!***********************************************************************


!***********************************************************************
SUBROUTINE MYSPLINE(x,N,alpha,beta,gamma2,delta)
!=======================================================================
!This subroutine finds arrays alpha, beta, gamma and delta describing
!a cubic spline approximation of an unknown function f(x) given as an
!array f(i)=f(x(i)) with i=1..N. The cubic spline approximation is:
!f(x)=a(i) + b(i)*t + c(i)*t^2 + d(i)*t^3  for x(i).LE.x.LE.x(i+1)
!and t = (x-x(i))/(x(i+1)-x(i)), i=1..N-1. Coefficients a,b,c,d are
!equal to:
!a(i) = alpha(i,1)*f(1) + alpha(i,2)*f(2) + ... + alpha(i,N)*f(N)
!and b,c,d analogously.    [Z.I., Dec. 1995]
!=======================================================================
  use common
  use interfaces
  IMPLICIT none
  !--- parameter
  integer :: N,nY
  double precision,allocatable :: x(:),alpha(:,:),beta(:,:),gamma2(:,:),delta(:,:)
  !--- local
  INTEGER i, j, dummy, Kron
  double precision,allocatable :: secnder(:,:), yaux(:), deraux(:)
  DOUBLE PRECISION y2at1, y2atN, D
  EXTERNAL Kron
  
  allocate(secnder(N,N))
  allocate(yaux(N))
  allocate(deraux(N))
  ! -----------------------------------------------------------------------
  ! generate second derivatives, secnder(j,l)
  DO j = 1, N
     DO dummy = 1, N
        IF (dummy.EQ.j) THEN
           yaux(dummy) = 1.0D+00
        ELSE
           yaux(dummy) = 0.0D+00
        END IF
     END DO
     y2at1 = (yaux(2)-yaux(1))/(x(2)-x(1))
     y2atN = (yaux(N)-yaux(N-1))/(x(N)-x(N-1))
     CALL SPLINE(x,yaux,N,y2at1,y2atN,deraux)
     DO i = 1, N
        secnder(i,j) =  deraux(i)
        ! secnder(i,j) = 0.0
     END DO
  END DO
  ! generate alpha, beta, gamma, delta
  DO i = 1, N-1
     D = (x(i+1) - x(i))*(x(i+1) - x(i)) / 6.0D+00
     DO j = 1, N
        alpha(i,j) = Kron(i,j)*1.0D+00
        beta(i,j) = Kron(i+1,j) - Kron(i,j)
        beta(i,j) = beta(i,j)-D*(2.0D+00*secnder(i,j)+secnder(i+1,j))
        gamma2(i,j) = 3.0D+00 * D * secnder(i,j)
        delta(i,j) = D*(secnder(i+1,j)-secnder(i,j))
     END DO
  END DO
  !-----------------------------------------------------------------------
  deallocate(secnder)
  deallocate(yaux)
  deallocate(deraux)
  RETURN
END SUBROUTINE MYSPLINE
!***********************************************************************

!***********************************************************************
INTEGER FUNCTION Kron(i1,i2)
!=======================================================================
!This function is Kronecker delta-function defined as:
!Kron(i1,i2) = 1 for i1=i2
!Kron(i1,i2) = 0 otherwise.[Z.I., Dec. 1995]
!=======================================================================
  IMPLICIT none
  INTEGER i1, i2
  !-----------------------------------------------------------------------
  IF (i1.EQ.i2) THEN
     Kron = 1
  ELSE
     Kron = 0
  END IF
  !-----------------------------------------------------------------------
  RETURN
END FUNCTION Kron
!***********************************************************************


!***********************************************************************
SUBROUTINE WEIGHTS(TAUaux,iP,iL,nZ,alpha,beta,gamma2,delta,wgp,wgm,nY)
!=======================================================================
!This subroutine calculates weights wgp(iZ,iY) and wgm(iZ,iY) for
!integrations:
!INT(S(w)*exp(sign*ETAzp(iP,iZ')/w^2)dETAzp(iP,iZ')]
!from ETAzp(iP,iZ) to ETAzp(iP,iZ+1), where w is local radius
!corresponding to TAU(iP,iZ'), and sign=1 for wgp and -1 for wgm.
!Integrals are evaluated as:
!INT = wg(iZ,1)*S(1) + wg(iZ,2)*S(2) + ... + wg(iZ,nY)*S(nY) with
!iZ=1..nZ-1. The method is based on approximation of S by cubic spline
!in radial optical depth given through matrices alpha, beta, gamma and
!delta (see MYSPLINE).                         [ZI,Dec'95;MN,Sep'97]
!=======================================================================
  use common
  use interfaces
  IMPLICIT none
  !---parameter
  integer iP,iL,nZ,nY
  double precision,allocatable :: TAUaux(:,:,:),alpha(:,:), beta(:,:),&
       gamma2(:,:),delta(:,:),wgp(:,:), wgm(:,:)
  !---local
  INTEGER iW, iZ, j
  DOUBLE PRECISION  waux
  double precision,allocatable :: K1p(:),K2p(:), K3p(:),K4p(:),&
       K1m(:),K2m(:),K3m(:),K4m(:)
  allocate(K1p(nY))
  allocate(K2p(nY))
  allocate(K3p(nY))
  allocate(K4p(nY))
  allocate(K1m(nY))
  allocate(K2m(nY))
  allocate(K3m(nY))
  allocate(K4m(nY))
  !-----------------------------------------------------------------------
  ! generate integrals of 'TAUr**n'
  CALL Kint4(TAUaux,iP,iL,nZ,K1p,K2p,K3p,K4p,K1m,K2m,K3m,K4m)
  ! loop over position on the line of sight
  DO iZ = 1, nZ
     iW = iYfirst(iP) + iZ - 1
     ! loop over radial position
     DO j = 1, nY
        IF (iZ.GT.1) THEN
           waux = alpha(iW-1,j)*K1p(iZ) + beta(iW-1,j)*K2p(iZ)
           wgp(iZ,j)=waux + gamma2(iW-1,j)*K3p(iZ)+delta(iW-1,j)*K4p(iZ)
!           if (wgp(iZ,j).lt.dynrange) wgp(iZ,j) = 0.0D0
        ELSE
           wgp(1,j) = 0.0D+00
        END IF
        IF (iZ.LT.nZ) THEN
           wgm(iZ,j) = alpha(iW,j)*K1m(iZ) + beta(iW,j)*K2m(iZ)
           wgm(iZ,j) = wgm(iZ,j)+gamma2(iW,j)*K3m(iZ)+delta(iW,j)*K4m(iZ)
!           if (wgm(iZ,j).lt.dynrange) wgm(iZ,j) = 0.0D0
        ELSE
           wgm(nZ,j) = 0.0D+00
        END IF
     END DO
  END DO
  !-----------------------------------------------------------------------
  deallocate(K1p)
  deallocate(K2p)
  deallocate(K3p)
  deallocate(K4p)
  deallocate(K1m)
  deallocate(K2m)
  deallocate(K3m)
  deallocate(K4m)
  RETURN
END SUBROUTINE WEIGHTS
!***********************************************************************
!!$
!!$
!***********************************************************************
SUBROUTINE Kint4(TAUaux,iP,iL,nZ,K1p,K2p,K3p,K4p,K1m,K2m,K3m,K4m)
!=======================================================================
!For given wavelength (iL) and impact parameter (iP), this subroutine
!calculates integrals Knp and Knm defined as:
!              Knp(iZ)=INT[PHIn(tz)*exp(tz)*dtz]
!from tz1=TAUaux(iL,iP,iZ) to tz2=TAUaux(iL,iP,iZ+1) and analogously for
!Km with exp(tz) replaced by exp(-tz). Function PHIn is defined as
!x**(n-1)/Yloc^2, where Ylo!is the local radius corresponding to tz,
!and x measures relative radial tau: x = (rt - tL)/(tR-tL). Here rt is
!the radial optical depth corresponding to tz and tL and tR are radial
!optical depths at the boundaries of the integration interval:
!tL = TAUaux(iL,1,iZ) = rt(iZ) and tR = TAUaux(iL,1,iZ+1) = rt(iZ+1).
!Integration is performed in z space by Romberg integration implemented
!in subroutine ROMBERG2 (slightly changed version of 'qromb' from Num.
!Recipes).          [ZI,Feb'96;MN,Sep'97]
!=======================================================================
  use common
  use interfaces
  IMPLICIT none
  !---parameter
  INTEGER iP, iL, nZ,iic
  DOUBLE PRECISION,allocatable :: TAUaux(:,:,:), K1p(:),K2p(:),K3p(:),&
       K4p(:), K1m(:), K2m(:), K3m(:), K4m(:)
  !---local
  INTEGER iZ, iW1, iLaux
  double precision :: Rresult(8), Kaux(8), deltrton(4),tRL,paux,&
       w1,w2,wL, delTAUzp, z1, z2
  !-----------------------------------------------------------------------
  paux = P(iP)
  iLaux = iL
  ! iLaux is needed to avoid compiler errors since it is in COMMON
  ! /phi2/ (here and in 'TWOfun'), while iL is transferred as a
  ! argument loop over positions on the line of sight
  DO iZ = 1, nZ-1
     ! index for the local radial position (left boundary)
     iW1 = iYfirst(iP) + iZ - 1
     ! radii at the boundaries
     wL = Y(iW1)
     IF (iZ.EQ.1) THEN
        if (paux.GT.1.0D+00) then
           w1 = paux
        else
           w1 = 1.0D+00
        end if
     ELSE
        w1 = Y(iW1)
     END IF
     w2 = Y(iW1+1)
     z1 = dsqrt(DABS(w1*w1 - paux*paux))
     z2 = dsqrt(w2*w2 - paux*paux)
     ! radial tau-difference at the bound., scaled to tot. opt. depth
     tRL = TAUaux(iL,1,iW1+1)-TAUaux(iL,1,iW1)
     ! auxiliary quantity aux/tRL**(n-1)
     deltrton(1) = TAUtot(iL)
     DO iic= 1, 3
        deltrton(iiC+1) = deltrton(iiC)/tRL
     END DO
     ! delTAUzp is needed in PHIn fun's
     delTAUzp = TAUaux(iL,iP,iZ+1)-TAUaux(iL,iP,iZ)
     ! integrate this step for all 8 cases
     CALL ROMBERG2(z1,z2,Rresult,w1,wl,iW1,iLaux,delTAUzp,paux)
     ! generate output values
     DO iic= 1, 4
        Kaux(iiC) = Rresult(iiC) * deltrton(iiC)
        Kaux(iiC+4) = Rresult(iiC+4) * deltrton(iiC)
     END DO
     K1m(iZ) = Kaux(1)
     K2m(iZ) = Kaux(2)
     K3m(iZ) = Kaux(3)
     K4m(iZ) = Kaux(4)
     K1p(iZ+1) = Kaux(5)
     K2p(iZ+1) = Kaux(6)
     K3p(iZ+1) = Kaux(7)
     K4p(iZ+1) = Kaux(8)
  END DO
  ! set undefined elements to 0
  K1m(nZ) = 0.0D+00
  K2m(nZ) = 0.0D+00
  K3m(nZ) = 0.0D+00
  K4m(nZ) = 0.0D+00
  K1p(1) = 0.0D+00
  K2p(1) = 0.0D+00
  K3p(1) = 0.0D+00
  K4p(1) = 0.0D+00
  !-----------------------------------------------------------------------
  RETURN
END SUBROUTINE Kint4
!***********************************************************************

!***********************************************************************
SUBROUTINE ROMBERG2(a,b,ss8,w1,wl,iW1,iLaux,delTAUzp,paux)
!=======================================================================
!This subroutine performs Romberg integration of 8 functions calculated
!in trapzd2 (by calling subroutine TWOFUN) on interval [a,b].
!The results are returned in ss8(1..8). Desired accuracy accRomb is
!user supplied and comes through COMMON /numerics/ read in from
!'numerics.inc'. This subroutine is based on slightly changed versions
!of 'qromb' and 'qromo' from Numerical Recipes.
!                        [MN & ZI,Aug'96]
!=======================================================================
  use common
  use interfaces
  IMPLICIT NONE
  INTEGER fconv(8),JMAX,JMAXP,K,KM, J, idone, kaux, iW1,iLaux,iic
  PARAMETER (JMAX=50, JMAXP=JMAX+1, K=5, KM=K-1)
  DOUBLE PRECISION ss, ss8(8), S2D(8,JMAXP), h(JMAXP), sjKM(JMAXP),&
       a, b, h0, EPS_romb, dss, s8(8), chk(8), w1,wl,delTAUzp,paux
!-----------------------------------------------------------------------
  EPS_romb = accRomb
  h0 = 0.0D+00
  h(1)=1.0D+00
  ! intialize convergence flags
  DO iic = 1, 8
     fconv(iiC) = 0
  END DO
  ! integrate until all 8 intergrals converge
  idone = 0
  j = 0
  DO WHILE(idone.NE.1.and.j.LE.JMAX)
     j = j + 1
     ! integrate with j division points
     call trapzd2(a,b,s8,j,w1,wl,iW1,iLaux,delTAUzp,paux)
     DO iic = 1, 8
        S2D(iiC,j) = S8(iiC)
     END DO
     ! check if any of 8 integrals has converged
     IF (j.ge.K) THEN
        idone = 1
        DO iic = 1, 8
           IF (fconv(iiC).EQ.0) THEN
              ! generate array for polint
              DO kaux = 1, j
                 sjKM(kaux) = S2D(iiC,kaux)
              END DO
              ! predict the integral for stepsize h->h0=0.0
              CALL polint(h(j-KM),sjKM(j-KM),K,h0,ss,dss)
              IF (dabs(dss).le.EPS_romb*dabs(ss)) THEN
                 SS8(iiC) = ss
                 fconv(iiC) = 1
              ELSE
                 chk(iiC) = dabs(dss)/dabs(ss)
              END IF
           END IF
           idone = idone*fconv(iiC)
        END DO
     END IF
     h(j+1)=0.25D+00*h(j)
  END DO
  IF (j.GE.jMAX) THEN
     write(*,*)' Reached the limiting number of steps in ROMBERG2'
     write(*,*)'You might want to change accRomb in the input file'
  END IF
  !-----------------------------------------------------------------------
  RETURN
END SUBROUTINE ROMBERG2
!***********************************************************************


!***********************************************************************
SUBROUTINE trapzd2(a,b,s,n,w1,wl,iW1,iLaux,delTAUzp,paux)
!=======================================================================
!This function integrates prescribed 8 functions from z=a to z=b with n
!divisions and stores the results to s(1..8). It is a heavily modified
!version of subroutine 'trapzd' (Num.Rec.'92).        [MN & ZI, Aug'96]
!=======================================================================
  IMPLICIT none
  INTEGER it,iC,i,n,j,iW1,iLaux
  DOUBLE PRECISION s(8),a,b,funcx(8),funca(8),funcb(8),del,summ(8),&
       tnm, x, ff, gp, gm, w1,wl,delTAUzp,paux
!-----------------------------------------------------------------------
  IF (n.eq.1) then
     ! calculate auxiliary functions at a and at b
     CALL TWOFUN(a,ff,gp,gm,w1,wl,iW1,iLaux,delTAUzp,paux)
     funca(1) =  gm
     funca(5) =  gp
     DO iC= 2, 4
        funca(iC) = funca(iC-1) * ff
        funca(4+iC) = funca(3+iC) * ff
     END DO
     CALL TWOFUN(b,ff,gp,gm,w1,wl,iW1,iLaux,delTAUzp,paux)
     funcb(1) =  gm
     funcb(5) =  gp
     DO iC= 2, 4
        funcb(iC) = funcb(iC-1) * ff
        funcb(4+iC) = funcb(3+iC) * ff
     END DO
     ! calculate integrals for all 8 functions
     DO i = 1, 8
        s(i) = 0.5D+00*(b-a)*(funca(i)+funcb(i))
     END DO
  ELSE
     it=2**(n-2)
     tnm=1.0D+00*(it)
     del=(b-a)/tnm
     x=a+0.5D+00*del
     DO i=1,8
        summ(i)=0.0D+00
     END DO
     ! calculate contributions of all 'it' divisions
     DO j = 1, it
        ! auxiliary functions at x
        CALL TWOFUN(x,ff,gp,gm,w1,wl,iW1,iLaux,delTAUzp,paux)
        ! generate (8) integrated functions at x
        funcx(1) = gm
        funcx(5) = gp
        DO iC= 2, 4
           funcx(iC) = funcx(iC-1) * ff
           funcx(4+iC) = funcx(3+iC) * ff
        END DO
        DO i=1,8
           summ(i)=summ(i)+funcx(i)
        END DO
        !        next x
        x=x+del
     END DO
     !      evaluate new value of the integral for all 8 cases
     DO i=1,8
        s(i)=0.5D+00*(s(i)+(b-a)*summ(i)/tnm)
     END DO
  END IF
  !-----------------------------------------------------------------------
  RETURN
END SUBROUTINE trapzd2
!***********************************************************************

!***********************************************************************
SUBROUTINE TWOFUN(z,ff,gp,gm,w1,wl,iW1,iLaux,delTAUzp,paux)
!=======================================================================
!This function evaluates auxiliary functions needed in trapzd2.
!               [MN & ZI,Aug'96; MN,Sep'97]
!=======================================================================
!-----------------------------------------------------------------------
  use common
  use interfaces
  implicit none
  integer iic
  DOUBLE PRECISION w,w1,wl,paux,z,auxw,delTAUzp,etaloc,ff,gm,gm1,gp,gp1,pp,&
       IntETA_matrix
  INTEGER iLaux,iW1


  ! local radius
  w = dsqrt(paux*paux + z*z)
  IF (w.LT.w1) w = w1
  ! find local value for ETA function
  etaloc= 0.0D+00
  auxw = 1.
  DO iiC= 1, 4
     etaloc= etaloc+ ETAcoef(iW1,iiC)*auxw
     auxw = auxw/w
  END DO
  ! ff, i.e. radial optical depth:
  pp = 0.0D+00
  ff = IntETA_matrix(pp,iW1,wL,w)*TAUtot(iLaux)
  ! g functions:
  gp1 = dexp(IntETA_matrix(paux,iW1,w1,w)*TAUtot(iLaux)-delTAUzp)
  gm1 = dexp(-IntETA_matrix(paux,iW1,w1,w)*TAUtot(iLaux))
  gp = etaloc/w/w * gp1
  gm = etaloc/w/w * gm1
  !-----------------------------------------------------------------------
  RETURN
END SUBROUTINE TWOFUN
!***********************************************************************

!***********************************************************************
DOUBLE PRECISION FUNCTION IntETA_matrix(p2,iW1,w1,w)
!=======================================================================
!This function calculates the integral over the normalized dens. prof.
!along the line of sight with impact parameter p and between the points
!corresponding to y=w1 and y=w. The method used is spline approximation
!for normalized density distribution ETA and subsequent integration
!performed analytically by MAPLE (these results are given through
!soubroutine Maple3).                         [ZI,Feb'96,MN,Aug'97]
!=======================================================================
  use common
  use interfaces
  IMPLICIT none
  INTEGER iW1,iic
  DOUBLE PRECISION  p2, w1, w, aux(4), z, z1, aux1(4)
  !-----------------------------------------------------------------------

  z = dsqrt(w*w-p2*p2)
  z1 = dsqrt(w1*w1-p2*p2)
  !    integrals calculated by MAPLE
  CALL Maple3(w,z,p2,aux)
  CALL Maple3(w1,z1,p2,aux1)
  DO iiC = 1, 4
     aux(iiC) = aux(iiC) - aux1(iiC)
  END DO
  IntETA_matrix = 0.0D+00
  DO iiC = 1, 4
     IntETA_matrix = IntETA_matrix + ETAcoef(iW1,iiC) * aux(iiC)
  END DO
  !-----------------------------------------------------------------------
  RETURN
END FUNCTION IntETA_matrix

!***********************************************************************
!!$
!!$!***********************************************************************
!!$SUBROUTINE ANALYSIS_matrix(model,error)
!!$!=======================================================================
!!$!This subroutine analyzes the solution. It finds the flux conservation
!!$!accuracy and evaluates many output quantites (e.g. QF(y), TAUF(y),Psi, F1
!!$!the rad.pressure force, dynamical quantities etc.)
!!$!This is with new additions acc. to IE'00           [ZI,Mar'96;MN,Mar'99]
!!$!=======================================================================
!!$  use common
!!$  IMPLICIT none
!!$  DOUBLE PRECISION qaux(npL),qaux2(npL), K1(npY), K2(npY), QpTd(npG,npY),&
!!$       QpStar(npY), mx, aux, C1, C2, C3, delta, Eps1, Fi, Gie2000, maxFerr,&
!!$       q_star,ugas_out, QUtot1, resaux, S4, tauV, Tei, Teo, xP, Yok(npY), &
!!$       Pok(npP), planck, Us(npL,npY)
!!$  INTEGER iL, iY, model, iP, error,denstyp
!!$  !-----------------------------------------------------------------------
!!$  !    make sure that grids correspond to accepted solution
!!$  nY = nYok
!!$  DO iY = 1, nY
!!$     Y(iY) = Yok(iY)
!!$  END DO
!!$  nP = nPok
!!$  DO iP = 1, nP
!!$     P(iP) = Pok(iP)
!!$  END DO
!!$  !-------------
!!$  ! spectrum (flux at the outer edge as a function of wavelength)
!!$  DO iL = 1, nL
!!$     Spectrum(iL) = dabs(ftot(iL,nY))
!!$     ! added in version dusty16.for - to prevent taking log from zero
!!$     ! in Spectral [MN]:
!!$     IF (Spectrum(iL).LE.1.0D-20) Spectrum(iL) = 1.0D-20
!!$  END DO
!!$  !-------------
!!$  ! analyze bolometri!flux error (1/2 of the max spread of fbol)
!!$  CALL FindErr(nY,fbol,maxFerr)
!!$  ! find the flux averaged optical depth, tauF(y)
!!$  IF (denstyp.NE.0) THEN
!!$     ! for spherical shell
!!$     tauF(1) = 0.0D+00
!!$     DO iY = 2, nY
!!$        ! generate auxiliary function for integration:
!!$        ! loop over iL (wavelength)
!!$        ! N.B. the definition: ETAzp(1,y) = taur(y)/tauT so that
!!$        ! tau(iL,iY) = TAUtot(iL)*ETAzp(1,iY)
!!$        DO iL = 1, nL
!!$           qaux(iL)=TAUtot(iL)*ETAzp(1,iY)*dabs(ftot(iL,iY))/lambda(iL)
!!$        END DO
!!$        CALL Simpson(npL,1,nL,lambda,qaux,resaux)
!!$        ! tauF(iY) = <tau(iL,iY)*ftot(iL,iY)>
!!$        tauF(iY) = resaux
!!$     END DO
!!$     ! for full RDW calculation redo tauF to be consistent with CalcEta
!!$     IF (RDW) THEN
!!$        ! generate ETA and its integral (normalization constant)
!!$        DO iY = 1, nY
!!$           K1(iY) = vrat(1,iY)/ugas(iY)/Y(iY)/Y(iY)
!!$        END DO
!!$        CALL SIMPSON(npY,1,nY,Y,K1,resaux)
!!$        ! find tauF
!!$        DO iY = 1, nY
!!$           K2(iY) = qF(iY)*K1(iY)/resaux
!!$           CALL SIMPSON(npY,1,iY,Y,K2,aux)
!!$           tauF(iY) = TAUfid*aux
!!$        END DO
!!$     END IF
!!$  ELSE
!!$     ! for slab
!!$     tauF(1) = 0.0D+00
!!$     DO iY = 1, nY
!!$        ! generate auxiliary function for integration:
!!$        ! loop over iL (wavelength)
!!$        DO iL = 1, nL
!!$           qaux(iL)=TAUslb(iL,iY)*dabs(fTot(iL,iY))/lambda(iL)
!!$           CALL Simpson(npL,1,nL,lambda,qaux,resaux)
!!$           tauF(iY) = resaux
!!$        END DO
!!$     END DO
!!$  END IF
!!$  !-------------
!!$  ! ratio of gravitational to radiation pressure force (isotropi!
!!$  ! scattering) per unit volume
!!$  ! s4 = (L4sol/Msol)/(4*Pi*G*c*rho_s)/1D-6;
!!$  ! rho_s=3000 kg.m-3, grain radius 'a' is in microns, aveV=4/3*Pi*<a^3>
!!$  IF(denstyp.NE.0) THEN
!!$     s4 = 1.925D+00 / (4.0D+00*Pi*Gconst*3.0D+008*3000.0D+00*1.0D-06)
!!$     ! in case of sigma's from a file aveV=1 (initialized in GetOptPr)
!!$     DO iY = 1, nY
!!$        DO iL = 1, nL
!!$           qaux(iL)=(SigmaA(1,iL)+SigmaS(1,iL))/aveV * &
!!$                dabs(ftot(iL,iY))/lambda(iL)
!!$        END DO
!!$        CALL Simpson(npL,1,nL,lambda,qaux,resaux)
!!$        rg(1,iY) = s4 * resaux / r_gd
!!$        ! If dust drift (dynamics case):
!!$        IF (RDW) rg(1,iY) = rg(1,iY)*vrat(1,iY)
!!$        IF (iY.EQ.1) THEN
!!$           Phi = resaux
!!$        END IF
!!$     END DO
!!$     ! the terminal value of the reddening profile, normalized to y=1
!!$     Phi = resaux / Phi
!!$  END IF
!!$  !-------------
!!$  ! find the Planck averaged absorption efficiencies
!!$  DO iY = 1, nY
!!$     ! generate auxiliary function for integration over wavelengths:
!!$     DO iL = 1, nL
!!$        qaux(iL) = SigmaA(1,iL) * Us(iL,iY) / lambda(iL)
!!$        xP = 14400.0D+00 / Td(1,iY) / lambda(iL)
!!$        qaux2(iL) = SigmaA(1,iL) * Planck(xP) / lambda (iL)
!!$     END DO
!!$     CALL Simpson(npL,1,nL,lambda,qaux,resaux)
!!$     QpStar(iY) = resaux
!!$     CALL Simpson(npL,1,nL,lambda,qaux2,resaux)
!!$     QpTd(1,iY) = resaux
!!$  END DO
!!$  !----------
!!$  ! find parameter Psi (see Ivezi!& Elitzur, 1996)
!!$  ! generate auxiliary function for integration:
!!$  ! loop over iL (wavelength)
!!$  DO iL = 1, nL
!!$     qaux(iL) = SigmaA(1,iL) * Utot(iL,1) / lambda (iL)
!!$  END DO
!!$  CALL Simpson(npL,1,nL,lambda,qaux,resaux)
!!$  QUtot1 = resaux
!!$  Psi = QUtot1 / QpTd(1,1)
!!$  ! for slab Psi is defined by the flux at normal ill.
!!$  IF (SLB) Psi = dabs(mu1)*QUtot1 / QpTd(1,1)
!!$  !-------------
!!$  IF(denstyp.NE.0) THEN
!!$     ! ratio r1/r* (see Ivezi!& Elitzur, 1996, eq. 27)
!!$     r1rs = 0.5D+00 * dsqrt(Psi) * (Tstar(1) / Td(1,1))**2.0D+00
!!$     IF(Left.eq.0) r1rs = 1.0D+00
!!$  END IF
!!$  !-------------
!!$  ! Find epsilon - the relative contribution of the diffuse radiation
!!$  DO iY = 1, nY
!!$     aux = QpStar(iY)/QpTd(1,iY)/Psi*(Td(1,1)/Td(1,iY))**4.0D+00
!!$     IF (SLB) THEN
!!$        aux = aux*dabs(mu1)
!!$     ELSE
!!$        aux = aux/ Y(iY)/Y(iY)
!!$     END IF
!!$     Eps(iY) = 1.0D+00 - aux
!!$  END DO
!!$  Eps1 = 1.0D+00 - QpStar(1) / QUtot1
!!$  ! store these parameters in the storage array
!!$  SmC(1,model) = Psi
!!$  SmC(2,model) = Eps1
!!$  SmC(3,model) = QpStar(1)
!!$  SmC(4,model) = QpTd(1,1)
!!$  SmC(5,model) = maxFerr
!!$  !-------------
!!$  ! additional output quantities
!!$  ! bolometri!flux at r1 (in W/m2)
!!$  IF (typEntry(1).EQ.1) THEN
!!$     ! The constant is 4*Sigma*1000^4 (2.27E5 = 4*5.67D-08*1000**4)
!!$     Fi = 2.27D+5 / Psi * (Tsub(1)/1000.0D+00)**4.0D+00
!!$  ELSE
!!$     IF(Left.eq.0) THEN
!!$        Fi = Ubol(1)*4.0D+00*sigma*(Y(nY)*Teo**2.0D+00)**2.0D+00
!!$     ELSE
!!$        Fi = sigma*Tei**4.0D+00
!!$     END IF
!!$  END IF
!!$  ! inner radius (in cm) in case it is not an input
!!$  ! 5.53E16 = sqrt(10^4*Lo/4/Pi)
!!$  IF (typEntry(1).NE.3) THEN
!!$     IF(SLB) THEN
!!$        ! r1 is found from Fi = L/(4*pi*r1^2). Since in sub Input
!!$        ! Fi=Fi*mu1, here the mu1 dependence has to be removed
!!$        Cr1 =  5.53D+16 / dsqrt(Fi/abs(mu1))
!!$     ELSE
!!$        Cr1 = 5.53D+16 / dsqrt(Fi)
!!$     END IF
!!$  END IF
!!$  IF (denstyp.NE.0) THEN
!!$     ! angular diameter of inner cavity if Fbol=1D-6 W/m2
!!$     theta1 = 412.6D+00 / dsqrt(Fi)
!!$     ! check if the pt.source assumption is still obeyed
!!$     ! (only for BB-type spectrum including EM-function)
!!$     IF(startyp(1).eq.1.OR.startyp(1).eq.2) THEN
!!$        mx = sqrt(sqrt(Fi/sigma))
!!$        Te_min = 2.0D+00 * DMAX1(Td(1,1), mx)
!!$     END IF
!!$  END IF
!!$  IF (SLB) THEN
!!$     ! Teff for the left illuminating source in slab geometry
!!$     ! Teff = (Fi/sigma)^0.25D+00
!!$     SmC(7,model) = Tei
!!$     IF (ksi.GT.0.) THEN
!!$        ! Teff for the right illuminating source in slab geometry
!!$        SmC(8,model) = SmC(7,model)*sqrt(sqrt(ksi))
!!$     ELSE
!!$        SmC(8,model) = 0.0D+00
!!$     END IF
!!$  END IF
!!$  ! calculate conversion constants for dynamics
!!$  IF (RDWA.OR.RDW) THEN
!!$     ! for analytical approximation
!!$     ! (otherwise I1,2,3 are found in Gammafun)
!!$     IF (RDWA) THEN
!!$        I1 = 2.0D+00 * (1.0D+00-pow)/(1.0D+00+pow)/tauF(nY)
!!$        I2 = I1
!!$        I3 = I1 * tauF(nY) / TAUfid
!!$        Gamma(nY) = 0.5D+00
!!$     END IF
!!$     ! terminal expansion velocity, full formula:
!!$     ugas_out = tauF(nY) * (1.0D+00-Gamma(nY)) / (1.0D+00-pow)
!!$     ! The coefficients come from the units conversion
!!$     C1 = 0.2845D+00*TAUfid*sqrt(Psi)/I2/(SigExfid/aveV)*&
!!$          1.0D+006/Td(1,1)/Td(1,1)
!!$     C2 = 2.040D+00*ugas_out
!!$     C3 = 6.628D+00*I3*SigExfid/aveV*Gamma(nY)/I1
!!$     ! from version 2.0 stellar mass is defined as the maximal stellar
!!$     ! mass which does not quench the wind; the calculation is done
!!$     ! with half that mass since any smaller mass will have no effect
!!$     ! on the radial velocity and density profile (see IE2000)
!!$     ! n.b. erroneous Gamma(nY) is removed
!!$     CM = 6.628D+00*I3*SigExfid/aveV/I1
!!$     ! new definitions for output
!!$     ! mass-loss rate in Msol/yr
!!$     CMdot = 1.0D-05 * sqrt(C1)
!!$     ! terminal expansion velocity in km/s
!!$     Cve = 10.0D+00* C2 / sqrt(C1)
!!$     ! *** this is conversion to the nomenclature as in IE2001
!!$     IF (denstyp.EQ.6) THEN
!!$        ! IF (RDW) THEN
!!$        ! size averaged extinction efficiency
!!$        QV = SigExfid / aveA
!!$        tauV = TAUfid
!!$        q_star = qF(1)
!!$        zeta1 = vrat(1,1)
!!$        G1 = Gamma(1)
!!$        Ginf = Gamma(nY)
!!$        IF (G1.GT.0.0D+00) THEN
!!$           Gie2000 = 1.0D+00 / zeta1 / G1
!!$           delta = 1.0D+00 / (Gie2000 - 1.0D+00)
!!$        ELSE
!!$           delta = 0.0D+00
!!$        END IF
!!$        PIrdw = tauV / QV
!!$        Prdw = dsqrt(2.D+00*PIrdw/I2/QV/q_star)
!!$        winf = ugas_out / QV / q_star
!!$     END IF
!!$  END IF
!!$  !-----------------------------------------------------------------------
!!$  RETURN
!!$END SUBROUTINE ANALYSIS_matrix
!!$!***********************************************************************
!!$
!!$! ***********************************************************************
!!$SUBROUTINE SetGrids_matrix(pstar,iPstar,error,TAU)
!!$! =======================================================================
!!$! Sets the Y and P grids based on GrayBody flux conservation.
!!$!                                                     [MN & ZI, July'96]
!!$! =======================================================================
!!$  use common
!!$  IMPLICIT none
!!$  INTEGER error, iPstar, consfl
!!$  DOUBLE PRECISION pstar, Ugb(npY), fgb(npY), albedo,&
!!$       aux, faccs, TAU, accur, delTAUin, delTAUsc, facc
!!$
!!$! -----------------------------------------------------------------------
!!$  ! store the default value for delTAUsc and facc
!!$  faccs = facc
!!$  delTAUin = delTAUsc
!!$  ! change the delTAUsc seed for the initial Y grid if TAU is large
!!$  IF (TAU.LT.1.0D+00) delTAUsc = delTAUin * 2.0D+00
!!$  IF (TAU.GE.1.0D+00.and.TAU.LT.5.0D+00)delTAUsc = delTAUin*1.5D+00
!!$  IF (TAU.EQ.5.0D+00) delTAUsc = delTAUin
!!$  IF (TAU.GT.5.0D+0.and.TAU.LT.10.0D+0) delTAUsc = delTAUin/1.2D+00
!!$  IF (TAU.GE.10.0D+0.and.TAU.LT.20.0D+0)delTAUsc = delTAUin/1.3D+00
!!$  ! The grid is set with TAU=min{TAUlim,TAUmax}, so the lines below are obsol
!!$  ! IF (TAU.GE.20.0.and.TAU.LT.30.0) delTAUsc = delTAUin / 1.4
!!$  ! IF (TAU.GE.30.0.and.TAU.LT.50.0) delTAUsc = delTAUin / 1.5
!!$  ! IF (TAU.GE.50.0) delTAUsc = delTAUin / 2.0
!!$  ! for steep density distributions (RDW, including analyt.approximation):
!!$  IF (RDWA.OR.RDW) delTAUsc = delTAUin / 1.2D+00
!!$  ! change the facc seed for the initial Y grid if Yout is very small
!!$  IF (Yout.LT.1000.0D+00) facc = dsqrt(faccs)
!!$  IF (Yout.LT.100.0D+00) facc = dsqrt(facc)
!!$  IF (Yout.LT.10.0D+00) facc = dsqrt(facc)
!!$  IF (Yout.LT.2.0D+00) facc = dsqrt(facc)
!!$  IF (Yout.LT.1.2D+00) facc = dsqrt(facc)
!!$  IF (Yout.LT.1.05D+00) facc = dsqrt(facc)
!!$
!!$  albedo = 1.0D+00
!!$  aux = 1.0D+00
!!$  ! generate initial grids
!!$  CALL Ygrid(pstar,iPstar,error)
!!$ !!$  ! increase the grid if large tau and external illumination only
!!$ !!$  IF(Left.eq.0.AND.taumax.ge.99.0D+00) CALL DblYgrid(error)
!!$  IF (error.NE.0) goto 101
!!$  CALL Pgrid(pstar,iPstar,error)
!!$  IF (error.NE.0) goto 101
!!$  IF (iX.GE.1) THEN
!!$     write(18,'(a24,i3)')' Y grid generated, nY =',nY
!!$     write(18,'(a24,i3)')'                   nP =',nP
!!$     write(18,'(a24,i3)')'                 Nins =',Nins
!!$     write(18,'(a24,i3)')'                 Ncav =',Ncav
!!$  END IF
!!$  ! solve for gray body (i.e. pure scattering)
!!$  CALL GrayBody(albedo,TAU,Ugb,fgb)
!!$  IF (iVerb.EQ.2) write(*,*) 'Done with GrayBody'
!!$  ! find the max deviation of fgb (FindRMS called with flag 1)
!!$  ! (for grid generation purpose aux is set to 1.)
!!$  CALL FindRMS(1,fgb,aux,accur,nY)
!!$  IF (iX.GE.1) THEN
!!$     IF (accur.GT.accuracy) THEN
!!$        write(18,'(a25)')' Grids need improvement:'
!!$        write(18,'(a29,1p,e10.3)') &
!!$             '                   fTot(nY):',fgb(nY)
!!$        write(18,'(a29,1p,e10.3)')'      Single wavelength TAU:',TAU
!!$        write(18,'(a29,1p,e10.3)') &
!!$             '          Required accuracy:',accuracy
!!$     END IF
!!$     write(18,'(a29,1p,e10.3)')' Single wavelength accuracy:',accur
!!$  END IF
!!$  IF(accur.GT.accuracy) THEN
!!$     ! ChkFlux checks the bolometric flux conservation for the given
!!$     ! grid and decreases the step if conservation is not satisfactory
!!$     consfl = 5
!!$     CALL ChkFlux(fgb,accuracy,consfl,error)
!!$     IF (error.NE.0) goto 101
!!$     ! consfl=5 means everything was fine in ChkFlux
!!$     IF (consfl.EQ.5) THEN
!!$        IF (iX.GE.1) write(18,'(a23,i3)')' Y grid improved, nY =',nY
!!$        ! generate new impact parameter grid
!!$        CALL Pgrid(pstar,iPstar,error)
!!$        ! if P grid is not OK end this model
!!$        IF (error.NE.0) goto 101
!!$     ELSE
!!$        IF (iX.GE.1) THEN
!!$           write(18,'(a59,i3)') &
!!$                ' Although single wavelength accuracy was not satisfactory,'
!!$           write(18,'(a56,i3)') &
!!$                ' Y grid could not be improved because npY is too small.'
!!$           write(18,'(a58,i3)') &
!!$                ' Continuing calculation with a hope that it will be fine.'
!!$        END IF
!!$     END IF
!!$  END IF
!!$  ! return the default value for facc
!!$101 facc = faccs
!!$  delTAUsc = delTAUin
!!$  ! -----------------------------------------------------------------------
!!$  RETURN
!!$END SUBROUTINE SetGrids_matrix
!!$! ***********************************************************************
!!$
!!$
!!$!***********************************************************************
!!$SUBROUTINE GrayBody(albedo,TAUgbTot,Ugb,fgb)
!!$!=======================================================================
!!$!This subroutine solves the gray body problem for albedo=1 (or
!!$!equivalently pure scattering) and scattering with absorption (but no
!!$!emission) for albedo<1, in a spherically symmetri!envelope. Total
!!$!optical depth is TAUtot, and density law is specified elsewhere.
!!$!This subroutine was designed to be a part of Dusty and to use already
!!$!existing subroutines as much as possible, so some parts might seem to
!!$!be a little awkward.                           [ZI,Jul'96;MN,Sep'97]
!!$!=======================================================================
!!$  use common
!!$  IMPLICIT none
!!$  DOUBLE PRECISION Us(npL,npY),fs(npL,npY),Em(npG,npL,npY),omega(npL,npY), &
!!$       Ugb(npY), fgb(npY), albedo, Dummy1(npL,npP,npY),&
!!$       Dummy2(npL,npP,npY), Dummy3(npL,npY), TAUgbTot, TAUstore,&
!!$       mat0(npL,npY,npY), mat1(npL,npY,npY), pGB
!!$  INTEGER iPGB, iY, nLstore, error
!!$
!!$!----------------------------------------------------------------------
!!$!    Values needed in this subroutine only
!!$  pGB = 0.0D+00
!!$  iPGB = 0
!!$  nLstore = nL
!!$  nL = 1
!!$  TAUstore = TAUtot(1)
!!$  TAUtot(1) = TAUgbTot
!!$  ! generate spline coefficients for ETA
!!$  CALL setupETA
!!$  ! evaluate ETAzp
!!$  CALL getETAzp(ETAzp)
!!$  ! generate some temporary arrays
!!$  DO iY = 1, nY
!!$     Us(1,iY) = dexp(-ETAzp(1,iY)*TAUgbTot)
!!$     fs(1,iY) = Us(1,iY)
!!$     Em(1,1,iY) = 0.0D+00
!!$     fde(1,iY) = 0.0D+00
!!$     omega(1,iY) = albedo
!!$  END DO
!!$  ! find radiative transfer matrices
!!$  CALL Matrix(pGB,iPGB,mat0,mat1,Dummy1,Dummy2)
!!$  ! solve for Utot
!!$  CALL Invert(1,mat0,Us,Em,omat,Utot,error)
!!$  ! calculate flux, ftot
!!$  CALL Multiply(1,npY,nY,npL,nL,mat1,Utot,omat,1,fs,ftot)
!!$  ! store to the output arrays
!!$  DO iY = 1, nY
!!$     Ugb(iY) = Utot(1,iY)
!!$     fgb(iY) = ftot(1,iY)
!!$  END DO
!!$  nL = nLstore
!!$  TAUtot(1) = TAUstore
!!$  !-----------------------------------------------------------------------
!!$  RETURN
!!$END SUBROUTINE GrayBody
!!$!***********************************************************************
!!$
!!$
!!$!***********************************************************************
!!$SUBROUTINE FindRMS(typ,X,val,accur,N)
!!$!=======================================================================
!!$!Finds relative deviations 'accur' of an array X(N) from a given value val.
!!$!For typ=1 accur is maximal deviation, and for typ=2 the rms deviation.
!!$!                                                        [ZI'95; MN'99]
!!$!=======================================================================
!!$  IMPLICIT NONE
!!$  INTEGER N, i, typ
!!$  DOUBLE PRECISION X(N), val, accur, ss, dev
!!$!-----------------------------------------------------------------------
!!$  IF (typ.EQ.1) THEN
!!$     accur = 0.0D+00
!!$     DO i = 1, N
!!$        dev = (X(i)-val)/val
!!$        IF (DABS(dev).GT.accur) accur = DABS(dev)
!!$     END DO
!!$  ELSE
!!$     ss = 0.0D+00
!!$     DO i = 1, N
!!$        dev = X(i)-val
!!$        ss = ss + dev*dev
!!$     END DO
!!$     accur = dsqrt(ss/N/(N-1.0D+00))
!!$  END IF
!!$  !-----------------------------------------------------------------------
!!$  RETURN
!!$END SUBROUTINE FindRMS
!!$!***********************************************************************
!!$
!!$
!***********************************************************************
SUBROUTINE ChkFlux(nY,nYprev,flux,tolern,consfl,iterEta)
!=======================================================================
!Checks the bolometri!flux conservation at any point of a given Ygrid.
!In case of nonconservation increases the number of points at certain
!places. The current criterion is increasing the flux difference from
!tolern to its maximum value.                         [MN & ZI,July'96]
!=======================================================================
  use common
  use interfaces
  IMPLICIT none
  !---parameter
  DOUBLE PRECISION,allocatable :: flux(:)
  DOUBLE PRECISION :: tolern
  INTEGER consfl, nY,nYprev,iterEta
  !---local
  integer iY,k,kins,flag,istop,iDm
  integer,allocatable :: iYins(:)
  DOUBLE PRECISION delTAUMax,devfac,devmax,ee,ff,ffold,fmax,Yloc,ETA,fmed
  double precision,allocatable :: EtaTemp(:),Yins(:)
  !--------------------------------------------------------------------
  allocate(iYins(nY))
  allocate(EtaTemp(nY))
  allocate(Yins(nY))
  ! save old grid and values of Eta (important for denstyp = 5 or 6)
  IF (denstyp.eq.3) THEN !3(RDW)
     DO iY = 1, nY
        Yprev(iY) = Y(iY)
        EtaTemp(iY) = ETAdiscr(iY)
     END DO
     nYprev = nY
  END IF
  IF(Right.gt.0) THEN
     ! Find fmed - the median value of the bol.flux
     ! (if there is an external source fbol < 1)
     ! CALL SLBmisc(flux,fmax,fmed,AveDev,RMS,nY)
     print*,'! CALL SLBmisc(flux,fmax,fmed,AveDev,RMS,nY)'
  ELSE
     fmed = 1.0D+00
  END IF
  error = 0
  kins = 0
  devmax = 0.0D+00
  ! maximal delTAU is no more than 2 times the average value
  delTAUmax = 2.0D+00*TAUtot(1)*ETAzp(1,nY)/nY
  ! maximal deviation from fmed
  DO iY = 2, nY
     IF (dabs(flux(iY)-fmed).GT.devmax) devmax = dabs(flux(iY)-fmed)
  END DO
  ff = 0.0D+00
  istop = 0
  devfac = 0.1D+00
  ! search for places to improve the grid
  DO WHILE (istop.NE.1)
     DO iY = 2, nY
        ffold = ff
        ff = dabs(flux(iY) - fmed)
        flag = 0
        ! if any of these criteria is satisfied insert a point:
        ! 1) if error is increasing too fast
        IF (abs(ff-ffold).GT.devfac*devmax) flag = 1
        ! 2) if delTAU is too large
        IF (TAUtot(1)*(ETAzp(1,iY)-ETAzp(1,iY-1)).GT. &
             delTAUmax) flag = 1
        IF(flag.EQ.1.AND.devmax.GE.tolern) THEN
           kins = kins + 1
           Yins(kins) = Y(iY-1)+0.5D+00*(Y(iY)-Y(iY-1))
           iYins(kins) = iY-1
        END IF
     END DO
     IF (devmax.LT.tolern.OR.devfac.LT.0.01D+00) THEN
        istop = 1
     ELSE
        IF (kins.GT.0) istop = 1
     END IF
     devfac = devfac / 2.0D+00
  END DO
  IF (kins.EQ.0) THEN
     IF (consfl.NE.5) consfl = 1
  ELSE
     ! Add all new points to Y(nY). This gives the new Y(nY+kins).
     ! However, check if npY is large enough to insert all points:
     IF ((nY+kins).GT.npY) THEN
        ! consfl.EQ.5 is a signal that Chkflux was called from SetGrids,
        ! in this case continue without inserting new points. If this is
        ! full problem then give it up.
        IF (consfl.NE.5) THEN
           consfl = 1
        ELSE
           consfl = 7
           goto 777
        END IF
        IF (iX.GE.1) THEN
           write(18,*)' ****************     WARNING   ******************'
           write(18,*)'  The new Y grid can not accomodate more points!'
           write(18,'(a,i3)')'   Specified accuracy would require',nY+kins
           write(18,'(a,i3,a)')'   points, while npY =',npY,'.'
           write(18,*)'  For the required accuracy npY must be increased,'
           write(18,*)'  (see the manual S3.5 Numerical Accuracy).'
           write(18,*)' *************************************************'
        END IF
        kins = npY - nY
        error = 2
     END IF
     DO k = 1, kins
        CALL SHIFT(Y,nY,nY+k-1,Yins(k),iYins(k)+k-1)
     END DO
  END IF
  ! new size of the Y grid
  nY = nY + kins
  ! intepolate ETAdiscr to new Y grid for denstyp = 5 or 6
  DO iY = 1, nY
     Yloc = Y(iY)
     IF (iterETA.GT.1) THEN
        CALL LinInter(nY,nYprev,Yprev,EtaTemp,Yloc,iDm,ee)
        ETAdiscr(iY) = ee
     ELSE
        ETAdiscr(iY) = ETA(Yloc)
     END IF
  END DO
  !-----------------------------------------------------------------------
777 deallocate(iYins)
  deallocate(EtaTemp)
  deallocate(Yins)
  RETURN
END SUBROUTINE ChkFlux
!***********************************************************************

!***********************************************************************
SUBROUTINE SHIFT(X,Nmax,N,Xins,i)
!=======================================================================
!Rearranges a vector X by inserting a new element Xins.    [MN, Aug'96]
!=======================================================================
  implicit none
  integer Nmax, N, i,j
  DOUBLE PRECISION X(Nmax),Xins
!-----------------------------------------------------------------------
  DO j = N+1, i+2, -1
     x(j) = x(j-1)
  END DO
  x(i+1) = xins
  !-----------------------------------------------------------------------
  RETURN
END SUBROUTINE SHIFT
!***********************************************************************
!!$
!!$!***********************************************************************
!!$SUBROUTINE Emission_matrix(geom,flag,nG,Uin,Emiss)
!!$!=======================================================================
!!$!This subroutine calculates emission term from the temperature and abund
!!$!arrays for flag=0, and adds U to it for flag=1.
!!$!                                                     [Z.I., Mar. 1996]
!!$!=======================================================================
!!$  use common
!!$  IMPLICIT none
!!$  INTEGER iL,iY,iG,nG, flag, geom
!!$  DOUBLE PRECISION TT,Emiss(npL,npY), EmiG, xP, Planck, Uin(npL,npY), &
!!$       Tei, Teo
!!$!-----------------------------------------------------------------------
!!$  Tei = Ji*4*pi/sigma
!!$  Teo = Jo*4*pi/sigma
!!$  ! first initialize Emiss
!!$  ! loop over wavelengths
!!$  DO iL = 1, nL
!!$     ! loop over radial coordinate
!!$     DO iY = 1, nY
!!$        Emiss(iL,iY) = 0.0D+00
!!$     END DO
!!$  END DO
!!$  ! calculate emission term for each component and add it to Emiss
!!$  ! loop over wavelengths
!!$  DO iL = 1, nL
!!$     ! loop over radial coordinate
!!$     DO iY = 1, nY
!!$        ! loop over grains
!!$        DO iG = 1, nG
!!$           xP = 14400.0D+00 / lambda(iL) / Td(iG,iY)
!!$           IF(geom.NE.0) THEN
!!$              IF(Left.eq.0) THEN
!!$                 TT = (Td(iG,iY)**4.0D+00)/&
!!$                      (0.25D+00*Tei**4.0D+00+Y(nY)*Y(nY)*Teo**4.0D+00)
!!$              ELSE
!!$                 TT = (Td(iG,iY)**4.0D+00)/&
!!$                      (0.25D+00*Tei**4.0D+00+Teo**4.0D+00)
!!$              END IF
!!$              TT = TT * Y(iY)**2.0D+00
!!$           ELSE
!!$              TT = 4.0D+00 * (Td(iG,iY)/Tei)**4.0D+00
!!$           END IF
!!$           EmiG = abund(iG,iY) * TT * Planck(xP)
!!$           ! add contribution for current grains
!!$           Emiss(iL,iY) = Emiss(iL,iY) + EmiG
!!$        END DO
!!$        ! if needed add Uin
!!$        IF (flag.EQ.1) THEN
!!$           Emiss(iL,iY) = Emiss(iL,iY) + Uin(iL,iY)
!!$        END IF
!!$        IF (Emiss(iL,iY).LT.dynrange*dynrange) Emiss(iL,iY) = 0.0D+00
!!$     END DO
!!$  END DO
!!$  !-----------------------------------------------------------------------
!!$  RETURN
!!$END SUBROUTINE Emission_matrix
!!$!***********************************************************************


! *************************************************************************
SUBROUTINE FindErr(nY,flux,maxFerr)
!========================================================================
!This subroutine finds maximum err in flux conservation for both
!spherical and slab case as (fmax-fmin)/(fmax+fmin)   [MN,Aug'99]
!=========================================================================
  use common
  use interfaces
  IMPLICIT none
  !--- parameter
  integer nY
  DOUBLE PRECISION maxFerr
  double precision,allocatable :: flux(:)
  !--- local
  INTEGER iY
  DOUBLE PRECISION fmin, fmax, aux, accFbol, tune_acc
!---------------------------------------------------------------------
  ! accFbol = 1e-3 of input flux
  ! accFbol = max(Ji,Jo)*4*pi*1.0d-03
  tune_acc = 1e-1
  if (slb) then 
     accFbol = min(Ji,Jo)*4*pi*accFlux*tune_acc / Jext(1)
  else if (sph) then
     ! in the spherical case the flux is only zero if there is zero no
     ! central source (Ji=0) therefor accFbok is changed!
     accFbol = max(Ji,Jo)*4*pi*accFlux*tune_acc / max(Jext(1),Jext(nY))
  endif
  ! Find the min and max of fbol values
  ! The abs and lower limit on fbol are protection for the case
  ! of completely symmetri! slab illumination. The lower limit
  ! is bound by the numerical accuracy of the flux calculation
  fmin = 1.e5
  fmax = 0.
  DO iY = 1, nY
     aux = flux(iY)
!     IF (ksi.eq.1.0) aux = dabs(aux)
     IF (dabs(aux).LE.accFbol) aux = accFbol
     IF(aux.LT.fmin) fmin = aux
     IF(aux.GT.fmax) fmax = aux
     !print*,iY,(fmax - dabs(fmin))/(fmax + dabs(fmin))
  END DO
  if (fmax.LE.0.) then
     ! bad solution; overall flux cannot be negative
     maxFerr = 1.
  else if ((fmax.eq.fmin).and.(fmax.eq.accFbol)) then
     maxFerr = accFlux*tune_acc
  else
     maxFerr = (fmax - dabs(fmin))/(fmax + dabs(fmin))
  end if
  ! ----------------------------------------------------------------------
  RETURN
END SUBROUTINE FindErr
! *************************************************************************
