!***********************************************************************
subroutine Kernel(path,lpath,tau,Nmodel)
!=======================================================================
! This subroutine generates solution for all optical depth steps starting
! with intial optical depth given in the input file. It prints output
! only for input models
!                                               [MN, May'10]
!=======================================================================
  use common
  implicit none
  integer iG, iY, iL, i,j, GridType, model, Nmodel, itereta, &
       iterfbol, fbolOK, istop, lpath, nY, nYprev, nP, nCav, nIns
  double precision ratio, delta, tau0, taumax
  double precision,allocatable :: tau(:)
  character*235 path
  logical initial
  !----------------------------------------------------------------------
  istop = 0
  model = 1
  i = 1
  j = 0
  ! loop over optical depths
  do while (istop.eq.0)
     ! intialize print and initial flag
     j = j + 1
     if (j.eq.1) then
        initial = .true.
     else
        initial = .false.
     end if
     tau0 = tau(model)
     taufid = tau0
     call GetTaumax(tau0,taumax)
     if (iVerb.gt.0)  write(*,'(a9,i4,a6,f12.4)') ' model = ',model,', tau=',tau0
     call OPPEN(model,path,lpath)
     if (iVerb.eq.2) write(*,*) ' going to Solve '
     ! solve radiative transfer for this particular optical depth
     if (sph_matrix) then 
        CALL Solve_matrix(model,taumax,nY,nYprev,itereta,nP,nCav,nIns,initial,delta,iterfbol,fbolOK)
     else
        call Solve(model,taumax,nY,nYprev,itereta,nP,nCav,nIns,initial,delta,iterfbol,fbolOK)
     end if
     ! old dustys way
     ! if flux is conserved, the solution is obtained. So write out the values
     ! to output files for specified models
     if (fbolOK.eq.1) then
        if (error.eq.0) then
           ! call Spectral(model) ! no more spectral props.
           ! if (iVerb.eq.2) write(*,*) 'Done with Spectral'
           call PrOut(nY,nP,nYprev,itereta,model,delta)
           print*, 'Done with prOut'
        else
           go to 10
        end if
        call CLLOSE(model,Nmodel)
        if (iVerb.eq.2) write(*,*) ' ----------- '
        if (iVerb.eq.2) write(*,*) '             '
        model = model + 1
     end if
     if ((model-1).eq.Nmodel.and.fbolOK.eq.1) then
        istop = 1
     end if
     ! end of loop over tau
  END DO
!-----------------------------------------------------------------------
10  return
end subroutine Kernel
!***********************************************************************

!***********************************************************************
subroutine GetTauMax(tau0,taumax)
!=======================================================================
! This subroutine finds tautot(nL) and its max. value taumax.
!                                                      [MN, May'10]
!=======================================================================
  use common
  implicit none
  !---parameter
  double precision tau0,taumax
  !---local variables
  integer iL
  double precision :: sigAfid, sigSfid !x-section fiduc. wavel.
  double precision :: tau_tmp
  double precision,allocatable :: faux1(:), faux2(:)
  !---------------------------------------------------------------------
  allocate(faux1(nL))
  allocate(faux2(nL))
  do iL = 1, nL
     faux1(iL) = sigmaA(nG+1,iL)
     faux2(iL) = sigmaS(nG+1,iL)
  end do
  if (lamfid.lt.lambda(1)) then
     write(12,*)' fiducial wavelength was too small.'
     write(12,'(a8,e9.3,a17)')' using ',lambda(1),' micron instead.'
  end if
  if (lamfid.gt.lambda(nL)) then
     write(12,*)' fiducial wavelength was too large.'
     write(12,'(a8,e9.3,a17)')' using ',lambda(nL),' micron instead.'
  end if
  call lininter(nL,nL,lambda,faux1,lamfid,iLfid,sigAfid)
  call lininter(nL,nL,lambda,faux2,lamfid,iLfid,sigSfid)
  ! Extinction efficiency at fiducial wavelength
  sigexfid = sigAfid + sigSfid
  taumax = 0.0d0
  do iL = 1, nL
     tautot(iL) = tau0*(SigmaA(nG,iL) + SigmaS(nG,iL))/SigExfid
     if (tautot(iL).ge.taumax) taumax = tautot(iL)
  end do
  !-------------------------------------------------------------------
  deallocate(faux1)
  deallocate(faux2)
  return
end subroutine GetTauMax
!***********************************************************************

!***********************************************************************
subroutine Solve(model,taumax,nY,nYprev,itereta,nP,nCav,nIns,initial,&
  delta,iterfbol,fbolOK)
!=======================================================================
! This subroutine solves the continuum radiative transfer problem in
! spherical and planar geometry.                       [Deka,'09, MN'09]
!=======================================================================
  use common
  implicit none
  INTERFACE
     subroutine Rad_Transf(initial,nY,nYprev,nP,itereta,pstar,y_incr,us,fs,emiss, &
     iterfbol,initTemp,T4_ext)
       use common
       logical, intent(in) :: initial,initTemp
       integer, intent(in) :: y_incr,iterfbol
       integer :: nY,nP,nYprev,itereta
       double precision pstar
       double precision,allocatable :: us(:,:), fs(:,:) 
       double precision,allocatable :: T4_ext(:)
       double precision,allocatable :: emiss(:,:,:)
     end subroutine Rad_Transf
     SUBROUTINE SPH_Int(nY,nP,fs)
       integer nY,nP
       double precision,allocatable :: fs(:,:)
     END SUBROUTINE SPH_Int
     subroutine SLBdiff(nY,flag,grid,T4_ext,em,fp,fm)
       !---parameter
       integer :: nY,flag
       double precision,allocatable :: grid(:,:),T4_ext(:),em(:,:,:),&
            fp(:,:),fm(:,:)
     end subroutine SLBdiff
     subroutine Bolom(q,qbol,nY)
       integer nY
       double precision, allocatable :: q(:,:), qbol(:)
     end subroutine Bolom
     subroutine add2(nY,flxs,flxe,fbsum)
       integer :: nY
       double precision, allocatable :: flxs(:,:),flxe(:,:),fbsum(:)
     end subroutine add2
     subroutine add(np1,nr1,np2,nr2,q1,q2,q3,qout)
       integer np1, nr1, np2, nr2
       double precision, allocatable :: q1(:,:), q2(:,:), q3(:,:),qout(:,:)
     end subroutine add
     subroutine Flux_Consv(nY,nYprev,Ncav,itereta,iterfbol, fbolom,fbol_em,fbol_sc,fbolOK,maxrat)
       integer nY,nYprev,itereta,fbolOK,iterfbol
       double precision :: maxrat
       double precision, allocatable :: fbolom(:),fbol_em(:),fbol_sc(:)
     end subroutine Flux_Consv
     subroutine SLBintensity(nY,em)
       integer :: nY
       double precision, allocatable :: em(:,:,:)
     end subroutine SLBintensity
     subroutine Analysis(nY,model,us,T4_ext,delta,maxrat)
       integer nY,model
       double precision :: delta,maxrat
       double precision, allocatable ::  us(:,:),T4_ext(:)
     end subroutine Analysis
     SUBROUTINE FindErr(nY,flux,maxFerr)
       integer nY
       DOUBLE PRECISION maxFerr
       double precision,allocatable :: flux(:)
     END SUBROUTINE FindErr
  END INTERFACE
  !--- parameter 
  integer :: model,iterfbol,fbolOK,itereta, nY, nP, nYprev, nCav, nIns
  double precision :: delta, taumax
  logical initial
  !--- local variables
  integer etaOK, nY_old, iPstar, y_incr, iY, iL, imu
  double precision taulim,pstar, maxFerr, maxFerr_old, maxrat
  double precision, allocatable :: fs(:,:),us(:,:),T4_ext(:), &
       emiss(:,:,:),fDebol(:),fDsbol(:), &
       fdsm(:,:),fdsp(:,:),fdem(:,:),fdep(:,:)
  logical initTemp
  !---------------------------------------------------------------------
  allocate(fs(nL,npY))
  allocate(us(nL,npY))
  allocate(T4_ext(npY))
  allocate(emiss(nG,nL,npY))
  allocate(fDebol(npY))
  allocate(fDsbol(npY))
  allocate(fdsm(nL,npY))
  allocate(fdsp(nL,npY))
  allocate(fdem(nL,npY))
  allocate(fdep(nL,npY))
  if (iX.ne.0) then
     call line(0,2,18)
     write(18,'(a7,i3,a20)') ' model ',model,'  RUN-TIME MESSAGES '
     call line(0,1,18)
     write(18,*)'===========  starting to solve  ==========='
  end if
  error = 0
  if(sph) then
     ! solve for spherical envelope:
     ! temporarily the star is approximated by a point source
     pstar = 0.0d0
     iPstar = 1
     ! select optical depth for the grid calculation
     ! based on dynamical range
     taulim = 0.5d0*dlog(1.0d0/dynrange)
     ! if actual maximal tau is smaller use that value
     if (taumax.lt.taulim) taulim = taumax
     ! counter over eta (for radiatively driven winds only)
     itereta = 0
     EtaOK = 0
  end if
  ! counter for iterations over bolometric flux conservation
  iterfbol = 0
  fbolOK = 0
  itereta = 0
  EtaOK = 0
  y_incr = 0
  nCav = 0
  nP = 0
  maxFerr = 0.0
  if (model.eq.1) then 
     initTemp = .true.
  else
     initTemp = .false.
  end if
  !------------ loop over bol.flux conservation ------------
  do while(fbolOK.eq.0)
     iterfbol = iterfbol + 1
     if (iX.ge.1) then
        write(18,*)'  ',iterfbol,' iteration over fbol'
     end if
     if (iVerb.eq.2) then
        write(*,'(a14,i3,a20)') ' In Solve: ',iterfbol,' iteration over fbol '
     end if
     ! set grids for the initial optical depth
     call SetGrids(pstar,iPstar,taumax,nY,nYprev,nP,nCav,nIns,initial,iterfbol,itereta)
     ! assign number of grid points to nY_old
     nY_old = nY
     if (iVerb.eq.2) write(*,'(a19,i3,a12)') '  Calculating with ', nY,' grid points.'
     if (iX.ge.1) then
        write(18,'(a19,i3,a12)') '  Calculating with ', nY,' grid points.'
     end if
     if (nY.gt.npY) then
        if (iVerb.eq.2) write(*,*)' npY needs to be increased!'
        if (iX.ge.1) then
           write(18,*) ' ********** message from solve *********'
           write(18,*) ' npY has to be increased in file userpar.inc'
           write(18,*) ' ***************************************'
           go to 999
        end if
     end if
     ! solve the radiative transfer problem
     call Rad_Transf(initial,nY,nYprev,nP,itereta,pstar,y_incr,us,fs,emiss,iterfbol,initTemp,T4_ext)
     if (iVerb.eq.2) write(*,*)' Done with radiative transfer. '
     ! calculate diffuse flux
     ! for slab
     if(slb) then
        ! find the diffuse scattered flux(fl=1 for scatt. and fl=0 is for emission)
        call SLBdiff(nY,1,TAUslb,T4_ext,emiss,fdsp,fdsm)
        ! find the diffuse emitted flux
        call SLBdiff(nY,0,TAUslb,T4_ext,emiss,fdep,fdem)
        ! find bolometric diffuse flux
        call add2(nY,fdsp,fdep,fpbol)
        call add2(nY,fdsm,fdem,fmbol)
        ! find calculated bolometric diffuse flux
        !  for sphere
     elseif(sph) then
        !!** added here, instead of before Analysis [Sep.1'10]
        call add(nY,nY,nL,nL,fs,fds,fde,ftot)
     end if
     ! Find accfbol as in old Dusty versions.
     call BOLOM(fs,fsbol,nY)
     call BOLOM(fde,fDebol,nY)
     call BOLOM(fds,fDsbol,nY)
     DO iY = 1, nY
        fbol(iY) = fDebol(iY)+fDsbol(iY)+fsbol(iY)
     END DO
!!$     DO iY = 1, nY
!!$        write(7777,*) Y(iY),fDebol(iY),fDsbol(iY),fsbol(iY),Td(1,iY)
!!$     END DO
!!$     close(7777)
     maxFerr_old = maxFerr
     call Finderr(nY,fbol,maxFerr)
     if ((maxFerr-maxFerr_old.gt.0.).and.(.not.initTemp)) then
        initTemp=.true.
        print*,'true'
     else 
        initTemp=.false.
        print*,'false'
     end if
     if(iVerb.eq.2) write(*,'(A,F10.3,A)') '  Achieved error in bolometric Flux:',maxFerr*100,'%'
     call Flux_Consv(nY,nYprev,Ncav,itereta,iterfbol,fbol,fDebol,fDsbol,fbolOK,maxrat)
     if (iVerb.eq.2) write(*,'(a,i3,a,i3,a,i3)') '  After Flux_Cons nY=', nY,' nP=',nP,' Ncav=',Ncav
     ! initialize flag, y_incr, which records whether there is an increase in
     ! y-grid points
     y_incr = 0
     ! if new number of grid points is greater than the old, then taux = 1
     if (nY.gt.nY_old) y_incr = 1
     ! if the grid size limit is reached error=2
     if (error.eq.2.and.iterfbol.eq.1) then
        ! if this is the first calculation end this model
        if(iX.ge.1) then
           write(18,*)' =========== IMPORTANT WARNING =========== '
           write(18,*)' The limit for grid size is already reached'
           write(18,*)' and flux conservation can not be improved.'
           write(18,*)' Treat all results with caution!'
           write(18,'(a20,1p,F6.2,a2)')'  Achieved accuracy:',maxrat*100., ' %'
        end if
        error = 0
        fbolOK = 2
     end if
     ! if this is a higher iteration use previous solution
     if (error.eq.2) then
        if (iX.ge.1.and.iterfbol.ge.1) then
           write(18,*)' ======= IMPORTANT WARNING ======== '
           write(18,*)' In trying to conserve fbol, '
           write(18,*)' the limit for grid sizes is reached.  '
           write(18,'(a20,1p,F6.2,a2)')'  Achieved accuracy:',maxrat*100., ' %'
           write(18,*)' Treat all results with caution!'
        end if
        error = 0
        fbolOK = 2
     end if
     ! if fbol not conserved try again with a finer grid
     if (fbolOK.eq.0 .and. iX.ge.1) then
        write(18,*)'  ******** MESSAGE from Solve ********'
        write(18,'(a20,1p,e12.3,a2)')'  Achieved accuracy:',maxrat*100., ' %'
        write(18,*)'  Trying again with finer grids'
     end if
  end do !  end of loop over flux conservion
  if ((denstyp.eq.3).and.sph) then !3(RDW)
     ! counter over eta (for radiatively driven winds only)
     ! do while (EtaOK.eq.0)
     itereta = itereta + 1
     if (iX.ne.0.and.(denstyp.eq.3)) then
        write(18,*)'----------------------------------------'
        write(18,*)' ',itereta,'. iteration over eta'
     end if
     if (iVerb.eq.2.and.(denstyp.eq.3)) write(*,*) ' ',itereta,'. iteration over eta'
     ! for winds check if eta has converged...
     if (((denstyp.eq.3)).and.fbolOK.ne.2) then
        ! ptr(2) is specified in input and controls converg. crit.
        if (ptr(2).lt.1.0e-6.and.itereta.gt.2)then
           EtaOK = 1
        else
           call Winds(nY,nYprev,EtaOK)
        end if
        if (itereta.gt.10.and.EtaOK.eq.0) then
           EtaOK = 2
           if (iX.ne.0) then
              write(18,*)' *********  WARNING  *********'
              write(18,*)' Could not converge on eta in '
              write(18,*)' 10 iterations.'
              write(18,*)' *****************************'
           end if
        end if
        ! ...or otherwise finish right away
     else
        EtaOK = 1
     end if
     ! if(EtaOK.eq.0) go to 10
  end if
  ! error=4 means npY not large enough for oblique illumination grid
  if (error.eq.4) then
     call msg(15)
     goto 999
  end if
  !    add diffuse flux to the transmitted flux to find total flux
  !!** sub Add has to be called before WINDS! ftot is needed there! [MN]
  !!**     call add(npY,nY,npL,nL,fs,fds,fde,ftot)
  call bolom(utot,ubol,nY)
  if (sph) then
     !     calculate intensity (at the outer edge) if required
     if (iC.ne.0) then
        if (iX.ne.0) write(18,*) 'Calculating intensities.'
        call sph_int(nY,nP,fs)
        if (iVerb.eq.2) write(*,*) ' Done with finding intensity '
     end if
     ! if needed convolve intensity with the psf
!!$     if (iPsf.ne.0) then
!!$        call convolve
!!$        if (iVerb.eq.2) write(*,*) ' Done with Convolve. '
!!$     end if
!!$     ! if needed find the visibility function
!!$     if (iV.ne.0) then
!!$        call visibili
!!$        if (iVerb.eq.2) write(*,*) 'Done with Visibili.'
!!$     end if
  else if(slb) then
     ! if emerging intensity is required (for left side ill.only):
     if(iC.gt.0) then
        if (iVerb.eq.2) write(*,*)' Calculating intensities. '
        call SLBintensity(nY,emiss)
        do iL = 1, nL
           do imu = 1, nmu
              SLBintm(imu,iL) = SLBintm(imu,iL)/(4.0d0*pi) ! needed to change later [MD]
              SLBintp(imu,iL) = SLBintp(imu,iL)/(4.0d0*pi) ! needed to change later
              !iauxr(iL) = SLBintm(imu,iL)/lambda(iL)
              !iauxl(iL) = SLBintp(imu,iL)/lambda(iL)
           end do
        end do
        if (iVerb.eq.2) write(*,*) ' Done with finding intensity. '
     end if
  end if       !end if for sphere or slab
  ! analyze the solution and calculate some auxiliary quantities
  CALL analysis(nY,model,us,T4_ext,delta,maxrat)
  if (iVerb.eq.2) write(*,*) 'Done with Analysis.'
  if (iX.ne.0) then
     write(18,'(a36,F5.1,a2)') '  Accuracy of computed bolom. flux: ',maxFerr*100.0,' %'
     write(18,*)'  ==== SOLVE successfully completed ====='
     write(18,*)'  ======================================='
  end if
  !---------------------------------------------------------------------
999 deallocate(fs)
  deallocate(us)
  deallocate(T4_ext)
  deallocate(emiss)
  deallocate(fDebol)
  deallocate(fDsbol)
  deallocate(fdsm)
  deallocate(fdsp)
  deallocate(fdem)
  deallocate(fdep)
 return
end subroutine solve
!***********************************************************************

!***********************************************************************
subroutine SetGrids(pstar,iPstar,taumax,nY,nYprev,nP,nCav,nIns,initial,&
     iterfbol,itereta)
!=======================================================================
! Sets the y and p grids for sphere and tau-grid for slab.
!                                                     [MN & ZI, July'96]
!=======================================================================
  use common
  implicit none
  integer iL,iY, iPstar, itereta, iterfbol,taux, nY,nYprev, nP, nCav, nIns
  double precision pstar,taumax
  double precision, allocatable :: tau(:)
  logical initial
  !-----------------------------------------------------------------------
  allocate(tau(npY))
  if (sph) then
     if (initial.and.iterfbol.eq.1) then
        ! generate y-grid
        call Ygrid(pstar,iPstar,itereta,iterfbol,taumax,nY,nYprev,nP,nCav,nIns)
        ! generate p and z grid
        call Pgrid(pstar,iPstar,nY,nP,nCav,nIns)
     else
        ! redefine p and z grid for the new y-grid
        IF ((TAUmax.GT.10000.0)) THEN
           ! extreme taus
           Ncav = nY*3
        ELSE IF ((TAUmax.GT.2000.0)) THEN
           ! huge taus
           Ncav = nY*2
        ELSE IF ((TAUmax.GT.400.0)) THEN
           ! large taus
           Ncav = (nY*2)/3
        ELSE
           Ncav = (nY*2)/3
        END IF
        ! limit Ncav
        ncav = (nY*2)/3
        if (Ncav.gt.100) Ncav = 100
        IF (iX.NE.0) write(18,'(a20,i3)')' Ncav increased to:',Ncav
        call Pgrid(pstar,iPstar,nY,nP,nCav,nIns)
     end if
  elseif (slb) then
     if (initial.and.iterfbol.eq.1) then
        nY = 15
        DO iY = 1, nY
           IF(iY.LE.5) THEN
              Y(iY) = 0.5*(exp(0.5*iY)-exp(0.5))
           ELSE IF(iY.GT.5.and.iY.LE.12) THEN
              Y(iY) = Y(iY-1)
           ELSE
              Y(iY) = Y(2+nY-iY)
           END IF
        END DO
     end if
     tau(1) = 0.0d0
     do iY = 2, nY
        tau(iY) = tau(iY-1) + Y(iY)
     end do
     do iL = 1, nL
        do iY = 1, nY
           TAUslb(iL,iY) = tautot(iL)*tau(iY)/tau(nY)
        end do
     end do
  end if
  if (iX.ge.1) then
     if (slb) then
        write(18,'(a25,i4)')'  Tau grid generated, nY=', nY
     elseif (sph) then
        write(18,'(a24,i3)')' y grid generated, nY = ',nY
        write(18,'(a24,i3)')'                   nP = ',nP
        write(18,'(a24,i3)')'                 Nins = ',Nins
        write(18,'(a24,i3)')'                 Ncav = ',Ncav
     end if
  end if
  deallocate(tau)
  !---------------------------------------------------------------------
  return
end subroutine setGrids
!***********************************************************************

!***********************************************************************
subroutine Ygrid(pstar,iPstar,itereta,iterfbol,taumax,nY,nYprev,nP,&
     nCav,nIns)
!=======================================================================
! This subroutine generates the radial (Y) grid. Yout is the relative
! thickness, Yout=rout/r1. pstar is the impact parameter for star (0<=p<=1).
! First few points are prescribed arbitrarily.  This subroutine calls
! function ETA which evaluates the normalized density profile.
!                                       [ZI, Nov'95; MN,Sep'99, Deka'08]
!=======================================================================
  use common
  implicit none
  integer itereta, iterfbol, nY,nYprev, nP, nCav, nIns
  integer i,j, jy,itr,istop,iPstar, iter, iYdummy, iY, ir,irmax
  double precision pstar, dyR, dyT, eta, aux, EtaTemp(npY), &
       ee, Yloc,  Yold, Ynew, y1, delY1, ymid,      &
       rat, dY, etamin, etamax,deltau, delY,tauloc2,    &
       tausc, nh, fac, dyF, EtaRat, TAUloc, delTausc, facc, &
       taumax
  external eta
!--------------------------------------------------------------------
  !!** these are initialized as in old Dusty's sub Input. [MN]
  delTAUsc = 0.1 
  facc = 2.0
  EtaRat = 4.0
  !!**
  if(denstyp.ne.5) then
     ! max number iter. over improving ratio of two Eta's
     irmax = 20
     ! save old grid and values of Eta (important for denstyp = 5 or 6)
     IF (nY.GT.0.AND.(denstyp.eq.3)) THEN
        DO iY = 1, nY
           Yprev(iY) = Y(iY)
           EtaTemp(iY) = ETAdiscr(iY)
        END DO
        nYprev = nY
     END IF
     y1 = 1.0
     iter = 0
101  error = 0
     iter = iter + 1
     ! resolve inner boundary in tau space
     ! from requiring TAU(2) = TAUmax*ETA(1)*(Y(2)-1) =~ 1
     Y(1) = y1
     IF (TAUmax*ETA(y1,nY,nYprev,itereta).GT.2000.0) THEN
        delY1 = 2.0 / TAUmax / ETA(y1,nY,nYprev,itereta)
     ELSE
        delY1 = 1.0 / TAUmax / ETA(y1,nY,nYprev,itereta)
     END IF
     ! if very thin generate at least about 10-15 pts.
     IF (delY1.GT.0.01/ETA(y1,nY,nYprev,itereta)) delY1 = 0.01/ETA(y1,nY,nYprev,itereta)
     ! do not push it over the limit of spline stability
     IF (delY1.LT.0.00005) delY1 = 0.00005
     i = 1
     istop = 0
     DO WHILE (istop.NE.1)
        i = i + 1
        Y(i) = Y(i-1) + delY1
        IF (Y(i).GT.facc*Y(i-1)) THEN
           Y(i) = facc * Y(i-1)
        ELSE
           delY1 = 2. * delY1
        END IF
        Ymid = dsqrt(Y(i)*Y(i-1))
        TAUloc = TAUmax * (Y(i)-1.0) * ETA(Ymid,nY,nYprev,itereta)
        IF (Y(i).GE.1.01.OR.TAUloc.GT.10.0) istop = 1
        ! in case of shell thinner than 1.01 comment the above line and
        ! uncomment the next line, unless it is a case of RDW at high TauV.
        ! These require more points near the origin and 1.01 is a better limit.
        ! IF (Y(i).GE.1.0001.OR.TAUloc.GT.10.0) istop = 1
     END DO
     Ynew = Y(i)
     ! some rule of thumb estimates for factor nh
     IF (TAUmax.GT.10000.0) THEN
        ! extreme taus
        nh = 15.0
        Ncav = 80
     ELSE IF (TAUmax.GT.2000.0) THEN
        ! huge taus
        nh = 10.0
        Ncav = 40
     ELSE IF (TAUmax.GT.400.0) THEN
        ! large taus
        nh = 8.0
        Ncav = 20
     ELSE
        ! normal' taus
        nh = 5.0
        Ncav = 10
     END IF
     ! very small taus
     IF (TAUmax.LT.10.0) nh = nh / 2.0
     IF (TAUmax.LT.1.0) nh = nh / 2.0
     ! empirically: ~1/r needs more points for small y:
     IF ((pow-1.4)*(pow-0.6).LE.0.0) nh = nh * 1.5
     ! same for for steep density distributions:
     IF (denstyp.eq.3.OR.denstyp.eq.4) nh = nh * 1.5
     tausc = ETA(Y(i),nY,nYprev,itereta)+ETA(Y(i-1),nY,nYprev,itereta)
     tausc = (Y(i)-Y(i-1)) * 0.5 * tausc
     fac = dexp(-dlog(tausc) / nh)
     istop = 0
     ! for broken power-laws
     IF (Ntr.GE.1) THEN
        itr = 1
        DO j = 1, i
           IF (Y(j).GE.Ytr(itr)) itr = itr + 1
        END DO
     END IF
     ! generate the rest of Y grid points
     DO WHILE (istop.NE.1)
        i = i + 1
        Yold = Ynew
        ! find maximal increase in Y allowed by the ratio facc
        dyR = Yold * (facc-1.)
        ! find maximal increase in Y allowed by delTausc
        dyT = delTausc / ETA(Yold,nY,nYprev,itereta)
        ! find maximal increase in Y allowed by the ratio of tausc
        dyF = tausc*(fac-1.) / ETA(Yold,nY,nYprev,itereta)
        ! find new Y
        ! Ynew = Yold + MIN(dyR,dyT,dyF)
        dY = MIN(dyR,dyT,dyF)
        ! Check if the max ratio btw. two Eta values is less than EtaRat
        ! and insert additional y-pts. where necessary. This prevents sharp
        ! drops in Utot(npL,npY) in case of steep Eta's [MN'99].
        DO ir = 1 , irmax
           Ynew = Yold + dY
           rat = ETA(Yold,nY,nYprev,itereta)/ETA(Ynew,nY,nYprev,itereta)
           IF(rat.GE.1./EtaRat .AND. rat.LE.EtaRat) goto 10
           dY = 0.5*dY
        END DO
        CALL MSG(16)
10      continue
        Y(i) = Ynew
        ! make sure that all transition points are included in the grid
        IF (Ntr.GE.1) THEN
           IF (Y(i).GE.Ytr(itr)) THEN
              Y(i) = Ytr(itr)
              Ynew = Y(i)
              itr = itr + 1
           END IF
        END IF
        aux = ETA(Ynew,nY,nYprev,itereta)+ETA(Yold,nY,nYprev,itereta)
        aux = (Ynew-Yold) * 0.5 * aux
        tausc = tausc + aux
        ! finish when Yout is reached
        IF (Ynew.GE.Yout) istop = 1
     END DO
     Y(i) = Yout
     nY = i
     ! insert additional penultimate point to avoid spline oscillations
     Y(nY+1) = Yout
     Y(nY) = dsqrt(Y(nY)*Y(nY-1))
     nY = nY + 1
     ! check that outer edge is well resolved in tau space
     ! (important for flat ETAs)
     istop = 0
     DO WHILE (istop.NE.1)
        IF ((Yout-Y(nY-1))*TAUmax*ETA(Yout,nY,nYprev,itereta).GT.1.0) THEN
           Y(nY+1) = Yout
           Y(nY) = dsqrt(Y(nY)*Y(nY-1))
           nY = nY + 1
        ELSE
           istop = 1
        END IF
     END DO
     ! check dynamical range of Eta to avoid nonphysical results or code errors [
     Etamax = 0.
     Etamin = 1.e+20
     DO iY = 1, nY
        IF(ETA(Y(iY),nY,nYprev,itereta).lt.Etamin) Etamin = ETA(Y(iY),nY,nYprev,itereta)
        IF(ETA(Y(iY),nY,nYprev,itereta).gt.Etamax) Etamax = ETA(Y(iY),nY,nYprev,itereta)
        IF (ETA(Y(iY),nY,nYprev,itereta).LT.1.e-12) THEN
           IF (iX.GT.0) THEN
              write(18,*)'      Y          ETA  '
              DO jY = 1, nY
                 write(18,'(1p,2e12.3)') Y(jY),ETA(Y(jY),nY,nYprev,itereta)
              END DO
           END IF
           CALL MSG(17)
           error = 6
           goto 102
        END IF
     END DO
     IF ((Etamin/Etamax).LT.1.e-12) THEN
        CALL MSG(18)
        error = 6
        goto 102
     END IF
     ! check that the Y grid is not too large
     IF (nY.GT.npY) THEN
        delTAUsc = delTAUsc * 1.5
        IF (iX.NE.0) THEN
           IF (iter.EQ.1) call line(0,2,18)
           write(18,'(a)') ' ************   WARNING   *******************'
           write(18,'(a46,i3)')   &
                ' Initial delTAUsc resulted in too many points:',nY
           write(18,'(a,i3,a)')' You need to increase npY in userpar.inc'
           write(18,*)'Multiplying delTAUsc by 1.5 and trying again'
           write(18,'(a14,1p,e10.3)')' New delTAUsc:',delTAUsc
        END IF
        IF (iter.LT.5) THEN
           goto 101
        ELSE
           IF (iX.NE.0) THEN
              write(18,'(a)') ' **********  GIVING UP  ******************'
              call line(0,2,18)
           END IF
           error = 2
           goto 102
        END IF
     END IF
     DO iY = 1, nY
        Yloc = Y(iY)
        IF (iterETA.GT.1) THEN
           CALL LinInter(nY,nYprev,Yprev,EtaTemp,Yloc,iYdummy,ee)
           ETAdiscr(iY) = ee
        ELSE
           ETAdiscr(iY) = ETA(Yloc,nY,nYprev,itereta)
        END IF
      END DO
   else
      ! take the radial grid as in density profile file [MN,Apr'02]
      nY = nYEtaf
      do iY = 1, nY
         Y(iY) = yEtaf(iY)
      end do
   end if
   !--------------------------------------------------------------------
102 return
end subroutine ygrid
!***********************************************************************

!***********************************************************************
subroutine Pgrid(pstar,iPstar,nY,nP,nCav,nIns)
!=======================================================================
! after having the Y grid, generate the P grid (impact parameters)
!                                                              [Deka'08]
!=======================================================================
  use common
  implicit none
  !---parameter
  double precision pstar
  integer ipstar,nY,nP,nCav,nIns
  !---local
  integer i,ii,k,iP,iz,iw,nZ,j, Naux, istop, NinsLoc
  double precision delP, eta
  external eta
  !--------------------------------------------------------------------
  error  = 0
  ! Ncav = # rays in cavitiy
  ! Ncav = 30
  p(1) = 0.0d0
  Naux = 1
  iYfirst(Naux) = 1
  YPequal(Naux) = 1
  ! grid within the cavity (points concentrated towards p=1, sqrt!)
  nPcav = Ncav + Naux
  do iP = Naux+1, nPcav
     P(iP) = dble(iP-1)/dble(nPcav)
     iYfirst(iP) = 1
     YPequal(iP) = 1
  end do
  ! insert two rays if the stellar disk size is finite
  if (pstar.gt.0.0d0) then
     call insert(pstar,iP,iPstar)
     do i = nPcav+1, iP
        iYfirst(i) = 1
        YPequal(i) = 1
     end do
     nPcav = iP
  else
     iP  = nPcav
     iPstar = 0
  endif
  do i = 1, nY-1
     ! Nins = # of rays per radial step
     if((Y(i+1) - Y(i)).gt.20.0d0) then
        Nins = 20 !2
     elseif((Y(i+1) - Y(i)).gt.10.0d0.and.(Y(i+1) - Y(i)).le.20.d0) then
        Nins = 15 !2
     elseif ((Y(i+1) - Y(i)).gt.1.5d0.and.(Y(i+1) - Y(i)).le.10.d0) then
        Nins = 10 !2
     elseif ((Y(i+1) - Y(i)).gt.1.0d0.and.(Y(i+1) - Y(i)).le.1.5d0) then
        Nins = 3
     elseif ((Y(i+1) - Y(i)).gt.0.5d0.and.(Y(i+1) - Y(i)).le.1.0d0) then
        Nins = 5 !7
     else
        Nins = 7 !10
     end if
     ! in old dusty nins is allways 2
     Nins = 2
     k = Nins
     Plast(i) = iP + 1
     delP = (Y(i+1) - Y(i))/dble(k)
     do j = 1, k
        iP = iP + 1
        iYfirst(iP) = i
        if (j.eq.1) then
           YPequal(iP) = 1
        else
           YPequal(iP) = 0
        end if
        P(iP) = Y(i) + (j-1)*delP
     end do
  end do
  ! number of points for the p grid
  nP = iP + 1
  P(nP) = Y(nY)
  iYfirst(nP) = nY
  YPequal(nP) = 1
  Plast(nY) = nP
  Nins = nP - nPcav + 1
  ! check that the p grid is not too large
  if (nP.gt.npP) then
     error = 2
     call line(0,2,12)
     write(12,'(a)')' this model terminated because needed accuracy'
     write(12,'(a,i6)')' results in too many points in p grid:',nP
     write(12,'(a,i4,a)')'          (see userpar.inc, npP =',npP,')'
     call line(0,2,12)
  end if
  !-----------------------------------------------------------------------
  return
end subroutine pgrid

!***********************************************************************
!!$
!!$! **********************************************************************
!!$subroutine SLBtau(initial,iterfbol,tautot,Y,nL,nY,TAUslb)
!!$! ======================================================================
!!$! It generates TAUslb(npL,npY) grid with given spacing delT, calculated
!!$! in sub SLBy. In the current version it is exp near the slab boundaries
!!$! and equidistant in the middle.                  [Deka,'09, MN,'98,'09]
!!$! ======================================================================
!!$  implicit none
!!$  integer iterfbol, nY, iY, nL, iL
!!$  integer npY, npP, npX, npL, npG, npR
!!$  include '../userpar.inc'
!!$  double precision TAUslb(npL,npY),tautot(npL),tau(npY),Y(npY)
!!$  logical initial
!!$! ----------------------------------------------------------------------
!!$  if (initial.and.iterfbol.eq.1) then
!!$     nY = 15
!!$     DO iY = 1, nY
!!$        IF(iY.LE.5) THEN
!!$           Y(iY) = 0.5*(exp(0.5*iY)-exp(0.5))
!!$        ELSE IF(iY.GT.5.and.iY.LE.12) THEN
!!$           Y(iY) = Y(iY-1)
!!$        ELSE
!!$           Y(iY) = Y(2+nY-iY)
!!$        END IF
!!$     END DO
!!$  end if
!!$  tau(1) = 0.0d0
!!$  do iY = 2, nY
!!$     tau(iY) = tau(iY-1) + Y(iY)
!!$  end do
!!$  do iL = 1, nL
!!$     do iY = 1, nY
!!$        TAUslb(iL,iY) = tautot(iL)*tau(iY)/tau(nY)
!!$     end do
!!$  end do
!!$!----------------------------------------------------------------------
!!$  return
!!$end subroutine SLBtau
!!$!***********************************************************************
!!$
!!$
!***********************************************************************
 double precision function ETA(Yy,nY,nYprev,itereta)
!=======================================================================
! This subroutine evaluates the normalized density profile. denstype is
! a type of density law: 1 and 2 - power law, 3 for exponential density law,
! 4 for radiatively driven winds (the gray-body approximation), 5,6 - RDW and
! 7 - d.d.from a file. pow is parameter describing the choosen density law:
! power for 1 and 2, v1/v8 for RDW (the ratio of expansion velocities at the
! inner and outer radii), sigma for 3 [i.e. rho = dexp(-(y/sigma)**2)].
! Yout is the relative thickness, Yout=rout/r1. Y is the radial position.
!                                                          [ZI'95; ZI'99]
!=======================================================================
   use common
   implicit none
   integer i, istop, iYdummy,nY,nYprev,itereta
   double precision Yy, c, intaux, powaux, prod, eps_loc, facteta
!-----------------------------------------------------------------------
   ! this is an adjustable value regulating the initial eta approximation
   ! for dynamics
   facteta = 0.5d0
   if (Yy.gt.Yout) Yy = Yout
   if (itereta.gt.1) then
      ! if this is iteration over eta (but not the first one) in the case
      ! of dynamical calculation for radiatively driven winds calculate
      ! eta by linear interpolation of eta from the previous iteration
      call lininter(nY,nYprev,Yprev,etadiscr,Yy,iYdummy,eta)
      ! otherwise use prescribed formulae for different cases
   else
      ! smooth power-law
      if (denstyp.eq.1 .and. Ntr.eq.0) then
         ! find normalization constant
         if (pow.ne.1.0d0) then
            c = (1.0d0 - Yout**(1.0d0-pow)) / (pow - 1.0d0)
         else
            c = dlog(Yout)
         endif
         if (Ntr.ge.1) then
            do i = 1, Ntr
               powaux = pow - ptr(i)
               if (powaux.ne.1.0d0) then
                  intaux =(1.0d0-Yout**(1.0d0-powaux))/(powaux-1.0d0)
               else
                  intaux = dlog(Yout)
               end if
               c = c + intaux / Ytr(i)**ptr(i)
            end do
         endif
         c = 1.0d0 / c
         ! calculate density
         if (Yy.ge.(1.0d0-1.0d-08)) then
            eta = c / Yy**pow
            if (Ntr.ge.1) then
               do i = 1, Ntr
                  eta = eta + c * Yy**(ptr(i)-pow) / Ytr(i)**ptr(i)
               end do
            end if
         else
            eta = 0.0d0
         endif
      end if
      ! broken power-law
      if (denstyp.eq.1.and.Ntr.gt.0) then
         Ytr(Ntr+1) = Yout
         ! find normalization constants
         if (pow.ne.1.0) then
            c = (1.0d0 - Ytr(1)**(1.0d0-pow)) / (pow - 1.0d0)
         else
            c = dlog(Ytr(1))
         endif
         if (Ntr.ge.1) then
            do i = 1, Ntr
               call doproduct(10,Ytr,ptr,pow,i,prod)
               if (ptr(i).ne.1.0d0) then
                  intaux=Ytr(i)**(1.0d0-ptr(i))-Ytr(i+1)**(1.0d0-ptr(i))
                  intaux=prod * intaux / (ptr(i) - 1.0d0)
               else
                  intaux = prod * dlog(Ytr(i+1)/Ytr(i))
               end if
               c = c + intaux
            end do
         endif
         c = 1.0d0 / c
         ! calculate density
         if (Yy.ge.1.0d0-1.0d-8) then
            if (Yy.le.Ytr(1)) then
               eta = c / Yy**pow
            else
               istop = 0
               i = 0
               do while (istop.ne.1)
                  i = i + 1
                  if (Yy.le.Ytr(i+1)) istop = 1
               end do
               call doproduct(10,Ytr,ptr,pow,i,prod)
               eta = c * prod / Yy**ptr(i)
            end if
         else
            eta = 0.0d0
         endif
      end if
      ! exponential law
      if (denstyp.eq.2) then
         if (Yy.ge.1.0d0-1.0d-8) then
            eta = (Yout-1.) * (1.0d0-dexp(-pow)) / pow
            eta = dexp(-pow*(Yy-1.0d0)/(Yout-1.0d0)) / eta
         else
            eta = 0.0d0
         end if
      end if
      ! radiatively driven winds (the gray-body approximation)
      if (denstyp.eq.4) then
         eps_loc = pow
         if (Yy.ge.1.0d0-1.0d-8) then
            eta = (1.0d0+eps_loc)/2.0d0/Yy/Yy/sqrt(1.0d0-(1.0d0-eps_loc*eps_loc)/Yy)
         else
            eta = 0.0d0
         endif
      end if
      ! radiatively driven winds (i.e. denstyp=3 or 6)
      if (denstyp.eq.3.or.denstyp.eq.6) then
         ! if this is the first iteration use analytic approximation
         if (itereta.lt.2) then
            ! for 3 eps_loc is pow, for 6 assume eps_loc=0.1
            if (denstyp.eq.3) then
               eps_loc = pow
            else
               eps_loc = 0.1d0
            end if
            if (Yy.ge.(1.0d0-1.0d-8)) then
               eta = (1.0d0+eps_loc)/2.0d0/Yy/Yy/sqrt(1.0d0-(1.0d0-eps_loc*eps_loc)/Yy)
               ! empirical improvement for the initial approximation
               ! good only for large optical depths, but small ones
               ! are fast anyway (zi, may99)
               if (Yy.le.2.0d0) then
                  eta = eta / (1.0d0 + facteta / Yy**10.0d0)
               end if
            else
               eta = 0.0d0
            endif
            ! or interpolate from the previous solution
         else
            call lininter(nY,nYprev,Yprev,etadiscr,Yy,iYdummy,eta)
         end if
      end if
      ! user specified function for eta
      if (denstyp.eq.5) then
         if (Yy.lt.yetaf(nYetaf)) then
            call lininter(nY,nYetaf,yetaf,etaf,Yy,iYdummy,eta)
         else
            eta = etaf(nYetaf)
         end if
      end if
      ! done
   end if
!-----------------------------------------------------------------------
   return
end function eta
!***********************************************************************
!!$
!***********************************************************************
subroutine Insert(pstar,iP,iPstar)
!=======================================================================
! This subroutine inserts two rays if the stellar disk is finite. The
! first ray, corresponding to the disk edge, has index iPstar, and the
! following ray with 1.0001 greater impact parameter has index iPstar+1.
! The only exception is if pstar>0.9999 when only a single ray is
! inserted.                                            [Z.I., Feb. 1996]
! =======================================================================
  use common
  implicit none
  integer iP, iPstar
  double precision pstar
! -----------------------------------------------------------------------
  if (pstar.ge.0.999999d0) then
     pstar = 0.999999d0
     iP = nPcav + 1
     iPstar = nPcav + 1
     P(iP) = pstar
  else
     if (pstar.ge.P(nPcav)) then
        if (pstar.eq.P(nPcav)) pstar = 1.00001d0*P(nPcav)
        P(nPcav+1) = pstar
        iPstar = nPcav + 1
        P(nPcav+2) = 1.00001d0*pstar
        iP = nPcav + 2
     else
        iPstar = 0
        iP = 0
        do while (iPstar.eq.0)
           iP = iP + 1
           if (P(iP).gt.pstar.and.iPstar.eq.0) iPstar = iP
        end do
        do iP = 1, nPcav-iPstar+1
           P(nPcav+3-iP) = P(nPcav+1-iP)
        end do
        if (pstar.eq.P(iPstar-1)) pstar = 1.00001d0*P(iPstar-1)
        P(iPstar) = pstar
        P(iPstar+1) = 1.00001d0*pstar
        iP = nPcav + 1
     end if
  end if
! -----------------------------------------------------------------------
  return
end subroutine Insert
!***********************************************************************

!********************************************************************
subroutine Rad_Transf(initial,nY,nYprev,nP,itereta,pstar,y_incr,us,fs,emiss, &
     iterfbol,initTemp,T4_ext)
!======================================================================
  use omp_lib
  use common
  implicit none
  INTERFACE
     subroutine find_Text(nY,T4_ext)
       integer nY
       double precision, allocatable :: T4_ext(:)
     end subroutine find_Text
     subroutine OccltMSG(us)
       double precision,allocatable :: us(:,:)
     end subroutine OccltMSG
     subroutine Find_Tran(pstar,nY,nP,T4_ext,us,fs)
       integer nY, nP
       double precision :: pstar
       double precision, allocatable :: T4_ext(:)
       double precision, allocatable :: fs(:,:),us(:,:)
     end subroutine Find_Tran
     subroutine Emission(nY,T4_ext,emiss,emiss_total)
       integer nY
       double precision, allocatable :: T4_ext(:),emiss(:,:,:),emiss_total(:,:)
     end subroutine Emission
     subroutine find_diffuse(nY,nP,initial,moment,iter,iterfbol,T4_ext,us,emiss)
       integer nY,nP,iter,iterfbol,moment
       logical initial
       double precision, allocatable :: T4_ext(:)
       double precision, allocatable :: us(:,:)
       double precision, allocatable :: emiss(:,:,:)
     end subroutine find_diffuse
     subroutine init_temp(nY,T4_ext,us)
       integer nY
       double precision,allocatable :: us(:,:),T4_ext(:)
     end subroutine init_temp
     subroutine find_temp(nY,T4_ext)
       integer :: nY
       double precision, allocatable :: T4_ext(:)
     end subroutine find_temp
     subroutine SPH_DIFF(flag,moment,nY,nP,initial,iter,iterfbol,T4_ext,emiss,us,vec2)
       integer nY,nP,iter,iterfbol,flag,moment
       double precision, allocatable :: T4_ext(:),emiss(:,:,:),us(:,:),vec2(:,:)
       logical initial
     end subroutine SPH_DIFF
  END INTERFACE
  logical, intent(in) :: initial,initTemp
  integer, intent(in) :: y_incr,iterfbol
  integer :: nY,nP,nYprev,itereta
  double precision pstar
  double precision,allocatable :: us(:,:), fs(:,:),ubol_old(:)
  double precision,allocatable :: T4_ext(:)
  double precision,allocatable :: emiss(:,:,:),emiss_total(:,:)
  !---- local variable
  double precision,allocatable :: utot_old2(:,:)
  integer :: itlim,conv,iter,iG,iL,iY,iY1, moment, thread_id,i,istop,iOut
  double precision aux1,maxerrT(max_threads),maxerrU(max_threads),m,n,JL,JR,xx
  external eta

  allocate(emiss_total(nL,npY))
  allocate(ubol_old(npY))
  allocate(utot_old2(nL,npY))

  !------------------------------------------------------------------------
  error = 0
  if(sph) then
     ! generate spline coefficients for ETA as in old Dusty [MN'Aug,10]
     CALL setupETA(nY,nYprev,itereta)
     ! evaluate ETAzp (carried in common)
     CALL getETAzp(nY,nP)
  end if
!!$  ! the tau-profile at the fiducious lambda (needed in prout)
!!$  if (slb) then
!!$     do iY = 1, nY
!!$        tr(iY) = TAUslb(iLfid,iY)/TAUslb(iLfid,nY)
!!$     end do
!!$  elseif(sph) then
!!$     do iY = 1, nY
!!$        tr(iY) = ETAzp(1,iY)/ETAzp(1,nY)
!!$     end do
!!$  end if
  ! generate albedo through the envelope
  call getOmega(nY)
  ! generate stellar spectrum
  call Find_Tran(pstar,nY,nP,T4_ext,us,fs)
  if(iVerb.eq.2) write(*,*)' Done with transmitted radiation.'
  ! issue a message in fname.out about the condition for neglecting
  ! occultation only if T1 is given in input:
  if(typentry(1).eq.5.and.sph) then
     if(iterfbol.eq.1.and.itereta.eq.1.and.right.eq.0) call OccltMSG(us)
  end if
  ! finish when file with the stellar spectrum is not available
  if (error.eq.3) goto 999
  ! in the case of first (lowest) optical depth,
  ! us is the intial approximation for utot(iL,iY) for the first iteration over Td
  ! Find initial approximation of Td for the case of first iteration over Fbol or flux error to large.
  !if ((initial.and.(iterfbol.eq.1)).or.(initTemp.and.(iterfbol.gt.2))) then
  if ((iterfbol.eq.1).or.(initTemp.and.(iterfbol.gt.2))) then
     call Init_Temp(nY,T4_ext,us)
     if(iVerb.eq.2) write(*,*)' Done with initial dust temperature.'
  end if
  do iY = 1,nY
     do iG = 1,nG
        Td_old(iG,iY) = Td(iG,iY)
     end do
  end do
  itlim = 2000
  conv = 0
  iter = 0
  !=== iterations over dust temperature =========
  do while (conv.eq.0.and.iter.le.itlim)
     iter = iter + 1
     !print*,iter
     ! find T_external for the new y-grid if T(1) given in input
     if (typentry(1).eq.5) call find_Text(nY,T4_ext)
     ! find emission term
     call Emission(nY,T4_ext,emiss,emiss_total)
     ! moment = 1 is for finding total energy density only
     moment = 1
     call Find_Diffuse(nY,nP,initial,moment,iter,iterfbol,T4_ext,us,emiss)
     call Find_Temp(nY,T4_ext)
     ! assign previus Td to Td_old
     do thread_id =1,max_threads
        maxerrU(thread_id) = 0.0d0  !just for info, iterations are over Td [MN]
        maxerrT(thread_id) = 0.0d0
     end do
     !$OMP PARALLEL DO PRIVATE(thread_id,iL,aux1,iG,m,n)
     do iY = 1,nY
        thread_id = omp_get_thread_num()+1
        do iL=1,nL
           if (abs(utot(iL,iY)).gt.dynrange) then 
              aux1 = abs(utot_old(iL,iY)-utot(iL,iY))/(abs(utot_old(iL,iY))+abs(utot(iL,iY)))
           else 
              aux1 = 0.0D0
           end if
           utot_old2(iL,iY) = utot_old(iL,iY)
           utot_old(iL,iY) = utot(iL,iY)
           if ((iter.gt.5000).and.(mod(iter,5).eq.0)) then 
              m = (utot_old(iL,iY)-utot_old2(iL,iY))
              if (m.gt.dynrange) then 
                 n = utot_old2(iL,iY)
                 utot_old(iL,iY) = m*(iter/3)+n
                 utot(iL,iY) = m*(1+iter/3)+n
                 !cut off everything smaller than dynrange**2 including negative values
                 if (utot_old(iL,iY).lt.dynrange*dynrange) utot_old(iL,iY) = 0.0D0
                 if (utot(iL,iY).lt.dynrange*dynrange) utot(iL,iY) = 0.0D0 
                 if (abs(utot(iL,iY)).gt.dynrange) then 
                    aux1 = abs(utot_old(iL,iY)-utot(iL,iY))/(abs(utot_old(iL,iY))+abs(utot(iL,iY)))
                 else
                    aux1 = 0.0D0
                 end if
              end if
           end if
           if (maxerrU(thread_id).lt.aux1) maxerrU(thread_id)=aux1
        end do
        do iG = 1,nG
           aux1 = dabs(Td_old(iG,iY) - Td(iG,iY))/Td(iG,iY)
           if (aux1.gt.maxerrT(thread_id)) maxerrT(thread_id) = aux1
           Td_old(iG,iY) = Td(iG,iY)
        end do
     end do
     !$OMP END PARALLEL DO
!!$     if ((iter.gt.5).and.(mod(iter,5).eq.0)) then 
!!$        call Find_Temp(nY,T4_ext)
!!$     end if
     !print*,maxval(maxerrU),Td(1,1),Td(1,nY)
      if(iVerb.eq.2) write(*,fmt='(a1)',advance='no') '.'
     if ((maxval(maxerrT).le.accTemp).and.(maxval(maxerrU).lt.(accFlux*9.e-1))) conv = 1
     if (iter.eq.itlim) print'(A,I6,A)','  !!! Reached iteration limit of ',itlim,' !!!!'
  enddo
   if(iVerb.eq.2) write(*,*) ' '
  !=== the end of iterations over Td ===
  if(iVerb.eq.2) then 
     write(*,'(A,I3,A)') '  Done with finding dust temperature after ',iter,' iterations'
     write(*,'(A,1PE9.3,A,1PE9.3)') '    errT: ',maxval(maxerrT),' errU: ',maxval(maxerrU)
  end if
  ! find T_external for the converged dust temperature
  if (typentry(1).eq.5) call find_Text(nY,T4_ext)
  ! find Jext, needed in PrOut [MN]
  do iY = 1, nY
     Jext(iY) = sigma/pi * T4_ext(iY)
  end do
  ! calculate the emission term using the converged Td
  call Emission(nY,T4_ext,emiss,emiss_total)
  ! calculate total energy density and diffuse flux using the converged Td
  moment = 2
  call Find_Diffuse(nY,nP,initial,moment,iter,iterfbol,T4_ext,us,emiss)
  if(iVerb.eq.2) write(*,*) ' Done with finding energy density and diffuse flux.'
  !-----------------------------------------------------------------
  ! Find the energy density profiles if required.
  ! They are normalized in PrOut [MN,11]
  IF (iJ.gt.0) THEN
     if (sph_matrix) then 
        print*,'  J output not available for matrix method!!!'
     else
        call SPH_diff(1,0,nY,nP,initial,iter,iterfbol,T4_ext,emiss,us,Ude)
        call SPH_diff(2,0,nY,nP,initial,iter,iterfbol,T4_ext,emiss,us,Uds)
     end if
     ! interpolate J-output(iOut) to Y(iOut)
      DO iOut = 1, nJOut
         ! bracket the needed wavelength
        istop = 0
        i = 0
        DO WHILE (istop.EQ.0)
          i = i + 1
          IF (Y(i).GT.YJOut(iOut)) istop = 1
          IF (i.EQ.nJout) istop = 1
        END DO
        ! interpolate intensity
        xx = (YJOut(iOut)-Y(i-1))/(Y(i)-Y(i-1))
        DO iL = 1, nL
          JL = Ude(iL,i-1) + Uds(iL,i-1)
          JR = Ude(iL,i) + Uds(iL,i)
          JOut(iL,iOut) = JL + xx*(JR - JL)
        END DO
      END DO
   END IF
999 deallocate(emiss_total)
  deallocate(ubol_old)
  deallocate(utot_old2)
  return
end subroutine Rad_Transf
!***********************************************************************

!!$!********************************************************************
!!$subroutine Rad_Transf(initial,nY,nYprev,nP,itereta,pstar,y_incr,us,fs,emiss, &
!!$     iterfbol,initTemp,T4_ext)
!!$!======================================================================
!!$  use omp_lib
!!$  use common
!!$  implicit none
!!$  INTERFACE
!!$     subroutine find_Text(nY,T4_ext)
!!$       integer nY
!!$       double precision, allocatable :: T4_ext(:)
!!$     end subroutine find_Text
!!$     subroutine OccltMSG(us)
!!$       double precision,allocatable :: us(:,:)
!!$     end subroutine OccltMSG
!!$     subroutine Find_Tran(pstar,nY,nP,T4_ext,us,fs)
!!$       integer nY, nP
!!$       double precision :: pstar
!!$       double precision, allocatable :: T4_ext(:)
!!$       double precision, allocatable :: fs(:,:),us(:,:)
!!$     end subroutine Find_Tran
!!$     subroutine Emission(nY,T4_ext,emiss,emiss_total)
!!$       integer nY
!!$       double precision, allocatable :: T4_ext(:),emiss(:,:,:),emiss_total(:,:)
!!$     end subroutine Emission
!!$     subroutine find_diffuse(nY,nP,initial,moment,iter,iterfbol,T4_ext,us,emiss)
!!$       integer nY,nP,iter,iterfbol,moment
!!$       logical initial
!!$       double precision, allocatable :: T4_ext(:)
!!$       double precision, allocatable :: us(:,:)
!!$       double precision, allocatable :: emiss(:,:,:)
!!$     end subroutine find_diffuse
!!$     subroutine init_temp(nY,T4_ext,us)
!!$       integer nY
!!$       double precision,allocatable :: us(:,:),T4_ext(:)
!!$     end subroutine init_temp
!!$     subroutine find_temp(nY,T4_ext)
!!$       integer :: nY
!!$       double precision, allocatable :: T4_ext(:)
!!$     end subroutine find_temp
!!$     subroutine SPH_DIFF(flag,moment,nY,nP,initial,iter,iterfbol,T4_ext,emiss,us,vec2)
!!$       integer nY,nP,iter,iterfbol,flag,moment
!!$       double precision, allocatable :: T4_ext(:),emiss(:,:,:),us(:,:),vec2(:,:)
!!$       logical initial
!!$     end subroutine SPH_DIFF
!!$  END INTERFACE
!!$  logical, intent(in) :: initial,initTemp
!!$  integer, intent(in) :: y_incr,iterfbol
!!$  integer :: nY,nP,nYprev,itereta
!!$  double precision pstar
!!$  double precision,allocatable :: us(:,:), fs(:,:),ubol_old(:)
!!$  double precision,allocatable :: T4_ext(:)
!!$  double precision,allocatable :: emiss(:,:,:),emiss_total(:,:)
!!$  !---- local variable
!!$  integer :: itlim,conv,iter,iG,iL,iY,iY1, moment, thread_id,i,istop,iOut
!!$  double precision aux1,maxerrT,maxerrU,m,n,JL,JR,xx
!!$  external eta
!!$
!!$  allocate(emiss_total(nL,npY))
!!$  allocate(ubol_old(npY))
!!$
!!$  !------------------------------------------------------------------------
!!$  error = 0
!!$  if(sph) then
!!$     ! generate spline coefficients for ETA as in old Dusty [MN'Aug,10]
!!$     CALL setupETA(nY,nYprev,itereta)
!!$     ! evaluate ETAzp (carried in common)
!!$     CALL getETAzp(nY,nP)
!!$  end if
!!$  ! generate albedo through the envelope
!!$  call getOmega(nY)
!!$  ! generate stellar spectrum
!!$  call Find_Tran(pstar,nY,nP,T4_ext,us,fs)
!!$  if(iVerb.eq.2) write(*,*)' Done with transmitted radiation.'
!!$  ! issue a message in fname.out about the condition for neglecting
!!$  ! occultation only if T1 is given in input:
!!$  if(typentry(1).eq.5.and.sph) then
!!$     if(iterfbol.eq.1.and.itereta.eq.1.and.right.eq.0) call OccltMSG(us)
!!$  end if
!!$  ! finish when file with the stellar spectrum is not available
!!$  if (error.eq.3) goto 999
!!$  ! in the case of first (lowest) optical depth,
!!$  ! us is the intial approximation for utot(iL,iY) for the first iteration over Td
!!$  ! Find initial approximation of Td for the case of first iteration over Fbol or flux error to large.
!!$  !if ((initial.and.(iterfbol.eq.1)).or.(initTemp.and.(iterfbol.gt.2))) then
!!$  if ((iterfbol.eq.1).or.(initTemp.and.(iterfbol.gt.2))) then
!!$     call Init_Temp(nY,T4_ext,us)
!!$     if(iVerb.eq.2) write(*,*)' Done with initial dust temperature.'
!!$  end if
!!$  do iY = 1,nY
!!$     do iG = 1,nG
!!$        Td_old(iG,iY) = Td(iG,iY)
!!$     end do
!!$  end do
!!$  itlim = 2000
!!$  conv = 0
!!$  iter = 0
!!$  !=== iterations over dust temperature =========
!!$  do while (conv.eq.0.and.iter.le.itlim)
!!$     iter = iter + 1
!!$     !print*,iter
!!$     ! find T_external for the new y-grid if T(1) given in input
!!$     if (typentry(1).eq.5) call find_Text(nY,T4_ext)
!!$     ! find emission term
!!$     call Emission(nY,T4_ext,emiss,emiss_total)
!!$     ! moment = 1 is for finding total energy density only
!!$     moment = 1
!!$     call Find_Diffuse(nY,nP,initial,moment,iter,iterfbol,T4_ext,us,emiss)
!!$     call Find_Temp(nY,T4_ext)
!!$     ! assign previus Td to Td_old
!!$     maxerrU = 0.0d0  
!!$     maxerrT = 0.0d0
!!$     do iY = 1,nY
!!$        do iL=1,nL
!!$           if (abs(utot(iL,iY)).gt.dynrange*dynrange) then 
!!$              aux1 = dabs((utot_old(iL,iY)-utot(iL,iY))/utot(iL,iY))
!!$           else 
!!$              aux1 = 0.0D0
!!$           end if
!!$           if (maxerrU.lt.aux1) then 
!!$              maxerrU=aux1
!!$              !print*,iL,iY,utot(iL,iY),utot_old(iL,iY),aux1
!!$           end if
!!$           utot_old(iL,iY) = utot(iL,iY)
!!$        end do
!!$        do iG = 1,nG
!!$           aux1 = dabs(Td_old(iG,iY) - Td(iG,iY))/Td(iG,iY)
!!$           if (aux1.gt.maxerrT) maxerrT = aux1
!!$           Td_old(iG,iY) = Td(iG,iY)
!!$        end do
!!$     end do
!!$     if(iVerb.eq.2) write(*,fmt='(a1)',advance='no') '.'
!!$     if ((maxerrT.le.accTemp).and.(maxerrU.lt.(accFlux*9.e-1))) conv = 1
!!$     if (iter.eq.itlim) print'(A,I6,A)','  !!! Reached iteration limit of ',itlim,' !!!!'
!!$  enddo
!!$   if(iVerb.eq.2) write(*,*) ' '
!!$  !=== the end of iterations over Td ===
!!$  if(iVerb.eq.2) then 
!!$     write(*,'(A,I3,A)') '  Done with finding dust temperature after ',iter,' iterations'
!!$     write(*,'(A,1PE9.3,A,1PE9.3)') '    errT: ',maxerrT,' errU: ',maxerrU
!!$  end if
!!$  ! find T_external for the converged dust temperature
!!$  if (typentry(1).eq.5) call find_Text(nY,T4_ext)
!!$  ! find Jext, needed in PrOut [MN]
!!$  do iY = 1, nY
!!$     Jext(iY) = sigma/pi * T4_ext(iY)
!!$  end do
!!$  ! calculate the emission term using the converged Td
!!$  call Emission(nY,T4_ext,emiss,emiss_total)
!!$  ! calculate total energy density and diffuse flux using the converged Td
!!$  moment = 2
!!$  call Find_Diffuse(nY,nP,initial,moment,iter,iterfbol,T4_ext,us,emiss)
!!$  if(iVerb.eq.2) write(*,*) ' Done with finding energy density and diffuse flux.'
!!$  !-----------------------------------------------------------------
!!$  ! Find the energy density profiles if required.
!!$  ! They are normalized in PrOut [MN,11]
!!$  IF (iJ.gt.0) THEN
!!$     if (sph_matrix) then 
!!$        print*,'  J output not available for matrix method!!!'
!!$     else
!!$        call SPH_diff(1,0,nY,nP,initial,iter,iterfbol,T4_ext,emiss,us,Ude)
!!$        call SPH_diff(2,0,nY,nP,initial,iter,iterfbol,T4_ext,emiss,us,Uds)
!!$     end if
!!$     ! interpolate J-output(iOut) to Y(iOut)
!!$      DO iOut = 1, nJOut
!!$         ! bracket the needed wavelength
!!$        istop = 0
!!$        i = 0
!!$        DO WHILE (istop.EQ.0)
!!$          i = i + 1
!!$          IF (Y(i).GT.YJOut(iOut)) istop = 1
!!$          IF (i.EQ.nJout) istop = 1
!!$        END DO
!!$        ! interpolate intensity
!!$        xx = (YJOut(iOut)-Y(i-1))/(Y(i)-Y(i-1))
!!$        DO iL = 1, nL
!!$          JL = Ude(iL,i-1) + Uds(iL,i-1)
!!$          JR = Ude(iL,i) + Uds(iL,i)
!!$          JOut(iL,iOut) = JL + xx*(JR - JL)
!!$        END DO
!!$      END DO
!!$   END IF
!!$999 deallocate(emiss_total)
!!$  deallocate(ubol_old)
!!$  return
!!$end subroutine Rad_Transf
!!$!***********************************************************************

!***********************************************************************
subroutine Find_Tran(pstar,nY,nP,T4_ext,us,fs)
!=======================================================================
! This subroutine generates the transmitted energy density and  flux for
! the slab and sphere. is=1 is for the left src, is=2 is the right src
!                                               [MN, Feb.'99, Deka, 2008]
!=======================================================================
  use common
  implicit none
  INTERFACE
     subroutine Bolom(q,qbol,nY)
       integer nY
       double precision, allocatable :: q(:,:), qbol(:)
     end subroutine Bolom
     subroutine SPH_ext_illum(m0,m1,m1p,m1m,nY,nP)
       integer nY, nP
       double precision,allocatable ::  m0(:,:), m1(:,:),m1p(:,:),m1m(:,:)
     end subroutine SPH_ext_illum
  END INTERFACE
  !--- parameter
  integer nY, nP
  double precision :: pstar
  double precision, allocatable :: T4_ext(:)
  double precision, allocatable :: fs(:,:),us(:,:)

  integer i, iY, iL, iz, is, nZ, nn
  double precision  arg, dyn2, eint2, eint3, x, res, &
       denom, expow,  tauaux(npL,npY), zeta, result1, eta
  double precision,allocatable ::  usL(:,:),usR(:,:),&
       m0(:,:), m1(:,:),m1p(:,:),m1m(:,:)
  external eta
!-----------------------------------------------------------------------

  allocate(usL(nL,nY))
  allocate(usR(nL,nY))
  allocate(m0(nL,nY))
  allocate(m1(nL,nY))
  allocate(m1p(nL,nY))
  allocate(m1m(nL,nY))

  dyn2 = 1.0d-30
  ! define T_external for typEntry(1).ne.5
  if (typEntry(1).ne.5) then
     ! for slab
     if(slb) then
        do iY = 1, nY
           if (left.eq.1.and.right.eq.0) then
              if (mu1.eq.-1.0d0) then
                 T4_ext(iY) = pi*Ji/(2.0d0*sigma)
              else
                 T4_ext(iY) = pi*Ji/sigma
              end if
           elseif (left.eq.0.and.right.eq.1) then
              if (mu2.eq.-1.0d0) then
                 T4_ext(iY) = pi*ksi*Ji/(2.0d0*sigma)
              else
                 T4_ext(iY) = pi*ksi*Ji/(sigma)
              end if
           elseif (left.eq.1.and.right.eq.1) then
              if (mu1.eq.-1.0d0.and.mu2.eq.-1.0d0) then
                 T4_ext(iY) = pi*(Ji + ksi*Ji)/(2.0d0*sigma)
              else
                 T4_ext(iY) = pi*(Ji + ksi*Ji)/sigma
              end if
           end if
        end do
     ! for sphere
     elseif(sph) then
        do iY = 1, nY
           if (left.eq.1.and.right.eq.0) then
              T4_ext(iY) =  (pi/sigma)*Ji/(Y(iY)**2.0d0)
           elseif (left.eq.0.and.right.eq.1) then
              T4_ext(iY) =  pi*Jo/sigma
           elseif (left.eq.1.and.right.eq.1) then
              T4_ext(iY) =  (pi/sigma)*(Ji/(Y(iY)**2.0d0) + Jo)
           end if
        end do
     end if
  else
  end if
  do iY = 1, nY
     ! loop over wavelengths
     do iL = 1, nL
        ! find stellar part of flux and en.density
        ! for slab
        if(slb) then
           ! for the source on the left
           ! for isotropic illumination (from a half-plane)
           if(mu1.eq.-1.0d0) then
              usL(iL,iY) = shpL(iL)*eint2(TAUslb(iL,iY))
              fsL(iL,iY) = 4.0d0*pi*shpL(iL)*eint3(TAUslb(iL,iY))
           else
              ! For directional illumination
              x = TAUslb(iL,iY)/mu1
              if(x.ge.50.0d0) then
                 fsL(iL,iY) = 0.0d0
                 usL(iL,iY) = 0.0d0
              else
                 usL(iL,iY) = shpL(iL)*dexp(-x)
                 fsL(iL,iY) = 4.0d0*pi*mu1*shpL(iL)*dexp(-x)
              end if
           end if
           ! For the second source on the right
           if(right.gt.0) then
              if(mu2.eq.-1.0d0) then
                 ! If isotropic illumination from the right
                 arg = TAUslb(iL,nY)-TAUslb(iL,iY)
                 usR(iL,iY) = ksi*shpR(iL)*eint2(arg)
                 fsR(iL,iY) = 4.0d0*pi*ksi*shpR(iL)*eint3(arg)
              else
                 !  For directional illumination from the right
                 x = (TAUslb(iL,nY)-TAUslb(iL,iY))/mu2
                 if(x.ge.50.0d0) then
                    fsR(iL,iY) = 0.0d0
                    usR(iL,iY) = 0.0d0
                 else
                    usR(iL,iY) = ksi*shpR(iL)*dexp(-x)
                    fsR(iL,iY) = 4.0d0*pi*mu2*ksi*shpR(iL)*dexp(-x)
                 end if
              end if
              ! end if for source on the right
           end if
           ! for sphere
        elseif(sph) then
           ! for the central (left) source
           if (left.eq.1) then
              ! effect of the star's finite size
              if (pstar.gt.0.0d0) then
                 zeta = 2.0d0*(1.0d0-sqrt(1.0d0-(pstar/Y(iY))**2.0d0))* &
                      (Y(iY)/pstar)**2.0d0
              else
                 zeta = 1.0d0
              end if
              expow = ETAzp(1,iY)*TAUtot(iL)
              if(expow.lt.50.0d0) then
                 ! en. denisity
                 usL(iL,iY) = shpL(iL)*zeta*exp(-expow)
                 ! flux
                 fsL(iL,iY) = 4.0d0*pi*shpL(iL)*exp(-expow)
              else
                 usL(iL,iY) = 0.0d0
                 fsL(iL,iY) = 0.0d0
              end if
           end if
        end if
     end do
  end do
  if(sph.and.right.eq.1) then
     call SPH_ext_illum(m0,m1,m1p,m1m,nY,nP)
     do iY = 1, nY
        do iL = 1, nL
           ! find the en.density of the ext.radiation
           usR(iL,iY) = shpR(iL)*m0(iL,iY)/2.0d0
           ! and the external flux
           fsR(iL,iY) = 2.0d0*pi*shpR(iL)*m1(iL,iY)
        end do
     end do
  end if
  !scaling with energy density
  do iY = 1, nY
     ! loop over wavelengths
     do iL = 1, nL
        if (slb) then
           fs(iL,iY) = (fsL(iL,iY) -  fsR(iL,iY))/(1.0d0 + ksi)
           us(iL,iY) = (usL(iL,iY) +  usR(iL,iY))/(1.0d0 + ksi)
           fsL(iL,iY) = fsL(iL,iY)/(1.0d0 + ksi)
           fsR(iL,iY) = fsR(iL,iY)/(1.0d0 + ksi)
           ! fs(iL,iY) =  (fsL(iL,iY) - fsR(iL,iY))
        elseif(sph) then
           if (left.eq.1.and.right.eq.0) then
              us(iL,iY) = usL(iL,iY)
              fs(iL,iY) = fsL(iL,iY)
           elseif(left.eq.0.and.right.eq.1) then
              !!** fs(iL,iY) = fsR(iL,iY)  a '-' sign is needed (see below, MN)
              fs(iL,iY) = -fsR(iL,iY)
              us(iL,iY) = usR(iL,iY)
           elseif(typentry(1).ne.5.and.left.eq.1.and.right.eq.1) then
              denom = Ji/(Y(iY)**2.0d0) + Jo
              us(iL,iY) = (Ji*usL(iL,iY)/(Y(iY)**2.0d0) + Jo*usR(iL,iY))/denom
              fsL(iL,iY) = Ji*fsL(iL,iY)/(Y(iY)**2.0d0)/denom
              fsR(iL,iY) = Jo*fsR(iL,iY)/denom
              fs(iL,iY) = fsL(iL,iY) - fsR(iL,iY)
           end if
        end if
     end do
  end do
  do iY = 1, nY
     ! loop over wavelengths
     do iL = 1, nL
        ! here only us needs limit from below; fs can be negative though
        if (us(iL,iY).lt.dyn2) us(iL,iY) = 0.0d0
     end do
  end do
  call bolom(fsL,fsLbol,nY)
  call bolom(fsR,fsRbol,nY)
  call bolom(fs,fsbol,nY)
  error = 0
  goto 999
998 write(12,*)' *** fatal error in dusty! *************************'
  write(12,*)' file with the spectral shape of external radiation:'
  write(12,*)' is missing or not properly formatted?!'
  write(12,*)' ***************************************************'
  error = 3
  !--------------------------------------------------------------------
  deallocate(usL)
  deallocate(usR)
  deallocate(m0)
  deallocate(m1)
  deallocate(m1p)
  deallocate(m1m)
999 return
end subroutine Find_Tran
!***********************************************************************
!!$
!**********************************************************************
double precision function eint1(x)
!======================================================================
! Needed for the slab geometry. It calculates the first exponential
! integral E1(x) by analytical f-la (13.13) from Abramovitz & Stegun(1994)
!                                                         [MN,Dec'97]
! ======================================================================
  implicit none
  integer i
  double precision ac(4),bc(4), cc(6), x, aux, poly, denom
  data ac/8.5733287401d0,18.0590169730d0,8.6347608925d0,0.2677737343d0/
  data bc/9.5733223454d0,25.6329561486d0,21.0996530827d0,3.9584969228d0/
  data cc/-0.57721566d0,0.99999193d0,-0.24991055d0,0.05519968d0,-0.00976004d0,0.00107857d0/
  ! ----------------------------------------------------------------------

  ! For x=1D-15, E1~30 (used below to limit the value at x=0);for x>1, E1<1D-8
  ! Two approximations are used, for x>1 and x<1, respectively
  if(x.lt.0.0d0) x = dabs(x)
  if (x.gt.1.0d0) then
     poly = 0.0d0
     denom = 0.0d0
     aux = 1.0d0
     do i = 1, 4
        poly = poly + ac(5-i)*aux
        denom = denom + bc(5-i)*aux
        aux = aux * x
     end do
     poly = poly + aux
     denom = denom + aux
     eint1 = poly/denom/x*dexp(-x)
  else
     ! if (x.gt.0.0d0.and.x.le.1.0d-15) x=1.0d-15
     poly = 0.0d0
     aux = 1.0d0
     do i = 1, 6
        poly = poly + cc(i)*aux
        aux = aux * x
     end do
     eint1 = poly - dlog(x)
  end if
  ! ----------------------------------------------------------------------
  return
end function eint1
!**********************************************************************

!**********************************************************************
double precision function eint2(x)
!======================================================================
! Needed for the slab geometry. it calculates the second exponential
! integral e2(x) by the recurrence f-la. (see abramovitz & stegun,1994)
!                                                         [MN,Dec'97]
! ======================================================================
  implicit none
  double precision x, eint1
  ! -------------------------------------------------------------------
  if(x.lt.0.0d0) x = dabs(x)
  if (x.lt.1.0d-15) then
     eint2 = 1.0d0
  else
     eint2 = dexp(-x) - x*eint1(x)
  end if
  ! --------------------------------------------------------------------
  return
end function eint2
!**********************************************************************

!**********************************************************************
double precision function eint3(x)
!======================================================================
! Needed for the slab geometry. It calculates the third exponential
! integral E3(x) by the recurrence f-la. (see Abramovitz & Stegun,1994)
!                                                        [MN,Dec'97]
! ======================================================================
  implicit none
  double precision x, eint2
  !---------------------------------------------------------------------
  if(x.le.0.0d0) x = dabs(x)
  if (x.lt.1.0d-15) then
     eint3 = 1.0d0/(2.0d0)
  else
     eint3 = (dexp(-x) - x*eint2(x))/(2.0d0)
  end if
  ! --------------------------------------------------------------------
  return
end function eint3
!**********************************************************************

!***********************************************************************
subroutine Bolom(q,qbol,nY)
!=======================================================================
! This subroutine integrates given radiation field, q (function of
! wavelength and radial coordinate), over wavelength. q is a matrix
! of physical size (npL,npY) [coming from paramet.inc] and real size
! (nL,nY) [coming from grids.inc], and qBol is an array of physical size
! (npY) and real size nY.                              [Z.I., Mar. 1996]
!=======================================================================
  use common
  implicit none
  INTERFACE
     subroutine Simpson(n,n1,n2,x,y,integral)
       integer n, n1, n2
       double precision integral
       double precision,allocatable ::  x(:), y(:)
     end subroutine Simpson
  END INTERFACE
  !---parameter
  integer nY
  double precision, allocatable :: q(:,:), qbol(:)
  !---local
  integer iL, iY
  double precision  resaux
  double precision, allocatable :: qaux(:)
  !---------------------------------------------------------------------
  ! loop over iY (radial coordinate)
  allocate(qaux(nL))
  do iY = 1, nY
     ! generate auxiliary function for integration
     ! loop over iL (wavelength)
     do iL = 1, nL
        qaux(iL) = q(iL,iY)/lambda(iL)
     end do
     call Simpson(nL,1,nL,lambda,qaux,resaux)
     qbol(iY) = resaux
  end do
  !---------------------------------------------------------------------
  deallocate(qaux)
  return
end subroutine Bolom
!***********************************************************************

!***********************************************************************
subroutine SPH_ext_illum(m0,m1,m1p,m1m,nY,nP)
!=======================================================================
! This subroutine finds the transmitted energy density and flux for sphere
! if there is external illumination.                          [Deka, 2008]
!=======================================================================
  use common
  implicit none
  !--- parameter
  integer nY, nP
  double precision,allocatable ::  m0(:,:), m1(:,:),m1p(:,:),m1m(:,:)
  !--- local variables
  integer i,ii, iP,iL, n,nn,nZ, iz, izz, iY,iNloc
  double precision eta, result1, result2, &
        x1, x2, p1, dyn2,p_loc,Yloc1,  expow1, expow2
  double precision,allocatable :: term1(:,:),term2(:,:),term_aux1(:), &
       tauaux(:), z(:,:),angle(:),xg(:),  wg(:)
  external eta

  allocate(term1(nL,nP))
  allocate(term2(nL,nP)) 
  allocate(term_aux1(nP))
  allocate(tauaux(nY)) 
  allocate(z(nP,nY))
  allocate(angle(nP))
!-----------------------------------------------------------------------
  dyn2 = 1.0d-30
  error = 0
  term1 = 0.0d0
  term2 = 0.0d0
  do iY = 1, nY
     do iP = 1, Plast(iY)
        ! upper limit for the counter of z position
        nZ = nY + 1 - iYfirst(iP)
        do iZ = 1, nZ
           tauaux(iZ) = ETAzp(iP,iZ)
        end do
        iZz  = iY + 1 - iYfirst(iP)
        ! loop over wavelengths
        do iL = 1, nL
           expow1 = tautot(iL)*(tauaux(nZ) + tauaux(iZz))
           expow2 = tautot(iL)*(tauaux(nZ) - tauaux(iZz))
           if (expow1.lt.50.0d0) then
              term1(iL,iP) = exp(-expow1)
           else
              term1(iL,iP) = 0.0d0
           end if
           if (expow2.lt.50.0d0) then
              term2(iL,iP) = exp(-expow2)
           else
              term2(iL,iP) = 0.0d0
           end if
        end do
     end do
     do iP = 1, Plast(iY)
        angle(iP) = asin(P(iP)/Y(iY))
     end do
     do iL = 1, nL
        term_aux1 = 0.0d0
        do iP = 1, Plast(iY)
           term_aux1(iP) = (term1(iL,iP) + term2(iL,iP))
        end do
        result1 = 0.0d0
        do iP = 1, Plast(iY) - 1
           result1 = result1 + abs(term_aux1(iP) + term_aux1(iP+1))* &
                abs(cos(angle(iP)) - cos(angle(iP+1)))/2.0d0
        end do
        m0(iL,iY) =  result1
     end do
     do iL = 1, nL
        term_aux1 = 0.0d0
        do iP = 1, Plast(iY)
           term_aux1(iP) = abs(term1(iL,iP))
        end do
        result1 = 0.0d0
        do iP = 1, Plast(iY)-1
           x1 = angle(iP)
           x2 = angle(iP+1)
           nn = 2*int(abs(x2 - x1)/5.0d0) + 11
           if(allocated(xg)) deallocate(xg)
           if(allocated(wg)) deallocate(wg)
           allocate(xg(nn))
           allocate(wg(nn))
           call gauleg(x1,x2,xg,wg,nn)
           do i = 1, nn
              p_loc = xg(i)
              call lininter(Plast(iY),Plast(iY),angle,term_aux1,p_loc,iNloc,Yloc1)
              result1 = result1 + Yloc1*wg(i)*sin(xg(i))*cos(xg(i))
           end do
        end do
        term_aux1 = 0.0d0
        do iP = 1, Plast(iY)
           term_aux1(iP) = abs(term2(iL,iP))
        end do
        result2 = 0.0d0
        do iP = 1, Plast(iY)-1
           x1 = angle(iP)
           x2 = angle(iP+1)
           nn = 2*int(abs(x2 - x1)/5.0d0) + 11
           if(allocated(xg)) deallocate(xg)
           if(allocated(wg)) deallocate(wg)
           allocate(xg(nn))
           allocate(wg(nn))
           call gauleg(x1,x2,xg,wg,nn)
           do i = 1, nn
              p_loc = xg(i)
              call lininter(Plast(iY),Plast(iY),angle,term_aux1,p_loc,iNloc,Yloc1)
              result2 = result2 + Yloc1*wg(i)*sin(xg(i))*cos(xg(i))
           end do
        end do
        m1(iL,iY) = abs(result1 - result2)
     end do
  end do
  !-----------------------------------------------------------------------
  deallocate(term1)
  deallocate(term2) 
  deallocate(term_aux1)
  deallocate(tauaux) 
  deallocate(z)
  deallocate(angle)
  return
end subroutine SPH_ext_illum
!***********************************************************************


!***********************************************************************
subroutine Emission(nY,T4_ext,emiss,emiss_total)
!=======================================================================
! This subroutine calculates emission term from the temperature and abund
! arrays for flag=0, and adds U to it for flag=1.
!                                                      [Z.I., Mar. 1996]
!=======================================================================
  use common
  implicit none
  ! ---- parameter
  integer nY
  double precision, allocatable :: T4_ext(:)
  double precision, allocatable :: emiss(:,:,:),emiss_total(:,:)
  ! ---- local variable 
  integer iG, iY, iL
  double precision  emig, tt, xP,Planck
  ! -------------------------------------------------------------------

  ! first initialize Emiss
  emiss = 0.0d0
  ! calculate emission term for each component and add it to emiss
  ! loop over wavelengths
  !$OMP PARALLEL DO &
  !$OMP PRIVATE(iL,iY,iG,xP,tt,emig)
  do iL = 1, nL
     ! loop over radial coordinate
     do iY = 1, nY
        ! loop over grains
        emiss_total(iL,iY) = 0.
        do iG = 1, nG
           if (destroyed(iG,iY).gt.0.0) then 
              xP = 14400.0d0/(lambda(iL)*Td(iG,iY))
              tt = (Td(iG,iY)**4.0d0) / T4_ext(iY)
              emig = abund(iG,iY)*tt*Planck(xP)
              if (emig.lt.dynrange) emig = 0.0d0
              ! add contribution for current grains
              emiss(iG,iL,iY) = emig
              if (emiss(iG,iL,iY).lt.dynrange) emiss(iG,iL,iY) = 0.0d0
              emiss_total(iL,iY)  = emiss(iG,iL,iY)
           else
              emiss_total(iL,iY)  = 0.0
           end if
        end do
     end do
  end do
  !$OMP END PARALLEL DO
  ! --------------------------------------------------------------------
  return
end subroutine Emission
!***********************************************************************

!***********************************************************************
subroutine Find_Diffuse(nY,nP,initial,moment,iter,iterfbol,T4_ext,us,emiss)
!=======================================================================
! This subroutine finds the diffuse en.density for slab and sphere,
! and the diffuse flux for sphere.                       [Deka,'08]
! Renamed some variables for easier reading:
! faux becomes Sfn_em, faux1 becomes Sfn_sc, removed b(npY),using Utot directly instead.
!                                                        [MN'10]
!=======================================================================
  use common
  implicit none
  INTERFACE
     subroutine SPH_DIFF(flag,moment,nY,nP,initial,iter,iterfbol,T4_ext,emiss,us,vec2)
       integer nY,nP,iter,iterfbol,flag,moment
       double precision, allocatable :: T4_ext(:),emiss(:,:,:),us(:,:),vec2(:,:)
       logical initial
     end subroutine SPH_DIFF
  END INTERFACE
  ! --- parameter
  integer nY,nP,iter,iterfbol,moment
  logical initial
  double precision, allocatable :: T4_ext(:)
  double precision, allocatable :: us(:,:)
  double precision, allocatable :: emiss(:,:,:)
  ! --- local variabels
  integer iL,iY,iP,j, iYaux, kronecker,flag, nZ,iG, moment_loc
  double precision eint2, &
       Sfn_sc(npY),sum1,sum2,dyn2, frac
  double precision,allocatable :: sph_em(:,:),sph_sc(:,:),tau(:),Sfn_em(:)
  external eint2
  !----------------------------------------------------------------------
  allocate(sph_em(nL,nY))
  allocate(sph_sc(nL,nY))
  allocate(tau(nY))
  allocate(Sfn_em(nY))
  !----------------------------------------------------------------------
  dyn2 = 1.0d-30
  error = 0
  Sfn_em(:) = 0.0d0
  Sfn_sc(:) = 0.0d0
  ! for slab calculate new energy density
  ! loop over wavelengths
  !------------for slab ----------------
  if(slb) then 
     do iL = 1, nL
        do iY = 1,nY
           tau(iY) = TAUslb(iL,iY)
           Sfn_em(iY) = 0
           Sfn_sc(iY) = 0
           do iG=1,nG
              frac = destroyed(iG,iY)*(sigmaA(iG,iL)+sigmaS(iG,iL))/(sigmaA(nG+1,iL)+sigmaS(nG+1,iL))
              Sfn_em(iY) = Sfn_em(iY) + frac*(1.0d0-omega(iG,iL))*emiss(iG,iL,iY)/2.0d0
              if (initial.and.iter.eq.1.and.iterfbol.eq.1) then
                 Sfn_sc(iY) = Sfn_sc(iY) + frac*omega(iG,iL)*us(iL,iY)/2.0d0
              else
                 Sfn_sc(iY) = Sfn_sc(iY) + frac*omega(iG,iL)*utot(iL,iY)/2.0d0
              end if
           end do
        end do
        ! integrate to get diffuse emission and scattering
        do iY = 1, nY
           sum1 = 0.0d0
           sum2 = 0.0d0
           do j = 1, nY-1
              sum1 = sum1 + 0.5d0*(Sfn_em(j)+Sfn_em(j+1))* &
                   dabs(eint2(tau(iY)-tau(j))-eint2(tau(iY)-tau(j+1)))
              sum2 = sum2 + 0.5d0*(Sfn_sc(j)+Sfn_sc(j+1))* &
                   dabs(eint2(tau(iY)-tau(j))-eint2(tau(iY)-tau(j+1)))
           end do
           !!** find Ude,Uds
           Ude(iL,iY) = sum1
           Uds(iL,iY) = sum2
        end do
        do iY = 1, nY
           utot(iL,iY) = us(iL,iY) + Ude(iL,iY) + Uds(iL,iY)
        end do
     end do  !end do over lambda
  !-----------for sphere ----------------
  elseif(sph) then 
     do moment_loc = 1, moment
        if (moment_loc.eq.1) then
           ! Find diffuse en.density
           call SPH_diff(0,0,nY,nP,initial,iter,iterfbol,T4_ext,emiss,us,sph_em)
           ! find total energy density
           do iY = 1, nY
              do  iL = 1, nL
                 !utot(iL,iY) = us(iL,iY) + sph_em(iL,iY) + sph_sc(iL,iY)
                 utot(iL,iY) = us(iL,iY) + sph_em(iL,iY)
                 !Ude(iL,iY) = sph_em(iL,iY)
                 !Uds(iL,iY) = sph_sc(iL,iY)
              end do
           end do
        elseif (moment_loc.eq.2) then
           ! Find diffuse fluxes
           call SPH_diff(1,1,nY,nP,initial,iter,iterfbol,T4_ext,emiss,us,fde)
           call SPH_diff(2,1,nY,nP,initial,iter,iterfbol,T4_ext,emiss,us,fds)
        end if
     end do
  end if
  !---------------------------------------------------------------------
  deallocate(sph_em)
  deallocate(sph_sc)
  deallocate(tau)
  deallocate(Sfn_em)
  return
end subroutine Find_Diffuse
!***********************************************************************

!***********************************************************************
subroutine Find_Temp(nY,T4_ext)
!=======================================================================
! This subroutine finds new temperature from Utot.
! Temperature is obtained by solving:
!   f2(y)*f1(y) - g(T) = 0
! where
!   f1(y) = Int(Qabs*Utot*dlambda)
!   f2(y) = T4_ext
!   g(T) = qP(Td)*Td**4
!   T_ext is the effective T of the illuminating radiation
!                                                [ZI'96, MN,'00, Deka, 2008]
!=======================================================================
  use common
  implicit none
  INTERFACE
     subroutine Simpson(n,n1,n2,x,y,integral)
       integer n, n1, n2
       double precision integral
       double precision,allocatable ::  x(:), y(:)
     end subroutine Simpson
     subroutine find_Text(nY,T4_ext)
       integer nY
       double precision, allocatable :: T4_ext(:)
     end subroutine find_Text
  END INTERFACE
  !--- parameter
  integer :: nY
  double precision, allocatable :: T4_ext(:)
  double precision, allocatable :: us(:,:)
  !--- local variables
  integer iG, iY, iL
  double precision  xP,Planck, qpt1,qu1, gg, fnum1,wgth
  double precision, allocatable :: fnum(:),ff(:)
  external Planck
  !---------------------------------------------------------------------
  allocate(fnum(nG))
  allocate(ff(nG))
  ! if T1 given in input:
  if(typentry(1).eq.5) call find_Text(nY,T4_ext)
  ! loop over grains
  do iG = 1, nG
     ! loop over radial positions (solving f1-f2*g(t)=0)
     !OMP PRALLEL DO FIRSTPRIVAT(iG) &
     !OMP PRIVATE(xP,wgth,fnum1,gg)
     do iY = 1, nY
      ! calculate f1 and f2
        fnum1 = 0.
        gg = 0.
        do iL = 1, nL
           fnum(iG) = sigmaA(iG,iL)*utot(iL,iY)/lambda(iL)
           xP = 14400.0d0/lambda(iL)/Td(iG,iY)
           ff(iG) = sigmaA(iG,iL)*Planck(xP)/ lambda(iL)
           !begin simpson
           if (iL.ne.1.and.iL.ne.nL) then
              wgth = 0.5d0*(lambda(iL+1)-lambda(iL-1))
           else
              if (iL.eq.1) wgth = 0.5d0*(lambda(1+1)-lambda(1))
              if (iL.eq.nL) wgth = 0.5d0*(lambda(nL)-lambda(nL-1))
           end if
           ! add contribution to the integral
           fnum1 = fnum1 + fnum(iG)*wgth
           gg = gg + ff(iG)*wgth
           !end simpson
        end do
!        call Simpson(nL,1,nL,lambda,fnum,fnum1)
!        call Simpson(nL,1,nL,lambda,ff,gg)
        if (ISNAN(Td(iG,iY))) print*,'Td NAN',fnum1,T4_ext(iY),gg
        !Td(iG,iY) = (fnum1*T4_ext(iY)/gg)**(1.0d0/4.0d0)
        Td(iG,iY) = (fnum1/gg)**(1.0d0/4.0d0)*T4_ext(iY)**(1.0d0/4.0d0)
        if (Td(iG,iY).gt.Tsub(iG)) then 
           destroyed(iG,iY) = 0.0
        else 
           destroyed(iG,iY) = 1.0
        end if
     end do
     !OMP END PRALLEL DO &
  end do
  !---------------------------------------------------------------------
  deallocate(fnum)
  deallocate(ff)
  return
end subroutine Find_Temp
!***********************************************************************

!***********************************************************************
subroutine Find_Text(nY,T4_ext)
!=======================================================================
  use common
  implicit none
  INTERFACE
     subroutine Simpson(n,n1,n2,x,y,integral)
       integer n, n1, n2
       double precision integral
       double precision,allocatable ::  x(:), y(:)
     end subroutine Simpson
  END INTERFACE
  !--- parameter
  integer :: nY
  double precision, allocatable :: T4_ext(:)
  !--- local
  integer iG,iY,iL
  double precision  xP,Planck, qPT1,qU1,tmp
  double precision, allocatable :: fnum(:),fdenum(:)
!-----------------------------------------------------------------------

  allocate(fnum(nL))
  allocate(fdenum(nL))
! loop over grains
  do iG = 1, nG
     do iL = 1, nL
        fnum(iL) = sigmaA(iG,iL)*utot(iL,1)/lambda(iL)
        xP = 14400.0d0/lambda(iL)/ Tinner(1)
        fdenum(iL) = sigmaA(iG,iL)*Planck(xP)/lambda(iL)
     end do
     call Simpson(npL,1,nL,lambda,fnum,qU1)
     call Simpson(npL,1,nL,lambda,fdenum,qPT1)
     if (sph) then
        do iY = 1, nY
           !T4_ext(iY) = (Tsub(1)**4.0d0)*(qPT1/qU1)/(Y(iY)**2.0d0)
           T4_ext(iY) = (Tinner(1)**4.0d0)*(qPT1/qU1)/(Y(iY)**2.0d0)
        end do
     elseif(slb) then 
        T4_ext = Tinner(1)**4.0d0*(qPT1/qU1)
     end if
  end do
!-----------------------------------------------------------------------
  if (typentry(1).eq.5) then
     Ji = sigma/pi*T4_ext(1)
     !Jo = ksi * Ji
  endif
  deallocate(fnum)
  deallocate(fdenum)
   return
end subroutine Find_Text
!***********************************************************************

!***********************************************************************
subroutine Find_Text_multi(nY,T4_ext)
!=======================================================================
  use common
  implicit none
  INTERFACE
     subroutine Simpson(n,n1,n2,x,y,integral)
       integer n, n1, n2
       double precision integral
       double precision,allocatable ::  x(:), y(:)
     end subroutine Simpson
  END INTERFACE
  !--- parameter
  integer :: nY
  double precision, allocatable :: T4_ext(:)
  !--- local variables
  integer iG,iY,iL
  double precision  xP,Planck, qPT1,qU1
  double precision, allocatable :: fnum(:),fdenum(:)
  !----------------------------------------------------------------------
  allocate(fnum(nL))
  allocate(fdenum(nL))
  ! set to fiducial Grain
  do iL = 1, nL
     fnum(iL) = sigmaA(ifidG,iL)*utot(iL,1)/lambda(iL)
!     print*,sigmaA(ifidG,iL),utot(iL,1)
     xP = 14400.0d0/lambda(iL)/ Tinner(ifidG)
     ! xP = 14400.0d0/lambda(iL)/ Tsub(ifidG)
     fdenum(iL) = sigmaA(ifidG,iL)*Planck(xP)/lambda(iL)
  end do
  call Simpson(nL,1,nL,lambda,fnum,qU1)
  call Simpson(nL,1,nL,lambda,fdenum,qPT1)
  if (sph) then
     do iY = 1, nY
        T4_ext(iY) = (Tinner(ifidG)**4.0d0)*(qPT1/qU1)/(Y(iY)**2.0d0)
        !T4_ext(iY) = (Tsub(ifidG)**4.0d0)*(qPT1/qU1)/(Y(iY)**2.0d0)
     end do
  elseif(slb) then
     do iY = 1, nY
        T4_ext(iY) = Tinner(ifidG)**4.0d0*(qPT1/qU1)
        !T4_ext(iY) = Tsub(ifidG)**4.0d0*(qPT1/qU1)
     end do
  end if
  do iG=1,nG
     if (iG.ne.ifidG) then
        do iL = 1, nL
           fnum(iL) = sigmaA(iG,iL)*utot(iL,1)/lambda(iL)
           xP = 14400.0d0/lambda(iL) / Tinner(ifidG)
           !xP = 14400.0d0/lambda(iL) / Tsub(ifidG)
           fdenum(iL) = sigmaA(iG,iL)*Planck(xP)/lambda(iL)
        end do
        call Simpson(nL,1,nL,lambda,fnum,qU1)
        call Simpson(nL,1,nL,lambda,fdenum,qPT1)
        Tinner(iG) = (T4_ext(1)*qU1/qPT1)**(1.0d0/4.0d0)
     end if
  end do
  if (typentry(1).eq.5) then
     Ji = sigma/pi*T4_ext(1)
     Jo = ksi * Ji
  endif
  !---------------------------------------------------------------------
  deallocate(fnum)
  deallocate(fdenum)
  return
end subroutine Find_Text_multi
!***********************************************************************

!***********************************************************************
subroutine SPH_diff(flag,moment,nY,nP,initial,iter,iterfbol,&
     T4_ext,emiss,us,vec2)
!=======================================================================
! This subroutine finds the diffuse part for both emission and scattering
! for spherical shell                                        [Deka, 2008]
!! If flag=0 vec=Em+Us If flag=1 for vec2=Em; flag=2 for vec2=Us
!!** Introduced integration of diffuse radiation with Nordlund, as in old Dusty.
!!** Restored ETAzp instead of the time-consuming GauLeg integration
!!** New commments added.                                   [MN, July'10]
!=======================================================================
  use omp_lib
  use common
  implicit none
  ! --- parameter
  integer nY,nP,iter,iterfbol,iZaux,flag,moment
  double precision, allocatable :: T4_ext(:),emiss(:,:,:),us(:,:),vec2(:,:)
  logical initial
  ! --- local variables 
  integer iP,iL,iZ,iZz,nZ,iY,iYy,iG,isum,first
  double precision, allocatable :: xN(:), yN(:),S_unscaled(:,:)
  double precision  P_integ, frac, wgth, wcor1,wcor2,&
       daux1,daux2,int1,int2,int3,tmp,eta,wgth1,tauZz,S_fun,&
       S_fun1,S_fun2,m,n,sum_yN,tautot_iL,Wsimp,Wcorr,wC1,wC2
  REAL*8 CUBIC_SPLINT
  external eta
  !-----------------------------------------------------------------------
  allocate(xN(nP))
  allocate(yN(nP))
  allocate(S_unscaled(nL,nY))
  first = 1
  !!** Prepcalcualte unscaled source function
  !$OMP PARALLEL DO PRIVATE(iY,iL,iG)
  do iY = 1, nY
     do iL=1,nL
        S_unscaled(iL,iY) = 0.0D0
        do iG = 1, nG
           frac = destroyed(iG,iY)*(sigmaA(iG,iL)+sigmaS(iG,iL))/(sigmaA(nG+1,iL)+sigmaS(nG+1,iL))
           if ((flag.eq.0).or.(flag.eq.1)) then 
              tmp = frac*(1.0d0-omega(iG,iL))
              S_unscaled(iL,iY) = S_unscaled(iL,iY) + tmp*emiss(iG,iL,iY)*T4_ext(iY)
           end if
           if ((flag.eq.0).or.(flag.eq.2)) then 
              tmp = frac*omega(iG,iL)
              if (initial.and.iter.eq.1) then
                 S_unscaled(iL,iY) = S_unscaled(iL,iY) + tmp*us(iL,iY)*T4_ext(iY)
              else
                 S_unscaled(iL,iY) = S_unscaled(iL,iY) + tmp*utot(iL,iY)*T4_ext(iY)
              end if
           end if
        end do
     end do
  end do
  !!** for each radial grid point calculate the integrals from Sec.4.1 in Blueprint:
  !$OMP END PARALLEL DO
  do iY = 1, nY
     do iP = 1, Plast(iY)
        xN(iP) = sqrt(1.0-(P(iP)/Y(iY)*P(iP)/Y(iY)))
     end do
     !$OMP PARALLEL DO FIRSTPRIVATE(iY) &
     !$OMP PRIVATE(frac,iYy,iP,iZz,nZ) &
     !$OMP PRIVATE(int1,int2,int3,wgth,wcor1,wcor2,daux1,daux2,iZ,iL,m,n,sum_yN,tautot_iL) &
     !$OMP PRIVATE(tmp,tauZz,S_fun,S_fun1,S_fun2,P_integ,Wsimp,Wcorr,wC1,wC2) 
     do iL = 1, nL
        sum_yN = 0.0D0
        tautot_iL = tautot(iL)
        P_integ = 0.0D0
        do iP = 1, Plast(iY)
           iZz  = iY + 1 - iYfirst(iP)  !this is for z in eq.(4.1.5)
           ! upper limit for the counter of z position
           nZ  = nY + 1 - iYfirst(iP)   !nZ is index for zmax=sqrt(Y**2-p**2) [MN]
           int1 = 0.
           int2 = 0.
           int3 = 0.
           tauZz = tautot_iL*ETAzp(iP,iZz)
           if ((nZ.gt.1)) then 
              do iZ = 1, nZ
                 iYy = iYfirst(iP) + iZ - 1
                 if((iZ.eq.1).and.(P(iP).gt.Y(iYy)).and.(P(iP).lt.Y(iYy+1))) then
                    S_fun1 = S_unscaled(iL,iYy+0)
                    S_fun2 = S_unscaled(iL,iYy+1)
                    m = (S_fun2-S_fun1) / (Y(iYy+1)-Y(iYy))
                    n = S_fun1 - m*(Y(iYy))
                    S_fun = m*P(iP)+n
                 else
                    S_fun =  S_unscaled(iL,iYy)
                 end if
                 if (abs(S_fun).gt.dynrange*dynrange) then
                    !call simpson_int(iP,iZ,iZz,nZ,tautot_iL,tauZz,S_fun,int1,int2,int3)
                    !Simpson wgth
                    if ((iZ.ne.1).and.(iZ.ne.nZ)) then 
                       daux1 = exp(-tautot_iL*ETAzp(iP,iZ-1))
                       daux2 = exp(-tautot_iL*ETAzp(iP,iZ+1))
                    ELSE
                       if (iZ.eq.1) then 
                          daux1 = exp(-tautot_iL*ETAzp(iP,iZ))
                          daux2 = exp(-tautot_iL*ETAzp(iP,iZ+1))
                       end if
                       if (iZ.eq.nZ) then 
                          daux1 = exp(-tautot_iL*ETAzp(iP,iZ-1))
                          daux2 = exp(-tautot_iL*ETAzp(iP,iZ))
                       end if
                    END IF
                    wgth = 0.5*(daux2-daux1)
                    int1 = int1 + S_fun*wgth
                    ! add contribution to the integral
                    if ((iZ.le.iZz).and.(iZz.gt.1)) then 
                       if ((iZ.ne.1).and.(iZ.ne.iZz)) then 
                          daux1 = exp(-abs(tauZz-tautot_iL*ETAzp(iP,iZ-1)))
                          daux2 = exp(-abs(tauZz-tautot_iL*ETAzp(iP,iZ+1)))
                       else
                          IF (iZ.eq.1) THEN
                             daux1 = exp(-abs(tauZz-tautot_iL*ETAzp(iP,iZ)))
                             daux2 = exp(-abs(tauZz-tautot_iL*ETAzp(iP,iZ+1)))
                          END IF
                          IF (iZ.eq.iZz) THEN  
                             daux1 = exp(-abs(tauZz-tautot_iL*ETAzp(iP,iZ-1)))
                             daux2 = exp(-abs(tauZz-tautot_iL*ETAzp(iP,iZ)))
                          END IF
                       end if
                       wgth = 0.5*(daux2-daux1)
                       int2 = int2 + S_fun*wgth
                    end if
                    if ((iZ.ge.iZz).and.(iZz.lt.nZ)) then 
                       if ((iZ.ne.iZz).and.(iZ.ne.nZ)) then 
                          daux1 = exp(-abs(tautot_iL*ETAzp(iP,iZ-1)-tauZz))
                          daux2 = exp(-abs(tautot_iL*ETAzp(iP,iZ+1)-tauZz))
                       else
                          if (iZ.eq.iZz) then 
                             daux1 = exp(-abs(tautot_iL*ETAzp(iP,iZ)-tauZz))
                             daux2 = exp(-abs(tautot_iL*ETAzp(iP,iZ+1)-tauZz))
                          end if
                          if (iZ.eq.nZ) then 
                             daux1 = exp(-abs(tautot_iL*ETAzp(iP,iZ-1)-tauZz))
                             daux2 = exp(-abs(tautot_iL*ETAzp(iP,iZ)-tauZz))
                          end if
                       end if
                       wgth = 0.5*(daux2-daux1)
                       int3 = int3 + S_fun*wgth
                    end if
                    !spline_simpson_int is not better than normal simpson
                    !call spline_simpson_int(iP,iZ,iZz,nZ,tautot_iL,tauZz,S_fun,int1,int2,int3)
                 end if
              end do
           end if
           ! generate intensity array for NORDLUND
           if (moment.eq.0) tmp = abs(int1)*exp(-tauZz) + abs(int2) + abs(int3)
           if (moment.eq.1) tmp = abs(abs(int1)*exp(-tauZz) + abs(int2) - abs(int3))
           !yN(iP) = tmp
           sum_yN = sum_yN + tmp
           ! Begin former Nordlund
           wSimp = 0.0
           wCorr = 0.0
           wC1 = 0.0
           wC2 = 0.0
           IF ((iP.GT.1).AND.(iP.LT.Plast(iY))) THEN
              wSimp = 0.5 * (xN(iP+1)-xN(iP-1))
           ELSE
              IF (iP.eq.1) wSimp = 0.5 * (xN(1+1)-xN(1))
              IF (iP.eq.Plast(iY)) wSimp = 0.5 * (xN(Plast(iY))-xN(Plast(iY)-1))
           END IF
           ! ... and then correction term for cubic spline (Nordlund, eq. III-14,
           ! second term and eq. III-16) (wC1 and wC2 are auxiliary quantities)
           IF (iP.GT.1+1) THEN
              wC1 = xN(iP) - 2.0*xN(iP-1) + xN(iP-2)
           ELSE
              IF (iP.EQ.1) wC1 =  first * (xN(1+1) - xN(1))
              IF (iP.EQ.1+1) wC1 = first * (xN(1) - xN(1+1))
           ENDIF
           IF (iP.LE.(Plast(iY)-2)) THEN
              wC2 = xN(iP+2) - 2.0*xN(ip+1) + xN(iP)
           ELSE
              IF (iP.EQ.(Plast(iY)-1)) wC2 = first * (xN(Plast(iY)-1) - xN(Plast(iY)))
              IF (iP.EQ.Plast(iY)) wC2 = first * (xN(Plast(iY)) - xN(Plast(iY)-1))
           ENDIF
           wCorr = (wC1 - wC2) / 12.
           ! add contribution to the integral
           if (moment.eq.0) P_integ = P_integ - tmp * (wSimp + wCorr)
           if (moment.eq.1) P_integ = P_integ - xN(iP) * tmp * (wSimp + wCorr)
           ! End Norlund
        end do
        if (sum_yN.gt.dynrange*dynrange) then
           !CALL NORDLUND(nY,nP,0,xN,yN,1,Plast(iY),0,P_integ)
           if (moment.eq.0) vec2(iL,iY) = 0.5*(P_integ)/T4_ext(iY)
           if (moment.eq.1) vec2(iL,iY) = 2*pi*(P_integ)/T4_ext(iY)
        else
           vec2(iL,iY) = 0.0
        end if
     end do  !end do over lambda
     !$OMP END PARALLEL DO
  end do !end do over radial grid Y(iY)  [MN]
  !----------------------------------------------------------------------
  deallocate(xN)
  deallocate(yN)
999 return
end subroutine SPH_diff
!***********************************************************************

subroutine simpson_int(iP,iZ,iZz,nZ,tautot_iL,tauZz,S_fun,int1,int2,int3)
  use common
  implicit none
  !---parameter
  integer :: iP,iZ,nZ,iZz
  double precision :: tautot_iL,int1,int2,int3,S_fun,tauZz
  !--- local
  double precision ::  daux1,daux2,wgth
  !Simpson wgth
  if ((iZ.ne.1).and.(iZ.ne.nZ)) then 
     daux1 = exp(-tautot_iL*ETAzp(iP,iZ-1))
     daux2 = exp(-tautot_iL*ETAzp(iP,iZ+1))
  ELSE
     if (iZ.eq.1) then 
        daux1 = exp(-tautot_iL*ETAzp(iP,iZ))
        daux2 = exp(-tautot_iL*ETAzp(iP,iZ+1))
     end if
     if (iZ.eq.nZ) then 
        daux1 = exp(-tautot_iL*ETAzp(iP,iZ-1))
        daux2 = exp(-tautot_iL*ETAzp(iP,iZ))
     end if
  END IF
  wgth = 0.5*(daux2-daux1)
  int1 = int1 + S_fun*wgth
  ! add contribution to the integral
  if ((iZ.le.iZz).and.(iZz.gt.1)) then 
     if ((iZ.ne.1).and.(iZ.ne.iZz)) then 
        daux1 = exp(-abs(tauZz-tautot_iL*ETAzp(iP,iZ-1)))
        daux2 = exp(-abs(tauZz-tautot_iL*ETAzp(iP,iZ+1)))
     else
        IF (iZ.eq.1) THEN
           daux1 = exp(-abs(tauZz-tautot_iL*ETAzp(iP,iZ)))
           daux2 = exp(-abs(tauZz-tautot_iL*ETAzp(iP,iZ+1)))
        END IF
        IF (iZ.eq.iZz) THEN  
           daux1 = exp(-abs(tauZz-tautot_iL*ETAzp(iP,iZ-1)))
           daux2 = exp(-abs(tauZz-tautot_iL*ETAzp(iP,iZ)))
        END IF
     end if
     wgth = 0.5*(daux2-daux1)
     int2 = int2 + S_fun*wgth
  end if
  if ((iZ.ge.iZz).and.(iZz.lt.nZ)) then 
     if ((iZ.ne.iZz).and.(iZ.ne.nZ)) then 
        daux1 = exp(-abs(tautot_iL*ETAzp(iP,iZ-1)-tauZz))
        daux2 = exp(-abs(tautot_iL*ETAzp(iP,iZ+1)-tauZz))
     else
        if (iZ.eq.iZz) then 
           daux1 = exp(-abs(tautot_iL*ETAzp(iP,iZ)-tauZz))
           daux2 = exp(-abs(tautot_iL*ETAzp(iP,iZ+1)-tauZz))
        end if
        if (iZ.eq.nZ) then 
           daux1 = exp(-abs(tautot_iL*ETAzp(iP,iZ-1)-tauZz))
           daux2 = exp(-abs(tautot_iL*ETAzp(iP,iZ)-tauZz))
        end if
     end if
     wgth = 0.5*(daux2-daux1)
     int3 = int3 + S_fun*wgth
  end if
end subroutine simpson_int

subroutine spline_simpson_int(iP,iZ,iZz,nZ,tautot_iL,tauZz,S_fun,int1,int2,int3)
  use common
  implicit none
  !---parameter
  integer :: iP,iZ,nZ,iZz
  double precision :: tautot_iL,int1,int2,int3,S_fun,tauZz
  !--- local
  double precision ::  wcor1,wcor2,&
       wgth,first,tau1,tau2,tau3,tau4,tau5
  first = 0.
  if (iZ.gt.2)    tau1 = tautot_iL*ETAzp(iP,iZ-2)
  if (iZ.gt.1)    tau2 = tautot_iL*ETAzp(iP,iZ-1)
  if (iZ.gt.0)    tau3 = tautot_iL*ETAzp(iP,iZ)
  if (iZ.lt.nZ)   tau4 = tautot_iL*ETAzp(iP,iZ+1)
  if (iZ.lt.nZ-1) tau5 = tautot_iL*ETAzp(iP,iZ+2)
  !---first integral
  if ((iZ.ne.1).and.(iZ.ne.nZ)) then 
     wgth = 0.5*(exp(-tau4)-exp(-tau2))
  ELSE
     IF (iZ.eq.1)  wgth = 0.5*(exp(-tau4)-exp(-tau3))
     IF (iZ.eq.nZ) wgth = 0.5*(exp(-tau3)-exp(-tau2))
  END IF
  ! spline correction
  if (iZ.gt.2) then 
     wcor1 = exp(-tau3)-2.0D0*exp(-tau2)+exp(-tau1)
  else
     if (iZ.eq.1) wcor1 = first*(exp(-tau4)-exp(-tau3))
     if (iZ.eq.2) wcor1 = first*(exp(-tau3)-exp(-tau4))
  end if
  if (iZ.le.nZ-2) then 
     wcor2 = exp(-tau5)-2.0D0*exp(-tau4)+exp(-tau3)
  else
     if (iZ.eq.nZ-1) wcor2 = first*(exp(-tau4)-exp(-tau5))
     if (iZ.eq.nZ)   wcor2 = first*(exp(-tau5)-exp(-tau4))
  end if
  wgth = wgth + (wcor1-wcor2)/12.
  int1 = int1 + S_fun*(wgth+(wcor1-wcor2)/12.)
  !----second integral
  if ((iZ.le.iZz).and.(iZz.gt.1)) then 
     if ((iZ.ne.1).and.(iZ.ne.iZz)) then 
        wgth = 0.5*(exp(-abs(tau4-tauZz))-exp(-abs(tau2-tauZz)))
     ELSE
        IF (iZ.eq.1)   wgth = 0.5*(exp(-abs(tau4-tauZz))-exp(-abs(tau3-tauZz)))
        IF (iZ.eq.iZz) wgth = 0.5*(exp(-abs(tau3-tauZz))-exp(-abs(tau2-tauZz)))
     END IF
     !spline correction 
     if (iZ.gt.2) then 
        wcor1 = exp(-abs(tau3-tauZz))-2.0D0*exp(-abs(tau2-tauZz))+exp(-abs(tau1-tauZz))
     else
        if (iZ.eq.1) wcor1 = first*(exp(-abs(tau4-tauZz))-exp(-abs(tau3-tauZz)))
        if (iZ.eq.2) wcor1 = first*(exp(-abs(tau3-tauZz))-exp(-abs(tau4-tauZz)))
     end if
     if (iZ.le.iZz-2) then 
        wcor2 = exp(-abs(tau5-tauZz))-2.0D0*exp(-abs(tau4-tauZz))+exp(-abs(tau3-tauZz))
     else
        if (iZ.eq.iZz-1) wcor2 = first*(exp(-abs(tau4-tauZz))-exp(-abs(tau5-tauZz)))
        if (iZ.eq.iZz)   wcor2 = first*(exp(-abs(tau5-tauZz))-exp(-abs(tau4-tauZz)))
     end if
     wgth = wgth + (wcor1-wcor2)/12.
     int2 = int2 + S_fun*wgth
  end if
  if ((iZ.ge.iZz).and.(iZz.lt.nZ)) then 
     if ((iZ.ne.1).and.(iZ.ne.iZz)) then 
        wgth = 0.5*(exp(-abs(tauZz-tau4))-exp(-abs(tauZz-tau2)))
     ELSE
        IF (iZ.eq.1)   wgth = 0.5*(exp(-abs(tauZz-tau4))-exp(-abs(tauZz-tau3)))
        IF (iZ.eq.iZz) wgth = 0.5*(exp(-abs(tauZz-tau3))-exp(-abs(tauZz-tau2)))
     END IF
     !spline correction
     if (iZ.gt.2) then 
        wcor1 = exp(-abs(tauZz-tau3))-2.0D0*exp(-abs(tauZz-tau2))+exp(-abs(tauZz-tau1))
     else
        if (iZ.eq.1) wcor1 = first*(exp(-abs(tauZz-tau4))-exp(-abs(tauZz-tau3)))
        if (iZ.eq.2) wcor1 = first*(exp(-abs(tauZz-tau3))-exp(-abs(tauZz-tau4)))
     end if
     if (iZ.le.iZz-2) then 
        wcor2 = exp(-abs(tau5-tauZz))-2.0D0*exp(-abs(tauZz-tau4))+exp(-abs(tauZz-tau3))
     else
        if (iZ.eq.iZz-1) wcor2 = first*(exp(-abs(tauZz-tau4))-exp(-abs(tauZz-tau5)))
        if (iZ.eq.iZz)   wcor2 = first*(exp(-abs(tauZz-tau5))-exp(-abs(tauZz-tau4)))
     end if
     wgth = wgth + (wcor1-wcor2)/12.
     int3 = int3 + S_fun*wgth
  end if


  return 
end subroutine spline_simpson_int

!***********************************************************************
SUBROUTINE NORDLUND(nY,nP,flag,x,f,N1,N2,m,intfdx)
!=======================================================================
! This subroutine calculates integral I(x**m*y(x)*dx). Both y and x are
! 1D arrays, y(i), x(i) with i=1,npP (npP comes from 'paramet.inc'). Lower
! and upper integration limits are x(N1) and x(N2), respectively. The
! method used is approximation of y(x) by a piecewise cubic spline (see
! Nordlund: Spherical Transfer with Single-Ray Approximation, in
! 'Methods in radiative transfer', ed. W. Kalkofen, Cambridge University
! Press, 1984). The resulting integral is sum of y(i)*w(i), i=N1,N2.
! Weights w(i) are determined from the array x(i). Here w(i) = wSimp(i)
! + wCorr(i).
! To improve accuracy of angular integration in radiative
! transfer, if flag=1 the contribution of the last Nanal (specified
! below) points is calculated in subroutine ANALINT which fits a
! function of the form: y = P(x) + d/sqrt(1-x*x), where P(x) is the
! polynomial of order Nanal-1, to these points and evaluates integral
! analytically.                                        [Z.I., Nov. 1995]
! =======================================================================
  use common
  IMPLICIT none
  !---parameter
  integer :: nY, nP, flag, N1, N2, m
  double precision :: intfdx
  double precision, allocatable :: x(:),f(:)
  !---local
  INTEGER i, N2n, Nanal, first
  DOUBLE PRECISION  wC1, wC2, am, xaux(4), faux(4), aux
  double precision, allocatable :: wSimp(:), wCorr(:)
  ! ------------------------------------------------------------------
  allocate(wSimp(nP))
  allocate(wCorr(nP))
  error = 0
  ! parameter 'first' selects choice for derivatives at boundary points.
  ! For first.EQ.0 derivatives are 0 and first*(f2-f1)/(x2-x1) otherwise.
  ! first=1 works much better for functions encountered here.
  first = 1
  ! number of points for analytic integration
  Nanal = 4
  ! do not include points for analytic integration
  IF ((flag.EQ.1).AND.(N2.GT.N1+Nanal)) THEN
     N2n = N2 - Nanal + 1
  ELSE
     N2n = N2
  END IF
  ! set integral to 0 and accumulate result in the loop
  intfdx = 0.0
  ! generate weighting factors, w(i), and integrate in the same loop
  DO i = N1, N2n
     ! first usual Simpson factors (Nordlund, eq. III-14, first term)...
     IF ((i.NE.N1).AND.(i.NE.N2n)) THEN
        wSimp(i) = 0.5 * (x(i+1)-x(i-1))
     ELSE
        IF (i.eq.N1) wSimp(i) = 0.5 * (x(N1+1)-x(N1))
        IF (i.eq.N2n) wSimp(i) = 0.5 * (x(N2n)-x(N2n-1))
     END IF
     ! ... and then correction term for cubic spline (Nordlund, eq. III-14,
     ! second term and eq. III-16) (wC1 and wC2 are auxiliary quantities)
     IF (i.GT.N1+1) THEN
        wC1 = x(i) - 2.0*x(i-1) + x(i-2)
     ELSE
        IF (i.EQ.N1) wC1 =  first * (x(N1+1) - x(N1))
        IF (i.EQ.N1+1) wC1 = first * (x(N1) - x(N1+1))
     ENDIF
     IF (i.LE.(N2n-2)) THEN
        wC2 = x(i+2) - 2.0*x(i+1) + x(i)
     ELSE
        IF (i.EQ.(N2n-1)) wC2 = first * (x(N2n-1) - x(N2n))
        IF (i.EQ.N2n) wC2 = first * (x(N2n) - x(N2n-1))
     ENDIF
     wCorr(i) = (wC1 - wC2) / 12.
     ! add contribution to the integral
     IF (m.EQ.0) THEN
        intfdx = intfdx + f(i) * (wSimp(i) + wCorr(i))
     ELSE IF(m.EQ.1) THEN
        intfdx = intfdx + x(i)*f(i)*(wSimp(i) + wCorr(i))
     ELSE IF(m.EQ.2) THEN
        intfdx = intfdx + x(i)*x(i)*f(i)*(wSimp(i) + wCorr(i))
     END IF
  END DO
  ! change the sign (x [i.e. mu] array is in descending order!!!)
  intfdx = -intfdx
  ! if flag=1 use analytic approximation for the last Nanal points
  IF ((flag.EQ.1).AND.(N2n.GT.N1+Nanal)) THEN
     ! generate auxiliary arrays for ANALINT
     DO i=1,Nanal
        xaux(i) = x(N2n+Nanal-i)
        faux(i) = f(N2n+Nanal-i)
     END DO
     ! calculate the contribution of the last Nanal points
     ! produce REAL copy of m
     am = 1.0*(m)
     CALL ANALINT(nY,Nanal,xaux,faux,am,aux)
     IF(error.NE.0) THEN
        RETURN
     END IF
     ! add the contribution of the last Nanal points
     intfdx = intfdx + aux
  END IF
  ! ----------------------------------------------------------------
  deallocate(wSimp)
  deallocate(wCorr)
  RETURN
END SUBROUTINE NORDLUND
!**********************************************************************
!!$
!***********************************************************************
SUBROUTINE setupETA(nY,nYprev,itereta)
!=======================================================================
! This subroutine finds spline coefficients ETAcoef such that
! the normalized density function ETA(Y(iY)) is:
! ETAcoef(iY,1)+ETAcoef(iY,2)/Y(iY)+...+ETAcoef(iY,2)/Y(iY)^3
! If spline approximation differs more than maxerr (see below) at the
! midpoint, then a straight line is used instead. (In case of wavelength
! depend. ETA, use ETAfun where any new dens. laws should be described).
! Coefficients ETAcoef are later used in getETAzp to calculate ETAzp.
!                                                [ZI, Feb'96; MN,Aug'97]
! =======================================================================
  use common
  implicit none
  INTEGER iY, nY,nYprev,itereta, iCoeff
  DOUBLE PRECISION coef(npY,4), ETA, maxerr, Ymid, Yinverse(npY), &
       ETAaux(npY), ETAmid(npY)
! -----------------------------------------------------------------------
  ! generate input function for SPLINE2
  DO iY = 1, nY
     Yinverse(iY) = 1. / Y(iY)
     ETAaux(iY) = ETA(Y(iY),nY,nYprev,itereta)
     IF (iY.LT.nY) THEN
        Ymid = dsqrt(Y(iY)*Y(iY+1))
        ETAmid(iY) = ETA(Ymid,nY,nYprev,itereta)
     END IF
  END DO
  ! calculate spline coefficients
  CALL SPLINE2(Yinverse,ETAaux,nY,coef)
  ! check and fix spline coefficients
  maxerr = 0.1
  ! RDW is initialized in Input
  CALL CHKSPLIN(Yinverse,ETAaux,ETAmid,nY,coef,maxerr)
  ! copy coefficients to the output array ETAcoef
  DO iY = 1, nY
     DO iCoeff = 1, 4
        ETAcoef(iY,iCoeff) = coef(iY,iCoeff)
     END DO
  END DO
! -----------------------------------------------------------------------
  RETURN
END SUBROUTINE setupETA
!***********************************************************************

! ***********************************************************************
SUBROUTINE CHKSPLIN(x,fun,funmid,N,coef,maxerr)
! ======================================================================
! This subroutine checks the spline coefficients coef(i,j):
! fun(x)=coef(i,1) + coef(i,2)*x + coef(i,3)*x^2 + coef(i,4)*x^3,
! for x(i).LE.x.LE.x(i+1) with i=1..N. Array funmid(1..N-1) contains the
! values of function fun at mid points defined as
! xmid(i)=SQRT(x(i)*x(i+1). If spline approximation produces error
! greater than maxerr, or funmid<0, a straight line is produced between
! x(i) and x(i+1).                                   [Z.I., Feb. 1995]
! ======================================================================
  use common
  IMPLICIT none
  INTEGER N, i, iCoeff
  DOUBLE PRECISION x(npY), fun(npY), funmid(npY), coef(npY,4),   &
       maxerr, error_spl, slope, xmid, funSpline, aux, power, yR, yL
  ! -------------------------------------------------------------------
  ! check the midpoints
  DO i = 1, N - 1
     xmid = dsqrt(x(i)*x(i+1))
     funSpline = 0.0
     DO iCoeff = 1,4
        IF (xmid.EQ.0.0.AND.iCoeff.EQ.1) THEN
           aux = 1.0
        ELSE
           aux = xmid**(float(iCoeff)-1.0)
        END IF
        funSpline = funSpline + coef(i,iCoeff)*aux
     END DO
     error_spl = DABS((funSpline-funmid(i))/funmid(i))
     ! check for the deviation at the midpoint
     IF (error_spl.GE.maxerr.OR.funSpline.LE.0.0) THEN
        slope = (fun(i+1) - fun(i)) / (x(i+1)-x(i))
        coef(i,1) = fun(i) - x(i) * slope
        coef(i,2) = slope
        coef(i,3) = 0.0
        coef(i,4) = 0.0
     END IF
     ! check for the logarithmic derivative (only for RDW)
     IF(denstyp.eq.3) THEN !denstyp=3 -> RDW
        yL = fun(i)
        yR = fun(i+1)
        IF (x(i)*x(i+1).GT.0.AND.yL*yR.GT.0) THEN
           power = log(yR/yL)/log(x(i+1)/x(i))
           IF (abs(power).GT.10.) THEN
              slope = (yR - yL) / (x(i+1)-x(i))
              coef(i,1) = yL - x(i) * slope
              coef(i,2) = slope
              coef(i,3) = 0.0
              coef(i,4) = 0.0
           END IF
        END IF
     END IF
  END DO
  ! --------------------------------------------------------------------
  RETURN
end subroutine CHKSPLIN
! ***********************************************************************
!!$
!***********************************************************************
SUBROUTINE getETAzp(nY,nP)
!=======================================================================
! This function calculates ETAzp(iP,iZ) along the line of sight with
! impact parameter P(iP) and iZ=1, nZ. Here iZ = 1 corresponds to z=0
! and iZ=nZ to the outer edge. Other grid points coincide with the
! radial grid. The method used is spline approximation for normalized
! density distribution ETA, with subsequent z-integration performed
! analytically in function IntETA
!                                               [ZI,Feb'95; MN,Aug'97]
! =======================================================================
  use common
  implicit none
  INTEGER iP, nZ, iZ, iW, nY, nP
  DOUBLE PRECISION IntETA, auxEta, w1, w2
  EXTERNAL IntEta
! -----------------------------------------------------------------------
  ! loop over impact parameters
  DO iP = 1, nP
     ! maximal number of points along tangential position, z
     nZ = nY + 1 - iYfirst(iP)
     ! starting values for z and ETAzp(iP,iZ)
     IF (P(iP).GE.1.0) THEN
        w2 = P(iP)
     ELSE
        w2 = 1.0
     END IF
     ! initialize ETAzp(iP,iZ)*TAUtot(iL)
     ETAzp(iP,1) = 0.0
     ! loop over z
     DO iZ = 2, nZ
        ! index for local radius, w2
        iW = iYfirst(iP) + iZ - 1
        ! limits for integration
        w1 = w2
        w2 = Y(iW)
        ! find next step in ETAzp
        auxEta = IntETA(P(iP),iW-1,w1,w2)
        ! add next step in ETAzp
        ETAzp(iP,iZ) = ETAzp(iP,iZ-1) + auxEta
     END DO
  END DO
! -----------------------------------------------------------------------
  RETURN
END SUBROUTINE getETAzp
!***********************************************************************

!***********************************************************************
DOUBLE PRECISION FUNCTION IntETA(paux,iW1,w1,w)
!=======================================================================
! This function calculates the integral over the normalized dens. prof.
! along the line of sight with impact parameter p and between the points
! corresponding to y=w1 and y=w. The method used is spline approximation
! for normalized density distribution ETA and subsequent integration
! performed analytically by MAPLE (these results are given through
! soubroutine Maple3).                         [ZI,Feb'96,MN,Aug'97]
! =======================================================================
  use common
  implicit none
  INTEGER iW1, iCoeff
  DOUBLE PRECISION  paux, w1, w, aux(4), z, z1, aux1(4)
! -----------------------------------------------------------------------
  z = dsqrt(w*w-paux*paux)
  z1 = dsqrt(w1*w1-paux*paux)
  ! integrals calculated by MAPLE
  CALL Maple3(w,z,paux,aux)
  CALL Maple3(w1,z1,paux,aux1)
  DO iCoeff = 1, 4
     aux(iCoeff) = aux(iCoeff) - aux1(iCoeff)
  END DO
  IntETA = 0.0
  DO iCoeff = 1, 4
     IntETA = IntETA + ETAcoef(iW1,iCoeff) * aux(iCoeff)
  END DO
! -----------------------------------------------------------------------
  RETURN
END FUNCTION IntETA
!***********************************************************************

!***********************************************************************
subroutine Init_Temp(nY,T4_ext,us)
!=======================================================================
! This subroutine calculates the initial approximation for the temperature.
! Temperature is obtained by solving:
!   g(T) = f2(y)*f1(y)
! where
!   g(T) = qP(Td)*Td**4
!   f1(y) = Int(Qabs*Utot*dlambda)
!   f2(y) = T4_ext
!   T_ext is the effective T of the illuminating radiation [Deka, July'08]
!
!!**   Minor editing to avoid having IF's inside do-loops [MN, Aug'10]
!=======================================================================
  use common
  implicit none
  INTERFACE
     subroutine Simpson(n,n1,n2,x,y,integral)
       integer n, n1, n2
       double precision integral
       double precision,allocatable ::  x(:), y(:)
     end subroutine Simpson
  END INTERFACE
  !--- Parameter
  integer nY
  double precision,allocatable :: us(:,:),T4_ext(:)
  !--- locale variables
  integer iL, iY, iG, iw
  double precision xP, Planck, fnum1,qP
  double precision,allocatable ::fnum(:),ff(:)

  allocate(fnum(nL))
  allocate(ff(nL))
  
  !--------------------------------------------------------------------------
  if(typentry(1).eq.5) then
     ! this is if Tinner(ifidG) given in input
     if(sph) then
        do iY = 1, nY
           T4_ext(iY) = Tinner(ifidG)**4.0d0/Y(iY)**2.0d0   !eq.(4.1.22)
        end do
     else if (slb) then
        ! for slab
        do iY = 1, nY
           T4_ext(iY) = Tinner(ifidG)**4.0d0
        end do
     end if
     ! first approximation for temperature
     do iG = 1, nG
        do iY = 1, nY
           Td(iG,iY) = T4_ext(iY)**(1.0d0/4.0d0)
        end do
     end do
  else
     do iG = 1, nG
        ! loop over radial positions
        do iY = 1, nY
           Td(iG,iY) = T4_ext(iY)**(1.0d0/4.0d0)
           ! calculate fnum1 and qP integrals
           do iL = 1, nL
              fnum(iL) = sigmaA(iG,iL)*us(iL,iY)/lambda(iL)
              xP = 14400.0d0/(lambda(iL)*Td(iG,iY))
              ff(iL) = sigmaA(iG,iL)*Planck(xP)/lambda(iL)
           end do
           call Simpson(nL,1,nL,lambda,fnum,fnum1)
           call Simpson(nL,1,nL,lambda,ff,qP)
           ! get initial temperature
           Td(iG,iY) = (T4_ext(iY)*fnum1/qP)**(1.0d0/4.0d0)
        end do
     end do
  end if
  !--------------------------------------------------------------------
  deallocate(fnum)
  deallocate(ff)
  return
end subroutine Init_Temp
!***********************************************************************
!!$
!***********************************************************************
subroutine Flux_Consv(nY,nYprev,Ncav,itereta,iterfbol,fbolom,fbol_em,fbol_sc,fbolOK,maxrat)
!=======================================================================
! Replaces the former SUBROUTINE ChkFlux(flux,tolern,consfl,error,ETAzp)
!
! This subroutine checks the bolometric diffuse flux conservation at any
! point of the grid. In case of nonconservation inserts a number of points
!  at certain places.                            [Deka'08, MN'99; ZI'96]
!=======================================================================
  use common
  implicit none
  INTERFACE
     SUBROUTINE FindErr(nY,flux,maxFerr)
       integer nY
       DOUBLE PRECISION maxFerr
       double precision,allocatable :: flux(:)
     END SUBROUTINE FindErr
  END INTERFACE
  !---parameter
  integer nY,nYprev,Ncav,fbolOK,itereta,iterfbol
  double precision :: maxrat
  double precision, allocatable :: fbolom(:),fbol_em(:),fbol_sc(:)
  !---local
  integer :: iY,flag, kins, istop,i_ins,n_ins,iG,iL,idm,k,j
  integer,allocatable :: iYins(:)
  double precision :: eta,deltaumax,temp_mean,ee,Yloc,tmp1,tmp2,avg_flux,&
       fact,ff,ffold,fmed,devfac,devmax
  double precision, allocatable :: ratio(:),tauaux(:),etatemp(:),Yins(:),tmp(:)
  real*8 :: median
  external eta
  !---------------------------------------------------------------------
  allocate(ratio(npY))
  allocate(tauaux(npY))
  allocate(etatemp(npY))
  allocate(Yins(npY))
  allocate(iYins(npY))
  flag= 0
  error = 0.0d0
  kins = 0
  maxrat = 0.0d0
  if(sph) then
     ! save old grid and values of Eta (important for denstyp = 5 or 6)
     ! for spherical case
     if (denstyp.eq.3) then !3(RDW)
        Yprev = y
        etatemp = etadiscr
        nYprev = nY
     end if
     do iY = 1, nY
        tauaux(iY) = TAUtot(1)*ETAzp(1,iY)
     end do
     !!** I am not sure if deltaumax has to be found as below: [MN]
     !  maximal deltau is no more than 4 times the average value (2 refinements)
     delTAUmax = 4.0d0*tauaux(nY)/nY
  elseif(slb) then
     delTAUmax = 4.0d0*tautot(1)/nY
  end if
  call FindErr(nY,fbolom,maxrat)
  IF (maxval(fbolom).lt.dynrange**2) maxrat = dynrange
  fmed = MEDIAN(fbolom,nY)
  devmax = 0.0D+00                                                                                  
  ff = 0.0D+00                                                                                      
  istop = 0                                                                                     
  devfac = 0.1D+00                                                                                  
  DO iY = 2, nY                                                                                 
     IF (dabs(fbolom(iY)-fmed).GT.devmax) devmax = dabs(fbolom(iY)-fmed)
  END DO
!!$  if (maxrat.gt.accFlux) then 
!!$     DO WHILE (istop.ne.1)
!!$        DO iY = 2, nY                     
!!$           !ffold = ff
!!$          ffold = dabs(fbolom(iY-1) - fmed)                                                          
!!$          ff = dabs(fbolom(iY) - fmed)
!!$          flag = 0                                                                                  
!!$          !if any of these criteria is satisfied insert a point:                                     
!!$          !1) if error is increasing too fast
!!$          if ((sph).and.(left.eq.1).and.(iY.eq.2).and.(iterfbol.le.5)) flag = 1
!!$          IF (abs(ff-ffold).GT.devfac*devmax) flag = 1                                              
!!$          !2) if delTAU is too large                                                                 
!!$          IF (TAUtot(1)*(ETAzp(1,iY)-ETAzp(1,iY-1)).GT.delTAUmax) flag = 1
!!$          IF(flag.EQ.1.AND.devmax.GE.accFlux) THEN                                                   
!!$            kins = kins + 1                                                                         
!!$            Yins(kins) = Y(iY-1)+0.5D+00*(Y(iY)-Y(iY-1))                                                
!!$            iYins(kins) = iY-1                                                                      
!!$          END IF                                                                                    
!!$        END DO                                                                                      
!!$        IF (devmax.LT.accFlux.OR.devfac.LT.0.01D+00) THEN                                                
!!$          istop = 1                                                                                 
!!$          ELSE                                                                                      
!!$          IF (kins.GT.0) istop = 1                                                                  
!!$        END IF                                                                                      
!!$        devfac = devfac / 2.0D+00                                                                       
!!$     END DO
!!$  endif


  if (maxrat.gt.accFlux) then 
!!$     kins = kins + 1
!!$     Yins(kins) = 0.5*(Y(2)+Y(1))
!!$     iYins(kins) = 1
     tmp2 = 0.0
     do iY = 2, nY
        tmp1 = dabs(fbolom(iY-1)-fbolom(iY))
        tmp2 = tmp2 + tmp1
        !if (tmp1.gt.tmp2) tmp2 = tmp1
     end do
     fact = -0.00001
     n_ins = nY
     do while (n_ins.gt.nY*0.3)
        fact = fact+0.00001
        n_ins = 0
        do iY = 2, nY
           tmp1 = dabs(fbolom(iY-1)-fbolom(iY))
           if  (tmp1.ge.fact*tmp2) n_ins = n_ins + 1
        end do
     end do
     !if # of added points is less than 5% of nY 
     ! or less than 2 additinal points
     if ((n_ins.lt.nY*0.05).or.(n_ins.lt.2)) fact = 0.0
     do iY = 2, nY
        tmp1 = dabs(fbolom(iY-1)-fbolom(iY))
        if  ((TAUtot(1)*(ETAzp(1,iY)-ETAzp(1,iY-1)).GT.delTAUmax).or.&
             (tmp1.ge.fact*tmp2).or.&
             ((iterfbol.lt.3).and.&
             ((((fbol_em(iY)/fbolom(iY)).gt.0.01).and.&
             ((fbol_em(iY)/fbolom(iY)).lt.0.99)).or. &
             (((fbol_sc(iY)/fbolom(iY)).gt.0.01).and.&
             ((fbol_sc(iY)/fbolom(iY)).lt.0.99))))) then
           !print*,iY,tmp1,tmp2,fact*tmp2
           n_ins = 1
           kins = kins + n_ins
           do i_ins = 1,n_ins
              Yins(kins-n_ins+i_ins) = Y(iY-1)+1.*i_ins/(1.*(n_ins+1))*(Y(iY)-Y(iY-1))
              iYins(kins-n_ins+i_ins) = iY-1
           end do
        endif
     enddo
  endif
  do iY=1,nY
     Yprev(iY) = Y(iY)
     do iG=1,nG
        Td_old(iG,iY) = Td(iG,iY)
     end do
     do iL=1,nL
        utot_old(iL,iY) = utot(iL,iY)
     end do
  end do
  IF (kins.eq.0) THEN
     fbolOK = 1
  ELSE
     ! add all new points to Y(nY). this gives the new Y(nY+kins).
     ! however, check if npY is large enough to insert all points:
     if ((nY+kins).gt.npY) then
        fbolOK = 1
        if (iX.ge.1) then
           write(18,*)' ****************     WARNING   ******************'
           write(18,*)'  The new Y-grid can not accomodate more points!'
           write(18,'(a,i5)')'   The specified accuracy would require',nY+kins
           write(18,'(a,i5,a)')'   points, while npY =',npY,'.'
           write(18,*)'  For the required accuracy npY must be increased,'
           write(18,*)'  (see the manual s3.5 numerical accuracy).'
           write(18,'(a37,F5.1,a2)')'   The currently achieved accuracy is ', maxrat*100.0, ' %'
           write(18,*)' *************************************************'
        end if
        !! kins = npY - nY     !!this is in the old code, but doesn't work here.
        !! error = 2           !!this is in the old code, but doesn't work here.
        go to 777
     else
        do k = 1, kins
           do j = nY+k-1+1, iYins(k)+k-1+2, -1
              Y(j) = Y(j-1)
           end do
           Y(iYins(k)+k-1+1) = Yins(k)
           do iG=1,nG
              temp_mean = (0.5*(Td(iG,iYins(k)+k-1)**4. + Td(iG,iYins(k)+k)**4.))**(1./4.)
              !temp_mean = max(Td(iG,iYins(k)+k-1) ,Td(iG,iYins(k)+k))
              !shiftIns(x,Nmax,n,xins,i)
              do j = nY+k-1+1, iYins(k)+k-1+2, -1
                 Td(iG,j) = Td(iG,j-1)
              end do
              Td(iG,iYins(k)+k-1+1) = temp_mean
              ! call shiftIns(Td(iG,:),npY,nY+k-1,temp_mean,iYins(k)+k-1)
           end do
           do iL=1,nL
              temp_mean = 0.5*(utot(iL,iYins(k)+k-1)+utot(iL,iYins(k)+k))
              !temp_mean = max(utot(iL,iYins(k)+k-1),utot(iL,iYins(k)+k))
              !shiftIns(x,Nmax,n,xins,i)
              do j = nY+k-1+1, iYins(k)+k-1+2, -1
                 utot(iL,j) = utot(iL,j-1)
              end do
              utot(iL,iYins(k)+k-1+1) = temp_mean
              ! call shiftIns(utot(iL,:),npY,nY+k-1,temp_mean,iYins(k)+k-1)
           end do
        end do
     end if
  END IF
  ! new size of the y grid
  nY = nY + kins
  ! intepolate etadiscr to new y grid for denstyp = 5 or 6
  if(sph) then
     do iY = 1, nY
        Yloc = Y(iY)
        if (itereta.gt.1) then
           call lininter(nY,nYprev,Yprev,etatemp,Yloc,idm,ee)
           etadiscr(iY) = ee
        else
           etadiscr(iY) = eta(Yloc,nY,nYprev,itereta)
        end if
     end do
  end if
!--------------------------------------------------------------------------
777 return
end subroutine Flux_Consv
!*********************************************************************
!!$
!*********************************************************************
subroutine SLBdiff(nY,flag,grid,T4_ext,em,fp,fm)
!=========================================================================
! Integration of U(lam,t)*E2|Tau-t| to get the diffuse flux.
! flag=1 is for scattered and flag=0 for emitted flux. [MN, Apr'98]
!=========================================================================
  use common
  implicit none
  !---parameter
  integer :: nY,flag
  double precision,allocatable :: grid(:,:),T4_ext(:),em(:,:,:),&
       fp(:,:),fm(:,:)
  !---local
  integer iL, iY, j, iG
  double precision frac, efact, fave, sum, eint3
  double precision, allocatable ::  tau(:),faux(:)
  external eint3
  allocate(tau(nY))
  allocate(faux(nY))
  !--------------------------------------------------------------------
  do iL = 1, nL
     !$OMP PARALLEL DO private(iY,iG,frac) firstprivate(iL)
     do iY = 1, nY
        tau(iY) = grid(iL,iY)
        faux(iY) = 0
        do iG=1,nG
           frac = destroyed(iG,iY)*(sigmaA(iG,iL)+sigmaS(iG,iL))/(sigmaA(nG+1,iL)+sigmaS(nG+1,iL))
           if (flag.eq.1) then
              faux(iY) = faux(iY) + frac*omega(iG,iL)*utot(iL,iY)
           else
              faux(iY) = faux(iY) + frac*(1.0d0-omega(iG,iL))*em(iG,iL,iY)
           end if
        end do
     end do
     !$OMP END PARALLEL DO
     ! find f(+) (in arg tau>t)
     !$OMP PARALLEL DO private(iY,j,efact,fave,sum) firstprivate(iL)
     do iY = 1, nY
        sum = 0.0d0
        do j = 1, iY-1
           efact = abs(eint3(tau(iY)-tau(j))-eint3(tau(iY)-tau(j+1)))
           fave = 0.5d0*(faux(j)+faux(j+1))
           sum = sum + fave*efact
        end do
        fp(iL,iY) = 2.0d0*pi*sum
        ! and f(-) (in arg tau<t)
        sum = 0.0d0
        do j = iY, nY-1
           efact = abs(eint3(tau(iY)-tau(j))-eint3(tau(iY)-tau(j+1)))
           fave = 0.5d0*(faux(j)+faux(j+1))
           sum = sum + fave*efact
        end do
        fm(iL,iY) = 2.0d0*pi*sum
     end do
     !$OMP END PARALLEL DO
     do iY = 1, nY
        if (flag.eq.1) then
           fds(iL,iY) = (fp(iL,iY) - fm(iL,iY))
        else 
           fde(iL,iY) = (fp(iL,iY) - fm(iL,iY))
        end if
     end do
     ! end of loop over iL
  end do
  deallocate(tau)
  deallocate(faux)
  return
  !--------------------------------------------------------------------
end subroutine SLBdiff
!*********************************************************************

!**********************************************************************
subroutine add2(nY,flxs,flxe,fbsum)
!======================================================================
! This subroutine is auxiliary for finding the bolometric
! diffuse flux.   [MN, May'99]
!======================================================================
  use common
  implicit none
  INTERFACE
     subroutine Bolom(q,qbol,nY)
       integer nY
       double precision, allocatable :: q(:,:), qbol(:)
     end subroutine Bolom
  END INTERFACE
  !---parameter
  integer :: nY
  double precision, allocatable :: flxs(:,:),flxe(:,:),fbsum(:)
  !---local
  integer iY
  double precision, allocatable :: flxsb(:),flxeb(:)
  !-------------------------------------------------------------------
  allocate(flxsb(nY))
  allocate(flxeb(nY))
  call bolom(flxs,flxsb,nY)
  call bolom(flxe,flxeb,nY)
  do iY = 1, nY
     fbsum(iY) = flxsb(iY) + flxeb(iY)
  end do
  !-------------------------------------------------------------------
  deallocate(flxsb)
  deallocate(flxeb)
  return
end subroutine add2
!**********************************************************************

!***********************************************************************
SUBROUTINE SPH_Int(nY,nP,fs)
!***********************************************************************
!!$! This is the former  SUBROUTINE FindInt(nG,ETAzp) [MN].
!!$! This subroutine finds the intensity distribution at outer edge and for
!!$! user specified wavelengths lamOut. It also evaluates the angular size
!!$! of the stellar disk and adds two impact parameters describing the star
!!$! to the P grid, thus producing bOut grid. All intensities are indeed
!!$! dimensionless quantities lambda*I_lambda/F1 where I_lambda is real
!!$! physical quantity defined as usual and F1 is the bolometric flux at
!!$! the dust sublimation radius, r1. For conversion to the physical value
!!$! lambda*I_lambda, I_lambda from the program has to be multiplied by F1.
!!$! F1 can obtained either as:
!!$!      1) F1 = 4*sigma*Tsub**4/Psi (IE96, eq. 15),
!!$! where Tsub is sublimation temperature and parameter Psi is given in
!!$! *.spp file; or as:
!!$!      2) F1 = Fbol/alpha1**2 (IE96, eq. 34)
!!$! where Fbol is the bolometric flux and alpha1 is the angular size of r1
!!$! at any particular distance from the envelope (i.e. both F1 and alpha1
!!$! correspond to observed quantities). Also note that
!!$!     INT(I_lambda(p)*2Pi*P*dP) = f_lambda
!!$! where I_lambda is the scaled quantity from the program, P is impact
!!$! parameter, and f_lambda is the spectral shape F_lambda/Fbol.
!!$!                                                      [Z.I., Aug. 1996]
!!$! =======================================================================
  use common
  implicit none
  INTERFACE
     subroutine Simpson(n,n1,n2,x,y,integral)
       integer n, n1, n2
       double precision integral
       double precision,allocatable ::  x(:), y(:)
     end subroutine Simpson
  END INTERFACE
  !---parameter
  integer nY,nP
  double precision,allocatable :: fs(:,:)
  !---local
  integer :: i,iL,iP,iY,iZ, iW,k, nZ, Nzpt,izloc,iLout,iLstop
  double precision :: xP,Planck, resaux, QUtot1, QpTsub,res,denum, &
       pst, numcorr, stelfact,w1, w2, ep1, z1, z2, delz, delTau, &
       lw12, pT, palb, UtotL, UtotR, pUtot, palf, tauzp1, zloc, &
       wloc,Tz, IntETA, alb, Utotloc, factaux, alfa, tauInf, &
       exterm, resint, ETAzpStar, Idboth, Idfront, xx, IntL, IntR
  double precision,allocatable :: qaux(:),qaux2(:),alpha(:,:),&
       fnum(:),fdenum(:),Istell(:),tauOut(:), Ids(:,:), Ide(:,:), &
       Istsc(:,:),Istem(:,:),tzp(:),Semis(:),Sscat(:),Sstem(:),&
       Sstsc(:)
  external Planck,IntETA
!!$  INTEGER iL,nG,iY,k,i,iLout,iLstop, iP, iW, nZ, Nzpt, iZ, izloc
!!$  DOUBLE PRECISION qaux(npL), alpha(1,npY), resaux, QUtot1,   &
!!$       QpTsub, xP, pst, stelfact, Istell(npL), Planck, z1, z2,    &
!!$       IntL, IntR, xx , alb, Ids(npL,npP), Ide(npL,npP),w1, w2,lw12, &
!!$       numcorr, ETAzpStar, qaux2(npL), fs(npL,npY), omega(npG,npL), &
!!$       delz, zloc, wloc, resint, pT, Tz, Idboth, tzp(100), pUtot, &
!!$       Semis(100), Sscat(100), IntETA, palb, palf, alfa, exterm,  &
!!$       Utotloc, Sstem(100), Sstsc(100), Istem(npL,100), Idfront,  &
!!$       Istsc(npL,100), delTau, factaux, UtotL, UtotR, ep1,        &
!!$       tauzp1, tauInf, fnum(npL), fdenum(npL), res, denum
!!$  EXTERNAL IntETA
!!$! -----------------------------------------------------------------------
  ! temporary
  allocate(qaux(nL))
  allocate(qaux2(nL))
  allocate(alpha(1,nY))
  allocate(fnum(nL))
  allocate(fdenum(nL))
  allocate(Istell(nL))
  allocate(tauOut(nL))
  allocate(Ids(nL,nP))
  allocate(Ide(nL,nP))
  allocate(Istsc(nL,nP))
  allocate(Istem(nL,nP))
  allocate(tzp(nP))
  allocate(Semis(nP))
  allocate(Sscat(nP))
  allocate(Sstem(nP))
  allocate(Sstsc(nP))
  IF (nG.GT.1.AND.iX.GE.1) THEN
     print*, ' FindInt should be fixed, nG>1 ! kernel.f90 LINE 3141'
     stop
  END IF
  ! find impact parameter tangential to the stellar disk
  ! first find the Planck averaged absorption efficiencies at Y=1
  ! print*,'still single grain line 3529'
  ! stop
  DO iL = 1, nL
     qaux(iL) = SigmaA(1,iL) * Utot(iL,1) / lambda (iL)
     xP = 14400.0 / Td(1,1) / lambda(iL)
     qaux2(iL) = SigmaA(1,iL) * Planck(xP) / lambda (iL)
  END DO
  CALL Simpson(nL,1,nL,lambda,qaux,resaux)
  QUtot1 = resaux
  CALL Simpson(npL,1,nL,lambda,qaux2,resaux)
  QpTsub = resaux
  ! parameter Psi (see Ivezic & Elitzur, 1996, eq. C4)
  Psi = QUtot1 / QpTsub
  alpha(1,1) = Psi
  ! ***alpha is a local array, to match the changes with Zeljko's old expressions
  DO iY = 2, nY
     ! calculate f1 and f2
     DO iL = 1, nL
        fnum(iL) = SigmaA(1,iL) * Utot(iL,iY) / lambda(iL)
     END DO
     CALL Simpson(nL,1,nL,lambda,fnum,res)
     ! calculate alpha
     DO iL = 1, nL
        xP = 14400.0 / lambda(iL) / Td(1,iY)
        fdenum(iL) = SigmaA(1,iL) * Planck(xP) / lambda(iL)
     END DO
     CALL Simpson(npL,1,nL,lambda,fdenum,denum)
     alpha(1,iY) = res /  denum
  END DO
!!$  ! *********
  ! ratio pst = rstar/rsub (see Ivezic & Elitzur, 1996, eq. 27)
  ! Added Apr.07 [MN] to take care of the case of no central source
  IF (Left.eq.0) THEN
     pst = 1.
  ELSE
     pst = 2.0 / dsqrt(Psi) * (Td(1,1) / Tstar(1))**2.0
  END IF
  ! this is if only central source is present
  IF (pst.GE.0.5.AND.Right.eq.0) THEN
     IF (iX.GE.1) THEN
        write(18,*)' FindInt: specified dust temperature at the '
        write(18,*)' inner radius results in r*/r1 >= 0.5: '
        write(18,*)'    r*/r1 =', pst
        write(18,*)' This violates some of Dusty`s assumptions'
        write(18,*)'  ------  Please consult the manual ------'
        write(18,*)'  ####  r*/r1 changed by hand to 0.5  ####'
     END IF
     pst = 0.5
  END IF
  stelfact = 1.0 / pst / pst / Pi
  ! generate bOut, i.e. insert two points such that
  ! bOut(k)=0.999*pst and bOut(k+1)=1.001*pst
  CALL GetbOut(nP,pst,k)
  ! correction for numerical errors in tau
  numcorr = 1. / TAUtot(1)
  ! loop over wavelengths
  DO iL = 1, nL
     ! stellar intensity, Istell (extinction already included)
     Istell(iL) = fs(iL,nY) * stelfact
     ! total optical depth along a line of sight
     tauOut(iL) = numcorr*TAUtot(iL)
  END DO
  ! generate diffuse intensities, Ide (emission) and Ids (scat)
  ! loop over wavelengths
  DO iL = 1, nL
     DO iP = 1, nP
        ! maximal number of points along tangential position, z
        nZ = nY + 1 - iYfirst(iP)
        ! starting value for local radius
        IF (P(iP).GE.1.0) THEN
           w2 = P(iP)
        ELSE
           w2 = 1.0
        END IF
        ! initialize intensities
        Ide(iL,iP) = 0.0
        Ids(iL,iP) = 0.0
        IF (iP.LE.k+1) THEN
           Istem(iL,iP) = 0.0
           Istsc(iL,iP) = 0.0
        END IF
        ! total optical depth along this impact parameter
        ep1 = ETAzp(iP,nZ)*TAUtot(iL)
        ! loop over z, i.e. steps over points crossing the y grid
        DO iZ = 2, nZ
           ! index for the ending local radius
           iW = iYfirst(iP) + iZ - 1
           ! local boundary radii
           w1 = w2
           w2 = Y(iW)
           ! corresponding displacements along a line of sight
           z1 = sqrt(abs(w1**2.-P(iP)**2.))
           z2 = sqrt(abs(w2**2.-P(iP)**2.))
           ! # of pts. for z integration, should increase with deltaTau
           ! it is messy because INT function which would do the job is
           ! not in F77 standard set
           Nzpt = 5
           delTau = (ETAzp(iP,iW)-ETAzp(iP,iW-1))*TAUtot(iL)
           IF (delTau.GT.1) Nzpt = 10
           IF (delTau.GT.5) Nzpt = 20
           IF (delTau.GT.10) Nzpt = 30
           IF (delTau.GT.20) Nzpt = 40
           IF (delTau.GT.50) Nzpt = 50
           delz = (z2-z1) / (Nzpt-1)
           ! powers for power-law interpolations between 2 y pts.
           lw12 = dlog(Y(iW-1)/Y(iW))
           ! for T
           pT = dlog(Td(1,iW)/Td(1,iW-1)) / lw12
           ! for albedo
           ! stop
           IF (omega(nG+1,iL).GT.0.0.AND.omega(nG+1,iL).GT.0.0) THEN
              palb = dlog(omega(nG+1,iL)/omega(nG+1,iL)) / lw12
           ELSE
              palb = 0.0
           END IF
           ! for Utot
           UtotL = Utot(iL,iW-1)
           UtotR = Utot(iL,iW)
           CALL ChkRange(dynrange,UtotL)
           CALL ChkRange(dynrange,UtotR)
           IF (UtotL.GT.0.0.AND.UtotR.GT.0) THEN
              pUtot = dlog(UtotR/UtotL) / lw12
           ELSE
              pUtot = 0.0
           END IF
           ! for alpha
           palf = dlog(alpha(1,iW)/alpha(1,iW-1)) / lw12
           ! tauzp between z=0 and z=z1
           tauzp1 = ETAzp(iP,iZ-1)*TAUtot(iL)
           ! integrate between adjacent grid points
           DO izloc = 1, Nzpt
              zloc = z1 + (izloc-1)*delz
              wloc = sqrt(zloc**2 + P(iP)**2)
              ! find local TAUzp(w(z))-TAUzp(w1=w(z1))
              tzp(izloc) = IntETA(P(iP),iW-1,w1,wloc)*TAUtot(iL)
              ! find Tz = T(zloc) = T(wloc), this works for single
              ! size grains only; for multigrain case one needs to
              ! get Semis by summation over all Td
              Tz = Td(1,iW-1) * (Y(iW-1)/wloc)**pT
              xP = 14400/lambda(iL)/Tz
              ! power-law interpolation for albedo
              alb = omega(nG+1,iL) * (Y(iW-1)/wloc)**palb
              ! power-law interpolation for Utot
              IF (UtotL.GT.0) THEN
                 UtotLoc = UtotL * (Y(iW-1)/wloc)**pUtot
              ELSE
                 UtotLoc = 0.0
              END IF
              CALL ChkRange(dynrange,UtotLoc)
              ! power-law interpolation for alpha
              alfa = alpha(1,iW-1) * (Y(iW-1)/wloc)**palf
              ! source functions (wloc**2 because D uses scaled quant.)
              factaux = 1 / wloc**2 / (4. * Pi)
              Semis(izloc) = (1-alb) * alfa * Planck(xP) * factaux
              Sscat(izloc) = alb * UtotLoc * factaux
              ! check for the dynamic range
              CALL ChkRange(dynrange,Semis(izloc))
              CALL ChkRange(dynrange,Sscat(izloc))
              ! optical depth from infinity along the line of sight
              tauInf = ep1 - tauzp1 - tzp(izloc)
              ! for a line of sight terminating on the star find
              ! contribution only from the front part of the envelope
              IF (iP.LE.k+1) THEN
                 IF (tauInf.LT.50) THEN
                    exterm = dexp(-tauInf)
                 ELSE
                    exterm = 0.0
                 END IF
                 Sstem(izloc) = Semis(izloc) * exterm
                 Sstsc(izloc) = Sscat(izloc) * exterm
              END IF
              ! otherwise take both the front and back contributions
              IF (tauInf.LT.50) THEN
                 exterm = dexp(-tauInf)+dexp(-tauInf-ep1)
              ELSE
                 exterm = 0.0
              END IF
              Semis(izloc) = Semis(izloc) * exterm
              Sscat(izloc) = Sscat(izloc) * exterm
              ! end of local loop over z
           END DO
           ! integrate and add contribution from this step
           CALL SIMPSON(100,1,Nzpt,tzp,Semis,resint)
           CALL ChkRange(dynrange,resint)
           Ide(iL,iP) = Ide(iL,iP) + resint
           CALL SIMPSON(100,1,Nzpt,tzp,Sscat,resint)
           CALL ChkRange(dynrange,resint)
           Ids(iL,iP) = Ids(iL,iP) + resint
           IF (iP.LE.k+1) THEN
              CALL SIMPSON(100,1,Nzpt,tzp,Sstem,resint)
              CALL ChkRange(dynrange,resint)
              Istem(iL,iP) = Istem(iL,iP) + resint
              CALL SIMPSON(100,1,Nzpt,tzp,Sstsc,resint)
              CALL ChkRange(dynrange,resint)
              Istsc(iL,iP) = Istsc(iL,iP) + resint
           END IF
           ! end of loop over z
        END DO
        ! end of loop over impact parameter, iP
     END DO
     ! end of loop over wavelengths, iL
  END DO
  ! add all intensities, Istell, Ide, Ids
  DO iL = 1, nL
     ! interpolate optical depth  at pstar
     IF (iL.EQ.iLfid) THEN
        ETAzpStar = (ETAzp(k,nY) - ETAzp(k-1,nY))
        ETAzpStar = ETAzpStar * (pst-P(k-1)) / (P(k)-P(k-1))
        ETAzpStar = ETAzp(k-1,nY) + ETAzpStar
     END IF
     ! find diffuse contribution at pstar (by linear interpolation)
     Idfront = Istsc(iL,k)+Istem(iL,k)-Istsc(iL,k-1)-Istem(iL,k-1)
     Idfront = Idfront * (pst-P(k-1)) / (P(k) - P(k-1))
     Idfront = Idfront + Istsc(iL,k-1) + Istem(iL,k-1)
     Idboth = Ids(iL,k) + Ide(iL,k) - Ids(iL,k-1) - Ide(iL,k-1)
     Idboth = Idboth * (pst-P(k-1)) / (P(k) - P(k-1))
     Idboth = Idboth + Ids(iL,k-1) + Ide(iL,k-1)
     ! first for p<pstar, all three contributions
     DO i = 1, k-1
        Intens(iL,i) = Istell(iL) + Istsc(iL,i) + Istem(iL,i)
        IF (iL.EQ.iLfid) tauZout(i) = ETAzp(i,nY)/ETAzp(1,nY)
     END DO
     ! barely on the stellar disk
     Intens(iL,k) = Istell(iL) + Idfront
     tauZout(k) = ETAzpStar/ETAzp(1,nY)
     ! barely off the stellar disk
     Intens(iL,k+1) = Idboth
     tauZout(k+1) = 2. * tauZout(k)
     ! all other p>pstar
     DO i = k, nP
        Intens(iL,i+2) = Ids(iL,i)+Ide(iL,i)
        IF (iL.EQ.iLfid) THEN
           nZ = nY + 1 - iYfirst(i)
           tauZout(i+2) = 2. * ETAzp(i,nZ)/ETAzp(1,nY)
        END IF
     END DO
  END DO
  ! check dynamic range
  DO iL = 1, nL
     DO i = 1, nP+2
        CALL ChkRange(dynrange,Intens(iL,i))
     END DO
  END DO
  ! now interpolate Intens(lambda) to lamOut
  DO iLout = 1, NlambdaOut
     ! bracket the needed wavelength
     iLstop = 0
     iL = 0
     DO WHILE (iLstop.EQ.0)
        iL = iL + 1
        IF (lambda(iL).GT.LambdaOut(iLout)) iLstop = 1
        IF (iL.EQ.nL) iLstop = 1
     END DO
     ! interpolate intensity
     xx = (LambdaOut(iLout)-lambda(iL-1))/(lambda(iL)-lambda(iL-1))
     DO i = 1, nP+2
        IntL = Intens(iL-1,i)
        IntR = Intens(iL,i)
        IntOut(iLout,i) = IntL + xx*(IntR - IntL)
        CALL ChkRange(dynrange,IntOut(iLout,i))
     END DO
  END DO
  ! -------------------------------------------------------------------
999 deallocate(qaux)
  deallocate(qaux2)
  deallocate(alpha)
  deallocate(fnum)
  deallocate(fdenum)
  deallocate(Istell)
  deallocate(tauOut)
  deallocate(Ids)
  deallocate(Ide)
  deallocate(Istsc)
  deallocate(Istem)
  deallocate(tzp)
  deallocate(Semis)
  deallocate(Sscat)
  deallocate(Sstem)
  deallocate(Sstsc)
  RETURN
end subroutine SPH_Int
! ***********************************************************************

! ***********************************************************************
  subroutine SLBintensity(nY,em)
! =======================================================================
  use common
  implicit none
  !--- parameter
  integer :: nY
  double precision, allocatable :: em(:,:,:)
  !--- local
  integer :: iL, iY, iG, imu
  double precision :: kron, idifp, idifm,res, Sexp, frac
  double precision, allocatable :: tau1(:)
  external Sexp
  ! --------------------------------------------------------------------
  allocate(tau1(nY))
  ! Loop over wavelengths
  do iL = 1, nL
     taut = TAUslb(iL,nY)
     ! loop over angles (i.e. muobs), read in input from 'slb_mugrid.dat'
     do imu = 1, nmu
        muobs = dcos(theta(imu))
        if(dabs(mu1-muobs).lt.1.0d-4) then
           kron = 1.0d0
        else
           kron = 0.0d0
        end if
        ! scaled diffuse intensity:
        idifm = 0.0d0
        idifp = 0.0d0
        tau1(1) = 0.0d0
        if ((sigmaA(nG+1,iL)+sigmaS(nG+1,iL)).gt.0.0) then 
           do iY = 2, nY
              tau1(iY) = TAUslb(iL,iY)
              Sfn = 0
              do iG=1,nG
                 frac = (sigmaA(iG,iL)+sigmaS(iG,iL))/(sigmaA(nG+1,iL)+sigmaS(nG+1,iL))
                 Sfn = Sfn + frac*(1.0d0-omega(iG,iL))*em(iG,iL,iY)
                 Sfn = Sfn + frac*omega(iG,iL)*utot(iL,iY)
              end do
              if (Sfn.gt.0.0) then 
                 ! transmit=1 for tau < t, transmit=0 for tau > t
                 transmit = 1
                 call romby(Sexp,tau1(iY-1),tau1(iY),res)
                 idifm = idifm + res
                 transmit = 0
                 call romby(Sexp,tau1(iY-1),tau1(iY),res)
                 idifp = idifp + res
              end if
           end do
        end if
        if(idifm.lt.1.d-20) idifm = 0.0d0
        if(idifp.lt.1.d-20) idifp = 0.0d0
        SLBintm(imu,iL) = idifm
        SLBintp(imu,iL) = idifp
        ! enddo over angles
     end do
     ! enddo over lambda
  end do
  deallocate(tau1)
! -----------------------------------------------------------------
  return
  end subroutine SLBintensity
! ***********************************************************************

!***********************************************************************
 subroutine Analysis(nY,model,us,T4_ext,delta,maxrat)
!=======================================================================
! This subroutine analyzes the solution. It finds the flux conservation
! accuracy and evaluates many output quantites (e.g. QF(y), TAUF(y),Psi, F1
! the rad.pressure force, dynamical quantities etc.)
! This is with new additions acc. to IE'00           [ZI,Mar'96;MN,Mar'99]
!=======================================================================
  use common
  implicit none
  INTERFACE
     subroutine Simpson(n,n1,n2,x,y,integral)
       integer n, n1, n2
       double precision integral
       double precision,allocatable ::  x(:), y(:)
     end subroutine Simpson
     SUBROUTINE FindErr(nY,flux,maxFerr)
       integer nY
       DOUBLE PRECISION maxFerr
       double precision,allocatable :: flux(:)
     END SUBROUTINE FindErr
  END INTERFACE
  !---parameter
  integer nY,model
  double precision :: delta,maxrat
  double precision, allocatable ::  us(:,:),T4_ext(:)
  !---local
  integer :: i,iL,iY,iG
  double precision eta, maxFerr, resaux, aux, s4, Planck, xp, qutot1, Eps1, &
       mx, theta1_loc, sig_22, tauV, Tc3, C1, C2, C3, ugas_out, Gie2000
  double precision, allocatable :: spectrum(:),qaux(:), K1(:), K2(:),&
       qpTd(:,:), qpstar(:), qaux2(:)
  external eta,Planck
  !-------------------------------------------------------------------

!!$  integer i, iL, iY, nn, model, iP, error, iG, nG
!!$  double precision eta, qpTd(npG,npY), qpstar(npY), &
!!$       qaux(npL), qaux2(npL), resaux, xP, Planck, qutot1,       &
!!$       eps1, aux, C1, C2, C3, theta1_loc, ugas_out, s4, mx,     &
!!$       tauV, Gie2000, K1(npY), K2(npY), tauaux(npL,npY),       &
!!$       delta, us(npL,npY), x1, x2, result1, T4_ext(npY), maxFerr, &
!!$       L4, Mo, Tc3, sig_22
!!$  external eta
  allocate(spectrum(nL))
  allocate(qaux(nL))
  allocate(qaux2(nL))
  allocate(K1(nY))
  allocate(K2(nY))
  allocate(qpTd(nG,nY))
  allocate(qpstar(nY))
  ! spectrum (flux at the outer edge as a function of wavelength)
  do iL = 1, nL
     spectrum(iL) = dabs(ftot(iL,nY))
     ! to prevent taking log from zero in spectral [MN]:
     if (spectrum(iL).le.1.0d-20) spectrum(iL) = 1.0d-20
  end do
  !-------------
  ! analyze bolometric flux error (1/2 of the max spread of fbol)
  CALL FindErr(nY,fbol,maxFerr)
  ! find the flux averaged optical depth, tauF(y)
  if (sph) then
     ! for spherical shell
     tauF(1) = 0.0
     DO iY = 2, nY
        ! generate auxiliary function for integration:
        ! loop over iL (wavelength)
        ! N.B. the definition: ETAzp(1,y) = taur(y)/tauT so that
        ! tau(iL,iY) = TAUtot(iL)*ETAzp(1,iY)
        DO iL = 1, nL
           qaux(iL)=TAUtot(iL)*ETAzp(1,iY)*dabs(ftot(iL,iY))/lambda(iL)
        END DO
        CALL Simpson(nL,1,nL,lambda,qaux,resaux)
        ! tauF(iY) = <tau(iL,iY)*ftot(iL,iY)>
        tauF(iY) = resaux
     END DO
     ! for full RDW calculation redo tauF to be consistent with CalcEta
     IF (denstyp.eq.3) THEN !3(RDW)
        ! generate ETA and its integral (normalization constant)
        DO iY = 1, nY
           K1(iY) = vrat(1,iY)/ugas(iY)/Y(iY)/Y(iY)
        END DO
        CALL SIMPSON(nY,1,nY,Y,K1,resaux)
        ! find tauF
        DO iY = 1, nY
           K2(iY) = qF(iY)*K1(iY)/resaux
           CALL SIMPSON(nY,1,iY,Y,K2,aux)
           tauF(iY) = TAUfid*aux
        END DO
     END IF
  elseif(slb) then
     ! for slab
     tauF(1) = 0.0d0
     do iY = 1, nY
        ! generate auxiliary function for integration:
        do iL = 1, nL
           qaux(iL) = TAUslb(iL,iY)*dabs(ftot(iL,iY))/(fsLbol(1)*lambda(iL))
           call Simpson(npL,1,nL,lambda,qaux,resaux)
           tauF(iY) = resaux
        end do
     end do
  end if
  ! ------------
  ! ratio of gravitational to radiation pressure force (isotropic scattering) per unit volume
  ! s4 = (L4sol/Msol)/(4*Pi*G*c*rho_s)/1e-6;
  ! rho_s=3000 kg.m-3, grain radius 'a' is in microns, aveV=4/3*Pi*<a^3>
  IF(sph) THEN
     s4 = 1.925 / (4.0*Pi*Gconst*3.0d08*3000.0*1.0D-06)
     ! in case of sigma's from a file aveV=1 (initialized in GetOptPr)
     DO iY = 1, nY
        do iG = 1,nG+1
           DO iL = 1, nL
              qaux(iL)=(SigmaA(iG,iL)+SigmaS(iG,iL))/aveV*dabs(ftot(iL,iY))/lambda(iL)
           END DO
           CALL Simpson(nL,1,nL,lambda,qaux,resaux)
           rg(iG,iY) = s4 * resaux / r_gd
           ! If dust drift (dynamics case):
           IF (denstyp.eq.3) rg(iG,iY) = rg(iG,iY)*vrat(1,iY) !3(RDW)
           IF ((iY.EQ.1).AND.(iG.eq.nG+1)) THEN
              Phi = resaux
           END IF
        END DO
     END DO
     ! the terminal value of the reddening profile, normalized to y=1
     Phi = resaux / Phi
  END IF
  !-------------
  ! Find the Planck averaged absorption efficiencies
  do iG = 1,nG
     DO iY = 1, nY
        ! generate auxiliary function for integration over wavelengths:
        DO iL = 1, nL
           qaux(iL) = SigmaA(iG,iL) * Us(iL,iY) / lambda(iL)
           xP = 14400.0 / Td(iG,iY) / lambda(iL)
           qaux2(iL) = SigmaA(iG,iL) * Planck(xP) / lambda (iL)
        END DO
        CALL Simpson(nL,1,nL,lambda,qaux,resaux)
        QpStar(iY) = resaux
        CALL Simpson(nL,1,nL,lambda,qaux2,resaux)
        QpTd(iG,iY) = resaux
     END DO
  end do
  ! ----------
  ! find parameter Psi (see Ivezic & Elitzur, 1996)
  ! generate auxiliary function for integration:
  ! loop over iL (wavelength)
  DO iL = 1, nL
     qaux(iL) = SigmaA(nG+1,iL) * Utot(iL,1) / lambda (iL)
  END DO
  CALL Simpson(nL,1,nL,lambda,qaux,resaux)
  QUtot1 = resaux
  Psi = QUtot1 / QpTd(1,1)
  !!**  added Psi0 from eq.(41) of IE'01 [MN]
  Psi0 = QpStar(1) / QpTd(1,1)
  ! for slab Psi is defined by the flux at normal ill.
  IF (SLB) Psi = dabs(mu1)*QUtot1 / QpTd(1,1)
  ! -------------
  IF(sph) THEN
     ! ratio r1/r* (see Ivezic & Elitzur, 1996, eq. 27)
     r1rs = 0.5 * dsqrt(Psi) * (Tstar(1) / Td(1,1))**2.0
     IF(Left.eq.0) r1rs = 1.0
  END IF
  ! -------------
  ! Find epsilon - the relative contribution of the diffuse radiation
  DO iY = 1, nY
     aux = QpStar(iY)/QpTd(1,iY)/Psi*(Td(1,1)/Td(1,iY))**4.
     IF (SLB) THEN
        aux = aux*dabs(mu1)
     ELSE
        aux = aux/ Y(iY)/Y(iY)
     END IF
     Eps(iY) = 1. - aux
  END DO
  Eps1 = 1.0 - QpStar(1) / QUtot1
  ! store these parameters in the storage array
  SmC(1,model) = Psi
  SmC(2,model) = Eps1
  SmC(3,model) = QpStar(1)
  SmC(4,model) = QpTd(1,1)
  !maxrat is the max err of calculated diffuse flux as in Blueprint.
  SmC(5,model) = maxrat
  ! inner radius (in cm) in case it is not an input
  ! 5.53e16 = sqrt(10^4*Lo/4/pi)
  if (typentry(1).ne.2) then
     if(slb) then
        Cr1 =  5.53d+16 / dsqrt(Ji*4*pi/abs(mu1))
     else
        Cr1 = 5.53d+16 / dsqrt(Ji*4*pi)
     end if
  end if
  if (sph) then
     ! angular diameter of inner cavity if Fbol=1d-6 W/m2
     theta1_loc = 412.6d0 / dsqrt(Ji*4*pi)
     ! check if the pt.source assumption is still obeyed
     ! (This is only for BB-type spectrum including Engelke-Marengo function)
     if(startyp(1).eq.1.or.startyp(1).eq.2) then
        mx = sqrt(sqrt(Ji*4*pi/sigma))
        Te_min = 2.0d0 * dmax1(Td(1,1), mx)
     end if
  end if
  ! calculate conversion constants for dynamics
  Tc3 = Td(1,1)/1000.0
  ! from eq.(5): 10^(-22)*sig_22 = 10^(-8)*aveA*nd/ng = aveA/r_gd; aveA is in micron^2
  ! Dusty uses SI units, so convert to m^2;
  ! rho_s = 3000 kg/m3, aveA/aveV = 1/a, where grain radius a is in microns
  sig_22 = 1.67d04 /r_gd/3000.0d00 * aveA/aveV   ! eq.(86) in IE'01
  QV = SigExfid / aveA
  tauV = TAUfid
  if (denstyp.eq.3) then !RDW
     ! ---- printout in fname.out according to IE2001 -------
     ! from eq.(43)
     CMdot = 1.0d-06 * 1.98*dsqrt(Qstar/sig_22) * Prdw*(Psi**0.25)/Tc3
     ! from eq.(44) where ugas(iY) = Qstar*w(iY) found in Winds
     Cve = (2.0d-05/3.) * ugas(nY)/CMdot
     ! from eq.(15) with Gamma = gmax
     CM = 45.8*Qstar*sig_22 / gmax
     write(39,*) '    '
     write(39,'(a7,1p,e12.3)') ' tauV/QV=', tauV/QV
     write(39,'(a7,1p,e12.3)') '    QV=', QV
     write(39,'(a7,1p,e12.3)') ' Qstar=', Qstar
  end if
  ! for analytical approximation ---------------
  if (denstyp.eq.4) then !(RDWA)
     I1_dyn = 2.0d0 * (1.0d0-pow)/(1.0d0+pow)/tauF(nY)
     I2_dyn = I1_dyn
     I3_dyn = I1_dyn * tauF(nY) / taufid
     Gamma(nY) = 0.5d0
     ! terminal expansion velocity, full formula:
     ugas_out = tauF(nY) * (1.-Gamma(nY)) / (1.-pow)
     ! The coefficients come from the units conversion
     C1 = 0.2845*TAUfid*sqrt(Psi)/I2_dyn/(SigExfid/aveV)/ Tc3
     C2 = 2.040*ugas_out
     C3 = 6.628*I3_dyn*SigExfid/aveV*Gamma(nY)/I1_dyn
     ! from version 2.0 stellar mass is defined as the maximal stellar
     ! mass which does not quench the wind; the calculation is done
     ! with half that mass since any smaller mass will have no effect
     ! on the radial velocity and density profile (see IE2001)
     CM = 6.628*I3_dyn*SigExfid/aveV/I1_dyn
     ! mass-loss rate in Msol/yr
     CMdot = 1.0E-05 * sqrt(C1)
     ! terminal expansion velocity in km/s
     Cve = 10.* C2 / sqrt(C1)
  end if
  ! *** this is conversion to the nomenclature as in IE2001
  IF (denstyp.eq.5) THEN !(RDWPR)
     ! size averaged extinction efficiency
      Qstar = qF(1)
      zeta1 = vrat(1,1)
      G1 = Gamma(1)
      Ginf = Gamma(nY)
      IF (G1.GT.0) THEN
          Gie2000 = 1 / zeta1 / G1
          delta = 1 / (Gie2000 - 1)
      ELSE
          delta = 0.0
      END IF
      PIrdw = tauV / QV
      Prdw = sqrt(2*PIrdw/I2_dyn/QV/Qstar)
      winf = ugas_out / QV / Qstar
  END IF
  deallocate(qpTd)
  deallocate(qpstar)
  deallocate(spectrum)
  deallocate(qaux)
  deallocate(qaux2)
  deallocate(K1)
  deallocate(K2)
  return
  !--------------------------------------------------------------------
end subroutine Analysis
!***********************************************************************
!!$
!!$! ***********************************************************************
!!$SUBROUTINE Visibili
!!$! =======================================================================
!!$! This subroutine finds visibility functions corresponding to IntOut.
!!$! The work horse is subroutine Visi2D, and this subroutine is used to
!!$! prepare everything.                                  [Z.I., Jan. 1997]
!!$! =======================================================================
!!$  use common
!!$  IMPLICIT none
!!$  INTEGER i, j, N1, N2
!!$  DOUBLE PRECISION  Visi(1000), Int1D(npP+2)
!!$! -----------------------------------------------------------------------
!!$  ! generate spatial frequency (q) grid
!!$  ! first N1 points up to qtheta1=1.22 (Rayleigh limit for a disk)
!!$  N1 = 80
!!$  ! first 2 points manually:
!!$  ! there must be 0!
!!$  qtheta1(1) = 0.0D+00
!!$  ! make sure the whole envelope is resolved
!!$  qtheta1(2) = 0.5D+00 / bOut(nP+2)
!!$  ! and the rest on logarithmic grid up to 1.22
!!$  DO i = 1, N1-2
!!$     qtheta1(i+2)=qtheta1(2)*(1.22D+00/qtheta1(2))**(i*1.0D+00/(N1-2))
!!$  END DO
!!$  ! envelope is well sampled, now to be sure that the star will be OK
!!$  ! for small taus add N2 points on a logarithmic grid up to 1.22/p*
!!$  N2 = 20
!!$  DO i = 1, N2
!!$     qtheta1(N1+i) = 1.22D+00 / bOut(2)**(i*1.0D+00/N2)
!!$  END DO
!!$  Nvisi = N1 + N2
!!$  ! find visibility wavelength by wavelength
!!$  DO j = 1, NlambdaOut
!!$     DO i = 1, nP+2
!!$        Int1D(i) = IntOut(j,i)
!!$        CALL ChkRange(dynrange,Int1D(i))
!!$        IF (Int1D(i).LT.dynrange) Int1D(i)=0.0D+00
!!$     END DO
!!$     CALL Visi2D(npP+2,nP+2,bOut,Int1D,1000,N1+N2,qtheta1,Visi)
!!$     ! copy 1D convolved visibility to Visib
!!$     DO i = 1, N1+N2
!!$        ! check dynamic range
!!$        CALL ChkRange(dynrange,Visi(i))
!!$        Visib(j,i) = Visi(i)
!!$     END DO
!!$  END DO
!!$! -----------------------------------------------------------------------
!!$  RETURN
!!$END SUBROUTINE Visibili
!!$! ***********************************************************************
!!$
!!$! ***********************************************************************
!!$SUBROUTINE Visi2D(NinMax,Nin,Xin,Yin,Noutmax,Nout,Xout,Yout_loc)
!!$! =======================================================================
!!$! This subroutine finds the visibility function (the spatial Fourier
!!$! transform of the intensity distribution) corresponding to the
!!$! intensity Yin(Xin[i]), i=1,Nin. Visibility, Yout, is evaluated at q
!!$! positions (spatial frequency) given in Xout[i], i=1,Nout. Maximum size
!!$! of Xin is NinMax, maximum size of Xout is NoutMax. The Bessel function
!!$! of the zeroth order is provided separately. The integration is done by
!!$! calling subroutine ROMBY (Bessel function is called from IMGFN).
!!$! Note:
!!$! The visibility function V(q) for a circularly symmetric intensity
!!$! I(x) is:
!!$!          V(q) = F(q)/F(0)
!!$! where Jo is the Bessel function of the zeroth order, and
!!$!          F(q) = Int[Jo(2Pi*q*x)*I(x)*2Pi*x*dx]
!!$! Note that F(0) is nothing more than flux. For more details see
!!$! Ivezic & Elitzur, 1996, MNRAS, 279, 1019 and ref. therein.
!!$!                                                      [Z.I., Jan. 1997]
!!$! =======================================================================
!!$  use common
!!$  IMPLICIT none
!!$  INTEGER NinMax, Nin, NoutMax, Nout, iq, iXin
!!$  DOUBLE PRECISION Xin(NinMax), Yin(NinMax), Xout(NoutMax),  &
!!$       Yout_loc(NoutMax),  F(1000), F0, int1, int2, A, B, imagfn
!!$  EXTERNAL imagfn
!!$! --------------------------------------------------------------------
!!$  ! loop over spatial frequency q (= Xout)
!!$  DO iq = 1, Nout
!!$     Cqtheta = Xout(iq)
!!$     F(iq) = 0.0D+00
!!$     ! loop over radial positions
!!$     !! DO iXin = 1, Nin - corrected after Zeljko's email from 7/29/09
!!$     DO iXin = 1, Nin-1
!!$        ! find F(q)
!!$        ftype = 2
!!$        Ckn = 1.0D+00
!!$        CALL ROMBY(imagfn,Xin(iXin),Xin(iXin+1),int1)
!!$        Ckn = 2.0D+00
!!$        CALL ROMBY(imagfn,Xin(iXin),Xin(iXin+1),int2)
!!$        ! contribution from this annulus (lin. approx. for intensity)
!!$        A = Xin(iXin+1)*Yin(iXin)-Xin(iXin)*Yin(iXin+1)
!!$        A = A /(Xin(iXin+1)-Xin(iXin))
!!$        B = (Yin(iXin+1)-Yin(iXin))/(Xin(iXin+1)-Xin(iXin))
!!$        F(iq) = F(iq) + A*int1 + B*int2
!!$     END DO
!!$  END DO
!!$  ! flux
!!$  F0 = F(1)
!!$  DO iq = 1, Nout
!!$     IF(F0.EQ.0.0D+00) THEN
!!$        Yout_loc(iq) = 0.0D+00
!!$     ELSE
!!$        Yout_loc(iq) = dabs(F(iq) / F0)
!!$     END IF
!!$  END DO
!!$!-----------------------------------------------------------------------
!!$  RETURN
!!$END SUBROUTINE Visi2D
!!$! ***********************************************************************
!!$
!**********************************************************************
double precision function Sexp(t)
!=======================================================================
! This is the function under the tau-integral; it is called from Romby.
! Here t = tau(iY); the flag 'transmit' is in 'SLBintens.inc'.
!======================================================================
  use common
  implicit none
  double precision t, arg, efact
  !-----------------------------------------------------------------
  if(transmit.eq.1) then
     ! for transmitted intensity
     arg = (taut-t)/muobs
  else
     ! for reflected intensity
     arg = t/muobs
  end if
  ! limits (in case of exp over/under flow)
  if(arg.gt.50.0d0) then
     efact = 0.0d0
  else
     efact = dexp(-arg)
  end if
  if(arg.lt.1.0d-6) efact = 1.0d0
  ! and finally the function under the integral:
  Sexp = Sfn * efact
  return
  !-----------------------------------------------------------------
end function Sexp
!***********************************************************************
!!$
!!$! ***********************************************************************
!!$DOUBLE PRECISION FUNCTION IMAGFN(Yy)
!!$! =======================================================================
!!$! This function evaluates auxiliary functions needed to produce
!!$! visibility curves and convolved images. It is called from the image
!!$! integration subroutine ROMBY.                        [Z.I., Jan. 1997]
!!$! =======================================================================
!!$  use common
!!$  IMPLICIT none
!!$  DOUBLE PRECISION x, Yy, PSFN, Bessel
!!$! -----------------------------------------------------------------------
!!$  IF (ftype.EQ.1) THEN
!!$     ! this part is for convolution
!!$     x = dsqrt(dabs(Cxout*Cxout+Yy*Yy-2.0D+00*Cxout*Yy*dcos(Cphi)))
!!$     imagfn = PSFN(x) * Yy**Ckn
!!$  ELSE
!!$     ! this part is for visibility
!!$     ! argument is Pi*q*Yy (not 2*Pi*q*Yy) to account for the fact that
!!$     ! theta1 is diameter rather than radius (so V is function of
!!$     ! q*theta1, like in IE, '96, MNRAS 279, 1019)
!!$     imagfn = Bessel(2.0D+00*dSIN(1.0D+00)*Cqtheta*Yy) * Yy**Ckn
!!$  END IF
!!$! -----------------------------------------------------------------------
!!$  RETURN
!!$END FUNCTION IMAGFN
!!$!***********************************************************************
!!$
!!$! ***********************************************************************
!!$SUBROUTINE Convolve
!!$! =======================================================================
!!$! This subroutine convolves intensity IntOut with the point spread
!!$! function to produce convolved images ConvInt. The work horse is
!!$! subroutine Conv2D, and this subroutine is used to prepare everything.
!!$!                                                      [Z.I., Jan. 1997]
!!$! Changed normalization of convolved images: instead of profiles normalized
!!$! at the center, now the convolved intensities are normalized by the area
!!$! A = 2*pi*Int{psf(x)*x}dx                             [MN, 2005]
!!$! =======================================================================
!!$  use common
!!$  IMPLICIT none
!!$  INTEGER i, j
!!$  DOUBLE PRECISION yang(npP+2), Youtang, deltaOff,    &
!!$       Int1D(npP+2), ConvS, Conv(npP+2), FWHM1max, FWHM2max, PSFN, &
!!$       PSFN1, psf_Y(1000), psf_X(1000), res
!!$!-----------------------------------------------------------------------
!!$  ! find the largest FWHMs
!!$  FWHM1max = FWHM1(1)
!!$  FWHM2max = FWHM2(1)
!!$  IF (psftype.LT.3) THEN
!!$     DO i = 1, NlambdaOut
!!$        IF (FWHM1(i).GT.FWHM1max) FWHM1max = FWHM1(i)
!!$        IF (FWHM2(i).GT.FWHM2max) FWHM2max = FWHM2(i)
!!$     END DO
!!$  END IF
!!$  ! scale angular coordinate to theta1
!!$  DO i = 1, nP+2
!!$     yang(i) = bOut(i) * Theta1 / 2.0D+00
!!$  END DO
!!$  ! generate off-set grid
!!$  Youtang = Y(nY) * Theta1
!!$  IF (Youtang.GT.FWHM1max.AND.Youtang.GT.FWHM2max) THEN
!!$     ! the envelope is well resolved, take impact parameter grid
!!$     Nconv = nP + 2
!!$     DO i = 1, Nconv
!!$        Offset(i) = yang(i)
!!$     END DO
!!$  ELSE
!!$     ! the envelope is not well resolved, take equidistant grid
!!$     ! to 2FWHM1max, i.e. image will be more or less the PSF itself
!!$     Nconv = 30
!!$     deltaOff = 2.0D+00 * FWHM1max / (Nconv-1)
!!$     IF (FWHM2max.GT.FWHM1max) THEN
!!$        deltaOff = 2.0D+00 *FWHM2max / (Nconv-1)
!!$     END IF
!!$     DO i = 1, Nconv
!!$        Offset(i) = deltaOff * 1.0D+00*(i-1)
!!$     END DO
!!$  END IF
!!$  ! !!!!!!!added normalization of psf, PSFN1(x) is the original fn
!!$  ! psfArea(iLambda) b/c the Gaussian option allows FWHM(iLambda)
!!$  write(12,'(a36)')'  lambdaOut(mic)   psfArea(arcsec^2)'
!!$  DO iLambda = 1, NlambdaOut
!!$     DO i = 1, Nconv
!!$        psf_X(i) = Offset(i)
!!$        psf_Y(i) = PSFN(Offset(i))
!!$        CALL ChkRange(dynrange,psf_Y(i))
!!$     END DO
!!$     ! CALL ScaleTo1(1000,Nconv,psf_Y)
!!$     CALL ScaletoArea(1000,Nconv,psf_X,psf_Y,res)
!!$     ! now for ea. lambda psf_Y is normalized to Area
!!$     psfArea(iLambda) = res
!!$     write(12,'(1p,2e15.3)') lambdaOut(iLambda), psfArea(iLambda)
!!$     j = iLambda
!!$     ! generate 1D intensity vector for subroutine Conv2D
!!$     ! take only diffuse emission, stellar contribution will
!!$     ! be added below (a shortcut to avoid inaccuracies or too many
!!$     ! points in Conv2D)
!!$     DO i = 1, nP+2
!!$        IF (i.LE.2) THEN
!!$           Int1D(i) = IntOut(j,3)
!!$        ELSE
!!$           Int1D(i) = IntOut(j,i)
!!$        END IF
!!$        CALL ChkRange(dynrange,Int1D(i))
!!$        IF (Int1D(i).LT.dynrange) Int1D(i)=0.0D+00
!!$     END DO
!!$     ! convolve
!!$     CALL Conv2D(npP+2,nP+2,yang,Int1D,1000,NConv,Offset,Conv)
!!$     DO i = 1, Nconv
!!$        CALL ChkRange(dynrange,Conv(i))
!!$        ConvInt(j,i) = Conv(i)
!!$     END DO
!!$     ! add stellar contribution
!!$     DO i = 1, nP+2
!!$        ConvS=2.0D+00*ASIN(1.0D+00)*(yang(2)**2.0D+00)*IntOut(j,1)*PSFN(Offset(i))
!!$        Conv(i) = Conv(i) + ConvS
!!$     END DO
!!$     ! scale to 1 at the center
!!$     ! CALL ScaleTo1(1000,Nconv,Conv)
!!$     ! copy 1D convolved intensity to ConvInt
!!$     ! ConvInt is normalized to area in Sub PrOut
!!$     DO i = 1, Nconv
!!$        CALL ChkRange(dynrange,Conv(i))
!!$        ConvInt(j,i) = Conv(i)
!!$     END DO
!!$     ! end do over iLambda - lambda for conv,images
!!$  END DO
!!$  write(12,*)' --------------------------------------------'
!!$  ! --------------------------------------------------------------------
!!$  RETURN
!!$END SUBROUTINE Convolve
!!$! ***********************************************************************
!!$
!!$! ***********************************************************************
!!$SUBROUTINE Conv2D(NinMax,Nin,Xin,Yin,Noutmax,Nout,Xout,Yout_loc)
!!$! =======================================================================
!!$! This subroutine convolves intensity Yin(Xin[i]), i=1,Nin with
!!$! the point spread function PSFN(x) (provided as a separate function).
!!$! It is assumed that both the intensity yin and PSFN(x) are circularly
!!$! symmetric functions of radial coordinate x, i.e., this subroutine
!!$! performs two-dimensional convolution. Convolved intensity, Yout, is
!!$! evaluated for very position Xout[i], i=1,Nout, as:
!!$!        Yout(Xout) = Int[Yin(Yloc)*PSF(xloc)*Yloc*dYloc*dphi]
!!$! where xloc = sqrt(Yloc **2+Xout**2-2*Yloc*Xout*cos(phi), with Yloc
!!$! and phi being dummy integration variables. Declared size of Xin is
!!$! NinMax, the one for Xout is NoutMax. The radial integration is done
!!$! using subroutine ROMBY and angular integration is done by using
!!$! Simpson rule.                                        [Z.I., Jan. 1997]
!!$! =======================================================================
!!$  use common
!!$  IMPLICIT none
!!$  INTEGER NinMax, Nin, NoutMax, Nout, iPhi, iXin, Nphi, iXOut
!!$  DOUBLE PRECISION Xin(NinMax), Yin(NinMax), Xout(NoutMax), A, B,  &
!!$       Yout_loc(NoutMax), dphi, phi_loc(1000), fphi(1000), int1, &
!!$       int2, imagfn
!!$  EXTERNAL imagfn
!!$!-----------------------------------------------------------------------
!!$  ! Parameters for integration:
!!$  ! number of angular points
!!$  Nphi = 9
!!$  ! step in angle phi
!!$  dphi = 2.0D+00*ASIN(1.0D+00) / (Nphi-1)
!!$  ! flag for imgfn
!!$  ftype = 1
!!$  ! Start integrations
!!$  ! loop over output positions
!!$  DO iXout = 1, Nout
!!$     Cxout = Xout(iXout)
!!$     ! loop over angular wedges (phi integration)
!!$     DO iPhi = 1, Nphi
!!$        phi_loc(iPhi) = dphi*1.0D+00*(iPhi-1)
!!$        Cphi = phi_loc(iPhi)
!!$        fphi(iPhi) = 0.0D+00
!!$        ! loop over input radial positions (radial integration)
!!$        DO iXin = 1, Nin-1
!!$           Ckn = 1.0D+00
!!$           CALL ROMBY(imagfn,Xin(iXin),Xin(iXin+1),int1)
!!$           Ckn = 2.0D+00
!!$           CALL ROMBY(imagfn,Xin(iXin),Xin(iXin+1),int2)
!!$           ! contribution from this annulus (lin. approx. for intensity)
!!$           A = Xin(iXin+1)*Yin(iXin) - Xin(iXin)*Yin(iXin+1)
!!$           A = A / (Xin(iXin+1)-Xin(iXin))
!!$           B = (Yin(iXin+1)-Yin(iXin)) / (Xin(iXin+1)-Xin(iXin))
!!$           fphi(iPhi) = fphi(iPhi) + A*int1 + B*int2
!!$        END DO
!!$     END DO
!!$     CALL Simpson(1000,1,Nphi,phi,fphi,Yout_loc(iXout))
!!$  END DO
!!$!-----------------------------------------------------------------------
!!$  RETURN
!!$END SUBROUTINE Conv2D
!!$! ***********************************************************************
!!$
! ***********************************************************************
DOUBLE PRECISION FUNCTION PSFN(x)
! =======================================================================
! This function evaluates the point spread function. For psftype.EQ.1
! the function is evaluated as a sum of two Gaussians, for psftype.EQ.3
! it is provided by user in a file. psftype and all other relevant
! parameters come from COMMON /psf/ and are initialized in subroutine
! INPUT.                                               [Z.I., Jan. 1997]
! =======================================================================
  use common
  IMPLICIT none
  DOUBLE PRECISION x
  INTEGER idummy
  ! -----------------------------------------------------------------------
  IF (psftype.LT.3) THEN
     psfn = dexp(-(1.665d0*x/FWHM1(iLambda))**2.0d0)
     IF (psftype.EQ.2)  &
          psfn = (psfn + kPSF(iLambda) *  &
          dexp(-(1.665d0*x/FWHM2(iLambda))**2.0d0))/ &
          (1.0d0+kPSF(iLambda))
  ELSE
     CALL LinInter(1000,Npsf,xpsf,ypsf,x,idummy,psfn)
  ENDIF
  ! -----------------------------------------------------------------------
  RETURN
END FUNCTION PSFN
! ***********************************************************************

!***********************************************************************
SUBROUTINE GetbOut(nP,pstar,k)
!=======================================================================
! This subroutine inserts two impact parameters corresponding to pstar,
! producing bOut(nP+2) from P(nP). The inserted elements are bOut(k) and
! bOut(k+1)                                            [Z.I., Aug. 1996]
! =======================================================================
  use common
  IMPLICIT none
  !---parameter
  integer nP
  double precision :: pstar
  !---local
  INTEGER k, kstop, i
! -----------------------------------------------------------------------
  k = 0
  kstop = 0
  DO WHILE (kstop.NE.1)
     k = k + 1
     bOut(k) = P(k)
     IF (1.001D+00*pstar.LE.P(k).OR.k.EQ.nP) kstop = 1
  END DO
  IF (0.999D+00*pstar.GT.P(k-1)) THEN
     bOut(k) = 0.999D+00*pstar
  ELSE
     bOut(k) = 0.5D+00*(P(k-1)+1.001D+00*pstar)
  END IF
  IF (1.001D+00*pstar.LT.P(k)) THEN
     bOut(k+1) = 1.001D+00*pstar
  ELSE
     bOut(k+1) = 0.5D+00*(P(k)+0.999D+00*pstar)
  END IF
  DO i = k, nP
     bOut(i+2) = P(i)
  END DO
  ! -----------------------------------------------------------------------
  RETURN
END SUBROUTINE GetbOut
!***********************************************************************

FUNCTION CUBIC_SPLINT(N, XI, FI, P2, X)
  IMPLICIT NONE
  INTEGER, INTENT (IN) :: N
  REAL*8, INTENT (IN), DIMENSION (N):: XI, FI, P2
  REAL*8, INTENT (IN) :: X
  REAL*8 :: CUBIC_SPLINT
  INTEGER klo,khi,k
  REAL*8 :: h,b,a,tmp
  INTEGER LOCAT

  klo=max(min(locat(XI,N,X),N-1),1)
  khi=klo+1
  h=XI(khi)-XI(klo)
  a=(XI(khi)-X)/h
  b=(X-XI(klo))/h
  tmp = a*FI(klo)+b*FI(khi)+((a*a*a-a)*P2(klo)+&
       (b*b*b-b)*P2(khi))*(h*h)/6.0D0
  !if ((tmp.ge.min(FI(klo),FI(khi))).and.(tmp.le.max(FI(klo),FI(khi)))) then 
     CUBIC_SPLINT= tmp
  !ELSE
  !   CUBIC_SPLINT=0.0
  !END IF
END FUNCTION CUBIC_SPLINT

function locat(v, n, x)
  !   Sei v ein geordneter Vektor der Lange n und x eine beliebige Zahl.
  !   Wenn  v(1) < v(2) < ... < v(n),  so liegt x  im halboffenen 
  !   Intervall  (v(j), v(j+1)].
  !   Falls   x <= v(1):   j = 0.   Falls   v(n) < x:   j = n
  !   Wenn  v(1) > v(2) > ... > v(n),  so liegt x  im halboffenen 
  !   Intervall  (v(j+1), v(j)].
  !   Falls   x > v(1):   j = 0.    Falls   v(n) >= x:  j = n
  implicit none
  real*8,intent(in)	:: v(n)
  integer,intent(in) 	:: n
  real*8,intent(in) 	:: x
  integer               :: locat
  integer		:: jlow,jup,jm
  jlow = 0
  jup  = n + 1
10 continue
  if(jup-jlow .gt. 1) then
     jm = (jup + jlow) / 2
     if( (v(n).gt.v(1) ) .eqv. ( x.gt.v(jm)) ) then
        jlow = jm
     else
        jup  = jm
     end if
     go to 10
  end if
  locat      = jlow
  return
end function locat

FUNCTION  Median(X, N)
  IMPLICIT  NONE
  REAL*8, DIMENSION(1:N), INTENT(IN) :: X
  INTEGER, INTENT(IN)                :: N
  REAL*8, DIMENSION(1:N)            :: Temp
  INTEGER                            :: i
  REAL*8 :: MEDIAN
  
  DO i = 1, N                       ! make a copy
     Temp(i) = X(i)
  END DO
  CALL  Sort(Temp, N)               ! sort the copy
  IF (MOD(N,2) == 0) THEN           ! compute the median
     Median = (Temp(N/2) + Temp(N/2+1)) / 2.0
  ELSE
     Median = Temp(N/2+1)
  END IF
END FUNCTION  Median
