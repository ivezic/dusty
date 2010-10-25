!!======================================================================
!! This is the current last version from August 2010, Dusty3.10.f90.
!! It is obtained from Dusty08.10.f90 by cleaning some unnecessary variables and comments,
!! fixing the RDW output and minor editing in Solve, Sub Flux_Cons and Sub SPH_ext_illum.
!!
!! All routines and functions are arranged in the following blocks:
!!    RADIATIVE TRANSFER CALCULATION
!!    GRIDS
!!    INPUT/OUTPUT
!!    ANALYSIS of the obtained solution
!!    RDW
!!    AUXILIARY MATH
!!  ----------
!!  More significant fixes of the version from 2009:
!!  1) p-integration of U and flux is done now by Nordlund integration, which greatly
!!     improves the accuracy (especially for the thermal bath test case).
!!
!!  2) Fixed the bug in Subroutine SPH_Diffuse.
!!     (the index of Utot() and Em() was iY instead of iYy, while Td had the dummy iYy
!!     in the same expression).
!!
!!  3) Fixed the intensities output. Before there were NaNs and 0.0s.
!!     Done by implementing relevant subroutines from old Dusty.
!!
!!  4) Set accConv for T-iterations to 1e-4.
!!
!!  5) Minor fixes:
!!    = Restored flag iPhys in 'dusty.inp'; when it is set to 1 the output fluxes and intensities
!!      are multiplied by the physical scale Jext(nY) in [W/m2] (see sub PrOut).
!!    = limval in sub PrOut was changed to 1e-30 but not used.
!!      Restored it to 1e-20, used as a lower limit (as in old Dusty).
!!    = Fixed the RDW output
!!
!!   General comments:
!!    - fluxes and en.densities are in lambda*f_lambda units, see sub Bolom for the integration.
!!
!!    - T_ext(npY) is Text^4 from the Blueprint;
!!
!!    - flag iInn is for internal use (search for iInn =), it gives additional output:
!!      '*.err' files are produced to trace convergence; in m-files Ubols, fbols profiles are added.
!!
!!    - Compared to the old Dusty there is a 4pi factor in fs(iL,iY) in the new code
!!      (see expressions in Table 5.1). This results in a 4pi factor difference in tauF and rg
!!      produced by new and old Dusty.
!!                                                         [Aug.2010, MN]


!!===========================================================================
  PROGRAM DUSTY
!============================================================================
! This program solves the continuum radiative transfer problem for a
! spherically symmetric envelope or for a plane-parallel slab.
! It is assumed that all input data are given in files named *.inp,
! and that the list of these files is given in a master input file 'dusty.inp'.
!
! This is Dusty version from summer 2010, with major changes in the numerical
! method, and greatly improved numerical stability at high optical depth.
! Input files for older Dusty versions are incompatible with this one.
! For details see the Manual.
!                              [Zeljko Ivezic, Maia Nenkova, Mridupawan Deka]
!============================================================================
  use common
  implicit none

  integer error, length, nG, model, Nmodel, GridType, io1, empty, lpath, &
       Nrec, lambdaOK,n_init_model,iterfbol,fbolOK
! Nrec is the max number of records for taugrid in a file
  parameter (Nrec = 1000)
  double precision value, tau1, tau2, tauIn(Nrec), pstar, us(npL,npY),delta, RDINP
  double precision, allocatable :: tau(:)
  character*12 version
  character*235 dustyinpfile,arg, path, apath, nameIn, nameOut, nameQ(npG), &
       nameNK(10), stdf(7), str
  logical UCASE,equal,initial, Lprint
!----------------------------------------------------------------------

!**************************
!*** ABOUT THIS VERSION ***
!**************************

  version= '3.10'

!********************  MAIN  *******************************

  equal = .true.
! First read lambda grid
  call ChkLambda(lambdaOK)
  if (lambdaOK.eq.0) then
     goto 999
  end if
! Open master input file; first determine whether user supplies custom
! DUSTY input file as the 1st argument on the command line. If yes,
! use it. Else revert to the default file ./dusty.inp
  call getarg(1,dustyinpfile)
  if (trim(dustyinpfile) == "") then
     write(*,*) "No input file name found on command line. Proceeding with default file ./dusty.inp"
     dustyinpfile = "dusty.inp"
  else
     write(*,*) "Found input file ", trim(dustyinpfile), " on on command line."
  endif
!  open(13,err=998,file='dusty.inp',status='old')
  open(13,err=998,file=trim(dustyinpfile),status='old')
  io1 = 0
! read the verbose mode
  iVerb = RDINP(Equal,13)
! read flag for spectra in W/m^2
  iPhys = RDINP(Equal,13)
! loop over input files
  do while (io1.ge.0)
! read a line from master input file using
100  read(13,'(a)',iostat=io1) apath
     if(io1.lt.0) then
        stop
     end if
     call clean(apath, path, lpath)
! if not eof and if line is not empty, or commented, proceed
     if (empty(path).ne.1) then
! get input/output file names
        call attach(path,lpath,'.inp',nameIn)
        call attach(path,lpath,'.out',nameOut)
! read input data
        call Input(nameIn,nG,nameOut,nameQ,nameNK,tau1,tau2,&
             tauIn,Nrec,GridType,Nmodel,error,version,stdf)
        if (iVerb.gt.0) write(*,'(a24,a235)') ' working on input file: ',nameIn
        if (iVerb.eq.2) write(*,*) 'Done with reading input'
! if an error reading files go to the next input file
! error=3 means some files are missing
        if (error.eq.3) goto 100
! get optical properties
        call getOptPr(nG,nameQ,nameNK,error,stdf)
! if an error reading files go to the next input file
        if (error.eq.3) goto 100
        if (iVerb.eq.2) write(*,*) 'Done with getOptPr'
        if(allocated(tau)) deallocate(tau)
        allocate(tau(Nmodel))
        if (error.eq.0) then
           call GetTau(nG,tau1,tau2,tauIn,Nrec,GridType,Nmodel,tau)
           if (iVerb.eq.2) write(*,*) 'Done with GetTau'
           call Kernel(nG,path,lpath,tauIn,tau,Nrec,Nmodel,GridType,error,Lprint)
        endif
        if (error.ne.0) go to 100
     end if
! end of the loop over input files
  end do
  if (iVerb.gt.0.and.Lprint) write(*,*) ' End of input files '
  close(13)
! end this run
  goto 999
! to execute if the master input file is missing
998 write(*,*)' *********** Fatal Error in DUSTY ********************'
  write(*,*)' * Problem finding input file ',dustyinpfile,'!? '
  write(*,*)' *****************************************************'
!----------------------------------------------------------------------

999 stop
end program DUSTY
!**********************************************************************


!!=======================================================================
! Routines and funcitons related to the RADIATIVE TRANSFER CALCULATION,
! arranged in alphabetical order.                        [MN, Aug.2010]
!!=======================================================================


!***********************************************************************
subroutine Emission(nG,T_ext,emiss)
!=======================================================================
! This subroutine calculates emission term from the temperature and abund
! arrays for flag=0, and adds U to it for flag=1.
!                                                      [Z.I., Mar. 1996]
!=======================================================================
  use common, only: nY,npY,nL,npL,lambda,Td, abund,dynrange
  implicit none

  integer iG, iY, iL, nG
  double precision  emiss(npL,npY),emig, tt, xP,Planck,T_ext(npY)
! -----------------------------------------------------------------------

! first initialize Emiss
  emiss = 0.0d0
! calculate emission term for each component and add it to emiss
! loop over wavelengths
  do iL = 1, nL
! loop over radial coordinate
     do iY = 1, nY
!   loop over grains
        do iG = 1, nG
           xP = 14400.0d0/(lambda(iL)*Td(iG,iY))
           tt = (Td(iG,iY)**4.0d0) / T_ext(iY)
           emig = abund(iG,iY)*tt*Planck(xP)
           if (emig.lt.dynrange*dynrange) emig = 0.0d0
!     add contribution for current grains
           emiss(iL,iY) = emiss(iL,iY) + emig
        end do
        if (emiss(iL,iY).lt.dynrange*dynrange) emiss(iL,iY) = 0.0d0
     end do
  end do
! -----------------------------------------------------------------------
  return
end subroutine Emission
!***********************************************************************

!***********************************************************************
subroutine Find_Diffuse(initial,iter,iterfbol,T_ext,us,em,omega,error)
!=======================================================================
! This subroutine finds the diffuse en.density for slab and sphere,
! and the diffuse flux for sphere.                       [Deka,'08]
! Renamed some variables for easier reading:
! faux becomes Sfn_em, faux1 becomes Sfn_sc, removed b(npY),using Utot directly instead.
!                                                        [MN'10]
!=======================================================================
  use common
  implicit none

  integer iL,iY,iP,j, iYaux, kronecker, error,iter,iterfbol,flag, nZ
  double precision us(npL,npY), em(npL,npY), omega(npG,npL),eint2, &
       sph_em(npL,npY),sph_sc(npL,npY),tau(npY),T_ext(npY),Sfn_em(npY),Sfn_sc(npY),sum1,sum2, dyn2
  logical initial
  external eint2
!-----------------------------------------------------------------------
  dyn2 = 1.0d-30
  error = 0
  Sfn_em = 0.0d0
  Sfn_sc = 0.0d0

! for slab calculate new energy density
! loop over wavelengths
  if(slb) then !------------for slab ----------------
     do iL = 1, nL
        do iY = 1,nY
           tau(iY) = TAUslb(iL,iY)
           Sfn_em(iY) = (1.0d0-omega(1,iL))*em(iL,iY)/2.0d0
           if (initial.and.iter.eq.1.and.iterfbol.eq.1) then
              Sfn_sc(iY) = omega(1,iL)*us(iL,iY)/2.0d0
           else
              Sfn_sc(iY) = omega(1,iL)*utot(iL,iY)/2.0d0
           end if
        end do
!    integrate to get diffuse emission and scattering
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
  elseif(sph) then !-----------for sphere ----------------
     do moment_loc = 1, moment
        if (moment_loc.eq.1) then
!      Find diffuse en.density
           call SPH_diff(1,initial,iter,iterfbol,T_ext,em,omega,sph_em)
	   call SPH_diff(2,initial,iter,iterfbol,T_ext,us,omega,sph_sc)
!      find total energy density
           do iY = 1, nY
              do  iL = 1, nL
                 utot(iL,iY) = us(iL,iY) + sph_em(iL,iY) + sph_sc(iL,iY)
                 Ude(iL,iY) = sph_em(iL,iY)
                 Uds(iL,iY) = sph_sc(iL,iY)
              end do
           end do
        elseif (moment_loc.eq.2) then
!      Find diffuse fluxes
           call SPH_diff(1,initial,iter,iterfbol,T_ext,em,omega,fde)
           call SPH_diff(2,initial,iter,iterfbol,T_ext,us,omega,fds)
        end if
     end do
  end if
!-----------------------------------------------------------------------
  return
end subroutine Find_Diffuse
!***********************************************************************

!***********************************************************************
subroutine Find_Temp(nG,T_ext)
!=======================================================================
! This subroutine finds new temperature from Utot.
! Temperature is obtained by solving:
!   f2(y)*f1(y) - g(T) = 0
! where
!   f1(y) = Int(Qabs*Utot*dlambda)
!   f2(y) = T_ext
!   g(T) = qP(Td)*Td**4
!   T_ext is the effective T of the illuminating radiation
!                                                [ZI'96, MN,'00, Deka, 2008]
!=======================================================================
  use common
  implicit none

  integer iG, iY, iL,nG
  double precision  fnum(npL),fdenum(npL),xP,Planck, qpt1,qu1,us(npL,npY),&
       T_ext(npY),gg, ff(npL), fnum1
!-----------------------------------------------------------------------

! loop over grains
  do iG = 1, nG
! if T1 given in input:
     if(typentry(1).eq.5) call find_Text(nG,T_ext)
! loop over radial positions (solving f1-f2*g(t)=0)
     do iY = 1, nY
! calculate f1 and f2
        do iL = 1, nL
           fnum(iL) = sigmaA(iG,iL)*utot(iL,iY)/lambda(iL)
           xP = 14400.0d0/lambda(iL)/Td(iG,iY)
           ff(iL) = sigmaA(iG,iL)*Planck(xP)/ lambda(iL)
        end do
        call Simpson(npL,1,nL,lambda,fnum,fnum1)
        call Simpson(npL,1,nL,lambda,ff,gg)
        Td(iG,iY) = (fnum1*T_ext(iY)/gg)**(1.0d0/4.0d0)
     end do
  end do
!-----------------------------------------------------------------------

  return
end subroutine Find_Temp
!***********************************************************************


!***********************************************************************
subroutine Find_Text(nG,T_ext)
!=======================================================================
  use common
  implicit none

  integer iG,iY,iL,nG
  double precision  fnum(npL),fdenum(npL),xP,Planck, qPT1,qU1,T_ext(npY)
!-----------------------------------------------------------------------

! loop over grains
  do iG = 1, nG
     do iL = 1, nL
        fnum(iL) = sigmaA(iG,iL)*utot(iL,1)/lambda(iL)
        xP = 14400.0d0/lambda(iL)/ Tsub(1)
        fdenum(iL) = sigmaA(iG,iL)*Planck(xP)/lambda(iL)
     end do
     call Simpson(npL,1,nL,lambda,fnum,qU1)
     call Simpson(npL,1,nL,lambda,fdenum,qPT1)
     if (sph) then
        do iY = 1, nY
           T_ext(iY) = (Tsub(1)**4.0d0)*(qPT1/qU1)/(Y(iY)**2.0d0)
        end do
     elseif(slb) then
        T_ext = Tsub(1)**4.0d0*(qPT1/qU1)
     end if
  end do
!-----------------------------------------------------------------------
  return
end subroutine Find_Text
!***********************************************************************


!***********************************************************************
subroutine Find_Tran(pstar,T_ext,us,fs)
!=======================================================================
! This subroutine generates the transmitted energy density and  flux for
! the slab and spere. is=1 is for the left src, is=2 is the right src
!                                               [MN, Feb.'99, Deka, 2008]
!=======================================================================
  use common
  implicit none

  integer i, iY, iL, iz, is, nZ, nn, error
  double precision  shpL(npL), shpR(npL), usL(npL,npY),arg, dyn2,usR(npL,npY), &
       eint2, eint3, x, res, &
        denom, expow, m0(npL,npY), m1(npL,npY),m1p(npL,npY), &
       m1m(npL,npY), tauaux(npL,npY), zeta, result1, eta
  double precision, intent(in) :: pstar
  double precision, intent(out) :: fs(npL,npY),us(npL,npY),T_ext(npY)
  external eta
!-----------------------------------------------------------------------

  dyn2 = 1.0d-30

! intialise
  usL = 0.0d0
  usR = 0.0d0
  us = 0.0d0
  fsL = 0.0d0
  fsR = 0.0d0
  fs = 0.0d0
  fsLbol = 0.0d0
  fsRbol = 0.0d0
  fsbol = 0.0d0

! left-side source is always present
  if(left.gt.0) then
     call getSpShape(shpL,1)
  else
     shpl = 0.0d0
  end if
  if(right.gt.0) then
     call getSpShape(shpR,2)
  else
     shpr = 0.0d0
  end if

! define T_external for typEntry(1).ne.5
  if (typEntry(1).ne.5) then
! for slab
     if(slb) then
        if (left.eq.1.and.right.eq.0) then
           if (mu1.eq.-1.0d0) then
              T_ext = pi*Ji/(2.0d0*sigma)
           else
              T_ext = pi*Ji/sigma
           end if
        elseif (left.eq.0.and.right.eq.1) then
           if (mu2.eq.-1.0d0) then
              T_ext = pi*ksi*Ji/(2.0d0*sigma)
           else
              T_ext = pi*ksi*Ji/(sigma)
           end if
        elseif (left.eq.1.and.right.eq.1) then
           if (mu1.eq.-1.0d0.and.mu2.eq.-1.0d0) then
              T_ext = pi*(Ji + ksi*Ji)/(2.0d0*sigma)
           else
              T_ext = pi*(Ji + ksi*Ji)/sigma
           end if
        end if
! for sphere
     elseif(sph) then
        do iY = 1, nY
           if (left.eq.1.and.right.eq.0) then
              T_ext(iY) =  (pi/sigma)*Ji/(Y(iY)**2.0d0)
           elseif (left.eq.0.and.right.eq.1) then
              T_ext(iY) =  pi*Jo/sigma
           elseif (left.eq.1.and.right.eq.1) then
              T_ext(iY) =  (pi/sigma)*(Ji/(Y(iY)**2.0d0) + Jo)
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
     call SPH_ext_illum(m0,m1,m1p,m1m)
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
!    fs(iL,iY) =  (fsL(iL,iY) - fsR(iL,iY))
        elseif(sph) then
           if (left.eq.1.and.right.eq.0) then
              us(iL,iY) = usL(iL,iY)
              fs(iL,iY) = fsL(iL,iY)
           elseif(left.eq.0.and.right.eq.1) then
!!**     fs(iL,iY) = fsR(iL,iY)  a '-' sign is needed (see below, MN)
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
  call bolom(fsL,fsLbol)
  call bolom(fsR,fsRbol)
  call bolom(fs,fsbol)
  error = 0
  goto 999

998 write(12,*)' *** fatal error in dusty! *************************'
  write(12,*)' file with the spectral shape of external radiation:'
  write(12,'(a2,a256)')'  ',namestar(1), namestar(2)
  write(12,*)' is missing or not properly formatted?!'
  write(12,*)' ***************************************************'
  error = 3
!-----------------------------------------------------------------------
999 return
end subroutine Find_Tran
!***********************************************************************
!***********************************************************************
subroutine Flux_Consv(flux1,flux2,fbolOK,error,Lprint)
!=======================================================================
! Replaces the former SUBROUTINE ChkFlux(flux,tolern,consfl,error,ETAzp)
!
! This subroutine checks the bolometric diffuse flux conservation at any
! point of the grid. In case of nonconservation inserts a number of points
!  at certain places.                            [Deka'08, MN'99; ZI'96]
!=======================================================================

  use common
  implicit none

  integer iYins(npY), k, kins, i, iY, idm, nn, flag, error, istop, fbolOK
  double precision  tauaux(npY),Yins(npY), flux1(npY), flux2(npY), &
       deltaumax, ratio(npY), etatemp(npY),ee,result1,Yloc, eta,x1, x2, tmp
  double precision, dimension(:), allocatable:: xg,  wg
  external eta
  LOGICAL Lprint
!--------------------------------------------------------------------------

  flag= 0
  error = 0.0d0
  kins = 0
  istop = 0
  ratio = 0
  maxrat = 0.0d0

  if(sph) then
! save old grid and values of Eta (important for denstyp = 5 or 6)
! for spherical case
   if (rdw) then
    Yprev = y
    etatemp = etadiscr
    nYprev = nY
   end if

   do iY = 1, nY
      tauaux(iY) = TAUtot(1)*ETAzp(1,iY)
   end do
!!** I am not sure if deltaumax has to be found as below: [MN]
!  maximal deltau is no more than 2 times the average value
   deltaumax = 2.0d0*tauaux(nY)/nY
  elseif(slb) then
     deltaumax = 2.0d0*tautot(1)/nY
  end if

! search for places to improve the grid
  do  iY = 2, nY
   if(sph) then
     if(Y(iY).gt.1.5d0) then
       if(left.eq.1) then
         ratio(iY) = abs(abs(flux1(iY))-abs(flux2(iY)))/(abs(flux1(iY))+abs(flux2(iY)))
       elseif(left.eq.0.and.right.eq.1.and.(flux1(iY).gt.0.05d0*pi)) then
         ratio(iY) = abs(abs(flux1(iY))-abs(flux2(iY)))/(abs(flux1(iY))+abs(flux2(iY)))
       end if
     end if
   elseif(slb) then
     ratio(iY) = abs(abs(flux1(iY))-abs(flux2(iY)))/(abs(flux1(iY))+abs(flux2(iY)))
   end if
    if (ratio(iY).gt.maxrat) maxrat = abs(ratio(iY))
  end do

! For very high optical depth (taumax > 200), change the accuracy to 20% [Deka'09]
!  if(taumax.ge.200.0d0) accuracy = 0.2d0
  if (taumax.ge.200.0d0.and.(2*nY-1).gt.npY) accuracy = 0.2d0

! if any of these criteria is satisfied insert a point:
  DO WHILE (istop.ne.1)
    if (maxrat.gt.accuracy) flag=1
    do iY = 2, nY
     if(sph) then
! if deltau is too large
      if ((tauaux(iY)-tauaux(iY-1)).gt.deltaumax) flag = 1
     elseif(slb) then
! if deltau is too large at the left edge:
      if(TAUslb(1,iY).lt.5.0d0) then
       if ((TAUslb(1,iY)-TAUslb(1,iY-1)).gt.deltaumax) flag = 1
      end if
     end if
     if(flag.eq.1.and.maxrat.gt.accuracy) then
      kins = kins + 1
      Yins(kins) = Y(iY-1)+0.5d0*(Y(iY)-Y(iY-1))
      iYins(kins) = iY-1
     end if
    enddo

   if (maxrat.le.accuracy) then
     istop = 1
   else
     if (kins.gt.0) istop = 1
   end if
! end do over places to improve the grid
  END DO

  IF (kins.eq.0) THEN
!!! why do we need the if statements for slb, sph, why not fbolOK=1
   if(sph) then
     fbolOK = 1
   elseif(slb) then
     fbolOK = 1
   end if
  ELSE
! add all new points to Y(nY). this gives the new Y(nY+kins).
! however, check if npY is large enough to insert all points:
   if ((nY+kins).gt.npY) then
    fbolOK = 1
    if (iX.ge.1.AND.Lprint) then
     write(18,*)' ****************     WARNING   ******************'
     write(18,*)'  The new Y-grid can not accomodate more points!'
     write(18,'(a,i5)')'   The specified accuracy would require',nY+kins
     write(18,'(a,i5,a)')'   points, while npY =',npY,'.'
     write(18,*)'  For the required accuracy npY must be increased,'
     write(18,*)'  (see the manual s3.5 numerical accuracy).'
     write(18,'(a37,F5.1,a2)')'   The currently achieved accuracy is ', maxrat*100.0, ' %'
     write(18,*)' *************************************************'
    end if
!!     kins = npY - nY     !!this is in the old code, but doesn't work here.
     iWARNING = iWARNING + 1
!!     error = 2           !!this is in the old code, but doesn't work here.
    go to 777
   else
     do k = 1, kins
       call shift(Y,npY,nY+k-1,Yins(k),iYins(k)+k-1)
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
     call lininter(npY,nYprev,Yprev,etatemp,Yloc,idm,ee)
     etadiscr(iY) = ee
    else
     etadiscr(iY) = eta(Yloc)
    end if
   end do
  end if
!--------------------------------------------------------------------------

777 return
end subroutine Flux_Consv
!*********************************************************************

!***********************************************************************
subroutine Init_Temp(nG,T_ext,us)
!=======================================================================
! This subroutine calculates the initial approximation for the temperature.
! Temperature is obtained by solving:
!   g(T) = f2(y)*f1(y)
! where
!   g(T) = qP(Td)*Td**4
!   f1(y) = Int(Qabs*Utot*dlambda)
!   f2(y) = T_ext
!   T_ext is the effective T of the illuminating radiation [Deka, July'08]
!
!!**   Minor editing to avoid having IF's inside do-loops [MN, Aug'10]
!=======================================================================

  use common
  implicit none

  integer iL, iY, iG, nG, iw
  double precision xP, Planck, us(npL,npY),T_ext(npY),fnum(npL),qP,ff(npL), fnum1
!--------------------------------------------------------------------------

  if(typentry(1).eq.5) then
!   this is if Tsub(1) given in input
    if(sph) then
      do iY = 1, nY
        T_ext(iY) = Tsub(1)**4.0d0/Y(iY)**2.0d0   !eq.(4.1.22)
      end do
    else
!     for slab
      do iY = 1, nY
        T_ext(iY) = Tsub(1)**4.0d0
      end do
    end if
!   first approximation for temperature
    do iG = 1, nG
      do iY = 1, nY
	     Td(iG,iY) = Tsub(1)
      end do
    end do

  else !if flux Fe1 is given in input in some form - at this point
!   T_ext(iY) is found in sub Find_Tran
    do iG = 1, nG
      do iY = 1, nY
         Td(iG,iY) = T_ext(iY)**(1.0d0/4.0d0)
      end do
    end do
  end if

 if(typentry(1).ne.5) then
   do iG = 1, nG
!   loop over radial positions
    do iY = 1, nY
!     calculate fnum1 and qP integrals
      do iL = 1, nL
        fnum(iL) = sigmaA(iG,iL)*us(iL,iY)/lambda(iL)
        xP = 14400.0d0/(lambda(iL)*Td(iG,iY))
        ff(iL) = sigmaA(iG,iL)*Planck(xP)/lambda(iL)
      end do
      call Simpson(npL,1,nL,lambda,fnum,fnum1)
      call Simpson(npL,1,nL,lambda,ff,qP)
!     get initial temperature
      Td(iG,iY) = (T_ext(iY)*fnum1/qP)**(1.0d0/4.0d0)
    end do
   end do
 end if

!-----------------------------------------------------------------------
  return
end subroutine Init_Temp
!***********************************************************************

!***********************************************************************
subroutine Kernel(nG,path,lpath,tauIn,tau,Nrec,Nmodel,GridType,error,Lprint)
!=======================================================================
! This subroutine generates solution for all optical depth steps starting
! with intial optical depth given in the input file. It prints output
! only for input models
!                                               [MN, May'10]
!=======================================================================
  use common
  implicit none

  integer iG, nG, iY, iL, i,j, GridType, model, Nmodel, Nrec, &
       iterfbol, fbolOK, istop, error, lpath
  double precision tauIn(Nrec),tau(Nmodel), ratio, delta, tau0
  character*235 path
  logical initial,Lprint
!-----------------------------------------------------------------------

! intialize print and initial flag
  Lprint = .false.
  initial = .false.

  istop = 0
!  init_tau = 1.0
  model = 1
  i = 1
  j = 0

! loop over optical depths
  do while (istop.eq.0)
! intialize print and initial flag
    Lprint = .false.
    initial = .false.
!   j is to select the small initial optical depth (set by user in .inp file)
    j = j + 1
    if (j.eq.1) then
       initial = .true.
       call GetTau0(tau0,1)
    else
       initial = .false.
    end if

    tau0 = tau(model)
    Lprint = .true.
    taufid = tau0
    
    print*, 'tau0=', tau0
    call GetTaumax(tau0,nG)

    if (iVerb.gt.0.and.Lprint)  write(*,'(a9,i4)') ' model = ',model
    if(tau0.eq.tau(model)) call OPPEN(model,path,lpath)
    if (iVerb.eq.2.and.Lprint) write(*,*) ' going to Solve '

!   solve radiative transfer for this particular optical depth
    call Solve(model,Lprint,initial,nG,error,delta,iterfbol,fbolOK)

!  if flux is conserved, the solution is obtained. So write out the values
!  to output files for specified models
    if (fbolOK.eq.1) then
      if(Lprint) then
        if (error.eq.0) then
          call Spectral(model)
          if (iVerb.eq.2) write(*,*) 'Done with Spectral'
          call PrOut(model,nG,delta)
        else
          go to 10
        end if
        call CLLOSE(error,model,Nmodel)
        if (iVerb.eq.2) write(*,*) ' ----------- '
        if (iVerb.eq.2) write(*,*) '             '
        model = model + 1
      end if
    end if
!    15 is the file with spectra
!     close(15)
!    16 is the file with radial properties
!     close(16)
!    17 is the file for images
!     close(17)
!    25 is the file for slab z-spectra
!     close(25)
!    18 is the message file
!     close(18)
!  if the last model is solved, then stop
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
double precision function Planck(x)
!=======================================================================
! This function evaluates the Planck function multiplied by wavelength
! and normalized by sigma*T^4/Pi.                      [Z.I., Mar. 1996]
! =======================================================================
  implicit none
  double precision x
! -----------------------------------------------------------------------

  if (x.gt.100.0d0) then
   Planck = 0.0d0
  else
   if (x.lt.0.00001d0) then
    Planck = 0.155d0*x**3.0d0
   else
    Planck = 0.155d0*x**4.0d0/(dexp(x) - 1.0d0)
   end if
  end if
! -----------------------------------------------------------------------
  return
end function Planck
!***********************************************************************

!********************************************************************
subroutine Rad_Transf(nG,Lprint,initial,pstar,y_incr,us,fs,em,omega, &
     iterfbol,T_ext)
!======================================================================

  use common
  implicit none
  integer i,iY,iY1,iL,iG,nn,itlim,imu,error,conv,iter,iPstar, iP, iZ, nZ
  integer, intent(in)::nG, y_incr,iterfbol
  double precision  em(npL,npY), tauaux(npY), pstar, result1, omega(npG,npL),&
       fs(npL,npY), us(npL,npY), T_ext(npY), T_old(nG,npY), u_old(npL,npY), &
       maxerrT,maxerrU, aux1,aux2, x1,x2, eta, &
       fDebol(npY), fDsbol(npY), Usbol(npY), Udebol(npY), Udsbol(npY), &
	   U_prev(npL,npY), Ubol_old(npY)
  double precision, dimension(:), allocatable:: xg, wg
  logical, intent(in) ::  initial,Lprint
  external eta
!--------------------------------------------------------------------------

  error = 0
  if(sph) then
!     generate spline coefficients for ETA as in old Dusty [MN'Aug,10]
      CALL setupETA
!     evaluate ETAzp (carried in common)
      CALL getETAzp
  end if

! the tau-profile at the fiducious lambda (needed in prout)
  if (slb) then
     do iY = 1, nY
        tr(iY) = TAUslb(iLfid,iY)/TAUslb(iLfid,nY)
     end do
  elseif(sph) then
     do iY = 1, nY
       tr(iY) = ETAzp(1,iY)/ETAzp(1,nY)
     end do
  end if

! generate albedo through the envelope
  call getOmega(nG,omega)
! generate stellar spectrum
  call Find_Tran(pstar,T_ext,us,fs)
  if(iVerb.eq.2.and.Lprint) write(*,*)' Done with transmitted radiation.'
! issue a message in fname.out about the condition for neglecting
! occultation only if T1 is given in input:
  if(typentry(1).eq.5.and.sph) then
   if(iterfbol.eq.1.and.itereta.eq.1.and.right.eq.0) call OccltMSG(us)
  end if
! finish when file with the stellar spectrum is not available
  if (error.eq.3) goto 999

! in the case of first (lowest) optical depth,
! us is the intial approximation for utot(iL,iY) for the first iteration over Td
! Find initial approximation of Td for the case of first (lowest) optical depth.
  if (initial.and.iterfbol.eq.1) then
   call Init_Temp(nG,T_ext,us)
   if(iVerb.eq.2.and.Lprint) write(*,*)' Done with initial dust temperature.'
  end if

!!**  added to check convergence, recorded in 'fname.err' (unit=38)[MN]
    IF(iInn.eq.1.AND.Lprint) write(38,'(2(a8,f5.0))') 'taufid0=',taufid0
    IF(iInn.eq.1.AND.Lprint) write(38,*)'iter    maxerrT    maxerrU'
!!**

  itlim = 50000
  conv = 0
  iter = 0
!=== iterations over dust temperature =========
  do while (conv.eq.0.and.iter.le.itlim)
   iter = iter + 1
!  if y-grid points are increased, linearly interpolate previous
!  dust temperature for the new y-grid
   if (y_incr.eq.1.and.iter.eq.1) then
    T_old = Td
    u_old = utot
    do iY = 1, nY, 2
     iY1 = int((iY + 1)/2)
     do iG = 1, nG
      Td(iG,iY) = T_old(iG,iY1)
     end do
     do iL = 1, nL
      utot(iL,iY) = u_old(iL,iY1)
     end do
    end do
    do  iY = 2, nY-1, 2
     do iG = 1, nG
      Td(iG,iY) = abs(Td(iG,iY+1) + Td(iG,iY-1))/2.0d0
     end do
     do iL = 1, nL
      utot(iL,iY) =  abs(utot(iL,iY+1) + utot(iL,iY-1))/2.0d0
     end do
    end do
   end if
!  find T_external for the new y-grid if T(1) given in input
   if (initial.and.iterfbol.ne.1.and.typentry(1).eq.5) then
    call find_Text(nG,T_ext)
   elseif (.not.initial.and.typentry(1).eq.5) then
    call find_Text(nG,T_ext)
   end if
!!** this is to check convergence on U [MN]
   U_prev = Utot

!  find emission term
   call Emission(nG,T_ext,em)
!  moment = 1 is for finding total energy density only
   moment = 1
   call Find_Diffuse(initial,iter,iterfbol,T_ext,us,em,omega,error)

!  assign previus Td to T_old
   T_old = Td
!  find Td
   call Find_Temp(nG,T_ext)
!  check convergence for dust temperature
   maxerrT = 0.0d0
   aux1 = 0.0d0
   do iG = 1,nG
    do iY = 1,nY
     aux1 = dabs(T_old(iG,iY) - Td(iG,iY))/Td(iG,iY)
     if (aux1.gt.maxerrT) maxerrT = aux1
    enddo
   end do
   if (maxerrT.le.accConv) conv = 1

!!**
   IF(iInn.eq.1 .AND. Lprint) THEN
      CALL Bolom(utot,Ubol)
      CALL Bolom(U_prev,Ubol_old)
      maxerrU = 0.0d0  !just for info, iterations are over Td [MN]
      do iY = 1,nY
        aux1 = dabs(Ubol_old(iY) - Ubol(iY)) / Ubol(iY)
        if (aux1.gt.maxerrU) maxerrU = aux1
      enddo
      write(38,'(i4,1p,2e12.3)') iter, maxerrT, maxerrU
   END IF
!!** -------------------------

  enddo
!=== the end of iterations over Td ===

  if(iVerb.eq.2.and.Lprint) write(*,*) ' Done with finding dust temperature. '
  if (iX.ge.1) then
   if (Lprint.and.iter.le.itlim) then
    write(18,*)' Convergence achieved, number of'
    write(18,'(a36,i4)') '  iterations over dust temperature: ',iter
   end if
  end if
! find T_external for the converged dust temperature
  if (typentry(1).eq.5) call find_Text(nG,T_ext)  !if T1 is given in input
! calculate the emission term using the converged Td
  call Emission(nG,T_ext,em)
! calculate total energy density and diffuse flux using the converged Td
  moment = 2
  call Find_Diffuse(initial,iter,iterfbol,T_ext,us,em,omega,error)
  if(iVerb.eq.2.and.Lprint) write(*,*) ' Done with finding energy density and diffuse flux.'
!-----------------------------------------------------------------

!!** additional printout [MN] vvvvvv
  IF(iInn.eq.1 .AND. Lprint) THEN
    write(18,'(a9,f5.0))') ' taufid0=',taufid0
    write(18,'(2(a9,i4))') 'iterfbol=',iterfbol,'     nY =', nY

    call Bolom(Utot,Ubol)
    call Bolom(Us,Usbol)
    call Bolom(Ude,Udebol)
    call Bolom(Uds,Udsbol)
!   fbols are printed from Solve (because slab fluxes are found there)
    IF (SLB) THEN
        write(18,*)'      tr       Ubol      Usbol      Udsbol     Udebol'
        DO iY = 1, nY
          write(18,'(1p,5E11.3)') tr(iY),Ubol(iY),Usbol(iY),Udsbol(iY),Udebol(iY)
        END DO
    ELSE
        write(18,*)'      Y       Ubol      Usbol      Udsbol     Udebol'
        DO iY = 1, nY
          write(18,'(1p,5E11.3)') Y(iY),Ubol(iY),Usbol(iY),Udsbol(iY),Udebol(iY)
        END DO
    END IF
  END IF
!!** ^^^^^^^^^^^^^^^^^^^^^

999 return
end subroutine Rad_Transf
!***********************************************************************

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
!-----------------------------------------------------------------

  return
end function Sexp
!***********************************************************************

!*********************************************************************
subroutine SLBdiff(flag,om,grid,T_ext,mat1,nL,nY,mat2,fp,fm)
!=========================================================================
! Integration of U(lam,t)*E2|Tau-t| to get the diffuse flux.
! flag=1 is for scattered and flag=0 for emitted flux. [MN, Apr'98]
!=========================================================================

  implicit none
  integer npY, npP, npX, npL, npG, npR
  include 'userpar.inc'
  parameter (npG=1)
  integer iL, iY, j, nL, nY, flag
  double precision mat1(npL,npY), mat2(npL,npY), grid(npL,npY), &
       om(npG,npL), fp(npL,npY), fm(npL,npY), tau(npY), &
       faux(npY), efact, fave, sum, eint3,pi,T_ext(npY)
!--------------------------------------------------------------------

pi=4.0d0*atan(1.0d0)

  do iL = 1, nL
   do iY = 1, nY
    tau(iY) = grid(iL,iY)
    if (flag.eq.1) then
     faux(iY) = om(1,iL)*mat1(iL,iY)
    else
     faux(iY) = (1.0d0-om(1,iL))*mat1(iL,iY)
    end if
   end do

! find f(+) (in arg tau>t)
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

   do iY = 1, nY
    mat2(iL,iY) = (fp(iL,iY) - fm(iL,iY))
   end do
! end of loop over iL
  end do
!--------------------------------------------------------------------
  return
end subroutine SLBdiff
!*********************************************************************


!***********************************************************************
subroutine Solve(model,Lprint,initial,nG,error,delta,iterfbol,fbolOK)
!=======================================================================
! This subroutine solves the continuum radiative transfer problem in
! spherical and planar geometry.                       [Deka,'09, MN'09]
!=======================================================================

  use common
  implicit none

  integer model, error, nG, iterfbol, fbolOK,grid,iY,iL,nY_old,y_incr,imu, &
       iPstar,EtaOK , iP, iZ, nZ
  double precision pstar,taulim, us(npL,npY),comp_fdiff_bol(npY),comp_fdiff_bol1(npY),calc_fdiff(npY), &
       delta, em(npL,npY), omega(npG,npL), iauxl(npL),iauxr(npL), &
       fdsp(npL,npY),fdsm(npL,npY), fdep(npL,npY),fdem(npL,npY),fs(npL,npY), &
       T_ext(npY), accfbol, fbol_max, fbol_min, aux, &
       Udbol(npY), Usbol(npY), fDebol(npY),fDsbol(npY), maxFerr
  logical initial, Lprint
!--------------------------------------------------------------------------

  if (Lprint.and.iX.ne.0) then
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
!------------- loop over bol.flux conservation -------------

10 do while(fbolOK.eq.0)
   iterfbol = iterfbol + 1
   if (Lprint.and.iX.ge.1) then
     write(18,*)'  ',iterfbol,' iteration over fbol'
   end if
   if (iVerb.eq.2.and.Lprint) write(*,'(a14,i3,a20)') &
        ' In Solve: ',iterfbol,' iteration over fbol '
! set grids for the initial optical depth
   call SetGrids(Lprint,pstar,iPstar,error,initial,iterfbol)
! assign number of grid points to nY_old
   nY_old = nY

   if (iVerb.eq.2) &
      write(*,'(a19,i3,a12)') '  Calculating with ', nY,' grid points.'
   if (iX.ge.1.and.Lprint) then
     write(18,'(a19,i3,a12)') '  Calculating with ', nY,' grid points.'
   end if
   if (nY.gt.npY) then
     if (iVerb.eq.2.and.Lprint) write(*,*)' npY needs to be increased!'
     if (Lprint.and.iX.ge.1) then
       write(18,*) ' ********** message from solve *********'
       write(18,*) ' npY has to be increased in file userpar.inc'
       write(18,*) ' ***************************************'
       iWarning = iWarning + 1
       go to 999
     end if
   end if

! solve the radiative transfer problem
   call Rad_Transf(nG,Lprint,initial,pstar,y_incr,us,fs,em,omega,iterfbol,T_ext)
   if (iVerb.eq.2.and.Lprint) write(*,*)' Done with radiative transfer. '

! calculate diffuse flux
! for slab
   if(slb) then
!    find the diffuse scattered flux(fl=1 for scatt. and fl=0 is for emission)
     call SLBdiff(1,omega,TAUslb,T_ext,utot,nL,nY,fds,fdsp,fdsm)
!    find the diffuse emitted flux
     call SLBdiff(0,omega,TAUslb,T_ext,em,nL,nY,fde,fdep,fdem)
!    find bolometric diffuse flux
     call add2(fds,fde,comp_fdiff_bol,nY)
     call add2(fdsp,fdep,fpbol,nY)
     call add2(fdsm,fdem,fmbol,nY)
!    find calculated bolometric diffuse flux
     do iY = 1, nY
       calc_fdiff(iY) = (fsLbol(1)-fsRbol(1)) - (fsLbol(iY)-fsRbol(iY))
       comp_fdiff_bol1(iY) = comp_fdiff_bol(iY) - comp_fdiff_bol(1)
     end do
     comp_fdiff_bol = comp_fdiff_bol1
!  for sphere
   elseif(sph) then
!!** added here, instead of before Analysis [Sep.1'10]
     call add(npY,nY,npL,nL,fs,fds,fde,ftot)
     call add2(fds,fde,comp_fdiff_bol,nY)
     do iY = 1, nY
       calc_fdiff(iY) = abs(fsbol(1) - fsbol(iY))
     end do
   end if

!!!Added by MN to find accfbol as in old Dusty versions.
     call BOLOM(fs,fsbol)
     call BOLOM(fde,fDebol)
     call BOLOM(fds,fDsbol)
     DO iY = 1, nY
        fbol(iY) = fDebol(iY)+fDsbol(iY)+fsbol(iY)
     END DO
     call Finderr(fbol,maxFerr)

     IF(iInn.eq.1 .AND. Lprint) THEN
        write(18,'(a20,1p,e12.3)') 'from Solve: maxFerr=',maxFerr
        IF(SLB) THEN
           write(18,*)'      tr       fbol      fsbol      fdsbol     fdebol'
           DO iY = 1, nY
              write(18,'(1p,5E11.3)') tr(iY),fbol(iY),fsbol(iY),fdsbol(iY),fdebol(iY)
           END DO
        ELSE
           write(18,*)'      Y       fbol      fsbol      fdsbol     fdebol'
           DO iY = 1, nY
              write(18,'(1p,5E11.3)') Y(iY),fbol(iY),fsbol(iY),fdsbol(iY),fdebol(iY)
           END DO
        END IF
     END IF
!!!----------------------

! check for flux conservation, and if there is no conservation
! increase number of grid points

   call Flux_Consv(calc_fdiff,comp_fdiff_bol,fbolOK,error,Lprint)

   if (iVerb.eq.2) &
      write(*,'(a20,i3)') ' After Flux_Cons nY=', nY

! initialize flag, y_incr, which records whether there is an increase in
! y-grid points
   y_incr = 0
! if new number of grid points is greater than the old, then taux = 1
   if (nY.gt.nY_old) y_incr = 1

! if the grid size limit is reached error=2
   if (error.eq.2.and.iterfbol.eq.1) then
! if this is the first calculation end this model
    if(Lprint.and.iX.ge.1) then
     write(18,*)' =========== IMPORTANT WARNING =========== '
     write(18,*)' The limit for grid size is already reached'
     write(18,*)' and flux conservation can not be improved.'
     write(18,*)' Treat all results with caution!'
     write(18,'(a20,1p,F6.2,a2)')'  Achieved accuracy:',maxrat*100., ' %'
    end if
    error = 0
    fbolOK = 2
    iWarning = iWarning + 1
   end if
! if this is a higher iteration use previous solution
   if (error.eq.2) then
    if (Lprint.and.iX.ge.1.and.iterfbol.ge.1) then
     write(18,*)' ======= IMPORTANT WARNING ======== '
     write(18,*)' In trying to conserve fbol, '
     write(18,*)' the limit for grid sizes is reached.  '
     write(18,'(a20,1p,F6.2,a2)')'  Achieved accuracy:',maxrat*100., ' %'
    write(18,*)' Treat all results with caution!'
    end if
    error = 0
    fbolOK = 2
    iWarning = iWarning + 1
   end if
! if fbol not conserved try again with a finer grid
   if (fbolOK.eq.0 .and. iX.ge.1.and.Lprint) then
    write(18,*)'  ******** MESSAGE from Solve ********'
     write(18,'(a20,1p,e12.3,a2)')'  Achieved accuracy:',maxrat*100., ' %'
    write(18,*)'  Trying again with finer grids'
   end if
! if could not conserve fbol in 5 trials give it up
   if (Lprint.and.fbolOK.eq.0 .and. iterfbol.gt.5) then
    if (iX.ge.1) then
     write(18,*)' **********  WARNING from Solve  **********'
     write(18,*)' Could not obtain the required accuracy '
     write(18,'(a4,i4,a8)')' in ',iterfbol,' trials.'
     write(18,'(a20,1p,e12.3,a2)')'  Achieved accuracy:',maxrat*100., ' %'
     write(18,*)' !!  Treat all results with caution  !!'
     write(18,*)' ******************************************'
    end if
    iWarning = iWarning + 1
    fbolOK = 2
    end if
!  end of loop over flux conservion
   end do

   if (rdw.and.sph) then
! counter over eta (for radiatively driven winds only)
! do while (EtaOK.eq.0)
    itereta = itereta + 1
    if (iX.ne.0.and.rdw) then
     write(18,*)'----------------------------------------'
     write(18,*)' ',itereta,'. iteration over eta'
    end if
    if (iVerb.eq.2.and.rdw) &
         write(*,*) ' ',itereta,'. iteration over eta'
! for winds check if eta has converged...
    if ((rdw).and.fbolOK.ne.2) then
! ptr(2) is specified in input and controls converg. crit.
     if (ptr(2).lt.1.0e-6.and.itereta.gt.2)then
      EtaOK = 1
     else
      call Winds(nG,EtaOK)
     end if
     if (itereta.gt.10.and.EtaOK.eq.0) then
      EtaOK = 2
      iWarning = iWarning + 1
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
    if(EtaOK.eq.0) go to 10
   end if

! error=4 means npY not large enough for oblique illumination grid
   if (error.eq.4) then
    call msg(15)
    iWarning = iWarning + 1
    goto 999
   end if

   if(Lprint) then
!    add diffuse flux to the transmitted flux to find total flux
!!** sub Add has to be called before WINDS! ftot is needed there! [MN]
!!**     call add(npY,nY,npL,nL,fs,fds,fde,ftot)
     call bolom(utot,ubol)
     if (sph) then
!     calculate intensity (at the outer edge) if required
      if (iC.ne.0) then
        if (iX.ne.0) write(18,*) 'Calculating intensities.'
        call sph_int(nG,omega,fs)
        if (iVerb.eq.2) write(*,*) ' Done with finding intensity '
      end if
!     if needed convolve intensity with the psf
      if (iPsf.ne.0) then
        call convolve
        if (iVerb.eq.2) write(*,*) ' Done with Convolve. '
      end if
!     if needed find the visibility function
      if (iV.ne.0) then
        call visibili
        if (iVerb.eq.2) write(*,*) 'Done with Visibili.'
      end if
     elseif(slb) then
!      if emerging intensity is required (for left side ill.only):
       if(iC.gt.0) then
         if (iVerb.eq.2) write(*,*)' Calculating intensities. '
         call SLBintensity(omega,em)
         do iL = 1, nL
          do imu = 1, nmu
            SLBintm(imu,iL) = SLBintm(imu,iL)/(4.0d0*pi) ! needed to change later [MD]
            SLBintp(imu,iL) = SLBintp(imu,iL)/(4.0d0*pi) ! needed to change later
            iauxr(iL) = SLBintm(imu,iL)/lambda(iL)
            iauxl(iL) = SLBintp(imu,iL)/lambda(iL)
          end do
        end do
        if (iVerb.eq.2) write(*,*) ' Done with finding intensity. '
       end if
     end if       !end if for sphere or slab

!   analyze the solution and calculate some auxiliary quantities
	call Analysis(model,error,us,T_ext,delta)
    if (iVerb.eq.2) write(*,*) 'Done with Analysis.'
    if (iX.ne.0) then
     write(18,'(a36,F5.1,a2)') '  Accuracy of computed diffuse flux:',maxrat*100.0,' %'
     IF (iInn.eq.1) write(18,'(a36,F5.1,a2)') '  Accuracy of computed bolom. flux: ',maxFerr*100.0,' %'
	 write(18,*)'  ==== SOLVE successfully completed ====='
	 write(18,*)'  ======================================='
    end if
   end if  !end if for printout
!-----------------------------------------------------------------------

999 return
end subroutine solve
!***********************************************************************

!***********************************************************************
subroutine SPH_diff(flag1,initial,iter,iterfbol,T_ext,vec1,omega,vec2)
!=======================================================================
! This subroutine finds the diffuse part for both emission and scattering
! for spherical shell                                        [Deka, 2008]
!! If flag1=1 for vec1=Em; flag1=2 for vec1=Us
!!** Introduced integration of diffuse radiation with Nordlund, as in old Dusty.
!!** Restored ETAzp instead of the time-consuming GauLeg integration
!!** New commments added.                                   [MN, July'10]
!=======================================================================

  use common
  implicit none

  integer iP,iL,iZ,iZz,nZ,iY,iYy,iter,iterfbol, iNloc, flag1,flagN, error
  double precision vec1(npL,npY), vec2(npL,npY), omega(npG,npL), Iplus1(npP,npL), &
       Iplus2(npP,npL), Iminus(npP,npL), T_ext(npY), aux2(npY), diff(npY),        &
	   S_fun(npY), func(npY), faux3(npY), xN(npP), yN(npP), S_loc, p_loc, expow1, &
	   result1, result2, res1

  logical initial
!-----------------------------------------------------------------------

  Iplus1 = 0.0d0
  Iplus2 = 0.0d0
  Iminus = 0.0d0
!!** for each radial grid point calculate the integrals from Sec.4.1 in Blueprint:
  do iY = 1, nY
   do iP = 1, Plast(iY)
     iZz  = iY + 1 - iYfirst(iP)   !this is for z in eq.(4.1.5)
!    upper limit for the counter of z position
     nZ  = nY + 1 - iYfirst(iP)   !nZ is index for zmax=sqrt(Y**2-p**2) [MN]

     do iL = 1, nL
	   do iYy = 1, nY
         if (flag1.eq.1) then
!          find the diffuse emission term, vec1=em
           S_fun(iYy) = (1.0d0-omega(1,iL))*vec1(iL,iYy)*T_ext(iYy)
         elseif (flag1.eq.2) then
!          find the scattering term, vec1=Us initially, or Utot afterwards
           if (initial.and.iter.eq.1.and.iterfbol.eq.1) then
             S_fun(iYy) = omega(1,iL)*vec1(iL,iYy)*T_ext(iYy)
           else
             S_fun(iYy) = omega(1,iL)*utot(iL,iYy)*T_ext(iYy)
           end if
         end if
       end do ! end do over dummy iYy

       if (P(iP).le.1.0d0) then
!        inside the cavity
         do iZ = 1, nZ
           iYy = iYfirst(iP) + iZ - 1
           diff(iZ) = S_fun(iYy)
         end do
       else
!        in the shell
         do iZ = 1, nZ-1
           iYy = iYfirst(iP) + iZ - 1
           if(iZ.eq.1.and.P(iP).gt.Y(iYy).and.P(iP).lt.Y(iYy+1)) then
             p_loc = P(iP)
             call LININTER(npY,nY,Y,S_fun,p_loc,iNloc,S_loc)
             diff(iZ) = abs(S_loc)
           else
             diff(iZ) = S_fun(iYy)
           end if
         end do
         diff(nZ) = S_fun(nY)
       end if

       do iZ = 1, nZ
         aux2(iZ) = tautot(iL)*ETAzp(iP,iZ)
       end do

! 1st term in the energy density or flux. See blueprint, Table 4.1, or eq.(4.1.5)
!     from z0-midpoint to the outer edge of the shell (on the left side)  [MN]
      DO iZ = 1, nZ
	    faux3(iZ) = exp(-aux2(iZ))
        func(iZ) =  diff(iZ)
      END DO
      CALL Simpson(npY,1,nZ,faux3,func,res1)
      Iplus1(iP,iL) = abs(res1)*exp(-aux2(iZz))
! 2nd term in the energy density or flux. See blueprint, Table 4.1, or eq.(4.1.5)
!     from z0-midpoint to the running z-point
      DO iZ = 1, iZz
	    expow1 = tautot(iL)*abs(ETAzp(iP,iZz) - ETAzp(iP,iZ)) !this is tau(z',z;p)
	    faux3(iZ) = exp(-expow1)
        func(iZ) =  diff(iZ)
      END DO
      CALL Simpson(npY,1,iZz,faux3,func,res1)
      Iplus2(iP,iL) = abs(res1)
! 3rd term in the energy density or flux. See blueprint, Table 4.1, or eq.(4.1.5)
!     from the running z-point to the outer edge of the shell [MN]
      DO iZ = iZz, nZ
	    expow1 = tautot(iL)*abs(ETAzp(iP,iZ) - ETAzp(iP,iZz)) !this is tau(z,z';p)
	    faux3(iZ) = exp(-expow1)
        func(iZ) =  diff(iZ)
      END DO
      CALL Simpson(npY,iZz,nZ,faux3,func,res1)
      Iminus(iP,iL) = abs(res1)

	 end do ! end loop over wavelengths
   end do ! end loop over impact parameters P=1..Plast(iY)


!  Find diffuse energy density (U)
   if (moment_loc.eq.1) then
!!**   xN is mu,the integration variable.
!!**   U ~ Int[yN*dmu], while flux ~ Int[yN*mu*dmu]. Sub Nordlund takes care of this difference.
!!**   When calling NORDLUND(flagN,xN,yN,N1,N2,m,intfdx,error)
!!**    flagN=1 is for analytic integration; m=0 for U and m=1 for flux [MN]
       do iL = 1, nL
            DO iP = 1, Plast(iY)
              xN(iP) = sqrt(1.0-(P(iP)/Y(iY)*P(iP)/Y(iY)))
            END DO
!           generate intensity array for NORDLUND
            DO iP = 1, Plast(iY)
              yN(iP) = abs(Iplus1(iP,iL) + Iplus2(iP,iL) + Iminus(iP,iL))
            END DO
            CALL NORDLUND(0,xN,yN,1,nPcav+1,0,result1,error)
!           flag for analytic integration outside cavity
!!            IF (iY.GT.6) THEN ! this was in the old Dusty
            IF (iY.GT.4) THEN  ! take it as 4 in case the grid has a few pts. only.
               flagN = 1
            ELSE
               flagN = 0
            ENDIF
!           angular integration outside cavity
            IF (iY.GT.1) THEN
               CALL NORDLUND(flagN,xN,yN,nPcav+1,Plast(iY),0,result2,error)
               IF (error.NE.0) GOTO 999
            ELSE
              result2 = 0.0
            END IF
!!**        result1 is for inside, result2 is for outside the cavity [MN]
            vec2(iL,iY)= 0.5*(result1 + result2)/T_ext(iY)
       end do  !end do over lambda
! ------------end of U integration ------------------

!  Find diffuse flux
   elseif(moment_loc.eq.2) then
        do iL = 1, nL
            DO iP = 1, Plast(iY)
              xN(iP) = sqrt(1.0-(P(iP)/Y(iY)*P(iP)/Y(iY)))
            END DO
!           generate intensity array for NORDLUND
            DO iP = 1, Plast(iY)
              yN(iP) = abs(Iplus1(iP,iL) + Iplus2(iP,iL) - Iminus(iP,iL))
            END DO
!           integration inside the cavity
            CALL NORDLUND(0,xN,yN,1,nPcav+1,1,result1,error)
!           flag for analytic integration outside cavity
!!            IF (iY.GT.6) THEN ! this was in the old Dusty
            IF (iY.GT.4) THEN
               flagN = 1
            ELSE
               flagN = 0
            ENDIF
!           angular integration outside cavity
            IF (iY.GT.1) THEN
               CALL NORDLUND(flagN,xN,yN,nPcav+1,Plast(iY),1,result2,error)
               IF (error.NE.0) GOTO 999
            ELSE
              result2 = 0.0
            END IF
!!**        result1 is for inside, result2 is for outside the cavity [MN]
            vec2(iL,iY) = 2.0d0*pi*abs(result1+result2)/T_ext(iY)
	   end do  !end do over lambda
   end if !end if for diffuse flux

  end do !end do over radial grid Y(iY)  [MN]

!-----------------------------------------------------------------------

999  return
end subroutine SPH_diff
!***********************************************************************

!***********************************************************************
subroutine SPH_ext_illum(m0,m1,m1p,m1m)
!=======================================================================
! This subroutine finds the transmitted energy density and flux for sphere
! if there is external illumination.                          [Deka, 2008]
!=======================================================================
  use common
  implicit none

  integer i,ii, iP,iL, n,nn,nZ, iz, izz, iY, error,iNloc
  double precision eta,tauaux(npY), z(np,nY), m0(npL,npY), m1(npL,npY), &
       m1m(npL,npY),m1p(npL,npY),result1, result2, term1(npL,npP),term2(npL,npP), &
       term_aux1(npP), x1, x2, p1, dyn2,p_loc,Yloc1, angle(npP), expow1, expow2
  double precision, dimension(:), allocatable:: xg,  wg
  external eta
!-----------------------------------------------------------------------
  dyn2 = 1.0d-30
  error = 0
  term1 = 0.0d0
  term2 = 0.0d0

  do iY = 1, nY
   do iP = 1, Plast(iY)
!   upper limit for the counter of z position
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

  return
end subroutine SPH_ext_illum
!***********************************************************************


!!=======================================================================
! Routines and functions related to the various GRIDS and to initializing
! of the optical properties.                         [MN, Aug.2010]
!!=======================================================================

!***********************************************************************
subroutine find_z(iP,nZ,z)
!=======================================================================
  use common, only: iYfirst,nP,Y, P
  implicit none

  integer iW,iY,iP,iZ, nZ
  double precision w, z(nP,nZ)
! =======================================================================
  z = 0.0d0
! maximal number of points along tangential position, z
  if (P(iP).ge.1.0d0) then
   z(iP,1) = 0.0d0
  else
   z(iP,1) = dsqrt(1.0d0-P(iP)**2.0d0)
  end if

! loop over z
  do iZ = 2, nZ
   iW = iYfirst(iP) + iZ - 1
   z(iP,iZ) = dsqrt(Y(iW)**2.0d0-P(iP)**2.0d0)
  end do
! -----------------------------------------------------------------------
  return
end subroutine find_z
!***********************************************************************

!***********************************************************************
subroutine getOmega(nG,omega)
!=======================================================================
! This subroutine generates albedo omega(iL,iY) from the abs/sca cross-
! sections and the component abundancies. This is temporary (trivial)
! version  for single size grains.                     [Z.I., Mar. 1996]
!
!!** Note that Omega(iG,iL) here is re-defined compared to the old Dusty. [MN]
!=======================================================================

  use common
  implicit none

  integer  iG, nG, iL, iY
  double precision omega(npG,npL)
!-----------------------------------------------------------------------

! generate overall albedo through the envelope
! ** this is for future multigrain code **
! ** for single grains it is trivial **
  do iG = 1, nG
! calculate albedo
   do iL = 1, nL
    omega(iG,iL) = sigmaS(iG,iL)/(sigmaA(iG,iL) + sigmaS(iG,iL))
   end do
! calculate relative abundances
   do iY = 1, nY
    abund(iG,iY) = 1.0d0
   end do
  end do
!----------------------------------------------------------------------
  return
end subroutine getOmega
!**********************************************************************


!**********************************************************************
subroutine getOptPr(nG,nameQ,nameNK,er,stdf)
!=====================================================================
! This subroutine calculates the absorption and scattering efficiences
! Qabs and Qsca in the wavelength range of the code or in case of
! user supplied efficiences reads them from a file.
!                                                 [ZI Mar96; MN Aug97]
!=====================================================================

  use common
  implicit none

  character*235 nameQ(npG), nameNK(10), fname, dummy*132
  integer iG, nG, io1, iL, nLin, iiLaux,Nprop,nA, iiA, iiC,iCuser,er, Nmax, npA
! Nmax is the number of records in user supplied file with opt.prop.
! and npA is the dimension of the array of grain sizes
  parameter (Nmax=10000, npA=100)
! parameter (Nmax=10000, npA=1000)
  double precision aa,bb,cc,lambdain(Nmax),Qain(Nmax),Qsin(Nmax), &
       n(npL),k(npL), aQabs(npA,npL),aQsca(npA,npL), amax,       &
       nsd(npA), a(npA), faux1(npA), faux2(npA), f(npA), int,    &
       ala(Nmax), sigAbs(npA,npL), sigSca(npA,npL), sizedist,    &
       aQa(Nmax), aQs(Nmax),  Cnorm, a3ave, a2ave,               &
       n_int(npL), k_int(npL)
  character stdf(7)*235
! ----------------------------------------------------------------
! this should never change
  nL = npL
  Nprop = 7
!----------------------------------------------------------------

  er = 0
! first check that the user supplied wavelength grid is
! monotonously increasing
  if (top.lt.3) then
! calculate efficiencies from n and k by mie theory
! generate the size array
   if (szds.gt.2) then
    amax = 5.0d0*a2
   else
    amax = a2
   end if
   if (dabs(a1-a2).le.1.0d-03) then
    nA = 1
   else
    nA =50
   end if
! build-up the array of sizes a(nA)
   call getSizes(npA,nA,a1,amax,a)
! evaluate the normalization constant for the size
! distribution nsd(nA)
   do iiA = 1, nA
    nsd(iiA) = sizedist(qsd,a(iiA),szds,a2)
   end do
   call powerint(npA,1,nA,a,nsd,Cnorm)
! find the average grain volume aveV and average grain  eff.
! area aveA (needed in dynamics)
   if(dabs(a1-a2).le.1.d-3) then
    aveV = 4.0d0/3.0d0*pi*a1**3
    aveA = pi*a1**2.0d0
   else
    do iiA = 1, nA
     faux1(iiA)=nsd(iiA)*a(iiA)**3.0d0
    end do
    call powerint(npA,1,nA,a,faux1,a3ave)
    aveV = 4.0d0/3.0d0*pi*a3ave/Cnorm
    do iiA = 1, nA
     faux1(iiA)=nsd(iiA)*a(iiA)**2.0d0
    end do
    call powerint(npA,1,nA,a,faux1,a2ave)
    aveA = pi*a2ave/Cnorm
   end if
!--  loop over supported components --
   do iiC= 1, Nprop
    f(iiC) = xC(iiC)
    fname = stdf(iiC)
    call getprop(npL,lambda,nL,fname,n,k,er)
    if (er.eq.3) goto 999
! calculate qabs and qsca for supported grains
    call mie(npL,nL,lambda,n,k,npA,nA,a,1,aQabs,aQsca)
! for each lambda integrate pi*a^2*qext with n(a)da
    do iL = 1, nL
     do iiA = 1, nA
      faux1(iiA)=nsd(iiA)*aQabs(iiA,iL)*pi*a(iiA)**2.0d0
      faux2(iiA)=nsd(iiA)*aQsca(iiA,iL)*pi*a(iiA)**2.0d0
     end do
     call powerint(npA,1,nA,a,faux1,int)
     sigAbs(iiC,iL) = int/Cnorm
     call powerint(npA,1,nA,a,faux2,int)
     sigSca(iiC,iL) = int/Cnorm
    end do
   end do
   if (top.eq.2) then
!--  loop over user supplied components --
    do iCuser = 1, nfiles
     iiC = Nprop + iCuser
     f(iiC) = xCuser(iCuser)
! read in optical properties
     fname = nameNK(iCuser)
     call getprop(npL,lambda,nL,fname,n,k,er)
     if (er.eq.3) goto 999
! calculate qabs and qsca
     call mie(npL,nL,lambda,n,k,npA,nA,a,1,aQabs,aQsca)
! for each lambda integrate pi*a^2*qext with n(a)da
     do iL = 1, nL
      do iiA = 1, nA
       faux1(iiA)=nsd(iiA)*aQabs(iiA,iL)*pi*a(iiA)**2.0d0
       faux2(iiA)=nsd(iiA)*aQsca(iiA,iL)*pi*a(iiA)**2.0d0
      end do
      call powerint(npA,1,nA,a,faux1,int)
      sigAbs(iiC,iL) = int/Cnorm
      call powerint(npA,1,nA,a,faux2,int)
      sigSca(iiC,iL) = int/Cnorm
     end do
    end do
   else
    nfiles = 0
   end if
! mix them together (syntetic grain model)
   do iL = 1, nL
    sigmaA(1,iL) = 0.0d0
    sigmaS(1,iL) = 0.0d0
    do iiC= 1, Nprop+nfiles
     sigmaA(1,iL) = sigmaA(1,iL) + f(iiC) * sigAbs(iiC,iL)
     sigmaS(1,iL) = sigmaS(1,iL) + f(iiC) * sigSca(iiC,iL)
    end do
   end do

  else
! this is for top.ge.3 - [Sigma/V] from a file
! initialize aveV and aveA for this case
   aveV = 1.0d0
   aveA = 1.0d0
! read in lambda grid and optical properties
   do iG = 1, nG
    open(1,err=998,file=nameQ(iG),status='old')
    read(1,'(a)',err=998)dummy
    read(1,'(a)',err=998)dummy
    read(1,'(a)',err=998)dummy
    iL = 0
    io1 = 0
    do while (io1.ge.0)
     read(1,*,end=900,err=998,iostat=io1) aa, bb, cc
     if (io1.ge.0) then
      iL = iL + 1
      lambdain(iL) = aa
      Qain(iL) = bb
      Qsin(iL) = cc
     end if
    end do
900 close(1)
    if (iL.lt.2) goto 998
    nLin = iL
! if input wavelengths in descending order turn them around
    if (lambdain(1).gt.lambdain(2)) then
     do iL = 1, nLin
      ala(iL) = lambdain(iL)
      aQa(iL) = Qain(iL)
      aQs(iL) = Qsin(iL)
     end do
     do iL = 1, nLin
      lambdain(iL) = ala(nLin+1-iL)
      Qain(iL) = aQa(nLin+1-iL)
      Qsin(iL) = aQs(nLin+1-iL)
     end do
    end if
! interpolate to dusty's wavelength grid
    do iL = 1, nL
     call powerinter(Nmax,nLin,lambdain,Qain,lambda(iL),iiLaux,aa)
     sigmaA(iG,iL) = aa
! call lininter(Nmax,nLin,lambdain,Qsin,lambda(iL),iiLaux,aa)
     call powerinter(Nmax,nLin,lambdain,Qsin,lambda(iL),iiLaux,aa)
     sigmaS(iG,iL) = aa
    end do
   end do
  end if
  goto 999
998 write(12,*)' ***  FATAL ERROR IN DUSTY ***********'
    write(12,*)' The file with optical properties:'
    write(12,'(a2,a70)')'  ',nameQ(iG)
    write(12,*)' is missing or not properly formatted?!'
    write(12,*)' **************************************'
! close(12)
  er = 3
!-----------------------------------------------------------------------

999 return
end subroutine getOptPr
!***********************************************************************


!***********************************************************************
subroutine GetProp(npL,lambda,nL,fname,en,ek,error)
!=======================================================================
! This subroutine reads optical properties en(i,j), ek(i,j) from file
! fname(nf), with i=nf, j=1..nLl(nf), and interpolates them onto
! wavelength grid lambda(1..nL)                        [z.i., mar. 1996]
! =======================================================================
  implicit none
  character*235 fname
  character*232 line
  integer i, nL, iLoc, iL, npL, io1, error, Nmax
! Nmax is the number of records in the user supplied file
  parameter (Nmax=10000)
  double precision en(npL), ek(npL), lambda(npL), pw(Nmax),   &
       pren(Nmax), pimn(Nmax), a(Nmax), b(Nmax), c(Nmax), aa, &
       bb, cc
!-----------------------------------------------------------------------

  error = 0
  open(2,err=998,file=fname,status='old')
! read in a header from the input file
  do i = 1, 7
   read(2,'(a)',err=998)line
  end do
! read in input data
  iL = 0
  io1 = 0
  do while (io1.ge.0)
   read(2,*,end=900,err=998,iostat=io1) aa, bb, cc
   if (io1.ge.0) then
    iL = iL + 1
    pw(iL) = aa
    pren(iL) = bb
    pimn(iL) = cc
   end if
  end do
900 close(2)
! if input wavelengths in descending order turn them around
  if (iL.lt.2) goto 998
  if (pw(1).gt.pw(2)) then
   do i = 1, iL
    a(i) = pw(i)
    b(i) = pren(i)
    c(i) = pimn(i)
   end do
   do i = 1, iL
    pw(i) = a(iL+1-i)
    pren(i) = b(iL+1-i)
    pimn(i) = c(iL+1-i)
   end do
  end if
! interpolate
  do i = 1, nL
   call lininter(Nmax,iL,pw,pren,lambda(i),iLoc,en(i))
   call lininter(Nmax,iL,pw,pimn,lambda(i),iLoc,ek(i))
  end do
  goto 999
998 write(12,*)' ***  fatal error in dusty  ***********'
  write(12,*)' file with optical properties:'
  write(12,'(a2,a70)')'  ',fname
  write(12,*)' is missing or not properly formatted?!'
  write(12,*)' **************************************'
! close(12)
  error = 3
!-----------------------------------------------------------------------
999 return
end subroutine GetProp
!***********************************************************************

!***********************************************************************
subroutine getSizes(nn,n,x1,x2,x)
!=======================================================================
! This subroutine generates an array x(i=1..N) of physical size NN,
! with N elements logarithmically spaced between x1 and x2.
!                                              [ZI,Aug'96;MN,Nov'97]
!=======================================================================

 implicit none
  integer nn, n, i
  double precision x(nn), x1, x2, fac, pw1, pw2
!-----------------------------------------------------------------------

  if (n.gt.1) then
   pw1 = 1.0d0/(n-1)
   fac = (x2/x1)**pw1
   do i = 1, n
    pw2 = 1.0d0*(i-1)
    x(i) = x1*fac **pw2
   end do
  else
   x(1) = x1
  end if
!-----------------------------------------------------------------------

  return
end subroutine getSizes
!***********************************************************************

!***********************************************************************
subroutine GetTau(nG,tau1,tau2,tauIn,Nrec,GridType,Nmodel,tau)
!=======================================================================
! This subroutine generates total optical depth TAUtot.
!                                                      [Z.I., Mar. 1996]
!=======================================================================

  use common
  implicit none

  integer model, Nmodel, nG, iL, GridType, Nrec
  double precision  q,TAUin(Nrec), tau(Nmodel), TAU1, TAU2
!-----------------------------------------------------------------------

  if (ng.gt.1) then
   write(18,*)'fix GetTau, ng>1 !'
   stop
  end if

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

!**********************************************************************
subroutine GetTau0(tau0,i)
!=======================================================================
! This subroutine finds tau0 - the running variable for optical depth of
! the set of models. For the first iteration this is init_tau from the .inp file.
!                                                      [MN, May'10]
!=======================================================================

  use common
  implicit none

  integer i
  double precision tau0, tauaux(npY), deltaumax
!-----------------------------------------------------------------------
!   find tau0 for each loop over increasing optical
    if (i.eq.1) then
!     for initial grid generation
      if(sph) tau0 = init_tau
      if(slb) tau0 = abs(init_tau*mu1)
!!!    else
!     for smaller tau-steps if solution not obtained
!!!      if(sph) tau0 = dabs(init_tau)*2.0d0**(1.0d0/dfloat(i))
!!!      if(slb) tau0 = dabs(init_tau*mu1)*2.0d0**(1.0d0/dfloat(i))
    end if
!-----------------------------------------------------------------------

  return
end subroutine GetTau0
!***********************************************************************


!***********************************************************************
subroutine GetTauMax(tau0,nG)
!=======================================================================
! This subroutine finds tautot(nL) and its max. value taumax.
!                                                      [MN, May'10]
!=======================================================================

  use common
  implicit none

  integer iL, nG
  double precision faux1(npL), faux2(npL), tau0, sigAfid, sigSfid
!-----------------------------------------------------------------------
   if (nG.gt.1) then
     write(18,*)'nG>1 is not implemented yet!'
     stop
   end if

   taumax = 0.0d0
   taufid0 = tau0
   do iL = 1, nL
    faux1(iL) = sigmaA(nG,iL)
    faux2(iL) = sigmaS(nG,iL)
   end do
   if (lamfid.lt.lambda(1)) then
    write(12,*)' fiducial wavelength was too small.'
    write(12,'(a8,e9.3,a17)')' using ',lambda(1),' micron instead.'
   end if
   if (lamfid.gt.lambda(nL)) then
    write(12,*)' fiducial wavelength was too large.'
    write(12,'(a8,e9.3,a17)')' using ',lambda(nL),' micron instead.'
   end if
   call lininter(npL,nL,lambda,faux1,lamfid,iLfid,sigAfid)
   call lininter(npL,nL,lambda,faux2,lamfid,iLfid,sigSfid)

! Extinction efficiency at fiducial wavelength
   sigexfid = sigAfid + sigSfid
   do iL = 1, nL
    tautot(iL) = taufid0*(SigmaA(nG,iL) + SigmaS(nG,iL))/SigExfid
    if (tautot(iL).ge.taumax) taumax = tautot(iL)
   end do
!-----------------------------------------------------------------------

  return
end subroutine GetTauMax
!***********************************************************************


!***********************************************************************
double precision function SizeDist(q,aa,sdtype,a0)
!=======================================================================
! This subroutine calculates size distribution n(a) for a=aa. the size
! distribution is mrn type n(a)~1/a**q for sdtype.le.2 and kmh type
! n(a)~dexp(-a/a0)/a**q otherwise
!                                                      [Z.I., Aug. 1996]
!=======================================================================

  implicit none
  integer sdtype
  double precision aa, a0, q
!-----------------------------------------------------------------------

  if (sdtype.le.2) then
   sizedist = 1.0d0/aa**q
  else
   sizedist = dexp(-aa/a0)/aa**q
  end if
!-----------------------------------------------------------------------

  return
end function SizeDist
!***********************************************************************

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


!***********************************************************************
subroutine SetGrids(Lprint,pstar,iPstar,error,initial,iterfbol)
!=======================================================================
! Sets the y and p grids for sphere and tau-grid for slab.
!                                                     [MN & ZI, July'96]
!=======================================================================
  use common
  implicit none

  integer iL,iY,error, iPstar, iterfbol,taux
  double precision pstar,tau
  logical initial, Lprint
!-----------------------------------------------------------------------

  if (sph) then
     if (initial.and.iterfbol.eq.1) then
! generate y-grid
        call Ygrid(pstar,iPstar,error)
! generate p and z grid
        call Pgrid(pstar,iPstar,error)
     else
! redefine p and z grid for the new y-grid
	call Pgrid(pstar,iPstar,error)
     end if
  elseif (slb) then
     call SLBtau(initial,iterfbol,tautot,y,nL,nY,TAUslb)
  end if
  if (Lprint.and.iX.ge.1) then
   if (slb) then
!!     write(18,'(a21)')'  Tau grid generated.'
     write(18,'(a25,i4)')'  Tau grid generated, nY=', nY
   elseif (sph) then
    write(18,'(a24,i3)')' y grid generated, nY = ',nY
    write(18,'(a24,i3)')'                   nP = ',nP
    write(18,'(a24,i3)')'                 Nins = ',Nins
    write(18,'(a24,i3)')'                 Ncav = ',Ncav
   end if
  end if
!-----------------------------------------------------------------------

  return
end subroutine setGrids
!***********************************************************************

!***********************************************************************
subroutine Ygrid(pstar,iPstar,error)
!=======================================================================
! This subroutine generates the radial (Y) grid. Yout is the relative
! thickness, Yout=rout/r1. pstar is the impact parameter for star (0<=p<=1).
! First few points are prescribed arbitrarily.  This subroutine calls
! function ETA which evaluates the normalized density profile.
!                                       [ZI, Nov'95; MN,Sep'99, Deka'08]
!=======================================================================

  use common
  implicit none
  integer i,j, jy,itr,istop,error,iPstar, iter, iYdummy, iY, ir,irmax
  double precision pstar, dyR, dyT, eta, aux, EtaTemp(npY), &
       ee, Yloc,  Yold, Ynew, y1, delY1, ymid,      &
       rat, dY, etamin, etamax,deltau, delY,tauloc2,    &
       tausc, nh, fac, dyF, EtaRat, TAUloc, delTausc, facc
  external eta
!-----------------------------------------------------------------------
!!** these are initialized as in old Dusty's sub Input. [MN]
    delTAUsc = 0.3
    facc = 2.0
    EtaRat = 4.0
!!**

  if(.not.fild) then
! save old grid and values of eta (important for denstyp = 5 or 6)
  if (nY.gt.0.and.(rdw)) then
   Yprev  = Y
   etatemp = etadiscr
   nYprev = nY
  end if

!     max number iter. over improving ratio of two Eta's
      irmax = 20
!     save old grid and values of Eta (important for denstyp = 5 or 6)
      IF (nY.GT.0.AND.(RDW)) THEN
        DO iY = 1, nY
          Yprev(iY) = Y(iY)
          EtaTemp(iY) = ETAdiscr(iY)
        END DO
        nYprev = nY
      END IF
      pstar = pstar
      iPstar = iPstar
      y1 = 1.0
      iter = 0
101   error = 0
      iter = iter + 1
!     resolve inner boundary in tau space
!     from requiring TAU(2) = TAUmax*ETA(1)*(Y(2)-1) =~ 1
      Y(1) = y1
      IF (TAUmax*ETA(y1).GT.2000.0) THEN
        delY1 = 2.0 / TAUmax / ETA(y1)
      ELSE
        delY1 = 1.0 / TAUmax / ETA(y1)
      END IF
!     if very thin generate at least about 10-15 pts.
      IF (delY1.GT.0.01/ETA(y1)) delY1 = 0.01/ETA(y1)
!     do not push it over the limit of spline stability
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
        TAUloc = TAUmax * (Y(i)-1.0) * ETA(Ymid)
        IF (Y(i).GE.1.01.OR.TAUloc.GT.10.0) istop = 1
!       in case of shell thinner than 1.01 comment the above line and
!       uncomment the next line, unless it is a case of RDW at high TauV.
!       These require more points near the origin and 1.01 is a better limit.
!       IF (Y(i).GE.1.0001.OR.TAUloc.GT.10.0) istop = 1
      END DO
      Ynew = Y(i)
!     some rule of thumb estimates for factor nh
      IF (TAUmax.GT.10000.0) THEN
!       extreme taus
        nh = 15.0
        Ncav = 80
      ELSE IF (TAUmax.GT.2000.0) THEN
!         huge taus
          nh = 10.0
          Ncav = 40
      ELSE IF (TAUmax.GT.400.0) THEN
!         large taus
          nh = 8.0
          Ncav = 20
      ELSE
!         normal' taus
          nh = 5.0
          Ncav = 10
      END IF
!     very small taus
      IF (TAUmax.LT.10.0) nh = nh / 2.0
      IF (TAUmax.LT.1.0) nh = nh / 2.0
!     empirically: ~1/r needs more points for small y:
      IF ((pow-1.4)*(pow-0.6).LE.0.0) nh = nh * 1.5
!     same for for steep density distributions:
      IF (RDWA.OR.RDW) nh = nh * 1.5
      tausc = ETA(Y(i))+ETA(Y(i-1))
      tausc = (Y(i)-Y(i-1)) * 0.5 * tausc
      fac = dexp(-dlog(tausc) / nh)
      istop = 0
!     for broken power-laws
      IF (Ntr.GE.1) THEN
        itr = 1
        DO j = 1, i
          IF (Y(j).GE.Ytr(itr)) itr = itr + 1
        END DO
      END IF
!     generate the rest of Y grid points
      DO WHILE (istop.NE.1)
        i = i + 1
        Yold = Ynew
!       find maximal increase in Y allowed by the ratio facc
        dyR = Yold * (facc-1.)
!       find maximal increase in Y allowed by delTausc
        dyT = delTausc / ETA(Yold)
!       find maximal increase in Y allowed by the ratio of tausc
        dyF = tausc*(fac-1.) / ETA(Yold)
!       find new Y
!        Ynew = Yold + MIN(dyR,dyT,dyF)
        dY = MIN(dyR,dyT,dyF)
!       Check if the max ratio btw. two Eta values is less than EtaRat
!       and insert additional y-pts. where necessary. This prevents sharp
!       drops in Utot(npL,npY) in case of steep Eta's [MN'99].
        DO ir = 1 , irmax
          Ynew = Yold + dY
          rat = ETA(Yold)/ETA(Ynew)
          IF(rat.GE.1./EtaRat .AND. rat.LE.EtaRat) goto 10
          dY = 0.5*dY
        END DO
        CALL MSG(16)
10      continue
        Y(i) = Ynew
!       make sure that all transition points are included in the grid
        IF (Ntr.GE.1) THEN
         IF (Y(i).GE.Ytr(itr)) THEN
           Y(i) = Ytr(itr)
           Ynew = Y(i)
           itr = itr + 1
         END IF
        END IF
        aux = ETA(Ynew)+ETA(Yold)
        aux = (Ynew-Yold) * 0.5 * aux
        tausc = tausc + aux
!       finish when Yout is reached
        IF (Ynew.GE.Yout) istop = 1
      END DO
      Y(i) = Yout
      nY = i
!     insert additional penultimate point to avoid spline oscillations
      Y(nY+1) = Yout
      Y(nY) = dsqrt(Y(nY)*Y(nY-1))
      nY = nY + 1
!     check that outer edge is well resolved in tau space
!     (important for flat ETAs)
      istop = 0
      DO WHILE (istop.NE.1)
        IF ((Yout-Y(nY-1))*TAUmax*ETA(Yout).GT.1.0) THEN
          Y(nY+1) = Yout
          Y(nY) = dsqrt(Y(nY)*Y(nY-1))
          nY = nY + 1
        ELSE
          istop = 1
        END IF
      END DO
!     check dynamical range of Eta to avoid nonphysical results or code errors [
      Etamax = 0.
      Etamin = 1.e+20
      DO iY = 1, nY
       IF(ETA(Y(iY)).lt.Etamin) Etamin = ETA(Y(iY))
       IF(ETA(Y(iY)).gt.Etamax) Etamax = ETA(Y(iY))
       IF (ETA(Y(iY)).LT.1.e-12) THEN
        IF (iX.GT.0) THEN
         write(18,*)'      Y          ETA  '
         DO jY = 1, nY
           write(18,'(1p,2e12.3)') Y(jY),ETA(Y(jY))
         END DO
        END IF
        CALL MSG(17)
        error = 6
        iERROR = iERROR + 1
        goto 102
       END IF
      END DO
      IF ((Etamin/Etamax).LT.1.e-12) THEN
       CALL MSG(18)
       error = 6
       iERROR = iERROR + 1
       goto 102
      END IF
!     check that the Y grid is not too large
      IF (nY.GT.npY) THEN
        delTAUsc = delTAUsc * 1.5
        iWARNING = iWARNING + 1
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
!     intepolate ETAdiscr to new Y grid (for RDW (denstyp=5 or 6))
!      write(18,*)' ***** From Ygrid *****'
!      write(18,*)'      Y      ETAdiscr      '
      DO iY = 1, nY
        Yloc = Y(iY)
        IF (iterETA.GT.1) THEN
          CALL LinInter(npY,nYprev,Yprev,EtaTemp,Yloc,iYdummy,ee)
          ETAdiscr(iY) = ee
        ELSE
          ETAdiscr(iY) = ETA(Yloc)
        END IF
!        write(18,'(1p,2e12.4)')Y(iY), ETAdiscr(iY)
      END DO
!      write(18,*)' *************************'


 elseif(fild) then
! take the radial grid as in density profile file [MN,Apr'02]
  nY = nYEtaf
  do iY = 1, nY
   Y(iY) = yEtaf(iY)
  end do
 end if
!-----------------------------------------------------------------------

102 return
end subroutine ygrid
!***********************************************************************

!***********************************************************************
subroutine Pgrid(pstar,iPstar,error)
!=======================================================================
! after having the Y grid, generate the P grid (impact parameters)
!                                                              [Deka'08]
!=======================================================================

  use common
  implicit none

  integer i,ii,k,iP,iz,iw,nZ,j,error, Naux, iPstar, istop, NinsLoc
  double precision pstar, delP, eta
  external eta
!-----------------------------------------------------------------------

  error  = 0
! Ncav = # rays in cavitiy
  Ncav = 30

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
    Nins = 2 !20
   elseif((Y(i+1) - Y(i)).gt.10.0d0.and.(Y(i+1) - Y(i)).le.20.d0) then
    Nins = 2 !15
   elseif ((Y(i+1) - Y(i)).gt.1.5d0.and.(Y(i+1) - Y(i)).le.10.d0) then
    Nins = 2 !10
   elseif ((Y(i+1) - Y(i)).gt.1.0d0.and.(Y(i+1) - Y(i)).le.1.5d0) then
    Nins = 3
   elseif ((Y(i+1) - Y(i)).gt.0.5d0.and.(Y(i+1) - Y(i)).le.1.0d0) then
    Nins = 7
   else
    Nins = 10
   end if
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
   iError = iError + 1
  end if
!-----------------------------------------------------------------------

  return
end subroutine pgrid
!***********************************************************************


! **********************************************************************
  subroutine SLBtau(initial,iterfbol,tautot,Y,nL,nY,TAUslb)
! ======================================================================
! It generates TAUslb(npL,npY) grid with given spacing delT, calculated
! in sub SLBy. In the current version it is exp near the slab boundaries
! and equidistant in the middle.                  [Deka,'09, MN,'98,'09]
! ======================================================================

  implicit none
  integer iterfbol, nY, iY, nL, iL
  integer npY, npP, npX, npL, npG, npR
  include 'userpar.inc'
  parameter (npG=1)
  double precision TAUslb(npL,npY),tautot(npL),tau(npY),Y(npY)
  logical initial
! ----------------------------------------------------------------------

  if (initial.and.iterfbol.eq.1) then
!!   nY = 6
!!   do iY = 1, nY
!!    if(iY.le.3) then
!!     Y(iY) = 0.5d0*(dexp(0.5d0*iY)-dexp(0.5d0))
!!    else if(iY.gt.3.and.iY.le.5) then
!!     Y(iY) = Y(iY-1)
!!    else
!!     Y(iY) = Y(2+nY-iY)
!!    end if
!!   end do
!  Use more points initially. [MN]
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
!     nY = 30
!     DO iY = 1, nY
!        IF(iY.LE.10) THEN
!           Y(iY) = 0.5*(exp(0.5*iY)-exp(0.5))
!        ELSE IF(iY.GT.10.and.iY.LE.22) THEN
!           Y(iY) = Y(iY-1)
!        ELSE
!           Y(iY) = Y(2+nY-iY)
!        END IF
!     END DO
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
!----------------------------------------------------------------------

  return
end subroutine SLBtau
!***********************************************************************

!***********************************************************************
SUBROUTINE getETAzp
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

      INTEGER iP, nZ, iZ, iW
      DOUBLE PRECISION IntETA, auxEta, w1, w2
      EXTERNAL IntEta
! -----------------------------------------------------------------------
!     loop over impact parameters
      DO iP = 1, nP
!       maximal number of points along tangential position, z
        nZ = nY + 1 - iYfirst(iP)
!       starting values for z and ETAzp(iP,iZ)
        IF (P(iP).GE.1.0) THEN
          w2 = P(iP)
          ELSE
          w2 = 1.0
        END IF
!       initialize ETAzp(iP,iZ)*TAUtot(iL)
        ETAzp(iP,1) = 0.0
!       loop over z
        DO iZ = 2, nZ
!         index for local radius, w2
          iW = iYfirst(iP) + iZ - 1
!         limits for integration
          w1 = w2
          w2 = Y(iW)
!           find next step in ETAzp
            auxEta = IntETA(P(iP),iW-1,w1,w2)
!           add next step in ETAzp
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
!     integrals calculated by MAPLE
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
      END
!***********************************************************************

!***********************************************************************
SUBROUTINE setupETA
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

      INTEGER iY, iCoeff
      DOUBLE PRECISION coef(npY,4), ETA, maxerr, Ymid, Yinverse(npY), &
            ETAaux(npY), ETAmid(npY)
! -----------------------------------------------------------------------
!       generate input function for SPLINE2
        DO iY = 1, nY
          Yinverse(iY) = 1. / Y(iY)
          ETAaux(iY) = ETA(Y(iY))
          IF (iY.LT.nY) THEN
            Ymid = dsqrt(Y(iY)*Y(iY+1))
            ETAmid(iY) = ETA(Ymid)
          END IF
        END DO
!       calculate spline coefficients
        CALL SPLINE2(Yinverse,ETAaux,nY,coef)
!       check and fix spline coefficients
        maxerr = 0.1
!       RDW is initialized in Input
        CALL CHKSPLIN(Yinverse,ETAaux,ETAmid,nY,coef,maxerr)
!       copy coefficients to the output array ETAcoef
        DO iY = 1, nY
          DO iCoeff = 1, 4
            ETAcoef(iY,iCoeff) = coef(iY,iCoeff)
          END DO
        END DO
! -----------------------------------------------------------------------
      RETURN
END SUBROUTINE setupETA
!***********************************************************************

!=======================================================================
! Routines and funcitons related to ANALYSIS of the obtained solution,
! arranged in alphabetical order.                        [MN, Aug.2010]
!=======================================================================

!***********************************************************************
 subroutine Analysis(model,error,us,T_ext,delta)
!=======================================================================
! This subroutine analyzes the solution. It finds the flux conservation
! accuracy and evaluates many output quantites (e.g. QF(y), TAUF(y),Psi, F1
! the rad.pressure force, dynamical quantities etc.)
! This is with new additions acc. to IE'00           [ZI,Mar'96;MN,Mar'99]
!=======================================================================

  use common
  implicit none

  integer i, iL, iY, nn, model, iP, error
  double precision eta, qpTd(npG,npY), qpstar(npY), psi,&
       qaux(npL), qaux2(npL), resaux, xP, Planck, qutot1,       &
       eps1, aux, C1, C2, C3, theta1_loc, ugas_out, s4, mx,     &
       tauV, Gie2000, K1(npY), K2(npY), tauaux(npL,npY),       &
       delta, us(npL,npY), x1, x2, result1, T_ext(npY), maxFerr
!!**  double precision, dimension(:), allocatable:: xg, wg
  external eta
!-----------------------------------------------------------------------

! spectrum (flux at the outer edge as a function of wavelength)
  do iL = 1, nL
   spectrum(iL) = dabs(ftot(iL,nY))
! to prevent taking log from zero in spectral [MN]:
   if (spectrum(iL).le.1.0d-20) spectrum(iL) = 1.0d-20
  end do
!-------------
! analyze bolometric flux error (1/2 of the max spread of fbol)
  CALL FindErr(fbol,maxFerr)

! find the flux averaged optical depth, tauF(y)
  if (sph) then
!   for spherical shell
    tauF(1) = 0.0
    DO iY = 2, nY
!    generate auxiliary function for integration:
!    loop over iL (wavelength)
!    N.B. the definition: ETAzp(1,y) = taur(y)/tauT so that
!    tau(iL,iY) = TAUtot(iL)*ETAzp(1,iY)
     DO iL = 1, nL
        qaux(iL)=TAUtot(iL)*ETAzp(1,iY)*dabs(ftot(iL,iY))/lambda(iL)
     END DO
     CALL Simpson(npL,1,nL,lambda,qaux,resaux)
!    tauF(iY) = <tau(iL,iY)*ftot(iL,iY)>
     tauF(iY) = resaux
    END DO
!   for full RDW calculation redo tauF to be consistent with CalcEta
    IF (RDW) THEN
!     generate ETA and its integral (normalization constant)
      DO iY = 1, nY
         K1(iY) = vrat(1,iY)/ugas(iY)/Y(iY)/Y(iY)
      END DO
      CALL SIMPSON(npY,1,nY,Y,K1,resaux)
!     find tauF
      DO iY = 1, nY
         K2(iY) = qF(iY)*K1(iY)/resaux
         CALL SIMPSON(npY,1,iY,Y,K2,aux)
         tauF(iY) = TAUfid*aux
      END DO
    END IF
  elseif(slb) then
!  for slab
   tauF(1) = 0.0d0
   do iY = 1, nY
!  generate auxiliary function for integration:
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
!   in case of sigma's from a file aveV=1 (initialized in GetOptPr)
    DO iY = 1, nY
      DO iL = 1, nL
       qaux(iL)=(SigmaA(1,iL)+SigmaS(1,iL))/aveV*dabs(ftot(iL,iY))/lambda(iL)
      END DO
      CALL Simpson(npL,1,nL,lambda,qaux,resaux)
      rg(1,iY) = s4 * resaux / r_gd
!     If dust drift (dynamics case):
      IF (RDW) rg(1,iY) = rg(1,iY)*vrat(1,iY)
       IF (iY.EQ.1) THEN
           Phi = resaux
       END IF
    END DO
!   the terminal value of the reddening profile, normalized to y=1
    Phi = resaux / Phi
  END IF
!-------------
! Find the Planck averaged absorption efficiencies
  DO iY = 1, nY
!    generate auxiliary function for integration over wavelengths:
     DO iL = 1, nL
       qaux(iL) = SigmaA(1,iL) * Us(iL,iY) / lambda(iL)
       xP = 14400.0 / Td(1,iY) / lambda(iL)
       qaux2(iL) = SigmaA(1,iL) * Planck(xP) / lambda (iL)
     END DO
     CALL Simpson(npL,1,nL,lambda,qaux,resaux)
     QpStar(iY) = resaux
     CALL Simpson(npL,1,nL,lambda,qaux2,resaux)
     QpTd(1,iY) = resaux
  END DO
! ----------
!     find parameter Psi (see Ivezic & Elitzur, 1996)
!     generate auxiliary function for integration:
!     loop over iL (wavelength)
      DO iL = 1, nL
        qaux(iL) = SigmaA(1,iL) * Utot(iL,1) / lambda (iL)
      END DO
      CALL Simpson(npL,1,nL,lambda,qaux,resaux)
      QUtot1 = resaux
      Psi = QUtot1 / QpTd(1,1)
!     for slab Psi is defined by the flux at normal ill.
      IF (SLB) Psi = dabs(mu1)*QUtot1 / QpTd(1,1)
! -------------
      IF(sph) THEN
!      ratio r1/r* (see Ivezic & Elitzur, 1996, eq. 27)
       r1rs = 0.5 * dsqrt(Psi) * (Tstar(1) / Td(1,1))**2.0
       IF(Left.eq.0) r1rs = 1.0
      END IF
!-------------
!     Find epsilon - the relative contribution of the diffuse radiation
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
!     store these parameters in the storage array
      SmC(1,model) = Psi
      SmC(2,model) = Eps1
      SmC(3,model) = QpStar(1)
      SmC(4,model) = QpTd(1,1)
!!**      SmC(5,model) = maxFerr  !maxFerr is the asymmetry of fbol
      SmC(5,model) = maxrat  !maxrat is the max err of calculated diffuse flux as in Blueprint.
!-------------
! additional output quantities
! bolometric flux at r1 (in W/m2)
  if (typentry(1).eq.5) then
   Fi = 4.0d0*sigma*T_ext(1)
  else
   if(left.eq.0) then
    Fi = 10.0d-20
   elseif(left.eq.1.and.typentry(1).ne.5) then
    Fi = Fint
   end if
  end if

! inner radius (in cm) in case it is not an input
! 5.53e16 = sqrt(10^4*Lo/4/pi)
  if (typentry(1).ne.2) then
   if(slb) then
! r1 is found from Fi = L/(4*pi*r1^2). since in sub input
! Fi=Fi*mu1, here the mu1 dependence has to be removed
    Cr1 =  5.53d+16 / dsqrt(Fi/abs(mu1))
   else
    Cr1 = 5.53d+16 / dsqrt(Fi)
   end if
  end if

  if (sph) then
!  angular diameter of inner cavity if Fbol=1d-6 W/m2
   theta1_loc = 412.6d0 / dsqrt(Fi)
!  check if the pt.source assumption is still obeyed
! (This is only for BB-type spectrum including Engelke-Marengo function)
   if(startyp(1).eq.1.or.startyp(1).eq.2) then
    mx = sqrt(sqrt(Fi/sigma))
    Te_min = 2.0d0 * dmax1(Td(1,1), mx)
   end if
  end if
  if (slb) then
! Teff for the left illuminating source in slab geometry
! Teff = (Fi/sigma)^0.25d0
   SmC(7,model) = (Fi/sigma)**0.25d0
   if (ksi.gt.0.) then
!  Teff for the right illuminating source in slab geometry
    SmC(8,model) = SmC(7,model)*sqrt(sqrt(ksi))
   else
    SmC(8,model) = 0.0d0
   end if
  end if
! calculate conversion constants for dynamics
!!  check the output, not tested for rdw case [mn, 08]
  if (rdwa.or.rdw) then
! for analytical approximation
! (otherwise i1,2,3 are found in gammafun)
!!**   if (rdwa) then
!!** we need these for both RDW and RDWA output! [MN]
    I1 = 2.0d0 * (1.0d0-pow)/(1.0d0+pow)/tauF(nY)
    I2 = I1
    I3 = I1 * tauF(nY) / taufid
    gamma(nY) = 0.5d0
!!**   end if

! terminal expansion velocity, full formula:
   ugas_out = tauF(nY) * (1.0d0-gamma(nY)) / (1.0d0-pow)
! the coefficients come from the units conversion
   C1 = 0.2845d0*taufid*sqrt(psi)/I2/(sigexfid/aveV)*1.0d06/Td(1,1)/Td(1,1)
   C2 = 2.040d0*ugas_out
   C3 = 6.628d0*I3*sigexfid/aveV*gamma(nY)/I1

! from version 2.0 stellar mass is defined as the maximal stellar
! mass which does not quench the wind; the calculation is done
! with half that mass since any smaller mass will have no effect
! on the radial velocity and density profile (see IE2000)
! N.B. erroneous gamma(nY) is removed
   CM = 6.628d0*I3*sigExfid/aveV/I1
!  new definitions for output
!  mass-loss rate in Msol/yr Td(iG,iY) = Td(iG,1)
   CMdot = 1.0d-05 * sqrt(C1)
! terminal expansion velocity in km/s
   Cve = 10.0d0* C2 / sqrt(C1)
!*** this is conversion to the nomenclature as in IE2001
  if (rdwpr) then
! if (rdw) then
! size averaged extinction efficiency
   QV = sigexfid / aveA
   tauV = taufid
   q_star = qF(1)
   zeta1 = vrat(1,1)
   G1 = gamma(1)
   Ginf = gamma(nY)
   if (g1.gt.0.0d0) then
    Gie2000 = 1.0d0 / zeta1 / G1
    delta = 1.0d0 / (Gie2000 - 1.0d0)
   else
    delta = 0.0d0
   end if
   PIrdw = tauV / QV
   Prdw = dsqrt(2.d0*PIrdw/I2/QV/q_star)
   winf = ugas_out / QV / q_star
  end if
 end if
!-----------------------------------------------------------------------

 return
end subroutine Analysis
!***********************************************************************


! ***********************************************************************
      SUBROUTINE Convolve
! =======================================================================
! This subroutine convolves intensity IntOut with the point spread
! function to produce convolved images ConvInt. The work horse is
! subroutine Conv2D, and this subroutine is used to prepare everything.
!                                                      [Z.I., Jan. 1997]
! Changed normalization of convolved images: instead of profiles normalized
! at the center, now the convolved intensities are normalized by the area
! A = 2*pi*Int{psf(x)*x}dx                             [MN, 2005]
! =======================================================================
        use common
      IMPLICIT none

      INTEGER i, j
      DOUBLE PRECISION yang(npP+2), Youtang, deltaOff,    &
           Int1D(npP+2), ConvS, Conv(npP+2), FWHM1max, FWHM2max, PSFN, &
           PSFN1, psf_Y(1000), psf_X(1000), res
! -----------------------------------------------------------------------
!     find the largest FWHMs
      FWHM1max = FWHM1(1)
      FWHM2max = FWHM2(1)
      IF (psftype.LT.3) THEN
        DO i = 1, NlambdaOut
          IF (FWHM1(i).GT.FWHM1max) FWHM1max = FWHM1(i)
          IF (FWHM2(i).GT.FWHM2max) FWHM2max = FWHM2(i)
        END DO
      END IF
!     scale angular coordinate to theta1
      DO i = 1, nP+2
        yang(i) = bOut(i) * Theta1 / 2.0D+00
      END DO
!     generate off-set grid
      Youtang = Y(nY) * Theta1
      IF (Youtang.GT.FWHM1max.AND.Youtang.GT.FWHM2max) THEN
!       the envelope is well resolved, take impact parameter grid
        Nconv = nP + 2
        DO i = 1, Nconv
          Offset(i) = yang(i)
        END DO
      ELSE
!       the envelope is not well resolved, take equidistant grid
!       to 2FWHM1max, i.e. image will be more or less the PSF itself
        Nconv = 30
        deltaOff = 2.0D+00 * FWHM1max / (Nconv-1)
        IF (FWHM2max.GT.FWHM1max) THEN
          deltaOff = 2.0D+00 *FWHM2max / (Nconv-1)
        END IF
        DO i = 1, Nconv
          Offset(i) = deltaOff * 1.0D+00*(i-1)
        END DO
      END IF

! !!!!!!!added normalization of psf, PSFN1(x) is the original fn
!     psfArea(iLambda) b/c the Gaussian option allows FWHM(iLambda)
      write(12,'(a36)')'  lambdaOut(mic)   psfArea(arcsec^2)'
      DO iLambda = 1, NlambdaOut
        DO i = 1, Nconv
          psf_X(i) = Offset(i)
          psf_Y(i) = PSFN(Offset(i))
          CALL ChkRange(dynrange,psf_Y(i))
        END DO
!        CALL ScaleTo1(1000,Nconv,psf_Y)
        CALL ScaletoArea(1000,Nconv,psf_X,psf_Y,res)
!       now for ea. lambda psf_Y is normalized to Area
        psfArea(iLambda) = res
        write(12,'(1p,2e15.3)') lambdaOut(iLambda), psfArea(iLambda)

        j = iLambda
!       generate 1D intensity vector for subroutine Conv2D
!       take only diffuse emission, stellar contribution will
!       be added below (a shortcut to avoid inaccuracies or too many
!       points in Conv2D)
        DO i = 1, nP+2
          IF (i.LE.2) THEN
            Int1D(i) = IntOut(j,3)
          ELSE
            Int1D(i) = IntOut(j,i)
          END IF
          CALL ChkRange(dynrange,Int1D(i))
          IF (Int1D(i).LT.dynrange) Int1D(i)=0.0D+00
        END DO
!       convolve
        CALL Conv2D(npP+2,nP+2,yang,Int1D,1000,NConv,Offset,Conv)

        DO i = 1, Nconv
          CALL ChkRange(dynrange,Conv(i))
          ConvInt(j,i) = Conv(i)
        END DO

!       add stellar contribution
        DO i = 1, nP+2
          ConvS=2.0D+00*ASIN(1.0D+00)*(yang(2)**2.0D+00)*IntOut(j,1)*PSFN(Offset(i))
          Conv(i) = Conv(i) + ConvS
        END DO
!       scale to 1 at the center
!        CALL ScaleTo1(1000,Nconv,Conv)
!       copy 1D convolved intensity to ConvInt
!       ConvInt is normalized to area in Sub PrOut
        DO i = 1, Nconv
          CALL ChkRange(dynrange,Conv(i))
          ConvInt(j,i) = Conv(i)
        END DO
!      end do over iLambda - lambda for conv,images
      END DO
      write(12,*)' --------------------------------------------'
!      STOP
! -----------------------------------------------------------------------
      RETURN
      END SUBROUTINE Convolve
! ***********************************************************************

! ***********************************************************************
      SUBROUTINE Conv2D(NinMax,Nin,Xin,Yin,Noutmax,Nout,Xout,Yout_loc)
! =======================================================================
! This subroutine convolves intensity Yin(Xin[i]), i=1,Nin with
! the point spread function PSFN(x) (provided as a separate function).
! It is assumed that both the intensity yin and PSFN(x) are circularly
! symmetric functions of radial coordinate x, i.e., this subroutine
! performs two-dimensional convolution. Convolved intensity, Yout, is
! evaluated for very position Xout[i], i=1,Nout, as:
!        Yout(Xout) = Int[Yin(Yloc)*PSF(xloc)*Yloc*dYloc*dphi]
! where xloc = sqrt(Yloc **2+Xout**2-2*Yloc*Xout*cos(phi), with Yloc
! and phi being dummy integration variables. Declared size of Xin is
! NinMax, the one for Xout is NoutMax. The radial integration is done
! using subroutine ROMBY and angular integration is done by using
! Simpson rule.                                        [Z.I., Jan. 1997]
! =======================================================================
        use common
      IMPLICIT none

      INTEGER NinMax, Nin, NoutMax, Nout, iPhi, iXin, Nphi, iXOut
      DOUBLE PRECISION Xin(NinMax), Yin(NinMax), Xout(NoutMax), A, B,  &
           Yout_loc(NoutMax), dphi, phi_loc(1000), fphi(1000), int1, &
           int2, imagfn
      EXTERNAL imagfn
! -----------------------------------------------------------------------
!     Parameters for integration:
!     number of angular points
      Nphi = 9
!     step in angle phi
      dphi = 2.0D+00*ASIN(1.0D+00) / (Nphi-1)
!     flag for imgfn
      ftype = 1
!     Start integrations
!     loop over output positions
      DO iXout = 1, Nout
        Cxout = Xout(iXout)
!       loop over angular wedges (phi integration)
        DO iPhi = 1, Nphi
          phi_loc(iPhi) = dphi*1.0D+00*(iPhi-1)
          Cphi = phi_loc(iPhi)
          fphi(iPhi) = 0.0D+00
!         loop over input radial positions (radial integration)
          DO iXin = 1, Nin-1
            Ckn = 1.0D+00
            CALL ROMBY(imagfn,Xin(iXin),Xin(iXin+1),int1)
            Ckn = 2.0D+00
            CALL ROMBY(imagfn,Xin(iXin),Xin(iXin+1),int2)
!           contribution from this annulus (lin. approx. for intensity)
            A = Xin(iXin+1)*Yin(iXin) - Xin(iXin)*Yin(iXin+1)
            A = A / (Xin(iXin+1)-Xin(iXin))
            B = (Yin(iXin+1)-Yin(iXin)) / (Xin(iXin+1)-Xin(iXin))
            fphi(iPhi) = fphi(iPhi) + A*int1 + B*int2
          END DO
        END DO
        CALL Simpson(1000,1,Nphi,phi,fphi,Yout_loc(iXout))
      END DO
! -----------------------------------------------------------------------
      RETURN
      END SUBROUTINE Conv2D
! ***********************************************************************


! *************************************************************************
      SUBROUTINE FindErr(flux,maxFerr)
! ========================================================================
! This subroutine finds maximum err in flux conservation for both
! spherical and slab case as (fmax-fmin)/(fmax+fmin)   [MN,Aug'99]
! =========================================================================
     use common
      IMPLICIT none

      INTEGER iY
      DOUBLE PRECISION flux(npY), maxFerr, fmin, fmax, aux, accFbol
! -----------------------------------------------------------------------
     accFbol = 1.0d-04

!     Find the min and max of fbol values
!     The abs and lower limit on fbol are protection for the case
!     of completely symmetri! slab illumination. The lower limit
!     is bound by the numerical accuracy of the flux calculation
        fmin = 1.e5
        fmax = 0.
        DO iY = 1, nY
           aux = flux(iY)
           IF (ksi.eq.1.0) aux = dabs(aux)
           IF (dabs(aux).LE.accFbol) aux = accFbol
           IF(aux.LT.fmin) fmin = aux
           IF(aux.GT.fmax) fmax = aux
        END DO
        if (fmax.LT.0.) then
!     bad solution; overall flux cannot be negative
            maxFerr = 1
        else
            maxFerr = (fmax - fmin)/(fmax + dabs(fmin))
        end if
! -----------------------------------------------------------------------
      RETURN
      END
! *************************************************************************

!***********************************************************************
      SUBROUTINE GetbOut(npP,nP,P,pstar,bOut,k)
!=======================================================================
! This subroutine inserts two impact parameters corresponding to pstar,
! producing bOut(nP+2) from P(nP). The inserted elements are bOut(k) and
! bOut(k+1)                                            [Z.I., Aug. 1996]
! =======================================================================
      IMPLICIT none
      INTEGER npP, nP, k, kstop, i
      DOUBLE PRECISION P(npP), bOut(npP+2), pstar
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

! ***********************************************************************
      DOUBLE PRECISION FUNCTION IMAGFN(Yy)
! =======================================================================
! This function evaluates auxiliary functions needed to produce
! visibility curves and convolved images. It is called from the image
! integration subroutine ROMBY.                        [Z.I., Jan. 1997]
! =======================================================================
     use common
      IMPLICIT none

      DOUBLE PRECISION x, Yy, PSFN, Bessel

! -----------------------------------------------------------------------
      IF (ftype.EQ.1) THEN
!       this part is for convolution
        x = dsqrt(dabs(Cxout*Cxout+Yy*Yy-2.0D+00*Cxout*Yy*dcos(Cphi)))
        imagfn = PSFN(x) * Yy**Ckn
      ELSE
!       this part is for visibility
!       argument is Pi*q*Yy (not 2*Pi*q*Yy) to account for the fact that
!       theta1 is diameter rather than radius (so V is function of
!       q*theta1, like in IE, '96, MNRAS 279, 1019)
        imagfn = Bessel(2.0D+00*dSIN(1.0D+00)*Cqtheta*Yy) * Yy**Ckn
      END IF
! -----------------------------------------------------------------------
      RETURN
      END
! ***********************************************************************

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
  write(18,*) 'Tsub(1)=', Tsub(1)
! write(18,*) '  lambda(iL)  SigmaA(1,iL) Planck(xP)   xP'
  do iL = 1, nL
   qaux(iL) = sigmaA(1,iL)*us(iL,1)/lambda(iL)
   xP = 14400.0d0/Tsub(1)/lambda(iL)
   qaux2(iL) = sigmaA(1,iL)*Planck(xP)/lambda (iL)
  end do

  call Simpson(npL,1,nL,lambda,qaux,res1)
  call Simpson(npL,1,nL,lambda,qaux2,res2)
! approximate psi for opt.thin case:
  psitn = res1/res2
  mx = Tsub(1)*sqrt(sqrt(4.0d0/psitn))
  if(Tsub(1).lt.mx) then
   Te_min_loc = 2.0d0*mx
  else
   Te_min_loc = 2.0d0*Tsub(1)
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


! *******************************************************************
  subroutine philam(alam,f,al,phi_loc)
! =====================================================================
! interpolates IRAS filters [f(4)] for a given wavelength alam

  implicit none
  integer i,im,istop
  double precision a,b,alam,al(4,7),phi_loc(4,7),f(4)
! -----------------------------------------------

  do i = 1, 4
   f(i) = 0.0d0
   im = 0
   istop = 0
   do while (istop.ne.1)
    im = im + 1
    if ( (alam-al(i,im))*(alam-al(i,im+1)).le.0.0d0) then
     a = (phi_loc(i,im+1)-phi_loc(i,im))/dlog10(al(i,im+1)/al(i,im))
     b = phi_loc(i,im) - a*log10(al(i,im))
     f(i) = a*dlog10(alam) + b
     if ( f(i).gt.1.0d0 ) f(i) = 1.0d0
    end if
    if (im.eq.6) istop = 1
   end do
  end do
! ------------------------------------------------------
  return
  end subroutine philam
! ***********************************************************************

! ***********************************************************************
  subroutine SpFeatur(model,spectr,charac)
! =======================================================================
! This subroutine calculates IRAS colors and other spectral quantities
! Filters data from Neugebauer et al, 1984, ApJ, 278, L1.
! Procedure described in Bedijn, 1987, A&A, 186, 136.  [Z.I., Mar. 1996]
! =======================================================================

  use common
  implicit none

  integer i, j, iaux, model
  double precision  ff1(7), ff2(7), f3(7), f4(7), w1_loc(7), w2(7), &
       w3(7),w4(7), phi_loc(4,7), al(4,7), cl(4), prz(4),        &
       tinf(4), flxy(9),wwav(9), spectr(npL), charac(11), wmid,  &
       flx1,flx2, faux, B98, B11, rat9818, beta813, beta1422,    &
       f12,f25, f60, f100, an98, f98c, f11c, an11
! -----------------------------------------------------------------------

! data for 4 IRAS filters
! wavelengths
  data w1_loc/7.55d+0,8.0d+0,10.3d+0,11.5d+0,13.4d+0,14.7d+0,15.5d+0/
  data w2/16.6d+0,22.5d+0,25.6d+0,26.8d+0,27.5d+0,29.3d+0,31.1d+0/
  data w3/30.5d+0,40.1d+0,40.2d+0,65.2d+0,74.2d+0,83.8d+0,83.9d+0/
  data w4/72.7d0,95.4d0,111.2d0,116.6d0,137.4d+0,137.5d0,137.6d0/
! transmittivities
  data ff1/0.000d0,0.618d0,0.940d0,0.750d0,1.022d0,0.906d0,0.000d0/
  data ff2/0.235d0,0.939d0,0.939d0,0.745d0,0.847d0,0.847d0,0.000d0/
  data f3/0.000d0,0.102d0,0.260d0,1.026d0,0.842d0,0.001d0,0.000d0/
  data f4/0.000d0,0.910d0,1.000d0,0.330d0,0.002d0,0.001d0,0.000d0/
! ------------------------------------------------------------
! initialization
  do i = 1, 4
   tinf(i) = 0.0d0
   cl(i) = 0.0d0
  end do
  do j = 1, 7
   al(1,j) = w1_loc(j)
   phi_loc(1,j) = ff1(j)
   al(2,j) = w2(j)
   phi_loc(2,j) = ff2(j)
   al(3,j) = w3(j)
   phi_loc(3,j) = f3(j)
   al(4,j) = w4(j)
   phi_loc(4,j) = f4(j)
  end do
! ------------------------------------------------------------------------
! first find IRAS colors
  do j = 2, nL
! middle wavelength
   wmid = 0.5d0*(lambda(j-1)+lambda(j))
! interpolate filters for wmid
   call philam(wmid,prz,al,phi_loc)
! convert spectrum to flambda
   flx1 = spectr(j-1) / lambda(j-1)
   flx2 = spectr(j) / lambda(j)
! add contribution to the integral (index is over filters)
   do i = 1, 4
    tinf(i) = tinf(i) + prz(i)*0.5d0*(flx2+flx1)*(lambda(j)-lambda(j-1))
    cl(i) = cl(i) + prz(i) * (lambda(j)-lambda(j-1))
   end do
  end do
  do i = 1, 4
   tinf(i) = tinf(i) / cl(i)
  end do
! spectrum corrected for IRAS filters
  f12 = tinf(1)*12.0d0
  f25 = tinf(2)*25.0d0
  f60 = tinf(3)*60.0d0
  f100 = tinf(4)*100.0d0
! now find  B98, B11, f98/f18, beta 8-13, beta 14-22
! find fluxes at all needed wavelengths (energy increases with index)
  data wwav/2.2, 8.0, 9.8, 11.3, 13.0, 14.0, 18.0, 22.0, 0.55/
  do j = 1, 9
   call lininter(npL,nL,lambda,spectr,wwav(j),iaux,faux)
   flxy(j) = faux
  end do
! the feature strength at 9.8 and 11.4 microns
  if((flxy(2)*flxy(5)*flxy(3)*flxy(4)).gt.0.0d0) then
   an98 = log(flxy(5)/flxy(2))/log(wwav(5)/wwav(2))
   f98c = flxy(2)*(wwav(3)/wwav(2))**an98
   B98 = log(flxy(3)/f98c)
   an11 = log(flxy(5)/flxy(3))/log(wwav(5)/wwav(3))
   f11c = flxy(3)*(wwav(4)/wwav(3))**an11
   B11 = log(flxy(4)/f11c)
  else
   B98 = 0.0d0
   B11 = 0.0d0
   if(iX.ge.1) write(18,*) ' no 10 micron feature strength, fluxes are 0.'
  end if
! ratio f9.8/f18
  if((flxy(2)*flxy(5)*flxy(7)).gt.0.0d0) then
   rat9818 = flxy(3)/flxy(7)*wwav(7)/wwav(3)
! beta 8-13 and beta 14-22 (see neugebauer)
   beta813 = dlog(flxy(5)/flxy(2))/dlog(13.0d0/8.0d0) - 1.0d0
  else
   rat9818 = 0.0d0
   beta813 = 0.0d0
   if(iX.ge.1) write(18,*) ' no 8-13 microns slope, fluxes are 0.'
  end if
  if((flxy(6)*flxy(8)).gt.0.0d0) then
   beta1422 = dlog(flxy(8)/flxy(6))/dlog(22.0d0/14.0d0)-1.0d0
  else
   beta1422 = 0.0d0
   if(iX.ge.1) write(18,*) ' no 14-22 microns slope, fluxes are 0.'
  end if
! store specchar to output array specchar
  charac(1) = B98
  charac(2) = B11
  charac(3) = rat9818
  charac(4) = beta813
  charac(5) = beta1422
  charac(6) = flxy(9)
  charac(7) = flxy(1)
  charac(8) = f12
  if((f12*f25).gt.0.0d0) then
   charac(9) = dlog10(25.0d0*f25/f12/12.0d0)
  else
   charac(9) = 0.0d0
  end if
  if((f12*f60).gt.0.0d0) then
   charac(10) = dlog10(60.0d0*f60/f12/12.0d0)
  else
   charac(10) = 0.0d0
  end if
  if((f60*f100).gt.0.0d0) then
   charac(11) = dlog10(f100*100.0d0/f60/60.0d0)
  else
   charac(11) = 0.0d0
  end if
! ------------------------------------------------------------------
  return
  end subroutine SpFeatur
! *******************************************************************

! ***********************************************************************
  subroutine Spectral(model)
! =======================================================================
! This subroutine finds the spectral features for spp and zpp files.
! It employs Sub SpFeatur                              [MN, Jan'99]
! =======================================================================

  use common
  implicit none

  integer i, iL, model, nchar
  parameter (nchar=11)
  double precision charac(nchar),spectr(npL)
! -----------------------------------------------------------------------

! find features for *.spp file
  if (slb) then
   do iL = 1, nL
    spectr(iL) = ftot(iL,nY) + ksi*fsR(iL,nY)
   end do
  else
   do iL = 1, nL
!  only the shell emission, external contribution not added
    spectr(iL) = fsL(iL,nY) + fde(iL,nY) + fds(iL,nY)
   end do
  end if
  call spfeatur(model,spectr,charac)
  do i = 1, nchar
   SpecChar(i,model) = charac(i)
  end do
! find features for *.zpp file for slab
  if (slb) then
   do iL = 1, nL
    spectr(iL) = dabs(ftot(iL,1) - fsL(iL,1))
   end do
   call spfeatur(model,spectr,charac)
   do i = 1, nchar
    SpecChar(i+11,model) = charac(i)
   end do
  end if
! -----------------------------------------------------------------------
  return
  end subroutine Spectral
! ***********************************************************************

! ***********************************************************************
  subroutine SLBintensity(omega,em)
! =======================================================================
  use common
  implicit none

  integer  iL, iY, imu
  double precision em(npL,npY),omega(npG,npL),tau1(npY),idifp, idifm,res, Sexp, Kron
  external Sexp
! -----------------------------------------------------------------------

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
    do iY = 2, nY
     tau1(iY) = TAUslb(iL,iY)
     Sfn=(1.0d0-omega(1,iL))*em(iL,iY)+omega(1,iL)*utot(iL,iY)
! transmit=1 for tau < t, transmit=0 for tau > t
     transmit = 1
     call romby(Sexp,tau1(iY-1),tau1(iY),res)
     idifm = idifm + res
     transmit = 0
     call romby(Sexp,tau1(iY-1),tau1(iY),res)
     idifp = idifp + res
    end do
    if(idifm.lt.1.d-20) idifm = 0.0d0
    if(idifp.lt.1.d-20) idifp = 0.0d0
    SLBintm(imu,iL) = idifm
    SLBintp(imu,iL) = idifp
! enddo over angles
   end do
! enddo over lambda
  end do
! -----------------------------------------------------------------

  return
  end subroutine SLBintensity
! ***********************************************************************


!***********************************************************************
SUBROUTINE SPH_Int(nG,omega,fs)
!***********************************************************************
! This is the former  SUBROUTINE FindInt(nG,ETAzp) [MN].
! This subroutine finds the intensity distribution at outer edge and for
! user specified wavelengths lamOut. It also evaluates the angular size
! of the stellar disk and adds two impact parameters describing the star
! to the P grid, thus producing bOut grid. All intensities are indeed
! dimensionless quantities lambda*I_lambda/F1 where I_lambda is real
! physical quantity defined as usual and F1 is the bolometric flux at
! the dust sublimation radius, r1. For conversion to the physical value
! lambda*I_lambda, I_lambda from the program has to be multiplied by F1.
! F1 can obtained either as:
!      1) F1 = 4*sigma*Tsub**4/Psi (IE96, eq. 15),
! where Tsub is sublimation temperature and parameter Psi is given in
! *.spp file; or as:
!      2) F1 = Fbol/alpha1**2 (IE96, eq. 34)
! where Fbol is the bolometric flux and alpha1 is the angular size of r1
! at any particular distance from the envelope (i.e. both F1 and alpha1
! correspond to observed quantities). Also note that
!     INT(I_lambda(p)*2Pi*P*dP) = f_lambda
! where I_lambda is the scaled quantity from the program, P is impact
! parameter, and f_lambda is the spectral shape F_lambda/Fbol.
!                                                      [Z.I., Aug. 1996]
! =======================================================================
  use common
  implicit none

      INTEGER iL,nG,iY,k,i,iLout,iLstop, iP, iW, nZ, Nzpt, iZ, izloc
      DOUBLE PRECISION qaux(npL), Psi, alpha(1,npY), resaux, QUtot1,   &
            QpTsub, xP, pst, stelfact, Istell(npL), Planck, z1, z2,    &
            IntL, IntR, xx , alb, Ids(npL,npP), Ide(npL,npP),w1, w2,lw12, &
            numcorr, ETAzpStar, qaux2(npL), fs(npL,npY), omega(1,npL), &
            delz, zloc, wloc, resint, pT, Tz, Idboth, tzp(100), pUtot, &
            Semis(100), Sscat(100), IntETA, palb, palf, alfa, exterm,  &
            Utotloc, Sstem(100), Sstsc(100), Istem(npL,100), Idfront,  &
            Istsc(npL,100), delTau, factaux, UtotL, UtotR, ep1,        &
            tauzp1, tauInf, fnum(npL), fdenum(npL), res, denum
EXTERNAL IntETA

! -----------------------------------------------------------------------
!     temporary
      IF (nG.GT.1.AND.iX.GE.1) THEN
        write(18,*)' FindInt should be fixed, nG>1 !'
        stop
      END IF
!     find impact parameter tangential to the stellar disk
!     first find the Planck averaged absorption efficiencies at Y=1
      DO iL = 1, nL
        qaux(iL) = SigmaA(1,iL) * Utot(iL,1) / lambda (iL)
        xP = 14400.0 / Td(1,1) / lambda(iL)
        qaux2(iL) = SigmaA(1,iL) * Planck(xP) / lambda (iL)
      END DO
      CALL Simpson(npL,1,nL,lambda,qaux,resaux)
      QUtot1 = resaux
      CALL Simpson(npL,1,nL,lambda,qaux2,resaux)
      QpTsub = resaux
!     parameter Psi (see Ivezic & Elitzur, 1996, eq. C4)
      Psi = QUtot1 / QpTsub
      alpha(1,1) = Psi
! ***alpha is a local array, to match the changes with Zeljko's old expressions
      DO iY = 2, nY
!       calculate f1 and f2
        DO iL = 1, nL
          fnum(iL) = SigmaA(1,iL) * Utot(iL,iY) / lambda(iL)
        END DO
        CALL Simpson(npL,1,nL,lambda,fnum,res)
!        calculate alpha
        DO iL = 1, nL
           xP = 14400.0 / lambda(iL) / Td(1,iY)
           fdenum(iL) = SigmaA(1,iL) * Planck(xP) / lambda(iL)
        END DO
        CALL Simpson(npL,1,nL,lambda,fdenum,denum)
        alpha(1,iY) = res /  denum
      END DO
! *********
!     ratio pst = rstar/rsub (see Ivezic & Elitzur, 1996, eq. 27)
!     Added Apr.07 [MN] to take care of the case of no central source
      IF (Left.eq.0) THEN
         pst = 1.
      ELSE
         pst = 2.0 / dsqrt(Psi) * (Td(1,1) / Tstar(1))**2.0
      END IF
!      IF (pst.GE.0.5) THEN
!     this is if only central source is present
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
!     generate bOut, i.e. insert two points such that
!     bOut(k)=0.999*pst and bOut(k+1)=1.001*pst
      CALL GetbOut(npP,nP,P,pst,bOut,k)
!     correction for numerical errors in tau
      numcorr = 1. / TAUtot(1)
!     loop over wavelengths
      DO iL = 1, nL
!        stellar intensity, Istell (extinction already included)
         Istell(iL) = fs(iL,nY) * stelfact
!        total optical depth along a line of sight
         tauOut(iL) = numcorr*TAUtot(iL)
      END DO

!     generate diffuse intensities, Ide (emission) and Ids (scat)
!     loop over wavelengths
      DO iL = 1, nL
        DO iP = 1, nP
!         maximal number of points along tangential position, z
          nZ = nY + 1 - iYfirst(iP)
!         starting value for local radius
          IF (P(iP).GE.1.0) THEN
            w2 = P(iP)
          ELSE
            w2 = 1.0
          END IF
!         initialize intensities
          Ide(iL,iP) = 0.0
          Ids(iL,iP) = 0.0
          IF (iP.LE.k+1) THEN
            Istem(iL,iP) = 0.0
            Istsc(iL,iP) = 0.0
          END IF
!         total optical depth along this impact parameter
           ep1 = ETAzp(iP,nZ)*TAUtot(iL)
!         loop over z, i.e. steps over points crossing the y grid
          DO iZ = 2, nZ
!          index for the ending local radius
            iW = iYfirst(iP) + iZ - 1
!          local boundary radii
            w1 = w2
            w2 = Y(iW)
!           corresponding displacements along a line of sight
            z1 = sqrt(abs(w1**2.-P(iP)**2.))
            z2 = sqrt(abs(w2**2.-P(iP)**2.))
!           # of pts. for z integration, should increase with deltaTau
!           it is messy because INT function which would do the job is
!           not in F77 standard set
            Nzpt = 5
            delTau = (ETAzp(iP,iW)-ETAzp(iP,iW-1))*TAUtot(iL)
            IF (delTau.GT.1) Nzpt = 10
            IF (delTau.GT.5) Nzpt = 20
            IF (delTau.GT.10) Nzpt = 30
            IF (delTau.GT.20) Nzpt = 40
            IF (delTau.GT.50) Nzpt = 50
            delz = (z2-z1) / (Nzpt-1)
!           powers for power-law interpolations between 2 y pts.
            lw12 = dlog(Y(iW-1)/Y(iW))
!           for T
            pT = dlog(Td(1,iW)/Td(1,iW-1)) / lw12
!           for albedo
            IF (omega(iL,iW-1).GT.0.0.AND.omega(iL,iW).GT.0.0) THEN
              palb = dlog(omega(iL,iW)/omega(iL,iW-1)) / lw12
            ELSE
              palb = 0.0
            END IF
!           for Utot
            UtotL = Utot(iL,iW-1)
            UtotR = Utot(iL,iW)
            CALL ChkRange(dynrange,UtotL)
            CALL ChkRange(dynrange,UtotR)
            IF (UtotL.GT.0.0.AND.UtotR.GT.0) THEN
              pUtot = dlog(UtotR/UtotL) / lw12
            ELSE
              pUtot = 0.0
            END IF
!           for alpha
            palf = dlog(alpha(1,iW)/alpha(1,iW-1)) / lw12
!           tauzp between z=0 and z=z1
            tauzp1 = ETAzp(iP,iZ-1)*TAUtot(iL)
!           integrate between adjacent grid points
            DO izloc = 1, Nzpt
              zloc = z1 + (izloc-1)*delz
              wloc = sqrt(zloc**2 + P(iP)**2)
!             find local TAUzp(w(z))-TAUzp(w1=w(z1))
              tzp(izloc) = IntETA(P(iP),iW-1,w1,wloc)*TAUtot(iL)
!             find Tz = T(zloc) = T(wloc), this works for single
!             size grains only; for multigrain case one needs to
!             get Semis by summation over all Td
              Tz = Td(1,iW-1) * (Y(iW-1)/wloc)**pT
              xP = 14400/lambda(iL)/Tz
!             power-law interpolation for albedo
              alb = omega(iL,iW-1) * (Y(iW-1)/wloc)**palb
!             power-law interpolation for Utot
              IF (UtotL.GT.0) THEN
                UtotLoc = UtotL * (Y(iW-1)/wloc)**pUtot
              ELSE
                UtotLoc = 0.0
              END IF
              CALL ChkRange(dynrange,UtotLoc)
!             power-law interpolation for alpha
              alfa = alpha(1,iW-1) * (Y(iW-1)/wloc)**palf
!             source functions (wloc**2 because D uses scaled quant.)
              factaux = 1 / wloc**2 / (4. * Pi)
              Semis(izloc) = (1-alb) * alfa * Planck(xP) * factaux
              Sscat(izloc) = alb * UtotLoc * factaux
!             check for the dynamic range
              CALL ChkRange(dynrange,Semis(izloc))
              CALL ChkRange(dynrange,Sscat(izloc))
!             optical depth from infinity along the line of sight
              tauInf = ep1 - tauzp1 - tzp(izloc)
!             for a line of sight terminating on the star find
!             contribution only from the front part of the envelope
              IF (iP.LE.k+1) THEN
                 IF (tauInf.LT.50) THEN
                   exterm = dexp(-tauInf)
                 ELSE
                   exterm = 0.0
                 END IF
                 Sstem(izloc) = Semis(izloc) * exterm
                 Sstsc(izloc) = Sscat(izloc) * exterm
              END IF
!             otherwise take both the front and back contributions
              IF (tauInf.LT.50) THEN
                 exterm = dexp(-tauInf)+dexp(-tauInf-ep1)
              ELSE
                 exterm = 0.0
              END IF
              Semis(izloc) = Semis(izloc) * exterm
              Sscat(izloc) = Sscat(izloc) * exterm
!             end of local loop over z
            END DO
!           integrate and add contribution from this step
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
!         end of loop over z
          END DO
!       end of loop over impact parameter, iP
        END DO
!     end of loop over wavelengths, iL
      END DO
!     add all intensities, Istell, Ide, Ids
      DO iL = 1, nL
!       interpolate optical depth  at pstar
        IF (iL.EQ.iLfid) THEN
          ETAzpStar = (ETAzp(k,nY) - ETAzp(k-1,nY))
          ETAzpStar = ETAzpStar * (pst-P(k-1)) / (P(k)-P(k-1))
          ETAzpStar = ETAzp(k-1,nY) + ETAzpStar
        END IF
!       find diffuse contribution at pstar (by linear interpolation)
        Idfront = Istsc(iL,k)+Istem(iL,k)-Istsc(iL,k-1)-Istem(iL,k-1)
        Idfront = Idfront * (pst-P(k-1)) / (P(k) - P(k-1))
        Idfront = Idfront + Istsc(iL,k-1) + Istem(iL,k-1)
        Idboth = Ids(iL,k) + Ide(iL,k) - Ids(iL,k-1) - Ide(iL,k-1)
        Idboth = Idboth * (pst-P(k-1)) / (P(k) - P(k-1))
        Idboth = Idboth + Ids(iL,k-1) + Ide(iL,k-1)
!       first for p<pstar, all three contributions
        DO i = 1, k-1
          Intens(iL,i) = Istell(iL) + Istsc(iL,i) + Istem(iL,i)
          IF (iL.EQ.iLfid) tauZout(i) = ETAzp(i,nY)/ETAzp(1,nY)
        END DO
!       barely on the stellar disk
        Intens(iL,k) = Istell(iL) + Idfront
        tauZout(k) = ETAzpStar/ETAzp(1,nY)
!       barely off the stellar disk
        Intens(iL,k+1) = Idboth
        tauZout(k+1) = 2. * tauZout(k)
!       all other p>pstar
        DO i = k, nP
          Intens(iL,i+2) = Ids(iL,i)+Ide(iL,i)
          IF (iL.EQ.iLfid) THEN
            nZ = nY + 1 - iYfirst(i)
            tauZout(i+2) = 2. * ETAzp(i,nZ)/ETAzp(1,nY)
          END IF
        END DO
      END DO
!     check dynamic range
      DO iL = 1, nL
        DO i = 1, nP+2
          CALL ChkRange(dynrange,Intens(iL,i))
        END DO
      END DO
!     now interpolate Intens(lambda) to lamOut
      DO iLout = 1, NlambdaOut
!       bracket the needed wavelength
        iLstop = 0
        iL = 0
        DO WHILE (iLstop.EQ.0)
          iL = iL + 1
          IF (lambda(iL).GT.LambdaOut(iLout)) iLstop = 1
          IF (iL.EQ.nL) iLstop = 1
        END DO
!       interpolate intensity
        xx = (LambdaOut(iLout)-lambda(iL-1))/(lambda(iL)-lambda(iL-1))
        DO i = 1, nP+2
          IntL = Intens(iL-1,i)
          IntR = Intens(iL,i)
          IntOut(iLout,i) = IntL + xx*(IntR - IntL)
          CALL ChkRange(dynrange,IntOut(iLout,i))
        END DO
      END DO
! -----------------------------------------------------------------------
999   RETURN
end subroutine SPH_Int
! ***********************************************************************

! ***********************************************************************
      SUBROUTINE Visibili
! =======================================================================
! This subroutine finds visibility functions corresponding to IntOut.
! The work horse is subroutine Visi2D, and this subroutine is used to
! prepare everything.                                  [Z.I., Jan. 1997]
! =======================================================================
        use common
      IMPLICIT none

      INTEGER i, j, N1, N2
      DOUBLE PRECISION  Visi(1000), Int1D(npP+2)
! -----------------------------------------------------------------------
!     generate spatial frequency (q) grid
!     first N1 points up to qtheta1=1.22 (Rayleigh limit for a disk)
      N1 = 80
!     first 2 points manually:
!     there must be 0!
      qtheta1(1) = 0.0D+00
!     make sure the whole envelope is resolved
      qtheta1(2) = 0.5D+00 / bOut(nP+2)
!     and the rest on logarithmic grid up to 1.22
      DO i = 1, N1-2
       qtheta1(i+2)=qtheta1(2)*(1.22D+00/qtheta1(2))**(i*1.0D+00/(N1-2))
      END DO
!     envelope is well sampled, now to be sure that the star will be OK
!     for small taus add N2 points on a logarithmic grid up to 1.22/p*
      N2 = 20
      DO i = 1, N2
        qtheta1(N1+i) = 1.22D+00 / bOut(2)**(i*1.0D+00/N2)
      END DO
      Nvisi = N1 + N2
!     find visibility wavelength by wavelength
      DO j = 1, NlambdaOut
        DO i = 1, nP+2
          Int1D(i) = IntOut(j,i)
          CALL ChkRange(dynrange,Int1D(i))
          IF (Int1D(i).LT.dynrange) Int1D(i)=0.0D+00
        END DO
        CALL Visi2D(npP+2,nP+2,bOut,Int1D,1000,N1+N2,qtheta1,Visi)
!       copy 1D convolved visibility to Visib
        DO i = 1, N1+N2
!         check dynamic range
          CALL ChkRange(dynrange,Visi(i))
          Visib(j,i) = Visi(i)
        END DO
      END DO
! -----------------------------------------------------------------------
      RETURN
      END
! ***********************************************************************

! ***********************************************************************
  SUBROUTINE Visi2D(NinMax,Nin,Xin,Yin,Noutmax,Nout,Xout,Yout_loc)
! =======================================================================
! This subroutine finds the visibility function (the spatial Fourier
! transform of the intensity distribution) corresponding to the
! intensity Yin(Xin[i]), i=1,Nin. Visibility, Yout, is evaluated at q
! positions (spatial frequency) given in Xout[i], i=1,Nout. Maximum size
! of Xin is NinMax, maximum size of Xout is NoutMax. The Bessel function
! of the zeroth order is provided separately. The integration is done by
! calling subroutine ROMBY (Bessel function is called from IMGFN).
! Note:
! The visibility function V(q) for a circularly symmetric intensity
! I(x) is:
!          V(q) = F(q)/F(0)
! where Jo is the Bessel function of the zeroth order, and
!          F(q) = Int[Jo(2Pi*q*x)*I(x)*2Pi*x*dx]
! Note that F(0) is nothing more than flux. For more details see
! Ivezic & Elitzur, 1996, MNRAS, 279, 1019 and ref. therein.
!                                                      [Z.I., Jan. 1997]
! =======================================================================
     use common
      IMPLICIT none

      INTEGER NinMax, Nin, NoutMax, Nout, iq, iXin
      DOUBLE PRECISION Xin(NinMax), Yin(NinMax), Xout(NoutMax),  &
           Yout_loc(NoutMax),  F(1000), F0, int1, int2, A, B, imagfn
      EXTERNAL imagfn
! -----------------------------------------------------------------------
!     loop over spatial frequency q (= Xout)
      DO iq = 1, Nout
        Cqtheta = Xout(iq)
        F(iq) = 0.0D+00
!       loop over radial positions
!!        DO iXin = 1, Nin - corrected after Zeljko's email from 7/29/09
        DO iXin = 1, Nin-1
!         find F(q)
          ftype = 2
          Ckn = 1.0D+00
          CALL ROMBY(imagfn,Xin(iXin),Xin(iXin+1),int1)
          Ckn = 2.0D+00
          CALL ROMBY(imagfn,Xin(iXin),Xin(iXin+1),int2)
!         contribution from this annulus (lin. approx. for intensity)
          A = Xin(iXin+1)*Yin(iXin)-Xin(iXin)*Yin(iXin+1)
          A = A /(Xin(iXin+1)-Xin(iXin))
          B = (Yin(iXin+1)-Yin(iXin))/(Xin(iXin+1)-Xin(iXin))
          F(iq) = F(iq) + A*int1 + B*int2
        END DO
      END DO
!     flux
      F0 = F(1)
      DO iq = 1, Nout
       IF(F0.EQ.0.0D+00) THEN
         Yout_loc(iq) = 0.0D+00
       ELSE
         Yout_loc(iq) = dabs(F(iq) / F0)
       END IF
      END DO
! -----------------------------------------------------------------------
      RETURN
      END SUBROUTINE Visi2D
! ***********************************************************************






! ***********************************************************************
      SUBROUTINE ChkConv(accuracy_loc,Aold,Anew,Aconv_loc)
! =======================================================================
! This subroutine checks convergence of an array A(nY) between values
! given in Aold and Anew. If the relative difference for EVERY element
! is smaller than accuracy, Aconv is assigned 1, otherwise 0.
!                                                      [Z.I., Jul. 1996]
! =======================================================================
        use common
      IMPLICIT none

      INTEGER iY, Aconv_loc
      DOUBLE PRECISION accuracy_loc, Aold(npY), Anew(npY), delta
! -----------------------------------------------------------------------
      Aconv_loc = 1
!     loop over radial positions
      DO iY = 1, nY
!       find relative difference
        delta = dabs(Anew(iY)-Aold(iY))
        IF (delta.GT.dabs(Anew(iY))*accuracy_loc) Aconv_loc = 0
      END DO
! -----------------------------------------------------------------------
      RETURN
      END
! ***********************************************************************


!***********************************************************************
 double precision function ETA(Yy)
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

   integer i, istop, iYdummy
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
  call lininter(npY,nYprev,Yprev,etadiscr,Yy,iYdummy,eta)
! otherwise use prescribed formulae for different cases
  else
! smooth power-law
  if (powd .and. Ntr.eq.0) then
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
  if (powd.and.Ntr.gt.0) then
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
  if (expd) then
   if (Yy.ge.1.0d0-1.0d-8) then
    eta = (Yout-1.) * (1.0d0-dexp(-pow)) / pow
    eta = dexp(-pow*(Yy-1.0d0)/(Yout-1.0d0)) / eta
   else
    eta = 0.0d0
   end if
  end if
! radiatively driven winds (the gray-body approximation)
  if (rdwa) then
   eps_loc = pow
   if (Yy.ge.1.0d0-1.0d-8) then
    eta = (1.0d0+eps_loc)/2.0d0/Yy/Yy/sqrt(1.0d0-(1.0d0-eps_loc*eps_loc)/Yy)
   else
    eta = 0.0d0
   endif
  end if
! radiatively driven winds (i.e. denstyp=3 or 6)
  if (rdw.or.rdwpr) then
! if this is the first iteration use analytic approximation
   if (itereta.lt.2) then
! for 3 eps_loc is pow, for 6 assume eps_loc=0.1
    if (rdw) then
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
    call lininter(npY,nYprev,Yprev,etadiscr,Yy,iYdummy,eta)
   end if
  end if
! user specified function for eta
  if (fild) then
   if (Yy.lt.yetaf(nYetaf)) then
    call lininter(npY,nYetaf,yetaf,etaf,Yy,iYdummy,eta)
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




!!=====================================================================
! Below are INPUT/OUTPUT related subroutines and functions
! arranged in alphabetical order.                       [MN, Aug,2010]
!!=====================================================================

!***********************************************************************
subroutine attach(root,length,ext,fname)
!=======================================================================
! Attaches extensions to the root cleaned by Clean
!=======================================================================
  implicit none
  integer i, length
  character*(*) root, ext, fname
! -----------------------------------------------------------------------

  do i = 1, len(fname)
   fname(i:i) = ' '
  end do
  fname(:length) = root(:length)
  fname(length + 1:) = ext
!-----------------------------------------------------------------------
  return
end subroutine attach
!***********************************************************************

!***********************************************************************
subroutine clean(strin, strout, length)
!=======================================================================
! Find meaningful part of strin without leading and trailing junk
! It is returned left-justified in StrOut, right-padded with blanks
! The number of meaningful characters is returned in Length.
! In case of any problems, StrOut is empty. This sub should be used to
! clean every input filename immediately after DUSTY reads it. [ME,'99]
!=======================================================================
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
!-----------------------------------------------------------------------
  return
end subroutine clean
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
  if(lambda(1).gt.0.01d0) then
   write(*,*)' *************** WARNING! ********************** '
   write(*,*)'  the shortest wavelength in lambda_grid.dat has '
   write(*,*)'  to be 0.01 microns. correct this and try again!'
   write(*,*)' *********************************************** '
   goto 999
  end if
  if(lambda(npL).lt.36000.0d0) then
   write(*,*)' *************** WARNING! ******************* '
   write(*,*)'  the longest wavelength in lambda_grid.dat   '
   write(*,*)'  has to be 3.6e4 um. correct this and try again!'
   write(*,*)' ******************************************** '
   goto 999
  end if
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

!***********************************************************************
SUBROUTINE ChkRange(dr,x)
!=======================================================================
! This subroutine checks if x is within the allowed range defined by
! dr<<1:
!         dr**2 < x < 1/dr**2
! If it is not then x = 0.0                            [Z.I., Jan. 1997]
!=======================================================================
      IMPLICIT none
      DOUBLE PRECISION x, dr
!-----------------------------------------------------------------------
      IF ((x-dr*dr)*(x-1.0d0/dr/dr).LT.0.0d0) THEN
        continue
      ELSE
!        continue
        x = 0.0d0
      END IF
!-----------------------------------------------------------------------
      RETURN
END SUBROUTINE ChkRange
!***********************************************************************

!***********************************************************************
subroutine CLLOSE(error,model,Nmodel)
! =======================================================================
! This subroutine closes output files.             [ZI,Feb'96; MN,Apr'99]
! =======================================================================
  use common
  implicit none

  character*235 su1, su2, s3, s4, txtt, txtf
  integer  error, model, Nmodel, im
! -----------------------------------------------------------------------

! close the default output file:
  if (iError.ne.0) then
   write(12,'(a42,i4)') ' There are some error messages for model:',model
   write(12,*) ' Please check m## file (if not produced then rerun)'
  end if
  if (iWarning.ne.0.and.iError.eq.0) then
   write(12,'(a36,i4)') ' There are some warnings for model:',model
   write(12,*)' Please check m## file (if not produced then rerun)'
  end if
  iCumm = iCumm + iError + iWarning
  if (model.eq.Nmodel.or.error.eq.3.or.error.eq.4) then
   if (error.ne.3) then
    su1=' ============================================================='
    su2='============================'
    if(rdw.or.rdwa.or.rdwpr) then
     write(12,'(a62,a28)')su1,su2
    else
     if(sph) then
      write(12,'(a62)')su1
     else
      su1=' ====================================================='
      if(ksi.gt.0) then
       su2='================='
       write(12,'(a53,a16)')su1,su2
      else
       su2='====='
       write(12,'(a53,a5)')su1,su2
      end if
     end if
    end if
    write(12,'(a23,1p,e8.1,a8)')'   (1) optical depth at',lamfid,' microns'
!  ----------  for slab output ----------------------------
    if(slb) then
! txtf = '  (2) bol.flux of the left-side source at the slab left boundary'
! txtt = '  (2) dust temperature at the slab left boundary'
     if(typentry(1).eq.5) then
      write(12,*) '  (2) bol.flux of the left-side source at the slab left boundary' ! txtf
     else
      write(12,*) '  (2) dust temperature at the slab left boundary' ! txtt
     end if
!     write(12,*)'  (3) f1=f/fi, where f is the overall bol.flux in the slab'

     if(typentry(1).eq.3) then
      write(12,*)'  (3) external bolometric flux at the slab left boundary'
     else
      write(12,*)'  (3) position of the left slab boundary for L = 1e4 Lo'
     end if
     write(12,*)'  (4) dust temperature at the right slab face'
     write(12,*)'  (5) effective temperature of the left source'
     if(ksi.gt.0.0d0) then
      write(12,*)'  (6) effective temperature of the right source'
!!**      write(12,*)'  (7) maximum error in flux conservation (%)'
      write(12,*)'  (7) maximum error in calculated diffuse flux (%)'
     else
!!**      write(12,*)'  (6) maximum error in flux conservation (%)'
      write(12,*)'  (6) maximum error in calculated diffuse flux (%)'
     end if
    else
!---------- for spherical shell ----------------------------
     if(typentry(1).eq.5) then
      write(12,*) '  (2) bolometric flux at the inner radius '
     else
      write(12,*) '  (2) dust temperature at the inner radius '
     end if
     if(typentry(1).eq.3) then
      write(12,*)'  (3) bolometric flux at the inner radius '
     else
      write(12,*)'  (3) inner radius for L = 1e4 Lo'
     end if
     write(12,*)'  (4) ratio of the inner to the stellar radius'
     write(12,*)'  (5) angular size (in arcsec) when Fbol=1e-6 W/m2'
     write(12,*)'  (6) dust temperature at the outer edge'
!!**     write(12,*)'  (7) maximum error in flux conservation (%)'
     write(12,*)'  (7) maximum error in calculated diffuse flux (%)'

     if(rdw.or.rdwa.or.rdwpr) then
      write(12,*)'  (8) mass-loss rate (in Mo/year)'
      write(12,*)'  (9) terminal outflow velocity (in km/s)'
      write(12,*)'  (10) upper limit of the stellar mass (in Msun)'
     end if
    end if
    write(12,*)'================================================='
    if(iCumm.eq.0) write(12,*)' Everything is ok for all models.'
    if(iSPP.ne.0) then
     if(slb.and.iSPP.eq.3) then
      write(12,*)' Tables with spectral properties are in files *.spp and *.zpp'
     else
      write(12,*)' Table with spectral properties is in file *.spp'
     end if
    end if
!---- private file ---
    if(rdwpr) then
     write(12,*)' Table with the wind properties is in file *.rdw'
    end if
!------------------------- spectra ------------------
    if(ia.ne.0) then
     if (ia.eq.1) then
      write(12,*)' All spectra are in file *.stb'
     else
      if (slb.and.ia.eq.3) then
       write(12,*)' Spectra are in files *.s## and *.z##'
      else
       write(12,*)' Spectra are in files *.s##'
      end if
     end if
    end if
!------------------------- images ------------------
    if (ic.ne.0) then
     if (abs(ic).eq.1) then
      if(slb) then
       write(12,*)' All intensities are in file *.itb'
      else
       write(12,*)' All imaging quantities are in file *.itb'
      end if
     else
      if(slb) then
       write(12,*)' Intensities are in files *.i##'
      else
       write(12,*)' Images are in files *.i##'
      end if
      if (iv.ne.0.and.abs(ic).eq.3) write(12,*)' Visibility curves are in files *.v##'
     end if
     if (ic.eq.-3.and.iPsf.ne.0) then
      write(12,*)' Convolved images are in files *.c##'
      if (psftype.lt.3) write(12,*)' Point spread functions are in file *.psf'
     end if
    end if
!------------------------- radial ------------------
    if(ib.ne.0) then
     if (ib.eq.1) then
      write(12,*)' All radial profiles are in file *.rtb'
     else
      write(12,*)' Radial profiles are in files *.r##'
     end if
    end if
!------------------------- messages ------------------
    if (iX.eq.1) write(12,*)' All run-time messages are in file *.mtb'
    if (iX.gt.1) write(12,*)' Run-time messages are in files *.m##'
   else
    write(12,*)' Ending calculations for this input file.'
   end if
  end if

  if (model.eq.Nmodel.or.error.eq.3.or.error.eq.4) then
   write(12,*) '========== the end =============================='
   close(12)
  end if

! table with spectral properties
  if (iSPP.ne.0) then
! in case of slab: add the zpp table after the spp (if desired)
   if(slb.and.model.eq.Nmodel.and.iSPP.ne.3) then
    s3='###   tau0      Psi      fV       fK       f12    C21  '
    s4=' C31   C43  b8-13 b14-22 B9.8 B11.4  R9.8-18  '
    write(19,'(a49)') '# ==============================================='
    write(19,'(a49)') '# properties of spectra from the slab left side  '
    write(19,'(a49)') '# -----------------------------------------------'
    write(19,'(a55,a46)')s3,s4
!!!    write(19,'(a100)') zline(model)
!!! the do loop over im is necessary, otherwise you don't append all previous zline's [MN]
    DO im = 1, Nmodel
       write(19,'(a100)') zline(im)
    END DO
!   (19) is the .zpp file; (24) is the .spp file
    if (model.eq.Nmodel) close(19)
    if (model.eq.Nmodel) close(24)
   end if
  end if

! close the psf file
  if (model.eq.1.or.error.eq.3) then
   if (iPsf.eq.1.and.psftype.lt.3) close(23)
  end if
! conditionally close the spectral files
  if(iA.eq.1) then
   if(model.eq.Nmodel) close(15)
  else
   close(15)
  end if
! conditionally close the radial files
  if(iB.eq.1) then
   if(model.eq.Nmodel) close(16)
  else
   close(16)
  end if
! conditionally close the imaging files
  if(abs(iC).eq.1) then
   if(model.eq.Nmodel) close(17)
  else
   close(17)
  end if
  if(iX.eq.1) then
   if(model.eq.Nmodel) close(18)
  else
   close(18)
  end if

  if (iPsf.ne.0) close(21)
  if (iV.ne.0) close(22)
!-----------------------------------------------------------------------

  return
end subroutine CLLOSE
!***********************************************************************


!***********************************************************************
double precision function EMfunc(lambda,teff,xSiO)
!***********************************************************************
! This is modeled after subroutine Engelke by M. Marengo. Here are his
! original comments:
!=================================================================
! This subroutine computes a modified black body spectrum using an
! "Engelke" function (see Engelke 1992, AJ 104, 1248):
! Bnu = Bnu(Tb) with Tb = 0.738*Teff*(1+79450/(lambda*Teff))**0.182
!
! Molecular SiO absorption is modelled from the alpha Tau spectrum
! of Cohen et al. 1992, AJ 104,2030 with a 5th order polinomial,
! and added to the modified bb.
!
! M. Marengo - mmarengo@cfa.harvard.edu - Sep 1998
!=================================================================
!
! This version makes use of the scaled quantities and Dusty's function
! Planck(x)                                                  [ZI, Feb 99]
!=======================================================================

  implicit none
  integer j
  double precision lambda, teff, xSiO, x, Planck, tenG, sioc(6), &
       lambda1, lambda2, sio8m, siof
!-----------------------------------------------------------------------
! SiO fit data from Massimo:
! Polinomial coeff for SiO absorption model (5th order),
! wavelength interval in which to apply the absorption
! and given absorption at 8 micron (to rescale for SiO)
  lambda1 =  7.8636d0
  lambda2 = 11.4280d0
  sio8m = 1.0701447d0
  sioc(1) = -300.43916d0
  sioc(2) =  149.32134d0
  sioc(3) =  -29.493280d0
  sioc(4) =    2.9067144d0
  sioc(5) =   -0.14304663d0
  sioc(6) =    0.0028134070d0
!-----------------------------------------------------------------------

! Engelke's effective temperature
  teng = 0.738d0*teff*(1.0d0 + 79450.0d0/(lambda*teff))**0.182d0
  x = 14400.0d0 / lambda / teng
  EMfunc = (teng/teff)**4.0d0 * Planck(x)
! If lambda is in SiO region, compute and apply the SiO absorption
  if ((lambda-lambda1)*(lambda-lambda2).lt.0.0d0) then
   siof = 0.0d0
   do j = 1, 6
    siof = siof + sioc(j) * lambda**(1.0d0*j-1)
   end do
   EMfunc = EMfunc / (1.0d0+ (siof-1)/(sio8m-1)*xSiO*0.01d0)
  end if
!-----------------------------------------------------------------------

  return
end function EMfunc
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
subroutine getSpShape(shp,is)
!=======================================================================
! Produces source sp. shape. This was part of subroutine star, separated
! here for clarity.  [MN]
!=======================================================================

  use common
  implicit none

  integer iY, iL, iLs, nLs, k, kstop, i, error, nLambdam, nis, is
!  nLambdam is the max number entries for a user supplied stellar spectrum
  parameter (nLambdam = 10000, nis = 2)
  double precision lambdas(nLambdam), Llamstar(nLambdam), shp(npL), &
       stellar(nLambdam), fl(100), fpl(npL), Lstar, EMfunc, bp, bb, &
       x, Planck,fplbol
!-----------------------------------------------------------------------

! is=1 for the enclosed source, is=2 for the external shell illumination
  if ((startyp(is).ge.4).and.(startyp(is).le.6)) then
   call readspectar(lambdas,Llamstar,Lstar,nLs,is,error)
! Generate dimensionless stellar spectrum
   do iLs = 1, nLs
    stellar(iLs) = lambdas(iLs)*Llamstar(iLs)/Lstar
   end do
  else
! if startyp.eq.3 generate power-law spectrum
   if (startyp(is).eq.3) then
    fl(1) = 1.0d0
    if (nLamtr(1).gt.1) then
     do i = 2, nLamtr(1)
      fl(i) = fl(i-1)*(lamtr(is,i-1)/lamtr(is,i))**klam(is,i-1)
     end do
    end if
    do iL = 1, nL
     if ((lambda(iL)-lamtr(is,1))*(lambda(iL)-lamtr(is,nLamtr(1)+1)).le.0.0d0) then
      kstop = 0
      k = 0
      do whiLe (kstop.eq.0)
       k = k + 1
       if (lambda(iL).ge.lamtr(is,k).and. &
            lambda(iL).le.lamtr(is,k+1)) then  ! This is Matt's correction for reading more than one powers:
        kstop = 1
        fpl(iL) = fl(k)*(lamtr(is,k)/lambda(iL))**klam(is,k)
        fpl(iL) = fpl(iL)/lambda(iL)
       end if
      end do
     else
      fpl(iL) = 0.0d0
     end if
    end do
   end if
  end if

  call Simpson(npL,1,nL,lambda,fpl,fplbol)

  do iY = 1, nY
! loop over wavelengths
   do iL = 1, nL
    if (startyp(is).eq.1) then
     bb = 0.0d0
     do k = 1, nbb(is)
      x = 14400.0d0/(lambda(iL)*Tbb(is,k))  ! hc/k = 14400 micron.Kelvin
      bb = bb + rellum(is,k)*Planck(x)
     end do
    else if (startyp(is).eq.2) then
     bb = EMfunc(lambda(iL),Tbb(is,1),xSiO)
    else if (startyp(is).eq.3) then
     bb = lambda(iL)*fpl(iL)/(fplbol)
!  for lambda longer than the longest entry in namestar
!  assume rayleigh-jeans tail
    else if (lambda(iL).gt.lambdas(nLs)) then
     bb = stellar(nLs) * (lambdas(nLs)/lambda(iL))**3.0d0
    else if (lambda(iL).lt.lambdas(1)) then
! if shorter than the shortest assume 0
     bb = 0.0d0
!  if within limits interpolate and startyp(is).ge.4
    else
     call powerinter(nLambdam,nLs,lambdas,stellar,lambda(iL),iLs,bp)
     bb = bp
    end if
    shp(iL) = bb
   end do
  end do
  error = 0
!-----------------------------------------------------------------------

 return
end subroutine getSpShape
!***********************************************************************

!***********************************************************************
subroutine Input(nameIn,nG,nameOut,nameQ,nameNK,tau1,tau2,tauIn, &
           Nrec,GridType,Nmodel,error,version,stdf)
!=======================================================================
! This subroutine reads input data from the file 'filename.inp'. It
! utilizes the function RDINP and subroutine RDINPS2 written by Moshe Elitzur.
!                                           [ZI,NOV'95; MN,JAN'00, MN'09]
!=======================================================================

  use common
  implicit none
  integer i, iG, nG, Nmodel, EtaOK, error,GridType, istop, nLs, Nrec, &
       ioverflw, Nmax, nLambdam, Nis, imu, geom, denstyp, ang_type, L
! Nmax is the size of user supplied eta file
! nLambdam is the max number entries for a user supplied stellar spectrum
  parameter (Nmax = 1000, nLambdam = 10000, nis = 2)
  double precision tau1, tau2, sum, a, b, xx(Nmax), e(Nmax), &
       aa(Nmax), bb(Nmax), tauIn(Nrec), ceta, x1, psf1, &
       lum, dist, var1, var2, var3, Lstar, th1, th2,      &
       lambdas(nLambdam), Llamstar(nLambdam), res, value, RDINP
  character lamstr(20)*72, strpow*72, strg*40, version*(*),stdf(7)*235
  character*(*) nameIn,nameOut, nameQ(npG),nameNK(10), namepsf*100, &
       nametau*100, anggrid*100,str*235
  logical Equal, noEqual, UCASE
!-----------------------------------------------------------------------

  UCASE = .true.
  Equal = .true.
  noEqual = .false.
  error = 0
  geom = 0

! Open output file
  open(12,file=nameOut,status='unknown')
  write(12,*)'==========================='
  write(12,*)' Output from program dusty '
  write(12,*)' version: ',version
  write(12,*)'==========================='
  write(12,*)' '
  write(12,*)' Input parameters from file: '
  write(12,'(2x,a140)')nameIn
  write(12,*)' '

! Open input file
  open(1,err=998,file=nameIn,status='old')
  rewind(1)

!********************************************
!** I. Geometry **
!********************************************
  call rdinps2(Equal,1,str,L,UCASE)
  if(str(1:L).eq.'SPHERE') then
   slb = .false.
   sph = .true.
   geom = 1
  elseif(str(1:L).eq.'SLAB') then
   slb = .true.
   sph = .false.
   geom = 0
  end if


!********************************************
!** II. Physical parameters **
!********************************************
! (1) Flags for presence of sources for slab these have the meaning
! of "left" and "right" source,
! illumination
  call rdinps2(Equal,1,str,L,UCASE)
  if(str(1:L).eq.'ON') then
   left = 1
  elseif(str(1:L).eq.'OFF') then
   left = 0
  end if

  if(left.eq.0) then
   if(slb) then
    call msg(23)
    left = 1
   end if
  end if

  call rdinps2(Equal,1,str,L,UCASE)
  if(str(1:L).eq.'ON') then
   right = 1
  elseif(str(1:L).eq.'OFF') then
   right = 0
  end if

! 1.1) CENTRAL SOURCE RADIATION (left-side source for slab)
  if(left.gt.0) then
   call inp_rad(error,1,nameIn)
   if(error.ne.0) goto 996
! typentry give the scale of input radiation
   call rdinps2(Equal,1,str,L,UCASE)
   if (str(1:L).eq.'FLUX') then
    typentry(1) = 1
   elseif (str(1:L).eq.'LUM_R1') then
    typentry(1) = 2
   elseif (str(1:L).eq.'ENERGY_DEN') then
    typentry(1) = 3
   elseif (str(1:L).eq.'DILUTN_FAC') then
    typentry(1) = 4
   elseif (str(1:L).eq.'T1') then
    typentry(1) = 5
   end if

! check if the entered value is acceptable
   if (typentry(1).lt.1.or.typentry(1).gt.5) then
    call msg(21)
    error = 1
    goto 999
   end if
   if (typentry(1).eq.1) then
    if (startyp(1).gt.3) then
!    for source spectrum in a file
!    get the scale of the input flux from the file
     call readspectar(lambdas,Llamstar,Lstar,nLs,1,error)
     dilutn = RDINP(Equal,1)
     Fint = dilutn*Lstar/pi
    else
!    typentry(1)=1: enter Fint, [W/m2]
!    N.B. in slab case Tei is calculated in this case
!    since the scaling flux is mu1*Fi,the local bol.flux Fint=L/(4*pi*r^2)
     Fint = RDINP(Equal,1)
    end if
    var1 = Fint
    Ji = Fint/(4.0d0*pi)
    Tei = (pi*Ji/sigma)**0.25d0
    var2 = Tei
   else if (typentry(1).eq.2) then
!   enter luminosity [in Lo] of the source and distance r1[cm] to the source
    Lum = RDINP(Equal,1)
    dist = RDINP(Equal,1)
! all units in dusty are in SI, so convert the input
    var1 = Lum
    var2 = dist
    Lum = Lum*3.862d+26
    dist = dist/100.0d0
    Fint = Lum/(4.0d0*pi*dist*dist)
    Ji = Fint/(4.0d0*pi)
    Cr1 = dist
    Tei = (Fint/sigma)**0.25d0
   else if (typentry(1).eq.3) then
    Ji = RDINP(Equal,1)
    var3 = Ji
   else if (typentry(1).eq.4) then
!   entry of dilution (normalization) factor
    dilutn = RDINP(Equal,1)
    if (startyp(1).gt.3) then
! get the scale of the input flux from the file
     call readspectar(lambdas,Llamstar,Lstar,nLs,1,error)
     Ji = dilutn*Lstar
    else
! for one Bbody Tstar=Tbb, for any other shape Dusty's default is Tstar=1e4 K.
     Ji = dilutn*sigma/pi*Tstar(1)**4.0d0
    end if
! var3 = dilutn
    var3 = Ji
    Fint = 4.0d0*pi*Ji
    Tei = (Fint/sigma)**0.25d0
   elseif (typentry(1).eq.5) then
!   enter dust temperature on inner boundary, T1[K]
    Tsub(1) = RDINP(Equal,1)
    var1 = Tsub(1)
! in this case Tei (or Fint) is determined from the Rad.Equillibrium condition
! at the first y-grid point. the ini.approximation is in
! InitTemp and Tei is adjusted in FindTemp for ea. Td-iteration.
   end if
  else
! if no central source
   Tei = 0.0d0
   typentry(1) = 0
! end if for presence of central(left) source
  end if

! for slab:
  if (slb) then
   write(12,'(a33)') ' Calculation in planar geometry:'
! find the kind of illumination
   call rdinps2(Equal,1,str,L,UCASE)
   if (str(1:L).eq.'DIRECTIONAL') then
   write(12,'(a41)') ' Directional illumination from the left.'
! enter incident theta_in:
! th1 the left illumination angle (in degrees) measured from the normal
   th1 = RDINP(Equal,1)
   elseif (str(1:L).eq.'ISOTROPIC') then
    th1 = -1.0d0
   end if
! for isotropic illumination
   if (th1.eq.-1.0d0) then
    write(12,'(a40)') ' Isotropic illumination from the left.'
! in this case th1=-1.0 is a flag for diffuse slab illumination,
! mu1 is set to -1.0 as a flag as well, not an actual value
    mu1 = -1.0d0
    if (typentry(1).eq.1.or.typentry(1).eq.2) then
     Ji = Fint/pi
    end if
! for directional illumination
   else
    write(12,'(a28,F5.2,a8)') ' Left illumination angle =',th1,' degrees'
    call chkangle(th1)
! convert to radians
    th1 = th1*pi/180.0d0
    mu1 = dcos(th1)
    if (typentry(1).eq.1.or.typentry(1).eq.2) then
     Ji = Fint/(4.0d0*pi)
    end if
   end if
   if (typentry(1).ge.2.and.th1.ne.-1.0d0) then
! for oblique illumination the input flux is Fint*mu1,
! scaling is with this flux. Since mu1 is read here
! Fint (resp. Tei) need recalculation:
! dabs(mu1) to accomodate the case with mu1=-1 (isotropic ill.)
! Fint = Fint*dabs(mu1)
    Tei = (pi*Ji/ sigma)**0.25d0
   end if
   var1 = Fint
   var2 = Tei
! even if no second source is supplied, mu2 needs a value
! of 1 to avoid crashing in the formulae in SLBStar [MN]
   mu2 = 1.0d0
   ksi = 0.0d0
  end if

! 1.2) external radiation (right-side source for slab)
  if (right.gt.0) then
   call inp_rad(error,2,nameIn)
   if (error.ne.0) goto 996
   if (slb) then
      call rdinps2(Equal,1,str,L,UCASE)
    if (str(1:L).eq.'DIRECTIONAL') then
       write(12,'(a42)') '  Directional illumination from the right.'
       th2 = RDINP(Equal,1)
    elseif (str(1:L).eq.'ISOTROPIC') then
       write(12,'(a41)') ' Isotropic illumination from the right.'
       th2 = -1.0d0
    end if
!  ksi is the relative bol.flux of the second source
    ksi = RDINP(Equal,1)
    if (ksi.lt.0.0) ksi = 0.0d0
    if (ksi.gt.1.0) ksi = 1.0d0
    write(12,'(a49,F5.2)') ' Relative bol.flux fraction of right source: R =',ksi
    if (th2.eq.-1.0d0) then
! in this case th2=-1.0 is a flag for diffuse slab illumination,
! mu2is set to -1.0 as a flag as well, not an actual value
     mu2 = -1.0d0
    else
     write(12,'(a29,F5.2,a8)') ' Right illumination angle =',th2,' degrees'
     call chkangle(th2)
! convert to radians
     th2 = th2*pi/180.0d0
     mu2 = dcos(th2)
    end if
   elseif(sph) then
! for sphere:
    call rdinps2(Equal,1,str,L,UCASE)
! typentry(2) can be 1(flux), 2(lum,r1) 3 (en_den), 4 Tsub(1)
    if (str(1:L).eq.'FLUX') then
     typentry(2) = 1
    elseif (str(1:L).eq.'LUM_R1') then
     typentry(2) = 2
    elseif (str(1:L).eq.'ENERGY_DEN') then
     typentry(2) = 3
    elseif (str(1:L).eq.'DILUTN_FAC') then
     typentry(2) = 4
    end if
    if (typentry(2).gt.4) then
     call msg(22)
     error = 1
     goto 999
    end if
    if (typentry(2).eq.1) then
     if (startyp(2).gt.3) then
      call readspectar(lambdas,Llamstar,Lstar,nLs,1,error)
      dilutn = RDINP(Equal,1)
      Fo = dilutn*Lstar/pi
     else
!  enter Fo, [W/m2]
       Fo = RDINP(Equal,1)
     end if
     Jo = Fo/pi
     Teo = (Fo/sigma)**0.25d0
    elseif(typentry(2).eq.2) then
! typEntry(2)=2: enter luminosity,[in Lo] of the source and
! distance r1,[cm] to the source
      lum = RDINP(Equal,1)
     dist = RDINP(Equal,1)
! all units in dusty are in SI, so convert the input
     lum = lum*3.862d+26
     dist = dist/100.0d+00
!  Cr1 is carried in common/dyn/
     Cr1 = dist
     Fo = lum/(4.0d+00*pi*dist*dist)
     Teo = (Fo/sigma)**0.25d0
     Jo = Fo/pi
    elseif (typentry(2).eq.3) then
     Jo = RDINP(Equal,1)
     Fo = pi*Jo
     Teo = (Fo/sigma)**0.25d0
    elseif(typentry(2).eq.4) then
! typentry(2)=4 for entry of dilution (normalization) factor
     dilutn = RDINP(Equal,1)
     if (startyp(2).gt.3) then
! get the scale of the input flux from the file
     call readspectar(lambdas,Llamstar,Lstar,nLs,2,error)
     Jo = dilutn*Lstar
    else
! for one blackbody Tstar=Tbb, for any other shape
! dusty's default is Tstar =1e4 K.
     Jo = dilutn*sigma/pi*Tstar(2)**4.0d0
    end if
!   var3 = dilutn
    Fo = pi*Ji
    Teo = (Fo/sigma)**0.25d0
   end if
   var3 = Jo
! end if for geometry
   end if
  else
! if no external radiation -- set to 0
   Teo = 0.0d0
   typentry(2) = 0
! end if for second source
  end if
!! line added for clarity in the .out file [MN]
  write(12,*) ' --------------------------------------------'

!=========  END READING OF SOURCE PARAMETERS ===================

! (2) DUST PROPERTIES
! # of different dust grains, to be used in a future version
  nG = 1
! 2.1 Chemical composition
! Type of optical properties
  call rdinps2(Equal,1,str,L,UCASE)
  if (str(1:L).eq.'COMMON_GRAIN') then
   top = 1
  elseif (str(1:L).eq.'COMMON_AND_ADDL_GRAIN') then
   top = 2
  elseif (str(1:L).eq.'TABULATED') then
   top = 3
  end if
  if (top.ne.1.and.top.ne.2.and.top.ne.3) then
   call msg(9)
   error = 1
   goto 999
  end if
! For top.lt.3 read in abundances for supported grains
  if (top.lt.3) then
     xC(1) = RDINP(Equal,1)
     if (xC(1).lt.0.0d0) xC(1) = 0.0d0
     sum = xC(1)
     do i = 2, 7
! Special care to be taken of graphite (1/3-2/3 rule):
       if (i.ne.5) then
          xC(i) = RDINP(noEqual,1)
          if (xC(i).lt.0.0d0) xC(i) = 0.0d0
!        i Equal 4 is data for graphite (parallel to c axis):
          if(i.eq.4) xC(i) = 1.0d0*xC(i)/3.0d0
      else
! graphite (perpendicular to c axis) :
          xC(i) = 2.0d0 * xC(i-1)
      end if
      sum = sum + xC(i)
     end do
  end if

! Assign supported dust filenames to stdf
  do i = 1,7
   if (i.eq.1) write(stdf(i),'(a)')"stnd_dust_lib/OssOdef.nk"
   if (i.eq.2) write(stdf(i),'(a)')"stnd_dust_lib/OssOrich.nk"
   if (i.eq.3) write(stdf(i),'(a)')"stnd_dust_lib/sil-dlee.nk"
   if (i.eq.4) write(stdf(i),'(a)')"stnd_dust_lib/gra-par-draine.nk"
   if (i.eq.5) write(stdf(i),'(a)')"stnd_dust_lib/gra-perp-draine.nk"
   if (i.eq.6) write(stdf(i),'(a)')"stnd_dust_lib/amC-hann.nk"
   if (i.eq.7) write(stdf(i),'(a)')"stnd_dust_lib/SiC-peg.nk"
  enddo
! user supplied n and k:
  if (top.eq.2) then
     nfiles = RDINP(Equal,1)
! File names
   strg = 'optical constants:'
   do i = 1, nfiles
    call filemsg(nameNK(i),strg)
   end do
   if(error.ne.0) goto 996
! Abundances
   xCuser(1) = RDINP(Equal,1)
   if (xCuser(1).lt.0.0d0) xCuser(1) = 0.0d0
   sum = sum + xCuser(1)
   if (nfiles.gt.1) then
    do i = 2, nfiles
     xCuser(i) = RDINP(noEqual,1)
     if (xCuser(i).lt.0.0d0) xCuser(i) = 0.0d0
     sum = sum + xCuser(i)
    end do
   end if
  end if
  if (top.lt.3) then
   if (sum.le.0.0d0) then
    call msg(5)
    error = 1
    goto 999
   end if
! Normalize abundances for supported grains:
   do i = 1, 7
    xC(i) = xC(i) / sum
   end do
! Normalize abundances for user supplied grains
   if (top.eq.2) then
    do i = 1, nfiles
     xCuser(i) = xCuser(i) / sum
    end do
   end if
  end if
! user supplied cross-sections:
  if (top.eq.3) then
! filename for qabs and qsca
   strg= 'abs. and scatt. cross-sections:'
   do iG = 1, nG
    call filemsg(nameQ(iG),strg)
   end do
  end if
! 2.2 Grain size distribution
  if (top.ne.3) then
! Type of size distribution
   call rdinps2(Equal,1,str,L,UCASE)
   if (str(1:L).eq.'MRN') then
    szds = 1
   elseif (str(1:L).eq.'MODIFIED_MRN') then
    szds = 2
   elseif (str(1:L).eq.'KMH') then
    szds = 3
   end if

   if (szds.ne.1.and.szds.ne.2.and.szds.ne.3) then
    call msg(10)
    error = 1
    goto 999
   end if
! Grain sizes
   if (szds.gt.1) then
      qsd = RDINP(Equal,1)
      a1 = RDINP(Equal,1)
    if (a1.le.0.0) a1 = 0.0001d0
      a2 = RDINP(Equal,1)
    if (szds.eq.2.and.a2.lt.a1) a2 = a1
   else
    qsd = 3.5d0
    a1 = 0.005d0
    a2 = 0.25d0
   end if
  end if
!=========  END READING DUST PROPERTIES ===================

! WriteOut prints all input data, read so far, in fname.out
! var1 is t1,fe1,luminosity or teff; var2 is r1; var3 is ext.rad. input
  call WriteOut(var1,var2,var3,nG,nameQ,nameNK)

! (3) Density distribution
! For sphere only:
  if(sph) then
   powd  = .false.
   expd  = .false.
   rdw   = .false.
   rdwa  = .false.
   fild  = .false.
   rdwpr = .false.
! Parameter describing eta function:
   call rdinps2(Equal,1,str,L,UCASE)
   if (str(1:L).eq.'POWD') then
    powd = .true.
    denstyp = 1
   elseif (str(1:L).eq.'EXPD') then
    expd = .true.
    denstyp = 2
! *** Winds ***
! denstyp.eq.3 is RDW with default values of v1/ve=0.2, GravCor=0.5
! denstyp.eq.6 is a private option with additional input for v1/ve and
! GravCor=max(Fgrav/Frad);
   elseif (str(1:L).eq.'RDW') then
    rdw = .true.
    denstyp = 3
! analytical (gray) approximation for rdw
  elseif (str(1:L).eq.'RDWA') then
    rdwa = .true.
    denstyp = 4
! file with user supplied density distribution
!!!   elseif (str(1:L).eq.'USR_SUPPLD') then
    elseif (str(1:L).eq.'USER_SUPPLIED') then
    fild = .true.
    denstyp = 5
!  private option for RDW with additional output
   elseif (str(1:L).eq.'RDWPR') then
    rdwpr = .true.
    denstyp = 6
   end if
! initialize EtaOK and Ntr
   EtaOK = 0
   Ntr = 0
! read parameters for each type of density distribution
! smooth or broken power laws
   if (powd) then
    EtaOK = 1
    Ntr = RDINP(Equal,1)
! changed definition
    Ntr = Ntr - 1
! read in transition radii
    if (Ntr.gt.0) then
      Ytr(1) = RDINP(Equal,1)
      if (Ntr.gt.1) then
        do i = 2, Ntr
          Ytr(i) = RDINP(NoEqual,1)
        end do
      end if
      Yout = RDINP(noEqual,1)
    else
! for smooth density power law
      Yout = RDINP(Equal,1)
    end if
    if (Yout.le.1.0d0) Yout = 1.001d0
! read in powers
    pow = RDINP(Equal,1)
    if (Ntr.gt.0) then
      do i = 1, Ntr
         ptr(i) = RDINP(NoEqual,1)
      end do
    end if
! print info to the output file
    if (Ntr.eq.0) then
     call getfs(pow,2,0,strpow)
     write(12,'(a38,a5)') ' density described by 1/r**k with k =',strpow
     write(12,'(a21,1p,e10.3)')'  relative thickness:',Yout
    else
     write(12,*)' density described by a broken power law:'
     write(12,*)'  power   Ytransition'
     write(12,*)'  -------------------'
     write(12,*)'              1.0'
     call getfs(pow,2,0,strpow)
     write(12,'(a2,a5)')'  ',strpow
     do i = 1, Ntr
      write(12,'(a10,1p,e10.3)')'          ',Ytr(i)
      call getfs(ptr(i),2,0,strpow)
      write(12,'(a2,a5)')'  ',strpow
     end do
     write(12,'(a10,1p,e10.3)')'          ',Yout
    end if
   end if

! exponential law
   if (expd) then
     EtaOK = 1
     Yout = RDINP(Equal,1)
     if (Yout.le.1.0d0) Yout = 1.001d0
       pow = RDINP(Equal,1)
       if (pow.le.0.0d0) then
       EtaOK = 0
     else
       write(12,*)' density described by exponential distribution'
       write(12,'(a21,1p,e10.3)')'               sigma:',pow
       write(12,'(a21,1p,e10.3)')'  relative thickness:',Yout
     end if
   end if
! default approximation and default numerics for rad. driven winds
   if (rdwa.or.rdw) then
     EtaOK = 1
     Yout = RDINP(Equal,1)
     if (Yout.le.1.0d0) Yout = 1.001d0
! ** default ** for epsilon = v1/ve = u1/ue:
     pow = 0.2d0
     if(rdw) then
! ** default ** for max(gravcor = fgrav/frad_press):
      ptr(1) = 0.5d0
! convergence criterion:
      ptr(2) = 1.0d0
! default linear version of the eq. for velocity
      ver = 1
     end if
     write(12,*)' Density for radiatively driven winds from'
     if (rdwa) then
       write(12,*)' Analytic approximation for gray dust.'
     else
       write(12,*)' Full dynamic calculation.'
     end if
     write(12,'(a21,1p,e10.3)')'  Relative thickness:',Yout
   end if
! full dynamical calculation for radiatively driven winds (private option)
! the user can specify parameters that have default values in denstyp=3
! user specified table for eta
   if(fild) then
    EtaOK = 1
    strg = 'Dust density distribution:'
    call filemsg(nameeta,strg)
    write(12,*)' Density distribution supplied from file:'
    write(12,'(2x,a100)') nameeta
    call prHeader(3,nameeta)
! read in the density
    open(26,err=997,file=nameeta,status='old')
! three lines in the header:
    do i = 1, 3
      read(26,*,err=997) strpow
    end do
    istop = 0
    i = 0
    do while (istop.ge.0)
     read(26,*,end=900,err=997,iostat=istop) a, b
     if (istop.ge.0) then
      i = i + 1
      xx(i) = a
      e(i) = b
      if (i.eq.1) x1 = xx(i)
      yetaf(i) = xx(i) / x1
     end if
    end do
900 close(26)
    nYetaf = i
    if (nYetaf.lt.2) goto 997
! if input positions in descending order turn them around
    if (yetaf(1).gt.yetaf(2)) then
     do i = 1, nYetaf
      aa(i) = yetaf(i)
      bb(i) = e(i)
     end do
     do i = 1, nYetaf
      yetaf(i) = aa(nYetaf+1-i)
      e(i) = bb(nYetaf+1-i)
     end do
    end if
! relative thickness
    Yout = yetaf(nYetaf)
    write(12,'(a21,1p,e10.3)')'  relative thickness:',Yout
    if (Yout.le.1.0d0) Yout = 1.001d0
! integrate and ...
    call Simpson(Nmax,1,nYetaf,yetaf,e,ceta)
! ... renormalize
    do i = 1, nYetaf
     etaf(i) = e(i) / ceta
    end do
   end if
! Done with the reading of density distribution
   if (EtaOK.ne.1) then
    call msg(3)
    error = 1
    goto 999
   end if
   write(12,*)' --------------------------------------------'
  end if
!=========  End reading density distribution ===================

! 4) Optical depth
! Grid type
  call rdinps2(Equal,1,str,L,UCASE)
  if (str(1:L).eq.'LINEAR') then
   GridType = 1
  elseif (str(1:L).eq.'LOGARITHMIC') then
   GridType = 2
  elseif (str(1:L).eq.'USER_SUPPLIED') then
   GridType = 3
  end if
  if (GridType.eq.3) then
! tau-grid from a file
    strg = 'user supplied tau-grid:'
    call filemsg(nametau,strg)
! read optical depths
    open(27,err=992,file=nametau,status='old')
! fiducial wavelength
! (the second argument of rdinp is the unit)
    lamfid = RDINP(Equal,27)
! number of models in the list
    Nmodel = RDINP(Equal,27)
    do i = 1, Nmodel
      read(27,*) tauIn(i)
    end do
902 close(27)
! Sort the tau-grid if there is more than one model:
    if(Nmodel.gt.1) then
      call sort(tauIn,Nmodel)
    end if
    tau1 = tauIn(1)
    if (tau1.le.0.0d0) tau1 = 0.0001d0
    tau2 = tauIn(Nmodel)
  else
! fiducial wavelength
      lamfid = RDINP(Equal,1)
! total optical depths at lamfid
      TAU1 = RDINP(Equal,1)
   if (tau1.le.0.0d0) tau1 = 0.0001d0
      TAU2 = RDINP(Equal,1)
   if (tau2.le.tau1) then
    tau2 = tau1
    Nmodel = 1
   end if
! read number of models
     Nmodel = RDINP(Equal,1)
! Nrec = 1000, initialized in MAIN
   if (Nmodel.gt.(Nrec-1)) Nmodel = Nrec-1
   if (Nmodel.lt.1) Nmodel = 1
  end if
  if (Nmodel.gt.1) then
   write(12,'(a19,1p,e8.1,a8)')' Optical depths at',lamfid, 'microns'
   write(12,'(a14,1p,e9.2,a3,e9.2)')' ranging from',tau1,' to',tau2
   if (GridType.eq.1) strg=' models with linear grid    '
   if (GridType.eq.2) strg=' models with logarithmic grid'
   if (GridType.eq.3) strg=' models with grid from file  '
   write(12,'(a1,i4,a)')' ', Nmodel, strg
   if (GridType.eq.3) write(12,'(a4,a70)')'    ',nametau
  else
   write(12,'(a18,1p,e8.1,a9,e9.2)')' Optical depth at',lamfid, ' microns:',tau1
  end if

! 5) disk
! for disk calculations make sure you have npX>1 in 'userpar.inc' !!
  if (npX.gt.1) then
   if (iVerb.ge.1) write(*,*) 'No disk option in this version.'
   goto 999
  end if

!********************************************
!** III. Numerical accuracy **
!********************************************
! accuracy for convergence (typical 0.0001)
!!** this is accConv for dust temperature
!!**  accConv = 10.0d-3  !! Too rough.
  accConv = 1.0d-4  !tests on July,8,2010  1e-4 and 1e-6 give very close results
! accuracy for flux conservation
  accuracy = RDINP(Equal,1)
  if (accuracy.le.0.0d0) accuracy = 0.02d0
! Protect against a very large value for accuracies
  if (accuracy.gt.0.25d0) accuracy = 0.25d0
! starting optical depth
  init_tau = RDINP(Equal,1)
! increment in optical depth
  dtau = RDINP(Equal,1)
! dynamical range
  dynrange = 1.0d-15
  if (accuracy.ge.0.1d0) then
   call getfs(accuracy*1000.0d0,0,1,strpow)
   write(12,'(a20,a3,a1)')' Required accuracy:',strpow,'%'
  else
   call getfs(accuracy*100.0d0,0,1,strpow)
   write(12,'(a20,a2,a1)')' Required accuracy:',strpow,'%'
  end if
  write(12,*)' --------------------------------------------'

!********************************************
!** IV. Output flags **
!********************************************
! Internal flag for additional miscellaneous output  [MN]:
! if iInn=1: print err.vs.iter in unt=38 (fname.err) for all models
! and additionally list scaled fbol(y) and ubol(y) in m-files.
  iInn = 1
! spectral properties
  iSPP = RDINP(Equal,1)
!  spectra
  iA = RDINP(Equal,1)
  iC = RDINP(Equal,1)
! images (intensity)
  if (iC.ne.0) then
    if (slb) then
! Read angular grid (this is theta_out) for slab intensity output.
! the output intensities are in units of lambda*I_lambda*cos(theta_out)/Fe
! where Fe=L/(4*pi*r^2), the local bolometric flux.
     ang_type = RDINP(Equal,1)
!    Create the grid depending on grid type
!    1-equidistant in theta, 2-equidistant in cos(theta), 3-from a file
     call input_slb_ang(ang_type)
!    Convert to radians
     do imu = 1, nmu
      theta(imu) = theta(imu)*pi/180.0d0
     end do
     iV = 0
     iPsf = 0
    else
! for spherical case
     NlambdaOut = RDINP(Equal,1)
     if (nLambdaOut.ge.1) then
      do i = 1, nLambdaOut
         LambdaOut(i) = RDINP(NoEqual,1)
!  make sure the wavelengths are inside dusty's range
        if (LambdaOut(i).le.0.01d0) LambdaOut(i) = 0.01d0
        if (LambdaOut(i).gt.36000.0d0) LambdaOut(i) = 36000.0d0
      end do
      ioverflw = 0
      do i = 1, nLambdaOut
        if (LambdaOut(i).lt.0.995d0) then
           call getfs(LambdaOut(i),2,1,lamstr(i))
        else
          if (LambdaOut(i).lt.9.95d0) then
            call getfs(LambdaOut(i),1,0,lamstr(i))
          else
             if (LambdaOut(i).lt.99.5d0) then
               call getfs(LambdaOut(i),0,0,lamstr(i))
             else
               call getfs(LambdaOut(i),0,1,lamstr(i))
             end if
          end if
        end if
        if (LambdaOut(i).gt.9999.5d0) then
          ioverflw = 1
          strpow = lamstr(i)
          strpow(4:4) = '*'
          strpow(5:5) = ' '
          lamstr(i) = strpow
        end if
      end do
     end if
    write(12,*)' Images requested for these wavelengths (mic)'
    write(12,'(a1,20a5)')' ',(lamstr(i),i=1,nLambdaOut)
    if (ioverflw.eq.1) write(12,*)'  *: in mm'
! Convolved images  (only for our use)
    if (iC.lt.0) then
     iPsf = 1
! iPsf = rdinp(Equal,1)
     if (iPsf.ne.0) then
      Theta1 = RDINP(Equal,1)
      write(12,'(a39,1p,e7.1)') ' Convolved images produced for theta1=',theta1
       psftype = RDINP(Equal,1)
      if (psftype.ne.1.and.psftype.ne.2.and.psftype.ne.3) goto 994
      if (psftype.lt.3) then
! Gaussians, read in parameters
! FWHM for the first component
       FWHM1(1) = RDINP(Equal,1)
       if (nLambdaOut.gt.1) then
        do i = 2, nLambdaOut
           FWHM1(i) = RDINP(NoEqual,1)
        end do
       end if
       if (psftype.eq.2) then
! Relative strength for the second component
        kPSF(1) = RDINP(Equal,1)
        if (nLambdaOut.gt.1) then
         do i = 2, nLambdaOut
           kPSF(i) = RDINP(NoEqual,1)
         end do
        end if
!       FWHM for the second component
        FWHM2(1) = RDINP(Equal,1)
        if (nLambdaOut.gt.1) then
         do i = 2, nLambdaOut
           FWHM2(i) = RDINP(NoEqual,1)
         end do
        end if
       end if
       write(12,*)' the point spread functions are gaussians'
      else
! user supplied psf
       strg = 'point spread function:'
       call filemsg(namepsf,strg)
       write(12,*)' the point spread function supplied from file'
       write(12,'(2x,a100)')namepsf
       open(28,err=995,file=namepsf,status='old')
! Three lines in the header:
       do i = 1, 3
        read(28,*,err=995)
       end do
       istop = 0
       i = 0
       do while (istop.ge.0)
        read(28,*,end=901,err=995,iostat=istop)a, b
        if (istop.ge.0) then
         i = i + 1
         if (i.eq.1) then
          psf1 = b
          if (a.ne.0.0d0) goto 995
         end if
         xpsf(i) = a
         ypsf(i) = b / psf1
        end if
       end do
901    close(28)
       Npsf = i
!      scale to 1 at the center. This is only to get FWHM here.
!      ypsf is normalized to area in Subroutine Convolve [MN]
       call scaleto1(1000,npsf,ypsf)
!      Find equivalent fwhm
       istop = 0
       i = 1
       do while (istop.eq.0)
        i = i + 1
        if (ypsf(i).le.0.5d0) istop = 1
       end do
!      Linear interpolation
       FWHM1(1) = (xpsf(i)-xpsf(i-1))/(ypsf(i)-ypsf(i-1))
       FWHM1(1) = (fwhm1(1)*(0.5d0-ypsf(i-1))+xpsf(i-1))*2.0d0
       FWHM2(1) = 0.0d0
       write(12,'(a18,1p,e8.1)')' equivalent FWHM:',FWHM1(1)
! end if for psf from a file
      end if
! end if for psf
     end if
! end if for convolved images
    end if

! visibility (only if the intensity is requested)
    iV = RDINP(Equal,1)
    if(iV.ne.0) iV = abs(iC)
! end if for geometry
   end if
   write(12,*)' --------------------------------------------'
  else
! if iC=0 set the other flags to 0 (just in case).
   iPsf = 0
   iV = 0
   write(12,*)' --------------------------------------------'
  end if
! radial quantities
  iB = RDINP(Equal,1)
! run-time messages
  iX = RDINP(Equal,1)
! *** DONE READING INPUT PARAMETERS ***
! if everything is ok, close the input file and finish
999 goto 996
! or in the case of err reading files...
920 write(12,*)' ***  FATAL ERROR IN DUSTY  *************'
  write(12,*)' File with user supplied angular grid:   '
  write(12,*)'     slab_ang_grid.dat                   '
  write(12,*)' is missing or not properly formatted?!  '
  write(12,*)' ****************************************'
  close(12)
  error = 3
992 write(12,*)' ***  FATAL ERROR IN DUSTY  *************'
  write(12,*)' File with user supplied TAU-grid:       '
  write(12,'(2x,a100)') nameTAU
  write(12,*)' is missing or not properly formatted?!  '
  write(12,*)' ****************************************'
!  close(12)
  error = 3
  goto 996
994 call MSG(12)
!  close(12)
  error = 3
  goto 996
995 write(12,*)' ***  FATAL ERROR IN DUSTY  *************'
  write(12,*)' File with the point spread function:    '
  write(12,'(a2,a100)')'  ', namePSF
  write(12,*)' is missing or not properly formatted?!  '
  write(12,*)' ****************************************'
!  close(12)
  error = 3
997 write(12,*)' ***  FATAL ERROR IN DUSTY  *************'
  write(12,*)' File with the dust density distribution:'
  write(12,'(2x,a100)') nameETA
  write(12,*)' is missing or not properly formatted?!  '
  write(12,*)' ****************************************'
!  close(12)
  error = 3
998 write(12,*)' ***  FATAL ERROR IN DUSTY  ****'
  write(12,*)' Input file:'
  write(12,'(2x,a100)') nameIn
  write(12,*)' is missing?!'
  write(12,*)' *******************************'
!  close(12)
  error = 3
!-----------------------------------------------------------------------
996  close(1)
  return
end subroutine Input
!***********************************************************************


!***********************************************************************
subroutine inp_rad(error,is,nameIn)
!=======================================================================
! This is the former SUBROUTINE InpStar(error,is,nameIn)
! This subroutine is for reading the input radiation parameters
!                                                              [MN,Mar'99]
!=======================================================================
  use common
  implicit none

  integer error, i, is, l
  double precision sum, tsum, value, RDINP
  character strg*40, str*235
  character*(*) nameIn
  logical Equal, noEqual, UCASE
!-------------------------------------------------------------------------
  UCASE = .true.
  Equal = .true.
  noEqual = .false.
  error = 0

! Generic temperature for cases other than bb or Engelke-Marengo shape.
  Tstar(is) = 10000.0d0
! Flag for the external spectrum
  call rdinps2(Equal,1,str,L,UCASE)
  if (str(1:L).eq.'BLACK_BODY') then
   startyp(is) = 1
  elseif(str(1:L).eq.'ENGELKE_MARENGO') then
   startyp(is) = 2
  elseif(str(1:L).eq.'POWER_LAW') then
   startyp(is) = 3
  elseif(str(1:L).eq.'FILE_LAMBDA_F_LAMBDA') then
   startyp(is) = 4
  elseif(str(1:L).eq.'FILE_F_LAMBDA') then
   startyp(is) = 5
  elseif(str(1:L).eq.'FILE_F_NU') then
   startyp(is) = 6
  end if

! (1) Black body(ies) for startyp=1
  if(startyp(is).eq.1) then
! Number of black bodies
   nBB(is) = RDINP(Equal,1)
! Stellar temperature(s)
   Tbb(is,1) = RDINP(Equal,1)
   if (Tbb(is,1).le.0.0d0) then
    call msg(8)
    error = 1
    goto 999
   end if
! Single black body
   if (nbb(is).eq.1) then
    Tstar(is) = Tbb(is,1)
!    relative luminosity
    rellum(is,1) = 1.0d0
   endif  !end if for one bb
! Multiple black bodies
   if (nbb(is).gt.1) then
    do i = 2, nbb(is)
     Tbb(is,i) = RDINP(NoEqual,1)
     if (Tbb(is,i).le.0.0d0) then
      call msg(8)
      error = 1
      goto 999
     end if
    end do
! Read in relative luminosities
    rellum(is,1) = RDINP(Equal,1)
    sum = rellum(is,1)
    do i = 2, nbb(is)
     rellum(is,i) = RDINP(NoEqual,1)
     sum = sum + rellum(is,i)
    end do
    if (sum.le.0.0d0) then
     call msg(7)
     error = 1
     goto 999
    end if
! Normalize
    tsum = 0.0d0
    do i = 1, nbb(is)
     rellum(is,i) = rellum(is,i)/sum
     tsum = tsum + rellum(is,i)*Tbb(is,i)**(4.0d0)
    end do
    Tstar(is) = (tsum)**(0.25d0)
   end if ! end if for multiple bb
  end if  ! end if for bb-type or startype 1

! (2) engelkd-marengo function for startyp=2
  if(startyp(is).eq.2) then
! Effective stellar temperature
    Tbb(is,1) = RDINP(Equal,1)
    Tstar(is) = Tbb(1,1)
! Depth of SiO abs.feature in %
    xSiO = RDINP(Equal,1)
   if (xSiO.le.0.0d0) xSiO = 0.0001d0
   if (xSiO.gt.100.0d0) xSiO = 100.0d0
  end if

! (3) Power-law(s) for startyp=3
  if(startyp(is).eq.3) then
! Number of transitions
   Nlamtr(is)= RDINP(Equal,1)
   if (nLamtr(is).gt.0) then
    lamtr(is,1) = RDINP(Equal,1)
    if (nLamtr(is).gt.1) then
     do i = 2, nLamtr(is)+1
      lamtr(is,i) = RDINP(NoEqual,1)
      if (lamtr(is,i).lt.lamtr(is,i-1)) then
       call msg(6)
       error = 1
       goto 999
      end if
     end do
    endif
    klam(is,1) = RDINP(Equal,1)
    if (nLamtr(is).gt.1) then
     do i = 2, nLamtr(is)
       klam(is,i) = RDINP(NoEqual,1)
     end do
    end if
   else
    startyp(is) = 1
    Tstar(is) = 10000.0d0
   end if
  end if

 if (startyp(is).ge.4.and.startyp(is).le.6) then
  strg = 'Spectral shape of external radiation:'
  call filemsg(namestar(is),strg)
 end if
!-----------------------------------------------------------------------
999 return
end subroutine inp_rad
!***********************************************************************


!***********************************************************************
subroutine input_slb_ang(ang_type)
!=======================================================================
! This subroutine reads the set of input or output illumination angles
! for slab case.                                           [MN, 2005]
! =======================================================================
  use common
  implicit none

  integer ang_type, error, imu, length
  double precision th_min, th_max, angstep, cth_min, cth_max, caux, value, RDINP
  character*70 anggrid, strg,str*235
  logical Equal, noEqual,UCASE
!-----------------------------------------------------------------------

  Equal = .true.
  noEqual = .false.

! if ang_type=1 (equidistant in theta, given min,max,step)
  if (ang_type.eq.1) then
    th_min = RDINP(Equal,1)
    call chkangle(th_min)
    th_max = RDINP(Equal,1)
    call chkangle(th_max)
    if (th_max.le.th_min) then
     th_max = th_min
     nmu = 1
    end if
! step equidistant in theta
    angstep = RDINP(Equal,1)
!   create the grid:
    imu = 1
    theta(1)=th_min
    do while(theta(imu).lt.th_max)
     theta(imu+1) = theta(imu) + angstep
     imu = imu+1
    end do
    nmu = imu
    theta(nmu) = th_max
  end if

! if ang_type=2 (equidistant in cos theta)
  if (ang_type.eq.2) then
    th_min = RDINP(Equal,1)
    call chkangle(th_min)
    th_max = RDINP(Equal,1)
    call chkangle(th_max)
    cth_min = dcos(th_min*pi/180.0d0)
    cth_max = dcos(th_max*pi/180.0d0)
!   Step, equidistant in cos(theta)
    AngStep = RDINP(Equal,1)
!   Create the grid:
    imu = 1
    caux = cth_min
    do while (caux.gt.0.0d0)
     theta(imu) = dacos(caux)*180.0d0/pi
     caux = caux - angstep
     imu = imu+1
    end do
    Nmu = imu
    theta(Nmu) = th_max
  end if

  if (ang_type.eq.3) then
! Angular grid from a file, angles in degrees
   call FileMSG(ANGgrid,strg)
   open(7,ERR=92,file=ANGgrid,STATUS='OLD')
   Nmu = RDINP(Equal,7)
   do imu = 1, Nmu
     read(7,*) theta(imu)
   end do
92 close(7)
  end if
  write(12,*)' Intensity requested for these theta_out(deg):'
  write(12,'(a1,8f7.1,/,x,8f7.1,/,x,10f7.1)')' ', (theta(imu), imu = 1, nmu)
  if (ang_type.eq.1) write(12,'(a34,f4.1)') '  equidistant in theta_out, step=', angstep
  if (ang_type.eq.2) write(12,'(a39,f4.1)') '  equidistant in cos(theta_out), step=', angstep
  if (ang_type.eq.3) write(12,'(a18,a70)') '  grid from file: ', anggrid
! -----------------------------------------------------------------------
  return
end subroutine input_slb_ang
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
      END
! ***********************************************************************

!***********************************************************************
subroutine MakeTable(Elems,rows,cols,unt)
! =======================================================================
!     This is an auxiliary subroutine for print out of tables
!     of Elems(Nrows,Ncols) in output unit 'unt'. Nrows = max{npL,npY}.
!     The array size is defined in Subroutine PrOut.         [MN, Mar'98]
! =======================================================================
      IMPLICIT NONE
      INTEGER rows, cols, unt, k, i
      DOUBLE PRECISION Elems(rows,cols)
! -----------------------------------------------------------------------
      DO k = 1, rows
        write(unt,'(1p,21E11.3)') (Elems(k,i),i=1,cols)
      END DO
! -----------------------------------------------------------------------
      RETURN
END subroutine MakeTable
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
subroutine OPPEN(model,rootname,length)
!=======================================================================
! This subroutine prints the results out.              [Z.I., Feb. 1996]
!=======================================================================

  use common
  implicit none
  character ch5*5, rootname*(*), fname*235
  character*72 header1, s3, s4
  integer model, length, i
!-----------------------------------------------------------------------

! set up the status indicators
  iError = 0
  iWarning = 0
  if (model.eq.1) iCumm = 0
! the following files pertain to all models and are open if model.eq.1
  if (model.eq.1) then
! the header to output file *.out is moved to prout [mn,sep'99]

! open file with spectral properties rootname.spp
   if (iSPP.ne.0) then
    call attach(rootname,length,'.spp',fname)
    open(19,file=fname,status='unknown')
    if(slb) then
     write(19,'(a49)') '# ==============================================='
     write(19,'(a49)') '# properties of spectra from the slab right side '
     write(19,'(a49)') '# -----------------------------------------------'
    else
     call line(1,2,19)
     write(19,'(a23)')'#  spectral properties '
     call line(1,1,19)
    end if
    s3='###   tau0      Psi      fV       fK       f12    C21  '
    s4=' C31   C43  b8-13 b14-22 B9.8 B11.4  r9.8-18  '
    write(19,'(a55,a46)')s3,s4
    if (slb.and.iSPP.eq.3) then
     call attach(rootname,length,'.zpp',fname)
     open(24,file=fname,status='unknown')
     write(24,'(a49)') '# ==============================================='
     write(24,'(a49)')'# properties of spectra from the slab left side  '
     write(24,'(a49)') '# -----------------------------------------------'
     s3='###   tau0      Psi      fV       fK       f12    C21  '
     s4=' C31   C43  b8-13 b14-22 B9.8 B11.4  R9.8-18  '
     write(24,'(a55,a46)')s3,s4
    end if
   end if
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
   if (abs(ic).eq.1) then
    call attach(rootname,length,'.itb',fname)
    open(17,file=fname,status='unknown')
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
!     open the disk files rootname.d##, rootname.e## and rootname.w##
!      if (iD.ge.1) then
!        write(ch5,'(a2,i3.3)') '.d', model
!        call attach(rootname,length,ch5,fname)
!        open(30,file=fname,status='unknown')
!        if (iD.gt.1) then
!          write(ch5,'(a2,i3.3)') '.e', model
!          call attach(rootname,length,ch5,fname)
!           open(31,file=fname,status='unknown')
!          write(ch5,'(a2,i3.3)') '.w', model
!          call attach(rootname,length,ch5,fname)
!           open(32,file=fname,status='unknown')
!        end if
!      end if
!-----------------------------------------------------------------------
  return
end subroutine OPPEN
!***********************************************************************


!***********************************************************************
      SUBROUTINE Pass_Header(iunit)
!     Get past Header lines of opened file on unit # iunit
!     The end of headers is marked by a line that contains just ">"
!=======================================================================
      implicit none

      integer iunit
      character*128 header
!-----------------------------------------------------------------------
      header = "XXX"
      do while(.not. (header(1:1).eq.">"))
         read(iunit, '(a)') header
      end do
!-----------------------------------------------------------------------
      return
end subroutine Pass_Header
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


!***********************************************************************
subroutine PrOut(model,nG,delta)
!=======================================================================
! This subroutine prints the results out.        [ZI,Feb'96; MN,Mar'99]
!=======================================================================

  use common
  implicit none
  integer iY, iL, i, model, j, unt, imu, nrows, ncols,nG
!  parameter (nrows=200, ncols=25)
  double precision psfn, psffunc(20,1000),eta,faux(npL), &
       tht1, xs, xds, xde, res, fnorm, var1, var2, dmax, &
       limval, GinfG1, omega(npG,npL), Jext(npY), delta
  double precision, allocatable::Elems(:,:)
  character*72 SC21,SC31,SC43,SB98,SB11,Sbet1,Sbet2,STemp,Serr,hdint, hdcon, &
       hdvis, s1, su1, s2, su2, tstr*10
  character*132 hdsp1,hdsp2,hdrslb1,hdrslb2,hdrsph1,hdrsph2,hdrdyn
!----------------------------------------------------------------------

  allocate(Elems(npL,8))
!FH 10/18/10 line Ji=Fi/(4*pi) added, because Ji not yet set
  Ji = Fi / (4*pi)
  IF(iPhys.eq.1) THEN
     DO iY = 1, nY
        IF (sph) THEN
           Jext(iY) = (Ji/Y(iY)**2.) + Jo
        ELSE
           Jext(iY) = Ji + Jo
        ENDIF
     END DO
  END IF

  res = 0.0d00
! this is the cut-off for printout of small values (in spectra)
  limval = 1.0d-20
! SmC(1..5, model) are found in Sub Analysis
  SmC(6,model) = Td(1,nY)
! SmC(7,model) and SmC(8,model) are Teff(L) and Teff(R), respectively,
! the effetive temperatures of the slab illumination sources, found from Fi in sub Analysis.
! The above-mentioned effective temperatures are not printed out anymore [MN].

! Zeljko's calculation of theta1, the ang. size (in arcsec) of the cavity for Fbol=1e-6 W/m2
  tht1 = 412.6d0/(dsqrt(Fi))
! colors
  call getfs(SpecChar(9,model),2,0,SC21)
  call getfs(SpecChar(10,model),2,0,SC31)
  call getfs(SpecChar(11,model),2,0,SC43)
! error in %
  if (SmC(5,model).lt.0.1d0) then
   call getfs(SmC(5,model)*100.0d0,0,0,Serr)
  else if (SmC(5,model).ge.0.1d0.and.SmC(5,model).lt.1.0d0) then
   call getfs(SmC(5,model)*100.0d0,0,1,Serr)
  else
   call getfs(SmC(5,model)*100.0d0,0,2,Serr)
  end if
! dust temperature at y=Y
  if (SmC(6,model).lt.99.5d0) then
   call getfs(SmC(6,model),0,0,STemp)
  else
   call getfs(SmC(6,model),0,1,STemp)
  end if
!--------------  overall parameters to *.out file -----------------------
! write header to output file *.out
  if (model.eq.1) then
   write(12,*)'         '
   write(12,*)' RESULTS:'
   write(12,*)' --------'
   if (slb) then
! ---------- slab output ----------------
    if(typentry(1).eq.5) then
     s1=' ###   Tau0    F_L(W/m2)   r1(cm)    Td(K)    Te_L(K)'
    else if (typentry(1).eq.3) then
     s1=' ###   Tau0     T1(K)    F_L(W/m2)   Td(K)    Te_L(K)'
    else
     s1=' ###   Tau0     T1(K)     r1(cm)     Td(K)    Te_L(K)'
    end if
    su1=' ###    1         2         3         4         5    '
    if(ksi.gt.0) then
     s2 ='   Te_R(K)  err'
     su2='     6       7 '
     write(12,'(a53,a15)') s1,s2
     write(12,'(a53,a15)') su1,su2
     su1=' ===================================================='
     su2='================'
     write(12,'(a53,a16)') su1,su2
    else
      s2='  err'
     su2='   6 '
     write(12,'(a53,a5)') s1,s2
     write(12,'(a53,a5)') su1,su2
     su1=' ===================================================='
     su2='====='
     write(12,'(a53,a6)')su1,su2
    end if
!------------- output for sphere -----------------------------
   elseif(sph) then
    if(typentry(1).eq.5) then
     s1= ' ###   tau0   Fi(W/m2)  r1(cm)    r1/rc   theta1   Td(Y)   err'
    else if (typentry(1).eq.2) then
     s1= ' ###   tau0     T1(K)   Fi(W/m2)  r1/rc   theta1   Td(Y)   err'
    else
     s1= ' ###   tau0     T1(K)    r1(cm)   r1/rc   theta1   Td(Y)   err'
    end if
    su1= ' ###     1        2        3        4        5       6      7 '
    if(rdwa.or.rdw) then
      s2='   Mdot      Ve       M> '
     su2='     8        9       10 '
     write(12,'(a62,a25)')s1,s2
     write(12,'(a62,a25)')su1,su2
     su1=' =============================================================='
     su2='============================'
     write(12,'(a63,a28)')su1,su2
! ** private rdw file **
     if (rdwpr) then
      s1= '###   tau0      tauF     Mdot      Ve       M>       '
      su1='###    1          2        3        4       5       6'
      s2= 'Ginf/G1   P    delta  d/sqrt(w1)  winf     Phi    zeta(1)'
      su2='        7        8        9        10       11       12'
      write(66,'(a53,a57)')s1,s2
      write(66,'(a53,a55)')su1,su2
     end if
    else
     write(12,'(a63)')s1
     write(12,'(a63)')su1
     su1=' ============================================================='
     write(12,'(a63)')su1
    end if
   end if
  end if
! print output tables for ea.model
  if(typentry(1).eq.5) then
! if T1 in input
   var1 = Fi
   var2 = Cr1
  else if (typentry(1).eq.2) then
! if L,r1 in input
   var1 = Td(1,1)
   var2 = Fi
  else
! if Fi or Tei in input
   var1 = Td(1,1)
   var2 = Cr1
  end if
!---------------- Output for slab: ---------------------------
  if(slb) then
   if (ksi.gt.0) then
! if second source on the right
     write(12,'(i4,1p,6e10.2,a3)') &
         model, taufid, var1, var2, Td(1,nY), SmC(7,model),SmC(8,model), Serr
   else
    write(12,'(i4,1p,5e10.2,a3)')  &
         model, taufid, var1, var2, Td(1,nY), SmC(7,model), Serr
   end if
!---------- for spherical shell ------------------------
  elseif(sph) then
   if (rdwa.or.rdw.or.rdwpr) then
    write(12,'(i4,1p,6e9.2,a1,a3,1p,3e9.2)') &
         model, taufid, var1, var2, r1rs, tht1,Td(1,nY),' ',Serr, CMdot, CVe, CM
   else
    write(12,'(i4,1p,6e9.2,a1,a3)') &
         model, taufid, var1, var2, r1rs, tht1,Td(1,nY),' ',Serr
   end if
!!!  if (left.eq.1) then
! If rc/r1 > 0.1% issue a warning about violation of point source assumption
!!!  I don't see the need for this. We have such a condition using Teff (see sub OccultMSG). [MN,2009]
!!!  This appears in almost all spherical runs.
!!!    if (r1rs.lt.1.0d+3) then
!!!     write(12,'(a86)') 'The user specified input radiation violates point source assumption (rc/r1<0.1%)'
!!!     write(12,'(a40)') 'Treat all solutions carefully'
!!!    end if
!!!   end if

   if (rdwpr) then
!** private rdw file **
    if (model.eq.1) then
     write(66,'(a11,1p,e9.3,a10,1p,e9.3,a12,1p,e9.3,a13,e9.3)')  &
          '###   qv = ',qv,', Qstar = ',q_star,', v1/vinf = ',pow, &
          ', (g/r)max = ',ptr(1)
    end if
    dmax = dsqrt(pow*winf)
    if (G1.gt.0) then
     GinfG1 =  Ginf / G1
    else
     GinfG1 = 0
    end if
    write(66,'(i4,1p,5e9.2,7e9.2)') &
         model, taufid, tauF(nY), CMdot, CVe, CM, &
         GinfG1, Prdw, delta, delta/dmax, winf, Phi, zeta1
   end if

   if(right.eq.0) then
    if (startyp(1).eq.1.or.startyp(1).eq.2) then
     if(Tstar(1).lt.Te_min) then
      call getfs(Tstar(1),0,1,Tstr)
      write(12,'(a50,a5,a5)') &
           ' ** WARNING: the input spectrum is a black-body at ',Tstr,' K **'
      call getfs(Te_min,0,1,Tstr)
      write(12,'(a50,a5,a5)') &
           ' *the point-source assumption requires min Teff of ',Tstr,' K **'
     end if
    end if
   end if
! end if for geometry
  end if

  call getfs(SpecChar(1,model),2,0,SB98)
  call getfs(SpecChar(2,model),2,0,SB11)
  call getfs(SpecChar(4,model),2,0,Sbet1)
  call getfs(SpecChar(5,model),2,0,Sbet2)
!--------------     spectral properties to *.spp file   -----------------------
  if (iSPP.ne.0) then
   write(19,'(i3,1p,5e9.2,7a6,e9.2)') model, taufid, SmC(1,model), &
        (SpecChar(i,model),i=6,8), SC21, SC31, SC43, Sbet1, Sbet2, &
        SB98, SB11, SpecChar(3,model)
  end if

!--------------   spectrum to *.s##  file   ------------------------
  if (iA.ne.0) then
   if(slb) then
    hdsp1 = '#  lambda     fRight       xAtt       xDs   '
   else
    hdsp1 = '#  lambda     fTot        xAtt       xDs   '
   end if
   hdsp2 = '    xDe       fInp      TauTot     albedo  '
   unt = 15
   call line(1,2,unt)
   if(slb) then
    write(unt,'(a7,i3,a8,f8.3,a36)')  &
         '# model',model,' taufid=',taufid,'   spectrum from the right slab side '
   else
    write(unt,'(a7,i3,a8,f8.3,a12)') '# model',model,' taufid=',taufid,'   spectrum '
   end if
   call line(1,1,unt)
   do iL = 1, nL
    if(slb) then
! the right-side spectra for slab: fsL + fds + fde
!     ftot(iL,nY) = ftot(iL,nY) + ksi*fsR(iL,nY)
! FH 10/18/10
     ftot(iL,nY) = fsL(iL,nY) + fde(iL,nY) + fds(iL,nY) + ksi*fsR(iL,nY)
!!**     if(ftot(iL,nY).LE.0.0) ftot(iL,nY) = dabs(ftot(iL,nY))
!!**    flux can be negative if flowing "right-to-left"
    else
! only the outward propagating flux
! fsL - for the central source; fsRp - for the external rad.
! only the shell sp.shape, external source not added
     ftot(iL,nY) =  fsL(iL,nY) + fde(iL,nY) + fds(iL,nY)
     if (dabs(ftot(iL,nY)).LT.limval) ftot(iL,nY) = 0.0d0
    end if
    faux(iL) = ftot(iL,nY)/lambda(iL)
   end do
   call Simpson(npL,1,nL,lambda,faux,res)
! normalization factor for output spectra
   fnorm = res
   call getOmega(nG,omega)
   do iL = 1, nL
    if (ftot(iL,nY).ne.0.0d0) then
     xs = fsL(iL,nY)/ftot(iL,nY)
     xds = fds(iL,nY)/ftot(iL,nY)
     xde = fde(iL,nY)/ftot(iL,nY)
    else
     xs = 0.0d0
     xds = 0.0d0
     xde = 0.0d0
    end if
!  no need to print negligible values
    if (xs.lt.limval) xs = 0.0d0
    if (xds.lt.limval) xds = 0.0d0
    if (xde.lt.limval) xde = 0.0d0
    if (fsL(iL,1).lt.limval) fsL(iL,1) = 0.0d0
! rescale ftot with the bolom flux
    ftot(iL,nY) = ftot(iL,nY)/fnorm
! if flag-selected, print ftot in [W/m2]
    IF (iPhys.eq.1) ftot(iL,nY) = ftot(iL,nY)*Jext(nY)*4*pi
    Elems(iL,1) = lambda(iL)
    Elems(iL,2) = ftot(iL,nY)
    Elems(iL,3) = xs
    Elems(iL,4) = xds
    Elems(iL,5) = xde
    if (slb) then
      Elems(iL,6) = fsL(iL,1)/fsLbol(1)
    else
! the input radiation is from central source, at y=1
! and/or from external, at y=Yout
! Elems(iL,6) = (fsL(iL,1) + fsrm(iL,nY))/fnorm
     Elems(iL,6) = fsL(iL,1) /fnorm
    end if
    Elems(iL,7) = tautot(iL)
    Elems(iL,8) = omega(1,iL)
   end do
!------ tabulate the spectra in the desired form ----------
   write(unt,'(2(a45))') hdsp1,hdsp2
   call maketable(Elems,npL,8,unt)
! spectra from the illuminated slab side (file *.z##)
   if (slb) then
    do iL = 1, nL
! the left-side spectra are: |R*fsR + fm|=|fsL-ftot|
!     ftot(iL,1) = dabs(fsL(iL,1) - ftot(iL,1))
! FH 10/19/10
     ftot(iL,1) =  dabs(fds(iL,1) - fde(iL,nY) - fds(iL,nY)) 
     faux(iL) = ftot(iL,1)/lambda(iL)
    end do
    call Simpson(npL,1,nL,lambda,faux,res)
    fnorm = res
    call getOmega(nG,omega)
    do iL = 1, nL
     if (ftot(iL,1).ne.0.0d0) then
      xs =  fsR(iL,1)/ftot(iL,1)
      xds = dabs(fds(iL,1)/ftot(iL,1))
      xde = dabs(fde(iL,1)/ftot(iL,1))
     else
      xs = 0.0d0
      xds = 0.0d0
      xde = 0.0d0
     end if
     if (xs.lt.limval) xs =0.0d0
     if (xds.lt.limval) xds =0.0d0
     if (xde.lt.limval) xde =0.0d0
     if(fsR(iL,nY).lt.limval) fsR(iL,nY) = 0.0d0
!     if(ftot(iL,1).lt.limval) ftot(iL,1) = 0.0d0
! rescale ftot with the bolom flux for z-spectra
     ftot(iL,1) = dabs(ftot(iL,1))/fnorm
!FH 10/18/10
     IF (iPhys.eq.1) ftot(iL,1) = ftot(iL,1)*Jext(1)*4*pi
     Elems(iL,1) = lambda(iL)
     Elems(iL,2) = ftot(iL,1)
     Elems(iL,3) = xs
     Elems(iL,4) = xds
     Elems(iL,5) = xde
     if (ksi.gt.0) then
      Elems(iL,6) = fsR(iL,nY)/fsRbol(1)
     else
      Elems(iL,6) = 0.0d0
     end if
     Elems(iL,7) = tautot(iL)
     Elems(iL,8) = omega(1,iL)
    end do
    if (iA.eq.3) unt=25
! append to the .s## file or write in a separate .z## file (if iA=3)
    call line(1,1,unt)
    write(unt,'(a7,i3,a8,f8.3,a36)') '# model',model,  &
         ' taufid=',taufid,'   spectrum from the left slab side'
    call line(1,1,unt)
    hdsp1 = '#  lambda      fLeft       xAtt       xDs   '
    write(unt,'(2(a45))') hdsp1,hdsp2
    call maketable(Elems,nL,8,unt)
   end if
  end if

! spectral properties for *.zpp file in slab case
  if (slb.and.iSPP.ne.0) then
   call getfs(SpecChar(12,model),2,0,SB98)
   call getfs(SpecChar(13,model),2,0,SB11)
   call getfs(SpecChar(15,model),2,0,Sbet1)
   call getfs(SpecChar(16,model),2,0,Sbet2)
   call getfs(SpecChar(20,model),2,0,SC21)
   call getfs(SpecChar(21,model),2,0,SC31)
   call getfs(SpecChar(22,model),2,0,SC43)

   if (iSPP.eq.3) then
! write spectral properties to *.zpp file
    write(24,'(i3,1p,5e9.2,7a6,e9.2)')  &
         model, taufid, SmC(1,model), (SpecChar(i,model),i=17,19), &
         SC21,SC31,SC43, Sbet1,Sbet2, SB98,SB11,SpecChar(14,model)
   else
! write into string zline. This will be appended to spp file in Cllose
    write(zline(model),'(i3,1p,5e9.2,7a6,e9.2)')  &
         model, taufid, SmC(1,model), (SpecChar(i,model),i=17,19), &
         SC21,SC31,SC43, Sbet1,Sbet2,SB98,SB11,SpecChar(14,model)
   end if
  end if

!-----------  radial quantities to *.r## (old *.bxx) file -------------
  if (iB.ne.0) then
   deallocate(Elems)
   allocate(Elems(nY,9))
!   hdrslb1= '#     t        tauF      epsilon       Td '
!   hdrsph1= '#     y         eta         t         tauF '
!   hdrsph2= '     epsilon      Td         rg '

   hdrslb1= '#     t        Td      epsilon       tauF '
   hdrslb2= '      fsbol      fRbol      fLbol '
   hdrsph1= '#     y         Td         eta         t '
   hdrsph2= '     tauF      epsilon        rg '
   hdrdyn= '         u        drift'
   unt = 16
   call line(1,2,unt)
   if(slb) then
    write(unt,'(a7,i3,a8,f8.3,a18)') '# model',model,' taufid=',taufid,'  spatial profiles'
   else
    write(unt,'(a7,i3,a8,f8.3,a18)') '# model',model,' taufid=',taufid,'  radial profiles '
   end if
   call line(1,1,unt)
!--------- for slab ---------
   if (slb) then
    do iY = 1, nY
     Elems(iY,1) = tr(iY)
     Elems(iY,2) = Td(1,iY)
!     Elems(iY,2) = tauF(iY) !rearranged Td as second column for plotting convenience [MN]
     Elems(iY,3) = eps(iY)
     Elems(iY,4) = tauF(iY)
     Elems(iY,5) = fsbol(iY)/fsbol(1)
     Elems(iY,6) = fpbol(iY)/fsbol(1)
     Elems(iY,7) = fmbol(iY)/fsbol(1)
    end do
    write(unt,'(a42,a34)') hdrslb1,hdrslb2
    call maketable(Elems,nY,7,unt)
!------  for spherical shell --------
   elseif(sph) then
    do iY = 1, nY
     Elems(iY,1) = Y(iY)
     Elems(iY,2) = Td(1,iY) !rearranged Td as 2nd column [MN]
     Elems(iY,3) = eta(Y(iY))
     Elems(iY,4) = tr(iY)
     Elems(iY,5) = tauF(iY)
     Elems(iY,6) = eps(iY)
     Elems(iY,7) = rg(1,iY)
     if (rdwpr) then
! redefine for private rdw (denstyp.eq.6) option
      Elems(iY,5) = gamma(iY)
      Elems(iY,7) = qF(iY)
     end if
    end do
! check values:
    do i = 1, 7
     do iY = 1, nY
      if(Elems(iY,i).lt.limval) Elems(iY,i) = 0.0d0
     end do
    end do
! with dynamics
    if (rdw) then
     do iY = 1, nY
      Elems(iY,8) = ugas(iY)/ugas(nY)
      Elems(iY,9) = vrat(1,iY)
     end do
! check values:
     do i = 8, 9
      do iY = 1, nY
       if(Elems(iY,i).lt.limval) Elems(iY,i) = 0.0d0
      end do
     end do
     write(unt,'(a42,a32,a23)') hdrsph1,hdrsph2,hdrdyn
     call maketable(Elems,nY,9,unt)
    else
     write(unt,'(a42,a32)') hdrsph1,hdrsph2
     call maketable(Elems,nY,7,unt)
    end if
! end if for geometry
   end if
! end if for the iB (radial) flag
  end if

!--------------   intensities to *.inn (old *.cxx) file  --------------
  if (abs(ic).ne.0) then
! slab intensity (found at the end of subroutine slbradt)
! theta(nmu) are the angles of output intensities
   if (slb) then
    deallocate(Elems)
    allocate(Elems(npL,nmu+2))

    hdint = '   lambda'
    unt = 17
    call line(1,2,unt)
    write(unt,'(a7,i3,a8,f8.3,a32)')'# model',model,' taufid=', &
         taufid,' transmitted i(theta)*cos(theta)'
    call line(1,1,unt)
    do iL = 1, nL
     Elems(iL,1) = lambda(iL)
     do imu = 1, nmu
      if(iPhys.eq.1) SLBintm(imu,iL) = SLBintm(imu,iL)*Jext(nY)
      Elems(iL,imu+1) = SLBintm(imu,iL)
     end do
     Elems(iL,nmu+2) = istR(iL)
    end do
! write(unt,'(a9,21f11.3)')hdint,(theta(imu),imu=1,nmu)
! printout angles in degrees
! write(unt,'(a9,37f11.1,a9)') hdint,
! &                    (theta(imu)*180.0d0/pi,imu=1,nmu),'     IstR'
    write(unt,'(a9,100f11.1)') hdint,(theta(imu)*180.0d0/pi,imu=1,nmu)
    call maketable(Elems,npL,nmu+1,unt)
! adding the column with stellar ints at the end of the table
!        call maketable(Elems,nL,nmu+2,unt)

    hdint = '   lambda'
    unt = 17
    call line(1,2,unt)
    write(unt,'(a7,i3,a8,f8.3,a32)')'# model',model, &
         ' taufid=',taufid,' reflected cos(theta)*i(theta)'
    call line(1,1,unt)
    do iL = 1, nL
     Elems(iL,1) = lambda(iL)
     do imu = 1, nmu
      if(iPhys.eq.1) SLBintp(imu,iL) = SLBintp(imu,iL)*Jext(1)
      Elems(iL,imu+1) = SLBintp(imu,iL)
     end do
    end do
!       write(unt,'(a9,21f11.3)')hdint,(theta(imu),imu=1,nmu)
!       printout angles in degrees
    write(unt,'(a9,99f11.1)')hdint,(theta(imu)*180.0d0/pi,imu=1,nmu)
    call maketable(Elems,npL,nmu+1,unt)
!------  for spherical shell --------
   elseif(sph) then
    deallocate(Elems)
    allocate(Elems(np+2,nLambdaOut+2))

    hdint = '#     b          t(b)'
    hdcon = '#   offset '
    hdvis = '#     q    '
    unt = 17
    call line(1,2,unt)
    write(unt,'(a7,i3,a8,f8.3,a14)') '# model',model,' taufid=',taufid,'   raw image  '
    call line(1,1,unt)
    do i = 1, nP+2
     Elems(i,1) = bOut(i)
     Elems(i,2) = tauZout(i)
     do j = 1, nLambdaOut
! check values:
      if(IntOut(j,i).ne.IntOut(j,i).or.IntOut(j,i).lt.limval) then
       IntOut(j,i) = 0.0d0
      end if
      Elems(i,j+2) = IntOut(j,i)
! we want intensity in Jy/arcsec^2
! this was the bug in intensity output for sphere,
! the missing 4piY^2 factor for intensity output [June 2006]
! Elems(i,j+2) = 7.83 * LambdaOut(j) * Fi * Elems(i,j+2)
      IF (iPhys.eq.1) THEN
        Elems(i,j+2) = 7.834d0*LambdaOut(j)*(Jext(nY)*4.0d0*pi*Yout**2.0d0)*Elems(i,j+2)
      ELSE
        Elems(i,j+2) = 7.834d0*LambdaOut(j)*(4.0d0*pi*Yout**2.0d0)*Elems(i,j+2)
      END IF
     end do
    end do
    write(unt,'(a21,20f11.2)')hdint,(LambdaOut(j),j=1,nLambdaOut)
    call maketable(Elems,np+2,nLambdaOut+2,unt)
   end if
  end if
  if (iC.lt.0) then
!---------  convolved images either add to .i## file or write in *.c## file --
   if(iC.eq.-3) unt = 21
   call line(1,2,unt)
   write(unt,'(a7,i3,a8,f8.3,a20)') '# model',model,' taufid=',taufid,'   convolved image  '
   call line(1,1,unt)

   if(allocated(Elems)) deallocate(Elems)
   allocate(Elems(nconv,nLambdaOut+1))

   do i = 1, nconv
    Elems(i,1) = offset(i)
    do j = 1, nLambdaOut
     if(convint(j,i).lt.limval) convint(j,i) = 0.0d0
     Elems(i,j+1) = convint(j,i)
    end do
   end do
   write(unt,'(a11,20f11.2)')hdcon,(LambdaOut(i),i=1,nLambdaOut)
   call maketable(Elems,nconv,nLambdaOut+1,unt)
   if (psftype.lt.3.and.model.eq.1) then
! wavelength dependent psfs, print them separately in *.psf
! first generate wavelength dependent psfs
    do j = 1, nLambdaOut
     iLambda = j
! added dec.04 [mn]
     do i = 1, nconv
      psffunc(j,i) = psfn(offset(i))
! norm.needs to be done here again (after call to psfn)
      psffunc(j,i) = psffunc(j,i)/psfarea(j)
! check dynamic range
      call chkrange(dynrange,psffunc(j,i))
     end do
    end do
! print them out
    do i = 1, nconv
     write(23,'(1p,e12.5,20e10.3)')offset(i),(psffunc(j,i),j=1,nLambdaOut)
    end do
   end if
  end if
!--------------  visibility curves to *.vnn file    ------------------------

  if (sph) then
   if (iV.ne.0) then
    if(abs(iC).eq.3) unt = 22
    call line(1,2,unt)
    write(unt,'(a7,i3,a8,f8.3,a14)') '# model',model, &
         ' taufid=',taufid,'  visibility  '
    call line(1,1,unt)

    if(allocated(Elems)) deallocate(Elems)
    allocate(Elems(nvisi,nLambdaOut+1))

    do i = 1, nvisi
     Elems(i,1) = qtheta1(i)
     do j = 1, nLambdaOut
      if(visib(j,i).lt.limval) visib(j,i) = 0.0d0
      Elems(i,j+1) = visib(j,i)
     end do
    end do
    write(unt,'(a11,20f11.2)')hdvis,(LambdaOut(i),i=1,nLambdaOut)
    call maketable(Elems,nvisi,nLambdaOut+1,unt)
   end if
  endif

 deallocate(Elems)
!-----------------------------------------------------------------------

 return
end subroutine PrOut
!***********************************************************************

! ***********************************************************************
! This function is taken from Moshe Elitzur.           [Z.I., Nov. 1995]
! =======================================================================
double precision function RDINP(Equal,iUnit)
! =======================================================================
!     Read lines, up to 232 long, from pre-opened unit iUnit and extract
!     all input numbers from them. When EQUAL is set, numeric input data
!     must be preceded by an equal sign.All non-numeric data and numbers
!      not preceded by = when EQUAL is on are ignored.RDINP = next number
!      encountered (after equal sign) and terminated by a nun-numeric
!      symbol (termination with blank is best). Commas and exponential
!      notation are allowed.  All text after % is ignored, as in TeX.
!      Lines with * in the first column are echoed to the output device.
!      The number is comprised of an actual VALUE, decimal part FRAC and
!      (integer) exponent PWR.  It has a SIGN, and the exponential power
!      has sign SIGNEX. Logical flag to the decimal part is set in
!      DECIMAL. The search is conducted between FIRST, which is
!     continuously increased, and LAST.  A new line is read when FIRST
!     exceeds LAST, and can be forced by calling with -iUnit.  Actual
!     extraction of numerical value is done in separate FUNCTION VAL.
! =======================================================================
      IMPLICIT None
      Integer iUnit,ind,First, Last
      DOUBLE PRECISION Value,Val,Frac,Pwr,Sign,Signex
      CHARACTER Card*(232),CR,prev,Term,Next
      Logical Equal,digit,minus,sgn,dot,E,decimal
      Save Card,First,Last
      DATA First/1/, Last/0/
! -----------------------------------------------------------------------
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
!     Start the search for the next number:
  1   RDINP  = 0.
      VALUE  = 0.
      FRAC   = 0.
      PWR    = 0.
      SIGN   = 1.
      SIGNEX = 1.
      Decimal = .False.
      If (first.gt.last) then
!     Time to get a new line
         READ (iUnit, '(A)' , END = 99) Card
         first = 1
         last = len(Card)
!        Find start of trailing junk:
         DO WHILE (Card(last:last).LE.' ')
          last = last - 1
          if (last.lt.first) goto 1
         END DO
         IF (Card(first:first).EQ.'*') WRITE (12,'(A)') Card(1:last)
         ind = Index(Card,'%')
         if (ind.gt.0) last = ind - 1
      End If

!     Get past the next '=' when the EQUAL flag is set
      If (Equal) then
        DO WHILE (Card(first:first).NE.'=')
          first = first + 1
          IF (first.GT.last) goto 1
        END DO
      End If
!     OK, start searching for the next digit:
!      Do While (.not.digit(Card(first:first))) - corrected by ME
       Do While (.not.(Card(first:first).ge.'0'.AND.Card(first:first).le.'9'))
          first = first + 1
          if (first.gt.last) goto 1
       End Do
!     Check if it is a negative or decimal number
      If (first.gt.1) then
         prev = Card(first-1:first-1)
!         necessary syntax changes, otherwise the code crushes on ncx[MN]
!         if (minus(prev)) sign = -1.
         if (prev.eq.'-') sign = -1.
!         if (dot(prev)) then
         if (prev.eq.'.') then
           decimal = .True.
           if (first.gt.2 .and.(Card(first-2:first-2).eq.'-')) sign = -1.
         end if
      End If
!     Extract the numerical value
      IF (.not.Decimal) Then
         Value = VAL(Card,first,last,decimal,Term)
!        Check for a possible decimal part.  Termination with '.E' or
!        '.e' is equivalent to the same ending without '.'
!         if (first.lt.last.and.dot(Term)) then
         if (first.lt.last.and.Term.eq.'.') then
            first = first + 1
            next = Card(first:first)
!            if (digit(next)) decimal = .true.
            if (next.GE.'0' .AND. next.LE.'9') decimal = .true.
!            if (E(next)) Term = 'E'
            if (next.eq.'E' .OR. next.EQ.'e') Term = 'E'
         end if
      END IF
!     Extract the decimal fraction, when it exists
      IF (Decimal) Frac = VAL(Card,first,last,decimal,Term)
!     An exponential may exist if any part terminated with 'E' or 'e'
!       IF (first.lt.last.and.E(term)) then
       IF (first.lt.last.AND.(term.eq.'E'.or.term.eq.'e')) then
         first = first + 1
         next = Card(first:first)
!         if (first.lt.last.and.sgn(next))then
         if (first.lt.last.AND.(next.EQ.'+'.or.next.EQ.'-')) then
            first = first + 1
!            if (minus(next)) Signex = -1.
            if (next.eq.'-') Signex = -1.
         end if
         decimal = .False.
         PWR = VAL(Card,first,last,decimal,Term)
      END IF
!     Finally, put the number together
      RDINP = Sign*(Value + Frac)*10**(Signex*PWR)
      Return

99    WRITE (12,'(3(1x,a,/))')                                           &
      ' ****TERMINATED. EOF reached by RDINP while looking for input. ', &
      ' *** Last line read:',Card
! -----------------------------------------------------------------------
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
! This function is taken from Moshe Elitzur            [Z.I., Nov. 1995]
!=======================================================================
double precision function Val(Card,first,last,decimal,term)
!=======================================================================
!  Extract numerical value from CARD, begining at position FIRST up
!  to the first non-digit encountered.  The terminating character is
!  returned in TERM, its position in FIRST. Commas are treated as
!  simple separators.
!=======================================================================
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
    do iG = 1, nG
     if (SLB) then
      write(12,'(a45,1p,e9.2,a2)')  &
           ' Dust temperature on the slab left boundary:', Tsub(iG),' K'
     else
      write(12,'(a41,1p,e9.2,a2)')  &
           ' Dust temperature on the inner boundary:', Tsub(iG),' K'
     end if
    end do
   else if (typEntry(1).eq.1) then
    if (slb) then
     if (startyp(1).gt.3) then
      write(12,'(a33,1p,e9.2,a6)') &
           ' Flux at the slab left boundary:', dilutn*var1,' W/m^2'
     else
      write(12,'(a33,1p,e9.2,a6)') &
           ' Flux at the slab left boundary:', var1,' W/m^2'
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


!!=======================================================================
! Routines and funcitons related to RDW (Radiatively driven winds).
!                                                    [MN, Aug,2010]
!!=======================================================================

!***********************************************************************
subroutine Winds(nG,EtaOK)
!=======================================================================
! This subroutine takes care of the interface between radiatively driven
! winds and radiative transfer.  It is entered
! after a radiative transfer calculation with given eta.
! This sub caclulates the reddening profile phi and passes it to the dynamics
! module, which returns the velocity and density profiles corresponding
! to the given phi. Convergence is achived when the eta returned from
! the dynamics calculation is the same as that used to produce phi.
!
! Notations follow EI01 (MNRAS 327, 403)
!=======================================================================

  use common
  implicit none

  integer nG, EtaOK, iY, iL, err
  double precision etaold(npY),acceta, resaux,   &
       faux(npL), Qstar, Qfid, phi1, localp, eps_loc, gammamax, &
       wscale, uacc, reddn(npY), w(npY), gmax
! the parameter ver determines the version of the velocity formal solution
! 1 for linear, 2 quadratic.  it is specified in input and carried in /dyn/
! the quadratic solution is from equation d1 in ei01. the linear is
! obtained similarly from the differential equation 24 by using dw^2 = 2wdw
! and dividing through by 2w
! -----------------------------------------------------------------------

  if (iX.ge.1) then
   write(18,*)' doing dynamics'
   if (ver.eq.1) then
    write(18,*)' linear version of velocity formal solution'
   elseif (ver.eq.2) then
    write(18,*)' quadratic version of velocity formal solution'
   else
    write(12,*)' **************************** '
    write(12,'(a,i3)')'  illegal input ver = ', ver
    write(12,*)'     ver must be 1 or 2!      '
    write(12,*)'       program stopped        '
    write(12,*)' **************************** '
    stop
   end if
  end if
  if(iVerb.eq.2) write(*,*)' doing dynamics'
! so far it works for nG=1 only:
  if (nG.gt.1) then
   write(12,*)' **************************** '
   write(12,*)' change dynamics sub to nG>1! '
   write(12,*)'       program stopped        '
   write(12,*)' **************************** '
   stop
  end if
!-----------------------------------------------------------------------
! assign input parameters to local variables
  gammamax = ptr(1)
  eps_loc = pow
! accuracy for velocity convergence same as for utot:
  uacc = accConv
! extinction efficiency at the fiducial wavelength
  Qfid = sigexfid/aveA
! calculate Qstar and the scale factor of w:
  do iL = 1, nL
   faux(iL) = (sigmaA(1,iL)+sigmaS(1,iL))*ftot(iL,1)/lambda(iL)
  end do
  call Simpson(npL,1,nL,lambda,faux,resaux)

! Qstar is from EI'01 equation 4, wscale is from equation 29
  Qstar = resaux / aveA
  wscale = taufid/Qfid
!
!-----------------------------------------------------------------------
! here's the eta that was used in the radiative transfer
  do iY = 1, nY
   etaold(iY) = etadiscr(iY)
  end do
! and here's the resulting reddening profile
  do iY = 1, nY
   do iL = 1, nL
    faux(iL) = (sigmaA(1,iL)+sigmaS(1,iL))*ftot(iL,iY)/lambda(iL)
   end do
   call Simpson(npL,1,nL,lambda,faux,resaux)
   if (iY.eq.1)    phi1 = resaux
   reddn(iY) = resaux / phi1
  end do
! now find new eta
  err = 0
  call dynamics(eps_loc, gammamax, wscale, reddn, w, vrat, gmax, &
       localp, etadiscr, y, nY, nG, uacc, err, iX, ver)
! and check convergence (ptr(2) is specified in input)
  acceta = ptr(2) * accuracy
  call chkconv(acceta,etaold,etadiscr,EtaOK)
  if (iX.ge.1) then
   write(18,'(2(a,1pe10.3))')'  p = ', localp, '  gmax = ', gmax
   write(18,*) '     y    ugas(new)   tauf      etaold    etanew    ratio'
   do iY = 1, nY
!**************************************
! for output compatibiLity; we can do away with qF and tauF, which
! have only nostalgic reasons.  EI01 never uses them
    ugas(iY) = Qstar*w(iY)
    qF(iY)   = (Qstar/Qfid)*reddn(iY)
    faux(iY) = qF(iY)*etadiscr(iY)
    call Simpson(npY,1,iY,y,faux,resaux)
    taufdyn(iY) = taufid*resaux
!**************************************
    acceta = etaold(iY) / etadiscr(iY)
    write(18,'(1p,6e10.3)')Y(iY), ugas(iY), taufdyn(iY),etaold(iY), etadiscr(iY),acceta
   end do
   if (EtaOK.eq.1) then
    write(18,*)' Convergence on eta achieved.'
   else
    write(18,*)' Convergence on eta not achieved.'
    write(18,*)' Going to the next iteration.'
   end if
  end if
! save y to Yprev and nY to nYprev
  do iY = 1, nY
   Yprev(iY) = Y(iY)
  end do
  nYprev = nY
! -----------------------------------------------------------------------
  return
end subroutine Winds
! ***********************************************************************


!***********************************************************************
subroutine Dynamics(eps_loc, f, uscale, phi_loc, u, zeta, gmax, &
     p, eta, y, nY, nG, acc, err, ipr, ver)
!=======================================================================
! Calculates the velocity structure of a radiatively driven
! wind given the reddening profile phi.  It returns the
! profiles w, eta and zeta, and the wind parameters P and gmax.
!
! This subroutine calculates the initial guess and controlls the
! convergence of the iterative procedure, with CalcVel actually
! calculating a new velocity profile from the previous one.
! Convergence is checked both for u(y) (= w(y)) and its derivative eta(y)
! because eta is crucial for the radiative transfer.
! *** This version works for single size grains only ***
!
! Implementing equations from EI01 (MNRAS 327, 403)
!=======================================================================

  implicit none
  integer npY, npP, npX, npL, npG, npR
  include 'userpar.inc'
  parameter (npG=1)
  integer nG, nY, iY, itr, etaconv, uconv, err, ipr, itmax, ver
  double precision eta(npY), etaold(npY), u(npY), uold(npY),  &
       phi_loc(npY), zeta(npY), Y(npY), eps_loc, f, acc, gmax, &
       uscale, n, p, wf, k, e1
  data   itmax/100/, k/0.4/
!  we may wish to control itmax and k as input parameters
!-----------------------------------------------------------------------
! initial approximation for u(y) from EI01, eq C6
! wf is from eq. 29 with epsilon correction (eq. C8)
  wf = (1.0d+00/(1.0d+00 - eps_loc))*phi_loc(nY)*uscale
! add a correction for the finite outer radius so that wf = u(nY):
  e1 = 1.0d+00 - eps_loc**(1.0d+00/k)
  wf = wf/(1.0d+00 - e1/Y(nY))**k
! and now calculate all u from eq. C6
  do iY = 1, nY
   uold(iY) = wf*(1.0d+00 - e1/Y(iY))**k
! initial eta is irrelevant; might as well use
! the one passed from radiative transfer:
   etaold(iY) = eta(iY)
  end do
! iterations until u and eta converge within acc
  do itr = 1, itmax
   call calcvel(eps_loc,f,uscale,phi_loc,uold,u,zeta,gmax,y,nY,ver)
   call calceta(u, zeta, eta, n, y, nY)
   p = dsqrt(uscale/n)  ! eq.(46) in IE'01

! check convergence of u and eta
   call chkconv(acc,uold,u,uconv)
   call chkconv(acc,etaold,eta,etaconv)
! convergence required for both u(y) and eta(y)
   err = 1 - etaconv * uconv
   if (err.ne.0) then
! did not converge, repeat the exercise...
    do iY =1, nY
     uold(iY) = u(iY)
     etaold(iY) = eta(iY)
    end do
   else
! we're done:
    if (ipr.ge.1) write(18,'(a35,i3)')' Number of iterations to converge:',itr
    return
   end if
  end do
! -----------------------------------------------------------------------
  return
end subroutine Dynamics
!***********************************************************************


! ***********************************************************************
      SUBROUTINE CalcVel(eps_loc, f, ws, phi_loc, wold, w, zeta, gmax,Y, nY, ver)
! =======================================================================
! Calculates the scaled gas velocity w(y) from wold(y), the previous
! velocity profile, and the given reddening profile phi
! The calculation follows the formal solution in Appendix D
! ver = 1 triggers the linear version, 2 the quadratic (eq. D1)
! All symbols are as defined there, ws = tauV/QV and P2 = P^2
! =======================================================================
      IMPLICIT none
      INTEGER npY, npP, npX, npL, npG, npR
      INCLUDE 'userpar.inc'
      PARAMETER (npG=1)
      INTEGER ver, nY, iY
      DOUBLE PRECISION eps_loc, f, wold(npY), w(npY), phi_loc(npY), &
           zeta(npG,npY), Y(npY), z(npY), zz(npY), gmax, g,        &
           ws, N, P2, ww1, aux, F1(npY), F2(npY)
! -----------------------------------------------------------------------
!     first get the drift profile zeta
      CALL CalcDrift(phi_loc, wold, zeta, nY)
!     then the normalization N (= EtaINT; eq D3)
      DO iY = 1, nY
        F1(iY) = zeta(1,iY)/(wold(iY)*Y(iY)*Y(iY))
      END DO
      CALL SIMPSON(npY,1,nY,Y, F1, N)
!     and finally P (eq. D3):
      P2 = 2.0D+00*ws/N
!     Now the two versions diverge:
      IF (ver.eq.2) THEN
!        Quadratic version (eq. D1):
!         Get the profile z = integral(zeta*phi/y^2) and gmax
!         (= Gamma_min), the maximum gravitational correction
!         as defined in eqs D4 and D5:
          F1(1) = 0.0D+00
          z(1)  = 0.0D+00
          gmax  = 1.0D+00/zeta(1,1)
          DO iY = 2, nY
             F1(iY) = zeta(1,iY)*phi_loc(iY)/(Y(iY)*Y(iY))
             CALL SIMPSON(npY,1,iY,Y, F1, z(iY))
             aux = (1.0D+00 - 1.0D+00/Y(iY))/z(iY)
             if (aux.gt.gmax) gmax = aux
          END DO
          g = f/gmax
!         All ready.  Calculate the new w1 (ww1 = w1^2):
          aux  = eps_loc*eps_loc
          ww1  =(aux/(1.0D+00-aux))*P2*(z(nY)-g*(1.0D+00-1.0D+00/Y(nY)))
          w(1) = dsqrt(ww1)
!         and now the rest of the profile:
          DO iY = 2, nY
             aux   = ww1 + P2*(z(iY)-g*(1.0D+00-1.0D+00/Y(iY)))
             w(iY) = dsqrt(aux)
          END DO
      ELSE
!        Linear version.
!         obtained from the differential equation 24
!         by using dw^2 = 2wdw and dividing through by 2w
!         Now we need the profiles z = integral(zeta*phi/w*y^2),
!         zz = integral(1/w*y^2) and then gmax = max(zz/z)
          F1(1) = 0.0D+00
          F2(1) = 0.0D+00
          z(1)  = 0.0D+00
          zz(1) = 0.0D+00
          gmax  = 1./zeta(1,1)
          DO iY = 2, nY
             F2(iY) = 1.0D+00/(wold(iY)*Y(iY)*Y(iY))
             F1(iY) = F2(iY)*zeta(1,iY)*phi_loc(iY)
             CALL SIMPSON(npY,1,iY,Y, F1, z(iY))
             CALL SIMPSON(npY,1,iY,Y, F2, zz(iY))
             aux = zz(iY)/z(iY)
             if (aux.gt.gmax) gmax = aux
          END DO
          g = f/gmax
!         All ready.  Calculate the new w1:
          w(1)  = (eps_loc/(1.0D+00-eps_loc))*0.5D+00*P2*(z(nY) - g*zz(nY))
!         and now the rest of the profile:
          DO iY = 2, nY
             w(iY) = w(1) + 0.5D+00*P2*(z(iY) - g*zz(iY))
          END DO
      END IF
! -----------------------------------------------------------------------
      RETURN
      END
! ***********************************************************************


! ***********************************************************************
      SUBROUTINE CalcETA(w, zeta, Eta, EtaINT, Y, nY)
! =======================================================================
! Calculates the dimensionless density profile ETA(y) from EI eq. 25 given
! the velocity profile w(y) and its corresponding drift zeta(y)
! =======================================================================
      IMPLICIT none
      INTEGER npY, npP, npX, npL, npG, npR
      INCLUDE 'userpar.inc'
      PARAMETER (npG=1)
      INTEGER iY, nY
      DOUBLE PRECISION w(npY), Eta(npY), zeta(npG,npY), Y(npY),EtaINT
! ======================================================================
      DO iY = 1, nY
        Eta(iY) = zeta(1,iY)/(w(iY)*Y(iY)*Y(iY))
      END DO
!     now normalize eta:
      CALL SIMPSON(npY,1,nY,Y,Eta,EtaINT)
      DO iY = 1, nY
        Eta(iY) = Eta(iY)/EtaINT
      END DO
! -----------------------------------------------------------------------
      RETURN
      END
! ***********************************************************************

! ***********************************************************************
      SUBROUTINE CalcDrift(phi_loc,w,zeta,nY)
! =======================================================================
! Calculates the drift profile zeta from EI01 eq. 24
! without the correction for sub-sonic drift (theta = 0).
! =======================================================================
      IMPLICIT none
      INTEGER npY, npP, npX, npL, npG, npR
      INCLUDE 'userpar.inc'
      PARAMETER (npG=1)
      INTEGER nY, iY
      DOUBLE PRECISION phi_loc(npY), w(npY), zeta(npG,npY)
! -----------------------------------------------------------------------
      DO iY = 1, nY
         zeta(1,iY) = 1.0D+00 / (1.0D+00 + dsqrt(phi_loc(iY)/w(iY)))
      END DO
! -----------------------------------------------------------------------
      RETURN
      END
! ***********************************************************************




!!=====================================================================
! Below are AUXILIARY MATH and miscellaneous subroutines and functions
! arranged in alphabetical order.                       [MN, Aug,2010]
!!=====================================================================

!**********************************************************************
subroutine add(np1,nr1,np2,nr2,q1,q2,q3,qout)
!======================================================================
! This subroutine evaluates the following expression:
! [qOut] = [q1] + [q2] + [q3]. qout, q1, q2 and q2 are matrices of
! physical size (np2,np1) and real size (nr2,nr1).     [Z.I., Nov. 1995]
! ======================================================================
  implicit none
  integer npY, npP, npX, npL, npG, npR
  include 'userpar.inc'
  parameter (npG=1)

  integer  np1, nr1, np2, nr2, i2, i1
  double precision  q1(np2,np1), q2(np2,np1), q3(np2,np1),qout(np2,np1)
! ----------------------------------------------------------------------

! loop over index 2
  do i2 = 1, nr2
! loop over index 1
   do i1 = 1, nr1
    qout(i2,i1) = q1(i2,i1) +  q2(i2,i1) + q3(i2,i1)
   end do
  end do
! ----------------------------------------------------------------------
  return
end subroutine add
!**********************************************************************

!**********************************************************************
subroutine add2(flxs,flxe,fbsum,nY)
!======================================================================
! This subroutine is auxiliary for finding the bolometric
! diffuse flux.   [MN, May'99]
!======================================================================
  implicit none
  integer npY, npP, npX, npL, npG, npR
  include 'userpar.inc'
  parameter (npG=1)

  integer nY, iY
  double precision flxs(npL,npY), flxe(npL,npY), flxsb(npY),flxeb(npY), fbsum(npY)
!----------------------------------------------------------------------

  call bolom(flxs,flxsb)
  call bolom(flxe,flxeb)
  do iY = 1, nY
   fbsum(iY) = flxsb(iY) + flxeb(iY)
  end do
!----------------------------------------------------------------------
  return
end subroutine add2
!**********************************************************************

!**********************************************************************
SUBROUTINE ANALINT(Nanal,xaux,yaux,m,aux,error)
!======================================================================
! This subroutine calculates integral I(x**m*y(x)*dx). Both y and x are
! 1D arrays, y(i), x(i) with i=1,Nanal. The method used is approximation
! of y(x) by y = P(x) + d/sqrt(1-x*x), where P(x) is the polynomial of
! order Nanal-1, and analytic evaluation of the integral. It is assumed
! that xaux(1)=0. Coefficients are determined from the set of Nanal
! linear equations and subsequent call to the linear system solver
! LINSYS.                                              [Z.I., Nov. 1995]
! ANALINT is called from Nordlund to evaluate analytically the contribution
! of Nanal grid points. [MN]
! =======================================================================
     use common
      IMPLICIT none

      INTEGER i, j, Nanal, error
      DOUBLE PRECISION xaux(Nanal),yaux(Nanal),coeff(4),A(npY,npY),m,aux,b
! ----------------------------------------------------------------------
      error = 0
!     generate matrix A and vector B
      DO i = 1, Nanal
        DO j = 1, Nanal-1
          IF (xaux(i).EQ.0.0.AND.j.EQ.1) THEN
            A(i,j) = 1.0
          ELSE
            A(i,j) = xaux(i)**(1.0*j-1.0)
          END IF
        END DO
        A(i,Nanal) = 1.0/sqrt(1.0-xaux(i)*xaux(i))
      END DO

!     solve for the coefficients
      CALL LINSYS(Nanal,A,yaux,coeff,error)
        IF(error.NE.0) THEN
         CALL MSG(19)
         iERROR = iERROR + 1
         RETURN
        END IF
!     upper limit for integration:
      b = xaux(Nanal)
!     evaluate m-dependent contribution of the last term
      IF (m.GT.0.1) THEN
        IF (m.GT.1.1) THEN
!         this is for m=2
          aux = 0.5*(DASIN(b)-b*sqrt(1.-b*b))
        ELSE
!         this is for m=1
          aux = 1.0 - sqrt(1.-b*b)
        ENDIF
      ELSE
!       this is for m=0
        aux = DASIN(b)
      ENDIF
      aux = aux * coeff(Nanal)
!     add contribution from the polynom
      DO i = 1, Nanal-1
        aux = aux + coeff(i) * (b**(m+1.0*i)) / (m+1.0*i)
      END DO
! -----------------------------------------------------------------------
999   RETURN
END SUBROUTINE ANALINT
!***********************************************************************

!***********************************************************************
! This subroutine obtained from prof. P. Menguc, Dept. of Mechanical
! Engineering, University of Kentucky.                 [Z.I., Aug. 1996]
!-----------------------------------------------------------------------
!_______________________________________________________________________
!
! Subroutine bhmie calculates amplitude scattering matrix elements
! & efficiencies for extinction, total scattering and bacscattering,
! for a given size parameter and relative refractive index
!_______________________________________________________________________

subroutine bhmie (x,refrel,nang,s1,s2,qext,qsca,qback)
!=======================================================================

  implicit none

  double precision amu(100),theta(100),pi(100),tau(100)
  double precision pi0(100),pi1(100)
  integer n,nmx,nang,nn,nstop,j,jj
  complex d(3000),y,refrel,xi,xi0,xi1,an,bn,s1(200),s2(200)
  double precision xstop,ymod,psi0,psi1,psi,dn,dx,x,qext,qsca,qback
  double precision rn,t,p,fn,dang,chi,chi0,chi1,apsi,apsi0
  double precision apsi1
  dx=x
  y=x*refrel

!___________________________________________________________________
! series terminated after nstop terms
!___________________________________________________________________

  xstop=x+4.0d0*x**0.3333d0 +2.0d0
  nstop=xstop
  ymod=abs(y)
  nmx=dmax1(xstop,ymod) + 15
! dang=1.570796327d0/float(nang-1)
  dang=1.570796327d0/dble(nang-1)
  do 555 j = 1,nang
! theta(j)= (float(j)-1.0d0)*dang
   theta(j)= (dble(j)-1.0d0)*dang
555 amu(j)=dcos(theta(j))
!__________________________________________________________________
! logarithmic derivative d(j) calculated by downward recurrence
! beginning with initial value 0.0 + i*0.0 at j = nmx
! __________________________________________________________________

   d(nmx)=cmplx(0.0,0.0)
   nn=nmx-1
   do 120 n=1,nn
    rn=nmx-n+1
    d(nmx-n)=(rn/y)-(1.0d0/(d(nmx-n+1)+rn/y))
120 continue
    do 666 j=1,nang
     pi0(j)=0.0d0
     pi1(j)=1.0d0
666 continue
     nn=2*nang-1
     do 777 j=1,nn
      s1(j)=cmplx(0.0,0.0)
      s2(j)=cmplx(0.0,0.0)
777 continue
!__________________________________________________________________
! riccati bessel functions with real argument x calculated by upward
! recurrence
!__________________________________________________________________

    psi0=dcos(dx)
    psi1=dsin(dx)
    chi0=-1.0d0*dsin(x)
    chi1=dcos(x)
    apsi0=psi0
    apsi1=psi1
    xi0=cmplx(apsi0,-chi0)
    xi1=cmplx(apsi1,-chi1)
    qsca=0.0d0
    n=1
200 dn=n
    rn=n
    fn=(2.0d0*rn+1.0d0)/(rn*(rn+1.0d0))
    psi=(2.0d0*dn-1.0d0)*psi1/dx-psi0
    apsi=psi
    chi=(2.0d0*rn-1.0d0)*chi1/x -  chi0
    xi = cmplx(apsi,-chi)
    an=(d(n)/refrel+rn/x)*apsi - apsi1
    an=an/((d(n)/refrel+rn/x)*xi - xi1)
    bn=(refrel *d(n)+rn/x)*apsi - apsi1
    bn=bn/((refrel*d(n)+rn/x)*xi - xi1)
    qsca=qsca+(2.0d0*rn+1.0d0)*(abs(an)*abs(an)+abs(bn)*abs(bn))
  do 789 j=1,nang
     jj=2*nang-j
     pi(j)=pi1(j)
     tau(j)=rn*amu(j)*pi(j) - (rn+1.)*pi0(j)
     p=(-1.0d0)**(n-1)
     s1(j)=s1(j)+fn*(an*pi(j)+bn*tau(j))
     t=(-1.0d0)**n
     s2(j)=s2(j) + fn*(an*tau(j)+bn*pi(j))
     if (j .eq. jj) go to 789
     s1(jj)=s1(jj) + fn*(an*pi(j)*p + bn*tau(j)*t)
     s2(jj)=s2(jj) + fn*(an*tau(j)*t + bn*pi(j)*p)
789 continue
     psi0=psi1
     psi1=psi
     apsi1=psi1
     chi0=chi1
     chi1=chi
     xi1=cmplx(apsi1,-chi1)
     n=n+1
     rn=n
    do 999 j=1,nang
     pi1(j)=((2.0d0*rn-1.0d0)/(rn-1.0d0))*amu(j)*pi(j)
     pi1(j)=pi1(j) - rn*pi0(j)/(rn-1.0d0)
     pi0(j) = pi(j)
999 continue
     if((n-1-nstop).ge.0)then
      go to 300
     else
      go to 200
     endif
300 qsca=(2.0d0/(x*x))*qsca
     qext=(4.0d0/(x*x))*real(s1(1))
     qback=(4.0d0/(x*x))*abs(s1(2*nang -1))*abs(s1(2*nang -1))
!-----------------------------------------------------------------------

   return
end subroutine bhmie
!***********************************************************************

!***********************************************************************
subroutine Bolom(q,qbol)
!=======================================================================
! This subroutine integrates given radiation field, q (function of
! wavelength and radial coordinate), over wavelength. q is a matrix
! of physical size (npL,npY) [coming from paramet.inc] and real size
! (nL,nY) [coming from grids.inc], and qBol is an array of physical size
! (npY) and real size nY.                              [Z.I., Mar. 1996]
!=======================================================================
  use common, only: npL, nL, lambda, nY, npY
  implicit none

  integer iL, iY
  double precision q(npL,npY), qaux(npL), qbol(npY), resaux
!-----------------------------------------------------------------------
! loop over iY (radial coordinate)
  do iY = 1, nY
! generate auxiliary function for integration
! loop over iL (wavelength)
   do iL = 1, nL
    qaux(iL) = q(iL,iY)/lambda(iL)
   end do
   call Simpson(npL,1,nL,lambda,qaux,resaux)
   qbol(iY) = resaux
  end do
!-----------------------------------------------------------------------
  return
end subroutine Bolom
!***********************************************************************

!***********************************************************************
DOUBLE PRECISION FUNCTION Bessel(x)
!=======================================================================
! This function evaluates the Bessel function of the zeroth kind.
! Formulae are from Abramowitz & Stegun.               [Z.I., Jan. 1997]
! =======================================================================
     use common
      IMPLICIT none

      INTEGER i
      DOUBLE PRECISION x, c(6)
! -----------------------------------------------------------------------
      c(1) = -2.2499997D+00
      c(2) =  1.2656208D+00
      c(3) = -0.3163866D+00
      c(4) =  0.0444479D+00
      c(5) = -0.0039444D+00
      c(6) =  0.00021D+00
      Bessel=0.0D+00
      IF (x.LE.3.0D+00)THEN
        DO i=1,6
          Bessel = Bessel + c(i)*(x/3.0D+00)**(2.0D+00*i)
        END DO
        Bessel = 1.0D+00 + Bessel
        ELSE
        Bessel = dsqrt(2.0D+00/Pi/x) * dcos(x-Pi/4.0D+00)
      ENDIF
! -----------------------------------------------------------------------
      RETURN
END FUNCTION Bessel
!***********************************************************************

!***********************************************************************
subroutine Chk_deltau(tauaux,deltaumax)
! This checks Dejan's condition  [MN,May'10]
!=======================================================================
  use common
  implicit none

  integer iY, i, nn
  double precision etatemp(npY), tauaux(npY), deltaumax, x1, x2, result1, eta
  double precision, dimension(:), allocatable:: xg,  wg
  external eta

!-----------------------------------------------------------------------
  if(sph) then
!     save old grid and values of Eta (important for denstyp = 5 or 6)
!     for spherical case
      if (rdw) then
        Yprev = Y
        ETAtemp = ETAdiscr
        nYprev = nY
!     end if for rdw
      end if

     tauaux= 0.0d0
     do iY = 2, nY
       x1 = Y(iY-1)
       x2 = Y(iY)
       nn = 2*int(abs(x2 - x1)/10.0d0) + 11
       if(allocated(xg)) deallocate(xg)
       if(allocated(wg)) deallocate(wg)
       allocate(xg(nn))
       allocate(wg(nn))
       call gauleg(x1,x2,xg,wg,nn)
       result1 = 0.0d0
       do i = 1, nn
         result1 = result1 + eta(xg(i))*wg(i)
       end do
       tauaux(iY) = tauaux(iY-1) + tautot(1)*result1
      end do
!   maximal deltau is no more than 2 times the average value
    deltaumax = 2.0d0*tauaux(nY)*result1/nY
  elseif(slb) then
!     maximal deltau is no more than 2 times the average value
      deltaumax = 2.0d0*tautot(1)/nY
  end if
!-----------------------------------------------------------------------
  return
end subroutine Chk_deltau
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
            maxerr, error, slope, xmid, funSpline, aux, power, yR, yL
! ---------------------------------------------------------
!       check the midpoints
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
          error = DABS((funSpline-funmid(i))/funmid(i))
!         check for the deviation at the midpoint
          IF (error.GE.maxerr.OR.funSpline.LE.0.0) THEN
            slope = (fun(i+1) - fun(i)) / (x(i+1)-x(i))
            coef(i,1) = fun(i) - x(i) * slope
            coef(i,2) = slope
            coef(i,3) = 0.0
            coef(i,4) = 0.0
          END IF
!         check for the logarithmic derivative (only for RDW)
          IF(RDW) THEN
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
! ---------------------------------------------------------
      RETURN
end subroutine CHKSPLIN
! ***********************************************************************


!***********************************************************************
subroutine doProduct(nn,yt,pt,p0,j,prd)
!=======================================================================
! This is an auxiliary subroutine which evaluates a messy expression
! needed to calculate normalization constants for a broken power law
! density.                                             [Z.I., Aug. 1996]
!=======================================================================
  implicit none

  integer nn, i, j
  double precision yt(nn), pt(nn), prd, p0
!-----------------------------------------------------------------------

  prd = yt(1)**(pt(1) - p0)
  if (j.gt.1) then
   do i = 2, j
    prd = prd * yt(i)**(pt(i) - pt(i-1))
   end do
  end if
!-----------------------------------------------------------------------

  return
end subroutine doProduct
!***********************************************************************

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
!   if (x.gt.0.0d0.and.x.le.1.0d-15) x=1.0d-15
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
! -------------------------------------------------------------------------

  if(x.lt.0.0d0) x = dabs(x)
  if (x.lt.1.0d-15) then
   eint2 = 1.0d0
  else
   eint2 = dexp(-x) - x*eint1(x)
  end if
! -------------------------------------------------------------------------
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
!-------------------------------------------------------------------------

  if(x.le.0.0d0) x = dabs(x)
  if (x.lt.1.0d-15) then
   eint3 = 1.0d0/(2.0d0)
  else
   eint3 = (dexp(-x) - x*eint2(x))/(2.0d0)
  end if
! -------------------------------------------------------------------------
  return
end function eint3
!**********************************************************************

!**********************************************************************
double precision function eint4(x)
!======================================================================
! Needed for the slab geometry. It calculates the fourth exponential
! integral E4(x) by the recurrence f-la. (see Abramovitz & Stegun,1994)
!                                                         [MN,Jan'97]
! ======================================================================
  implicit none
  double precision x, eint3
! -------------------------------------------------------------------------

  if(x.lt.0.0d0) x = dabs(x)
  if (x.lt.1.0d-15) then
   eint4 = 1.0d0/3.0d0
  else
   eint4 = (dexp(-x)-x*eint3(x))/3.0d0
  end if
!-------------------------------------------------------------------------
  return
end function eint4
!**********************************************************************

!**********************************************************************
double precision function eint5(x)
!======================================================================
! Needed for the slab geometry. It calculates the fifth exponential
! integral E5(x) by the recurrence f-la. (see Abramovitz & Stegun,1994)
!                                                         [MN,Jan'97]
! ======================================================================
  implicit none
  double precision x, eint4
! -----------------------------------------------------------------------

  if(x.lt.0.0d0) x = dabs(x)
  if (x.lt.1.0d-15) then
   eint5 = 1./(4.0d0)
  else
   eint5 = (dexp(-x)-x*eint4(x))/(4.0d0)
  end if
! -----------------------------------------------------------------------
  return
end function eint5
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
! -----------------------------------------------------------------------
      l = LEN(line)
      EMPTY = 1
      iTeX = 0
      DO i = 1, l
        ch = line(i:i)
        IF(EMPTY.EQ.1.AND.ch.EQ.'%') iTeX = 1
         IF (ch.NE.' ') EMPTY = 0
      END DO
      IF (iTeX.EQ.1) EMPTY = 1
! -----------------------------------------------------------------------
      RETURN
      END
!***********************************************************************

!***********************************************************************
subroutine gauleg(x1,x2,xg,wg,n)
!=====================================================================
  implicit none

  integer i,m,n,j
  double precision x1,x2,xm,xl,eps,delj,p,eta
  double precision xg(n),wg(n),sum,ff,p1,p2,p3,z1,z,pp
  parameter (eps=1.0d-14)
! ---------------------------------------------------------------------

  xg = 0.0d0
  wg = 0.0d0
  m = int((n+1)/2)
  xm = 0.5d0*(x2+x1)
  xl = 0.5d0*(x2-x1)
  do 12 i = 1, m
   z = cos(3.1415926535898d0*(dble(i) - 0.25d0)/(dble(n) + 0.5d0))
1  continue
   p1 = 1.0d0
   p2 = 0.0d0
   do 11 j = 1, n
    p3 = p2
    p2 = p1
    p1 = ((2.0d0*dble(j)-1.0d0)*z*p2-(j-1.0d0)*p3)/dble(j)
11  continue
    pp = dble(n)*(z*p1 - p2)/(z*z - 1.0d0)
    z1 = z
    z = z1-p1/pp
    if (abs(z-z1).gt.eps) go to 1
    xg(i) = xm-xl*z
    xg(n+1-i) = xm+xl*z
    wg(i) = 2.d0*xl/((1.d0-z*z)*pp*pp)
    wg(n+1-i) = wg(i)
12  continue
!-----------------------------------------------------------------------

    return
end subroutine gauleg
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
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
  return
end subroutine getfs
!***********************************************************************

!***********************************************************************
integer function Kron(i1,i2)
!=======================================================================
! This function is Kronecker delta-function defined as:
! Kron(i1,i2) = 1 for i1=i2
! Kron(i1,i2) = 0 otherwise.                           [Z.I., Dec. 1995]
!=======================================================================

   implicit none
   integer i1, i2
!-----------------------------------------------------------------------
   if (i1.eq.i2) then
    kron = 1
   else
    kron = 0
   end if
!-----------------------------------------------------------------------
   return
end function Kron
!***********************************************************************

!***********************************************************************
subroutine LinInter(nn,n,x,y,xloc,iNloc,Yloc)
!=======================================================================
! This subroutine performs linear interpolation for y(x) such that
! Yloc = y(xloc). It is assumed that x is monotonously increasing.
!                                                      [Z.I., Mar. 1996]
!=======================================================================

  implicit none
  integer nn, n, i, istop, iNloc
  double precision x(nn), y(nn), xloc, Yloc
!-----------------------------------------------------------------------

  if (n.gt.1) then
   if ((x(1)-xloc)*(x(n)-xloc).le.0.0d0) then
    istop = 0
    i = 1
    do while (istop.ne.1)
     i = i + 1
     if (i.gt.n) stop 'lininter ???'
     if (x(i).ge.xloc) then
      istop = 1
      iNloc = i
      Yloc = y(i-1) + (y(i)-y(i-1))/(x(i)-x(i-1))*(xloc-x(i-1))
     end if
    end do
   else
    if (xloc.le.x(1)) Yloc = y(1)
    if (xloc.ge.x(n)) Yloc = y(n)
   end if
  else
   Yloc = y(1)
  end if
!-----------------------------------------------------------------------

  return
end subroutine LinInter
!***********************************************************************

!***********************************************************************
      SUBROUTINE LINSYS(Nreal,A,B,X,error)
!=======================================================================
! This subroutine solves the set of linear equations [A]*[X] = [B] for
! X [A(k,1)*X(1)+A(k,2)*X(2)+...+A(k,Nreal)*X(Nreal) = B(k), k=1,Nreal).
! The real size of matrix A is Nreal x Nreal and its physical dimension
! is npY x npY, where npY comes from INCLUDE 'userpar.inc'. Both vectors
! B and X have real lengths Nreal. The set is solved by calls to LUDCMP
! and LUBKSB and the solution is improved subsequently by a call to
! MPROVE. These three subroutines are taken from Numerical Recipes.
!                                                      [Z.I., Nov. 1995]
! =======================================================================
     use common
      IMPLICIT none

      INTEGER Nreal, indx(npY), i, j, error
      DOUBLE PRECISION A(npY,npY), B(npY), X(npY), &
              A1c(npY,npY), B1(npY), A2c(npY,npY), B2(npY), d
! -----------------------------------------------------------------------
      error = 0
! generate DOUBLE PRECISION copies of A and B (two copies because they
! are changed in LUDCMP and LUBKSB, but still needed for MPROVE)
      DO i = 1, Nreal
        B1(i) = B(i)
        B2(i) = B(i)
        DO j = 1, Nreal
           A1c(i,j) = A(i,j)
           A2c(i,j) = A(i,j)
        END DO
      END DO
!     solve the system
      CALL LUDCMP(A1c,Nreal,npY,indx,d,error)
      IF (error.NE.0) RETURN
      CALL LUBKSB(A1c,Nreal,npY,indx,B1)
!     improve the solution (saved in B)
      CALL MPROVE(A2c,A1c,Nreal,npY,indx,B2,B1)
!     copy the improved solution to output vector X
      DO i = 1, Nreal
        X(i) = B1(i)
      END DO
! -----------------------------------------------------------------------
      RETURN
      END SUBROUTINE LINSYS
!***********************************************************************

! ***********************************************************************
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
! =======================================================================
      DIMENSION INDX(NP)
      DOUBLE PRECISION A(NP,NP),B(NP)
! -------------------------------------------------------------------
      II=0
      DO 12 I=1,N
      LL=INDX(I)
      SUM=B(LL)
      B(LL)=B(I)
      IF (II.NE.0)THEN
        DO 11 J=II,I-1
          SUM=SUM-A(I,J)*B(J)
11        CONTINUE
      ELSE IF (SUM.NE.0.) THEN
        II=I
      ENDIF
      B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
      SUM=B(I)
      IF(I.LT.N)THEN
        DO 13 J=I+1,N
          SUM=SUM-A(I,J)*B(J)
13        CONTINUE
      ENDIF
      B(I)=SUM/A(I,I)
14    CONTINUE
! -------------------------------------------------------------------
      RETURN
      END
! ***********************************************************************

! ***********************************************************************
      SUBROUTINE LUDCMP(A,N,NP,INDX,D,error)
! =======================================================================
      PARAMETER (NMAX=10000,TINY=1.0E-20)
      DIMENSION INDX(NP)
      INTEGER error
      DOUBLE PRECISION A(NP,NP),VV(NMAX), D, SUM
! -------------------------------------------------------------------
      error = 0
      D = 1.
      DO I = 1, N
       AAMAX=0.
       DO J = 1, N
        IF (DABS(A(I,J)).GT.AAMAX) AAMAX=DABS(A(I,J))
       END DO
!       IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
       IF (AAMAX.EQ.0.) THEN
        error = 5
        RETURN
       ENDIF
       VV(I)=1./AAMAX
      END DO
      DO J = 1 , N
       IF (J.GT.1) THEN
        DO I = 1, J-1
          SUM=A(I,J)
          IF (I.GT.1)THEN
            DO K = 1, I-1
             SUM=SUM-A(I,K)*A(K,J)
            END DO
            A(I,J)=SUM
          ENDIF
        END DO
       ENDIF
       AAMAX=0.
       DO I = J, N
        SUM=A(I,J)
        IF (J.GT.1)THEN
          DO K = 1, J-1
            SUM=SUM-A(I,K)*A(K,J)
          END DO
          A(I,J)=SUM
        ENDIF
        DUM=VV(I)*DABS(SUM)
        IF (DUM.GE.AAMAX) THEN
          IMAX=I
          AAMAX=DUM
        ENDIF
       END DO
       IF (J.NE.IMAX)THEN
        DO K = 1, N
          DUM=A(IMAX,K)
          A(IMAX,K)=A(J,K)
          A(J,K)=DUM
        END DO
        D=-D
        VV(IMAX)=VV(J)
       ENDIF
       INDX(J)=IMAX
       IF(J.NE.N)THEN
        IF(A(J,J).EQ.0.)A(J,J)=TINY
        DUM=1./A(J,J)
        DO I = J+1, N
          A(I,J)=A(I,J)*DUM
        END DO
       ENDIF
      END DO
      IF(A(N,N).EQ.0.)A(N,N)=TINY
! -------------------------------------------------------------------
      RETURN
      END
! ***********************************************************************

! ***********************************************************************
      SUBROUTINE Maple3(w,z,p,MpInt)
! =======================================================================
! This function calculates indefinite integral:
!    MpInt(iC) = INT(w^(2-iC) / sqrt(w^2-p^2) * dw), for iC=1,2,3,4.
!                                                      [Z.I., Apr. 1996]
! =======================================================================
      IMPLICIT none
      DOUBLE PRECISION w, z, p, MpInt(4)
! -----------------------------------------------------------------------
!     integrals
      MpInt(1) = z
      MpInt(2) = dlog(w+z)
      IF (p.GT.0.0) THEN
        MpInt(3) = dacos(p/w)/p
        MpInt(4) = z/w/p/p
        ELSE
        MpInt(3) = -1.0 / w
        MpInt(4) = -0.5 / w / w
      END IF
! -----------------------------------------------------------------------
      RETURN
      END
! ***********************************************************************


!***********************************************************************
subroutine Mie(npL,nL,lambda,ere,eim,npA,na,a,ng1,qabs,qsca)
!=======================================================================
! This subroutine calculates qabs and qsca for a given diffractive
! index ere, eim, wavelength lambda and size a. here, lambda is an
! array (1..nL), ere and eim are given on this array, a is an array
! of sizes (1..na). qabs and qsca are arrays (ng1..ng1+na,nL), i.e. for
! each wavelength lambda, qabs and qsca are evaluated for na different
! sizes. the numbering, however, does not start from 1, but rather from
! ng1.                                                [Z.I., Aug. 1996]
!=======================================================================

  implicit none
  integer npL, nL, npA, na, ng1, iL, ia
  double precision lambda(npL), ere(npL), eim(npL), a(npA),   &
       qabs(npA,npL), qsca(npA,npL)
  double precision xx, qex, qsc, qback
  complex refrel, s1(200), s2(200)
!-----------------------------------------------------------------------

! loop over wavelengths
  do iL = 1, nL
! complex index of refraction
   refrel = cmplx(ere(iL),eim(iL))
! loop over sizes
   do ia = 1, na
! size parameter
    xx=2.0d0*3.14159265d0*a(ia)/lambda(iL)
! if size parameter xx>100 use xx=100 (geometrical optics)
    if (xx.gt.100.0) xx = 100.0d0
! calculate efficiencies
    call bhmie(xx,refrel,2,s1,s2,qex,qsc,qback)
! store the result
    qabs(ng1+ia-1,iL) = qex - qsc
    qsca(ng1+ia-1,iL) = qsc
   end do
  end do
!-----------------------------------------------------------------------
  return
end subroutine Mie
!***********************************************************************


! ***********************************************************************
      SUBROUTINE MPROVE(A,ALUD,N,NP,INDX,B,X)
! =======================================================================
      PARAMETER (NMAX=10000)
      DIMENSION INDX(N)
      DOUBLE PRECISION SDP,A(NP,NP),ALUD(NP,NP),B(N),X(N),R(NMAX)
! -----------------------------------------------------------------------
      DO i = 1, N
       SDP = -B(i)
       DO j = 1, N
         SDP = SDP + A(i,j)*X(j)
       END DO
       R(i) = SDP
      END DO
      CALL LUBKSB(ALUD,N,NP,INDX,R)
      DO i = 1, N
       X(i) = X(i) - R(i)
      END DO
! -----------------------------------------------------------------------
      RETURN
      END
! ***********************************************************************

!***********************************************************************
SUBROUTINE NORDLUND(flag,x,f,N1,N2,m,intfdx,error)
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

      INTEGER i, flag, N1, N2, N2n, Nanal, m, first, error
      DOUBLE PRECISION x(npP), f(npP), wSimp(npP), wCorr(npP), wC1, wC2,  &
             am, intfdx, xaux(4), faux(4), aux
! -----------------------------------------------------------------------
      error = 0
! parameter 'first' selects choice for derivatives at boundary points.
! For first.EQ.0 derivatives are 0 and first*(f2-f1)/(x2-x1) otherwise.
! first=1 works much better for functions encountered here.
      first = 1
! number of points for analytic integration
      Nanal = 4
! do not include points for analytic integration
      IF (flag.EQ.1.AND.N2.GT.N1+Nanal) THEN
        N2n = N2 - Nanal + 1
      ELSE
        N2n = N2
      END IF
! set integral to 0 and accumulate result in the loop
      intfdx = 0.0
! generate weighting factors, w(i), and integrate in the same loop
      DO i = N1, N2n
! first usual Simpson factors (Nordlund, eq. III-14, first term)...
        IF (i.NE.N1.AND.i.NE.N2n) THEN
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
!       add contribution to the integral
        IF (m.EQ.0) THEN
          intfdx = intfdx + f(i) * (wSimp(i) + wCorr(i))
         ELSE IF(m.EQ.1) THEN
            intfdx = intfdx + x(i)*f(i)*(wSimp(i) + wCorr(i))
         ELSE IF(m.EQ.2) THEN
            intfdx = intfdx + x(i)*x(i)*f(i)*(wSimp(i) + wCorr(i))
        END IF
      END DO
!     change the sign (x [i.e. mu] array is in descending order!!!)
      intfdx = -intfdx
!     if flag=1 use analytic approximation for the last Nanal points
      IF (flag.EQ.1.AND.N2n.GT.N1+Nanal) THEN
!       generate auxiliary arrays for ANALINT
        DO i=1,Nanal
          xaux(i) = x(N2n+Nanal-i)
          faux(i) = f(N2n+Nanal-i)
        END DO
!     calculate the contribution of the last Nanal points
!       produce REAL copy of m
        am = 1.0*(m)
        CALL ANALINT(Nanal,xaux,faux,am,aux,error)
        IF(error.NE.0) THEN
          RETURN
        END IF
!       add the contribution of the last Nanal points
        intfdx = intfdx + aux
      END IF
! -----------------------------------------------------------------------
      RETURN
END SUBROUTINE NORDLUND
!**********************************************************************


!***********************************************************************
subroutine polint(xa,ya,n,x,y,dy)
! For polinomial interpolation, used in Subroutine Romby.
!=======================================================================
  implicit none
  integer n,Nmax
  double precision dy,x,y,xa(n),ya(n)
  parameter (Nmax=1000)
  integer i,m,ns
  double precision den,dif,dift,ho,hp,w,c(Nmax),d(Nmax)
!-----------------------------------------------------------------------
  c = 0.0d0
  d = 0.0d0
  ns=1
  dif=dabs(x-xa(1))
  do 11 i=1,n
   dift=dabs(x-xa(i))
   if (dift.lt.dif) then
    ns=i
    dif=dift
   endif
   c(i)=ya(i)
   d(i)=ya(i)
11 continue
   y=ya(ns)
   ns=ns-1
   do 13 m=1,n-1
    do 12 i=1,n-m
     ho=xa(i)-x
     hp=xa(i+m)-x
     w=c(i+1)-d(i)
     den=ho-hp
     if(den.eq.0.0d0) then
      write(6,'(a)') 'failure in polint'
      stop
     endif
     den=w/den
     d(i)=hp*den
     c(i)=ho*den
12   continue
     if (2*ns.lt.n-m)then
      dy=c(ns+1)
     else
      dy=d(ns)
      ns=ns-1
     endif
     y=y+dy
13  continue
!-----------------------------------------------------------------------
  return
 end subroutine polint
!***********************************************************************

!***********************************************************************
subroutine PowerInt(n,n1,n2,x,y,integral)
!=======================================================================
! This subroutine calculates integral I(y(x)*dx). Both y and x are
! 1D arrays, y(i), x(i) with i=1,N (declared with NN). Lower and upper
! integration limits are x(N1) and x(N2), respectively. The method used
! is a power-law approximation for y(x) between any two points .
! (This subroutine is used for integration over size distribution) [ZI,'96]
!=======================================================================

  implicit none
  integer i, n, n1, n2
  double precision x(n), y(n), integral, pow, c, delint
!-----------------------------------------------------------------------
! set integral to 0 and accumulate result in the loop
  integral = 0.0d0
! calculate weight, wgth, and integrate in the same loop
  if (n2.gt.n1) then
   do i = n1, n2-1
    pow = dlog(Y(i+1)/Y(i)) / dlog(x(i+1)/x(i))
    c = Y(i) / x(i)**pow
    delint=(x(i+1)**(pow+1.0d+0)-x(i)**(pow+1.0d+0))*c/(pow+1.0d+0)
! add contribution to the integral
    integral = integral + delint
   end do
  else
   integral = 0.0d0
! this was in case of single size grains
! integral = Y(1)
  end if
!-----------------------------------------------------------------------

  return
end subroutine PowerInt
!***********************************************************************

!***********************************************************************
subroutine PowerInter(nn,n,x,y,xloc,iNloc,Yloc)
!=======================================================================
! This subroutine performs power law interpolation for y(x) such that
! Yloc = y(xloc). It is assumed that x is monotonously increasing.
! [based on sub LinInter by ZI'96, modified for power law interp. by MN'03]
!=======================================================================

  implicit none
  integer nn, n, i, istop, iNloc
  double precision x(nn), y(nn), xloc, Yloc, pow
! -----------------------------------------------------------------------

  if (n.gt.1) then
   if ((x(1)-xloc)*(x(n)-xloc).le.0.0d0) then
     istop = 0
     i = 1
     do while (istop.ne.1)
      i = i + 1
      if (i.gt.n) stop 'powinter ???'
      if (x(i).ge.xloc) then
        istop = 1
        iNloc = i
        if ((y(i)*y(i-1)).gt.0.0d0) then
          pow = dlog(y(i)/y(i-1))/dlog(x(i)/x(i-1))
          Yloc = y(i-1)*((xloc/x(i-1))**pow)
        else
          Yloc = 0.0d0
        end if
      end if
     end do
   else
     if (xloc.le.x(1)) Yloc = y(1)
     if (xloc.ge.x(n)) Yloc = y(n)
   end if
  else
    Yloc = y(1)
  end if
!-----------------------------------------------------------------------
  return
end subroutine PowerInter
!***********************************************************************

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
      END
! ***********************************************************************


!***********************************************************************
subroutine ROMBY(fnc,a,b,ss)
!=======================================================================
! This subroutine performs Romberg integration of function func on
! interval [a,b]. The result is returned in ss. Desired accuracy is set
! to 0.002.                                            [Z.I., Feb. 1996]
! =======================================================================
      IMPLICIT NONE
      INTEGER JMAX,JMAXP,K,KM, J
      PARAMETER (JMAX=30, JMAXP=JMAX+1, K=3, KM=K-1)
      DOUBLE PRECISION a,b,fnc,ss,EPS_loc, aux, dss,h(JMAXP),s(JMAXP)
      EXTERNAL fnc
! -----------------------------------------------------------------------
      EPS_loc = 0.002d0
      h(1)=1.0d0
      do 11 j=1,JMAX
        call trapzd(fnc,a,b,s(j),j)
        if (j.ge.K) then
          aux = 0.0d0
          call polint(h(j-KM),s(j-KM),K,aux,ss,dss)
          IF (dabs(dss).le.EPS_loc*dabs(ss)) RETURN
        endif
        s(j+1)=s(j)
        h(j+1)=0.25d0*h(j)
11    continue
! -----------------------------------------------------------------------
      RETURN
END subroutine ROMBY
!***********************************************************************

!***********************************************************************
subroutine scaleto1(Nmax,n,y)
!=======================================================================
! This subroutine scales vector Y such that Y(1) = 1.0
!                                                      [Z.I., Jan. 1997]
!=======================================================================
  implicit none

  integer Nmax, N, i
  double precision Y(Nmax), Scale
!-----------------------------------------------------------------------
  Scale = Y(1)
  do i = 1, N
   Y(i) = Y(i) / Scale
  end do
!-----------------------------------------------------------------------
  return
end subroutine scaleto1
!***********************************************************************

! ***********************************************************************
      SUBROUTINE ScaletoArea(Nmax,N,X,Y,Area)
! =======================================================================
! This subroutine scales a function Y(x) by the area A=Int{2*Pi Y(x)xdx}.
! X and Y are 1D arrays. (Used for PSF normalization.)       [MN, Sep'04]
! =======================================================================
      IMPLICIT none
      INTEGER Nmax, N, i
      DOUBLE PRECISION Y(Nmax),X(Nmax),Fn(Nmax),Area,Pi
! -----------------------------------------------------------------------
      Pi = 2.0D+00*ASIN(1.0)
      DO i = 1, N
	  Fn(i) = Y(i)*X(i)
      END DO
!     Integrate:
      Area = 0.0D+00
      DO i = 1, N-1
        Area = Area + 0.5D+00*(Fn(i+1)+Fn(i))*(X(i+1)-X(i))
	END DO
	Area = 2.0D+00*Pi*Area
!     Normalize:
      DO i = 1, N
        Y(i) = Y(i) / Area
      END DO
! -----------------------------------------------------------------------
      RETURN
      END
! ***********************************************************************

!***********************************************************************
subroutine shift(x,Nmax,n,xins,i)
!=======================================================================
! Rearranges a vector X by inserting a new element Xins.    [MN, Aug'96]
! =======================================================================
  implicit none
  integer Nmax, n, i,j
  double precision x(Nmax),xins
! -----------------------------------------------------------------------

  do j = n+1, i+2, -1
   x(j) = x(j-1)
  end do
  x(i+1) = xins
! -----------------------------------------------------------------------
  return
end subroutine shift
!***********************************************************************

!***********************************************************************
subroutine Simpson(n,n1,n2,x,y,integral)
!=======================================================================
! This subroutine calculates integral I(y(x)*dx). Both y and x are
! 1D arrays, y(i), x(i) with i=1,N (declared with NN). Lower and upper
! integration limits are x(N1) and x(N2), respectively. The method used
! is Simpson (trapezoid) approximation. The resulting integral is sum of
! y(i)*wgth, i=N1,N2.                                  [Z.I., Mar. 1996]
! =======================================================================
  implicit none
  integer i, n, n1, n2
  double precision x(n), y(n), wgth, integral, dyn2
! -----------------------------------------------------------------------
  dyn2 = 0.0d0

! set integral to 0 and accumulate result in the loop
  integral = 0.0d0
! calculate weight, wgth, and integrate in the same loop
  if (n2.gt.n1) then
   do i = n1, n2
! weigths
    if (i.ne.n1.and.i.ne.n2) then
     wgth = 0.5d0*(x(i+1)-x(i-1))
    else
     if (i.eq.n1) wgth = 0.5d0*(x(n1+1)-x(n1))
     if (i.eq.n2) wgth = 0.5d0*(x(n2)-x(n2-1))
    end if
! add contribution to the integral
    integral = integral + y(i)*wgth
   end do
  else
   integral = 0.0d0
  end if
! -----------------------------------------------------------------------
  return
end subroutine Simpson
!***********************************************************************

!***********************************************************************
subroutine sort(ra,n)
!=======================================================================
  implicit none

  integer i,j,l,n,ir
  double precision ra(n),rra
!-----------------------------------------------------------------------
  l=n/2+1
  ir=n
10 continue
  if(l.gt.1)then
   l=l-1
   rra=ra(l)
  else
   rra=ra(ir)
   ra(ir)=ra(1)
   ir=ir-1
   if(ir.eq.1)then
    ra(1)=rra
    return
   endif
  endif
  i=l
  j=l+l
20 if(j.le.ir)then
   if(j.lt.ir)then
    if(ra(j).lt.ra(j+1))j=j+1
   endif
   if(rra.lt.ra(j))then
    ra(i)=ra(j)
    i=j
    j=j+j
   else
    j=ir+1
   endif
   go to 20
  endif
  ra(i)=rra
  go to 10
! -----------------------------------------------------------------------
  return
end subroutine sort
! ***********************************************************************

!***********************************************************************
SUBROUTINE Spline(x,y,n,yp1,ypn,y2)
!=======================================================================
      INTEGER n,NMAX
      DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)
      INTEGER i,k
      DOUBLE PRECISION p,qn,sig,un,u(NMAX)
! -----------------------------------------------------------------------
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i)) &
        -(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
! -----------------------------------------------------------------------
  return
end subroutine Spline
!***********************************************************************


!***********************************************************************
SUBROUTINE SPLINE2(x,fun,N,coef)
! =======================================================================
! This subroutine finds coefficients coef(i,j) such that
! fun(x)=coef(i,1) + coef(i,2)*x + coef(i,3)*x^2 + coef(i,4)*x^3
! for x(i).LE.x.LE.x(i+1) is a cubic spline approximation of fun(x),
! with i=1..N.                                         [Z.I., Feb. 1995]
! =======================================================================
   use common
      IMPLICIT none

      INTEGER N, i
      DOUBLE PRECISION x(npY), coef(npY,4), secnder(npY), y2at1, y2atN, &
           Dd, xL, xR, dR, dL, fun(npY), fL, fR
! -----------------------------------------------------------------------
!     find second derivative, secnder
        y2at1 = (fun(2)-fun(1))/(x(2)-x(1))
        y2atN = (fun(N)-fun(N-1))/(x(N)-x(N-1))
        CALL SPLINE(x,fun,N,y2at1,y2atN,secnder)
!     generate coef(i,j), j=1,2,3,4
      DO i = 1, N-1
        Dd = x(i+1) - x(i)
        xL = x(i)
        xR = x(i+1)
        dL = secnder(i)
        dR = secnder(i+1)
        fL = fun(i)
        fR = fun(i+1)
        coef(i,1) = (xR*fL-xL*fR)/Dd + dL*xR*Dd/6.*((xR/Dd)**2.-1.)
        coef(i,1) = coef(i,1) - dR*xL*Dd/6. *((xL/Dd)**2.-1.)
        coef(i,2) = (fR-fL)/Dd + dL*Dd/6.*(1.-3.*(xR/Dd)**2.)
        coef(i,2) = coef(i,2) - dR*Dd/6.*(1.-3.*(xL/Dd)**2.)
        coef(i,3) = (dL*xR-dR*xL)/Dd/2.
        coef(i,4) = (dR-dL)/6./Dd
      END DO
! -----------------------------------------------------------------------
      RETURN
END SUBROUTINE SPLINE2
!***********************************************************************

!***********************************************************************
SUBROUTINE trapzd(func,a,b,s,n)
! =======================================================================
      IMPLICIT NONE
      INTEGER n
      DOUBLE PRECISION a,b,s,func
      EXTERNAL func
      INTEGER it,j
      DOUBLE PRECISION del,sum,tnm,x
! -------------------------------------------------------------------------
      IF (n.eq.1) THEN
        s=0.5d0*(b-a)*(func(a)+func(b))
      ELSE
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5d0*del
        sum=0.
        DO j = 1, it
          sum=sum+func(x)
          x=x+del
        END DO
        s=0.5d0*(s+(b-a)*sum/tnm)
      END IF
! -------------------------------------------------------------------------
      RETURN
END SUBROUTINE trapzd
!***********************************************************************

!***********************************************************************
!************************  THE END *************************************
!***********************************************************************
