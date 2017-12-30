!DUSTY RELEASE
Module common
  implicit none
  double precision pi, sigma, Gconst, r_gd, clight, mprot
  parameter (pi=3.141592653589793116)   ! pi
  parameter (sigma=5.67D-08)            ! Stefan-Boltzman constant
  parameter (Gconst= 6.67D-11)          ! gravitational constant
  parameter (r_gd = 2.0D+02)            ! gas to dust ratio
  parameter (clight = 3.0D+08)          ! speed of light
  parameter (mprot = 1.67D-27)          ! proton mass
  integer max_threads
  parameter (max_threads = 16)          ! max number of threads if parallel
  integer NREC
  PARAMETER (NREC = 1000)
  double precision :: dynrange          ! dynamical range (1d-20 .. 1d20)
  PARAMETER (dynrange = 1.d-30)
  double precision :: accRomb
  PARAMETER (accRomb = 1e-6)            ! 1e-5 neeeded for tauv=1000 1/r2 in solve_matrix
  integer nOutput                       ! number of output columns
  PARAMETER (nOutput=20)
  integer npY,npP,npL,npG,npR
  !PARAMETER (npY = 1000)
  !PARAMETER (npP = 20*npY+20)
  PARAMETER (npL = 119)
  !PARAMETER (npG = 16)
  !PARAMETER (npR = 90)
  include "userpar.inc"
  double precision,allocatable :: ETAdiscr(:)
  integer,allocatable ::  iYfirst(:), YPequal(:), Plast(:)
  integer left,right

  integer :: Ntr,nYetaF,nPcav
  integer ::  transmit !<----needed in slb intensity .... bad style!!!
  double precision :: muobs,taut,Sfn !<----needed in slb intensity .... bad style!!!
  integer :: startyp(2) 
  double precision Pow
  double precision Yout
  double precision aveV,aveA,sigExfid
  double precision accConv,accuracy

  double precision,allocatable :: Gamma(:)
  double precision,allocatable :: destroyed(:,:)
  double precision,allocatable :: ptr(:)
  double precision,allocatable :: Ytr(:)
  double precision,allocatable :: yetaf(:)
  double precision,allocatable :: etaf(:)
  double precision,allocatable :: ETAzp(:,:)
  double precision,allocatable :: abund(:,:)
  double precision,allocatable :: omega(:,:)
  double precision,allocatable :: TAUslb(:,:)
  double precision,allocatable :: SLBIntm(:,:)
  double precision,allocatable :: SLBIntp(:,:)
  double precision,allocatable :: TAUtot(:)
  double precision,allocatable :: Jext(:)
  double precision,allocatable :: fsL(:,:), fsR(:,:), fsLbol(:), fsRbol(:),fsbol(:)
  double precision,allocatable :: fds(:,:)
  double precision,allocatable :: fde(:,:)
  double precision,allocatable :: ftot(:,:)
  double precision,allocatable :: Utot(:,:)
  double precision,allocatable :: Utot_old(:,:)
  double precision,allocatable :: Ude(:,:)
  double precision,allocatable :: Uds(:,:)
  double precision,allocatable :: Td(:,:)
  double precision,allocatable :: Td_old(:,:)
  double precision,allocatable :: tauF(:)
  double precision,allocatable :: RPr(:)
  double precision,allocatable :: Eps(:)
  double precision,allocatable :: ugas(:)
  double precision,allocatable :: vrat(:,:)
  double precision,allocatable :: IntOut(:,:)
  double precision,allocatable :: IstR(:)
  double precision,allocatable :: bOut(:)
  double precision,allocatable :: tauZout(:)
  double precision,allocatable :: fbol(:)
  double precision,allocatable :: fpbol(:)
  double precision,allocatable :: fmbol(:)
  double precision,allocatable :: ubol(:)
  double precision,allocatable :: qF(:)
  double precision,allocatable :: rg(:,:)
  double precision,allocatable :: Intens(:,:)
  double precision,allocatable :: tauFdyn(:)
  double precision :: SmC(30,99)
  double precision :: gmax

  character*4 :: version
  parameter (version='4.00-$Rev$')
  integer :: error,warning,iverb
  double precision :: Tstar(2) ! needed to calculate intensities ???
  double precision, allocatable :: lambda(:)  ! lambda grid dusty
  double precision, allocatable :: shpL(:)    ! shape left side ilumination
  double precision, allocatable :: shpR(:)    ! shape right side ilumination
  double precision,allocatable  :: Tsub(:)    ! sublimation temperature for each grain
  double precision,allocatable  :: Tinner(:)  ! Temperature at the inner boundary for fiducial grain
  double precision,allocatable  :: Y(:)       ! Y - grid
  double precision,allocatable  :: Yprev(:)   ! previous Y - grid
  double precision,allocatable  :: P(:)       ! P - grid
  double precision :: TAUin(Nrec) ! user specified input tau grid
  integer :: nL                 ! size of lambda grid
  integer :: nG                 ! number of grain types
  logical :: slb,sph,sph_matrix     ! geometry
  ! remove TypEntry ???? --- not necessary ??? 
  integer :: TypEntry(2)        ! type of ilumination for (0-left and 1-right)
  integer :: denstyp            ! density type  1(POWD) 2(EXPD) 3(RDW) 4(RDWA) 5(USER_SUPPLIED) 6(RDWPR)
  double precision ksi          ! the ratio of the right/left bol. fluxes (<1) for slab
  double precision Ji,Jo        ! input mean energy density - scaling in dusty
  
  double precision mu1,mu2      ! cosine of ilumination anlge slab case
  integer :: iFidG              ! id number of fiducial grain
  double precision :: lamfid    ! fiducial wavelength [micron]
  integer :: iLfid              ! index fiducial wavelength
  double precision,allocatable  :: SigmaA(:,:) ! absorbtion crosssection
  double precision,allocatable  :: SigmaS(:,:) ! scattering crosssection
  double precision :: tauFid    ! optical depth at fiducial wavelength
!  integer :: ver                ! RDW variable
  double precision :: accFlux   ! flux accuracy
  double precision :: accTemp   ! temperature accuracy ((1.+accFlux)**(1./4.)-1.)*0.1

  double precision :: psi,psi0   ! psi as defined in IE97
  double precision :: ETAcoef(npY,4)
  ! output parameter - should be a structure and not in global ... search if it is used elseweher
  double precision CMdot, CM, Cve, Cr1, G1, Ginf, Phi, Prdw, QV, Qstar, r1rs, winf, &
       zeta1, I1_dyn, I2_dyn, I3_dyn, PIrdw
  integer psftype, Npsf, iLambda
  integer :: iA,iB,iC,iX,iPSF,iV,NlambdaOut,Nconv,nMu,Nvisi,iD,iJ,nJOut
  double precision,allocatable :: theta(:)
  double precision :: lambdaOut(nOutput), FWHM1(nOutput), FWHM2(nOutput),kPSF(nOutput), &
       xpsf(1000), ypsf(1000), psfArea(nOutput), YJout(nOutput), theta1
  double precision ConvInt(nOutput,1000), Visib(nOutput,1000),     &
       Offset(1000), qtheta1(1000), Te_min, JOut(npL,nOutput)


!!$  ! ===========================================================================
!!$  ! Commons with various parameters related to numerical accuracy.
!!$  ! All are initialized in subroutine INPUT.                [ZI'95;MN'99; MD'07]
!!$  ! Former COMMON /mumerics/
!!$  ! ===========================================================================
!!$  ! accConv - desired accuracy for energy.density convergence, by default it is
!!$  !           set to accConv = accuracy/500.
!!$  ! accFbol - desired accuracy for bol.flux convergence, by default
!!$  !           accFbol = 10.*accConv
!!$  ! accRomb - desired accuracy for numerical integration in ROMBERG
!!$  ! EtaRat  - limit on the ratio of density profile for any two consecutive
!!$  !           radial grid points, used in Ygrid.
!!$  ! delTAUs! - max. increase of scaled TAU/TAUtot, used in Ygrid
!!$  ! fac! - max. increase in the ratio of two y, used in Ygrid
!!$  ! Ncav - number of p-rays inside the cavity
!!$  ! Nins - number of p-rays inserted per radial step
!!$  ! ---------------------------------------------------------------------------
!!$  integer Ncav, Nins
!!$  double precision accuracy, accConv, dtau, init_tau, dynrange, accRomb
!!$  parameter (accRomb=1e-4)
!!$  ! ===========================================================================
!!$  ! Common statements for the grids. /grids1/ is for integers, /grids2/ for reals.
!!$  !                                                      [ZI'95-97; MN'99; MD'07]
!!$  ! =============================================================================
!!$  ! npY - max size for the radial (Y) grid
!!$  ! npP - max size for the impact parameter (P) and angular (mu) grids
!!$  ! npX - max size for the x-grid x(npX) in the optional disk calculation
!!$  ! npL - max size for the wavelength grid (as defined in 'lambda_grid.dat')
!!$  ! npG - max size for the number of grains in future MG version
!!$  ! npR - max size for the output inclination angle for slab (in *.i### files)
!!$  ! nY, nP, nL - the actual sizes of the grids: nY, nP are determined in
!!$  !         subroutines Ygrid and Pgrid, respectively, while nL = npL.
!!$  ! nPcav - the number of rays in the impact parameter grid, passing through
!!$  !            the inner cavity. They are calculated in subroutine Pgrid.
!!$  ! Yprev,nYprev - used in case of iterations over ETA (RDW) to keep Y and nY
!!$  !            from previous ETA iteration. Used in ChkFlux, Ygrid and Winds.
!!$  ! bOut(npP+2) - the impact parameter grid with 2 additional rays to take care
!!$  !            of pstar. Found in GetbOut, used in FindInt, Convolve and Visibili.
!!$  ! ----------------------------------------------------------------------------
!!$  integer  npY,npP,npX, npL, npG, npR
!!$  include '../userpar.inc'
!!$  !  parameter (npG=1)
!!$  integer nY, nYprev, nP, nPcav, nL
!!$  double precision Y(npY), Yprev(npY),P(npP),bOut(npP+2),
!!$  ! ===========================================================================
!!$  ! Common statement and variable definitions for array iYfirst(iP) which shows
!!$  ! the starting radial position for the line of sight with impact parameter P(iP),
!!$  ! and array YPequal(iP) which is 1 for P(iP).EQ.Y(iYfirst(iP)).
!!$  !                                          .                            [ZI,'95]
!!$  ! -----------------------------------------------------------------------------
!!$  integer iYfirst(npP), YPequal(npP), Plast(npY)
!!$  ! =============================================================================
!!$  ! Common statements for variables describing the dust density law ETA.
!!$  ! /dens1/ is for integers, /dens2/ - for character,  /dens3/ - for reals,
!!$  ! /dens4/ - for logical variables.
!!$  !                                                       [ZI'95-97; MN'99; MD'07]
!!$  ! ============================================================================
!!$  ! denstyp - flag for density type, entered as integer in 'fname.inp':
!!$  !           Flag values: 1-power law, 2-exponential density law,
!!$  !           3-RDW (radiatively driven winds), 4- analytical (gray-body) RDW
!!$  !           approximation, 5-ETA from a file. Internally, denstyp.eq.0 is for slab
!!$  !           (Private option for full RDW is denstyp.eq.6)
!!$  !!
!!$  !! FOR IMPROVEMENT: no need to carry both logical and denstyp in common.
!!$  !! denstyp will be local for Sub INPUT then only Logical flags will be carried.
!!$  !!
!!$  ! Ntr     - number transition points for broken power low
!!$  ! iterETA - counter over ETA iterations (in case of RDW, denstyp 5 or 6)
!!$  ! nYEtaf  - # pts. for ETA from file nameETA
!!$  ! pow     - real parameter for the chosen density law: power for 1, sigma for 2
!!$  !           [rho = exp(-(y/sigma)**2)], v1/ve for 3 (the ratio of expansion
!!$  !           velocities at the inner and outer radii).
!!$  ! ptr(10), Ytr(10) - powers and scaled transitional radii for denstyp=1
!!$  ! ETAcoef(npY,4) - the four coefficients for spline approximation of ETA
!!$  ! ETAdiscr(npY) - needed in ETA interpolation
!!$  ! yEtaf(1000),Etaf(1000) - the user supplied table for Eta(y) in file nameETA
!!$  ! Yout - the relative thickness, Yout=rout/r1.
!!$  ! -----------------------------------------------------------------------------
!!$  integer Ntr, iterETA, nYEtaf
!!$  character*235 nameETA
!!$  logical POWD, EXPD, RDW, RDWA, FILD, RDWPR, SLB, SPH
!!$  double precision pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10)
!!$  double precision ETAdiscr(npY), yEtaf(1000), Etaf(1000)
!!$  ! ===========================================================================
!!$  ! This file contains the common statement for variables describing the
!!$  ! spectrum of the illuminating source(s).            [ZI,'96, MN,'99,'00]
!!$  ! Former COMMON /source/
!!$  ! =============================================================================
!!$  !  Source index: 1- central source for sphere or left source for slab;
!!$  !                2- external source for sph or right source for slab;
!!$  !  Left, Right - flags for presence of sources
!!$  !
!!$  !  These are initialized in Input or InpStar:
!!$  !  typEntry(2) - type of boundary condition - Tsub(1) or flux in some form:
!!$  !                flux Fi, luminosity and distance, or Teff from the flux.
!!$  !  startyp(2) - type of spectral shape ,
!!$  !      nBB(2) - # of black bodies with temperatures Tbb(2,10) and relative
!!$  !                luminosities rellum(2,10);
!!$  !  Nlamtr(2)   - # transitional pts. for broken power law for input spectrum
!!$  !  lamtr(2,101)-transition wavelength; klam(3,100)- powers of the power law..
!!$  !
!!$  !  nameStar(2) - file with stellar spectrum
!!$  !  xSiO - depth in % of the SiO feature for Engelke-Marengo function
!!$  !  r1rs - the ratio of inner shell radius/stellar radius, found in Analysis
!!$  !
!!$  !  mu1, mu2 - the cosines of the left and right illumination angles in slab
!!$  !  ksi - the ratio of the right/left bol. fluxes (<1) for slab
!!$  !
!!$  !  chi - relative contribution of the external illumination of spherical shell
!!$  ! -----------------------------------------------------------------------------
!!$  integer startyp(2), Nlamtr(2), nBB(2), , Left, Right
!!$  character nameStar(2)*235
!!$  double precision  lamtr(2,101), klam(2,100), &
!!$       Tbb(2,10), rellum(2,10), Tstar(2), mu1, ksi, mu2, xSiO, &
!!$       r1rs, chi, dilutn
!!$  ! ===========================================================================
!!$  ! Commons for variables related to dust optical properties.
!!$  !                                                      [ZI'95-97; MN'99; MD'07]
!!$  ! Former COMMON /optprop/
!!$  ! =============================================================================
!!$  !  The following variables are initialized in Subroutine Input:
!!$  !  Nfiles -    number of user supplied .nk files
!!$  !  qsd,a1,a2 - parameters of the built-in size distribution functions
!!$  !  szds -      type of size distribution function
!!$  !  top -       opt.properties index ("top" means "type of properties")
!!$  !  xC(10), xCuser(10) - fractional number abundances for built-in and
!!$  !              user supplied components, respectively
!!$  !  Tsub(npG) - sublimation temperature for a given grain type
!!$  !
!!$  !  The following variables are determined in Subroutine GetTau:
!!$  !  iLfid -     the index of the fiducial wavelength lamfid, initilized in
!!$  !              subroutine GetTau
!!$  !  TAUfid -    opt. depth at lamfid
!!$  !  TAUmax -    max value of TAUtot(npL), used in Ygrid and ChkFlux
!!$  !  TAUtot(npL)- total opt.depth
!!$  !
!!$  !  These are determined in Subroutine GetOptPr
!!$  !  SigmaA(npG+1,npL), SigmaS(npG+1,npL) - absorp. and scatt. cross sections
!!$  !  SigExfid - the extinction cross section at lamfid
!!$  !  aveV -     for a single grain size it is the single grain volume,
!!$  !             otherwise it is the volume averaged over size distribution.
!!$  !  aveA -     for a single grain size it is the eff. grain area (Pi*a^2),
!!$  !             otherwise it is the effective area averaged over the size
!!$  !             distribution function
!!$  !
!!$  ! The following arrays carry an extra y-dependence as a preparation for the
!!$  ! mutigrain code, and are determined in subroutine GetOmega:
!!$  !  abund(npG,npY) - relative abundances of grain types, currently 1
!!$  !  omega(npG+1,npL) - scattering albedo
!!$  ! -----------------------------------------------------------------------------
!!$  integer iLfid, ifidG, szds, top, Nfiles,noprint
!!$  double precision TAUtot(npL),SigmaA(npG+1,npL), SigmaS(npG+1,npL),&
!!$       Tsub(npG), abund(npG,npY), TAUmax, xC(10), xCuser(10), &
!!$       SigExfid, TAUfid, taufid0,lamfid, qsd, a1, a2, aveV, aveA
!!$  ! ===========================================================================
!!$  ! Common statement for dynamical quantites (radiatively driven winds).
!!$  !                                                                [ZI'96; MN'99]
!!$  ! Former COMMON /dyn/
!!$  ! =============================================================================
!!$  ! Gamma(npY) - the ratio of gravitational to the radiative force on gas.
!!$  ! qF(npY) - effective flux averaged extinction efficiency
!!$  ! ugas(npY) - the gas velocity profile scaled by its terminal velocity
!!$  !             (at Yout)
!!$  ! vrat(npG,npY) - the ratio v(y)/vd(y), v is the gas velocity, vd is the
!!$  !                 dust velocity (different for each dust component)
!!$  ! I1,I2,I3,CMdot,Cve,CM,Cr1 - conversion constants for dynamics (found
!!$  !                in Analysis)
!!$  ! G1, Ginf, Prdw, delta, w1, Phi, PI - private quantities
!!$  ! -----------------------------------------------------------------------------
!!$  integer ver
!!$  double precision ugas(npY), qF(npY), vrat(npG,npY), Gamma(npY),  &
!!$       I1, I2, I3, CMdot, Cve, CM, Cr1, G1, Ginf, Prdw, gmax,      &
!!$       winf, Phi, PIrdw, QV, Qstar, zeta1, tauFdyn(npY)
!!$  ! =============================================================================
!!$  ! Related to dynamics calculation (RDW) [ZI, ME,2002]
!!$  ! -----------------------------------------------------------------------------
!!$  integer nW
!!$  double precision wav(npL), Eff(npL), f1, f2
!!$  ! =============================================================================
!!$  ! This file contains the common statement for the various solution arrays.
!!$  ! =============================================================================
!!$  !  Genrally, U is energy density, f is flux profile lambda*f_lambda, scaled by
!!$  !  the bolometri! flux. In array names 'de' stands for dust emission component,
!!$  !  'ds' for dust scattering, 's' is for star (or source).
!!$  !                                                     [ZI'95; MN'97,'00; MD'07]
!!$  !  Former COMMON /solution/ 
!!$  ! -----------------------------------------------------------------------------  
!!$  ! tauOut(npL) - total optical depth along a line of sight, used in sph_int
!!$  integer nYok, nPok, moment_loc,moment
!!$  double precision fde(npL,npY), fds(npL,npY), Utot(npL,npY), ftot(npL,npY), Td(npG,npY), &
!!$       Ubol(npY), fbol(npY), Spectrum(npL), SmC(30,99), tauF(npY), tr(npY), rg(npG,npY), &
!!$       Intens(npL,npP+2), IntOut(20,npP+2), tauZout(npP+2), tauOut(npL), Eps(npY), &
!!$       fsL(npL,npY), fsR(npL,npY), fsLbol(npY), fsRbol(npY),fsbol(npY), RPr(npY), Jext(npY), &
!!$       Ude(npL,npY), Uds(npL,npY), ETAzp(npP,npY), Ji, Jo, Psi, Psi0, RPr1
!!$  ! =============================================================================
!!$  ! Common statement for flags controling production of output files.
!!$  ! Here are also arrays related to intensity output from spherical shell.
!!$  !          .                                         [ZI'95; MN'97,'00; MD'07]
!!$  ! Former COMMON /output/
!!$  ! =============================================================================
!!$  !  Flags set in the master 'dusty.inp':
!!$  !   iVerb - for additional screen printout
!!$  !  Flags set by the user in each 'fname.inp':
!!$  !   iOUT  - used in Cllose, now is obsolete ;
!!$  !    see if iOut can be removed from the code.
!!$  !  iSPP  - for spectral properites file production
!!$  !  iA - for spectra; iB - for radial profiles; iC - for images,
!!$  !  iV - for visibility files; iPSF - for convolved images (private option);
!!$  !  iJ - for energy density profiles at user-selected radii
!!$  !  iX - for message files;
!!$  !  iInn - private option, set from inside Sub Input, it is for additional
!!$  !         printout in file 'fname.err'. Good for checking convergence.
!!$  !
!!$  !  Te_min - calculated in OccultMSG, to warn user for min required source Teff
!!$  !         (this is in case of BB type spectrum)
!!$  !!
!!$  !! see if Te_min can be removed.from common
!!$  !!
!!$  !  The following variables are used in imaging subroutines:
!!$  !  NlambdaOut - number of user required wavelengths (up to 20)
!!$  !  LambdaOut(20) - their wavelengths;  Visib(20,1000) - visibility results found
!!$  !  in subroutine Visibili; Offset(Nconv) - used for convolved images;
!!$  !  ConvInt(20,1000) - convolved intensity
!!$  !  nJOut - user-selected number of radii for J-output
!!$  !  YJOut(10) - the Y(iY) radii for J-output
!!$  ! -----------------------------------------------------------------------------
!!$  integer iVerb,iPhys,iA,iB,iC,iX,iInn,iPSF,iV,NlambdaOut,Nconv,Nvisi,iD,iJ,nJOut
!!$  character*100 zline(999)
!!$  double precision LambdaOut(20), ConvInt(20,1000), Visib(20,1000),     &
!!$       Offset(1000), qtheta1(1000), Te_min, YJOut(10), JOut(npL,10)
!!$  ! =============================================================================
!!$  !  Common flags for status.
!!$  ! =============================================================================
!!$  !  iWARNING - counter for number of warnings issued by the code
!!$  !  iERROR - counter for number of errors encountered by the code
!!$  !  iCUMM - counter for any warnings or errors. If 0, message issued that all OK.
!!$  !  Former COMMON /status/
!!$  ! -----------------------------------------------------------------------------
!!$  integer iWARNING, iERROR, iCUMM
!!$  ! =============================================================================
!!$  ! Common statement for variables used in imaging subroutines.       [ZI, '97]
!!$  ! Quantites describing the user supplied point spread function.
!!$  ! Added psfArea, used in normalization of convolved images [MN, Sep.'04]
!!$  !  Former COMMON /psf1/ and /psf2/
!!$  ! -----------------------------------------------------------------------------
!!$  integer psftype, Npsf, iLambda
!!$  double precision kPSF(20), FWHM1(20), FWHM2(20), Theta1, &
!!$       xpsf(1000), ypsf(1000), psfArea(20)
!!$  ! =============================================================================
!!$  ! Related to imaging                                                   [ZI,'97]
!!$  ! -----------------------------------------------------------------------------
!!$  integer ftype
!!$  double precision Ckn, Cxout, Cphi, Cqtheta
!!$  ! =============================================================================
!!$  ! Common statement for slab intensity.                            [MN, '99-00]
!!$  ! Used in Subroutines SLBRadT and SLBIntens and Real function Sexp(t).
!!$  ! =============================================================================
!!$  !  muobs - element of mu-array, needed in SLBIntens()
!!$  !  tauT = TAUslb(iL,nY) is the total opt.depth (lambda dependent)
!!$  !  Sfn - the local value of the source function               [MN, Jan'00]
!!$  ! -----------------------------------------------------------------------------
!!$  integer Nmu, transmit
!!$  double precision theta(npR), muobs, tauT, Sfn
!!$  ! =============================================================================
!!$  ! Common statement for slab calculation                           [MN, '99-00]
!!$  ! -----------------------------------------------------------------------------
!!$  
!!$  double precision TAUslb(npL,npY), fpbol(npY), fmbol(npY), fmed, &
!!$       SLBIntm(npR,npL), SLBIntp(npR,npL), IstR(npL),maxrat
!!$
!!$  !=============================================================================
end module common
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
!***********************************************************************
subroutine Input(nameIn,nameOut,tau1,tau2,GridType,Nmodel)
!======================================================================
! This subroutine reads input data from the file 'filename.inp'. It
! utilizes the function RDINP and subroutine RDINPS2 
! written by Moshe Elitzur.
!                                        [ZI,NOV'95; MN,JAN'00, MN'09]
!======================================================================
  use common
  implicit none
  INTERFACE
     subroutine Simpson(n,n1,n2,x,y,integral)
       integer n, n1, n2
       double precision integral
       double precision,allocatable ::  x(:), y(:)
     end subroutine Simpson
     subroutine getOptPr(nameQ,nameNK,er,stdf,top,szds,qsd,a1,a2,nFiles,xC,XCuser)
       integer er,top,szds,nFiles
       character*235,allocatable,nameQ(:)
       character*235 nameNK(10),stdf(7)
       double precision :: qsd,a1,a2,xC(10),xCuser(10)
     end subroutine getOptPr
     subroutine inp_rad(is,shp,spec_scale,styp)
       integer :: is,styp
       double precision :: spec_scale
       double precision,allocatable :: shp(:)
     end subroutine inp_rad
  END INTERFACE
  logical ::  Equal, noEqual, UCASE, composite
  character(len=235) :: stdf(7), str, nameIn, nameOut, nameEta,& 
       nameTau, strg, namepsf
  character(len=72) :: strpow,lamstr(nOutput)
  character(len=235),allocatable :: nameNK(:),nameQ(:)
  integer :: i, istop, GridType,Nmodel,L,top,iG,iFiles, &
       nFiles,szds, EtaOK, ang_type, imu, ioverflw
  double precision :: a,b,tau1,tau2,Lum,dist,RDINP,spec_scale, &
       dilutn,th1,th2,xC(10),xCuser(10),sum,qsd,a1,a2,&
       x1, ceta, Fi, Fo, psf1,Tinner_fidG
  double precision, allocatable ::aa(:),bb(:),xx(:),e(:) 
  !-----------------------------------------------------------------------
  UCASE = .true.
  Equal = .true.
  noEqual = .false.
  error = 0
  ! Open output file
  open(12,file=nameOut,status='unknown')
  write(12,*)'==========================='
  write(12,*)' Output from program dusty '
  write(12,*)' version: ',version
  write(12,*)'==========================='
  write(12,*)' '
  write(12,*)' Input parameters from file: '
  write(12,*) nameIn
  ! Open input file
  open(1,file=nameIn,status='old')
  call skip_header(1)
  !********************************************
  !** I. Geometry **
  !********************************************
  call rdinps2(Equal,1,12,str,L,UCASE)
  if(str(1:L).eq.'SPHERE') then
   slb = .false.
   sph = .true.
   sph_matrix = .false.
  elseif(str(1:L).eq.'SPHERE_MATRIX') then
   slb = .false.
   sph = .true.
   sph_matrix = .true.
  elseif(str(1:L).eq.'SLAB') then
   slb = .true.
   sph = .false.
   sph_matrix = .false.
  end if
  ! initialize sources (start without any source)
  left = 0
  right = 0
  Ji = 0.
  Jo = 0.
  ksi = 0.
  startyp(1) = 0
  startyp(2) = 0
  !********************************************
  !** II. Physical parameters **
  !********************************************
  ! (1) Flags for presence of sources for slab these have the meaning
  ! of "left" and "right" source,
  ! illumination
!!$  call rdinps2(Equal,1,12,str,L,UCASE)
!!$  if(str(1:L).eq.'ON') left = 1
!!$  if(str(1:L).eq.'OFF') left = 0
!!$  if(left.eq.0) then
!!$     if(slb) then
!!$        call msg(23)
!!$        left = 1
!!$     end if
!!$  end if
!!$  call rdinps2(Equal,1,12,str,L,UCASE)
!!$  if(str(1:L).eq.'ON') right = 1
!!$  if(str(1:L).eq.'OFF') right = 0
!!$  if (slb.and.(left.eq.0.and.right.eq.1)) then
!!$     left=1
!!$     right=0
!!$  endif
  !FOR SPHERE 
  if (sph) then
     call rdinps2(Equal,1,12,str,L,UCASE)
     if(str(1:L).eq.'ON') left = 1
     if(str(1:L).eq.'OFF') left = 0
     if (left.gt.0) then
        write(12,*) ' Central source spectrum described by'
        if (allocated(shpL)) deallocate(shpL)
        allocate(shpL(nL))
        shpL = 0
        call inp_rad(1,shpL,spec_scale,startyp(1))
        !typentry give the scale of input radiation
        call rdinps2(Equal,1,12,str,L,UCASE)
        if (str(1:L).eq.'FLUX') TypEntry(1) = 1
        if (str(1:L).eq.'LUM_R1') TypEntry(1) = 2
        if (str(1:L).eq.'ENERGY_DEN') TypEntry(1) = 3
        if (str(1:L).eq.'DILUTN_FAC') TypEntry(1) = 4
        if (str(1:L).eq.'T1') TypEntry(1) = 5
        !check if the entered value is acceptable
        if (TypEntry(1).lt.1.or.TypEntry(1).gt.5) then
           call msg(21)
           error = 1
        end if
        if (typentry(1).eq.1) then 
           Fi = RDINP(Equal,1,12)
           write(12,*) ' Flux at the inner boundary:', Fi,' W/m^2'
           Ji = Fi / (4.0d0 * pi)
        end if
        if (typentry(1).eq.2) then
           !enter luminosity [in Lo] of the source and distance r1[cm] to the source
           Lum = RDINP(Equal,1,12)
           dist = RDINP(Equal,1,12)
           !all units in dusty are in SI, so convert the input
           Lum = Lum*3.862d+26
           dist = dist/100.0d0
           write(12,*) ' Source luminosity ', Lum,' Lo and distance ', dist, ' m'
           Ji = Lum/(4.0d0*pi*dist*dist)/(4.0d0*pi)
        endif
        if (typentry(1).eq.3) then 
           Ji = RDINP(Equal,1,12)
        end if
        if (typentry(1).eq.4) then
           !entry of dilution (normalization) factor
           dilutn = RDINP(Equal,1,12)
           Ji = dilutn*spec_scale/(4*pi)
        endif
        if (typentry(1).eq.5) then
           !enter dust temperature on inner boundary, T1[K]
           Tinner_fidG = RDINP(Equal,1,12)
           if (right.eq.1) then 
              print*,' !!!!!Error!!!!!'
              print*,' Input of Temperature at the inner boundery not possible for two side illumination:'
              print*,' !!!!!Error!!!!!'
           end if
           write(12,*) ' Dust temperature on the inner boundary:', Tinner_fidG,' K'
        end if
     else
        typentry(1) = 0
     endif
     call rdinps2(Equal,1,12,str,L,UCASE)
     if(str(1:L).eq.'ON') right = 1
     if(str(1:L).eq.'OFF') right = 0
     if (right.gt.0) then
        write(12,*) ' External source spectrum described by'
        if (allocated(shpR)) deallocate(shpR)
        allocate(shpR(nL))
        shpR = 0
        call inp_rad(2,shpR,spec_scale,startyp(2))
        !typentry give the scale of input radiation
        if (left.eq.0) then
           call rdinps2(Equal,1,12,str,L,UCASE)
           if (str(1:L).eq.'ENERGY_DEN') typentry(2) = 3
           if (str(1:L).eq.'DILUTN_FAC') typentry(2) = 4
           !check if the entered value is acceptable
           if (typentry(2).lt.3.or.typentry(2).gt.4) then
              call msg(21)
              error = 1
              print*,'MSG(21)'
           end if
           if (typentry(2).eq.3) Jo = RDINP(Equal,1,12)
           if (typentry(2).eq.4) then
              !entry of dilution (normalization) factor
              dilutn = RDINP(Equal,1,12)
              Jo = dilutn*spec_scale/(pi)
           endif
           if (typentry(2).eq.5) then
              print*,' !!!!!Error!!!!!'
              print*,' Input of Temperature at the inner boundery not possible for two side illumination:'
              print*,' !!!!!Error!!!!!'
           end if
        else
           ! ksi is the relative bol.flux of the second source
           ksi = RDINP(Equal,1,12)
           if (ksi.lt.0.0) ksi = 0.0d0
           if (ksi.gt.1.0) ksi = 1.0d0
           Jo = ksi*Ji
           write(12,'(a49,F5.2)') ' Relativ bol.flux fraction of right source: R =',ksi
        end if
     end if
     if ((left.eq.1).and.(right.eq.1)) print*,"ERROR internal and external source is not supported in the spherical case!!!!!"
  end if
  !FOR SLAB
  if (slb) then
     call rdinps2(Equal,1,12,str,L,UCASE)
     if(str(1:L).eq.'ON') left = 1
     if(str(1:L).eq.'OFF') left = 0
     if (left.gt.0) then
        write(12,*) ' Left-side source spectrum described by'
        if (allocated(shpL)) deallocate(shpL)
        allocate(shpL(nL))
        shpL = 0
        call inp_rad(1,shpL,spec_scale,startyp(1))
        ! typentry give the scale of input radiation
        call rdinps2(Equal,1,12,str,L,UCASE)
        if (str(1:L).eq.'FLUX') typentry(1) = 1
        if (str(1:L).eq.'ENERGY_DEN') typentry(1) = 3
        if (str(1:L).eq.'DILUTN_FAC') typentry(1) = 4
        if (str(1:L).eq.'T1') typentry(1) = 5
        if ((typentry(1).lt.3.or.typentry(1).gt.5).and.(typentry(1).ne.1)) then
           call msg(21)
           error = 1
           print*,'error msg(21) see out file'
        end if
        if (typentry(1).eq.1) then 
           Fi = RDINP(Equal,1,12)
           write(12,*) ' Flux at the slab left boundary:', Fi,' W/m^2'
           Ji = Fi / (4.d0 * pi)
        end if
        if (typentry(1).eq.3) then 
           Ji = RDINP(Equal,1,12)
        end if
        if (typentry(1).eq.4) then
           ! entry of dilution (normalization) factor
           dilutn = RDINP(Equal,1,12)
           Ji = dilutn*spec_scale/pi
        endif
        if (typentry(1).eq.5) then
           ! enter dust temperature on inner boundary, T1[K]
           Tinner_fidG = RDINP(Equal,1,12)
           write(12,*) ' Dust temperature on the inner boundary:', Tinner_fidG,' K'
        end if
        write(12,'(a33)') ' Calculation in planar geometry:'
        !find the kind of illumination
        call rdinps2(Equal,1,12,str,L,UCASE)
        if (str(1:L).eq.'DIRECTIONAL') then
           write(12,'(a41)') ' Directional illumination from the left.'
           ! enter incident theta_in:
           ! th1 the left illumination angle (in degrees) measured from the normal
           th1 = RDINP(Equal,1,12)
           call chkangle(th1)
           th1 = th1*pi/180.0d0
           mu1 = dcos(th1)
        elseif (str(1:L).eq.'ISOTROPIC') then
           th1 = -1.0d0
           write(12,'(a40)') ' Isotropic illumination from the left.'
           mu1 = -1.0d0
        end if
     endif
     call rdinps2(Equal,1,12,str,L,UCASE)
     if(str(1:L).eq.'ON') right = 1
     if(str(1:L).eq.'OFF') right = 0
     if (right.gt.0) then
        write(12,*) ' Right-side source spectrum described by'
        if (allocated(shpR)) deallocate(shpR)
        allocate(shpR(nL))
        shpR = 0
        call inp_rad(2,shpR,spec_scale,startyp(2))
        call rdinps2(Equal,1,12,str,L,UCASE)
        if (str(1:L).eq.'DIRECTIONAL') then
           write(12,'(a41)') ' Directional illumination from the right.'
           ! enter incident theta_in:
           ! th1 the left illumination angle (in degrees) measured from the normal
           th2 = RDINP(Equal,1,12)
           call chkangle(th2)
           th2 = th2*pi/180.0d0
           mu2 = dcos(th2)
        elseif (str(1:L).eq.'ISOTROPIC') then
           th2 = -1.0d0
           write(12,'(a40)') ' Isotropic illumination from the right.'
           mu2 = -1.0d0
        end if
        ! ksi is the relative bol.flux of the second source
        ksi = RDINP(Equal,1,12)
        if (ksi.lt.0.0) ksi = 0.0d0
        if (ksi.gt.1.0) ksi = 1.0d0
        write(12,'(a49,F5.2)') ' Relative bol.flux fraction of right source: R =',ksi
        Jo = Ji*ksi
     endif
     ! Sab case isotropic ilumination:
     ! The input flux is the half flux of the sphere and therfore 
     ! F = pi*J -> J = F/pi instead of J = F/(4pi)
     if (typentry(1).eq.1) then
        if (th1.eq.-1.0d0) Ji = Ji*4
        if (th2.eq.-1.0d0) Jo = Jo*4
     endif
  endif
  write(12,*) ' --------------------------------------------'
  !=========  END READING OF SOURCE PARAMETERS ===================
  ! (2) DUST PROPERTIES
  ! # of different dust grains, to be used in a future version
  nG = 1
  ! 2.1 Chemical composition
  ! Type of optical properties
  call rdinps2(Equal,1,12,str,L,UCASE)
  ifidG = 1
  if (str(1:L).eq.'COMMON_GRAIN') then
     top = 1
     nG = 6
     composite = .false.
  elseif (str(1:L).eq.'COMMON_GRAIN_COMPOSITE') then
     top = 1
     composite = .true.
     nG = 1
  elseif (str(1:L).eq.'COMMON_AND_ADDL_GRAIN') then
     top = 2
     composite = .false.
     nG = 6
  elseif (str(1:L).eq.'COMMON_AND_ADDL_GRAIN_COMPOSITE') then
     top = 2
     composite = .true.
     nG = 1
  elseif (str(1:L).eq.'TABULATED') then
     top = 3
     nG = 1
  elseif (str(1:L).eq.'TABULATED_MULTI') then
     top = 3
     nG = RDINP(Equal,1,12)
  end if
  if (top.ne.1.and.top.ne.2.and.top.ne.3.and.top.ne.4) then
     call msg(9)
     error = 1
     print*,'error msg(9)'
  end if
  ! For top.lt.3 read in abundances for supported grains
  if (top.lt.3) then
     xC(1) = RDINP(Equal,1,12)
     if (xC(1).lt.0.0d0) xC(1) = 0.0d0
     if ((xC(1).eq.0.0d0).and.(nG.gt.1)) nG = nG - 1
     sum = xC(1)
     do iG = 2, 7
        ! Special care to be taken of graphite (1/3-2/3 rule):
        if (iG.ne.5) then
           xC(iG) = RDINP(noEqual,1,12)
           if (xC(iG).lt.0.0d0) xC(iG) = 0.0d0
           if ((xC(iG).eq.0.0d0).and.(composite .eqv. .false.)) then 
              nG = nG - 1
           endif
           ! i Equal 4 is data for graphite (parallel to c axis):
           if(iG.eq.4) xC(iG) = 1.0d0*xC(iG)/3.0d0
        else
           ! graphite (perpendicular to c axis) :
           xC(iG) = 2.0d0 * xC(iG-1)
        end if
        sum = sum + xC(iG)
     end do
  end if
  ! Assign supported dust filenames to stdf
  do iG = 1,7
     if (iG.eq.1) write(stdf(iG),'(a)')"data/stnd_dust_lib/OssOdef.nk"
     if (iG.eq.2) write(stdf(iG),'(a)')"data/stnd_dust_lib/OssOrich.nk"
     if (iG.eq.3) write(stdf(iG),'(a)')"data/stnd_dust_lib/sil-dlee.nk"
     if (iG.eq.4) write(stdf(iG),'(a)')"data/stnd_dust_lib/gra-par-draine.nk"
     if (iG.eq.5) write(stdf(iG),'(a)')"data/stnd_dust_lib/gra-perp-draine.nk"
     if (iG.eq.6) write(stdf(iG),'(a)')"data/stnd_dust_lib/amC-hann.nk"
     if (iG.eq.7) write(stdf(iG),'(a)')"data/stnd_dust_lib/SiC-peg.nk"
  enddo
  ! user supplied n and k:
  if (top.eq.2) then
     nFiles = RDINP(Equal,1,12)
     if (composite .eqv. .false.) nG = nfiles + nG
     ! File names
     allocate(nameNK(nFiles))
     do iFiles = 1, nFiles
        call filemsg(nameNK(iFiles),'optical constants:')
     end do
     if(error.ne.0) goto 996
     ! Abundances
     xCuser(1) = RDINP(Equal,1,12)
     if (xCuser(1).lt.0.0d0) xCuser(1) = 0.0d0
     sum = sum + xCuser(1)
     if (nfiles.gt.1) then
        do iFiles = 2, nfiles
           xCuser(iFiles) = RDINP(noEqual,1,12)
           if (xCuser(iFiles).lt.0.0d0) xCuser(iFiles) = 0.0d0
           sum = sum + xCuser(iFiles)
        end do
     end if
  end if
  if (top.lt.3) then
     if (sum.le.0.0d0) then
        call msg(5)
        error = 1
        print*,'msg(5)'
     end if
     ! Normalize abundances for supported grains:
     do iG = 1, 7
        xC(iG) = xC(iG) / sum
     end do
     ! Normalize abundances for user supplied grains
     if (top.eq.2) then
        do iFiles = 1, nfiles
           xCuser(iFiles) = xCuser(iFiles) / sum
        end do
     end if
  end if
  ! user supplied cross-sections:
  allocate(nameQ(nG))
  if (top.eq.3) then
     ! filename for qabs and qsca
     do iG = 1, nG
        call filemsg(nameQ(iG),'abs. and scatt. cross-sections:')
     end do
     if (nG.gt.1) then
        xC(1) = RDINP(Equal,1,12)
        do iG = 2, nG
           xC(iG) = RDINP(noEqual,1,12)
        end do
     else 
        xC(1) = 1.0
     end if
  end if
  ! 2.1 Grain size distribution
  if (top.ne.3) then
     ! Type of size distribution
     call rdinps2(Equal,1,12,str,L,UCASE)
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
        print*,'error msg(10)'
     end if
     ! Grain sizes
     if (szds.gt.1) then
        qsd = RDINP(Equal,1,12)
        a1 = RDINP(Equal,1,12)
        if (a1.le.0.0) a1 = 0.0001d0
        a2 = RDINP(Equal,1,12)
        if (szds.eq.2.and.a2.lt.a1) a2 = a1
     else
        qsd = 3.5d0
        a1 = 0.005d0
        a2 = 0.25d0
     end if
  end if
  ! 2.2 Temperatures
  if (allocated(Tsub)) deallocate(Tsub)
  allocate(Tsub(nG))
  Tsub = 0
  if (allocated(Tinner)) deallocate(Tinner)
  allocate(Tinner(nG))
  Tinner = 0
  Tsub(1) = RDINP(Equal,1,12)
  print'(A,I3,A,F12.3)',' Grain:',1,' Sublimation Temperature:',Tsub(1)
  if (nG.gt.1) then
     do iG = 2, nG
        Tsub(iG) = RDINP(noEqual,1,12)
        if (Tsub(iG).gt.Tsub(ifidG)) ifidG = iG
        print'(A,I3,A,F12.3)',' Grain:',iG,' Sublimation Temperature:',Tsub(iG)
     end do
  end if
!!$  if (Tinner_fidG.gt.Tsub(ifidG)) then 
!!$     Tinner_fidG = Tsub(ifidG)
!!$     write(12,*) ' *** Warning***'
!!$     write(12,*) ' Inner boundary temperature larger than sublimation temperature'
!!$     write(12,*) ' -> reduced to ',Tinner_fidG
!!$     write(12,*) ' *** Warning***'
!!$     write(6,*) ' *** Warning***'
!!$     write(6,*) ' Inner boundary temperature larger than sublimation temperature'
!!$     write(6,*) ' -> reduced to ',Tinner_fidG
!!$     write(6,*) ' *** Warning***'
!!$  end if
  if (typentry(1).eq.5) then 
     Tinner(ifidG) = Tinner_fidG
     print'(a,i3,a,f8.2)',' Inner Boundary Temperature of fiducial Grain(',ifidG,')=',Tinner(ifidG)
  end if
  if (allocated(SigmaA)) deallocate(SigmaA)
  allocate(SigmaA(nG+1,nL))
  SigmaA = 0
  if (allocated(SigmaS)) deallocate(SigmaS)
  allocate(SigmaS(nG+1,nL))
  SigmaS = 0
  call getOptPr(nameQ,nameNK,error,stdf,top,szds,qsd,a1,a2,nFiles,xC,XCuser)
  IF (iVerb.ge.2) print*,'Done with getOptPr'
  !=========  END READING DUST PROPERTIES ===================
  ! WriteOut prints all input data, read so far, in fname.out
  ! var1 is t1,fe1,luminosity or teff; var2 is r1; var3 is ext.rad. input
  ! call WriteOut(nameQ,nameNK,var1,var2,var3,a1,a2, dilutn,left,right)
  ! (3) Density distribution
  ! For sphere only:
  if(sph) then
     ! Parameter describing eta function:
     call rdinps2(Equal,1,12,str,L,UCASE)
     if (str(1:L).eq.'POWD') then
        denstyp = 1
     elseif (str(1:L).eq.'EXPD') then
        denstyp = 2
     elseif (str(1:L).eq.'RDW') then
        ! *** Winds ***
        ! denstyp.eq.3 is RDW with default values of v1/ve=0.2, GravCor=0.5
        ! denstyp.eq.6 is a private option with additional input for v1/ve and
        ! GravCor=max(Fgrav/Frad);
        denstyp = 3
     elseif (str(1:L).eq.'RDWA') then
        ! analytical (gray) approximation for rdw
        denstyp = 4
     elseif (str(1:L).eq.'USER_SUPPLIED') then
        ! file with user supplied density distribution
        denstyp = 5
     elseif (str(1:L).eq.'RDWPR') then
        ! private option for RDW with additional output
        denstyp = 6
     end if
     ! initialize EtaOK and Ntr
     EtaOK = 0
     Ntr = 0
     ! read parameters for each type of density distribution
     ! smooth or broken power laws
     if (denstyp.eq.1) then
        EtaOK = 1
        Ntr = RDINP(Equal,1,12)
        ! changed definition
        Ntr = Ntr - 1
        ! read in transition radii
        if (Ntr.gt.0) then
           if (allocated(Ytr)) deallocate(Ytr)
           allocate(Ytr(Ntr))
           Ytr = 0
           Ytr(1) = RDINP(Equal,1,12)
           if (Ntr.gt.1) then
              do i = 2, Ntr
                 Ytr(i) = RDINP(NoEqual,1,12)
              end do
           end if
           Yout = RDINP(noEqual,1,12)
        else
           ! for smooth density power law
           Yout = RDINP(Equal,1,12)
        end if
        if (Yout.le.1.0d0) Yout = 1.001d0
        ! read in powers
        pow = RDINP(Equal,1,12)
        if (Ntr.gt.0) then
           if (allocated(ptr)) deallocate(ptr)
           allocate(ptr(Ntr))
           ptr = 0
           do i = 1, Ntr
              ptr(i) = RDINP(NoEqual,1,12)
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
     if (denstyp.eq.2) then
        EtaOK = 1
        Yout = RDINP(Equal,1,12)
        if (Yout.le.1.0d0) Yout = 1.001d0
        pow = RDINP(Equal,1,12)
        if (pow.le.0.0d0) then
           EtaOK = 0
        else
           write(12,*)' density described by exponential distribution'
           write(12,'(a21,1p,e10.3)')'               sigma:',pow
           write(12,'(a21,1p,e10.3)')'  relative thickness:',Yout
        end if
     end if
     ! default approximation and default numerics for rad. driven winds
     if (denstyp.eq.3.or.denstyp.eq.4) then
        EtaOK = 1
        Yout = RDINP(Equal,1,12)
        if (Yout.le.1.0d0) Yout = 1.001d0
        ! ** default ** for epsilon = v1/ve = u1/ue:
        pow = 0.2d0
        if(denstyp.eq.3) then !RDW
           ! ** default ** for max(gravcor = fgrav/frad_press):
           ptr(1) = 0.5d0
           ! convergence criterion:
           ptr(2) = 1.0d0
           ! default linear version of the eq. for velocity
           ! ver = 1
        end if
        write(12,*)' Density for radiatively driven winds from'
        if (denstyp.eq.4) then !RDWA
           write(12,*)' Analytic approximation for gray dust.'
        else
           write(12,*)' Full dynamic calculation.'
        end if
        write(12,'(a21,1p,e10.3)')'  Relative thickness:',Yout
     end if
     ! full dynamical calculation for radiatively driven winds (private option)
     ! the user can specify parameters that have default values in denstyp=3
     ! user specified table for eta
     if(denstyp.eq.5) then
        EtaOK = 1
        call filemsg(nameeta,'Dust density distribution:')
        write(12,*)' Density distribution supplied from file:'
        write(12,'(2x,a100)') nameeta
        call prHeader(3,nameeta)
        ! read in the density
        open(26,err=997,file=nameeta,status='old')
        call skip_header(26)
        ! three lines in the header:
        do i = 1, 3
           read(26,*,err=997) strpow
        end do
        istop = 0
        i = 0
        do while (istop.ge.0) 
           read(26,*,iostat=istop) a, b
           i = i + 1
        end do
        nYetaf = i - 1
        if (allocated(xx)) deallocate(xx)
        allocate(xx(nYetaf))
        xx = 0
        if (allocated(aa)) deallocate(aa)
        allocate(aa(nYetaf))
        aa = 0
        if (allocated(bb)) deallocate(bb)
        allocate(bb(nYetaf))
        bb = 0
        if (allocated(e)) deallocate(e)
        allocate(e(nYetaf))
        e = 0
        if (allocated(yetaf)) deallocate(yetaf)
        allocate(yetaf(nYetaf))
        yetaf = 0
        if (allocated(etaf)) deallocate(etaf)
        allocate(etaf(nYetaf))
        etaf = 0
        rewind(26)
        ! # Read header again  --- **FH** need function to read header !!!!
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
900     close(26)
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
        call Simpson(nYetaf,1,nYetaf,yetaf,e,ceta)
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
  call rdinps2(Equal,1,12,str,L,UCASE)
  if (str(1:L).eq.'LINEAR') then
     GridType = 1
  elseif (str(1:L).eq.'LOGARITHMIC') then
     GridType = 2
  elseif (str(1:L).eq.'USER_SUPPLIED') then
     GridType = 3
  end if
  if (GridType.eq.3) then
     ! tau-grid from a file
     call filemsg(nametau,'user supplied tau-grid:')
     ! read optical depths
     open(27,err=992,file=nametau,status='old')
     call skip_header(27)
     ! fiducial wavelength
     ! (the second argument of rdinp is the unit)
     lamfid = RDINP(Equal,27,12)
     ! number of models in the list
     Nmodel = RDINP(Equal,27,12)
     do i = 1, Nmodel
        read(27,*) tauIn(i)
     end do
902  close(27)
     ! Sort the tau-grid if there is more than one model:
     if(Nmodel.gt.1) then
        call sort(tauIn,Nmodel)
     end if
     tau1 = tauIn(1)
     if (tau1.le.0.0d0) tau1 = 0.0001d0
     tau2 = tauIn(Nmodel)
  else
     ! fiducial wavelength
     lamfid = RDINP(Equal,1,12)
     ! total optical depths at lamfid
     TAU1 = RDINP(Equal,1,12)
     if (tau1.le.0.0d0) tau1 = 0.0001d0
     TAU2 = RDINP(Equal,1,12)
     if (tau2.le.tau1) then
        tau2 = tau1
        Nmodel = 1
     end if
     ! read number of models
     Nmodel = RDINP(Equal,1,12)
     ! Nrec = 1000, set in common
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
  !********************************************
  !** III. Numerical accuracy **
  !********************************************
  ! accuracy for convergence (typical 0.0001)
  ! accuracy for flux conservation
  accFlux = RDINP(Equal,1,12)
  accTemp = ((1.+accFlux)**(1./4.)-1.)*1e-1
!RN This converges better, but still not great.  accTemp = ((1.+accFlux)**(1./4.)-1.)*1.e-3
  if ((accFlux.le.0.0d0).or.(accFlux.gt.0.75)) then 
     print*,' Problem with specified Flux accuracy !!!'
  end if
  write(12,'(A,F12.2,A)') '   Required Flux accuracy:',accFlux*100.,'%'
  write(12,'(A,F12.2,A)') '   Required Temp accuracy:',accTemp*100.,'%' 
  write(12,*)' --------------------------------------------'
  !********************************************
  !** IV. Output flags **
  !********************************************
  ! Internal flag for additional miscellaneous output  [MN]:
  !  spectra
  iA = RDINP(Equal,1,12)
  iC = RDINP(Equal,1,12)
  ! images (intensity)
  if (iC.ne.0) then
     if (slb) then
        ! Read angular grid (this is theta_out) for slab intensity output.
        ! the output intensities are in units of lambda*I_lambda*cos(theta_out)/Fe
        ! where Fe=L/(4*pi*r^2), the local bolometric flux.
        ang_type = RDINP(Equal,1,12)
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
        NlambdaOut = RDINP(Equal,1,12)
        if (nLambdaOut.ge.1) then
           do i = 1, nLambdaOut
              LambdaOut(i) = RDINP(NoEqual,1,12)
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
           ! iPsf = rdinp(Equal,1,12)
           if (iPsf.ne.0) then
              Theta1 = RDINP(Equal,1,12)
              write(12,'(a39,1p,e7.1)') ' Convolved images produced for theta1=',theta1
              psftype = RDINP(Equal,1,12)
              if (psftype.ne.1.and.psftype.ne.2.and.psftype.ne.3) goto 994
              if (psftype.lt.3) then
                 ! Gaussians, read in parameters
                 ! FWHM for the first component
                 FWHM1(1) = RDINP(Equal,1,12)
                 if (nLambdaOut.gt.1) then
                    do i = 2, nLambdaOut
                       FWHM1(i) = RDINP(NoEqual,1,12)
                    end do
                 end if
                 if (psftype.eq.2) then
                    ! Relative strength for the second component
                    kPSF(1) = RDINP(Equal,1,12)
                    if (nLambdaOut.gt.1) then
                       do i = 2, nLambdaOut
                          kPSF(i) = RDINP(NoEqual,1,12)
                       end do
                    end if
                    !       FWHM for the second component
                    FWHM2(1) = RDINP(Equal,1,12)
                    if (nLambdaOut.gt.1) then
                       do i = 2, nLambdaOut
                          FWHM2(i) = RDINP(NoEqual,1,12)
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
                 call skip_header(28)
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
901              close(28)
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
        ! visibility not calculated in this version ...
        !iV = RDINP(Equal,1,12)
        !if(iV.ne.0) iV = abs(iC)
        ! end if for geometry
     end if
     write(12,*)' --------------------------------------------'
  else
     ! if iC=0 set the other flags to 0 (just in case).
     iPsf = 0
     iV = 0
     write(12,*)' --------------------------------------------'
  end if
  ! ---- added printout of lam*J_lam/J for sphere [MN'10] ------------
  if(SPH) then
     iJ = RDINP(Equal,1,12)
     if(iJ.GT.0) then
        nJOut = RDINP(Equal,1,12)
        if (nJOut.ge.1) then
           do i = 1, nJOut
              YJOut(i) = RDINP(NoEqual,1,12)
              !        make sure the radii are inside Dusty's range
              if (YJOut(i).le.1.0) YJOut(i) = 1.0
              if (YJOut(i).gt.Yout) YJOut(i) = Yout
           end do
        end if
        write(12,*)' En.density profile requested for these y:'
        write(12,'(a1,1p,10e12.3)')' ',(YJOut(i),i=1,nJOut)
        write(12,*)' --------------------------------------------'
     end if
  end if
  ! radial quantities
  iB = RDINP(Equal,1,12)
  ! run-time messages
  iX = RDINP(Equal,1,12)
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
    goto 996
998 write(12,*)' ***  FATAL ERROR IN DUSTY  ****'
    write(12,*)' Input file:'
    write(12,'(2x,a100)') nameIn
    write(12,*)' is missing?!'
    write(12,*)' *******************************'
    !  close(12)
    error = 3
    !-----------------------------------------------------------------------
996 close(1)
    close(66)
    return
end subroutine Input
!***********************************************************************

!***********************************************************************
subroutine inp_rad(is,shp,spec_scale,styp)
!=======================================================================
! This is the former SUBROUTINE InpStar(error,is,nameIn)
! This subroutine is for reading the input radiation parameters
!                                                              [MN,Mar'99]
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
  integer :: is,styp
  double precision :: spec_scale
  double precision,allocatable :: shp(:)
  !---local
  integer i,iL,iLs,nLs,k, l, nBB,ios1,filetype,nLamtr,kstop
  double precision sum, value, RDINP, xSiO, bb, x, planck
  double precision a,b, EMfunc,fplbol,Tsum
  double precision,allocatable :: Tbb(:),rellum(:),lambda_s(:),shp_s(:),&
       tmp_sort1(:),tmp_sort2(:),lamTr(:),klam(:),fl(:),fpl(:)
  character strg*40, str*235,filename*235,line*235
  logical Equal, noEqual, UCASE
  !-------------------------------------------------------------------------
  UCASE = .true.
  Equal = .true.
  noEqual = .false.
  error = 0

  ! Flag for the external spectrum
  call rdinps2(Equal,1,12,str,L,UCASE)
  ! help variable if input file is present
  filetype = 0
  if (str(1:L).eq.'BLACK_BODY') then 
     ! Number of black bodies
     styp=1
     nBB = RDINP(Equal,1,12)
     ! Stellar temperature(s)
     if (allocated(Tbb)) deallocate(Tbb)
     allocate(Tbb(nBB))
     Tbb = 0
     if (allocated(rellum)) deallocate(rellum)
     allocate(rellum(nBB))
     rellum = 0
     Tbb(1) = RDINP(Equal,1,12)
     if (Tbb(1).le.0.0d0) then
        call msg(8)
        error = 1
        print*,'MSG(8)'
     end if
     ! Single black body
     if (nbb.eq.1) then
        !relative luminosity
        rellum(1) = 1.0d0
        Tstar(is) = Tbb(1)
        write(12,'(A,F12.2,A)') '   a single Black Body with temperature', Tbb(1),' K'
     endif  !end if for one bb
     ! Multiple black bodies
     if (nbb.gt.1) then
        do i = 2, nbb
           Tbb(i) = RDINP(NoEqual,1,12)
           if (Tbb(i).le.0.0d0) then
              call msg(8)
              error = 1
              print*,'MSG(8)'
           end if
        end do
        ! Read in relative luminosities
        rellum(1) = RDINP(Equal,1,12)
        sum = rellum(1)
        do i = 2, nbb
           rellum(i) = RDINP(NoEqual,1,12)
           sum = sum + rellum(i)
        end do
        if (sum.le.0.0d0) then
           call msg(7)
           error = 1
           print*,'MSG(7)'
        end if
        ! Normalize
        do i = 1, nbb
           rellum(i) = rellum(i)/sum
           Tsum = Tsum + rellum(i)*Tbb(i)**4.0D+00                                                 
        end do
        Tsum= Tsum**0.25D+00                                                                          
        Tstar(is) = Tsum                                                                          
        write(12,'(a2,i2,a13)')'  ', nBB,' black bodies'
        write(12,'(a28)')'  with temperatures (in K):'
        write(12,'(2x,1p,10e10.3)')(Tbb(i),i=1,nBB)
        write(12,'(a43)')'  and relative luminosities, respectively:'
        write(12,'(1p,10e10.1)')(rellum(i),i=1,nBB)
     end if ! end if for multiple bb
     do iL = 1, nL
        bb = 0.0d0
        do k = 1, nbb
           x = 14400.0d0/(lambda(iL)*Tbb(k))  ! hc/k = 14400 micron.Kelvin
           bb = bb + rellum(k)*Planck(x)
        end do
        shp(iL) = bb
     end do
     call Simpson(nL,1,nL,lambda,shp,spec_scale)
     spec_scale = 0.0D0
     do k = 1, nbb
        spec_scale = spec_scale + rellum(k)*sigma*(Tbb(k)**4.0D0)
     end do
     deallocate(Tbb)
     deallocate(rellum)
  else if (str(1:L).eq.'ENGELKE_MARENGO') then
     styp=2
     ! Effective stellar temperature
     nBB = 1
     if (allocated(Tbb)) deallocate(Tbb)
     allocate(Tbb(nBB))
     Tbb = 0
     if (allocated(rellum)) deallocate(rellum)
     allocate(rellum(nBB))
     rellum = 0 
     Tbb(1) = RDINP(Equal,1,12)
     Tstar(is) = Tbb(1)
     ! Depth of SiO abs.feature in %
     xSiO = RDINP(Equal,1,12)
     if (xSiO.le.0.0d0) xSiO = 0.0001d0
     if (xSiO.gt.100.0d0) xSiO = 100.0d0
     do iL=1,nL
        shp = EMfunc(lambda(iL),Tbb(1),xSiO)
        write(12,*) '   ENGELKE_MARENGO function'
        write(12,'(a13,1p,e10.3,a16)')' with Teff =',Tbb(1), ' K and depth of'
        write(12,'(a30,F6.1,a2)')' the SiO absorption feature =', xSiO,' %'
     end do
     deallocate(Tbb)
     deallocate(rellum)
  else if (str(1:L).eq.'POWER_LAW') then
     styp=3
     ! Number of transitions
     nLamTr= RDINP(Equal,1,12)
     if (allocated(lamtr)) deallocate(lamtr)
     allocate(lamTr(nLamTr+1))
     lamTr = 0
     if (allocated(klam)) deallocate(klam)
     allocate(klam(nLamTr+1))
     klam = 0
     if (nLamTr.gt.0) then
        lamtr(1) = RDINP(Equal,1,12)
        if (nLamTr.gt.1) then
           do i = 2, nLamtr+1
              lamtr(i) = RDINP(NoEqual,1,12)
              if (lamtr(i).lt.lamtr(i-1)) then
                 call msg(6)
                 error = 1
                 print*,'MSG(6)'
              end if
           end do
        endif
        klam(1) = RDINP(Equal,1,12)
        if (nLamtr.gt.1) then
           do i = 2, nLamtr
              klam(i) = RDINP(NoEqual,1,12)
           end do
        end if
     else
        print*,'something wrong with powerlaw input'
     end if
     if (allocated(fl)) deallocate(fl)
     allocate(fl(nLamTr))
     fl = 0 
     if (allocated(fpl)) deallocate(fpl)
     allocate(fpl(nL))
     fpl = 0
     fl(1) = 1.0d0
     if (nLamtr.gt.1) then
        do i = 2, nLamtr
           fl(i) = fl(i-1)*(lamtr(i-1)/lamtr(i))**klam(i-1)
        end do
     end if
     do iL = 1, nL
        if ((lambda(iL)-lamtr(1))*(lambda(iL)-lamtr(nLamtr+1)).le.0.0d0) then
           kstop = 0
           k = 0
           do whiLe (kstop.eq.0)
              k = k + 1
              ! This is Matt's correction for reading more than one powers:
              if (lambda(iL).ge.lamtr(k).and. &
                   lambda(iL).le.lamtr(k+1)) then
                 kstop = 1
                 fpl(iL) = fl(k)*(lamtr(k)/lambda(iL))**klam(k)
                 fpl(iL) = fpl(iL)/lambda(iL)
              end if
           end do
        else
           fpl(iL) = 0.0d0
        end if
     end do
     call Simpson(nL,1,nL,lambda,fpl,fplbol)
     do iL = 1,nL
        shp(iL) = lambda(iL)*fpl(iL)/(fplbol)
     enddo
     if (Nlamtr.gt.0) then
        write(12,*) '   Power law with:'
        write(12,*)'    lambda      k'
        do i = 1, Nlamtr
           write(12,'(1x,1p,e10.3)')lamtr(i)
           write(12,'(11x,1p,e10.3)')klam(i)
        end do
        write(12,'(1x,1p,e10.3)')lamtr(Nlamtr+1)
     else
        write(12,*)' Input data for the source spectrum is not good.'
        write(12,*)' Changed to a 10000 K black body'
     end if
     spec_scale = 1
     deallocate(fl)
     deallocate(fpl)
     deallocate(lamTr)
     deallocate(klam)
  else if (str(1:L).eq.'FILE_LAMBDA_F_LAMBDA') then 
     write(12,*)'    Spectrum supplied from file (lambda_F_lambda):'
     styp=4
     filetype = 1
  else if (str(1:L).eq.'FILE_F_LAMBDA') then 
     write(12,*)'    Spectrum supplied from file (F_lambda):'
     styp=5
     filetype = 2
  else if (str(1:L).eq.'FILE_F_NU') then 
     write(12,*)'    Spectrum supplied from file (F_nu):'
     styp=6
     filetype = 3
  end if
  if (filetype.gt.0) then 
     strg = 'Spectral shape of external radiation:'
     call filemsg(filename,strg)
     write(12,*)'   ',filename
     open(3,file=filename,status='old')
     call skip_header(3)
     iLs=0
     ios1=0
     do while (ios1.ge.0)
        read(3,'(a235)',iostat=ios1) line
        iLs = iLs + 1
     end do
     ! 3 header lines and one line to much from loop
     nLs = iLs-1
     if (nLs.lt.2) then
        print*,'error while reading file',filename
     endif
     rewind(3)
     call skip_header(3)
     ios1 = 0
     iLs = 0
     allocate(lambda_s(nLs))
     lambda_s = 0
     allocate(shp_s(nLs))
     shp_s = 0
     allocate(tmp_sort1(nLs))
     tmp_sort1 = 0
     allocate(tmp_sort2(nLs))
     tmp_sort2 = 0
     do iLs = 1,nLs
        read(3,*,end=900,iostat=ios1) a, b
        if(ios1.ge.0) then
           lambda_s(iLs) = a
           ! shp_s is always f_lambda
           ! if filetype.eq.1 then file gives lambda*f_lambda
           if (filetype.eq.1) shp_s(iLs) = b/a
           ! if filetype.eq.2 then file gives f_lambda
           if (filetype.eq.2) shp_s(iLs) = b
           ! if filetype.eq.3 then file gives lnu=lambda**2*f_lambda
           if (filetype.eq.3) shp_s(iLs) = b/(a**2)
        end if
     end do
900  close(3)
     ! if input wavelengths in descending order turn them around
     if (lambda_s(1).gt.lambda_s(2)) then
        do iLs = 1, nLs
           tmp_sort1(iLs) = lambda_s(iLs)
           tmp_sort2(iLs) = shp_s(iLs)
        end do
        do iLs = 1, nLs
           lambda_s(iLs) = tmp_sort1(nLs+1-iLs)
           shp_s(iLs)    = tmp_sort2(nLs+1-iLS)
        end do
     end if
     ! interpolate to dusty grid
     do iL=1,nL
        call powerinter(nLs,nLs,lambda_s,shp_s,lambda(iL),iLs,shp(iL))
     end do

     call Simpson(nL,1,nL,lambda,shp,spec_scale)
     do iL=1,nL
        shp(iL) = shp(iL)/spec_scale*lambda(iL)
     end do

     error = 0
     deallocate(lambda_s)
     deallocate(shp_s)
     deallocate(tmp_sort1)
     deallocate(tmp_sort2)
  end if
  ! get total scale
  return
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
  
  integer ang_type, imu, length
  double precision th_min, th_max, angstep, cth_min, cth_max, caux, value, RDINP
  character*70 anggrid, strg,str*235
  logical Equal, noEqual,UCASE
  !----------------------------------------------------------------------
  Equal = .true.
  noEqual = .false.
  
  ! if ang_type=1 (equidistant in theta, given min,max,step)
  if (ang_type.eq.1) then
     th_min = RDINP(Equal,1,12)
     call chkangle(th_min)
     th_max = RDINP(Equal,1,12)
     call chkangle(th_max)
     if (th_max.le.th_min) then
        th_max = th_min
        nmu = 1
     end if
     ! step equidistant in theta
     angstep = RDINP(Equal,1,12)
     ! create the grid:
     imu = 1
     if (allocated(theta)) deallocate(theta)
     allocate(theta(int((th_max-th_min)/angstep)+1))
     theta = 0
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
     th_min = RDINP(Equal,1,12)
     call chkangle(th_min)
     th_max = RDINP(Equal,1,12)
     call chkangle(th_max)
     cth_min = dcos(th_min*pi/180.0d0)
     cth_max = dcos(th_max*pi/180.0d0)
     ! Step, equidistant in cos(theta)
     AngStep = RDINP(Equal,1,12)
     ! Create the grid:
     imu = 1
     if (allocated(theta)) deallocate(theta)
     allocate(theta(int((cth_max-cth_min)/AngStep)+1))
     theta = 0
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
     call skip_header(7)
     Nmu = RDINP(Equal,7,12)
     do imu = 1, Nmu
        read(7,*) theta(imu)
     end do
92   close(7)
  end if
  write(12,*)' Intensity requested for these theta_out(deg):'
  write(12,'(a1,8f7.1,/,x,8f7.1,/,x,10f7.1)')' ', (theta(imu), imu = 1, nmu)
  if (ang_type.eq.1) write(12,'(a34,f4.1)') '  equidistant in theta_out, step=', angstep
  if (ang_type.eq.2) write(12,'(a39,f4.1)') '  equidistant in cos(theta_out), step=', angstep
  if (ang_type.eq.3) write(12,'(a18,a70)') '  grid from file: ', anggrid
  return
end subroutine input_slb_ang
!***********************************************************************

!$
!***********************************************************************
subroutine getOmega(nY)
!=======================================================================
! This subroutine generates albedo omega(iL,iY) from the abs/sca cross-
! sections and the component abundancies. This is temporary (trivial)
! version  for single size grains.                     [Z.I., Mar. 1996]
!
!!** Note that Omega(iG,iL) here is re-defined compared to the old Dusty. [MN]
!=======================================================================
  use common
  implicit none
  integer  iG, iL, iY, nY
  double precision,allocatable :: ext(:),sca(:)
  !----------------------------------------------------------------------
  ! generate overall albedo through the envelope
  ! ** this is for future multigrain code **
  ! ** for single grains it is trivial **
  allocate(ext(nL))
  ext = 0
  allocate(sca(nL))
  sca = 0
  do iL = 1, nL
     ext(iL) = 0
     sca(iL) = 0
     do iG = 1, nG
        ! calculate albedo
        ext(iL) = ext(iL) + (sigmaA(iG,iL) + sigmaS(iG,iL))
        sca(iL) = sca(iL) + sigmaS(iG,iL)
        omega(iG,iL) = sigmaS(iG,iL)/(sigmaA(iG,iL) + sigmaS(iG,iL))
     end do
     omega(nG+1,iL) = sca(iL) / ext(iL)
  end do
  ! calculate relative abundances
  do iG = 1, nG
     do iY = 1, nY
        abund(iG,iY) = 1.0d0
     end do
  end do
  deallocate(ext)
  deallocate(sca)
  !--------------------------------------------------------------------
  return
end subroutine getOmega
!**********************************************************************

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
  !---------------------------------------------------------------------
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
  !--------------------------------------------------------------------
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
  !---------------------------------------------------------------------
  return
end function EMfunc
!***********************************************************************

!***********************************************************************
subroutine PrOut(nY,nP,nYprev,itereta,model,delta)
!=======================================================================
! This subroutine prints the results out.        [ZI,Feb'96; MN,Mar'99]
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
  integer :: model,nY,nP,nYprev,itereta
  double precision :: delta
  !---local variables
  integer :: i, j, iLV, iG, iL, iY, unt, imu, iOut, iNloc, tsub_reached
  double precision, allocatable::Elems(:,:),ftotL(:),ftotR(:),faux(:),sigma_tmp(:)
  double precision :: sigmaVs,sigmaVa,sigmaVe,Y_loc,J_loc,Jbol(10),temp_sub_reached
  double precision :: FbolL, FbolR, FbolIL,FbolIR,res, xAttTotL,&
       xAttTotR,xDsTotL,xDsTotR,xDeTotL,xDeTotR,temp1,temp2, &
       fnormL, fnormR, limval, tht1, dmax, GinfG1, xs, xde, xds, tr, &
       eta,psffunc(nOutput,1000),psfn
  character*120 STemp,Serr,hdint, hdcon,hdvis, s1, su1, s2, su2, tstr*10
  character*132 hdsp1,hdsp2,hdrslb1,hdrslb2,hdrsph1,hdrsph2,hdrdyn
  character*255 crossfilename
  CHARACTER*30 FMT
  external eta,psfn
  !---------------------------------------------------------------------
  allocate(ftotL(nL))
  ftotL = 0
  allocate(ftotR(nL))
  ftotR = 0
  allocate(faux(nL))
  faux = 0
  allocate(sigma_tmp(nL))
  sigma_tmp = 0
  if(allocated(Elems)) deallocate(Elems)
  if (nG.gt.1) allocate(Elems(nL,3+2*nG))
  if (nG.eq.1) allocate(Elems(nL,3))
  Elems = 0
  do iL=1,nL
     sigma_tmp(iL) = sigmaS(nG+1,iL)
  end do
  call lininter(nL,nL,lambda,sigma_tmp,lamfid,iLV,sigmaVs)
  do iL=1,nL
     sigma_tmp(iL) = sigmaA(nG+1,iL)
  end do
  call lininter(nL,nL,lambda,sigma_tmp,lamfid,iLV,sigmaVa)
  sigmaVe = sigmaVa+sigmaVs
  Elems(:,1) = lambda(:)
  Elems(:,2) = SigmaA(nG+1,:)/(sigmaVa+sigmaVs)
  Elems(:,3) = SigmaS(nG+1,:)/(sigmaVa+sigmaVs)
  if (model.eq.1) then 
     write(855,'(A,f7.2,A,e12.5)') '# total extinction at',lamfid,' micron <V>=',sigmaVa + sigmaVs
     hdrdyn = '#  lambda    <abs>/<V>  <sca>/<V>'
     if (nG.gt.1) then
        do iG=1,nG
           hdrdyn(34+(iG-1)*22:34+(iG-0)*22) = '  <abs>/<V>  <sca>/<V>'
           write(hdrdyn(38+(iG-1)*22:39+(iG-1)*22),'(i2.2)') iG
           write(hdrdyn(49+(iG-1)*22:50+(iG-1)*22),'(i2.2)') iG
           Elems(:,2+2*iG) = SigmaA(iG,:)/(sigmaVa+sigmaVs)
           Elems(:,3+2*iG) = SigmaS(iG,:)/(sigmaVa+sigmaVs)
        end do
        write(855,'(A)') hdrdyn(:34+nG*22)
        call maketable(Elems,nL,3+2*nG,855)
     else
        write(855,*) hdrdyn
        call maketable(Elems,nL,3,855)
     end if
     close(855)
  end if
  do iL=1,nL
     faux(iL) = fsL(iL,1)/lambda(iL)
  end do
  call Simpson(nL,1,nL,lambda,faux,FbolIL)
  do iL=1,nL
     faux(iL) = fsR(iL,nY)/lambda(iL)
  end do
  call Simpson(nL,1,nL,lambda,faux,FbolIR)
  FbolIL=FbolIL*Jext(1)
  FbolIR=FbolIR*Jext(nY)
  if(allocated(Elems)) deallocate(Elems)
  allocate(Elems(nL,8))
  Elems = 0
  !  find the bolometric fluxes at the boundaries [MN]
  !** FH changed to find ftot everywhere
  do iY = 1, nY
     do iL = 1, nL
        ! the emerging spectra for sphere (or right-side spectra for slab)
        ftot(iL,iY) = fsL(iL,iY) + fde(iL,iY) + fds(iL,iY) - ksi*fsR(iL,iY)
        if (abs(ftot(iL,iY)).lt.dynrange) ftot(iL,iY) = 0.
        faux(iL) = ftot(iL,nY)/lambda(iL)
     end do
  enddo
  do iL = 1,nL
     ftotL(iL) = fde(iL,1) + fds(iL,1) - ksi*fsR(iL,1)
     ftotR(iL) = fsL(iL,nY) + fde(iL,nY) + fds(iL,nY)
     if (abs(ftotL(iL)).lt.dynrange) ftotL(iL) = 0.
     if (abs(ftotR(iL)).lt.dynrange) ftotR(iL) = 0.
  end do
  do iL=1,nL
     faux(iL) = ftotL(iL)/lambda(iL)
  end do
  call Simpson(nL,1,nL,lambda,faux,temp1)
  do iL=1,nL
     faux(iL) = fsR(iL,1)/lambda(iL)
  end do
  call Simpson(nL,1,nL,lambda,faux,temp2)
  xAttTotL = abs(temp2/temp1)
  do iL=1,nL
     faux(iL) = fde(iL,1)/lambda(iL)
  end do
  call Simpson(nL,1,nL,lambda,faux,temp2)
  xDeTotL = abs(temp2/temp1)
  do iL=1,nL
     faux(iL) = fds(iL,1)/lambda(iL)
  end do
  call Simpson(nL,1,nL,lambda,faux,temp2)
  xDsTotL = abs(temp2/temp1)
  do iL=1,nL
     faux(iL) = ftotR(iL)/lambda(iL)
  end do
  call Simpson(nL,1,nL,lambda,faux,temp1)
  do iL=1,nL
     faux(iL) = fsL(iL,nY)/lambda(iL)
  end do
  call Simpson(nL,1,nL,lambda,faux,temp2)
  xAttTotR = abs(temp2/temp1)
  do iL=1,nL
     faux(iL) = fde(iL,nY)/lambda(iL)
  end do
  call Simpson(nL,1,nL,lambda,faux,temp2)
  xDeTotR = abs(temp2/temp1)
  do iL=1,nL
     faux(iL) = fds(iL,nY)/lambda(iL)
  end do
  call Simpson(nL,1,nL,lambda,faux,temp2)
  xDsTotR = abs(temp2/temp1)
  ! normalization factor for output spectra
  do iL=1,nL
     faux(iL) = ftotR(iL)/lambda(iL)
  end do
  call Simpson(nL,1,nL,lambda,faux,fnormR)
  do iL=1,nL
     faux(iL) = ftotL(iL)/lambda(iL)
  end do
  call Simpson(nL,1,nL,lambda,faux,fnormL)
  ! the emerging bolometric flux
  FbolR = fnormR * Jext(nY)
  if (slb) FbolL = fnormL * Jext(1)
  ! calculation of radiation pressure
  ! nG + 1 contains the sum of all sigma(iG) 1<=iG<=nG
  do iY=1,nY
     do iL=1,nL
        faux(iL) = (sigmaS(nG+1,iL)+sigmaA(nG+1,iL))*ftot(iL,iY)/lambda(iL)
     end do
     call Simpson(nL,1,nL,lambda,faux,temp1)
     RPr(iY) = temp1/(4*pi*clight*mprot*Gconst)*1.0D4*3.84e26/1.988e30*5.0D-26*Jext(iY)/sigmaVe
  enddo
  res = 0.0d00
  ! this is the cut-off for printout of small values (in spectra)
  limval = 1.0d-20
  ! Zeljko's calculation of theta1, the ang. size (in arcsec) of the cavity for Fbol=1e-6 W/m2
  tht1 = 412.6d0/(dsqrt(Ji*4*pi))
  ! error in %
  if (SmC(5,model).lt.0.1d0) then
     call getfs(SmC(5,model)*100.0d0,0,0,Serr)
  else if (SmC(5,model).ge.0.1d0.and.SmC(5,model).lt.1.0d0) then
     call getfs(SmC(5,model)*100.0d0,0,1,Serr)
  else
     call getfs(SmC(5,model)*100.0d0,0,2,Serr)
  end if
  !--------------  overall parameters to *.out file -----------------------
  ! write header to output file *.out
  if (model.eq.1) then
     write(12,*)'         '
     write(12,*)' RESULTS:'
     write(12,*)' --------'
     if (slb) then
        !    slab output
        s1=' ###    Tau0    Psi/Psi0     FiL      FiR       FbolL    FbolR     r1(cm)    TdL(K)    TdR(K)    RPr(1)   e(%)'
        su1=' ###      1        2          3        4          5        6         7          8         9       10      11'
        write(12,'(a)') s1
        write(12,'(a)') su1
        write(12,'(a)') &
             ' ============================================================================================================='
        !  output for sphere
     elseif(sph) then
        s1= ' ###   tau0   Psi/Psi0 Fi(W/m2)  r1(cm)   r1/rc    theta1   T1(K)    Td(K)    RPr(1)  e(%)'
        su1= ' ###     1       2        3        4        5        6        7        8        9      10'
        if((denstyp.eq.3).or.(denstyp.eq.4)) then ! 3(RDW) 4(RDWA)
           s2='  Mdot      Ve       M> '
           su2='   11       12       13 '
           write(12,'(a,a)') s1,s2
           write(12,'(a,a)')su1,su2
           write(12,'(a)') &
                ' =======================================================&
                =============================================================='
           ! **  private rdw file **
           if (denstyp.eq.6) then ! 6(RDWPR)
              s1= '###   tau0      tauF     Mdot      Ve       M>       '
              su1='###    1          2        3        4       5       6'
              s2= 'Ginf/G1   P    delta  d/sqrt(w1)  winf     Phi    zeta(1)'
              su2='        7        8        9        10       11       12'
              write(66,'(a53,a57)') s1,s2
              write(66,'(a53,a55)')su1,su2
           end if
        else
           write(12,'(a)') s1
           write(12,'(a)') su1
           write(12,'(a)') &
                ' ========================================================================================'
        end if
     end if !end if for sphere
  end if
  ! print output tables for ea.model
  !---------------- Output for slab: ---------------------------
  tsub_reached = 0
  do iY=1,nY
     do iG=1,nG
        if (Td(iG,iY).gt.(1+2*accTemp)*Tsub(1)) then 
           tsub_reached = 1
           if (temp_sub_reached.lt.Td(iG,iY)) then 
              temp_sub_reached = Td(iG,iY)
           end if
        end if
     end do
  end do
  if (tsub_reached.eq.1) then 
     if (iX.eq.1) then 
        write(18,*) ' ***Warning***'
        write(18,*) ' dust temperature is higher than sublimation temperature'
        write(18,*) ' ***Warning***'
     end if
     write(6,*) ' ***Warning***'
     write(6,*) ' dust temperature is higher than sublimation temperature'
     write(6,*) ' ***Warning***'
  end if
  if(slb) then
     if (tsub_reached.eq.0) then 
        write(12,'(i4,1p,10e10.2,a3)') model, taufid, Psi/Psi0,FbolIL, FbolIR, FbolL, FbolR, Cr1, Td(1,1), Td(1,nY), RPr(1), Serr
     else
        write(12,'(i4,1p,10e10.2,a3,a,1e11.3,a,1e11.3,a)') model, taufid, Psi/Psi0,FbolIL, FbolIR, FbolL, &
             FbolR, Cr1, Td(1,1), Td(1,nY), RPr(1), Serr,' #WARNING: dust temperature (',temp_sub_reached,&
             'K) is higher than sublimation temperature (',Tsub(1),'K)'
     endif
     !---------- for spherical shell ------------------------------
  elseif(sph) then
     if ((denstyp.eq.3).or.(denstyp.eq.4).or.(denstyp.eq.6)) then ! 3(RDW) 4(RDWA) 6(RDWPR)
        if (tsub_reached.eq.0) then 
           write(12,'(i4,1p,9e9.2,a1,a3,a1,1p,3e9.2)') &
                model, taufid, Psi/Psi0, Ji*4*pi, Cr1, r1rs, tht1, Td(1,1), Td(1,nY), RPr(1),' ',Serr,' ',CMdot, CVe, CM
        else
           write(12,'(i4,1p,9e9.2,a1,a3,a1,1p,3e9.2,a,1e11.3,a,1e11.3,a)') &
                model, taufid, Psi/Psi0, Ji*4*pi, Cr1, r1rs, tht1, Td(1,1), Td(1,nY), RPr(1),' ',Serr,' ',CMdot, CVe, CM,&
                ' #WARNING: dust temperature (',temp_sub_reached,'K) is higher than sublimation temperature (',Tsub(1),'K)'
        end if
     else
        if (tsub_reached.eq.0) then 
           write(12,'(i4,1p,9e9.2,a1,a3)') &
                model, taufid, Psi/Psi0, Ji*4*pi, Cr1, r1rs, tht1, Td(1,1), Td(1,nY), RPr(1),' ',Serr
        else
           write(12,'(i4,1p,9e9.2,a1,a3,a,1e11.3,a,1e11.3,a)') &
                model, taufid, Psi/Psi0, Ji*4*pi, Cr1, r1rs, tht1, Td(1,1), Td(1,nY), RPr(1),' ',Serr,&
                ' #WARNING: dust temperature (',temp_sub_reached,'K) is higher than sublimation temperature (',Tsub(1),'K)'
        endif
     end if
     if ((denstyp.eq.6)) then ! 6(RDWPR)
        !** private rdw file **
        if (model.eq.1) then
           write(66,'(a11,1p,e9.3,a10,1p,e9.3,a12,1p,e9.3,a13,e9.3)')  &
                '###   qv = ',qv,', Qstar = ',Qstar,', v1/vinf = ',pow, &
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
  !--------------   spectrum to *.s##  file   ------------------------
  if (iA.ne.0) then
     unt = 15
     call line(1,2,unt)
     if(slb) then
        write(unt,'(a7,i3,a8,f8.3,a36)')'# model',model,' taufid=',taufid,'  spectrum from the right slab side'
        write(unt,'(a13,1p,e9.2)') '# Fbol[W/m2]=',FbolR
     else
        write(unt,'(a7,i3,a8,f8.3,a10)') '# model',model,' taufid=',taufid,'  spectrum'
        write(unt,'(a13,1p,e9.2)') '# Fbol[W/m2]=',FbolR
     end if
     call line(1,1,unt)
     call getOmega(nY)
     do iL = 1, nL
        if ((ftot(iL,nY).ne.0.0d0).and.(ftotR(iL).ne.0.)) then
           xs = fsL(iL,nY)/ftotR(iL)
           xds = fds(iL,nY)/ftotR(iL)
           xde = fde(iL,nY)/ftotR(iL)
        else
           xs = 0.0d0
           xds = 0.0d0
           xde = 0.0d0
        end if
        !  no need to print negligible values
        if (dabs(xs).lt.limval) xs = 0.0d0
        if (dabs(xds).lt.limval) xds = 0.0d0
        if (dabs(xde).lt.limval) xde = 0.0d0
        if (dabs(fsL(iL,1)).lt.limval) fsL(iL,1) = 0.0d0
        !   Printing normalized spectral shapes. Bol. flux values are in the headers. [MN]
        if (dabs(ftot(iL,nY)).lt.limval) ftot(iL,nY) = 0.0d0
        Elems(iL,1) = lambda(iL)
        if (abs(fnormR).gt.dynrange) then  
           Elems(iL,2) = ftotR(iL)/fnormR
        else
           Elems(iL,2) = 0.0
        end if
        Elems(iL,3) = xs
        Elems(iL,4) = xds
        Elems(iL,5) = xde
        if (abs(fsLbol(1)).gt.dynrange) then  
           Elems(iL,6) = fsL(iL,1)/fsLbol(1)
        else 
           Elems(iL,6) = 0.0
        end if
        Elems(iL,7) = tautot(iL)
        Elems(iL,8) = omega(nG+1,iL)
     end do
     !------ tabulate the spectra in the desired form ----------
     if(slb) then
        hdsp1 = '#   lambda     fRight     xAtt       xDs        xDe        fInp_L     TauTot     albedo'
     else
        hdsp1 = '#   lambda     fTot       xAtt       xDs        xDe        fInp       TauTot     albedo'
     end if
     write(unt,'(A90)') hdsp1
     if(slb) then
        write(hdsp1,'(A,1p,E10.3,A,E10.3,A,E10.3,A,E10.3,A,E10.3,A)')  '    -1      ',FbolR,' ',xAttTotR,' ',&
             xDsTotR,' ',xDeTotR,' ',FbolIL,'     -1.       -1.'
     else
        write(hdsp1,'(A,E9.3,A,E9.3,A)')  '   -1.       ',FbolR,'                                   ',FbolIL,'      -1.         -1.'
     end if
     write(unt,'(A90)') hdsp1
     call maketable(Elems,nL,8,unt)
     !  spectra from the left (illuminated) slab side (file *.z##)
     if (slb) then
        call getOmega(nY)
        do iL = 1, nL
           if ((ftot(iL,1).ne.0.0d0).and.(ftotL(iL).ne.0.)) then
              xs =  fsR(iL,1)/ftotL(iL)
              xds = fds(iL,1)/ftotL(iL)
              xde = fde(iL,1)/ftotL(iL)
           else
              xs = 0.0d0
              xds = 0.0d0
              xde = 0.0d0
           end if
           if (dabs(xs).lt.limval) xs =0.0d0
           if (dabs(xds).lt.limval) xds =0.0d0
           if (dabs(xde).lt.limval) xde =0.0d0
           if (dabs(fsR(iL,nY)).lt.limval) fsR(iL,nY) = 0.0d0
           ! rescale ftot with the bolom flux for z-spectra
           if (dabs(ftot(iL,1)).lt.limval) ftot(iL,1) = 0.0d0
           Elems(iL,1) = lambda(iL)
           Elems(iL,2) = ftotL(iL)/fnormL
           Elems(iL,3) = xs
           Elems(iL,4) = xds
           Elems(iL,5) = xde
           if ((ksi.gt.0).and.(abs(fsRbol(nY)).gt.dynrange)) then
              Elems(iL,6) = ksi*fsR(iL,nY)/fsRbol(nY)
           else
              Elems(iL,6) = 0.0d0
           end if
           Elems(iL,7) = tautot(iL)
           Elems(iL,8) = omega(nG+1,iL)
        end do
        if (iA.eq.3) unt=25
        ! append to the .s## file or write in a separate .z## file (if iA=3)
        call line(1,1,unt)
        write(unt,'(a7,i3,a8,f8.3,a33)') '# model',model,' taufid=',taufid,' spectrum from the left slab side'
        write(unt,'(a13,1p,e9.2)') '# Fbol[W/m2]=',FbolL
        call line(1,1,unt)
        write(unt,'(a)')'#   lambda     fLeft      xAtt       xDs        xDe        fInp_R     TauTot     albedo'
        write(hdsp1,'(A,1p,E10.3,A,E10.3,A,E10.3,A,E10.3,A,E10.3,A)')  '    -1      ',&
             FbolL,' ',xAttTotL,' ',xDsTotL,' ',xDeTotL,' ',FbolIR,'     -1.       -1.'
        write(unt,'(A)') hdsp1
        call maketable(Elems,nL,8,unt)
     end if
  end if
  !-----------  radial quantities to *.r## (old *.bxx) file -------------
  if (iB.ne.0) then
     if(allocated(Elems)) deallocate(Elems)
     allocate(Elems(nY,9+nG))
     Elems = 0
     hdrslb1= '#     t        epsilon     tauF       RPr   '
     do iG = 1,nG
        write(hdrslb2(1+(iG-1)*11:1+(iG)*11),'(a7,i2.2,a2)'), '    Td(',iG,') '
     end do
     hdrsph1= '#     y         eta         t        tauF      epsilon      RPr'
     do iG = 1,nG
        write(hdrsph2(1+(iG-1)*11:1+(iG)*11),'(a7,i2.2,a2)'), '    Td(',iG,') '
     end do
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
           tr = TAUslb(iLfid,iY)/TAUslb(iLfid,nY)
           Elems(iY,1) = tr
           Elems(iY,2) = eps(iY)
           Elems(iY,3) = tauF(iY)
           if (abs(RPr(1)).gt.dynrange) then 
              Elems(iY,4) = RPr(iY)/RPr(1)
           else
              Elems(iY,4) = 0.0
           end if
           do iG=1,nG
              if (destroyed(iG,iY).gt.0) then 
                 Elems(iY,4+iG) = Td(iG,iY)
              else
                 Elems(iY,4+iG) = 0.0
              end if
           end do
        end do
        WRITE(FMT,'("(a43,a", I0, ")")') nG*11
        write(unt,FMT) hdrslb1,hdrslb2
        call maketable(Elems,nY,4+nG,unt)
        !------  for spherical shell --------
     elseif(sph) then
        do iY = 1, nY
           tr = ETAzp(1,iY)/ETAzp(1,nY)
           Elems(iY,1) = Y(iY)
           Elems(iY,2) = eta(Y(iY),nY,nYprev,itereta)
           Elems(iY,3) = tr
           Elems(iY,4) = tauF(iY)
           Elems(iY,5) = eps(iY)
           if (abs(RPr(1)).gt.dynrange) then 
              Elems(iY,6) = RPr(iY)/RPr(1)
           else 
              Elems(iY,6) = 0.0
           end if
           do iG=1,nG
              if (destroyed(iG,iY).gt.0) then 
                 Elems(iY,6+iG) = Td(iG,iY)
              else
                 Elems(iY,6+iG) = 0.0
              end if
           end do
           !     Elems(iY,8) = rg(1,iY)*Jext(iY)
           !     if (rdwpr) then
           ! redefine for private rdw (denstyp.eq.6) option
           !      Elems(iY,8) = rg(1,iY)
           !      Elems(iY,9) = gamma(iY)
           !      Elems(iY,10) = qF(iY)
           !     end if
        end do
        ! check values:
        do i = 1, 7
           do iY = 1, nY
              if(Elems(iY,i).lt.limval) Elems(iY,i) = 0.0d0
           end do
        end do
        ! with dynamics
        if (denstyp.eq.3) then ! 3(RDW)
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
           WRITE(FMT,'("(a65,a23,a", I0, ")")') nG*11
           write(unt,FMT) hdrsph1,hdrdyn,hdrsph2
           call maketable(Elems,nY,8+nG,unt)
        else
           WRITE(FMT,'("(a65,a", I0, ")")') nG*11
           write(unt,FMT) hdrsph1,hdrsph2
           call maketable(Elems,nY,6+nG,unt)
        end if
        ! end if for geometry
     end if
     ! end if for the iB (radial) flag
  end if
  !--------------   intensities to *.inn (old *.cxx) file  --------------
  if (abs(iC).ne.0) then
     ! slab intensity (found at the end of subroutine slbradt)
     ! theta(nmu) are the angles of output intensities
     if (slb) then
        if(allocated(Elems)) deallocate(Elems)
        allocate(Elems(nL,nmu+2))
        Elems = 0
        hdint = '#  lambda'
        unt = 17
        call line(1,2,unt)
        write(unt,'(a7,i3,a8,f8.3,a32)')'# model',model,' taufid=',taufid,' transmitted i(theta)*cos(theta)'
        call line(1,1,unt)
        do iL = 1, nL
           Elems(iL,1) = lambda(iL)
           do imu = 1, nmu
              !if(iPhys.eq.1) SLBintm(imu,iL) = SLBintm(imu,iL)*Jext(nY)
              !4pi comes from slbintp since it is divided by 4pi need to be changed!!
              Elems(iL,imu+1) = SLBintm(imu,iL)*Jext(1)*4*pi
           end do
           Elems(iL,nmu+2) = istR(iL)
        end do
        ! write(unt,'(a9,21f11.3)')hdint,(theta(imu),imu=1,nmu)
        ! printout angles in degrees
        ! write(unt,'(a9,37f11.1,a9)') hdint,
        ! &                    (theta(imu)*180.0d0/pi,imu=1,nmu),'     IstR'
        write(unt,'(a9,100f11.1)') hdint,(theta(imu)*180.0d0/pi,imu=1,nmu)
        call maketable(Elems,nL,nmu+1,unt)
        ! adding the column with stellar ints at the end of the table
        !  call maketable(Elems,nL,nmu+2,unt)
        hdint = '#  lambda'
        unt = 17
        call line(1,2,unt)
        write(unt,'(a7,i3,a8,f8.3,a32)')'# model',model,' taufid=',taufid,' reflected cos(theta)*i(theta)'
        call line(1,1,unt)
        do iL = 1, nL
           Elems(iL,1) = lambda(iL)
           do imu = 1, nmu
              !if(iPhys.eq.1) SLBintp(imu,iL) = SLBintp(imu,iL)*Jext(1)
              !4pi comes from slbintp since it is divided by 4pi need to be changed!!
              Elems(iL,imu+1) = SLBintp(imu,iL)*Jext(nY)*4*pi !4pi comes from slbintp
           end do
        end do
        !write(unt,'(a9,21f11.3)')hdint,(theta(imu),imu=1,nmu)
        !printout angles in degrees
        write(unt,'(a9,99f11.1)')hdint,(theta(imu)*180.0d0/pi,imu=1,nmu)
        call maketable(Elems,nL,nmu+1,unt)
     !------  for spherical shell --------
     elseif(sph) then
        if(allocated(Elems)) deallocate(Elems)
        allocate(Elems(np+2,nLambdaOut+2))
        Elems = 0
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
              !IF (iPhys.eq.1) THEN <--- iphys allways 1
              Elems(i,j+2) = 7.834d0*LambdaOut(j)*(Jext(nY)*4.0d0*pi*Yout**2.0d0)*Elems(i,j+2)
              ! ELSE
              !   Elems(i,j+2) = 7.834d0*LambdaOut(j)*(4.0d0*pi*Yout**2.0d0)*Elems(i,j+2)
              ! END IF
           end do
        end do
        write(unt,'(a21,20f11.2)')hdint,(LambdaOut(j),j=1,nLambdaOut)
        call maketable(Elems,nP+2,nLambdaOut+2,unt)
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
     Elems = 0
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
        Elems = 0
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
  !----------  energy density profiles to *jnn file  -------------------------
  if(sph) then
     if(allocated(Elems)) deallocate(Elems)
     allocate(Elems(nL,nJOut+1))
     Elems = 0
     if(iJ.gt.0) then
        unt = 19
        call line(1,2,unt)
        write(unt,'(a7,i3,a8,f8.3,a20)') '# model',model,' taufid=',taufid,'   Energy density   '
        call line(1,1,unt)
        !   Jext is found in the beginning of sub PrOut.
        !   Find the scale of en. density in [W/m2] for the required YJout(iOut) [MN]
        do iOut = 1, nJOut
           y_loc = YJOut(iOut)
           call LININTER(nY,nY,Y,Jext,y_loc,iNloc,J_loc)
           Jbol(iOut) = J_loc
        end do
        !   write the scale in the header
        write(unt,'(a11,1p,10e11.3)') 'Jbol[W/m2]=',(Jbol(iOut),iOut=1,nJout)
        write(unt,'(a11,10f11.2)')'      Y =  ',(YJOut(iOut),iOut=1,nJout)
        write(unt,'(a11)')'   lambda  '
        !   normalize en. density profiles
        do iOut = 1, nJOut
           do iL = 1, nL
              if(JOut(iL,iOut).lt.limval) JOut(iL,iOut) = 0.0d0
              faux(iL) = JOut(iL,iOut) / lambda(iL)
           end do
           call Simpson(nL,1,nL,lambda,faux,res)
           fnormR = res
           do iL = 1, nL
              Elems(iL,1) = lambda(iL)
              Elems(iL,iOut+1) = JOut(iL,iOut) / fnormR
           end do
        end do
        call maketable(Elems,nL,nJout+1,unt)
     end if
  end if
  if(allocated(Elems)) deallocate(Elems)
  !--------------------------------------------------------------------
  deallocate(ftotL)
  deallocate(ftotR)
  deallocate(faux)
  deallocate(sigma_tmp)
  return
end subroutine PrOut
!***********************************************************************

!***********************************************************************
subroutine CLLOSE(model,Nmodel)
  ! =======================================================================
  ! This subroutine closes output files.             [ZI,Feb'96; MN,Apr'99]
  ! =======================================================================
  use common
  implicit none
  character*235 su1, su2, s3, s4, txtt, txtf
  integer model, Nmodel, im
  ! ----------------------------------------------------------------------
  ! close the default output file:
  if (Error.ne.0) then
     write(12,'(a42,i4)') ' There are some error messages for model:',model
     write(12,*) ' Please check m## file (if not produced then rerun)'
     if (iverb.gt.0) print*, 'There are some error messages for model:',model,' (check m## file)'
  end if
  if (Warning.ne.0.and.Error.eq.0) then
     write(12,'(a36,i4)') ' There are some warnings for model:',model
     write(12,*)' Please check m## file (if not produced then rerun)'
     if (iverb.gt.0) print*, 'There are some warnings for model:',model,' (check m## file)'
  end if
  if (model.eq.Nmodel.or.error.eq.3.or.error.eq.4) then
     if (error.ne.3) then
        if((denstyp.eq.3).or.(denstyp.eq.4)) then
           write(12,'(a)')  ' ====================================='
        else
           if(sph) then
              write(12,'(a)')  ' ====================================='
           else
              write(12,'(a)')  ' ====================================='
           end if
        end if
        write(12,'(a22,1p,e8.1,a8)')'  (1) optical depth at',lamfid,' microns'
        write(12,'(a69,1p,e9.2)')'  (2) Psi as defined by eq.14 in IE97 with optically thin value Psi0=', Psi0
        if(slb) then
           !  ----------  for slab output ----------------------------
           write(12,'(a)')'  (3) input bol.flux (in W/m2) of the left-side source at the left slab boundary'
           write(12,'(a)')'  (4) input bol.flux (in W/m2) of the right-side source at the right slab boundary'
           write(12,'(a)')'  (5) bolometric flux (in W/m2) at the left slab boundary'
           write(12,'(a)')'  (6) bolometric flux (in W/m2) at the right slab boundary'
           write(12,'(a)')'  (7) position of the left slab boundary for L = 1e4 Lo'
           write(12,'(a)')'  (8) dust temperature at the left slab boundary'
           write(12,'(a)')'  (9) dust temperature at the right slab boundary'
           write(12,'(a)')'  (10) radiation pressure on left boundary; see manual for units'
           write(12,'(a)')'  (11) maximum error in flux conservation (Fmax-Fmin)/(Fmax+Fmin)'
        else
           !---------- for spherical shell ----------------------------
           write(12,'(a)')'  (3) bolometric flux at the inner radius '
           write(12,'(a)')'  (4) inner radius for L = 1e4 Lo'
           write(12,'(a)')'  (5) ratio of the inner to the stellar radius'
           write(12,'(a)')'  (6) angular size (in arcsec) when Fbol=1e-6 W/m2'
           write(12,'(a)')'  (7) dust temperature at the inner radius '
           write(12,'(a)')'  (8) dust temperature at the outer edge'
           write(12,'(a)')'  (9) radiation pressure on inner boundary; see manual for units'
           write(12,'(a)')' (10) maximum error in flux conservation (Fmax-Fmin)/(Fmax+Fmin)'
           if((denstyp.eq.3).or.(denstyp.eq.4)) then
              write(12,'(a)')' (11) mass-loss rate (in Mo/year)'
              write(12,'(a)')' (12) terminal outflow velocity (in km/s)'
              write(12,'(a)')' (13) upper limit of the stellar mass (in Mo)'
           end if
        end if
        write(12,*)' ================================================'
        if((error+warning).eq.0) write(12,*)' Everything is ok for all models.'
        !---- private file ---
!        if(rdwpr) write(12,*)' Table with the wind properties is in file *.rdw'
        !------------------------- spectra ------------------
        if(iA.ne.0) then
           if (iA.eq.1) then
              write(12,*)' All spectra are in file *.stb'
           else
              if (slb.and.iA.eq.3) then
                 write(12,*)' Spectra are in files *.s## and *.z##'
              else
                 write(12,*)' Spectra are in files *.s##'
              end if
           end if
        end if
        !------------------------- images ------------------
        if (iC.ne.0) then
           if (abs(iC).eq.1) then
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
              if (iV.ne.0.and.abs(ic).eq.3) write(12,*)' Visibility curves are in files *.v##'
           end if
           if (iC.eq.-3.and.iPsf.ne.0) then
              write(12,*)' Convolved images are in files *.c##'
              if (psftype.lt.3) write(12,*)' Point spread functions are in file *.psf'
           end if
        end if
        !----------------------- en.density ------------------
        if(sph) then
           if(iJ.ne.0) then
              if (iJ.eq.1) then
                 write(12,*)' All energy density profiles are in file *.jtb'
              else
                 write(12,*)' Energy density profiles are in files *.j##'
              end if
           end if
        end if
        !------------------------- radial ------------------
        if(iB.ne.0) then
           if (iB.eq.1) then
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
     end if  !end if for error.ne.3
  end if  !end if for models
  if (model.eq.Nmodel.or.error.eq.3.or.error.eq.4) then
     write(12,*) '========== the end =============================='
     close(12)
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
  ! conditionally close the energy density files
  if(iJ.eq.1) then
     if(model.eq.Nmodel) close(19)
  else
     close(19)
  end if
  if (iPsf.ne.0) close(21)
  if (iV.ne.0) close(22)
!-----------------------------------------------------------------------
  return
end subroutine CLLOSE
!***********************************************************************

!***********************************************************************
subroutine MakeTable(Elems,rows,cols,unt)
! =======================================================================
!     This is an auxiliary subroutine for print out of tables
!     of Elems(Nrows,Ncols) in output unit 'unt'. Nrows = max{nL,npY}.
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
subroutine fileMSG(fname,strg)
!=======================================================================
! Prints a message in *.out file in case of error opening the user
! supplied files.
!=======================================================================
  implicit none
  character aux*235, strg*(*), fname*(*)
  integer length, empty
  !---------------------------------------------------------------------
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
  close(12)
  !--------------------------------------------------------------------
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
  Error = 0
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
  end if
  !------------------------------------------------------------
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
  faux1 = 0
  allocate(faux2(nL))
  faux2 = 0
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
  fs = 0
  allocate(us(nL,npY))
  us = 0
  allocate(T4_ext(npY))
  T4_ext = 0
  allocate(emiss(nG,nL,npY))
  Emiss = 0
  allocate(fDebol(npY))
  fDebol = 0
  allocate(fDsbol(npY))
  fDsbol = 0
  allocate(fdsm(nL,npY))
  fdsm = 0
  allocate(fdsp(nL,npY))
  fdsp = 0
  allocate(fdem(nL,npY))
  fdem = 0
  allocate(fdep(nL,npY))
  fdep = 0
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
     if (iVerb.eq.2) print'(a,i4,a,i4,a,i4)', '  Calculating with nY=', nY,' nP=',nP,' nCav=',nCav
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
     else 
        initTemp=.false.
     end if
     if(iVerb.eq.2) write(*,'(A,F10.3,A)') '  Achieved error in bolometric Flux:',maxFerr*100,'%'
     call Flux_Consv(nY,nYprev,Ncav,itereta,iterfbol,fbol,fDebol,fDsbol,fbolOK,maxrat)
     if (iVerb.eq.2) write(*,'(a,i3,a,i3,a,i3)') '  After Flux_Cons nY=', nY
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
  tau = 0
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
  integer ipstar,nY,nP,nCav,nIns,nIns1,nIns2
  !---local
  integer i,ii,k,iP,iz,iw,nZ,j, Naux, istop, NinsLoc
  double precision delP, eta,tmp1,tmp2
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
     tmp1 = etadiscr(i)/etadiscr(i+1)
     tmp2 = etadiscr(i+1)/etadiscr(i)
     !Max 33% difference in gradient
     nIns1 = int((max(tmp1,tmp2)-1.)*3.)
     !Linear P grid with size close to Y - grid
     nIns2 = int((Y(i+1)-Y(i))/Y(nY)*1.5*nY) 
     nIns = max(nIns1,nIns2)
     if (Nins>9) Nins = 9
     if (Nins<1) Nins = 1
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
  emiss_total = 0
  allocate(ubol_old(npY))
  ubol_old = 0
  allocate(utot_old2(nL,npY))
  utot_old2 = 0
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
!  if ((iterfbol.eq.1).or.(initTemp.and.(iterfbol.gt.2))) then
!     call Init_Temp(nY,T4_ext,us)
!     if(iVerb.eq.2) write(*,*)' Done with initial dust temperature.'
!  end if
  call Init_Temp(nY,T4_ext,us)
  if(iVerb.eq.2) write(*,*)' Done with initial dust temperature.'
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
      if ((maxval(maxerrT).le.accTemp*9e-1).and.(maxval(maxerrU).lt.(accFlux*9.e-1))) conv = 1
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
  usL = 0
  allocate(usR(nL,nY))
  usR = 0
  allocate(m0(nL,nY))
  m0 = 0
  allocate(m1(nL,nY))
  m1 = 0
  allocate(m1p(nL,nY))
  m1p = 0
  allocate(m1m(nL,nY))
  m1m = 0

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
  qaux = 0
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
  term1 = 0
  allocate(term2(nL,nP)) 
  term2 = 0
  allocate(term_aux1(nP))
  term_aux1 = 0
  allocate(tauaux(nY)) 
  tauaux = 0
  allocate(z(nP,nY))
  z = 0
  allocate(angle(nP))
  angle = 0
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
           xg = 0
           allocate(wg(nn))
           wg = 0
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
           xg = 0
           allocate(wg(nn))
           wg = 0
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
  sph_em = 0
  allocate(sph_sc(nL,nY))
  sph_sc = 0
  allocate(tau(nY))
  tau = 0 
  allocate(Sfn_em(nY))
  Sfn_em = 0
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
  fnum = 0
  allocate(ff(nG))
  ff = 0
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
        if (Td(iG,iY).gt.(1.+2*accTemp)*Tsub(iG)) then 
           destroyed(iG,iY) = 1.0
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
  fnum = 0
  allocate(fdenum(nL))
  fdenum = 0
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
  fnum = 0
  allocate(fdenum(nL))
  fdenum = 0
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
  xN = 0
  allocate(yN(nP))
  yN = 0
  allocate(S_unscaled(nL,nY))
  S_unscaled = 0
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
  wSimp = 0
  allocate(wCorr(nP))
  wCorr = 0
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
  fnum = 0
  allocate(ff(nL))
  ff = 0
  
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
  ratio = 0
  allocate(tauaux(npY))
  tauaux = 0
  allocate(etatemp(npY))
  etatemp = 0
  allocate(Yins(npY))
  Yins = 0
  allocate(iYins(npY))
  iYins = 0
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
  DO iY = 1, nY-1
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


!!$  if (maxrat.gt.accFlux) then 
!!$     tmp2 = 0.0
!!$     do iY = 2, nY
!!$        tmp1 = dabs(fbolom(iY-1)-fbolom(iY))
!!$        tmp2 = tmp2 + tmp1
!!$        !if (tmp1.gt.tmp2) tmp2 = tmp1
!!$     end do
!!$     fact = -0.00001
!!$     n_ins = nY
!!$     do while (n_ins.gt.nY*0.3)
!!$        fact = fact+0.00001
!!$        n_ins = 0
!!$        do iY = 2, nY
!!$           tmp1 = dabs(fbolom(iY-1)-fbolom(iY))
!!$           if  (tmp1.ge.fact*tmp2) n_ins = n_ins + 1
!!$        end do
!!$     end do
!!$     !if # of added points is less than 5% of nY 
!!$     ! or less than 2 additinal points
!!$     if ((n_ins.lt.nY*0.05).or.(n_ins.lt.2)) fact = 0.0
!!$     do iY = 2, nY
!!$        tmp1 = dabs(fbolom(iY-1)-fbolom(iY))
!!$        if  ((TAUtot(1)*(ETAzp(1,iY)-ETAzp(1,iY-1)).GT.delTAUmax).or.&
!!$             (tmp1.ge.fact*tmp2).or.&
!!$             ((iterfbol.lt.3).and.&
!!$             ((((fbol_em(iY)/fbolom(iY)).gt.0.01).and.&
!!$             ((fbol_em(iY)/fbolom(iY)).lt.0.99)).or. &
!!$             (((fbol_sc(iY)/fbolom(iY)).gt.0.01).and.&
!!$             ((fbol_sc(iY)/fbolom(iY)).lt.0.99))))) then
!!$           !print*,iY,tmp1,tmp2,fact*tmp2
!!$           n_ins = 1
!!$           kins = kins + n_ins
!!$           do i_ins = 1,n_ins
!!$              Yins(kins-n_ins+i_ins) = Y(iY-1)+1.*i_ins/(1.*(n_ins+1))*(Y(iY)-Y(iY-1))
!!$              iYins(kins-n_ins+i_ins) = iY-1
!!$           end do
!!$        endif
!!$     enddo
!!$  endif

  if (maxrat.gt.accFlux) then 
     n_ins = 0
     fact=0.0
     do while (n_ins.lt.nY*0.5)
        fact = fact+devmax*0.00001
        n_ins = 0
        do iY = 2, nY
           tmp1 = dabs(fbolom(iY)-fmed)
           if  (tmp1.ge.(devmax-fact)) n_ins = n_ins + 1
        end do
     end do
     do iY = 2, nY
        tmp1 = dabs(fbolom(iY)-fmed)
        if (tmp1.ge.(devmax-fact)) then
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
  tau  = 0
  allocate(faux(nY))
  faux  = 0 
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
  flxsb = 0
  allocate(flxeb(nY))
  flxeb = 0
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
  qaux = 0
  allocate(qaux2(nL))
  qaux2 = 0
  allocate(alpha(1,nY))
  alpha = 0
  allocate(fnum(nL))
  fnum = 0 
  allocate(fdenum(nL))
  fdenum = 0
  allocate(Istell(nL))
  Istell = 0 
  allocate(tauOut(nL))
  tauOut = 0
  allocate(Ids(nL,nP))
  Ids = 0
  allocate(Ide(nL,nP))
  Ide = 0
  allocate(Istsc(nL,nP))
  Istsc = 0
  allocate(Istem(nL,nP))
  Istem = 0
  allocate(tzp(nP))
  tzp = 0
  allocate(Semis(nP))
  Semis = 0 
  allocate(Sscat(nP))
  Sscat = 0
  allocate(Sstem(nP))
  Sstem = 0
  allocate(Sstsc(nP))
  Sstsc = 0
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
  tau1 = 0 
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
  spectrum = 0
  allocate(qaux(nL))
  qaux = 0 
  allocate(qaux2(nL))
  qaux2 = 0 
  allocate(K1(nY))
  K1 = 0
  allocate(K2(nY))
  K2 = 0
  allocate(qpTd(nG,nY))
  qpTd = 0
  allocate(qpstar(nY))
  QpStar = 0
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
  goto 10
  ! ---------------------------------------------------------------------
end subroutine sort
!***********************************************************************

!***********************************************************************
double precision function Planck(x)
!=======================================================================
! This function evaluates the Planck function multiplied by wavelength
! and normalized by sigma*T^4/Pi.                      [Z.I., Mar. 1996]
! =======================================================================
  implicit none
  double precision x
  ! ---------------------------------------------------------------------
  if (x.gt.100.0d0) then
     Planck = 0.0d0
  else
     if (x.lt.0.00001d0) then
        Planck = 0.155d0*x**3.0d0
     else
        Planck = 0.155d0*x**4.0d0/(dexp(x) - 1.0d0)
     end if
  end if
  ! ---------------------------------------------------------------------
  return
end function Planck
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
  !---parameter
  integer n, n1, n2
  double precision integral
  double precision,allocatable ::  x(:), y(:)
  !---locale
  integer i
  double precision wgth, dyn2
  ! ---------------------------------------------------------------------
  dyn2 = 0.0d0
  ! set integral to 0 and accumulate result in the loop
  integral = 0.0d0
  ! calculate weight, wgth, and integrate in the same loop
  if ((n2-n1).gt.100) then
     !$OMP PARALLEL DO reduction(+:integral) private(i,wgth)
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
     !$OMP END PARALLEL DO
  else  if (n2.gt.n1) then
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
  ! --------------------------------------------------------------------
  return
end subroutine Simpson
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
  ! --------------------------------------------------------------------
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
  !---------------------------------------------------------------------

  return
end subroutine PowerInt
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
  !---------------------------------------------------------------------
  Scale = Y(1)
  do i = 1, N
     Y(i) = Y(i) / Scale
  end do
  !---------------------------------------------------------------------
  return
end subroutine scaleto1
!***********************************************************************

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
  ! ---------------------------------------------------------------------
  EPS_loc = 0.002d0
  h(1)=1.0d0
  do j=1,JMAX
     call trapzd(fnc,a,b,s(j),j)
     if (j.ge.K) then
        aux = 0.0d0
        call polint(h(j-KM),s(j-KM),K,aux,ss,dss)
        IF (dabs(dss).le.EPS_loc*dabs(ss)) RETURN
     endif
     s(j+1)=s(j)
     h(j+1)=0.25d0*h(j)
  end do
  ! --------------------------------------------------------------------
  RETURN
END subroutine ROMBY
!***********************************************************************
!!$
!!$
!!$!***********************************************************************
!!$SUBROUTINE ScaletoArea(Nmax,N,X,Y,Area)
!!$! =======================================================================
!!$! This subroutine scales a function Y(x) by the area A=Int{2*Pi Y(x)xdx}.
!!$! X and Y are 1D arrays. (Used for PSF normalization.)       [MN, Sep'04]
!!$! =======================================================================
!!$  IMPLICIT none
!!$  INTEGER Nmax, N, i
!!$  DOUBLE PRECISION Y(Nmax),X(Nmax),Fn(Nmax),Area,Pi
!!$  !-----------------------------------------------------------------------
!!$  Pi = 2.0D+00*ASIN(1.0)
!!$  DO i = 1, N
!!$     Fn(i) = Y(i)*X(i)
!!$  END DO
!!$  ! Integrate:
!!$  Area = 0.0D+00
!!$  DO i = 1, N-1
!!$     Area = Area + 0.5D+00*(Fn(i+1)+Fn(i))*(X(i+1)-X(i))
!!$  END DO
!!$  Area = 2.0D+00*Pi*Area
!!$  ! Normalize:
!!$  DO i = 1, N
!!$     Y(i) = Y(i) / Area
!!$  END DO
!!$  ! -----------------------------------------------------------------------
!!$  RETURN
!!$END SUBROUTINE ScaletoArea
!!$! ***********************************************************************
!!$
!!$!***********************************************************************
!!$subroutine shiftIns(x,Nmax,n,xins,i)
!!$!=======================================================================
!!$! Rearranges a vector X by inserting a new element Xins.    [MN, Aug'96]
!!$! =======================================================================
!!$  use interfaces
!!$  implicit none
!!$  !---parameter
!!$  integer :: Nmax,n,i
!!$  double precision :: xins
!!$  double precision,allocatable :: x(:)
!!$  !---local
!!$  integer j
!!$  ! ---------------------------------------------------------------------
!!$  do j = n+1, i+2, -1
!!$     x(j) = x(j-1)
!!$  end do
!!$  x(i+1) = xins
!!$  ! -----------------------------------------------------------------------
!!$  return
!!$end subroutine shiftIns
!***********************************************************************
!!$
!!$
!***********************************************************************
SUBROUTINE Spline(x,y,n,yp1,ypn,y2)
!=======================================================================
  INTEGER n,NMAX
  DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)
  PARAMETER (NMAX=500)
  INTEGER i,k
  DOUBLE PRECISION p,qn,sig,un,u(NMAX)
  ! --------------------------------------------------------------------
  if (yp1.gt..99e30) then
     y2(1)=0.
     u(1)=0.
  else
     y2(1)=-0.5
     u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  endif
  do i=2,n-1
     sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
     p=sig*y2(i-1)+2.
     y2(i)=(sig-1.)/p
     u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i)) &
          -(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  end do
  if (ypn.gt..99e30) then
     qn=0.
     un=0.
  else
     qn=0.5
     un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  endif
  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
  do k=n-1,1,-1
     y2(k)=y2(k)*y2(k+1)+u(k)
  end do
  ! --------------------------------------------------------------------
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
  ! find second derivative, secnder
  y2at1 = (fun(2)-fun(1))/(x(2)-x(1))
  y2atN = (fun(N)-fun(N-1))/(x(N)-x(N-1))
  CALL SPLINE(x,fun,N,y2at1,y2atN,secnder)
  ! generate coef(i,j), j=1,2,3,4
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
  ! ---------------------------------------------------------------------
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
  ! ----------------------------------------------------------------------
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
subroutine polint(xa,ya,n,x,y,dy)
  ! For polinomial interpolation, used in Subroutine Romby.
  ! ====================================================================
  implicit none
  integer n,Nmax
  double precision dy,x,y,xa(n),ya(n)
  parameter (Nmax=1000)
  integer i,m,ns
  double precision den,dif,dift,ho,hp,w,c(Nmax),d(Nmax)
  !---------------------------------------------------------------------
  c = 0.0d0
  d = 0.0d0
  ns=1
  dif=dabs(x-xa(1))
  do i=1,n
     dift=dabs(x-xa(i))
     if (dift.lt.dif) then
        ns=i
        dif=dift
     endif
     c(i)=ya(i)
     d(i)=ya(i)
  end do
  y=ya(ns)
  ns=ns-1
  do m=1,n-1
     do i=1,n-m
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
     end do
     if (2*ns.lt.n-m)then
        dy=c(ns+1)
     else
        dy=d(ns)
        ns=ns-1
     endif
     y=y+dy
  end do
!-----------------------------------------------------------------------
  return
end subroutine polint
!***********************************************************************
!!$
!!$
!!$
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
  !---------------------------------------------------------------------
  prd = yt(1)**(pt(1) - p0)
  if (j.gt.1) then
     do i = 2, j
        prd = prd * yt(i)**(pt(i) - pt(i-1))
     end do
  end if
  !---------------------------------------------------------------------
  return
end subroutine doProduct
!***********************************************************************
!!$
!!$
!***********************************************************************
subroutine gauleg(x1,x2,xg,wg,n)
!=====================================================================
  implicit none

  integer i,m,n,j
  double precision x1,x2,xm,xl,eps,delj,p,eta
  double precision xg(n),wg(n),sum,ff,p1,p2,p3,z1,z,pp
  parameter (eps=1.0d-14)
  ! -------------------------------------------------------------------
  xg = 0.0d0
  wg = 0.0d0
  m = int((n+1)/2)
  xm = 0.5d0*(x2+x1)
  xl = 0.5d0*(x2-x1)
  do i = 1, m
     z = cos(3.1415926535898d0*(dble(i) - 0.25d0)/(dble(n) + 0.5d0))
1    continue
     p1 = 1.0d0
     p2 = 0.0d0
     do j = 1, n
        p3 = p2
        p2 = p1
        p1 = ((2.0d0*dble(j)-1.0d0)*z*p2-(j-1.0d0)*p3)/dble(j)
     end do
     pp = dble(n)*(z*p1 - p2)/(z*z - 1.0d0)
     z1 = z
     z = z1-p1/pp
     if (abs(z-z1).gt.eps) go to 1
     xg(i) = xm-xl*z
     xg(n+1-i) = xm+xl*z
     wg(i) = 2.d0*xl/((1.d0-z*z)*pp*pp)
     wg(n+1-i) = wg(i)
  end do
  !-------------------------------------------------------------------
  return
end subroutine gauleg
!***********************************************************************

!**********************************************************************
SUBROUTINE ANALINT(nY,Nanal,xaux,yaux,m,aux)

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
  INTERFACE
     SUBROUTINE LINSYS(Nreal,A,B,X)
       integer Nreal
       DOUBLE PRECISION,allocatable :: A(:,:), B(:), X(:)
     END SUBROUTINE LINSYS
  END INTERFACE
  !---parameter
  integer nY, Nanal
  double precision :: xaux(Nanal),yaux(Nanal)
  double precision :: m,aux
  !---local
  INTEGER i, j
  DOUBLE PRECISION b
  double precision,allocatable :: A(:,:),xaux_tmp(:),yaux_tmp(:),coeff(:)
  ! ---------------------------------------------------------------------
  allocate(A(Nanal,Nanal))
  A = 0
  allocate(xaux_tmp(Nanal))
  xaux_tmp = 0
  allocate(yaux_tmp(Nanal))
  yaux_tmp = 0
  allocate(coeff(Nanal))
  coeff = 0
  do i=1,Nanal
     xaux_tmp(i) = xaux(i)
     yaux_tmp(i) = yaux(i)
  end do
  error = 0
  ! generate matrix A and vector B
  DO i = 1, Nanal
     DO j = 1, Nanal-1
        IF (xaux_tmp(i).EQ.0.0.AND.j.EQ.1) THEN
           A(i,j) = 1.0
        ELSE
           A(i,j) = xaux_tmp(i)**(1.0*j-1.0)
        END IF
     END DO
     A(i,Nanal) = 1.0/sqrt(1.0-xaux_tmp(i)*xaux_tmp(i))
  END DO
  ! solve for the coefficients
  CALL LINSYS(Nanal,A,yaux_tmp,coeff)
  IF(error.NE.0) THEN
     CALL MSG(19)
     print*,"MSG(19)"
     RETURN
  END IF
  ! upper limit for integration:
  b = xaux_tmp(Nanal)
  ! evaluate m-dependent contribution of the last term
  IF (m.GT.0.1) THEN
     IF (m.GT.1.1) THEN
        ! this is for m=2
        aux = 0.5*(DASIN(b)-b*sqrt(1.-b*b))
     ELSE
        ! this is for m=1
        aux = 1.0 - sqrt(1.-b*b)
     ENDIF
  ELSE
     ! this is for m=0
     aux = DASIN(b)
  ENDIF
  aux = aux * coeff(Nanal)
  ! add contribution from the polynom
  DO i = 1, Nanal-1
     aux = aux + coeff(i) * (b**(m+1.0*i)) / (m+1.0*i)
  END DO
! -----------------------------------------------------------------------
  do i=1,Nanal
     xaux(i) = xaux_tmp(i)
     yaux(i) = yaux_tmp(i)
  end do
999 deallocate(A)
  deallocate(xaux_tmp)
  deallocate(yaux_tmp)
  deallocate(coeff)
  RETURN
END SUBROUTINE ANALINT
!***********************************************************************

! ***********************************************************************
SUBROUTINE ChkConv(nY,accuracy_loc,Aold,Anew,Aconv_loc)
! =======================================================================
! This subroutine checks convergence of an array A(nY) between values
! given in Aold and Anew. If the relative difference for EVERY element
! is smaller than accuracy, Aconv is assigned 1, otherwise 0.
!                                                      [Z.I., Jul. 1996]
! =======================================================================
  use common
  IMPLICIT none
  ! --- parameter
  INTEGER Aconv_loc,nY
  DOUBLE PRECISION accuracy_loc 
  DOUBLE PRECISION,allocatable :: Aold(:), Anew(:) 
  ! --- local
  INTEGER iY
  DOUBLE PRECISION delta
! -----------------------------------------------------------------------
  Aconv_loc = 1
  ! loop over radial positions
  DO iY = 1, nY
     ! find relative difference
     delta = dabs(Anew(iY)-Aold(iY))
     IF (delta.GT.dabs(Anew(iY))*accuracy_loc) Aconv_loc = 0
  END DO
  ! -----------------------------------------------------------------------
  RETURN
END SUBROUTINE ChkConv
! ***********************************************************************

!***********************************************************************
SUBROUTINE LINSYS(Nreal,A,B,X)
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
  !--- parameter
  integer nY,Nreal
  DOUBLE PRECISION,allocatable :: A(:,:), B(:), X(:)
  !--- local
  INTEGER i, j
  integer ,allocatable :: indx(:)
  DOUBLE PRECISION d
  double precision,allocatable :: A1c(:,:), B1(:), A2c(:,:), B2(:)
  INTERFACE
     SUBROUTINE LUDCMP(A,N,NP,INDX,D)
       integer N,NP
       double precision :: D
       integer, allocatable :: indx(:)
       double precision,allocatable :: A(:,:),B(:)
     END SUBROUTINE LUDCMP
     SUBROUTINE LUBKSB(A,N,NP,INDX,B)
       integer N,NP
       integer, allocatable :: indx(:)
       double precision,allocatable :: A(:,:),B(:)
     END SUBROUTINE LUBKSB
     SUBROUTINE MPROVE(A,ALUD,N,NP,INDX,B,X)
       integer :: n,np
       integer,allocatable :: INDX(:)
       double precision,allocatable :: A(:,:),ALUD(:,:),B(:),X(:)
     END SUBROUTINE MPROVE
  END INTERFACE
  ! ---------------------------------------------------------------------
  allocate(indx(Nreal))
  indx = 0
  allocate(A1c(Nreal,Nreal))
  A1c = 0
  allocate(B1(Nreal))
  B1 = 0
  allocate(A2c(Nreal,Nreal))
  A2c = 0
  allocate(B2(Nreal))
  B2 = 0
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
  ! solve the system
  CALL LUDCMP(A1c,Nreal,nY,indx,d)
  IF (error.NE.0) RETURN
  CALL LUBKSB(A1c,Nreal,nY,indx,B1)
  ! improve the solution (saved in B)
  CALL MPROVE(A2c,A1c,Nreal,nY,indx,B2,B1)
  ! copy the improved solution to output vector X
  DO i = 1, Nreal
     X(i) = B1(i)
  END DO
  ! --------------------------------------------------------------------
  deallocate(indx)
  deallocate(A1c)
  deallocate(B1)
  deallocate(A2c)
  deallocate(B2)
  RETURN
END SUBROUTINE LINSYS
!***********************************************************************

! ***********************************************************************
SUBROUTINE LUBKSB(A,N,NP,INDX,B)
  IMPLICIT none
  !---parameter
  integer N,NP
  integer, allocatable :: indx(:)
  double precision,allocatable :: A(:,:),B(:)
  !---local
  integer :: i,j,ii,ll
  double precision :: sum
!!$  DIMENSION INDX(NP)
!!$  DOUBLE PRECISION A(NP,NP),B(NP)
  ! -------------------------------------------------------------------
  II=0
  !impossible to parallelize since B(J) needs to is changed and used!
  DO I=1,N
     LL=INDX(I)
     SUM=B(LL)
     B(LL)=B(I)
     IF (II.NE.0)THEN
        DO J=II,I-1
           SUM=SUM-A(I,J)*B(J)
        END DO
     ELSE IF (SUM.NE.0.) THEN
        II=I
     ENDIF
     B(I)=SUM
  END DO
  !impossible to parallelize since B(J) needs to is changed and used!
  DO I=N,1,-1
     SUM=B(I)
     IF(I.LT.N)THEN
        DO J=I+1,N
           SUM=SUM-A(I,J)*B(J)
        END DO
     ENDIF
     B(I)=SUM/A(I,I)
  END DO
  ! -------------------------------------------------------------------
  RETURN
END SUBROUTINE LUBKSB
! ***********************************************************************

! ***********************************************************************
SUBROUTINE LUDCMP(A,N,NP,INDX,D)
  IMPLICIT none
  !---parameter
  integer N,NP
  double precision :: D
  integer, allocatable :: indx(:)
  double precision,allocatable :: A(:,:),B(:)
  !---local
  integer :: i,k,j,imax,error
  DOUBLE PRECISION TINY
  PARAMETER (TINY=1.D-20)
  DOUBLE PRECISION,allocatable :: VV(:)
  double precision :: SUM,aamax,DUM
  ! ------------------------------------------------------------------
  allocate(VV(N))
  VV = 0
  error = 0
  D = 1.
  DO I = 1, N
     AAMAX=0.
     DO J = 1, N
        IF (DABS(A(I,J)).GT.AAMAX) AAMAX=DABS(A(I,J))
     END DO
     ! IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
     IF (AAMAX.EQ.0.) THEN
        error = 5
        RETURN
     ENDIF
     VV(I)=1./AAMAX
  END DO
!!  !$OMP PARALLEL DO private(J,I,SUM,DUM,IMAX,AAMAX,D)
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
!!  !$OMP END PARALLEL DO
  IF(A(N,N).EQ.0.)A(N,N)=TINY
  !------------------------------------------------------------------
  deallocate(VV)
  RETURN
END SUBROUTINE LUDCMP
! ***********************************************************************

!***********************************************************************
SUBROUTINE MPROVE(A,ALUD,N,NP,INDX,B,X)
  IMPLICIT none
  INTERFACE
     SUBROUTINE LUBKSB(A,N,NP,INDX,B)
       integer N,NP
       integer, allocatable :: indx(:)
       double precision,allocatable :: A(:,:),B(:)
     END SUBROUTINE LUBKSB
  END INTERFACE
  !---parameter
  integer :: n,np
  integer,allocatable :: INDX(:)
  double precision,allocatable :: A(:,:),ALUD(:,:),B(:),X(:)
  !---local
  integer i,j
  DOUBLE PRECISION SDP
  double precision,allocatable :: R(:)
  ! ---------------------------------------------------------------------
  allocate(R(N))
  R = 0
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
  deallocate(R)
  RETURN
END SUBROUTINE MPROVE
!***********************************************************************

! ***********************************************************************
SUBROUTINE Maple3(w,z,p,MpInt)
! =====================================================================
! This function calculates indefinite integral:
!    MpInt(iC) = INT(w^(2-iC) / sqrt(w^2-p^2) * dw), for iC=1,2,3,4.
!                                                     [Z.I., Apr. 1996]
! =====================================================================
  IMPLICIT none
  DOUBLE PRECISION w, z, p, MpInt(4)
  ! ---------------------------------------------------------------------
  ! integrals
  MpInt(1) = z
  MpInt(2) = dlog(w+z)
  IF (p.GT.0.0) THEN
     MpInt(3) = dacos(p/w)/p
     MpInt(4) = z/w/p/p
  ELSE
     MpInt(3) = -1.0 / w
     MpInt(4) = -0.5 / w / w
  END IF
  ! ---------------------------------------------------------------------
  RETURN
END SUBROUTINE Maple3
!***********************************************************************

!**********************************************************************
subroutine add(np1,nr1,np2,nr2,q1,q2,q3,qout) !only needed in matrix method
!======================================================================
! This subroutine evaluates the following expression:
! [qOut] = [q1] + [q2] + [q3]. qout, q1, q2 and q2 are matrices of
! physical size (np2,np1) and real size (nr2,nr1).     [Z.I., Nov. 1995]
! ======================================================================
  implicit none
  !---parameter
  integer np1, nr1, np2, nr2
  double precision, allocatable :: q1(:,:), q2(:,:), q3(:,:),qout(:,:)
  !---local
  integer  i2, i1
  !--------------------------------------------------------------------
  ! loop over index 2
  do i2 = 1, nr2
     ! loop over index 1
     do i1 = 1, nr1
        qout(i2,i1) = q1(i2,i1) +  q2(i2,i1) + q3(i2,i1)
     end do
  end do
  !--------------------------------------------------------------------
  return
end subroutine add
!**********************************************************************
!!$
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
     ! continue
     x = 0.0d0
  END IF
!-----------------------------------------------------------------------
  RETURN
END SUBROUTINE ChkRange
!***********************************************************************
!!$
!!$!***********************************************************************
!!$DOUBLE PRECISION FUNCTION Bessel(x)
!!$!=======================================================================
!!$! This function evaluates the Bessel function of the zeroth kind.
!!$! Formulae are from Abramowitz & Stegun.               [Z.I., Jan. 1997]
!!$! =======================================================================
!!$  use common
!!$  IMPLICIT none
!!$  INTEGER i
!!$  DOUBLE PRECISION x, c(6)
!!$  ! ---------------------------------------------------------------------
!!$  c(1) = -2.2499997D+00
!!$  c(2) =  1.2656208D+00
!!$  c(3) = -0.3163866D+00
!!$  c(4) =  0.0444479D+00
!!$  c(5) = -0.0039444D+00
!!$  c(6) =  0.00021D+00
!!$  Bessel=0.0D+00
!!$  IF (x.LE.3.0D+00)THEN
!!$     DO i=1,6
!!$        Bessel = Bessel + c(i)*(x/3.0D+00)**(2.0D+00*i)
!!$     END DO
!!$     Bessel = 1.0D+00 + Bessel
!!$  ELSE
!!$     Bessel = dsqrt(2.0D+00/Pi/x) * dcos(x-Pi/4.0D+00)
!!$  ENDIF
!!$  ! --------------------------------------------------------------------
!!$  RETURN
!!$END FUNCTION Bessel
!!$!***********************************************************************

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
  !---------------------------------------------------------------------
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
  !---------------------------------------------------------------------
  return
end subroutine LinInter
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
  call skip_header(28)
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
  if (msgno.eq.1) then
   write(18,*)' ************  WARNING  *************'
   write(18,*)' Temperature calculation in FindTemp'
   write(18,*)' achieved the limit of 500 iterations'
  end if
  if (msgno.eq.2) then
   write(18,*)' ************  WARNING  **************'
   write(18,*)' Energy density iterations in radtransf'
   write(18,*)' achieved the limit of 1000 iterations'
  end if
  if (msgno.eq.3) then
   write(12,*)' **********  INPUT ERROR ***********'
   write(12,*)' * Denstyp is not between 0 and 5! *'
   write(12,*)' * Check input file and try again.  *'
   write(12,*)' ***********************************'
  end if
  if (msgno.eq.4) then
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
  if (msgno.eq.14) then
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
  end if
  if (msgno.eq.17) then
   write(12,*)' *****************  WARNING  ********************'
   write(12,*)'  Eta is too steep and reaches values less than  '
   write(12,*)'  1e-12. Try decreasing the outer radius Yout.   '
   write(12,*)'  (see manual,3.3.3)                             '
   write(12,*)' ************************************************'
  end if
  if (msgno.eq.18) then
   write(12,*)' ************  WARNING  ************************ '
   write(12,*)'  The dynamical range of eta is more than 1e-12. '
   write(12,*)'  The outer radius Yout must be decreased so that'
   write(12,*)'  eta does not go below 1e-12 (see manual,3.3.3) '
   write(12,*)' *********************************************** '
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
subroutine OccltMSG(us)
!=======================================================================
! Prints a message informing the user about the min Teff required to
! neglect occultation by the central source.
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
  double precision,allocatable :: us(:,:)
  !---local
  character*10 tstrg
  integer iL
  double precision,allocatable :: qaux(:),qaux2(:)
  double precision :: res1, res2, Te_min_loc, mx, psitn, Planck, xp
  !----------------------------------------------------------------------
  allocate(qaux(nL))
  qaux = 0
  allocate(qaux2(nL))
  qaux2 = 0
  ! Estimate min Teff required to neglect occultation (eq.(5) in Manual):
  write(18,*) 'Tsub(',ifidG,')=', Tsub(ifidG)
  ! write(18,*) '  lambda(iL)  SigmaA(1,iL) Planck(xP)   xP'
  do iL = 1, nL
     qaux(iL) = sigmaA(ifidG,iL)*us(iL,1)/lambda(iL)
     xP = 14400.0d0/Tsub(ifidG)/lambda(iL)
     qaux2(iL) = sigmaA(1,iL)*Planck(xP)/lambda (iL)
  end do
  call Simpson(nL,1,nL,lambda,qaux,res1)
  call Simpson(nL,1,nL,lambda,qaux2,res2)
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
  deallocate(qaux)
  deallocate(qaux2)
  return
end subroutine OccltMSG
!***********************************************************************
integer function omp_get_thread_num()
  omp_get_thread_num = 0
  return 
end function omp_get_thread_num
!**********************************************************************
subroutine getOptPr(nameQ,nameNK,er,stdf,top,szds,qsd,a1,a2,nFiles,xC,XCuser)
!=====================================================================
! This subroutine calculates the absorption and scattering efficiences
! Qabs and Qsca in the wavelength range of the code or in case of
! user supplied efficiences reads them from a file.
!                                                 [ZI Mar96; MN Aug97]
!=====================================================================
  use common
  implicit none
  !---parameter
  integer er,top,szds,nFiles
  character*235,allocatable,nameQ(:)
  character*235 nameNK(10),stdf(7)
  double precision :: qsd,a1,a2,xC(10),xCuser(10)
  !---local
  character*235 fname, dummy
  integer iG, io1, iL, nLin, iiLaux,Nprop,nA, iiA, iiC,iCuser, Nmax, npA
  ! Nmax is the number of records in user supplied file with opt.prop.
  ! and npA is the dimension of the array of grain sizes
  parameter (Nmax=10000, npA=100)
  ! parameter (Nmax=10000, npA=1000)
  double precision,allocatable:: n(:),k(:),aQabs(:,:),aQsca(:,:), &
       sigAbs(:,:), sigSca(:,:),n_int(:), k_int(:)
  double precision aa,bb,cc,lambdain(Nmax),Qain(Nmax),Qsin(Nmax), &
       amax, nsd(npA), a(npA), faux1(npA), faux2(npA), f(npA), int,&
       ala(Nmax),  sizedist, aQa(Nmax), aQs(Nmax),  Cnorm, a3ave, a2ave
  !-----------------------------------------------------------------
  ! this should never change
  Nprop = 7
  !----------------------------------------------------------------
  allocate(n(nL))
  n = 0
  allocate(k(nL))
  k = 0
  allocate(aQabs(npA,nL))
  aQabs = 0
  allocate(aQsca(npA,nL))
  aQsca = 0
  allocate(SigAbs(npA,nL))
  SigAbs = 0
  allocate(SigSca(npA,nL))
  SigSca = 0
  allocate(n_int(nL))
  n_int = 0
  allocate(k_int(nL))
  k_int = 0
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
     if (1.eq.nA) then
        print*,"Single grain size setting cnorm to nsd(1)"
        cnorm = nsd(1)
     endif
     ! find the average grain volume aveV and average grain  eff.
     ! area aveA (needed in dynamics)
     if(dabs(a1-a2).le.1.d-3) then
        aveV = 4.0d0/3.0d0*pi*a1**3
        print*,'aveV:',aveV,'a1:',a1
        aveA = pi*a1**2.0d0
     else
        do iiA = 1, nA
           faux1(iiA)=nsd(iiA)*a(iiA)**3.0d0
        end do
        call powerint(npA,1,nA,a,faux1,a3ave)
        if (1.eq.nA) a3ave=faux1(1) !Single grain size
        aveV = 4.0d0/3.0d0*pi*a3ave/Cnorm
        print*,'aveV:',aveV,'a3ave/Cnorm:',a3ave/Cnorm
        do iiA = 1, nA
           faux1(iiA)=nsd(iiA)*a(iiA)**2.0d0
        end do
        call powerint(npA,1,nA,a,faux1,a2ave)
        if (1.eq.nA) a2ave=faux1(1) !Single grain size
        aveA = pi*a2ave/Cnorm
     end if
     !--  loop over supported components --
     do iiC= 1, Nprop
        f(iiC) = xC(iiC)
        fname = stdf(iiC)
        call getprop(nL,lambda,nL,fname,n,k,er)
        if (er.eq.3) goto 999
        ! calculate qabs and qsca for supported grains
        call mie(nL,nL,lambda,n,k,npA,nA,a,1,aQabs,aQsca)
        ! for each lambda integrate pi*a^2*qext with n(a)da
        do iL = 1, nL
           do iiA = 1, nA
              faux1(iiA)=nsd(iiA)*aQabs(iiA,iL)*pi*a(iiA)**2.0d0
              faux2(iiA)=nsd(iiA)*aQsca(iiA,iL)*pi*a(iiA)**2.0d0
           end do
           call powerint(npA,1,nA,a,faux1,int) 
           if (1.eq.nA) int=faux1(1) !Single grain size
           sigAbs(iiC,iL) = int/Cnorm
           call powerint(npA,1,nA,a,faux2,int) 
           if (1.eq.nA) int=faux2(1) !Single grain size
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
           call getprop(nL,lambda,nL,fname,n,k,er)
           if (er.eq.3) goto 999
           ! calculate qabs and qsca
           call mie(nL,nL,lambda,n,k,npA,nA,a,1,aQabs,aQsca)
           ! for each lambda integrate pi*a^2*qext with n(a)da
           do iL = 1, nL
              do iiA = 1, nA
                 faux1(iiA)=nsd(iiA)*aQabs(iiA,iL)*pi*a(iiA)**2.0d0
                 faux2(iiA)=nsd(iiA)*aQsca(iiA,iL)*pi*a(iiA)**2.0d0
              end do
              call powerint(npA,1,nA,a,faux1,int)
              if (1.eq.nA) int=faux1(1) !Single grain size
              sigAbs(iiC,iL) = int/Cnorm
              call powerint(npA,1,nA,a,faux2,int)
              if (1.eq.nA) int=faux2(1) !Single grain size
              sigSca(iiC,iL) = int/Cnorm
           end do
        end do
     else
        nfiles = 0
     end if
     ! mix them together and store in sigmaA(1,*) sigmaS(1,*)
     ! as well as storing them individually in sigmaA(iG+1,*) sigmaS(iG,*)
     do iL = 1, nL
        sigmaA(nG+1,iL) = 0.0d0
        sigmaS(nG+1,iL) = 0.0d0
        iG = 1
        do iiC= 1, Nprop+nfiles
           if ((f(iic).gt.0.0).and.(iic.ne.5)) then
              sigmaA(nG+1,iL) = sigmaA(nG+1,iL) + f(iiC) * sigAbs(iiC,iL)
              sigmaS(nG+1,iL) = sigmaS(nG+1,iL) + f(iiC) * sigSca(iiC,iL)
              sigmaA(iG,iL) = f(iiC) * sigAbs(iiC,iL)
              sigmaS(iG,iL) = f(iiC) * sigSca(iiC,iL)
              !parallel perpendicular graphites
              if (iic.eq.4) then
                 sigmaA(iG,iL) =  sigmaA(iG,iL) + f(iiC+1) * sigAbs(iiC+1,iL)
                 sigmaS(iG,iL) =  sigmaS(iG,iL) + f(iiC+1) * sigSca(iiC+1,iL)
                 sigmaA(nG+1,iL) = sigmaA(nG+1,iL) + f(iiC+1) * sigAbs(iiC+1,iL)
                 sigmaS(nG+1,iL) = sigmaS(nG+1,iL) + f(iiC+1) * sigSca(iiC+1,iL)
              endif
              if (nG.ne.1) iG = iG +1
           end if
        end do
     end do
     if (nG.eq.1) then
        sigmaA(1,:) = sigmaA(2,:)
        sigmaS(1,:) = sigmaS(2,:)
        print*,'COMPOSITE GRAIN'
     end if
     print*,'-------------------------------'
  else
     ! this is for top.ge.3 - [Sigma/V] from a file
     ! initialize aveV and aveA for this case
     aveV = 1.0d0
     aveA = 1.0d0
     ! read in lambda grid and optical properties
     do iL = 1, nL
        sigmaA(nG+1,iL) = 0.0d0
        sigmaS(nG+1,iL) = 0.0d0
     end do
     do iG = 1, nG
        open(123,err=998,file=nameQ(iG),status='old')
        call skip_header(123)
        read(123,'(a)',err=998)dummy
        read(123,'(a)',err=998)dummy
        read(123,'(a)',err=998)dummy
        iL = 0
        io1 = 0
        do while (io1.ge.0)
           read(123,*,end=900,err=998,iostat=io1) aa, bb, cc
           if (io1.ge.0) then
              iL = iL + 1
              lambdain(iL) = aa
              Qain(iL) = bb
              Qsin(iL) = cc
           end if
        end do
900     close(123)
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
           sigmaA(iG,iL) = xC(iG) * aa
           ! call lininter(Nmax,nLin,lambdain,Qsin,lambda(iL),iiLaux,aa)
           call powerinter(Nmax,nLin,lambdain,Qsin,lambda(iL),iiLaux,aa)
           sigmaS(iG,iL) = xC(iG) * aa
           sigmaA(nG+1,iL) = sigmaA(nG+1,iL) + sigmaA(iG,iL)
           sigmaS(nG+1,iL) = sigmaS(nG+1,iL) + sigmaS(iG,iL)
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
999 deallocate(n)
  deallocate(k)
  deallocate(aQabs)
  deallocate(aQsca)
  deallocate(SigAbs)
  deallocate(SigSca)
  deallocate(n_int)
  deallocate(k_int)
  return
end subroutine getOptPr
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
  !---------------------------------------------------------------------
  return
end subroutine getSizes
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
  call skip_header(2)
  ! read in a header from the input file
  !do i = 1, 7
  !   read(2,'(a)',err=998)line
  !end do
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
  error = 3
!-----------------------------------------------------------------------
999 return
end subroutine GetProp
!***********************************************************************

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
  !---------------------------------------------------------------------

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




!****************************************************************************!
!                                                                           !
! Replacement of the old RDINP. Simpler operation -- the number is          !
! obtained straight by reading from internal unit instead of the old        !
! method which required parsing and evaluating. Also, introduces            !
! the input variable outUnit, the unit number to which this                 !
! function will write output. In the old version this was hard coded        !
! as unit 12, which was really a bug.                                       !
!                                                    [M.E., March 2012]     !
!                                                                           !
!===========================================================================!
double precision function RDINP(Equal, inUnit, outUnit)
! ==========================================================================!
!  Read lines, up to 255 long, from pre-opened unit inUnit and extract      !
!  all input numbers from them. When EQUAL is set, numeric input data       !
!  must be preceded by an equal sign. All non-numeric data and numbers      !
!  not preceded by = are ignored when EQUAL is on.                          !
!                                                                           !
!  RDINP = next number encountered (after equal sign) and terminated blank. !
!  The blank can be optionally preceded by comma.                           !
!  Input numbers can be written in any FORTRAN allowable format.            !
!  In addition, comma separation is allowed in the input numbers.           !
!                                                                           !
!  All text after % (or !) is ignored, as in TeX (or F90).                  !
!  Lines with * in the first column are echoed to output on pre-opened      !
!  unit outUnit.                                                            !
!                                                                           !
!  The search is conducted between FIRST, which is                          !
!  continuously increased, and LAST.  A new line is read when FIRST         !
!  exceeds LAST, and can be forced by calling with -iUnit.                  !
!                                                                           !
!===========================================================================!
   IMPLICIT None
   Integer :: inUnit, outUnit, ind, First = 1, Last = 0
   Logical, intent(in) :: Equal
   CHARACTER(255) Card, no_comma
   Save Card, First, Last
! -----------------------------------------------------------------------
!
   IF (inUnit.lt.0) Then                         ! force a new line
      First = Last + 1
      inUnit = -inUnit
   END IF

   DO
      If (first > last) then                     ! Time to get a new line
         READ (inUnit, '(A)' , END = 99) Card
         if (len_trim(Card) == 0) cycle          ! ignore empty lines
                                                 ! Echo to output lines that start with * 
         IF (Card(1:1) == '*') &
             WRITE (outUnit,'(A)') TRIM(Card(1:len_trim(card)))
         Card = trim(no_comma(Card))             ! remove commas to allow comma-separated numbers
         first = 1
         last = len_Trim(Card)
         ind = MAX(Index(Card,'%'), Index(Card,'!'))
         if (ind.gt.0) last = ind - 1            ! Everything after % and ! is ignored 
      End If

     !Get past the next '=' when the EQUAL flag is set
      If (Equal) then
        DO WHILE ((Card(first:first) /= '=') .and. (first <= last) )
          first = first + 1
        END DO
      End If
      first = first + 1
      IF (first > last) cycle

     !Find start of next number; necessary when the EQUAL flag is off
      Do While (.not. (Card(first:first).ge.'0'.AND.Card(first:first).le.'9') &
                .and. (first <= last))
          first = first + 1
      End Do
      if (first > last) cycle

     !OK; time to get the number
      READ(card(first - 1:last), *, ERR = 98) RDINP 
      
     !and move past its end
      Do While (Card(first:first) /= ' ')
          first = first + 1
      End Do
      
      return
   End DO

98 WRITE (outUnit,'(3(1x,a,/))')                               &
   ' ****ERROR. RDINP could not read a proper number from', &
   Card(first - 1:last),																			 &
   ' ****Number should be preceded and terminated by spaces'
   RETURN
99 WRITE (outUnit,'(3(1x,a,/))')                                      &
   ' ****TERMINATED. EOF reached by RDINP while looking for input. ', &
   ' *** Last line read:', Card
   RETURN
end function RDINP
! ***********************************************************************


!***********************************************************************
!=======================================================================
function no_comma(str)
!=======================================================================
! Remove all commas from string str. This enables input to RDINP of 
! comma separated numbers
!=======================================================================
   character(len = *), intent(in) :: str
   character(255) temp, no_comma
   integer l, k
   
   ! First create a blank string:      
   do l = 1, len(temp)
      temp(l:l) = ' '  
   end do 

   ! Now fill the blank string with the characters
   ! of str, skipping the commas:
   k = 1
   do l = 1, len_trim(str)
      if (str(l:l) /= ',') then
         temp(k:k) = str(l:l)
         k = k + 1
      end if
   end do

   no_comma = trim(temp)
   return
   
end function no_comma
!***********************************************************************



!***********************************************************************
subroutine rdinps2(equal,iUnit,outUnit,str,LENGTH,UCASE)
!=======================================================================
!     2nd version. Returns also length of meaningful string; if UCASE flag
!     is set, the returned string is in upper case to avoid keyboard entry 
!     problems.
!     Added outUnit for output; original version was hard coded to write
!     to unit 6---a bug.
!=======================================================================
      integer i,iUnit,outUnit,ind,first,last,next,length
      character card*(255),chr,cr
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
        if(chr(first) .eq. '*') write(outUnit,'(a)') card
        last = len(card)
        ind = MAX(index(card,'%'), index(card,'!'))
        if(ind .gt. 0) last = ind - 1
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

99    write(outUnit,'(3(1x,a,/))')  &
     ' Terminated. EOF reached by rdinp while looking for input. ', &
     ' last line read:',card
!-----------------------------------------------------------------------
      return
end subroutine rdinps2
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

! ***********************************************************************
SUBROUTINE solve_matrix(model,taumax,nY,nYprev,itereta,nP,nCav,nIns,initial,deviat,iterfbol,fbolOK)
! =======================================================================
! This subroutine solves the continuum radiative transfer problem for a
! spherically symmetric envelope.                      [Z.I., Nov. 1995]
  ! =====================================================================
  use common
  IMPLICIT none
  INTERFACE
     SUBROUTINE RADTRANSF_matrix(pstar,iPstar,nY,nYprev,nP,nCav,nIns,TAUlim,FbolOK,initial,deviat,&
          iterFbol,iterEta,model,us,u_old,fs,T4_ext,emiss)
       integer iPstar, FbolOK,iterFbol,iterEta,model,nY,nYprev,nP,nCav,nIns
       double precision :: pstar,TAUlim,deviat
       double precision,allocatable :: us(:,:),u_old(:,:),fs(:,:),T4_ext(:),emiss(:,:,:)
       logical initial
     END SUBROUTINE RADTRANSF_matrix
     subroutine Analysis(nY,model,us,T4_ext,delta,maxrat)
       integer nY,model
       double precision :: delta,maxrat
       double precision, allocatable ::  us(:,:),T4_ext(:)
     end subroutine Analysis
  END INTERFACE
  !--- parameter 
  integer :: model,iterfbol,fbolOK, itereta, nY, nP, nYprev, nCav, nIns
  double precision :: deviat, taumax
  logical initial
  !--- local variables
  integer iPstar, EtaOK, iY
  double precision pstar,TAUlim, delta
  double precision, allocatable :: fs(:,:),u_old(:,:),us(:,:),T4_ext(:)
  double precision, allocatable :: emiss(:,:,:)
  !----------------------------------------------------------------------
  allocate(fs(nL,npY))
  fs = 0
  allocate(us(nL,npY))
  us = 0
  allocate(T4_ext(npY))
  T4_ext = 0 
  allocate(emiss(nG,nL,npY))
  Emiss = 0
  allocate(u_old(nL,npY))
  u_old = 0
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
             ! ptr(2) is specified in INPUT and controls converg. crit.
             IF (ptr(2).LT.1.0D-6.AND.iterETA.GT.2)THEN
                EtaOK = 1
             ELSE
                CALL WINDS(nY,nYprev,EtaOK)
             END IF
             IF (iterETA.GT.10.AND.EtaOK.EQ.0) THEN
                EtaOK = 2
                IF (iX.NE.0) THEN
                   write(18,*)' *********  WARNING  *********'
                   write(18,*)' Could not converge on ETA in '
                   write(18,*)' 10 iterations.'
                   write(18,*)' *********************************'
                END IF
             END IF
             ! ...or otherwise finish right away
          ELSE
             EtaOK = 1
          END IF
        ! end of loop over ETA
     END DO
  ELSE
     ! solve for slab case
     PRINT*,'Slab case is not implemented for the matrix method!'
     STOP
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
  IMPLICIT none
  INTERFACE
     subroutine matrix(pstar,iPstar,m0,m1,mifront,miback,nP,nY,nPok,nYok,T4_ext)
       integer nP,nY,nPok,nYok,iPstar
       double precision :: pstar 
       double precision,allocatable :: m0(:,:,:), m1(:,:,:), mifront(:,:,:), &
            miback(:,:,:),T4_ext(:)
     end subroutine matrix
     subroutine invert(nY,mat,Us,Em,Uold,omat)
       integer nY
       double precision,allocatable :: Us(:,:), Uold(:,:), Em(:,:,:),&
            mat(:,:,:),omat(:,:)
     end subroutine invert
     subroutine lambda_iter(nY,mat,Us,Em,Uold,omat)
       integer nY
       double precision,allocatable :: Us(:,:), Uold(:,:), Em(:,:,:),&
            mat(:,:,:),omat(:,:)
     end subroutine lambda_iter
     subroutine Emission(nY,T4_ext,emiss,emiss_total)
       integer nY
       double precision, allocatable :: T4_ext(:),emiss(:,:,:),emiss_total(:,:)
     end subroutine Emission
     subroutine find_Text(nY,T4_ext)
       integer nY
       double precision, allocatable :: T4_ext(:)
     end subroutine find_Text
     subroutine find_temp(nY,T4_ext)
       integer :: nY
       double precision, allocatable :: T4_ext(:)
     end subroutine find_temp
     subroutine init_temp(nY,T4_ext,us)
       integer nY
       double precision,allocatable :: us(:,:),T4_ext(:)
     end subroutine init_temp
     SUBROUTINE MULTIPLY(type,np1,nr1,np2,nr2,mat,vec1,omat,flag,q1,q2)
       integer :: type,np1,nr1,np2,nr2,flag
       double precision, allocatable :: mat(:,:,:), vec1(:,:), omat(:,:), &
            q1(:,:), q2(:,:)
     END SUBROUTINE MULTIPLY
     subroutine Bolom(q,qbol,nY)
       integer nY
       double precision, allocatable :: q(:,:), qbol(:)
     end subroutine Bolom
     SUBROUTINE Converg1(nY,Aold,Anew,Aconv,dmax)
       integer nY,Aconv
       double precision :: dmax
       double precision,allocatable ::  Aold(:), Anew(:)
     end SUBROUTINE Converg1
     SUBROUTINE Converg2(nY,Aold,Anew,Aconv,dmax)
       integer nY,Aconv
       double precision :: dmax
       double precision,allocatable :: Aold(:,:), Anew(:,:)
     end SUBROUTINE Converg2
     SUBROUTINE SPH_Int(nY,nP,fs)
       integer nY,nP
       double precision,allocatable :: fs(:,:)
     END SUBROUTINE SPH_Int
     SUBROUTINE FindErr(nY,flux,maxFerr)
       integer nY
       DOUBLE PRECISION maxFerr
       double precision,allocatable :: flux(:)
     END SUBROUTINE FindErr
     subroutine Find_Tran(pstar,nY,nP,T4_ext,us,fs)
       integer nY, nP
       double precision :: pstar
       double precision, allocatable :: T4_ext(:)
       double precision, allocatable :: fs(:,:),us(:,:)
     end subroutine Find_Tran
     subroutine add(np1,nr1,np2,nr2,q1,q2,q3,qout)
       integer np1, nr1, np2, nr2
       double precision, allocatable :: q1(:,:), q2(:,:), q3(:,:),qout(:,:)
     end subroutine add
     subroutine OccltMSG(us)
       double precision,allocatable :: us(:,:)
     end subroutine OccltMSG
     subroutine Flux_Consv(nY,nYprev,Ncav,itereta,iterfbol, fbolom,fbol_em,fbol_sc,fbolOK,maxrat)
       integer nY,nYprev,itereta,fbolOK,iterfbol
       double precision :: maxrat
       double precision, allocatable :: fbolom(:),fbol_em(:),fbol_sc(:)
     end subroutine Flux_Consv
  END INTERFACE
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
 
  if (allocated(omat)) deallocate(omat)
  allocate(omat(nL,npY))
  omat = 0
  if (allocated(emiss_total)) deallocate(emiss_total)
  allocate(emiss_total(nL,npY))
  emiss_total = 0
  if (allocated(UbolChck)) deallocate(UbolChck)
  allocate(UbolChck(npY))
  UbolChck = 0
  if (allocated(Uchck)) deallocate(Uchck)
  allocate(Uchck(nL,npY))
  Uchck = 0
  if (allocated(fbolold)) deallocate(fbolold)
  allocate(fbolold(npY))
  fbolold = 0
  if (allocated(fdebol)) deallocate(fdebol)
  allocate(fdebol(npY))
  fdebol = 0
  if (allocated(fdsbol)) deallocate(fdsbol)
  allocate(fdsbol(npY))
  fdsbol = 0
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
  mat0 = 0
  allocate(mat1(nL,nY,nY))
  mat1 = 0
  allocate(mifront(nL,(20*nY+20),nY))
  mifront = 0
  allocate(miback(nL,(20*nY+20),nY))
  miback = 0
  ! generate spline coefficients for ETA
  CALL setupETA(nY,nYprev,itereta)
  ! evaluate ETAzp
  CALL getETAzp(nY,nP)
  ! generate albedo through the envelope
  CALL getOmega(nY)
  do iY=1,nY
     do iL=1,nL
          omat(iL,iY) = SigmaS(nG+1,iL) / (SigmaA(nG+1,iL) + SigmaS(nG+1,iL))        
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
  ! === Iterations over dust temperature =========
  DO WHILE (Conv.EQ.0.AND.iter.LE.itlim)
     iter = iter + 1
     if (typentry(1).eq.5) call find_Text(nY,T4_ext)
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
     ! CALL Bolom(Utot,Ubol,nY)
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
        IF (abs(maxFerr).LE.accFlux) THEN
           BolConv = 1
        ELSE
           BolConv = 0
        END IF
        ! total criterion for convergence: Utot must converge, and ftot
        ! must either converge or have the required accuracy
        IF (Uconv*(Fconv+BolConv).GT.0) Conv = 1
     END IF
     ! --------------------------------------
  END DO
  ! === The End of Iterations over Td ===
  if (typentry(1).eq.5) call find_Text(nY,T4_ext)
  do iY = 1, nY
     Jext(iY) = sigma/pi * T4_ext(iY)
  end do
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
!!$  DO iY = 1, nY
!!$     write(7777,*) Y(iY),fDebol(iY),fDsbol(iY),fsbol(iY),Td(1,iY)
!!$  END DO
!!$  close(7777)
  ! check whether, and how well, is bolometric flux conserved
  CALL FindErr(nY,fbol,maxFerr)
  ! added in ver.2.06
  IF (maxFerr.LT.accFlux) FbolOK = 1
  IF (iVerb.EQ.2) write(*,*)' Achieved Flux accuracy:',maxFerr
  deviat = maxFerr
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
  IF(iC.NE.0) THEN
     IF (iX.NE.0) write(18,*) 'Calculating intensities'
     call sph_int(nY,nP,fs)
     IF (iVerb.EQ.2) write(*,*) 'Done with SPH_INT(FindInt)'
  END IF
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
  IMPLICIT none
  INTERFACE
     SUBROUTINE LINSYS(Nreal,A,B,X)
       integer Nreal
       DOUBLE PRECISION,allocatable :: A(:,:), B(:), X(:)
     END SUBROUTINE LINSYS
  END INTERFACE
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
  B = 0
  allocate(A(nY,nY))
  A = 0
  allocate(X(nY))
  X = 0
  error = 0
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
  INTERFACE
     SUBROUTINE MYSPLINE(x,N,alpha,beta,gamma2,delta)
       integer :: N
       double precision,allocatable :: x(:),alpha(:,:),beta(:,:),gamma2(:,:),delta(:,:)
     END SUBROUTINE MYSPLINE
     SUBROUTINE WEIGHTS(TAUaux,iP,iL,nZ,alpha,beta,gamma2,delta,wgp,wgm,nY)
       integer iP,iL,nZ,nY
       double precision,allocatable :: TAUaux(:,:,:),alpha(:,:), beta(:,:),&
            gamma2(:,:),delta(:,:),wgp(:,:), wgm(:,:)
     end SUBROUTINE WEIGHTS
     SUBROUTINE NORDLUND(nY,nP,flag,x,f,N1,N2,m,intfdx)
       integer :: nY,nP, flag, N1, N2, m
       double precision :: intfdx
       double precision, allocatable :: x(:),f(:)
     END SUBROUTINE NORDLUND
  END INTERFACE
  !---------------------------------------------------------------------
  allocate(nZ(nP))
  nZ = 0
  allocate(haux(nP))
  haux = 0 
  allocate(TAUaux(nL,nP,nY))
  TAUaux = 0
  allocate(Tplus(nP,nY,nY))
  Tplus = 0
  allocate(Tminus(nP,nY,nY))
  Tminus = 0 
  allocate(xN(nP))
  xN = 0
  allocate(yN(nP))
  yN = 0
  allocate(TAUr(nY))
  TAUr = 0
  allocate(wm(nY))
  wm = 0
  allocate(wmT(nY))
  wmt = 0
  allocate(wp(npY))
  wp = 0
  allocate(alpha(nY,nY))
  alpha = 0 
  allocate(beta(nY,nY))
  beta = 0
  allocate(gamma2(nY,nY))
  gamma2 = 0
  allocate(delta(nY,nY))
  delta = 0
  allocate(wgmatp(nY,nY))
  wgmatp = 0
  allocate(wgmatm(nY,nY))
  wgmatm = 0
  allocate(Yok(nY))
  Yok = 0
  allocate(Pok(nP))
  Pok = 0
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
  secnder = 0
  allocate(yaux(N))
  yaux = 0 
  allocate(deraux(N))
  deraux = 0
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
  INTERFACE
     SUBROUTINE Kint4(TAUaux,iP,iL,nZ,K1p,K2p,K3p,K4p,K1m,K2m,K3m,K4m)
       INTEGER iP, iL, nZ
       DOUBLE PRECISION,allocatable :: TAUaux(:,:,:), K1p(:),K2p(:),K3p(:),&
            K4p(:), K1m(:), K2m(:), K3m(:), K4m(:)
     end SUBROUTINE Kint4
  END INTERFACE
  allocate(K1p(nY))
  K1p = 0
  allocate(K2p(nY))
  K2p = 0
  allocate(K3p(nY))
  K3p = 0
  allocate(K4p(nY))
  K4p = 0
  allocate(K1m(nY))
  K1m = 0
  allocate(K2m(nY))
  K2m = 0
  allocate(K3m(nY))
  K3m = 0
  allocate(K4m(nY))
  K4m = 0
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
  iYins = 0 
  allocate(EtaTemp(nY))
  EtaTemp = 0
  allocate(Yins(nY))
  Yins = 0
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
        ETAdiscr(iY) = ETA(Yloc,nY,nYprev,itereta)
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
SUBROUTINE FindErr_old(nY,flux,maxFerr)
!========================================================================
!This subroutine finds maximum err in flux conservation for both
!spherical and slab case as (fmax-fmin)/(fmax+fmin)   [MN,Aug'99]
!=========================================================================
  use common
  IMPLICIT none
  !--- parameter
  integer nY
  DOUBLE PRECISION maxFerr
  double precision,allocatable :: flux(:)
  !--- local
  INTEGER iY
  DOUBLE PRECISION fmin, fmax, aux, accFbol, tune_acc
  !---------------------------------------------------------------------
  ! tune_acc defines the accuracy in the case of zero flux
  ! tune_acc = 1/2*accFbol*Jext(nY) 
  tune_acc = 6e-1
  if ((slb).and.(left.eq.1).and.(right.eq.1)) then 
     accFbol = min(Ji,Jo)*4*pi*accFlux*tune_acc / Jext(1)
  else if ((sph).and.(right.eq.1)) then
     ! in the spherical case the flux is only zero if there is zero no
     ! central source (Ji=0) therefor accFbol is changed!
     accFbol = max(Ji,Jo)*4*pi*accFlux*tune_acc/ max(Jext(1),Jext(nY))
  else
     accFbol = 0.
  endif
  ! Find the min and max of fbol values
  ! The abs and lower limit on fbol are protection for the case
  ! of completely symmetri! slab illumination. The lower limit
  ! is bound by the numerical accuracy of the flux calculation
  fmin = 1.e30
  fmax = 0.
  DO iY = 1, nY
     aux = flux(iY)
     IF (ksi.eq.1.0) aux = dabs(aux)
     IF (dabs(aux).LE.accFbol) aux = accFbol
     IF(aux.LT.fmin) fmin = aux
     IF(aux.GT.fmax) fmax = aux
     !print*,iY,(fmax - dabs(fmin))/(fmax + dabs(fmin)),flux(iY),fmax,fmin
  END DO
  if (fmax.LE.0.) then
     ! bad solution; overall flux cannot be negative
     maxFerr = 1.
  else if ((fmax.eq.fmin).and.(fmax.eq.accFbol)) then
     maxFerr = accFlux*tune_acc
  else
     maxFerr = dabs(fmax - dabs(fmin))/(fmax + dabs(fmin))
  end if
  ! ----------------------------------------------------------------------
  RETURN
END SUBROUTINE FindErr_old
! *************************************************************************

SUBROUTINE FindErr(nY,flux,maxFerr)
!========================================================================
!This subroutine finds maximum err in flux conservation for both
!spherical and slab case as (fmax-fmin)   [MN,Aug'99]
!=========================================================================
  use common
  IMPLICIT none
  !--- parameter
  integer nY
  DOUBLE PRECISION maxFerr
  double precision,allocatable :: flux(:)
  !--- local
  INTEGER iY
  DOUBLE PRECISION fmin, fmax, aux, accFbol, tune_acc

  !---------------------------------------------------------------------
  tune_acc = 1.
  ! Find the min and max of fbol values
  fmin = minval(flux(:nY))
  fmax = maxval(flux(:nY))
  maxFerr = (fmax-fmin)/(4*pi) !/(4*pi*maxval(Jext(:nY)))
  RETURN
END SUBROUTINE FindErr
! *************************************************************************
! ***********************************************************************
SUBROUTINE WINDS(nY,nYprev,EtaOK)
! =======================================================================
! This subroutine takes care of the interface between radiatively driven
! winds and radiative transfer.  It is entered
! after a radiative transfer calculation with given eta.
! This sub caclulates the reddening profile phi and passes it to the dynamics
! module, which returns the velocity and density profiles corresponding
! to the given phi. Convergence is achived when the eta returned from
! the dynamics calculation is the same as that used to produce phi.
!
! Notations follow EI01 (MNRAS 327, 403)
! =======================================================================
  use common
  implicit none
  INTERFACE
     subroutine Simpson(n,n1,n2,x,y,integral)
       integer n, n1, n2
       double precision integral
       double precision,allocatable ::  x(:), y(:)
     end subroutine Simpson
     SUBROUTINE DYNAMICS(eps_loc, f, uScale, phi_loc, u, zeta,  &
          P_dyn, eta, nY, acc, err, ver)
       integer nY, ver, err
       double precision eps_loc, f, uScale, acc, P_dyn
       double precision, allocatable :: eta(:), zeta(:,:), u(:), phi_loc(:)
     END SUBROUTINE DYNAMICS
     SUBROUTINE ChkConv(nY,accuracy_loc,Aold,Anew,Aconv_loc)
       INTEGER Aconv_loc,nY
       DOUBLE PRECISION accuracy_loc 
       DOUBLE PRECISION,allocatable:: Aold(:), Anew(:) 
     END SUBROUTINE CHKCONV
  END INTERFACE
  ! --- parameter
  integer EtaOK, nYprev, nY
  ! -- local
  integer iY, iL, err, ver
  double precision, allocatable :: ETAold(:),faux(:), reddn(:), w(:)
  double precision acceta, resaux,  Qfid, phi1, localP, eps_loc, GammaMax, wScale, uacc

  ! the parameter ver determines the version of the velocity formal solution
  ! 1 for linear, 2 quadratic.  it is specified in input and carried in /dyn/
  ! the quadratic solution is from equation d1 in ei01. the linear is
  ! obtained similarly from the differential equation 24 by using dw^2 = 2wdw
  ! and dividing through by 2w
  ! -----------------------------------------------------------------------

  allocate(ETAold(npY))
  ETAold = 0
  allocate(faux(nL))
  faux = 0
  allocate(reddn(npY))
  reddn = 0
  allocate(w(npY))
  w = 0
  ver = 1
  IF (iX.GE.1) THEN
     write(18,*)' Doing Dynamics'
     IF (ver.EQ.1) THEN
        write(18,*)' Linear version of velocity formal solution'
     ELSEIF (ver.EQ.2) THEN
        write(18,*)' Quadratic version of velocity formal solution'
     ELSE
        write(12,*)' **************************** '
        write(12,'(a,i3)')'  Illegal Input ver = ', ver
        write(12,*)'     ver must be 1 or 2!      '
        write(12,*)'       PROGRAM STOPPED        '
        write(12,*)' **************************** '
        stop
     END IF
  END IF
  IF(iVerb.EQ.2) write(*,*)' Doing Dynamics'
  ! so far it works for nG=1 only:
  IF (nG.GT.1) THEN
     write(12,*)' **************************** '
     write(12,*)' Change dynamics sub to nG>1! '
     write(12,*)'       PROGRAM STOPPED        '
     write(12,*)' **************************** '
     stop
  END IF
  ! -----------------------------------------------------------------------
  ! assign input parameters to local variables
  GammaMax = ptr(1)
  eps_loc = pow
  ! accuracy for velocity convergence same as for Utot:
  uacc = accConv
  ! extinction efficiency at the fiducial wavelength
  Qfid = SigExfid/aveA
  ! calculate Qstar and the scale factor of w:
  DO iL = 1, nL
     Faux(iL) = (SigmaA(1,iL)+SigmaS(1,iL))*ftot(iL,1)/lambda(iL)
  END DO
  CALL Simpson(nL,1,nL,lambda,Faux,resaux)
  !
  ! Qstar is from EI01 equation 4, wScale is from equation 29
  !
  Qstar = resaux / aveA
  wScale = TAUfid/Qfid
  !
  ! -----------------------------------------------------------------------
  ! Here's the eta that was used in the radiative transfer
  DO iY = 1, nY
     ETAold(iY) = ETAdiscr(iY)
  END DO
  ! and here's the resulting reddening profile
  DO iY = 1, nY
     DO iL = 1, nL
        Faux(iL) = (SigmaA(1,iL)+SigmaS(1,iL))*ftot(iL,iY)/lambda(iL)
        print*,'still single grain L7370'
        stop
     END DO
     CALL Simpson(nL,1,nL,lambda,Faux,resaux)
     if (iY.eq.1)    phi1 = resaux
     reddn(iY) = resaux / phi1  ! eq.(3);reddn(iY)=phi_loc in Dynamics
  END DO
  ! Now Find new ETA
  err = 0
  CALL DYNAMICS(eps_loc, GammaMax, wScale, reddn, w, vrat,  &
       localP, ETAdiscr, nY, uacc, err, ver)
  !!** Prdw=P is stored in common /dyn/, needed in Analysis [MN]
  Prdw = localP
  ! and check convergence (ptr(2) is specified in INPUT)
  accETA = ptr(2) * accuracy
  CALL ChkConv(nY,accETA,ETAold,ETAdiscr,EtaOK)
  IF (iX.GE.1) THEN
     write(18,'(2(a,1pe10.3))') '  P = ', localP, '  gmax = ', gmax
     write(18,*) '     Y    ugas(new)   tauF      ETAold    ETAnew    ratio'
     DO iY = 1, nY
        !**************************************
        ! for output compatibility; we can do away with qF and tauF, which
        ! have only nostalgic reasons.  EI01 never uses them
        ! ugas is in common /dyn/ and used in Analysis [MN]
        ugas(iY) = Qstar*w(iY)
        qF(iY)   = (Qstar/Qfid)*reddn(iY)
        Faux(iY) = qF(iY)*ETAdiscr(iY)
        CALL SIMPSON(nY,1,iY,Y,Faux,resaux)
        tauFdyn(iY) = TAUfid*resaux
        !**************************************
        accETA = ETAold(iY) / ETAdiscr(iY)
        write(18,'(1p,6e10.3)')Y(iY), ugas(iY), tauFdyn(iY),   &
             ETAold(iY), ETAdiscr(iY),accETA
     END DO
     IF (EtaOK.EQ.1) THEN
        write(18,*)' Convergence on Eta achieved'
     ELSE
        write(18,*)' Convergence on Eta not achieved.'
        write(18,*)' Going to the next iteration.'
     END IF
  END IF
  ! save Y to Yprev and nY to nYprev
  DO iY = 1, nY
     Yprev(iY) = Y(iY)
  END DO
  nYprev = nY
  ! ----------------------------------------------------------------------
  deallocate(ETAold)
  deallocate(faux)
  deallocate(reddn)
  deallocate(w)
  RETURN
END subroutine Winds
! ***********************************************************************

! ***********************************************************************
SUBROUTINE DYNAMICS(eps_loc, f, uScale, phi_loc, u, zeta,  &
     P_dyn, eta, nY, acc, err, ver)
! =======================================================================
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
  use common
  implicit none
  INTERFACE
     SUBROUTINE CalcVel(eps_loc, f, ws, phi_loc, wold, w, zeta, nY, ver)
       INTEGER ver, nY
       DOUBLE PRECISION eps_loc,f,ws
       double precision, allocatable :: wold(:), w(:), phi_loc(:),zeta(:,:)
     END SUBROUTINE CalcVel
     SUBROUTINE CalcETA(w, zeta, Eta, EtaINT, nY)
       INTEGER  nY
       DOUBLE PRECISION EtaINT
       DOUBLE PRECISION,allocatable :: w(:), Eta(:), zeta(:,:)
     END SUBROUTINE CalcETA
     SUBROUTINE ChkConv(nY,accuracy_loc,Aold,Anew,Aconv_loc)
       INTEGER Aconv_loc,nY
       DOUBLE PRECISION accuracy_loc 
       DOUBLE PRECISION,allocatable:: Aold(:), Anew(:) 
     END SUBROUTINE CHKCONV
  END INTERFACE
  ! --- parameter
  integer nY, ver, err
  double precision eps_loc, f, uScale, acc, P_dyn
  double precision, allocatable :: eta(:), zeta(:,:), u(:), phi_loc(:)
  ! --- local
  integer iY, itr, ETAconv, uconv, itmax
  double precision  N, wf, k, e1
  double precision, allocatable :: etaold(:), uold(:)
  data   itMax/100/, k/0.4/

  allocate(etaold(npY))
  etaold = 0
  allocate(uold(npY))
  uold = 0
  allocate(zeta(npG,npY))
  zeta = 0

  !  we may wish to control itMax and k as input parameters
  ! -----------------------------------------------------------------------
  !  for information: phi_loc(nY) = reddn(nY)     [MN]
  !                 zeta(npG,npY) = vrat(npG,npY)
  ! Initial approximation for u(y) from EI01, eq C6
  ! wf is from eq. 29 with epsilon correction (eq. C8)
  wf = (1.0d+00/(1.0d+00 - eps_loc))*phi_loc(nY)*uScale
  ! add a correction for the finite outer radius so that wf = u(nY):
  e1 = 1.0d+00 - eps_loc**(1.0d+00/k)
  wf = wf/(1.0d+00 - e1/Y(nY))**k
  ! and now calculate all u from eq. C6
  DO iY = 1, nY
     uold(iY) = wf*(1.0d+00 - e1/Y(iY))**k
     ! initial eta is irrelevant; might as well use
     ! the one passed from radiative transfer:
     ETAold(iY) = eta(iY)
  END DO
  ! ITERATIONS until u and eta converge within acc
  DO itr = 1, itMax
     Call CalcVel(eps_loc,f,uScale,phi_loc,uold,u,zeta,nY,ver)
     CALL CalcETA(u, zeta, eta, N, nY)
     ! here N=EtaINT found in CalcEta
     P_dyn = dsqrt(uScale/N)       !eq.(46) in IE'01
     ! check convergence of u and Eta
     CALL ChkConv(nY,acc,uold,u,uconv)
     CALL ChkConv(nY,acc,ETAold,eta,ETAconv)
     ! convergence required for both u(y) and ETA(y)
     err = 1 - ETAconv * uconv
     IF (err.NE.0) THEN
        ! did not converge, repeat the exercise...
        DO iY =1, nY
           uold(iY) = u(iY)
           ETAold(iY) = eta(iY)
        END DO
     ELSE
        ! we're done:
        IF (iX.GE.1) write(18,'(a35,i3)')' Number of iterations to converge:',itr
        RETURN
     END IF
  END DO
  ! -----------------------------------------------------------------------
  return
  deallocate(eta)
  deallocate(etaold)
  deallocate(u)
  deallocate(uold)
  deallocate(phi_loc)
  deallocate(zeta)
end subroutine Dynamics
! ***********************************************************************

! ***********************************************************************
SUBROUTINE CalcVel(eps_loc, f, ws, phi_loc, wold, w, zeta, nY, ver)
! =======================================================================
! Calculates the scaled gas velocity w(y) from wold(y), the previous
! velocity profile, and the given reddening profile phi
! The calculation follows the formal solution in Appendix D
! ver = 1 triggers the linear version, 2 the quadratic (eq. D1)
! All symbols are as defined there, ws = tauV/QV and P2 = P^2
! =======================================================================
  use common
  IMPLICIT none
  INTERFACE
     subroutine Simpson(n,n1,n2,x,y,integral)
       integer n, n1, n2
       double precision integral
       double precision,allocatable ::  x(:), y(:)
     end subroutine Simpson
  END INTERFACE
  ! --- parameter
  INTEGER ver, nY
  DOUBLE PRECISION eps_loc,f,ws
  double precision, allocatable :: wold(:), w(:), phi_loc(:),zeta(:,:)
  ! --- local
  INTEGER iY
  DOUBLE PRECISION  g, N, P2, ww1, aux
  double precision, allocatable :: z(:), zz(:), F1(:), F2(:) 
  ! -----------------------------------------------------------------------
  
  allocate(z(npY))
  z = 0
  allocate(zz(npY))
  zz = 0
  allocate(F1(npY))
  F1 = 0
  allocate(F2(npY))
  F2 = 0
  ! first get the drift profile zeta
  CALL CalcDrift(phi_loc, wold, zeta, nY)
  ! then the normalization N (= EtaINT; eq D3)
  DO iY = 1, nY
     F1(iY) = zeta(1,iY)/(wold(iY)*Y(iY)*Y(iY))
  END DO
  CALL SIMPSON(nY,1,nY,Y, F1, N)
  ! and finally P (eq. D3):
  P2 = 2.0D+00*ws/N
  ! Now the two versions diverge:
  IF (ver.eq.2) THEN
     ! Quadratic version (eq. D1):
     ! Get the profile z = integral(zeta*phi/y^2) and gmax
     ! (= Gamma_min), the maximum gravitational correction
     ! as defined in eqs D4 and D5:
     F1(1) = 0.0D+00
     z(1)  = 0.0D+00
     gmax  = 1.0D+00/zeta(1,1)
     DO iY = 2, nY
        F1(iY) = zeta(1,iY)*phi_loc(iY)/(Y(iY)*Y(iY))
        CALL SIMPSON(nY,1,iY,Y, F1, z(iY))
        aux = (1.0D+00 - 1.0D+00/Y(iY))/z(iY)
        if (aux.gt.gmax) gmax = aux
     END DO
     g = f/gmax
     ! All ready.  Calculate the new w1 (ww1 = w1^2):
     aux  = eps_loc*eps_loc
     ww1  =(aux/(1.0D+00-aux))*P2*(z(nY)-g*(1.0D+00-1.0D+00/Y(nY)))
     w(1) = dsqrt(ww1)
     ! and now the rest of the profile:
     DO iY = 2, nY
        aux   = ww1 + P2*(z(iY)-g*(1.0D+00-1.0D+00/Y(iY)))
        w(iY) = dsqrt(aux)
     END DO
  ELSE
     ! Linear version.
     ! obtained from the differential equation 24
     ! by using dw^2 = 2wdw and dividing through by 2w
     ! Now we need the profiles z = integral(zeta*phi/w*y^2),
     ! zz = integral(1/w*y^2) and then gmax = max(zz/z)
     F1(1) = 0.0D+00
     F2(1) = 0.0D+00
     z(1)  = 0.0D+00
     zz(1) = 0.0D+00
     gmax  = 1./zeta(1,1)
     DO iY = 2, nY
        F2(iY) = 1.0D+00/(wold(iY)*Y(iY)*Y(iY))
        F1(iY) = F2(iY)*zeta(1,iY)*phi_loc(iY)
        CALL SIMPSON(nY,1,iY,Y, F1, z(iY))
        CALL SIMPSON(nY,1,iY,Y, F2, zz(iY))
        aux = zz(iY)/z(iY)
        if (aux.gt.gmax) gmax = aux
     END DO
     g = f/gmax
     ! All ready.  Calculate the new w1:
     w(1)  = (eps_loc/(1.0D+00-eps_loc))*0.5D+00*P2*(z(nY) - g*zz(nY))
     ! and now the rest of the profile:
     DO iY = 2, nY
        w(iY) = w(1) + 0.5D+00*P2*(z(iY) - g*zz(iY))
     END DO
  END IF
  ! -----------------------------------------------------------------------
  deallocate(z)
  deallocate(zz)
  deallocate(F1)
  deallocate(F2)
  RETURN
END SUBROUTINE CalcVel
! ***********************************************************************

! ***********************************************************************
SUBROUTINE CalcETA(w, zeta, Eta, EtaINT, nY)
! =======================================================================
! Calculates the dimensionless density profile ETA(y) from EI eq. 25 given
! the velocity profile w(y) and its corresponding drift zeta(y)
! =======================================================================
  use common
  IMPLICIT none
  INTERFACE
     subroutine Simpson(n,n1,n2,x,y,integral)
       integer n, n1, n2
       double precision integral
       double precision,allocatable ::  x(:), y(:)
     end subroutine Simpson
  END INTERFACE
  ! --- parameter
  INTEGER  nY
  DOUBLE PRECISION EtaINT
  DOUBLE PRECISION,allocatable :: w(:), Eta(:), zeta(:,:)
  ! --- local
  INTEGER iY
  ! =====================================================================
  DO iY = 1, nY
     Eta(iY) = zeta(1,iY)/(w(iY)*Y(iY)*Y(iY))
  END DO
  ! now normalize eta:
  CALL SIMPSON(nY,1,nY,Y,Eta,EtaINT)
  DO iY = 1, nY
     Eta(iY) = Eta(iY)/EtaINT
  END DO
! -----------------------------------------------------------------------
  RETURN
END SUBROUTINE CalcETA
! ***********************************************************************

! ***********************************************************************
SUBROUTINE CalcDrift(phi_loc,w,zeta,nY)
! =======================================================================
! Calculates the drift profile zeta from EI01 eq. 24
! without the correction for sub-sonic drift (theta = 0).
! =======================================================================
  use common
  IMPLICIT none
  !INTEGER npY, npP, npX, npL, npG, npR
  INTEGER nY, iY
  DOUBLE PRECISION phi_loc(npY), w(npY), zeta(npG,npY)
! -----------------------------------------------------------------------
  DO iY = 1, nY
     zeta(1,iY) = 1.0D+00 / (1.0D+00 + dsqrt(phi_loc(iY)/w(iY)))
  END DO
! -----------------------------------------------------------------------
  RETURN
END SUBROUTINE CalcDrift
! ***********************************************************************
