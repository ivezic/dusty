!! this is common.f90 for version Dusty3.10 from August 2010 [MN]
!! Make sure you compile with the correct common file.

module common
  implicit none
  
  double precision pi, sigma, Gconst, r_gd, clight, mprot
  parameter (pi=3.141592653589793116)
  parameter (sigma=5.67D-08)
  parameter (Gconst= 6.67D-11)
  parameter (r_gd = 2.0D+02)
  parameter (clight = 3.0D+08)
  parameter (mprot = 1.67D-27)
! =============================================================================
  

! =============================================================================
! Commons with various parameters related to numerical accuracy.
! All are initialized in subroutine INPUT.                [ZI'95;MN'99; MD'07]
! Former COMMON /mumerics/
! =============================================================================
! accConv - desired accuracy for energy.density convergence, by default it is
!           set to accConv = accuracy/500.
! accFbol - desired accuracy for bol.flux convergence, by default
!           accFbol = 10.*accConv
! accRomb - desired accuracy for numerical integration in ROMBERG
! EtaRat  - limit on the ratio of density profile for any two consecutive
!           radial grid points, used in Ygrid.
! delTAUs! - max. increase of scaled TAU/TAUtot, used in Ygrid
! fac! - max. increase in the ratio of two y, used in Ygrid
! Ncav - number of p-rays inside the cavity
! Nins - number of p-rays inserted per radial step
! -----------------------------------------------------------------------------

  integer Ncav, Nins
  double precision accuracy, accConv, dtau, init_tau, dynrange
! =============================================================================
  
! =============================================================================
! Common statements for the grids. /grids1/ is for integers, /grids2/ for reals.
!                                                      [ZI'95-97; MN'99; MD'07]
! =============================================================================
! npY - max size for the radial (Y) grid
! npP - max size for the impact parameter (P) and angular (mu) grids
! npX - max size for the x-grid x(npX) in the optional disk calculation
! npL - max size for the wavelength grid (as defined in 'lambda_grid.dat')
! npG - max size for the number of grains in future MG version
! npR - max size for the output inclination angle for slab (in *.i### files)
! nY, nP, nL - the actual sizes of the grids: nY, nP are determined in
!         subroutines Ygrid and Pgrid, respectively, while nL = npL.
! nPcav - the number of rays in the impact parameter grid, passing through
!            the inner cavity. They are calculated in subroutine Pgrid.
! Yprev,nYprev - used in case of iterations over ETA (RDW) to keep Y and nY
!            from previous ETA iteration. Used in ChkFlux, Ygrid and Winds.
! bOut(npP+2) - the impact parameter grid with 2 additional rays to take care
!            of pstar. Found in GetbOut, used in FindInt, Convolve and Visibili.
! -----------------------------------------------------------------------------

  integer  npY,npP,npX, npL, npG, npR
  include 'userpar.inc'
  parameter (npG=1)
  integer nY, nYprev, nP, nPcav, nL
  double precision Y(npY), Yprev(npY),P(npP),bOut(npP+2),lambda(npL)
! =============================================================================
  
! =============================================================================
! Common statement and variable definitions for array iYfirst(iP) which shows
! the starting radial position for the line of sight with impact parameter P(iP),
! and array YPequal(iP) which is 1 for P(iP).EQ.Y(iYfirst(iP)).
!                                          .                            [ZI,'95]
! -----------------------------------------------------------------------------

  integer iYfirst(npP), YPequal(npP), Plast(npY)
! =============================================================================
  
! =============================================================================
! Common statements for variables describing the dust density law ETA.
! /dens1/ is for integers, /dens2/ - for character,  /dens3/ - for reals,
! /dens4/ - for logical variables.
!                                                       [ZI'95-97; MN'99; MD'07]
! =============================================================================
! denstyp - flag for density type, entered as integer in 'fname.inp':
!           Flag values: 1-power law, 2-exponential density law,
!           3-RDW (radiatively driven winds), 4- analytical (gray-body) RDW
!           approximation, 5-ETA from a file. Internally, denstyp.eq.0 is for slab
!           (Private option for full RDW is denstyp.eq.6)
!!
!! FOR IMPROVEMENT: no need to carry both logical and denstyp in common.
!! denstyp will be local for Sub INPUT then only Logical flags will be carried.
!!
! Ntr     - number transition points for broken power low
! iterETA - counter over ETA iterations (in case of RDW, denstyp 5 or 6)
! nYEtaf  - # pts. for ETA from file nameETA
! pow     - real parameter for the chosen density law: power for 1, sigma for 2
!           [rho = exp(-(y/sigma)**2)], v1/ve for 3 (the ratio of expansion
!           velocities at the inner and outer radii).
! ptr(10), Ytr(10) - powers and scaled transitional radii for denstyp=1
! ETAcoef(npY,4) - the four coefficients for spline approximation of ETA
! ETAdiscr(npY) - needed in ETA interpolation
! yEtaf(1000),Etaf(1000) - the user supplied table for Eta(y) in file nameETA
! Yout - the relative thickness, Yout=rout/r1.
! -----------------------------------------------------------------------------
  
  integer Ntr, iterETA, nYEtaf
  character*235 nameETA
  logical POWD, EXPD, RDW, RDWA, FILD, RDWPR, SLB, SPH
  double precision pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),   &
       ETAdiscr(npY), yEtaf(1000), Etaf(1000)
! =============================================================================

! =============================================================================
! This file contains the common statement for variables describing the
! spectrum of the illuminating source(s).            [ZI,'96, MN,'99,'00]
! Former COMMON /source/
! =============================================================================
!  Source index: 1- central source for sphere or left source for slab;
!                2- external source for sph or right source for slab;
!  Left, Right - flags for presence of sources
!
!  These are initialized in Input or InpStar:
!  typEntry(2) - type of boundary condition - Tsub(1) or flux in some form:
!                flux Fi, luminosity and distance, or Teff from the flux.
!  startyp(2) - type of spectral shape ,
!      nBB(2) - # of black bodies with temperatures Tbb(2,10) and relative
!                luminosities rellum(2,10);
!  Nlamtr(2)   - # transitional pts. for broken power law for input spectrum
!  lamtr(2,101)-transition wavelength; klam(3,100)- powers of the power law..
!
!  nameStar(2) - file with stellar spectrum
!  xSiO - depth in % of the SiO feature for Engelke-Marengo function
!  r1rs - the ratio of inner shell radius/stellar radius, found in Analysis
!
!  mu1, mu2 - the cosines of the left and right illumination angles in slab
!  ksi - the ratio of the right/left bol. fluxes (<1) for slab
!
!  chi - relative contribution of the external illumination of spherical shell
! -----------------------------------------------------------------------------

  integer startyp(2), Nlamtr(2), nBB(2), typEntry(2), Left, Right
  character nameStar(2)*235
  double precision  lamtr(2,101), klam(2,100), &
       Tbb(2,10), rellum(2,10), Tstar(2), mu1, ksi, mu2, xSiO, &
       r1rs, chi, dilutn
! =============================================================================
  
! =============================================================================
! Commons for variables related to dust optical properties.
!                                                      [ZI'95-97; MN'99; MD'07]
! Former COMMON /optprop/
! =============================================================================
!  The following variables are initialized in Subroutine Input:
!  Nfiles -    number of user supplied .nk files
!  qsd,a1,a2 - parameters of the built-in size distribution functions
!  szds -      type of size distribution function
!  top -       opt.properties index ("top" means "type of properties")
!  xC(10), xCuser(10) - fractional number abundances for built-in and
!              user supplied components, respectively
!  Tsub(npG) - sublimation temperature for a given grain type
!
!  The following variables are determined in Subroutine GetTau:
!  iLfid -     the index of the fiducial wavelength lamfid, initilized in
!              subroutine GetTau
!  TAUfid -    opt. depth at lamfid
!  TAUmax -    max value of TAUtot(npL), used in Ygrid and ChkFlux
!  TAUtot(npL)- total opt.depth
!
!  These are determined in Subroutine GetOptPr
!  SigmaA(npG,npL), SigmaS(npG,npL) - absorp. and scatt. cross sections
!  SigExfid - the extinction cross section at lamfid
!  aveV -     for a single grain size it is the single grain volume,
!             otherwise it is the volume averaged over size distribution.
!  aveA -     for a single grain size it is the eff. grain area (Pi*a^2),
!             otherwise it is the effective area averaged over the size
!             distribution function
!
! The following arrays carry an extra y-dependence as a preparation for the
! mutigrain code, and are determined in subroutine GetOmega:
!  abund(npG,npY) - relative abundances of grain types, currently 1
!  omega(npL,npY) - scattering albedo
! -----------------------------------------------------------------------------
  
  integer iLfid, szds, top, Nfiles,noprint
  double precision TAUtot(npL),SigmaA(npG,npL), SigmaS(npG,npL), &
       Tsub(npG), abund(npG,npY), TAUmax, xC(10), xCuser(10), &
       SigExfid, TAUfid, taufid0,lamfid, qsd, a1, a2, aveV, aveA
! =============================================================================


! =============================================================================
! Common statement for dynamical quantites (radiatively driven winds).
!                                                                [ZI'96; MN'99]
! Former COMMON /dyn/
! =============================================================================
! Gamma(npY) - the ratio of gravitational to the radiative force on gas.
! qF(npY) - effective flux averaged extinction efficiency
! ugas(npY) - the gas velocity profile scaled by its terminal velocity
!             (at Yout)
! vrat(npG,npY) - the ratio v(y)/vd(y), v is the gas velocity, vd is the
!                 dust velocity (different for each dust component)
! I1,I2,I3,CMdot,Cve,CM,Cr1 - conversion constants for dynamics (found
!                in Analysis)
! G1, Ginf, Prdw, delta, w1, Phi, PI - private quantities
! -----------------------------------------------------------------------------
  
  integer ver
  double precision ugas(npY), qF(npY), vrat(npG,npY), Gamma(npY),  &
       I1, I2, I3, CMdot, Cve, CM, Cr1, G1, Ginf, Prdw, gmax,      &
       winf, Phi, PIrdw, QV, Qstar, zeta1, tauFdyn(npY)
! =============================================================================

! =============================================================================
! Related to dynamics calculation (RDW) [ZI, ME,2002]
! -----------------------------------------------------------------------------
  
  integer nW
  double precision wav(npL), Eff(npL), f1, f2
! =============================================================================
  

! =============================================================================
! This file contains the common statement for the various solution arrays.
! =============================================================================
!  Genrally, U is energy density, f is flux profile lambda*f_lambda, scaled by
!  the bolometri! flux. In array names 'de' stands for dust emission component,
!  'ds' for dust scattering, 's' is for star (or source).
!                                                     [ZI'95; MN'97,'00; MD'07]
!  Former COMMON /solution/ 
! -----------------------------------------------------------------------------  
! tauOut(npL) - total optical depth along a line of sight, used in sph_int
  integer nYok, nPok, moment_loc,moment
  double precision fde(npL,npY), fds(npL,npY), Utot(npL,npY), ftot(npL,npY), Td(npG,npY), &
       Ubol(npY), fbol(npY), Spectrum(npL), SmC(30,99), tauF(npY), tr(npY), rg(npG,npY), &
       Intens(npL,npP+2), IntOut(20,npP+2), tauZout(npP+2), tauOut(npL), Eps(npY), &
       fsL(npL,npY), fsR(npL,npY), fsLbol(npY), fsRbol(npY),fsbol(npY), RPr(npY), Jext(npY), &
       Ude(npL,npY), Uds(npL,npY), ETAzp(npP,npY), Ji, Jo, Psi, Psi0, RPr1
  
! =============================================================================


! =============================================================================
! Common statement for flags controling production of output files.
! Here are also arrays related to intensity output from spherical shell.
!          .                                         [ZI'95; MN'97,'00; MD'07]
! Former COMMON /output/
! =============================================================================
!  Flags set in the master 'dusty.inp':
!   iVerb - for additional screen printout
!  Flags set by the user in each 'fname.inp':
!   iOUT  - used in Cllose, now is obsolete ;
!    see if iOut can be removed from the code.
!  iSPP  - for spectral properites file production
!  iA - for spectra; iB - for radial profiles; iC - for images,
!  iV - for visibility files; iPSF - for convolved images (private option);
!  iJ - for energy density profiles at user-selected radii
!  iX - for message files;
!  iInn - private option, set from inside Sub Input, it is for additional
!         printout in file 'fname.err'. Good for checking convergence.
!
!  Te_min - calculated in OccultMSG, to warn user for min required source Teff
!         (this is in case of BB type spectrum)
!!
!! see if Te_min can be removed.from common
!!
!  The following variables are used in imaging subroutines:
!  NlambdaOut - number of user required wavelengths (up to 20)
!  LambdaOut(20) - their wavelengths;  Visib(20,1000) - visibility results found
!  in subroutine Visibili; Offset(Nconv) - used for convolved images;
!  ConvInt(20,1000) - convolved intensity
!  nJOut - user-selected number of radii for J-output
!  YJOut(10) - the Y(iY) radii for J-output
! -----------------------------------------------------------------------------
 
  integer iVerb,iPhys,iA,iB,iC,iX,iInn,iPSF,iV,NlambdaOut,Nconv,Nvisi,iD,iJ,nJOut
  character*100 zline(999)
  double precision LambdaOut(20), ConvInt(20,1000), Visib(20,1000),     &
       Offset(1000), qtheta1(1000), Te_min, YJOut(10), JOut(npL,10)
! =============================================================================
  
  
! =============================================================================
!  Common flags for status.
! =============================================================================
!  iWARNING - counter for number of warnings issued by the code
!  iERROR - counter for number of errors encountered by the code
!  iCUMM - counter for any warnings or errors. If 0, message issued that all OK.
!  Former COMMON /status/
! -----------------------------------------------------------------------------

  integer iWARNING, iERROR, iCUMM
! =============================================================================
  
  
! =============================================================================
! Common statement for variables used in imaging subroutines.       [ZI, '97]
! Quantites describing the user supplied point spread function.
! Added psfArea, used in normalization of convolved images [MN, Sep.'04]
!  Former COMMON /psf1/ and /psf2/
! -----------------------------------------------------------------------------
  
  integer psftype, Npsf, iLambda
  double precision kPSF(20), FWHM1(20), FWHM2(20), Theta1, &
       xpsf(1000), ypsf(1000), psfArea(20)
! =============================================================================
  
  
  
! =============================================================================
! Related to imaging                                                   [ZI,'97]
! -----------------------------------------------------------------------------

  integer ftype
  double precision Ckn, Cxout, Cphi, Cqtheta
! =============================================================================
  
! =============================================================================
! Common statement for slab intensity.                            [MN, '99-00]
! Used in Subroutines SLBRadT and SLBIntens and Real function Sexp(t).
! =============================================================================
!  muobs - element of mu-array, needed in SLBIntens()
!  tauT = TAUslb(iL,nY) is the total opt.depth (lambda dependent)
!  Sfn - the local value of the source function               [MN, Jan'00]
! -----------------------------------------------------------------------------

  integer Nmu, transmit
  double precision theta(npR), muobs, tauT, Sfn
! =============================================================================
  
! =============================================================================
! Common statement for slab calculation                           [MN, '99-00]
! -----------------------------------------------------------------------------

  double precision TAUslb(npL,npY), fpbol(npY), fmbol(npY), fmed, &
       SLBIntm(npR,npL), SLBIntp(npR,npL), IstR(npL),maxrat
!=============================================================================

end module common
