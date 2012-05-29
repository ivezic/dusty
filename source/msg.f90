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
