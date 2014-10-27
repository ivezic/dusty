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


