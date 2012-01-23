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
  logical ::  Equal, noEqual, UCASE
  character(len=235) :: stdf(7),str,nameIn,nameOut
  character(len=235),allocatable :: nameNK(:),nameQ(:)
  double precision :: tau1,tau2,Lum,dist,RDINP,spec_scale, &
       dilutn,Tstar(2),th1,th2,xC(10),xCuser(10),sum,qsd,a1,a2,&
       var1,var2,var3
  integer :: GridType,Nmodel,L,right,left,top,iG,iFiles,nFiles,szds

  interface
     subroutine inp_rad(shp,spec_scale)
       double precision :: shp(:),spec_scale
     end subroutine inp_rad
  end interface
  
!!$  integer i, iG, nG, Nmodel, EtaOK,GridType, istop, nLs, Nrec, &
!!$       ioverflw, Nmax, nLambdam, Nis, imu, denstyp, ang_type, L
!!$  ! Nmax is the size of user supplied eta file
!!$  ! nLambdam is the max number entries for a user supplied stellar spectrum
!!$  parameter (Nmax = 1000, nLambdam = 10000, nis = 2)
!!$  double precision tau1, tau2, sum, a, b, xx(Nmax), e(Nmax), &
!!$       aa(Nmax), bb(Nmax), tauIn(Nrec), ceta, x1, psf1, &
!!$       lum, dist, var1, var2, var3, spec_scale, th1, th2,      &
!!$       lambdas(nLambdam), Llamstar(nLambdam), res, value, RDINP
!!$  character lamstr(20)*72, strpow*72, strg*40, version*(*),stdf(7)*235
!!$  character*(*) nameIn,nameOut, nameQ(npG),nameNK(10), namepsf*100, &
!!$       nametau*100, anggrid*100,str*235
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
  write(12,'(2x,a140)')nameIn
  write(12,*)' '
  ! Open input file
  open(1,file=nameIn,status='old')
  rewind(1)
  !********************************************
  !** I. Geometry **
  !********************************************
  call rdinps2(Equal,1,str,L,UCASE)
  if(str(1:L).eq.'SPHERE') then
   slb = .false.
   sph = .true.
  elseif(str(1:L).eq.'SLAB') then
   slb = .true.
   sph = .false.
  end if
  !********************************************
  !** II. Physical parameters **
  !********************************************
  ! (1) Flags for presence of sources for slab these have the meaning
  ! of "left" and "right" source,
  ! illumination
  call rdinps2(Equal,1,str,L,UCASE)
  if(str(1:L).eq.'ON') left = 1
  if(str(1:L).eq.'OFF') left = 0
  if(left.eq.0) then
     if(slb) then
        call msg(23)
        left = 1
     end if
  end if
  call rdinps2(Equal,1,str,L,UCASE)
  if(str(1:L).eq.'ON') right = 1
  if(str(1:L).eq.'OFF') right = 0
  if (slb.and.(left.eq.0.and.right.eq.1)) then
     left=1
     right=0
  endif
  !FOR SPHERE 
  if (sph) then
     if (left.gt.0) then
        allocate(shpL(nL))
        call inp_rad(shpL,spec_scale)
        !typentry give the scale of input radiation
        call rdinps2(Equal,1,str,L,UCASE)
        if (str(1:L).eq.'FLUX') TypEntry(1) = 1
        if (str(1:L).eq.'LUM_R1') TypEntry(1) = 2
        if (str(1:L).eq.'ENERGY_DEN') TypEntry(1) = 3
        if (str(1:L).eq.'DILUTN_FAC') TypEntry(1) = 4
        if (str(1:L).eq.'T1') TypEntry(1) = 5
        !check if the entered value is acceptable
        if (TypEntry(1).lt.1.or.TypEntry(1).gt.5) then
           call msg(21)
           error = 1
           stop
        end if
        if (typentry(1).eq.2) then
           !enter luminosity [in Lo] of the source and distance r1[cm] to the source
           Lum = RDINP(Equal,1)
           dist = RDINP(Equal,1)
           !all units in dusty are in SI, so convert the input
           Lum = Lum*3.862d+26
           dist = dist/100.0d0
           Ji = Lum/(4.0d0*pi*dist*dist)/(4.0d0*pi)
!           Cr1 = dist
        endif
        if (typentry(1).eq.3) Ji = RDINP(Equal,1)
        if (typentry(1).eq.4) then
           !entry of dilution (normalization) factor
           dilutn = RDINP(Equal,1)
           Ji = dilutn*spec_scale/(4*pi)
        endif
        if (typentry(1).eq.5) then
           !enter dust temperature on inner boundary, T1[K]
           Tinner = RDINP(Equal,1)
        end if
     else
        typentry(1) = 0
     endif
     if (right.gt.0) then
        allocate(shpR(nL))
        call inp_rad(shpR,spec_scale)
        !typentry give the scale of input radiation
        call rdinps2(Equal,1,str,L,UCASE)
        if (left.eq.0) then
           if (str(1:L).eq.'ENERGY_DEN') typentry(2) = 3
           if (str(1:L).eq.'DILUTN_FAC') typentry(2) = 4
           !check if the entered value is acceptable
           if (typentry(2).lt.3.or.typentry(2).gt.4) then
              call msg(21)
              error = 1
              stop
           end if
           if (typentry(2).eq.3) Jo = RDINP(Equal,1)
           if (typentry(2).eq.4) then
              !entry of dilution (normalization) factor
              dilutn = RDINP(Equal,1)
              Jo = dilutn*spec_scale/(pi)
           endif
        else
           ! ksi is the relative bol.flux of the second source
           ksi = RDINP(Equal,1)
           if (ksi.lt.0.0) ksi = 0.0d0
           if (ksi.gt.1.0) ksi = 1.0d0
           write(12,'(a49,F5.2)') ' Relative bol.flux fraction of right source: R =',ksi
        end if
     end if
  end if
  !FOR SLAB
  if (slb) then
     if (left.gt.0) then
        allocate(shpL(nL))
        call inp_rad(shpL,spec_scale)
        ! typentry give the scale of input radiation
        call rdinps2(Equal,1,str,L,UCASE)
        if (str(1:L).eq.'FLUX') typentry(1) = 1
        if (str(1:L).eq.'ENERGY_DEN') typentry(1) = 3
        if (str(1:L).eq.'DILUTN_FAC') typentry(1) = 4
        if (str(1:L).eq.'T1') typentry(1) = 5
        if ((typentry(1).lt.3.or.typentry(1).gt.5).and.(typentry(1).ne.1)) then
           call msg(21)
           error = 1
           print*,'error msg(21) see out file'
           stop
        end if
        if (typentry(1).eq.1) Ji = RDINP(Equal,1) / 4.d0 / pi
        if (typentry(1).eq.3) Ji = RDINP(Equal,1)
        if (typentry(1).eq.4) then
           ! entry of dilution (normalization) factor
           dilutn = RDINP(Equal,1)
           Ji = dilutn*spec_scale/pi
        endif
        ! var3 = dilutn
!        if (typentry(1).eq.5) then
!           ! enter dust temperature on inner boundary, T1[K]
!           Tsub(1) = RDINP(Equal,1)
!        end if
        write(12,'(a33)') ' Calculation in planar geometry:'
        !find the kind of illumination
        call rdinps2(Equal,1,str,L,UCASE)
        if (str(1:L).eq.'DIRECTIONAL') then
           write(12,'(a41)') ' Directional illumination from the left.'
           ! enter incident theta_in:
           ! th1 the left illumination angle (in degrees) measured from the normal
           th1 = RDINP(Equal,1)
           call chkangle(th1)
           th1 = th1*pi/180.0d0
           mu1 = dcos(th1)
        elseif (str(1:L).eq.'ISOTROPIC') then
           th1 = -1.0d0
           write(12,'(a40)') ' Isotropic illumination from the left.'
           mu1 = -1.0d0
        end if
     endif
     if (right.gt.0) then
        allocate(shpR(nL))
        call inp_rad(shpR,spec_scale)
        call rdinps2(Equal,1,str,L,UCASE)
        if (str(1:L).eq.'DIRECTIONAL') then
           write(12,'(a41)') ' Directional illumination from the right.'
           ! enter incident theta_in:
           ! th1 the left illumination angle (in degrees) measured from the normal
           th2 = RDINP(Equal,1)
           call chkangle(th2)
           th2 = th2*pi/180.0d0
           mu2 = dcos(th2)
        elseif (str(1:L).eq.'ISOTROPIC') then
           th2 = -1.0d0
           write(12,'(a40)') ' Isotropic illumination from the right.'
           mu2 = -1.0d0
        end if
        ! ksi is the relative bol.flux of the second source
        ksi = RDINP(Equal,1)
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
  call rdinps2(Equal,1,str,L,UCASE)
  ifidG = 1
  if (str(1:L).eq.'COMMON_GRAIN') then
     top = 1
     nG = 6
  elseif (str(1:L).eq.'COMMON_GRAIN_COMPOSITE') then
     top = 1
     nG = 1
  elseif (str(1:L).eq.'COMMON_AND_ADDL_GRAIN') then
     top = 2
  elseif (str(1:L).eq.'COMMON_AND_ADDL_GRAIN_COMPOSITE') then
     top = 4
  elseif (str(1:L).eq.'TABULATED') then
     top = 3
     nG = 1
  end if
  if (top.ne.1.and.top.ne.2.and.top.ne.3.and.top.ne.4) then
     call msg(9)
     error = 1
     print*,'error msg(9)'
     stop
  end if
  ! For top.lt.3 read in abundances for supported grains
  if (top.lt.3) then
     xC(1) = RDINP(Equal,1)
     if (xC(1).lt.0.0d0) xC(1) = 0.0d0
     if ((xC(1).eq.0.0d0).and.(nG.gt.1)) nG = nG - 1
     sum = xC(1)
     do iG = 2, 7
        ! Special care to be taken of graphite (1/3-2/3 rule):
        if (iG.ne.5) then
           xC(iG) = RDINP(noEqual,1)
           if (xC(iG).lt.0.0d0) xC(iG) = 0.0d0
           if ((xC(iG).eq.0.0d0).and.(nG.gt.1)) nG = nG - 1
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
     if (iG.eq.1) write(stdf(iG),'(a)')"stnd_dust_lib/OssOdef.nk"
     if (iG.eq.2) write(stdf(iG),'(a)')"stnd_dust_lib/OssOrich.nk"
     if (iG.eq.3) write(stdf(iG),'(a)')"stnd_dust_lib/sil-dlee.nk"
     if (iG.eq.4) write(stdf(iG),'(a)')"stnd_dust_lib/gra-par-draine.nk"
     if (iG.eq.5) write(stdf(iG),'(a)')"stnd_dust_lib/gra-perp-draine.nk"
     if (iG.eq.6) write(stdf(iG),'(a)')"stnd_dust_lib/amC-hann.nk"
     if (iG.eq.7) write(stdf(iG),'(a)')"stnd_dust_lib/SiC-peg.nk"
  enddo
  ! user supplied n and k:
  if ((top.eq.2).or.(top.eq.4)) then
     nFiles = RDINP(Equal,1)
     if (top.eq.2) nG = nfiles + 7
     if (top.eq.4) nG = 1
     ! File names
     allocate(nameNK(nFiles))
     do iFiles = 1, nFiles
        call filemsg(nameNK(iFiles),'optical constants:')
     end do
     if(error.ne.0) goto 996
     ! Abundances
     xCuser(1) = RDINP(Equal,1)
     if (xCuser(1).lt.0.0d0) xCuser(1) = 0.0d0
     sum = sum + xCuser(1)
     if (nfiles.gt.1) then
        do iFiles = 2, nfiles
           xCuser(iFiles) = RDINP(noEqual,1)
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
        stop
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
  end if
  ! 2.1 Temperatures
  allocate(Tsub(nG))
  Tsub(1) = RDINP(Equal,1)
  print'(A,I3,A,F12.3)',' Grain:',1,' Sublimation Temperature:',Tsub(1)
  if (nG.gt.1) then
     do iG = 2, nG
        Tsub(iG) = RDINP(noEqual,1)
        if (Tsub(iG).gt.Tsub(ifidG)) ifidG = iG
        print*,'Grain:',iG,'Sublimation Temperature:',Tsub(iG)
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
        print*,'error msg(10)'
        stop
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
  allocate(SigmaA(nG+1,nL))
  allocate(SigmaS(nG+1,nL))
  call getOptPr(nameQ,nameNK,error,stdf,top,szds,qsd,a1,a2,nFiles,xC,XCuser)

  !=========  END READING DUST PROPERTIES ===================
  ! WriteOut prints all input data, read so far, in fname.out
  ! var1 is t1,fe1,luminosity or teff; var2 is r1; var3 is ext.rad. input
  ! call WriteOut(var1,var2,var3,nG,nameQ,nameNK)
!!$  ! (3) Density distribution
!!$  ! For sphere only:
!!$  if(sph) then
!!$     powd  = .false.
!!$     expd  = .false.
!!$     rdw   = .false.
!!$     rdwa  = .false.
!!$     fild  = .false.
!!$     rdwpr = .false.
!!$     ! Parameter describing eta function:
!!$     call rdinps2(Equal,1,str,L,UCASE)
!!$     if (str(1:L).eq.'POWD') then
!!$        powd = .true.
!!$        denstyp = 1
!!$     elseif (str(1:L).eq.'EXPD') then
!!$        expd = .true.
!!$        denstyp = 2
!!$        ! *** Winds ***
!!$        ! denstyp.eq.3 is RDW with default values of v1/ve=0.2, GravCor=0.5
!!$        ! denstyp.eq.6 is a private option with additional input for v1/ve and
!!$        ! GravCor=max(Fgrav/Frad);
!!$     elseif (str(1:L).eq.'RDW') then
!!$        rdw = .true.
!!$        denstyp = 3
!!$        ! analytical (gray) approximation for rdw
!!$     elseif (str(1:L).eq.'RDWA') then
!!$        rdwa = .true.
!!$        denstyp = 4
!!$        ! file with user supplied density distribution
!!$     elseif (str(1:L).eq.'USER_SUPPLIED') then
!!$        fild = .true.
!!$        denstyp = 5
!!$        ! private option for RDW with additional output
!!$     elseif (str(1:L).eq.'RDWPR') then
!!$        rdwpr = .true.
!!$        denstyp = 6
!!$     end if
!!$     ! initialize EtaOK and Ntr
!!$     EtaOK = 0
!!$     Ntr = 0
!!$     ! read parameters for each type of density distribution
!!$     ! smooth or broken power laws
!!$     if (powd) then
!!$        EtaOK = 1
!!$        Ntr = RDINP(Equal,1)
!!$        ! changed definition
!!$        Ntr = Ntr - 1
!!$        ! read in transition radii
!!$        if (Ntr.gt.0) then
!!$           Ytr(1) = RDINP(Equal,1)
!!$           if (Ntr.gt.1) then
!!$              do i = 2, Ntr
!!$                 Ytr(i) = RDINP(NoEqual,1)
!!$              end do
!!$           end if
!!$           Yout = RDINP(noEqual,1)
!!$        else
!!$           ! for smooth density power law
!!$           Yout = RDINP(Equal,1)
!!$        end if
!!$        if (Yout.le.1.0d0) Yout = 1.001d0
!!$        ! read in powers
!!$        pow = RDINP(Equal,1)
!!$        if (Ntr.gt.0) then
!!$           do i = 1, Ntr
!!$              ptr(i) = RDINP(NoEqual,1)
!!$           end do
!!$        end if
!!$        ! print info to the output file
!!$        if (Ntr.eq.0) then
!!$           call getfs(pow,2,0,strpow)
!!$           write(12,'(a38,a5)') ' density described by 1/r**k with k =',strpow
!!$           write(12,'(a21,1p,e10.3)')'  relative thickness:',Yout
!!$        else
!!$           write(12,*)' density described by a broken power law:'
!!$           write(12,*)'  power   Ytransition'
!!$           write(12,*)'  -------------------'
!!$           write(12,*)'              1.0'
!!$           call getfs(pow,2,0,strpow)
!!$           write(12,'(a2,a5)')'  ',strpow
!!$           do i = 1, Ntr
!!$              write(12,'(a10,1p,e10.3)')'          ',Ytr(i)
!!$              call getfs(ptr(i),2,0,strpow)
!!$              write(12,'(a2,a5)')'  ',strpow
!!$           end do
!!$           write(12,'(a10,1p,e10.3)')'          ',Yout
!!$        end if
!!$     end if
!!$     ! exponential law
!!$     if (expd) then
!!$        EtaOK = 1
!!$        Yout = RDINP(Equal,1)
!!$        if (Yout.le.1.0d0) Yout = 1.001d0
!!$        pow = RDINP(Equal,1)
!!$        if (pow.le.0.0d0) then
!!$           EtaOK = 0
!!$        else
!!$           write(12,*)' density described by exponential distribution'
!!$           write(12,'(a21,1p,e10.3)')'               sigma:',pow
!!$           write(12,'(a21,1p,e10.3)')'  relative thickness:',Yout
!!$        end if
!!$     end if
!!$     ! default approximation and default numerics for rad. driven winds
!!$     if (rdwa.or.rdw) then
!!$        EtaOK = 1
!!$        Yout = RDINP(Equal,1)
!!$        if (Yout.le.1.0d0) Yout = 1.001d0
!!$        ! ** default ** for epsilon = v1/ve = u1/ue:
!!$        pow = 0.2d0
!!$        if(rdw) then
!!$           ! ** default ** for max(gravcor = fgrav/frad_press):
!!$           ptr(1) = 0.5d0
!!$           ! convergence criterion:
!!$           ptr(2) = 1.0d0
!!$           ! default linear version of the eq. for velocity
!!$           ver = 1
!!$        end if
!!$        write(12,*)' Density for radiatively driven winds from'
!!$        if (rdwa) then
!!$           write(12,*)' Analytic approximation for gray dust.'
!!$        else
!!$           write(12,*)' Full dynamic calculation.'
!!$        end if
!!$        write(12,'(a21,1p,e10.3)')'  Relative thickness:',Yout
!!$     end if
!!$     ! full dynamical calculation for radiatively driven winds (private option)
!!$     ! the user can specify parameters that have default values in denstyp=3
!!$     ! user specified table for eta
!!$     if(fild) then
!!$        EtaOK = 1
!!$        strg = 'Dust density distribution:'
!!$        call filemsg(nameeta,strg)
!!$        write(12,*)' Density distribution supplied from file:'
!!$        write(12,'(2x,a100)') nameeta
!!$        call prHeader(3,nameeta)
!!$        ! read in the density
!!$        open(26,err=997,file=nameeta,status='old')
!!$        ! three lines in the header:
!!$        do i = 1, 3
!!$           read(26,*,err=997) strpow
!!$        end do
!!$        istop = 0
!!$        i = 0
!!$        do while (istop.ge.0)
!!$           read(26,*,end=900,err=997,iostat=istop) a, b
!!$           if (istop.ge.0) then
!!$              i = i + 1
!!$              xx(i) = a
!!$              e(i) = b
!!$              if (i.eq.1) x1 = xx(i)
!!$              yetaf(i) = xx(i) / x1
!!$           end if
!!$        end do
!!$900     close(26)
!!$        nYetaf = i
!!$        if (nYetaf.lt.2) goto 997
!!$        ! if input positions in descending order turn them around
!!$        if (yetaf(1).gt.yetaf(2)) then
!!$           do i = 1, nYetaf
!!$              aa(i) = yetaf(i)
!!$              bb(i) = e(i)
!!$           end do
!!$           do i = 1, nYetaf
!!$              yetaf(i) = aa(nYetaf+1-i)
!!$              e(i) = bb(nYetaf+1-i)
!!$           end do
!!$        end if
!!$        ! relative thickness
!!$        Yout = yetaf(nYetaf)
!!$        write(12,'(a21,1p,e10.3)')'  relative thickness:',Yout
!!$        if (Yout.le.1.0d0) Yout = 1.001d0
!!$        ! integrate and ...
!!$        call Simpson(Nmax,1,nYetaf,yetaf,e,ceta)
!!$        ! ... renormalize
!!$        do i = 1, nYetaf
!!$           etaf(i) = e(i) / ceta
!!$        end do
!!$     end if
!!$     ! Done with the reading of density distribution
!!$     if (EtaOK.ne.1) then
!!$        call msg(3)
!!$        error = 1
!!$        goto 999
!!$     end if
!!$     write(12,*)' --------------------------------------------'
!!$  end if
!!$  !=========  End reading density distribution ===================
!!$  ! 4) Optical depth
!!$  ! Grid type
!!$  call rdinps2(Equal,1,str,L,UCASE)
!!$  if (str(1:L).eq.'LINEAR') then
!!$     GridType = 1
!!$  elseif (str(1:L).eq.'LOGARITHMIC') then
!!$     GridType = 2
!!$  elseif (str(1:L).eq.'USER_SUPPLIED') then
!!$     GridType = 3
!!$  end if
!!$  if (GridType.eq.3) then
!!$     ! tau-grid from a file
!!$     strg = 'user supplied tau-grid:'
!!$     call filemsg(nametau,strg)
!!$     ! read optical depths
!!$     open(27,err=992,file=nametau,status='old')
!!$     ! fiducial wavelength
!!$     ! (the second argument of rdinp is the unit)
!!$     lamfid = RDINP(Equal,27)
!!$     ! number of models in the list
!!$     Nmodel = RDINP(Equal,27)
!!$     do i = 1, Nmodel
!!$        read(27,*) tauIn(i)
!!$     end do
!!$902  close(27)
!!$     ! Sort the tau-grid if there is more than one model:
!!$     if(Nmodel.gt.1) then
!!$        call sort(tauIn,Nmodel)
!!$     end if
!!$     tau1 = tauIn(1)
!!$     if (tau1.le.0.0d0) tau1 = 0.0001d0
!!$     tau2 = tauIn(Nmodel)
!!$  else
!!$     ! fiducial wavelength
!!$     lamfid = RDINP(Equal,1)
!!$     ! total optical depths at lamfid
!!$     TAU1 = RDINP(Equal,1)
!!$     if (tau1.le.0.0d0) tau1 = 0.0001d0
!!$     TAU2 = RDINP(Equal,1)
!!$     if (tau2.le.tau1) then
!!$        tau2 = tau1
!!$        Nmodel = 1
!!$     end if
!!$     ! read number of models
!!$     Nmodel = RDINP(Equal,1)
!!$     ! Nrec = 1000, initialized in MAIN
!!$     if (Nmodel.gt.(Nrec-1)) Nmodel = Nrec-1
!!$     if (Nmodel.lt.1) Nmodel = 1
!!$  end if
!!$  if (Nmodel.gt.1) then
!!$     write(12,'(a19,1p,e8.1,a8)')' Optical depths at',lamfid, 'microns'
!!$     write(12,'(a14,1p,e9.2,a3,e9.2)')' ranging from',tau1,' to',tau2
!!$     if (GridType.eq.1) strg=' models with linear grid    '
!!$     if (GridType.eq.2) strg=' models with logarithmic grid'
!!$     if (GridType.eq.3) strg=' models with grid from file  '
!!$     write(12,'(a1,i4,a)')' ', Nmodel, strg
!!$     if (GridType.eq.3) write(12,'(a4,a70)')'    ',nametau
!!$  else
!!$     write(12,'(a18,1p,e8.1,a9,e9.2)')' Optical depth at',lamfid, ' microns:',tau1
!!$  end if
!!$  ! 5) disk
!!$  ! for disk calculations make sure you have npX>1 in 'userpar.inc' !!
!!$  if (npX.gt.1) then
!!$     if (iVerb.ge.1) write(*,*) 'No disk option in this version.'
!!$     goto 999
!!$  end if
!!$  !********************************************
!!$  !** III. Numerical accuracy **
!!$  !********************************************
!!$  ! accuracy for convergence (typical 0.0001)
!!$  !!** this is accConv for dust temperature
!!$  !!**  accConv = 10.0d-3  !! Too rough.
!!$  !  accConv = 1.0d-4  !tests on July,8,2010  1e-4 and 1e-6 give very close results
!!$  ! accuracy for flux conservation
!!$  accuracy = RDINP(Equal,1)
!!$  accConv = RDINP(Equal,1)
!!$  accConv = ((1.+accuracy)**(1./4.)-1.)*0.1
!!$  if (accuracy.le.0.0d0) accuracy = 0.02d0
!!$  ! Protect against a very large value for accuracies
!!$  if (accuracy.gt.0.25d0) accuracy = 0.25d0
!!$  ! starting optical depth
!!$  init_tau = RDINP(Equal,1)
!!$  ! increment in optical depth
!!$  dtau = RDINP(Equal,1)
!!$  ! dynamical range
!!$  dynrange = 1.0d-15
!!$  if (accuracy.ge.0.1d0) then
!!$     call getfs(accuracy*1000.0d0,0,1,strpow)
!!$     write(12,'(a20,a3,a1)')' Required accuracy:',strpow,'%'
!!$  else
!!$     call getfs(accuracy*100.0d0,0,1,strpow)
!!$     write(12,'(a20,a2,a1)')' Required accuracy:',strpow,'%'
!!$  end if
!!$  !  write(12,*)' --------------------------------------------'
!!$  !********************************************
!!$  !** IV. Output flags **
!!$  !********************************************
!!$  ! Internal flag for additional miscellaneous output  [MN]:
!!$  ! if iInn=1: print err.vs.iter in unt=38 (fname.err) for all models
!!$  ! and additionally list scaled fbol(y) and ubol(y) in m-files.
!!$  iInn = 1
!!$  iPhys = 0
!!$  !  spectra
!!$  iA = RDINP(Equal,1)
!!$  iC = RDINP(Equal,1)
!!$  ! images (intensity)
!!$  if (iC.ne.0) then
!!$     if (slb) then
!!$        ! Read angular grid (this is theta_out) for slab intensity output.
!!$        ! the output intensities are in units of lambda*I_lambda*cos(theta_out)/Fe
!!$        ! where Fe=L/(4*pi*r^2), the local bolometric flux.
!!$        ang_type = RDINP(Equal,1)
!!$        !    Create the grid depending on grid type
!!$        !    1-equidistant in theta, 2-equidistant in cos(theta), 3-from a file
!!$        call input_slb_ang(ang_type)
!!$        !    Convert to radians
!!$        do imu = 1, nmu
!!$           theta(imu) = theta(imu)*pi/180.0d0
!!$        end do
!!$        iV = 0
!!$        iPsf = 0
!!$     else
!!$        ! for spherical case
!!$        NlambdaOut = RDINP(Equal,1)
!!$        if (nLambdaOut.ge.1) then
!!$           do i = 1, nLambdaOut
!!$              LambdaOut(i) = RDINP(NoEqual,1)
!!$              !  make sure the wavelengths are inside dusty's range
!!$              if (LambdaOut(i).le.0.01d0) LambdaOut(i) = 0.01d0
!!$              if (LambdaOut(i).gt.36000.0d0) LambdaOut(i) = 36000.0d0
!!$           end do
!!$           ioverflw = 0
!!$           do i = 1, nLambdaOut
!!$              if (LambdaOut(i).lt.0.995d0) then
!!$                 call getfs(LambdaOut(i),2,1,lamstr(i))
!!$              else
!!$                 if (LambdaOut(i).lt.9.95d0) then
!!$                    call getfs(LambdaOut(i),1,0,lamstr(i))
!!$                 else
!!$                    if (LambdaOut(i).lt.99.5d0) then
!!$                       call getfs(LambdaOut(i),0,0,lamstr(i))
!!$                    else
!!$                       call getfs(LambdaOut(i),0,1,lamstr(i))
!!$                    end if
!!$                 end if
!!$              end if
!!$              if (LambdaOut(i).gt.9999.5d0) then
!!$                 ioverflw = 1
!!$                 strpow = lamstr(i)
!!$                 strpow(4:4) = '*'
!!$                 strpow(5:5) = ' '
!!$                 lamstr(i) = strpow
!!$              end if
!!$           end do
!!$        end if
!!$        write(12,*)' Images requested for these wavelengths (mic)'
!!$        write(12,'(a1,20a5)')' ',(lamstr(i),i=1,nLambdaOut)
!!$        if (ioverflw.eq.1) write(12,*)'  *: in mm'
!!$        ! Convolved images  (only for our use)
!!$        if (iC.lt.0) then
!!$           iPsf = 1
!!$           ! iPsf = rdinp(Equal,1)
!!$           if (iPsf.ne.0) then
!!$              Theta1 = RDINP(Equal,1)
!!$              write(12,'(a39,1p,e7.1)') ' Convolved images produced for theta1=',theta1
!!$              psftype = RDINP(Equal,1)
!!$              if (psftype.ne.1.and.psftype.ne.2.and.psftype.ne.3) goto 994
!!$              if (psftype.lt.3) then
!!$                 ! Gaussians, read in parameters
!!$                 ! FWHM for the first component
!!$                 FWHM1(1) = RDINP(Equal,1)
!!$                 if (nLambdaOut.gt.1) then
!!$                    do i = 2, nLambdaOut
!!$                       FWHM1(i) = RDINP(NoEqual,1)
!!$                    end do
!!$                 end if
!!$                 if (psftype.eq.2) then
!!$                    ! Relative strength for the second component
!!$                    kPSF(1) = RDINP(Equal,1)
!!$                    if (nLambdaOut.gt.1) then
!!$                       do i = 2, nLambdaOut
!!$                          kPSF(i) = RDINP(NoEqual,1)
!!$                       end do
!!$                    end if
!!$                    !       FWHM for the second component
!!$                    FWHM2(1) = RDINP(Equal,1)
!!$                    if (nLambdaOut.gt.1) then
!!$                       do i = 2, nLambdaOut
!!$                          FWHM2(i) = RDINP(NoEqual,1)
!!$                       end do
!!$                    end if
!!$                 end if
!!$                 write(12,*)' the point spread functions are gaussians'
!!$              else
!!$                 ! user supplied psf
!!$                 strg = 'point spread function:'
!!$                 call filemsg(namepsf,strg)
!!$                 write(12,*)' the point spread function supplied from file'
!!$                 write(12,'(2x,a100)')namepsf
!!$                 open(28,err=995,file=namepsf,status='old')
!!$                 ! Three lines in the header:
!!$                 do i = 1, 3
!!$                    read(28,*,err=995)
!!$                 end do
!!$                 istop = 0
!!$                 i = 0
!!$                 do while (istop.ge.0)
!!$                    read(28,*,end=901,err=995,iostat=istop)a, b
!!$                    if (istop.ge.0) then
!!$                       i = i + 1
!!$                       if (i.eq.1) then
!!$                          psf1 = b
!!$                          if (a.ne.0.0d0) goto 995
!!$                       end if
!!$                       xpsf(i) = a
!!$                       ypsf(i) = b / psf1
!!$                    end if
!!$                 end do
!!$901              close(28)
!!$                 Npsf = i
!!$                 !      scale to 1 at the center. This is only to get FWHM here.
!!$                 !      ypsf is normalized to area in Subroutine Convolve [MN]
!!$                 call scaleto1(1000,npsf,ypsf)
!!$                 !      Find equivalent fwhm
!!$                 istop = 0
!!$                 i = 1
!!$                 do while (istop.eq.0)
!!$                    i = i + 1
!!$                    if (ypsf(i).le.0.5d0) istop = 1
!!$                 end do
!!$                 !      Linear interpolation
!!$                 FWHM1(1) = (xpsf(i)-xpsf(i-1))/(ypsf(i)-ypsf(i-1))
!!$                 FWHM1(1) = (fwhm1(1)*(0.5d0-ypsf(i-1))+xpsf(i-1))*2.0d0
!!$                 FWHM2(1) = 0.0d0
!!$                 write(12,'(a18,1p,e8.1)')' equivalent FWHM:',FWHM1(1)
!!$                 ! end if for psf from a file
!!$              end if
!!$              ! end if for psf
!!$           end if
!!$           ! end if for convolved images
!!$        end if
!!$        ! visibility (only if the intensity is requested)
!!$        iV = RDINP(Equal,1)
!!$        if(iV.ne.0) iV = abs(iC)
!!$        ! end if for geometry
!!$     end if
!!$     write(12,*)' --------------------------------------------'
!!$  else
!!$     ! if iC=0 set the other flags to 0 (just in case).
!!$     iPsf = 0
!!$     iV = 0
!!$     write(12,*)' --------------------------------------------'
!!$  end if
!!$
!!$  ! ---- added printout of lam*J_lam/J for sphere [MN'10] ------------
!!$  if(SPH) then
!!$     iJ = RDINP(Equal,1)
!!$     if(iJ.GT.0) then
!!$        nJOut = RDINP(Equal,1)
!!$        if (nJOut.ge.1) then
!!$           do i = 1, nJOut
!!$              YJOut(i) = RDINP(NoEqual,1)
!!$              !        make sure the radii are inside Dusty's range
!!$              if (YJOut(i).le.1.0) YJOut(i) = 1.0
!!$              if (YJOut(i).gt.Yout) YJOut(i) = Yout
!!$           end do
!!$        end if
!!$        write(12,*)' En.density profile requested for these y:'
!!$        write(12,'(a1,1p,10e12.3)')' ',(YJOut(i),i=1,nJOut)
!!$        write(12,*)' --------------------------------------------'
!!$     end if
!!$  end if
!!$  ! radial quantities
!!$  iB = RDINP(Equal,1)
!!$  ! run-time messages
!!$  iX = RDINP(Equal,1)
!!$  ! *** DONE READING INPUT PARAMETERS ***
!!$  ! if everything is ok, close the input file and finish
!!$999 goto 996
!!$  ! or in the case of err reading files...
!!$920 write(12,*)' ***  FATAL ERROR IN DUSTY  *************'
!!$    write(12,*)' File with user supplied angular grid:   '
!!$    write(12,*)'     slab_ang_grid.dat                   '
!!$    write(12,*)' is missing or not properly formatted?!  '
!!$    write(12,*)' ****************************************'
!!$    close(12)
!!$    error = 3
!!$992 write(12,*)' ***  FATAL ERROR IN DUSTY  *************'
!!$    write(12,*)' File with user supplied TAU-grid:       '
!!$    write(12,'(2x,a100)') nameTAU
!!$    write(12,*)' is missing or not properly formatted?!  '
!!$    write(12,*)' ****************************************'
!!$    !  close(12)
!!$    error = 3
!!$    goto 996
!!$994 call MSG(12)
!!$    !  close(12)
!!$    error = 3
!!$    goto 996
!!$995 write(12,*)' ***  FATAL ERROR IN DUSTY  *************'
!!$    write(12,*)' File with the point spread function:    '
!!$    write(12,'(a2,a100)')'  ', namePSF
!!$    write(12,*)' is missing or not properly formatted?!  '
!!$    write(12,*)' ****************************************'
!!$    !  close(12)
!!$    error = 3
!!$997 write(12,*)' ***  FATAL ERROR IN DUSTY  *************'
!!$    write(12,*)' File with the dust density distribution:'
!!$    write(12,'(2x,a100)') nameETA
!!$    write(12,*)' is missing or not properly formatted?!  '
!!$    write(12,*)' ****************************************'
!!$    !  close(12)
!!$    error = 3
!!$998 write(12,*)' ***  FATAL ERROR IN DUSTY  ****'
!!$    write(12,*)' Input file:'
!!$    write(12,'(2x,a100)') nameIn
!!$    write(12,*)' is missing?!'
!!$    write(12,*)' *******************************'
!!$    !  close(12)
!!$    error = 3
!!$    !-----------------------------------------------------------------------
996 close(1)
    return
end subroutine Input
!***********************************************************************

!***********************************************************************
subroutine inp_rad(shp,spec_scale)
!=======================================================================
! This is the former SUBROUTINE InpStar(error,is,nameIn)
! This subroutine is for reading the input radiation parameters
!                                                              [MN,Mar'99]
!=======================================================================
  use common
  implicit none

  integer i,iL,iLs,nLs,k, l, nBB,ios1,filetype,nLamtr,kstop
  double precision shp(:)
  double precision sum, tsum, value, RDINP, xSiO, bb, spec_scale, x, planck
  double precision a,b, EMfunc,fplbol
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
  call rdinps2(Equal,1,str,L,UCASE)
  ! help variable if input file is present
  filetype = 0
  if (str(1:L).eq.'BLACK_BODY') then 
     ! Number of black bodies
     nBB = RDINP(Equal,1)
     ! Stellar temperature(s)
     allocate(Tbb(nBB))
     allocate(rellum(nBB))
     Tbb(1) = RDINP(Equal,1)
     if (Tbb(1).le.0.0d0) then
        call msg(8)
        error = 1
        stop
     end if
     ! Single black body
     if (nbb.eq.1) then
        !relative luminosity
        rellum(1) = 1.0d0
     endif  !end if for one bb
     ! Multiple black bodies
     if (nbb.gt.1) then
        do i = 2, nbb
           Tbb(i) = RDINP(NoEqual,1)
           if (Tbb(i).le.0.0d0) then
              call msg(8)
              error = 1
              stop
           end if
        end do
        ! Read in relative luminosities
        rellum(1) = RDINP(Equal,1)
        sum = rellum(1)
        do i = 2, nbb
           rellum(i) = RDINP(NoEqual,1)
           sum = sum + rellum(i)
        end do
        if (sum.le.0.0d0) then
           call msg(7)
           error = 1
           stop
        end if
        ! Normalize
        tsum = 0.0d0
        do i = 1, nbb
           rellum(i) = rellum(i)/sum
           tsum = tsum + rellum(i)*Tbb(i)**(4.0d0)
        end do
     end if ! end if for multiple bb
     do iL = 1, nL
        bb = 0.0d0
        do k = 1, nbb
           x = 14400.0d0/(lambda(iL)*Tbb(k))  ! hc/k = 14400 micron.Kelvin
           bb = bb + rellum(k)*Planck(x)
        end do
        shp(iL) = bb
     end do
     deallocate(Tbb)
     deallocate(rellum)
  else if (str(1:L).eq.'ENGELKE_MARENGO') then
     ! Effective stellar temperature
     nBB = 1
     allocate(Tbb(nBB))
     allocate(rellum(nBB))
     Tbb(1) = RDINP(Equal,1)
     ! Depth of SiO abs.feature in %
     xSiO = RDINP(Equal,1)
     if (xSiO.le.0.0d0) xSiO = 0.0001d0
     if (xSiO.gt.100.0d0) xSiO = 100.0d0
     do iL=1,nL
        shp = EMfunc(lambda(iL),Tbb(1),xSiO)
     end do
  else if (str(1:L).eq.'POWER_LAW') then
     ! Number of transitions
     nLamTr= RDINP(Equal,1)
     allocate(lamTr(nLamTr))
     allocate(klam(nLamTr))
     if (nLamTr.gt.0) then
        lamtr(1) = RDINP(Equal,1)
        if (nLamTr.gt.1) then
           do i = 2, nLamtr+1
              lamtr(i) = RDINP(NoEqual,1)
              if (lamtr(i).lt.lamtr(i-1)) then
                 call msg(6)
                 error = 1
                 stop
              end if
           end do
        endif
        klam(1) = RDINP(Equal,1)
        if (nLamtr.gt.1) then
           do i = 2, nLamtr
              klam(i) = RDINP(NoEqual,1)
           end do
        end if
     else
        print*,'something wrong with powerlaw input'
        stop
     end if
     allocate(fl(nLamTr))
     allocate(fpl(nL))
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
     deallocate(fl)
     deallocate(fpl)
     deallocate(lamTr)
     deallocate(klam)
  else if (str(1:L).eq.'FILE_LAMBDA_F_LAMBDA') then 
     filetype = 1
  else if (str(1:L).eq.'FILE_F_LAMBDA') then 
     filetype = 2
  else if (str(1:L).eq.'FILE_F_NU') then 
     filetype = 3
  end if

  if (filetype.gt.0) then 
     strg = 'Spectral shape of external radiation:'
     call filemsg(filename,strg)
     open(3,file=filename,status='old')
     rewind(3)
     iLs=0
     ios1=0
     do while (ios1.ge.0)
        read(3,'(a235)',iostat=ios1) line
        iLs = iLs + 1
     end do
     ! 3 header lines and one line to much from loop
     nLs = iLs-1-3
     if (nLs.lt.2) then
        print*,'error while reading file',filename
        stop
     endif
     rewind(3)
     do iLs = 1, 3
        read(3,'(a235)') line
     end do
     ios1 = 0
     iLs = 0
     allocate(lambda_s(nLs))
     allocate(shp_s(nLs))
     allocate(tmp_sort1(nLs))
     allocate(tmp_sort2(nLs))
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
        call powerinter(nLs,nLs,lambda_s,shp_s,lambda(iL),iL,shp(iL))
     end do
     error = 0
     deallocate(lambda_s)
     deallocate(shp_s)
     deallocate(tmp_sort1)
     deallocate(tmp_sort2)
  end if
  ! get total scale
  call Simpson(nL,1,nL,lambda,shp,spec_scale)
  return
end subroutine inp_rad

!***********************************************************************
!!$
!!$
!!$!***********************************************************************
!!$subroutine input_slb_ang(ang_type)
!!$!=======================================================================
!!$! This subroutine reads the set of input or output illumination angles
!!$! for slab case.                                           [MN, 2005]
!!$! =======================================================================
!!$  use common
!!$  implicit none
!!$
!!$  integer ang_type, error, imu, length
!!$  double precision th_min, th_max, angstep, cth_min, cth_max, caux, value, RDINP
!!$  character*70 anggrid, strg,str*235
!!$  logical Equal, noEqual,UCASE
!!$!-----------------------------------------------------------------------
!!$
!!$  Equal = .true.
!!$  noEqual = .false.
!!$
!!$! if ang_type=1 (equidistant in theta, given min,max,step)
!!$  if (ang_type.eq.1) then
!!$    th_min = RDINP(Equal,1)
!!$    call chkangle(th_min)
!!$    th_max = RDINP(Equal,1)
!!$    call chkangle(th_max)
!!$    if (th_max.le.th_min) then
!!$     th_max = th_min
!!$     nmu = 1
!!$    end if
!!$! step equidistant in theta
!!$    angstep = RDINP(Equal,1)
!!$!   create the grid:
!!$    imu = 1
!!$    theta(1)=th_min
!!$    do while(theta(imu).lt.th_max)
!!$     theta(imu+1) = theta(imu) + angstep
!!$     imu = imu+1
!!$    end do
!!$    nmu = imu
!!$    theta(nmu) = th_max
!!$  end if
!!$
!!$! if ang_type=2 (equidistant in cos theta)
!!$  if (ang_type.eq.2) then
!!$    th_min = RDINP(Equal,1)
!!$    call chkangle(th_min)
!!$    th_max = RDINP(Equal,1)
!!$    call chkangle(th_max)
!!$    cth_min = dcos(th_min*pi/180.0d0)
!!$    cth_max = dcos(th_max*pi/180.0d0)
!!$!   Step, equidistant in cos(theta)
!!$    AngStep = RDINP(Equal,1)
!!$!   Create the grid:
!!$    imu = 1
!!$    caux = cth_min
!!$    do while (caux.gt.0.0d0)
!!$     theta(imu) = dacos(caux)*180.0d0/pi
!!$     caux = caux - angstep
!!$     imu = imu+1
!!$    end do
!!$    Nmu = imu
!!$    theta(Nmu) = th_max
!!$  end if
!!$
!!$  if (ang_type.eq.3) then
!!$! Angular grid from a file, angles in degrees
!!$   call FileMSG(ANGgrid,strg)
!!$   open(7,ERR=92,file=ANGgrid,STATUS='OLD')
!!$   Nmu = RDINP(Equal,7)
!!$   do imu = 1, Nmu
!!$     read(7,*) theta(imu)
!!$   end do
!!$92 close(7)
!!$  end if
!!$  write(12,*)' Intensity requested for these theta_out(deg):'
!!$  write(12,'(a1,8f7.1,/,x,8f7.1,/,x,10f7.1)')' ', (theta(imu), imu = 1, nmu)
!!$  if (ang_type.eq.1) write(12,'(a34,f4.1)') '  equidistant in theta_out, step=', angstep
!!$  if (ang_type.eq.2) write(12,'(a39,f4.1)') '  equidistant in cos(theta_out), step=', angstep
!!$  if (ang_type.eq.3) write(12,'(a18,a70)') '  grid from file: ', anggrid
!!$! -----------------------------------------------------------------------
!!$  return
!!$end subroutine input_slb_ang
!!$!***********************************************************************
!!$
!!$!***********************************************************************
!!$subroutine GetTau(nG,tau1,tau2,tauIn,Nrec,GridType,Nmodel,tau)
!!$!=======================================================================
!!$! This subroutine generates total optical depth TAUtot.
!!$!                                                      [Z.I., Mar. 1996]
!!$!=======================================================================
!!$  use common
!!$  implicit none
!!$  integer model, Nmodel, nG, iL, GridType, Nrec
!!$  double precision  q,TAUin(Nrec), tau(Nmodel), TAU1, TAU2
!!$!-----------------------------------------------------------------------
!!$  taumax = 0.0d0
!!$  do model = 1, Nmodel
!!$     if(GridType.eq.3) then
!!$        tau(model) = tauIn(model)
!!$     else
!!$        ! Calculate taufid for given model
!!$        if (model.eq.1) then
!!$           tau(model) = tau1
!!$        else
!!$           if (model.eq.Nmodel) then
!!$              tau(model) = tau2
!!$           else
!!$              if(GridType.eq.1) then
!!$                 q =  (tau2 - tau1)/(Nmodel-1.0d0)
!!$                 tau(model) = tau1 + q*(model-1.0d0)
!!$              else
!!$                 q = dexp(dlog(tau2/tau1)/(Nmodel-1.0d0))
!!$                 tau(model) = tau1*q**(model-1.0d0)
!!$              end if
!!$           end if
!!$        end if
!!$     end if
!!$  end do
!!$!-----------------------------------------------------------------------
!!$  return
!!$end subroutine GetTau
!!$!***********************************************************************
!!$
!!$!***********************************************************************
!!$subroutine getOmega(nG,omega)
!!$!=======================================================================
!!$! This subroutine generates albedo omega(iL,iY) from the abs/sca cross-
!!$! sections and the component abundancies. This is temporary (trivial)
!!$! version  for single size grains.                     [Z.I., Mar. 1996]
!!$!
!!$!!** Note that Omega(iG,iL) here is re-defined compared to the old Dusty. [MN]
!!$!=======================================================================
!!$  use common
!!$  implicit none
!!$  integer  iG, nG, iL, iY
!!$  double precision omega(npG+1,npL),ext(npL),sca(nPL)
!!$!-----------------------------------------------------------------------
!!$
!!$  ! generate overall albedo through the envelope
!!$  ! ** this is for future multigrain code **
!!$  ! ** for single grains it is trivial **
!!$  do iL = 1, nL
!!$     ext(iL) = 0
!!$     sca(iL) = 0
!!$     do iG = 1, nG
!!$        ! calculate albedo
!!$        ext(iL) = ext(iL) + (sigmaA(iG,iL) + sigmaS(iG,iL))
!!$        sca(iL) = sca(iL) + sigmaS(iG,iL)
!!$        omega(iG,iL) = sigmaS(iG,iL)/(sigmaA(iG,iL) + sigmaS(iG,iL))
!!$     end do
!!$     omega(nG+1,iL) = sca(iL) / ext(iL)
!!$  end do
!!$  ! calculate relative abundances
!!$  do iG = 1, nG
!!$     do iY = 1, nY
!!$        abund(iG,iY) = 1.0d0
!!$     end do
!!$  end do
!!$  !--------------------------------------------------------------------
!!$  return
!!$end subroutine getOmega
!!$!**********************************************************************
!!$
!!$!***********************************************************************
!!$subroutine getSpShape(shp,is)
!!$!=======================================================================
!!$! Produces source sp. shape. This was part of subroutine star, separated
!!$! here for clarity.  [MN]
!!$!=======================================================================
!!$
!!$  use common
!!$  implicit none
!!$
!!$  integer iY, iL, iLs, nLs, k, kstop, i, error, nLambdam, nis, is
!!$!  nLambdam is the max number entries for a user supplied stellar spectrum
!!$  parameter (nLambdam = 10000, nis = 2)
!!$  double precision lambdas(nLambdam), Llamstar(nLambdam), shp(npL), &
!!$       stellar(nLambdam), fl(100), fpl(npL), Lstar, EMfunc, bp, bb, &
!!$       x, Planck,fplbol
!!$!-----------------------------------------------------------------------
!!$
!!$  !set nLs=1 to avoid out of bounce
!!$  nLs = 1
!!$  ! is=1 for the enclosed source, is=2 for the external shell illumination
!!$  if ((startyp(is).ge.4).and.(startyp(is).le.6)) then
!!$     call readspectar(lambdas,Llamstar,Lstar,nLs,is,error)
!!$     ! Generate dimensionless stellar spectrum
!!$     do iLs = 1, nLs
!!$        stellar(iLs) = lambdas(iLs)*Llamstar(iLs)/Lstar
!!$     end do
!!$  else
!!$     ! if startyp.eq.3 generate power-law spectrum
!!$     if (startyp(is).eq.3) then
!!$        fl(1) = 1.0d0
!!$        if (nLamtr(1).gt.1) then
!!$           do i = 2, nLamtr(1)
!!$              fl(i) = fl(i-1)*(lamtr(is,i-1)/lamtr(is,i))**klam(is,i-1)
!!$           end do
!!$        end if
!!$        do iL = 1, nL
!!$           if ((lambda(iL)-lamtr(is,1))*(lambda(iL)-lamtr(is,nLamtr(1)+1)).le.0.0d0) then
!!$              kstop = 0
!!$              k = 0
!!$              do whiLe (kstop.eq.0)
!!$                 k = k + 1
!!$                 ! This is Matt's correction for reading more than one powers:
!!$                 if (lambda(iL).ge.lamtr(is,k).and. &
!!$                      lambda(iL).le.lamtr(is,k+1)) then
!!$                      kstop = 1
!!$                      fpl(iL) = fl(k)*(lamtr(is,k)/lambda(iL))**klam(is,k)
!!$                      fpl(iL) = fpl(iL)/lambda(iL)
!!$                 end if
!!$              end do
!!$           else
!!$              fpl(iL) = 0.0d0
!!$           end if
!!$        end do
!!$     end if
!!$  end if
!!$  call Simpson(npL,1,nL,lambda,fpl,fplbol)
!!$  ! loop over wavelengths
!!$  do iL = 1, nL
!!$     if (startyp(is).eq.1) then
!!$        bb = 0.0d0
!!$        do k = 1, nbb(is)
!!$           x = 14400.0d0/(lambda(iL)*Tbb(is,k))  ! hc/k = 14400 micron.Kelvin
!!$           bb = bb + rellum(is,k)*Planck(x)
!!$        end do
!!$     else if (startyp(is).eq.2) then
!!$        bb = EMfunc(lambda(iL),Tbb(is,1),xSiO)
!!$     else if (startyp(is).eq.3) then
!!$        bb = lambda(iL)*fpl(iL)/(fplbol)
!!$        !  for lambda longer than the longest entry in namestar
!!$        !  assume rayleigh-jeans tail
!!$     else if (lambda(iL).gt.lambdas(nLs)) then
!!$        bb = stellar(nLs) * (lambdas(nLs)/lambda(iL))**3.0d0
!!$     else if (lambda(iL).lt.lambdas(1)) then
!!$        ! if shorter than the shortest assume 0
!!$        bb = 0.0d0
!!$        !  if within limits interpolate and startyp(is).ge.4
!!$     else
!!$        call powerinter(nLambdam,nLs,lambdas,stellar,lambda(iL),iLs,bp)
!!$        bb = bp
!!$     end if
!!$     shp(iL) = bb
!!$  end do
!!$  error = 0
!!$  return
!!$  !---------------------------------------------------------------------
!!$end subroutine getSpShape
!!$!***********************************************************************
!!$
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
!!$
!!$!***********************************************************************
!!$subroutine PrOut(model,nG,delta)
!!$!=======================================================================
!!$! This subroutine prints the results out.        [ZI,Feb'96; MN,Mar'99]
!!$!=======================================================================
!!$
!!$  use common
!!$  implicit none
!!$  integer iG, iY, iL, i, model, j, unt, imu, nrows, ncols,nG, iOut, iNloc , iLV
!!$!  parameter (nrows=200, ncols=25)
!!$  double precision psfn, psffunc(20,1000),eta, faux(npL), omega(npG+1,npL),ftotL(npL),ftotR(npL),  &
!!$       tht1, xs, xds, xde, res, fnormL, fnormR, dmax, limval, GinfG1, delta, &
!!$       y_loc, J_loc, Jbol(10), FbolL, FbolR, FbolIL,FbolIR,xAttTotL,xAttTotR,xDsTotL,xDsTotR,xDeTotL,xDeTotR,temp1,temp2
!!$  double precision, allocatable::Elems(:,:)
!!$  character*120 STemp,Serr,hdint, hdcon,hdvis, s1, su1, s2, su2, tstr*10
!!$  character*132 hdsp1,hdsp2,hdrslb1,hdrslb2,hdrsph1,hdrsph2,hdrdyn
!!$  character*255 crossfilename
!!$  double precision sigmaVs,sigmaVa,sigmaVe
!!$!----------------------------------------------------------------------
!!$
!!$  dynrange = 1.d-15
!!$  if(allocated(Elems)) deallocate(Elems)
!!$  if (nG.gt.1) allocate(Elems(npL,3+2*nG))
!!$  if (nG.eq.1) allocate(Elems(npL,3))
!!$  call lininter(npL,nL,lambda,sigmaS(nG+1,:),lamfid,iLV,sigmaVs)
!!$  call lininter(npL,nL,lambda,sigmaA(nG+1,:),lamfid,iLV,sigmaVa)
!!$  Elems(:,1) = lambda(:)
!!$  Elems(:,2) = SigmaA(nG+1,:)/(sigmaVa+sigmaVs)
!!$  Elems(:,3) = SigmaS(nG+1,:)/(sigmaVa+sigmaVs)
!!$  write(855,'(A,f7.2,A,e12.5)') '# total extinction at',lamfid,' micron <V>=',sigmaVa + sigmaVs
!!$  hdrdyn = '#  lambda    <abs>/<V>  <sca>/<V>'
!!$  if (nG.gt.1) then
!!$     do iG=1,nG
!!$        hdrdyn(34+(iG-1)*22:34+(iG-0)*22) = '  <abs>/<V>  <sca>/<V>'
!!$        write(hdrdyn(38+(iG-1)*22:39+(iG-1)*22),'(i2.2)') iG
!!$        write(hdrdyn(49+(iG-1)*22:50+(iG-1)*22),'(i2.2)') iG
!!$        Elems(:,2+2*iG) = SigmaA(iG,:)/(sigmaVa+sigmaVs)
!!$        Elems(:,3+2*iG) = SigmaS(iG,:)/(sigmaVa+sigmaVs)
!!$     end do
!!$     write(855,'(A)') hdrdyn(:34+nG*22)
!!$     call maketable(Elems,npL,3+2*nG,855)
!!$  else
!!$     write(855,*) hdrdyn
!!$     call maketable(Elems,npL,3,855)
!!$  end if
!!$  close(855)
!!$
!!$  call Simpson(npL,1,nL,lambda,fsL(:,1)/lambda,FbolIL)
!!$  call Simpson(npL,1,nL,lambda,fsR(:,nY)/lambda,FbolIR)
!!$  FbolIL=FbolIL*Jext(1)
!!$  FbolIR=FbolIR*Jext(nY)
!!$  if(allocated(Elems)) deallocate(Elems)
!!$  allocate(Elems(npL,8))
!!$
!!$  !  find the bolometric fluxes at the boundaries [MN]
!!$  !** FH changed to find ftot everywhere
!!$  do iY = 1, nY
!!$     do iL = 1, nL
!!$        ! the emerging spectra for sphere (or right-side spectra for slab)
!!$        ftot(iL,iY) = fsL(iL,iY) + fde(iL,iY) + fds(iL,iY) - ksi*fsR(iL,iY)
!!$        if (abs(ftot(iL,iY)).lt.dynrange) ftot(iL,iY) = 0.
!!$        faux(iL) = ftot(iL,nY)/lambda(iL)
!!$     end do
!!$  enddo
!!$  do iL = 1,nL
!!$     ftotL(iL) = fde(iL,1) + fds(iL,1) - ksi*fsR(iL,1)
!!$     ftotR(iL) = fsL(iL,nY) + fde(iL,nY) + fds(iL,nY)
!!$     if (abs(ftotL(iL)).lt.dynrange) ftotL(iL) = 0.
!!$     if (abs(ftotR(iL)).lt.dynrange) ftotR(iL) = 0.
!!$  end do
!!$  call Simpson(npL,1,nL,lambda,faux,res)
!!$  call Simpson(npL,1,nL,lambda,ftotL/lambda,temp1)
!!$  call Simpson(npL,1,nL,lambda,fsR(:,1)/lambda(:),temp2)
!!$  xAttTotL = abs(temp2/temp1)
!!$  call Simpson(npL,1,nL,lambda,fde(:,1)/lambda(:),temp2)
!!$  xDeTotL = abs(temp2/temp1)
!!$  call Simpson(npL,1,nL,lambda,fds(:,1)/lambda(:),temp2)
!!$  xDsTotL = abs(temp2/temp1)
!!$  call Simpson(npL,1,nL,lambda,ftotR/lambda,temp1)
!!$  call Simpson(npL,1,nL,lambda,fsL(:,nY)/lambda(:),temp2)
!!$  xAttTotR = abs(temp2/temp1)
!!$  call Simpson(npL,1,nL,lambda,fde(:,nY)/lambda(:),temp2)
!!$  xDeTotR = abs(temp2/temp1)
!!$  call Simpson(npL,1,nL,lambda,fds(:,nY)/lambda(:),temp2)
!!$  xDsTotR = abs(temp2/temp1)
!!$
!!$  ! normalization factor for output spectra
!!$  call Simpson(npL,1,nL,lambda,ftotR/lambda,fnormR)
!!$  call Simpson(npL,1,nL,lambda,ftotL/lambda,fnormL)
!!$  ! the emerging bolometric flux
!!$  FbolR = fnormR * Jext(nY)
!!$  if (slb) FbolL = fnormL * Jext(1)
!!$
!!$  ! calculation of radiation pressure
!!$  ! nG + 1 contains the sum of all sigma(iG) 1<=iG<=nG
!!$  call lininter(npL,nL,lambda,sigmaS(nG+1,:),lamfid,iLV,sigmaVs)
!!$  call lininter(npL,nL,lambda,sigmaA(nG+1,:),lamfid,iLV,sigmaVa)
!!$  sigmaVe = sigmaVa+sigmaVs
!!$  do iY=1,nY
!!$     call Simpson(npL,1,nL,lambda,(sigmaS(nG+1,:)+sigmaA(nG+1,:))*ftot(:,iY)/lambda(:),temp1)
!!$     RPr(iY) = temp1/(4*pi*clight*mprot*Gconst)*1.0D4*3.84e26/1.988e30*5.0D-26*Jext(iY)/sigmaVe
!!$  enddo
!!$
!!$  res = 0.0d00
!!$! this is the cut-off for printout of small values (in spectra)
!!$  limval = 1.0d-20
!!$! Zeljko's calculation of theta1, the ang. size (in arcsec) of the cavity for Fbol=1e-6 W/m2
!!$  tht1 = 412.6d0/(dsqrt(Ji*4*pi))
!!$! error in %
!!$  if (SmC(5,model).lt.0.1d0) then
!!$   call getfs(SmC(5,model)*100.0d0,0,0,Serr)
!!$  else if (SmC(5,model).ge.0.1d0.and.SmC(5,model).lt.1.0d0) then
!!$   call getfs(SmC(5,model)*100.0d0,0,1,Serr)
!!$  else
!!$   call getfs(SmC(5,model)*100.0d0,0,2,Serr)
!!$  end if
!!$!--------------  overall parameters to *.out file -----------------------
!!$! write header to output file *.out
!!$  if (model.eq.1) then
!!$   write(12,*)'         '
!!$   write(12,*)' RESULTS:'
!!$   write(12,*)' --------'
!!$   if (slb) then
!!$!    slab output
!!$      s1=' ###   Tau0   Psi/Psi0    FiL     FiR      FbolL   FbolR    r1(cm)   TdL(K)   TdR(K)   RPr(1)  e(%)'
!!$     su1=' ###     1       2         3       4         5       6        7         8        9      10     11'
!!$     write(12,'(a)') s1
!!$     write(12,'(a)') su1
!!$     write(12,'(a)') &
!!$         ' ==================================================================================================='
!!$!  output for sphere
!!$   elseif(sph) then
!!$      s1= ' ###   tau0   Psi/Psi0 Fi(W/m2)  r1(cm)   r1/rc    theta1   T1(K)    Td(K)    RPr(1)  e(%)'
!!$     su1= ' ###     1       2        3        4        5        6        7        8        9      10'
!!$     if(rdwa.or.rdw) then
!!$      s2='  Mdot      Ve       M> '
!!$     su2='   11       12       13 '
!!$     write(12,'(a,a)') s1,s2
!!$     write(12,'(a,a)')su1,su2
!!$     write(12,'(a)') &
!!$     ' ====================================================================================================================='
!!$! **  private rdw file **
!!$      if (rdwpr) then
!!$       s1= '###   tau0      tauF     Mdot      Ve       M>       '
!!$       su1='###    1          2        3        4       5       6'
!!$       s2= 'Ginf/G1   P    delta  d/sqrt(w1)  winf     Phi    zeta(1)'
!!$       su2='        7        8        9        10       11       12'
!!$       write(66,'(a53,a57)') s1,s2
!!$       write(66,'(a53,a55)')su1,su2
!!$      end if
!!$     else
!!$      write(12,'(a)') s1
!!$      write(12,'(a)') su1
!!$      write(12,'(a)') &
!!$         ' ========================================================================================'
!!$     end if
!!$   end if !end if for sphere
!!$  end if
!!$! print output tables for ea.model
!!$!---------------- Output for slab: ---------------------------
!!$  if(slb) then
!!$     write(12,'(i4,1p,10e9.2,a3)') model, taufid, Psi/Psi0,FbolIL, FbolIR, FbolL, FbolR, Cr1, Td(1,1), Td(1,nY), RPr(1), Serr
!!$!---------- for spherical shell ------------------------------
!!$  elseif(sph) then
!!$     if (rdwa.or.rdw.or.rdwpr) then
!!$        write(12,'(i4,1p,9e9.2,a1,a3,a1,1p,3e9.2)') &
!!$             model, taufid, Psi/Psi0, Ji*4*pi, Cr1, r1rs, tht1, Td(1,1), Td(1,nY), RPr(1),' ',Serr,' ',CMdot, CVe, CM
!!$     else
!!$        write(12,'(i4,1p,9e9.2,a1,a3)') &
!!$             model, taufid, Psi/Psi0, Ji*4*pi, Cr1, r1rs, tht1, Td(1,1), Td(1,nY), RPr(1),' ',Serr
!!$     end if
!!$     if (rdwpr) then
!!$        !** private rdw file **
!!$        if (model.eq.1) then
!!$           write(66,'(a11,1p,e9.3,a10,1p,e9.3,a12,1p,e9.3,a13,e9.3)')  &
!!$                '###   qv = ',qv,', Qstar = ',Qstar,', v1/vinf = ',pow, &
!!$                ', (g/r)max = ',ptr(1)
!!$        end if
!!$        dmax = dsqrt(pow*winf)
!!$        if (G1.gt.0) then
!!$           GinfG1 =  Ginf / G1
!!$        else
!!$           GinfG1 = 0
!!$        end if
!!$        write(66,'(i4,1p,5e9.2,7e9.2)') &
!!$             model, taufid, tauF(nY), CMdot, CVe, CM, &
!!$             GinfG1, Prdw, delta, delta/dmax, winf, Phi, zeta1
!!$     end if
!!$     if(right.eq.0) then
!!$        if (startyp(1).eq.1.or.startyp(1).eq.2) then
!!$           if(Tstar(1).lt.Te_min) then
!!$              call getfs(Tstar(1),0,1,Tstr)
!!$              write(12,'(a50,a5,a5)') &
!!$                   ' ** WARNING: the input spectrum is a black-body at ',Tstr,' K **'
!!$              call getfs(Te_min,0,1,Tstr)
!!$              write(12,'(a50,a5,a5)') &
!!$                   ' *the point-source assumption requires min Teff of ',Tstr,' K **'
!!$           end if
!!$        end if
!!$     end if
!!$     ! end if for geometry
!!$  end if
!!$  !--------------   spectrum to *.s##  file   ------------------------
!!$  if (iA.ne.0) then
!!$   unt = 15
!!$   call line(1,2,unt)
!!$   if(slb) then
!!$    write(unt,'(a7,i3,a8,f8.3,a36)')'# model',model,' taufid=',taufid,'  spectrum from the right slab side'
!!$    write(unt,'(a13,1p,e9.2)') '# Fbol[W/m2]=',FbolR
!!$   else
!!$    write(unt,'(a7,i3,a8,f8.3,a10)') '# model',model,' taufid=',taufid,'  spectrum'
!!$    write(unt,'(a13,1p,e9.2)') '# Fbol[W/m2]=',FbolR
!!$   end if
!!$   call line(1,1,unt)
!!$   call getOmega(nG,omega)
!!$   do iL = 1, nL
!!$    if (ftot(iL,nY).ne.0.0d0) then
!!$     xs = fsL(iL,nY)/ftotR(iL)
!!$     xds = fds(iL,nY)/ftotR(iL)
!!$     xde = fde(iL,nY)/ftotR(iL)
!!$    else
!!$     xs = 0.0d0
!!$     xds = 0.0d0
!!$     xde = 0.0d0
!!$    end if
!!$!  no need to print negligible values
!!$    if (dabs(xs).lt.limval) xs = 0.0d0
!!$    if (dabs(xds).lt.limval) xds = 0.0d0
!!$    if (dabs(xde).lt.limval) xde = 0.0d0
!!$    if (dabs(fsL(iL,1)).lt.limval) fsL(iL,1) = 0.0d0
!!$!   Printing normalized spectral shapes. Bol. flux values are in the headers. [MN]
!!$    if (dabs(ftot(iL,nY)).lt.limval) ftot(iL,nY) = 0.0d0
!!$    Elems(iL,1) = lambda(iL)
!!$    Elems(iL,2) = ftotR(iL)/fnormR
!!$    Elems(iL,3) = xs
!!$    Elems(iL,4) = xds
!!$    Elems(iL,5) = xde
!!$    Elems(iL,6) = fsL(iL,1)/fsLbol(1)
!!$    Elems(iL,7) = tautot(iL)
!!$    Elems(iL,8) = omega(nG+1,iL)
!!$   end do
!!$!------ tabulate the spectra in the desired form ----------
!!$   if(slb) then
!!$    hdsp1 = '#   lambda     fRight     xAtt       xDs        xDe        fInp_L     TauTot     albedo'
!!$   else
!!$    hdsp1 = '#   lambda     fTot       xAtt       xDs        xDe        fInp       TauTot     albedo'
!!$   end if
!!$   write(unt,'(A90)') hdsp1
!!$   if(slb) then
!!$      write(hdsp1,'(A,1p,E10.3,A,E10.3,A,E10.3,A,E10.3,A,E10.3,A)')  '    -1      ',FbolR,' ',xAttTotR,' ',&
!!$           xDsTotR,' ',xDeTotR,' ',FbolIL,'     -1.       -1.'
!!$   else
!!$      write(hdsp1,'(A,E9.3,A,E9.3,A)')  '   -1.       ',FbolR,'                                   ',FbolIL,'      -1.         -1.'
!!$   end if
!!$   write(unt,'(A90)') hdsp1
!!$   call maketable(Elems,npL,8,unt)
!!$!  spectra from the left (illuminated) slab side (file *.z##)
!!$   if (slb) then
!!$    call getOmega(nG,omega)
!!$    do iL = 1, nL
!!$       if (ftot(iL,1).ne.0.0d0) then
!!$          xs =  fsR(iL,1)/ftotL(iL)
!!$          xds = fds(iL,1)/ftotL(iL)
!!$          xde = fde(iL,1)/ftotL(iL)
!!$       else
!!$          xs = 0.0d0
!!$          xds = 0.0d0
!!$          xde = 0.0d0
!!$       end if
!!$       if (dabs(xs).lt.limval) xs =0.0d0
!!$       if (dabs(xds).lt.limval) xds =0.0d0
!!$       if (dabs(xde).lt.limval) xde =0.0d0
!!$       if (dabs(fsR(iL,nY)).lt.limval) fsR(iL,nY) = 0.0d0
!!$       ! rescale ftot with the bolom flux for z-spectra
!!$       if (dabs(ftot(iL,1)).lt.limval) ftot(iL,1) = 0.0d0
!!$       Elems(iL,1) = lambda(iL)
!!$       Elems(iL,2) = ftotL(iL)/fnormL
!!$       Elems(iL,3) = xs
!!$       Elems(iL,4) = xds
!!$       Elems(iL,5) = xde
!!$       if (ksi.gt.0) then
!!$          Elems(iL,6) = ksi*fsR(iL,nY)/fsRbol(nY)
!!$       else
!!$          Elems(iL,6) = 0.0d0
!!$       end if
!!$       Elems(iL,7) = tautot(iL)
!!$       Elems(iL,8) = omega(nG+1,iL)
!!$    end do
!!$    if (iA.eq.3) unt=25
!!$    ! append to the .s## file or write in a separate .z## file (if iA=3)
!!$    call line(1,1,unt)
!!$    write(unt,'(a7,i3,a8,f8.3,a33)') '# model',model,' taufid=',taufid,' spectrum from the left slab side'
!!$    write(unt,'(a13,1p,e9.2)') '# Fbol[W/m2]=',FbolL
!!$    call line(1,1,unt)
!!$    write(unt,'(a)')'#   lambda     fLeft      xAtt       xDs        xDe        fInp_R     TauTot     albedo'
!!$    write(hdsp1,'(A,1p,E10.3,A,E10.3,A,E10.3,A,E10.3,A,E10.3,A)')  '    -1      ',&
!!$         FbolL,' ',xAttTotL,' ',xDsTotL,' ',xDeTotL,' ',FbolIR,'     -1.       -1.'
!!$    write(unt,'(A)') hdsp1
!!$    call maketable(Elems,nL,8,unt)
!!$   end if
!!$  end if
!!$
!!$!-----------  radial quantities to *.r## (old *.bxx) file -------------
!!$  if (iB.ne.0) then
!!$   if(allocated(Elems)) deallocate(Elems)
!!$   allocate(Elems(nY,9+nG))
!!$   hdrslb1= '#     t        epsilon     tauF       RPr   '
!!$   do iG = 1,nG
!!$      write(hdrslb2(1+(iG-1)*11:1+(iG)*11),'(a7,i2.2,a2)'), '    Td(',iG,') '
!!$   end do
!!$   hdrsph1= '#     y         Td         eta         t '
!!$   hdrsph2= '     tauF      epsilon        RPr'
!!$   hdrdyn= '         u        drift'
!!$   unt = 16
!!$   call line(1,2,unt)
!!$   if(slb) then
!!$    write(unt,'(a7,i3,a8,f8.3,a18)') '# model',model,' taufid=',taufid,'  spatial profiles'
!!$   else
!!$    write(unt,'(a7,i3,a8,f8.3,a18)') '# model',model,' taufid=',taufid,'  radial profiles '
!!$   end if
!!$   call line(1,1,unt)
!!$!--------- for slab ---------
!!$   if (slb) then
!!$    do iY = 1, nY
!!$     Elems(iY,1) = tr(iY)
!!$     Elems(iY,2) = eps(iY)
!!$     Elems(iY,3) = tauF(iY)
!!$     Elems(iY,4) = RPr(iY)/RPr(1)
!!$     do iG=1,nG
!!$        Elems(iY,4+iG) = Td(iG,iY)
!!$     end do
!!$    end do
!!$    write(unt,'(a43,a)') hdrslb1,hdrslb2
!!$    call maketable(Elems,nY,4+nG,unt)
!!$!------  for spherical shell --------
!!$   elseif(sph) then
!!$    do iY = 1, nY
!!$     Elems(iY,1) = Y(iY)
!!$     Elems(iY,2) = Td(1,iY)
!!$     Elems(iY,3) = eta(Y(iY))
!!$     Elems(iY,4) = tr(iY)
!!$     Elems(iY,5) = tauF(iY)
!!$     Elems(iY,6) = eps(iY)
!!$     if (RPr(1).ne.0) Elems(iY,7) = RPr(iY)/RPr(1)
!!$!     Elems(iY,8) = rg(1,iY)*Jext(iY)
!!$!     if (rdwpr) then
!!$! redefine for private rdw (denstyp.eq.6) option
!!$!      Elems(iY,8) = rg(1,iY)
!!$!      Elems(iY,9) = gamma(iY)
!!$!      Elems(iY,10) = qF(iY)
!!$!     end if
!!$    end do
!!$! check values:
!!$    do i = 1, 7
!!$     do iY = 1, nY
!!$      if(Elems(iY,i).lt.limval) Elems(iY,i) = 0.0d0
!!$     end do
!!$    end do
!!$! with dynamics
!!$    if (rdw) then
!!$     do iY = 1, nY
!!$        Elems(iY,8) = ugas(iY)/ugas(nY)
!!$        Elems(iY,9) = vrat(1,iY)
!!$     end do
!!$! check values:
!!$     do i = 8, 9
!!$      do iY = 1, nY
!!$       if(Elems(iY,i).lt.limval) Elems(iY,i) = 0.0d0
!!$      end do
!!$     end do
!!$     write(unt,'(a42,a42,a23)') hdrsph1,hdrsph2,hdrdyn
!!$     call maketable(Elems,nY,9,unt)
!!$    else
!!$     write(unt,'(a42,a42)') hdrsph1,hdrsph2
!!$     call maketable(Elems,nY,7,unt)
!!$    end if
!!$! end if for geometry
!!$   end if
!!$! end if for the iB (radial) flag
!!$  end if
!!$
!!$!--------------   intensities to *.inn (old *.cxx) file  --------------
!!$  if (abs(iC).ne.0) then
!!$! slab intensity (found at the end of subroutine slbradt)
!!$! theta(nmu) are the angles of output intensities
!!$     if (slb) then
!!$        if(allocated(Elems)) deallocate(Elems)
!!$        allocate(Elems(npL,nmu+2))
!!$        hdint = '   lambda'
!!$        unt = 17
!!$        call line(1,2,unt)
!!$        write(unt,'(a7,i3,a8,f8.3,a32)')'# model',model,' taufid=',taufid,' transmitted i(theta)*cos(theta)'
!!$        call line(1,1,unt)
!!$        do iL = 1, nL
!!$           Elems(iL,1) = lambda(iL)
!!$           do imu = 1, nmu
!!$              !if(iPhys.eq.1) SLBintm(imu,iL) = SLBintm(imu,iL)*Jext(nY)
!!$              !4pi comes from slbintp since it is divided by 4pi need to be changed!!
!!$              Elems(iL,imu+1) = SLBintm(imu,iL)*Jext(1)*4*pi
!!$           end do
!!$           Elems(iL,nmu+2) = istR(iL)
!!$        end do
!!$        ! write(unt,'(a9,21f11.3)')hdint,(theta(imu),imu=1,nmu)
!!$        ! printout angles in degrees
!!$        ! write(unt,'(a9,37f11.1,a9)') hdint,
!!$        ! &                    (theta(imu)*180.0d0/pi,imu=1,nmu),'     IstR'
!!$        write(unt,'(a9,100f11.1)') hdint,(theta(imu)*180.0d0/pi,imu=1,nmu)
!!$        call maketable(Elems,npL,nmu+1,unt)
!!$        ! adding the column with stellar ints at the end of the table
!!$        !  call maketable(Elems,nL,nmu+2,unt)
!!$        hdint = '   lambda'
!!$        unt = 17
!!$        call line(1,2,unt)
!!$        write(unt,'(a7,i3,a8,f8.3,a32)')'# model',model,' taufid=',taufid,' reflected cos(theta)*i(theta)'
!!$        call line(1,1,unt)
!!$        do iL = 1, nL
!!$           Elems(iL,1) = lambda(iL)
!!$           do imu = 1, nmu
!!$!              if(iPhys.eq.1) SLBintp(imu,iL) = SLBintp(imu,iL)*Jext(1)
!!$               !4pi comes from slbintp since it is divided by 4pi need to be changed!!
!!$
!!$              Elems(iL,imu+1) = SLBintp(imu,iL)*Jext(nY)*4*pi !4pi comes from slbintp
!!$           end do
!!$        end do
!!$        !write(unt,'(a9,21f11.3)')hdint,(theta(imu),imu=1,nmu)
!!$        !printout angles in degrees
!!$        write(unt,'(a9,99f11.1)')hdint,(theta(imu)*180.0d0/pi,imu=1,nmu)
!!$        call maketable(Elems,npL,nmu+1,unt)
!!$     !------  for spherical shell --------
!!$     elseif(sph) then
!!$        if(allocated(Elems)) deallocate(Elems)
!!$        allocate(Elems(np+2,nLambdaOut+2))
!!$        hdint = '#     b          t(b)'
!!$        hdcon = '#   offset '
!!$        hdvis = '#     q    '
!!$        unt = 17
!!$        call line(1,2,unt)
!!$        write(unt,'(a7,i3,a8,f8.3,a14)') '# model',model,' taufid=',taufid,'   raw image  '
!!$        call line(1,1,unt)
!!$        do i = 1, nP+2
!!$           Elems(i,1) = bOut(i)
!!$           Elems(i,2) = tauZout(i)
!!$           do j = 1, nLambdaOut
!!$              ! check values:
!!$              if(IntOut(j,i).ne.IntOut(j,i).or.IntOut(j,i).lt.limval) then
!!$                 IntOut(j,i) = 0.0d0
!!$              end if
!!$              Elems(i,j+2) = IntOut(j,i)
!!$              ! we want intensity in Jy/arcsec^2
!!$              ! this was the bug in intensity output for sphere,
!!$              ! the missing 4piY^2 factor for intensity output [June 2006]
!!$              ! Elems(i,j+2) = 7.83 * LambdaOut(j) * Fi * Elems(i,j+2)
!!$              IF (iPhys.eq.1) THEN
!!$                 Elems(i,j+2) = 7.834d0*LambdaOut(j)*(Jext(nY)*4.0d0*pi*Yout**2.0d0)*Elems(i,j+2)
!!$              ELSE
!!$                 Elems(i,j+2) = 7.834d0*LambdaOut(j)*(4.0d0*pi*Yout**2.0d0)*Elems(i,j+2)
!!$              END IF
!!$           end do
!!$        end do
!!$        write(unt,'(a21,20f11.2)')hdint,(LambdaOut(j),j=1,nLambdaOut)
!!$        call maketable(Elems,nP+2,nLambdaOut+2,unt)
!!$     end if
!!$  end if
!!$  if (iC.lt.0) then
!!$     !---------  convolved images either add to .i## file or write in *.c## file --
!!$     if(iC.eq.-3) unt = 21
!!$     call line(1,2,unt)
!!$     write(unt,'(a7,i3,a8,f8.3,a20)') '# model',model,' taufid=',taufid,'   convolved image  '
!!$     call line(1,1,unt)
!!$     if(allocated(Elems)) deallocate(Elems)
!!$     allocate(Elems(nconv,nLambdaOut+1))
!!$     do i = 1, nconv
!!$        Elems(i,1) = offset(i)
!!$        do j = 1, nLambdaOut
!!$           if(convint(j,i).lt.limval) convint(j,i) = 0.0d0
!!$           Elems(i,j+1) = convint(j,i)
!!$        end do
!!$     end do
!!$     write(unt,'(a11,20f11.2)')hdcon,(LambdaOut(i),i=1,nLambdaOut)
!!$     call maketable(Elems,nconv,nLambdaOut+1,unt)
!!$     if (psftype.lt.3.and.model.eq.1) then
!!$        ! wavelength dependent psfs, print them separately in *.psf
!!$        ! first generate wavelength dependent psfs
!!$        do j = 1, nLambdaOut
!!$           iLambda = j
!!$           ! added dec.04 [mn]
!!$           do i = 1, nconv
!!$              psffunc(j,i) = psfn(offset(i))
!!$              ! norm.needs to be done here again (after call to psfn)
!!$              psffunc(j,i) = psffunc(j,i)/psfarea(j)
!!$              ! check dynamic range
!!$              call chkrange(dynrange,psffunc(j,i))
!!$           end do
!!$        end do
!!$        ! print them out
!!$        do i = 1, nconv
!!$           write(23,'(1p,e12.5,20e10.3)')offset(i),(psffunc(j,i),j=1,nLambdaOut)
!!$        end do
!!$     end if
!!$  end if
!!$  !--------------  visibility curves to *.vnn file    ------------------------
!!$  if (sph) then
!!$   if (iV.ne.0) then
!!$    if(abs(iC).eq.3) unt = 22
!!$    call line(1,2,unt)
!!$    write(unt,'(a7,i3,a8,f8.3,a14)') '# model',model, &
!!$         ' taufid=',taufid,'  visibility  '
!!$    call line(1,1,unt)
!!$
!!$    if(allocated(Elems)) deallocate(Elems)
!!$    allocate(Elems(nvisi,nLambdaOut+1))
!!$
!!$    do i = 1, nvisi
!!$     Elems(i,1) = qtheta1(i)
!!$     do j = 1, nLambdaOut
!!$      if(visib(j,i).lt.limval) visib(j,i) = 0.0d0
!!$      Elems(i,j+1) = visib(j,i)
!!$     end do
!!$    end do
!!$    write(unt,'(a11,20f11.2)')hdvis,(LambdaOut(i),i=1,nLambdaOut)
!!$    call maketable(Elems,nvisi,nLambdaOut+1,unt)
!!$   end if
!!$  endif
!!$!----------  energy density profiles to *jnn file  -------------------------
!!$  if(sph) then
!!$   if(allocated(Elems)) deallocate(Elems)
!!$   allocate(Elems(nL,nJOut+1))
!!$
!!$   if(iJ.gt.0) then
!!$    unt = 19
!!$    call line(1,2,unt)
!!$    write(unt,'(a7,i3,a8,f8.3,a20)') '# model',model,' taufid=',taufid,'   Energy density   '
!!$    call line(1,1,unt)
!!$!   Jext is found in the beginning of sub PrOut.
!!$!   Find the scale of en. density in [W/m2] for the required YJout(iOut) [MN]
!!$    do iOut = 1, nJOut
!!$      y_loc = YJOut(iOut)
!!$      call LININTER(npY,nY,Y,Jext,y_loc,iNloc,J_loc)
!!$      Jbol(iOut) = J_loc
!!$    end do
!!$!   write the scale in the header
!!$    write(unt,'(a11,1p,10e11.3)')'Jbol[W/m2]=',(Jbol(iOut),iOut=1,nJout)
!!$    write(unt,'(a11,10f11.2)')'      Y =  ',(YJOut(iOut),iOut=1,nJout)
!!$    write(unt,'(a11)')'   lambda  '
!!$
!!$!   normalize en. density profiles
!!$    do iOut = 1, nJOut
!!$      do iL = 1, nL
!!$        if(JOut(iL,iOut).lt.limval) JOut(iL,iOut) = 0.0d0
!!$        faux(iL) = JOut(iL,iOut) / lambda(iL)
!!$      end do
!!$      call Simpson(npL,1,nL,lambda,faux,res)
!!$      fnormR = res
!!$      do iL = 1, nL
!!$        Elems(iL,1) = lambda(iL)
!!$        Elems(iL,iOut+1) = JOut(iL,iOut) / fnormR
!!$      end do
!!$    end do
!!$
!!$    call maketable(Elems,nL,nJout+1,unt)
!!$   end if
!!$  end if
!!$
!!$ if(allocated(Elems)) deallocate(Elems)
!!$!-----------------------------------------------------------------------
!!$
!!$ return
!!$end subroutine PrOut
!!$!***********************************************************************
!!$
!!$
!!$!***********************************************************************
!!$subroutine CLLOSE(error,model,Nmodel)
!!$! =======================================================================
!!$! This subroutine closes output files.             [ZI,Feb'96; MN,Apr'99]
!!$! =======================================================================
!!$  use common
!!$  implicit none
!!$  character*235 su1, su2, s3, s4, txtt, txtf
!!$  integer  error, model, Nmodel, im
!!$! -----------------------------------------------------------------------
!!$  ! close the default output file:
!!$  if (iError.ne.0) then
!!$     write(12,'(a42,i4)') ' There are some error messages for model:',model
!!$     write(12,*) ' Please check m## file (if not produced then rerun)'
!!$     if (iverb.gt.0) print*, 'There are some error messages for model:',model,' (check m## file)'
!!$  end if
!!$  if (iWarning.ne.0.and.iError.eq.0) then
!!$     write(12,'(a36,i4)') ' There are some warnings for model:',model
!!$     write(12,*)' Please check m## file (if not produced then rerun)'
!!$     if (iverb.gt.0) print*, 'There are some warnings for model:',model,' (check m## file)'
!!$  end if
!!$  iCumm = iCumm + iError + iWarning
!!$  if (model.eq.Nmodel.or.error.eq.3.or.error.eq.4) then
!!$     if (error.ne.3) then
!!$        if(rdw.or.rdwa.or.rdwpr) then
!!$           write(12,'(a)')  ' ====================================='
!!$        else
!!$           if(sph) then
!!$              write(12,'(a)')  ' ====================================='
!!$           else
!!$              write(12,'(a)')  ' ====================================='
!!$           end if
!!$        end if
!!$        write(12,'(a22,1p,e8.1,a8)')'  (1) optical depth at',lamfid,' microns'
!!$        write(12,'(a69,1p,e9.2)')'  (2) Psi as defined by eq.14 in IE97 with optically thin value Psi0=', Psi0
!!$        if(slb) then
!!$           !  ----------  for slab output ----------------------------
!!$           write(12,'(a)')'  (3) input bol.flux (in W/m2) of the left-side source at the left slab boundary'
!!$           write(12,'(a)')'  (4) input bol.flux (in W/m2) of the right-side source at the right slab boundary'
!!$           write(12,'(a)')'  (5) bolometric flux (in W/m2) at the left slab boundary'
!!$           write(12,'(a)')'  (6) bolometric flux (in W/m2) at the right slab boundary'
!!$           write(12,'(a)')'  (7) position of the left slab boundary for L = 1e4 Lo'
!!$           write(12,'(a)')'  (8) dust temperature at the left slab boundary'
!!$           write(12,'(a)')'  (9) dust temperature at the right slab boundary'
!!$           write(12,'(a)')'  (10) radiation pressure on left boundary; see manual for units'
!!$           write(12,'(a)')'  (11) maximum error in flux conservation (Fmax-Fmin)/(Fmax+Fmin)'
!!$        else
!!$           !---------- for spherical shell ----------------------------
!!$           write(12,'(a)')'  (3) bolometric flux at the inner radius '
!!$           write(12,'(a)')'  (4) inner radius for L = 1e4 Lo'
!!$           write(12,'(a)')'  (5) ratio of the inner to the stellar radius'
!!$           write(12,'(a)')'  (6) angular size (in arcsec) when Fbol=1e-6 W/m2'
!!$           write(12,'(a)')'  (7) dust temperature at the inner radius '
!!$           write(12,'(a)')'  (8) dust temperature at the outer edge'
!!$           write(12,'(a)')'  (9) radiation pressure on inner boundary; see manual for units'
!!$           write(12,'(a)')' (10) maximum error in flux conservation (Fmax-Fmin)/(Fmax+Fmin)'
!!$           if(rdw.or.rdwa.or.rdwpr) then
!!$              write(12,'(a)')' (11) mass-loss rate (in Mo/year)'
!!$              write(12,'(a)')' (12) terminal outflow velocity (in km/s)'
!!$              write(12,'(a)')' (13) upper limit of the stellar mass (in Mo)'
!!$           end if
!!$        end if
!!$        write(12,*)' ================================================'
!!$        if(iCumm.eq.0) write(12,*)' Everything is ok for all models.'
!!$        !---- private file ---
!!$        if(rdwpr) write(12,*)' Table with the wind properties is in file *.rdw'
!!$        !------------------------- spectra ------------------
!!$        if(iA.ne.0) then
!!$           if (iA.eq.1) then
!!$              write(12,*)' All spectra are in file *.stb'
!!$           else
!!$              if (slb.and.iA.eq.3) then
!!$                 write(12,*)' Spectra are in files *.s## and *.z##'
!!$              else
!!$                 write(12,*)' Spectra are in files *.s##'
!!$              end if
!!$           end if
!!$        end if
!!$        !------------------------- images ------------------
!!$        if (iC.ne.0) then
!!$           if (abs(iC).eq.1) then
!!$              if(slb) then
!!$                 write(12,*)' All intensities are in file *.itb'
!!$              else
!!$                 write(12,*)' All imaging quantities are in file *.itb'
!!$              end if
!!$           else
!!$              if(slb) then
!!$                 write(12,*)' Intensities are in files *.i##'
!!$              else
!!$                 write(12,*)' Images are in files *.i##'
!!$              end if
!!$              if (iV.ne.0.and.abs(ic).eq.3) write(12,*)' Visibility curves are in files *.v##'
!!$           end if
!!$           if (iC.eq.-3.and.iPsf.ne.0) then
!!$              write(12,*)' Convolved images are in files *.c##'
!!$              if (psftype.lt.3) write(12,*)' Point spread functions are in file *.psf'
!!$           end if
!!$        end if
!!$        !----------------------- en.density ------------------
!!$        if(sph) then
!!$           if(iJ.ne.0) then
!!$              if (iJ.eq.1) then
!!$                 write(12,*)' All energy density profiles are in file *.jtb'
!!$              else
!!$                 write(12,*)' Energy density profiles are in files *.j##'
!!$              end if
!!$           end if
!!$        end if
!!$        !------------------------- radial ------------------
!!$        if(iB.ne.0) then
!!$           if (iB.eq.1) then
!!$              write(12,*)' All radial profiles are in file *.rtb'
!!$           else
!!$              write(12,*)' Radial profiles are in files *.r##'
!!$           end if
!!$        end if
!!$        !------------------------- messages ------------------
!!$        if (iX.eq.1) write(12,*)' All run-time messages are in file *.mtb'
!!$        if (iX.gt.1) write(12,*)' Run-time messages are in files *.m##'
!!$     else
!!$        write(12,*)' Ending calculations for this input file.'
!!$     end if  !end if for error.ne.3
!!$  end if  !end if for models
!!$  if (model.eq.Nmodel.or.error.eq.3.or.error.eq.4) then
!!$     write(12,*) '========== the end =============================='
!!$     close(12)
!!$  end if
!!$  ! close the psf file
!!$  if (model.eq.1.or.error.eq.3) then
!!$     if (iPsf.eq.1.and.psftype.lt.3) close(23)
!!$  end if
!!$  ! conditionally close the spectral files
!!$  if(iA.eq.1) then
!!$     if(model.eq.Nmodel) close(15)
!!$  else
!!$     close(15)
!!$  end if
!!$  ! conditionally close the radial files
!!$  if(iB.eq.1) then
!!$     if(model.eq.Nmodel) close(16)
!!$  else
!!$     close(16)
!!$  end if
!!$  ! conditionally close the imaging files
!!$  if(abs(iC).eq.1) then
!!$     if(model.eq.Nmodel) close(17)
!!$  else
!!$     close(17)
!!$  end if
!!$  if(iX.eq.1) then
!!$     if(model.eq.Nmodel) close(18)
!!$  else
!!$     close(18)
!!$  end if
!!$  ! conditionally close the energy density files
!!$  if(iJ.eq.1) then
!!$     if(model.eq.Nmodel) close(19)
!!$  else
!!$     close(19)
!!$  end if
!!$  if (iPsf.ne.0) close(21)
!!$  if (iV.ne.0) close(22)
!!$!-----------------------------------------------------------------------
!!$  return
!!$end subroutine CLLOSE
!!$!***********************************************************************
!!$
!!$!***********************************************************************
!!$subroutine MakeTable(Elems,rows,cols,unt)
!!$! =======================================================================
!!$!     This is an auxiliary subroutine for print out of tables
!!$!     of Elems(Nrows,Ncols) in output unit 'unt'. Nrows = max{npL,npY}.
!!$!     The array size is defined in Subroutine PrOut.         [MN, Mar'98]
!!$! =======================================================================
!!$      IMPLICIT NONE
!!$      INTEGER rows, cols, unt, k, i
!!$      DOUBLE PRECISION Elems(rows,cols)
!!$! -----------------------------------------------------------------------
!!$      DO k = 1, rows
!!$        write(unt,'(1p,21E11.3)') (Elems(k,i),i=1,cols)
!!$      END DO
!!$! -----------------------------------------------------------------------
!!$      RETURN
!!$END subroutine MakeTable
!!$!***********************************************************************

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
