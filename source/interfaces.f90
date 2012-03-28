module interfaces
  INTERFACE
     subroutine Rad_Transf(initial,nY,nYprev,nP,itereta,pstar,y_incr,us,fs,emiss, &
     iterfbol,T4_ext)
       use common
       logical, intent(in) :: initial
       integer, intent(in) :: y_incr,iterfbol
       integer :: nY,nP,nYprev,itereta
       double precision pstar
       double precision,allocatable :: us(:,:), fs(:,:) 
       double precision,allocatable :: T4_ext(:)
       double precision,allocatable :: emiss(:,:,:)
     end subroutine Rad_Transf
     SUBROUTINE RADTRANSF_matrix(pstar,iPstar,nY,nYprev,nP,nCav,nIns,TAUlim,FbolOK,initial,deviat,&
          iterFbol,iterEta,model,us,u_old,fs,T4_ext,emiss)
       integer iPstar, FbolOK,iterFbol,iterEta,model,nY,nYprev,nP,nCav,nIns
       double precision :: pstar,TAUlim,deviat
       double precision,allocatable :: us(:,:),u_old(:,:),fs(:,:),T4_ext(:),emiss(:,:,:)
       logical initial
     END SUBROUTINE RADTRANSF_matrix
     subroutine Analysis(nY,model,us,T4_ext,delta)
       integer nY,model
       double precision :: delta
       double precision, allocatable ::  us(:,:),T4_ext(:)
     end subroutine Analysis
     subroutine Find_Tran(pstar,nY,nP,T4_ext,us,fs)
       integer nY, nP
       double precision :: pstar
       double precision, allocatable :: T4_ext(:)
       double precision, allocatable :: fs(:,:),us(:,:)
     end subroutine Find_Tran
     subroutine find_Text(nY,T4_ext)
       integer nY
       double precision, allocatable :: T4_ext(:)
     end subroutine find_Text
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
     subroutine SPH_ext_illum(m0,m1,m1p,m1m,nY,nP)
       integer nY, nP
       double precision,allocatable ::  m0(:,:), m1(:,:),m1p(:,:),m1m(:,:)
     end subroutine SPH_ext_illum
     subroutine SPH_DIFF(nY,nP,flag1,moment_loc,initial,iter,iterfbol,T4_ext,emiss,us,vec2)
       integer nY,nP,flag1,moment_loc,iter,iterfbol
       double precision, allocatable :: T4_ext(:),emiss(:,:,:),us(:,:),vec2(:,:)
       logical initial
     end subroutine SPH_DIFF
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
     subroutine inp_rad(shp,spec_scale,startyp)
       integer :: startyp
       double precision :: spec_scale
       double precision,allocatable :: shp(:)
     end subroutine inp_rad
     subroutine CHKFlux(nY,nYprev,flux,tolern,consfl,iterEta)
       DOUBLE PRECISION,allocatable :: flux(:)
       DOUBLE PRECISION :: tolern
       INTEGER consfl,nY,nYprev,iterEta
     end subroutine CHKFlux
     subroutine invert(nY,mat,Us,Em,Uold,omat)
       integer nY
       double precision,allocatable :: Us(:,:), Uold(:,:), Em(:,:,:),&
            mat(:,:,:),omat(:,:)
     end subroutine invert
     subroutine matrix(pstar,iPstar,m0,m1,mifront,miback,nP,nY,nPok,nYok)
       integer nP,nY,nPok,nYok,iPstar
       double precision :: pstar 
       double precision,allocatable :: m0(:,:,:), m1(:,:,:), mifront(:,:,:), &
            miback(:,:,:)
     end subroutine matrix
     subroutine add(np1,nr1,np2,nr2,q1,q2,q3,qout)
       integer np1, nr1, np2, nr2
       double precision, allocatable :: q1(:,:), q2(:,:), q3(:,:),qout(:,:)
     end subroutine add
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
     SUBROUTINE ChkBolom(nY,qbol,accur,dev,FbolOK)
       integer nY,FbolOK
       double precision :: accur,dev
       DOUBLE PRECISION,allocatable ::  qBol(:) 
     END SUBROUTINE ChkBolom
     subroutine Bolom(q,qbol,nY)
       integer nY
       double precision, allocatable :: q(:,:), qbol(:)
     end subroutine Bolom
     SUBROUTINE LINSYS(Nreal,A,B,X)
       integer Nreal
       DOUBLE PRECISION,allocatable :: A(:,:), B(:), X(:)
     END SUBROUTINE LINSYS
     SUBROUTINE MYSPLINE(x,N,alpha,beta,gamma2,delta)
       integer :: N
       double precision,allocatable :: x(:),alpha(:,:),beta(:,:),gamma2(:,:),delta(:,:)
     end SUBROUTINE MYSPLINE
     SUBROUTINE WEIGHTS(TAUaux,iP,iL,nZ,alpha,beta,gamma2,delta,wgp,wgm,nY)
       integer iP,iL,nZ,nY
       double precision,allocatable :: TAUaux(:,:,:),alpha(:,:), beta(:,:),&
            gamma2(:,:),delta(:,:),wgp(:,:), wgm(:,:)
     end SUBROUTINE WEIGHTS
     SUBROUTINE Kint4(TAUaux,iP,iL,nZ,K1p,K2p,K3p,K4p,K1m,K2m,K3m,K4m)
       INTEGER iP, iL, nZ
       DOUBLE PRECISION,allocatable :: TAUaux(:,:,:), K1p(:),K2p(:),K3p(:),&
            K4p(:), K1m(:), K2m(:), K3m(:), K4m(:)
     end SUBROUTINE Kint4
     SUBROUTINE LUBKSB(A,N,NP,INDX,B)
       integer N,NP
       integer, allocatable :: indx(:)
       double precision,allocatable :: A(:,:),B(:)
     END SUBROUTINE LUBKSB
     SUBROUTINE LUDCMP(A,N,NP,INDX,D)
       integer N,NP
       double precision :: D
       integer, allocatable :: indx(:)
       double precision,allocatable :: A(:,:),B(:)
     END SUBROUTINE LUDCMP
     SUBROUTINE MPROVE(A,ALUD,N,NP,INDX,B,X)
       integer :: n,np
       integer,allocatable :: INDX(:)
       double precision,allocatable :: A(:,:),ALUD(:,:),B(:),X(:)
     END SUBROUTINE MPROVE
     SUBROUTINE NORDLUND(nY,nP,flag,x,f,N1,N2,m,intfdx)
       integer :: nY,nP, flag, N1, N2, m
       double precision :: intfdx
       double precision, allocatable :: x(:),f(:)
     END SUBROUTINE NORDLUND
     subroutine add2(nY,flxs,flxe,fbsum)
       integer :: nY
       double precision, allocatable :: flxs(:,:),flxe(:,:),fbsum(:)
     end subroutine add2
     subroutine Flux_Consv(nY,nYprev,Ncav,itereta, flux1,flux2,fbolOK,maxrat)
       integer nY,nYprev,itereta,fbolOK
       double precision :: maxrat
       double precision, allocatable :: flux1(:),flux2(:)
     end subroutine Flux_Consv
     subroutine shiftIns(x,Nmax,n,xins,i)
       integer :: Nmax,n,i
       double precision :: xins
       double precision,allocatable :: x(:)
     end subroutine shiftIns
     subroutine getOptPr(nameQ,nameNK,er,stdf,top,szds,qsd,a1,a2,nFiles,xC,XCuser)
       integer er,top,szds,nFiles
       character*235,allocatable,nameQ(:)
       character*235 nameNK(10),stdf(7)
       double precision :: qsd,a1,a2,xC(10),xCuser(10)
     end subroutine getOptPr
     SUBROUTINE MULTIPLY(type,np1,nr1,np2,nr2,mat,vec1,omat,flag,q1,q2)
       integer :: type,np1,nr1,np2,nr2,flag
       double precision, allocatable :: mat(:,:,:), vec1(:,:), omat(:,:), &
            q1(:,:), q2(:,:)
     END SUBROUTINE MULTIPLY
     subroutine SLBdiff(nY,flag,grid,T4_ext,em,fp,fm)
       !---parameter
       integer :: nY,flag
       double precision,allocatable :: grid(:,:),T4_ext(:),em(:,:,:),&
            fp(:,:),fm(:,:)
     end subroutine SLBdiff
     SUBROUTINE SPH_Int(nY,nP,fs)
       integer nY,nP
       double precision,allocatable :: fs(:,:)
     END SUBROUTINE SPH_Int
     SUBROUTINE GetbOut(nP,pstar,k)
       integer nP
       double precision :: pstar
     END SUBROUTINE GetbOut
  END INTERFACE
end module interfaces
