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
  integer i, n, n1, n2
  double precision x(n), y(n), wgth, integral, dyn2
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



!!$!***********************************************************************
!!$subroutine ROMBY(fnc,a,b,ss)
!!$!=======================================================================
!!$! This subroutine performs Romberg integration of function func on
!!$! interval [a,b]. The result is returned in ss. Desired accuracy is set
!!$! to 0.002.                                            [Z.I., Feb. 1996]
!!$! =======================================================================
!!$  IMPLICIT NONE
!!$  INTEGER JMAX,JMAXP,K,KM, J
!!$  PARAMETER (JMAX=30, JMAXP=JMAX+1, K=3, KM=K-1)
!!$  DOUBLE PRECISION a,b,fnc,ss,EPS_loc, aux, dss,h(JMAXP),s(JMAXP)
!!$  EXTERNAL fnc
!!$  ! ---------------------------------------------------------------------
!!$  EPS_loc = 0.002d0
!!$  h(1)=1.0d0
!!$  do j=1,JMAX
!!$     call trapzd(fnc,a,b,s(j),j)
!!$     if (j.ge.K) then
!!$        aux = 0.0d0
!!$        call polint(h(j-KM),s(j-KM),K,aux,ss,dss)
!!$        IF (dabs(dss).le.EPS_loc*dabs(ss)) RETURN
!!$     endif
!!$     s(j+1)=s(j)
!!$     h(j+1)=0.25d0*h(j)
!!$  end do
!!$  ! --------------------------------------------------------------------
!!$  RETURN
!!$END subroutine ROMBY
!!$!***********************************************************************
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
!!$  implicit none
!!$  integer Nmax, n, i,j
!!$  double precision x(Nmax),xins
!!$  ! ---------------------------------------------------------------------
!!$  do j = n+1, i+2, -1
!!$     x(j) = x(j-1)
!!$  end do
!!$  x(i+1) = xins
!!$  ! -----------------------------------------------------------------------
!!$  return
!!$end subroutine shiftIns
!!$!***********************************************************************
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
!!$
!!$!***********************************************************************
!!$SUBROUTINE trapzd(func,a,b,s,n)
!!$! =======================================================================
!!$  IMPLICIT NONE
!!$  INTEGER n
!!$  DOUBLE PRECISION a,b,s,func
!!$  EXTERNAL func
!!$  INTEGER it,j
!!$  DOUBLE PRECISION del,sum,tnm,x
!!$  ! ----------------------------------------------------------------------
!!$  IF (n.eq.1) THEN
!!$     s=0.5d0*(b-a)*(func(a)+func(b))
!!$  ELSE
!!$     it=2**(n-2)
!!$     tnm=it
!!$     del=(b-a)/tnm
!!$     x=a+0.5d0*del
!!$     sum=0.
!!$     DO j = 1, it
!!$        sum=sum+func(x)
!!$        x=x+del
!!$     END DO
!!$     s=0.5d0*(s+(b-a)*sum/tnm)
!!$  END IF
!!$  ! -------------------------------------------------------------------------
!!$  RETURN
!!$END SUBROUTINE trapzd
!!$!***********************************************************************
!!$
!!$!***********************************************************************
!!$subroutine polint(xa,ya,n,x,y,dy)
!!$  ! For polinomial interpolation, used in Subroutine Romby.
!!$  ! ====================================================================
!!$  implicit none
!!$  integer n,Nmax
!!$  double precision dy,x,y,xa(n),ya(n)
!!$  parameter (Nmax=1000)
!!$  integer i,m,ns
!!$  double precision den,dif,dift,ho,hp,w,c(Nmax),d(Nmax)
!!$  !---------------------------------------------------------------------
!!$  c = 0.0d0
!!$  d = 0.0d0
!!$  ns=1
!!$  dif=dabs(x-xa(1))
!!$  do i=1,n
!!$     dift=dabs(x-xa(i))
!!$     if (dift.lt.dif) then
!!$        ns=i
!!$        dif=dift
!!$     endif
!!$     c(i)=ya(i)
!!$     d(i)=ya(i)
!!$  end do
!!$  y=ya(ns)
!!$  ns=ns-1
!!$  do m=1,n-1
!!$     do i=1,n-m
!!$        ho=xa(i)-x
!!$        hp=xa(i+m)-x
!!$        w=c(i+1)-d(i)
!!$        den=ho-hp
!!$        if(den.eq.0.0d0) then
!!$           write(6,'(a)') 'failure in polint'
!!$           stop
!!$        endif
!!$        den=w/den
!!$        d(i)=hp*den
!!$        c(i)=ho*den
!!$     end do
!!$     if (2*ns.lt.n-m)then
!!$        dy=c(ns+1)
!!$     else
!!$        dy=d(ns)
!!$        ns=ns-1
!!$     endif
!!$     y=y+dy
!!$  end do
!!$!-----------------------------------------------------------------------
!!$  return
!!$end subroutine polint
!!$!***********************************************************************
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
!!$!***********************************************************************
!!$subroutine gauleg(x1,x2,xg,wg,n)
!!$!=====================================================================
!!$  implicit none
!!$
!!$  integer i,m,n,j
!!$  double precision x1,x2,xm,xl,eps,delj,p,eta
!!$  double precision xg(n),wg(n),sum,ff,p1,p2,p3,z1,z,pp
!!$  parameter (eps=1.0d-14)
!!$  ! -------------------------------------------------------------------
!!$  xg = 0.0d0
!!$  wg = 0.0d0
!!$  m = int((n+1)/2)
!!$  xm = 0.5d0*(x2+x1)
!!$  xl = 0.5d0*(x2-x1)
!!$  do i = 1, m
!!$     z = cos(3.1415926535898d0*(dble(i) - 0.25d0)/(dble(n) + 0.5d0))
!!$1    continue
!!$     p1 = 1.0d0
!!$     p2 = 0.0d0
!!$     do j = 1, n
!!$        p3 = p2
!!$        p2 = p1
!!$        p1 = ((2.0d0*dble(j)-1.0d0)*z*p2-(j-1.0d0)*p3)/dble(j)
!!$     end do
!!$     pp = dble(n)*(z*p1 - p2)/(z*z - 1.0d0)
!!$     z1 = z
!!$     z = z1-p1/pp
!!$     if (abs(z-z1).gt.eps) go to 1
!!$     xg(i) = xm-xl*z
!!$     xg(n+1-i) = xm+xl*z
!!$     wg(i) = 2.d0*xl/((1.d0-z*z)*pp*pp)
!!$     wg(n+1-i) = wg(i)
!!$  end do
!!$  !-------------------------------------------------------------------
!!$  return
!!$end subroutine gauleg
!!$!***********************************************************************
!!$
!!$!**********************************************************************
!!$SUBROUTINE ANALINT(Nanal,xaux,yaux,m,aux,error)
!!$!======================================================================
!!$! This subroutine calculates integral I(x**m*y(x)*dx). Both y and x are
!!$! 1D arrays, y(i), x(i) with i=1,Nanal. The method used is approximation
!!$! of y(x) by y = P(x) + d/sqrt(1-x*x), where P(x) is the polynomial of
!!$! order Nanal-1, and analytic evaluation of the integral. It is assumed
!!$! that xaux(1)=0. Coefficients are determined from the set of Nanal
!!$! linear equations and subsequent call to the linear system solver
!!$! LINSYS.                                              [Z.I., Nov. 1995]
!!$! ANALINT is called from Nordlund to evaluate analytically the contribution
!!$! of Nanal grid points. [MN]
!!$! =======================================================================
!!$  use common
!!$  IMPLICIT none
!!$
!!$  INTEGER i, j, Nanal, error
!!$  DOUBLE PRECISION xaux(Nanal),yaux(Nanal),coeff(Nanal),A(npY,npY),m,aux,b
!!$  ! ---------------------------------------------------------------------
!!$  error = 0
!!$  ! generate matrix A and vector B
!!$  DO i = 1, Nanal
!!$     DO j = 1, Nanal-1
!!$        IF (xaux(i).EQ.0.0.AND.j.EQ.1) THEN
!!$           A(i,j) = 1.0
!!$        ELSE
!!$           A(i,j) = xaux(i)**(1.0*j-1.0)
!!$        END IF
!!$     END DO
!!$     A(i,Nanal) = 1.0/sqrt(1.0-xaux(i)*xaux(i))
!!$  END DO
!!$  ! solve for the coefficients
!!$  CALL LINSYS(Nanal,A,yaux,coeff,error)
!!$  IF(error.NE.0) THEN
!!$     CALL MSG(19)
!!$     iERROR = iERROR + 1
!!$     RETURN
!!$  END IF
!!$  ! upper limit for integration:
!!$  b = xaux(Nanal)
!!$  ! evaluate m-dependent contribution of the last term
!!$  IF (m.GT.0.1) THEN
!!$     IF (m.GT.1.1) THEN
!!$        ! this is for m=2
!!$        aux = 0.5*(DASIN(b)-b*sqrt(1.-b*b))
!!$     ELSE
!!$        ! this is for m=1
!!$        aux = 1.0 - sqrt(1.-b*b)
!!$     ENDIF
!!$  ELSE
!!$     ! this is for m=0
!!$     aux = DASIN(b)
!!$  ENDIF
!!$  aux = aux * coeff(Nanal)
!!$  ! add contribution from the polynom
!!$  DO i = 1, Nanal-1
!!$     aux = aux + coeff(i) * (b**(m+1.0*i)) / (m+1.0*i)
!!$  END DO
!!$! -----------------------------------------------------------------------
!!$999   RETURN
!!$END SUBROUTINE ANALINT
!!$!***********************************************************************
!!$
!!$! ***********************************************************************
!!$SUBROUTINE ChkConv(accuracy_loc,Aold,Anew,Aconv_loc)
!!$! =======================================================================
!!$! This subroutine checks convergence of an array A(nY) between values
!!$! given in Aold and Anew. If the relative difference for EVERY element
!!$! is smaller than accuracy, Aconv is assigned 1, otherwise 0.
!!$!                                                      [Z.I., Jul. 1996]
!!$! =======================================================================
!!$  use common
!!$  IMPLICIT none
!!$  INTEGER iY, Aconv_loc
!!$  DOUBLE PRECISION accuracy_loc, Aold(npY), Anew(npY), delta
!!$! -----------------------------------------------------------------------
!!$  Aconv_loc = 1
!!$  ! loop over radial positions
!!$  DO iY = 1, nY
!!$     ! find relative difference
!!$     delta = dabs(Anew(iY)-Aold(iY))
!!$     IF (delta.GT.dabs(Anew(iY))*accuracy_loc) Aconv_loc = 0
!!$  END DO
!!$  ! -----------------------------------------------------------------------
!!$  RETURN
!!$END SUBROUTINE ChkConv
!!$! ***********************************************************************
!!$
!!$!***********************************************************************
!!$SUBROUTINE LINSYS(Nreal,A,B,X,error)
!!$!=======================================================================
!!$! This subroutine solves the set of linear equations [A]*[X] = [B] for
!!$! X [A(k,1)*X(1)+A(k,2)*X(2)+...+A(k,Nreal)*X(Nreal) = B(k), k=1,Nreal).
!!$! The real size of matrix A is Nreal x Nreal and its physical dimension
!!$! is npY x npY, where npY comes from INCLUDE 'userpar.inc'. Both vectors
!!$! B and X have real lengths Nreal. The set is solved by calls to LUDCMP
!!$! and LUBKSB and the solution is improved subsequently by a call to
!!$! MPROVE. These three subroutines are taken from Numerical Recipes.
!!$!                                                      [Z.I., Nov. 1995]
!!$! =======================================================================
!!$  use common
!!$  IMPLICIT none
!!$  
!!$  INTEGER Nreal, indx(npY), i, j, error
!!$  DOUBLE PRECISION A(npY,npY), B(npY), X(npY), &
!!$              A1c(npY,npY), B1(npY), A2c(npY,npY), B2(npY), d
!!$  ! ---------------------------------------------------------------------
!!$  error = 0
!!$  ! generate DOUBLE PRECISION copies of A and B (two copies because they
!!$  ! are changed in LUDCMP and LUBKSB, but still needed for MPROVE)
!!$  DO i = 1, Nreal
!!$     B1(i) = B(i)
!!$     B2(i) = B(i)
!!$     DO j = 1, Nreal
!!$        A1c(i,j) = A(i,j)
!!$        A2c(i,j) = A(i,j)
!!$     END DO
!!$  END DO
!!$  ! solve the system
!!$  CALL LUDCMP(A1c,Nreal,npY,indx,d,error)
!!$  IF (error.NE.0) RETURN
!!$  CALL LUBKSB(A1c,Nreal,npY,indx,B1)
!!$  ! improve the solution (saved in B)
!!$  CALL MPROVE(A2c,A1c,Nreal,npY,indx,B2,B1)
!!$  ! copy the improved solution to output vector X
!!$  DO i = 1, Nreal
!!$     X(i) = B1(i)
!!$  END DO
!!$  ! --------------------------------------------------------------------
!!$  RETURN
!!$END SUBROUTINE LINSYS
!!$!***********************************************************************
!!$
!!$! ***********************************************************************
!!$SUBROUTINE LUBKSB(A,N,NP,INDX,B)
!!$  ! =====================================================================
!!$  DIMENSION INDX(NP)
!!$  DOUBLE PRECISION A(NP,NP),B(NP)
!!$  ! -------------------------------------------------------------------
!!$  II=0
!!$  !impossible to parallelize since B(J) needs to is changed and used!
!!$  DO I=1,N
!!$     LL=INDX(I)
!!$     SUM=B(LL)
!!$     B(LL)=B(I)
!!$     IF (II.NE.0)THEN
!!$        DO J=II,I-1
!!$           SUM=SUM-A(I,J)*B(J)
!!$        END DO
!!$     ELSE IF (SUM.NE.0.) THEN
!!$        II=I
!!$     ENDIF
!!$     B(I)=SUM
!!$  END DO
!!$  !impossible to parallelize since B(J) needs to is changed and used!
!!$  DO I=N,1,-1
!!$     SUM=B(I)
!!$     IF(I.LT.N)THEN
!!$        DO J=I+1,N
!!$           SUM=SUM-A(I,J)*B(J)
!!$        END DO
!!$     ENDIF
!!$     B(I)=SUM/A(I,I)
!!$  END DO
!!$  ! -------------------------------------------------------------------
!!$  RETURN
!!$END SUBROUTINE LUBKSB
!!$! ***********************************************************************
!!$
!!$! ***********************************************************************
!!$SUBROUTINE LUDCMP(A,N,NP,INDX,D,error)
!!$  ! =====================================================================
!!$  PARAMETER (NMAX=10000,TINY=1.0E-20)
!!$  DIMENSION INDX(NP)
!!$  INTEGER error
!!$  DOUBLE PRECISION A(NP,NP),VV(NMAX), D, SUM
!!$  ! ------------------------------------------------------------------
!!$  error = 0
!!$  D = 1.
!!$  DO I = 1, N
!!$     AAMAX=0.
!!$     DO J = 1, N
!!$        IF (DABS(A(I,J)).GT.AAMAX) AAMAX=DABS(A(I,J))
!!$     END DO
!!$     ! IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
!!$     IF (AAMAX.EQ.0.) THEN
!!$        error = 5
!!$        RETURN
!!$     ENDIF
!!$     VV(I)=1./AAMAX
!!$  END DO
!!$!!  !$OMP PARALLEL DO private(J,I,SUM,DUM,IMAX,AAMAX,D)
!!$  DO J = 1 , N
!!$     IF (J.GT.1) THEN
!!$        DO I = 1, J-1
!!$           SUM=A(I,J)
!!$           IF (I.GT.1)THEN
!!$              DO K = 1, I-1
!!$                 SUM=SUM-A(I,K)*A(K,J)
!!$              END DO
!!$              A(I,J)=SUM
!!$           ENDIF
!!$        END DO
!!$     ENDIF
!!$     AAMAX=0.
!!$     DO I = J, N
!!$        SUM=A(I,J)
!!$        IF (J.GT.1)THEN
!!$           DO K = 1, J-1
!!$              SUM=SUM-A(I,K)*A(K,J)
!!$           END DO
!!$           A(I,J)=SUM
!!$        ENDIF
!!$        DUM=VV(I)*DABS(SUM)
!!$        IF (DUM.GE.AAMAX) THEN
!!$           IMAX=I
!!$           AAMAX=DUM
!!$        ENDIF
!!$     END DO
!!$     IF (J.NE.IMAX)THEN
!!$        DO K = 1, N
!!$           DUM=A(IMAX,K)
!!$           A(IMAX,K)=A(J,K)
!!$           A(J,K)=DUM
!!$        END DO
!!$        D=-D
!!$        VV(IMAX)=VV(J)
!!$     ENDIF
!!$     INDX(J)=IMAX
!!$     IF(J.NE.N)THEN
!!$        IF(A(J,J).EQ.0.)A(J,J)=TINY
!!$        DUM=1./A(J,J)
!!$        DO I = J+1, N
!!$           A(I,J)=A(I,J)*DUM
!!$        END DO
!!$     ENDIF
!!$  END DO
!!$!!  !$OMP END PARALLEL DO
!!$  IF(A(N,N).EQ.0.)A(N,N)=TINY
!!$  !------------------------------------------------------------------
!!$  RETURN
!!$END SUBROUTINE LUDCMP
!!$! ***********************************************************************
!!$
!!$!***********************************************************************
!!$SUBROUTINE MPROVE(A,ALUD,N,NP,INDX,B,X)
!!$!=======================================================================
!!$  PARAMETER (NMAX=10000)
!!$  DIMENSION INDX(N)
!!$  DOUBLE PRECISION SDP,A(NP,NP),ALUD(NP,NP),B(N),X(N),R(NMAX)
!!$  ! ---------------------------------------------------------------------
!!$  DO i = 1, N
!!$     SDP = -B(i)
!!$     DO j = 1, N
!!$        SDP = SDP + A(i,j)*X(j)
!!$     END DO
!!$     R(i) = SDP
!!$  END DO
!!$  CALL LUBKSB(ALUD,N,NP,INDX,R)
!!$  DO i = 1, N
!!$     X(i) = X(i) - R(i)
!!$  END DO
!!$  ! -----------------------------------------------------------------------
!!$  RETURN
!!$END SUBROUTINE MPROVE
!!$!***********************************************************************
!!$
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
!!$
!!$!**********************************************************************
!!$subroutine add(np1,nr1,np2,nr2,q1,q2,q3,qout)
!!$!======================================================================
!!$! This subroutine evaluates the following expression:
!!$! [qOut] = [q1] + [q2] + [q3]. qout, q1, q2 and q2 are matrices of
!!$! physical size (np2,np1) and real size (nr2,nr1).     [Z.I., Nov. 1995]
!!$! ======================================================================
!!$  implicit none
!!$  integer npY, npP, npX, npL, npG, npR
!!$  include '../userpar.inc'
!!$
!!$  integer  np1, nr1, np2, nr2, i2, i1
!!$  double precision  q1(np2,np1), q2(np2,np1), q3(np2,np1),qout(np2,np1)
!!$! ----------------------------------------------------------------------
!!$
!!$! loop over index 2
!!$  do i2 = 1, nr2
!!$! loop over index 1
!!$   do i1 = 1, nr1
!!$    qout(i2,i1) = q1(i2,i1) +  q2(i2,i1) + q3(i2,i1)
!!$   end do
!!$  end do
!!$! ----------------------------------------------------------------------
!!$  return
!!$end subroutine add
!!$!**********************************************************************
!!$
!!$!***********************************************************************
!!$SUBROUTINE ChkRange(dr,x)
!!$!=======================================================================
!!$! This subroutine checks if x is within the allowed range defined by
!!$! dr<<1:
!!$!         dr**2 < x < 1/dr**2
!!$! If it is not then x = 0.0                            [Z.I., Jan. 1997]
!!$!=======================================================================
!!$  IMPLICIT none
!!$  DOUBLE PRECISION x, dr
!!$!-----------------------------------------------------------------------
!!$  IF ((x-dr*dr)*(x-1.0d0/dr/dr).LT.0.0d0) THEN
!!$     continue
!!$  ELSE
!!$     ! continue
!!$     x = 0.0d0
!!$  END IF
!!$!-----------------------------------------------------------------------
!!$  RETURN
!!$END SUBROUTINE ChkRange
!!$!***********************************************************************
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
