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
