!c ======================================================================
!c This file defines parameters which control the sizes of Dusty's arrays.
!c The values of these parameters need to be changed only in this file.
!c                                                          [ZI,MN,'99]
!c ======================================================================
!c npY - max size for the radial (Y) grid
!c npP - max size for the impact parameter (P) and angular (mu) grids
!c npR - max size for the output inclination angle tables given in the *.i### files 

       PARAMETER (npR=90)
       PARAMETER (npG=10)

!c     spatial grids definitions:
!c     a safe limit which works even on old machines and PCs

!c     PARAMETER (npY=128)

!c     for a typical workstation  

       PARAMETER (npY=512)

!c     do not change these coefficients!
       PARAMETER (npP=20*npY+20)  
