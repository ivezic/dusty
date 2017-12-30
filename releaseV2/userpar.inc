c ======================================================================
c This file defines parameters which control the sizes of Dusty's arrays.
c The values of these parameters need to be changed only in this file.
c This file MUST stay with DUSTY's source code in the same directory.
c ======================================================================
c npL - size for the wavelength grid (as defined in lambda.inc)
c npY - size for the radial (Y) grid
c npP - size for the impact parameter (P) and angular (mu) grids

c The npY and npP parameters below are set for the default case of 
c spherical geometry. For plane-parallel models with large tauV (>~100), 
c comment the 2 lines for spherical case and uncomment the line for slab
c ----------------------------------------------------------------------
c     this parameter must be the same as the number of wavelength 
c     entries in file lambda.inc
      PARAMETER (npL=105)
c     spatial grids definitions:
       
c     ****** spherical case ******
c     a safe limit which works even on old machines and PCs
      PARAMETER (npY=40)
c     for a typical workstation  
c      PARAMETER (npY=65)
c     do not change these coefficients!
      PARAMETER (npP=4*npY+95)  

c     ****** slab case ******
c      PARAMETER (npY=120, npP=1)                      
c ======================================================================
