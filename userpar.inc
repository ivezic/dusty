!c ======================================================================
!c This file defines parameters which control the sizes of Dusty's arrays.
!c The values of these parameters need to be changed only in this file.
!c                                                          [ZI,MN,'99]
!c ======================================================================
!c npL - max size for the wavelength grid (as defined in lambda_grid.dat)
!c npY - max size for the radial (Y) grid
!c npP - max size for the impact parameter (P) and angular (mu) grids
!c npX - max size for the x-grid x(npX) in the optional disk calculation 
!c npR - max size for the output inclination angle tables given in the *.i### files 

!c The npY and npP parameters below are set for the default case of 
!c spherical geometry. For plane-parallel models with large tauV (>~100), 
!c comment the 2 lines for spherical case and uncomment the line for slab
!c ----------------------------------------------------------------------
!c     this parameter must be the same as the number of wavelength 
!c     entries in file lambda_grid.dat
!c      PARAMETER (npL=105)
!c      PARAMETER (npL=118)
      PARAMETER (npL=125)
      PARAMETER (npR=90)
      PARAMETER (npG=10)

!c     spatial grids definitions:
       
!c     ****** spherical case ******
!c     a safe limit which works even on old machines and PCs
!c      PARAMETER (npY=55)
!c     for a typical workstation  
      PARAMETER (npY=970)
!c     do not change these coefficients!
      PARAMETER (npP=20*npY+20)  
!c     do not change except if added disk: 
      PARAMETER (npX=1)   
!c     choose for opt.thick, geom.thin disk calculation:
!c      PARAMETER (npX=2*npY)  

!c     ****** slab case ******
!c      PARAMETER (npY=90, npP=1, npX=1)                      
!c ======================================================================
