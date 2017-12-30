                                                                                
c**********************************************************************         
c     This is the block data with optical constants for the supported           
c     grain types. It has to be at the beginning of the program.                
c     The data was compiled from different sources by Z. Ivezic (1996).         
c                                                          [MN, Apr'98]         
c =====================================================================         
      BLOCK DATA                                                                
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER npLnk                                                             
      PARAMETER (npLnk=98)                                                      
      DOUBLE PRECISION n_sil_ow(npLnk),k_sil_ow(npLnk),n_sil_oc(npLnk),         
     &       k_sil_oc(npLnk), n_sil_dl(npLnk), k_sil_dl(npLnk),                 
     &       n_amc_hn(npLnk), k_amc_hn(npLnk), n_sic_pg(npLnk),                 
     &       k_sic_pg(npLnk), n_gr1_dl(npLnk), k_gr1_dl(npLnk),                 
     &       n_gr2_dl(npLnk), k_gr2_dl(npLnk), lam_nk(npLnk)                    
      COMMON /nkdat/ n_sil_ow, k_sil_ow, n_sil_oc, k_sil_oc,                    
     &               n_sil_dl, k_sil_dl, n_amc_hn, k_amc_hn,                    
     &               n_sic_pg, k_sic_pg, n_gr1_dl, k_gr1_dl,                    
     &               n_gr2_dl, k_gr2_dl, lam_nk                                 
      DATA lam_nk/1.00E-02,2.00E-02,3.00E-02,4.00E-02,5.00E-02,                 
     &   6.00E-02,8.00E-02,1.00E-01,1.20E-01,1.50E-01,2.00E-01,2.50E-01,        
     &   3.00E-01,3.60E-01,4.40E-01,5.50E-01,7.00E-01,8.50E-01,1.00E+00,        
     &   1.15E+00,1.30E+00,1.70E+00,2.00E+00,2.20E+00,2.70E+00,3.00E+00,        
     &   3.50E+00,4.00E+00,4.50E+00,5.00E+00,5.50E+00,6.00E+00,6.50E+00,        
     &   7.00E+00,7.50E+00,8.00E+00,8.50E+00,9.00E+00,9.40E+00,9.55E+00,        
     &   9.70E+00,9.85E+00,1.00E+01,1.05E+01,1.10E+01,1.13E+01,1.16E+01,        
     &   1.20E+01,1.25E+01,1.30E+01,1.35E+01,1.40E+01,1.45E+01,1.50E+01,        
     &   1.60E+01,1.70E+01,1.80E+01,1.90E+01,2.00E+01,2.20E+01,2.30E+01,        
     &   2.40E+01,2.50E+01,2.60E+01,2.70E+01,2.80E+01,3.00E+01,3.50E+01,        
     &   4.00E+01,4.50E+01,5.00E+01,5.50E+01,6.00E+01,6.50E+01,7.00E+01,        
     &   7.50E+01,8.00E+01,8.50E+01,9.00E+01,9.50E+01,1.00E+02,1.05E+02,        
     &   1.10E+02,1.20E+02,1.30E+02,1.40E+02,1.50E+02,2.00E+02,2.50E+02,        
     &   3.00E+02,4.00E+02,5.00E+02,7.00E+02,1.30E+03,4.00E+03,1.30E+04,        
     &   2.00E+04,3.60E+04/                                                     
      DATA n_sil_ow/1.81E+00,1.81E+00,1.81E+00,1.81E+00,1.81E+00,               
     &   1.81E+00,1.81E+00,1.81E+00,1.81E+00,1.81E+00,1.81E+00,1.81E+00,        
     &   1.81E+00,1.81E+00,1.82E+00,1.84E+00,1.85E+00,1.85E+00,1.85E+00,        
     &   1.85E+00,1.85E+00,1.87E+00,1.88E+00,1.88E+00,1.89E+00,1.89E+00,        
     &   1.88E+00,1.87E+00,1.87E+00,1.85E+00,1.83E+00,1.81E+00,1.79E+00,        
     &   1.75E+00,1.71E+00,1.64E+00,1.55E+00,1.45E+00,1.40E+00,1.43E+00,        
     &   1.44E+00,1.46E+00,1.51E+00,1.73E+00,1.91E+00,2.00E+00,2.11E+00,        
     &   2.21E+00,2.21E+00,2.16E+00,2.09E+00,2.06E+00,2.03E+00,1.98E+00,        
     &   1.91E+00,1.90E+00,1.96E+00,2.04E+00,2.12E+00,2.27E+00,2.34E+00,        
     &   2.36E+00,2.38E+00,2.40E+00,2.42E+00,2.43E+00,2.45E+00,2.60E+00,        
     &   2.66E+00,2.70E+00,2.74E+00,2.76E+00,2.77E+00,2.79E+00,2.80E+00,        
     &   2.81E+00,2.82E+00,2.83E+00,2.85E+00,2.86E+00,2.87E+00,2.87E+00,        
     &   2.87E+00,2.87E+00,2.88E+00,2.88E+00,2.88E+00,2.89E+00,2.90E+00,        
     &   2.90E+00,2.90E+00,2.90E+00,2.91E+00,2.91E+00,2.91E+00,2.91E+00,        
     &   2.91E+00,2.91E+00/                                                     
      DATA k_sil_ow/9.99E-02,9.99E-02,9.99E-02,9.99E-02,9.99E-02,               
     &   9.99E-02,9.99E-02,9.99E-02,9.99E-02,9.99E-02,9.99E-02,9.99E-02,        
     &   9.99E-02,9.99E-02,9.08E-02,7.17E-02,6.03E-02,5.64E-02,5.46E-02,        
     &   6.12E-02,6.75E-02,7.23E-02,7.33E-02,6.96E-02,6.19E-02,5.88E-02,        
     &   5.48E-02,5.24E-02,5.07E-02,5.14E-02,5.24E-02,5.50E-02,5.79E-02,        
     &   6.36E-02,7.07E-02,9.03E-02,1.31E-01,2.63E-01,4.31E-01,4.96E-01,        
     &   5.60E-01,6.30E-01,7.12E-01,8.10E-01,7.94E-01,7.86E-01,7.68E-01,        
     &   6.13E-01,4.93E-01,3.98E-01,3.74E-01,3.69E-01,3.66E-01,3.72E-01,        
     &   4.39E-01,5.63E-01,6.68E-01,7.23E-01,7.51E-01,7.55E-01,7.34E-01,        
     &   7.10E-01,6.85E-01,6.61E-01,6.37E-01,6.44E-01,6.58E-01,6.28E-01,        
     &   5.78E-01,5.23E-01,4.69E-01,4.34E-01,4.13E-01,3.92E-01,3.70E-01,        
     &   3.49E-01,3.28E-01,3.06E-01,2.85E-01,2.64E-01,2.42E-01,2.37E-01,        
     &   2.32E-01,2.22E-01,2.12E-01,2.02E-01,1.92E-01,1.42E-01,1.00E-01,        
     &   9.09E-02,7.24E-02,5.39E-02,3.79E-02,2.15E-02,7.26E-03,2.39E-03,        
     &   2.39E-03,2.39E-03/                                                     
      DATA n_sil_oc/1.77E+00,1.77E+00,1.77E+00,1.77E+00,1.77E+00,               
     &   1.77E+00,1.77E+00,1.77E+00,1.77E+00,1.77E+00,1.77E+00,1.77E+00,        
     &   1.77E+00,1.77E+00,1.79E+00,1.82E+00,1.83E+00,1.82E+00,1.81E+00,        
     &   1.81E+00,1.81E+00,1.84E+00,1.85E+00,1.86E+00,1.87E+00,1.88E+00,        
     &   1.88E+00,1.87E+00,1.86E+00,1.85E+00,1.83E+00,1.80E+00,1.78E+00,        
     &   1.74E+00,1.69E+00,1.62E+00,1.51E+00,1.39E+00,1.34E+00,1.36E+00,        
     &   1.37E+00,1.39E+00,1.44E+00,1.67E+00,1.86E+00,1.96E+00,2.08E+00,        
     &   2.22E+00,2.23E+00,2.17E+00,2.09E+00,2.04E+00,2.01E+00,1.92E+00,        
     &   1.78E+00,1.73E+00,1.80E+00,1.98E+00,2.17E+00,2.40E+00,2.48E+00,        
     &   2.53E+00,2.58E+00,2.63E+00,2.69E+00,2.70E+00,2.74E+00,2.90E+00,        
     &   2.95E+00,2.98E+00,3.01E+00,3.02E+00,3.03E+00,3.03E+00,3.03E+00,        
     &   3.04E+00,3.04E+00,3.05E+00,3.05E+00,3.05E+00,3.06E+00,3.06E+00,        
     &   3.06E+00,3.06E+00,3.06E+00,3.06E+00,3.06E+00,3.07E+00,3.07E+00,        
     &   3.07E+00,3.07E+00,3.08E+00,3.08E+00,3.08E+00,3.08E+00,3.08E+00,        
     &   3.08E+00,3.08E+00/                                                     
      DATA k_sil_oc/8.95E-02,8.95E-02,8.95E-02,8.95E-02,8.95E-02,               
     &   8.95E-02,8.95E-02,8.95E-02,8.95E-02,8.95E-02,8.95E-02,8.95E-02,        
     &   8.95E-02,8.95E-02,8.23E-02,6.64E-02,5.09E-02,4.46E-02,4.69E-02,        
     &   6.37E-02,7.89E-02,9.62E-02,1.03E-01,1.00E-01,9.27E-02,8.97E-02,        
     &   8.50E-02,8.14E-02,7.85E-02,7.80E-02,7.77E-02,8.06E-02,8.40E-02,        
     &   8.91E-02,9.51E-02,1.10E-01,1.48E-01,2.95E-01,4.78E-01,5.49E-01,        
     &   6.19E-01,6.96E-01,7.87E-01,9.11E-01,9.10E-01,9.14E-01,9.07E-01,        
     &   7.40E-01,5.94E-01,4.74E-01,4.39E-01,4.28E-01,4.20E-01,4.35E-01,        
     &   5.35E-01,7.36E-01,9.52E-01,1.10E+00,1.13E+00,1.07E+00,1.04E+00,        
     &   1.01E+00,0.97E+00,0.94E+00,8.98E-01,8.77E-01,8.42E-01,7.52E-01,        
     &   6.60E-01,5.68E-01,4.75E-01,4.25E-01,4.04E-01,3.82E-01,3.60E-01,        
     &   3.39E-01,3.17E-01,2.96E-01,2.74E-01,2.52E-01,2.31E-01,2.26E-01,        
     &   2.21E-01,2.12E-01,2.02E-01,1.92E-01,1.83E-01,1.35E-01,9.44E-02,        
     &   8.56E-02,6.82E-02,5.07E-02,3.57E-02,2.02E-02,6.82E-03,2.29E-03,        
     &   2.29E-03,2.29E-03/                                                     
      DATA n_sil_dl/8.66E-01,8.66E-01,8.66E-01,8.25E-01,7.81E-01,               
     &   8.90E-01,1.29E+00,1.60E+00,1.87E+00,2.26E+00,1.93E+00,1.80E+00,        
     &   1.76E+00,1.74E+00,1.73E+00,1.72E+00,1.72E+00,1.71E+00,1.71E+00,        
     &   1.71E+00,1.71E+00,1.71E+00,1.71E+00,1.71E+00,1.70E+00,1.70E+00,        
     &   1.69E+00,1.68E+00,1.66E+00,1.64E+00,1.60E+00,1.57E+00,1.52E+00,        
     &   1.47E+00,1.34E+00,1.21E+00,1.16E+00,1.11E+00,1.20E+00,1.24E+00,        
     &   1.30E+00,1.35E+00,1.39E+00,1.57E+00,1.75E+00,1.82E+00,1.90E+00,        
     &   2.00E+00,2.04E+00,2.09E+00,2.04E+00,2.00E+00,1.91E+00,1.82E+00,        
     &   1.73E+00,1.69E+00,1.71E+00,1.78E+00,1.89E+00,2.04E+00,2.11E+00,        
     &   2.18E+00,2.26E+00,2.30E+00,2.35E+00,2.40E+00,2.50E+00,2.63E+00,        
     &   2.76E+00,2.89E+00,3.02E+00,3.07E+00,3.13E+00,3.18E+00,3.23E+00,        
     &   3.25E+00,3.27E+00,3.29E+00,3.31E+00,3.33E+00,3.35E+00,3.35E+00,        
     &   3.36E+00,3.37E+00,3.38E+00,3.39E+00,3.40E+00,3.41E+00,3.42E+00,        
     &   3.42E+00,3.43E+00,3.43E+00,3.43E+00,3.43E+00,3.43E+00,3.43E+00,        
     &   3.43E+00,3.43E+00/                                                     
      DATA k_sil_dl/1.39E-01,1.39E-01,1.39E-01,2.55E-01,4.29E-01,               
     &   6.73E-01,8.78E-01,9.26E-01,7.22E-01,5.30E-01,5.32E-02,2.77E-02,        
     &   2.84E-02,2.88E-02,2.91E-02,2.94E-02,2.97E-02,3.00E-02,3.03E-02,        
     &   3.06E-02,3.09E-02,3.21E-02,3.31E-02,3.39E-02,3.60E-02,3.72E-02,        
     &   3.94E-02,4.11E-02,4.25E-02,4.40E-02,4.72E-02,5.05E-02,5.36E-02,        
     &   5.66E-02,1.14E-01,1.71E-01,3.68E-01,5.66E-01,7.65E-01,8.29E-01,        
     &   8.71E-01,8.97E-01,9.24E-01,9.63E-01,1.00E+00,9.65E-01,9.28E-01,        
     &   8.78E-01,7.71E-01,6.63E-01,5.70E-01,4.77E-01,4.83E-01,4.88E-01,        
     &   5.83E-01,7.14E-01,8.55E-01,9.76E-01,1.05E+00,1.08E+00,1.10E+00,        
     &   1.11E+00,1.13E+00,1.12E+00,1.12E+00,1.12E+00,1.11E+00,1.06E+00,        
     &   1.01E+00,9.66E-01,9.19E-01,8.66E-01,8.13E-01,7.59E-01,7.06E-01,        
     &   6.72E-01,6.38E-01,6.03E-01,5.69E-01,5.35E-01,5.00E-01,4.83E-01,        
     &   4.67E-01,4.33E-01,3.99E-01,3.65E-01,3.32E-01,2.48E-01,2.06E-01,        
     &   1.65E-01,1.32E-01,9.87E-02,7.04E-02,3.83E-02,2.46E-02,2.46E-02,        
     &   2.46E-02,2.46E-02/                                                     
      DATA n_amc_hn/8.40E-01,8.40E-01,8.40E-01,8.40E-01,7.40E-01,               
     &   6.90E-01,9.30E-01,1.53E+00,1.74E+00,1.55E+00,1.22E+00,1.40E+00,        
     &   1.60E+00,1.71E+00,1.78E+00,1.85E+00,1.94E+00,2.01E+00,2.11E+00,        
     &   2.21E+00,2.31E+00,2.50E+00,2.63E+00,2.70E+00,2.82E+00,2.86E+00,        
     &   2.95E+00,3.03E+00,3.04E+00,3.04E+00,3.09E+00,3.15E+00,3.16E+00,        
     &   3.16E+00,3.25E+00,3.35E+00,3.38E+00,3.42E+00,3.47E+00,3.49E+00,        
     &   3.51E+00,3.53E+00,3.55E+00,3.59E+00,3.63E+00,3.65E+00,3.67E+00,        
     &   3.70E+00,3.73E+00,3.75E+00,3.79E+00,3.84E+00,3.86E+00,3.88E+00,        
     &   3.96E+00,3.99E+00,4.10E+00,4.14E+00,4.18E+00,4.26E+00,4.31E+00,        
     &   4.34E+00,4.36E+00,4.42E+00,4.47E+00,4.50E+00,4.57E+00,4.77E+00,        
     &   4.94E+00,5.11E+00,5.25E+00,5.32E+00,5.44E+00,5.49E+00,5.67E+00,        
     &   5.71E+00,5.85E+00,5.90E+00,5.99E+00,5.94E+00,6.32E+00,6.70E+00,        
     &   6.79E+00,6.99E+00,7.17E+00,7.37E+00,7.55E+00,8.50E+00,9.32E+00,        
     &   1.01E+01,1.15E+01,1.25E+01,1.41E+01,1.65E+01,1.65E+01,1.65E+01,        
     &   1.65E+01,1.65E+01/                                                     
      DATA k_amc_hn/1.08E-01,1.08E-01,1.08E-01,1.08E-01,1.77E-01,               
     &   3.80E-01,9.00E-01,8.40E-01,5.60E-01,1.77E-01,3.21E-01,7.40E-01,        
     &   7.20E-01,6.86E-01,6.70E-01,6.95E-01,7.70E-01,8.25E-01,9.00E-01,        
     &   9.38E-01,9.60E-01,9.95E-01,1.02E+00,1.01E+00,9.96E-01,9.90E-01,        
     &   1.01E+00,1.04E+00,1.04E+00,1.03E+00,1.09E+00,1.15E+00,1.20E+00,        
     &   1.25E+00,1.34E+00,1.42E+00,1.44E+00,1.47E+00,1.50E+00,1.51E+00,        
     &   1.52E+00,1.53E+00,1.54E+00,1.57E+00,1.60E+00,1.61E+00,1.62E+00,        
     &   1.63E+00,1.66E+00,1.69E+00,1.71E+00,1.74E+00,1.76E+00,1.78E+00,        
     &   1.83E+00,1.89E+00,1.95E+00,1.95E+00,1.98E+00,2.05E+00,2.08E+00,        
     &   2.13E+00,2.19E+00,2.23E+00,2.27E+00,2.32E+00,2.40E+00,2.60E+00,        
     &   2.77E+00,2.92E+00,3.00E+00,3.18E+00,3.31E+00,3.50E+00,3.70E+00,        
     &   3.80E+00,4.00E+00,4.15E+00,4.30E+00,4.48E+00,4.59E+00,4.70E+00,        
     &   4.79E+00,4.97E+00,5.15E+00,5.33E+00,5.51E+00,6.41E+00,7.17E+00,        
     &   7.92E+00,9.43E+00,1.04E+01,1.25E+01,1.41E+01,1.41E+01,1.41E+01,        
     &   1.41E+01,1.41E+01/                                                     
      DATA n_sic_pg/7.06E-01,7.06E-01,7.06E-01,7.06E-01,7.06E-01,               
     &   7.06E-01,7.06E-01,7.06E-01,9.69E-01,2.07E+00,4.29E+00,5.16E+00,        
     &   4.75E+00,3.42E+00,2.59E+00,2.55E+00,2.51E+00,2.50E+00,2.47E+00,        
     &   2.49E+00,2.50E+00,2.51E+00,2.55E+00,2.58E+00,2.65E+00,2.69E+00,        
     &   2.73E+00,2.76E+00,2.76E+00,2.77E+00,2.74E+00,2.71E+00,2.67E+00,        
     &   2.63E+00,2.57E+00,2.51E+00,2.42E+00,2.33E+00,2.17E+00,2.11E+00,        
     &   2.04E+00,1.96E+00,1.88E+00,1.56E+00,1.25E+00,1.47E+00,1.69E+00,        
     &   1.99E+00,3.07E+00,4.14E+00,4.21E+00,4.28E+00,4.11E+00,3.93E+00,        
     &   3.71E+00,3.59E+00,3.52E+00,3.49E+00,3.46E+00,3.43E+00,3.42E+00,        
     &   3.40E+00,3.39E+00,3.39E+00,3.39E+00,3.38E+00,3.38E+00,3.38E+00,        
     &   3.39E+00,3.39E+00,3.39E+00,3.40E+00,3.41E+00,3.41E+00,3.42E+00,        
     &   3.42E+00,3.43E+00,3.43E+00,3.44E+00,3.44E+00,3.45E+00,3.45E+00,        
     &   3.46E+00,3.46E+00,3.47E+00,3.48E+00,3.49E+00,3.51E+00,3.52E+00,        
     &   3.52E+00,3.52E+00,3.52E+00,3.52E+00,3.52E+00,3.52E+00,3.52E+00,        
     &   3.52E+00,3.52E+00/                                                     
      DATA k_sic_pg/1.53E+00,1.53E+00,1.53E+00,1.53E+00,1.53E+00,               
     &   1.53E+00,1.53E+00,1.53E+00,1.43E+00,1.28E+00,1.06E+00,3.50E-01,        
     &   2.50E-01,1.62E-01,1.04E-01,1.13E-01,1.21E-01,1.30E-01,1.44E-01,        
     &   1.62E-01,1.80E-01,2.17E-01,2.63E-01,2.76E-01,2.77E-01,2.78E-01,        
     &   2.51E-01,2.23E-01,1.88E-01,1.53E-01,1.22E-01,9.08E-02,8.07E-02,        
     &   7.06E-02,6.69E-02,6.31E-02,6.59E-02,6.87E-02,9.61E-02,1.05E-01,        
     &   1.09E-01,1.24E-01,1.38E-01,5.43E-01,9.48E-01,1.40E+00,1.85E+00,        
     &   2.45E+00,2.45E+00,2.45E+00,1.69E+00,9.19E-01,7.00E-01,4.80E-01,        
     &   3.62E-01,2.90E-01,2.80E-01,2.82E-01,2.70E-01,2.57E-01,2.51E-01,        
     &   2.45E-01,2.38E-01,2.35E-01,2.32E-01,2.29E-01,2.23E-01,2.13E-01,        
     &   2.04E-01,1.94E-01,1.84E-01,1.80E-01,1.75E-01,1.71E-01,1.66E-01,        
     &   1.63E-01,1.60E-01,1.57E-01,1.54E-01,1.51E-01,1.48E-01,1.46E-01,        
     &   1.43E-01,1.39E-01,1.34E-01,1.29E-01,1.25E-01,1.01E-01,9.22E-02,        
     &   8.30E-02,8.30E-02,8.30E-02,8.30E-02,8.30E-02,8.30E-02,8.30E-02,        
     &   8.30E-02,8.30E-02/                                                     
      DATA n_gr1_dl/9.81E-01,9.25E-01,8.48E-01,7.72E-01,7.23E-01,               
     &   7.65E-01,9.91E-01,1.39E+00,2.57E+00,1.94E+00,1.61E+00,1.56E+00,        
     &   1.93E+00,2.17E+00,2.35E+00,2.34E+00,2.29E+00,2.25E+00,2.23E+00,        
     &   2.21E+00,2.20E+00,2.19E+00,2.18E+00,2.17E+00,2.17E+00,2.16E+00,        
     &   2.16E+00,2.15E+00,2.15E+00,2.14E+00,2.13E+00,2.13E+00,2.12E+00,        
     &   2.11E+00,2.10E+00,2.09E+00,2.08E+00,2.07E+00,2.06E+00,2.06E+00,        
     &   2.05E+00,2.05E+00,2.04E+00,2.03E+00,2.01E+00,1.98E+00,2.05E+00,        
     &   2.01E+00,1.99E+00,1.97E+00,1.96E+00,1.95E+00,1.93E+00,1.92E+00,        
     &   1.89E+00,1.87E+00,1.84E+00,1.82E+00,1.80E+00,1.76E+00,1.74E+00,        
     &   1.73E+00,1.72E+00,1.72E+00,1.71E+00,1.71E+00,1.71E+00,1.76E+00,        
     &   1.83E+00,1.92E+00,2.02E+00,2.11E+00,2.21E+00,2.30E+00,2.40E+00,        
     &   2.49E+00,2.57E+00,2.66E+00,2.74E+00,2.82E+00,2.90E+00,2.97E+00,        
     &   3.05E+00,3.19E+00,3.33E+00,3.46E+00,3.58E+00,4.15E+00,4.65E+00,        
     &   5.10E+00,5.90E+00,6.59E+00,7.80E+00,9.33E+00,9.33E+00,9.33E+00,        
     &   9.33E+00,9.33E+00/                                                     
      DATA k_gr1_dl/3.08E-03,2.78E-02,9.28E-02,2.09E-01,3.70E-01,               
     &   5.39E-01,9.54E-01,1.65E+00,6.55E-01,1.94E-01,2.75E-01,6.38E-01,        
     &   8.16E-01,6.57E-01,4.79E-01,2.55E-01,1.41E-01,9.41E-02,6.96E-02,        
     &   5.51E-02,4.60E-02,3.33E-02,2.87E-02,2.68E-02,2.42E-02,2.35E-02,        
     &   2.32E-02,2.36E-02,2.44E-02,2.57E-02,2.72E-02,2.91E-02,3.13E-02,        
     &   3.41E-02,3.74E-02,4.14E-02,4.61E-02,5.16E-02,5.65E-02,5.85E-02,        
     &   6.05E-02,6.26E-02,6.47E-02,7.23E-02,8.11E-02,8.83E-02,9.49E-02,        
     &   9.92E-02,1.10E-01,1.22E-01,1.34E-01,1.47E-01,1.61E-01,1.75E-01,        
     &   2.07E-01,2.42E-01,2.79E-01,3.20E-01,3.63E-01,4.55E-01,5.04E-01,        
     &   5.55E-01,6.06E-01,6.59E-01,7.11E-01,7.64E-01,8.68E-01,1.11E+00,        
     &   1.33E+00,1.52E+00,1.69E+00,1.85E+00,1.99E+00,2.11E+00,2.23E+00,        
     &   2.34E+00,2.44E+00,2.54E+00,2.63E+00,2.72E+00,2.81E+00,2.89E+00,        
     &   2.97E+00,3.12E+00,3.27E+00,3.40E+00,3.53E+00,4.12E+00,4.62E+00,        
     &   5.08E+00,5.88E+00,6.58E+00,7.79E+00,9.32E+00,9.32E+00,9.32E+00,        
     &   9.32E+00,9.32E+00/                                                     
      DATA n_gr2_dl/9.80E-01,9.20E-01,8.39E-01,7.64E-01,5.55E-01,               
     &   5.43E-01,1.01E+00,2.08E+00,1.87E+00,1.27E+00,6.88E-01,1.25E+00,        
     &   2.39E+00,2.58E+00,2.66E+00,2.74E+00,2.87E+00,3.03E+00,3.19E+00,        
     &   3.34E+00,3.47E+00,3.76E+00,3.94E+00,4.05E+00,4.29E+00,4.43E+00,        
     &   4.64E+00,4.82E+00,4.99E+00,5.14E+00,5.29E+00,5.43E+00,5.57E+00,        
     &   5.69E+00,5.80E+00,5.90E+00,6.01E+00,6.11E+00,6.20E+00,6.22E+00,        
     &   6.25E+00,6.28E+00,6.31E+00,6.39E+00,6.47E+00,6.52E+00,6.56E+00,        
     &   6.62E+00,6.67E+00,6.74E+00,6.81E+00,6.89E+00,6.96E+00,7.03E+00,        
     &   7.15E+00,7.26E+00,7.36E+00,7.46E+00,7.53E+00,7.83E+00,7.99E+00,        
     &   8.13E+00,8.25E+00,8.47E+00,8.73E+00,9.01E+00,9.60E+00,1.14E+01,        
     &   1.35E+01,1.59E+01,1.88E+01,2.13E+01,2.33E+01,2.48E+01,2.59E+01,        
     &   2.66E+01,2.70E+01,2.70E+01,2.66E+01,2.62E+01,2.56E+01,2.46E+01,        
     &   2.34E+01,2.05E+01,1.75E+01,1.49E+01,1.33E+01,1.33E+01,1.68E+01,        
     &   2.11E+01,3.02E+01,3.89E+01,5.47E+01,7.40E+01,7.40E+01,7.40E+01,        
     &   7.40E+01,7.40E+01/                                                     
      DATA k_gr2_dl/3.09E-03,2.79E-02,9.38E-02,1.36E-01,3.62E-01,               
     &   7.20E-01,1.48E+00,1.24E+00,3.38E-01,2.77E-01,1.09E+00,2.27E+00,        
     &   2.08E+00,1.67E+00,1.56E+00,1.56E+00,1.68E+00,1.82E+00,1.94E+00,        
     &   2.03E+00,2.12E+00,2.33E+00,2.48E+00,2.58E+00,2.82E+00,2.96E+00,        
     &   3.18E+00,3.39E+00,3.60E+00,3.79E+00,4.00E+00,4.19E+00,4.37E+00,        
     &   4.55E+00,4.72E+00,4.91E+00,5.09E+00,5.27E+00,5.41E+00,5.46E+00,        
     &   5.50E+00,5.55E+00,5.60E+00,5.77E+00,5.94E+00,6.04E+00,6.14E+00,        
     &   6.27E+00,6.44E+00,6.63E+00,6.82E+00,7.00E+00,7.17E+00,7.35E+00,        
     &   7.69E+00,8.05E+00,8.42E+00,8.80E+00,9.24E+00,1.01E+01,1.05E+01,        
     &   1.10E+01,1.15E+01,1.20E+01,1.24E+01,1.29E+01,1.38E+01,1.59E+01,        
     &   1.76E+01,1.90E+01,1.95E+01,1.92E+01,1.83E+01,1.73E+01,1.62E+01,        
     &   1.50E+01,1.38E+01,1.26E+01,1.16E+01,1.08E+01,9.98E+00,9.24E+00,        
     &   8.73E+00,8.70E+00,9.94E+00,1.24E+01,1.54E+01,2.94E+01,3.95E+01,        
     &   4.78E+01,6.08E+01,7.09E+01,8.61E+01,1.03E+02,1.03E+02,1.03E+02,        
     &   1.03E+02,1.03E+02/                                                     
                                                                                
      END                                                                       
c =====================================================================         
                                                                                
c **********************************************************************        
      program Dusty                                                             
c =====================================================================         
c This program solves the continuum radiative transfer problem for a            
c spherically symmetric envelope or for a plane-parallel slab. Only this        
c file needs to be compiled. It is assumed that all input data are given        
c in files named *.inp, and that the list of these files is given in a          
c master input file dusty.inp. For details see the Manual.                      
c                                                       [Z.I. and M.N.]         
c ======================================================================        
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER nYok, nPok                                                        
      DOUBLE PRECISION                                                          
     &       Ude(npL,npY), Uds(npL,npY), Us(npL,npY), fde(npL,npY),             
     &       fds(npL,npY), fs(npL,npY), Utot(npL,npY), ftot(npL,npY),           
     &       Td(npG,npY), Ubol(npY), fbol(npY), Uchck(npL,npY),                 
     &       UbolChck(npY), Spectrum(npL), SpecChar(30,99), tauF(npY),          
     &       Intens(npL,npP+2), IntOut(20,npP+2), SmC(30,99),                   
     &       Yok(npY), Pok(npP), tauOut(npL), tauZout(npP+2), F1,               
     &       tr(npY), rg(npG,npY), fsL(npL,npY), fsR(npL,npY), Eps(npY)         
      COMMON /solution/ Ude, Uds, Us, fde, fds, fs, Utot, ftot, Td,             
     &       Ubol, fbol, Uchck, UbolChck, Spectrum, SpecChar, tauF,             
     &       Intens, IntOut, SmC, Yok, Pok, tauOut, tauZout, F1,                
     &       tr, rg, fsL, fsR, Eps, nYok, nPok                                  
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      CHARACTER*3 version                                                       
      CHARACTER*230 path, apath                                                 
      CHARACTER*235 nameIn, nameOut, nameQ(npG), nameNK(10)                     
      INTEGER error, nG, model, Nmodel, GridType, io1, Empty, lpath,            
     &        iL, Nrec, Nlam                                                    
c     Nrec is the max number of records for TAUgrid in a file                   
      PARAMETER (Nrec = 1000)                                                   
      DOUBLE PRECISION ETAzp(npP,npY), TAUin(Nrec), TAU1, TAU2, RDINP           
      LOGICAL Equal                                                             
      Equal = .True.                                                            
c ----------------------------------------------------------------------        
c     **************************                                                
c     *** ABOUT THIS VERSION ***                                                
c     **************************                                                
      version= '2.01'                                                            
c     Updated versions of Dusty(2.0) with minor changes and bug fixes
c     start from 2.01.

c     Version (2.0) is the first public release. The code has been              
c     significantly improved in terms of speed and I/O options. All             
c     suggestions of the users of version(1.0) had been taken into              
c     consideration. Finished Oct.'99                                           
                                                                                
c     Version 1.0 is a beta version sent to a few people before the             
c     first public release.  Finished Nov,'96.                                  
c **********************************************************************        
c     *** MAIN ***                                                              
c     first open the file with lambda grid and check that the grid satisfies    
c     certain conditions (the wavelengths are in microns) :                     
      open(4, file='lambda_grid.dat', status = 'OLD')                           
      Nlam = RDINP(Equal,4)                                                     
      IF (Nlam.NE.npL) THEN                                                     
       write(*,*)' *************** A BIG ERROR !!! ***************** '          
       write(*,*)'  The number of wavelengths in lambda_grid.dat is  '          
       write(*,*)'  not equal to the specified npL in userpar.inc    '          
       write(*,*)'  Make sure the numbers are the same, recompile    '          
       write(*,*)'  and try again.                                   '          
       write(*,*)' ************************************************* '          
       goto 999                                                                 
      END IF                                                                    
c     Initialize lambda array                                                   
      read(4,*,END=99) (lambda(iL), iL = 1, npL)                                
 99   close(4)                                                                 
      CALL sort(lambda,npL)                                                     
c     Check the ends of the lambda grid :                                       
      IF(lambda(1).GT.0.01) THEN                                                
       write(*,*)' *************** WARNING! ********************** '            
       write(*,*)'  The shortest wavelength in lambda_grid.dat has '            
       write(*,*)'  to be 0.01 microns. Correct this and try again!'            
       write(*,*)' *********************************************** '            
       goto 999                                                                 
      END IF                                                                    
      IF(lambda(npL).LT.36000.) THEN                                            
       write(*,*)' *************** WARNING! ******************* '               
       write(*,*)'  The longest wavelength in lambda_grid.dat   '               
       write(*,*)'  has to be 36 mm. Correct this and try again!'               
       write(*,*)' ******************************************** '               
       goto 999                                                                 
      END IF                                                                    
c     Check the resolution:                                                     
      DO iL = 2, npL                                                            
        IF (lambda(iL)/lambda(iL-1).GT.1.51) THEN                               
          write(*,*)' ***************** WARNING!  *******************'          
          write(*,*)' The ratio of two consecutive wavelengths in the'          
          write(*,*)' grid has to be no bigger than 1.5. You have    '          
          write(*,'(2(a4,1p,e8.2))') '    ',lambda(iL)/lambda(iL-1),            
     &                                ' at ', lambda(iL)                        
          write(*,*)' Correct this and try again!                    '          
          write(*,*)' ***********************************************'          
          goto 999                                                              
        END IF                                                                  
      END DO                                                                    
c     open master input file dusty.inp                                          
      open(13,ERR=998,file='dusty.inp',STATUS='OLD')                            
      io1 = 0                                                                   
c     loop over input files                                                     
      DO WHILE (io1.GE.0)                                                       
c       read a line from master input file using                                
100     read(13,'(a)',iostat=io1) apath                                         
        if(io1.lt.0) then                                                       
         stop                                                                   
        end if                                                                  
        CALL Clean(apath, path, lpath)                                          
c       if not EOF and if line is not empty, or commented, proceed              
        IF (Empty(path).NE.1) THEN                                              
c         get input/output file names                                           
          CALL ATTACH(path,lpath,'.inp',nameIn)                                 
          CALL ATTACH(path,lpath,'.out',nameOut)                                
c         read input data                                                       
          CALL Input(nameIn,nG,nameOut,nameQ,nameNK,                            
     &               TAU1,TAU2,TAUin,Nrec,GridType,Nmodel,error,version)        
           IF (iVerb.GT.0)                                                      
     &          write(*,'(a24,a80)') ' Working on Input File: ',nameIn          
           IF (iVerb.EQ.2) write(*,*) 'Done with Reading Input'                 
c         if an error reading files go to the next input file                   
c         error=3 means some files are missing                                  
          IF (error.EQ.3) goto 100                                              
c         get optical properties                                                
          CALL getOptPr(nG,nameQ,nameNK,error)                                  
c         if an error reading files go to the next input file                   
          IF (error.EQ.3) goto 100                                              
          IF (iVerb.EQ.2) write(*,*) 'Done with getOptPr'                       
c         solve for every model                                                 
          DO model = 1, Nmodel                                                  
           IF (iVerb.GT.0)  write(*,'(a9,i4)') ' model = ',model                
           IF (error.EQ.0) THEN                                                 
c            open output files                                                  
             CALL Oppen(model,path,lpath)                                       
c            calculate optical depth for current model                          
             CALL GetTau(model,nG,TAU1,TAU2,TAUin,Nrec,GridType,Nmodel)         
             IF (iVerb.EQ.2)                                                    
     &          write(*,*) 'Done with GetTau. Going to Solve'                   
c            solve                                                              
             CALL Solve(model,nG,error,ETAzp)                                   
             IF (error.EQ.0) THEN                                               
c              calculate spectral characteristics                               
               CALL Spectral(model,denstyp,nL,Lambda)                           
c              write results out                                                
               CALL PrOut(model)                                                
             ELSE                                                               
              goto 100                                                          
             END IF                                                             
c            close output files                                                 
             CALL Cllose(error,model,Nmodel)                                    
           END IF                                                               
           IF (iVerb.EQ.2) write(*,*) ' ----------- '                           
           IF (iVerb.EQ.2) write(*,*) '    '                                    
c         end of the loop over models                                           
          END DO                                                                
        END IF                                                                  
c     end of the loop over input files                                          
      END DO                                                                    
      IF (iVerb.GT.0) write(*,*) ' End of Input Files '                         
      close(13)                                                                 
c     end this run                                                              
      goto 999                                                                  
c     to execute if the master input file is missing                            
998   write(*,*)' *********** FATAL ERROR IN DUSTY ***********'                 
      write(*,*)' * Master input file dusty.inp is missing!? *'                 
      write(*,*)' ********************************************'                 
c ----------------------------------------------------------------------        
999   STOP                                                                      
      END                                                                       
c **********************************************************************        
                                                                                
c =======================================================================       
c This is the include file with I/O subroutines              [MN, Apr'99]       
c =======================================================================       
C     Table of Contents                                                         
C                                                                               
C     CLLOSE                                                                    
C     INPSTAR                                                                   
C     INPUT                                                                     
C     OPPEN                                                                     
C     PROUT                                                                     
C     RDINP                                                                     
C     VAL                                                                       
c =======================================================================       
                                                                                
c ***********************************************************************       
      SUBROUTINE CLLOSE(error,model,Nmodel)                                     
c =======================================================================       
c This subroutine closes output files.             [ZI,Feb'96; MN,Apr'99]       
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER startyp(2), Nlamtr(2), nBB(2)                                     
      CHARACTER nameStar(2)*235                                                 
      DOUBLE PRECISION Tstar, lamtr(2,101), klam(2,100), Tbb(2,10),             
     &                 rellum(2,10), mu1, ksi, mu2, xSiO, r1rs                  
      COMMON /source/ Tstar, lamtr, klam, Tbb, rellum, mu1, ksi, mu2,           
     &                xSiO, r1rs, startyp, Nlamtr, nBB, nameStar                
      INTEGER psftype, Npsf, iLambda                                            
      DOUBLE PRECISION kPSF(20), FWHM1(20), FWHM2(20), theta1,                  
     &       xpsf(1000), ypsf(1000)                                             
      COMMON /psf1/ iLambda, psftype, Npsf                                      
      COMMON /psf2/ kPSF, FWHM1, FWHM2, Theta1,                                 
     &       xpsf, ypsf                                                         
      CHARACTER*72 su1, su2, s3, s4                                             
      INTEGER  error, model, Nmodel                                             
c -----------------------------------------------------------------------       
c     Close the default output file:                                            
      IF (iERROR.NE.0.AND.iOUT.EQ.1) THEN                                       
        write(12,'(a42,i4)')                                                    
     &                 ' There are some error messages for model:',model        
        write(12,*)                                                             
     &         ' Please check m## file (if not produced then rerun)'            
      END IF                                                                    
      IF (iWARNING.NE.0.AND.iOUT.EQ.1.AND.iERROR.EQ.0) THEN                     
        write(12,'(a36,i4)')                                                    
     &  ' There are some warnings for model:',model                             
        write(12,*)                                                             
     &         ' Please check m## file (if not produced then rerun)'            
      END IF                                                                    
      iCUMM = iCUMM + iERROR + iWARNING                                         
      IF (model.EQ.Nmodel.OR.error.EQ.3.OR.error.EQ.4) THEN                     
      IF (error.NE.3) THEN                                                      
       su1=' =========================================================='        
       su2='============================'                                       
        IF(denstyp.eq.4.OR.RDW) THEN                                            
          write(12,'(a59,a28)')su1,su2                                          
        ELSE                                                                    
         IF(denstyp.NE.0) THEN                                                  
           write(12,'(a59)')su1                                                 
         ELSE                                                                   
          su1=' ====================================================='          
          IF(ksi.GT.0) THEN                                                     
           su2='=============='                                                 
           write(12,'(a54,a14)')su1,su2                                         
          ELSE                                                                  
           su2='====='                                                          
           write(12,'(a54,a5)')su1,su2                                          
          END IF                                                                
         END IF                                                                 
        END IF                                                                  
        write(12,'(a23,1p,e8.1,a8)')'   (1) Optical depth at',lamfid,           
     &                              ' microns'                                  
        IF(denstyp.eq.0) THEN                                                   
c  ----------  for slab output ----------------------------                     
         write(12,*)                                                            
     &'  (2) Bol.flux of the left-side source at the slab left boundary'        
         write(12,*)                                                            
     & '  (3) f1=F/Fe1, where F is the overall bol.flux in the slab'            
         write(12,*)                                                            
     & '  (4) Position of the left slab boundary for L=1E4 Lsun'                
         write(12,*)                                                            
     & '  (5) Dust temperature at the right slab face'                          
         write(12,*)                                                            
     & '  (6) Effective temperature of the left source (in K)'                  
         IF(ksi.GT.0) THEN                                                      
           write(12,*)                                                          
     & '  (7) Effective temperature of the right source (in K)'                 
           write(12,*)'  (8) Maximum error in flux conservation (%)'            
         ELSE                                                                   
           write(12,*)'  (7) Maximum error in flux conservation (%)'            
         END IF                                                                 
        ELSE                                                                    
c    ---------- for spherical shell ----------------------------                
         write(12,*)'  (2) Bolometric flux at the inner radius '                
         write(12,*)'  (3) Inner radius for L=1E4 Lsun'                         
         write(12,*)'  (4) Ratio of the inner to the stellar radius'            
         write(12,*)'  (5) Angular size (in arcsec) when Fbol=1E-6 W/m2'        
         write(12,*)'  (6) Dust temperature at the outer edge (in K)'           
         write(12,*)'  (7) Maximum error in flux conservation (%)'              
         IF(denstyp.eq.4.OR.RDW) THEN                                           
          write(12,*)'  (8) Mass-loss rate (in Msun/yr)'                        
          write(12,*)'  (9) Terminal outflow velocity (in km/s)'                
          write(12,*)'  (10) Upper limit of the stellar mass (in Msun)'         
         END IF                                                                 
        END IF                                                                  
        write(12,*)'================================================='          
        IF(iCUMM.EQ.0) write(12,*)' Everything is OK for all models'            
        IF(iSPP.NE.0) THEN                                                      
         IF(denstyp.eq.0.AND.iSPP.eq.3) THEN                                    
          write(12,*)                                                           
     & ' Tables with spectral properties are in files *.spp and *.zpp'          
         ELSE                                                                   
          write(12,*)' Table with spectral properties is in file *.spp'         
         END IF                                                                 
        END IF                                                                  
       IF(iA.NE.0) THEN                                                         
        IF (iA.EQ.1) THEN                                                       
          write(12,*)' All spectra are in file *.stb'                           
        ELSE                                                                    
         IF (denstyp.EQ.0.AND.iA.EQ.3) THEN                                     
           write(12,*)' Spectra are in files *.s## and *.z##'                   
         ELSE                                                                   
           write(12,*)' Spectra are in files *.s##'                             
         END IF                                                                 
        END IF                                                                  
       END IF                                                                   
       IF (denstyp.NE.0.AND.iC.NE.0) THEN                                       
        IF (abs(iC).EQ.1) THEN                                                  
          write(12,*)' All imaging quantities are in file *.itb'                
        ELSE                                                                    
          write(12,*)' Images are in files *.i##'                               
          IF (iV.NE.0.AND.abs(iC).EQ.3)                                         
     &      write(12,*)' Visibility curves are in files *.v##'                  
        END IF                                                                  
        IF (iC.EQ.-3.AND.iPSF.NE.0) THEN                                        
          write(12,*)' Convolved images are in files *.c##'                     
          IF (psftype.LT.3)                                                     
     &    write(12,*)' Point spread functions are in file *.psf'                
        END IF                                                                  
       END IF                                                                   
       IF(iB.NE.0) THEN                                                         
        IF (iB.EQ.1) THEN                                                       
          write(12,*)' All radial profiles are in file *.rtb'                   
        ELSE                                                                    
          write(12,*)' Radial profiles are in files *.r##'                      
        END IF                                                                  
       END IF                                                                   
       IF (iX.EQ.1)                                                             
     &     write(12,*)' All run-time messages are in file *.mtb'                
       IF (iX.GT.1)                                                             
     &     write(12,*)' Run-time messages are in files *.m##'                   
      ELSE                                                                      
        write(12,*)' Ending calculations for this input file'                   
      END IF                                                                    
      END IF                                                                    
                                                                                
      IF (model.EQ.Nmodel.OR.error.EQ.3.OR.error.EQ.4) THEN                     
        write(12,*)                                                             
     &    '========== THE END =============================='                   
      END IF                                                                    
                                                                                
c     Table with spectral properties                                            
      IF (iSPP.NE.0) THEN                                                       
c       In case of slab: add the zpp table after the spp (if desired)           
        IF(denstyp.eq.0.AND.model.eq.Nmodel.AND.iSPP.NE.3) THEN                 
          s3='###   tau0      Psi      fV       fK       f12    C21  '          
          s4=' C31   C43  b8-13 b14-22 B9.8 B11.4  R9.8-18  '                   
          write(19,'(a49)')                                                     
     &              '# ==============================================='         
          write(19,'(a49)')                                                     
     &              '# Properties of Spectra from the slab left side  '         
          write(19,'(a49)')                                                     
     &              '# -----------------------------------------------'         
          write(19,'(a55,a46)')s3,s4                                            
          DO model = 1, Nmodel                                                  
            write(19,'(a100)') zline(model)                                     
          END DO                                                                
         Close(19)                                                              
         Close(24)                                                              
        END IF                                                                  
      END IF                                                                    
      IF (model.EQ.1.OR.error.EQ.3) THEN                                        
        IF (iPSF.EQ.1.AND.psftype.LT.3) Close(23)                               
      END IF                                                                    
c     conditionally close the spectral files                                    
      IF(iA.EQ.1) THEN                                                          
       IF(model.EQ.Nmodel) close(15)                                            
      ELSE                                                                      
       close(15)                                                                
      END IF                                                                    
c     conditionally close the radial files                                      
      IF(iB.EQ.1) THEN                                                          
       IF(model.EQ.Nmodel) close(16)                                            
      ELSE                                                                      
       close(16)                                                                
      END IF                                                                    
c     conditionally close the imaging files                                     
      IF(abs(iC).EQ.1) THEN                                                     
       IF(model.EQ.Nmodel) close(17)                                            
      ELSE                                                                      
       close(17)                                                                
      END IF                                                                    
      IF(iX.EQ.1) THEN                                                          
       IF(model.EQ.Nmodel) close(18)                                            
      ELSE                                                                      
       close(18)                                                                
      END IF                                                                    
                                                                                
      IF (iPSF.NE.0) Close(21)                                                  
      IF (iV.NE.0) Close(22)                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE InpStar(error,is,nameIn)                                       
c =======================================================================       
c     This subroutine is for reading the stellar input parameters               
c                                                             [MN,Mar'99]       
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER startyp(2), Nlamtr(2), nBB(2)                                     
      CHARACTER nameStar(2)*235                                                 
      DOUBLE PRECISION Tstar, lamtr(2,101), klam(2,100), Tbb(2,10),             
     &                 rellum(2,10), mu1, ksi, mu2, xSiO, r1rs                  
      COMMON /source/ Tstar, lamtr, klam, Tbb, rellum, mu1, ksi, mu2,           
     &                xSiO, r1rs, startyp, Nlamtr, nBB, nameStar                
      CHARACTER strg*40, nameIn*(*)                                             
      INTEGER error, i, is                                                      
      DOUBLE PRECISION RDINP, sum                                               
      LOGICAL Equal, NoEqual                                                    
      Equal = .True.                                                            
      NoEqual = .False.                                                         
c -------------------------                                                     
      error = 0                                                                 
c     generic temperature for cases other than BB or Engelke-Marengo shape.     
      Tstar = 10000.0                                                           
c     flag for the external spectrum                                            
      startyp(is) = RDINP(Equal,1)                                              
c     if stellar flag is not btw.1 and 6 stop:                                  
      IF(startyp(is).LT.1 .OR. startyp(is).GT.6) THEN                           
       CALL MSG(11)                                                             
       error = 3                                                                
       goto 999                                                                 
      END IF                                                                    
c     #1: black body(ies) for startyp=1                                         
      IF (startyp(is).EQ.1) THEN                                                
c       number of black bodies                                                  
        nBB(is) = RDINP(Equal,1)                                                
c       stellar temperature(s)                                                  
        Tbb(is,1) = RDINP(Equal,1)                                              
        IF (nBB(is).GT.1) THEN                                                  
c       multiple black bodies                                                   
          Tstar = 0.0                                                           
          DO i = 2, nBB(is)                                                     
            Tbb(is,i) = RDINP(NoEqual,1)                                        
            IF (Tbb(is,i).LE.0.0) THEN                                          
              CALL MSG(8)                                                       
              error = 1                                                         
              GOTO 999                                                          
            END IF                                                              
          END DO                                                                
c         read in relative luminosities                                         
          rellum(is,1) = RDINP(Equal,1)                                         
          sum = rellum(is,1)                                                    
          DO i = 2, nBB(is)                                                     
            rellum(is,i) = RDINP(NoEqual,1)                                     
            sum = sum + rellum(is,i)                                            
          END DO                                                                
          IF (sum.LE.0.0) THEN                                                  
             CALL MSG(7)                                                        
             error = 1                                                          
             GOTO 999                                                           
          END IF                                                                
c         normalize to 1                                                        
          DO i = 1, nBB(is)                                                     
            rellum(is,i) = rellum(is,i) / sum                                   
            Tstar = Tstar + rellum(is,i)*Tbb(is,i)**4                           
          END DO                                                                
          Tstar= Tstar**0.25                                                    
        ELSE                                                                    
c         this is for a single black body - Tbb(is,1)                           
          IF(is.eq.1) THEN                                                      
           Tstar = Tbb(1,1)                                                     
           rellum(1,1) = 1.0                                                    
          ELSE                                                                  
           Tstar = Tbb(2,1)                                                     
           rellum(2,1) = 1.0                                                    
          END IF                                                                
          IF (Tbb(is,1).LE.0.0) THEN                                            
             CALL MSG(8)                                                        
             error = 1                                                          
             GOTO 999                                                           
          ELSE                                                                  
c        end if for first or second source                                      
         END IF                                                                 
c       end if for one or many BB                                               
        END IF                                                                  
c     end if for BB-type                                                        
      END IF                                                                    
                                                                                
c     #2: Engelke-Marengo function for startyp=2                                
      IF (startyp(is).EQ.2) THEN                                                
c        effective stellar temperature                                          
         Tbb(is,1) = RDINP(Equal,1)                                             
         Tstar = Tbb(is,1)                                                      
c        depth of SiO abs.feature in %                                          
         xSiO = RDINP(Equal,1)                                                  
         IF (xSiO.LE.0.0) xSiO = 0.0001                                         
         IF (xSiO.GT.100.0) xSiO = 100.0                                        
      END IF                                                                    
                                                                                
c     #3: power-law(s) for startyp=3                                            
      IF (startyp(is).EQ.3) THEN                                                
c       number of transitions                                                   
        Nlamtr(is)= RDINP(Equal,1)                                              
        IF (Nlamtr(is).GT.0) THEN                                               
          lamtr(is,1) = RDINP(Equal,1)                                          
          DO i = 2, Nlamtr(is)+1                                                
            lamtr(is,i) = RDINP(NoEqual,1)                                      
            IF (lamtr(is,i).LT.lamtr(is,i-1)) THEN                              
             CALL MSG(6)                                                        
             error = 1                                                          
             GOTO 999                                                           
            END IF                                                              
          END DO                                                                
          klam(is,1) = RDINP(Equal,1)                                           
          IF (Nlamtr(is).GT.1) THEN                                             
            DO i = 2, Nlamtr(is)                                                
              klam(is,i) = RDINP(NoEqual,1)                                     
            END DO                                                              
          END IF                                                                
        ELSE                                                                    
          startyp(is) = 1                                                       
          Tstar = 10000.0                                                       
        END IF                                                                  
      END IF                                                                    
                                                                                
c     spectrum from a file for startyp=4,5,6                                    
c     startyp.EQ.4 - file gives lambda*L_lambda                                 
c     startyp.EQ.5 - file gives L_lambda                                        
c     startyp.EQ.6 - file gives Lnu = lambda**2*L_lambda                        
      IF (startyp(is).EQ.4.OR.startyp(is).EQ.5.OR.startyp(is).EQ.6) THEN        
        strg = 'spectral shape of external radiation:'                          
        CALL FileMSG(nameStar(is),strg)                                         
      END IF                                                                    
                                                                                
999   RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE INPUT(nameIn,nG,nameOut,nameQ,nameNK,TAU1,TAU2,TAUin,          
     &                 Nrec,GridType,Nmodel,error,version)                      
c =======================================================================       
c This subroutine reads input data from the file 'filename'. It                 
c utilizes the subroutine RDINP written by Moshe Elitzur.                       
c                                                      [Z.I., Nov. 1995]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER startyp(2), Nlamtr(2), nBB(2)                                     
      CHARACTER nameStar(2)*235                                                 
      DOUBLE PRECISION Tstar, lamtr(2,101), klam(2,100), Tbb(2,10),             
     &                 rellum(2,10), mu1, ksi, mu2, xSiO, r1rs                  
      COMMON /source/ Tstar, lamtr, klam, Tbb, rellum, mu1, ksi, mu2,           
     &                xSiO, r1rs, startyp, Nlamtr, nBB, nameStar                
      INTEGER psftype, Npsf, iLambda                                            
      DOUBLE PRECISION kPSF(20), FWHM1(20), FWHM2(20), theta1,                  
     &       xpsf(1000), ypsf(1000)                                             
      COMMON /psf1/ iLambda, psftype, Npsf                                      
      COMMON /psf2/ kPSF, FWHM1, FWHM2, Theta1,                                 
     &       xpsf, ypsf                                                         
      CHARACTER LamStr(20)*72, strpow*72, strg*40, version*3                    
      CHARACTER*(*) nameIn, nameOut, nameQ(npG), nameNK(10), namePSF*70,        
     &         nameTAU*70                                                       
      INTEGER iG, nG, Nmodel, i, EtaOK, error, GridType, istop,                 
     &        ioverflw, Nrec, Nmax                                              
c     Nmax is the size of user supplied ETA file                                
      PARAMETER (Nmax = 1000)                                                   
      DOUBLE PRECISION RDINP, TAU1, TAU2, sum, a, b, x(Nmax), e(Nmax),          
     &       aa(Nmax), bb(Nmax), TAUin(Nrec), Ceta, x1, psf1                    
      LOGICAL Equal, NoEqual                                                    
      Equal = .True.                                                            
      NoEqual = .False.                                                         
c -----------------------------------------------------------------------       
      error = 0                                                                 
c     open output file                                                          
      open(12,file=nameOut,STATUS='UNKNOWN')                                    
      write(12,*)'==========================='                                  
      write(12,*)' Output from program Dusty '                                  
      write(12,*)' Version: ',version                                           
      write(12,*)'==========================='                                  
      write(12,*)' '                                                            
      write(12,*)' INPUT parameters from file:'                                 
      write(12,'(2x,a70)')nameIn                                                
      write(12,*)' '                                                            
c     open input file                                                           
      open(1,ERR=998,file=nameIn,STATUS='OLD')                                  
      rewind(1)                                                                 
                                                                                
c     *************************                                                 
c     ** PHYSICAL PARAMETERS **                                                 
c     *************************                                                 
                                                                                
c     1) EXTERNAL RADIATION                                                     
      CALL InpStar(error,1,nameIn)                                              
      IF(error.NE.0) goto 996                                                   
                                                                                
c     2) DUST PROPERTIES                                                        
c     # of different dust grains, to be used in a future version                
c     nG = RDINP(Equal,1)                                                       
      nG = 1                                                                    
c     2.1 Chemical Composition                                                  
c     type of optical properties                                                
      top = RDINP(Equal,1)                                                      
      IF (top.NE.1.AND.top.NE.2.AND.top.NE.3) THEN                              
        CALL MSG(9)                                                             
        error = 1                                                               
        GOTO 999                                                                
      END IF                                                                    
c     for top.LT.3 read in abundances for supported grains                      
      IF (top.LT.3) THEN                                                        
        xC(1) = RDINP(Equal,1)                                                  
        IF (xC(1).LT.0.0) xC(1) = 0.0                                           
        sum = xC(1)                                                             
        DO i = 2, 7                                                             
c         special care to be taken of graphite (1/3-2/3 rule):                  
          IF (i.NE.5) THEN                                                      
            xC(i) = RDINP(NoEqual,1)                                            
            IF (xC(i).LT.0.0) xC(i) = 0.0                                       
c           graphite (perpendicular to c axis) :                                
            IF (i.EQ.4) xC(i) = 2.*xC(i)/3.                                     
          ELSE                                                                  
c           graphite (parallel to c axis) :                                     
            xC(i) = 0.5 * xC(i-1)                                               
          END IF                                                                
          sum = sum + xC(i)                                                     
        END DO                                                                  
      END IF                                                                    
c     user supplied n and k:                                                    
      IF (top.EQ.2) THEN                                                        
        Nfiles = RDINP(Equal,1)                                                 
c       file names                                                              
        strg = 'optical constants:'                                             
        DO i = 1, Nfiles                                                        
          CALL FileMSG(nameNK(i),strg)                                          
        END DO                                                                  
        IF(error.NE.0) goto 996                                                 
c       abundances                                                              
        xCuser(1) = RDINP(Equal,1)                                              
        IF (xCuser(1).LT.0.0) xCuser(1) = 0.0                                   
        sum = sum + xCuser(1)                                                   
        IF (Nfiles.GT.1) THEN                                                   
          DO i = 2, Nfiles                                                      
            xCuser(i) = RDINP(NoEqual,1)                                        
            IF (xCuser(i).LT.0.0) xCuser(i) = 0.0                               
            sum = sum + xCuser(i)                                               
          END DO                                                                
        END IF                                                                  
      END IF                                                                    
      IF (top.LT.3) THEN                                                        
        IF (sum.LE.0.0) THEN                                                    
          CALL MSG(5)                                                           
          error = 1                                                             
          GOTO 999                                                              
        END IF                                                                  
c       normalize abundances for supported grains:                              
        DO i = 1, 7                                                             
          xC(i) = xC(i) / sum                                                   
        END DO                                                                  
c       normalize abundances for user supplied grains                           
        IF (top.EQ.2) THEN                                                      
          DO i = 1, Nfiles                                                      
            xCuser(i) = xCuser(i) / sum                                         
          END DO                                                                
        END IF                                                                  
      END IF                                                                    
c     user supplied cross-sections:                                             
      IF (top.EQ.3) THEN                                                        
c       filename for Qabs and Qsca                                              
        strg= 'abs. and scatt. cross-sections:'                                 
        DO iG = 1, nG                                                           
          CALL FileMSG(nameQ(iG),strg)                                          
        END DO                                                                  
      END IF                                                                    
                                                                                
c     2.2 Grain size distribution                                               
      IF (top.NE.3) THEN                                                        
c       type of size distribution                                               
        szds = RDINP(Equal,1)                                                   
        IF (szds.NE.1.AND.szds.NE.2.AND.szds.NE.3) THEN                         
          CALL MSG(10)                                                          
          error = 1                                                             
          GOTO 999                                                              
        END IF                                                                  
c       grain sizes                                                             
        IF (szds.GT.1) THEN                                                     
          qsd = RDINP(Equal,1)                                                  
          a1 = RDINP(Equal,1)                                                   
          IF (a1.LE.0.0) a1 = 0.0001                                            
          a2 = RDINP(Equal,1)                                                   
          IF (szds.EQ.2.AND.a2.LT.a1) a2 = a1                                   
        ELSE                                                                    
          qsd = 3.5                                                             
          a1 = 0.005                                                            
          a2 = 0.25                                                             
        END IF                                                                  
      END IF                                                                    
                                                                                
c     2.3 Dust temperature on inner boundary                                    
      DO iG = 1, nG                                                             
        Tsub(iG) = RDINP(Equal,1)                                               
      END DO                                                                    
                                                                                
c     3) DENSITY DISTRIBUTION                                                   
                                                                                
c     parameter describing eta function:                                        
      denstyp = RDINP(Equal,1)                                                  
c     WriteOut prints all input data, read so far, in fname.out 	               
c     It has to be placed here b/c geometry (spherical or slab) is not determine
c     until denstyp is read [MN,Sep'99].                                        
      CALL WriteOut(1,nG,nameQ,nameNK)                                          
c     The following transfers are caused by the change of order                 
c     of density types [Dec'96]                                                 
      IF (denstyp.EQ.1) THEN                                                    
        denstyp = 2                                                             
        ELSE                                                                    
        IF (denstyp.EQ.2) THEN                                                  
          denstyp = 3                                                           
          ELSE                                                                  
          IF (denstyp.EQ.3) THEN                                                
            denstyp = 5                                                         
            ELSE                                                                
            IF (denstyp.EQ.4) THEN                                              
              denstyp = 4                                                       
              ELSE                                                              
              IF (denstyp.EQ.5) denstyp = 7                                     
            END IF                                                              
          END IF                                                                
        END IF                                                                  
      END IF                                                                    
c     initialize the logical variable RDW (in 'denstyp.inc', used in ChkSplin)  
      IF(denstyp.eq.5.OR.denstyp.eq.6) THEN                                     
        RDW = .true.                                                            
      ELSE                                                                      
        RDW = .false.		  	                                                      
      END IF                                                                    
                                                                                
      EtaOK = 0                                                                 
      Ntr = 0                                                                   
c     Read parameters for each type of density distribution                     
c     slab geometry                                                             
      IF (denstyp.EQ.0) THEN                                                    
        EtaOK = 1                                                               
        write(12,'(a33)') ' Calculation in planar geometry:'                    
        mu1 = RDINP(Equal,1)                                                    
        IF (mu1.GT.1.0) mu1 = 1.0                                               
        IF (mu1.LT.0.0) mu1 = -1.0                                              
        write(12,'(a35,1p,e11.3)')                                              
     &                ' cos of left illumination angle = ',mu1                  
        ksi = RDINP(Equal,1)                                                    
        IF (ksi.LT.0.0) ksi = 0.0                                               
        IF (ksi.GT.1.0) THEN                                                    
          ksi = 1.0                                                             
        END IF                                                                  
        write(12,'(a6,1p,e11.3)')'  R = ',ksi                                   
        IF(ksi.GT.0.) THEN                                                      
c         in case of additional source on the right                             
          mu2 = RDINP(Equal,1)                                                  
          IF (mu2.GT.1.0) mu2 = 1.0                                             
          IF (mu2.LT.0.0) mu2 = -1.0                                            
          write(12,'(a36,1p,e10.3)')                                            
     &                ' cos of right illumination angle = ',mu2                 
          CALL InpStar(error,2,nameIn)                                          
          CALL WriteOut(2,nG,nameQ,nameNK)                                      
        ELSE                                                                    
c         even if no second source is supplied, mu2 needs a value               
c         of 1 to avoid crashing in the formulae in SLBStar [MN]                
          mu2 = 1.                                                              
        END IF                                                                  
      ELSE                                                                      
        ksi = 0.                                                                
      END IF                                                                    
c     smooth or broken power laws                                               
      IF (denstyp.EQ.1.OR.denstyp.EQ.2) THEN                                    
        EtaOK = 1                                                               
        Ntr = RDINP(Equal,1)                                                    
c       changed definition                                                      
        Ntr = Ntr - 1                                                           
c       read in transition radii                                                
        IF (Ntr.GT.0) THEN                                                      
          Ytr(1) = RDINP(Equal,1)                                               
          IF (Ntr.GT.1) THEN                                                    
            DO i = 2, Ntr                                                       
              Ytr(i) = RDINP(NoEqual,1)                                         
            END DO                                                              
          END IF                                                                
          Yout = RDINP(NoEqual,1)                                               
        ELSE                                                                    
          Yout = RDINP(Equal,1)                                                 
        END IF                                                                  
        IF (Yout.LE.1.0) Yout = 1.001                                           
c       read in powers                                                          
        pow = RDINP(Equal,1)                                                    
        IF (Ntr.GT.0) THEN                                                      
          DO i = 1, Ntr                                                         
            ptr(i) = RDINP(NoEqual,1)                                           
          END DO                                                                
        END IF                                                                  
c       print info to the output file                                           
        IF (Ntr.EQ.0) THEN                                                      
          CALL getfs(pow,2,0,strpow)                                            
          write(12,'(a38,a5)')                                                  
     &      ' Density described by 1/r**k with k =',strpow                      
          write(12,'(a21,1p,e10.3)')'  Relative thickness:',Yout                
        ELSE                                                                    
          IF (denstyp.EQ.1) THEN                                                
            write(12,*)' Density described by smooth power law'                 
          ELSE                                                                  
            write(12,*)' Density described by a broken power law:'              
          END IF                                                                
          write(12,*)'  power   Ytransition'                                    
          write(12,*)'  -------------------'                                    
          write(12,*)'              1.0'                                        
          CALL getfs(pow,2,0,strpow)                                            
          write(12,'(a2,a5)')'  ',strpow                                        
          DO i = 1, Ntr                                                         
            write(12,'(a10,1p,e10.3)')'          ',Ytr(i)                       
            CALL getfs(ptr(i),2,0,strpow)                                       
            write(12,'(a2,a5)')'  ',strpow                                      
          END DO                                                                
          write(12,'(a10,1p,e10.3)')'          ',Yout                           
        END IF                                                                  
      END IF                                                                    
c     exponential law                                                           
      IF (denstyp.EQ.3) THEN                                                    
        EtaOK = 1                                                               
        Yout = RDINP(Equal,1)                                                   
        IF (Yout.LE.1.0) Yout = 1.001                                           
        pow = RDINP(Equal,1)                                                    
        IF (pow.LE.0.0) THEN                                                    
          EtaOK = 0                                                             
          ELSE                                                                  
          write(12,*)' Density described by exponential distribution'           
          write(12,'(a21,1p,e10.3)')'               Sigma:',pow                 
          write(12,'(a21,1p,e10.3)')'  Relative thickness:',Yout                
        END IF                                                                  
      END IF                                                                    
c     default approximation and default numerics for rad. driven winds          
      IF (denstyp.EQ.4.OR.denstyp.EQ.5) THEN                                    
        EtaOK = 1                                                               
        Yout = RDINP(Equal,1)                                                   
        IF (Yout.LE.1.0) Yout = 1.001                                           
c       ** DEFAULT ** for epsilon = v1/ve = u1/ue:                              
        pow = 0.2                                                               
        IF(denstyp.EQ.5) THEN                                                   
c         ** DEFAULT ** for Max(GravCor = Fgrav/Frad_press):                    
          ptr(1) = 0.5                                                          
c         convergence criterion:                                                
          ptr(2) = 1.0                                                          
        END IF                                                                  
        write(12,*)' Density for radiatively driven winds from'                 
        IF (denstyp.EQ.4) THEN                                                  
          write(12,*)' analytic approximation for gray dust.'                   
        ELSE                                                                    
          write(12,*)' full dynamic calculation.'                               
        END IF                                                                  
        write(12,'(a21,1p,e10.3)')'  Relative thickness:',Yout                  
      END IF                                                                    
c     full dynamical calculation for radiatively driven winds                   
      IF (denstyp.EQ.6) THEN                                                    
        EtaOK = 1                                                               
c       Yout:                                                                   
        Yout = RDINP(Equal,1)                                                   
        IF (Yout.LE.1.0) Yout = 1.001                                           
c       u1:                                                                     
        pow = RDINP(Equal,1)                                                    
c       GravCor:                                                                
        ptr(1) = RDINP(Equal,1)                                                 
c       convergence criterion:                                                  
        ptr(2) = 1.0                                                            
c       default limits                                                          
        IF (ptr(1).GT.1.0) ptr(1) = 0.9999                                      
        IF (ptr(1).LT.0.0) ptr(1) = 0.0                                         
        write(12,*)' Density for radiatively driven winds:'                     
        write(12,*)' Full dynamic calculation (denstyp.EQ.6),'                  
        write(12,'(a6,1p,e10.3,a20,e10.3)')                                     
     &                   '  u1 =',pow,' Gamma^{-1}_{max} =',ptr(1)              
      END IF                                                                    
c     user specified table for ETA                                              
      IF(denstyp.EQ.7) THEN                                                     
        EtaOK = 1                                                               
        strg = 'dust density distribution:'                                     
        CALL FileMSG(nameETA,strg)                                              
        write(12,*)' Density distribution supplied from file:'                  
        write(12,'(2x,a70)') nameETA                                            
c       read in the density                                                     
        open(31,ERR=997,file=nameETA,STATUS='OLD')                              
c       three lines in the header:                                              
        DO i = 1, 3                                                             
          read(31,*,ERR=997)                                                    
        END DO                                                                  
        istop = 0                                                               
        i = 0                                                                   
        DO WHILE (istop.GE.0)                                                   
          read(31,*,END=900,ERR=997,iostat=istop)a, b                           
          IF (istop.GE.0) THEN                                                  
            i = i + 1                                                           
            x(i) = a                                                            
            e(i) = b                                                            
            IF (i.EQ.1) x1 = x(i)                                               
            yEta7(i) = x(i) / x1                                                
          END IF                                                                
        END DO                                                                  
900     close(31)                                                               
        Eta7OK = 7                                                              
        nYEta7 = i                                                              
        IF (nYEta7.LT.2) goto 997                                               
c       if input positions in descending order turn them around                 
        IF (yEta7(1).GT.yEta7(2)) THEN                                          
          DO i = 1, nYEta7                                                      
            aa(i) = yEta7(i)                                                    
            bb(i) = e(i)                                                        
          END DO                                                                
          DO i = 1, nYEta7                                                      
            yEta7(i) = aa(nYEta7+1-i)                                           
            e(i) = bb(nYEta7+1-i)                                               
          END DO                                                                
        END IF                                                                  
c       relative thickness                                                      
        Yout = yEta7(nYEta7)                                                    
        write(12,'(a21,1p,e10.3)')'  Relative thickness:',Yout                  
        IF (Yout.LE.1.0) Yout = 1.001                                           
c       integrate and ...                                                       
        CALL SIMPSON(Nmax,1,nYEta7,yEta7,e,Ceta)                                
c       ... renormalize                                                         
        DO i = 1, nYEta7                                                        
          Eta7(i) = e(i) / Ceta                                                 
        END DO                                                                  
      END IF                                                                    
c     done with the reading of density distribution                             
      IF (EtaOK.NE.1) THEN                                                      
        CALL MSG(3)                                                             
        error = 1                                                               
        GOTO 999                                                                
      END IF                                                                    
      write(12,*)' --------------------------------------------'                
                                                                                
c     4) OPTICAL DEPTH                                                          
c     grid type                                                                 
      GridType = RDINP(Equal,1)                                                 
      IF (GridType.EQ.3) THEN                                                   
c     TAU-grid from a file                                                      
        strg = 'user supplied TAU-grid:'                                        
        CALL FileMSG(nameTAU,strg)                                              
c       read optical depths                                                     
        open(32,ERR=992,file=nameTAU,STATUS='OLD')                              
c       fiducial wavelength                                                     
c       (the second argument of RDINP is the unit)                              
        lamfid = RDINP(Equal,32)                                                
        i = 0                                                                   
        DO WHILE (i.LT.Nrec)                                                    
          read(32,*,END=902,ERR=992) a                                          
          i = i + 1                                                             
          TAUin(i) = a                                                          
        END DO                                                                  
902     close(32)                                                               
        Nmodel = i                                                              
c       sort the tau-grid                                                       
        CALL Sort(TAUin,Nmodel)                                                 
        TAU1 = TAUin(1)                                                         
        IF (TAU1.LE.0.0) TAU1 = 0.0001                                          
        TAU2 = TAUin(Nmodel)                                                    
      ELSE                                                                      
c      fiducial wavelength                                                      
       lamfid = RDINP(Equal,1)                                                  
c      total optical depths at lamfid                                           
       TAU1 = RDINP(Equal,1)                                                    
       IF (TAU1.LE.0.0) TAU1 = 0.0001                                           
       TAU2 = RDINP(Equal,1)                                                    
       IF (TAU2.LE.TAU1) THEN                                                   
         TAU2 = TAU1                                                            
         Nmodel = 1                                                             
       END IF                                                                   
c      read number of models                                                    
       Nmodel = RDINP(Equal,1)                                                  
       IF (Nmodel.GT.(Nrec-1)) Nmodel = Nrec-1                                  
       IF (Nmodel.LT.1) Nmodel = 1                                              
      END IF                                                                    
      IF (Nmodel.GT.1) THEN                                                     
        write(12,'(a19,1p,e8.1,a8)')' Optical depths at',lamfid,                
     &                              'microns'                                   
        write(12,'(a14,1p,e9.2,a3,e9.2)')' ranging from',TAU1,' to',TAU2        
        IF (GridType.EQ.1) strg=' models with linear grid    '                  
        IF (GridType.EQ.2) strg=' models with logarithmic grid'                 
        IF (GridType.EQ.3) strg=' models with grid from file  '                 
        write(12,'(a1,i4,a)')' ', Nmodel, strg                                  
        IF (GridType.EQ.3) write(12,'(a4,a70)')'    ',nameTAU                   
      ELSE                                                                      
        write(12,'(a18,1p,e8.1,a9,e9.2)')' Optical depth at',lamfid,            
     &                                  ' microns:',TAU1                        
      END IF                                                                    
                                                                                
c     ************************                                                  
c     ** NUMERICAL ACCURACY **                                                  
c     ************************                                                  
c     some of the parameters were originally left to user to                    
c     specify them, most of them are now unchangeable defaults                  
c     max. increase of scaled TAU/TAUtot                                        
c      delTAUsc = RDINP(Equal,1)                                                
c      IF (delTAUsc.LE.0.0) delTAUsc = 0.3                                      
      delTAUsc = 0.3                                                            
c     max. increase in the ratio of two y                                       
c      facc = RDINP(Equal,1)                                                    
c      IF (facc.LE.0.0) facc = 2.0                                              
      facc = 2.0                                                                
c     max allowed ratio of Eta for two consecutive pts.                         
      EtaRat = 4.0                                                              
c     # of rays per radial step                                                 
c      Nins = RDINP(Equal,1)                                                    
c      IF (Nins.LE.0) Nins = 2                                                  
      Nins = 2                                                                  
c     accuracy for numerical integration in ROMBERG                             
c      accRomb = RDINP(Equal,1)                                                 
c      IF (accRomb.LE.0.0) accRomb = 0.0001                                     
c      IF (accRomb.GT.0.005) accRomb = 0.005                                    
      accRomb = 0.0001                                                          
c     accuracy for flux conservation                                            
      accuracy = RDINP(Equal,1)                                                 
      IF (accuracy.LE.0.0) accuracy = 0.05                                      
c     protect against a very large value for accuracies                         
      IF (accuracy.GT.0.25) accuracy = 0.25                                     
c     accuracy for convergence (typical 0.0001)                                 
      accConv = accuracy / 500.0                                                
      accFbol = 10.*accConv                                                     
c     dynamical range                                                           
c      dynrange = RDINP(Equal,1)                                                
c      IF (dynrange.LE.0.0) dynrange = 1.0e-14                                  
c      IF (dynrange.GT.0.0001) dynrange = 1.0e-4                                
      dynrange = 1.0e-15                                                        
      IF (accuracy.GE.0.1) THEN                                                 
        CALL getfs(accuracy*100,0,1,strpow)                                     
        write(12,'(a20,a3,a1)')' Required accuracy:',strpow,'%'                 
      ELSE                                                                      
        CALL getfs(accuracy*100,0,1,strpow)                                     
        write(12,'(a20,a2,a1)')' Required accuracy:',strpow,'%'                 
      END IF                                                                    
      write(12,*)' --------------------------------------------'                
c     ******************                                                        
c     ** OUTPUT FLAGS **                                                        
c     ******************                                                        
c     A flag for additional output for us only [MN]:                            
c     If iInn=1: print err.vs.iter in unt=38 (fname.err) for all models         
c     and additionally list scaled fbol(y) and Ubol(y) in m-files.              
c     iInn replaces iLTS in 'output.inc'                                        
      iInn = 0                                                                  
c     these two by default                                                      
      iSUM = 1                                                                  
      iOUT = 1                                                                  
      iVerb = RDINP(Equal,1)                                                    
c     spectral properties                                                       
      iSPP = RDINP(Equal,1)                                                     
c     spectra                                                                   
      iA = RDINP(Equal,1)                                                       
c     images (intensity) for spherical only                                     
      IF(denstyp.NE.0) THEN                                                     
       iC = RDINP(Equal,1)                                                      
       IF (iC.NE.0) THEN                                                        
        NlambdaOut = RDINP(Equal,1)                                             
        IF (NlambdaOut.GE.1) THEN                                               
          DO i = 1, NlambdaOut                                                  
            LambdaOut(i) = RDINP(NoEqual,1)                                     
c           make sure the wavelengths are inside Dusty's range                  
            IF (LambdaOut(i).LE.0.01) LambdaOut(i) = 0.01                       
            IF (LambdaOut(i).GT.36000) LambdaOut(i) = 36000                     
          END DO                                                                
          ioverflw = 0                                                          
          DO i = 1, NlambdaOut                                                  
            IF (LambdaOut(i).LT.0.995) THEN                                     
              CALL getfs(LambdaOut(i),2,1,LamStr(i))                            
            ELSE                                                                
              IF (LambdaOut(i).LT.9.95) THEN                                    
               CALL getfs(LambdaOut(i),1,0,LamStr(i))                           
              ELSE                                                              
                IF (LambdaOut(i).LT.99.5) THEN                                  
                 CALL getfs(LambdaOut(i),0,0,LamStr(i))                         
                ELSE                                                            
                 CALL getfs(LambdaOut(i),0,1,LamStr(i))                         
c                  IDINT truncates the numbers                                  
c                  prn = IDINT(LambdaOut(i))                                    
c                  write (LamStr(i),'(I6)') prn                                 
                END IF                                                          
              END IF                                                            
            END IF                                                              
            IF (LambdaOut(i).GT.9999.5) THEN                                    
              ioverflw = 1                                                      
              strpow = LamStr(i)                                                
              strpow(4:4) = '*'                                                 
              strpow(5:5) = ' '                                                 
              LamStr(i) = strpow                                                
            END IF                                                              
          END DO                                                                
        END IF                                                                  
        write(12,*)' Images requested for these wavelengths (mic)'              
        write(12,'(a1,20a5)')' ',(LamStr(i),i=1,NlambdaOut)                     
        IF (ioverflw.EQ.1) write(12,*)'  *: in mm'                              
c       convolved images  (only for our use)                                    
        IF (iC.LT.0) THEN                                                       
         iPSF = RDINP(Equal,1)                                                  
         IF (iPSF.NE.0) THEN                                                    
           Theta1 = RDINP(Equal,1)                                              
           write(12,'(a39,1p,e7.1)')                                            
     &           ' Convolved images produced for theta1=',Theta1                
           psftype = RDINP(Equal,1)                                             
           IF (psftype.NE.1.AND.psftype.NE.2.AND.psftype.NE.3) goto 994         
           IF (psftype.LT.3) THEN                                               
c            Gaussians, read in parameters                                      
c            FWHM for the first component                                       
             FWHM1(1) = RDINP(Equal,1)                                          
             IF (NlambdaOut.GT.1) THEN                                          
               DO i = 2, NlambdaOut                                             
                FWHM1(i) = RDINP(NoEqual,1)                                     
               END DO                                                           
             END IF                                                             
             IF (psftype.EQ.2) THEN                                             
c             relative strength for the second component                        
              kPSF(1) = RDINP(Equal,1)                                          
              IF (NlambdaOut.GT.1) THEN                                         
                DO i = 2, NlambdaOut                                            
                  kPSF(i) = RDINP(NoEqual,1)                                    
                END DO                                                          
              END IF                                                            
c             FWHM for the second component                                     
              FWHM2(1) = RDINP(Equal,1)                                         
              IF (NlambdaOut.GT.1) THEN                                         
                DO i = 2, NlambdaOut                                            
                  FWHM2(i) = RDINP(NoEqual,1)                                   
                END DO                                                          
              END IF                                                            
             END IF                                                             
             write(12,*)' The point spread functions are Gaussians'             
           ELSE                                                                 
c           user supplied PSF                                                   
            strg = 'point spread function:'                                     
            CALL FileMSG(namePSF,strg)                                          
            write(12,*)' The point spread function supplied from file'          
            write(12,'(2x,a70)')namePSF                                         
            open(31,ERR=995,file=namePSF,STATUS='OLD')                          
c           three lines in the header:                                          
            DO i = 1, 3                                                         
              read(31,*,ERR=995)                                                
            END DO                                                              
            istop = 0                                                           
            i = 0                                                               
            DO WHILE (istop.GE.0)                                               
              read(31,*,END=901,ERR=995,iostat=istop)a, b                       
              IF (istop.GE.0) THEN                                              
                i = i + 1                                                       
                IF (i.EQ.1) THEN                                                
                  psf1 = b                                                      
                  IF (a.NE.0.0) goto 995                                        
                END IF                                                          
                xpsf(i) = a                                                     
                ypsf(i) = b / psf1                                              
              END IF                                                            
            END DO                                                              
901         close(31)                                                           
            Npsf = i                                                            
c           scale to 1 at the center                                            
            CALL ScaleTo1(1000,Npsf,ypsf)                                       
c           find equivalent FWHM                                                
            istop = 0                                                           
            i = 1                                                               
            DO WHILE (istop.EQ.0)                                               
              i = i + 1                                                         
              IF (ypsf(i).LE.0.5) istop = 1                                     
            END DO                                                              
c           linear interpolation                                                
            FWHM1(1) = (xpsf(i)-xpsf(i-1))/(ypsf(i)-ypsf(i-1))                  
            FWHM1(1) = (FWHM1(1)*(0.5-ypsf(i-1))+xpsf(i-1))*2                   
            FWHM2(1) = 0.0                                                      
            write(12,'(a18,1p,e8.1)')' Equivalent FWHM:',FWHM1(1)               
           END IF                                                               
         END IF                                                                 
        END IF                                                                  
c       visibility (only if the intensity is requested)                         
        iV = RDINP(Equal,1)                                                     
        IF(iV.NE.0) iV = iC                                                     
        write(12,*)' --------------------------------------------'              
       ELSE                                                                     
        iPSF = 0                                                                
        iV = 0                                                                  
       END IF                                                                   
      ELSE                                                                      
       iC = 0                                                                   
      END IF                                                                    
      write(12,*)' '                                                            
c     radial quantities                                                         
      iB = RDINP(Equal,1)                                                       
c     run-time messages                                                         
      iX = RDINP(Equal,1)                                                       
      iINP = 0                                                                  
c     *** DONE ***                                                              
c     if everything is OK, close the input file and finish                      
999   goto 996                                                                  
c     or in the case of err reading files...                                    
992   write(12,*)' ***  FATAL ERROR IN DUSTY  *************'                    
      write(12,*)' File with user supplied TAU-grid:'                           
      write(12,'(2x,a70)') nameTAU                                              
      write(12,*)' is not properly formatted?!'                                 
      write(12,*)' ****************************************'                    
      close(12)                                                                 
      error = 3                                                                 
      goto 996                                                                  
994   CALL MSG(12)                                                              
      close(12)                                                                 
      error = 3                                                                 
      goto 996                                                                  
995   write(12,*)' ***  FATAL ERROR IN DUSTY  *************'                    
      write(12,*)' File with the point spread function:'                        
      write(12,'(a2,a70)')'  ', namePSF                                         
      write(12,*)' is not properly formatted?!'                                 
      write(12,*)' ****************************************'                    
      close(12)                                                                 
      error = 3                                                                 
997   write(12,*)' ***  FATAL ERROR IN DUSTY  *************'                    
      write(12,*)' File with the dust density distribution:'                    
      write(12,'(2x,a70)') nameETA                                              
      write(12,*)' is not properly formatted?!'                                 
      write(12,*)' ****************************************'                    
      close(12)                                                                 
      error = 3                                                                 
998   write(12,*)' ***  FATAL ERROR IN DUSTY  ****'                             
      write(12,*)' Input file:'                                                 
      write(12,'(2x,a70)') nameIn                                               
      write(12,*)' is missing?!'                                                
      write(12,*)' *******************************'                             
      close(12)                                                                 
      error = 3                                                                 
c -----------------------------------------------------------------------       
996   close(1)                                                                  
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE OPPEN(model,RootName,length)                                   
c =======================================================================       
c This subroutine prints the results out.              [Z.I., Feb. 1996]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER startyp(2), Nlamtr(2), nBB(2)                                     
      CHARACTER nameStar(2)*235                                                 
      DOUBLE PRECISION Tstar, lamtr(2,101), klam(2,100), Tbb(2,10),             
     &                 rellum(2,10), mu1, ksi, mu2, xSiO, r1rs                  
      COMMON /source/ Tstar, lamtr, klam, Tbb, rellum, mu1, ksi, mu2,           
     &                xSiO, r1rs, startyp, Nlamtr, nBB, nameStar                
      INTEGER psftype, Npsf, iLambda                                            
      DOUBLE PRECISION kPSF(20), FWHM1(20), FWHM2(20), theta1,                  
     &       xpsf(1000), ypsf(1000)                                             
      COMMON /psf1/ iLambda, psftype, Npsf                                      
      COMMON /psf2/ kPSF, FWHM1, FWHM2, Theta1,                                 
     &       xpsf, ypsf                                                         
      CHARACTER ch5*5, RootName*(*), fname*235                                  
      CHARACTER*72 header1, s3, s4                                              
      INTEGER model, length, i                                                  
c -----------------------------------------------------------------------       
c     set up the status indicators                                              
      iERROR = 0                                                                
      iWARNING = 0                                                              
      IF (model.EQ.1) iCUMM = 0                                                 
c     The following files pertain to ALL models and are open if model.EQ.1      
      IF (model.EQ.1) THEN                                                      
c      the header to output file *.OUT is moved to PrOut [MN,Sep'99]            
                                                                                
c      open file with spectral properties RootName.SPP                          
       IF (iSPP.NE.0) THEN                                                      
        CALL Attach(RootName,length,'.spp',fname)                               
        open(19,file=fname,STATUS='UNKNOWN')                                    
        IF(denstyp.EQ.0) THEN                                                   
          write(19,'(a49)')                                                     
     &              '# ==============================================='         
          write(19,'(a49)')                                                     
     &              '# Properties of Spectra from the slab right side '         
          write(19,'(a49)')                                                     
     &              '# -----------------------------------------------'         
        ELSE                                                                    
          call line(1,2,19)                                                     
          write(19,'(a23)')'#  Spectral Properties '                            
          call line(1,1,19)                                                     
        END IF                                                                  
         s3='###   tau0      Psi      fV       fK       f12    C21  '           
         s4=' C31   C43  b8-13 b14-22 B9.8 B11.4  R9.8-18  '                    
         write(19,'(a55,a46)')s3,s4                                             
        IF (denstyp.EQ.0.AND.iSPP.EQ.3) THEN                                    
        CALL Attach(RootName,length,'.zpp',fname)                               
          open(24,file=fname,STATUS='UNKNOWN')                                  
          write(24,'(a49)')                                                     
     &              '# ==============================================='         
          write(24,'(a49)')                                                     
     &              '# Properties of Spectra from the slab left side  '         
          write(24,'(a49)')                                                     
     &              '# -----------------------------------------------'         
         s3='###   tau0      Psi      fV       fK       f12    C21  '           
         s4=' C31   C43  b8-13 b14-22 B9.8 B11.4  R9.8-18  '                    
           write(24,'(a55,a46)')s3,s4                                           
        END IF                                                                  
       END IF                                                                   
c       open file for point spread function                                     
        IF (iPSF.NE.0.AND.psftype.LT.3) THEN                                    
c         wavelength dependent PSF are also printed out                         
          CALL Attach(RootName,length,'.psf',fname)                             
          open(23,file=fname,STATUS='UNKNOWN')                                  
          header1 = '    Offset'                                                
          write(23,'(a10,20f10.2)')header1,(lambdaOut(i),i=1,NlambdaOut)        
        END IF                                                                  
c       Spectra for all models in one file '*.stb' if flag=1                    
        IF(iA.eq.1) THEN                                                        
         CALL Attach(RootName,length,'.stb',fname)                              
         open(15,file=fname,STATUS='UNKNOWN')                                   
        END IF                                                                  
c       All radial profiles in one file '*.rtb' if flag=1                       
        IF(iB.eq.1) THEN                                                        
         CALL Attach(RootName,length,'.rtb',fname)                              
         open(16,file=fname,STATUS='UNKNOWN')                                   
        END IF                                                                  
c       All imaging files in  '*.itb' if flag=1                                 
        IF (abs(iC).eq.1) THEN                                                  
         CALL Attach(RootName,length,'.itb',fname)                              
         open(17,file=fname,STATUS='UNKNOWN')                                   
        END IF                                                                  
        IF(iX.eq.1) THEN                                                        
         CALL Attach(RootName,length,'.mtb',fname)                              
         open(18,file=fname,STATUS='UNKNOWN')                                   
        END IF                                                                  
        IF(iInn.eq.1) THEN                                                      
         CALL Attach(RootName,length,'.err',fname)                              
         open(38,file=fname,STATUS='UNKNOWN')                                   
        END IF                                                                  
c     end if for model=1                                                        
      END IF                                                                    
c -------------------------------------------------------------                 
c     the following files are open for EVERY model                              
                                                                                
c     (the headers for .s## and .r## files are moved to PrOut, MN)              
c     open the spectrum file RootName.s##  (## = model number)                  
      IF(iA.GT.1) THEN                                                          
       write(ch5,'(a2,I3.3)') '.s', model                                       
       CALL Attach(RootName,length,ch5,fname)                                   
       open(15,file=fname,STATUS='UNKNOWN')                                     
       IF(denstyp.EQ.0) THEN                                                    
         IF(iA.eq.3) THEN                                                       
           write(ch5,'(a2,I3.3)') '.z', model                                   
           CALL Attach(RootName,length,ch5,fname)                               
           open(25,file=fname,STATUS='UNKNOWN')                                 
         END IF                                                                 
       END IF                                                                   
      END IF                                                                    
c     open the file RootName.r## (y-dependent quantities)                       
      IF(iB.GT.1) THEN                                                          
        write(ch5,'(a2,I3.3)') '.r', model                                      
        CALL Attach(RootName,length,ch5,fname)                                  
        open(16,file=fname,STATUS='UNKNOWN')                                    
      END IF                                                                    
c     open the file RootName.i## (surface brightness)                           
      IF(abs(iC).GT.1) THEN                                                     
        write(ch5,'(a2,I3.3)') '.i', model                                      
        CALL Attach(RootName,length,ch5,fname)                                  
        open(17,file=fname,STATUS='UNKNOWN')                                    
      END IF                                                                    
c     open the file RootName.c## (convolved images) for flag iC<0               
c      if iC=-3 - in a separate file fname.c##                                  
      IF(iC.EQ.-3.AND.iPSF.GT.0) THEN                                           
        write(ch5,'(a2,I3.3)') '.c', model                                      
        CALL Attach(RootName,length,ch5,fname)                                  
        open(21,file=fname,STATUS='UNKNOWN')                                    
      END IF                                                                    
c     open the file RootName.v## (visibility curves)                            
      IF((abs(iC).EQ.3).AND.(iV.GT.0)) THEN                                     
        write(ch5,'(a2,I3.3)') '.v', model                                      
        CALL Attach(RootName,length,ch5,fname)                                  
        open(22,file=fname,STATUS='UNKNOWN')                                    
      END IF                                                                    
c     open the output file RootName.m##                                         
      IF (iX.GT.1) THEN                                                         
        write(ch5,'(a2,I3.3)') '.m', model                                      
        CALL Attach(RootName,length,ch5,fname)                                  
        open(18,file=fname,STATUS='UNKNOWN')                                    
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE PROUT(model)                                                   
c =======================================================================       
c This subroutine prints the results out.        [ZI,Feb'96; MN,Mar'99]         
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      DOUBLE PRECISION TAUslb(npL,npY),fsbol(npY),fpbol(npY),fmbol(npY),        
     &                 AveDev, RMS, maxFerr, fmax, fmed                         
      COMMON /slab/ TAUslb, fsbol, fpbol, fmbol, AveDev, RMS, maxFerr,          
     &                 fmax, fmed                                               
      INTEGER nYok, nPok                                                        
      DOUBLE PRECISION                                                          
     &       Ude(npL,npY), Uds(npL,npY), Us(npL,npY), fde(npL,npY),             
     &       fds(npL,npY), fs(npL,npY), Utot(npL,npY), ftot(npL,npY),           
     &       Td(npG,npY), Ubol(npY), fbol(npY), Uchck(npL,npY),                 
     &       UbolChck(npY), Spectrum(npL), SpecChar(30,99), tauF(npY),          
     &       Intens(npL,npP+2), IntOut(20,npP+2), SmC(30,99),                   
     &       Yok(npY), Pok(npP), tauOut(npL), tauZout(npP+2), F1,               
     &       tr(npY), rg(npG,npY), fsL(npL,npY), fsR(npL,npY), Eps(npY)         
      COMMON /solution/ Ude, Uds, Us, fde, fds, fs, Utot, ftot, Td,             
     &       Ubol, fbol, Uchck, UbolChck, Spectrum, SpecChar, tauF,             
     &       Intens, IntOut, SmC, Yok, Pok, tauOut, tauZout, F1,                
     &       tr, rg, fsL, fsR, Eps, nYok, nPok                                  
      INTEGER startyp(2), Nlamtr(2), nBB(2)                                     
      CHARACTER nameStar(2)*235                                                 
      DOUBLE PRECISION Tstar, lamtr(2,101), klam(2,100), Tbb(2,10),             
     &                 rellum(2,10), mu1, ksi, mu2, xSiO, r1rs                  
      COMMON /source/ Tstar, lamtr, klam, Tbb, rellum, mu1, ksi, mu2,           
     &                xSiO, r1rs, startyp, Nlamtr, nBB, nameStar                
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      DOUBLE PRECISION ugas(npY), qF(npY), vrat(npG,npY), Gamma(npY),           
     &       I1, I2, I3, CMdot, Cve, CM, Cr1                                    
      COMMON /dyn/ ugas, qF, vrat, Gamma, I1, I2, I3, CMdot, Cve, CM,           
     &       Cr1                                                                
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER psftype, Npsf, iLambda                                            
      DOUBLE PRECISION kPSF(20), FWHM1(20), FWHM2(20), theta1,                  
     &       xpsf(1000), ypsf(1000)                                             
      COMMON /psf1/ iLambda, psftype, Npsf                                      
      COMMON /psf2/ kPSF, FWHM1, FWHM2, Theta1,                                 
     &       xpsf, ypsf                                                         
      CHARACTER*72 SC21,SC31,SC43,SB98,SB11,Sbet1,Sbet2,STemp,Serr,             
     &             hdint, hdcon, hdvis, s1, su1, s2, su2, Tstr*10               
      CHARACTER*132 hdsp1,hdsp2,hdrslb1,hdrslb2,hdrsph1,hdrsph2,hdrdyn          
      INTEGER iY, iL, i, model, j, unt                                          
      DOUBLE PRECISION PSFN, PSFfunc(20,1000), Elems(25,200),  ETA,             
     &        faux(npL), fsRbol(npY), tht1, xs, xds, xde, g, res, fnorm         
c ----------------------------------------------------------------------        
c     SmC(1..5, model) are found in ANALYSIS                                    
      SmC(6,model) = Td(1,nY)                                                   
c     SmC(7,model) and SmC(8,model) are Teff(L) and Teff(R), respectively,      
c     the Teff of the slab illumination sources, found from F1 in Analysis      
      tht1 = 412.6/(dsqrt(F1))                                                  
c     colors                                                                    
      CALL getFS(SpecChar(9,model),2,0,SC21)                                    
      CALL getFS(SpecChar(10,model),2,0,SC31)                                   
      CALL getFS(SpecChar(11,model),2,0,SC43)                                   
c     error in %                                                                
      IF (SmC(5,model).LT.0.1) THEN                                             
        CALL getFS(SmC(5,model)*100,0,0,Serr)                                   
      ELSE IF (SmC(5,model).LT.1.0) THEN                                        
        CALL getFS(SmC(5,model)*100,0,1,Serr)                                   
      ELSE                                                                      
        CALL getFS(SmC(5,model)*100,0,2,Serr)                                   
      END IF                                                                    
c     dust temperature at y=Y                                                   
      IF (SmC(6,model).LT.99.5) THEN                                            
        CALL getFS(SmC(6,model),0,0,STemp)                                      
      ELSE                                                                      
        CALL getFS(SmC(6,model),0,1,STemp)                                      
      END IF                                                                    
                                                                                
c --------------  overall parameters to *.OUT file -----------------------      
c     write header to output file *.OUT                                         
      IF (model.EQ.1) THEN                                                      
       write(12,*)'         '                                                   
       write(12,*)' RESULTS:'                                                   
       write(12,*)' --------'                                                   
       IF (denstyp.EQ.0) THEN                                                   
c     ---------- slab output ----------------                                   
        s1=' ###   tau0    Fe1(W/m2)   f1     r1(cm)  Td(K)  Te(L)'             
       su1=' ###     1        2         3       4       5      6  '             
        IF(ksi.GT.0) THEN                                                       
          s2='   Te(R)   err'                                                   
         su2='     7      8 '                                                   
         write(12,'(a54,a14)')s1,s2                                             
         write(12,'(a54,a14)')su1,su2                                           
         su1=' ====================================================='           
         su2='=============='                                                   
         write(12,'(a54,a14)')su1,su2                                           
        ELSE                                                                    
          s2='  err '                                                           
         su2='   7  '                                                           
         write(12,'(a54,a6)')s1,s2                                              
         write(12,'(a54,a6)')su1,su2                                            
         su1=' ====================================================='           
         su2='====='                                                            
         write(12,'(a54,a6)')su1,su2                                            
        END IF                                                                  
       ELSE                                                                     
c     ------------- output for sphere -----------------------------             
        s1=' ###   tau0   F1(W/m2)  r1(cm)    r1/rc   theta1  Td(Y) err'        
       su1=' ###     1        2        3        4        5      6    7 '        
        IF(denstyp.eq.4.OR.RDW) THEN                                            
        s2='   Mdot      Ve       M> '                                          
       su2='     8        9       10 '                                          
         write(12,'(a59,a25)')s1,s2                                             
         write(12,'(a59,a25)')su1,su2                                           
       su1=' =========================================================='        
       su2='============================'                                       
         write(12,'(a59,a28)')su1,su2                                           
        ELSE                                                                    
         write(12,'(a59)')s1                                                    
         write(12,'(a59)')su1                                                   
       su1=' =========================================================='        
         write(12,'(a59)')su1                                                   
        END IF                                                                  
       END IF                                                                   
      END IF                                                                    
c     print output tables                                                       
c     for slab:                                                                 
      IF(denstyp.EQ.0) THEN                                                     
      IF (fmed.LE.accFbol) fmed = 0.                                            
        IF(ksi.NE.0.) THEN                                                      
          CALL Bolom(fsR,fsRbol)                                                
          g = ksi*fsRbol(1) + fmbol(1)                                          
        ELSE                                                                    
          g = fmbol(1)                                                          
        END IF                                                                  
        IF (ksi.GT.0) THEN                                                      
         IF (SmC(6,model).GE.999.5) THEN                                        
          write(12,'(i4,1p,2E9.2,E10.2,E9.2,a5,1p,2E9.2,a4)')                   
     &    model, TAUfid, F1, fmed, Cr1, STemp, SmC(7,model),                    
     &                                         SmC(8,model), Serr               
         ELSE                                                                   
          write(12,'(i4,1p,2E9.2,E10.2,E9.2,a1,a4,1p,2E9.2,a4)')                
     &    model, TAUfid, F1, fmed, Cr1, ' ',STemp, SmC(7,model),                
     &                                             SmC(8,model), Serr           
         END IF                                                                 
        ELSE                                                                    
         IF (SmC(6,model).GE.999.5) THEN                                        
          write(12,'(i4,1p,2E9.2,E10.2,E9.2,a5,1p,E9.2,a3)')                    
     &    model, TAUfid, F1, fmed, Cr1, STemp, SmC(7,model), Serr               
         ELSE                                                                   
          write(12,'(i4,1p,2E9.2,E10.2,E9.2,a1,a4,1p,E9.2,a3)')                 
     &    model, TAUfid, F1, fmed, Cr1,' ', STemp, SmC(7,model), Serr           
         END IF                                                                 
        END IF                                                                  
                                                                                
      ELSE                                                                      
c     for spherical shell                                                       
        IF (SmC(6,model).LT.9.5) THEN                                           
          IF (denstyp.eq.4.OR.RDW) THEN                                         
           write(12,'(i4,1p,5E9.2,a2,a3,a1,a3,1p,3E9.2)')                       
     &     model, TAUfid, F1, Cr1, r1rs, tht1,'  ', STemp,' ',Serr,             
     &                                    CMdot, CVe, CM                        
          ELSE                                                                  
           write(12,'(i4,1p,5E9.2,a2,a4,a3)')                                   
     &     model, TAUfid, F1, Cr1, r1rs, tht1,'  ', STemp, Serr                 
          END IF                                                                
        ELSE                                                                    
         IF (denstyp.eq.4.OR.RDW) THEN                                          
          write(12,'(i4,1p,5E9.2,2(a1,a4),1p,3E9.2)')                           
     &    model, TAUfid, F1, Cr1, r1rs, tht1,' ',STemp,' ', Serr,               
     &                                    CMdot, CVe, CM                        
         ELSE                                                                   
          write(12,'(i4,1p,5E9.2,2(a1,a4))')                                    
     &    model, TAUfid, F1, Cr1, r1rs, tht1,' ',STemp,' ', Serr                
         END IF                                                                 
        END IF                                                                  
                                                                                
       IF (startyp(1).EQ.1.OR.startyp(1).EQ.2) THEN                             
        IF(Tstar.LT.Te_min) THEN                                                
          CALL GetFS(Tstar,0,1,Tstr)                                            
          write(12,'(a50,a5,a5)')                                               
     &' ** WARNING: the input spectrum is a black-body at ',Tstr,' K **'        
          CALL GetFS(Te_min,0,1,Tstr)                                           
          write(12,'(a50,a5,a5)')                                               
     &' *The point-source assumption requires min Teff of ',Tstr,' K **'        
        END IF                                                                  
       END IF                                                                   
c     end if for geometry                                                       
      END IF                                                                    
                                                                                
      CALL getFS(SpecChar(1,model),2,0,SB98)                                    
      CALL getFS(SpecChar(2,model),2,0,SB11)                                    
      CALL getFS(SpecChar(4,model),2,0,Sbet1)                                   
      CALL getFS(SpecChar(5,model),2,0,Sbet2)                                   
c --------------     spectral properties to *.SPP file   -----------------------
      IF (iSPP.NE.0) THEN                                                       
       write(19,'(i3,1p,5E9.2,7a6,E9.2)') model, TAUfid, SmC(1,model),          
     &      (SpecChar(i,model),i=6,8), SC21, SC31, SC43, Sbet1, Sbet2,          
     &       SB98, SB11, SpecChar(3,model)                                      
      END IF                                                                    
                                                                                
c --------------   spectrum to *.s## (old *.Axx) file   ------------------------
      IF (iA.NE.0) THEN                                                         
       IF(denstyp.EQ.0) THEN                                                    
        hdsp1 = '#  lambda     fRight       xAtt       xDs   '                  
       ELSE                                                                     
        hdsp1 = '#  lambda      fTot        xAtt       xDs   '                  
       END IF                                                                   
        hdsp2 = '     xDe       fInp       tauT      albedo  '                  
        unt = 15                                                                
        call line(1,2,unt)                                                      
        IF(denstyp.EQ.0) THEN                                                   
          write(unt,'(a8,i4,a36)')                                              
     &         '#  model',model,'   SPECTRUM from the right slab side '         
        ELSE                                                                    
          write(unt,'(a8,i4,a12)') '#  model',model,'   SPECTRUM '              
        END IF                                                                  
        call line(1,1,unt)                                                      
        DO iL = 1, nL                                                           
c         The right-side spectra are: ftot + ksi*fsR = fsL + fds + fde          
          ftot(iL,nYok) = ftot(iL,nYok) + ksi*fsR(iL,nYok)                      
          faux(iL) = ftot(iL,nYok)/lambda(iL)                                   
        END DO                                                                  
        CALL Simpson(npL,1,nL,lambda,faux,res)                                  
c       Normalization factor for slab s-spectra                                 
        fnorm = res	                                                            
        DO iL = 1, nL                                                           
           IF (ftot(iL,nYok).NE.0.0) THEN                                       
             xs = fsL(iL,nYok)/ftot(iL,nYok)                                    
             xds = fds(iL,nYok)/ftot(iL,nYok)                                   
             xde = fde(iL,nYok)/ftot(iL,nYok)                                   
           ELSE                                                                 
             xs = 0.0                                                           
             xds = 0.0                                                          
             xde = 0.0                                                          
           END IF                                                               
c          change introduced for version 11a - no need to print negligible value
           IF (ftot(iL,nYok).LT.1.0E-20) ftot(iL,nYok) = 0.0                    
           IF (xs.LT.1.0E-20) xs = 0.0                                          
           IF (xds.LT.1.0E-20) xds = 0.0                                        
           IF (xde.LT.1.0E-20) xde = 0.0                                        
           IF (fsL(iL,1).LT.1.0E-20) fsL(iL,1) = 0.0                            
c          for slab rescale fTot with the bolom flux                            
           IF (denstyp.eq.0) ftot(iL,nYok) = ftot(iL,nYok)/fnorm                
           Elems(1,iL) = lambda(iL)                                             
           Elems(2,iL) = fTot(iL,nYok)                                          
           Elems(3,iL) = xs                                                     
           Elems(4,iL) = xds                                                    
           Elems(5,iL) = xde                                                    
           Elems(6,iL) = fsL(iL,1)                                              
           Elems(7,iL) = TAUtot(iL)                                             
           Elems(8,iL) = omega(iL,1)                                            
        END DO                                                                  
c     ------ Tabulate the spectra in the desired form ----------                
        write(unt,'(2(a45))') hdsp1,hdsp2                                       
        CALL MakeTable(Elems,nL,8,unt)                                          
c       spectra from the illuminated slab side (file *.z##)                     
        IF (denstyp.EQ.0) THEN                                                  
           DO iL = 1, nL                                                        
c           The left-side spectra are: |R*fsR + fm|=|fsL-ftot|                  
            ftot(iL,1) = dabs(fsL(iL,1) - ftot(iL,1))                           
c           normalization of z-spectra for slab                                 
            faux(iL) = ftot(iL,1)/lambda(iL)                                    
           END DO                                                               
           CALL Simpson(npL,1,nL,lambda,faux,res)                               
           fnorm = res	                                                         
           DO iL = 1, nL                                                        
             IF (ftot(iL,1).NE.0.0) THEN                                        
                xs =  ksi*fsR(iL,1)/ftot(iL,1)                                  
                xds = dabs(fds(iL,1)/ftot(iL,1))                                
                xde = dabs(fde(iL,1)/ftot(iL,1))                                
             ELSE                                                               
                xs = 0.0                                                        
                xds = 0.0                                                       
                xde = 0.0                                                       
             END IF                                                             
             IF (xs.LT.1.0E-20) xs =0.0                                         
             IF (xds.LT.1.0E-20) xds =0.0                                       
             IF (xde.LT.1.0E-20) xde =0.0                                       
             IF(fsR(iL,nYok).LT.1.0E-20) fsR(iL,nYok) = 0.0                     
             IF(ftot(iL,1).LT.1.0E-20) ftot(iL,1) = 0.0                         
c            rescale fTot with the bolom flux for z-spectra                     
             ftot(iL,1) = ftot(iL,1)/fnorm                                      
             Elems(1,iL) = lambda(iL)                                           
             Elems(2,iL) = fTot(iL,1)                                           
             Elems(3,iL) = xs                                                   
             Elems(4,iL) = xds                                                  
             Elems(5,iL) = xde                                                  
             Elems(6,iL) = ksi*fsR(iL,nYok)                                     
             Elems(7,iL) = TAUtot(iL)                                           
             Elems(8,iL) = omega(iL,1)                                          
           END DO                                                               
           IF (iA.EQ.3) unt=25                                                  
c          append to the .s## file or write in a separate .z## file (if iA=3)   
           call line(1,1,unt)                                                   
           write(unt,'(a8,i4,a36)')                                             
     &         '#  model',model,'   SPECTRUM from the left slab side'           
           call line(1,1,unt)                                                   
           hdsp1 = '#  lambda      fLeft       xAtt       xDs   '               
           write(unt,'(2(a45))') hdsp1,hdsp2                                    
           CALL MakeTable(Elems,nL,8,unt)                                       
        END IF                                                                  
      END IF                                                                    
                                                                                
c     spectral properties for *.ZPP file in slab case                           
      IF (denstyp.eq.0.AND.iSPP.NE.0) THEN                                      
        CALL getFS(SpecChar(12,model),2,0,SB98)                                 
        CALL getFS(SpecChar(13,model),2,0,SB11)                                 
        CALL getFS(SpecChar(15,model),2,0,Sbet1)                                
        CALL getFS(SpecChar(16,model),2,0,Sbet2)                                
        CALL getFS(SpecChar(20,model),2,0,SC21)                                 
        CALL getFS(SpecChar(21,model),2,0,SC31)                                 
        CALL getFS(SpecChar(22,model),2,0,SC43)                                 
c       write spectral properties to *.ZPP file                                 
         IF (iSPP.EQ.3) THEN                                                    
            write(24,'(i3,1p,5E9.2,7a6,E9.2)')                                  
     &      model, TAUfid, SmC(1,model), (SpecChar(i,model),i=17,19),           
     &      SC21,SC31,SC43, Sbet1,Sbet2, SB98,SB11,SpecChar(14,model)           
         ELSE                                                                   
            write(zline(model),'(i3,1p,5E9.2,7a6,E9.2)')                        
     &      model, TAUfid, SmC(1,model), (SpecChar(i,model),i=17,19),           
     &      SC21,SC31,SC43, Sbet1,Sbet2,SB98,SB11,SpecChar(14,model)            
         END IF                                                                 
      END IF                                                                    
                                                                                
c -----------  y-dependent quantities to *.r## (old *.Bxx) file -------------   
      IF (iB.NE.0) THEN                                                         
       hdrslb1= '#     t        tauF      epsilon       Td '                    
       hdrslb2= '      febol      fRbol      fLbol '                            
       hdrsph1= '#     y         eta         t         tauF '                   
       hdrsph2= '     epsilon      Td         rg '                              
        hdrdyn= '         u        drift'                                       
       unt = 16                                                                 
       call line(1,2,unt)                                                       
       IF(denstyp.EQ.0) THEN                                                    
          write(unt,'(a8,i4,a18)') '#  model',model,'  SPATIAL PROFILES'        
       ELSE                                                                     
          write(unt,'(a8,i4,a18)') '#  model',model,'  RADIAL PROFILES '        
       END IF                                                                   
       call line(1,1,unt)                                                       
c    --------- for slab ---------                                               
       IF (denstyp.EQ.0) THEN                                                   
         DO iY = 1, nYok                                                        
           Elems(1,iY) = tr(iY)                                                 
           Elems(2,iY) = tauF(iY)                                               
           Elems(3,iY) = Eps(iY)                                                
           Elems(4,iY) = Td(1,iY)                                               
           Elems(5,iY) = fsbol(iY)                                              
           Elems(6,iY) = fpbol(iY)                                              
           Elems(7,iY) = fmbol(iY)                                              
         END DO                                                                 
c      check for nonsense values:                                               
         DO i = 1, 8                                                            
          DO iY = 1, nYok                                                       
           IF (Elems(i,iY).NE.Elems(i,iY)) THEN                                 
              Elems(i,iY) = 0.0                                                 
           END IF                                                               
          END DO                                                                
         END DO                                                                 
         write(unt,'(a42,a34)') hdrslb1,hdrslb2                                 
         CALL MakeTable(Elems,nYok,7,unt)                                       
       ELSE                                                                     
c     ------  for spherical shell --------                                      
        DO iY = 1, nYok                                                         
          Elems(1,iY) = Yok(iY)                                                 
          Elems(2,iY) = ETA(Y(iY))                                              
          Elems(3,iY) = tr(iY)                                                  
          Elems(4,iY) = tauF(iY)                                                
          Elems(5,iY) = Eps(iY)                                                 
          Elems(6,iY) = Td(1,iY)                                                
          Elems(7,iY) = rg(1,iY)                                                
        END DO                                                                  
c       check for nonsense values:                                              
        DO i = 1, 7                                                             
         DO iY = 1, nYok                                                        
          IF(Elems(i,iY).NE.Elems(i,iY).OR.Elems(i,iY).LT.1.e-20) THEN          
            Elems(i,iY) = 0.0                                                   
          END IF                                                                
         END DO                                                                 
        END DO                                                                  
c       with dynamics                                                           
        IF (RDW) THEN                                  
         DO iY = 1, nYok                                                        
           Elems(8,iY) = ugas(iY)/ugas(nYok)                                    
           Elems(9,iY) = vrat(1,iY)                                             
         END DO                                                                 
c        check values:                                                          
         DO i = 8, 9                                                            
          DO iY = 1, nYok                                                       
           IF(Elems(i,iY).NE.Elems(i,iY).OR.Elems(i,iY).LT.1.e-20) THEN         
             Elems(i,iY) = 0.0                                                  
           END IF                                                               
          END DO                                                                
         END DO                                                                 
         write(unt,'(a42,a32,a23)') hdrsph1,hdrsph2,hdrdyn                      
         CALL MakeTable(Elems,nYok,9,unt)                                       
        ELSE                                                                    
         write(unt,'(a42,a32)') hdrsph1,hdrsph2                                 
         CALL MakeTable(Elems,nYok,7,unt)                                       
        END IF                                                                  
c      end if for geometry                                                      
       END IF                                                                   
c     end if for the iB (radial) flag                                           
      END IF                                                                    
                                                                                
c --------------   intensities to *.inn (old *.Cxx) file  --------------        
      IF (abs(iC).NE.0) THEN                                                    
        hdint = '#     b          t(b)'                                         
        hdcon = '#   Offset '                                                   
        hdvis = '#     q    '                                                   
        unt = 17                                                                
        CALL LINE(1,2,unt)                                                      
        write(unt,'(a8,i4,a14)') '#  model',model,'   RAW IMAGE  '              
        CALL LINE(1,1,unt)                                                      
        DO i = 1, nPok+2                                                        
          Elems(1,i) = bOut(i)                                                  
          Elems(2,i) = tauZout(i)                                               
          DO j = 1, NLambdaOut                                                  
c           check values:                                                       
            IF(IntOut(j,i).NE.IntOut(j,i).OR.IntOut(j,i).LT.1.e-20) THEN        
               IntOut(j,i) = 0.0                                                
            END IF                                                              
            Elems(j+2,i) = IntOut(j,i)                                          
c           we want intensity in Jy/arcsec2                                     
            Elems(j+2,i) = 7.83 * lambdaOut(j) * F1 * Elems(j+2,i)              
          END DO                                                                
        END DO                                                                  
        write(unt,'(a21,20f11.2)')hdint,(lambdaOut(j),j=1,NlambdaOut)           
        CALL MakeTable(Elems,nPok+2,NlambdaOut+2,unt)                           
      END IF                                                                    
      IF (iC.LT.0) THEN                                                         
c ---------  convolved images either add to .i## file or write in *.c## file -- 
        IF(iC.EQ.-3) unt = 21                                                   
        call line(1,2,unt)                                                      
        write(unt,'(a8,i4,a20)') '#  model',model,'   CONVOLVED IMAGE  '        
        call line(1,1,unt)                                                      
        DO i = 1, Nconv                                                         
         Elems(1,i) = Offset(i)                                                 
         DO j = 1, NLambdaOut                                                   
c        check values:                                                          
         IF(ConvInt(j,i).NE.ConvInt(j,i).OR.ConvInt(j,i).LT.1.e-20) THEN        
            ConvInt(j,i) = 0.0                                                  
         END IF                                                                 
         Elems(j+1,i) = ConvInt(j,i)                                            
         END DO                                                                 
        END DO                                                                  
        write(unt,'(a11,20f11.2)')hdcon,(lambdaOut(i),i=1,NlambdaOut)           
        CALL MakeTable(Elems,Nconv,NLambdaOut+1,unt)                            
        IF (psftype.LT.3.AND.model.EQ.1) THEN                                   
c         Wavelength dependent PSFs, print them separately in *.psf             
c         first generate wavelength dependent PSFs                              
          DO j = 1, NlambdaOut                                                  
            iLambda = j                                                         
            DO i = 1, Nconv                                                     
              PSFfunc(j,i) = PSFN(Offset(i))                                    
c             check dynamic range                                               
              CALL CHKRANGE(dynrange,PSFfunc(j,i))                              
            END DO                                                              
          END DO                                                                
c         print them out                                                        
          DO i = 1, Nconv                                                       
           write(23,'(1p,e12.5,20e10.3)')Offset(i),                             
     &                                 (PSFfunc(j,i),j=1,NlambdaOut)            
          END DO                                                                
        END IF                                                                  
      END IF                                                                    
c --------------     visibility curves to *.vnn file    ------------------------
      IF (iV.NE.0) THEN                                                         
        IF(abs(iC).EQ.3) unt = 22                                               
        call line(1,2,unt)                                                      
        write(unt,'(a8,i4,a14)') '#  model',model,'  VISIBILITY  '              
        call line(1,1,unt)                                                      
        DO i = 1, Nvisi                                                         
          Elems(1,i) = qtheta1(i)                                               
          DO j = 1, NLambdaOut                                                  
c          check for nonsense values:                                           
           IF(Visib(j,i).NE.Visib(j,i).OR.Visib(j,i).LT.1.e-20) THEN            
              Visib(j,i) = 0.0                                                  
           END IF                                                               
           Elems(j+1,i) = Visib(j,i)                                            
          END DO                                                                
        END DO                                                                  
        write(unt,'(a11,20f11.2)')hdvis,(lambdaOut(i),i=1,NlambdaOut)           
        CALL MakeTable(Elems,Nvisi,NLambdaOut+1,unt)                            
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
c This function is taken from Moshe Elitzur.           [Z.I., Nov. 1995]        
c =======================================================================       
      DOUBLE PRECISION FUNCTION RDINP(Equal,IUNIT)                              
c =======================================================================       
C     Read lines, up to 232 long, from pre-opened unit IUNIT and extract        
C     all input numbers from them. When EQUAL is set, numeric input data        
C     must be preceded by an equal sign.All non-numeric data and numbers        
C     not preceded by = when EQUAL is on are ignored.RDINP = next number        
C     encountered (after equal sign) and terminated by a nun-numeric            
C     symbol (termination with blank is best). Commas and exponential           
C     notation are allowed.  All text after % is ignored, as in TeX.            
C     Lines with * in the first column are echoed to the output device.         
C     The number is comprised of an actual VALUE, decimal part FRAC and         
C     (integer) exponent PWR.  It has a SIGN, and the exponential power         
C     has sign SIGNEX. Logical flag to the decimal part is set in               
C     DECIMAL. The search is conducted between FIRST, which is                  
C     continuously increased, and LAST.  A new line is read when FIRST          
C     exceeds LAST, and can be forced by calling with -IUNIT.  Actual           
C     extraction of numerical value is done in separate FUNCTION VAL.           
c =======================================================================       
      IMPLICIT None                                                             
      Integer IUnit,Ind,First, Last                                             
      DOUBLE PRECISION Value,Val,Frac,Pwr,Sign,Signex                           
      CHARACTER Card*(232),CR,prev,Term,Next                                    
      Logical Equal,digit,minus,sgn,dot,E,decimal                               
      Save Card,First,Last                                                      
      DATA First/1/, Last/0/                                                    
c -----------------------------------------------------------------------       
C FUNCTION STATEMENTS                                                           
      digit(CR) = CR.GE.'0' .AND. CR.LE.'9'                                     
      minus(CR) = CR.EQ.'-'                                                     
      sgn(CR)   = CR.EQ.'+' .OR. CR.EQ.'-'                                      
      dot(CR)   = CR.EQ.'.'                                                     
      E(CR)     = CR.EQ.'E' .OR. CR.EQ.'e'                                      
C                                                                               
      IF (IUnit.lt.0) Then                                                      
         First = Last + 1                                                       
         IUnit = -IUnit                                                         
      END IF                                                                    
c     Start the search for the next number:                                     
  1   RDINP  = 0.                                                               
      VALUE  = 0.                                                               
      FRAC   = 0.                                                               
      PWR    = 0.                                                               
      SIGN   = 1.                                                               
      SIGNEX = 1.                                                               
      Decimal = .False.                                                         
      If (first.gt.last) then                                                   
c     Time to get a new line                                                    
         READ (IUNIT, '(A)' , END = 99) Card                                    
         first = 1                                                              
         last = len(Card)                                                       
c        Find start of trailing junk:                                           
         DO WHILE (Card(last:last).LE.' ')                                      
          last = last - 1                                                       
          if (last.lt.first) goto 1                                             
         END DO	                                                                
         IF (Card(first:first).EQ.'*') WRITE (12,'(A)') Card(1:last)            
         ind = Index(Card,'%')                                                  
         if (ind.gt.0) last = ind - 1                                           
      End If                                                                    
                                                                                
c     Get past the next '=' when the EQUAL flag is set                          
      If (Equal) then                                                           
        DO WHILE (Card(first:first).NE.'=')                                     
          first = first + 1                                                     
          IF (first.GT.last) goto 1                                             
        END DO                                                                  
      End If                                                                    
c     OK, start searching for the next digit:                                   
      Do While (.not.digit(Card(first:first)))                       		         
	    first = first + 1                                                          
		if (first.gt.last) goto 1                                                     
	End Do                                                                         
c     Check if it is a negative or decimal number                               
      If (first.gt.1) then                                                      
         prev = Card(first-1:first-1)                                           
         if (minus(prev)) sign = -1.                                            
         if (dot(prev)) then                                                    
           decimal = .True.                                                     
           if (first.gt.2 .and. minus(Card(first-2:first-2))) sign = -1.        
         end if                                                                 
      End If                                                                    
c     Extract the numerical value                                               
      IF (.not.Decimal) Then                                                    
         Value = VAL(Card,first,last,decimal,Term)                              
c        Check for a possible decimal part.  Termination with '.E' or           
c        '.e' is equivalent to the same ending without '.'                      
         if (first.lt.last.and.dot(Term)) then                                  
            first = first + 1                                                   
            next = Card(first:first)                                            
            if (digit(next)) decimal = .true.                                   
            if (E(next)) Term = 'E'                                             
         end if                                                                 
      END IF                                                                    
c     Extract the decimal fraction, when it exists                              
      IF (Decimal) Frac = VAL(Card,first,last,decimal,Term)                     
c     An exponential may exist if any part terminated with 'E' or 'e'           
      IF (first.lt.last.and.E(term)) then                                       
         first = first + 1                                                      
         next = Card(first:first)                                               
         if (first.lt.last.and.sgn(next))then                                   
            first = first + 1                                                   
            if (minus(next)) Signex = -1.                                       
         end if                                                                 
         decimal = .False.                                                      
         PWR = VAL(Card,first,last,decimal,Term)                                
      END IF                                                                    
c     Finally, put the number together                                          
      RDINP = Sign*(Value + Frac)*10**(Signex*PWR)                              
      Return                                                                    
                                                                                
99    WRITE (12,'(3(1x,a,/))')                                                  
     *' ****TERMINATED. EOF reached by RDINP while looking for input. ',        
     *' *** Last line read:',Card                                               
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
                                                                                
c ***********************************************************************       
c This function is taken from Moshe Elitzur            [Z.I., Nov. 1995]        
c =======================================================================       
      DOUBLE PRECISION FUNCTION VAL(Card,first,last,decimal,Term)               
c     Extract numerical value from CARD, begining at position FIRST up          
c     to the first non-digit encountered.  The terminating character is         
c     returned in TERM, its position in FIRST. Commas are treated as            
c     simple separators.                                                        
c =======================================================================       
      IMPLICIT None                                                             
      Character Card*(*), Term, CH                                              
      Logical Decimal, digit                                                    
      Integer IVAL, first, last, first0                                         
      Real pwr                                                                  
c -----------------------------------------------------------------------       
C     FUNCTION STATEMENTS                                                       
      IVAL(CH) = ICHAR(CH) - ICHAR('0')                                         
      digit(CH) = CH.GE.'0' .AND. CH.LE.'9'                                     
c                                                                               
      VAL = 0.                                                                  
      pwr = 1.                                                                  
      first0 = first                                                            
      DO 10 first = first0, last                                                
         Term = Card(first:first)                                               
         if (Term.eq.',') goto 10                                               
         if (.not.digit(Term)) Return                                           
         If (decimal) then                                                      
            pwr = pwr*0.1                                                       
            Val = Val + pwr*IVAL(Term)                                          
         Else                                                                   
            Val = 10.*Val + IVAL(Term)                                          
         End If                                                                 
  10  CONTINUE                                                                  
      Term = ' '                                                                
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
c =========================================================================     
c     This is are the core subroutines for radiative transfer calculation.      
c                                                             [MN, Mar'99]      
c =========================================================================     
C     Table of Contents                                                         
C                                                                               
C     EMISSION                                                                  
C     GFUN                                                                      
C     FINDTEMP                                                                  
C     GRAYBODY                                                                  
C     INITTEMP                                                                  
C     INSERT                                                                    
C     INTETA                                                                    
C     INVERT                                                                    
C     KINT4                                                                     
C     MATRIX                                                                    
C     NORDLUND                                                                  
C     OCCLTMSG                                                                  
C     PGRID                                                                     
C     PLANCK                                                                    
C     RADTRANSF                                                                 
C     SETGRIDS                                                                  
C     SOLVE                                                                     
C     STAR                                                                      
C     TRAPZD2                                                                   
C     TWOFUN                                                                    
C     WEIGHTS                                                                   
C     YGRID                                                                     
c =========================================================================     
                                                                                
c ***********************************************************************       
      SUBROUTINE Emission(flag,nG,Td,alpha,abund,Uin,Emiss)                     
c =======================================================================       
c This subroutine calculates emission term from the temperature, abund          
c and alpha arrays for flag=0, and adds U to it for flag=1.                     
c                                                      [Z.I., Mar. 1996]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      INTEGER iG, iY, iL, nG, flag                                              
      DOUBLE PRECISION Uin(npL,npY), Td(npG,npY), alpha(npG,npY),               
     &       Emiss(npL,npY), abund(npG,npY), EmiG, xP, Planck                   
c -----------------------------------------------------------------------       
c     first initialize Emiss                                                    
c     loop over wavelengths                                                     
      DO iL = 1, nL                                                             
c       loop over radial coordinate                                             
        DO iY = 1, nY                                                           
          Emiss(iL,iY) = 0.0                                                    
        END DO                                                                  
      END DO                                                                    
c     calculate emission term for each component and add it to Emiss            
c     loop over wavelengths                                                     
      DO iL = 1, nL                                                             
c       loop over radial coordinate                                             
        DO iY = 1, nY                                                           
c         loop over grains                                                      
          DO iG = 1, nG                                                         
            xP = 14400.0 / lambda(iL) / Td(iG,iY)                               
            EmiG = abund(iG,iY) * alpha(iG,iY) * Planck(xP)                     
c           add contribution for current grains                                 
            Emiss(iL,iY) = Emiss(iL,iY) + EmiG                                  
          END DO                                                                
c         if needed add Uin                                                     
          IF (flag.EQ.1) THEN                                                   
            Emiss(iL,iY) = Emiss(iL,iY) + Uin(iL,iY)                            
          END IF                                                                
          IF (Emiss(iL,iY).LT.dynrange*dynrange) Emiss(iL,iY) = 0.0             
        END DO                                                                  
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      DOUBLE PRECISION FUNCTION Gfun(Tt)                                        
c =======================================================================       
c This is an auxiliary function used in determining the dust temperature        
c                                                      [Z.I., Nov. 1995]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER iW, nW                                                            
      DOUBLE PRECISION Tt, wav(npL), Eff(npL), ff(npL), xP, gg, f1, f2,         
     &       Planck                                                             
      COMMON /gfunction/ wav, Eff, f1, f2, nW                                   
c -----------------------------------------------------------------------       
      DO iW = 1, nW                                                             
        xP = 14400.0 / wav(iW) / Tt                                             
        ff(iW) = Eff(iW) * Planck(xP) / wav(iW)                                 
      END DO                                                                    
      CALL Simpson(npL,1,nW,wav,ff,gg)                                          
      Gfun = f1 - f2 * gg * Tt**4.0                                             
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE FindTemp(geom,Utot,nG,Td,alpha)                                
c =======================================================================       
c This subroutine finds new temperature and alpha arrays from Utot.             
c Temperature is obtained by solving:                                           
c   f1(y) - f2(y)*g(T) = 0                                                      
c where                                                                         
c   f1(y) = Int(Qabs*Utot*dlambda)                                              
c   f2(y) = alpha(y=1)*y**2/Tsub**4                                             
c    g(T) = Td**4*Int(Qabs*Planck(Td)*dlambda)                                  
c alpha(y)= Int(Qabs*Utot*dlambda)/Int(Qabs*Planck(Td)*dlambda)                 
c The equation is solved by Ridder's method (subroutine from Num.Rec)           
c                                                      [Z.I., Oct. 1996]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER iG, iY, iL,nG, nW, succes, geom                                   
      DOUBLE PRECISION Utot(npL,npY), Td(npG,npY), alpha(npG,npY),              
     &       fnum(npL), fdenum(npL), num, denum, xP, Tin, aa, Planck,           
     &       wav(npL), f1, f2, gfun, x1, x2,  Eff(npL), Zriddr                  
      EXTERNAL gfun                                                             
      COMMON /gfunction/ wav, Eff, f1, f2, nW                                   
c -----------------------------------------------------------------------       
      nW = nL                                                                   
c     loop over grains                                                          
      DO iG = 1, nG                                                             
c       evaluate alpha(1)                                                       
        DO iL = 1, nL                                                           
          fnum(iL) = SigmaA(iG,iL) * Utot(iL,1) / lambda(iL)                    
          xP = 14400.0 / lambda(iL) / Tsub(iG)                                  
          fdenum(iL) = SigmaA(iG,iL) * Planck(xP) / lambda(iL)                  
          wav(iL) = lambda(iL)                                                  
          Eff(iL) = SigmaA(iG,iL)                                               
        END DO                                                                  
        CALL Simpson(npL,1,nL,lambda,fnum,num)                                  
        CALL Simpson(npL,1,nL,lambda,fdenum,denum)                              
        alpha(iG,1) = num /  denum                                              
        aa = alpha(iG,1) / Tsub(iG)**4.0                                        
        Td(iG,1) = Tsub(iG)                                                     
c       loop over radial positions (solving f1-f2*g(T)=0)                       
        DO iY = 2, nY                                                           
c         calculate f1 and f2                                                   
          DO iL = 1, nL                                                         
            fnum(iL) = SigmaA(iG,iL) * Utot(iL,iY) / lambda(iL)                 
          END DO                                                                
          CALL Simpson(npL,1,nL,lambda,fnum,f1)                                 
          IF (geom.EQ.0) THEN                                                   
c           This is for slab geometry ...                                       
            f2 = aa                                                             
          ELSE                                                                  
c          ... and for spherical                                                
             f2 = aa * Y(iY)**2.0                                               
          END IF                                                                
c         estimate range for temperature                                        
          x1 = Td(iG,iY)                                                        
          x2 = 1.01*Td(iG,iY)                                                   
c         make sure that correct solution is bracketted                         
          CALL Zbrac(gfun,x1,x2,100,succes)                                     
          IF (succes.NE.1) THEN                                                 
            CALL Msg(4)                                                         
            ELSE                                                                
c           save the old value of Td                                            
            Tin = Td(iG,iY)                                                     
c           get new temperature                                                 
            Td(iG,iY) = Zriddr(gfun,x1,x2,500,accConv)                          
c           calculate alpha                                                     
            DO iL = 1, nL                                                       
              xP = 14400.0 / lambda(iL) / Td(iG,iY)                             
              fdenum(iL) = SigmaA(iG,iL) * Planck(xP) / lambda(iL)              
            END DO                                                              
            CALL Simpson(npL,1,nL,lambda,fdenum,denum)                          
            alpha(iG,iY) = f1 /  denum                                          
          END IF                                                                
        END DO                                                                  
      END DO                                                                    
c -----------------------------------------------------------------------       
999   RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE GrayBody(albedo,TAUgbTot,Ugb,fgb)                              
c =======================================================================       
c This subroutine solves the gray body problem for albedo=1 (or                 
c equivalently pure scattering) and scattering with absorption (but no          
c emission) for albedo<1, in a spherically symmetric envelope. Total            
c optical depth is TAUtot, and density law is specified elsewhere.              
c This subroutine was designed to be a part of Dusty and to use already         
c existing subroutines as much as possible, so some parts might seem to         
c be a little awkward.                           [ZI,Jul'96;MN,Sep'97]          
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER nYok, nPok                                                        
      DOUBLE PRECISION                                                          
     &       Ude(npL,npY), Uds(npL,npY), Us(npL,npY), fde(npL,npY),             
     &       fds(npL,npY), fs(npL,npY), Utot(npL,npY), ftot(npL,npY),           
     &       Td(npG,npY), Ubol(npY), fbol(npY), Uchck(npL,npY),                 
     &       UbolChck(npY), Spectrum(npL), SpecChar(30,99), tauF(npY),          
     &       Intens(npL,npP+2), IntOut(20,npP+2), SmC(30,99),                   
     &       Yok(npY), Pok(npP), tauOut(npL), tauZout(npP+2), F1,               
     &       tr(npY), rg(npG,npY), fsL(npL,npY), fsR(npL,npY), Eps(npY)         
      COMMON /solution/ Ude, Uds, Us, fde, fds, fs, Utot, ftot, Td,             
     &       Ubol, fbol, Uchck, UbolChck, Spectrum, SpecChar, tauF,             
     &       Intens, IntOut, SmC, Yok, Pok, tauOut, tauZout, F1,                
     &       tr, rg, fsL, fsR, Eps, nYok, nPok                                  
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      INTEGER iPGB, iY, nLstore, error                                          
      DOUBLE PRECISION ETAzp(npP,npY), mat0(npL,npY,npY),                       
     &       mat1(npL,npY,npY), Em(npL,npY), Ugb(npY), fgb(npY),                
     &       Dummy1(npL,npP,npY), Dummy2(npL,npP,npY), Dummy3(npL,npY),         
     &       albedo, pGB, TAUgbTot, TAUstore                                    
c ----------------------------------------------------------------------        
c     Values needed in this subroutine only                                     
      pGB = 0.0                                                                 
      iPGB = 0                                                                  
      nLstore = nL                                                              
      nL = 1                                                                    
      TAUstore = TAUtot(1)                                                      
      TAUtot(1) = TAUgbTot                                                      
c     generate spline coefficients for ETA                                      
      CALL setupETA                                                             
c     evaluate ETAzp                                                            
      CALL getETAzp(ETAzp)                                                      
c     generate some temporary arrays                                            
      DO iY = 1, nY                                                             
        Us(1,iY) = dexp(-ETAzp(1,iY)*TAUgbTot)                                  
        fs(1,iY) = Us(1,iY)                                                     
        Em(1,iY) = 0.0                                                          
        fde(1,iY) = 0.0                                                         
        omega(1,iY) = albedo                                                    
      END DO                                                                    
c     find radiative transfer matrices                                          
      CALL Matrix(ETAzp,pGB,iPGB,mat0,mat1,Dummy1,Dummy2)                       
c     solve for Utot                                                            
      CALL Invert(nY,nL,mat0,Us,Em,omega,Utot,Dummy3,error)                     
c     calculate flux, ftot                                                      
      CALL Multiply(1,npY,nY,npL,nL,mat1,Utot,omega,1,fs,ftot)                  
c     store to the output arrays                                                
      DO iY = 1, nY                                                             
        Ugb(iY) = Utot(1,iY)                                                    
        fgb(iY) = ftot(1,iY)                                                    
      END DO                                                                    
      nL = nLstore                                                              
      TAUtot(1) = TAUstore                                                      
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE InitTemp(ETAzp,nG,alpha)                                       
c =======================================================================       
c This subroutine initializes the first approximations for temperature          
c and alpha arrays. It uses approximations given by eqs. B5 and B7 from         
c IE96. Alpha(y) is defined as the ratio of Qabs averaged with the total        
c energy density at y and Qabs Planck averaged with T(y). That is, Psi          
c is alpha(y=1) and (T/Tsub)^4 = alpha(y) / alpha(1).  [Z.I., Jul. 1996]        
c                                                                               
c *** This version works for single grains ONLY ***                             
c                                                                               
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER nYok, nPok                                                        
      DOUBLE PRECISION                                                          
     &       Ude(npL,npY), Uds(npL,npY), Us(npL,npY), fde(npL,npY),             
     &       fds(npL,npY), fs(npL,npY), Utot(npL,npY), ftot(npL,npY),           
     &       Td(npG,npY), Ubol(npY), fbol(npY), Uchck(npL,npY),                 
     &       UbolChck(npY), Spectrum(npL), SpecChar(30,99), tauF(npY),          
     &       Intens(npL,npP+2), IntOut(20,npP+2), SmC(30,99),                   
     &       Yok(npY), Pok(npP), tauOut(npL), tauZout(npP+2), F1,               
     &       tr(npY), rg(npG,npY), fsL(npL,npY), fsR(npL,npY), Eps(npY)         
      COMMON /solution/ Ude, Uds, Us, fde, fds, fs, Utot, ftot, Td,             
     &       Ubol, fbol, Uchck, UbolChck, Spectrum, SpecChar, tauF,             
     &       Intens, IntOut, SmC, Yok, Pok, tauOut, tauZout, F1,                
     &       tr, rg, fsL, fsR, Eps, nYok, nPok                                  
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      DOUBLE PRECISION TAUslb(npL,npY),fsbol(npY),fpbol(npY),fmbol(npY),        
     &                 AveDev, RMS, maxFerr, fmax, fmed                         
      COMMON /slab/ TAUslb, fsbol, fpbol, fmbol, AveDev, RMS, maxFerr,          
     &                 fmax, fmed                                               
      INTEGER iL, iY, iG, nG, istop, iW                                         
      DOUBLE PRECISION ETAzp(npP,npY), xP, Planck, resaux, aux(npY),            
     &       faux(npL), fstar(npY), Qfstar(npY),  fstarQ(npY),                  
     &       alpha(npG,npY), oldalpha(npY), Sigext(npL), QP(npY),               
     &       QF(npY),EtaQF, ETA, delta                                          
c -----------------------------------------------------------------------       
c     Lambda integral of fs=f_e*exp(-TAU) -> fstar                              
      DO iY = 1, nY                                                             
c       generate auxiliary function for lambda integration:                     
        DO iL = 1, nL                                                           
          faux(iL) = fs(iL,iY) / lambda (iL)                                    
        END DO                                                                  
        CALL Simpson(npL,1,nL,lambda,faux,resaux)                               
        fstar(iY) = resaux                                                      
      END DO                                                                    
c     loop over grains                                                          
      DO iG = 1, nG                                                             
c       first approximation for temperature                                     
        DO iY = 1, nY                                                           
          Td(iG,iY) = Tsub(iG) / Y(iY)**0.4                                     
          alpha(iG,iY) = 1.0                                                    
        END DO                                                                  
c       generate the extinction cross-section Sigext                            
        DO iL = 1, nL                                                           
          Sigext(iL) = SigmaA(iG,iL)+SigmaS(iG,iL)                              
        END DO                                                                  
c       Lambda integral of Sigext*fs -> Qfstar in IE97                          
        DO iY = 1, nY                                                           
c         generate auxiliary function for lambda integration:                   
          DO iL = 1, nL                                                         
            faux(iL) = Sigext(iL) * fs(iL,iY) / lambda (iL)                     
          END DO                                                                
          CALL Simpson(npL,1,nL,lambda,faux,resaux)                             
          Qfstar(iY) = resaux                                                   
        END DO                                                                  
c       iterate until Td and Psi converge (i.e. until alpha does)               
        istop = 0                                                               
        DO WHILE (istop.NE.1)                                                   
          istop = 1                                                             
c         find Planck average of Sigext and calculate QF (eq. B5)               
          DO iY = 1, nY                                                         
c           generate auxiliary function for lambda integration:                 
            DO iL = 1, nL                                                       
              xP = 14400.0 / Td(iG,iY) / lambda(iL)                             
              faux(iL) = Sigext(iL) * Planck(xP) / lambda (iL)                  
            END DO                                                              
            CALL Simpson(npL,1,nL,lambda,faux,resaux)                           
c           Planck average of Sigext                                            
            QP(iY) = resaux                                                     
c           calculate intetgral next to 1/y2 in eq. B7 -> fstarQ                
            DO iL = 1, nL                                                       
              faux(iL) = fs(iL,iY)*(1.-SigmaA(iG,iL)/QP(iY))/lambda(iL)         
            END DO                                                              
            CALL Simpson(npL,1,nL,lambda,faux,resaux)                           
            fstarQ(iY) = resaux                                                 
c           calculate QF (eq. B5)                                               
            QF(iY) = Qfstar(iY) + QP(iY) * (1.0-fstar(iY))                      
c           store current alpha                                                 
            oldalpha(iY) = alpha(iG,iY)                                         
          END DO                                                                
c         for each y calculate new alpha from eq. B7                            
c         two y loops are needed because of integral QF*ETA                     
          DO iY = 1, nY                                                         
c           first evaluate integral of QF with ETA...                           
            DO iW = iY, nY                                                      
              aux(iW) = QF(iW) * ETA(Y(iW)) / Y(iW)**2.0                        
            END DO                                                              
            CALL Simpson(npY,iY,nY,Y,aux,EtaQF)                                 
c           ... and alpha                                                       
            alpha(iG,iY) = 3.*TAUfid*EtaQF + (1-fstarQ(iY))                     
c           calculate temperature                                               
            Td(iG,iY) = Tsub(iG) * (alpha(iG,iY)/alpha(iG,1))**0.25             
            Td(iG,iY) = Td(iG,iY) / dsqrt(Y(iY))                                
            delta = DABS((alpha(iG,iY)-oldalpha(iY))/alpha(iG,iY))              
            IF (delta.GT.accConv) istop = 0                                     
          END DO                                                                
        END DO                                                                  
c       end of iterations                                                       
      END DO                                                                    
c     end of loop over grains (iG)                                              
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE INSERT(pstar,iP,iPstar)                                        
c =======================================================================       
c This subroutine inserts two rays if the stellar disk is finite. The           
c first ray, corresponding to the disk edge, has index iPstar, and the          
c following ray with 1.0001 greater impact parameter has index iPstar+1.        
c The only exception is if pstar>0.9999 when only a single ray is               
c inserted.                                            [Z.I., Feb. 1996]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER iP, iPstar                                                        
      DOUBLE PRECISION pstar                                                    
c -----------------------------------------------------------------------       
      IF (pstar.GE.0.999999) THEN                                               
        pstar = 0.999999                                                        
        iP = nPcav + 1                                                          
        iPstar = nPcav + 1                                                      
        P(iP) = pstar                                                           
        ELSE                                                                    
        IF (pstar.GE.P(nPcav)) THEN                                             
          IF (pstar.EQ.P(nPcav)) pstar = 1.00001*P(nPcav)                       
          P(nPcav+1) = pstar                                                    
          iPstar = nPcav + 1                                                    
          P(nPcav+2) = 1.00001*pstar                                            
          iP = nPcav + 2                                                        
        ELSE                                                                    
          iPstar = 0                                                            
          iP = 0                                                                
          DO WHILE (iPstar.EQ.0)                                                
            iP = iP + 1                                                         
            IF (P(iP).GT.pstar.AND.iPstar.EQ.0) iPstar = iP                     
          END DO                                                                
          DO iP = 1, nPcav-iPstar+1                                             
            P(nPcav+3-iP) = P(nPcav+1-iP)                                       
          END DO                                                                
          IF (pstar.EQ.P(iPstar-1)) pstar = 1.00001*P(iPstar-1)                 
          P(iPstar) = pstar                                                     
          P(iPstar+1) = 1.00001*pstar                                           
          iP = nPcav + 1                                                        
        END IF                                                                  
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      DOUBLE PRECISION FUNCTION IntETA(p,iW1,w1,w)                              
c =======================================================================       
c This function calculates the integral over the normalized dens. prof.         
c along the line of sight with impact parameter p and between the points        
c corresponding to y=w1 and y=w. The method used is spline approximation        
c for normalized density distribution ETA and subsequent integration            
c performed analytically by MAPLE (these results are given through              
c soubroutine Maple3).                         [ZI,Feb'96,MN,Aug'97]            
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER iW1, iC                                                           
      DOUBLE PRECISION  p, w1, w, aux(4), z, z1, aux1(4)                        
c -----------------------------------------------------------------------       
      z = dsqrt(w*w-p*p)                                                        
      z1 = dsqrt(w1*w1-p*p)                                                     
c     integrals calculated by MAPLE                                             
      CALL Maple3(w,z,p,aux)                                                    
      CALL Maple3(w1,z1,p,aux1)                                                 
      DO iC = 1, 4                                                              
        aux(iC) = aux(iC) - aux1(iC)                                            
      END DO                                                                    
      IntETA = 0.0                                                              
      DO iC = 1, 4                                                              
        IntETA = IntETA + ETAcoef(iW1,iC) * aux(iC)                             
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE INVERT(nY,nL,mat,Us,Em,omega,Utot,Uold,error)                  
c =======================================================================       
c This subroutine solves the linear system                                      
c [Utot] = [Us+Em] + [mat0]*[omega*Utot] by calling LINSYS subroutine.          
c                                                      [Z.I., Nov. 1995]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER nY, iL,nL,iY, iYaux, Kronecker, error                             
      DOUBLE PRECISION mat(npL,npY,npY), Us(npL,npY), Em(npL,npY),              
     &       Utot(npL,npY), Uold(npL,npY),B(npY), A(npY,npY), X(npY),           
     &       omega(npL,npY)                                                     
c -----------------------------------------------------------------------       
      error = 0                                                                 
c     first copy Utot to Uold                                                   
      DO iL = 1, nL                                                             
        DO iY = 1, nY                                                           
          Uold(iL,iY) = Utot(iL,iY)                                             
        END DO                                                                  
      END DO                                                                    
c     calculate new energy density                                              
c     loop over wavelengths                                                     
      DO iL = 1, nL                                                             
c       generate the vector of free coefficients, B, and matrix A               
        DO iY = 1, nY                                                           
          B(iY) = Us(iL,iY)                                                     
          DO iYaux = 1, nY                                                      
            Kronecker = 0                                                       
            IF (iY.EQ.iYaux) Kronecker = 1                                      
            B(iY) = B(iY) +                                                     
     &                (1.-omega(iL,iYaux))*Em(iL,iYaux)*mat(iL,iY,iYaux)        
            A(iY,iYaux) = Kronecker - omega(iL,iYaux) * mat(iL,iY,iYaux)        
          END DO                                                                
        END DO                                                                  
                                                                                
c       solve the system                                                        
        CALL LINSYS(nY,A,B,X,error)                                             
        IF(error.NE.0) THEN                                                     
         CALL MSG(20)                                                           
         iERROR = iERROR + 1                                                    
         RETURN                                                                 
        END IF                                                                  
c       store the result                                                        
        DO iY = 1, nY                                                           
          IF (X(iY).GE.dynrange*dynrange) THEN                                  
            Utot(iL,iY) = X(iY)                                                 
          ELSE                                                                  
            Utot(iL,iY) = 0.0                                                   
          END IF                                                                
        END DO                                                                  
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE Kint4(TAUaux,iP,iL,nZ,K1p,K2p,K3p,K4p,K1m,K2m,K3m,K4m)         
c =======================================================================       
c For given wavelength (iL) and impact parameter (iP), this subroutine          
c calculates integrals Knp and Knm defined as:                                  
c               Knp(iZ)=INT[PHIn(tz)*exp(tz)*dtz]                               
c from tz1=TAUaux(iL,iP,iZ) to tz2=TAUaux(iL,iP,iZ+1) and analogously for       
c Km with exp(tz) replaced by exp(-tz). Function PHIn is defined as             
c x**(n-1)/Yloc^2, where Yloc is the local radius corresponding to tz,          
c and x measures relative radial tau: x = (rt - tL)/(tR-tL). Here rt is         
c the radial optical depth corresponding to tz and tL and tR are radial         
c optical depths at the boundaries of the integration interval:                 
c tL = TAUaux(iL,1,iZ) = rt(iZ) and tR = TAUaux(iL,1,iZ+1) = rt(iZ+1).          
c Integration is performed in z space by Romberg integration implemented        
c in subroutine ROMBERG2 (slightly changed version of 'qromb' from Num.         
c Recipes).                                     [ZI,Feb'96;MN,Sep'97]           
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER iYfirst, YPequal, Plast                                           
      DIMENSION iYfirst(npP), YPequal(npP), Plast(npY)                          
      COMMON /Yfirst/ iYfirst, YPequal, Plast                                   
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER iP, iL, nZ, iZ, iW1, iLaux, iC                                    
      DOUBLE PRECISION TAUaux(npL,npP,npY), K1p(npY),K2p(npY),K3p(npY),         
     &    K4p(npY), K1m(npY), K2m(npY), K3m(npY), K4m(npY), result(8),          
     &    Kaux(8), deltrton(4),tRL,paux, w1,w2,wL, delTAUzp, z1, z2             
      COMMON /phi2/ paux, w1, wL, delTAUzp, iLaux, iW1                          
c -----------------------------------------------------------------------       
      paux = P(iP)                                                              
      iLaux = iL                                                                
c     iLaux is needed to avoid compiler errors since it is in COMMON            
c     /phi2/ (here and in 'TWOfun'), while iL is transferred as a               
c     argument loop over positions on the line of sight                         
      DO iZ = 1, nZ-1                                                           
c       index for the local radial position (left boundary)                     
        iW1 = iYfirst(iP) + iZ - 1                                              
c       radii at the boundaries                                                 
        wL = Y(iW1)                                                             
        IF (iZ.EQ.1) THEN                                                       
          if (paux.GT.1.0) then                                                 
            w1 = paux                                                           
            else                                                                
            w1 = 1.0                                                            
          end if                                                                
          ELSE                                                                  
          w1 = Y(iW1)                                                           
        END IF                                                                  
        w2 = Y(iW1+1)                                                           
        z1 = dsqrt(DABS(w1*w1 - paux*paux))                                     
        z2 = dsqrt(w2*w2 - paux*paux)                                           
c       radial tau-difference at the bound., scaled to tot. opt. depth          
        tRL = TAUaux(iL,1,iW1+1)-TAUaux(iL,1,iW1)                               
c       auxiliary quantity aux/tRL**(n-1)                                       
        deltrton(1) = TAUtot(iL)                                                
        DO iC = 1, 3                                                            
          deltrton(iC+1) = deltrton(iC)/tRL                                     
        END DO                                                                  
c       delTAUzp is needed in PHIn fun's                                        
        delTAUzp = TAUaux(iL,iP,iZ+1)-TAUaux(iL,iP,iZ)                          
c       integrate this step for all 8 cases                                     
        CALL ROMBERG2(z1,z2,result)                                             
c       generate output values                                                  
        DO iC = 1, 4                                                            
           Kaux(iC) = result(iC) * deltrton(iC)                                 
           Kaux(iC+4) = result(iC+4) * deltrton(iC)                             
        END DO                                                                  
        K1m(iZ) = Kaux(1)                                                       
        K2m(iZ) = Kaux(2)                                                       
        K3m(iZ) = Kaux(3)                                                       
        K4m(iZ) = Kaux(4)                                                       
        K1p(iZ+1) = Kaux(5)                                                     
        K2p(iZ+1) = Kaux(6)                                                     
        K3p(iZ+1) = Kaux(7)                                                     
        K4p(iZ+1) = Kaux(8)                                                     
      END DO                                                                    
c     set undefined elements to 0                                               
      K1m(nZ) = 0.0                                                             
      K2m(nZ) = 0.0                                                             
      K3m(nZ) = 0.0                                                             
      K4m(nZ) = 0.0                                                             
      K1p(1) = 0.0                                                              
      K2p(1) = 0.0                                                              
      K3p(1) = 0.0                                                              
      K4p(1) = 0.0                                                              
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE matrix(ETAzp,pstar,iPstar,m0,m1,mifront,miback)                
c =======================================================================       
c This subroutine evaluates radiative transfer matrix for spherically           
c symmetric envelope. Here m is the order of moment (0 for energy dens.,        
c 1 for flux, 2 for pressure etc.), ETAzp is array of optical depths            
c along the line of sight and mat is radiative transfer matrix.                 
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER m                                                                 
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER iYfirst, YPequal, Plast                                           
      DIMENSION iYfirst(npP), YPequal(npP), Plast(npY)                          
      COMMON /Yfirst/ iYfirst, YPequal, Plast                                   
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER nYok, nPok                                                        
      DOUBLE PRECISION                                                          
     &       Ude(npL,npY), Uds(npL,npY), Us(npL,npY), fde(npL,npY),             
     &       fds(npL,npY), fs(npL,npY), Utot(npL,npY), ftot(npL,npY),           
     &       Td(npG,npY), Ubol(npY), fbol(npY), Uchck(npL,npY),                 
     &       UbolChck(npY), Spectrum(npL), SpecChar(30,99), tauF(npY),          
     &       Intens(npL,npP+2), IntOut(20,npP+2), SmC(30,99),                   
     &       Yok(npY), Pok(npP), tauOut(npL), tauZout(npP+2), F1,               
     &       tr(npY), rg(npG,npY), fsL(npL,npY), fsR(npL,npY), Eps(npY)         
      COMMON /solution/ Ude, Uds, Us, fde, fds, fs, Utot, ftot, Td,             
     &       Ubol, fbol, Uchck, UbolChck, Spectrum, SpecChar, tauF,             
     &       Intens, IntOut, SmC, Yok, Pok, tauOut, tauZout, F1,                
     &       tr, rg, fsL, fsR, Eps, nYok, nPok                                  
      INTEGER iP, iL, nZ(npP), iZ, iW,  iY, flag, im, iPstar, jZ, error         
      DOUBLE PRECISION ETAzp(npP,npY),TAUaux(npL,npP,npY), haux(npP),H,         
     &       m0(npL,npY,npY), m1(npL,npY,npY), pstar, addplus,addminus,         
     &       Tplus(npP,npY,npY),Tminus(npP,npY,npY), xN(npP),yN(npP),           
     &       result1,result2,resaux,TAUr(npY), wm(npY),wmT(npY),wp(npY),        
     &       alpha(npY,npY),beta(npY,npY),gamma(npY,npY),delta(npY,npY),        
     &       wgmatp(npY,npY), wgmatm(npY,npY), mifront(npL,npP,npY),            
     &       miback(npL,npP,npY), fact, faux                                    
c -----------------------------------------------------------------------       
      error = 0                                                                 
c     generate auxiliary arrays haux & nZ                                       
      DO iP = 1, nP                                                             
c       parameter alowing for a line of sight terminating on the star           
c       H(x1,x2) is the step function.                                          
        haux(iP) = H(P(iP),pstar)                                               
c       upper limit for the counter of z position                               
        nZ(iP) = nY + 1 - iYfirst(iP)                                           
      END DO                                                                    
c     Using the local array TAUaux to avoid multiple calculations of the        
c     product                                                                   
      DO iL = 1, nL                                                             
       DO iP = 1, nP                                                            
        DO iY = 1, nY                                                           
           TAUaux(iL,iP,iY) = ETAzp(iP,iY)*TAUtot(iL)                           
        END DO                                                                  
       END DO                                                                   
      END DO                                                                    
c     -- evaluate matrix elements --                                            
c     loop over wavelengths                                                     
      DO iL = 1, nL                                                             
c       radial optical depths                                                   
        DO iY = 1, nY                                                           
          TAUr(iY) = ETAzp(1,iY)*TAUtot(iL)                                     
        END DO                                                                  
c       auxiliary arrays for given TAUr                                         
        CALL MYSPLINE(TAUr,nY,alpha,beta,gamma,delta)                           
c       loop over impact parameters                                             
        DO iP = 1, nP-1                                                         
c         set T-s to 0                                                          
          DO iY = 1, nY                                                         
            DO iW = 1, nY                                                       
              Tplus(iP,iY,iW) = 0.0                                             
              Tminus(iP,iY,iW) = 0.0                                            
            END DO                                                              
          END DO                                                                
c         generate weights matrices                                             
          CALL WEIGHTS(TAUaux,iP,iL,nZ(iP),nY,alpha,beta,gamma,delta,           
     &                 wgmatp,wgmatm)                                           
c         first position on the line of sight                                   
          IF (YPequal(iP).EQ.1) THEN                                            
            iY = iYfirst(iP)                                                    
            DO iW = 1, nY                                                       
c             cummulative weights for parts II & III                            
              wmT(iW) = 0.0                                                     
              DO jZ = 1, nZ(iP)-1                                               
                fact = dexp(-TAUaux(iL,iP,jZ))                                  
                wmT(iW) = wmT(iW) + fact * wgmatm(jZ,iW)                        
              END DO                                                            
              Tplus(iP,iY,iW) = Tplus(iP,iY,iW) + haux(iP)*wmT(iW)              
              Tminus(iP,iY,iW) = Tminus(iP,iY,iW) + wmT(iW)                     
            END DO                                                              
          END IF                                                                
c         loop over positions on the line of sight                              
          DO iZ = 2, nZ(iP)                                                     
c           increase index for radial position                                  
            iY = iYfirst(iP) + iZ - 1                                           
c           generate weights for this position                                  
            DO iW = 1, nY                                                       
              wp(iW) = 0.0                                                      
              wm(iW) = 0.0                                                      
              wmT(iW) = 0.0                                                     
c             part I                                                            
              DO jZ = 2, iZ                                                     
                fact = dexp(TAUaux(iL,iP,jZ)-TAUaux(iL,iP,iZ))                  
                wp(iW) = wp(iW) + fact * wgmatp(jZ,iW)                          
              END DO                                                            
c             part II & III                                                     
              DO jZ = 1, nZ(iP)-1                                               
                fact = dexp(-(TAUaux(iL,iP,iZ)+TAUaux(iL,iP,jZ)))               
                wmT(iW) = wmT(iW) + fact * wgmatm(jZ,iW)                        
              END DO                                                            
c             part IV                                                           
              IF (iZ.LT.nZ(iP)) THEN                                            
                DO jZ = iZ, nZ(iP)-1                                            
                  fact = dexp(-(TAUaux(iL,iP,jZ)-TAUaux(iL,iP,iZ)))             
                  wm(iW) = wm(iW) + fact * wgmatm(jZ,iW)                        
                END DO                                                          
              ELSE                                                              
                wm(iW) = 0.0                                                    
              END IF                                                            
c             add contribution from this step                                   
              addplus = wp(iW) + haux(iP)*wmT(iW)                               
              Tplus(iP,iY,iW) = Tplus(iP,iY,iW) + addplus                       
              addminus = wm(iW)                                                 
              Tminus(iP,iY,iW) = Tminus(iP,iY,iW) + addminus                    
            END DO                                                              
c         end of loop over iZ                                                   
          END DO                                                                
c       end of the impact parameter loop, iP                                    
        END DO                                                                  
c       add points on the edge                                                  
        DO iW = 1, nY                                                           
          Tplus(nP,nY,iW) = 0.0                                                 
          Tminus(nP,nY,iW) = 0.0                                                
        END DO                                                                  
c     ============================                                              
c       find mat(iL,iY,iW) -> angular (mu) integration                          
c       loop over moments (without calculation of rad.pressure)                 
        DO im = 1, 2                                                            
          m = im - 1                                                            
c         loop over radial positions                                            
          DO iY = 1, nY                                                         
c           generate mu arrray                                                  
            DO iP = 1, Plast(iY)                                                
              xN(iP) = sqrt(1.0-(P(iP)/Y(iY)*P(iP)/Y(iY)))                      
            END DO                                                              
c           loop over local (radial) positions                                  
            DO iW = 1, nY                                                       
c             generate intensity array for NORDLUND                             
              DO iP = 1, Plast(iY)                                              
c               'faux' is a representation of (-1)**m                           
                 faux = 1.0 - 2.0*MOD(m,2)                                      
                 yN(iP) = Tplus(iP,iY,iW) + faux*Tminus(iP,iY,iW)               
c               store matrix elements to evaluate intensity (*1/4Pi)            
                IF (im.EQ.1.AND.iY.EQ.nY) THEN                                  
                  mifront(iL,iP,iW) = 0.0795775 * Tplus(iP,iY,iW)               
                  miback(iL,iP,iW) = 0.0795775 * Tminus(iP,iY,iW)               
                END IF                                                          
              END DO                                                            
c             angular integration inside cavity                                 
              IF (pstar.GT.0.0) THEN                                            
                CALL NORDLUND(0,xN,yN,1,iPstar,m,resaux,error)                  
                IF (error.NE.0) GOTO 999                                        
                IF (nPcav.GT.iPstar)                                            
     &          CALL NORDLUND(0,xN,yN,iPstar+1,nPcav+1,m,result1,error)         
                IF (error.NE.0) GOTO 999                                        
                result1 = result1 + resaux                                      
              ELSE                                                              
                CALL NORDLUND(0,xN,yN,1,nPcav+1,m,result1,error)                
                IF (error.NE.0) GOTO 999                                        
              END IF                                                            
c             flag for analytic integration outside cavity                      
              IF (iY.GT.6) THEN                                                 
                flag = 1                                                        
              ELSE                                                              
                flag = 0                                                        
              ENDIF                                                             
c             angular integration outside cavity                                
              IF (iY.GT.1) THEN                                                 
                CALL NORDLUND(flag,xN,yN,nPcav+1,Plast(iY),                     
     &                                      m,result2,error)                    
                IF (error.NE.0) GOTO 999                                        
              ELSE                                                              
                result2 = 0.0                                                   
              END IF                                                            
c             store current matrix element                                      
              IF (m.EQ.0)                                                       
     &        m0(iL,iY,iW) = 0.5*Y(iY)*Y(iY)*(result1 + result2)                
              IF (m.EQ.1)                                                       
     &        m1(iL,iY,iW) = 0.5*Y(iY)*Y(iY)*(result1 + result2)                
            END DO                                                              
          END DO                                                                
        END DO                                                                  
c     =============================                                             
c     end of loop over wavelengths                                              
      END DO                                                                    
c     save Y and P grids to Yok and Pok, they are needed for analysis           
c     in cases when requirement for finer grids cannot be satisfied and         
c     previous solution is used for output                                      
      nYok = nY                                                                 
      DO iY = 1, nY                                                             
        Yok(iY) = Y(iY)                                                         
      END DO                                                                    
      nPok = nP                                                                 
      DO iP = 1, nP                                                             
        Pok(iP) = P(iP)                                                         
      END DO                                                                    
c -----------------------------------------------------------------------       
999   RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE NORDLUND(flag,x,y,N1,N2,m,intydx,error)                        
c =======================================================================       
c This subroutine calculates integral I(x**m*y(x)*dx). Both y and x are         
c 1D arrays, y(i), x(i) with i=1,nM (nM comes from 'paramslb.inc'). Lower       
c and upper integration limits are x(N1) and x(N2), respectively. The           
c method used is approximation of y(x) by a piecewise cubic spline (see         
c Nordlund: Spherical Transfer with Single-Ray Approximation, in                
c 'Methods in radiative transfer', ed. W. Kalkofen, Cambridge University        
c Press, 1984). The resulting integral is sum of y(i)*w(i), i=N1,N2.            
c Weights w(i) are determined from the array x(i). Here w(i) = wSimp(i)         
c + wCorr(i). To improve accuracy of angular integration in radiative           
c transfer, if flag=1 the contribution of the last Nanal (specified             
c below) points is calculated in subroutine ANALINT which fits a                
c function of the form: y = P(x) + d/sqrt(1-x*x), where P(x) is the             
c polynomial of order Nanal-1, to these points and evaluates integral           
c analytically.                                        [Z.I., Nov. 1995]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER i, flag, N1, N2, N2n, Nanal, m, first, error                      
      DOUBLE PRECISION x(npP), y(npP), wSimp(npP), wCorr(npP), wC1, wC2,        
     &        am, intydx, xaux(4), yaux(4), aux                                 
c -----------------------------------------------------------------------       
      error = 0                                                                 
c parameter 'first' selects choice for derivatives at boundary points.          
c For first.EQ.0 derivatives are 0 and first*(y2-y1)/(x2-x1) otherwise.         
c first=1 works much better for functions encountered here.                     
      first = 1                                                                 
c number of points for analytic integration                                     
      Nanal = 4                                                                 
c do not include points for analytic integration                                
      IF (flag.EQ.1.AND.N2.GT.N1+Nanal) THEN                                    
        N2n = N2 - Nanal + 1                                                    
      ELSE                                                                      
        N2n = N2                                                                
      END IF                                                                    
c set integral to 0 and accumulate result in the loop                           
      intydx = 0.0                                                              
c generate weighting factors, w(i), and integrate in the same loop              
      DO i = N1, N2n                                                            
c first usual Simpson factors (Nordlund, eq. III-14, first term)...             
        IF (i.NE.N1.AND.i.NE.N2n) THEN                                          
          wSimp(i) = 0.5 * (x(i+1)-x(i-1))                                      
        ELSE                                                                    
          IF (i.eq.N1) wSimp(i) = 0.5 * (x(N1+1)-x(N1))                         
          IF (i.eq.N2n) wSimp(i) = 0.5 * (x(N2n)-x(N2n-1))                      
        END IF                                                                  
c ... and then correction term for cubic spline (Nordlund, eq. III-14,          
c second term and eq. III-16) (wC1 and wC2 are auxiliary quantities)            
        IF (i.GT.N1+1) THEN                                                     
          wC1 = x(i) - 2.0*x(i-1) + x(i-2)                                      
        ELSE                                                                    
          IF (i.EQ.N1) wC1 =  first * (x(N1+1) - x(N1))                         
          IF (i.EQ.N1+1) wC1 = first * (x(N1) - x(N1+1))                        
        ENDIF                                                                   
        IF (i.LE.(N2n-2)) THEN                                                  
          wC2 = x(i+2) - 2.0*x(i+1) + x(i)                                      
        ELSE                                                                    
          IF (i.EQ.(N2n-1)) wC2 = first * (x(N2n-1) - x(N2n))                   
          IF (i.EQ.N2n) wC2 = first * (x(N2n) - x(N2n-1))                       
        ENDIF                                                                   
        wCorr(i) = (wC1 - wC2) / 12.                                            
c       add contribution to the integral                                        
        IF (m.EQ.0) THEN                                                        
          intydx = intydx + y(i) * (wSimp(i) + wCorr(i))                        
         ELSE IF(m.EQ.1) THEN                                                   
            intydx = intydx + x(i)*y(i)*(wSimp(i) + wCorr(i))                   
         ELSE IF(m.EQ.2) THEN                                                   
            intydx = intydx + x(i)*x(i)*y(i)*(wSimp(i) + wCorr(i))              
        END IF                                                                  
      END DO                                                                    
c     change the sign (x [i.e. mu] array is in descending order!!!)             
      intydx = -intydx                                                          
c     if flag=1 use analytic approximation for the last Nanal points            
      IF (flag.EQ.1.AND.N2n.GT.N1+Nanal) THEN                                   
c       generate auxiliary arrays for ANALINT                                   
        DO i=1,Nanal                                                            
          xaux(i) = x(N2n+Nanal-i)                                              
          yaux(i) = y(N2n+Nanal-i)                                              
        END DO                                                                  
c     calculate the contribution of the last Nanal points                       
c       produce REAL copy of m                                                  
        am = 1.0*(m)                                                            
        CALL ANALINT(Nanal,xaux,yaux,am,aux,error)                              
        IF(error.NE.0) RETURN                                                   
c       add the contribution of the last Nanal points                           
        intydx = intydx + aux                                                   
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE OccltMSG                                                       
c =======================================================================       
c     Prints a message informing the user about the min Teff required to        
c     neglect occultation by the central source.                                
c =======================================================================       
      IMPLICIT NONE                                                             
      CHARACTER*10 Tstrg                                                        
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER nYok, nPok                                                        
      DOUBLE PRECISION                                                          
     &       Ude(npL,npY), Uds(npL,npY), Us(npL,npY), fde(npL,npY),             
     &       fds(npL,npY), fs(npL,npY), Utot(npL,npY), ftot(npL,npY),           
     &       Td(npG,npY), Ubol(npY), fbol(npY), Uchck(npL,npY),                 
     &       UbolChck(npY), Spectrum(npL), SpecChar(30,99), tauF(npY),          
     &       Intens(npL,npP+2), IntOut(20,npP+2), SmC(30,99),                   
     &       Yok(npY), Pok(npP), tauOut(npL), tauZout(npP+2), F1,               
     &       tr(npY), rg(npG,npY), fsL(npL,npY), fsR(npL,npY), Eps(npY)         
      COMMON /solution/ Ude, Uds, Us, fde, fds, fs, Utot, ftot, Td,             
     &       Ubol, fbol, Uchck, UbolChck, Spectrum, SpecChar, tauF,             
     &       Intens, IntOut, SmC, Yok, Pok, tauOut, tauZout, F1,                
     &       tr, rg, fsL, fsR, Eps, nYok, nPok                                  
      INTEGER iL                                                                
      DOUBLE PRECISION qaux(npL), qaux2(npL), res1, res2, Te_min,               
     &        mx, Psitn, Planck, xP                                             
c -----------------------------------------------------------------------       
c     Estimate min Teff required to neglect occultation (eq.(5) in Manual):     
      DO iL = 1, nL                                                             
        qaux(iL) = SigmaA(1,iL) * Us(iL,1) / lambda(iL)                         
        xP = 14400.0 / Tsub(1) / lambda(iL)                                     
        qaux2(iL) = SigmaA(1,iL) * Planck(xP) / lambda (iL)                     
      END DO                                                                    
      CALL Simpson(npL,1,nL,lambda,qaux,res1)                                   
      CALL Simpson(npL,1,nL,lambda,qaux2,res2)                                  
c     Approximate Psi for opt.thin case:                                        
      Psitn = res1/res2                                                         
      mx = Tsub(1)*sqrt(sqrt(4./Psitn))                                         
      IF(Tsub(1).LT.mx) THEN                                                    
       Te_min = 2. * mx                                                         
      ELSE                                                                      
       Te_min = 2. * Tsub(1)                                                    
      END IF                                                                    
      CALL GetFS(Te_min,0,1,Tstrg)                                              
      write(12,*)                                                               
     & ' ====================================================  '                
      write(12,*)                                                               
     & ' For compliance with the point-source assumption, the'                  
      write(12,*)                                                               
     & ' following results should only be applied to sources '                  
      write(12,'(a37,a5,a3)')                                                   
     & '  whose effective temperature exceeds ',Tstrg, ' K.'                    
      write(12,*)                                                               
     & ' ===================================================='                  
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE Pgrid(pstar,iPstar,error)                                      
c =======================================================================       
c     After having the Y grid, generate the P grid (impact parameters)          
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iYfirst, YPequal, Plast                                           
      DIMENSION iYfirst(npP), YPequal(npP), Plast(npY)                          
      COMMON /Yfirst/ iYfirst, YPequal, Plast                                   
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER i, iP, j, error, Naux, iPstar, NinsLoc                            
      DOUBLE PRECISION pstar, delP                                              
c -----------------------------------------------------------------------       
      error = 0                                                                 
      P(1) = 0.0                                                                
      Naux = 1                                                                  
      iYfirst(Naux) = 1                                                         
      YPequal(Naux) = 1                                                         
c     grid within the cavity (points concentrated towards p=1, SQRT!)           
      nPcav = Ncav + Naux                                                       
      DO iP = Naux+1, nPcav                                                     
        P(iP) = P(Naux) + (1.0-P(Naux))*sqrt((iP-Naux)/(Ncav+1.0))              
        iYfirst(iP) = 1                                                         
        YPequal(iP) = 1                                                         
      END DO                                                                    
c     insert two rays if the stellar disk size is finite                        
      IF (pstar.GT.0.0) THEN                                                    
        CALL INSERT(pstar,iP,iPstar)                                            
        DO i = nPcav+1, iP                                                      
          iYfirst(i) = 1                                                        
          YPequal(i) = 1                                                        
        END DO                                                                  
        nPcav = iP                                                              
        ELSE                                                                    
        iP = nPcav                                                              
        iPstar = 0                                                              
      ENDIF                                                                     
c     impact parameters for the envelope (tangential to the radial grid)        
      DO i =1, nY-1                                                             
        Plast(i) = iP + 1                                                       
c       points close to the inner cavity should have more impact param.         
        NinsLoc = Nins                                                          
c       if tau is large add a few more points to the impact param. grid         
        IF (TAUmax.GT.10) THEN                                                  
          IF (i.LE.5) NinsLoc = Nins + 1                                        
        END IF                                                                  
        IF (TAUmax.GT.100) THEN                                                 
          IF (i.LE.10) NinsLoc = Nins + 1                                       
          IF (i.LE.5) NinsLoc = Nins + 2                                        
        END IF                                                                  
c       increment in P                                                          
        delP = (Y(i+1) - Y(i)) / (NinsLoc+1.0)                                  
        DO j = 1, NinsLoc+1                                                     
          iP = iP + 1                                                           
          iYfirst(iP) = i                                                       
          IF (j.eq.1) THEN                                                      
            YPequal(iP) = 1                                                     
            ELSE                                                                
            YPequal(iP) = 0                                                     
          END IF                                                                
          P(iP) = Y(i) + (j-1)*delP                                             
        END DO                                                                  
      END DO                                                                    
c     number of points for the P grid                                           
      nP = iP + 1                                                               
      P(nP) = Y(nY)                                                             
      iYfirst(nP) = nY                                                          
      YPequal(nP) = 1                                                           
      Plast(nY) = nP                                                            
c     check that the P grid is not too large                                    
      IF (nP.GT.npP) THEN                                                       
        error = 2                                                               
        CALL LINE(0,2,12)                                                       
        write(12,'(a)')' This model terminated because needed accuracy'         
        write(12,'(a,i6)')' results in too many points in P grid:',             
     &                    nP                                                    
        write(12,'(a,i4,a)')'          (See PARAMET.inc, npP =',npP,')'         
        CALL LINE(0,2,12)                                                       
        iERROR = iERROR + 1                                                     
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      DOUBLE PRECISION FUNCTION Planck(x)                                       
c =======================================================================       
c This function evaluates the Planck function multiplied by wavelength          
c and normalized by sigma*T^4/Pi.                      [Z.I., Mar. 1996]        
c =======================================================================       
      IMPLICIT none                                                             
      DOUBLE PRECISION x                                                        
c -----------------------------------------------------------------------       
      IF (x.GT.100) THEN                                                        
        Planck = 0.0                                                            
        ELSE                                                                    
        IF (x.LT.0.00001) THEN                                                  
          Planck = 0.154 * x**3.0                                               
          ELSE                                                                  
          Planck = 0.154 * x**4.0 / (dexp(x) - 1)                               
        END IF                                                                  
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE RADTRANSF(pstar,iPstar,TAUlim,nG,ETAzp,FbolOK,deviat,          
     &                     error,iterFbol,model)                                
c =======================================================================       
c This subroutine solves the continuum radiative transfer problem for a         
c spherically symmetric envelope.                      [Z.I., Nov. 1995]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER nYok, nPok                                                        
      DOUBLE PRECISION                                                          
     &       Ude(npL,npY), Uds(npL,npY), Us(npL,npY), fde(npL,npY),             
     &       fds(npL,npY), fs(npL,npY), Utot(npL,npY), ftot(npL,npY),           
     &       Td(npG,npY), Ubol(npY), fbol(npY), Uchck(npL,npY),                 
     &       UbolChck(npY), Spectrum(npL), SpecChar(30,99), tauF(npY),          
     &       Intens(npL,npP+2), IntOut(20,npP+2), SmC(30,99),                   
     &       Yok(npY), Pok(npP), tauOut(npL), tauZout(npP+2), F1,               
     &       tr(npY), rg(npG,npY), fsL(npL,npY), fsR(npL,npY), Eps(npY)         
      COMMON /solution/ Ude, Uds, Us, fde, fds, fs, Utot, ftot, Td,             
     &       Ubol, fbol, Uchck, UbolChck, Spectrum, SpecChar, tauF,             
     &       Intens, IntOut, SmC, Yok, Pok, tauOut, tauZout, F1,                
     &       tr, rg, fsL, fsR, Eps, nYok, nPok                                  
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER startyp(2), Nlamtr(2), nBB(2)                                     
      CHARACTER nameStar(2)*235                                                 
      DOUBLE PRECISION Tstar, lamtr(2,101), klam(2,100), Tbb(2,10),             
     &                 rellum(2,10), mu1, ksi, mu2, xSiO, r1rs                  
      COMMON /source/ Tstar, lamtr, klam, Tbb, rellum, mu1, ksi, mu2,           
     &                xSiO, r1rs, startyp, Nlamtr, nBB, nameStar                
      INTEGER error, iPstar, nG, Conv, iter, FbolOK, iterFbol, iaux, iY,        
     &        Fconv, Uconv, BolConv, model, itnum, itlim                        
      DOUBLE PRECISION pstar, ETAzp(npP,npY),mat0(npL,npY,npY), maxFerr,        
     &       mat1(npL,npY,npY), Em(npL,npY),  alpha(npG,npY), deviat,           
     &       Uold(npL,npY), fbolold(npY), mifront(npL,npP,npY), dmaxU,          
     &       miback(npL,npP,npY), fsbol(npY), fdbol(npY), TAUlim, dmaxF         
c -----------------------------------------------------------------------       
c     generate, or improve, or do not touch the Y and P grids                   
      IF (iterETA.EQ.1.OR.iterFbol.GT.1) THEN                                   
        IF (iterFbol.EQ.1) THEN                                                 
c         first time generate grids                                             
          CALL SetGrids(pstar,iPstar,error,TAUlim)                              
          IF (error.NE.0) goto 999                                              
          IF (iVerb.EQ.2) write(*,*) 'Done with SetGrids'                       
        ELSE                                                                    
c         or improve the grid from previous iteration                           
          CALL ChkFlux(fBol,accuracy,iaux,error,ETAzp)                          
          IF (error.NE.0) goto 999                                              
c         generate new impact parameter grid                                    
c         increase the number of rays through the cavity                        
          IF (Ncav.LT.80) THEN                                                  
            Ncav = 2 * Ncav                                                     
            IF (iX.NE.0) write(18,'(a20,i3)')' Ncav increased to:',Ncav         
          END IF                                                                
c         increase the number of rays per y-grid interval                       
          IF (iterFbol.EQ.3.AND.Nins.EQ.2) THEN                                 
            Nins = Nins + 1                                                     
            IF (iX.NE.0) write(18,'(a20,i3)')' Nins increased to:',Nins         
          END IF                                                                
          CALL Pgrid(pstar,iPstar,error)                                        
c         if P grid is not OK end this model                                    
          IF (error.NE.0) goto 999                                              
          IF (iX.NE.0) THEN                                                     
            write(18,'(a23,i3)')' Y grid improved, nY =',nY                     
            write(18,'(a23,i3)')'                  nP =',nP                     
            write(18,'(a23,i3)')'                Nins =',Nins                   
            write(18,'(a23,i3)')'                Ncav =',Ncav                   
          END IF                                                                
        END IF                                                                  
      ELSE                                                                      
        IF (iX.NE.0) write(18,*)' Using same Y and P grids'                     
      END IF                                                                    
c     generate spline coefficients for ETA                                      
      CALL setupETA                                                             
c     evaluate ETAzp                                                            
      CALL getETAzp(ETAzp)                                                      
c     generate albedo through the envelope                                      
      CALL getOmega(nG)                                                         
c     generate stellar moments                                                  
      CALL Star(pstar,ETAzp,error)                                              
      IF (iVerb.EQ.2) write(*,*) 'Done with Star'                               
c     issue a message in fname.out about the condition for neglecting           
c     occultation:                                                              
      IF(model.eq.1.AND.iterFbol.eq.1.AND.iterEta.eq.1) CALL OccltMSG           
c     generate the first approximations for Td and alpha                        
      CALL InitTemp(ETAzp,nG,alpha)                                             
      IF (iVerb.EQ.2) write(*,*) 'Done with InitTemp'                           
c     find radiative transfer matrices                                          
      IF (iX.NE.0) write(18,*)' Calculating weight matrices'                    
      IF (iVerb.EQ.2) write(*,*) 'Calculating weight matrices'                  
      CALL Matrix(ETAzp,pstar,iPstar,mat0,mat1,mifront,miback)                  
      Conv = 0                                                                  
      iter = 0                                                                  
c     itlim is an upper limit on number iterations                              
      itlim = 10000                                                             
      IF (iX.NE.0) write(18,*)' Weight matrices OK, calculating Tdust'          
      IF (iVerb.EQ.2) write(*,*)' Weight matrices OK, calculating Tdust'        
      IF (iInn.eq.1) THEN                                                       
        write(38,'(a8,i5)') '    nY= ',nY                                       
        write(38,*) '    iter   maxFerr     dmaxU       dmaxF'                  
      END IF                                                                    
c     === Iterations over dust temperature =========                            
      DO WHILE (Conv.EQ.0.AND.iter.LE.itlim)                                    
        iter = iter + 1                                                         
c       find emission term                                                      
        CALL Emission(0,nG,Td,alpha,abund,Us,Em)                                
c       solve for Utot                                                          
        CALL Invert(nY,nL,mat0,Us,Em,omega,Utot,Uold,error)                     
        IF(error.NE.0) goto 999                                                 
c       find new Td and alpha, and check convergence                            
        CALL FindTemp(1,Utot,nG,Td,alpha)                                       
c       --------------------------------------                                  
c       every itnum-th iteration check convergence:                             
          IF (iter.GT.80) THEN                                                  
             itnum = 10                                                         
            ELSE                                                                
             itnum = 6                                                          
          END IF                                                                
c       first find 'old' flux (i.e. in the previous iteration)                  
        IF (MOD(iter+1,itnum).EQ.0) THEN                                        
          CALL Multiply(1,npY,nY,npL,nL,mat1,Utot,omega,0,fs,fds)               
          CALL Multiply(0,npY,nY,npL,nL,mat1,Em,omega,0,fs,fde)                 
          CALL Add(npY,nY,npL,nL,fs,fds,fde,ftot)                               
c         find bolometric flux                                                  
          CALL Bolom(ftot,fbolold)                                              
        END IF                                                                  
        IF (MOD(iter,itnum).EQ.0) THEN                                          
c         first calculate total flux                                            
          CALL Multiply(1,npY,nY,npL,nL,mat1,Utot,omega,0,fs,fds)               
          CALL Multiply(0,npY,nY,npL,nL,mat1,Em,omega,0,fs,fde)                 
          CALL Add(npY,nY,npL,nL,fs,fds,fde,ftot)                               
c         find bolometric flux                                                  
          CALL Bolom(ftot,fbol)                                                 
c         check convergence of bolometric flux                                  
          CALL Converg1(nY,accFbol,dynrange,fbolold,fbol,Fconv,dmaxF)           
c         check convergence of energy density                                   
          CALL Converg2(nY,nL,accConv,dynrange,Uold,Utot,Uconv,dmaxU)           
c         find maximal fbol error                                               
          CALL FindErr(fbol,maxFerr,nY)                                         
c ------  printout of errors and convergence with iter.(inner flag): -------    
          IF(iInn.EQ.1) THEN                                                    
            write(38,'(i7,1p,3e12.4)') iter,maxFerr,dmaxU,dmaxF                 
          END IF                                                                
c --------------------------------------------------------------                
          IF (maxFerr.LE.accuracy) THEN                                         
            BolConv = 1                                                         
          ELSE                                                                  
            BolConv = 0                                                         
          END IF                                                                
c         total criterion for convergence: Utot must converge, and ftot         
c         must either converge or have the required accuracy                    
          IF (Uconv*(Fconv+BolConv).GT.0) Conv = 1                              
        END IF                                                                  
c       --------------------------------------                                  
      END DO                                                                    
c     === The End of Iterations over Td ===                                     
      IF (iX.NE.0) THEN                                                         
        IF (iter.LT.itlim) write(18,*)                                          
     &     ' Convergence achieved, number of'                                   
        write(18,'(a34,i4)')                                                    
     &     ' iterations over energy density: ',iter                             
        write(18,'(a30,1p,e8.1)')                                               
     &     ' Flux conservation OK within:',maxFerr                              
        IF (iter.GE.itlim) THEN                                                 
          CALL MSG(2)                                                           
          iWARNING = iWARNING + 1                                               
        END IF                                                                  
      END IF                                                                    
c     calculate the emission term for the converged Td                          
      CALL Emission(0,nG,Td,alpha,abund,Us,Em)                                  
c     calculate flux                                                            
      CALL Multiply(1,npY,nY,npL,nL,mat1,Utot,omega,0,fs,fds)                   
      CALL Multiply(0,npY,nY,npL,nL,mat1,Em,omega,0,fs,fde)                     
      CALL Add(npY,nY,npL,nL,fs,fds,fde,ftot)                                   
      CALL Bolom(ftot,fbol)                                                     
c     check whether, and how well, is bolometric flux conserved                 
      CALL ChkBolom(fbol,accuracy,deviat,FbolOK)                                
c     calculate additional output quantities                                    
c     1) energy densities                                                       
      CALL Multiply(1,npY,nY,npL,nL,mat0,Utot,omega,0,Us,Uds)                   
      CALL Multiply(0,npY,nY,npL,nL,mat0,Em,omega,0,Us,Ude)                     
      CALL Add(npY,nY,npL,nL,Us,Uds,Ude,Uchck)                                  
      CALL Bolom(Utot,Ubol)                                                     
      CALL Bolom(Uchck,UbolChck)                                                
c     2) scaled radial optical depth, tr                                        
      DO iY = 1, nY                                                             
        tr(iY) = ETAzp(1,iY) / ETAzp(1,nY)                                      
      END DO                                                                    
c     3) calculate intensity (at the outer edge) if required                    
      IF(iC.NE.0) THEN                                                          
       CALL FindInt(nG,alpha,ETAzp)                                             
       IF (iVerb.EQ.2) write(*,*) 'Done with FindInt'                           
      END IF                                                                    
c     if needed convolve intensity with the PSF                                 
      IF (iPSF.NE.0) THEN                                                       
        CALL Convolve(IntOut)                                                   
        IF (iVerb.EQ.2) write(*,*) 'Done with Convolve'                         
      END IF                                                                    
c     if needed find the visibility function                                    
      IF (iV.NE.0) THEN                                                         
        CALL Visibili(IntOut)                                                   
        IF (iVerb.EQ.2) write(*,*) 'Done with Visibili'                         
      END IF                                                                    
c ============ if the inner flag iInn=1:  =========                             
      IF(iX.GE.1 .AND. iInn.EQ.1) THEN                                          
        CALL Bolom(fs,fsbol)                                                    
        CALL ADD2(fds,fde,fdbol,nY)                                             
        write(18,'(a11,1p,E11.3)')'   TAUfid =',TAUfid                          
        write(18,'(a11,1p,E11.3)')'  MaxFerr =',maxFerr                         
        write(18,*)                                                             
     &'     tr      fbol       fsbol      fdbol       Ubol  '                   
       DO iY = 1, nY                                                            
        write(18,'(1p,5E11.3)') tr(iY), fbol(iY), fsbol(iY),                    
     &                          fdbol(iY), Ubol(iY)                             
       END DO                                                                   
      END IF                                                                    
c =====================                                                         
c -----------------------------------------------------------------------       
999   RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE SetGrids(pstar,iPstar,error,TAU)                               
c =======================================================================       
c Sets the Y and P grids based on GrayBody flux conservation.                   
c                                                     [MN & ZI, July'96]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER nYok, nPok                                                        
      DOUBLE PRECISION                                                          
     &       Ude(npL,npY), Uds(npL,npY), Us(npL,npY), fde(npL,npY),             
     &       fds(npL,npY), fs(npL,npY), Utot(npL,npY), ftot(npL,npY),           
     &       Td(npG,npY), Ubol(npY), fbol(npY), Uchck(npL,npY),                 
     &       UbolChck(npY), Spectrum(npL), SpecChar(30,99), tauF(npY),          
     &       Intens(npL,npP+2), IntOut(20,npP+2), SmC(30,99),                   
     &       Yok(npY), Pok(npP), tauOut(npL), tauZout(npP+2), F1,               
     &       tr(npY), rg(npG,npY), fsL(npL,npY), fsR(npL,npY), Eps(npY)         
      COMMON /solution/ Ude, Uds, Us, fde, fds, fs, Utot, ftot, Td,             
     &       Ubol, fbol, Uchck, UbolChck, Spectrum, SpecChar, tauF,             
     &       Intens, IntOut, SmC, Yok, Pok, tauOut, tauZout, F1,                
     &       tr, rg, fsL, fsR, Eps, nYok, nPok                                  
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER error, iPstar, dummy                                              
      DOUBLE PRECISION pstar, Ugb(npY), fgb(npY), ETAzp(npP,npY),albedo,        
     &       aux, faccs, TAU, accur, delTAUin                                   
c -----------------------------------------------------------------------       
c     store the default value for delTAUsc and facc                             
      faccs = facc                                                              
      delTAUin = delTAUsc                                                       
c     change the delTAUsc seed for the initial Y grid if TAU is large           
      IF (TAU.LT.1.0) delTAUsc = delTAUin * 2.0                                 
      IF (TAU.GE.1.0.and.TAU.LT.5.0) delTAUsc = delTAUin * 1.5                  
      IF (TAU.EQ.5.0) delTAUsc = delTAUin                                       
      IF (TAU.GT.5.0.and.TAU.LT.10.0) delTAUsc = delTAUin / 1.2                 
      IF (TAU.GE.10.0.and.TAU.LT.20.0) delTAUsc = delTAUin / 1.3                
c      The grid is set with TAU=min{TAUlim,TAUmax}, so the lines below are obsol
c      IF (TAU.GE.20.0.and.TAU.LT.30.0) delTAUsc = delTAUin / 1.4               
c      IF (TAU.GE.30.0.and.TAU.LT.50.0) delTAUsc = delTAUin / 1.5               
c      IF (TAU.GE.50.0) delTAUsc = delTAUin / 2.0                               
c     for steep density distributions (RDW, including analyt.approximation):    
      IF (denstyp.EQ.4.OR.RDW) delTAUsc = delTAUin / 1.2               
c     change the facc seed for the initial Y grid if Yout is very small         
      IF (Yout.LT.1000.0) facc = dsqrt(faccs)                                   
      IF (Yout.LT.100.0) facc = dsqrt(facc)                                     
      IF (Yout.LT.10.0) facc = dsqrt(facc)                                      
      IF (Yout.LT.2.0) facc = dsqrt(facc)                                       
      IF (Yout.LT.1.2) facc = dsqrt(facc)                                       
      IF (Yout.LT.1.05) facc = dsqrt(facc)                                      
                                                                                
      albedo = 1.0                                                              
      aux = 1.0                                                                 
c     generate initial grids                                                    
      CALL Ygrid(pstar,iPstar,error)                                            
      IF (error.NE.0) goto 101                                                  
      CALL Pgrid(pstar,iPstar,error)                                            
      IF (error.NE.0) goto 101                                                  
      IF (iX.GE.1) THEN                                                         
         write(18,'(a24,i3)')' Y grid generated, nY =',nY                       
         write(18,'(a24,i3)')'                   nP =',nP                       
         write(18,'(a24,i3)')'                 Nins =',Nins                     
         write(18,'(a24,i3)')'                 Ncav =',Ncav                     
      END IF                                                                    
c     solve for gray body (i.e. pure scattering)                                
      CALL GrayBody(albedo,TAU,Ugb,fgb)                                         
      IF (iVerb.EQ.2) write(*,*) 'Done with GrayBody'                           
c     find the max deviation of fgb (FindRMS called with flag 1)                
      CALL FindRMS(1,fgb,aux,accur,nY)                                          
      IF (iX.GE.1) THEN                                                         
        IF (accur.GT.accuracy) THEN                                             
          write(18,'(a25)')' Grids need improvement:'                           
          write(18,'(a29,1p,e10.3)')                                            
     &                          '                   fTot(nY):',fgb(nY)          
          write(18,'(a29,1p,e10.3)')'      Single wavelength TAU:',TAU          
          write(18,'(a29,1p,e10.3)')                                            
     &                         '          Required accuracy:',accuracy          
        END IF                                                                  
        write(18,'(a29,1p,e10.3)')' Single wavelength accuracy:',accur          
      END IF                                                                    
      IF(accur.GT.accuracy) THEN                                                
c       ChkFlux checks the bolometric flux conservation for the given           
c       grid and decreases the step if conservation is not satisfactory         
        dummy = 5                                                               
        CALL ChkFlux(fgb,accuracy,dummy,error,ETAzp)                            
        IF (error.NE.0) goto 101                                                
c       dummy=5 means everything was fine in ChkFlux                            
        IF (dummy.EQ.5) THEN                                                    
          IF (iX.GE.1) write(18,'(a23,i3)')' Y grid improved, nY =',nY          
c         generate new impact parameter grid                                    
          CALL Pgrid(pstar,iPstar,error)                                        
c         if P grid is not OK end this model                                    
          IF (error.NE.0) goto 101                                              
        ELSE                                                                    
          IF (iX.GE.1) THEN                                                     
            write(18,'(a59,i3)')                                                
     &      ' Although single wavelength accuracy was not satisfactory,'        
            write(18,'(a56,i3)')                                                
     &      ' Y grid could not be improved because npY is too small.'           
            write(18,'(a58,i3)')                                                
     &      ' Continuing calculation with a hope that it will be fine.'         
          END IF                                                                
        END IF                                                                  
      END IF                                                                    
c     return the default value for facc                                         
101    facc = faccs                                                             
       delTAUsc = delTAUin                                                      
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE solve(model,nG,error,ETAzp)                                    
c =======================================================================       
c This subroutine solves the continuum radiative transfer problem for a         
c spherically symmetric envelope.                      [Z.I., Nov. 1995]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER nYok, nPok                                                        
      DOUBLE PRECISION                                                          
     &       Ude(npL,npY), Uds(npL,npY), Us(npL,npY), fde(npL,npY),             
     &       fds(npL,npY), fs(npL,npY), Utot(npL,npY), ftot(npL,npY),           
     &       Td(npG,npY), Ubol(npY), fbol(npY), Uchck(npL,npY),                 
     &       UbolChck(npY), Spectrum(npL), SpecChar(30,99), tauF(npY),          
     &       Intens(npL,npP+2), IntOut(20,npP+2), SmC(30,99),                   
     &       Yok(npY), Pok(npP), tauOut(npL), tauZout(npP+2), F1,               
     &       tr(npY), rg(npG,npY), fsL(npL,npY), fsR(npL,npY), Eps(npY)         
      COMMON /solution/ Ude, Uds, Us, fde, fds, fs, Utot, ftot, Td,             
     &       Ubol, fbol, Uchck, UbolChck, Spectrum, SpecChar, tauF,             
     &       Intens, IntOut, SmC, Yok, Pok, tauOut, tauZout, F1,                
     &       tr, rg, fsL, fsR, Eps, nYok, nPok                                  
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER model,error, iPstar, iterFbol, nG, FbolOK, EtaOK, iY              
      DOUBLE PRECISION pstar, ETAzp(npP,npY), TAUlim, deviat                    
c -----------------------------------------------------------------------       
      IF(iInn.eq.1) THEN                                                        
        write(38,*)'============================================='              
        write(38,'(a7,i5)') ' model= ',model                                    
        write(38,*)'============================================='              
      END IF                                                                    
      IF (iX.NE.0) THEN                                                         
       CALL LINE(0,2,18)                                                        
       write(18,'(a7,i3,a20)') ' model ',model,'  RUN-TIME MESSAGES '           
       CALL LINE(0,1,18)                                                        
       write(18,*)'===========  STARTING SOLVE  ==========='                    
      END IF                                                                    
      error = 0                                                                 
      IF(denstyp.NE.0) THEN                                                     
c     Solve for spherical envelope:                                             
c      temporarily the star is approximated by a point source                   
       pstar = 0.0                                                              
       iPstar = 1                                                               
c      select optical depth for the grid calculation                            
c      based on dynamical range                                                 
       TAUlim = 0.5*log(1.0/dynrange)                                           
c      if actual maximal tau is smaller use that value                          
       IF (TAUmax.LT.TAUlim) TAUlim = TAUmax                                    
c      counter over ETA (for radiatively driven winds only)                     
       iterETA = 0                                                              
       EtaOK = 0                                                                
c      iterations over ETA                                                      
       DO WHILE (EtaOK.EQ.0)                                                    
         iterETA = iterETA + 1                                                  
         IF (iX.NE.0.AND.RDW) THEN                                              
           write(18,*)'----------------------------------------'                
           write(18,*)' ',iterETA,'. iteration over ETA'                        
         END IF                                                                 
         IF (iVerb.EQ.2.AND.RDW)                                                
     &          write(*,*) ' ',iterETA,'. iteration over ETA'                   
c        counter for iterations over bolometric flux conservation               
         iterFbol = 0                                                           
         FbolOK = 0                                                             
         DO WHILE (FbolOK.EQ.0)                                                 
           iterFbol = iterFbol + 1                                              
           IF (iX.NE.0) THEN                                                    
            write(18,*)'  ',iterFbol,'. iteration over Fbol'                    
           END IF                                                               
           IF (iVerb.EQ.2) write(*,*) iterFbol,'. iteration over Fbol'          
c          solve the radiative transfer problem                                 
           CALL RADTRANSF(pstar,iPstar,TAUlim,nG,ETAzp,FbolOK,deviat,           
     &                   error,iterFbol,model)                                  
           IF (iVerb.EQ.2) write(*,*) 'Done with RadTransf'                     
c          error.EQ.3 : file with stellar spectrum not available                
           IF (error.EQ.3) goto 999                                             
c          error.EQ.5 : Singular matrix                                         
           IF (error.EQ.5) goto 999                                             
c          error.EQ.6 : Eta exceeds limitations                                 
           IF (error.EQ.6) goto 999                                             
c          error.EQ.2 : P grid was not produced                                 
           IF (error.EQ.2.AND.iterFbol.EQ.1.AND.iterETA.EQ.1) THEN              
c           if this is the first calculation end this model                     
            iERROR = iERROR + 1                                                 
            goto 999                                                            
           ELSE                                                                 
c           if this is a higher iteration use previous solution                 
            IF (error.EQ.2) THEN                                                
              IF (iX.NE.0.AND.iterFbol.GT.1) THEN                               
              write(18,*)' ======= IMPORTANT WARNING ======== '                 
              write(18,*)' In trying to conserve Fbol reached'                  
              write(18,*)' the limit for grid sizes. Flux is '                  
              write(18,'(a,1p,e9.3)')'  conserved to within ', deviat           
              write(18,*)' Treat all results with caution!'                     
              END IF                                                            
              IF (iX.NE.0.AND.iterFbol.EQ.1) THEN                               
              write(18,*)' ======== IMPORTANT  WARNING ======== '               
              write(18,*)' In trying to converge on ETA reached'                
              write(18,*)' the limit for grid sizes. Flux is '                  
              write(18,'(a,1p,e9.3)')'  conserved to within ', deviat           
              write(18,*)' Treat all results with caution!'                     
              END IF                                                            
              error = 0                                                         
              FbolOK = 2                                                        
              iWARNING = iWARNING + 1                                           
            END IF                                                              
           END IF                                                               
c          just in case...                                                      
           IF (error.EQ.1) THEN                                                 
            IF (iX.NE.0) THEN                                                   
              write(18,*)' *********  FATAL ERROR  *********'                   
              write(18,*)' * Something was seriously wrong *'                   
              write(18,*)' * Contact Z. Ivezic, M. Elitzur *'                   
              write(18,*)' *********************************'                   
            END IF                                                              
            iERROR = iERROR + 1                                                 
            goto 999                                                            
           END IF                                                               
c          if Fbol not conserved try again with a finer grid                    
           IF (FbolOK.EQ.0.AND.iterFbol.LT.10.AND.iX.NE.0) THEN                 
            write(18,*)'  ******** MESSAGE from SOLVE ********'                 
            write(18,*)'  Full solution does not conserve Fbol'                 
            write(18,*)'       Y       TAU/TAUtot        fbol'                  
            DO iY =1, nY                                                        
              write(18,'(1p,3e13.4)')Y(iY),                                     
     &                  ETAzp(1,iY)/ETAzp(1,nY),fbol(iY)                        
            END DO                                                              
            write(18,*)'  Trying again with finer grids'                        
           END IF                                                               
c          if could not conserve Fbol in 10 trials give it up                   
           IF (FbolOK.EQ.0.AND.iterFbol.GE.10) THEN                             
            IF (RDW) THEN                                                       
            IF (iX.NE.0) THEN                                                   
            write(18,*)' **********  WARNING from SOLVE  **********'            
            write(18,*)' Could not obtain required accuracy in Fbol'            
            write(18,'(a26,1p,e10.3)')'  The achieved accuracy is:',            
     &                                 deviat                                   
            write(18,*)' Will try to converge on the dynamics, but '            
            write(18,*)' treat all results with caution !!         '            
            write(18,*)' If accuracy<=0.01, or TAUmax>1000, this   '            
            write(18,*)' code probably cannot do it. Otherwise,    '            
            write(18,*)' please contact Z. Ivezic or M. Elitzur    '            
            write(18,*)' ******************************************'            
            END IF                                                              
            iWARNING = iWARNING + 1                                             
            FbolOK = 1                                                          
            ELSE                                                                
            IF (iX.NE.0) THEN                                                   
            write(18,*)' **********  WARNING from SOLVE  **********'            
            write(18,*)' Could not obtain required accuracy in Fbol'            
            write(18,'(a26,1p,e10.3)')'  The achieved accuracy is:',            
     &                                 deviat                                   
            write(18,*)' !!!!  Treat all results with caution  !!!!'            
            write(18,*)' If accuracy<=0.01, or TAUmax>1000, this   '            
            write(18,*)' code probably cannot do it. Otherwise,    '            
            write(18,*)' please contact Z. Ivezic or M. Elitzur    '            
            write(18,*)' ******************************************'            
            END IF                                                              
            iWARNING = iWARNING + 1                                             
            FbolOK = 2                                                          
            END IF                                                              
           END IF                                                               
c        end of loop over flux conservation                                     
         END DO                                                                 
c        for winds check if ETA has converged...                                
         IF ((RDW).AND.FbolOK.NE.2) THEN                                        
c         ptr(2) is specified in INPUT and controls converg. crit.              
          IF (ptr(2).LT.1.0e-6.AND.iterETA.GT.2)THEN                            
            EtaOK = 1                                                           
          ELSE                                                                  
            CALL WINDS(nG,EtaOK,ETAzp,ftot)                                     
          END IF                                                                
          IF (iterETA.GT.10.AND.EtaOK.EQ.0) THEN                                
            EtaOK = 2                                                           
            iWARNING = iWARNING + 1                                             
            IF (iX.NE.0) THEN                                                   
              write(18,*)' *********  WARNING  *********'                       
              write(18,*)' Could not converge on ETA in '                       
              write(18,*)' 10 iterations.'                                      
              write(18,*)' *********************************'                   
            END IF                                                              
          END IF                                                                
c         ...or otherwise finish right away                                     
         ELSE                                                                   
          EtaOK = 1                                                             
         END IF                                                                 
c      end of loop over ETA                                                     
       END DO                                                                   
      ELSE                                                                      
c      solve for slab case                                                      
       CALL SLBsolve(model,nG,error)                                            
c      error=4 means npY not large enough for oblique illumination grid         
       IF (error.eq.4) THEN                                                     
        CALL MSG(15)                                                            
        iWARNING = iWARNING + 1                                                 
        goto 999                                                                
       END IF                                                                   
      END IF                                                                    
c     analyze the solution and calculate some auxiliary quantities              
      CALL Analysis(model,ETAzp,error)                                          
      IF (iVerb.EQ.2) write(*,*) 'Done with Analysis'                           
      IF (iX.NE.0) THEN                                                         
        write(18,*)' ==== SOLVE successfully completed ====='                   
        write(18,*)' ======================================='                   
      END IF                                                                    
c -----------------------------------------------------------------------       
999   RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE STAR(pstar,ETAzp,error)                                        
c =======================================================================       
c This subroutine generates the stellar moments   [ZI, Nov'95; MN,Mar'99]       
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      INTEGER nYok, nPok                                                        
      DOUBLE PRECISION                                                          
     &       Ude(npL,npY), Uds(npL,npY), Us(npL,npY), fde(npL,npY),             
     &       fds(npL,npY), fs(npL,npY), Utot(npL,npY), ftot(npL,npY),           
     &       Td(npG,npY), Ubol(npY), fbol(npY), Uchck(npL,npY),                 
     &       UbolChck(npY), Spectrum(npL), SpecChar(30,99), tauF(npY),          
     &       Intens(npL,npP+2), IntOut(20,npP+2), SmC(30,99),                   
     &       Yok(npY), Pok(npP), tauOut(npL), tauZout(npP+2), F1,               
     &       tr(npY), rg(npG,npY), fsL(npL,npY), fsR(npL,npY), Eps(npY)         
      COMMON /solution/ Ude, Uds, Us, fde, fds, fs, Utot, ftot, Td,             
     &       Ubol, fbol, Uchck, UbolChck, Spectrum, SpecChar, tauF,             
     &       Intens, IntOut, SmC, Yok, Pok, tauOut, tauZout, F1,                
     &       tr, rg, fsL, fsR, Eps, nYok, nPok                                  
      INTEGER startyp(2), Nlamtr(2), nBB(2)                                     
      CHARACTER nameStar(2)*235                                                 
      DOUBLE PRECISION Tstar, lamtr(2,101), klam(2,100), Tbb(2,10),             
     &                 rellum(2,10), mu1, ksi, mu2, xSiO, r1rs                  
      COMMON /source/ Tstar, lamtr, klam, Tbb, rellum, mu1, ksi, mu2,           
     &                xSiO, r1rs, startyp, Nlamtr, nBB, nameStar                
      DOUBLE PRECISION TAUslb(npL,npY),fsbol(npY),fpbol(npY),fmbol(npY),        
     &                 AveDev, RMS, maxFerr, fmax, fmed                         
      COMMON /slab/ TAUslb, fsbol, fpbol, fmbol, AveDev, RMS, maxFerr,          
     &                 fmax, fmed                                               
      CHARACTER*235 line                                                        
      INTEGER iY, iL, ios1, iLs, nLs, k, kstop, i, error, Nlambdam              
c     Nlambdam is the max number entries for a user supplied stellar spectrum   
      PARAMETER (Nlambdam = 10000)                                              
      DOUBLE PRECISION pstar, ETAzp(npP,npY), zeta, Bb, x, Planck, a, b,        
     &        lambdaS(Nlambdam), Llamstar(Nlambdam), Stellar(Nlambdam),         
     &        Lstar, fL(100), fpl(npL), llS(Nlambdam), lS(Nlambdam),            
     &        expow, dyn2, EMfunc                                               
c -----------------------------------------------------------------------       
c     for  startyp.GE.4 stellar spectrum is read from the file 'nameStar'       
c     it is unit=3 (1 is the input file)                                        
      IF ((startyp(1).GE.4).AND.(startyp(1).LE.6)) THEN                         
        open(3,ERR=998,file=nameStar(1),STATUS='OLD')                           
        rewind(3)                                                               
        read(3,'(a235)',ERR=998) line                                           
        read(3,'(a235)',ERR=998) line                                           
        read(3,'(a235)',ERR=998) line                                           
        ios1 = 0                                                                
        iLs = 0                                                                 
        DO WHILE (ios1.ge.0)                                                    
          read(3,*,END=900,ERR=998,iostat=ios1)a,b                              
          IF (ios1.ge.0) THEN                                                   
            iLs = iLs + 1                                                       
            lambdaS(iLs) = a                                                    
            IF (a.LE.0.0) goto 998                                              
c           it is assumed that Llamstar is L_Lambda, but...                     
c           if startyp.EQ.4 then file gives lambda*L_lambda                     
            IF (startyp(1).EQ.4) Llamstar(iLs) = b / a                          
c           if startyp.EQ.5 then file gives L_lambda                            
            IF (startyp(1).EQ.5) Llamstar(iLs) = b                              
c           if startyp.EQ.6 then file gives Lnu=lambda**2*L_lambda              
            IF (startyp(1).EQ.6) Llamstar(iLs) = b / a / a                      
           END IF                                                               
        END DO                                                                  
900     close(3)                                                                
        IF (iLs.LT.2) goto 998                                                  
        nLs = iLs                                                               
c       if input wavelengths in descending order turn them around               
        IF (lambdaS(1).GT.lambdaS(2)) THEN                                      
          DO iLs = 1, nLs                                                       
            llS(iLs) = lambdaS(iLs)                                             
            lS(iLs) = Llamstar(iLs)                                             
          END DO                                                                
          DO iLs = 1, nLs                                                       
            lambdaS(iLs) = llS(nLs+1-iLs)                                       
            Llamstar(iLs) = lS(nLs+1-iLs)                                       
          END DO                                                                
        END IF                                                                  
c       normalize stellar spectrum                                              
        CALL Simpson(Nlambdam,1,nLs,lambdaS,Llamstar,Lstar)                     
c       generate dimensionless stellar spectrum                                 
        DO iLs = 1, nLs                                                         
          Stellar(iLs) = lambdaS(iLs) * Llamstar(iLs) / Lstar                   
        END DO                                                                  
        ELSE                                                                    
c       if startyp.EQ.3 generate power-law spectrum                             
        IF (startyp(1).EQ.3) THEN                                               
          fL(1) = 1.0                                                           
          IF (Nlamtr(1).GT.1) THEN                                              
            DO i = 2, Nlamtr(1)                                                 
              fL(i) = fL(i-1) * (lamtr(1,i-1)/lamtr(1,i))**klam(1,i-1)          
            END DO                                                              
          END IF                                                                
          DO iL = 1, nL                                                         
            IF ((lambda(iL)-lamtr(1,1))*(lambda(iL)-                            
     &                              lamtr(1,Nlamtr(1)+1)).LE.0.0) THEN          
              kstop = 0                                                         
              k = 0                                                             
              DO WHILE (kstop.EQ.0)                                             
                k = k + 1                                                       
                IF (lambda(iL).GE.lamtr(1,k)) THEN                              
                  kstop = 1                                                     
                  fpl(iL) = fL(k) * (lamtr(1,k)/lambda(iL))**klam(1,k)          
                END IF                                                          
              END DO                                                            
            ELSE                                                                
              fpl(iL) = 0.0                                                     
            END IF                                                              
          END DO                                                                
        END IF                                                                  
      END IF                                                                    
      DO iY = 1, nY                                                             
c       effect of the star's finite size                                        
        IF (pstar.GT.0.0) THEN                                                  
          zeta = 2.*(1.-sqrt(1.-(pstar/Y(iY))**2.))*(Y(iY)/pstar)**2.           
        ELSE                                                                    
          zeta = 1.0                                                            
        END IF                                                                  
c       loop over wavelengths                                                   
        DO iL = 1, nL                                                           
          IF (startyp(1).EQ.1) THEN                                             
            Bb = 0.0                                                            
            DO k = 1, nBB(1)                                                    
              x = 14400.0 / lambda(iL) / Tbb(1,k)                               
              Bb = Bb + rellum(1,k)*Planck(x)                                   
            END DO                                                              
          ELSE IF (startyp(1).EQ.2) THEN                                        
              Bb = EMfunc(lambda(iL),Tbb(1,1),xSiO)                             
          ELSE IF (startyp(1).EQ.3) THEN                                        
              Bb = fpl(iL)                                                      
          ELSE IF (lambda(iL).GT.lambdaS(nLs)) THEN                             
c             for lambda longer than the longest entry in nameStar              
c             assume Rayleigh-Jeans tail                                        
              Bb = Stellar(nLs) * (lambdaS(nLs)/lambda(iL))**3.                 
          ELSE IF (lambda(iL).LT.lambdaS(1)) THEN                               
c             if shorter than the shortest assume 0                             
              Bb = 0.0                                                          
          ELSE                                                                  
c             if within limits interpolate                                      
              CALL LinInter(Nlambdam,nLs,lambdaS,Stellar,                       
     &                                             lambda(iL),iLs,Bb)           
          END IF                                                                
c   ----------  done with stellar spectrum ---------------                      
c         stellar part of en.density and flux                                   
          expow = ETAzp(1,iY)*TAUtot(iL)                                        
          Us(iL,iY) = Bb * zeta * dexp(-expow)                                  
c         flux                                                                  
          fs(iL,iY) = Bb * dexp(-expow)                                         
c         need to carry fsL to be compatible with the slab case                 
          fsL(iL,iY) = fs(iL,iY)                                                
          fsR(iL,iY) = 0.0                                                      
          dyn2 = dynrange*dynrange                                              
          IF (Us(iL,iY).LT.dyn2) Us(iL,iY) = 0.0                                
          IF (fs(iL,iY).LT.dyn2.AND.denstyp.NE.0) fs(iL,iY) = 0.0               
        END DO                                                                  
      END DO                                                                    
c     normalize stellar quantities with the stellar bolometric flux             
      CALL Bolom(fsL,fsbol)                                                     
      DO iY = 1, nY                                                             
        DO iL = 1, nL                                                           
           Us(iL,iY) = Us(iL,iY) / fsbol(1)                                     
           fs(iL,iY) = fs(iL,iY) / fsbol(1)                                     
           fsL(iL,iY) = fsL(iL,iY) / fsbol(1)                                   
        END DO                                                                  
      END DO                                                                    
      error = 0                                                                 
      goto 999                                                                  
998   write(12,*)' *** FATAL ERROR IN DUSTY! *************************'         
      write(12,*)' File with the spectral shape of external radiation:'         
      write(12,'(a2,a70)')'  ',nameStar(1)                                      
      write(12,*)' is missing or not properly formatted?!'                      
      write(12,*)' ***************************************************'         
      error = 3                                                                 
c -----------------------------------------------------------------------       
999   RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE trapzd2(a,b,s,n)                                               
c =======================================================================       
c This function integrates prescribed 8 functions from z=a to z=b with n        
c divisions and stores the results to s(1..8). It is a heavily modified         
c version of subroutine 'trapzd' (Num.Rec.'92).        [MN & ZI, Aug'96]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER it,iC,i,n,j                                                       
      DOUBLE PRECISION s(8),a,b,funcx(8),funca(8),funcb(8),del,sum(8),          
     &       tnm, x, ff, gp, gm                                                 
c -----------------------------------------------------------------------       
      IF (n.eq.1) then                                                          
c       calculate auxiliary functions at a and at b                             
        CALL TWOFUN(a,ff,gp,gm)                                                 
        funca(1) =  gm                                                          
        funca(5) =  gp                                                          
        DO iC = 2, 4                                                            
            funca(iC) = funca(iC-1) * ff                                        
            funca(4+iC) = funca(3+iC) * ff                                      
        END DO                                                                  
        CALL TWOFUN(b,ff,gp,gm)                                                 
        funcb(1) =  gm                                                          
        funcb(5) =  gp                                                          
        DO iC = 2, 4                                                            
          funcb(iC) = funcb(iC-1) * ff                                          
          funcb(4+iC) = funcb(3+iC) * ff                                        
        END DO                                                                  
c       calculate integrals for all 8 functions                                 
        DO i = 1, 8                                                             
          s(i) = 0.5*(b-a)*(funca(i)+funcb(i))                                  
        END DO                                                                  
      ELSE                                                                      
        it=2**(n-2)                                                             
        tnm=1.0*(it)                                                            
        del=(b-a)/tnm                                                           
        x=a+0.5*del                                                             
        DO i=1,8                                                                
          sum(i)=0.0                                                            
        END DO                                                                  
c       calculate contributions of all 'it' divisions                           
        DO j = 1, it                                                            
c         auxiliary functions at x                                              
          CALL TWOFUN(x,ff,gp,gm)                                               
c         generate (8) integrated functions at x                                
          funcx(1) = gm                                                         
          funcx(5) = gp                                                         
          DO iC = 2, 4                                                          
            funcx(iC) = funcx(iC-1) * ff                                        
            funcx(4+iC) = funcx(3+iC) * ff                                      
          END DO                                                                
          DO i=1,8                                                              
            sum(i)=sum(i)+funcx(i)                                              
          END DO                                                                
c         next x                                                                
          x=x+del                                                               
        END DO                                                                  
c       evaluate new value of the integral for all 8 cases                      
        DO i=1,8                                                                
          s(i)=0.5*(s(i)+(b-a)*sum(i)/tnm)                                      
        END DO                                                                  
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE TWOFUN(z,ff,gp,gm)                                             
c =======================================================================       
c This function evaluates auxiliary functions needed in trapzd2.                
c                                           [MN & ZI,Aug'96; MN,Sep'97]         
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER iLaux, iW1, iC                                                    
      DOUBLE PRECISION paux, w1,wL,etaloc, w, IntETA, pp, delTAUzp, z,          
     &        auxw, ff, gp, gm, gp1,gm1                                         
      COMMON /phi2/ paux, w1, wL, delTAUzp, iLaux, iW1                          
c -----------------------------------------------------------------------       
c     local radius                                                              
      w = dsqrt(paux*paux + z*z)                                                
      IF (w.LT.w1) w = w1                                                       
c     find local value for ETA function                                         
      etaloc = 0.0                                                              
      auxw = 1.                                                                 
      DO iC = 1, 4                                                              
        etaloc = etaloc + ETAcoef(iW1,iC)*auxw                                  
        auxw = auxw/w                                                           
      END DO                                                                    
c     ff, i.e. radial optical depth:                                            
      pp = 0.0                                                                  
      ff = IntETA(pp,iW1,wL,w)*TAUtot(iLaux)                                    
c     g functions:                                                              
      gp1 = dexp(IntETA(paux,iW1,w1,w)*TAUtot(iLaux)-delTAUzp)                  
      gm1 = dexp(-IntETA(paux,iW1,w1,w)*TAUtot(iLaux))                          
      gp = etaloc/w/w * gp1                                                     
      gm = etaloc/w/w * gm1                                                     
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE WEIGHTS(TAUaux,iP,iL,nZ,nY,alpha,beta,gamma,delta,             
     &                   wgp,wgm)                                               
c =======================================================================       
c This subroutine calculates weights wgp(iZ,iY) and wgm(iZ,iY) for              
c integrations:                                                                 
c INT(S(w)*exp(sign*ETAzp(iP,iZ')/w^2)dETAzp(iP,iZ')]                           
c from ETAzp(iP,iZ) to ETAzp(iP,iZ+1), where w is local radius                  
c corresponding to TAU(iP,iZ'), and sign=1 for wgp and -1 for wgm.              
c Integrals are evaluated as:                                                   
c INT = wg(iZ,1)*S(1) + wg(iZ,2)*S(2) + ... + wg(iZ,nY)*S(nY) with              
c iZ=1..nZ-1. The method is based on approximation of S by cubic spline         
c in radial optical depth given through matrices alpha, beta, gamma and         
c delta (see MYSPLINE).                         [ZI,Dec'95;MN,Sep'97]           
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER iYfirst, YPequal, Plast                                           
      DIMENSION iYfirst(npP), YPequal(npP), Plast(npY)                          
      COMMON /Yfirst/ iYfirst, YPequal, Plast                                   
      INTEGER iP, iL, nZ, nY, iW, iZ, j                                         
      DOUBLE PRECISION alpha(npY,npY), beta(npY,npY), gamma(npY,npY),           
     &        delta(npY,npY), TAUaux(npL,npP,npY), K1p(npY),K2p(npY),           
     &        K3p(npY),K4p(npY), K1m(npY),K2m(npY),K3m(npY),K4m(npY),           
     &        wgp(npY,npY), wgm(npY,npY), waux                                  
c -----------------------------------------------------------------------       
c     generate integrals of 'TAUr**n'                                           
      CALL Kint4(TAUaux,iP,iL,nZ,K1p,K2p,K3p,K4p,K1m,K2m,K3m,K4m)               
c     loop over position on the line of sight                                   
      DO iZ = 1, nZ                                                             
        iW = iYfirst(iP) + iZ - 1                                               
c       loop over radial position                                               
        DO j = 1, nY                                                            
         IF (iZ.GT.1) THEN                                                      
          waux = alpha(iW-1,j)*K1p(iZ) + beta(iW-1,j)*K2p(iZ)                   
          wgp(iZ,j)=waux + gamma(iW-1,j)*K3p(iZ)+delta(iW-1,j)*K4p(iZ)          
         ELSE                                                                   
          wgp(1,j) = 0.0                                                        
         END IF                                                                 
         IF (iZ.LT.nZ) THEN                                                     
          wgm(iZ,j) = alpha(iW,j)*K1m(iZ) + beta(iW,j)*K2m(iZ)                  
          wgm(iZ,j) = wgm(iZ,j)+gamma(iW,j)*K3m(iZ)+delta(iW,j)*K4m(iZ)         
         ELSE                                                                   
          wgm(nZ,j) = 0.0                                                       
         END IF                                                                 
        END DO                                                                  
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE Ygrid(pstar,iPstar,error)                                      
c =======================================================================       
c This subroutine generates the radial (Y) grid.                                
c Ncav is desired number of rays through the inner cavity (<=20), denstype      
c is type of density law (see input file). Yout is the relative thickness,      
c Yout=rout/r1. delTausc is maximum increase of the scaled optical depth        
c tau/tauTot between two successive radial grid points and facc is the          
c maximum ratio of their y coordinates. EtaRat is initialized in SUB Input      
c and limits the ratio of two successive Eta(iY), in order to prevent           
c nonphysical results in case of steep density distribution (pow>2). pstar is   
c the impact parameter for star (0<=p<=1). First few points are prescribed      
c arbitrarily. error is set to 1 if desired delTausc and facc result in too     
c many points (> npY from userpar.inc). This subroutine calls function ETA      
c which evaluates the normalized density profile.    [ZI, Nov'95; MN,Sep'99]    
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iYfirst, YPequal, Plast                                           
      DIMENSION iYfirst(npP), YPequal(npP), Plast(npY)                          
      COMMON /Yfirst/ iYfirst, YPequal, Plast                                   
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER i, istop,error,iPstar, iter, iYdummy, iY, itr, j, ir, jY,         
     &        irmax                                                             
      DOUBLE PRECISION pstar, dyR, dyT, ETA, aux, tausc, EtaTemp(npY),          
     &       ee, Yloc, nh, fac, dyF, yold, ynew, y1, delY1, Ymid,               
     &       TAUloc, rat, dY, Etamin, Etamax                                    
c -----------------------------------------------------------------------       
c     max number iter. over improving ratio of two Eta's                        
      irmax = 20                                                                
c     save old grid and values of Eta (important for denstyp = 5 or 6)          
      IF (nY.GT.0.AND.(RDW)) THEN                                               
        DO iY = 1, nY                                                           
          Yprev(iY) = Y(iY)                                                     
          EtaTemp(iY) = ETAdiscr(iY)                                            
        END DO                                                                  
        nYprev = nY                                                             
      END IF                                                                    
      pstar = pstar                                                             
      iPstar = iPstar                                                           
      y1 = 1.0                                                                  
      iter = 0                                                                  
101   error = 0                                                                 
      iter = iter + 1                                                           
c     resolve inner boundary in tau space                                       
c     from requiring TAU(2) = TAUmax*ETA(1)*(Y(2)-1) =~ 1                       
      Y(1) = y1                                                                 
      IF (TAUmax*ETA(y1).GT.2000.0) THEN                                        
        delY1 = 2.0 / TAUmax / ETA(y1)                                          
      ELSE                                                                      
        delY1 = 1.0 / TAUmax / ETA(y1)                                          
      END IF                                                                    
c     if very thin generate at least about 10-15 pts.                           
      IF (delY1.GT.0.01/ETA(y1)) delY1 = 0.01/ETA(y1)                           
c     do not push it over the limit of spline stability                         
      IF (delY1.LT.0.00005) delY1 = 0.00005                                     
      i = 1                                                                     
      istop = 0                                                                 
      DO WHILE (istop.NE.1)                                                     
        i = i + 1                                                               
        Y(i) = Y(i-1) + delY1                                                   
        IF (Y(i).GT.facc*Y(i-1)) THEN                                           
          Y(i) = facc * Y(i-1)                                                  
        ELSE                                                                    
          delY1 = 2. * delY1                                                    
        END IF                                                                  
        Ymid = dsqrt(Y(i)*Y(i-1))                                               
        TAUloc = TAUmax * (Y(i)-1.0) * ETA(Ymid)                                
        IF (Y(i).GE.1.01.OR.TAUloc.GT.10.0) istop = 1                           
c       in case of shell thinner than 1.01 comment the above line and           
c       uncomment the next line, unless it is a case of RDW at high TauV.       
c       These require more points near the origin and 1.01 is a better limit.   
c       IF (Y(i).GE.1.0001.OR.TAUloc.GT.10.0) istop = 1                         
      END DO                                                                    
      Ynew = Y(i)                                                               
c     some rule of thumb estimates for factor nh                                
      IF (TAUmax.GT.10000.0) THEN                                               
c       extreme taus                                                            
        nh = 15.0                                                               
        Ncav = 80                                                               
      ELSE IF (TAUmax.GT.2000.0) THEN                                           
c         huge taus                                                             
          nh = 10.0                                                             
          Ncav = 40                                                             
      ELSE IF (TAUmax.GT.400.0) THEN                                            
c         large taus                                                            
          nh = 8.0                                                              
          Ncav = 20                                                             
      ELSE                                                                      
c         normal' taus                                                          
          nh = 5.0                                                              
          Ncav = 10                                                             
      END IF                                                                    
c     very small taus                                                           
      IF (TAUmax.LT.10.0) nh = nh / 2.0                                         
      IF (TAUmax.LT.1.0) nh = nh / 2.0                                          
c     empirically: ~1/r needs more points for small y:                          
      IF ((pow-1.4)*(pow-0.6).LE.0.0) nh = nh * 1.5                             
c     same for for steep density distributions:                                 
      IF (denstyp.eq.4.OR.RDW) nh = nh * 1.5                           
      tausc = ETA(Y(i))+ETA(Y(i-1))                                             
      tausc = (Y(i)-Y(i-1)) * 0.5 * tausc                                       
      fac = dexp(-dlog(tausc) / nh)                                             
      istop = 0                                                                 
c     for broken power-laws                                                     
      IF (Ntr.GE.1) THEN                                                        
        itr = 1                                                                 
        DO j = 1, i                                                             
          IF (Y(j).GE.Ytr(itr)) itr = itr + 1                                   
        END DO                                                                  
      END IF                                                                    
c     generate the rest of Y grid points                                        
      DO WHILE (istop.NE.1)                                                     
        i = i + 1                                                               
        Yold = Ynew                                                             
c       find maximal increase in Y allowed by the ratio facc                    
        dyR = Yold * (facc-1.)                                                  
c       find maximal increase in Y allowed by delTausc                          
        dyT = delTausc / ETA(Yold)                                              
c       find maximal increase in Y allowed by the ratio of tausc                
        dyF = tausc*(fac-1.) / ETA(Yold)                                        
c       find new Y                                                              
c        Ynew = Yold + MIN(dyR,dyT,dyF)                                         
        dY = MIN(dyR,dyT,dyF)                                                   
c       Check if the max ratio btw. two Eta values is less than EtaRat          
c       and insert additional y-pts. where necessary. This prevents sharp       
c       drops in Utot(npL,npY) in case of steep Eta's [MN'99].                  
        DO ir = 1 , irmax                                                       
          Ynew = Yold + dY                                                      
          rat = ETA(Yold)/ETA(Ynew)                                             
          IF(rat.GE.1./EtaRat .AND. rat.LE.EtaRat) goto 10                      
          dY = 0.5*dY                                                           
        END DO                                                                  
        CALL MSG(16)                                                            
10      continue                                                                
        Y(i) = Ynew                                                             
c       make sure that all transition points are included in the grid           
        IF (Ntr.GE.1) THEN                                                      
         IF (Y(i).GE.Ytr(itr)) THEN                                             
           Y(i) = Ytr(itr)                                                      
           Ynew = Y(i)                                                          
           itr = itr + 1                                                        
         END IF                                                                 
        END IF                                                                  
        aux = ETA(Ynew)+ETA(Yold)                                               
        aux = (Ynew-Yold) * 0.5 * aux                                           
        tausc = tausc + aux                                                     
c       finish when Yout is reached                                             
        IF (Ynew.GE.Yout) istop = 1                                             
      END DO                                                                    
      Y(i) = Yout                                                               
      nY = i                                                                    
c     insert additional penultimate point to avoid spline oscillations          
      Y(nY+1) = Yout                                                            
      Y(nY) = dsqrt(Y(nY)*Y(nY-1))                                              
      nY = nY + 1                                                               
c     check that outer edge is well resolved in tau space                       
c     (important for flat ETAs)                                                 
      istop = 0                                                                 
      DO WHILE (istop.NE.1)                                                     
        IF ((Yout-Y(nY-1))*TAUmax*ETA(Yout).GT.1.0) THEN                        
          Y(nY+1) = Yout                                                        
          Y(nY) = dsqrt(Y(nY)*Y(nY-1))                                          
          nY = nY + 1                                                           
        ELSE                                                                    
          istop = 1                                                             
        END IF                                                                  
      END DO                                                                    
c     check dynamical range of Eta to avoid nonphysical results or code errors [
      Etamax = 0.                                                               
      Etamin = 1.e+20                                                           
      DO iY = 1, nY                                                             
       IF(ETA(Y(iY)).lt.Etamin) Etamin = ETA(Y(iY))                             
       IF(ETA(Y(iY)).gt.Etamax) Etamax = ETA(Y(iY))                             
       IF (ETA(Y(iY)).LT.1.e-12) THEN                                           
        IF (iX.GT.0) THEN                                                       
         write(18,*)'      Y          ETA  '                                    
         DO jY = 1, nY                                                          
           write(18,'(1p,2e12.3)') Y(jY),ETA(Y(jY))                             
         END DO                                                                 
        END IF                                                                  
        CALL MSG(17)                                                            
        error = 6                                                               
        iERROR = iERROR + 1                                                     
        goto 102                                                                
       END IF                                                                   
      END DO                                                                    
      IF ((Etamin/Etamax).LT.1.e-12) THEN                                       
       CALL MSG(18)                                                             
       error = 6                                                                
       iERROR = iERROR + 1                                                      
       goto 102                                                                 
      END IF                                                                    
c     check that the Y grid is not too large                                    
      IF (nY.GT.npY) THEN                                                       
        delTAUsc = delTAUsc * 1.5                                               
        iWARNING = iWARNING + 1                                                 
        IF (iX.NE.0) THEN                                                       
        IF (iter.EQ.1) call line(0,2,18)                                        
        write(18,'(a)')                                                         
     &  ' ****************   WARNING   *******************'                     
        write(18,'(a46,i3)')                                                    
     &             ' Initial delTAUsc resulted in too many points:',nY          
        write(18,'(a,i3,a)')' You need to increase npY in userpar.inc'          
        write(18,*)'Multiplying delTAUsc by 1.5 and trying again'               
        write(18,'(a14,1p,e10.3)')' New delTAUsc:',delTAUsc                     
        END IF                                                                  
        IF (iter.LT.5) THEN                                                     
          goto 101                                                              
        ELSE                                                                    
          IF (iX.NE.0) THEN                                                     
            write(18,'(a)')                                                     
     &            ' ****************  GIVING UP  *******************'           
            call line(0,2,18)                                                   
          END IF                                                                
          error = 2                                                             
          goto 102                                                              
        END IF                                                                  
      END IF                                                                    
c     intepolate ETAdiscr to new Y grid (for RDW (denstyp=5 or 6))              
c      write(18,*)' ***** From Ygrid *****'                                     
c      write(18,*)'      Y      ETAdiscr      '                                 
      DO iY = 1, nY                                                             
        Yloc = Y(iY)                                                            
        IF (iterETA.GT.1) THEN                                                  
          CALL LinInter(npY,nYprev,Yprev,EtaTemp,Yloc,iYdummy,ee)               
          ETAdiscr(iY) = ee                                                     
        ELSE                                                                    
          ETAdiscr(iY) = ETA(Yloc)                                              
        END IF                                                                  
c        write(18,'(1p,2e12.4)')Y(iY), ETAdiscr(iY)                             
      END DO                                                                    
c      write(18,*)' *************************'                                  
c -----------------------------------------------------------------------       
102   RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c =======================================================================       
c     These are the subroutines for a full dynamical calculation of             
c     radiatively driven winds.                              [MN, Mar'99]       
c =======================================================================       
C     Table of Contents                                                         
C                                                                               
C     CALCETA                                                                   
C     DYNAMICS                                                                  
C     EMFUNC                                                                    
C     GAMMAFUN                                                                  
C     UFUN                                                                      
C     VRATFUN                                                                   
C     WINDS                                                                     
c =======================================================================       
                                                                                
c ***********************************************************************       
      SUBROUTINE CalcETA(Y,qF,u,vrat,Eta,tauF,nY)                               
c =======================================================================       
c Calculates the dimensionless density profile ETA(y) and tauF(y) when          
c taking the drift into account                         [ZI & MN, Aug'96]       
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER iY, nY                                                            
      DOUBLE PRECISION qF(npY), Y(npY), u(npY), Eta(npY), tauF(npY),            
     &       F1(npY), F2(npY), vrat(npY), EtaINT, INT                           
c -----------------------------------------------------------------------       
c     generate ETA and its integral (normalization constant)                    
      DO iY = 1, nY                                                             
        F1(iY) = vrat(iY)/u(iY)/Y(iY)/Y(iY)                                     
      END DO                                                                    
      CALL SIMPSON(npY,1,nY,Y,F1,EtaINT)                                        
c     find tauF                                                                 
      DO iY = 1, nY                                                             
        Eta(iY) = F1(iY)/EtaINT                                                 
        F2(iY) = qF(iY)*ETA(iY)                                                 
        CALL SIMPSON(npY,1,iY,Y,F2,INT)                                         
        tauF(iY) = TAUfid*INT                                                   
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE DYNAMICS(nG,ETAzp,ftot)                                        
c =======================================================================       
c This subroutine finds the velocity structure of a radiatively driven          
c wind.                                                                         
c *** This version works for single size grains only ***                        
c                                                       [ZI & MN, Aug'96]       
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      DOUBLE PRECISION ugas(npY), qF(npY), vrat(npG,npY), Gamma(npY),           
     &       I1, I2, I3, CMdot, Cve, CM, Cr1                                    
      COMMON /dyn/ ugas, qF, vrat, Gamma, I1, I2, I3, CMdot, Cve, CM,           
     &       Cr1                                                                
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER nG, iL, iY, iconv, itr, ETAconv, uconv                            
      DOUBLE PRECISION ETAold(npY), u(npY), uold(npY), resaux,tauF(npY),        
     &       Faux(npL), u1, eps, uacc, GammaMax, ETAzp(npP,npY),                
     &       ftot(npL,npY),Sigaux                                               
c -----------------------------------------------------------------------       
c     Temporary local rescaling (to reconcile with Zeljko's f-lae which         
c     were for <sigma>/<V>) [MN]                                                
      Sigaux = SigExfid/aveV                                                    
      IF (iX.GE.1) THEN                                                         
         write(18,*)' Doing Dynamics'                                           
      END IF                                                                    
      IF(iVerb.EQ.2) write(*,*)' Doing Dynamics'                                
c     so far it works for nG=1 only:                                            
      IF (nG.GT.1) THEN                                                         
        write(12,*)' **************************** '                             
        write(12,*)' Change dynamics sub to nG>1! '                             
        write(12,*)'       PROGRAM STOPPED        '                             
        write(12,*)' **************************** '                             
        stop                                                                    
      END IF                                                                    
c     accuracy for velocity convergence same as for Utot:                       
      uacc = accConv                                                            
c     find the flux averaged extinction efficiency, qF(y), and flux             
c     averaged optical depth, tauF(y) (with old ETA, i.e. it can be             
c     determined by using ETAzp, see below)                                     
      DO iY = 1, nY                                                             
c       generate auxiliary function for lambda integration:                     
        DO iL = 1, nL                                                           
          Faux(iL) = (SigmaA(1,iL)+SigmaS(1,iL))*ftot(iL,iY)/lambda(iL)         
        END DO                                                                  
        CALL Simpson(npL,1,nL,lambda,Faux,resaux)                               
c       qF scaled with Qext at lamfid                                           
        qF(iY) = resaux / SigExfid                                              
c       tauF                                                                    
        IF (iY.EQ.1) THEN                                                       
          tauF(iY) = 0.0                                                        
        ELSE                                                                    
          DO iL = 1, nL                                                         
           Faux(iL)=ETAzp(1,iY)*TAUtot(iL)*ftot(iL,iY)/lambda(iL)               
          END DO                                                                
          CALL Simpson(npL,1,nL,lambda,Faux,resaux)                             
          tauF(iY) = resaux                                                     
        END IF                                                                  
      END DO                                                                    
c     assign input parameters to local variables                                
      GammaMax = ptr(1)                                                         
      IF (denstyp.EQ.5) THEN                                                    
        eps = pow                                                               
        u1 =  tauF(nY)*eps*(1.-0.5*GammaMax)/(1.-eps)                           
        ELSE                                                                    
        u1 = pow                                                                
      END IF                                                                    
c     Initial approximation for u(y)                                            
      DO iY = 1, nY                                                             
        uold(iY) = u1 + tauF(iY)*(1.-0.5*GammaMax)                              
      END DO                                                                    
c     Initial approximation for v/vd                                            
      CALL vratFun(Sigaux,qF,uold,vrat,nY)                                      
c     Initial approximation for Gamma                                           
      CALL GammaFun(GammaMax,qF,uold,vrat,Y,Gamma,nY,I1,I2,I3)                  
      iconv = 0                                                                 
      itr = 0                                                                   
c     ITERATIONS until u and eta converge within uacc                           
      DO WHILE(iconv.NE.1)                                                      
        itr = itr + 1                                                           
c       Find the new value of u(y)                                              
        IF (denstyp.EQ.5)                                                       
     &      u1 = tauF(nY)*eps*(1.-Gamma(nY))/(1.-eps)                           
        CALL uFun(u1,tauF,Gamma,u,nY)                                           
c       Find the new ETAdiscr(y) and tauF(y)                                    
        CALL CalcETA(Y,qF,u,vrat,ETAdiscr,tauF,nY)                              
c       generate new v/vd                                                       
        CALL vratFun(Sigaux,qF,u,vrat,nY)                                       
c       generate new Gamma                                                      
        CALL GammaFun(GammaMax,qF,u,vrat,Y,Gamma,nY,I1,I2,I3)                   
c       check convergence of u and Eta                                          
        IF (itr.GT.1) THEN                                                      
          CALL ChkConv(nY,uacc,uold,u,uconv)                                    
          CALL ChkConv(nY,uacc,ETAold,ETAdiscr,ETAconv)                         
c         convergence required for both u(y) and ETA(y)                         
          iconv = ETAconv * uconv                                               
        END IF                                                                  
        IF (iconv.NE.1) THEN                                                    
          DO iY =1, nY                                                          
            uold(iY) = u(iY)                                                    
            ETAold(iY) = ETAdiscr(iY)                                           
          END DO                                                                
          IF (itr.GE.100) iconv = 1                                             
        ELSE                                                                    
          DO iY = 1, nY                                                         
            ugas(iY) = u(iY)                                                    
          END DO                                                                
        END IF                                                                  
      END DO                                                                    
      IF (iX.GE.1)                                                              
     &  write(18,'(a35,i3)')' Number of iterations to converge:',itr            
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      DOUBLE PRECISION FUNCTION EMfunc(lambda,Teff,xSiO)                        
c This is modeled after subroutine engelke by M. Marengo. Here are his          
c original comments:                                                            
c     =================================================================         
c     This subroutine computes a modified black body spectrum using an          
c     "Engelke" function (see Engelke 1992, AJ 104, 1248):                      
c     Bnu = Bnu(Tb) with Tb = 0.738*Teff*(1+79450/(lambda*Teff))**0.182         
c                                                                               
c     Molecular SiO absorption is modelled from the alpha Tau spectrum          
c     of Cohen et al. 1992, AJ 104,2030 with a 5th order polinomial,            
c     and added to the modified bb.                                             
c                                                                               
c     M. Marengo - mmarengo@cfa.harvard.edu - Sep 1998                          
c     =================================================================         
c                                                                               
c                                                                               
c This version makes use of the scaled quantities and Dusty's function          
c Planck(x)                                                  [ZI, Feb 99]       
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER j                                                                 
      DOUBLE PRECISION lambda, Teff, xSiO, x, Planck, Teng, SiOc(6),            
     &                 lambda1, lambda2, SiO8m, SiOf                            
c -----------------------------------------------------------------------       
c  SiO fit data from Massimo:                                                   
c     Polinomial coeff for SiO absorption model (5th order),                    
c     wavelength interval in which to apply the absorption                      
c     and given absorption at 8 micron (to rescale for SiO)                     
      lambda1 =  7.8636                                                         
      lambda2 = 11.4280                                                         
      SiO8m = 1.0701447                                                         
      SiOc(1) = -300.43916                                                      
      SiOc(2) =  149.32134                                                      
      SiOc(3) =  -29.493280                                                     
      SiOc(4) =    2.9067144                                                    
      SiOc(5) =   -0.14304663                                                   
      SiOc(6) =    0.0028134070                                                 
c -----------------------------------------------------------------------       
c     Engelke's effective temperature                                           
      Teng = 0.738 * Teff * (1.0 + 79450.0/(lambda*Teff))**0.182                
      x = 14400.0 / lambda / Teng                                               
      EMfunc = (Teng/Teff)**4 * Planck(x)                                       
c     If lambda is in SiO region, compute and apply the SiO absorption          
      IF ((lambda-lambda1)*(lambda-lambda2).LT.0) THEN                          
          SiOf = 0.0                                                            
          DO j = 1, 6                                                           
             SiOf = SiOf + SiOc(j) * lambda**(1.0*j-1)                          
          END DO                                                                
          EMfunc = EMfunc / (1.0 + (SiOf-1)/(SiO8m-1)*xSiO*0.01)                
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE  GammaFun(GammaMax,qF,u,vrat,Y,Gamma,nY,Int1,Int2,Int3)        
c =======================================================================       
c This subroutine finds the ratio of the gravitaional to the radiative          
c force for the whole envelope Gamma^(-1). Note that this quantity is           
c called Gamma in the code.                                  [MN, Aug'96]       
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER jY, nY                                                            
      DOUBLE PRECISION qF(npY), u(npY), vrat(npY), Gamma(npY), Y(npY),          
     &       Int2, Int3, Faux1(npY), Faux2(npY), Faux3(npY), GammaMax,          
     &       C, Int1                                                            
c -----------------------------------------------------------------------       
c     generate integrals                                                        
      DO jY = 1, nY                                                             
        Faux1(jY) = 1.0/u(jY)/Y(jY)/Y(jY)                                       
        Faux2(jY) = vrat(jY)*Faux1(jY)                                          
        Faux3(jY) = qF(jY)*Faux2(jY)                                            
        CALL SIMPSON(npY,1,jY,Y,Faux1,Int1)                                     
        CALL SIMPSON(npY,1,jY,Y,Faux2,Int2)                                     
        CALL SIMPSON(npY,1,jY,Y,Faux3,Int3)                                     
        IF (jY.GT.1) THEN                                                       
          Gamma(jY) = Int1 / Int3                                               
          ELSE                                                                  
          Gamma(jY) = 1.0 / vrat(1) / qF(1)                                     
        END IF                                                                  
      END DO                                                                    
c     find normalization constants                                              
      CALL FindMax(npY,1,nY,Gamma,C)                                            
      C = GammaMax / C                                                          
c     normalize Gamma                                                           
      DO jY = 1, nY                                                             
        Gamma(jY) = C * Gamma(jY)                                               
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE uFun(u1,tauF,Gamma,u,nY)                                       
c =======================================================================       
c Calculates the scaled gas velocity u(y) with given u1=u(1), tauF(y),          
c and gravitational correction term Gamma^(-1) (named Gamma in the code).       
c                                                            [MN, Aug'96]       
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY,iY                                                             
      DOUBLE PRECISION u(npY), Gamma(npY), tauF(npY), u1                        
c -----------------------------------------------------------------------       
      DO iY = 1, nY                                                             
         u(iY) = u1 + tauF(iY) * (1.0-Gamma(iY))                                
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE vratFun(Qfid,qF,u,vrat,nY)                                     
c =======================================================================       
c Calculates the drift correction v(y)/vd(y).                [MN, Aug'96]       
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, iY                                                            
      DOUBLE PRECISION Qfid, qF(npY), u(npY), vrat(npG,npY)                     
c -----------------------------------------------------------------------       
      DO iY = 1, nY                                                             
         vrat(1,iY) = 1.0 / (1.0 + dsqrt(Qfid*qF(iY)/u(iY)))                    
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE WINDS(nG,EtaOK,ETAzp,ftot)                                     
c =======================================================================       
c This subroutine takes care of the interface between radiatively driven        
c winds and radiative transfer.                              [ZI, Aug'95]       
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER nG, EtaOK, iY                                                     
      DOUBLE PRECISION ETAold(npY), ETAzp(npP,npY), ftot(npL,npY),              
     &       accETA                                                             
c -----------------------------------------------------------------------       
      DO iY = 1, nY                                                             
        ETAold(iY) = ETAdiscr(iY)                                               
      END DO                                                                    
c     find new ETA (ETAdiscr in density.inc) based on current flux ftot         
      CALL DYNAMICS(nG,ETAzp,ftot)                                              
c     check convergence (ptr(2) is specified in INPUT)                          
      accETA = ptr(2) * accuracy                                                
      CALL ChkConv(nY,accETA,ETAold,ETAdiscr,EtaOK)                             
      IF (iX.GE.1) THEN                                                         
        write(18,*)'     Y         ETAold      ETAnew      ratio'               
        DO iY = 1, nY                                                           
          accETA = ETAold(iY) / ETAdiscr(iY)                                    
          write(18,'(1p,4e12.5)')Y(iY),ETAold(iY),ETAdiscr(iY),accETA           
        END DO                                                                  
        IF (EtaOK.EQ.1) THEN                                                    
          write(18,*)' Convergence on Eta achieved'                             
        ELSE                                                                    
          write(18,*)' Convergence on Eta not achieved.'                        
          write(18,*)' Going to the next iteration.'                            
        END IF                                                                  
      END IF                                                                    
c     save Y to Yprev and nY to nYprev                                          
      DO iY = 1, nY                                                             
         Yprev(iY) = Y(iY)                                                      
      END DO                                                                    
      nYprev = nY                                                               
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
                                                                                
c =========================================================================     
c          Below are the subroutines for calculation in slab geometry           
c                     arranged in alphabetical order                            
c =========================================================================     
C     Table of Contents                                                         
C                                                                               
C     SLBANALYT                                                                 
C     SLBDIFF                                                                   
C     EINT1                                                                     
C     EINT2                                                                     
C     EINT3                                                                     
C     EINT4                                                                     
C     EINT5                                                                     
C     SLBACC                                                                    
C     SLBGRAY                                                                   
C     SLBINIT                                                                   
C     SLBMAT                                                                    
C     SLBMISC                                                                   
C     SLBRADT                                                                   
C     SLBSOLVE                                                                  
C     SLBSTAR                                                                   
C     SLBTAU                                                                    
C     SLBTRACE                                                                  
C     SLBY                                                                      
c =========================================================================     
                                                                                
c *************************************************************************     
      SUBROUTINE SLBAnalyt(TAU,L1,L2,L3,L4,N,NN)                                
c =========================================================================     
c     Finds the integrals of the moments of the first exponential function E1   
c     Lk(Tau)=Int{t^(k-1) E1|t-Tau|}, k=1..4. The analytical expressions used   
c     are obtained integrating Lk by parts.                     [MN,Nov'98]     
c =========================================================================     
      IMPLICIT NONE                                                             
      INTEGER i, j, N, NN                                                       
      DOUBLE PRECISION L1(NN,NN), L2(NN,NN), L3(NN,NN), L4(NN,NN), del,         
     &          TAU(NN), arg, arg1, Eint2, Eint3, Eint4, Eint5, Lt1,            
     &          Lt2, Lt3, Lt4, E2a1, E3a1, E4a1                                 
c ----------------------------------------------------------------------        
      DO i = 1, N                                                               
       DO j = 1, N                                                              
        L1(i,j) = 0.                                                            
        L2(i,j) = 0.                                                            
        L3(i,j) = 0.                                                            
        L4(i,j) = 0.                                                            
       END DO                                                                   
      END DO                                                                    
c     Calculate Lk(i)=Int_t^t1 {x**(k-1)E1|x-TAU(i)| dx}                        
      DO j = 1, N                                                               
        DO i = 1, N-1                                                           
          del = TAU(i+1)-TAU(i)                                                 
          arg1 = dabs(TAU(j)-TAU(i+1))                                          
          arg = dabs(TAU(j)-TAU(i))                                             
                                                                                
          E2a1 = Eint2(arg1)                                                    
          E3a1 = Eint3(arg1)                                                    
          E4a1 = Eint4(arg1)                                                    
          Lt1 = E2a1-Eint2(arg)                                                 
          Lt2 =-(E3a1-Eint3(arg))/del                                           
          Lt3 = E2a1 + 2./del/del*(E4a1-Eint4(arg))                             
          Lt4 =-3./del*E3a1 -6./del/del/del*                                    
     *                     (Eint5(arg1)-Eint5(arg))                             
          IF(i.LT.j) THEN                                                       
             L1(i,j) = 0.5*Lt1                                                  
             L2(i,j) = 0.5*(Lt2 + E2a1)                                         
             L3(i,j) = 0.5*(Lt3 - 2./del*E3a1)                                  
             L4(i,j) = 0.5*(Lt4 + E2a1 + 6./del/del*E4a1)                       
           ELSE                                                                 
            L1(i,j) = -0.5*Lt1                                                  
            L2(i,j) = 0.5*(Lt2 - E2a1)                                          
            L3(i,j) = 0.5*(-Lt3 - 2./del*E3a1 )                                 
            L4(i,j) = 0.5*(Lt4 - E2a1 - 6./del/del*E4a1)                        
          END IF                                                                
        END DO                                                                  
      END DO                                                                    
                                                                                
      DO i =1 , N                                                               
        DO j = 1, N                                                             
          IF(L1(i,j).LT.1.e-20) L1(i,j) = 0.                                    
          IF(L2(i,j).LT.1.e-20) L2(i,j) = 0.                                    
          IF(L3(i,j).LT.1.e-20) L3(i,j) = 0.                                    
          IF(L4(i,j).LT.1.e-20) L4(i,j) = 0.                                    
        END DO                                                                  
      END DO                                                                    
c ---------------------------------------------------------------------         
      RETURN                                                                    
      END                                                                       
c *********************************************************************         
                                                                                
                                                                                
c *********************************************************************         
      SUBROUTINE SLBdifF(flag,om,grid,mat1,nL,nY,mat2,fp,fm)                    
c =========================================================================     
c     Integration of U(lam,t)*E2|Tau-t| to get the diffuse flux.                
c     flag=1 is for scattered and flag=0 for emitted flux. [MN, Apr'98]         
c =========================================================================     
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER iL, iY, j, nL, nY, flag                                           
      DOUBLE PRECISION mat1(npL,npY), mat2(npL,npY), grid(npL,npY),             
     &             om(npL,npY), fp(npL,npY), fm(npL,npY), TAU(npY),             
     &             faux(npY), Efact, fave, SUM, Eint3                           
c --------------------------------------------------------------------          
      DO iL = 1, nL                                                             
        DO iY = 1, nY                                                           
          TAU(iY) = grid(iL,iY)                                                 
          IF (flag.EQ.1) THEN                                                   
             faux(iY) = om(iL,iY)*mat1(iL,iY)                                   
          ELSE                                                                  
             faux(iY) = (1.0-om(iL,iY))*mat1(iL,iY)                             
          END IF                                                                
        END DO                                                                  
                                                                                
c       Find f(+) (in arg Tau>t)                                                
        DO iY = 1, nY                                                           
          SUM = 0.                                                              
          DO j = 1, iY-1                                                        
           Efact = dabs(Eint3(TAU(iY)-TAU(j))-Eint3(TAU(iY)-TAU(j+1)))          
           fave = 0.5*(faux(j)+faux(j+1))                                       
           SUM = SUM + fave*Efact                                               
          END DO                                                                
          fp(iL,iY) = 0.5*SUM                                                   
c       and f(-) (in arg Tau>t)                                                 
          SUM = 0.                                                              
          DO j = iY, nY-1                                                       
           Efact = dabs(Eint3(TAU(iY)-TAU(j))-Eint3(TAU(iY)-TAU(j+1)))          
           fave = 0.5*(faux(j)+faux(j+1))                                       
           SUM = SUM + fave*Efact                                               
          END DO                                                                
          fm(iL,iY) = 0.5*SUM                                                   
        END DO                                                                  
                                                                                
        DO iY = 1, nY                                                           
           mat2(iL,iY) = fp(iL,iY) - fm(iL,iY)                                  
        END DO                                                                  
c     End of loop over iL                                                       
      END DO                                                                    
c --------------------------------------------------------------------          
      RETURN                                                                    
      END                                                                       
c *********************************************************************         
                                                                                
                                                                                
c **********************************************************************        
      DOUBLE PRECISION FUNCTION Eint1(x)                                        
c ======================================================================        
c Needed for the slab geometry. It calculates the first exponential             
c integral E1(x) by analytical f-la (13.13) from Abramovitz & Stegun(1994)      
c                                                         [MN,Dec'97]           
c ======================================================================        
      IMPLICIT none                                                             
      INTEGER i                                                                 
      DOUBLE PRECISION AC(4),BC(4), CC(6), x, aux, poly, denom                  
      DATA AC/8.5733287401,18.0590169730,8.6347608925,0.2677737343/             
      DATA BC/9.5733223454,25.6329561486,21.0996530827,3.9584969228/            
      DATA CC/-0.57721566,0.99999193,-0.24991055,0.05519968,-0.00976004,        
     &        0.00107857/                                                       
c ----------------------------------------------------------------------        
c  For x=1E-15, E1~30 (used below to limit the value at x=0);for x>1, E1<1e-8   
c  Two approximations are used, for x>1 and x<1, respectively                   
      IF (x.GT.1.0) THEN                                                        
         poly = 0.0                                                             
         denom = 0.0                                                            
         aux = 1.0                                                              
         DO i = 1, 4                                                            
           poly = poly + AC(5-i)*aux                                            
           denom = denom + BC(5-i)*aux                                          
           aux = aux * x                                                        
         END DO                                                                 
           poly = poly + aux                                                    
           denom = denom + aux                                                  
         Eint1 = poly/denom/x*dexp(-x)                                          
      ELSE                                                                      
         IF (x.LE.1.E-15) x=1.0E-15                                             
         poly = 0.0                                                             
         aux = 1.0                                                              
         DO i = 1, 6                                                            
           poly = poly + CC(i)*aux                                              
           aux = aux * x                                                        
         END DO                                                                 
         Eint1 = poly - dlog(x)                                                 
      END IF                                                                    
c ----------------------------------------------------------------------        
      RETURN                                                                    
      END                                                                       
c **********************************************************************        
                                                                                
c **********************************************************************        
      DOUBLE PRECISION FUNCTION Eint2(x)                                        
c ======================================================================        
c Needed for the slab geometry. It calculates the second exponential            
c integral E2(x) by the recurrence f-la. (see Abramovitz & Stegun,1994)         
c                                                         [MN,Dec'97]           
c ======================================================================        
      IMPLICIT none                                                             
      DOUBLE PRECISION x, Eint1                                                 
c -------------------------------------------------------------------------     
       IF(x.LT.0.) x=dabs(x)                                                    
       IF (x.LT.1.0D-15) THEN                                                   
         Eint2 = 1.0                                                            
        ELSE                                                                    
         Eint2 = dexp(-x) - x*Eint1(x)                                          
       END IF                                                                   
c -------------------------------------------------------------------------     
      RETURN                                                                    
      END                                                                       
c **********************************************************************        
                                                                                
c **********************************************************************        
      DOUBLE PRECISION FUNCTION Eint3(x)                                        
c ======================================================================        
c Needed for the slab geometry. It calculates the third exponential             
c integral E3(x) by the recurrence f-la. (see Abramovitz & Stegun,1994)         
c                                                        [MN,Dec'97]            
c ======================================================================        
      IMPLICIT none                                                             
      DOUBLE PRECISION x, Eint1                                                 
c -------------------------------------------------------------------------     
       IF(x.LT.0.) x=dabs(x)                                                    
       IF (x.LT.1.0D-15) THEN                                                   
         Eint3 = 0.5                                                            
        ELSE                                                                    
         Eint3 = 0.5*((1.0-x)*dexp(-x)+x*x*Eint1(x))                            
       END IF                                                                   
c -------------------------------------------------------------------------     
      RETURN                                                                    
      END                                                                       
c **********************************************************************        
                                                                                
                                                                                
c **********************************************************************        
      DOUBLE PRECISION FUNCTION Eint4(x)                                        
c ======================================================================        
c Needed for the slab geometry. It calculates the fourth exponential            
c integral E4(x) by the recurrence f-la. (see Abramovitz & Stegun,1994)         
c                                                         [MN,Jan'97]           
c ======================================================================        
      IMPLICIT none                                                             
      DOUBLE PRECISION x, Eint3                                                 
c -------------------------------------------------------------------------     
       IF(x.LT.0.) x=dabs(x)                                                    
       IF (x.LT.1.0D-15) THEN                                                   
         Eint4 = 1.0d00/3.0d00                                                  
        ELSE                                                                    
         Eint4 = (dexp(-x)-x*Eint3(x))/3.0d00                                   
       END IF                                                                   
c -------------------------------------------------------------------------     
      RETURN                                                                    
      END                                                                       
c **********************************************************************        
                                                                                
                                                                                
c **********************************************************************        
      DOUBLE PRECISION FUNCTION Eint5(x)                                        
c ======================================================================        
c Needed for the slab geometry. It calculates the fifth exponential             
c integral E5(x) by the recurrence f-la. (see Abramovitz & Stegun,1994)         
c                                                         [MN,Jan'97]           
c ======================================================================        
      IMPLICIT none                                                             
      DOUBLE PRECISION x, Eint4                                                 
c -------------------------------------------------------------------------     
        IF(x.LT.0.) x=dabs(x)                                                   
       IF (x.LT.1.0d-15) THEN                                                   
         Eint5 = 0.25d00                                                        
        ELSE                                                                    
         Eint5 = 0.25d00*(dexp(-x)-x*Eint4(x))                                  
       END IF                                                                   
c -------------------------------------------------------------------------     
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE SLBacc(flux,accuracy,devmax,FbolOK,error)                      
c =======================================================================       
c This is SUB ChkFlux modified to be used for slab calculation.                 
c Checks the bolometric flux conservation at any point of the slab grid.        
c In case of nonconservation inserts a number of points at certain              
c places.                                             [MN,99; ZI'96]            
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      DOUBLE PRECISION TAUslb(npL,npY),fsbol(npY),fpbol(npY),fmbol(npY),        
     &                 AveDev, RMS, maxFerr, fmax, fmed                         
      COMMON /slab/ TAUslb, fsbol, fpbol, fmbol, AveDev, RMS, maxFerr,          
     &                 fmax, fmed                                               
      INTEGER iYins(npY), k, kins, iY, flag, error, istop, FbolOK               
      DOUBLE PRECISION Yins(npY), flux(npY), delTAUmax, devmax, ratio,          
     &        accuracy, devfac, ff, ffold                                       
c --------------------------------------------------------------------------    
      error = 0                                                                 
      kins = 0                                                                  
      devmax = 0.0                                                              
c     maximal delTAU is no more than 2 times the average value                  
      delTAUmax = 2.0*TAUtot(1)/nY                                              
c     find maximal relative deviation of fbol:                                  
      DO iY = 2, nY                                                             
        ratio = (dabs(flux(iY))-dabs(fmed))/dabs(fmed)                          
        IF (dabs(ratio).GT.devmax) devmax = dabs(ratio)                         
      END DO                                                                    
      ff = 0.0                                                                  
      istop = 0                                                                 
      devfac = 0.1                                                              
c     search for places to improve the grid                                     
      DO WHILE (istop.NE.1)                                                     
        DO iY = 2, nY                                                           
          ratio = (dabs(flux(iY))-dabs(fmed))/dabs(fmed)                        
          ffold = abs(flux(iY-1)/fmed-1.0)                                      
          ff = abs(flux(iY)/fmed-1.0)                                           
          flag = 0                                                              
c         if any of these criteria is satisfied insert a point:                 
          IF(ff.GT.accuracy) flag=1                                             
c         1) if error is increasing too fast                                    
          IF (abs(ff-ffold).GT.devfac*devmax) flag = 1                          
c         2) if delTAU is too large at the left edge:                           
          IF(TAUslb(1,iY).LT.5.) THEN                                           
           IF ((TAUslb(1,iY)-TAUslb(1,iY-1)).GT.                                
     &                             delTAUmax) flag = 1                          
          END IF                                                                
          IF(flag.EQ.1.AND.devmax.GE.accuracy) THEN                             
            kins = kins + 1                                                     
            Yins(kins) = Y(iY-1)+0.5*(Y(iY)-Y(iY-1))                            
            iYins(kins) = iY-1                                                  
          END IF                                                                
        END DO                                                                  
        IF (devmax.LT.accuracy.OR.devfac.LT.0.01) THEN                          
          istop = 1                                                             
        ELSE                                                                    
          IF (kins.GT.0) istop = 1                                              
        END IF                                                                  
        devfac = devfac / 2.0                                                   
      END DO                                                                    
                                                                                
      IF (kins.EQ.0) THEN                                                       
         FbolOK = 1                                                             
        ELSE                                                                    
c       Add all new points to Y(nY). This gives the new Y(nY+kins).             
c       However, check if npY is large enough to insert all points:             
        IF ((nY+kins).GT.npY) THEN                                              
         IF (iX.GE.1) THEN                                                      
         write(18,*)' ****************     WARNING   ******************'        
         write(18,*)'  The new Y grid can not accomodate more points!'          
         write(18,'(a,i5)')'   Specified accuracy would require',nY+kins        
         write(18,'(a,i5,a)')'   points, while npY =',npY,'.'                   
         write(18,*)'  For the required accuracy npY must be increased,'        
         write(18,*)'  (see the manual S3.5 Numerical Accuracy).'               
         write(18,*)' *************************************************'        
         END IF                                                                 
         kins = npY - nY                                                        
         iWARNING = iWARNING + 1                                                
         error = 2                                                              
        END IF                                                                  
        DO k = 1, kins                                                          
          CALL SHIFT(Y,npY,nY+k-1,Yins(k),iYins(k)+k-1)                         
        END DO                                                                  
      END IF                                                                    
c     new size of the Y grid                                                    
      nY = nY + kins                                                            
c --------------------------------------------------------------------------    
777   RETURN                                                                    
      END                                                                       
c *********************************************************************         
                                                                                
c *********************************************************************         
      SUBROUTINE SLBmisc(fbol,fmax,fmed,AveDev,RMS,nY)                          
c =======================================================================       
c     Finds some additional quantities: fave-the average and fmed- the          
c     median value of fbol in the slab, average and RMS deviations of           
c     fbol from fmed.                                     [MN,Aug'98]           
c =======================================================================       
      IMPLICIT NONE                                                             
      INTEGER iY, nY, imid, iup, idn                                            
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      DOUBLE PRECISION fbol(npY), fsort(npY), fave, fmed, AveDev,               
     &                 err, RMS, rnY, fmax, aux                                 
c     Find the median of fbol values distribution                               
      DO iY = 1, nY                                                             
        aux = DMAX1(dabs(fbol(iY)), accFbol)                                    
        fbol(iY) = aux                                                          
      END DO                                                                    
      DO iY = 1, nY                                                             
        fsort(iY) = fbol(iY)                                                    
      END DO                                                                    
      CALL sort(fsort,nY)                                                       
      IF (MOD(nY,2).NE.0) THEN                                                  
c       if nY odd:                                                              
        imid = (nY+1)/2                                                         
        fmed = fsort(imid)                                                      
      ELSE                                                                      
c       if nY even:                                                             
        iup = nY/2 + 1                                                          
        idn = nY/2                                                              
        fmed = 0.5*(fsort(iup) + fsort(idn))                                    
      END IF                                                                    
c     In sub PrOut: if(fmed.LE.accFbol) fmed =0                                 
      rnY = 1.0*nY                                                              
      fave = 0.                                                                 
      fmax = 0.                                                                 
c     find the max |fbol| (carried in slab.inc, needed in SLBacc)               
      DO iY = 1, nY                                                             
        fave = fave + fbol(iY)                                                  
        aux = dabs(fbol(iY))                                                    
        IF(aux.GT.fmax) fmax = aux                                              
      END DO                                                                    
      fave = fave / rnY                                                         
c     AveDev is the average relative deviation from the median                  
      AveDev = 0.                                                               
      RMS = 0.0                                                                 
      DO iY = 1, nY                                                             
        aux = dabs(fbol(iY))                                                    
c        err = dabs(fmed-aux) /fmax                                             
        err = dabs(fmed-aux) /fmed                                              
        AveDev = AveDev + err                                                   
        RMS = RMS + err*err                                                     
      END DO                                                                    
      RMS = sqrt(RMS/rnY/(rnY-1.))                                              
      AveDev = AveDev/rnY                                                       
c --------------------------------------------------------------------          
      RETURN                                                                    
      END                                                                       
c ********************************************************************          
                                                                                
c ********************************************************************          
      SUBROUTINE SLBgray(model,nG,error)                                        
c =======================================================================       
c     Solves the gray slab problem. For single lambda calculation               
c     comment the call to Spectral in main.for and the calls to FindErr         
c     in Analysis.for                                      [MN, May'98]         
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      DOUBLE PRECISION TAUslb(npL,npY),fsbol(npY),fpbol(npY),fmbol(npY),        
     &                 AveDev, RMS, maxFerr, fmax, fmed                         
      COMMON /slab/ TAUslb, fsbol, fpbol, fmbol, AveDev, RMS, maxFerr,          
     &                 fmax, fmed                                               
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER nYok, nPok                                                        
      DOUBLE PRECISION                                                          
     &       Ude(npL,npY), Uds(npL,npY), Us(npL,npY), fde(npL,npY),             
     &       fds(npL,npY), fs(npL,npY), Utot(npL,npY), ftot(npL,npY),           
     &       Td(npG,npY), Ubol(npY), fbol(npY), Uchck(npL,npY),                 
     &       UbolChck(npY), Spectrum(npL), SpecChar(30,99), tauF(npY),          
     &       Intens(npL,npP+2), IntOut(20,npP+2), SmC(30,99),                   
     &       Yok(npY), Pok(npP), tauOut(npL), tauZout(npP+2), F1,               
     &       tr(npY), rg(npG,npY), fsL(npL,npY), fsR(npL,npY), Eps(npY)         
      COMMON /solution/ Ude, Uds, Us, fde, fds, fs, Utot, ftot, Td,             
     &       Ubol, fbol, Uchck, UbolChck, Spectrum, SpecChar, tauF,             
     &       Intens, IntOut, SmC, Yok, Pok, tauOut, tauZout, F1,                
     &       tr, rg, fsL, fsR, Eps, nYok, nPok                                  
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER startyp(2), Nlamtr(2), nBB(2)                                     
      CHARACTER nameStar(2)*235                                                 
      DOUBLE PRECISION Tstar, lamtr(2,101), klam(2,100), Tbb(2,10),             
     &                 rellum(2,10), mu1, ksi, mu2, xSiO, r1rs                  
      COMMON /source/ Tstar, lamtr, klam, Tbb, rellum, mu1, ksi, mu2,           
     &                xSiO, r1rs, startyp, Nlamtr, nBB, nameStar                
      INTEGER iY, nG                                                            
      INTEGER error, nLst, model                                                
      character tt*58	                                                          
      DOUBLE PRECISION mat0(npL,npY,npY), Em(npL,npY), Uold(npL,npY),           
     &         fdsp(npL,npY),fdsm(npL,npY), fdep(npL,npY),fdem(npL,npY),        
     &         fp(npY), fm(npY), fd(npY)                                        
c --------------------------------------------------------------------------    
      nLst = nL                                                                 
      nL = 1                                                                    
c     generate some temporary arrays                                            
      DO iY = 1, nY                                                             
c       to prevent exp underflow:                                               
        IF(TAUslb(1,iY).GE.50.) THEN                                            
          Us(1,iY) = 0.                                                         
        ELSE                                                                    
          Us(1,iY) = dexp(-TAUslb(1,iY))                                        
        END IF                                                                  
        fs(1,iY) = Us(1,iY)                                                     
        Em(1,iY) = 0.                                                           
        omega(1,iY) = 1.                                                        
      END DO                                                                    
c     find radiative transfer matrices                                          
      IF (iX.GE.1) write(18,*)' Calculating matrix, pure scatt, iL=1'           
      CALL SLBMat(TAUslb,mat0)                                                  
c     solve for Utot                                                            
      CALL INVERT(nY,nL,mat0,Us,Em,omega,Utot,Uold,error)                       
c      find new Td                                                              
       DO iY = 1, nY                                                            
c       For a single lambda run:                                                
        Td(1,iY) = Tsub(1)*dsqrt(dsqrt(Utot(1,iY)/Utot(1,1)))                   
       END DO                                                                   
c      Find the diffuse scattered flux(fl=1 for scatt. and fl=0 is for emission)
       CALL SLBdifF(1,omega,TAUslb,Utot,nL,nY,fds,fdsp,fdsm)                    
c      Find the diffuse emitted flux                                            
       CALL SLBdifF(0,omega,TAUslb,Utot,nL,nY,fde,fdep,fdem)                    
       CALL Add(npY,nY,npL,nL,fs,fds,fde,ftot)                                  
c      overall conservation of flux - min/max err                               
       DO iY = 1, nY                                                            
        fbol(iY) = fTot(1,iY)                                                   
        fsbol(iY) = fs(1,iY)                                                    
        fp(iY) = fdsp(1,iY)+fdep(1,iY)                                          
        fm(iY) = fdsm(1,iY)+fdem(1,iY)                                          
        fd(iY) = fp(iY)-fm(iY)                                                  
       END DO                                                                   
c     FindErr calculates the err acc. to min/max values of fbol                 
      CALL FindErr(fbol,maxFerr,nY)                                             
      CALL SLBmisc(fbol,fmax,fmed,AveDev,RMS,nY)                                
c --------------------------------------------------------------                
c     The slab Tau-profile at the fiducious lambda (needed in PrOut)            
      DO iY = 1, nY                                                             
        tr(iY) = TAUslb(1,iY)                                                   
      END DO                                                                    
      IF(iInn.EQ.1) THEN                                                        
       write(18,*) '----- Single Lambda Case (iL=1) -----'                      
       write(18,'(a10,1p,E11.3)')'    tauT =', TAUslb(1,nY)                     
        write(18,'(a11,1p,E11.3)')'     fmed =',fmed                            
        write(18,'(a11,1p,E11.3)')'  RMS err =',RMS                             
        write(18,'(a11,1p,E11.3)')'  MAX err =',maxFerr                         
                                                                                
       tt = '     tr       fTot(i)      fs(i)      fd(i)     fp(i)   '          
       write(18,'(a58,a30)') tt , '  fm(i)     Utot(i)      Td(i)'              
       DO iY = 1, nY                                                            
        write(18,'(1p,8E11.3)') tr(iY), fbol(iY), fsbol(iY), fd(iY),            
     &                          fp(iY), fm(iY), Utot(1,iY), Td(1,iY)            
       END DO                                                                   
      END IF                                                                    
c     this is for running many models with single lambda case                   
      SmC(5,model) = maxFerr                                                    
      nL = nLst                                                                 
c --------------------------------------------------------------------------    
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE SLBIniT(nG,alpha)                                              
c =======================================================================       
c This subroutine calculates the initial approximation for the temperature      
c profile and alpha array. It is based on the analogous subroutine InitTemp     
c for spherical shell written by ZI, Jul'96.                [MN, Dec'97]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER nYok, nPok                                                        
      DOUBLE PRECISION                                                          
     &       Ude(npL,npY), Uds(npL,npY), Us(npL,npY), fde(npL,npY),             
     &       fds(npL,npY), fs(npL,npY), Utot(npL,npY), ftot(npL,npY),           
     &       Td(npG,npY), Ubol(npY), fbol(npY), Uchck(npL,npY),                 
     &       UbolChck(npY), Spectrum(npL), SpecChar(30,99), tauF(npY),          
     &       Intens(npL,npP+2), IntOut(20,npP+2), SmC(30,99),                   
     &       Yok(npY), Pok(npP), tauOut(npL), tauZout(npP+2), F1,               
     &       tr(npY), rg(npG,npY), fsL(npL,npY), fsR(npL,npY), Eps(npY)         
      COMMON /solution/ Ude, Uds, Us, fde, fds, fs, Utot, ftot, Td,             
     &       Ubol, fbol, Uchck, UbolChck, Spectrum, SpecChar, tauF,             
     &       Intens, IntOut, SmC, Yok, Pok, tauOut, tauZout, F1,                
     &       tr, rg, fsL, fsR, Eps, nYok, nPok                                  
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      DOUBLE PRECISION TAUslb(npL,npY),fsbol(npY),fpbol(npY),fmbol(npY),        
     &                 AveDev, RMS, maxFerr, fmax, fmed                         
      COMMON /slab/ TAUslb, fsbol, fpbol, fmbol, AveDev, RMS, maxFerr,          
     &                 fmax, fmed                                               
      INTEGER iL, iY, iG, nG, istop, iW                                         
      DOUBLE PRECISION Sigext(npL), aux(npY), faux(npL), fstar(npY),            
     &       Qfstar(npY), fstarQ(npY), QP(npY), QF(npY), oldalpha(npY),         
     &       alpha(npG,npY), IntQF, xP, Planck, resaux, delta, FovrF1,          
     &       TAU(npY)                                                           
c --------------------------------------------------------------------------    
c     Lambda integral of fs=f_e*exp(-TAU) -> fstar                              
      CALL Bolom(fs,fstar)                                                      
c     loop over grains                                                          
      DO iG = 1, nG                                                             
c       first approximation for temperature                                     
        DO iY = 1, nY                                                           
          Td(iG,iY) = Tsub(iG)                                                  
          alpha(iG,iY) = 1.0                                                    
        END DO                                                                  
c       generate the extinction cross-section Sigext                            
        DO iL = 1, nL                                                           
          Sigext(iL) = SigmaA(iG,iL)+SigmaS(iG,iL)                              
        END DO                                                                  
c       Lambda integral of Sigext*fs -> Qfstar in IE97                          
        DO iY = 1, nY                                                           
c         generate auxiliary function for lambda integration:                   
          DO iL = 1, nL                                                         
            faux(iL) = Sigext(iL) *fs(iL,iY)/lambda (iL)                        
          END DO                                                                
          CALL Simpson(npL,1,nL,lambda,faux,resaux)                             
          Qfstar(iY) = resaux                                                   
        END DO                                                                  
c      The slab geom. needs initial approximation for F/F1; Psi is alpha(nG,1)  
c      (The solution is not sensitive to this approximation, MN)                
        IF(TAUmax.LT.500.) THEN                                                 
           FovrF1 = 1.0                                                         
         ELSE                                                                   
           FovrF1 = 1. - 0.25*alpha(1,1)                                        
        END IF                                                                  
c       iterate until Td and Psi converge (i.e. until alpha does)               
        istop = 0                                                               
        DO WHILE (istop.NE.1)                                                   
          istop = 1                                                             
          DO iY = 1, nY                                                         
c           generate auxiliary function for lambda integration:                 
            DO iL = 1, nL                                                       
              xP = 14400.0 / Td(iG,iY) / lambda(iL)                             
              faux(iL) = Sigext(iL) * Planck(xP) / lambda (iL)                  
            END DO                                                              
            CALL Simpson(npL,1,nL,lambda,faux,resaux)                           
c           Planck average of Sigext                                            
            QP(iY) = resaux                                                     
c           calculate QF for slab (IV in the notes, the equivalent of eq.B5)    
            QF(iY) = Qfstar(iY)/FovrF1+QP(iY)*(1.0-fstar(iY)/FovrF1)            
c           Find the second term in III (in the notes)                          
            DO iL = 1, nL                                                       
              faux(iL) = fs(iL,iY)*(1-SigmaA(iG,iL)/QP(iY))/lambda(iL)          
            END DO                                                              
            CALL Simpson(npL,1,nL,lambda,faux,resaux)                           
            fstarQ(iY) = resaux                                                 
c           store current alpha                                                 
            oldalpha(iY) = alpha(iG,iY)                                         
          END DO                                                                
c         for each y calculate new alpha from III (the equivalent of eq.B7)     
          DO iY = 1, nY                                                         
            TAU(iY) = TAUslb(iLfid,iY)                                          
          END DO                                                                
          DO iY = 1, nY                                                         
            DO iW = iY, nY                                                      
              aux(iW) = FovrF1*QF(iY)                                           
            END DO                                                              
            CALL Simpson(npY,iY,nY,TAU,aux,IntQF)                               
            alpha(iG,iY) = 3.*IntQF - fstarQ(iY)                                
            alpha(iG,iY) = oldalpha(iY)                                         
c           calculate temperature                                               
            Td(iG,iY) = Tsub(iG) * (alpha(iG,iY)/alpha(iG,1))**0.25             
            delta = DABS((alpha(iG,iY)-oldalpha(iY))/alpha(iG,iY))              
            IF (delta.GT.accConv) istop = 0                                     
          END DO                                                                
c       end of iterations                                                       
        END DO                                                                  
c     end of loop over grains (iG)                                              
      END DO                                                                    
c --------------------------------------------------------------------------    
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c **********************************************************************        
      SUBROUTINE SLBmat(TAUslb,mat0)                                            
c ======================================================================        
c This subroutine evaluates the radiative transfer matrix for slab geometry.    
c TAUslb is the array of optical depths along the line of sight. [MN, Dec'97]   
c ======================================================================        
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER nYok, nPok                                                        
      DOUBLE PRECISION                                                          
     &       Ude(npL,npY), Uds(npL,npY), Us(npL,npY), fde(npL,npY),             
     &       fds(npL,npY), fs(npL,npY), Utot(npL,npY), ftot(npL,npY),           
     &       Td(npG,npY), Ubol(npY), fbol(npY), Uchck(npL,npY),                 
     &       UbolChck(npY), Spectrum(npL), SpecChar(30,99), tauF(npY),          
     &       Intens(npL,npP+2), IntOut(20,npP+2), SmC(30,99),                   
     &       Yok(npY), Pok(npP), tauOut(npL), tauZout(npP+2), F1,               
     &       tr(npY), rg(npG,npY), fsL(npL,npY), fsR(npL,npY), Eps(npY)         
      COMMON /solution/ Ude, Uds, Us, fde, fds, fs, Utot, ftot, Td,             
     &       Ubol, fbol, Uchck, UbolChck, Spectrum, SpecChar, tauF,             
     &       Intens, IntOut, SmC, Yok, Pok, tauOut, tauZout, F1,                
     &       tr, rg, fsL, fsR, Eps, nYok, nPok                                  
      INTEGER iL, iY, i, j, k                                                   
      DOUBLE PRECISION TAUslb(npL,npY), TAUr(npY), mat0(npL,npY,npY),           
     &         L1(npY,npY),L2(npY,npY),L3(npY,npY),L4(npY,npY),                 
     &         alpha(npY,npY), beta(npY,npY), gamma(npY,npY),                   
     &         delta(npY,npY), SUM                                              
c -------------------------------------------------------------------------     
      DO iL = 1, nL                                                             
       DO j = 1, nY                                                             
        DO k = 1, nY                                                            
          mat0(iL,j,k) = 0.0                                                    
        END DO                                                                  
       END DO                                                                   
      END DO                                                                    
c     -- evaluate matrix elements --                                            
c     loop over wavelengths                                                     
      DO iL = 1, nL                                                             
c      'Myspline' needs a 1D array                                              
        DO iY = 1, nY                                                           
          TAUr(iY) = TAUslb(iL,iY)                                              
        END DO                                                                  
c       Get the spline coefficients of the source function                      
        CALL MYSPLINE(TAUr,nY,alpha,beta,gamma,delta)                           
        CALL SLBAnalyt(TAUr,L1,L2,L3,L4,nY,npY)                                 
c       Matrix element for U:                                                   
        DO k = 1, nY                                                            
          DO j = 1, nY                                                          
            SUM = 0.                                                            
            DO i = 1, nY                                                        
             SUM = SUM + L1(i,k)*alpha(i,j)+L2(i,k)*beta(i,j)                   
     &                 + L3(i,k)*gamma(i,j)+L4(i,k)*delta(i,j)                  
            END DO                                                              
            mat0(iL,k,j) = SUM                                                  
          END DO                                                                
        END DO                                                                  
c     end of loop over wavelengths                                              
      END DO                                                                    
c     save Y to Yok; needed for analysis in cases when requirement for          
c     finer grids cannot be satisfied and previous solution is used for         
c     output                                                                    
      nYok = nY                                                                 
      DO iY = 1, nY                                                             
        Yok(iY) = Y(iY)                                                         
      END DO                                                                    
c ----------------------------------------------------------------------        
      RETURN                                                                    
      END                                                                       
c ********************************************************************          
                                                                                
                                                                                
c ********************************************************************          
      SUBROUTINE SLBRadT(nG,error)                                              
c ======================================================================        
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER nYok, nPok                                                        
      DOUBLE PRECISION                                                          
     &       Ude(npL,npY), Uds(npL,npY), Us(npL,npY), fde(npL,npY),             
     &       fds(npL,npY), fs(npL,npY), Utot(npL,npY), ftot(npL,npY),           
     &       Td(npG,npY), Ubol(npY), fbol(npY), Uchck(npL,npY),                 
     &       UbolChck(npY), Spectrum(npL), SpecChar(30,99), tauF(npY),          
     &       Intens(npL,npP+2), IntOut(20,npP+2), SmC(30,99),                   
     &       Yok(npY), Pok(npP), tauOut(npL), tauZout(npP+2), F1,               
     &       tr(npY), rg(npG,npY), fsL(npL,npY), fsR(npL,npY), Eps(npY)         
      COMMON /solution/ Ude, Uds, Us, fde, fds, fs, Utot, ftot, Td,             
     &       Ubol, fbol, Uchck, UbolChck, Spectrum, SpecChar, tauF,             
     &       Intens, IntOut, SmC, Yok, Pok, tauOut, tauZout, F1,                
     &       tr, rg, fsL, fsR, Eps, nYok, nPok                                  
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER startyp(2), Nlamtr(2), nBB(2)                                     
      CHARACTER nameStar(2)*235                                                 
      DOUBLE PRECISION Tstar, lamtr(2,101), klam(2,100), Tbb(2,10),             
     &                 rellum(2,10), mu1, ksi, mu2, xSiO, r1rs                  
      COMMON /source/ Tstar, lamtr, klam, Tbb, rellum, mu1, ksi, mu2,           
     &                xSiO, r1rs, startyp, Nlamtr, nBB, nameStar                
      DOUBLE PRECISION TAUslb(npL,npY),fsbol(npY),fpbol(npY),fmbol(npY),        
     &                 AveDev, RMS, maxFerr, fmax, fmed                         
      COMMON /slab/ TAUslb, fsbol, fpbol, fmbol, AveDev, RMS, maxFerr,          
     &                 fmax, fmed                                               
      INTEGER iY, nG, itnum, itlim                                              
      INTEGER error, Conv, iter, Fconv, Uconv, BolConv, m14                     
      DOUBLE PRECISION mat0(npL,npY,npY), Em(npL,npY), alpha(npG,npY),          
     &    fdsp(npL,npY),fdsm(npL,npY),fdep(npL,npY),fdem(npL,npY),dmaxU,        
     &    fbolold(npY),Uold(npL,npY), dmaxF, olderr, oldmaxU                    
c --------------------------------------------------------------------------    
      IF (iInn.eq.1) THEN                                                       
        write(38,'(a8,i5)') '    nY= ',nY                                       
        write(38,*) '    iter   maxFerr     dmaxU       dmaxF'                  
      END IF                                                                    
c     The slab Tau-profile at the fiducious lambda (needed in PrOut)            
      DO iY = 1, nY                                                             
        tr(iY) = TAUslb(iLfid,iY)/TAUslb(iLfid,nY)                              
      END DO                                                                    
c     generate stellar moments                                                  
      CALL SLBStar(error)                                                       
c     generate albedo through the envelope                                      
      CALL getOmega(nG)                                                         
c     finish when file with the stellar spectrum is not available               
      IF (error.EQ.3) goto 999                                                  
c     generate the first approximations for Td and alpha                        
      CALL SLBiniT(nG,alpha)                                                    
c     find radiative transfer matrices                                          
      IF (iX.GE.1) write(18,*)' Calculating weight matrices'                    
      CALL SLBmat(TAUslb,mat0)                                                  
      IF (iX.GE.1) THEN                                                         
        write(18,*)' Weight matrices OK, calculating Tdust'                     
      END IF                                                                    
      itlim = 9000                                                              
      itnum = 4                                                                 
      Conv = 0                                                                  
      iter = 0                                                                  
      oldmaxU = 0.                                                              
      olderr = 0.                                                               
      m14 = 0                                                                   
c     === Iterations over dust temperature =========                            
      DO WHILE (Conv.EQ.0.AND.iter.LE.itlim)                                    
        iter = iter + 1                                                         
c       find emission term                                                      
        CALL Emission(0,nG,Td,alpha,abund,Us,Em)                                
c       solve for Utot                                                          
        CALL Invert(nY,nL,mat0,Us,Em,omega,Utot,Uold,error)                     
c       find new Td and alpha, and check convergence                            
        CALL FindTemp(0,Utot,nG,Td,alpha)                                       
c       --------------------------------------                                  
c       every itnum-th iteration check convergence:                             
c       first find 'old' flux (i.e. in the previous iteration)                  
        IF (MOD(iter+1,itnum).EQ.0) THEN                                        
c         Find the diffuse scattered flux(fl=1 for scatt. and fl=0 is for emissi
          CALL SLBdifF(1,omega,TAUslb,Utot,nL,nY,fds,fdsp,fdsm)                 
c         Find the diffuse emitted flux                                         
          CALL SLBdifF(0,omega,TAUslb,Em,nL,nY,fde,fdep,fdem)                   
          CALL Add(npY,nY,npL,nL,fs,fds,fde,ftot)                               
c         find bolometric flux                                                  
          CALL Bolom(ftot,fbolold)                                              
        END IF                                                                  
        IF (MOD(iter,itnum).EQ.0) THEN                                          
c         Find the diffuse scattered flux                                       
          CALL SLBdifF(1,omega,TAUslb,Utot,nL,nY,fds,fdsp,fdsm)                 
c         Find the diffuse emitted flux                                         
          CALL SLBdifF(0,omega,TAUslb,Em,nL,nY,fde,fdep,fdem)                   
c         Add them to the stellar flux to find total flux                       
          CALL Add(npY,nY,npL,nL,fs,fds,fde,ftot)                               
c         find bolometric flux                                                  
          CALL Bolom(ftot,fbol)                                                 
c         check convergence of bolometric flux                                  
          CALL Converg1(nY,accFbol,dynrange,fbolold,fbol,Fconv,dmaxF)           
c         check convergence of energy density                                   
          CALL Converg2(nY,nL,accConv,dynrange,Uold,Utot,Uconv,dmaxU)           
c         FindErr calculates the err acc. to min/max values                     
          CALL FindErr(fbol,maxFerr,nY)                                         
c         flag m14 is to prevent unnecessary looping if convergence             
c         is extremely slow and err does not improve.                           
          IF(iter.GE.1000.AND.                                                  
     &           dabs(dmaxU-oldmaxU).LE.1.0d-5.AND.                             
     &             dabs(maxFerr-olderr).LE.1.0d-4) THEN                         
            m14 = 1                                                             
            Conv = 1                                                            
          END IF                                                                
          olderr = maxFerr                                                      
          oldmaxU = dmaxU                                                       
c ------  printout of errors and convergence with iter.(inner flag): ----       
          IF(iInn.EQ.1) THEN                                                    
            write(38,'(i7,1p,3e12.4)') iter, maxFerr, dmaxU, dmaxF              
          END IF                                                                
c ------------------------------------------------------------------------      
          IF (maxFerr.LE.accuracy) THEN                                         
            BolConv = 1                                                         
          ELSE                                                                  
            BolConv = 0                                                         
          END IF                                                                
c         total criterion for convergence: Utot must converge, and ftot         
c         must either converge or have the required accuracy                    
          IF (Uconv*(Fconv+BolConv).GT.0) Conv = 1                              
        END IF                                                                  
      END DO                                                                    
c     === The End of Iterations over Td ===                                     
      IF (iX.GE.1) THEN                                                         
        IF (iter.le.itlim) THEN                                                 
          write(18,*)' Convergence achieved, number of'                         
          write(18,'(a34,i4)')                                                  
     &      ' iterations over energy density: ',iter                            
        ELSE                                                                    
          write(18,'(a38,i5,a6)')                                               
     &    '  No convergence on energy density in',iter,' iter.'                 
          iWARNING = iWARNING + 1                                               
        END IF                                                                  
        write(18,'(a25,1p,E9.2)')                                               
     &    '  Max error in bol.flux:',maxFerr                                    
      END IF                                                                    
c     calculate the emission term for the converged Td                          
      CALL Emission(0,nG,Td,alpha,abund,Us,Em)                                  
c     calculate flux                                                            
      CALL SLBdifF(1,omega,TAUslb,Utot,nL,nY,fds,fdsp,fdsm)                     
      CALL SLBdifF(0,omega,TAUslb,Em,nL,nY,fde,fdep,fdem)                       
      CALL Add(npY,nY,npL,nL,fs,fds,fde,ftot)                                   
      CALL Add2(fdsp,fdep,fpbol,nY)                                             
      CALL Add2(fdsm,fdem,fmbol,nY)                                             
      CALL Bolom(ftot,fbol)                                                     
      CALL Bolom(fs,fsbol)                                                      
c     Find the err acc. to min/max values of fbol:                              
      CALL FindErr(fbol,maxFerr,nY)                                             
c     find average value of fbol and some misceleneous quantities               
      CALL SLBmisc(fbol,fmax,fmed,AveDev,RMS,nY)                                
c     calculate additional output quantities                                    
      CALL Multiply(1,npY,nY,npL,nL,mat0,Utot,omega,0,Us,Uds)                   
      CALL Multiply(0,npY,nY,npL,nL,mat0,Em,omega,0,Us,Ude)                     
      CALL Add(npY,nY,npL,nL,Us,Uds,Ude,Uchck)                                  
      CALL Bolom(Utot,Ubol)                                                     
      CALL Bolom(Uchck,UbolChck)                                                
c     fmed is printed in .out file as f1 = fmed = F/Fe1                         
c ============ if the inner flag iInn=1:  =========                             
      IF(iX.GE.1 .AND. iInn.EQ.1) THEN                                          
        write(18,'(a11,1p,E11.3)')'   TAUfid =',TAUfid                          
        write(18,'(a11,1p,E11.3)')'     fmed =',fmed                            
        write(18,'(a11,1p,E11.3)')'     fmax =',fmax                            
        write(18,'(a11,1p,E11.3)')'  MAX err =',maxFerr                         
        write(18,*)                                                             
     &'     tr      fbol       fsbol      fdbol       Ubol     fpbol'           
        DO iY = 1, nY                                                           
         write(18,'(1p,6E11.3)') tr(iY), fbol(iY), fsbol(iY),                   
     &               (fpbol(iY)-fmbol(iY)),Ubol(iY),fpbol(iY)                   
        END DO                                                                  
      END IF                                                                    
c ===================================================                           
c -----------------------------------------------------------------------       
999   RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE SLBsolve(model,nG,error)                                       
c =======================================================================       
c This subroutine solves the continuum radiative transfer problem in            
c planar geometry.                                      [MN, Dec.'97]           
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      DOUBLE PRECISION TAUslb(npL,npY),fsbol(npY),fpbol(npY),fmbol(npY),        
     &                 AveDev, RMS, maxFerr, fmax, fmed                         
      COMMON /slab/ TAUslb, fsbol, fpbol, fmbol, AveDev, RMS, maxFerr,          
     &                 fmax, fmed                                               
      INTEGER nYok, nPok                                                        
      DOUBLE PRECISION                                                          
     &       Ude(npL,npY), Uds(npL,npY), Us(npL,npY), fde(npL,npY),             
     &       fds(npL,npY), fs(npL,npY), Utot(npL,npY), ftot(npL,npY),           
     &       Td(npG,npY), Ubol(npY), fbol(npY), Uchck(npL,npY),                 
     &       UbolChck(npY), Spectrum(npL), SpecChar(30,99), tauF(npY),          
     &       Intens(npL,npP+2), IntOut(20,npP+2), SmC(30,99),                   
     &       Yok(npY), Pok(npP), tauOut(npL), tauZout(npP+2), F1,               
     &       tr(npY), rg(npG,npY), fsL(npL,npY), fsR(npL,npY), Eps(npY)         
      COMMON /solution/ Ude, Uds, Us, fde, fds, fs, Utot, ftot, Td,             
     &       Ubol, fbol, Uchck, UbolChck, Spectrum, SpecChar, tauF,             
     &       Intens, IntOut, SmC, Yok, Pok, tauOut, tauZout, F1,                
     &       tr, rg, fsL, fsR, Eps, nYok, nPok                                  
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER startyp(2), Nlamtr(2), nBB(2)                                     
      CHARACTER nameStar(2)*235                                                 
      DOUBLE PRECISION Tstar, lamtr(2,101), klam(2,100), Tbb(2,10),             
     &                 rellum(2,10), mu1, ksi, mu2, xSiO, r1rs                  
      COMMON /source/ Tstar, lamtr, klam, Tbb, rellum, mu1, ksi, mu2,           
     &                xSiO, r1rs, startyp, Nlamtr, nBB, nameStar                
      INTEGER model, error, nG, iterFbol, FbolOK, grid                          
      DOUBLE PRECISION devmax, oldFerr, mu, h, k                                
c --------------------------------------------------------------------------    
c     counter for iterations over bolometric flux conservation                  
      iterFbol = 0                                                              
      FbolOK = 0                                                                
      oldFerr = 100.                                                            
c     make the grid for the smaller mu                                          
      IF (dabs(mu1).LE.dabs(mu2)) THEN                                          
        mu = dabs(mu1)                                                          
      ELSE                                                                      
        mu = dabs(mu2)                                                          
      END IF                                                                    
      IF((mu.LT.0.1).AND.(TAUmax.GE.20.)) THEN                                  
        grid = 1                                                                
c      oblique illumination grid (grid=1) with initial parameters               
       IF(npY.LT.50) THEN                                                       
        iWARNING = iWARNING + 1                                                 
        error = 4                                                               
        RETURN                                                                  
       ELSE                                                                     
        nY = 50                                                                 
       END IF                                                                   
        h = 0.5                                                                 
        k = 5.0                                                                 
        CALL SLBtrace(TAUmax,TAUTot,TAUslb,Y,mu,h,k,nL,nY)                      
      ELSE                                                                      
        grid = 0                                                                
c       generate the increment for the initial TAU-grid (grid=0)                
        CALL SLBY(TAUmax,Y,nY)                                                  
      END IF                                                                    
c ------------- Loop over bol.flux conservation -------------                   
      DO WHILE (FbolOK.EQ.0)                                                    
       iterFbol = iterFbol + 1                                                  
        IF((iterFbol.GE.4).AND.(grid.EQ.1).AND.(iX.GE.1)) THEN                  
          write(18,*)'  Can not improve the oblique'                            
          write(18,*)'  illumination grid anymore  '                            
          iWARNING = iWARNING + 1                                               
          FbolOK = 2                                                            
          RETURN                                                                
        END IF                                                                  
        IF (iX.GE.1) THEN                                                       
          write(18,*)'  ',iterFbol,' iteration over Fbol'                       
        END IF                                                                  
        IF (iVerb.EQ.2) write(*,'(a14,i3,a20)')                                 
     &     ' In SLBsolve: ',iterFbol,' iteration over Fbol'                     
c       if the illumination angle is too large choose the proper grid           
        IF(grid.NE.1) THEN                                                      
          CALL SLBtau(TAUTot,TAUslb,Y,nL,nY)                                    
        END IF                                                                  
        IF (iX.GE.1) THEN                                                       
          write(18,'(a28,i3,a12)') '  Calculation for slab with ',              
     &                          nY,' grid points'                               
        END IF                                                                  
        IF (nY.GT.npY) THEN                                                     
         IF (iVerb.EQ.2) write(*,*)' npY needs to be increased!'                
         IF (iX.GE.1) THEN                                                      
          write(18,*) ' ********** MESSAGE from SOLVE *********'                
          write(18,*) ' npY has to be increased. Please, choose'                
          write(18,*) ' parameters for slab in file userpar.inc'                
          write(18,*) ' ***************************************'                
          iWARNING = iWARNING + 1                                               
          RETURN                                                                
         END IF                                                                 
        END IF                                                                  
c       solve the radiative transfer problem                                    
        CALL SLBRadT(nG,error)                                                  
         IF (iVerb.EQ.2) write(*,*)'Done with SLBRadT'                          
c       for tests of single lambda case:                                        
c       CALL SLBgray(model,nG,error)                                            
        IF(maxFerr.LE.accuracy) THEN                                            
          FbolOK = 1                                                            
        ELSE                                                                    
         IF(grid.EQ.1) THEN                                                     
           IF(iterFbol.EQ.1) THEN                                               
c            next try with 104pt grid:                                          
             h = h/1.5                                                          
             k = k*1.5                                                          
             nY = 104                                                           
             CALL SLBtrace(TAUmax,TAUTot,TAUslb,Y,mu,h,k,nL,nY)                 
           END IF                                                               
           IF(iterFbol.EQ.2) THEN                                               
c            last attempt with 160pt grid:                                      
             h = 0.25                                                           
             k = 7.5                                                            
             nY = 160                                                           
             CALL SLBtrace(TAUmax,TAUTot,TAUslb,Y,mu,h,k,nL,nY)                 
           END IF                                                               
         ELSE                                                                   
           CALL SLBacc(fbol,accuracy,devmax,FbolOK,error)                       
           IF(iterFbol.GT.1.AND.maxFerr.GT.oldFerr) THEN                        
            IF(iX.GE.1) THEN                                                    
             write(18,*)' ================ WARNING ================= '          
             write(18,*)' Error is increasing in spite of grid'                 
             write(18,*)' improvement. Stopping iterations over Fbol.'          
            END IF                                                              
            error = 0                                                           
            FbolOK = 2                                                          
            iWARNING = iWARNING + 1                                             
           END IF                                                               
           oldFerr = maxFerr                                                    
         END IF                                                                 
c         if the grid size limit is reached error=2                             
          IF (error.EQ.2.AND.iterFbol.EQ.1) THEN                                
c          if this is the first calculation end this model                      
           IF(iX.GE.1) THEN                                                     
            write(18,*)' =========== IMPORTANT WARNING =========== '            
            write(18,*)' The limit for grid size is already reached'            
            write(18,*)' and flux conservation can not be improved'             
            write(18,'(a,1p,e9.2)')'  Max deviation of Fbol ',maxFerr           
            write(18,*)' Treat all results with caution!'                       
           END IF                                                               
           error = 0                                                            
           FbolOK = 2                                                           
           iWARNING = iWARNING + 1                                              
          END IF                                                                
c         if this is a higher iteration use previous solution                   
          IF (error.EQ.2) THEN                                                  
            IF (iX.GE.1.AND.iterFbol.GT.1) THEN                                 
              write(18,*)' ======= IMPORTANT WARNING ======== '                 
              write(18,*)' In trying to conserve Fbol reached'                  
              write(18,*)' the limit for grid sizes.  '                         
              write(18,'(a,1p,e9.2)')'  Max deviation of Fbol:',maxFerr         
              write(18,*)' Treat all results with caution!'                     
            END IF                                                              
            error = 0                                                           
            FbolOK = 2                                                          
            iWARNING = iWARNING + 1                                             
          END IF                                                                
c         if Fbol not conserved try again with a finer grid                     
          IF (FbolOK.EQ.0 .AND. iX.GE.1) THEN                                   
            write(18,*)'  ******** MESSAGE from SOLVE ********'                 
            write(18,'(a,1p,e9.2)')                                             
     &                '  Max deviation of Fbol:', maxFerr                       
            write(18,*)'  Trying again with finer grids'                        
          END IF                                                                
c        if could not conserve Fbol in 10 trials give it up                     
         IF (FbolOK.EQ.0 .AND. iterFbol.GE.10) THEN                             
          IF (iX.GE.1) THEN                                                     
          write(18,*)' **********  WARNING from SOLVE  **********'              
          write(18,*)' Could not obtain required accuracy in 10 trials.'        
          write(18,'(a,1p,e9.2)')'  Max deviation of Fbol:',maxFerr             
          write(18,*)' !!!!  Treat all results with caution  !!!!'              
          write(18,*)' ******************************************'              
          END IF                                                                
          iWARNING = iWARNING + 1                                               
          FbolOK = 2                                                            
         END IF                                                                 
        END IF                                                                  
c     end of loop over flux conservation                                        
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE SLBStar(error)                                                 
c =======================================================================       
c This subroutine generates the stellar moments in case of slab geometry.       
c                                                          [MN, Feb.'99]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      INTEGER nYok, nPok                                                        
      DOUBLE PRECISION                                                          
     &       Ude(npL,npY), Uds(npL,npY), Us(npL,npY), fde(npL,npY),             
     &       fds(npL,npY), fs(npL,npY), Utot(npL,npY), ftot(npL,npY),           
     &       Td(npG,npY), Ubol(npY), fbol(npY), Uchck(npL,npY),                 
     &       UbolChck(npY), Spectrum(npL), SpecChar(30,99), tauF(npY),          
     &       Intens(npL,npP+2), IntOut(20,npP+2), SmC(30,99),                   
     &       Yok(npY), Pok(npP), tauOut(npL), tauZout(npP+2), F1,               
     &       tr(npY), rg(npG,npY), fsL(npL,npY), fsR(npL,npY), Eps(npY)         
      COMMON /solution/ Ude, Uds, Us, fde, fds, fs, Utot, ftot, Td,             
     &       Ubol, fbol, Uchck, UbolChck, Spectrum, SpecChar, tauF,             
     &       Intens, IntOut, SmC, Yok, Pok, tauOut, tauZout, F1,                
     &       tr, rg, fsL, fsR, Eps, nYok, nPok                                  
      INTEGER startyp(2), Nlamtr(2), nBB(2)                                     
      CHARACTER nameStar(2)*235                                                 
      DOUBLE PRECISION Tstar, lamtr(2,101), klam(2,100), Tbb(2,10),             
     &                 rellum(2,10), mu1, ksi, mu2, xSiO, r1rs                  
      COMMON /source/ Tstar, lamtr, klam, Tbb, rellum, mu1, ksi, mu2,           
     &                xSiO, r1rs, startyp, Nlamtr, nBB, nameStar                
      DOUBLE PRECISION TAUslb(npL,npY),fsbol(npY),fpbol(npY),fmbol(npY),        
     &                 AveDev, RMS, maxFerr, fmax, fmed                         
      COMMON /slab/ TAUslb, fsbol, fpbol, fmbol, AveDev, RMS, maxFerr,          
     &                 fmax, fmed                                               
      CHARACTER*235 line                                                        
      INTEGER iY,iL,ios1,iLs,nLs, k,kstop,i,is, error,isn, Nlambdam             
c     Nlambdam is the max number entries for a user supplied stellar spectrum   
      PARAMETER (Nlambdam = 10000)	                       	                     
      DOUBLE PRECISION lambdaS(Nlambdam), Llamstar(Nlambdam), Lstar,            
     &       Stellar(Nlambdam), llS(Nlambdam), lS(Nlambdam), fsrc(npL),         
     &       fpl(npL), Bb(2), arg, EMfunc, a, b, x, Planck, BP,                 
     &       dyn2, UsL, UsR, Eint2, Eint3                                       
c -----------------------------------------------------------------------       
      IF(ksi.GT.0) THEN                                                         
        isn = 2                                                                 
      ELSE                                                                      
        isn = 1                                                                 
      END IF                                                                    
      dyn2 = dynrange*dynrange                                                  
      DO is = 1, isn                                                            
c      if startyp.GE.4 stellar spectrum is read from file 'nameStar'            
       IF (startyp(is).GE.4.AND.startyp(is).LE.6) THEN                          
        open(3,ERR=998,file=nameStar(is),STATUS='OLD')                          
        rewind(3)                                                               
        read(3,'(a235)',ERR=998) line                                           
        read(3,'(a235)',ERR=998) line                                           
        read(3,'(a235)',ERR=998) line                                           
        ios1 = 0                                                                
        iLs = 0                                                                 
        DO WHILE (ios1.ge.0)                                                    
          read(3,*,END=900,ERR=998,iostat=ios1) a,b                             
          IF (ios1.ge.0) THEN                                                   
            iLs = iLs + 1                                                       
            lambdaS(iLs) = a                                                    
            IF (a.LE.0.0) goto 998                                              
c           it is assumed that Llamstar is L_Lambda, but...                     
c           if startyp.EQ.3 then file gives lambda*L_lambda                     
            IF (startyp(is).EQ.4) Llamstar(iLs) = b / a                         
c           if startyp.EQ.4 then file gives L_lambda                            
            IF (startyp(is).EQ.5) Llamstar(iLs) = b                             
c           if startyp.EQ.5 then file gives Lnu=lambda**2*L_lambda              
            IF (startyp(is).EQ.6) Llamstar(iLs) = b / a / a                     
           END IF                                                               
        END DO                                                                  
900     close(3)                                                                
        IF (iLs.LT.2) goto 998                                                  
        nLs = iLs                                                               
c       if input wavelengths in descending order turn them around               
        IF (lambdaS(1).GT.lambdaS(2)) THEN                                      
          DO iLs = 1, nLs                                                       
            llS(iLs) = lambdaS(iLs)                                             
            lS(iLs) = Llamstar(iLs)                                             
          END DO                                                                
          DO iLs = 1, nLs                                                       
            lambdaS(iLs) = llS(nLs+1-iLs)                                       
            Llamstar(iLs) = lS(nLs+1-iLs)                                       
          END DO                                                                
        END IF                                                                  
c       normalize stellar spectrum                                              
        CALL Simpson(Nlambdam,1,nLs,lambdaS,Llamstar,Lstar)                     
c       generate dimensionless stellar spectrum                                 
        DO iLs = 1, nLs                                                         
          Stellar(iLs) = lambdaS(iLs) * Llamstar(iLs) / Lstar                   
        END DO                                                                  
       ELSE                                                                     
c       if startyp.EQ.3 generate power-law spectrum                             
        IF (startyp(is).EQ.3) THEN                                              
          fsrc(is) = 1.0                                                        
          IF (Nlamtr(is).GT.1) THEN                                             
            DO i = 2, Nlamtr(is)                                                
              fsrc(i) = fsrc(i-1)*                                              
     &                        (lamtr(is,i-1)/lamtr(is,i))**klam(is,i-1)         
            END DO                                                              
          END IF                                                                
          DO iL = 1, nL                                                         
            IF((lambda(iL)-lamtr(is,1))*(lambda(iL)-                            
     &                             lamtr(is,Nlamtr(is)+1)).LE.0.0) THEN         
              kstop = 0                                                         
              k = 0                                                             
              DO WHILE (kstop.EQ.0)                                             
                k = k + 1                                                       
                IF (lambda(iL).GE.lamtr(is,k)) THEN                             
                  kstop = 1                                                     
                  fpl(iL) = fsrc(k)*(lamtr(is,k)/lambda(iL))**klam(is,k)        
                END IF                                                          
              END DO                                                            
            ELSE                                                                
              fpl(iL) = 0.0                                                     
            END IF                                                              
          END DO                                                                
        END IF                                                                  
       END IF                                                                   
       DO iY = 1, nY                                                            
c       loop over wavelengths                                                   
        DO iL = 1, nL                                                           
         IF (startyp(is).EQ.1) THEN                                             
            Bb(is) = 0.0                                                        
            DO k = 1, nBB(is)                                                   
              x = 14400.0 / lambda(iL) / Tbb(is,k)                              
              Bb(is) = Bb(is) + rellum(is,k)*Planck(x)                          
            END DO                                                              
         ELSE IF (startyp(is).EQ.2) THEN                                        
            Bb(is) = EMfunc(lambda(iL),Tbb(1,1),xSiO)                           
         ELSE IF (startyp(is).EQ.3) THEN                                        
            Bb(is) = fpl(iL)                                                    
         ELSE IF (lambda(iL).GT.lambdaS(nLs)) THEN                              
c           for lambda longer than the longest entry in nameStar                
c           assume Rayleigh-Jeans tail                                          
            Bb(is) = Stellar(nLs)*(lambdaS(nLs)/lambda(iL))**3.                 
         ELSE IF (lambda(iL).LT.lambdaS(1)) THEN                                
c           if shorter than the shortest assume 0                               
            Bb(is) = 0.0                                                        
         ELSE                                                                   
c          from a file: if within limits interpolate                            
           CALL LinInter(Nlambdam,nLs,lambdaS,Stellar,lambda(iL),iLs,BP)        
           Bb(is) = BP                                                          
         END IF                                                                 
c   ----------  done with stellar spectrum ---------------                      
c        stellar part of flux and en.density                                    
         IF (is.EQ.1) THEN                                                      
           IF((mu1+1.).LE.1.0e-5) THEN                                          
             fsL(iL,iY) = Bb(1)*Eint3(TAUslb(iL,iY))                            
           ELSE                                                                 
             x = TAUslb(iL,iY) / mu1                                            
             IF(x.GE.50.) THEN                                                  
               fsL(iL,iY) = 0.                                                  
             ELSE                                                               
               fsL(iL,iY) = Bb(1)*exp(-x)                                       
             END IF                                                             
           END IF                                                               
         END IF                                                                 
         IF (is.EQ.2) THEN                                                      
           IF((mu2+1.).LE.1.0e-5) THEN                                          
            arg = TAUslb(iL,nY)-TAUslb(iL,iY)                                   
            fsR(iL,iY) = Bb(2)*Eint3(arg)                                       
           ELSE                                                                 
             x = (TAUslb(iL,nY)-TAUslb(iL,iY)) / mu2                            
             IF(x.GE.50.) THEN                                                  
                 fsR(iL,iY) = 0.                                                
             ELSE                                                               
               fsR(iL,iY) = Bb(2)*exp(-x)                                       
             END IF                                                             
           END IF                                                               
         END IF                                                                 
        END DO                                                                  
       END DO                                                                   
c     end do over 'is' - the counter for sources                                
      END DO                                                                    
                                                                                
      DO iY = 1, nY                                                             
c       loop over wavelengths                                                   
        DO iL = 1, nL                                                           
          fs(iL,iY) = fsL(iL,iY) - ksi * fsR(iL,iY)                             
c         find en.density if diffuse illumination on the left (mu1=-1)          
          IF((mu1+1.).LE.1.0e-5) THEN                                           
c           to prevent exp. underflow                                           
            IF (TAUslb(iL,iY).GE.300.) THEN                                     
             UsL = 0.                                                           
            ELSE                                                                
             UsL = fsL(iL,iY)/Eint3(TAUslb(iL,iY))*Eint2(TAUslb(iL,iY))         
            END IF                                                              
          ELSE                                                                  
            UsL = fsL(iL,iY)/mu1                                                
          END IF                                                                
c         if diffuse illumination from the right                                
          IF((mu2+1.).LE.1.0e-5) THEN                                           
            arg = TAUslb(iL,nY)-TAUslb(iL,iY)                                   
c           to prevent exp. underflow                                           
            IF (arg.GE.300.) THEN                                               
             UsR = 0.                                                           
            ELSE                                                                
             UsR = fsR(iL,iY)/Eint3(arg)*Eint2(arg)                             
            END IF                                                              
          ELSE                                                                  
            UsR = fsR(iL,iY)/mu2                                                
          END IF                                                                
           Us(iL,iY) = UsL + ksi * UsR                                          
c         here only Us needs limit from below; fs can be negative though!       
          IF (Us(iL,iY).LT.dyn2) Us(iL,iY) = 0.0                                
        END DO                                                                  
      END DO                                                                    
                                                                                
c     normalize stellar quantities with the stellar bolometric flux             
      CALL Bolom(fsL,fsbol)                                                     
c     fsbol(1) is always non-zero                                               
      DO iY = 1, nY                                                             
        DO iL = 1, nL                                                           
           Us(iL,iY) = Us(iL,iY) / fsbol(1)                                     
           fs(iL,iY) = fs(iL,iY) / fsbol(1)                                     
           fsL(iL,iY) = fsL(iL,iY) / fsbol(1)                                   
           fsR(iL,iY) = fsR(iL,iY) / fsbol(1)                                   
        END DO                                                                  
      END DO                                                                    
      error = 0                                                                 
      goto 999                                                                  
998   write(12,*)' *** FATAL ERROR IN DUSTY! *************************'         
      write(12,*)' File with the spectral shape of external radiation:'         
      write(12,'(a2,a70)')'  ',nameStar(is)                                     
      write(12,*)' is missing or not properly formatted?!'                      
      write(12,*)' ***************************************************'         
      error = 3                                                                 
c -----------------------------------------------------------------------       
999   RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c **********************************************************************        
      SUBROUTINE SLBTau(TauTot,TAUslb,delT,nL,nY)                               
c ======================================================================        
c It generates TAUslb(npL,npY) grid with given spacing delT, calculated         
c in sub SLBy. In the current version it is exp near the slab faces and         
c equidistant in the middle.			             [MN,Sep'98]                         
c ======================================================================        
      IMPLICIT none                                                             
      INTEGER nL, nY, iY, iL                                                    
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      DOUBLE PRECISION TAUslb(npL,npY),TauTot(npL),TAU(npY),delT(npY)           
c ----------------------------------------------------------------------        
       TAU(1) = 0.                                                              
       DO iY = 2, nY                                                            
         TAU(iY) = TAU(iY-1) + delT(iY)                                         
       END DO                                                                   
      DO iL = 1, nL                                                             
        DO iY = 1, nY                                                           
          TAUslb(iL,iY) = TauTot(iL)*TAU(iY)/TAU(nY)                            
        END DO                                                                  
      END DO                                                                    
c ----------------------------------------------------------------------        
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE SLBtrace(TAUmax,TAUTot,TAUslb,Y,mu,h,k,nL,nY)                  
c ======================================================================        
c     This subroutine generates the TAU-grid for slab in the case of            
c     grazing angle of incidence.                            [MN,Jan'99]        
c ======================================================================        
      IMPLICIT NONE                                                             
      INTEGER nL, nY, iY, iL, Nthick, Nskin, iadj                               
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      DOUBLE PRECISION TAUslb(npL,npY), TAUTot(npL), TAU(npY), Y(npY),          
     &              fsort(npY), TAUmax, mu, h, k, step                          
c ----------------------------------------------------------------------        
c     Number of pts. in the skin layer with depth=min{mu1,mu2}:                 
      Nskin = k/h                                                               
c     Counter for the end of the layer adjacent to the skin layer               
      iadj = 2*Nskin                                                            
c     Number of pts. in the opt.thick inner part:                               
      Nthick = nY - 4*Nskin                                                     
      iY = 1                                                                    
      TAU(1) = 0.0                                                              
      TAU(nY) = 1.0                                                             
c     resolve the surface boundary layer, where tau = 0 .. k*mu                 
      DO WHILE(iY.LT.nY/2)                                                      
      iY = iY + 1                                                               
        IF((TAU(iY).LE.k*mu).AND.(iY.LE.Nskin)) THEN                            
          TAU(iY) = TAU(1) + (iY-1)*h*mu/TAUmax                                 
          TAU(nY-iY+1) = TAU(nY) - TAU(iY)                                      
        ELSE                                                                    
c         resolve the adjacent layer                                            
          IF(iY.LE.iadj) THEN                                                   
            TAU(iY) = TAU(iY-1) + h/TAUmax                                      
            TAU(nY-iY+1) = TAU(nY) - TAU(iY)                                    
          ELSE                                                                  
c           for the opt.thick inner part:                                       
            step = (TAU(nY-iadj+1)-TAU(iadj)) / (Nthick+1)                      
            TAU(iY) = TAU(iY-1) + step                                          
            TAU(nY-iY+1) = TAU(nY) - TAU(iY)                                    
          END IF                                                                
        END IF                                                                  
      END DO                                                                    
c     sorting (to avoid overlapping or mismatched pieces)                       
      DO iY = 1, nY                                                             
        fsort(iY) = TAU(iY)                                                     
      END DO                                                                    
      CALL sort(fsort,nY)                                                       
      DO iY = 1, nY                                                             
        TAU(iY) = fsort(iY)                                                     
        Y(iY) = TAU(iY)                                                         
      END DO                                                                    
      DO iL = 1, nL                                                             
        DO iY = 1, nY                                                           
          TAUslb(iL,iY) = TAUTot(iL)*TAU(iY)                                    
        END DO                                                                  
      END DO                                                                    
c ----------------------------------------------------------------------        
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE SLBy(TAUmax,dTAU,nY)                                           
c =======================================================================       
c     This subroutine generates the increment for the initial TAU-grid.         
c     It is equidistant in TAU in the middle and equidistant in log(TAU)        
c     near the two faces. It has 15 pts for small (tauV.LE.1) and 30pts         
c     for large total opt.depth (tauV.LE.1000). In some cases of tauV.GE.500.   
c     80pts grid is better to start with. If this is the case, comment          
c     the larger value of lim and uncomment lim=500, which will switch to       
c     80pt grid. Be sure you have selected the line for slab calculation        
c     in 'userpar.inc'.                                      [MN,Sep'98]        
c =======================================================================       
      IMPLICIT NONE                                                             
      INTEGER nY, iY                                                            
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      DOUBLE PRECISION dTAU(npY), TAUmax                                        
c ----------------------------------------------------------------------        
      IF(TAUmax.LE.10.) THEN                                                    
          nY = 15                                                               
          DO iY = 1, nY                                                         
           IF(iY.LE.5) THEN                                                     
             dTAU(iY) = 0.5*(exp(0.5*iY)-exp(0.5))                              
           ELSE IF(iY.GT.5.and.iY.LE.12) THEN                                   
             dTAU(iY) = dTAU(iY-1)                                              
           ELSE                                                                 
             dTAU(iY) = dTAU(2+nY-iY)                                           
           END IF                                                               
          END DO                                                                
      ELSE                                                                      
       IF(TAUmax.LE.10000.) THEN                                                
         nY = 30                                                                
         DO iY = 1, nY                                                          
           IF(iY.LE.10) THEN                                                    
             dTAU(iY) = 0.5*(exp(0.5*iY)-exp(0.5))                              
            ELSE IF(iY.GT.10.and.iY.LE.22) THEN                                 
             dTAU(iY) = dTAU(iY-1)                                              
            ELSE                                                                
             dTAU(iY) = dTAU(2+nY-iY)                                           
           END IF                                                               
        END DO                                                                  
       ELSE                                                                     
        nY = 80                                                                 
        DO iY = 1, nY                                                           
          IF(iY.LE.11) THEN                                                     
             dTAU(iY) = 0.5*(exp(0.5*iY)-exp(0.5))                              
          ELSE IF(iY.GT.11.and.iY.LE.71) THEN                                   
             dTAU(iY) = dTAU(iY-1)                                              
          ELSE                                                                  
             dTAU(iY) = dTAU(2+nY-iY)                                           
          END IF                                                                
        END DO                                                                  
       END IF                                                                   
      END IF                                                                    
c ----------------------------------------------------------------------        
      RETURN                                                                    
      END                                                                       
c *********************************************************************         
                                                                                
c =========================================================================     
c     This is are the subroutines analyzing the results of the radiative        
c     transfer calculation.                                   [MN, Mar'99]      
c =========================================================================     
C     Table of Contents                                                         
C                                                                               
C     ANALYSIS                                                                  
C     ASSPROP                                                                   
C     BOLOM                                                                     
C     CHKBOLOM                                                                  
C     CHKCONV                                                                   
C     CHKFLUX                                                                   
C     CONVERG1                                                                  
C     CONVERG2                                                                  
C     CONVOLVE                                                                  
C     CONV2D                                                                    
C     ETA                                                                       
C     ETAFUN                                                                    
C     FINDINT                                                                   
C     GETBOUT                                                                   
C     GETOMEGA                                                                  
C     GETOPTPR                                                                  
C     GETPROP                                                                   
C     GETSIZES                                                                  
C     GETTAU                                                                    
C     GETETAZP                                                                  
C     IMAGFN                                                                    
C     MIE                                                                       
C     PSFN                                                                      
C     SETUPETA                                                                  
C     SIZEDIST                                                                  
C     SPECTRAL                                                                  
C     SPFEATUR                                                                  
C     PHILAM                                                                    
C     VISIBILI                                                                  
C     VISI2D                                                                    
c =======================================================================       
                                                                                
c ***********************************************************************       
      SUBROUTINE ANALYSIS(model,ETAzp,error)                                    
c =======================================================================       
c This subroutine analyzes the solution. It finds the flux conservation         
c accuracy and evaluates many output quantites like QF(y), TAUF(y),Psi, F1,     
c the rad.pressure force, dynamical quantities etc.                             
c                                                   [ZI,Mar'96;MN,Mar'99]       
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      DOUBLE PRECISION TAUslb(npL,npY),fsbol(npY),fpbol(npY),fmbol(npY),        
     &                 AveDev, RMS, maxFerr, fmax, fmed                         
      COMMON /slab/ TAUslb, fsbol, fpbol, fmbol, AveDev, RMS, maxFerr,          
     &                 fmax, fmed                                               
      INTEGER nYok, nPok                                                        
      DOUBLE PRECISION                                                          
     &       Ude(npL,npY), Uds(npL,npY), Us(npL,npY), fde(npL,npY),             
     &       fds(npL,npY), fs(npL,npY), Utot(npL,npY), ftot(npL,npY),           
     &       Td(npG,npY), Ubol(npY), fbol(npY), Uchck(npL,npY),                 
     &       UbolChck(npY), Spectrum(npL), SpecChar(30,99), tauF(npY),          
     &       Intens(npL,npP+2), IntOut(20,npP+2), SmC(30,99),                   
     &       Yok(npY), Pok(npP), tauOut(npL), tauZout(npP+2), F1,               
     &       tr(npY), rg(npG,npY), fsL(npL,npY), fsR(npL,npY), Eps(npY)         
      COMMON /solution/ Ude, Uds, Us, fde, fds, fs, Utot, ftot, Td,             
     &       Ubol, fbol, Uchck, UbolChck, Spectrum, SpecChar, tauF,             
     &       Intens, IntOut, SmC, Yok, Pok, tauOut, tauZout, F1,                
     &       tr, rg, fsL, fsR, Eps, nYok, nPok                                  
      INTEGER startyp(2), Nlamtr(2), nBB(2)                                     
      CHARACTER nameStar(2)*235                                                 
      DOUBLE PRECISION Tstar, lamtr(2,101), klam(2,100), Tbb(2,10),             
     &                 rellum(2,10), mu1, ksi, mu2, xSiO, r1rs                  
      COMMON /source/ Tstar, lamtr, klam, Tbb, rellum, mu1, ksi, mu2,           
     &                xSiO, r1rs, startyp, Nlamtr, nBB, nameStar                
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      DOUBLE PRECISION ugas(npY), qF(npY), vrat(npG,npY), Gamma(npY),           
     &       I1, I2, I3, CMdot, Cve, CM, Cr1                                    
      COMMON /dyn/ ugas, qF, vrat, Gamma, I1, I2, I3, CMdot, Cve, CM,           
     &       Cr1                                                                
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER iL, iY, model, iP, error                                          
      DOUBLE PRECISION ETAzp(npP,npY), QpTd(npG,npY), QpStar(npY), Psi,         
     &        qaux(npL), qaux2(npL), resaux, xP, Planck, QUtot1, r_gd,          
     &        Eps1, aux, deltauF, C1,C2,C3, theta1, ugas_out, s4,Pi,mx          
c -----------------------------------------------------------------------       
      Pi = 2.0*ASIN(1.0)                                                        
      r_gd = 200.0                                                              
c     make sure that grids correspond to accepted solution                      
      nY = nYok                                                                 
      DO iY = 1, nY                                                             
        Y(iY) = Yok(iY)                                                         
      END DO                                                                    
      nP = nPok                                                                 
      DO iP = 1, nP                                                             
        P(iP) = Pok(iP)                                                         
      END DO                                                                    
c -------------                                                                 
c     spectrum (flux at the outer edge as a function of wavelength)             
      DO iL = 1, nL                                                             
        Spectrum(iL) = dabs(ftot(iL,nY))                                        
c       added in version dusty16.for - to prevent taking log from zero          
c       in Spectral [MN]:                                                       
        IF (Spectrum(iL).LE.1.0D-20) Spectrum(iL) = 1.0D-20                     
      END DO                                                                    
c -------------                                                                 
c     analyze bolometric flux error (1/2 of the max spread of fbol)             
      CALL FindErr(fbol,maxFerr,nY)                                             
c -------------                                                                 
c     find the flux averaged optical depth, tauF(y)                             
      IF (denstyp.NE.0) THEN                                                    
c      for spherical shell                                                      
       tauF(1) = 0.0                                                            
       DO iY = 2, nY                                                            
c        generate auxiliary function for integration:                           
c        loop over iL (wavelength)                                              
         DO iL = 1, nL                                                          
          qaux(iL)=TAUtot(iL)*ETAzp(1,iY)*dabs(ftot(iL,iY))/lambda(iL)          
         END DO                                                                 
         CALL Simpson(npL,1,nL,lambda,qaux,resaux)                              
         aux = resaux / ETAzp(1,iY)                                             
         deltauF = aux * (ETAzp(1,iY) - ETAzp(1,iY-1))                          
         tauF(iY) = tauF(iY-1) + deltauF                                        
       END DO                                                                   
      ELSE                                                                      
c      for slab                                                                 
       tauF(1) = 0.0                                                            
       DO iY = 1, nY                                                            
c        generate auxiliary function for integration:                           
c        loop over iL (wavelength)                                              
         DO iL = 1, nL                                                          
           qaux(iL)=TAUslb(iL,iY)*dabs(fTot(iL,iY))/lambda(iL)                  
           CALL Simpson(npL,1,nL,lambda,qaux,resaux)                            
           tauF(iY) = resaux                                                    
         END DO                                                                 
       END DO                                                                   
      END IF                                                                    
c -------------                                                                 
c     ratio of gravitational to radiation pressure force (isotropic             
c     scattering) per unit volume                                               
c     s4 = (L4sol/Msol)/(4*Pi*G*c*rho_s)/1e-6;                                  
c     rho_s=3000 kg.m-3, grain radius 'a' is in microns, aveV=4/3*Pi*<a^3>      
      IF(denstyp.NE.0) THEN                                                     
      s4 = 1.925 / (4.0*Pi*6.67d-11*3.0d08*3000.0*1.0D-06)                      
c      in case of sigma's from a file aveV=1 (initialized in GetOptPr)
       DO iY = 1, nY                                                            
        DO iL = 1, nL                                                           
          qaux(iL)=(SigmaA(1,iL)+SigmaS(1,iL))/aveV *                           
     &                             dabs(ftot(iL,iY))/lambda(iL)                 
        END DO                                                                  
        CALL Simpson(npL,1,nL,lambda,qaux,resaux)                               
        rg(1,iY) = s4 * resaux / r_gd                                           
c       If dust drift (dynamics case):                                          
        IF (RDW) rg(1,iY) = rg(1,iY)*vrat(1,iY)                                 
       END DO                                                                   
      END IF                                                                    
c -------------                                                                 
c     find the Planck averaged absorption efficiencies                          
      DO iY = 1, nY                                                             
c     generate auxiliary function for integration over wavelengths:             
       DO iL = 1, nL                                                            
         qaux(iL) = SigmaA(1,iL) * Us(iL,iY) / lambda(iL)                       
         xP = 14400.0 / Td(1,iY) / lambda(iL)                                   
         qaux2(iL) = SigmaA(1,iL) * Planck(xP) / lambda (iL)                    
       END DO                                                                   
       CALL Simpson(npL,1,nL,lambda,qaux,resaux)                                
       QpStar(iY) = resaux                                                      
       CALL Simpson(npL,1,nL,lambda,qaux2,resaux)                               
       QpTd(1,iY) = resaux                                                      
      END DO                                                                    
c ----------                                                                    
c     find parameter Psi (see Ivezic & Elitzur, 1996)                           
c     generate auxiliary function for integration:                              
c     loop over iL (wavelength)                                                 
      DO iL = 1, nL                                                             
        qaux(iL) = SigmaA(1,iL) * Utot(iL,1) / lambda (iL)                      
      END DO                                                                    
      CALL Simpson(npL,1,nL,lambda,qaux,resaux)                                 
      QUtot1 = resaux                                                           
      Psi = QUtot1 / QpTd(1,1)                                                  
c -------------                                                                 
      IF(denstyp.NE.0) THEN                                                     
c      ratio r1/r* (see Ivezic & Elitzur, 1996, eq. 27)                         
       r1rs = 0.5 * dsqrt(Psi) * (Tstar / Tsub(1))**2.0                         
      END IF                                                                    
c -------------                                                                 
c     Find epsilon - the relative contribution of the diffuse radiation         
      DO iY = 1, nY                                                             
       aux = QpStar(iY)/QpTd(1,iY)/Psi*(Tsub(1)/Td(1,iY))**4.                   
       IF (denstyp.NE.0) aux = aux/ Y(iY)/Y(iY)                                 
       Eps(iY) = 1. - aux                                                       
      END DO                                                                    
      Eps1 = 1.0 - QpStar(1) / QUtot1                                           
c     store these parameters in the storage array                               
      SmC(1,model) = Psi                                                        
      SmC(2,model) = Eps1                                                       
      SmC(3,model) = QpStar(1)                                                  
      SmC(4,model) = QpTd(1,1)                                                  
      SmC(5,model) = maxFerr                                                    
c -------------                                                                 
c     additional output quantities                                              
c     bolometric flux at r1 (in W/m2)                                           
c     The constant is 4*Sigma*T^4 (2.27E5 = 4*5.67E-08*(1000**4))               
      F1 = 2.27E5 / Psi * (Tsub(1)/1000.)**4.0                                  
c     inner radius (in cm)                                                      
      Cr1 = 5.53E16 / dsqrt(F1)                                                 
      IF (denstyp.NE.0) THEN                                                    
c      angular diameter of inner cavity if Fbol=1E-6 W/m2                       
       theta1 = 412.6 / dsqrt(F1)                                               
c      check if the pt.source assumption is still obeyed                        
c      (only for BB-type spectrum including EM-function)
       IF(startyp(1).eq.1.OR.startyp(1).eq.2) THEN                              
         mx = sqrt(sqrt(F1/5.67E-08))                                           
         Te_min = 2. * DMAX1(Tsub(1), mx)                                       
       END IF                                                                   
      END IF                                                                    
      IF (denstyp.EQ.0) THEN                                                    
c       Teff for the left illuminating source in slab geometry                  
c       Teff = (F1/sigma)^0.25                                                  
        SmC(7,model) = sqrt(sqrt(F1/5.67E-08))                                  
        IF (ksi.GT.0.) THEN                                                     
c       Teff for the right illuminating source in slab geometry                 
           SmC(8,model) = SmC(7,model)*sqrt(sqrt(ksi))                          
        ELSE                                                                    
           SmC(8,model) = 0.                                                    
        END IF                                                                  
      END IF                                                                    
c     calculate conversion constants for dynamics                               
      IF (denstyp.eq.4.OR.RDW) THEN                                    
        IF (denstyp.EQ.4) THEN                                                  
          I1 = 2. * (1.-pow)/(1.+pow)/tauF(nY)                                  
          I2 = I1                                                               
          I3 = I1 * tauF(nY) / TAUfid                                           
          Gamma(nY) = 0.5                                                       
        END IF                                                                  
c       terminal expansion velocity, full formula:                              
        ugas_out = tauF(nY) * (1.-Gamma(nY)) / (1.-pow)                         
c       The coefficients come from the units conversion                         
        C1 = 0.2845*TAUfid*sqrt(Psi)/I2/(SigExfid/aveV)*                        
     &                                        1.0D06/Tsub(1)/Tsub(1)            
        C2 = 2.040*ugas_out                                                     
        C3 = 6.628*I3*SigExfid/aveV*Gamma(nY)/I1                                
c       from version 2.0 stellar mass is defined as the maximal stellar         
c       mass which does not quench the wind; the calculation is done            
c       with half that mass since any smaller mass will have no effect          
c       on the radial velocity and density profile (see IE '99)                 
c       n.b. Gamma(nY) is removed                                               
        C3 = 6.628*I3*SigExfid/aveV/I1                                          
c       new definitions for output                                              
c       mass-loss rate in Msol/yr                                               
        CMdot = 1.0E-05 * sqrt(C1)                                              
c       terminal expansion velocity in km/s                                     
        Cve = 10.* C2 / sqrt(C1)                                                
c       stellar mass                                                            
        CM = C3                                                                 
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE AssProp(NN,nin,kin,nout,kout)                                  
c =======================================================================       
c This subroutine copies arrays nin and kin to nout and kout,                   
c respectively.                                        [Z.I., Nov. 1996]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER NN, i                                                             
      DOUBLE PRECISION nin(NN), kin(NN), nout(NN), kout(NN)                     
c -----------------------------------------------------------------------       
      DO i = 1, NN                                                              
        nout(i) = nin(i)                                                        
        kout(i) = kin(i)                                                        
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE BOLOM(q,qBol)                                                  
c =======================================================================       
c This subroutine integrates given radiation field, q (function of              
c wavelength and radial coordinate), over wavelength. q is a matrix             
c of physical size (npL,npY) [coming from paramet.inc] and real size            
c (nL,nY) [coming from grids.inc], and qBol is an array of physical size        
c (npY) and real size nY.                              [Z.I., Mar. 1996]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER iL, iY                                                            
      DOUBLE PRECISION q(npL,npY), qaux(npL), qBol(npY), resaux                 
c -----------------------------------------------------------------------       
c    loop over iY (radial coordinate)                                           
      DO iY = 1, nY                                                             
c       generate auxiliary function for integration                             
c       loop over iL (wavelength)                                               
        DO iL = 1, nL                                                           
          qaux(iL) = q(iL,iY) / lambda (iL)                                     
        END DO                                                                  
        CALL Simpson(npL,1,nL,lambda,qaux,resaux)                               
        qBol(iY) = resaux                                                       
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE ChkBolom(qbol,accur,dev,FbolOK)                                
c =======================================================================       
c This subroutine checks if any element of qbol(i), i=1,nY differs for          
c more than accuracy from 1.0. If so FbolOK = 0, otherwise FbolOK = 1.          
c dev is maximal deviation from 1.0.                   [Z.I., Mar. 1996]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER iY, FbolOK                                                        
      DOUBLE PRECISION qBol(npY), accur, dev                                    
c -----------------------------------------------------------------------       
      FbolOK = 1                                                                
      dev = 0.0                                                                 
c     loop over iY (radial coordinate)                                          
      DO iY = 1, nY                                                             
        IF (abs(1.0-qBol(iY)).GT.accur) FbolOK = 0                              
        IF (abs(1.0-qBol(iY)).GT.dev) dev = abs(1.0-qBol(iY))                   
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE ChkConv(nY,accuracy,Aold,Anew,Aconv)                           
c =======================================================================       
c This subroutine checks convergence of an array A(nY) between values           
c given in Aold and Anew. If the relative difference for EVERY element          
c is smaller than accuracy, Aconv is assigned 1, otherwise 0.                   
c                                                      [Z.I., Jul. 1996]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, iY, Aconv                                                     
      DOUBLE PRECISION accuracy, Aold(npY), Anew(npY), delta                    
c -----------------------------------------------------------------------       
      Aconv = 1                                                                 
c     loop over radial positions                                                
      DO iY = 1, nY                                                             
c       find relative difference                                                
        delta = dabs(Anew(iY)-Aold(iY))                                         
        IF (delta.GT.dabs(Anew(iY))*accuracy) Aconv = 0                         
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE ChkFlux(flux,tolern,consfl,error,ETAzp)                        
c =======================================================================       
c Checks the bolometric flux conservation at any point of a given Ygrid.        
c In case of nonconservation increases the number of points at certain          
c places. The current criterion is increasing the flux difference from          
c tolern to its maximum value.                         [MN & ZI,July'96]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER iYins(npY), k, kins, iY, consfl, flag, error, istop, iDm          
      DOUBLE PRECISION Yins(npY),flux(npY), tolern, ff, Yloc,                   
     &       delTAUMax, devfac, devmax, ffold, EtaTemp(npY),ee,                 
     &       ETA, ETAzp(npP,npY)                                                
c -----------------------------------------------------------------------       
c     save old grid and values of Eta (important for denstyp = 5 or 6)          
      IF (RDW) THEN                                                             
        DO iY = 1, nY                                                           
          Yprev(iY) = Y(iY)                                                     
          EtaTemp(iY) = ETAdiscr(iY)                                            
        END DO                                                                  
        nYprev = nY                                                             
      END IF                                                                    
      error = 0                                                                 
      kins = 0                                                                  
      devmax = 0.0                                                              
c     maximal delTAU is no more than 2 times the average value                  
      delTAUmax = 2.0*TAUtot(1)*ETAzp(1,nY)/nY                                  
c     maximal deviation from 1.0                                                
      DO iY = 2, nY                                                             
        IF (abs(flux(iY)-1.0).GT.devmax) devmax = abs(flux(iY)-1.0)             
      END DO                                                                    
      ff = 0.0                                                                  
      istop = 0                                                                 
      devfac = 0.1                                                              
c     search for places to improve the grid                                     
      DO WHILE (istop.NE.1)                                                     
        DO iY = 2, nY                                                           
          ffold = ff                                                            
          ff = abs(flux(iY)-1.0)                                                
          flag = 0                                                              
c         if any of these criteria is satisfied insert a point:                 
c         1) if error is increasing too fast                                    
          IF (abs(ff-ffold).GT.devfac*devmax) flag = 1                          
c         2) if delTAU is too large                                             
          IF (TAUtot(1)*(ETAzp(1,iY)-ETAzp(1,iY-1)).GT.                         
     &                             delTAUmax) flag = 1                          
          IF(flag.EQ.1.AND.devmax.GE.tolern) THEN                               
            kins = kins + 1                                                     
            Yins(kins) = Y(iY-1)+0.5*(Y(iY)-Y(iY-1))                            
            iYins(kins) = iY-1                                                  
          END IF                                                                
        END DO                                                                  
        IF (devmax.LT.tolern.OR.devfac.LT.0.01) THEN                            
          istop = 1                                                             
          ELSE                                                                  
          IF (kins.GT.0) istop = 1                                              
        END IF                                                                  
        devfac = devfac / 2.0                                                   
      END DO                                                                    
      IF (kins.EQ.0) THEN                                                       
        IF (consfl.NE.5) consfl = 1                                             
        ELSE                                                                    
c       Add all new points to Y(nY). This gives the new Y(nY+kins).             
c       However, check if npY is large enough to insert all points:             
        IF ((nY+kins).GT.npY) THEN                                              
c        consfl.EQ.5 is a signal that Chkflux was called from SetGrids,         
c        in this case continue without inserting new points. If this is         
c        full problem then give it up.                                          
         IF (consfl.NE.5) THEN                                                  
           consfl = 1                                                           
           ELSE                                                                 
           consfl = 7                                                           
           goto 777                                                             
         END IF                                                                 
         IF (iX.GE.1) THEN                                                      
         write(18,*)' ****************     WARNING   ******************'        
         write(18,*)'  The new Y grid can not accomodate more points!'          
         write(18,'(a,i3)')'   Specified accuracy would require',nY+kins        
         write(18,'(a,i3,a)')'   points, while npY =',npY,'.'                   
         write(18,*)'  For the required accuracy npY must be increased,'        
         write(18,*)'  (see the manual S3.5 Numerical Accuracy).'               
         write(18,*)' *************************************************'        
         END IF                                                                 
         kins = npY - nY                                                        
         iWARNING = iWARNING + 1                                                
         error = 2                                                              
        END IF                                                                  
        DO k = 1, kins                                                          
          CALL SHIFT(Y,npY,nY+k-1,Yins(k),iYins(k)+k-1)                         
        END DO                                                                  
      END IF                                                                    
c     new size of the Y grid                                                    
      nY = nY + kins                                                            
c     intepolate ETAdiscr to new Y grid for denstyp = 5 or 6                    
      DO iY = 1, nY                                                             
        Yloc = Y(iY)                                                            
        IF (iterETA.GT.1) THEN                                                  
          CALL LinInter(npY,nYprev,Yprev,EtaTemp,Yloc,iDm,ee)                   
          ETAdiscr(iY) = ee                                                     
          ELSE                                                                  
          ETAdiscr(iY) = ETA(Yloc)                                              
        END IF                                                                  
      END DO                                                                    
c -----------------------------------------------------------------------       
777   RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE Converg1(nY,accuracy,dynrange,Aold,Anew,Aconv,dmax)            
c =======================================================================       
c This subroutine checks convergence of an array A(nL,nY) between values        
c given in Aold and Anew, when the values are larger than dynrange. If          
c the maximum relative difference is smaller than the required accuracy,        
c Aconv is assigned 1, otherwise 0.              [Z.I.Jul 96;M.N.Apr.97]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, iY, Aconv                                                     
      DOUBLE PRECISION accuracy, dynrange, Aold(npY), Anew(npY), delta,         
     &       dmax                                                               
c -----------------------------------------------------------------------       
      Aconv = 1                                                                 
      dmax = 0.0                                                                
c     loop over radial positions                                                
      DO iY = 1, nY                                                             
c       do it only for elements larger than dynrange                            
        IF (Anew(iY).GE.dynrange) THEN                                          
c           find relative difference                                            
            delta = dabs((Anew(iY)-Aold(iY))/Anew(iY))                          
            IF (delta.GT.dmax) dmax = delta                                     
        END IF                                                                  
      END DO                                                                    
      IF (dmax.GT.accuracy) Aconv = 0                                           
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE Converg2(nY,nL,accuracy,dynrange,Aold,Anew,Aconv,dmax)         
c =======================================================================       
c This subroutine checks convergence of an array A(nL,nY) between values        
c given in Aold and Anew, when the values are larger than dynrange. If          
c the maximum relative difference is smaller than required accuracy,            
c Aconv is assigned 1, otherwise 0.             [Z.I.Jul 96; M.N.Apr.97]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nL, iY, iL, Aconv                                             
      DOUBLE PRECISION accuracy, dynrange, Aold(npL,npY), Anew(npL,npY),        
     &       delta, dmax                                                        
c -----------------------------------------------------------------------       
      Aconv = 1                                                                 
      dmax = 0.0                                                                
c     loop over wavelengths                                                     
      DO iL = 1, nL                                                             
c       loop over radial positions                                              
        DO iY = 1, nY                                                           
c         do it only for elements larger than dynrange                          
          IF (Anew(iL,iY).GE.dynrange) THEN                                     
c           find relative difference                                            
            delta = dabs((Anew(iL,iY)-Aold(iL,iY))/Anew(iL,iY))                 
            IF (delta.GT.dmax) dmax = delta                                     
          END IF                                                                
        END DO                                                                  
      END DO                                                                    
      IF (dmax.GT.accuracy) Aconv = 0                                           
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE Convolve(IntOut)                                               
c =======================================================================       
c This subroutine convolves intensity IntOut with the point spread              
c function to produce convolved images ConvInt. The work horse is               
c subroutine Conv2D, and this subroutine is used to prepare everything.         
c                                                      [Z.I., Jan. 1997]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER psftype, Npsf, iLambda                                            
      DOUBLE PRECISION kPSF(20), FWHM1(20), FWHM2(20), theta1,                  
     &       xpsf(1000), ypsf(1000)                                             
      COMMON /psf1/ iLambda, psftype, Npsf                                      
      COMMON /psf2/ kPSF, FWHM1, FWHM2, Theta1,                                 
     &       xpsf, ypsf                                                         
      INTEGER i, j                                                              
      DOUBLE PRECISION IntOut(20,npP+2), yang(npP+2), Youtang, deltaOff,        
     &       Int1D(npP+2),  ConvS, Conv(npP+2), FWHM1max, FWHM2max, PSFN        
c -----------------------------------------------------------------------       
c     find the largest FWHMs                                                    
      FWHM1max = FWHM1(1)                                                       
      FWHM2max = FWHM2(1)                                                       
      IF (psftype.LT.3) THEN                                                    
        DO i = 1, NlambdaOut                                                    
          IF (FWHM1(i).GT.FWHM1max) FWHM1max = FWHM1(i)                         
          IF (FWHM2(i).GT.FWHM2max) FWHM2max = FWHM2(i)                         
        END DO                                                                  
      END IF                                                                    
c     scale angular coordinate to theta1                                        
      DO i = 1, nP+2                                                            
        yang(i) = bOut(i) * Theta1 / 2.                                         
      END DO                                                                    
c     generate off-set grid                                                     
      Youtang = Y(nY) * Theta1                                                  
      IF (Youtang.GT.FWHM1max.AND.Youtang.GT.FWHM2max) THEN                     
c       the envelope is well resolved, take impact parameter grid               
        Nconv = nP + 2                                                          
        DO i = 1, Nconv                                                         
          Offset(i) = yang(i)                                                   
        END DO                                                                  
      ELSE                                                                      
c       the envelope is not well resolved, take equidistant grid                
c       to 2FWHM1max, i.e. image will be more or less the PSF itself            
        Nconv = 30                                                              
        deltaOff = 2.0 * FWHM1max / (Nconv-1)                                   
        IF (FWHM2max.GT.FWHM1max) THEN                                          
          deltaOff = 2.0 *FWHM2max / (Nconv-1)                                  
        END IF                                                                  
        DO i = 1, Nconv                                                         
          Offset(i) = deltaOff * 1.0*(i-1)                                      
        END DO                                                                  
      END IF                                                                    
c     convolve intensity wavelength by wavelength                               
      DO j = 1, NlambdaOut                                                      
c       needed to specify wavelength in PSFN                                    
        iLambda = j                                                             
c       generate 1D intensity vector for subrutine Conv2D                       
c       take only diffuse emission, stellar contribution will                   
c       be added below (a shortcut to avoid inaccuracies or too many            
c       points in Conv2D)                                                       
        DO i = 1, nP+2                                                          
          IF (i.LE.2) THEN                                                      
            Int1D(i) = IntOut(j,3)                                              
          ELSE                                                                  
            Int1D(i) = IntOut(j,i)                                              
          END IF                                                                
          CALL CHKRANGE(dynrange,Int1D(i))                                      
          IF (Int1D(i).LT.dynrange) Int1D(i)=0.0                                
        END DO                                                                  
c       convolve                                                                
        CALL Conv2D(npP+2,nP+2,yang,Int1D,1000,NConv,Offset,Conv)               
c       add stellar contribution                                                
        DO i = 1, nP+2                                                          
          ConvS=2.*ASIN(1.0)*(yang(2)**2.)*IntOut(j,1)*PSFN(Offset(i))          
          Conv(i) = Conv(i) + ConvS                                             
        END DO                                                                  
c       scale to 1 at the center                                                
        CALL ScaleTo1(1000,Nconv,Conv)                                          
c       copy 1D convolved intensity to ConvInt                                  
        DO i = 1, Nconv                                                         
          CALL CHKRANGE(dynrange,Conv(i))                                       
          ConvInt(j,i) = Conv(i)                                                
        END DO                                                                  
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE Conv2D(NinMax,Nin,Xin,Yin,Noutmax,Nout,Xout,Yout)              
c =======================================================================       
c This subroutine convolves intensity Yin(Xin[i]), i=1,Nin with                 
c the point spread function PSFN(x) (provided as a separate function).          
c It is assumed that both the intensity yin and PSFN(x) are circularly          
c symmetric functions of radial coordinate x, i.e., this subroutine             
c performs two-dimensional convolution. Convolved intensity, Yout, is           
c evaluated for very position Xout[i], i=1,Nout, as:                            
c        Yout(Xout) = Int[Yin(yloc)*PSF(xloc)*yloc*dyloc*dphi]                  
c where xloc = sqrt(yloc **2+Xout**2-2*yloc*Xout*cos(phi), with yloc            
c and phi being dummy integration variables. Declared size of Xin is            
c NinMax, the one for Xout is NoutMax. The radial integration is done           
c using subroutine ROMBY and angular integration is done by using               
c Simpson rule.                                        [Z.I., Jan. 1997]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER NinMax, Nin, NoutMax, Nout, iphi, iXin, Nphi, iXOut               
      DOUBLE PRECISION Xin(NinMax), Yin(NinMax), Xout(NoutMax), A, B,           
     &       Yout(NoutMax), dphi, phi(1000), fphi(1000), int1, int2,            
     &       imagfn                                                             
      INTEGER psftype, Npsf, iLambda                                            
      DOUBLE PRECISION kPSF(20), FWHM1(20), FWHM2(20), theta1,                  
     &       xpsf(1000), ypsf(1000)                                             
      COMMON /psf1/ iLambda, psftype, Npsf                                      
      COMMON /psf2/ kPSF, FWHM1, FWHM2, Theta1,                                 
     &       xpsf, ypsf                                                         
      INTEGER ftype                                                             
      DOUBLE PRECISION Ckn, Cxout, Cphi, Cqtheta                                
      COMMON /imfn1/ ftype                                                      
      COMMON /imfn2/ Ckn, Cxout, Cphi, Cqtheta                                  
      EXTERNAL imagfn                                                           
c -----------------------------------------------------------------------       
c     Parameters for integration:                                               
c     number of angular points                                                  
      Nphi = 9                                                                  
c     step in angle phi                                                         
      dphi = 2.0*ASIN(1.0) / (Nphi-1)                                           
c     flag for imgfn                                                            
      ftype = 1                                                                 
c     Start integrations                                                        
c     loop over output positions                                                
      DO iXout = 1, Nout                                                        
        Cxout = Xout(iXout)                                                     
c       loop over angular wedges (phi integration)                              
        DO iphi = 1, Nphi                                                       
          phi(iphi) = dphi*1.0*(iphi-1)                                         
          Cphi = phi(iphi)                                                      
          fphi(iphi) = 0.0                                                      
c         loop over input radial positions (radial integration)                 
          DO iXin = 1, Nin-1                                                    
            Ckn = 1.0                                                           
            CALL ROMBY(imagfn,Xin(iXin),Xin(iXin+1),int1)                       
            Ckn = 2.0                                                           
            CALL ROMBY(imagfn,Xin(iXin),Xin(iXin+1),int2)                       
c           contribution from this annulus (lin. approx. for intensity)         
            A = Xin(iXin+1)*Yin(iXin) - Xin(iXin)*Yin(iXin+1)                   
            A = A / (Xin(iXin+1)-Xin(iXin))                                     
            B = (Yin(iXin+1)-Yin(iXin)) / (Xin(iXin+1)-Xin(iXin))               
            fphi(iphi) = fphi(iphi) + A*int1 + B*int2                           
          END DO                                                                
        END DO                                                                  
        CALL Simpson(1000,1,Nphi,phi,fphi,Yout(iXout))                          
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      DOUBLE PRECISION FUNCTION ETA(Yy)                                         
c =======================================================================       
c This subroutine evaluates the normalized density profile. denstype is         
c a type of density law: 1 and 2 - power law, 3 for exponential density law,    
c 4 for radiatively driven winds (the gray-body approximation), 5,6 - RDW and   
c 7 - d.d.from a file. pow is parameter describing the choosen density law:     
c power for 1 and 2, v1/v8 for RDW (the ratio of expansion velocities at the    
c inner and outer radii), sigma for 3 [i.e. rho = dexp(-(y/sigma)**2)].         
c Yout is the relative thickness, Yout=rout/r1. Y is the radial position.       
c                                                          [ZI'95; ZI'99]       
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER i, istop, iYdummy                                                 
      DOUBLE PRECISION Yy, C, IntAux, powaux, Prod, eps, factEta                
c -----------------------------------------------------------------------       
c     this is an adjustable value regulating the initial Eta approximation for d
      factEta = 0.5                                                             
      IF (Yy.GT.Yout) Yy = Yout                                                 
      IF (iterETA.GT.1) THEN                                                    
c     if this is iteration over ETA (but not the first one) in the case         
c     of dynamical calculation for radiatively driven winds calculate           
c     ETA by linear interpolation of ETA from the previous iteration            
      CALL LinInter(npY,nYprev,Yprev,ETAdiscr,Yy,iYdummy,ETA)                   
c     otherwise use prescribed formulae for different cases                     
      ELSE                                                                      
c     smooth power-law                                                          
      IF (denstyp.EQ.1) THEN                                                    
c       find normalization constant                                             
        IF (pow.NE.1.0) THEN                                                    
          C = (1.0 - Yout**(1.0-pow)) / (pow - 1.0)                             
          ELSE                                                                  
          C = dlog(Yout)                                                        
        ENDIF                                                                   
        IF (Ntr.GE.1) THEN                                                      
          DO i = 1, Ntr                                                         
            powaux = pow - ptr(i)                                               
            IF (powaux.NE.1.0) THEN                                             
              IntAux = (1.0 - Yout**(1.0-powaux)) / (powaux - 1.0)              
            ELSE                                                                
              IntAux = dlog(Yout)                                               
            END IF                                                              
            C = C + IntAux / Ytr(i)**ptr(i)                                     
          END DO                                                                
        ENDIF                                                                   
        C = 1.0 / C                                                             
c       calculate density                                                       
        IF (Yy.GE.(1.0-1.0E-8)) THEN                                            
          ETA = C / Yy**pow                                                     
          IF (Ntr.GE.1) THEN                                                    
            DO i = 1, Ntr                                                       
              ETA = ETA + C * Yy**(ptr(i)-pow) / Ytr(i)**ptr(i)                 
            END DO                                                              
          END IF                                                                
          ELSE                                                                  
          ETA = 0.0                                                             
        ENDIF                                                                   
      END IF                                                                    
c     broken power-law                                                          
      IF (denstyp.EQ.2) THEN                                                    
        Ytr(Ntr+1) = Yout                                                       
c       find normalization constants                                            
        IF (pow.NE.1.0) THEN                                                    
          C = (1.0 - Ytr(1)**(1.0-pow)) / (pow - 1.0)                           
          ELSE                                                                  
          C = dlog(Ytr(1))                                                      
        ENDIF                                                                   
        IF (Ntr.GE.1) THEN                                                      
          DO i = 1, Ntr                                                         
            CALL DoProduct(10,Ytr,ptr,pow,i,Prod)                                 
            IF (ptr(i).NE.1.0) THEN                                             
              IntAux = Ytr(i)**(1.0-ptr(i)) - Ytr(i+1)**(1.0-ptr(i))            
              IntAux = Prod * IntAux / (ptr(i) - 1.0)                           
              ELSE                                                              
              IntAux = Prod * dlog(Ytr(i+1)/Ytr(i))                             
            END IF                                                              
            C = C + IntAux                                                      
          END DO                                                                
        ENDIF                                                                   
        C = 1.0 / C                                                             
c       calculate density                                                       
        IF (Yy.GE.1.0-1.0E-8) THEN                                              
          IF (Yy.LE.Ytr(1)) THEN                                                
            ETA = C / Yy**pow                                                   
          ELSE                                                                  
            istop = 0                                                           
            i = 0                                                               
            DO WHILE (istop.NE.1)                                               
              i = i + 1                                                         
              IF (Yy.LE.Ytr(i+1)) istop = 1                                     
            END DO                                                              
            CALL DoProduct(10,Ytr,ptr,pow,i,Prod)                                 
            ETA = C * Prod / Yy**ptr(i)                                         
          END IF                                                                
        ELSE                                                                    
          ETA = 0.0                                                             
        ENDIF                                                                   
      END IF                                                                    
c     exponential law                                                           
      IF (denstyp.EQ.3) THEN                                                    
        IF (Yy.GE.1.0-1.0E-8) THEN                                              
          ETA = (Yout-1.) * (1.-exp(-pow)) / pow                                
          ETA = exp(-pow*(Yy-1.)/(Yout-1.)) / ETA                               
        ELSE                                                                    
          ETA = 0.0                                                             
        END IF                                                                  
      END IF                                                                    
c     radiatively driven winds (the gray-body approximation)                    
      IF (denstyp.EQ.4) THEN                                                    
        eps = pow                                                               
        IF (Yy.GE.1.0-1.0E-8) THEN                                              
          ETA = (1.+eps)/2./Yy/Yy/sqrt(1.-(1.-eps*eps)/Yy)                      
        ELSE                                                                    
          ETA = 0.0                                                             
        ENDIF                                                                   
      END IF                                                                    
c     radiatively driven winds (full calculation)                               
       IF (RDW) THEN                                                            
c       if this is the first iteration use analytic approximation               
        IF (iterETA.LT.2) THEN                                                  
c         for 5 eps is pow, for 6 assume eps=0.1                                
          IF (denstyp.EQ.5) THEN                                                
            eps = pow                                                           
          ELSE                                                                  
            eps = 0.1                                                           
          END IF                                                                
          IF (Yy.GE.(1.0-1.0E-8)) THEN                                          
            ETA = (1.+eps)/2./Yy/Yy/sqrt(1.-(1.-eps*eps)/Yy)                    
c           empirical improvement for the initial approximation                 
c           good only for large optical depths, but small ones                  
c           are fast anyway (ZI, May99)                                         
            IF (Yy.LE.2) THEN                                                   
               ETA = ETA / (1. + factEta / Yy**10.)                             
            END IF                                                              
          ELSE                                                                  
            ETA = 0.0                                                           
          ENDIF                                                                 
c       or interpolate from the previous solution                               
        ELSE                                                                    
          CALL LinInter(npY,nYprev,Yprev,ETAdiscr,Yy,iYdummy,ETA)               
        END IF                                                                  
      END IF                                                                    
c     user specified function                                                   
      IF (denstyp.EQ.7) THEN                                                    
        IF (Yy.LT.yEta7(nYEta7)) THEN                                           
          CALL LinInter(npY,nYEta7,yEta7,Eta7,Yy,iYdummy,ETA)                   
        ELSE                                                                    
          ETA = Eta7(nYEta7)                                                    
        END IF                                                                  
      END IF                                                                    
c     done                                                                      
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
                                                                                
c ***********************************************************************       
      DOUBLE PRECISION FUNCTION ETAfun(iL,Y)                                    
c =======================================================================       
c This subroutine evaluates the normalized density profile as a function        
c of position and wavelength (in multigrain case). The MAIN purpose of          
c this function is to provide connection between the overall, prescribed        
c density distribution, ETA, and wavelength depedent function ETAfun.           
c It is NOT FINISHED, this is a trivial case, works only for single             
c grains.                                              [Z.I., Nov. 1995]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER iL                                                                
      DOUBLE PRECISION Y, eta                                                   
c -----------------------------------------------------------------------       
c     this is temporary, works only for single component grains                 
c     I need to add subroutines/functions to take care of temperature           
c     dependent and thus variable condensation radii for different              
c     components                                                                
      ETAfun = ETA(Y)                                                           
c     this is to avoid warning when compiling:                                  
      iL = iL                                                                   
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE FindInt(nG,alpha,ETAzp)                                        
c =======================================================================       
c This subroutine finds the intensity distribution at outer edge and for        
c user specified wavelengths lamOut. It also evaluates the angular size         
c of the stellar disk and adds two impact parameters describing the star        
c to the P grid, thus producing bOut grid. All intensities are indeed           
c dimensionless quantities lambda*I_lambda/F1 where I_lambda is real            
c physical quantity defined as usual and F1 is the bolometric flux at           
c the dust sublimation radius, r1. For conversion to the physical value         
c lambda*I_lambda, I_lambda from the program has to be multiplied by F1.        
c F1 can obtained either as:                                                    
c      1) F1 = 4*sigma*Tsub**4/Psi (IE96, eq. 15),                              
c where Tsub is sublimation temperature and parameter Psi is given in           
c *.SUM file; or as:                                                            
c      2) F1 = Fbol/alpha1**2 (IE96, eq. 34)                                    
c where Fbol is the bolometric flux and alpha1 is the angular size of r1        
c at any particular distance from the envelope (i.e. both F1 and alpha1         
c correspond to observed quantities). Also note that                            
c     INT(I_lambda(p)*2Pi*P*dP) = f_lambda                                      
c where I_lambda is the scaled quantity from the program, P is impact           
c parameter, and f_lambda is the spectral shape F_lambda/Fbol.                  
c                                                      [Z.I., Aug. 1996]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER nYok, nPok                                                        
      DOUBLE PRECISION                                                          
     &       Ude(npL,npY), Uds(npL,npY), Us(npL,npY), fde(npL,npY),             
     &       fds(npL,npY), fs(npL,npY), Utot(npL,npY), ftot(npL,npY),           
     &       Td(npG,npY), Ubol(npY), fbol(npY), Uchck(npL,npY),                 
     &       UbolChck(npY), Spectrum(npL), SpecChar(30,99), tauF(npY),          
     &       Intens(npL,npP+2), IntOut(20,npP+2), SmC(30,99),                   
     &       Yok(npY), Pok(npP), tauOut(npL), tauZout(npP+2), F1,               
     &       tr(npY), rg(npG,npY), fsL(npL,npY), fsR(npL,npY), Eps(npY)         
      COMMON /solution/ Ude, Uds, Us, fde, fds, fs, Utot, ftot, Td,             
     &       Ubol, fbol, Uchck, UbolChck, Spectrum, SpecChar, tauF,             
     &       Intens, IntOut, SmC, Yok, Pok, tauOut, tauZout, F1,                
     &       tr, rg, fsL, fsR, Eps, nYok, nPok                                  
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER iYfirst, YPequal, Plast                                           
      DIMENSION iYfirst(npP), YPequal(npP), Plast(npY)                          
      COMMON /Yfirst/ iYfirst, YPequal, Plast                                   
      INTEGER startyp(2), Nlamtr(2), nBB(2)                                     
      CHARACTER nameStar(2)*235                                                 
      DOUBLE PRECISION Tstar, lamtr(2,101), klam(2,100), Tbb(2,10),             
     &                 rellum(2,10), mu1, ksi, mu2, xSiO, r1rs                  
      COMMON /source/ Tstar, lamtr, klam, Tbb, rellum, mu1, ksi, mu2,           
     &                xSiO, r1rs, startyp, Nlamtr, nBB, nameStar                
      INTEGER iL, nG, k, i, iLout, iLstop, iP, iW, nZ, Nzpt, iZ, izloc          
      DOUBLE PRECISION   qaux(npL), Psi,alpha(npG,npY), resaux, QUtot1,         
     &       QpTsub, xP, pst, stelfact, Istell(npL), Planck, w1,                
     &       IntL, IntR, xx , alb, Ids(npL,npP), Ide(npL,npP),w2, lw12,         
     &       ETAzp(npP,npY), numcorr, ETAzpStar, qaux2(npL), z1, z2,            
     &       delz, zloc, wloc, resint, pT, Tz, Idboth, tzp(100), pUtot,         
     &       Semis(100), Sscat(100), IntETA, palb, palf, alfa, exterm,          
     &       Utotloc, Sstem(100), Sstsc(100), Istem(npL,100), Idfront,          
     &       Istsc(npL,100), delTau, factaux, UtotL, UtotR, ep1,                
     &       tauzp1, tauInf                                                     
c -----------------------------------------------------------------------       
c     temporary                                                                 
      IF (nG.GT.1.AND.iX.GE.1) THEN                                             
        write(18,*)' FindInt should be fixed, nG>1 !'                           
        stop                                                                    
      END IF                                                                    
c     find impact parameter tangential to the stellar disk                      
c     first find the Planck averaged absorption efficiencies at Y=1             
      DO iL = 1, nL                                                             
        qaux(iL) = SigmaA(1,iL) * Utot(iL,1) / lambda (iL)                      
        xP = 14400.0 / Tsub(1) / lambda(iL)                                     
        qaux2(iL) = SigmaA(1,iL) * Planck(xP) / lambda (iL)                     
      END DO                                                                    
      CALL Simpson(npL,1,nL,lambda,qaux,resaux)                                 
      QUtot1 = resaux                                                           
      CALL Simpson(npL,1,nL,lambda,qaux2,resaux)                                
      QpTsub = resaux                                                           
c     parameter Psi (see Ivezic & Elitzur, 1996, eq. C4)                        
      Psi = QUtot1 / QpTsub                                                     
c     ratio pst = rstar/rsub (see Ivezic & Elitzur, 1996, eq. 27)               
      pst = 2.0 / dsqrt(Psi) * (Tsub(1) / Tstar)**2.0                           
      IF (pst.GE.0.5) THEN                                                      
         IF (iX.GE.1) THEN                                                      
            write(18,*)' FindInt: specified dust temperature at the '           
            write(18,*)' inner radius results in r*/r1 >= 0.5: '                
            write(18,*)'    r*/r1 =', pst                                       
            write(18,*)' This violates some of Dusty`s assumptions'             
            write(18,*)'  ------  Please consult the manual ------'             
            write(18,*)'  ####  r*/r1 changed by hand to 0.5  ####'             
         END IF                                                                 
         pst = 0.5                                                              
      END IF                                                                    
      stelfact = 1.0 / pst / pst / 3.141593                                     
c     generate bOut, i.e. insert two points such that                           
c     bOut(k)=0.999*pst and bOut(k+1)=1.001*pst                                 
      CALL GetbOut(npP,nP,P,pst,bOut,k)                                         
c     correction for numerical errors in tau                                    
      numcorr = 1. / TAUtot(1)                                                  
c     loop over wavelengths                                                     
      DO iL = 1, nL                                                             
c        stellar intensity, Istell (extinction already included)                
         Istell(iL) = fs(iL,nY) * stelfact                                      
c        total optical depth along a line of sight                              
         tauOut(iL) = numcorr*TAUtot(iL)                                        
      END DO                                                                    
c     generate diffuse intensities, Ide (emission) and Ids (scat)               
c     loop over wavelengths                                                     
      DO iL = 1, nL                                                             
        DO iP = 1, nP                                                           
c         maximal number of points along tangential position, z                 
          nZ = nY + 1 - iYfirst(iP)                                             
c         starting value for local radius                                       
          IF (P(iP).GE.1.0) THEN                                                
            w2 = P(iP)                                                          
          ELSE                                                                  
            w2 = 1.0                                                            
          END IF                                                                
c         initialize intensities                                                
          Ide(iL,iP) = 0.0                                                      
          Ids(iL,iP) = 0.0                                                      
          IF (iP.LE.k+1) THEN                                                   
            Istem(iL,iP) = 0.0                                                  
            Istsc(iL,iP) = 0.0                                                  
          END IF                                                                
c         total optical depth along this impact parameter                       
          ep1 = ETAzp(iP,nY+1-iYfirst(iP))*TAUtot(iL)                           
c         loop over z, i.e. steps over points crossing the y grid               
          DO iZ = 2, nZ                                                         
c           index for the ending local radius                                   
            iW = iYfirst(iP) + iZ - 1                                           
c           local boundary radii                                                
            w1 = w2                                                             
            w2 = Y(iW)                                                          
c           corresponding displacements along a line of sight                   
            z1 = sqrt(abs(w1**2.-P(iP)**2.))                                    
            z2 = sqrt(abs(w2**2.-P(iP)**2.))                                    
c           # of pts. for z integration, should increase with deltaTau          
c           it is messy because INT function which would do the job is          
c           not in F77 standard set                                             
            Nzpt = 5                                                            
            delTau = (ETAzp(iP,iW)-ETAzp(iP,iW-1))*TAUtot(iL)                   
            IF (delTau.GT.1) Nzpt = 10                                          
            IF (delTau.GT.5) Nzpt = 20                                          
            IF (delTau.GT.10) Nzpt = 30                                         
            IF (delTau.GT.20) Nzpt = 40                                         
            IF (delTau.GT.50) Nzpt = 50                                         
            delz = (z2-z1) / (Nzpt-1)                                           
c           powers for power-law interpolations between 2 y pts.                
            lw12 = dlog(Y(iW-1)/Y(iW))                                          
c           for T                                                               
            pT = dlog(Td(1,iW)/Td(1,iW-1)) / lw12                               
c           for albedo                                                          
            IF (omega(iL,iW-1).GT.0.0.AND.omega(iL,iW).GT.0.0) THEN             
              palb = dlog(omega(iL,iW)/omega(iL,iW-1)) / lw12                   
            ELSE                                                                
              palb = 0.0                                                        
            END IF                                                              
c           for Utot                                                            
            UtotL = Utot(iL,iW-1)                                               
            UtotR = Utot(iL,iW)                                                 
            CALL CHKRANGE(dynrange,UtotL)                                       
            CALL CHKRANGE(dynrange,UtotR)                                       
            IF (UtotL.GT.0.0.AND.UtotR.GT.0) THEN                               
              pUtot = dlog(UtotR/UtotL) / lw12                                  
            ELSE                                                                
              pUtot = 0.0                                                       
            END IF                                                              
c           for alpha                                                           
            palf = dlog(alpha(1,iW)/alpha(1,iW-1)) / lw12                       
c           tauzp between z=0 and z=z1                                          
            tauzp1 = ETAzp(iP,iZ-1)*TAUtot(iL)                                  
c           integrate between adjacent grid points                              
            DO izloc = 1, Nzpt                                                  
              zloc = z1 + (izloc-1)*delz                                        
              wloc = sqrt(zloc**2 + P(iP)**2)                                   
c             find local TAUzp(w(z))-TAUzp(w1=w(z1))                            
              tzp(izloc) = IntETA(P(iP),iW-1,w1,wloc)*TAUtot(iL)                
c             find Tz = T(zloc) = T(wloc), this works for single                
c             size grains only; for multigrain case one needs to                
c             get Semis by summation over all Td                                
              Tz = Td(1,iW-1) * (Y(iW-1)/wloc)**pT                              
              xP = 14400/lambda(iL)/Tz                                          
c             power-law interpolation for albedo                                
              alb = omega(iL,iW-1) * (Y(iW-1)/wloc)**palb                       
c             power-law interpolation for Utot                                  
              IF (UtotL.GT.0) THEN                                              
                UtotLoc = UtotL * (Y(iW-1)/wloc)**pUtot                         
              ELSE                                                              
                UtotLoc = 0.0                                                   
              END IF                                                            
              CALL CHKRANGE(dynrange,UtotLoc)                                   
c             power-law interpolation for alpha                                 
              alfa = alpha(1,iW-1) * (Y(iW-1)/wloc)**palf                       
c             source functions (wloc**2 because D uses scaled quant.)           
              factaux = 1 / wloc**2 / (4 * 3.14159)                             
              Semis(izloc) = (1-alb) * alfa * Planck(xP) * factaux              
              Sscat(izloc) = alb * UtotLoc * factaux                            
c             check for the dynamic range                                       
              CALL CHKRANGE(dynrange,Semis(izloc))                              
              CALL CHKRANGE(dynrange,Sscat(izloc))                              
c             optical depth from infinity along the line of sight               
              tauInf = ep1 - tauzp1 - tzp(izloc)                                
c             for a line of sight terminating on the star find                  
c             contribution only from the front part of the envelope             
              IF (iP.LE.k+1) THEN                                               
                 IF (tauInf.LT.50) THEN                                         
                   exterm = dexp(-tauInf)                                       
                 ELSE                                                           
                   exterm = 0.0                                                 
                 END IF                                                         
                 Sstem(izloc) = Semis(izloc) * exterm                           
                 Sstsc(izloc) = Sscat(izloc) * exterm                           
              END IF                                                            
c             otherwise take both the front and back contributions              
              IF (tauInf.LT.50) THEN                                            
                 exterm = dexp(-tauInf)+dexp(-tauInf-ep1)                       
              ELSE                                                              
                 exterm = 0.0                                                   
              END IF                                                            
              Semis(izloc) = Semis(izloc) * exterm                              
              Sscat(izloc) = Sscat(izloc) * exterm                              
c             end of local loop over z                                          
            END DO                                                              
c           integrate and add contribution from this step                       
            CALL SIMPSON(100,1,Nzpt,tzp,Semis,resint)                           
            CALL CHKRANGE(dynrange,resint)                                      
            Ide(iL,iP) = Ide(iL,iP) + resint                                    
            CALL SIMPSON(100,1,Nzpt,tzp,Sscat,resint)                           
            CALL CHKRANGE(dynrange,resint)                                      
            Ids(iL,iP) = Ids(iL,iP) + resint                                    
            IF (iP.LE.k+1) THEN                                                 
              CALL SIMPSON(100,1,Nzpt,tzp,Sstem,resint)                         
              CALL CHKRANGE(dynrange,resint)                                    
              Istem(iL,iP) = Istem(iL,iP) + resint                              
              CALL SIMPSON(100,1,Nzpt,tzp,Sstsc,resint)                         
              CALL CHKRANGE(dynrange,resint)                                    
              Istsc(iL,iP) = Istsc(iL,iP) + resint                              
            END IF                                                              
c         end of loop over z                                                    
          END DO                                                                
c       end of loop over impact parameter, iP                                   
        END DO                                                                  
c     end of loop over wavelengths, iL                                          
      END DO                                                                    
c     add all intensities, Istell, Ide, Ids                                     
      DO iL = 1, nL                                                             
c       interpolate optical depth  at pstar                                     
        IF (iL.EQ.iLfid) THEN                                                   
          ETAzpStar = (ETAzp(k,nY) - ETAzp(k-1,nY))                             
          ETAzpStar = ETAzpStar * (pst-P(k-1)) / (P(k) - P(k-1))                
          ETAzpStar = ETAzp(k-1,nY) + ETAzpStar                                 
        END IF                                                                  
c       find diffuse contribution at pstar (by linear interpolation)            
        Idfront = Istsc(iL,k)+Istem(iL,k)-Istsc(iL,k-1)-Istem(iL,k-1)           
        Idfront = Idfront * (pst-P(k-1)) / (P(k) - P(k-1))                      
        Idfront = Idfront + Istsc(iL,k-1) + Istem(iL,k-1)                       
        Idboth = Ids(iL,k) + Ide(iL,k) - Ids(iL,k-1) - Ide(iL,k-1)              
        Idboth = Idboth * (pst-P(k-1)) / (P(k) - P(k-1))                        
        Idboth = Idboth + Ids(iL,k-1) + Ide(iL,k-1)                             
c       first for p<pstar, all three contributions                              
        DO i = 1, k-1                                                           
          Intens(iL,i) = Istell(iL) + Istsc(iL,i) + Istem(iL,i)                 
          IF (iL.EQ.iLfid)                                                      
     &        tauZout(i) = ETAzp(i,nY)/ETAzp(1,nY)                              
        END DO                                                                  
c       barely on the stellar disk                                              
        Intens(iL,k) = Istell(iL) + Idfront                                     
        tauZout(k) = ETAzpStar/ETAzp(1,nY)                                      
c       barely off the stellar disk                                             
        Intens(iL,k+1) = Idboth                                                 
        tauZout(k+1) = 2. * tauZout(k)                                          
c       all other p>pstar                                                       
        DO i = k, nP                                                            
          Intens(iL,i+2) = Ids(iL,i)+Ide(iL,i)                                  
          IF (iL.EQ.iLfid) THEN                                                 
            nZ = nY + 1 - iYfirst(i)                                            
            tauZout(i+2) = 2. * ETAzp(i,nZ)/ETAzp(1,nY)                         
          END IF                                                                
        END DO                                                                  
      END DO                                                                    
c     check dynamic range                                                       
      DO iL = 1, nL                                                             
        DO i = 1, nP+2                                                          
          CALL CHKRANGE(dynrange,Intens(iL,i))                                  
        END DO                                                                  
      END DO                                                                    
c     now interpolate Intens(lambda) to lamOut                                  
      DO iLout = 1, NlambdaOut                                                  
c       bracket the needed wavelength                                           
        iLstop = 0                                                              
        iL = 0                                                                  
        DO WHILE (iLstop.EQ.0)                                                  
          iL = iL + 1                                                           
          IF (lambda(iL).GT.LambdaOut(iLout)) iLstop = 1                        
          IF (iL.EQ.nL) iLstop = 1                                              
        END DO                                                                  
c       interpolate intensity                                                   
        xx = (LambdaOut(iLout)-lambda(iL-1))/(lambda(iL)-lambda(iL-1))          
        DO i = 1, nP+2                                                          
          IntL = Intens(iL-1,i)                                                 
          IntR = Intens(iL,i)                                                   
          IntOut(iLout,i) = IntL + xx*(IntR - IntL)                             
          CALL CHKRANGE(dynrange,IntOut(iLout,i))                               
        END DO                                                                  
      END DO                                                                    
c -----------------------------------------------------------------------       
999   RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE GetbOut(npP,nP,P,pstar,bOut,k)                                 
c =======================================================================       
c This subroutine inserts two impact parameters corresponding to pstar,         
c producing bOut(nP+2) from P(nP). The inserted elements are bOut(k) and        
c bOut(k+1)                                            [Z.I., Aug. 1996]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npP, nP, k, kstop, i                                              
      DOUBLE PRECISION P(npP), bOut(npP+2), pstar                               
c -----------------------------------------------------------------------       
      k = 0                                                                     
      kstop = 0                                                                 
      DO WHILE (kstop.NE.1)                                                     
        k = k + 1                                                               
        bOut(k) = P(k)                                                          
        IF (1.001*pstar.LE.P(k).OR.k.EQ.nP) kstop = 1                           
      END DO                                                                    
      IF (0.999*pstar.GT.P(k-1)) THEN                                           
        bOut(k) = 0.999*pstar                                                   
        ELSE                                                                    
        bOut(k) = 0.5*(P(k-1)+1.001*pstar)                                      
      END IF                                                                    
      IF (1.001*pstar.LT.P(k)) THEN                                             
        bOut(k+1) = 1.001*pstar                                                 
        ELSE                                                                    
        bOut(k+1) = 0.5*(P(k)+0.999*pstar)                                      
      END IF                                                                    
      DO i = k, nP                                                              
        bOut(i+2) = P(i)                                                        
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE getOmega(nG)                                                   
c =======================================================================       
c This subroutine generates albedo omega(iL,iY) from the abs/sca cross-         
c sections and the component abundancies. This is temporary (trivial)           
c version  for single size grains.                     [Z.I., Mar. 1996]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER  iG, nG, iL, iY                                                   
c -----------------------------------------------------------------------       
c     generate overall albedo through the envelope                              
c      ** this is for future multigrain code **                                 
c      ** for single grains it is trivial **                                    
      DO iY = 1, nY                                                             
c       calculate albedo                                                        
        DO iL = 1, nL                                                           
          omega(iL,iY) = SigmaS(nG,iL) / (SigmaA(nG,iL) + SigmaS(nG,iL))        
        END DO                                                                  
c       calculate relative abundances                                           
        DO iG = 1, nG                                                           
          abund(iG,iY) = 1.0                                                    
        END DO                                                                  
      END DO                                                                    
c ----------------------------------------------------------------------        
      RETURN                                                                    
      END                                                                       
c **********************************************************************        
                                                                                
c **********************************************************************        
      SUBROUTINE getOptPr(nG,nameQ,nameNK,er)                                   
c =====================================================================         
c This subroutine calculates the absorption and scattering efficiences          
c Qabs and Qsca in the wavelength range of the code or in case of               
c user supplied efficiences reads them from a file.                             
c                                                 [ZI Mar96; MN Aug97]          
c =====================================================================         
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER npLnk                                                             
      PARAMETER (npLnk=98)                                                      
      DOUBLE PRECISION n_sil_ow(npLnk),k_sil_ow(npLnk),n_sil_oc(npLnk),         
     &       k_sil_oc(npLnk), n_sil_dl(npLnk), k_sil_dl(npLnk),                 
     &       n_amc_hn(npLnk), k_amc_hn(npLnk), n_sic_pg(npLnk),                 
     &       k_sic_pg(npLnk), n_gr1_dl(npLnk), k_gr1_dl(npLnk),                 
     &       n_gr2_dl(npLnk), k_gr2_dl(npLnk), lam_nk(npLnk)                    
      COMMON /nkdat/ n_sil_ow, k_sil_ow, n_sil_oc, k_sil_oc,                    
     &               n_sil_dl, k_sil_dl, n_amc_hn, k_amc_hn,                    
     &               n_sic_pg, k_sic_pg, n_gr1_dl, k_gr1_dl,                    
     &               n_gr2_dl, k_gr2_dl, lam_nk                                 
      CHARACTER*235 nameQ(npG), nameNK(10), Fname, dummy*132                    
      INTEGER iG, nG, io1, iL, nLin, iLaux, Nprop, Na, iA, iC, iCuser,          
     &        er, Nmax, npA                                                     
c     Nmax is the number of records in user supplied file with opt.prop.        
c     and npA is the dimension of the array of grain sizes                      
      PARAMETER (Nmax=10000, npA=100)                                           
      DOUBLE PRECISION aa,bb,cc,lambdain(Nmax),Qain(Nmax),Qsin(Nmax),           
     &       n(npL),k(npL), aQabs(npA,npL),aQsca(npA,npL), amax, Pi,            
     &       nsd(npA), a(npA), faux1(npA), faux2(npA), f(npA), int,             
     &       ala(Nmax), SigAbs(npA,npL), SigSca(npA,npL), SizeDist,             
     &       aQa(Nmax), aQs(Nmax),  Cnorm, a3ave,                               
     &       n_int(npL), k_int(npL)                                             
c ----------------------------------------------------------------              
c     this should never change                                                  
      nL = npL                                                                  
      Nprop = 7                                                                 
c ----------------------------------------------------------------              
      Pi = 2.0*ASIN(1.0)                                                        
      er = 0                                                                    
c     first check that the user supplied wavelength grid is                     
c     monotonously increasing                                                   
      IF (top.LT.3) THEN                                                        
c       calculate efficiencies from n and k by Mie theory                       
c       generate the size array                                                 
        IF (szds.GT.2) THEN                                                     
          amax = 5.0*a2                                                         
        ELSE                                                                    
          amax = a2                                                             
        END IF                                                                  
        IF (dabs(a1-a2).LE.1.d-3) THEN                                          
          nA = 1                                                                
        ELSE                                                                    
          nA =50                                                                
        END IF                                                                  
c       Build-up the array of sizes a(nA)                                       
        CALL GETsizes(npA,nA,a1,amax,a)                                         
c       evaluate the normalization constant for the size                        
c       distribution nsd(nA)                                                    
        DO iA = 1, nA                                                           
          nsd(iA) = SizeDist(qsd,a(iA),szds,a2)                                 
        END DO                                                                  
        CALL PowerInt(npA,1,nA,a,nsd,Cnorm)                                     
c       find the average grain volume aveV (needed in dynamics)                 
        IF(dabs(a1-a2).LE.1.d-3) THEN                                           
          aveV = 4./3.*Pi*a1**3                                                 
        ELSE                                                                    
           DO iA = 1, nA                                                        
             faux1(iA)=nsd(iA)*a(iA)**3                                         
           END DO                                                               
           CALL PowerInt(npA,1,nA,a,faux1,a3ave)                                
           aveV = 4./3.*Pi*a3ave/Cnorm                                          
        END IF                                                                  
c      --  LOOP OVER SUPPORTED COMPONENTS --                                    
        DO iC= 1, Nprop                                                         
          f(iC) = xC(iC)                                                        
c         assign optical properties                                             
          IF (iC.EQ.1) CALL AssProp(npLnk,n_sil_ow,k_sil_ow,n,k)                
          IF (iC.EQ.2) CALL AssProp(npLnk,n_sil_oc,k_sil_oc,n,k)                
          IF (iC.EQ.3) CALL AssProp(npLnk,n_sil_dl,k_sil_dl,n,k)                
          IF (iC.EQ.4) CALL AssProp(npLnk,n_gr1_dl,k_gr1_dl,n,k)                
          IF (iC.EQ.5) CALL AssProp(npLnk,n_gr2_dl,k_gr2_dl,n,k)                
          IF (iC.EQ.6) CALL AssProp(npLnk,n_amc_hn,k_amc_hn,n,k)                
          IF (iC.EQ.7) CALL AssProp(npLnk,n_sic_pg,k_sic_pg,n,k)                
c         interpolate from opt. prop. grid to working grid                      
          DO iL = 1, nL                                                         
             CALL LinInter(npLnk,npLnk,lam_nk,n,lambda(iL),iLaux,aa)            
             n_int(iL) = aa                                                     
             CALL LinInter(npLnk,npLnk,lam_nk,k,lambda(iL),iLaux,aa)            
             k_int(iL) = aa                                                     
          END DO                                                                
c         calculate Qabs and Qsca (on working grid, i.e. lambda)                
          CALL MIE(npL,nL,lambda,n_int,k_int,npA,nA,a,1,aQabs,aQsca)            
c         for each lambda integrate Pi*a^2*Qext with n(a)da                     
          DO iL = 1, nL                                                         
            DO iA = 1, nA                                                       
              faux1(iA)=nsd(iA)*aQabs(iA,iL)*Pi*a(iA)**2                        
              faux2(iA)=nsd(iA)*aQsca(iA,iL)*Pi*a(iA)**2                        
            END DO                                                              
            CALL PowerInt(npA,1,nA,a,faux1,int)                                 
            sigAbs(iC,iL) = int/Cnorm                                           
            CALL PowerInt(npA,1,nA,a,faux2,int)                                 
            sigSca(iC,iL) = int/Cnorm                                           
          END DO                                                                
        END DO                                                                  
        IF (top.EQ.2) THEN                                                      
c         --  LOOP OVER USER SUPPLIED COMPONENTS --                             
          DO iCuser = 1, Nfiles                                                 
            iC = Nprop + iCuser                                                 
            f(iC) = xCuser(iCuser)                                              
c           read in optical properties                                          
            Fname = nameNK(iCuser)                                              
            CALL GetProp(npL,lambda,nL,Fname,n,k,er)                            
            IF (er.EQ.3) goto 999                                               
c           calculate Qabs and Qsca                                             
            CALL MIE(npL,nL,lambda,n,k,npA,nA,a,1,aQabs,aQsca)                  
c           for each lambda integrate Pi*a^2*Qext with n(a)da                   
            DO iL = 1, nL                                                       
              DO iA = 1, nA                                                     
                faux1(iA)=nsd(iA)*aQabs(iA,iL)*Pi*a(iA)**2                      
                faux2(iA)=nsd(iA)*aQsca(iA,iL)*Pi*a(iA)**2                      
              END DO                                                            
              CALL PowerInt(npA,1,nA,a,faux1,int)                               
              sigAbs(iC,iL) = int/Cnorm                                         
              CALL PowerInt(npA,1,nA,a,faux2,int)                               
              sigSca(iC,iL) = int/Cnorm                                         
            END DO                                                              
          END DO                                                                
        ELSE                                                                    
          Nfiles = 0                                                            
        END IF                                                                  
c       mix them together (syntetic grain model)                                
        DO iL = 1, nL                                                           
          SigmaA(1,iL) = 0.0                                                    
          SigmaS(1,iL) = 0.0                                                    
          DO iC= 1, Nprop+Nfiles                                                
            SigmaA(1,iL) = SigmaA(1,iL) + f(iC) * sigAbs(iC,iL)                 
            SigmaS(1,iL) = SigmaS(1,iL) + f(iC) * sigSca(iC,iL)                 
          END DO                                                                
        END DO                                                                  
      ELSE                                                                      
c     this is for top.GE.3                                                      
c      initialize aveV for this case
       aveV = 1.
c       read in lambda grid and optical properties                              
        DO iG = 1, nG                                                           
          open(1,ERR=998,file=nameQ(iG),STATUS='OLD')                           
          read(1,'(a)',ERR=998)dummy                                            
          read(1,'(a)',ERR=998)dummy                                            
          read(1,'(a)',ERR=998)dummy                                            
          iL = 0                                                                
          io1 = 0                                                               
          DO WHILE (io1.GE.0)                                                   
            read(1,*,END=900,ERR=998,iostat=io1) aa, bb, cc                     
            IF (io1.GE.0) THEN                                                  
              iL = iL + 1                                                       
              lambdain(iL) = aa                                                 
              Qain(iL) = bb                                                     
              Qsin(iL) = cc                                                     
            END IF                                                              
          END DO                                                                
900       close(1)                                                              
          IF (iL.LT.2) goto 998                                                 
          nLin = iL                                                             
c         if input wavelengths in descending order turn them around             
          IF (lambdain(1).GT.lambdain(2)) THEN                                  
            DO iL = 1, nLin                                                     
              ala(iL) = lambdain(iL)                                            
              aQa(iL) = Qain(iL)                                                
              aQs(iL) = Qsin(iL)                                                
            END DO                                                              
            DO iL = 1, nLin                                                     
              lambdain(iL) = ala(nLin+1-iL)                                     
              Qain(iL) = aQa(nLin+1-iL)                                         
              Qsin(iL) = aQs(nLin+1-iL)                                         
            END DO                                                              
          END IF                                                                
c         interpolate to Dusty's wavelength grid                                
          DO iL = 1, nL                                                         
            CALL LinInter(Nmax,nLin,lambdain,Qain,lambda(iL),iLaux,aa)          
            SigmaA(iG,iL) = aa                                                  
            CALL LinInter(Nmax,nLin,lambdain,Qsin,lambda(iL),iLaux,aa)          
            SigmaS(iG,iL) = aa                                                  
          END DO                                                                
        END DO                                                                  
      END IF                                                                    
      goto 999                                                                  
998   write(12,*)' ***  FATAL ERROR IN DUSTY  ***********'                      
      write(12,*)' File with optical properties:'                               
      write(12,'(a2,a70)')'  ',nameQ(iG)                                        
      write(12,*)' is missing or not properly formatted?!'                      
      write(12,*)' **************************************'                      
      close(12)                                                                 
      er = 3                                                                    
c -----------------------------------------------------------------------       
999   RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE GetProp(npL,lambda,nL,Fname,en,ek,error)                       
c =======================================================================       
c This subroutine reads optical properties en(i,j), ek(i,j) from file           
c fname(Nf), with i=Nf, j=1..NLL(Nf), and interpolates them onto                
c wavelength grid lambda(1..nL)                        [Z.I., Mar. 1996]        
c =======================================================================       
      IMPLICIT none                                                             
      CHARACTER*235 Fname                                                       
      CHARACTER*232 line                                                        
      INTEGER i, nL, iloc, iL, npL, io1, error, Nmax                            
c     Nmax is the number of records in the user supplied file                   
      PARAMETER (Nmax=10000)	                                                   
      DOUBLE PRECISION en(npL), ek(npL), lambda(npL), pw(Nmax),                 
     &       pren(Nmax), pimn(Nmax), a(Nmax), b(Nmax), c(Nmax), aa,             
     &       bb, cc                                                             
c -----------------------------------------------------------------------       
      error = 0                                                                 
      open(2,ERR=998,file=Fname,STATUS='OLD')                                   
c     read in a header from the input file                                      
      DO i = 1, 7                                                               
        read(2,'(a)',ERR=998)line                                               
      END DO                                                                    
c     read in input data                                                        
      iL = 0                                                                    
      io1 = 0                                                                   
      DO WHILE (io1.GE.0)                                                       
        read(2,*,END=900,ERR=998,iostat=io1) aa, bb, cc                         
        IF (io1.GE.0) THEN                                                      
          iL = iL + 1                                                           
          pw(iL) = aa                                                           
          pren(iL) = bb                                                         
          pimn(iL) = cc                                                         
        END IF                                                                  
      END DO                                                                    
900   close(2)                                                                  
      IF (iL.LT.2) goto 998                                                     
c     if input wavelengths in descending order turn them around                 
      IF (pw(1).GT.pw(2)) THEN                                                  
        DO i = 1, iL                                                            
          a(i) = pw(i)                                                          
          b(i) = pren(i)                                                        
          c(i) = pimn(i)                                                        
        END DO                                                                  
        DO i = 1, iL                                                            
          pw(i) = a(iL+1-i)                                                     
          pren(i) = b(iL+1-i)                                                   
          pimn(i) = c(iL+1-i)                                                   
        END DO                                                                  
      END IF                                                                    
c     interpolate                                                               
      DO i = 1, nL                                                              
        CALL LinInter(Nmax,iL,pw,pren,lambda(i),iloc,en(i))                     
        CALL LinInter(Nmax,iL,pw,pimn,lambda(i),iloc,ek(i))                     
      END DO                                                                    
      goto 999                                                                  
998   write(12,*)' ***  FATAL ERROR IN DUSTY  ***********'                      
      write(12,*)' File with optical properties:'                               
      write(12,'(a2,a70)')'  ',Fname                                            
      write(12,*)' is missing or not properly formatted?!'                      
      write(12,*)' **************************************'                      
      close(12)                                                                 
      error = 3                                                                 
c -----------------------------------------------------------------------       
999   RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE GETsizes(NN,N,x1,x2,x)                                         
c =======================================================================       
c This subroutine generates an array x(i=1..N) of physical size NN,             
c with N elements logarithmically spaced between x1 and x2.                     
c                                              [ZI,Aug'96;MN,Nov'97]            
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER NN, N, i                                                          
      DOUBLE PRECISION x(NN), x1, x2, fac, pw1, pw2                             
c -----------------------------------------------------------------------       
      IF (N.GT.1) THEN                                                          
        pw1 = 1.0/(N-1)                                                         
        fac = (x2/x1)**pw1                                                      
        DO i = 1, N                                                             
          pw2 = 1.0*(i-1)                                                       
          x(i) = x1*fac **pw2                                                   
        END DO                                                                  
      ELSE                                                                      
        x(1) = x1                                                               
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE getTau(model,nG,TAU1,TAU2,TAUin,Nrec,GridType,Nmodel)          
c =======================================================================       
c This subroutine generates total optical depth TAUtot.                         
c                                                      [Z.I., Mar. 1996]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      INTEGER model, Nmodel, nG, iL, GridType, Nrec                             
      DOUBLE PRECISION faux1(npL), faux2(npL), SigAfid, SigSfid, q,             
     &       TAUin(Nrec), TAU1, TAU2                                            
c -----------------------------------------------------------------------       
      IF (nG.GT.1) THEN                                                         
        write(18,*)'Fix getTAU, nG>1 !'                                         
        stop                                                                    
      END IF                                                                    
      TAUmax = 0.0                                                              
      IF(GridType.EQ.3) THEN                                                    
        TAUfid = TAUin(model)                                                   
      ELSE                                                                      
c      calculate TAUfid for given model                                         
       IF (model.EQ.1) THEN                                                     
        TAUfid = TAU1                                                           
       ELSE                                                                     
        IF (model.EQ.Nmodel) THEN                                               
          TAUfid = TAU2                                                         
        ELSE                                                                    
          IF (GridType.EQ.1) THEN                                               
            q =  (TAU2 - TAU1)/(Nmodel-1.0)                                     
            TAUfid = TAU1 + q*(model-1.0)                                       
          ELSE                                                                  
            q = dexp(log(TAU2/TAU1)/(Nmodel-1.0))                               
            TAUfid = TAU1 * q**(model-1.0)                                      
          END IF                                                                
        END IF                                                                  
       END IF                                                                   
      END IF                                                                    
c     generate TAUtot and find TAUmax                                           
      DO iL = 1, nL                                                             
        faux1(iL) = SigmaA(nG,iL)                                               
        faux2(iL) = SigmaS(nG,iL)                                               
      END DO                                                                    
      IF (lamfid.LT.lambda(1)) THEN                                             
        write(12,*)' Fiducial wavelength was too small.'                        
        write(12,'(a8,e9.3,a17)')' Using ',lambda(1),' micron instead.'         
      END IF                                                                    
      IF (lamfid.GT.lambda(nL)) THEN                                            
        write(12,*)' Fiducial wavelength was too large.'                        
        write(12,'(a8,e9.3,a17)')' Using ',lambda(nL),' micron instead.'        
      END IF                                                                    
      CALL LinInter(npL,nL,lambda,faux1,lamfid,iLfid,SigAfid)                   
      CALL LinInter(npL,nL,lambda,faux2,lamfid,iLfid,SigSfid)                   
c     extinction efficiency at fiducial wavelength                              
      SigExfid = SigAfid + SigSfid                                              
      DO iL = 1, nL                                                             
        TAUtot(iL) = TAUfid*(SigmaA(nG,iL) + SigmaS(nG,iL)) / SigExfid          
        IF (TAUtot(iL).GE.TAUmax) TAUmax = TAUtot(iL)                           
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE getETAzp(ETAzp)                                                
c =======================================================================       
c This function calculates ETAzp(iP,iZ) along the line of sight with            
c impact parameter P(iP) and iZ=1, nZ. Here iZ = 1 corresponds to z=0           
c and iZ=nZ to the outer edge. Other grid points coincide with the              
c radial grid. The method used is spline approximation for normalized           
c density distribution ETA, with subsequent z-integration performed             
c analytically in function IntETA                                               
c                                               [ZI,Feb'95; MN,Aug'97]          
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iYfirst, YPequal, Plast                                           
      DIMENSION iYfirst(npP), YPequal(npP), Plast(npY)                          
      COMMON /Yfirst/ iYfirst, YPequal, Plast                                   
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER iP, nZ, iZ, iW                                                    
      DOUBLE PRECISION ETAzp(npP,npY), IntETA, auxEta, w1, w2                   
c -----------------------------------------------------------------------       
c     loop over impact parameters                                               
      DO iP = 1, nP                                                             
c       maximal number of points along tangential position, z                   
        nZ = nY + 1 - iYfirst(iP)                                               
c       starting values for z and ETAzp(iP,iZ)                                  
        IF (P(iP).GE.1.0) THEN                                                  
          w2 = P(iP)                                                            
          ELSE                                                                  
          w2 = 1.0                                                              
        END IF                                                                  
c       initialize ETAzp(iP,iZ)*TAUtot(iL)                                      
        ETAzp(iP,1) = 0.0                                                       
c       loop over z                                                             
        DO iZ = 2, nZ                                                           
c         index for local radius, w2                                            
          iW = iYfirst(iP) + iZ - 1                                             
c         limits for integration                                                
          w1 = w2                                                               
          w2 = Y(iW)                                                            
c           find next step in ETAzp                                             
            auxEta = IntETA(P(iP),iW-1,w1,w2)                                   
c           add next step in ETAzp                                              
            ETAzp(iP,iZ) = ETAzp(iP,iZ-1) + auxEta                              
        END DO                                                                  
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      DOUBLE PRECISION FUNCTION IMAGFN(y)                                       
c =======================================================================       
c This function evaluates auxiliary functions needed to produce                 
c visibility curves and convolved images. It is called from the image           
c integration subroutine ROMBY.                        [Z.I., Jan. 1997]        
c =======================================================================       
      IMPLICIT none                                                             
      DOUBLE PRECISION x, y, PSFN, Bessel                                       
      INTEGER ftype                                                             
      DOUBLE PRECISION Ckn, Cxout, Cphi, Cqtheta                                
      COMMON /imfn1/ ftype                                                      
      COMMON /imfn2/ Ckn, Cxout, Cphi, Cqtheta                                  
c -----------------------------------------------------------------------       
      IF (ftype.EQ.1) THEN                                                      
c       this part is for convolution                                            
        x = sqrt(abs(Cxout*Cxout+y*y-2.*Cxout*y*dcos(Cphi)))                    
        imagfn = PSFN(x) * y**Ckn                                               
      ELSE                                                                      
c       this part is for visibility                                             
c       argument is Pi*q*y (not 2*Pi*q*y) to account for the fact that          
c       theta1 is diameter rather than radius (so V is function of              
c       q*theta1, like in IE, '96, MNRAS 279, 1019)                             
        imagfn = Bessel(2.*ASIN(1.0)*Cqtheta*y) * y**Ckn                        
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE MIE(npL,nL,lambda,Ere,Eim,npA,nA,a,nG1,Qabs,Qsca)              
c =======================================================================       
c This subroutine calculates Qabs and Qsca for a given diffractive              
c index Ere, Eim, wavelength lambda and size a. Here, lambda is an              
c array (1..nL), Ere and Eim are given on this array, a is an array             
c of sizes (1..nA). Qabs and Qsca are arrays (nG1..nG1+nA,nL), i.e. for         
c each wavelength lambda, Qabs and Qsca are evaluated for nA different          
c sizes. The numbering, however, does not start from 1, but rather from         
c nG1.                                                [Z.I., Aug. 1996]         
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npL, nL, npA, nA, nG1, iL, iA                                     
      DOUBLE PRECISION lambda(npL), Ere(npL), Eim(npL), a(npA),                 
     &       Qabs(npA,npL), Qsca(npA,npL)                                       
      REAL xx, Qex, Qsc, Qback                                                  
      COMPLEX refrel, s1(200), s2(200)                                          
c -----------------------------------------------------------------------       
c     loop over wavelengths                                                     
      DO iL = 1, nL                                                             
c       complex index of refraction                                             
        refrel = cmplx(Ere(iL),Eim(iL))                                         
c       loop over sizes                                                         
        DO iA = 1, nA                                                           
c         size parameter                                                        
          xx=2.0*3.14159265*a(iA)/lambda(iL)                                    
c         if size parameter xx>100 use xx=100 (geometrical optics)              
          IF (xx.GT.100.0) xx = 100.0                                           
c         calculate efficiencies                                                
          CALL bhmie(xx,refrel,2,s1,s2,Qex,Qsc,Qback)                           
c         store the result                                                      
          Qabs(nG1+iA-1,iL) = Qex - Qsc                                         
          Qsca(nG1+iA-1,iL) = Qsc                                               
        END DO                                                                  
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      DOUBLE PRECISION FUNCTION PSFN(x)                                         
c =======================================================================       
c This function evaluates the point spread function. For psftype.EQ.1           
c the function is evaluated as a sum of two Gaussians, for psftype.EQ.2         
c it is provided by user in a file. psftype and all other relevant              
c parameters come from COMMON /psf/ and are initialized in subroutine           
c INPUT.                                               [Z.I., Jan. 1997]        
c =======================================================================       
      IMPLICIT none                                                             
      DOUBLE PRECISION x                                                        
      INTEGER idummy                                                            
      INTEGER psftype, Npsf, iLambda                                            
      DOUBLE PRECISION kPSF(20), FWHM1(20), FWHM2(20), theta1,                  
     &       xpsf(1000), ypsf(1000)                                             
      COMMON /psf1/ iLambda, psftype, Npsf                                      
      COMMON /psf2/ kPSF, FWHM1, FWHM2, Theta1,                                 
     &       xpsf, ypsf                                                         
c -----------------------------------------------------------------------       
      IF (psftype.LT.3) THEN                                                    
        psfn = dexp(-(1.665*x/FWHM1(iLambda))**2.)                              
        IF (psftype.EQ.2)                                                       
     &    psfn = (psfn + kPSF(iLambda) *                                        
     &           dexp(-(1.665*x/FWHM2(iLambda))**2.))/(1.+kPSF(iLambda))        
        ELSE                                                                    
        CALL LinInter(1000,Npsf,xpsf,ypsf,x,idummy,psfn)                        
      ENDIF                                                                     
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE setupETA                                                       
c =======================================================================       
c This subroutine finds spline coefficients ETAcoef (defined in file            
c 'density.inc') such that normalized density function ETA(Y(iY)) is:           
c ETAcoef(iY,1)+ETAcoef(iY,2)/Y(iY)+...+ETAcoef(iY,2)/Y(iY)^3                   
c If spline approximation differs more than maxerr (see below) at the           
c midpoint, then a straight line is used instead. (In case of wavelength        
c depend. ETA, use ETAfun where any new dens. laws should be described).        
c Coefficients ETAcoef are later used in getETAzp to calculate ETAzp.           
c                                                [ZI, Feb'96; MN,Aug'97]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iY, iC                                                            
      DOUBLE PRECISION coef(npY,4), ETA, maxerr, Ymid, Yinverse(npY),           
     &       ETAaux(npY), ETAmid(npY)                                           
c -----------------------------------------------------------------------       
c       generate input function for SPLINE2                                     
        DO iY = 1, nY                                                           
          Yinverse(iY) = 1. / Y(iY)                                             
          ETAaux(iY) = ETA(Y(iY))                                               
          IF (iY.LT.nY) THEN                                                    
            Ymid = dsqrt(Y(iY)*Y(iY+1))                                         
            ETAmid(iY) = ETA(Ymid)                                              
          END IF                                                                
        END DO                                                                  
c       calculate spline coefficients                                           
        CALL SPLINE2(Yinverse,ETAaux,nY,coef)                                   
c       check and fix spline coefficients                                       
        maxerr = 0.1                                                            
c       RDW is initialized in Input   	                                         
        CALL CHKSPLIN(Yinverse,ETAaux,ETAmid,nY,coef,maxerr,RDW)                
c       copy coefficients to the output array ETAcoef                           
        DO iY = 1, nY                                                           
          DO iC = 1, 4                                                          
            ETAcoef(iY,iC) = coef(iY,iC)                                        
          END DO                                                                
        END DO                                                                  
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      DOUBLE PRECISION FUNCTION SizeDist(q,aa,sdtype,a0)                        
c =======================================================================       
c This subroutine calculates size distribution n(a) for a=aa. The size          
c distribution is MRN type n(a)~1/a**q for sdtype.LE.2 and KMH type             
c n(a)~dexp(-a/a0)/a**q otherwise                                               
c                                                      [Z.I., Aug. 1996]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER sdtype                                                            
      DOUBLE PRECISION aa, a0, q                                                
c -----------------------------------------------------------------------       
      IF (sdtype.LE.2) THEN                                                     
        SizeDist = 1.0/aa**q                                                    
        ELSE                                                                    
        SizeDist = dexp(-aa/a0)/aa**q                                           
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE Spectral(model,denstyp,nL,Lambda)                              
c =======================================================================       
c     This subroutine finds the spectral features for spp and zpp files.        
c     It employs SpFeatur().		       		     [MN, Jan'99]                        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nYok, nPok                                                        
      DOUBLE PRECISION                                                          
     &       Ude(npL,npY), Uds(npL,npY), Us(npL,npY), fde(npL,npY),             
     &       fds(npL,npY), fs(npL,npY), Utot(npL,npY), ftot(npL,npY),           
     &       Td(npG,npY), Ubol(npY), fbol(npY), Uchck(npL,npY),                 
     &       UbolChck(npY), Spectrum(npL), SpecChar(30,99), tauF(npY),          
     &       Intens(npL,npP+2), IntOut(20,npP+2), SmC(30,99),                   
     &       Yok(npY), Pok(npP), tauOut(npL), tauZout(npP+2), F1,               
     &       tr(npY), rg(npG,npY), fsL(npL,npY), fsR(npL,npY), Eps(npY)         
      COMMON /solution/ Ude, Uds, Us, fde, fds, fs, Utot, ftot, Td,             
     &       Ubol, fbol, Uchck, UbolChck, Spectrum, SpecChar, tauF,             
     &       Intens, IntOut, SmC, Yok, Pok, tauOut, tauZout, F1,                
     &       tr, rg, fsL, fsR, Eps, nYok, nPok                                  
      INTEGER startyp(2), Nlamtr(2), nBB(2)                                     
      CHARACTER nameStar(2)*235                                                 
      DOUBLE PRECISION Tstar, lamtr(2,101), klam(2,100), Tbb(2,10),             
     &                 rellum(2,10), mu1, ksi, mu2, xSiO, r1rs                  
      COMMON /source/ Tstar, lamtr, klam, Tbb, rellum, mu1, ksi, mu2,           
     &                xSiO, r1rs, startyp, Nlamtr, nBB, nameStar                
      INTEGER i, iL, model, nL, denstyp, Nchar                                  
      PARAMETER (Nchar=11)                                                      
      DOUBLE PRECISION Charac(Nchar),Spectr(npL),Lambda(npL)                    
c -----------------------------------------------------------------------       
c     find features for *.spp file                                              
      IF (denstyp.EQ.0) THEN                                                    
         DO iL = 1, nL                                                          
            Spectr(iL) = ftot(iL,nYok) + ksi*fsR(iL,nYok)                       
       END DO                                                                   
      ELSE                                                                      
         DO iL = 1, nL                                                          
            Spectr(iL) = Spectrum(iL)                                           
       END DO                                                                   
      END IF                                                                    
      CALL SpFeatur(model,nL,Lambda,Spectr,Charac)                              
      DO i = 1, Nchar                                                           
       SpecChar(i,model) = Charac(i)                                            
      END DO                                                                    
c     find features for *.zpp file                                              
      IF (denstyp.EQ.0) THEN                                                    
         DO iL = 1, nL                                                          
            Spectr(iL) = dabs(ftot(iL,1) - fsL(iL,1))                           
         END DO                                                                 
       CALL SpFeatur(model,nL,Lambda,Spectr,Charac)                             
         DO i = 1, Nchar                                                        
            SpecChar(i+11,model) = Charac(i)                                    
         END DO                                                                 
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE SpFeatur(model,nL,Lambda,Spectr,Charac)                        
c =======================================================================       
c This subroutine calculates IRAS colors and other spectral quantities          
c Filters data from Neugebauer et al, 1984, ApJ, 278, L1.                       
c Procedure described in Bedijn, 1987, A&A, 186, 136.  [Z.I., Mar. 1996]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER i, j, iaux, model, nL                                             
      DOUBLE PRECISION  f1(7), f2(7), f3(7), f4(7), w1(7), w2(7), w3(7),        
     &        w4(7), phi(4,7), al(4,7), cl(4), prz(4), tinf(4), flxy(9),        
     &        wav(9), lambda(npL), Spectr(npL), Charac(11), wmid, flx1,         
     &        flx2, faux, B98, B11, rat9818, beta813, beta1422, f12,            
     &        f25, f60, f100, an98, f98c, f11c, an11                            
c -----------------------------------------------------------------------       
c     Data for 4 IRAS filters                                                   
c     wavelengths                                                               
       DATA w1/7.55, 8.0, 10.3, 11.5, 13.4, 14.7, 15.5/                         
       DATA w2/16.6, 22.5, 25.6, 26.8, 27.5, 29.3, 31.1/                        
       DATA w3/30.5, 40.1, 40.2, 65.2, 74.2, 83.8, 83.9/                        
       DATA w4/72.7, 95.4, 111.2, 116.6, 137.4, 137.5, 137.6/                   
c     transmittivities                                                          
       DATA f1/0.000, 0.618, 0.940, 0.750, 1.022, 0.906, 0.000/                 
       DATA f2/0.235, 0.939, 0.939, 0.745, 0.847, 0.847, 0.000/                 
       DATA f3/0.000, 0.102, 0.260, 1.026, 0.842, 0.001, 0.000/                 
       DATA f4/0.000, 0.910, 1.000, 0.330, 0.002, 0.001, 0.000/                 
c     ------------------------------------------------------------              
c     initialization                                                            
       DO i = 1, 4                                                              
         tinf(i) = 0.0                                                          
         cl(i) = 0.0                                                            
       END DO                                                                   
       DO j = 1, 7                                                              
         al(1,j) = w1(j)                                                        
         phi(1,j) = f1(j)                                                       
         al(2,j) = w2(j)                                                        
         phi(2,j) = f2(j)                                                       
         al(3,j) = w3(j)                                                        
         phi(3,j) = f3(j)                                                       
         al(4,j) = w4(j)                                                        
         phi(4,j) = f4(j)                                                       
       END DO                                                                   
c ------------------------------------------------------------------------      
c     first find IRAS colors                                                    
      DO j = 2, nL                                                              
c       middle wavelength                                                       
        wmid = 0.5*(Lambda(j-1)+Lambda(j))                                      
c       interpolate filters for wmid                                            
        CALL PHILAM(wmid,prz,al,phi)                                            
c       convert Spectrum to Flambda                                             
        flx1 = Spectr(j-1) / Lambda(j-1)                                        
        flx2 = Spectr(j) / Lambda(j)                                            
c       add contribution to the integral (index is over filters)                
        DO i = 1, 4                                                             
          tinf(i) = tinf(i) + prz(i)*0.5*(flx2+flx1)*                           
     &                                    (Lambda(j)-Lambda(j-1))               
          cl(i) = cl(i) + prz(i) * (Lambda(j)-Lambda(j-1))                      
        END DO                                                                  
      END DO                                                                    
      DO i = 1, 4                                                               
        tinf(i) = tinf(i) / cl(i)                                               
      END DO                                                                    
c     Spectrum corrected for IRAS filters                                       
      f12 = tinf(1)*12.0                                                        
      f25 = tinf(2)*25.0                                                        
      f60 = tinf(3)*60.0                                                        
      f100 = tinf(4)*100.0                                                      
c     now find  B98, B11, F98/F18, beta 8-13, beta 14-22                        
c     find fluxes at all needed wavelengths (energy increases with index)       
      DATA wav/2.2, 8.0, 9.8, 11.3, 13.0, 14.0, 18.0, 22.0, 0.55/               
      DO j = 1, 9                                                               
        CALL LinInter(npL,nL,lambda,Spectr,wav(j),iaux,faux)                    
        flxy(j) = faux                                                          
      END DO                                                                    
c     the feature strength at 9.8 and 11.4 microns                              
      an98 = log(flxy(5)/flxy(2))/log(wav(5)/wav(2))                            
      f98c = flxy(2)*(wav(3)/wav(2))**an98                                      
      B98 = log(flxy(3)/f98c)                                                   
      an11 = log(flxy(5)/flxy(3))/log(wav(5)/wav(3))                            
      f11c = flxy(3)*(wav(4)/wav(3))**an11                                      
      B11 = log(flxy(4)/f11c)                                                   
c     ratio F9.8/F18                                                            
      rat9818 = flxy(3)/flxy(7)*wav(7)/wav(3)                                   
c     beta 8-13 and beta 14-22 (see Neugebauer)                                 
      beta813 = log(flxy(5)/flxy(2))/log(13.0/8.0) - 1.0                        
      beta1422 = log(flxy(8)/flxy(6))/log(22.0/14.0) - 1.0                      
c     store SpecChar to output array SpecChar                                   
      Charac(1) = B98                                                           
      Charac(2) = B11                                                           
      Charac(3) = rat9818                                                       
      Charac(4) = beta813                                                       
      Charac(5) = beta1422                                                      
      Charac(6) = flxy(9)                                                       
      Charac(7) = flxy(1)                                                       
      Charac(8) = f12                                                           
      Charac(9) = log10(25.0*f25/f12/12.0)                                      
      Charac(10) = log10(60.0*f60/f12/12.0)                                     
      Charac(11) = log10(f100*100.0/f60/60.0)                                   
c -------------------------------------------------------------------------     
      RETURN                                                                    
      END                                                                       
c *******************************************************************           
                                                                                
c *******************************************************************           
      SUBROUTINE PHILAM(Alam,F,Al,Phi)                                          
c     interpolates IRAS filters [f(4)] for a given wavelength alam              
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                        
      DIMENSION Phi(4,7), Al(4,7), F(4)                                         
c -----------------------------------------------                               
      DO i = 1, 4                                                               
        F(i) = 0.0                                                              
        im = 0                                                                  
        istop = 0                                                               
        DO WHILE (istop.NE.1)                                                   
              im = im + 1                                                       
            IF ( (Alam-Al(i,im))*(Alam-Al(i,im+1)).LE.0) THEN                   
               a = (Phi(i,im+1)-Phi(i,im))                                      
     &           /log10(Al(i,im+1)/Al(i,im))                                    
               b = Phi(i,im) - a*log10(Al(i,im))                                
               F(i) = a*log10(Alam) + b                                         
            IF ( F(i).GT.1.0 ) F(i) = 1.0                                       
              END IF                                                            
            IF (im.EQ.6) istop = 1                                              
        END DO                                                                  
      END DO                                                                    
c ------------------------------------------------------                        
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE Visibili(IntOut)                                               
c =======================================================================       
c This subroutine finds visibility functions corresponding to IntOut.           
c The work horse is subroutine Visi2D, and this subroutine is used to           
c prepare everything.                                  [Z.I., Jan. 1997]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER i, j, N1, N2                                                      
      DOUBLE PRECISION  IntOut(20,npP+2), Visi(1000), Int1D(npP+2)              
c -----------------------------------------------------------------------       
c     generate spatial frequency (q) grid                                       
c     first N1 points up to qtheta1=1.22 (Rayleigh limit for a disk)            
      N1 = 80                                                                   
c     first 2 points manually:                                                  
c     there must be 0!                                                          
      qtheta1(1) = 0.0                                                          
c     make sure the whole envelope is resolved                                  
      qtheta1(2) = 0.5 / bOut(nP+2)                                             
c     and the rest on logarithmic grid up to 1.22                               
      DO i = 1, N1-2                                                            
        qtheta1(i+2) = qtheta1(2) * (1.22/qtheta1(2))**(i*1.0/(N1-2))           
      END DO                                                                    
c     envelope is well sampled, now to be sure that the star will be OK         
c     for small taus add N2 points on a logarithmic grid up to 1.22/p*          
      N2 = 20                                                                   
      DO i = 1, N2                                                              
        qtheta1(N1+i) = 1.22 / bOut(2)**(i*1.0/N2)                              
      END DO                                                                    
      Nvisi = N1 + N2                                                           
c     find visibility wavelength by wavelength                                  
      DO j = 1, NlambdaOut                                                      
        DO i = 1, nP+2                                                          
          Int1D(i) = IntOut(j,i)                                                
          CALL CHKRANGE(dynrange,Int1D(i))                                      
          IF (Int1D(i).LT.dynrange) Int1D(i)=0.0                                
        END DO                                                                  
        CALL Visi2D(npP+2,nP+2,bOut,Int1D,1000,N1+N2,qtheta1,Visi)              
c       copy 1D convolved visibility to Visib                                   
        DO i = 1, N1+N2                                                         
c         check dynamic range                                                   
          CALL CHKRANGE(dynrange,Visi(i))                                       
          Visib(j,i) = Visi(i)                                                  
        END DO                                                                  
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE Visi2D(NinMax,Nin,Xin,Yin,Noutmax,Nout,Xout,Yout)              
c =======================================================================       
c This subroutine finds the visibility function (the spatial Fourier            
c transform of the intensity distribution) corresponding to the                 
c intensity Yin(Xin[i]), i=1,Nin. Visibility, Yout, is evaluated at q           
c positions (spatial frequency) given in Xout[i], i=1,Nout. Maximum size        
c of Xin is NinMax, maximum size of Xout is NoutMax. The Bessel function        
c of the zeroth order is provided separately. The integration is done by        
c calling subroutine ROMBY (Bessel function is called from IMGFN).              
c Note:                                                                         
c The visibility function V(q) for a circularly symmetric intensity             
c I(x) is:                                                                      
c          V(q) = F(q)/F(0)                                                     
c where Jo is the Bessel function of the zeroth order, and                      
c          F(q) = Int[Jo(2Pi*q*x)*I(x)*2Pi*x*dx]                                
c Note that F(0) is nothing more than flux. For more details see                
c Ivezic & Elitzur, 1996, MNRAS, 279, 1019 and ref. therein.                    
c                                                      [Z.I., Jan. 1997]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER NinMax, Nin, NoutMax, Nout, iq, iXin                              
      DOUBLE PRECISION Xin(NinMax), Yin(NinMax), Xout(NoutMax),                 
     &       Yout(NoutMax),  F(1000), F0, int1, int2, A, B, imagfn              
      INTEGER ftype                                                             
      DOUBLE PRECISION Ckn, Cxout, Cphi, Cqtheta                                
      COMMON /imfn1/ ftype                                                      
      COMMON /imfn2/ Ckn, Cxout, Cphi, Cqtheta                                  
      EXTERNAL imagfn                                                           
c -----------------------------------------------------------------------       
c     loop over spatial frequency q (= Xout)                                    
      DO iq = 1, Nout                                                           
        Cqtheta = Xout(iq)                                                      
        F(iq) = 0.0                                                             
c       loop over radial positions                                              
        DO iXin = 1, Nin                                                        
c         find F(q)                                                             
          ftype = 2                                                             
          Ckn = 1.0                                                             
          CALL ROMBY(imagfn,Xin(iXin),Xin(iXin+1),int1)                         
          Ckn = 2.0                                                             
          CALL ROMBY(imagfn,Xin(iXin),Xin(iXin+1),int2)                         
c         contribution from this annulus (lin. approx. for intensity)           
          A = Xin(iXin+1)*Yin(iXin)-Xin(iXin)*Yin(iXin+1)                       
          A = A /(Xin(iXin+1)-Xin(iXin))                                        
          B = (Yin(iXin+1)-Yin(iXin))/(Xin(iXin+1)-Xin(iXin))                   
          F(iq) = F(iq) + A*int1 + B*int2                                       
        END DO                                                                  
      END DO                                                                    
c     flux                                                                      
      F0 = F(1)                                                                 
      DO iq = 1, Nout                                                           
       IF(F0.EQ.0.) THEN                                                        
         Yout(iq) = 0.                                                          
       ELSE                                                                     
         Yout(iq) = dabs(F(iq) / F0)                                            
       END IF                                                                   
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
c ======================================================================        
c This is the file with auxiliary math and number-to-string conversion          
c subroutines.                                             [MN, Mar,99]         
c ======================================================================        
C     Table of Contents                                                         
C                                                                               
C     ADD                                                                       
C     ADD2                                                                      
C     ANALINT                                                                   
C     Attach                                                                    
C     BESSEL                                                                    
C     BHMIE                                                                     
C     CHKRANGE                                                                  
C     CHKSPLIN                                                                  
C     Clean                                                                     
C     EMPTY                                                                     
C     FileMSG                                                                   
C     FindErr                                                                   
C     FINDMAX                                                                   
C     FINDRMS                                                                   
C     GETFS                                                                     
C     H                                                                         
C     KRON                                                                      
C     LINE                                                                      
C     LININTER                                                                  
C     LINSYS                                                                    
C     LUBKSB                                                                    
C     LUDCMP                                                                    
C     MAKETABLE                                                                 
C     MAPLE3                                                                    
C     MIDSQL                                                                    
C     MPROVE                                                                    
C     MSG                                                                       
C     MULTIPLY                                                                  
C     MULTIP2                                                                   
C     MYSPLINE                                                                  
C     POLINT                                                                    
C     POWERINT                                                                  
C     PRODUCT                                                                   
C     ROMBERG2                                                                  
C     ROMBY                                                                     
C     SCALETO1                                                                  
C     SHIFT                                                                     
C     SIMPSON                                                                   
C     SORT                                                                      
C     SPLINE                                                                    
C     SPLINE2                                                                   
C     SPLINT                                                                    
C     TRAPZD                                                                    
C     WriteOut                                                                  
C     ZBRAC                                                                     
C     ZRIDDR                                                                    
c ======================================================================        
                                                                                
                                                                                
c **********************************************************************        
      SUBROUTINE ADD(np1,nr1,np2,nr2,q1,q2,q3,qOut)                             
c ======================================================================        
c This subroutine evaluates the following expression:                           
c [qOut] = [q1] + [q2] + [q3]. qout, q1, q2 and q2 are matrices of              
c physical size (np2,np1) and real size (nr2,nr1).     [Z.I., Nov. 1995]        
c ======================================================================        
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER  np1, nr1, np2, nr2, i2, i1                                       
      DOUBLE PRECISION  q1(np2,np1), q2(np2,np1), q3(np2,np1),                  
     &       qOut(np2,np1)                                                      
c ----------------------------------------------------------------------        
c     loop over index 2                                                         
      DO i2 = 1, nr2                                                            
c       loop over index 1                                                       
        DO i1 = 1, nr1                                                          
          qOut(i2,i1) = q1(i2,i1) +  q2(i2,i1) + q3(i2,i1)                      
        END DO                                                                  
      END DO                                                                    
c ----------------------------------------------------------------------        
      RETURN                                                                    
      END                                                                       
c **********************************************************************        
                                                                                
c **********************************************************************        
      SUBROUTINE ADD2(flxS,flxE,fBSum,nY)                                       
c ======================================================================        
c This subroutine is auxiliary for finding the bolometric                       
c components of the scattered and emitted diffuse flux.   [MN, May'99]          
c ======================================================================        
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, iY                                                            
      DOUBLE PRECISION flxS(npL,npY), flxE(npL,npY), flxSB(npY),                
     &          flxEB(npY), fBSum(npY)                                          
c ----------------------------------------------------------------------        
      CALL Bolom(flxS,flxSB)                                                    
      CALL Bolom(flxE,flxEB)                                                    
      DO iY = 1, nY                                                             
        fBSum(iY) = flxSB(iY) + flxEB(iY)                                       
      END DO                                                                    
c ----------------------------------------------------------------------        
      RETURN                                                                    
      END                                                                       
c **********************************************************************        
                                                                                
c **********************************************************************        
      SUBROUTINE ANALINT(Nanal,xaux,yaux,m,aux,error)                           
c ======================================================================        
c This subroutine calculates integral I(x**m*y(x)*dx). Both y and x are         
c 1D arrays, y(i), x(i) with i=1,Nanal. The method used is approximation        
c of y(x) by y = P(x) + d/sqrt(1-x*x), where P(x) is the polynomial of          
c order Nanal-1, and analytic evaluation of the integral. It is assumed         
c that xaux(1)=0. Coefficients are determined from the set of Nanal             
c linear equations and subsequent call to the linear system solver              
c LINSYS.                                              [Z.I., Nov. 1995]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER iWARNING, iERROR, iCUMM                                           
      COMMON /status/ iWARNING, iERROR, iCUMM                                   
      INTEGER i, j, Nanal, error                                                
      DOUBLE PRECISION xaux(Nanal), yaux(Nanal), coeff(4), A(npY,npY),          
     &               m, aux, b                                                  
c ----------------------------------------------------------------------        
      error = 0                                                                 
c     generate matrix A and vector B                                            
      DO i = 1, Nanal                                                           
        DO j = 1, Nanal-1                                                       
          IF (xaux(i).EQ.0.0.AND.j.EQ.1) THEN                                   
            A(i,j) = 1.0                                                        
          ELSE                                                                  
            A(i,j) = xaux(i)**(1.0*j-1.0)                                       
          END IF                                                                
        END DO                                                                  
        A(i,Nanal) = 1.0/sqrt(1.0-xaux(i)*xaux(i))                              
      END DO                                                                    
                                                                                
c     solve for the coefficients                                                
      CALL LINSYS(Nanal,A,yaux,coeff,error)                                     
        IF(error.NE.0) THEN                                                     
         CALL MSG(19)                                                           
         iERROR = iERROR + 1                                                    
         RETURN                                                                 
        END IF                                                                  
c     upper limit for integration:                                              
      b = xaux(Nanal)                                                           
c     evaluate m-dependent contribution of the last term                        
      IF (m.GT.0.1) THEN                                                        
        IF (m.GT.1.1) THEN                                                      
c         this is for m=2                                                       
          aux = 0.5*(DASIN(b)-b*sqrt(1.-b*b))                                   
        ELSE                                                                    
c         this is for m=1                                                       
          aux = 1.0 - sqrt(1.-b*b)                                              
        ENDIF                                                                   
      ELSE                                                                      
c       this is for m=0                                                         
        aux = DASIN(b)                                                          
      ENDIF                                                                     
      aux = aux * coeff(Nanal)                                                  
c     add contribution from the polynom                                         
      DO i = 1, Nanal-1                                                         
        aux = aux + coeff(i) * (b**(m+1.0*i)) / (m+1.0*i)                       
      END DO                                                                    
c -----------------------------------------------------------------------       
999   RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE Attach(root,length,ext,fname)                                  
c     Attaches extensions to the root cleaned by Clean                          
c =======================================================================       
      CHARACTER*(*) root, ext, fname                                            
      INTEGER i, length                                                         
c -----------------------------------------------------------------------       
      DO i = 1, LEN(fname)                                                      
        fname(i:i) = ' '                                                        
      END DO                                                                    
      fname(:length) = root(:length)                                            
      fname(length + 1:) = ext                                                  
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE Clean(StrIn, StrOut, Length)                                   
c =======================================================================       
c     Find meaningful part of StrIn without leading and trailing junk           
c     It is returned left-justified in StrOut, right-padded with blanks         
c     The number of meaningful characters is returned in Length.                
c     In case of any problems, StrOut is empty. This sub should be used to      
c     clean every input filename immediately after DUSTY reads it. [ME,'99]     
c =======================================================================       
      CHARACTER*(*) StrIn, StrOut                                               
      INTEGER i, first, last, Length                                            
c -----------------------------------------------------------------------       
      DO i = 1, LEN(StrOut)                                                     
        StrOut(i:i) = ' '                                                       
      END DO                                                                    
      first = 1                                                                 
      last = LEN(StrIn)                                                         
      If (first.gt.last) return                                                 
c     Find end of leading junk:                                                 
      DO WHILE (StrIn(first:first).LE.' ')                                      
       first = first + 1                                                        
       if (first.gt.last) return                                                
      END DO                                                                    
c     Find start of trailing junk:                                              
      DO WHILE (StrIn(last:last).LE.' ')                                        
       last = last - 1                                                          
       if (last.lt.first) return                                                
      END DO                                                                    
c     Now trim all junk:                                                        
      StrOut = StrIn(first:last)                                                
      Length = last - first + 1                                                 
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      DOUBLE PRECISION FUNCTION Bessel(x)                                       
c =======================================================================       
c This function evaluates the Bessel function of the zeroth kind.               
c Formulae are from Abramowitz & Stegun.               [Z.I., Jan. 1997]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER i                                                                 
      DOUBLE PRECISION x, c(6), Pi                                              
c -----------------------------------------------------------------------       
      Pi = 2.0*ASIN(1.0)                                                        
      c(1) = -2.2499997                                                         
      c(2) =  1.2656208                                                         
      c(3) = -0.3163866                                                         
      c(4) =  0.0444479                                                         
      c(5) = -0.0039444                                                         
      c(6) =  0.00021                                                           
      Bessel=0.0                                                                
      IF (x.LE.3.0)THEN                                                         
        DO i=1,6                                                                
          Bessel = Bessel + c(i)*(x/3.0)**(2.0*i)                               
        END DO                                                                  
        Bessel = 1.0 + Bessel                                                   
        ELSE                                                                    
        Bessel = sqrt(2.0/Pi/x) * dcos(x-Pi/4.0)                                
      ENDIF                                                                     
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
c This subroutine obtained from prof. P. Menguc, Dept. of Mechanical            
c Eng., University of Kentucky.                        [Z.I., Aug. 1996]        
c -----------------------------------------------------------------------       
c     __________________________________________________________________        
c                                                                               
c     SUBROUTINE BHMIE CALCULATES AMPLITUDE SCATTERING MATRIX ELEMENTS          
C     & EFFICIENCIES FOR EXTINCTION, TOTAL SCATTERING AND BACSCATTERING,        
C     FOR A GIVEN SIZE PARAMETER AND RELATIVE REFRACTIVE INDEX                  
C     __________________________________________________________________        
C                                                                               
      subroutine bhmie (x,refrel,nang,s1,s2,qext,qsca,qback)                    
      dimension amu(100),theta(100),pi(100),tau(100),pi0(100),pi1(100)          
      complex d(3000),y,refrel,xi,xi0,xi1,an,bn,s1(200),s2(200)                 
      double precision psi0,psi1,psi,dn,dx                                      
      dx=x                                                                      
      y=x*refrel                                                                
c     ___________________________________________________________________       
c     series terminated after nstop terms                                       
c     ___________________________________________________________________       
      xstop=x+4.*x**.3333 +2.0                                                  
      nstop=xstop                                                               
      ymod=cabs(y)                                                              
      nmx=amax1(xstop,ymod) + 15                                                
      dang=1.570796327/float(nang-1)                                            
      do 555 j = 1,nang                                                         
      theta(j)= (float(j)-1.)*dang                                              
555   amu(j)=cos(theta(j))                                                      
c     __________________________________________________________________        
c     logarithmic derivative d(j) calculated by downward recurrence             
c     beginning with initial value 0.0 + i*0.0 at j = nmx                       
c     __________________________________________________________________        
c                                                                               
      d(nmx)=cmplx(0.0,0.0)                                                     
      nn=nmx-1                                                                  
      do 120 n=1,nn                                                             
      rn=nmx-n+1                                                                
      d(nmx-n)=(rn/y)-(1./(d(nmx-n+1)+rn/y))                                    
120   continue                                                                  
      do 666 j=1,nang                                                           
      pi0(j)=0.0                                                                
      pi1(j)=1.0                                                                
666   continue                                                                  
      nn=2*nang-1                                                               
      do 777 j=1,nn                                                             
      s1(j)=cmplx(0.0,0.0)                                                      
      s2(j)=cmplx(0.0,0.0)                                                      
777   continue                                                                  
c     __________________________________________________________________        
c     riccati bessel functions with real argument x calculated by upward        
c     recurrence                                                                
c     __________________________________________________________________        
c                                                                               
      psi0=cos(dx)                                                              
      psi1=sin(dx)                                                              
      chi0=-sin(x)                                                              
      chi1=cos(x)                                                               
      apsi0=psi0                                                                
      apsi1=psi1                                                                
      xi0=cmplx(apsi0,-chi0)                                                    
      xi1=cmplx(apsi1,-chi1)                                                    
      qsca=0.0                                                                  
      n=1                                                                       
200   dn=n                                                                      
      rn=n                                                                      
      fn=(2.*rn+1.)/(rn*(rn+1.))                                                
      psi=(2.*dn-1.)*psi1/dx-psi0                                               
      apsi=psi                                                                  
      chi=(2.*rn-1.)*chi1/x -  chi0                                             
      xi = cmplx(apsi,-chi)                                                     
      an=(d(n)/refrel+rn/x)*apsi - apsi1                                        
      an=an/((d(n)/refrel+rn/x)*xi - xi1)                                       
      bn=(refrel *d(n)+rn/x)*apsi - apsi1                                       
      bn=bn/((refrel*d(n)+rn/x)*xi - xi1)                                       
      qsca=qsca+(2.*rn+1.)*(cabs(an)*cabs(an)+cabs(bn)*cabs(bn))                
      do 789 j=1,nang                                                           
      jj=2*nang-j                                                               
      pi(j)=pi1(j)                                                              
      tau(j)=rn*amu(j)*pi(j) - (rn+1.)*pi0(j)                                   
      p=(-1.)**(n-1)                                                            
      s1(j)=s1(j)+fn*(an*pi(j)+bn*tau(j))                                       
      t=(-1.)**n                                                                
      s2(j)=s2(j) + fn*(an*tau(j)+bn*pi(j))                                     
      if (j .eq. jj) go to 789                                                  
      s1(jj)=s1(jj) + fn*(an*pi(j)*p + bn*tau(j)*t)                             
      s2(jj)=s2(jj) + fn*(an*tau(j)*t + bn*pi(j)*p)                             
789   continue                                                                  
      psi0=psi1                                                                 
      psi1=psi                                                                  
      apsi1=psi1                                                                
      chi0=chi1                                                                 
      chi1=chi                                                                  
      xi1=cmplx(apsi1,-chi1)                                                    
      n=n+1                                                                     
      rn=n                                                                      
      do 999 j=1,nang                                                           
      pi1(j)=((2.*rn-1.)/(rn-1.))*amu(j)*pi(j)                                  
      pi1(j)=pi1(j) - rn*pi0(j)/(rn-1.)                                         
      pi0(j) = pi(j)                                                            
999   continue                                                                  
      if (n-1-nstop) 200, 300, 300                                              
300   qsca=(2./(x*x))*qsca                                                      
      qext=(4./(x*x))*real(s1(1))                                               
      qback=(4./(x*x))*cabs(s1(2*nang -1))*cabs(s1(2*nang -1))                  
      return                                                                    
      end                                                                       
c ***********************************************************************       
                                                                                
c***********************************************************************        
      SUBROUTINE CHKRANGE(dr,x)                                                 
c=======================================================================        
c This subroutine checks if x is within the allowed range defined by            
c dr<<1:                                                                        
c         dr**2 < x < 1/dr**2                                                   
c If it is not then x = 0.0                            [Z.I., Jan. 1997]        
c=======================================================================        
      IMPLICIT none                                                             
      DOUBLE PRECISION x, dr                                                    
c-----------------------------------------------------------------------        
      IF ((x-dr*dr)*(x-1./dr/dr).LT.0.0) THEN                                   
        continue                                                                
      ELSE                                                                      
c        continue                                                               
        x = 0.0                                                                 
      END IF                                                                    
c-----------------------------------------------------------------------        
      RETURN                                                                    
      END                                                                       
c***********************************************************************        
                                                                                
c ***********************************************************************       
      SUBROUTINE CHKSPLIN(x,fun,funmid,N,coef,maxerr,RDW)                       
c ======================================================================        
c This subroutine checks the spline coefficients coef(i,j):                     
c fun(x)=coef(i,1) + coef(i,2)*x + coef(i,3)*x^2 + coef(i,4)*x^3,               
c for x(i).LE.x.LE.x(i+1) with i=1..N. Array funmid(1..N-1) contains the        
c values of function fun at mid points defined as                               
c xmid(i)=SQRT(x(i)*x(i+1). If spline approximation produces error              
c greater than maxerr, or funmid<0, a straight line is produced between         
c x(i) and x(i+1).                                                              
c ======================================================================        
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER N, i, iC                                                          
      DOUBLE PRECISION x(npY), fun(npY), funmid(npY), coef(npY,4),              
     &       maxerr, error, slope, xmid, funSpline, aux, power, yR, yL          
      LOGICAL RDW                                                               
c ---------------------------------------------------------                     
c       check the midpoints                                                     
        DO i = 1, N - 1                                                         
          xmid = dsqrt(x(i)*x(i+1))                                             
          funSpline = 0.0                                                       
          DO iC=1,4                                                             
            IF (xmid.EQ.0.0.AND.iC.EQ.1) THEN                                   
              aux = 1.0                                                         
              ELSE                                                              
              aux = xmid**(float(iC)-1.0)                                       
            END IF                                                              
            funSpline = funSpline + coef(i,iC)*aux                              
          END DO                                                                
          error = DABS((funSpline-funmid(i))/funmid(i))                         
c         check for the deviation at the midpoint                               
          IF (error.GE.maxerr.OR.funSpline.LE.0.0) THEN                         
            slope = (fun(i+1) - fun(i)) / (x(i+1)-x(i))                         
            coef(i,1) = fun(i) - x(i) * slope                                   
            coef(i,2) = slope                                                   
            coef(i,3) = 0.0                                                     
            coef(i,4) = 0.0                                                     
          END IF                                                                
c         check for the logarithmic derivative (only for RDW)                   
          IF(RDW) THEN                                                          
            yL = fun(i)                                                         
            yR = fun(i+1)                                                       
            IF (x(i)*x(i+1).GT.0.AND.yL*yR.GT.0) THEN                           
              power = log(yR/yL)/log(x(i+1)/x(i))                               
              IF (abs(power).GT.10.) THEN                                       
                slope = (yR - yL) / (x(i+1)-x(i))                               
                coef(i,1) = yL - x(i) * slope                                   
                coef(i,2) = slope                                               
                coef(i,3) = 0.0                                                 
                coef(i,4) = 0.0                                                 
              END IF                                                            
            END IF                                                              
          END IF                                                                
        END DO                                                                  
c ---------------------------------------------------------                     
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE CHKSPLINold(x,fun,funmid,N,coef,maxerr)                        
c =======================================================================       
c This subroutine checks the spline coefficients coef(i,j):                     
c fun(x)=coef(i,1) + coef(i,2)*x + coef(i,3)*x^2 + coef(i,4)*x^3,               
c for x(i).LE.x.LE.x(i+1) with i=1..N. Array funmid(1..N-1) contains the        
c values of function fun at mid points defined as                               
c xmid(i)=SQRT(x(i)*x(i+1). If spline approximation produces error              
c greater than maxerr, or funmid<0, a straight line is produced between         
c x(i) and x(i+1).                                     [Z.I., Feb. 1995]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER N, i, iC                                                          
      DOUBLE PRECISION x(npY), fun(npY), funmid(npY), coef(npY,4), aux,         
     &     funSpline, maxerr, error, slope, xmid                                
c -----------------------------------------------------------------------       
c       check the midpoints                                                     
        DO i = 1, N - 1                                                         
          xmid = dsqrt(x(i)*x(i+1))                                             
          funSpline = 0.0                                                       
          DO iC=1,4                                                             
            IF (xmid.EQ.0.0.AND.iC.EQ.1) THEN                                   
              aux = 1.0                                                         
              ELSE                                                              
              aux = xmid**(float(iC)-1.0)                                       
            END IF                                                              
            funSpline = funSpline + coef(i,iC)*aux                              
          END DO                                                                
          error = DABS((funSpline-funmid(i))/funmid(i))                         
          IF (error.GE.maxerr.OR.funSpline.LE.0.0) THEN                         
            slope = (fun(i+1) - fun(i)) / (x(i+1)-x(i))                         
            coef(i,1) = fun(i) - x(i) * slope                                   
            coef(i,2) = slope                                                   
            coef(i,3) = 0.0                                                     
            coef(i,4) = 0.0                                                     
          END IF                                                                
        END DO                                                                  
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      INTEGER FUNCTION EMPTY(line)                                              
c =======================================================================       
c This function is 1 if string 'line' is empty, or if it contains only          
c '%', and 0 otherwise.                                                         
c                                                      [Z.I., Nov. 1996]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER i, iTeX, l                                                        
      CHARACTER ch                                                              
      CHARACTER*(*) line                                                        
c -----------------------------------------------------------------------       
      l = LEN(line)                                                             
      EMPTY = 1                                                                 
      iTeX = 0                                                                  
      DO i = 1, l                                                               
        ch = line(i:i)                                                          
        IF(EMPTY.EQ.1.AND.ch.EQ.'%') iTeX = 1                                   
         IF (ch.NE.' ') EMPTY = 0                                               
      END DO                                                                    
      IF (iTeX.EQ.1) EMPTY = 1                                                  
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE FileMSG(fname,strg)                                            
c =======================================================================       
c     Prints a message in *.out file in case of error opening the user          
c     supplied files.                                                           
c =======================================================================       
      IMPLICIT NONE                                                             
      CHARACTER aux*230, strg*(*), fname*(*)                                    
      INTEGER length, Empty                                                     
c -----------------------------------------------------------------------       
1     read(1,'(a)') aux                                                         
      IF (Empty(aux).EQ.1) goto 1                                               
      CALL Clean(aux,fname,length)                                              
                                                                                
      open(10, ERR=100, FILE=fname, STATUS='OLD') 	                             
      close(10)                                                                 
	RETURN                                                                         
                                                                                
100   write(12,*)' *** FATAL ERROR IN DUSTY! **************************'        
      write(12,*)' File with the ',strg                                         
      write(12,'(a2,a)')'  ',fname                                              
      write(12,*)' is missing ?!'                                               
      write(12,*)' ****************************************************'        
      close(12)                                                                 
c -----------------------------------------------------------------------       
      STOP                                                                      
      END                                                                       
c ***********************************************************************       
                                                                                
c *************************************************************************     
      SUBROUTINE FindErr(fbol,maxFerr,nY)                                       
c ========================================================================      
c This subroutine finds maximum err in flux conservation for both               
c spherical and slab case as (fmax-fmin)/(fmax+fmin)   [MN,Aug'99]              
c =========================================================================     
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      INTEGER startyp(2), Nlamtr(2), nBB(2)                                     
      CHARACTER nameStar(2)*235                                                 
      DOUBLE PRECISION Tstar, lamtr(2,101), klam(2,100), Tbb(2,10),             
     &                 rellum(2,10), mu1, ksi, mu2, xSiO, r1rs                  
      COMMON /source/ Tstar, lamtr, klam, Tbb, rellum, mu1, ksi, mu2,           
     &                xSiO, r1rs, startyp, Nlamtr, nBB, nameStar                
      INTEGER nY, iY                                                            
      DOUBLE PRECISION fbol(npY), maxFerr, fmin, fmax, aux                      
c -----------------------------------------------------------------------       
c     Find the min and max of fbol values                                       
c     The abs and lower limit on fbol are protection for the case               
c     of completely symmetric slab illumination. The lower limit                
c     is bound by the numerical accuracy of the flux calculation                
        fmin = 1.e5                                                             
        fmax = 0.                                                               
        DO iY = 1, nY                                                           
           aux = fbol(iY)                                                       
           IF (ksi.eq.1.0) aux = dabs(aux)                                      
           IF (dabs(aux).LE.accFbol) aux = accFbol                              
           IF(aux.LT.fmin) fmin = aux                                           
           IF(aux.GT.fmax) fmax = aux                                           
        END DO                                                                  
        if (fmax.LT.0.) then                                                    
c     bad solution; overall flux cannot be negative                             
            maxFerr = 1                                                         
        else                                                                    
            maxFerr = (fmax - fmin)/(fmax + dabs(fmin))                         
        end if                                                                  
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c *************************************************************************     
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE FindMax(NN,i1,i2,A,Amax)                                       
c =======================================================================       
c This subroutine finds maximum values, Amax, of an array A(nY) between         
c values A(i1) and A(i2).                              [Z.I., Jul. 1996]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER NN, i1, i2, i                                                     
      DOUBLE PRECISION A(NN), Amax                                              
c -----------------------------------------------------------------------       
      Amax = A(i1)                                                              
c     loop over radial positions                                                
      DO i = i1, i2                                                             
        IF (A(i).GT.Amax) Amax = A(i)                                           
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c =======================================================================       
                                                                                
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE FindRMS(typ,X,val,accur,N)                                     
c =======================================================================       
c Finds relative deviations 'accur' of an array X(N) from a given value val.    
c For typ=1 accur is maximal deviation, and for typ=2 the rms deviation.        
c                                                         [ZI'95; MN'99]        
c =======================================================================       
      IMPLICIT NONE                                                             
      INTEGER N, i, typ                                                         
      DOUBLE PRECISION X(N), val, accur, ss, dev                                
c -----------------------------------------------------------------------       
      IF (typ.EQ.1) THEN                                                        
        accur = 0.0                                                             
        DO i = 1, N                                                             
          dev = (X(i)-val)/val                                                  
          IF (DABS(dev).GT.accur) accur = DABS(dev)                             
        END DO                                                                  
      ELSE                                                                      
        ss = 0.0                                                                
        DO i = 1, N                                                             
          dev = X(i)-val                                                        
          ss = ss + dev*dev                                                     
        END DO                                                                  
        accur = sqrt(ss/N/(N-1.))                                               
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE GetFS(xx,nm,flag,str)                                          
c =======================================================================       
c This subroutine writes number xx to a string str according to a format        
c f?.nm. Here ? stands for the number of needed places. A blank is              
c inserted at the beginning, and for flag.NE.1 another one if number is         
c positive. If xx<0 second blank is replaced by '-'. For example, for           
c flag=0 and xx = -0.1234E+02, calling this subroutine with nm=1 will           
c result in str = ' -12.3', while xx = 0.0123 with nm=3 gives '  0.012'.        
c If flag=1 minus will be ignored, for example xx = -0.1234E+02 and nm=1        
c will result in str = ' 12.3',                        [Z.I., Nov. 1996]        
c =======================================================================       
      IMPLICIT none                                                             
      CHARACTER ch                                                              
      CHARACTER*(*) str                                                         
      INTEGER  flag, nm, db, i, d(20), j, k, dnmp1                              
      DOUBLE PRECISION xx, x, rest                                              
c -----------------------------------------------------------------------       
      DO i = 1, len(str)                                                        
        str(i:i) = ' '                                                          
      END DO                                                                    
                                                                                
      x = xx                                                                    
      str(1:1) = ' '                                                            
      i = 2                                                                     
      IF (flag.NE.1) THEN                                                       
         IF (x.LT.0.0) THEN                                                     
           str(i:i) = '-'                                                       
         ELSE                                                                   
           str(i:i) = ' '                                                       
         END IF                                                                 
         i = i + 1                                                              
      END IF                                                                    
      IF (x.LT.0.0) x = -x                                                      
c     first check if x will have to be rounded up                               
c     find (nm+1)-th decimal digit                                              
      dnmp1 = int(x*10.**(nm+1)-int(x*10.**nm)*10.)                             
      IF (dnmp1.GE.5) x = x + 1./10.0**nm                                       
      IF (x.GE.1.0) THEN                                                        
c       number of digits before the decimal sign                                
        db = int(log10(x)) + 1                                                  
c       copy all these digits to str                                            
        DO j = 1, db                                                            
          rest = x                                                              
          IF (j.GT.1) THEN                                                      
            DO k = 1, j-1                                                       
              rest = rest - d(k)*10.**(db-k)                                    
            END DO                                                              
          END IF                                                                
          d(j) = int(rest/10.**(db-j))                                          
          write(ch,'(i1)')d(j)                                                  
          str(i:i) = ch                                                         
          i = i + 1                                                             
        END DO                                                                  
        rest = rest - d(db)                                                     
        IF (nm.GT.0) THEN                                                       
          str(i:i) = '.'                                                        
          i = i + 1                                                             
        END IF                                                                  
      ELSE                                                                      
        str(i:i) = '0'                                                          
        i = i + 1                                                               
        IF (nm.GT.0) THEN                                                       
          str(i:i) = '.'                                                        
          i = i + 1                                                             
        END IF                                                                  
        rest = x                                                                
      END IF                                                                    
c     now copy all nm remaining decimal digits to str                           
      IF (nm.GT.0) THEN                                                         
        DO j = 1, nm                                                            
          d(j) = int(rest*10.**j)                                               
          IF (j.GT.1) THEN                                                      
            DO k = 1, j-1                                                       
              d(j)=d(j)-int(d(k)*10.**(j-k))                                    
            END DO                                                              
          END IF                                                                
          write(ch,'(i1)')d(j)                                                  
          str(i:i) = ch                                                         
          i = i + 1                                                             
        END DO                                                                  
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      DOUBLE PRECISION FUNCTION H(x1,x2)                                        
c =======================================================================       
c This function calculates the step function: H=1 for x1 >= x2 and H=0          
c for x1 < x2.                                         [Z.I., Nov. 1995]        
c =======================================================================       
      IMPLICIT none                                                             
      DOUBLE PRECISION x1, x2                                                   
c -----------------------------------------------------------------------       
      IF (x1.GE.x2) THEN                                                        
        H = 1.0                                                                 
      ELSE                                                                      
        H = 0.0                                                                 
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      INTEGER FUNCTION Kron(i1,i2)                                              
c =======================================================================       
c This function is Kronecker delta-function defined as:                         
c Kron(i1,i2) = 1 for i1=i2                                                     
c Kron(i1,i2) = 0 otherwise.                           [Z.I., Dec. 1995]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER i1, i2                                                            
c -----------------------------------------------------------------------       
      IF (i1.EQ.i2) THEN                                                        
        Kron = 1                                                                
      ELSE                                                                      
        Kron = 0                                                                
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE LINE(com,typ,unt)                                              
c =======================================================================       
c This subroutine writes a line into file open as unt. For type = 1             
c the line is '---', and for type = 2 '==='.If com=1 a comment sign # is        
c added in the beginning (this is when line is used in file headers)            
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER com, typ, unt                                                     
c -----------------------------------------------------------------------       
      IF(typ.EQ.1) THEN                                                         
       IF(com.eq.1) THEN                                                        
        write(unt,'(a50)')                                                      
     &   '# ------------------------------------------------'                   
       ELSE                                                                     
        write(unt,*)'--------------------------------------------------'        
       END IF                                                                   
      ELSE                                                                      
       IF(com.eq.1) THEN                                                        
        write(unt,'(a50)')                                                      
     &   '# ================================================'                   
       ELSE                                                                     
        write(unt,*)'=================================================='        
       END IF                                                                   
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE LinInter(NN,N,x,y,xloc,iNloc,yloc)                             
c =======================================================================       
c This subroutine performs linear interpolation for y(x) such that              
c yloc = y(xloc). It is assumed that x is monotonously increasing.              
c                                                      [Z.I., Mar. 1996]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER NN, N, i, istop, iNloc                                            
      DOUBLE PRECISION x(NN), y(NN), xloc, yloc                                 
c -----------------------------------------------------------------------       
      IF (N.GT.1) THEN                                                          
        IF ((x(1)-xloc)*(x(N)-xloc).LE.0.0) THEN                                
          istop = 0                                                             
          i = 1                                                                 
          DO WHILE (istop.NE.1)                                                 
            i = i + 1                                                           
            IF (i.GT.N) stop 'LinInter ???'                                     
            IF (x(i).GE.xloc) THEN                                              
              istop = 1                                                         
              iNloc = i                                                         
              yloc = y(i-1) + (y(i)-y(i-1))/(x(i)-x(i-1))*(xloc-x(i-1))         
            END IF                                                              
          END DO                                                                
          ELSE                                                                  
          IF (xloc.LE.x(1)) yloc = y(1)                                         
          IF (xloc.GE.x(N)) yloc = y(N)                                         
        END IF                                                                  
        ELSE                                                                    
        yloc = y(1)                                                             
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE LINSYS(Nreal,A,B,X,error)                                      
c =======================================================================       
c This subroutine solves the set of linear equations [A]*[X] = [B] for          
c X [A(k,1)*X(1)+A(k,2)*X(2)+...+A(k,Nreal)*X(Nreal) = B(k), k=1,Nreal).        
c The real size of matrix A is Nreal x Nreal and its physical dimension         
c is npY x npY, where npY comes from INCLUDE 'userpar.inc'. Both vectors        
c B and X have real lengths Nreal. The set is solved by calls to LUDCMP         
c and LUBKSB and the solution is improved subsequently by a call to             
c MPROVE. These three subroutines are taken from Numerical Recipes.             
c                                                      [Z.I., Nov. 1995]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER Nreal, indx(npY), i, j, error                                     
      DOUBLE PRECISION A(npY,npY), B(npY), X(npY)                               
      DOUBLE PRECISION A1(npY,npY), B1(npY), A2(npY,npY), B2(npY), d            
c -----------------------------------------------------------------------       
      error = 0                                                                 
c generate DOUBLE PRECISION copies of A and B (two copies because they          
c are changed in LUDCMP and LUBKSB, but still needed for MPROVE)                
      DO i = 1, Nreal                                                           
        B1(i) = B(i)                                                            
        B2(i) = B(i)                                                            
        DO j = 1, Nreal                                                         
           A1(i,j) = A(i,j)                                                     
           A2(i,j) = A(i,j)                                                     
        END DO                                                                  
      END DO                                                                    
c     solve the system                                                          
      CALL LUDCMP(A1,Nreal,npY,indx,d,error)                                    
      IF (error.NE.0) RETURN                                                    
      CALL LUBKSB(A1,Nreal,npY,indx,B1)                                         
c     improve the solution (saved in B)                                         
      CALL MPROVE(A2,A1,Nreal,npY,indx,B2,B1)                                   
c     copy the improved solution to output vector X                             
      DO i = 1, Nreal                                                           
        X(i) = B1(i)                                                            
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)                                          
c =======================================================================       
      DIMENSION INDX(NP)                                                        
      DOUBLE PRECISION A(NP,NP),B(NP)                                           
c -------------------------------------------------------------------           
      II=0                                                                      
      DO 12 I=1,N                                                               
      LL=INDX(I)                                                                
      SUM=B(LL)                                                                 
      B(LL)=B(I)                                                                
      IF (II.NE.0)THEN                                                          
        DO 11 J=II,I-1                                                          
          SUM=SUM-A(I,J)*B(J)                                                   
11        CONTINUE                                                              
      ELSE IF (SUM.NE.0.) THEN                                                  
        II=I                                                                    
      ENDIF                                                                     
      B(I)=SUM                                                                  
12    CONTINUE                                                                  
      DO 14 I=N,1,-1                                                            
      SUM=B(I)                                                                  
      IF(I.LT.N)THEN                                                            
        DO 13 J=I+1,N                                                           
          SUM=SUM-A(I,J)*B(J)                                                   
13        CONTINUE                                                              
      ENDIF                                                                     
      B(I)=SUM/A(I,I)                                                           
14    CONTINUE                                                                  
c -------------------------------------------------------------------           
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE LUDCMP(A,N,NP,INDX,D,error)                                    
c =======================================================================       
      PARAMETER (NMAX=10000,TINY=1.0E-20)                                       
      DIMENSION INDX(NP)                                                        
      INTEGER error                                                             
      DOUBLE PRECISION A(NP,NP),VV(NMAX), D, SUM                                
c -------------------------------------------------------------------           
      error = 0                                                                 
      D = 1.                                                                    
      DO I = 1, N                                                               
       AAMAX=0.                                                                 
       DO J = 1, N                                                              
        IF (DABS(A(I,J)).GT.AAMAX) AAMAX=DABS(A(I,J))                           
       END DO                                                                   
c       IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'                               
       IF (AAMAX.EQ.0.) THEN                                                    
        error = 5                                                               
        RETURN                                                                  
       ENDIF                                                                    
       VV(I)=1./AAMAX                                                           
      END DO                                                                    
      DO J = 1 , N                                                              
       IF (J.GT.1) THEN                                                         
        DO I = 1, J-1                                                           
          SUM=A(I,J)                                                            
          IF (I.GT.1)THEN                                                       
            DO K = 1, I-1                                                       
             SUM=SUM-A(I,K)*A(K,J)                                              
            END DO                                                              
            A(I,J)=SUM                                                          
          ENDIF                                                                 
        END DO                                                                  
       ENDIF                                                                    
       AAMAX=0.                                                                 
       DO I = J, N                                                              
        SUM=A(I,J)                                                              
        IF (J.GT.1)THEN                                                         
          DO K = 1, J-1                                                         
            SUM=SUM-A(I,K)*A(K,J)                                               
          END DO                                                                
          A(I,J)=SUM                                                            
        ENDIF                                                                   
        DUM=VV(I)*DABS(SUM)                                                     
        IF (DUM.GE.AAMAX) THEN                                                  
          IMAX=I                                                                
          AAMAX=DUM                                                             
        ENDIF                                                                   
       END DO                                                                   
       IF (J.NE.IMAX)THEN                                                       
        DO K = 1, N                                                             
          DUM=A(IMAX,K)                                                         
          A(IMAX,K)=A(J,K)                                                      
          A(J,K)=DUM                                                            
        END DO                                                                  
        D=-D                                                                    
        VV(IMAX)=VV(J)                                                          
       ENDIF                                                                    
       INDX(J)=IMAX                                                             
       IF(J.NE.N)THEN                                                           
        IF(A(J,J).EQ.0.)A(J,J)=TINY                                             
        DUM=1./A(J,J)                                                           
        DO I = J+1, N                                                           
          A(I,J)=A(I,J)*DUM                                                     
        END DO                                                                  
       ENDIF                                                                    
      END DO                                                                    
      IF(A(N,N).EQ.0.)A(N,N)=TINY                                               
c -------------------------------------------------------------------           
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE MakeTable(Elems,rows,cols,unt)                                 
c =======================================================================       
c     This is an auxiliary subroutine for print out of tables                   
c     of Elems(cols,rows) in output unit 'unt'. rows = max{npL,npY}.            
c     This array is defined in PrOut as well and if you change its size         
c     make sure you do this in both places.                 [MN, Mar'98]        
c =======================================================================       
      IMPLICIT NONE                                                             
      INTEGER rows, cols, unt, k, i                                             
      DOUBLE PRECISION Elems(25,200)                                            
c -----------------------------------------------------------------------       
      DO i = 1, rows                                                            
        write(unt,'(1p,21E11.3)') (Elems(k,i),k=1,cols)                         
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE Maple3(w,z,p,MpInt)                                            
c =======================================================================       
c This function calculates indefinite integral:                                 
c    MpInt(iC) = INT(w^(2-iC) / sqrt(w^2-p^2) * dw), for iC=1,2,3,4.            
c                                                      [Z.I., Apr. 1996]        
c =======================================================================       
      IMPLICIT none                                                             
      DOUBLE PRECISION w, z, p, MpInt(4)                                        
c -----------------------------------------------------------------------       
c     integrals                                                                 
      MpInt(1) = z                                                              
      MpInt(2) = dlog(w+z)                                                      
      IF (p.GT.0.0) THEN                                                        
        MpInt(3) = dacos(p/w)/p                                                 
        MpInt(4) = z/w/p/p                                                      
        ELSE                                                                    
        MpInt(3) = -1.0 / w                                                     
        MpInt(4) = -0.5 / w / w                                                 
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE midsql(funk,aa,bb,s,n)                                         
c =======================================================================       
      INTEGER n                                                                 
      DOUBLE PRECISION aa,bb,s,funk                                             
      EXTERNAL funk                                                             
      INTEGER it,j                                                              
      DOUBLE PRECISION ddel,del,sum,tnm,x,func,a,b                              
c -----------------------------------------------------------------------       
      func(x)=2.*x*funk(aa+x**2)                                                
      b=dsqrt(bb-aa)                                                            
      a=0.                                                                      
      if (n.eq.1) then                                                          
        s=(b-a)*func(0.5*(a+b))                                                 
      else                                                                      
        it=3**(n-2)                                                             
        tnm=it                                                                  
        del=(b-a)/(3.*tnm)                                                      
        ddel=del+del                                                            
        x=a+0.5*del                                                             
        sum=0.                                                                  
        do 11 j=1,it                                                            
          sum=sum+func(x)                                                       
          x=x+ddel                                                              
          sum=sum+func(x)                                                       
          x=x+del                                                               
11      continue                                                                
        s=(s+(b-a)*sum/tnm)/3.                                                  
      endif                                                                     
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE MPROVE(A,ALUD,N,NP,INDX,B,X)                                   
c =======================================================================       
      PARAMETER (NMAX=10000)                                                    
      DIMENSION INDX(N)                                                         
      DOUBLE PRECISION SDP,A(NP,NP),ALUD(NP,NP),B(N),X(N),R(NMAX)               
c -----------------------------------------------------------------------       
      DO i = 1, N                                                               
       SDP = -B(i)                                                              
       DO j = 1, N                                                              
         SDP = SDP + A(i,j)*X(j)                                                
       END DO                                                                   
       R(i) = SDP                                                               
      END DO                                                                    
      CALL LUBKSB(ALUD,N,NP,INDX,R)                                             
      DO i = 1, N                                                               
       X(i) = X(i) - R(i)                                                       
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE MSG(msgno)                                                     
c =======================================================================       
c This subroutine writes runtime messages to auxiliary file fname.m##           
c or to the output file fname.out.             [ZI,Feb'96; MN,Jul'99]           
c =======================================================================       
      IMPLICIT none                                                             
      CHARACTER*100 zline(999)                                                  
      INTEGER iINP, iSUM, iOUT, iVerb, iSPP, iA, iB, iC, iX, NlambdaOut,        
     &         iInn, iPsf, iV, Nconv, Nvisi                                     
      DOUBLE PRECISION LambdaOut(20), ConvInt(20,1000), Visib(20,1000),         
     &       Offset(1000), qtheta1(1000), Te_min                                
      COMMON /output/ LambdaOut, ConvInt, Visib, Offset, qtheta1,               
     &      Te_min, iPSF, NlambdaOut, iINP, iSUM, iOUT, iVerb, iSPP,            
     &      iA, iB, iC, iX, iInn, iV, Nconv, Nvisi, zline                       
      INTEGER  msgno                                                            
c -----------------------------------------------------------------------       
      IF (msgno.EQ.1.AND.iX.GT.0) THEN                                          
       write(18,*)' ************  WARNING  ************'                        
       write(18,*)' Temperature calculation in FindTemp'                        
       write(18,*)' achieved the limit of 500 iterations'                       
      END IF                                                                    
      IF (msgno.EQ.2.AND.iX.GT.0) THEN                                          
       write(18,*)' ************  WARNING  ************'                        
       write(18,*)' Energy density iterations in RADTRANSF'                     
       write(18,*)' achieved the limit of 10000 iterations'                     
      END IF                                                                    
      IF (msgno.EQ.3) THEN                                                      
       write(12,*)' ********** A BIG ERROR ***********'                         
       write(12,*)' * Denstyp is not between 1 and 7!*'                         
       write(12,*)' * Check input file and try again *'                         
       write(12,*)' ***********************************'                        
      END IF                                                                    
      IF (msgno.EQ.4.AND.iX.GT.0) THEN                                          
       write(18,*)' ************  WARNING  *****************'                   
       write(18,*)' Could not bracket in Zbrac (in sub FindTemp)'               
       write(18,*)' Something might be wrong in your input.     '               
      END IF                                                                    
      IF (msgno.EQ.5) THEN                                                      
       write(12,*)' ********** A BIG ERROR ***********'                         
       write(12,*)' *  All abundances must be >= 0!  *'                         
       write(12,*)' * Check input file and try again *'                         
       write(12,*)' ***********************************'                        
      END IF                                                                    
      IF (msgno.EQ.6) THEN                                                      
       write(12,*)' ********** A BIG ERROR ***********'                         
       write(12,*)' * Wavelengths for the power-law  *'                         
       write(12,*)' * spectrum must be ascending!    *'                         
       write(12,*)' * Check input file and try again *'                         
       write(12,*)' ***********************************'                        
      END IF                                                                    
      IF (msgno.EQ.7) THEN                                                      
       write(12,*)' ********** A BIG ERROR ***********'                         
       write(12,*)' * Relative luminosities must add *'                         
       write(12,*)' * up to a number >0!!!           *'                         
       write(12,*)' * Check input file and try again *'                         
       write(12,*)' ***********************************'                        
      END IF                                                                    
      IF (msgno.EQ.8) THEN                                                      
       write(12,*)' ********** A BIG ERROR ***********'                         
       write(12,*)' *    A black body temperature    *'                         
       write(12,*)' *        should be > 0 !!!       *'                         
       write(12,*)' * Check input file and try again *'                         
       write(12,*)' ***********************************'                        
      END IF                                                                    
      IF (msgno.EQ.9) THEN                                                      
       write(12,*)' ********** A BIG ERROR ***********'                         
       write(12,*)' * Flag for optical properties    *'                         
       write(12,*)' * should be between 1 and 3!!!   *'                         
       write(12,*)' * Check input file and try again *'                         
       write(12,*)' ***********************************'                        
      END IF                                                                    
      IF (msgno.EQ.10) THEN                                                     
       write(12,*)' ********** A BIG ERROR ***********'                         
       write(12,*)' * Flag for size distribution     *'                         
       write(12,*)' * should be between 1 and 3!!!   *'                         
       write(12,*)' * Check input file and try again *'                         
       write(12,*)' ***********************************'                        
      END IF                                                                    
      IF (msgno.EQ.11) THEN                                                     
       write(12,*)' ********** A BIG ERROR ***********'                         
       write(12,*)' * The flag for external spectrum *'                         
       write(12,*)' * should be between 1 and 6 !!!  *'                         
       write(12,*)' * Check input file and try again *'                         
       write(12,*)' ***********************************'                        
      END IF                                                                    
      IF (msgno.EQ.12) THEN                                                     
       write(12,*)' ***  FATAL ERROR IN DUSTY  **********'                      
       write(12,*)' Only three types of the point spread '                      
       write(12,*)' function are allowed: 1, 2 or 3 !!!  '                      
       write(12,*)' Check input file and try again       '                      
       write(12,*)' ***********************************'                        
      END IF                                                                    
c     msg 14 is not called in this version.                                     
      IF (msgno.EQ.14.AND.iX.GT.0) THEN                                         
       write(18,*)' ******** MESSAGE FROM SLBSolve *******'                     
       write(18,*)' Convergence on en.density is too slow.'                     
       write(18,*)' If the accuracy is not reached yet    '                     
       write(18,*)' will increase grid size and try again.'                     
       write(18,*)' **************************************'                     
      END IF                                                                    
      IF (msgno.EQ.15) THEN                                                     
       write(12,*) ' **************** WARNING ******************'               
       write(12,*) '  NO calculation in this case. Parameter npY'               
       write(12,*) '  needs to be at least 50. Use of the slab  '               
       write(12,*) '  parameters is suggested (see userpar.inc) '               
       write(12,*) ' *******************************************'               
      END IF                                                                    
      IF (msgno.EQ.16) THEN                                                     
       write(12,*)' ****************  WARNING  ******************'              
       write(12,*)'  The density profile Eta is too steep and the'              
       write(12,*)'  code can not handle this. Try decreasing the'              
       write(12,*)'  outer radius Y (see Manual, 3.3.3).         '              
       write(12,*)' *********************************************'              
       IF(iX.GT.0) THEN                                                         
       write(18,*)' ****************  WARNING  ******************'              
       write(18,*)'  The density profile Eta is too steep and the'              
       write(18,*)'  code can not handle this. Try decreasing the'              
       write(18,*)'  outer radius Y (see Manual, 3.3.3).         '              
       write(18,*)' *********************************************'              
       END IF                                                                   
      END IF                                                                    
      IF (msgno.EQ.17) THEN                                                     
       write(12,*)' *****************  WARNING  ********************'           
       write(12,*)'  Eta is too steep and reaches values less than  '           
       write(12,*)'  1e-12. Try decreasing the outer radius Y.      '           
       write(12,*)'  (see Manual,3.3.3)                             '           
       write(12,*)' ************************************************'           
       IF(iX.GT.0) THEN                                                         
       write(18,*)' *****************  WARNING  ********************'           
       write(18,*)'  Eta is too steep and reaches values less than  '           
       write(18,*)'  1e-12. Try decreasing the outer radius Y.      '           
       write(18,*)'  (see Manual,3.3.3)                             '           
       write(18,*)' ************************************************'           
       END IF                                                                   
      END IF                                                                    
      IF (msgno.EQ.18) THEN                                                     
       write(12,*)' ************  WARNING  ************************ '           
       write(12,*)'  Dynamical range of Eta more than 1.E-12.       '           
       write(12,*)'  The outer radius Y must be decreased so that   '           
       write(12,*)'  Eta does not go below 1.E-12 (see Manual,3.3.3)'           
       write(12,*)' *********************************************** '           
       IF(iX.GT.0) THEN                                                         
       write(18,*)' ************  WARNING  ************************ '           
       write(18,*)'  Dynamical range of Eta more than 1.E-12.       '           
       write(18,*)'  The outer radius Y must be decreased so that   '           
       write(18,*)'  Eta does not go below 1.E-12 (see Manual,3.3.3)'           
       write(18,*)' *********************************************** '           
       END IF                                                                   
      END IF                                                                    
      IF (msgno.EQ.19) THEN                                                     
       write(12,*)' ************ A BIG ERROR!!!************* '                  
       write(12,*)'  Singular matrix in LUDCMP when called   '                  
       write(12,*)'  from ANALINT. Stopping the calculation. '                  
       write(12,*)' **************************************** '                  
      END IF                                                                    
      IF (msgno.EQ.20) THEN                                                     
       write(12,*)' ************ A BIG ERROR!!! ************ '                  
       write(12,*)'  Singular matrix in LUDCMP when called   '                  
       write(12,*)'  from INVERT. Stopping the calculation.  '                  
       write(12,*)' **************************************** '                  
      END IF                                                                    
c -----------------------------------------------------------------------       
101   RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE MULTIPLY(type,np1,nr1,np2,nr2,mat,vec1,omat,flag,q1,q2)        
c =======================================================================       
c This subroutine evaluates the following expression:                           
c [q2] = flag*[q1] + [mat]*[tt*vec1]. Here tt is [omat] for type=1 and          
c 1-[omat] for type=2. mat is matrix of physical size (np2,np1,np1) and         
c real size (nr2,nr1,nr1). omat, vec1, q1 and q2 are matrices of                
c physical size (np2,np1) and real size (nr2,nr1).     [Z.I., Nov. 1995]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      INTEGER type, np1, nr1, np2, nr2, flag, i2, i1, idum                      
      DOUBLE PRECISION mat(np2,np1,np1), vec1(np2,np1), omat(np2,np1),          
     &       aux, q1(np2,np1), q2(np2,np1)                                      
c -----------------------------------------------------------------------       
c     loop over index 2                                                         
      DO i2 = 1, nr2                                                            
c       loop over index 1                                                       
        DO i1 = 1, nr1                                                          
          q2(i2,i1) = flag * q1(i2,i1)                                          
c         loop over dummy index (multiplication)                                
          DO idum = 1, nr1                                                      
            IF (type.EQ.1) THEN                                                 
              aux = omat(i2,idum)                                               
              ELSE                                                              
              aux = 1.0 - omat(i2,idum)                                         
            END IF                                                              
            q2(i2,i1) = q2(i2,i1) + mat(i2,i1,idum)*aux*vec1(i2,idum)           
          END DO                                                                
          IF (q2(i2,i1).LT.dynrange*dynrange) q2(i2,i1) = 0.0                   
        END DO                                                                  
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE MULTIP2(type,np1,nr1,np2,nr2,nr3,np3,mat,vec1,omat,q1)         
c =======================================================================       
c This subroutine evaluates the following expression:                           
c [q1] = [mat]*[tt*vec1] / 4Pi. Here tt is [omat] for type=1 and                
c 1-[omat] for type=2. mat is matrix of physical size (np2,np3,np1) and         
c real size (nr2,nr3,nr1). omat and vec1 are matrices of physical size          
c (np2,np1) and real size (nr2,nr1). q1 is a matrix of physical size            
c (np2,np3) and real size (nr2,nr3)                                             
c 1, 2 and 3 correspond to nY, nL and nP.              [Z.I., Nov. 1995]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER nY, nYprev, nP, nL, nPcav                                         
      DOUBLE PRECISION Y(npY), Yprev(npY), P(npP), lambda(npL),                 
     &       bOut(npP+2)                                                        
      COMMON /grids1/ nY, nYprev, nP, nPcav, nL                                 
      COMMON /grids2/ Y, Yprev, P, lambda, bOut                                 
      INTEGER type, np1, nr1, np2, nr2, np3, nr3, i2, i3, idum                  
      DOUBLE PRECISION mat(np2,np3,np1), vec1(np2,np1), omat(np2,np1),          
     &       aux, q1(np2,np3)                                                   
c -----------------------------------------------------------------------       
c     loop over index 2 (wavelength)                                            
      DO i2 = 1, nr2                                                            
c       loop over index 3 (impact parameter)                                    
        DO i3 = 1, nr3                                                          
          q1(i2,i3) = 0.0                                                       
c         loop over dummy index (multiplication)                                
          DO idum = 1, nr1                                                      
            IF (type.EQ.1) THEN                                                 
              aux = omat(i2,idum)                                               
            ELSE                                                                
              aux = 1.0 - omat(i2,idum)                                         
            END IF                                                              
            q1(i2,i3) = q1(i2,i3) + mat(i2,i3,idum)*aux*vec1(i2,idum)           
          END DO                                                                
        END DO                                                                  
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE MYSPLINE(x,N,alpha,beta,gamma,delta)                           
c =======================================================================       
c This subroutine finds arrays alpha, beta, gamma and delta describing          
c a cubic spline approximation of an unknown function f(x) given as an          
c array f(i)=f(x(i)) with i=1..N. The cubic spline approximation is:            
c f(x)=a(i) + b(i)*t + c(i)*t^2 + d(i)*t^3  for x(i).LE.x.LE.x(i+1)             
c and t = (x-x(i))/(x(i+1)-x(i)), i=1..N-1. Coefficients a,b,c,d are            
c equal to:                                                                     
c a(i) = alpha(i,1)*f(1) + alpha(i,2)*f(2) + ... + alpha(i,N)*f(N)              
c and b,c,d analogously.                               [Z.I., Dec. 1995]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER N, i, j, dummy, Kron                                              
      DOUBLE PRECISION x(npY), alpha(npY,npY), beta(npY,npY),                   
     &       delta(npY,npY), secnder(npY,npY), yaux(npY), deraux(npY),          
     &       y2at1, y2atN, D, gamma(npY,npY)                                    
      EXTERNAL Kron                                                             
c -----------------------------------------------------------------------       
c     generate second derivatives, secnder(j,l)                                 
      DO j = 1, N                                                               
        DO dummy = 1, N                                                         
          IF (dummy.EQ.j) THEN                                                  
            yaux(dummy) = 1.0                                                   
            ELSE                                                                
            yaux(dummy) = 0.0                                                   
          END IF                                                                
        END DO                                                                  
        y2at1 = (yaux(2)-yaux(1))/(x(2)-x(1))                                   
        y2atN = (yaux(N)-yaux(N-1))/(x(N)-x(N-1))                               
        CALL SPLINE(x,yaux,N,y2at1,y2atN,deraux)                                
        DO i = 1, N                                                             
          secnder(i,j) =  deraux(i)                                             
c          secnder(i,j) = 0.0                                                   
        END DO                                                                  
      END DO                                                                    
c     generate alpha, beta, gamma, delta                                        
      DO i = 1, N-1                                                             
        D = (x(i+1) - x(i))*(x(i+1) - x(i)) / 6.0                               
        DO j = 1, N                                                             
          alpha(i,j) = Kron(i,j)*1.0                                            
          beta(i,j) = Kron(i+1,j) - Kron(i,j)                                   
          beta(i,j) = beta(i,j) - D*(2.*secnder(i,j)+secnder(i+1,j))            
          gamma(i,j) = 3. * D * secnder(i,j)                                    
          delta(i,j) = D*(secnder(i+1,j)-secnder(i,j))                          
        END DO                                                                  
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE polint(xa,ya,n,x,y,dy)                                         
c =======================================================================       
      INTEGER n,NMAX                                                            
      DOUBLE PRECISION dy,x,y,xa(n),ya(n)                                       
      PARAMETER (NMAX=10)                                                       
      INTEGER i,m,ns                                                            
      DOUBLE PRECISION den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)                     
c -----------------------------------------------------------------------       
      ns=1                                                                      
      dif=DABS(x-xa(1))                                                         
      do 11 i=1,n                                                               
        dift=DABS(x-xa(i))                                                      
        if (dift.lt.dif) then                                                   
          ns=i                                                                  
          dif=dift                                                              
        endif                                                                   
        c(i)=ya(i)                                                              
        d(i)=ya(i)                                                              
11    continue                                                                  
      y=ya(ns)                                                                  
      ns=ns-1                                                                   
      do 13 m=1,n-1                                                             
        do 12 i=1,n-m                                                           
          ho=xa(i)-x                                                            
          hp=xa(i+m)-x                                                          
          w=c(i+1)-d(i)                                                         
          den=ho-hp                                                             
          if(den.eq.0.)pause 'failure in polint'                                
          den=w/den                                                             
          d(i)=hp*den                                                           
          c(i)=ho*den                                                           
12      continue                                                                
        if (2*ns.lt.n-m)then                                                    
          dy=c(ns+1)                                                            
        else                                                                    
          dy=d(ns)                                                              
          ns=ns-1                                                               
        endif                                                                   
        y=y+dy                                                                  
13    continue                                                                  
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE PowerInt(N,N1,N2,x,y,integral)                                 
c =======================================================================       
c This subroutine calculates integral I(y(x)*dx). Both y and x are              
c 1D arrays, y(i), x(i) with i=1,N (declared with NN). Lower and upper          
c integration limits are x(N1) and x(N2), respectively. The method used         
c is a power-law approximation for y(x) between any two points .                
c                                                      [Z.I., Mar. 1996]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER i, N, N1, N2                                                      
      DOUBLE PRECISION x(N), y(N), integral, pow, C, delint                     
c -----------------------------------------------------------------------       
c     set integral to 0 and accumulate result in the loop                       
      integral = 0.0                                                            
c     calculate weight, wgth, and integrate in the same loop                    
      IF (N2.GT.N1) THEN                                                        
        DO i = N1, N2-1                                                         
          pow = dlog(y(i+1)/y(i)) / dlog(x(i+1)/x(i))                           
          C = y(i) / x(i)**pow                                                  
          delint = (x(i+1)**(pow+1)-x(i)**(pow+1.))*C/(pow+1.)                  
c         add contribution to the integral                                      
          integral = integral + delint                                          
        END DO                                                                  
        ELSE                                                                    
c        integral = 0.0                                                         
        integral = y(1)                                                         
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE DoProduct(NN,Yt,pt,p0,j,Prd)                                     
c =======================================================================       
c This is an auxiliary subroutine which evaluates a messy expression            
c needed to calculate normalization constants for a broken power law            
c density.                                             [Z.I., Aug. 1996]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER NN, i, j                                                          
      DOUBLE PRECISION Yt(NN), pt(NN), Prd, p0                                  
c -----------------------------------------------------------------------       
      Prd = Yt(1)**(pt(1) - p0)                                                 
      IF (j.GT.1) THEN                                                          
        DO i = 2, j                                                             
          Prd = Prd * Yt(i)**(pt(i) - pt(i-1))                                  
        END DO                                                                  
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE ROMBERG2(a,b,ss8)                                              
c =======================================================================       
c This subroutine performs Romberg integration of 8 functions calculated        
c in trapzd2 (by calling subroutine TWOFUN) on interval [a,b].                  
c The results are returned in ss8(1..8). Desired accuracy accRomb is            
c user supplied and comes through COMMON /numerics/ read in from                
c 'numerics.inc'. This subroutine is based on slightly changed versions         
c of 'qromb' and 'qromo' from Numerical Recipes.                                
c                                                    [MN & ZI,Aug'96]           
c =======================================================================       
      INTEGER Ncav, Nins                                                        
      DOUBLE PRECISION accRomb, accuracy, accConv, delTAUsc, facc,              
     &                 dynrange, EtaRat, accFbol                                
      COMMON /numerics/ accRomb, accuracy, accConv, delTAUsc, facc,             
     &                  dynrange, EtaRat, accFbol, Ncav, Nins                   
      INTEGER fconv(8),JMAX,JMAXP,K,KM, J, iC, idone, kaux                      
      PARAMETER (JMAX=50, JMAXP=JMAX+1, K=5, KM=K-1)                            
      DOUBLE PRECISION ss, ss8(8), S2D(8,JMAXP), h(JMAXP), sjKM(JMAXP),         
     &                 a, b, EPS, h0, dss, s8(8), chk(8)                        
c -----------------------------------------------------------------------       
      EPS = accRomb                                                             
      h0 = 0.0                                                                  
      h(1)=1.0                                                                  
c     intialize convergence flags                                               
      DO iC = 1, 8                                                              
         fconv(iC) = 0                                                          
      END DO                                                                    
c     integrate until all 8 intergrals converge                                 
      idone = 0                                                                 
      j = 0                                                                     
      DO WHILE(idone.NE.1.and.j.LE.JMAX)                                        
        j = j + 1                                                               
c       integrate with j division points                                        
        call trapzd2(a,b,s8,j)                                                  
        DO iC = 1, 8                                                            
           S2D(iC,j) = S8(iC)                                                   
        END DO                                                                  
c       check if any of 8 integrals has converged                               
        IF (j.ge.K) THEN                                                        
           idone = 1                                                            
           DO iC = 1, 8                                                         
             IF (fconv(iC).EQ.0) THEN                                           
c              generate array for polint                                        
               DO kaux = 1, j                                                   
                 sjKM(kaux) = S2D(iC,kaux)                                      
               END DO                                                           
c              predict the integral for stepsize h->h0=0.0                      
               CALL polint(h(j-KM),sjKM(j-KM),K,h0,ss,dss)                      
               IF (dabs(dss).le.EPS*dabs(ss)) THEN                              
                 SS8(iC) = ss                                                   
                 fconv(iC) = 1                                                  
               ELSE                                                             
                 chk(iC) = dabs(dss)/dabs(ss)                                   
               END IF                                                           
             END IF                                                             
             idone = idone*fconv(iC)                                            
           END DO                                                               
        END IF                                                                  
        h(j+1)=0.25*h(j)                                                        
      END DO                                                                    
      IF (j.GE.jMAX) THEN                                                       
        write(*,*)' Reached the limiting number of steps in ROMBERG2'           
        write(*,*)'You might want to change accRomb in the input file'          
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE ROMBY(fnc,a,b,ss)                                              
c =======================================================================       
c This subroutine performs Romberg integration of function func on              
c interval [a,b]. The result is returned in ss. Desired accuracy is set         
c to 0.002.                                            [Z.I., Feb. 1996]        
c =======================================================================       
      INTEGER JMAX,JMAXP,K,KM, J                                                
      PARAMETER (JMAX=30, JMAXP=JMAX+1, K=3, KM=K-1)                            
      DOUBLE PRECISION a,b,fnc,ss,EPS, aux, dss,h(JMAXP),s(JMAXP)               
      EXTERNAL fnc                                                              
c -----------------------------------------------------------------------       
      EPS = 0.002                                                               
      h(1)=1.                                                                   
      do 11 j=1,JMAX                                                            
        call trapzd(fnc,a,b,s(j),j)                                             
        if (j.ge.K) then                                                        
          aux = 0.0                                                             
          call polint(h(j-KM),s(j-KM),K,aux,ss,dss)                             
          IF (dabs(dss).le.EPS*dabs(ss)) RETURN                                 
        endif                                                                   
        s(j+1)=s(j)                                                             
        h(j+1)=0.25*h(j)                                                        
11    continue                                                                  
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE ScaleTo1(Nmax,N,Y)                                             
c =======================================================================       
c This subroutine scales vector Y such that Y(1) = 1.0                          
c                                                      [Z.I., Jan. 1997]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER Nmax, N, i                                                        
      DOUBLE PRECISION Y(Nmax), Scale                                           
c -----------------------------------------------------------------------       
      Scale = Y(1)                                                              
      DO i = 1, N                                                               
        Y(i) = Y(i) / Scale                                                     
      END DO                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE SHIFT(X,Nmax,N,Xins,i)                                         
c =======================================================================       
c Rearranges a vector X by inserting a new element Xins.    [MN, Aug'96]        
c =======================================================================       
      implicit none                                                             
      integer Nmax, N, i,j                                                      
      DOUBLE PRECISION X(Nmax),Xins                                             
c -----------------------------------------------------------------------       
      DO j = N+1, i+2, -1                                                       
        x(j) = x(j-1)                                                           
      END DO                                                                    
      x(i+1) = xins                                                             
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE SIMPSON(N,N1,N2,x,y,integral)                                  
c =======================================================================       
c This subroutine calculates integral I(y(x)*dx). Both y and x are              
c 1D arrays, y(i), x(i) with i=1,N (declared with NN). Lower and upper          
c integration limits are x(N1) and x(N2), respectively. The method used         
c is Simpson (trapezoid) approximation. The resulting integral is sum of        
c y(i)*wgth, i=N1,N2.                                  [Z.I., Mar. 1996]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER i, N, N1, N2                                                      
      DOUBLE PRECISION x(N), y(N), wgth, integral                               
c -----------------------------------------------------------------------       
c     set integral to 0 and accumulate result in the loop                       
      integral = 0.0                                                            
c     calculate weight, wgth, and integrate in the same loop                    
      IF (N2.GT.N1) THEN                                                        
        DO i = N1, N2                                                           
c         weigths                                                               
          IF (i.NE.N1.AND.i.NE.N2) THEN                                         
            wgth = 0.5 * (x(i+1)-x(i-1))                                        
          ELSE                                                                  
            IF (i.eq.N1) wgth = 0.5 * (x(N1+1)-x(N1))                           
            IF (i.eq.N2) wgth = 0.5 * (x(N2)-x(N2-1))                           
          END IF                                                                
c         add contribution to the integral                                      
          integral = integral + y(i) * wgth                                     
        END DO                                                                  
      ELSE                                                                      
        integral = 0.0                                                          
      END IF                                                                    
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE SORT(RA,N)                                                     
c =======================================================================       
      INTEGER N                                                                 
      DOUBLE PRECISION RA(N)                                                    
c -----------------------------------------------------------------------       
      L=N/2+1                                                                   
      IR=N                                                                      
10    CONTINUE                                                                  
        IF(L.GT.1)THEN                                                          
          L=L-1                                                                 
          RRA=RA(L)                                                             
        ELSE                                                                    
          RRA=RA(IR)                                                            
          RA(IR)=RA(1)                                                          
          IR=IR-1                                                               
          IF(IR.EQ.1)THEN                                                       
            RA(1)=RRA                                                           
            RETURN                                                              
          ENDIF                                                                 
        ENDIF                                                                   
        I=L                                                                     
        J=L+L                                                                   
20      IF(J.LE.IR)THEN                                                         
          IF(J.LT.IR)THEN                                                       
            IF(RA(J).LT.RA(J+1))J=J+1                                           
          ENDIF                                                                 
          IF(RRA.LT.RA(J))THEN                                                  
            RA(I)=RA(J)                                                         
            I=J                                                                 
            J=J+J                                                               
          ELSE                                                                  
            J=IR+1                                                              
          ENDIF                                                                 
        GO TO 20                                                                
        ENDIF                                                                   
        RA(I)=RRA                                                               
      GO TO 10                                                                  
c -----------------------------------------------------------------------       
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE Spline(x,y,n,yp1,ypn,y2)                                       
c =======================================================================       
      INTEGER n,NMAX                                                            
      DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)                                  
      PARAMETER (NMAX=500)                                                      
      INTEGER i,k                                                               
      DOUBLE PRECISION p,qn,sig,un,u(NMAX)                                      
c -----------------------------------------------------------------------       
      if (yp1.gt..99e30) then                                                   
        y2(1)=0.                                                                
        u(1)=0.                                                                 
      else                                                                      
        y2(1)=-0.5                                                              
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)                     
      endif                                                                     
      do 11 i=2,n-1                                                             
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))                                       
        p=sig*y2(i-1)+2.                                                        
        y2(i)=(sig-1.)/p                                                        
        u(i)=(6.*((y(i+1)-y(i))/(x(i+                                           
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*                
     *u(i-1))/p                                                                 
11    continue                                                                  
      if (ypn.gt..99e30) then                                                   
        qn=0.                                                                   
        un=0.                                                                   
      else                                                                      
        qn=0.5                                                                  
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))                 
      endif                                                                     
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)                                      
      do 12 k=n-1,1,-1                                                          
        y2(k)=y2(k)*y2(k+1)+u(k)                                                
12    continue                                                                  
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
c ***********************************************************************       
      SUBROUTINE SPLINE2(x,fun,N,coef)                                          
c =======================================================================       
c This subroutine finds coefficients coef(i,j) such that                        
c fun(x)=coef(i,1) + coef(i,2)*x + coef(i,3)*x^2 + coef(i,4)*x^3                
c for x(i).LE.x.LE.x(i+1) is a cubic spline approximation of fun(x),            
c with i=1..N.                                         [Z.I., Feb. 1995]        
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER N, i                                                              
      DOUBLE PRECISION x(npY), coef(npY,4), secnder(npY), y2at1, y2atN,         
     &       Dd, xL, xR, dR, dL, fun(npY), fL, fR                               
c -----------------------------------------------------------------------       
c     find second derivative, secnder                                           
        y2at1 = (fun(2)-fun(1))/(x(2)-x(1))                                     
        y2atN = (fun(N)-fun(N-1))/(x(N)-x(N-1))                                 
        CALL SPLINE(x,fun,N,y2at1,y2atN,secnder)                                
c     generate coef(i,j), j=1,2,3,4                                             
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
c -----------------------------------------------------------------------       
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE splint(xa,ya,y2a,n,x,y)                                        
c =======================================================================       
      INTEGER n                                                                 
      DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)                                   
      INTEGER k,khi,klo                                                         
      DOUBLE PRECISION a,b,h                                                    
c -------------------------------------------------------------------------     
      klo=1                                                                     
      khi=n                                                                     
1     if (khi-klo.gt.1) then                                                    
        k=(khi+klo)/2                                                           
        if(xa(k).gt.x)then                                                      
          khi=k                                                                 
        else                                                                    
          klo=k                                                                 
        endif                                                                   
      goto 1                                                                    
      endif                                                                     
      h=xa(khi)-xa(klo)                                                         
      if (h.eq.0.) pause 'bad xa input in splint'                               
      a=(xa(khi)-x)/h                                                           
      b=(x-xa(klo))/h                                                           
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**          
     *2)/6.                                                                     
c -------------------------------------------------------------------------     
      return                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE trapzd(func,a,b,s,n)                                           
c =======================================================================       
      INTEGER n                                                                 
      DOUBLE PRECISION a,b,s,func                                               
      EXTERNAL func                                                             
      INTEGER it,j                                                              
      DOUBLE PRECISION del,sum,tnm,x                                            
c -------------------------------------------------------------------------     
      IF (n.eq.1) THEN                                                          
        s=0.5*(b-a)*(func(a)+func(b))                                           
      ELSE                                                                      
        it=2**(n-2)                                                             
        tnm=it                                                                  
        del=(b-a)/tnm                                                           
        x=a+0.5*del                                                             
        sum=0.                                                                  
        DO j = 1, it                                                            
          sum=sum+func(x)                                                       
          x=x+del                                                               
        END DO                                                                  
        s=0.5*(s+(b-a)*sum/tnm)                                                 
      END IF                                                                    
c -------------------------------------------------------------------------     
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE WriteOut(is,nG,nameQ,nameNK)                                   
c =======================================================================       
      IMPLICIT none                                                             
      INTEGER npY, npP, npL, npG                                                
      INCLUDE 'userpar.inc'                                                     
      PARAMETER (npG=1)                                                         
      INTEGER denstyp, Ntr, iterETA, Eta7OK, nYEta7                             
      CHARACTER*235 nameETA                                                     
      LOGICAL RDW                                                               
      DOUBLE PRECISION pow, Yout, ETAcoef(npY,4), ptr(10), Ytr(10),             
     &       ETAdiscr(npY), yEta7(1000), Eta7(1000)                             
      COMMON /density1/ denstyp, Ntr, iterETA, nYEta7, Eta7OK, nameETA          
      COMMON /density2/ pow, Yout, ETAcoef, Ytr, ptr, ETAdiscr, Eta7,           
     &                  yEta7, RDW                                              
      INTEGER iLfid, szds, top, Nfiles                                          
      DOUBLE PRECISION TAUtot(npL), SigmaA(npG,npL), SigmaS(npG,npL),           
     &     omega(npL,npY), Tsub(npG), abund(npG,npY), TAUmax,                   
     &     xC(10), xCuser(10), SigExfid, TAUfid, lamfid, qsd,                   
     &     a1, a2, aveV                                                         
      COMMON /optprop/ TAUtot, SigmaA, SigmaS, omega, Tsub, abund,              
     &            TAUmax, xC, xCuser, SigExfid, TAUfid, lamfid, qsd,            
     &            a1, a2, aveV, iLfid, szds, top, Nfiles                        
      INTEGER startyp(2), Nlamtr(2), nBB(2)                                     
      CHARACTER nameStar(2)*235                                                 
      DOUBLE PRECISION Tstar, lamtr(2,101), klam(2,100), Tbb(2,10),             
     &                 rellum(2,10), mu1, ksi, mu2, xSiO, r1rs                  
      COMMON /source/ Tstar, lamtr, klam, Tbb, rellum, mu1, ksi, mu2,           
     &                xSiO, r1rs, startyp, Nlamtr, nBB, nameStar                
      INTEGER is, iG, nG, i, length                                             
      CHARACTER*72 strpow, aux, src, chaux*3                                    
      CHARACTER*(*) nameQ(npG), nameNK(10)                                      
c -------------------------------------------------------------------------     
      IF (denstyp.eq.0) THEN                                                    
       IF (is.eq.1) THEN                                                        
        src = 'Left-side source spectrum described by'                          
       ELSE                                                                     
        src = 'Right-side source spectrum described by'                         
       END IF                                                                   
      ELSE                                                                      
        src = 'Central source spectrum described by'                            
      END IF                                                                    
      CALL Clean(src, aux, length)                                              
                                                                                
c     #1: black body(ies) for startyp=1                                         
      IF (startyp(is).EQ.1) THEN                                                
       IF (nBB(is).GT.1) THEN                                                   
         CALL ATTACH(aux, length, ' ', src)                                     
c        multiple black bodies                                                  
         write(12,'(a2,a37,i2,a13)')'  ', src, nBB(is),' black bodies'          
         write(12,'(a27)')' with temperatures (in K):'                          
         write(12,'(2x,1p,10e10.3)')(Tbb(is,i),i=1,nBB(is))                     
         write(12,'(a42)')' and relative luminosities, respectively:'           
         write(12,'(1p,10e10.1)')(rellum(is,i),i=1,nBB(is))                     
       ELSE                                                                     
c       for a single black body:                                                
        CALL ATTACH(aux,length,' a black body',src)                             
         write(12,'(a2,a)') '  ',src                                            
         IF (Tstar.LT.9999.999) THEN                                            
           CALL getfs(Tstar,0,1,strpow)                                         
           write(12,'(a19,a5,a2)')' with temperature:',strpow,' K'              
         ELSE                                                                   
           CALL getfs(Tstar,0,1,strpow)                                         
           write(12,'(a19,a6,a2)')' with temperature:',strpow,' K'              
         END IF                                                                 
       END IF                                                                   
      END IF                                                                    
                                                                                
c     #2: Engelke-Marengo function for startyp=2                                
      IF (startyp(is).EQ.2) THEN                                                
         CALL ATTACH(aux, length,' Engelke-Marengo function', src)              
         write(12,'(a2,a)') '  ',src                                            
         CALL getfs(Tbb(is,1),0,1,strpow)                                       
         write(12,'(a13,a5,a16)')' with Teff =',strpow,' K and depth of'        
         write(12,'(a30,F6.1,a2)')' the SiO absorption feature =',              
     &                              xSiO,' %'                                   
      END IF                                                                    
                                                                                
c     #3: power-law(s) for startyp=3                                            
      IF (startyp(is).EQ.3) THEN                                                
       IF (Nlamtr(is).GT.0) THEN                                                
         CALL ATTACH(aux,length,' power law:',src)                              
         write(12,'(a2,a)') '  ',src                                            
         write(12,*)'    lambda      k'                                         
         DO i = 1, Nlamtr(is)                                                   
           write(12,'(1x,1p,e10.3)')lamtr(is,i)                                 
           write(12,'(11x,1p,e10.3)')klam(is,i)                                 
         END DO                                                                 
         write(12,'(1x,1p,e10.3)')lamtr(is,Nlamtr(is)+1)                        
       ELSE                                                                     
         write(12,*)                                                            
     &      ' Input data for the source spectrum is not good.'                  
         write(12,*)' Changed to a 10000 K black body'                          
       END IF                                                                   
      END IF                                                                    
                                                                                
c     spectrum from a file for startyp=4,5,6                                    
      IF (startyp(is).GE.4.AND.startyp(is).LE.6) THEN                           
        write(12,*)' Stellar spectrum supplied from file:'                      
        write(12,'(a2,a70)') '  ',nameStar(is)                                  
      END IF                                                                    
      IF(is.eq.1)                                                               
     &  write(12,*)' --------------------------------------------'              
                                                                                
      IF(is.eq.1) THEN                                                          
c      2) DUST PROPERTIES                                                       
c      2.1 Chemical Composition                                                 
        write(12,*)' Abundances for supported grains:'                          
        write(12,*)' Sil-Ow Sil-Oc Sil-DL grf-DL amC-Hn SiC-Pg'                 
        write(12,'(6f7.3)')(xC(i),i=1,3),xC(4)+xC(5),(xC(i),i=6,7)              
        IF (top.EQ.2) THEN                                                      
          write(12,*)' Abundances for user supplied grains:'                    
          write(12,'(i6,9i7)')(i,i=1,Nfiles)                                    
          write(12,'(10f7.3)')(xCuser(i),i=1,Nfiles)                            
          write(12,*)' User supplied n and k from:'                             
          DO i = 1, Nfiles                                                      
            write(12,'(a2,i1,a2,a70)')'  ',i,') ',nameNK(i)                     
          END DO                                                                
        END IF                                                                  
c      user supplied cross-sections:                                            
       IF (top.EQ.3) THEN                                                       
        DO iG = 1, nG                                                           
          write(12,*)' Optical properties from file:'                           
          write(12,'(a2,a70)')'  ',nameQ(iG)                                    
        END DO                                                                  
       END IF                                                                   
c      2.2 Grain size distribution                                              
       IF (top.NE.3) THEN                                                       
         IF (szds.EQ.3) THEN                                                    
          chaux = 'KMH'                                                         
         ELSE                                                                   
          chaux = 'MRN'                                                         
         END IF                                                                 
         write(12,'(a2,a3,a19)')'  ',chaux,'size distribution:'                 
         CALL getfs(qsd,1,0,strpow)                                             
         write(12,'(a15,a5)')'      Power q:',strpow                            
         write(12,'(a15,1p,e9.2,a8)')                                           
     &                         ' Minimal size:',a1,' microns'                   
         IF (szds.EQ.3) THEN                                                    
            write(12,'(a22,1p,e9.2,a8)')                                        
     &                            ' Characteristic size:',a2,' microns'         
         ELSE                                                                   
           write(12,'(a15,1p,e9.2,a8)')' Maximal size:',a2,' microns'           
         END IF                                                                 
       END IF                                                                   
       write(12,*)' --------------------------------------------'               
c      2.3 Dust temperature on inner boundary                                   
       DO iG = 1, nG                                                            
        CALL getfs(Tsub(iG),0,1,strpow)                                         
        IF (denstyp.eq.0) THEN                                                  
          write(12,'(a45,a5,a2)')                                               
     &      ' Dust temperature on the slab left boundary:', strpow,' K'         
        ELSE                                                                    
          CALL getfs(Tsub(iG),0,1,strpow)                                       
          write(12,'(a41,a5,a2)')                                               
     &      ' Dust temperature on the inner boundary:', strpow,' K'             
        END IF                                                                  
       END DO                                                                   
       write(12,*)' --------------------------------------------'               
      END IF                                                                    
c -------------------------------------------------------------------------     
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      SUBROUTINE zbrac(func,x1,x2,Ntry,succes)                                  
c =======================================================================       
      INTEGER NTRY, succes                                                      
      DOUBLE PRECISION x1,x2,func,FACTOR                                        
      EXTERNAL func                                                             
      PARAMETER (FACTOR=1.6)                                                    
      INTEGER j                                                                 
      DOUBLE PRECISION f1,f2                                                    
c -------------------------------------------------------------------------     
      IF(x1.eq.x2) PAUSE 'you have to guess an initial range in zbrac'          
      f1=func(x1)                                                               
      f2=func(x2)                                                               
      succes=1                                                                  
      DO j = 1, NTRY                                                            
        IF(f1*f2.lt.0.) RETURN                                                  
        IF(DABS(f1).lt.DABS(f2)) THEN                                           
          x1=x1+FACTOR*(x1-x2)                                                  
c         f1=func(x1)                                                           
c     IF's added to prevent breaking of slab case on DEC [MN,Jun'99]            
          IF (x1.LT.(-100.)) THEN                                               
            f1 = 5.0e5                                                          
          ELSE                                                                  
            f1=func(x1)                                                         
          END IF                                                                
        ELSE                                                                    
          x2=x2+FACTOR*(x2-x1)                                                  
          f2=func(x2)                                                           
        END IF                                                                  
      END DO                                                                    
      succes=0                                                                  
c -------------------------------------------------------------------------     
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
c ***********************************************************************       
      DOUBLE PRECISION FUNCTION zriddr(func,x1,x2,MAXIT,xacc)                   
c =======================================================================       
      INTEGER MAXIT                                                             
      DOUBLE PRECISION x1,x2,xacc,func,UNUSED                                   
      PARAMETER (UNUSED=-1.11E30)                                               
      EXTERNAL func                                                             
      INTEGER j                                                                 
c     aux will be used as argument for SIGN function few lines below.           
c     It needs to be REAL to conform to FORTRAN90                               
      REAL aux                                                                  
      DOUBLE PRECISION fh,fl,fm,fnew,s,xh,xl,xm,xnew                            
c -------------------------------------------------------------------------     
c      fl=func(x1)                                                              
c     IF's added to prevent breaking of slab case on DEC [MN,Jun'99]            
      IF (x1.LT.(-100.)) THEN                                                   
        fl = 5.0e5                                                              
      ELSE                                                                      
        fl=func(x1)                                                             
      END IF                                                                    
      fh=func(x2)                                                               
      if((fl.gt.0..and.fh.lt.0.).or.(fl.lt.0..and.fh.gt.0.))then                
        xl=x1                                                                   
        xh=x2                                                                   
        zriddr=UNUSED                                                           
        do 11 j=1,MAXIT                                                         
          xm=0.5*(xl+xh)                                                        
          fm=func(xm)                                                           
          s=sqrt(fm**2.-fl*fh)                                                  
          if(s.eq.0.)return                                                     
          aux = fl-fh                                                           
          xnew=xm+(xm-xl)*(sign(1.,aux)*fm/s)                                   
          if (DABS(xnew-zriddr).le.xacc) return                                 
          zriddr=xnew                                                           
          fnew=func(zriddr)                                                     
          if (fnew.eq.0.) return                                                
          if(dsign(fm,fnew).ne.fm) then                                         
            xl=xm                                                               
            fl=fm                                                               
            xh=zriddr                                                           
            fh=fnew                                                             
          else if(dsign(fl,fnew).ne.fl) then                                    
            xh=zriddr                                                           
            fh=fnew                                                             
          else if(dsign(fh,fnew).ne.fh) then                                    
            xl=zriddr                                                           
            fl=fnew                                                             
          else                                                                  
            pause 'never get here in zriddr'                                    
          endif                                                                 
          if(dabs(xh-xl).le.xacc) return                                        
11      continue                                                                
        pause 'zriddr exceed maximum iterations'                                
      else if (fl.eq.0.) then                                                   
        zriddr=x1                                                               
      else if (fh.eq.0.) then                                                   
        zriddr=x2                                                               
      else                                                                      
        pause 'root must be bracketed in zriddr'                                
      endif                                                                     
c -------------------------------------------------------------------------     
      RETURN                                                                    
      END                                                                       
c ***********************************************************************       
                                                                                
                                                                                
                                                                                
                                                                                
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
c                             THE END                                           
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
                                                                                
                                                                                
                                                                                
                                                                                
