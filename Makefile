all:
	gfortran -O3 -o dusty common3.10.f90 Dusty3.10.f90
#	gfortran -g -o dusty common3.10.f90 Dusty3.10.f90

openmp:
	gfortran -O3 -fopenmp -o dusty common3.10.f90 Dusty3.10.f90