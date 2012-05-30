F90=gfortran -O3 

openmp:
	$(F90) -lgomp -fopenmp -o dusty Dusty4.f90

single:
	$(F90) -o dusty Dusty4.f90

