F90=gfortran -O3 -fbounds-check
#F90=gfortran -pg -fbounds-check
#F90=/opt/ekopath-4.0.10/bin/pathf90 -O3 -openmp
#F90=/opt/ekopath-4.0.10/bin/pathf90 -O3

all:
	$(F90) -o dusty common3.10.f90 Dusty3.10.f90 

openmp:
	$(F90) -fopenmp -o dusty common3.10.f90 Dusty3.10.f90
