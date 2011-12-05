#F90=gfortran -fbounds-check
F90=gfortran -O3 -fbounds-check -lgomp
#F90=gfortran -pg -fbounds-check
#F90=/opt/ekopath-4.0.10/bin/pathf90 -O3 -openmp
#F90=/opt/ekopath-4.0.10/bin/pathf90 -O3

all:
	$(F90) -o dusty common3.10.f90 Dusty3.10.f90 source/input.f90 source/math.f90
openmp:
	$(F90) -fopenmp -o dusty common3.10.f90 Dusty3.10.f90
new:
	$(F90) -o dusty source_new/common.f90 source_new/dusty.f90 source_new/misc.f90 source_new/inout.f90 source_new/math.f90 source_new/optprop.f90 source_new/kernel.f90 source_new/kernel_matrix.f90 source_new/winds.f90

new_openmp:
	$(F90) -fopenmp -o dusty source_new/common.f90 source_new/dusty.f90 source_new/misc.f90 source_new/inout.f90 source_new/math.f90 source_new/optprop.f90 source_new/kernel.f90 source_new/kernel_matrix.f90 source_new/winds.f90