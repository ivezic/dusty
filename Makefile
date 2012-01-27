#F90=gfortran -fbounds-check
#F90=gfortran -O3 -fbounds-check -lgomp
F90=gfortran -pg -fbounds-check -lgomp
#F90=/opt/ekopath-4.0.10/bin/pathf90 -O3 -openmp
#F90=/opt/ekopath-4.0.10/bin/pathf90 -O3

#all:
#	$(F90) -o dusty common3.10.f90 Dusty3.10.f90 source_old/input.f90 source_old/math.f90
#openmp:
#	$(F90) -fopenmp -o dusty common3.10.f90 Dusty3.10.f90
new:
	$(F90) -o dusty source/common.f90 source/dusty.f90 source/rdinp.f90 source/math.f90 source/misc.f90 source/inout.f90 source/msg.f90 source/optprop.f90 source/kernel.f90 #source/kernel_matrix.f90 source/winds.f90


#new_openmp:
#	$(F90) -fopenmp -o dusty source/common.f90 source/dusty.f90 source/misc.f90 source/inout.f90 source/math.f90 source/optprop.f90 source/kernel.f90 source/kernel_matrix.f90 source/winds.f90