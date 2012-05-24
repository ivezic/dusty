F90=gfortran -O3 -fbounds-check -lgomp

openmp:
	$(F90) -fopenmp -o dusty Dusty4.f90
#source/interfaces.f90  source/common.f90 source/dusty.f90 source/rdinp.f90 source/math.f90 source/misc.f90 source/inout.f90 source/msg.f90 source/optprop.f90 source/kernel.f90 source/solve_matrix.f90 # source/winds.f90

single:
#	$(F90) -o dusty Dusty4.f90
#source/interfaces.f90 source/common.f90 source/nonopenmp.f90 source/dusty.f90 source/rdinp.f90 source/math.f90 source/misc.f90 source/inout.f90 source/msg.f90 source/optprop.f90 source/kernel.f90 source/solve_matrix.f90 # source/winds.f90

