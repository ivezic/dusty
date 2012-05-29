#!/bin/bash

rm release/* -rf
mkdir release/dusty
echo '!DUSTY RELEASE' > release/dusty/dusty.f90
cat source/common.f90 >> release/dusty/dusty.f90
cat source/dusty.f90 >> release/dusty/dusty.f90
cat source/inout.f90 >> release/dusty/dusty.f90
cat source/kernel.f90 >> release/dusty/dusty.f90
cat source/math.f90 >> release/dusty/dusty.f90
cat source/misc.f90 >> release/dusty/dusty.f90
cat source/msg.f90 >> release/dusty/dusty.f90
cat source/nonopenmp.f90 >> release/dusty/dusty.f90
cat source/optprop.f90 >> release/dusty/dusty.f90
cat source/rdinp.f90 >> release/dusty/dusty.f90
cat source/solve_matrix.f90 >> release/dusty/dusty.f90
cat source/winds.f90 >> release/dusty/dusty.f90 

mkdir release/dusty/docs
cp data release/dusty/ -rf
cp examples release/dusty/ -rf

echo 'all:' > release/dusty/Makefile
echo '\t gfortran -O3 -lgmp -fopenmp -o dusty dusty.f90' >> release/dusty/Makefile
echo 'gfortran -O3 -lgomp -fopenmp -o dusty.exe dusty.f90' > release/dusty/compile.bat

tar -cf release/dusty.tar release/dusty
