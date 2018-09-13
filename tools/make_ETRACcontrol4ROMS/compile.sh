#!/bin/csh

setenv EXECUTABLE $1
setenv FORTRAN_COMPILER ifort
setenv NC_CONFIG /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/intel2016.4/netcdf-fortran/4.4.4/bin/nc-config

module load netcdf-fortran/4.4.4

make clean
make

exit
