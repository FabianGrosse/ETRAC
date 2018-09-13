#!/bin/csh

setenv EXECUTABLE ETRAC
setenv FORTRAN_COMPILER ifort
setenv BUILD_DIR ../Build
setenv SOURCE_DIR ./src

setenv NC_CONFIG /cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/intel2016.4/netcdf-fortran/4.4.4/bin/nc-config

module load netcdf-fortran/4.4.4

set RUN_DIR=$PWD
mkdir -p $BUILD_DIR

cd $SOURCE_DIR
echo $PWD

make clean
make

cd $RUN_DIR
cp $SOURCE_DIR/$EXECUTABLE $RUN_DIR/$BUILD_DIR/

exit
