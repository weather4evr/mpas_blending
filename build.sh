#! /usr/bin/env bash
#
# Author: Larissa Reames CIWRO/NOAA/NSSL/FRDD

set -eux

target=${1:-"NULL"}
compiler=${compiler:-"intel"}
echo $target, $compiler
if [[ "$target" == "linux.*" || "$target" == "macosx.*" ]]; then
    unset -f module
    set +x
    source ./modulefiles/build.$target > /dev/null
    set -x
else
    set +x
   #source ./machine-setup.sh
   #module use ./modulefiles
   #module load build.$target.$compiler.lua
   #module list
    set -x
fi

target="derecho"

if [[ "$target" == "cheyenne" ]] ; then # $target set in ./machine-setup.sh
   module purge
   module restore modules_intel19 # from user schwartz
fi
if [[ "$target" == "derecho" ]] ; then
   module --force purge
#  module restore default # default is from user schwartz
#  module swap hdf5/1.12.2 hdf5-mpi/1.12.2
#  module swap netcdf      netcdf-mpi
   module load ncarenv/23.09 craype/2.7.23 intel-oneapi/2023.2.1 ncarcompilers/1.0.0
   module load cray-mpich/8.1.27 hdf5-mpi/1.12.2 netcdf-mpi/4.9.2 esmf/8.6.0 cmake/3.26.3
   module load grib-util/1.2.4 wgrib2/3.1.1 mkl/2023.2.0
   module load parallel-netcdf/1.12.3
   module load parallelio/2.6.2
   module list
fi

if [[ "$target" == "hera" || "$target" == "orion" || "$target" == "wcoss2" ]]; then
  #CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=../ -DEMC_EXEC_DIR=ON -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Debug"
   CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=../ -DEMC_EXEC_DIR=ON -DBUILD_TESTING=OFF"
   #CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=../ -DEMC_EXEC_DIR=ON -DENABLE_DOCS=ON -DBUILD_TESTING=ON"
else
  #CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=../ -DEMC_EXEC_DIR=ON -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Debug"
  CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=../ -DEMC_EXEC_DIR=ON -DBUILD_TESTING=OFF"
fi

rm -fr ./build
mkdir ./build && cd ./build

cmake .. ${CMAKE_FLAGS}

make -j 8 VERBOSE=1
make install

#make test
#ctest -I 4,5

exit
