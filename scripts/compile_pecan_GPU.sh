#!/bin/sh

module load gcc/8.2.0 openmpi/4.0.1 cuda
export PATH=/hrtc/apps/devtools/cmake/x86_64/3.23.2/bin/:$PATH

export GEOSX_TPL_DIR=/data/gpfs/Users/j0436735/travis-deployments/GPU/GEOSX_TPL-224-965-1650e48
export GEOS_DIR=/data/gpfs/Users/j0436735/travis-deployments/GPU/GEOSX-7d9c0a7
export SEM_TPL_ROOT_DIR=${GEOSX_TPL_DIR}

#export CAMP_DIR=${SEM_TPL_ROOT_DIR}/raja/lib/cmake/camp
export CAMP_DIR=/data/gpfs/Users/j0535952/work2023/SEMCode_2023/TPL/camp/install-cuda/lib/cmake/camp
export adiak_DIR=${SEM_TPL_ROOT_DIR}/adiak/lib/cmake/adiak/

export RAJA_DIR=${SEM_TPL_ROOT_DIR}/raja/lib/cmake/raja
export CHAI_DIR=${SEM_TPL_ROOT_DIR}/chai/share/chai/cmake/
export CHAI_INCLUDE_DIR=${SEM_TPL_ROOT_DIR}/chai/include/
export UMPIRE_DIR=${SEM_TPL_ROOT_DIR}/chai/lib/cmake/umpire/
export LVARRAY_DIR=${GEOS_DIR}/share/lvarray/cmake/
export LVARRAY_INCLUDE_DIR=${GEOS_DIR}/include/
export CALIPER_DIR=${SEM_TPL_ROOT_DIR}/caliper/share/cmake/caliper
export caliper_INCLUDE_DIR=${SEM_TPL_ROOT_DIR}/caliper/include
export LD_LIBRARY_PATH=${SEM_TPL_ROOT_DIR}/caliper/lib64/:$LD_LIBRARY_PATH

#export KOKKOS_DIR=/data/gpfs/Users/j0535952/work2023/SEMCode_2023/TPL/kokkos/install-openmp/lib64/cmake/Kokkos/
