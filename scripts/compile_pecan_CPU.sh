#!/bin/sh

module load gcc/8.2.0 openmpi/4.0.1
export PATH=/hrtc/apps/devtools/cmake/x86_64/3.23.2/bin/:$PATH


export SEM_TPL_ROOT_DIR=/data/gpfs/Users/j0436735/travis-deployments/CPU/GEOSX_TPL-220-932-6d57932/
export RAJA_DIR=/data/gpfs/Users/j0436735/travis-deployments/CPU/GEOSX_TPL-220-932-6d57932/raja/lib/cmake/raja
export CAMP_DIR=/data/gpfs/Users/j0535952/work2023/GEOSX_2023/GEOSX_TPL/TPL_20230507/thirdPartyLibs/install-pecan-CPU-release/raja/lib/cmake/camp
export CHAI_DIR=/data/gpfs/Users/j0436735/travis-deployments/CPU/GEOSX_TPL-220-932-6d57932/chai/share/chai/cmake/
export CHAI_INCLUDE_DIR=/data/gpfs/Users/j0436735/travis-deployments/CPU/GEOSX_TPL-220-932-6d57932/chai/include/
export UMPIRE_DIR=/data/gpfs/Users/j0436735/travis-deployments/CPU/GEOSX_TPL-220-932-6d57932/chai/lib/cmake/umpire/
export LVARRAY_DIR=/data/gpfs/Users/j0436735/travis-deployments/CPU/GEOSX-c149a2e/share/lvarray/cmake/
export LVARRAY_INCLUDE_DIR=/data/gpfs/Users/j0436735/travis-deployments/CPU/GEOSX-c149a2e/include/
export CALIPER_DIR=/data/gpfs/Users/j0535952/work2023/SEMCode_2023/TPL/Caliper/install/share/cmake/caliper
export KOKKOS_DIR=/data/gpfs/Users/j0535952/work2023/SEMCode_2023/TPL/kokkos/install-openmp/lib64/cmake/Kokkos/
export LD_LIBRARY_PATH=/data/gpfs/Users/j0535952/work2023/SEMCode_2023/TPL/Caliper/install/lib64/:$LD_LIBRARY_PATH
