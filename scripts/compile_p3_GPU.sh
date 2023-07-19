#!/bin/sh

export MODULEPATH=/data_local/appli_local/MTS/GEOSX/modulefiles/:$MODULEPATH
module load tpl/gcc8.4.1-ompi4.1.2-stable.1-release
module load cmake/3.21.4 gcc/8.4.1 cuda/11.0.3 ompi/4.1.2 openblas/0.3.18 python4geosx/3.8.5-gcc-8.4.1-ompi-4.1.2

export SEM_TPL_ROOT_DIR=${GEOSX_TPL_DIR}
export GEOS_DIR=/data_local/appli_local/MTS/GEOSX/GEOSX/origin/develop-22-f45f894/install-pangea3-gcc8.4.1-openmpi-4.1.2-release/

export CAMP_DIR=/appli_RD/JIEMENG/SEMproxy/camp
export adiak_DIR=${SEM_TPL_ROOT_DIR}/adiak/lib/cmake/adiak/
export RAJA_DIR=${SEM_TPL_ROOT_DIR}/raja/lib/cmake/raja
export RAJA_INCLUDE_DIR=${SEM_TPL_ROOT_DIR}/raja/include/
export CHAI_DIR=/appli_RD/JIEMENG/ProgrammingModels/CHAI/install_p3/share/chai/cmake/
export CHAI_INCLUDE_DIR=/appli_RD/JIEMENG/ProgrammingModels/CHAI/install_p3/include/
export UMPIRE_DIR=/appli_RD/JIEMENG/ProgrammingModels/CHAI/install_p3/lib/cmake/umpire
export LVARRAY_DIR=/appli_RD/JIEMENG/ProgrammingModels/LvArray/install-p3-release-cuda/share/lvarray/cmake/
export LVARRAY_INCLUDE_DIR=/appli_RD/JIEMENG/ProgrammingModels/LvArray/install-p3-release-cuda/include
export CALIPER_DIR=${SEM_TPL_ROOT_DIR}/caliper/share/cmake/caliper
export caliper_INCLUDE_DIR=${SEM_TPL_ROOT_DIR}/caliper/include
export LD_LIBRARY_PATH=${SEM_TPL_ROOT_DIR}/caliper/lib64/:$LD_LIBRARY_PATH
