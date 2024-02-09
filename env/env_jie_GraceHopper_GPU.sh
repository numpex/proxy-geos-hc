module use /shared/data1/Projects/CSE_HPC/apps/test/gheval01_test/modules/
module load gcc cmake cuda/12.3.107 openmpi-gcc/4.1.6/cuda.12.3

export SEM_TPL_ROOT_DIR=/shared/data1/Users/j0535952/work2024/proxys_2024/tpl4ProxyApp/installTPL_GraceHopper_GPU
export RAJA_DIR=${SEM_TPL_ROOT_DIR}/raja/lib/cmake/raja
export CHAI_DIR=${SEM_TPL_ROOT_DIR}/chai/share/chai/cmake/
export CAMP_DIR=${SEM_TPL_ROOT_DIR}/raja/lib/cmake/camp
export UMPIRE_DIR=${SEM_TPL_ROOT_DIR}/chai/lib/cmake/umpire/
export KOKKOS_DIR=${SEM_TPL_ROOT_DIR}/kokkos/lib64/cmake/Kokkos
export KOKKOS_INCLUDE_DIR=${SEM_TPL_ROOT_DIR}/kokkos/include
export CUDA_ROOT=/hrtc/apps/cuda/12.3.107/aarch64/rocky9/


# don't forget to change CMakeLists.txt 
# option (X86_cypress "Compilation on Cypress Cluster" OFF)
# ption (ARM "ARM aarch64 architecture: e.g. grace CPU" ON)

# for kokkos: cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_KOKKOS=ON -DENABLE_CUDA=ON .. ; make 
# for raja: cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_RAJA=ON -DENABLE_CUDA=ON .. ; make 
