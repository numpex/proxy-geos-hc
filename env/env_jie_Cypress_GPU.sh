module load gcc/8.3.1 cuda/11.2.152 openmpi-gcc/4.1.0/cuda.11.2 cmake

export SEM_TPL_ROOT_DIR=/shared/data1/Users/j0535952/work2023/SEMCode_2023/tpl4ProxyApp/installTPL_Cypress_GPU
export RAJA_DIR=${SEM_TPL_ROOT_DIR}/raja/lib/cmake/raja
export CHAI_DIR=${SEM_TPL_ROOT_DIR}/chai/share/chai/cmake/
export CAMP_DIR=${SEM_TPL_ROOT_DIR}/raja/lib/cmake/camp
export UMPIRE_DIR=${SEM_TPL_ROOT_DIR}/chai/lib/cmake/umpire/
export KOKKOS_DIR=${SEM_TPL_ROOT_DIR}/kokkos/lib64/cmake/Kokkos
export KOKKOS_INCLUDE_DIR=${SEM_TPL_ROOT_DIR}/kokkos/include
export CUDA_ROOT=/hrtc/apps/cuda/11.2.152/x86_64/

