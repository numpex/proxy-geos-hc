#/ enable OpenMP for Kokkos

export SEM_TPL_ROOT_DIR=/appli_RD/JIEMENG/SEMproxy/proxyAppSEM_GPU/installTPL5
export RAJA_DIR=${SEM_TPL_ROOT_DIR}/raja/lib/cmake/raja
export CHAI_DIR=${SEM_TPL_ROOT_DIR}/chai/share/chai/cmake/
export CAMP_DIR=${SEM_TPL_ROOT_DIR}/raja/lib/cmake/camp
export UMPIRE_DIR=${SEM_TPL_ROOT_DIR}/chai/lib/cmake/umpire/
export KOKKOS_DIR=${SEM_TPL_ROOT_DIR}/kokkos/lib64/cmake/Kokkos
export KOKKOS_INCLUDE_DIR=${SEM_TPL_ROOT_DIR}/kokkos/include
