module use /shared/data1/Projects/CSE_HPC/apps/modules/amd
module load cmake llvm rocm/6.0.0

export SEM_TPL_ROOT_DIR=/shared/data1/Users/j0535952/work2024/proxys_2024/TPL/tpl4ProxyApp/installTPL_Cypress_mi100GPU
export RAJA_DIR=${SEM_TPL_ROOT_DIR}/raja/lib/cmake/raja
export CHAI_DIR=${SEM_TPL_ROOT_DIR}/chai/share/chai/cmake/
export CAMP_DIR=${SEM_TPL_ROOT_DIR}/raja/lib/cmake/camp
export UMPIRE_DIR=${SEM_TPL_ROOT_DIR}/chai/lib/cmake/umpire/
#export KOKKOS_DIR=${SEM_TPL_ROOT_DIR}/kokkos/lib64/cmake/Kokkos
#export KOKKOS_INCLUDE_DIR=${SEM_TPL_ROOT_DIR}/kokkos/include

export KOKKOS_DIR=/shared/data1/Users/j0535952/work2024/ProgrammingModels/kokkos/install_mi100_hipcc-6.0.0-gfx908/lib64/cmake/Kokkos
export KOKKOS_INCLUDE_DIR=/shared/data1/Users/j0535952/work2024/ProgrammingModels/kokkos/install_mi100_hipcc-6.0.0-gfx908/include


# don't forget to change CMakeLists.txt 
# option (X86_cypress "Compilation on Cypress Cluster" ON)

# for kokkos: cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_KOKKOS=ON -DENABLE_HIP=ON .. ; make install
# for raja: cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_RAJA=ON -DENABLE_HIP=ON .. ; make install

