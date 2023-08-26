#export MODULEPATH=/data_local/appli_local/MTS/GEOSX/modulefiles/:$MODULEPATH
#module load cmake/3.21.4 gcc/8.4.1 cuda/11.0.3 ompi/4.1.2 openblas/0.3.18 python4geosx/3.8.5-gcc-8.4.1-ompi-4.1.2

export SEM_TPL_ROOT_DIR=~/src/geosx/codes/proxyApp/buildInstall/install
export RAJA_DIR=${SEM_TPL_ROOT_DIR}/raja/lib/cmake/raja
export CAMP_DIR=${SEM_TPL_ROOT_DIR}/raja/lib/cmake/camp
export CHAI_DIR=${SEM_TPL_ROOT_DIR}/chai/share/chai/cmake/
export UMPIRE_DIR=${SEM_TPL_ROOT_DIR}/chai/lib/cmake/umpire/
kokkos_install_dir=~/src/geosx/codes/proxyApp/buildInstall/install/kokkos
export KOKKOS_DIR=$kokkos_install_dir/lib64/cmake/Kokkos
export KOKKOS_INCLUDE_DIR=$kokkos_install_dir/include

