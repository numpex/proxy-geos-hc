# Welcome to the  proxyApp project!

The proxyApp project  collects a suite of simple codes representing real applications.
It is intended to be a standard tool for evaluating and comparing the performance of different high-performance computing (HPC) systems, particularly those used for scientific simulations.


# Actual applications 

Current implementation of the proxyApp includes SEM ( Spectral finite Element Methods) and FD ( Finite Differences methods) to solve 2nd order acoustic wave equation in 2D and 3D spaces:  
* The SEM proxy applicaton is a benchmark designed to simulate wave propagation using the spectral element method (SEM), which is a Galerkin-based finite element method for solving partial differential equations.  
* The FD proxy applicaton is a benchmark designed to simulate wave propagation using finite differences stecnils operators for solving partial differential equations.

One of the key features of the SEM and FD proxy benchmarks are their adaptability to different programming models and HPC architectures. This makes them a useful proxy applications for advancing the state of the art in high-performance computing. In addition to their technical capabilities, they are designed to be easy to build and use. This makes them accessible to a wide range of users, from researchers to developers.

# What Programming Models users can select to run  proxy applications?

The programming models included in the current sem proxy implementations include:  
* OpenMP [https://www.openmp.org/] to parallelizing for loops  
* RAJA [https://raja.readthedocs.io/en/develop/]  
* KOKKOS [https://kokkos.github.io/kokkos-core-wiki/]  

RAJA and Kokkos must first be compiled and installed  in `TPL4ProxyApp` repository - described at **step 1**, below.

# What data containers users can select to run SEM proxy?

The data containers included in the current sem proxy implementations include:  
* LvArray [https://lvarray.readthedocs.io/en/latest/]  
* C++ std::vector  

# Quick Start to compile and install

First consider referring to the page on the prerequisites needed for [the general case](./README_prerequisites_general.md); while a specific configuration on a [Ubuntu 24.04 LTS](./README_prerequisites_specific.md) is provided.

## Step 1: [Building the third-party libraries](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl)

## Step 2: Compile and Install ProxyApp
Continue by exporting the path to the previously compiled TPLs required for linking and making the dependencies. This can be achieved by editing and sourcing the script `env/env.sh` in the `PROXY-GEOSX-HC` repository.  
Create a build folder `buildProxyApp` from where the executable will be built and installed  
```
  cmake -DCMAKE_INSTALL_PREFIX=.  CUDA_KOKKOS_RAJA_SETUP /path/to/PROXY-GEOSX-HC 
  make; make install
```

The ```CUDA_KOKKOS_RAJA_SETUP``` is a combination of the three options `USE_KOKKOS`, `ENABLE_CUDA` and `ENABLE_OPENMP` and is discussed below.  

### Available configuration options CUDA_KOKKOS_RAJA_SETUP
The default compilation (without any specification for `CUDA_KOKKOS_RAJA_SETUP`) is sequential mode with std::vector implementation. 
So you will get an executable named "sem_SEQUENTIAL.exe" and "fd_SEQUENTIAL.exe" in your installation directory. To enable shared memory or GPU-accelerated parallelization, the following options are available. 

#### OPEN_MP
in the case of OPENMP std::vector container is used,  CUDA_KOKKOS_RAJA_SETUP = `-DUSE_OMP=ON`

#### RAJA + OPEN_MP + CUDA
in the case of RAJA Lvarray container is used. `CUDA_KOKKOS_RAJA_SETUP` is set as:  
* on Nvidia GPUs: `-DUSE_RAJA=ON -DENABLE_CUDA=ON`  
* on AMD GPUs: `-DUSE_RAJA=ON -DENABLE_HIP=ON`   
* on Nvidia Grace-Hopper: `-DUSE_RAJA=ON -DENABLE_CUDA=ON -DARM=ON`  

#### KOKKOS + OPEN_MP + CUDA
in the case of KOKKOS which provides its own container, `CUDA_KOKKOS_RAJA_SETUP` is set as:  
* on Nvidia GPUs: `-DUSE_KOKKOS=ON -DENABLE_CUDA=ON`  
* on AMD GPUs: `-DUSE_KOKKOS=ON -DENABLE_HIP=ON`   
* on Nvidia Grace-Hopper: `-DUSE_KOKKOS=ON -DENABLE_CUDA=ON -DARM=ON`  

## Step 3: run the executable 
The executable is located in the `buildProxyApp/bin` folder and can be run as follow: for the default option
```
   path/to/buildProxyApp/bin/proxyName_SEQUENTIAL.exe (with proxyName: sem or fd)
```

# Reporting issues and things to be improved 
This ProxyApp will evolves with the aim of easing the deployment and portability on various HPC machines. The specific observations and things to be improved are reported [here](./README_SPECIFIC_and_TODO.md).