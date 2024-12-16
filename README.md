# Welcome to the  ProxyApp project!

The ProxyApp project  collects a suite of simple codes representing real applications.
It is intended to be a standard tool for evaluating and comparing the performance of different high-performance computing (HPC) systems, particularly those used for scientific simulations.


# Actual applications 

Current implementation of the proxyApp includes SEM ( Spectral finite Element Methods) and FD ( Finite Differences methods) to solve 2nd order acoustic wave equation in 2D and 3D spaces:  
* The SEM proxy applicaton is a benchmark designed to simulate wave propagation using the spectral element method (SEM), which is a Galerkin-based finite element method for solving partial differential equations.  
* The FD proxy applicaton is a benchmark designed to simulate wave propagation using finite differences stencils operators for solving partial differential equations.  

One of the key features of the SEM and FD proxy benchmarks are their adaptability to different programming models and HPC architectures. This makes them a useful proxy applications for advancing the state of the art in high-performance computing. In addition to their technical capabilities, they are designed to be easy to build and use, and therefore accessible to a wide range of users, from researchers to developers.

# What Programming Models and data containers ?

- The programming models available in the current ProxyApp implementations include:  
    * OpenMP [https://www.openmp.org/] to parallelizing for loops  
    * RAJA [https://raja.readthedocs.io/en/develop/]  
    * KOKKOS [https://kokkos.github.io/kokkos-core-wiki/]  
    
    RAJA, KOKKOS and other third-party libraries must first be compiled and installed  from the [third-party libraries repository](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl) - described at [**step 1**, below](#Quick Start to compile and install).  

- The data containers availbable in the current ProxyApp implementations include:  
    * LvArray [https://lvarray.readthedocs.io/en/latest/]  
    * C++ std::vector  

# Quick Start to compile and install

First consider referring to the page on the [prerequisites](./INSTALL_PREREQUISITES.md) needed.

## Step 1: [Build the Third-Party Libraries](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl)

## Step 2: Compile and Install ProxyApp

There are some environment variables of prototype `{LIBRARY}_DIR` which are used in the CMake configuring the ProxyApp.  Continue by exporting the paths to the third-party libraries required for linking and making the dependencies.  
This can be achieved by editing and sourcing the `<path-to-ProxyApp/env>/env.sh` script. Note that `PROXY-GEOSX-HC` is the default name for the `ProxyApp` folder when running the `git clone` command for the main branch of the `ProxyApp` repository.  
  
Create a build folder `buildProxyApp` from where the executable will be built and installed  
```
 cd <path-to-buildProxyApp>  
 cmake  -DCMAKE_BUILD_TYPE=RELEASE CUDA_KOKKOS_RAJA_SETUP -C config.cmake -B . -DCMAKE_INSTALL_PREFIX=<path-to-installProxyApp> -S <path-to-ProxyApp>
 make && make install
```

For consistency, one should use the same `config.cmake` script as the one used when [building the third-party-libraries](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl/), that is `<path-to-third-party-lib/configs>/config.cmake`. The ```CUDA_KOKKOS_RAJA_SETUP``` is discussed below. 
 

### Configuration options CUDA_KOKKOS_RAJA_SETUP

 The ```CUDA_KOKKOS_RAJA_SETUP``` is a combination of three options which are used for :  
- enabling a shared memory programming model `ENABLE_OPENMP=ON`  
- specifying whether building or not a GPU-accelerated application - depending on the vendor through `[ENABLE_CUDA|ENABLE_HIP|ARM]=ON`    
- considering performance portability accross various plateforms, with the use of KOKKOS or RAJA library `[USE_KOKKOS|USE_RAJA]=ON`. 
This enables abstractions either with respect to the parallel programming model or the data container and the corresponding layout.     

#### DEFAULT
The default compilation (without any specification for `CUDA_KOKKOS_RAJA_SETUP`) is the sequential mode with std::vector implementation. 

#### OPEN_MP
In the case of OPENMP std::vector container is used, and  `CUDA_KOKKOS_RAJA_SETUP = -DUSE_OMP=ON`.

#### RAJA + OPEN_MP + CUDA
In the case of RAJA Lvarray container is used, and `CUDA_KOKKOS_RAJA_SETUP` is set as:  
* on Nvidia GPUs: `-DUSE_RAJA=ON -DENABLE_CUDA=ON`  
* on AMD GPUs: `-DUSE_RAJA=ON -DENABLE_HIP=ON`   
* on Nvidia Grace-Hopper: `-DUSE_RAJA=ON -DENABLE_CUDA=ON -DARM=ON`  

#### KOKKOS + OPEN_MP + CUDA
In the case of KOKKOS which provides its own container, `CUDA_KOKKOS_RAJA_SETUP` is set as:  
* on Nvidia GPUs: `-DUSE_KOKKOS=ON -DENABLE_CUDA=ON`   
* on AMD GPUs: `-DUSE_KOKKOS=ON -DENABLE_HIP=ON`   
* on Nvidia Grace-Hopper: `-DUSE_KOKKOS=ON -DENABLE_CUDA=ON -DARM=ON`  

## Step 3: Run the executable 
The executables are installed in `<path-to-installProxyApp/bin>` and can be run as follow:   
```
<path-to-installProxyApp/bin>/{proxyName}_{SETTINGFLAG}.exe (with proxyName: sem or fd)
```
The argument `SETTINGFLAG` tag is  `CUDA_KOKKOS_RAJA_SETUP` dependent. For the for the default option `SETTINGFLAG=SEQUENTIAL`. 

# Reporting issues and things to be improved 
This ProxyApp will evolve with the aim of easing the deployment and portability on various HPC machines. The specific observations and things to be improved are reported [here](./NOTES_ISSUES.md).