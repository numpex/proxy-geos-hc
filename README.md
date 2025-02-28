# Welcome to the  ProxyApp Project!

The proxyApp project  collects a suite of simple codes representing real applications.
It is intended to be a standard tool for evaluating and comparing the performance of different high-performance computing (HPC) systems, particularly those used for scientific simulations.


# Actual applications 

Current implementation of the proxyApp includes SEM (Spectral finite Element Methods) and FD (Finite Differences methods) to solve 2nd order acoustic wave equation in 2D and 3D spaces:  
* The SEM proxy application is a benchmark designed to simulate wave propagation using the spectral element method (SEM), which is a Galerkin-based finite element method for solving partial differential equations.  
* The FD proxy applicaton is a benchmark designed to simulate wave propagation using finite differences stencils operators for solving partial differential equations.  

One of the key features of the SEM and FD proxy benchmarks are their adaptability to different programming models and HPC architectures. This makes them a useful proxy applications for advancing the state of the art in high-performance computing. In addition to their technical capabilities, they are designed to be easy to build and use, and therefore accessible to a wide range of users, from researchers to developers.

# What Programming Models and data containers ?

- The programming models available in the current proxyApp implementations include:  
    * OpenMP [https://www.openmp.org/] to parallelizing for loops  
    * RAJA [https://raja.readthedocs.io/en/develop/]  
    * KOKKOS [https://kokkos.github.io/kokkos-core-wiki/]  
    
    RAJA, KOKKOS and other Third-Party Libraries (TPLs) must first be compiled and installed  from the [third-party libraries repository](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl) - [**step 1**, below](#quick-start-to-compile-and-install).  

- The data containers availbable in the current proxyApp implementations include:   
    * LvArray [https://lvarray.readthedocs.io/en/latest/]  
    * C++ std::vector  

# Quick Start to compile and install

First consider referring to the page on the [prerequisites](./INSTALL_PREREQUISITES.md) needed.  

As a convention, the angle brackets `<variable>` are used as placeholder for *variable* or *option*.     

## Start by getting the source codes 
Using the following Git commands
```
git clone --recursive https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc   
git clone --recursive https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl 
```
will  create two folders `proxy-geos-hc` and `proxy-geos-hc_tpl`. The `--recursive` option allows to ship the relevant submodules: [BLT](https://github.com/LLNL/blt) for both repositories and  [LvArray](https://github.com/GEOS-DEV/LvArray) specifically for  `proxy-geos-hc`.    
## Step 1: [Build the Third-Party Libraries](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl)

## Step 2: Build and Install the ProxyApp

 1. Consider [exporting the environment variables defined at Step 1](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl#step-2-some-environment-variables-for-the-build). They are required for the config file `proxy-geos-hc/configs/config_proxy-app.cmake`, which serves as a wrapper for the config file (`config_<machine's name>.cmake`) that has been used to pre-load the cache when building the TPLs.  
3. Generate the Makefile and build the executable by running the following command lines 
```
cd proxy-geos-hc  
cmake  -DCMAKE_BUILD_TYPE=RELEASE <KOKKOS_RAJA_OMP> -DGUIX_INSTALLED_TPL=<BOOL> -C configs/config_proxy-app.cmake -B build -DCMAKE_INSTALL_PREFIX=install -S .
cd build  
make && make install
```
This will build and install the executable in the folder `build`. The configuration option `KOKKOS_RAJA_OMP` is discussed below.    
 

### Configuration option KOKKOS_RAJA_OMP

 The `KOKKOS_RAJA_OMP` is used to specify which programming model and portability enabling library is used. The available options include RAJA and KOKKOS. This enables cross-platform seamingless and abstractions either with respect to the parallel programming model or the data container and the corresponding layout. In the current proxyApp, Lvarray container is used for RAJA while  KOKKOS provides its own container. In cases where neither RAJA nor KOKKOS is used, std::vector container is used.  
 Some ready-to-use command lines for each of these configurations are provided [here - HowToInstall](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl/-/blob/develop/howToInstall.md?ref_type=heads).   
#### 1. DEFAULT option
The default option (without any specification for `KOKKOS_RAJA_OMP`) is sequential mode. [**Not supported at the moment**](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc/-/issues/8).  

#### 2. OpenMP
To use OMP, set `CUDA_KOKKOS_RAJA_OMP` as `-DUSE_OMP=ON`, for a shared-memory parallelization mode. [**Not  supported at the moment**](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc/-/issues/8).  

#### 3. RAJA  with OPENMP and GPU
To use RAJA, set `KOKKOS_RAJA_OMP` as `-DUSE_RAJA=ON`. This option is only valid when the OpenMP and GPU features are enabled in `proxy-geos-hc_tpls/configs/config_<machine's name>.cmake` - See [What Programming Models for the TPLs](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl#step-11-programming-models-enabled-for-the-tpls).  

#### 4. KOKKOS with OPENMP and GPU
To use KOKKOS, set `KOKKOS_RAJA_OMP` as `-DUSE_KOKKOS=ON`. This option is compatible with any combination of programming models. When none of the programming models is enabled, it is equivalent to a serial or sequential mode.   
## Step 3: Run the executable 
The executables are installed in the `proxy-geos-hc/install/bin`folder. The corresponding names have a specific prototype, which accounts of several inputs, and they can be run as follows:   
```
proxy-geos-hc/install/bin/<proxyName>_<LIB>_<HostModel>_<DEVICE>.exe 
```
1.  `proxyName: sem or fd` since the executables are installed for both FD and SEM solvers  
2. `LIB` is used as a label identifying the name of the abstraction enabling library. It is only relevant when KOKKOS or RAJA is used.  
3. `HostModel` is the tag for the programming model enabled on the host. When using OMP or KOKKOS, RAJA with  OMP enabled: `HostModel=OMP`, otherwise for the default option `HostModel=SEQUENTIAL`      
4. The tag of the device `DEVICE` is considered when a GPU-acceleration has been specified.   

For example, if KOKKOS is used and OMP enabled in addition to a GPU acceleration on a Nvidia `RTX2000` device, the following two executables will be installed `fd_Kokkos_OMP_RTX2000`, `sem_Kokkos_OMP_RTX2000`.  
 
# Tips and tricks
Some tips and tricks addressing common problems that you may encountered are reported [here](./TIPS_AND_TRICKS.md).
