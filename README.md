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
    
- The data containers availbable in the current proxyApp implementations include:   
    * LvArray [https://lvarray.readthedocs.io/en/latest/]  
    * C++ std::vector  

# How to compile and install

First consider referring to the page on the [prerequisites](./INSTALL_PREREQUISITES.md) needed.  

As a convention, the angle brackets `<variable>` are used as placeholder for *variable* or *option*.     

## Start by getting the ProxyApp source codes 
Using the following Git command
```
git clone --recursive https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc   
```
will  create the folder `proxy-geos-hc`. The `--recursive` option allows to ship the relevant submodules: [BLT](https://github.com/LLNL/blt) and  [LvArray](https://github.com/GEOS-DEV/LvArray).    

### TPLs install options 
KOKKOS, RAJA and the other Third-Party Libraries (TPLs) can be compiled and installed either from mirror as described  within the [third-party libraries repository](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl) or using a Package Manager (Guix or Spack).   
- In the former case ![TPL_MIRROR](https://img.shields.io/badge/TPL_from-Mirror-green) you shall consider [Getting the TPLs source code](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl/-/tree/dev_docs/guix?ref_type=heads#getting-the-tpls-source-codes) to create the folder `proxy-geos-hc_tpl` which provides with some tarballs of TPLs specific versions.  
- In the other case ![TPL_GUIX](https://img.shields.io/badge/TPL_from-Guix-blue), installing the TPLs with a Package Manager offers the possibility to build the ProxyApp within a containerized development environment with all the dependencies. 

This differentitation in the TPLs install options is carried through the CMake flag `BUILD_FROM_TPLMIRROR` which is `ON` by default ![TPL_MIRROR](https://img.shields.io/badge/TPL_from-Mirror-green).   

## Step 0. Environment variables for the CMake builds 
At most, two environment variables are required to configure the CMake when building either the TPLs from mirror or the [ProxyApp](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc#step-2-build-and-install-the-proxyapp).   
   
 1. Edit the script `proxy-geos-hc/env_var.sh` with the right arguments for the environment variables:   
   - `config_proxy` *the filename of the config file to be used to pre-load the CMake cache*. Some examples are provided in the folder `proxy-geos-hc/configs`   
   - `install_tpl_dir` *the absolute path of the install directory where the libraries are to be installed* ![TPL_MIRROR](https://img.shields.io/badge/TPL_from-Mirror-green)   
  2. Source the script `proxy-geos-hc/env_var.sh` to export these variables.  

The next step involves setting up the configuration file `$config_proxy` to be used during the CMake builds. The compilers and related CMake variables are relevant for set up only when building from mirror ![TPL_MIRROR](https://img.shields.io/badge/TPL_from-Mirror-green). Alternatively, when installing <!-- ![TPL_GUIX](https://img.shields.io/badge/TPL_from-Guix-blue) --> from a package manager, e.g. ![TPL_GUIX](https://img.shields.io/badge/TPL_from-Guix-blue),  with the aim to build the proxyApp within a containerized development environment, the CMake variables for the compilers are not required and we rely on CMake for a conistent and automatic set up. 

## Step 1. [Edit the configuration file for the CMake build](./Notes_config_setting.md)
## Step 2.  Installing the TPLs dependencies
#### [From mirror](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl/-/tree/dev_docs/guix?ref_type=heads)  ![TPL_MIRROR](https://img.shields.io/badge/TPL_from-Mirror-green)
#### [From Guix Package Manager](./Install_TPLs_Guix.md) ![TPL_GUIX](https://img.shields.io/badge/TPL_from-Guix-blue)

## Step 3. Build and Install the ProxyApp
 1. Export the environment variables by sourcing the `env_var.sh` file, as instructed at [Step 0](##step-0.-Environment-variables-for-the-CMake-builds). They are required to set up the config file `proxy-geos-hc/configs/config_proxy-app.cmake` which serves as a wrapper for the config `${config_proxy}` file that must be used to pre-load the CMake cache for the build. It also sets the relevant libraries paths in case where the ![TPL_MIRROR](https://img.shields.io/badge/TPL_from-Mirror-green) option is used.  
3. Generate the Makefile and build the executable by running the following command lines 
```
cd proxy-geos-hc  
cmake  -DCMAKE_BUILD_TYPE=RELEASE <KOKKOS_RAJA_OMP> -DBUILD_FROM_TPLMIRROR=<BOOL> -C configs/config_proxy-app.cmake -B ./build -DCMAKE_INSTALL_PREFIX=./install -S .
cd build  
make && make install
```
This will build the binaries and install the executables respectively in the folders `build` and `install/bin`. <!--The boolean `BOOL`, for the argument `BUILD_FROM_TPLMIRROR`, is used to specify whether the TPLs have been installed from mirror ![TPL_MIRROR](https://img.shields.io/badge/TPL_from-Mirror-green) or using a Package Manager ![TPL_GUIX](https://img.shields.io/badge/TPL_from-Guix-blue).--> The configuration option `KOKKOS_RAJA_OMP` is discussed below.    

### Configuration option KOKKOS_RAJA_OMP

 The `KOKKOS_RAJA_OMP` is used to specify which programming model and portability enabling library is used. The available options include RAJA and KOKKOS. This enables cross-platform seamingless and abstractions either with respect to the parallel programming model or the data container and the corresponding layout. In the current proxyApp, Lvarray container is used for RAJA while  KOKKOS provides its own container. In cases where neither RAJA nor KOKKOS is used, `std::vector` container is used.  
 Some ready-to-use command lines for each of these configurations are provided [here](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl/-/blob/develop/howToInstall.md?ref_type=heads).   
#### 1. DEFAULT option
The default option (without any specification for `KOKKOS_RAJA_OMP`) is sequential mode. [**Not supported at the moment**](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc/-/issues/8).  

#### 2. OpenMP
To use OMP, set `CUDA_KOKKOS_RAJA_OMP` as `-DUSE_OMP=ON`, for a shared-memory parallelization mode. [**Not  supported at the moment**](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc/-/issues/8).  

#### 3. RAJA  with OPENMP and GPU
To use RAJA, set `KOKKOS_RAJA_OMP` as `-DUSE_RAJA=ON`. This option is only valid when the OpenMP and GPU features are enabled in `proxy-geos-hc_tpls/configs/config_<machine's name>.cmake` - See [Programming Models Enabled for the TPLs](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc/-/blob/dev_docs/guix/Notes_config_setting.md?ref_type=heads#1-programming-models-enabled-for-the-tpls).  

#### 4. KOKKOS with OPENMP and GPU
To use KOKKOS, set `KOKKOS_RAJA_OMP` as `-DUSE_KOKKOS=ON`. This option is compatible with any combination of programming models. When none of the programming models is enabled, it is equivalent to a serial or sequential mode.   

## Step 4. Run the executable 
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
