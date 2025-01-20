# Welcome to the  ProxyApp project!

The ProxyApp project  collects a suite of simple codes representing real applications.
It is intended to be a standard tool for evaluating and comparing the performance of different high-performance computing (HPC) systems, particularly those used for scientific simulations.


# Actual applications 

Current implementation of the proxyApp includes SEM ( Spectral finite Element Methods) and FD ( Finite Differences methods) to solve 2nd order acoustic wave equation in 2D and 3D spaces:  
* The SEM proxy application is a benchmark designed to simulate wave propagation using the spectral element method (SEM), which is a Galerkin-based finite element method for solving partial differential equations.  
* The FD proxy applicaton is a benchmark designed to simulate wave propagation using finite differences stencils operators for solving partial differential equations.  

One of the key features of the SEM and FD proxy benchmarks are their adaptability to different programming models and HPC architectures. This makes them a useful proxy applications for advancing the state of the art in high-performance computing. In addition to their technical capabilities, they are designed to be easy to build and use, and therefore accessible to a wide range of users, from researchers to developers.

# What Programming Models and data containers ?

- The programming models available in the current ProxyApp implementations include:  
    * OpenMP [https://www.openmp.org/] to parallelizing for loops  
    * RAJA [https://raja.readthedocs.io/en/develop/]  
    * KOKKOS [https://kokkos.github.io/kokkos-core-wiki/]  
    
    RAJA, KOKKOS and other third-party libraries must first be compiled and installed  from the [Third-Party Libraries repository](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl) - described at [**step 1**, below](#quick-start-to-compile-and-install).  

- The data containers availbable in the current ProxyApp implementations include:   
    * LvArray [https://lvarray.readthedocs.io/en/latest/]  
    * C++ std::vector  

# Quick Start to compile and install

First consider referring to the page on the [prerequisites](./INSTALL_PREREQUISITES.md) needed.  
## Start by getting the source codes 
Using the following Git commands
```
git clone --recursive https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc   
git clone --recursive https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl 
```
will  create two folders `proxy-geos-hc` and `proxy-geos-hc_tpl`. The `--recursive` option allows to ship the relevant submodules: [BLT](https://github.com/LLNL/blt) for both repositories and  [LvArray](https://github.com/GEOS-DEV/LvArray) specifically for  `proxy-geos-hc`.    
## Step 1: [Build the Third-Party Libraries](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl)

## Step 2: Build and Install the ProxyApp
 
 1. Please make sure that the `proxy_tpl_dir` and `install_tpl_folder` are exported, by sourcing the script `source proxy-geos-hc_tpl/env_var.sh`. These variables are required to set the directory paths of the TPLs and to generate the build files.  
2. Generate the Makefile and build the executable by running the following commandlines. *Be aware that for consistency, the `config_proxy-app.cmake`  file must include the same `config_tpls.cmake` file which has been used when building the [Third-Party Libraries](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl/).*    

```
cd proxy-geos-hc  
cmake  -DCMAKE_BUILD_TYPE=RELEASE <KOKKOS_RAJA_SETUP> -C configs/config_proxy-app.cmake -B build -DCMAKE_INSTALL_PREFIX=install -S .
cd build  
make 
```
This will build and install the executable in the folder `install`.  `<KOKKOS_RAJA_SETUP>` is the placeholder for the configuration option discussed below.   
 

### Configuration options KOKKOS_RAJA_SETUP

 The `KOKKOS_RAJA_SETUP` is used to specify which model programming and portability enabling library is used. The available options include RAJA and KOKKOS. This enables cross-platform seamingless and abstractions either with respect to the parallel programming model or the data container and the corresponding layout. In the current proxy-app, Lvarray container is used for RAJA while  KOKKOS provides its own container. By default, without any specification for `KOKKOS_RAJA_SETUP`, std::vector container is used.  
#### 1. DEFAULT option
The default option (without any specification for `KOKKOS_RAJA_SETUP`) **is relevant for the sequential or a shared memory parallelization mode**. For the latest, one could set `KOKKOS_RAJA_SETUP` as `-DUSE_OMP=ON`.

#### 2. RAJA  with OPENMP and GPU
To use RAJA, set `KOKKOS_RAJA_SETUP` as `-DUSE_RAJA=ON`. This option is only compatible when the OpenMP (on the host) and GPU features are enabled in `proxy-geos-hc_tpls/configs/config_models.cmake`.  

#### 3. KOKKOS with OPENMP and GPU
To use KOKKOS, set `KOKKOS_RAJA_SETUP` as `-DUSE_KOKKOS=ON`. This option is compatible with any combination of programming models. When none of the programming models is enabled, it is equivalent to a serial or sequential mode. 
## Step 3: Run the executable 
The executables are installed in `proxy-geos-hc/install/bin`folder  and can be run as follow:   
```
proxy-geos-hc/install/bin/{proxyName}_{LIB-MODELS}.exe (with proxyName: sem or fd)
```
The tag `LIB-MODELS` is  `KOKKOS_RAJA_SETUP` and enabled programming models dependent. The first part `LIB` is used as a label identifying the name of the portability enabling library (`Kokkos`, `Raja` or empty for the default configuration). It is suffixed by a tag related to the enabled model on the host and the `CUDA_ARCH` flag of the device if a GPU-acceleration is required.   As an example, if KOKKOS is used and a shared-memory parallelization enabled on the host in addition to a  GPU accelaration on the device,  `LIB-MODELS=Kokkos-hOMP_d<CUDA_ARCH>`.  

# Tips and tricks
Some tips and tricks addressing common problems that you may encountered are reported [here](./TIPS_AND_TRICKS.md).