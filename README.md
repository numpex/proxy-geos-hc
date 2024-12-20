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
    
    RAJA, KOKKOS and other third-party libraries must first be compiled and installed  from the [third-party libraries repository](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl) - described at [**step 1**, below](#quick-start-to-compile-and-install).  

- The data containers availbable in the current ProxyApp implementations include:   
    * LvArray [https://lvarray.readthedocs.io/en/latest/]  
    * C++ std::vector  

# Quick Start to compile and install

First consider referring to the page on the [prerequisites](./INSTALL_PREREQUISITES.md) needed.  

As a convention, we use the generic notation `{ARGUMENT}` for an `argument` that may take several values. The angle brackets `<CONFIG_OPTIONS>` are used as placeholder for configuration options.      

## Start by getting the source codes 
Using the following Git commands
```
git clone https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc   
git clone https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl 
```

will  create two folders `proxy-geos-hc` and `proxy-geos-hc_tpl`.   
## Step 1: [Build the Third-Party Libraries](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl)

## Step 2: Build and Install the ProxyApp
 
 1. Edit the script `proxy-geos-hc/env/env.sh` with the right paths.  This script contains some environment variables of prototype `{LIBRARY}_DIR` which are used in the CMake configuring the ProxyApp; `{LIBRARY}` being a generic notation for any of the third-party libraries, for instance RAJA, KOKKOS.  
2. Source the script: `source proxy-geos-hc/env/env.sh `  
 
3. Generate the Makefile and build the executable by running the following commandline. *For consistency, make sure that the same `config.cmake` file is used when building both [third-party-libraries](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl/) and the current ProxyApp.*   

```
cmake  -DCMAKE_BUILD_TYPE=RELEASE <CUDA_KOKKOS_RAJA_SETUP> -C proxy-geos-hc_tpl/configs/config.cmake -B build -DCMAKE_INSTALL_PREFIX=install -S proxy-geos-hc
cd build  
make 
```
This will build and install the executable in the folder `install`.  
The ```CUDA_KOKKOS_RAJA_SETUP``` is discussed below. 
 

### Configuration options CUDA_KOKKOS_RAJA_SETUP

 The ```CUDA_KOKKOS_RAJA_SETUP``` is a combination of three options which are used for :  
 - considering model programming and portability enabling library such as KOKKOS or RAJA.  This enables abstractions either with respect to the parallel programming model or the data container and the corresponding layout.   In the case of RAJA Lvarray container is used while  KOKKOS provides its own container. By default std::vector container is used.  
- enabling or not a shared memory programming model `ENABLE_OPENMP=ON`  
- specifying whether building or not a GPU-accelerated application - depending on the vendor through `[ENABLE_CUDA|ENABLE_HIP|ARM]=ON`    

Below are some examples of usual configurations.  
#### SEQUENTIAL
The default compilation (without any specification for `CUDA_KOKKOS_RAJA_SETUP`) is the sequential mode. 

#### OPEN_MP
`CUDA_KOKKOS_RAJA_SETUP = -DUSE_OMP=ON`.

####  RAJA  with OPENMP and GPU
`CUDA_KOKKOS_RAJA_SETUP` is set as:  
- on Nvidia GPUs: `-DUSE_RAJA=ON -DENABLE_CUDA=ON -DENABLE_OPENMP=ON`  
- on AMD GPUs: `-DUSE_RAJA=ON -DENABLE_HIP=ON -DENABLE_OPENMP=ON`   
- on Nvidia Grace-Hopper: `-DUSE_RAJA=ON -DENABLE_CUDA=ON -DARM=ON -DENABLE_OPENMP=ON`  

#### KOKKOS with OPENMP and GPU
`CUDA_KOKKOS_RAJA_SETUP` is set as:  
- on Nvidia GPUs: `-DUSE_KOKKOS=ON -DENABLE_CUDA=ON`   
- on AMD GPUs: `-DUSE_KOKKOS=ON -DENABLE_HIP=ON`   
- on Nvidia Grace-Hopper: `-DUSE_KOKKOS=ON -DENABLE_CUDA=ON -DARM=ON -DENABLE_OPENMP`  

## Step 3: Run the executable 
The executables are installed in `install/bin`folder  and can be run as follow:   
```
install/bin/{proxyName}_{SETTINGFLAG}.exe (with proxyName: sem or fd)
```
The argument `SETTINGFLAG` tag is  `CUDA_KOKKOS_RAJA_SETUP` dependent. It is used as a label identifying the name of the portability enabling library used (`KOKKOS` or `RAJA`) . For the default option `SETTINGFLAG=SEQUENTIAL` when neither of KOKOSS nor RAJA is used. 

# Tips and tricks
Some tips and tricks addressing common problems that you may encountered are reported [here](./TIPS_AND_TRICKS.md).