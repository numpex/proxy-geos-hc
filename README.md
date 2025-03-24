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

First consider referring to the page on the [prerequisites](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc/-/blob/dev_docs/guix/readme-docs/INSTALL_PREREQUISITES.md?ref_type=heads#key-aspects-and-prerequisites--for-proxyapp) needed.  

As a convention, the angle brackets `<variable>` are used as placeholder for *variable* or *option*.     

## Getting the ProxyApp source code
Using the following Git command
```
git clone --recursive https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc   
```
will  create the folder `proxy-geos-hc`. The `--recursive` option allows to ship the relevant submodules: [BLT](https://github.com/LLNL/blt) and  [LvArray](https://github.com/GEOS-DEV/LvArray).    

## Getting Third-Party Libraries 
KOKKOS, RAJA and the other Third-Party Libraries (TPLs) needed by the proxyApp can be obtained by two different methods:
- From a mirror repository gathering tarballs of all TPLs in consistent versions
- Using a Functional Package Manager : currently only Guix is supported, Spack will be suppported latter on.

Labels ![TPL_MIRROR](https://img.shields.io/badge/TPL_from-Mirror-green) and ![TPL_GUIX](https://img.shields.io/badge/TPL_from-Guix-blue) will be used throughout the document to highlight when some specific handling is needed depending on the TPL installation method.

The CMake flag `BUILD_FROM_TPLMIRROR` is used to transmit the information on the actual TPL installation method to the build system.
Default value is `ON`, corresponding to ![TPL_MIRROR](https://img.shields.io/badge/TPL_from-Mirror-green).   

## Environment variables
The following environment variables are expected during the CMake build process:
- `config_proxy`:  the filename of the config file to be used (see [Step 1](#step-1-set-up-the-configuration-file))
- `install_tpl_dir`: the absolute path of the directory where the TPL libraries are installed (see [step 2](#step-2--install-the-tpls-dependencies)) in case of  ![TPL_MIRROR](https://img.shields.io/badge/TPL_from-Mirror-green) method. Not needed with ![TPL_GUIX](https://img.shields.io/badge/TPL_from-Guix-blue).     

The script `proxy-geos-hc/env_var.sh` can be used to keep a record of the relevant values for these environments variables.


## Step 1. Set up the configuration file
The configuration file allows to specify:
- which accelerations are enabled, depending on the underlying hardware capabilities
- the paths for compilers & tools, in case of ![TPL_MIRROR](https://img.shields.io/badge/TPL_from-Mirror-green) method.


1. Follow the [instructions for editing the configuration file](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc/-/blob/dev_docs/guix/readme-docs/Notes_config_setting.md?ref_type=heads#configuration-for-the-build-process). Examples are provided in the folder `proxy-geos-hc/configs`.
2. Make sure to update the `config_proxy` environment variable with the name of your specific configuration file.

## Step 2.  Install the TPLs dependencies
-  ![TPL_MIRROR](https://img.shields.io/badge/TPL_from-Mirror-green): 
   1. Apply the following [instructions](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl/-/tree/dev_docs/guix?ref_type=heads) 
   2. Make sure to update the `install_tpl_dir` environment variable.
- ![TPL_GUIX](https://img.shields.io/badge/TPL_from-Guix-blue):
   1. Apply the following [instructions](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc/-/blob/dev_docs/guix/readme-docs/Install_TPLs_Guix.md?ref_type=heads#getting-guix-package-manager) 


## Step 3. Build and Install the ProxyApp
1. Export the environment variables, e.g. by sourcing the `env_var.sh` file  
2. Generate the Makefile and build the executable by running the following command lines 
```
cd proxy-geos-hc  
cmake  -DCMAKE_BUILD_TYPE=RELEASE <KOKKOS_RAJA_OMP> -DBUILD_FROM_TPLMIRROR=<BOOL> -C configs/config_proxy-app.cmake -B ./build -DCMAKE_INSTALL_PREFIX=./install -S .
cd build  
make && make install
```
This will build the binaries and install the executables respectively in the folders `build` and `install/bin`. 

The `-DBUILD_FROM_TPLMIRROR=<BOOL>` option must be adjusted according to the chosen TPL installation method.

The configuration option `KOKKOS_RAJA_OMP` is discussed below.    

### Configuration option KOKKOS_RAJA_OMP

 The `KOKKOS_RAJA_OMP` is used to specify which programming model and portability enabling library is used. The available options include RAJA and KOKKOS. This enables cross-platform seamingless and abstractions either with respect to the parallel programming model or the data container and the corresponding layout. In the current proxyApp, Lvarray container is used for RAJA while  KOKKOS provides its own container. In cases where neither RAJA nor KOKKOS is used, `std::vector` container is used.  
 Some ready-to-use command lines for each of these configurations are provided [here](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl/-/blob/develop/howToInstall.md?ref_type=heads).   
#### 1. DEFAULT option
The default option (without any specification for `KOKKOS_RAJA_OMP`) is sequential mode. [**Not supported at the moment**](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc/-/issues/8).  

#### 2. OpenMP
To use OMP, set `KOKKOS_RAJA_OMP` as `-DUSE_OMP=ON`, for a shared-memory parallelization mode. [**Not  supported at the moment**](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc/-/issues/8).  

#### 3. RAJA  with OPENMP and GPU
To use RAJA, set `KOKKOS_RAJA_OMP` as `-DUSE_RAJA=ON`. This option is only valid when the OpenMP and GPU features are enabled in `${config_proxy}` - See [Programming Models Enabled for the TPLs](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc/-/blob/dev_docs/guix/readme-docs/Notes_config_setting.md?ref_type=heads#configuration-for-the-build-process).  

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
 
You may encounter path issues for the TPLs libraries. The following lines fix the issue:
```
 export LD_LIBRARY_PATH=/path/to/proxy-geos-hc_tpl/installTPLs/chai/lib/:$LD_LIBRARY_PATH
 export LD_LIBRARY_PATH=/path/to/proxy-geos-hc_tpl/installTPLs/kokkos/lib/:$LD_LIBRARY_PATH
 export LD_LIBRARY_PATH=/path/to/proxy-geos-hc_tpl/installTPLs/raja/lib/:$LD_LIBRARY_PATH
```
# Tips and tricks
Some tips and tricks addressing common problems that you may encountered are reported [here](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc/-/blob/dev_docs/guix/readme-docs/TIPS_AND_TRICKS.md?ref_type=heads#tips-and-tricks).
