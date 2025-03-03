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
    
    RAJA, KOKKOS and other Third-Party Libraries (TPLs) can be compiled and installed: either from source as described  within the [third-party libraries repository](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl) or using a Package Manager (Guix or Spack).  

- The data containers availbable in the current proxyApp implementations include:   
    * LvArray [https://lvarray.readthedocs.io/en/latest/]  
    * C++ std::vector  

# Quick Start to compile and install

First consider referring to the page on the [prerequisites](./INSTALL_PREREQUISITES.md) needed.  

As a convention, the angle brackets `<variable>` are used as placeholder for *variable* or *option*.     

## Start by getting the ProxyApp source codes 
Using the following Git command
```
git clone --recursive https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc   
```
will  create the folder `proxy-geos-hc`. The `--recursive` option allows to ship the relevant submodules: [BLT](https://github.com/LLNL/blt) and  [LvArray](https://github.com/GEOS-DEV/LvArray).    

The ProxyApp depends on some third-party libraries which can be installed using one of the two options described below. They include a build from source and a build using the Guix package manager. For the former, using the following Git command
```   
git clone --recursive https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl 
```
will fetch the third-party libraries source codes in the folder `proxy-geos-hc_tpl`.  

## Environment variables for the ProxyApp build 
Some environment variables are required to configure the CMake when building the [ProxyApp](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc#step-2-build-and-install-the-proxyapp). Their choice is specific to the method used for installing the TPLs. The [table](####environment-variables-in-env_var.sh-for-the-cMake-build) below summarizes the set up for the two build options.   

#### Environment variables in env_var.sh for the CMake build

| Build Option       |    proxy_config_root      | config_proxy | Sourcing required for| other environment variables|  
|: ---------- :|: ------------------- :| :-----------:| :------------------- :| :------------------- :|  
| From source   |  `proxy-geos-hc_tpl` | `config_<machine's name>.cmake` |TPLs and ProxyApp|`install_tpl` and `build_tpl` |   
| Using Guix   | `proxy-geos-hc` |`config_<machine's name>_guix.cmake` |ProxyApp only|  |   
   
 1. Edit the script `proxy-geos-hc/env_var.sh` with the right arguments for the environment variables, as discribed in the [table](####environment-variables-in-env_var.sh-for-the-cMake-build):   
   - `proxy_config_root` *the root path of the `configs` folder*  
   - `config_proxy` *the name of the config file to be used to pre-load the cache when building the ProxyApp, and [the TPLs when building from source](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl/)*.  
   - `build_tpl` and `install_tpl` *the name of its binary and install directories where the libraries are to be built and installed*, in the case where the TPLs are built from source   
  2. Source the script `source proxy-geos-hc/env_var.sh` to export these variables.  

The next step involves setting up the configuration file to be used during the cmake build of the proxyApp. They are also required when building [ the TPLs from source](###from-source). 

## Step 1. [Edit the configuration file for the build process](./Notes_config_setting.md)


## Step 2.  Installing the TPLs dependencies

Use one of the two options described below to install the TPLs

### From Source

The build and install of the TPLs from code proceeds as follows:  
- Source the env_var.sh script  
- [From source Build of the Third-Party Libraries](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl)  

### Using a Package Manager

The alternative option involves installing the TPLs using a Package Manager and is described on the following page 
#### [Guix package manager Build of the Third-Party Libraries](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc/-/issues/3)
It allows to build the ProxyApp within a containerized build environment with all the TPLs dependencies.   

## Step 3. Build and Install the ProxyApp

 1. Consider [exporting the environment variables](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc##environment-variables-for-the-proxyapp-build)  that are required, by sourcing the `env_var.sh` file. They are required for the config file `proxy-geos-hc/configs/config_proxy-app.cmake`. It serves as a wrapper for the config file `${proxy_config_root}/configs/${config_proxy}` that would have been used to pre-load the cache when building the TPLs from source.  
3. Generate the Makefile and build the executable by running the following command lines 
```
cd proxy-geos-hc  
cmake  -DCMAKE_BUILD_TYPE=RELEASE <KOKKOS_RAJA_OMP> -DGUIX_INSTALLED_TPL=<BOOL> -C configs/config_proxy-app.cmake -B ${build_tpl} -DCMAKE_INSTALL_PREFIX=${install_tpl} -S .
cd ${build_tpl}  
make && make install
```
This will build and install the executable in the folder `build`. The boolean `BOOL`, for the argument `GUIX_INSTALLED_TPL`, is used to specify whether the TPLs have been installed from source or using a Package Manager. The configuration option `KOKKOS_RAJA_OMP` is discussed below.    

### Configuration option KOKKOS_RAJA_OMP

 The `KOKKOS_RAJA_OMP` is used to specify which programming model and portability enabling library is used. The available options include RAJA and KOKKOS. This enables cross-platform seamingless and abstractions either with respect to the parallel programming model or the data container and the corresponding layout. In the current proxyApp, Lvarray container is used for RAJA while  KOKKOS provides its own container. In cases where neither RAJA nor KOKKOS is used, std::vector container is used.  
 Some ready-to-use command lines for each of these configurations are provided [here - HowToInstall](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl/-/blob/develop/howToInstall.md?ref_type=heads).   
#### 1. DEFAULT option
The default option (without any specification for `KOKKOS_RAJA_OMP`) is sequential mode. [**Not supported at the moment**](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc/-/issues/8).  

#### 2. OpenMP
To use OMP, set `CUDA_KOKKOS_RAJA_OMP` as `-DUSE_OMP=ON`, for a shared-memory parallelization mode. [**Not  supported at the moment**](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc/-/issues/8).  

#### 3. RAJA  with OPENMP and GPU
To use RAJA, set `KOKKOS_RAJA_OMP` as `-DUSE_RAJA=ON`. This option is only valid when the OpenMP and GPU features are enabled in `proxy-geos-hc_tpls/configs/config_<machine's name>.cmake` - See [What Programming Models for the TPLs](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl##11programming-models-enabled-for-the-tpls).  

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
