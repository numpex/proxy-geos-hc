# Configuration files for the build process

## Overview
In order to ensure consistency and to avoid redundancy, we have centralized all configuration data into an unique set of configuration files that are used to pre-load the cache during the CMake build process for **both** the proxyApp **and** the Third-Party Libraries in case of ![TPL_MIRROR](https://img.shields.io/badge/TPL_from-Mirror-green).

The folder `configs` in the Proxyapp repo includes :
1. A collection of machine-dependent configuration files, whose names follow the `config_<machine's name>.cmake` convention.
2. A file named `config_core.cmake` which contains generic configuration ; this file must be included by all machine-dependent configuration files
3. A file named `config-proxy-app.cmake` which contains all the specific configuration variables needed only in the case of the Proxyapp build - this file will automatically include the machine-dependent configuration file whose name is defined by the environment variable `config_proxy` (see [here](README.md#Environment-variables)).

## Creating your machine-dependent file

Have a look on the existing `config_<machine's name>.cmake` files and initiate a copy of the one that is closer to your architecture. Then adapt it according to your needs.

Those files are made of 4 distinct parts, as described in the table below:

| #Part       | Description         | Need to adjust | 
| ---------- | ------------------- | -----------|
| 1   |  specifies the enabled accelerations | Very likely |  
| 2   | platform specific parameters |Very likely|  
| 3 | TPLs config options for profiling (Caliper, OMPT) | Not likely |  
| 4  | includes the `config_core` which sets most CMake variables for the build |Very unlikely|    

### Part 1. Accelerations enabled
Edit this part to specify the accelerations that are supported by the underlying hardware architecture. You can enable a GPU-accelerated code and a shared memory parallelization on the host.  To do so,  
1. set `ENABLE_OPENMP=ON`,   
2. for the GPU programming models - depending on the vendor, set   
- on Nvidia GPUs: `ENABLE_CUDA=ON`   
- on AMD GPUs: `ENABLE_HIP=ON`   
- on Nvidia Grace-Hopper: `ENABLE_CUDA=ON ARM=ON`  

**Important note** ![TPL_GUIX](https://img.shields.io/badge/TPL_from-Guix-blue): when installing the TPLs with a package manager, make sure to select TPLs packages whose enabled accelerations are consistent with those specified in the part 1 of the `config_<machine's name>.cmake`.   

### Part 2. Platform-specific parameters

Edit this part to set the platform-specific variables.  It emboddies various variables which are used to specify the paths for compilers (for instance gcc, g++, gfortran, mpicc, nvcc, hipcc), for libraries and other compute options:   
- The path to the compilers on the host is typically `/usr/bin`  
- The variable `CUDA_ARCH` or `COMP_ARCH` is used to specify the architecture of the device. The related variables must be set in accordance to the GPU acceleration specified in [Part 1](#part-1-accelerations-enabled). We refer to the following webpage for a mapping between various [GPU microarchitectures and their corresponding flags or compute capabilities](https://kokkos.org/kokkos-core-wiki/keywords.html#gpu-architectures)  
- Some compilation flags for the targetted architecture (`mtune` and `mcpu`) are set through the variable `CPU_TUNE_FLAG`.   

Note that specification of the paths is only needed with the ![TPL_MIRROR](https://img.shields.io/badge/TPL_from-Mirror-green) method. In case of installation with package managers, these paths will be found automatically. This is the reason why all these settings are encompassed in a `if(BUILD_FROM_TPLMIRROR)`.

