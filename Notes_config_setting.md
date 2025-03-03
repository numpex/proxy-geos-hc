# Configuration for the build process

The current step follows the export of the [environment variables](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc##environment-variables-for-the-proxyapp-build). One of them defines the root path to the config folder. The folder `${proxy_config_root}/configs` contains the config files for the build process. Apart of the `config_core` file, there are several configuration files which are provided for various architectures and used to pre-load the cache when generating the CMake files at the [build stage](###build-the-third-party-libraries). They can be used as example to adapt the configuration file to any different architecture.   

The template of these machine-specific configuration files whose name follows the prototype `config_<machine's name>.cmake` includes four parts. They are described in the [table below](###table-of-the-config-files) and sorted in decreasing  order of setup likelihood.   

#### Constitutive parts of the `config_<machine's name>.cmake` 
| #Part       | Description         | Setup likelihood | 
| ---------- | ------------------- | -----------|
| 1   |  specifies the enabled programming models | Very high |  
| 2   | sets some platform specific parameters |Very high|  
| 3 | sets some TPLs config options for profiling (Caliper, OMPT) | Low |  
| 4  | includes the `config_core` which sets most CMake variables for the build |Very low|    

## 1. Programming models enabled for the TPLs
Edit the Part 1 of your `config_<machine's name>.cmake` file to specify the programming models which are considered in the TPLs. You can enable a GPU-accelerated code and a shared memory parallelization on the host.  To do so,  
1. set `ENABLE_OPENMP=ON`,   
2. for the GPU programming models - depending on the vendor, set   
- on Nvidia GPUs: ` ENABLE_CUDA=ON`   
- on AMD GPUs: `ENABLE_HIP=ON`   
- on Nvidia Grace-Hopper: `ENABLE_CUDA=ON ARM=ON`  

When installing the TPLs with a package manager, make sure to use some TPLs packages whose enabled programming models are consistent with those specified at Part 1 of the `config_<machine's name>.cmake`.   

## 2. Platform-specific parameters
Edit the Part 2 of your `config_<machine's name>.cmake` to set the platform-specific variables.  It emboddies various variables which are used to specify the compilers (for instance gcc, g++, gfortran, mpicc, nvcc, hipcc) and the compute options:   
- The root to the compilers for the programming model enabled on the host is typically `/usr/bin`  
- The variable `CUDA_ARCH` or `COMP_ARCH` is used to specify the architecture of the device. The related variables must be set in respect to the GPU programming model enabled and specified at [Part 1](##1.-programming models enabled for the TPLs). We refer to the following webpage for a mapping between various [GPU microarchitectures and their corresponding flags or compute capabilities](https://kokkos.org/kokkos-core-wiki/keywords.html#gpu-architectures)  
- Some compilation flags for the targetted architecture (`mtune` and `mcpu`) are set through the variable `CPU_TUNE_FLAG`.   

The host machine dependent variables and paths must be set and are only required when installing from source. This is done within the *conditional if* scopes on `BUILD_FROM_TPLMIRROR`. 
