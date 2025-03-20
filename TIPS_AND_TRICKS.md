## Tips and tricks

We report below some specific observations related to the deployment of the [ProxyApp](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc/) - and the associated [third-party libraries](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl), suggestions addressing the issues. 

### Unavailable library
It occurs that a product or a required specific version can be unavailable. This issue can be addressed by installing the product, either from source or using a package manager (Spack or Guix).   

A commonplace recommandation before any installing  is related to consider using [Environment Modules](https://modules.sourceforge.net/) package. 
#### Loading module
Indeed, on some machines, it is possible to modify dynamically the user environment throughout some modulefiles which may have been installed by the system administrator. First check if the `<product>`  is available on the machine by running
```
module avail <product>
```
and then load it
```
module load <product>
```
#### Installing from source 
Here is a typical case, on a given Grid5000 frontend machine where the CMake is out-of-date as the minimum required version for building the ProxyApp is 3.22.1. This has been addressed by installing from source.  

Another case involved the unavailability of CUDA features preventing from [getting installed the third-party libraries]((https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl) ). One way to address this issue involves installing a local CUDA Toolkit on the frontend machine. It proceeds as follows  

- downloading the specific CUDA release from the [Cuda-Toolkit-Archive](https://developer.nvidia.com/cuda-toolkit-archive) webpage   
 - installing by running the ```.run``` executable file with the following options: ```--silent, --toolkit``` while specifying the ```--toolkitpath``` and the ```--defaultroot```.  
 
By doing so, the next issue was related to [KOKKOS](https://kokkos.org/) which had failed to detect automatically the GPU architecture since the CUDA driver is not available. 
Using the option `-Dkokkos_ARCH_{VALUE}=ON` in the `CMakeLists.txt` during the building stage has allowed to proceed further without being on a GPU-capable computing machine.   For instance, on a Nvidia Ampere architecture with "compute capability 8.0", the argument must be set as ```VALUE=AMPERE80```.   Please refer to [Kokkoss - GPU Architectures](https://kokkos.org/kokkos-core-wiki/keywords.html) for in-depth description for other architectures.  

### Development issues under AMD architectures

Here are some issues encountered during the deployment of the [ProxyApp](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc/) and the linked [third-party libraries](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl) on AMD GPUs, especially on Grid5000:

- The [ProxyApp](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc/) source code uses some recent features of [KOKKOS](https://kokkos.org/) (version &ge; 4.0) but the installed version of ROCm and the HIP compiler is too old and [KOKKOS](https://kokkos.org/) can't support it.
- The build of the [TPLs](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc_tpl) using [RAJA](https://raja.readthedocs.io/en/develop/#) also needs a more recent version of ROCm.
- The build using the flag ```ENABLE_OPENMP=ON``` needs an additional package called `openmp-extras` that is included in the ROCm environment. Nevertheless, the installed ROCm version on Grid5000 (4.5.0) doesn't have this package which can only be used as a module. The version provided in the module (5.2) is different from the ROCm installation, leading to errors (even on a simple hello-world test including HIP and OpenMP).

These issues could be overcome by deploying a more recent version of ROCm. On Grid5000, this would be possible by deploying a more recent OS (currently Debian 11 and ROCm needs Debian 12). This could be done by using the kadeploy solution. This solution has not been explored yet.

Currently the best solution to be able to build and run the [ProxyApp](https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc/) on AMD architecture is to use the solution provided by guix which clearly simplify the deployment.
