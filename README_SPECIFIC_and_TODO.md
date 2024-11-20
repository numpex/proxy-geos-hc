## Why this document
This document aims to report some specific observations related to the deployment of the ProxyApp and things that can be improved.

### On Grid5000 

It occurs that the CMake version available in the build server can be outdated, as the minimum required version for the Proxy-app is 3.22.1.  This issue has been addressed by installing a recent CMake from a Git repository or using Guix.  


### Things to Improve 

We have encountered the issue of the unavailability of CUDA features on the Grid5000 frontend machine. Therefore, one needs to access the compute machine to utilize the CUDA features. The associated drawback is that, for now, the compilation is performed at the expense of compute resources. 
 
> This issue has been adressed by installing a local CUDA toolkit on the frontend machine. It proceeds as follows  

> > - downloading the ```.run``` file of the specific CUDA from the [Cuda-Toolkit-Archive](https://developer.nvidia.com/cuda-toolkit-archive) webpage   

>> - installing by running the ```.run``` executable file with the following options: ```--silent, --toolkit``` while specifying the ```--toolkitpath``` and the ```--defaultroot```.  

> By doing so, the next issue was related to [KOKKOS](https://kokkos.org/) which requires the GPU architecture, detected automatically by default. Using the option ```-Dkokkos_ARCH_{VALUE}=ON``` in the ```CMakeLists.txt``` building the Kokkos third-party libraries allows to proceed further without being on a GPU-capable computing machine. For instance, on a Nvidia Ampere architecture with "compute capability 8.0", the argument ```VALUE=AMPERE80```. Refer to [Kokkoss - GPU Architectures](https://kokkos.org/kokkos-core-wiki/keywords.html) for in-depth description.  

Another issue we observed is related to how to build executables or libraries targeting a specific and possibly different compute machine architecture.  
