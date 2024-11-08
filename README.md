## Welcome to the  proxyApp project!

The proxyApp project  collects a suite of simple codes representing real applications.
It is intended to be a standard tool for evaluating and comparing the performance of different high-performance computing (HPC) systems, particularly those used for scientific simulations.


## Actual applications 

Current implementation of the proxyApp includes SEM ( Spectral finite Element Methods) and FD ( Finite Differences methods) to solve 2nd order acoustic wave equation in 2D and 3D spaces:
* The SEM proxy applicaton is a benchmark designed to simulate wave propagation using the spectral element method (SEM), which is a Galerkin-based finite element method for solving partial differential equations.
* The FD proxy applicaton is a benchmark designed to simulate wave propagation using finite differences stecnils operators for solving partial differential equations.

One of the key features of the SEM and FD proxy benchmarks are their adaptability to different programming models and HPC architectures. This makes them a useful proxy applications for advancing the state of the art in high-performance computing. In addition to their technical capabilities, they are designed to be easy to build and use. This makes them accessible to a wide range of users, from researchers to developers.

## What Programming Models users can select to run  proxy applications?

The programming models included in the current sem proxy implementations include:
* OpenMP [https://www.openmp.org/] to parallelizing for loops
* RAJA [https://raja.readthedocs.io/en/develop/]
* KOKKOS [https://kokkos.github.io/kokkos-core-wiki/]

RAJA and Kokkos must first be compiled and installed  in TPL4ProxyApp repository - described at **step 1**.

## What data containers users can select to run SEM proxy?

The data containers included in the current sem proxy implementations include:
* LvArray [https://lvarray.readthedocs.io/en/latest/]
* C++ std::vector

## Quick Start to compile and install:

First consider referring to the page on the prerequisites needed for [the general case](./README_preinstall.md); while a specific configuration on a [Ubuntu 24.04 LTS](./README_preinstallnew_1.md) is provided.

### Step 1: Building the third-party libraries

#### Environment Variables for the Compilers  

Start by exporting some environment variables to access the compilers (gcc, g++, mpicc, nvcc, gfortran) 
> A typical setting is provided in ```env_var.sh ```, which should be adapted and sourced (*source /path/to/env_var.sh*). 

Most of the compilers, except for nvcc, should typically be located in ```/usr/bin```.


#### Compile the Third-Party Libraries 

The environment variables set up are used in the ```config.cmake``` to configure the required CMake files.  

##### Preconfigure: Set the GPU Architecture  

>In the ```config.cmake```, set the GPU architecture variable ```CUDA_ARCH``` appropriately to enable CUDA. 

> -This setting must also be changed in the file ```BLTOptions.cmake```
>>located in the directories ```/path/to/proxy-geos-hc_tpl/cmake/blt/cmake/``` and ```/path/to/proxy-geos-hc/blt/cmake/```   

> -For consistency, the same setting must be used in the main code of the ProxyApp when building the executable, 

>> specifically in the ```CMakeLists.txt``` file at ```/path/to/proxy-geos-hc/src/CMakeLists.txt.```  

> We refer to the following webpage for a mapping between various GPU microarchitecture and their architecture flags or compute capabilities: [GPU microarchitecture and associated flags](https://kokkos.org/kokkos-core-wiki/keywords.html). 

It may be relevant to consider adding a suffix to the gcc/g++ compiler using the variable ```CC_VERSION``` - for instance, with ```set(CC_VERSION "-11")``` if the targetted compiler is gcc-11
  
##### Configure and Build  

1- From the folder of the third-party libraries, run the configuration script:  
>python3 ./scripts/config-build.py --hostconfig=./configs/config.cmake --buildpath=../buildTPL --installpath=../buildTPL/installTPL/ --buildtype=Release 

to generate the ```CMakeFile```/```Makefile``` in the build directory (hereafter labeled ```buildTPL```). 

Note that in this example the buildTPL folder and the third-party libraries folder are at the same tree-level.  
Moreover, the ```buildTPL``` folder contains a subfolder installTPL for the installation.  

2- From the build directory, make the build: 

> make  

[](*It is worth mentioning that it is possible to specify the number of processors while running the script by adding the argument ```-DNUM_PROC=xx```.  
This argument will be passed to the subsequent make command.*) 

### Step 2: compile and install proxyApp
Continue by exporting the path to the previously compiled TPLs required for linking and making the dependencies 
>This can be achieved by editing and sourcing the script ```env/env.sh``` in the PROXY-GEOSX-HC repository. 

1- Create a build folder from where the executable will be built and installed  
```
  mkdir buildProxyApp  
  cd buildProxyApp  
  cmake -DCMAKE_INSTALL_PREFIX=../install  CUDA_KOKKOS_RAJA_SETUP /path/to/PROXY-GEOSX-HC 
  make; make install
```

The ```CUDA_KOKKOS_RAJA_SETUP``` is a combination of the three options USE_KOKKOS, ENABLE_CUDA and ENABLE_OPENMP and is discussed below.  

### Step 3: run the executable, for example:

```
   install/bin/proxyName_SEQUENTIAL.exe ( with proxyName: sem or fd)
```

## Available configuration at step 2:
The default compilation is sequential mode with std::vector implementation. 
So you will get an executable named "sem_SEQUENTIAL.exe" and FDTDSEQUENTIAL.exe in your installation directory.

### OPEN_MP
in the case of OPENMP std::vector container is used.
```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_OMP=ON  ..  
   make; make install
```
### RAJA + OPEN_MP + CUDA (on Nvidia GPUs)
in the case of RAJA Lvarray container is used.
```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_RAJA=ON -DENABLE_CUDA=ON  ..  
   make; make install
 
```

### KOKKOS + OPEN_MP + CUDA (on Nvidia GPUs)
KOKKOS provides its own data container.
```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_KOKKOS=ON -DENABLE_CUDA=ON  ..  
   make install
```

### RAJA + OPEN_MP + HIP (on AMD GPUs)
```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_RAJA=ON -DENABLE_HIP=ON  ..  
   make; make install
```

### KOKKOS + OPEN_MP + HIP (on AMD GPUs)
```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_KOKKOS=ON -DENABLE_HIP=ON  ..  
   make install
 
```

### RAJA + OPEN_MP + CUDA + ARM (on Nvidia Grace-Hopper)
```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_RAJA=ON -DENABLE_CUDA=ON -DARM=ON ..  
   make; make install
```

### KOKKOS + OPEN_MP + CUDA + ARM (on Nvidia Grace-Hopper)
```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_KOKKOS=ON -DENABLE_CUDA=ON  -DARM=ON ..  
   make install
```
