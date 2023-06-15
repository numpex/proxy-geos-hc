## Welcome to the SEM proxy project!

The SEM proxy application is a benchmark designed to simulate wave propagation using the spectral element method (SEM), which is a Galerkin-based finite element method for solving partial differential equations. The benchmark is intended to be a standard tool for evaluating and comparing the performance of different high-performance computing (HPC) systems, particularly those used for scientific simulations.

## Why SEM proxy?

One of the key features of the SEM proxy benchmark is its adaptability to different programming models and HPC architectures. This makes it a useful proxy application for advancing the state of the art in high-performance computing. In addition to its technical capabilities, the SEM proxy benchmark is also designed to be easy to build and use. This makes it accessible to a wide range of users, from researchers to developers.

## What Programming Models users can select to run SEM proxy?

The programming models included in the current sem proxy implementations include:
* OpenMP [https://www.openmp.org/] to parallelizing for loops
* RAJA [https://raja.readthedocs.io/en/develop/]
* KOKKOS [https://kokkos.github.io/kokkos-core-wiki/]

## What data containers users can select to run SEM proxy?

The data containers included in the current sem proxy implementations include:
* LvArray [https://lvarray.readthedocs.io/en/latest/]
* C++ std::vector

## Quick Start to compile and install:

### Step 1: compile and install proxyAppSEM

```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=../install -DSEM_USE_VECTOR=ON ..  
   make install
```

The default compilation is sequential mode with std::vector implementation. 
So you will get an executable named "sem_Sequential_VECTOR.exe" in your installatin directory.

### Step 2: run the executable, for example:

```
   install/bin/sem_Sequential_VECTOR.exe
```

## Option: to utilize OMP + std::vector

### Step 1: compile and install proxyAppSEM

```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=../install -DSEM_USE_OMP=ON -DSEM_USE_VECTOR=ON ..  
   make install
```

### Step 2: run the executable, for example:

```
   install/bin/sem_OMP_VECTOR.exe
```

## Option: to utilize RAJA + std::vector

### Step 1: install RAJA:
```
   git clone --recursive https://github.com/llnl/raja.git
   cd raja
   mkdir build && cd build
   cmake -DCMAKE_INSTALL_PREFIX=../install -DRAJA_ENABLE_TESTS=Off -DENABLE_OPENMP=On ..
   make
   make install
```
 
### Step 2: setup your environment variables:

```
   export RAJA_DIR={your Raja installation directory}/lib/cmake/raja
```

### Step 3: compile and install proxyAppSEM

```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=../install -DSEM_USE_RAJA=ON -DSEM_USE_VECTOR=ON ..  
   make install
```

### Step 4: run the executable, for example:

```
   install/bin/sem_Raja_VECTOR.exe
```

## Option: to utilize KOKKOS + std::vector

### Step 1: install KOKKOS:

```
   git clone --recursive https://github.com/kokkos/kokkos.git
   cd kokkos
   mkdir build && cd build
   cmake -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_CXX_COMPILER=g++ -DKokkos_ENABLE_OPENMP=On -DKokkos_ENABLE_TESTS=Off ..
   make install
```

### Step 2: setup your environment variables:

```
   export KOKKOS_DIR={your kokkos installation directory}/lib64/cmake/Kokkos/
```

### Step 3: compile and install proxyAppSEM

```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=../install -DSEM_USE_KOKKOS=ON -DSEM_USE_VECTOR=ON ..  
   make install
```

### Step 4: run the executable, for example:

```
   install/bin/sem_Kokkos_VECTOR.exe
```
## Option: to utilize LvArray

In order to intall LvArray you need to install camp and raja (please refer to install RAJA section above) first. 

### Step 1: install camp and LvArray:

```
   git clone https://github.com/LLNL/camp.git
   cd camp
   mkdir build && cd build
   cmake -DCMAKE_INSTALL_PREFIX=<path to install location> ..
   make install
```
The simplest way to build LvArray is to define a host-configs file, an example could be find in host-configs/corigpu-gcc.cmake.
```
   git clone https://github.com/GEOS-DEV/LvArray.git
   cd LvArray
   python ./scripts/config-build.py -hc host-configs/<your configure file> -bt Release
   cd build-pecan-base-release
   make -j 32; make install
```
### Step 2: setup your environment variables:

```
   export LVARRAY_DIR={your LvArray installation directory}/share/lvarray/cmake/
```

### Step 3: compile and install proxyAppSEM

```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=../install -DSEM_USE_LVARRAY=ON ..  
   make install
```

### Step 4: run the executable, for example:

```
   install/bin/sem_Sequential_LVARRAY.exe
```

## Option: to utilize OMP + LvArray
### Step 1: install LvArray:
Please find installation LvArray in the above section.
### Step 2: setup your environment variables:

```
   export LVARRAY_DIR={your LvArray installation directory}/share/lvarray/cmake/
```

### Step 3: compile and install proxyAppSEM

```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=../install -DSEM_USE_OMP=ON -DSEM_USE_LVARRAY=ON ..  
   make install
```

### Step 4: run the executable, for example:

```
   install/bin/sem_OMP_LVARRAY.exe
```
## Option: to utilize RAJA + LvArray
### Step 1: install RAJA and LvArray:
Please find how to install RAJA and LvArray in the above sections.
### Step 2: setup your environment variables:

```
   export RAJA_DIR={your Raja installation directory}/lib/cmake/raja
   export LVARRAY_DIR={your LvArray installation directory}/share/lvarray/cmake/
```

### Step 3: compile and install proxyAppSEM

```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=../install -DSEM_USE_RAJA=ON -DSEM_USE_LVARRAY=ON ..  
   make install
```

### Step 4: run the executable, for example:

```
   install/bin/sem_RAJA_LVARRAY.exe
```

## Option: to utilize KOKKOS + LvArray
### Step 1: install KOKKOS and LvArray:
Please find how to install KOKKOS and LvArray in the above sections.
### Step 2: setup your environment variables:

```
   export KOKKOS_DIR={your kokkos installation directory}/lib64/cmake/Kokkos/
   export LVARRAY_DIR={your LvArray installation directory}/share/lvarray/cmake/
```

### Step 3: compile and install proxyAppSEM

```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=../install -DSEM_USE_KOKKOS=ON -DSEM_USE_LVARRAY=ON ..  
   make install
```

### Step 4: run the executable, for example:

```
   install/bin/sem_KOKKOS_LVARRAY.exe
```
## Option: to utilize Caliper as a profiler to output timing information

### Step 1: install CALIPER:

```
   git clone https://github.com/LLNL/Caliper.git
   cd Caliper
   mkdir build && cd build
   cmake -DCMAKE_INSTALL_PREFIX=<path to install location> ..
   make
   make install
```
### Step 2: setup your environment variables:
 
```
   export CALIPER_DIR={your caliper installation directory}/share/cmake/caliper
```

### Step 3: compile and install proxyAppSEM
```
   using any combinations of programming models and data containeres as above
```

### Step 4: run the executable, for example:

```
   CALI_CONFIG=runtime-report <path_to_bin>/<sem executable>
```
