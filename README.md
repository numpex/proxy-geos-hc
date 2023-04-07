# Welcome to the SEM proxy project!

The SEM proxy application is a benchmark designed to simulate wave propagation using the spectral element method (SEM), which is a Galerkin-based finite element method for solving partial differential equations. The benchmark is intended to be a standard tool for evaluating and comparing the performance of different high-performance computing (HPC) systems, particularly those used for scientific simulations.

# Why SEM proxy?

One of the key features of the SEM proxy benchmark is its adaptability to different programming models and HPC architectures. This makes it a useful proxy application for advancing the state of the art in high-performance computing. In addition to its technical capabilities, the SEM proxy benchmark is also designed to be easy to build and use. This makes it accessible to a wide range of users, from researchers to developers.


# Quick Start to compile and install:

## Step 1: in the root path of proxyAppSEM:

```
   edit config.h and set options:

        set (SOLVER "solver" CACHE PATH "" FORCE)
        options for solver are:
             sequentialVector
             ompVector
```

## Step 2: compile and install proxyAppSEM

```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=.. ..
   make install
```

## Step 3: run the executable, for example:

```
   <path_to_bin>/sem.exe
```


# Optional: to utilize CALIPER

## Step 1: install CALIPER:

```
   git clone https://github.com/LLNL/Caliper.git
   cd Caliper
   mkdir build && cd build
   cmake -DCMAKE_INSTALL_PREFIX=<path to install location> ..
   make
   make install
```

## Step 2: setup your environment variables:
 
###   *if you have the permission to change .bashrc, add the following in your .bashrc file:*

```
   caliper_install_dir=your caliper installation directory
   export CALIPER_DIR=$caliper_intall_dir/share/cmake/caliper
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$caliper_install_dir/lib
```
###   *if not, run the following commands:*

```
   export CALIPER_DIR={your caliper installation directory}/share/cmake/caliper
```

## Step 3: compile and install proxyAppSEM

```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=.. ..
   make install
```

## Step 4: run the executable, for example:

```
   CALI_CONFIG=runtime-report <path_to_bin>/sem.exe
```


# Optional: to utilize RAJA


## Step 1: install RAJA:

```
   git clone --recursive https://github.com/llnl/raja.git
   cd raja
   mkdir build && cd build
   cmake -DCMAKE_INSTALL_PREFIX=<path to install location> ..
   make
   make install
```

## Step 2: setup your environment variables:
 
###   *if you have the permission to change .bashrc, add the following in your .bashrc file:*

```
   raja_install_dir=your raja installation directory
   export RAJA_DIR=$raja_intall_dir/share/cmake/raja
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$raja_install_dir/lib
```

###   *if not, run the following commands:*

```
   export RAJA_DIR={your raja installation directory}/lib/cmake/raja
```

## Step 3: compile and install proxyAppSEM

```
   mkdir ./build
   cd build
   cmake -DCMAKE_INSTALL_PREFIX=.. ..
   make install
```

## Step 4: run the executable, for example:

```
   <path_to_bin>/sem.exe
```

