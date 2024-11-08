## Our Objective: Why This README ?

This README aims to provide information with the prerequisites for compiling and running the Proxy-App provided by Henri CALANDRA [Proxy-Geos-HC](https://github.com/proxySem).  This app is derived from [GEOSX](https://geosx-geosx.readthedocs-hosted.com/en/latest/) - an open-source multi-physics simulation code. 

## Key Aspects and Prerequisites  

### Package Manager 

Install your package manager of choice, such as "synaptic", if it's not already installed on the laptop. For instance, use the command ```apt install synaptic```. 

### Check and Install the Minimal Prerequisites

Given our objective, we could refer to the following webpage to have an insight on the prerequisites to make things compile and run smoothly [Geosx prerequisites](https://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/buildGuide/Prerequisites.html#prerequisites)  

#### Nvidia Compiler for CUDA GPU-Accelerated Implementation. 

It may be natively installed on the laptop. Check with ```nvcc --version```. 


#### GCC/G++ compiler. 

check with ```gcc/g++ --version```. 


#### Gfortran compiler.

Most developments are written in C++, but some residual implementations in Fortran may require the gfortran compiler. 
Check with ```gfortran --version```.  

#### Install MPICH

It provides the necessary development kit for MPI programming.  

#### Git and Git-Lfs 

Git-Lfs is useful for managing large files. 


### Building, Linking, and testing based on CMake 

The libraries and executable are built using [CMake](https://cmake.org/). There are multiple dependencies which are embedded in a cross-platform 
friendly manner using BLT.

#### BLT: a CMake-based Foundation for Building, Linking, and Testing Large-Scale HPC Applications 

For a broader overview, please refer to the [BLT webpage](https://github.com/LLNL/blt?tab=readme-ov-file). BLT is used (as submodule) embedded in a CMake folder to provide dependencies (libraries or executables). Referring to the previous webpage will help you cover the prerequisites and achieve the complete setup for smooth compilation and execution. 

Please check for the following requisites:  

- compilers: g++, gcc, clang  
- HPC programming: MPICH, OpenMP  
- Documentation: Doxygen, Sphinx 
- Other: Astyle, ClangFormat, cmake-format, Uncrustify, Yapf 
- Code quality clang-query, clang-tidy, Cppcheck 