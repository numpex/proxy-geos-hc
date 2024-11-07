## Our Objective: Why This README ?

This README aims to provide information for configuring the workstation of any new (Exa-DI) member.  
Let us say that our objective is to compile and run the Proxy-App provided by Henri CALANDRA [Proxy-Geos-HC](https://github.com/proxySem).  
This app is derived from [GEOSX](https://geosx-geosx.readthedocs-hosted.com/en/latest/) - an open-source multi-physics simulation code. 

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

The libraries and executable are built using CMake https://cmake.org/. There are multiple dependencies which are embedded in a cross-platform 
friendly manner using BLT.

#### BLT: a CMake-based Foundation for Building, Linking, and Testing Large-Scale HPC Applications 

For a broader overview, please refer to the [BLT webpage](https://github.com/LLNL/blt?tab=readme-ov-file). BLT is used (as submodule) embedded in a CMake folder to provide dependencies (libraries or executables). Referring to the previous webpage will help you cover the prerequisites and achieve the complete setup for smooth compilation and execution. 

Please check for the following requisites:  

- compilers: g++, gcc, clang  
- HPC programming: MPICH, OpenMP  
- Documentation: Doxygen, Sphinx 
- Other: Astyle, ClangFormat, cmake-format, Uncrustify, Yapf 
- Code quality clang-query, clang-tidy, Cppcheck 

## Specific Case: Ubuntu 24.04 LTS on a DELL Precision 5490 with an Intel CORE ULTRA 7 Processor  

1- apt-get install synaptic  
2- Installation of  plocate (useful to find the location of some objects, e.g gcc, nvcc)  
3- Git: git, git-lfs  
4- Installation of CMake  
5- gfortran (libgfortran)  
6- Some Libraries: intel-mkl to have blas, lapack, omp and other..  
7- Debugger: valgrind  
8- MPI: mpich  
9- Compiler clang and related std libraries  
10- OpenMP: libopenmpi  
11- Python: python3-h5py / python3-mpi4py  
12- xml2-utils  
13- cppcheck, tinyxml2  
14- uncrustify  
15- clang-tidy-18 / clang-tidy *(may require a symbolic link)*  
16- libastyle3, doxygen, Astyle  
17- python3-yapf  
18- cmake-format, clang-format *(may require a symbolic link)*  
19- sphinx, flex  
20- bison  
21- Specific version of compilers: gcc-11, g++-11  
22- libgtk-3-0  

### Setting Symbolic Links 

Once done, you may need to create symbolic link for some binaries that may not be accessible due to the default ```PATH``` setting or because their names include the version numbers.  

For example, create a symbolic link using: 

```ln -s path/to/original_target path/to/symbolic/link```  

This trick may concern a symbolic link to ```clang-tidy-18```, named ```clang-tidy```. 