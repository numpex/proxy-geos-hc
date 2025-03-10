## Key Aspects and Prerequisites  for ProxyApp
  
We provide below with information on the prerequisites for compiling and running the [ProxyApp]((https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc)).  This ProxyApp is derived from [GEOSX](https://geosx-geosx.readthedocs-hosted.com/en/latest/) - an open-source multi-physics simulation code.  

Please consider the aspects below:

- Install your package manager of choice, such as "synaptic", if it's not already installed on the build server. For instance, use the command ```apt install synaptic```. 

- Check and Install the minimal prerequisites. Please refer to the following [GEOSX prerequisites webpage](https://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/buildGuide/Prerequisites.html#prerequisites)  to have an insight on some related prerequisites. These include, among others  
    - CMake build system generator  
    - build tools (GNU make or ninja on Linux, XCode on MacOS) 
    - python3 to run the python script    
    - Compatible C/C++ and Fortran compilers
    - CUDA toolkit providing nvcc compiler for GPU-accelerated implementation  
    - MPI runtime and compilers  (if building with MPI).  
    - Git and Git-lfs. 


- Additionnal prerequisites for Building, Linking and Testing (BLT) -based on CMake    
The libraries and executable are built in a cross-platform friendly manner using BLT which provides a  [CMake](https://cmake.org/)-based foundation for building, linking, and testing large-scale HPC applications.  
For a broader overview, please refer to the [BLT webpage](https://github.com/LLNL/blt?tab=readme-ov-file).  Doing so will help in the completeness of the prerequisites setup needed for a smooth compilation and execution. 

##  Summary of the prerequisites

- #### For building and compiling
|Build & debugger | Source code | Python & related                                     | Compilers |Libraries |
|:---------------:|:---------------:|:-------------------------------------------------------:|:---------------------------------------:|:---------------:|
|cmake, make, valgrind | git, git-lfs  |  python3, python3-h5py, python3-mpi4py, python3-yapf |clang, gcc, g++, gfortran (libgfortran)|blas, lapack, omp, libopenmpi, mpich, cuda|

- #### Quality check, documentation and other
|code quality| Code style | XML validation (GEOSX)| Documentation |Other |
|:---------------:|:-------------------------------:|:---------------:|:---------------:|:---------------:|
|clang-tidy, cppcheck, clang-query|libastyle3,cmake-format, clang-format , uncrustify| tinyxml2, xml2-utils |doxygen, sphinx, bison |flex, libgtk-3-0|

It is worth to mention that some installed binaries could be suffixed with their version number or located in folder not referenced in the `PATH` environment variable. Therefore some symbolic links might be useful.   

This trick may concern a symbolic link to `clang-tidy-18` with the name `clang-tidy`  and another for `clang-query`. For the former,  the symbolic link is typically made running
```
ln -s /lib/llvm-18/bin/clang-tidy   /usr/bin/clang-tidy
```  