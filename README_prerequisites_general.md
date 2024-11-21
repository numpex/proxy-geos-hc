## Why This README ?
  
This README aims to provide information with the prerequisites for compiling and running the [ProxyApp]((https://gitlab.inria.fr/numpex-pc5/wp2-co-design/proxy-geos-hc)) provided by Henri CALANDRA.  This App is derived from [GEOSX](https://geosx-geosx.readthedocs-hosted.com/en/latest/) - an open-source multi-physics simulation code. 

## Key Aspects and Prerequisites  

- Package Manager: Install your package manager of choice, such as "synaptic", if it's not already installed on the build server. For instance, use the command ```apt install synaptic```. 

- Check and Install the minimal prerequisites. Please refer to the following [GEOSX prerequisites webpage](https://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/buildGuide/Prerequisites.html#prerequisites)  to have an insight on some related prerequisites. These include, among others  
    - CMake build system generator  
    - build tools (GNU make or ninja on Linux, XCode on MacOS) 
    - python3 to run the python script    
    - Compatible C/C++ and Fortran compilers
    - CUDA toolkit providing nvcc compiler for GPU-accelerated Implementation  
    - MPI runtime and compilers  (if building with MPI).  
    - Git and Git-lfs. 


- #### Additionnal prerequisites for Building, Linking and testing based on CMake 
The libraries and executable are built in a cross-pwlatform friendly manner using BLT which provides a  [CMake](https://cmake.org/)-based Foundation for Building, Linking, and Testing Large-Scale HPC Applications. For a broader overview, please refer to the [BLT webpage](https://github.com/LLNL/blt?tab=readme-ov-file).  Referring to the aforementioned [BLT webpage](https://github.com/LLNL/blt?tab=readme-ov-file) will help you  complete the prerequisites setup for smooth compilation and execution. 

##  Summary of the prerequisites

|Build & debugger | Source code | Python & related| Compilers |Libraries |
|:---------------:|:---------------:|:---------------:|:---------------:|:---------------:|
|cmake, make, valgrind | git, git-lfs  |  python3, python3-h5py, python3-mpi4py, python3-yapf |clang, gcc, g++, gfortran (libgfortran)|blas, lapack, omp, libopenmpi, mpich, cuda|
|**code quality**  |**Code style**|**XML validation (GEOSX)**|**Documentation**|**Other**|
|clang-tidy, cppcheck, clang-query|libastyle3,cmake-format, clang-format , uncrustify| tinyxml2, xml2-utils |doxygen, sphinx, bison |flex, libgtk-3-0|

It is worth to mention that some installed binaries could be suffixed with their version number or located in folder not referenced in the `PATH` environment variable. Therefore some symbolic links might be useful.   
 The following command creates a symbolic link : 
```
ln -s path/to/original_target path/to/symbolic/link
```  
This trick may concern a symbolic link to `clang-tidy-18` with the name `clang-tidy`  and another for `clang-query`. 