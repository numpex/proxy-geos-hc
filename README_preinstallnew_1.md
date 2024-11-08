## Specific Case of a Ubuntu 24.04 laptop  

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