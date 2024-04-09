#!/bin/sh

# to compile on cypress

source env/env_jie_Cypress_GPU.sh

mkdir build; cd build

# sem
rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install .. ; make 
rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install -DCMAKE_BUILD_TYPE=RelWithDebInfo .. ; make 

rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_OMP=ON .. ; make 
rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_RAJA=ON -DENABLE_CUDA=ON .. ; make 
rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_KOKKOS=ON -DENABLE_CUDA=ON .. ; make 

rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_RAJA=ON -DENABLE_HIP=ON .. ; make 
rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_KOKKOS=ON -DENABLE_HIP=ON .. ; make 

# fd

rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_OMP=ON .. ; make 
rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_RAJA=ON -DENABLE_CUDA=ON .. ; make 
rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_KOKKOS=ON -DENABLE_CUDA=ON .. ; make 

rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_RAJA=ON -DENABLE_HIP=ON .. ; make 
rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install -DUSE_KOKKOS=ON -DENABLE_HIP=ON .. ; make 
