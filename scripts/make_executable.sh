#!/bin/sh

mkdir build; cd build

 rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install .. -DSEM_USE_OMP=ON -DSEM_USE_LVARRAY=ON; make; make install
 rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install .. -DSEM_USE_OMP=ON -DSEM_USE_VECTOR=ON; make; make install
 rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install .. -DSEM_USE_RAJA=ON -DSEM_USE_LVARRAY=ON; make; make install
 rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install .. -DSEM_USE_RAJA=ON -DSEM_USE_VECTOR=ON; make; make install
 rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install .. -DSEM_USE_KOKKOS=ON -DSEM_USE_LVARRAY=ON; make; make install
 rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install .. -DSEM_USE_KOKKOS=ON -DSEM_USE_VECTOR=ON; make; make install
 rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install .. -DSEM_USE_LVARRAY=ON; make; make install
 rm -rf ../build/* ; cmake -DCMAKE_INSTALL_PREFIX=../install .. -DSEM_USE_VECTOR=ON; make; make install
