//************************************************************************
//  SEM proxy application v.0.0.1
//
//  main.cpp: this main file is simply a driver
//************************************************************************

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <Kokkos_Macros.hpp>
#include <Kokkos_Core.hpp>

#define MemSpace Kokkos::CudaUVMSpace
//#define MemSpace Kokkos::HostSpace


using ExecSpace = MemSpace::execution_space;
using range_policy = Kokkos::RangePolicy<ExecSpace>;

typedef Kokkos::View<float**, Kokkos::LayoutLeft, MemSpace>   arrayReal;
template<class T>
  T allocateVector(int n1)
  {
     std::cout<<"allocate vector of size "<<n1<<std::endl;
     T vect("v",n1);
    return vect;
  }
template<class T>
  T allocateArray2D(int n1, int n2)
  {
    std::cout<<"allocate array of size "<<n1<<", "<<n2<<std::endl;
    T array("a",n1, n2);
    return array;
  }
              

KOKKOS_INLINE_FUNCTION  int  compute(int const n1, int const n2, arrayReal const & arrayViewIn, arrayReal const & arrayViewOut)
{
  for (int i=0; i<n1;i++)
    {
      for ( int j=0; j<n2;j++)
      {
        arrayViewOut(i,j)=sqrt(2*i*j)*arrayViewIn(i,j);
      }     
    }
  return(0);
}

/*
inline  arrayReal   compute1(int const n1, int const n2, arrayReal const arrayViewIn)
{
  arrayReal array(n1,n2);
  for (int i=0; i<n1;i++)
    {
      for ( int j=0; j<n2;j++)
      {
        array(i,j)=sqrt(2*i*j)*arrayViewIn(i,j);
      }     
    }
  return(array);
}
*/

int main( int argc, char *argv[] )
{
  Kokkos::initialize(argc,argv);  
  const int nIter=10000;
  const int n1=1000;
  const int n2=1000;

  //arrayReal arrayViewOut("data",n1,n2);
  //arrayReal arrayViewIn("data",n1,n2);
  arrayReal arrayViewIn;
  arrayViewIn=allocateArray2D<arrayReal>(n1,n2);
  arrayReal arrayViewOut;
  arrayViewOut=allocateArray2D<arrayReal>(n1,n2);

  //arrayReal::HostMirror h_arrayViewIn=Kokkos::create_mirror_view(arrayViewIn);
  //arrayReal::HostMirror h_arrayViewOut=Kokkos::create_mirror_view(arrayViewOut);

  for (int i=0; i<n1;i++)
  {
    for ( int j=0; j<n2;j++)
    {
      arrayViewIn(i,j)=sqrt(2*i*j);
      arrayViewOut(i,j)=0;
    }     
  }

  //Kokkos::deep_copy(arrayViewIn, h_arrayViewIn);
  //Kokkos::deep_copy(arrayViewOut, h_arrayViewOut);
   // Timer products.
  Kokkos::Timer timer;
  std::cout<<"start loop1 \n";
  Kokkos::parallel_for("loop", nIter , KOKKOS_LAMBDA  ( int e )
    {   for (int i=0; i<n1;i++)
        {
           for ( int j=0; j<n2;j++)
           {
              arrayViewOut(i,j)=sqrt(2*i*j)*arrayViewIn(i,j);
           }
        }
   });
  double time=timer.seconds();
  std::cout << "Elapsed Time Kokkos Cuda loop : "<<time <<" seconds.\n"<<std::endl;
  auto h_a=Kokkos::create_mirror_view(arrayViewOut);
  Kokkos::deep_copy(arrayViewOut, h_a);
  std::cout << h_a(500,500)<<std::endl;


  // test call subroutine
  // -------------------------------
  // Timer products.
  Kokkos::Timer timer1;
  Kokkos::parallel_for("routine", nIter , KOKKOS_LAMBDA ( int e)
  {
     int i=compute(n1,n2,arrayViewIn,arrayViewOut) ;
  });
  double time1=timer1.seconds();
  std::cout << "Elapsed Time Cuda function call  loop : "<<time1 <<" seconds.\n"<<std::endl;
  auto h_b=Kokkos::create_mirror_view(arrayViewOut);
  Kokkos::deep_copy(arrayViewOut, h_b);
  std::cout << h_b(500,500)<<std::endl;

  Kokkos::finalize();
  return (0);
}
