//************************************************************************
//  SEM proxy application v.0.0.1
//
//  main.cpp: this main file is simply a driver
//************************************************************************

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>

#include "RAJA/RAJA.hpp"
#include "Macros.hpp"
#include "Array.hpp"
#include "MallocBuffer.hpp"
#include "ChaiBuffer.hpp"

// test lvArray
using arrayRealMB=LvArray::Array< float,
                                 2,
                                 camp::idx_seq< 1,0 >,
                                 std::ptrdiff_t,
                                 LvArray::MallocBuffer >;
using arrayRealViewMB=LvArray::ArrayView< float,
                                         2,
                                         0,
                                         std::ptrdiff_t,
                                         LvArray::MallocBuffer >;

using arrayRealCB=LvArray::Array< float,
                                 2,
                                 camp::idx_seq< 1,0 >,
                                 std::ptrdiff_t,
                                 LvArray::ChaiBuffer >;
using arrayRealViewCB=LvArray::ArrayView< float,
                                         2,
                                         0,
                                         std::ptrdiff_t,
                                         LvArray::ChaiBuffer >;

              
// test MallocBuffer MallocBuffer 
// -------------------------------
arrayRealMB computeMBMB(int const n1, int const n2,arrayRealMB arrayIn)
{
  arrayRealMB array(n1,n2);
  for (int i=0; i<n1;i++)
    {
      for ( int j=0; j<n2;j++)
      {
        array[i][j]=sqrt(2*i*j)*arrayIn[i][j];
      }     
    }
  return (array);  
}

// test ChaiBuffer ChaiBuffer 
// -------------------------------
arrayRealCB computeCBCB(int const n1, int const n2, arrayRealViewCB &arrayView)
{
  arrayRealCB array(n1,n2);
  for (int i=0; i<n1;i++)
    {
      for ( int j=0; j<n2;j++)
      {
        array[i][j]=sqrt(2*i*j)*arrayView[i][j];
      }     
    }
  return (array);  
}

// test MallocBuffer ChaiBuffer 
// -------------------------------
arrayRealMB computeMBCB(int const n1, int const n2, arrayRealViewCB &arrayView)
{
  arrayRealMB array(n1,n2);
  for (int i=0; i<n1;i++)
    {
      for ( int j=0; j<n2;j++)
      {
        array[i][j]=sqrt(2*i*j)*arrayView[i][j];
      }     
    }
  return (array);  
}


int main( int argc, char *argv[] )
{
  
  
  const int nIter=1000000;
  const int n1=100;
  const int n2=100;

  arrayRealCB arrayToMove(n1,n2);
  arrayRealMB arrayIn(n1,n2);
  for (int i=0; i<n1;i++)
  {
    for ( int j=0; j<n2;j++)
    {
      arrayToMove[i][j]=sqrt(2*i*j);
      arrayIn[i][j]=sqrt(2*i*j);
    }     
  }
  arrayRealViewCB arrayView=arrayToMove.toView();

  // test MallocBuffer MallocBuffer OMP
  // -------------------------------
  std::chrono::time_point< std::chrono::system_clock > startTimeMB = std::chrono::system_clock::now();
  #pragma omp parallel for
  for (int e=0; e<nIter; e++)
  {
      arrayRealMB array=computeMBMB(n1,n2,arrayIn) ; 
      if (e==nIter-1)
        std::cout <<"Result OMP MBMB "<<array[n1/2][n2/2]<<std::endl;   
  }
  std::cout << "Elapsed Time OMP loop MBMB : "<<std::chrono::duration_cast< std::chrono::milliseconds >
            ( std::chrono::system_clock::now() - startTimeMB ).count() / 1000.0 <<" seconds.\n"<<std::endl;
  
  // test MallocBuffer ChaiBuffer OMP
  // -------------------------------
  std::chrono::time_point< std::chrono::system_clock > startTimeMBCB = std::chrono::system_clock::now();
  #pragma omp parallel for
  for (int e=0; e<nIter; e++)
  {
      arrayRealMB array=computeMBCB(n1,n2,arrayView) ;
      if (e==nIter-1)
        std::cout <<"Result OMP MBCB "<<array[n1/2][n2/2]<<std::endl;
  }
  std::cout << "Elapsed Time OMP loop MBCB : "<<std::chrono::duration_cast< std::chrono::milliseconds >
          ( std::chrono::system_clock::now() - startTimeMBCB ).count() / 1000.0 <<" seconds.\n"<<std::endl;

  // test ChaiBuffer ChaiBuffer OMP
  // -------------------------------
  std::chrono::time_point< std::chrono::system_clock > startTimeCB = std::chrono::system_clock::now();
  #pragma omp parallel for
  for (int e=0; e<nIter; e++)
  {
      arrayRealCB array=computeCBCB(n1,n2,arrayView) ;
      if (e==nIter-1)
        std::cout <<"Result OMP CBCB "<<array[n1/2][n2/2]<<std::endl;
  }
  std::cout << "Elapsed Time OMP loop CBCB : "<<std::chrono::duration_cast< std::chrono::milliseconds >
          ( std::chrono::system_clock::now() - startTimeCB ).count() / 1000.0 <<" seconds.\n"<<std::endl;


  // test MallocBuffer MallocBuffer RAJA
  // -------------------------------
  std::chrono::time_point< std::chrono::system_clock > startTimeRAJAMBMB = std::chrono::system_clock::now();
  RAJA::forall< RAJA::omp_parallel_for_exec >( RAJA::RangeSegment( 0, nIter ), [=,&arrayIn] ( int e )
  {
      arrayRealMB array=computeMBMB(n1,n2,arrayIn) ;
      if (e==nIter-1)
        std::cout <<"Result RAJA MBMB "<<array[n1/2][n2/2]<<std::endl;
  });
  std::cout << "Elapsed Time RAJA loop MBCB : "<<std::chrono::duration_cast< std::chrono::milliseconds >
          ( std::chrono::system_clock::now() - startTimeRAJAMBMB ).count() / 1000.0 <<" seconds.\n"<<std::endl;

  // test MallocBuffer ChaiBuffer RAJA
  // -------------------------------
  std::chrono::time_point< std::chrono::system_clock > startTimeRAJAMBCB = std::chrono::system_clock::now();
  RAJA::forall< RAJA::omp_parallel_for_exec >( RAJA::RangeSegment( 0, nIter ), [=,&arrayView] ( int e )
  {
      arrayRealMB array=computeMBCB(n1,n2,arrayView) ;
      if (e==nIter-1)
        std::cout <<"Result RAJA MBCB "<<array[n1/2][n2/2]<<std::endl;
  });
  std::cout << "Elapsed Time RAJA loop MBCB : "<<std::chrono::duration_cast< std::chrono::milliseconds >
          ( std::chrono::system_clock::now() - startTimeRAJAMBCB ).count() / 1000.0 <<" seconds.\n"<<std::endl;

  // test ChaiBuffer ChaiBuffer RAJA
  // -------------------------------
  std::chrono::time_point< std::chrono::system_clock > startTimeRAJACBCB = std::chrono::system_clock::now();
  RAJA::forall< RAJA::omp_parallel_for_exec >( RAJA::RangeSegment( 0, nIter ), [=,&arrayView] ( int e )
  {
      arrayRealCB array=computeCBCB(n1,n2,arrayView) ;
      if (e==nIter-1)
        std::cout <<"Result RAJA CBCB "<<array[n1/2][n2/2]<<std::endl;
  });
  std::cout << "Elapsed Time RAJA loop CBCB : "<<std::chrono::duration_cast< std::chrono::milliseconds >
          ( std::chrono::system_clock::now() - startTimeRAJACBCB ).count() / 1000.0 <<" seconds.\n"<<std::endl;
}
