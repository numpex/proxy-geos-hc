#ifndef DATATYPE_HPP_
#define DATATYPE_HPP_

#include "commonConfig.hpp"
#include <iostream>

using namespace std;

#ifdef SEM_USE_VECTOR
#include <vector>
template< class T > class Array2D
{
public:
  Array2D( int numRows, int numCols ): data( numRows, std::vector< T >( numCols )) {}
  Array2D(): data( 0, std::vector< T >( 0 )) {}

  std::vector< T > & operator[]( int index )
  {
    return data[index];
  }

private:
  std::vector< std::vector< T > > data;
};
using vectorInt=vector< int >;
using vectorReal=vector< float >;
using vectorDouble=vector< double >;
using arrayInt=Array2D< int >;
using arrayReal=Array2D< float >;
using arrayDouble=Array2D< double >;

#ifdef SEM_USE_RAJA
#include "RAJA/RAJA.hpp"
#endif

#endif //SEM_USE_VECTOR

#ifdef SEM_USE_LVARRAY

#include "RAJA/RAJA.hpp"
#include "LvArray/Array.hpp"
#include "LvArray/MallocBuffer.hpp"

// Create an 1D array of integers.
using vectorInt=LvArray::Array< int,
                                1,
                                camp::idx_seq< 0 >,
                                std::ptrdiff_t,
                                LvArray::MallocBuffer >;
using vectorReal=LvArray::Array< float,
                                 1,
                                 camp::idx_seq< 0 >,
                                 std::ptrdiff_t,
                                 LvArray::MallocBuffer >;
using vectorDouble=LvArray::Array< double,
                                   1,
                                   camp::idx_seq< 0 >,
                                   std::ptrdiff_t,
                                   LvArray::MallocBuffer >;
using arrayInt=LvArray::Array< int,
                               2,
                               camp::idx_seq< 0, 1 >,
                               std::ptrdiff_t,
                               LvArray::MallocBuffer >;
using arrayReal=LvArray::Array< float,
                                2,
                                camp::idx_seq< 0, 1 >,
                                std::ptrdiff_t,
                                LvArray::MallocBuffer >;
using arrayDouble=LvArray::Array< double,
                                  2,
                                  camp::idx_seq< 0, 1 >,
                                  std::ptrdiff_t,
                                  LvArray::MallocBuffer >;


#endif //SEM_USE_LVARRAY

#ifdef SEM_USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

#endif //DATATYPE_HPP_
