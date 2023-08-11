#ifndef DATATYPE_HPP_
#define DATATYPE_HPP_

#include "commonConfig.hpp"
#include <iostream>

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
  using vectorInt=std::vector< int >;
  using vectorReal=std::vector< float >;
  using vectorDouble=std::vector< double >;
  using arrayInt=Array2D< int >;
  using arrayReal=Array2D< float >;
  using arrayDouble=Array2D< double >;

  #ifdef SEM_USE_RAJA
    #include "RAJA/RAJA.hpp"
  #endif

#endif //SEM_USE_VECTOR

#ifdef SEM_USE_LVARRAY

  #include "RAJA/RAJA.hpp"
  #include "Array.hpp"
  #ifdef SEM_USE_CUDA
    #include "ChaiBuffer.hpp"
    // Create an 1D array of integers.
    using vectorInt=LvArray::Array< int,
                                    1,
                                    camp::idx_seq< 0 >,
                                    std::ptrdiff_t,
                                    LvArray::ChaiBuffer >;
    using vectorReal=LvArray::Array< float,
                                     1,
                                     camp::idx_seq< 0 >,
                                     std::ptrdiff_t,
                                     LvArray::ChaiBuffer >;
    using vectorDouble=LvArray::Array< double,
                                       1,
                                       camp::idx_seq< 0 >,
                                       std::ptrdiff_t,
                                       LvArray::ChaiBuffer >;
    using arrayInt=LvArray::Array< int,
                                   2,
                                   camp::idx_seq< 0, 1 >,
                                   std::ptrdiff_t,
                                   LvArray::ChaiBuffer >;
    using arrayReal=LvArray::Array< float,
                                    2,
                                    camp::idx_seq< 0, 1 >,
                                    std::ptrdiff_t,
                                    LvArray::ChaiBuffer >;
    using arrayDouble=LvArray::Array< double,
                                      2,
                                      camp::idx_seq< 0, 1 >,
                                      std::ptrdiff_t,
                                      LvArray::ChaiBuffer >;
  #else
    #include "MallocBuffer.hpp"
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
  #endif
#endif //SEM_USE_LVARRAY

#ifdef SEM_USE_KOKKOS
  #include <Kokkos_Core.hpp>
  #ifdef SEM_USE_KOKKOSVECTOR
    #include <Kokkos_Vector.hpp>
    template< class T > class Array2D
    {
      public:
        Array2D( int numRows, int numCols ): data( numRows, Kokkos::vector< T >( numCols )) {}
        Array2D(): data( 0, Kokkos::vector< T >( 0 )) {}

        Kokkos::vector< T > & operator[]( int index )
      {
      return data[index];
      }

      private:
      Kokkos::vector< Kokkos::vector< T > > data;
    };
    using vectorInt=Kokkos::vector< int >;
    using vectorReal=Kokkos::vector< float >;
    using vectorDouble=Kokkos::vector< double >;
    using arrayInt=Array2D< int >;
    using arrayReal=Array2D< float >;
    using arrayDouble=Array2D< double >;
  #endif //SEM_USE_KOKKOSVECTOR
#endif //SEM_USE_KOKKOS

#endif //DATATYPE_HPP_
