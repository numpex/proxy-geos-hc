#ifndef DATATYPE_HPP_
#define DATATYPE_HPP_

#include "commonConfig.hpp"
#include <iostream>

#ifdef SEM_USE_VECTOR
  #include <vector>
  template<class T> class Vector 
  {
  public:
    Vector(int numRows) : data(numRows) {}
    Vector() : data(0) {}

    T & operator()(int index)
    {
        return data[index];
    }
    T & operator[](int index)
    {
        return data[index];
    }
    int size(){return this->size();}

  private:
    std::vector<T> data ;
  };

  template< class T > class Array2D
  {
  public:
    Array2D( int numRows, int numCols ): data( numRows, std::vector< T >( numCols,0 )) {}
    Array2D(): data( 0, std::vector< T >( 0 )) {}

    std::vector< T > & operator[]( int index ){return data[index];}
    T& operator()(size_t row, size_t col) {return data[row][col];}

  private:
    std::vector< std::vector< T > > data;
  };
  using vectorInt=std::vector< int >;
  using vectorReal=std::vector< float >;
  using vectorDouble=std::vector< double >;
  //using vectorInt=Vector< int >;
  //using vectorReal=Vector< float >;
  //using vectorDouble=Vector< double >;
  using arrayInt=Array2D< int >;
  using arrayReal=Array2D< float >;
  using arrayDouble=Array2D< double >;

  #ifdef SEM_USE_RAJA
    #include "RAJA/RAJA.hpp"
  #endif

  template<class T>
  T allocateVector(int n1)
  {
     std::cout<<"allocate vector of size "<<n1<<std::endl;
     T vect(n1);
    return vect;
  }
  template<class T>
  T allocateArray2D(int n1, int n2)
  {
    std::cout<<"allocate array of size "<<n1<<", "<<n2<<std::endl;
    T array(n1, n2);
    return array;
  }

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

  template<class T>
  T allocateVector(int n1)
  {
     std::cout<<"allocate vector of size "<<n1<<std::endl;
     T vect(n1);
    return vect;
  }
  template<class T>
  T allocateArray2D(int n1, int n2)
  {
    std::cout<<"allocate array of size "<<n1<<", "<<n2<<std::endl;
    T array(n1, n2);
    return array;
  }
#endif //SEM_USE_LVARRAY

#ifdef SEM_USE_KOKKOS
  #include <Kokkos_Core.hpp>
  //#define MemSpace Kokkos::CudaUVMSpace
  #define MemSpace Kokkos::HostSpace
  using ExecSpace = MemSpace::execution_space;
  using range_policy = Kokkos::RangePolicy<ExecSpace>;

  typedef Kokkos::View<int*,     Kokkos::LayoutLeft, MemSpace> vectorInt;
  typedef Kokkos::View<float*,   Kokkos::LayoutLeft, MemSpace> vectorReal;
  typedef Kokkos::View<double*,  Kokkos::LayoutLeft, MemSpace> vectorDouble;

  typedef Kokkos::View<int**,    Kokkos::LayoutLeft, MemSpace> arrayInt;
  typedef Kokkos::View<float**,  Kokkos::LayoutLeft, MemSpace> arrayReal;
  typedef Kokkos::View<double**, Kokkos::LayoutLeft, MemSpace> arrayDouble;

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
#endif //SEM_USE_KOKKOS

#endif //DATATYPE_HPP_
