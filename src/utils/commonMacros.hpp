#include "commonConfig.hpp"

// define Macros for function type
#if defined (USE_RAJA)
  #define PROXY_HOST_DEVICE LVARRAY_HOST_DEVICE
#elif defined (USE_KOKKOS)
  #define PROXY_HOST_DEVICE KOKKOS_FUNCTION 
#else
  #define PROXY_HOST_DEVICE
#endif

#if defined (USE_RAJA)
  #define LOOPHEAD(Range, Iterator)\
    RAJA::forall< deviceExecPolicy >( RAJA::RangeSegment( 0, Range),  [=] LVARRAY_HOST_DEVICE  ( int Iterator ) {
  #define LOOPEND   });
#elif defined (USE_KOKKOS)
  #define LOOPHEAD(Range, Iterator) Kokkos::parallel_for( Range, KOKKOS_CLASS_LAMBDA ( const int Iterator ){
  #define LOOPEND   });
#elif defined (USE_OMP)
  #define LOOPHEAD(Range, Iterator)\
    _Pragma("omp parallel for")\
    for( int Iterator=0; Iterator<Range; Iterator++ ){
  #define LOOPEND   }
// sequential case
#else 
  #define LOOPHEAD(Range, Iterator)\
    for( int Iterator=0; Iterator<Range; Iterator++ ){
  #define LOOPEND   }
#endif

// define atomic add operation
#if defined (USE_RAJA)
  #define ATOMICADD(ADD1,ADD2) RAJA::atomicAdd< deviceAtomicPolicy >(&ADD1,ADD2);
  #define FENCE\
  if(timeStep%100==0)\
     RAJA::forall< RAJA::seq_exec>( RAJA::RangeSegment( 0, myMeshinfo.numberOfNodes ), [pnGlobal] LVARRAY_HOST_DEVICE ( int i ){});
#elif defined (USE_KOKKOS)
  #define ATOMICADD(ADD1,ADD2) Kokkos::atomic_add(&ADD1,ADD2)
  #define FENCE Kokkos::fence();
#else
  #define ATOMICADD(ADD1,ADD2) ADD1+=ADD2
  #define FENCE 
#endif

#if defined (USE_RAJA)
  #define ARRAY_DOUBLE_VIEW arrayDoubleView 
  #define ARRAY_REAL_VIEW arrayRealView
  #define ARRAY_INT_VIEW arrayIntView
  #define VECTOR_DOUBLE_VIEW vectorDoubleView
  #define VECTOR_INT_VIEW vectorIntView
#else
  #define ARRAY_DOUBLE_VIEW arrayDouble
  #define ARRAY_REAL_VIEW arrayReal
  #define ARRAY_INT_VIEW arrayInt
  #define VECTOR_DOUBLE_VIEW vectorDouble
  #define VECTOR_INT_VIEW vectorInt
#endif

#if defined (USE_KOKKOS)
  #define KOKKOSNAME "v",
#else
  #define KOKKOSNAME 
#endif
