#include "commonConfig.hpp"

// define Macros for function type
#if defined (USE_RAJA)
  #define PROXY_HOST_DEVICE LVARRAY_HOST_DEVICE
  #define SOLVER solverRaja
#elif defined (USE_KOKKOS)
  #define PROXY_HOST_DEVICE KOKKOS_FUNCTION 
  #define SOLVER solverKokkos
#elif defined (USE_OMP)
  #define PROXY_HOST_DEVICE
  #define SOLVER solverOMP
#else
  #define PROXY_HOST_DEVICE
  #define SOLVER solverSEQUENTIAL
#endif

// define Macros test SEM proxy 2D case 
#if defined (SEM2D)
  #define DIMENSION  2
  #define ROW 36
  #define COL 4
  #define ZEROED2D 0
// test SEM proxy 3D case
#else
  #define DIMENSION  3
  #define ROW 64
  #define COL 6
  #define ZEROED2D 1
#endif

// define Macros for loops
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

// create views only for RAJA
#if defined (USE_RAJA)
  #define CREATEVIEWS \
  vectorRealView massMatrixGlobal=h_massMatrixGlobal.toView(); \
  vectorRealView yGlobal=h_yGlobal.toView(); \
  vectorRealView ShGlobal=h_ShGlobal.toView(); \
  arrayIntView globalNodesList=h_globalNodesList.toView(); \
  vectorRealView model=h_model.toView(); \
  vectorIntView listOfInteriorNodes=h_listOfInteriorNodes.toView(); \
  vectorIntView listOfBoundaryNodes=h_listOfBoundaryNodes.toView(); \
  arrayIntView faceInfos=h_faceInfos.toView(); \
  arrayIntView localFaceNodeToGlobalFaceNode=h_localFaceNodeToGlobalFaceNode.toView(); \
  vectorDoubleView weights=h_weights.toView(); \
  arrayRealView globalNodesCoords=h_globalNodesCoords.toView(); \
  arrayDoubleView derivativeBasisFunction1D=h_derivativeBasisFunction1D.toView(); 
#else
  #define CREATEVIEWS
#endif
