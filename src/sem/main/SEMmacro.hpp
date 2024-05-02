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
  #define FENCE myPnGlobal=pnGlobal;
#endif

// create views only for RAJA
#if defined (USE_RAJA)
  #define CREATEVIEWS \
  arrayRealView rhsTerm=myRHSTerm.toView();\
  arrayRealView pnGlobal=myPnGlobal.toView();\
  vectorIntView rhsElement=myRhsElement.toView();\
  vectorRealView massMatrixGlobal=this->massMatrixGlobal.toView(); \
  vectorRealView yGlobal=this->yGlobal.toView(); \
  vectorRealView ShGlobal=this->ShGlobal.toView(); \
  arrayIntView globalNodesList=this->globalNodesList.toView(); \
  vectorRealView model=this->model.toView(); \
  vectorIntView listOfInteriorNodes=this->listOfInteriorNodes.toView(); \
  vectorIntView listOfBoundaryNodes=this->listOfBoundaryNodes.toView(); \
  arrayIntView faceInfos=this->faceInfos.toView(); \
  arrayIntView localFaceNodeToGlobalFaceNode=this->localFaceNodeToGlobalFaceNode.toView(); \
  vectorDoubleView weights=this->weights.toView(); \
  arrayRealView globalNodesCoords=this->globalNodesCoords.toView(); \
  arrayDoubleView derivativeBasisFunction1D=this->derivativeBasisFunction1D.toView(); 
#else
  #define CREATEVIEWS \
  arrayReal rhsTerm=myRHSTerm; \
  arrayReal pnGlobal=myPnGlobal; \
  vectorInt rhsElement=myRhsElement;
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
