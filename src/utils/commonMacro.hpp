#include "commonConfig.hpp"

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

// test SEM proxy 2D case 
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
