  #ifndef QKGL_HPP_
  #define QKGL_HPP_

  #include    <iostream>
  #include    <vector>
  #include    <cmath>
  #include    "dataType.hpp"

//using    namespace std;

/**
 * This class is the basis class for the hexahedron finite element cells with
 * shape functions defined on Gauss-Lobatto quadrature points.
 */
namespace FE
{

  class QkGL
  {
  private:
    int order;
  public:

   #ifdef USE_RAJA
   LVARRAY_HOST_DEVICE QkGL();
   #elif defined USE_KOKKOS
   KOKKOS_FUNCTION QkGL();
   #else
   QkGL();
   #endif

  #ifdef USE_RAJA
  LVARRAY_HOST_DEVICE ~QkGL();
  #elif defined USE_KOKKOS
  KOKKOS_FUNCTION ~QkGL();
  #else
  ~QkGL();
  #endif

  // get Gauss Lobatto quadrature points
  #ifdef USE_RAJA
  void gaussLobattoQuadraturePoints( int order, vectorDouble const & quadraturePoints ) const;
  #elif defined USE_KOKKOS
  void gaussLobattoQuadraturePoints( int order, vectorDouble const & quadraturePoints ) const;
  #else
  void gaussLobattoQuadraturePoints( int order, vectorDouble & quadraturePoints ) const;
  #endif

  // get Gauss Lobatto quadrature weights
  #ifdef USE_RAJA
  void  gaussLobattoQuadratureWeights( int order, vectorDouble const & weights ) const;
  #elif defined USE_KOKKOS
  void  gaussLobattoQuadratureWeights( int order, vectorDouble const & weights ) const;
  #else
  void  gaussLobattoQuadratureWeights( int order, vectorDouble & weights ) const;
  #endif

  // compute  1d shape Functions and derivatives
  std::vector<double>  shapeFunction1D( int order, double xi ) const;
  std::vector<double>  derivativeShapeFunction1D( int order, double xi ) const;

  // get  1d shape Functions and derivatives for all quadrature points
  #ifdef USE_RAJA
  void  getBasisFunction1D( int order, vectorDouble const & quadraturePoints ,arrayDouble const & basisFunction1D) const;
  #elif defined USE_KOKKOS
  void  getBasisFunction1D( int order, vectorDouble const & quadraturePoints ,arrayDouble const & basisFunction1D) const;
  #else
  void  getBasisFunction1D( int order, vectorDouble & quadraturePoints ,arrayDouble & basisFunction1D) const;
  #endif

  #ifdef USE_RAJA
  void getDerivativeBasisFunction1D( int order, vectorDouble const & quadraturePoints, 
                                                arrayDouble const & derivativeBasisFunction1D ) const;
  #elif defined USE_KOKKOS
  void getDerivativeBasisFunction1D( int order, vectorDouble const & quadraturePoints, 
                                                arrayDouble const & derivativeBasisFunction1D ) const;
  #else
  void getDerivativeBasisFunction1D( int order, vectorDouble & quadraturePoints, 
                                                arrayDouble & derivativeBasisFunction1D ) const;
  #endif

  // compute 2D gauss-lobatto weights
  #ifdef USE_RAJA
  void getGaussLobattoWeights2D( vectorDouble const & quadraturePoints,
                                           vectorDouble const & weights,
                                           vectorDouble const & W )const;
  #elif defined USE_KOKKOS
  void getGaussLobattoWeights2D( vectorDouble const & quadraturePoints,
                                           vectorDouble const & weights,
                                           vectorDouble const & W )const;
  #else
  void getGaussLobattoWeights2D( vectorDouble & quadraturePoints,
                                           vectorDouble & weights,
                                           vectorDouble & W )const;
  #endif

  // compute 3D gauss-lobatto weights
  #ifdef USE_RAJA
  void getGaussLobattoWeights3D( vectorDouble const & quadraturePoints,
                                           vectorDouble const & weights,
                                           vectorDouble const & W )const;
  #elif defined USE_KOKKOS
  void getGaussLobattoWeights3D( vectorDouble const & quadraturePoints,
                                           vectorDouble const & weights,
                                           vectorDouble const & W )const;
  #else
  void getGaussLobattoWeights3D( vectorDouble & quadraturePoints,
                                           vectorDouble & weights,
                                           vectorDouble & W )const;
  #endif

  // get  2d shape Functions  for all quadrature points
  #ifdef USE_RAJA
  void getBasisFunction2D( vectorDouble const & quadraturePoints,
                           arrayDouble const & a,
                           arrayDouble const & b,
                           arrayDouble const & c) const;
  #elif defined USE_KOKKOS
  void getBasisFunction2D( vectorDouble const & quadraturePoints,
                           arrayDouble const & a,
                           arrayDouble const & b,
                           arrayDouble const & c) const;
  #else
  void getBasisFunction2D( vectorDouble & quadraturePoints,
                           arrayDouble & a,
                           arrayDouble & b,
                           arrayDouble & c) const;
  #endif

  //2D
  // compute B and M
  #ifdef USE_RAJA
  LVARRAY_HOST_DEVICE int computeB( const int & elementNumber,
                                    const int & order,
                                    arrayIntView     const & nodesList,
                                    arrayRealView    const & nodesCoords,
                                    vectorDoubleView const & weights2D,
                                    arrayDoubleView  const & dPhi,
                                    double massMatrixLocal[],
                                    double   B[][4] ) const;
  #elif defined USE_KOKKOS
  KOKKOS_FUNCTION int computeB( const int & elementNumber,
                                const int & order,
                                arrayInt     const & nodesList,
                                arrayReal    const & nodesCoords,
                                vectorDouble const & weights2D,
                                arrayDouble  const & dPhi,
                                double massMatrixLocal[],
                                double   B[][4]) const;
  #else
  int computeB( const int & elementNumber,
                const int & order,
                arrayInt  & nodesList,
                arrayReal & nodesCoords,
                vectorDouble & weights2D,
                arrayDouble  & dPhi,
                double massMatrixLocal[],
                double   B[][4]) const;
  #endif

  //3D
  // compute B and M
  #ifdef USE_RAJA
  LVARRAY_HOST_DEVICE int computeB( const int & elementNumber,
                                    const int & order,
                                    arrayIntView     const & nodesList,
                                    arrayRealView    const & nodesCoords,
                                    vectorDoubleView const & weights3D,
                                    arrayDoubleView  const & dPhi,
                                    double massMatrixLocal[],
                                    double   B[][6] ) const;
  #elif defined USE_KOKKOS
  KOKKOS_FUNCTION int computeB( const int & elementNumber,
                                const int & order,
                                arrayInt     const & nodesList,
                                arrayReal    const & nodesCoords,
                                vectorDouble const & weights3D,
                                arrayDouble  const & dPhi,
                                double massMatrixLocal[],
                                double   B[][6]) const;
  #else
  int computeB( const int & elementNumber,
                const int & order,
                arrayInt  & nodesList,
                arrayReal & nodesCoords,
                vectorDouble & weights3D,
                arrayDouble  & dPhi,
                double massMatrixLocal[],
                double   B[][6]) const;
  #endif

  // 2D
  // compute the matrix $R_{i,j}=\int_{K}{\nabla{\phi_i}.\nabla{\phi_j}dx}$
  // Marc Durufle Formulae
  #ifdef USE_RAJA
  LVARRAY_HOST_DEVICE int  gradPhiGradPhi( const int & nPointsPerElement,
                                           const int & order,
                                           vectorDoubleView const & weights2D,
                                           double const B[][4],
                                           arrayDoubleView const & dPhi,
                                           double  R[][36] ) const;
  #elif defined USE_KOKKOS
  KOKKOS_FUNCTION int  gradPhiGradPhi( const int & nPointsPerElement,
                                       const int & order,
                                       vectorDouble const & weights2D,
                                       double const B[][4],
                                       arrayDouble const & dPhi,
                                       double  R[][36] ) const;
  #else
  int  gradPhiGradPhi( const int & nPointsPerElement,
                       const int & order,
                       vectorDouble  & weights2D,
                       double const B[][4],
                       arrayDouble  & dPhi,
                       double  R[][36]) const;
  #endif

  // 3D
  // compute the matrix $R_{i,j}=\int_{K}{\nabla{\phi_i}.\nabla{\phi_j}dx}$
  // Marc Durufle Formulae
  #ifdef USE_RAJA
  LVARRAY_HOST_DEVICE int  gradPhiGradPhi( const int & nPointsPerElement,
                                           const int & order,
                                           vectorDoubleView const & weights3D,
                                           double const B[][6],
                                           arrayDoubleView const & dPhi,
                                           double  R[][8] ) const;
  #elif defined USE_KOKKOS
  KOKKOS_FUNCTION int  gradPhiGradPhi( const int & nPointsPerElement,
                                       const int & order,
                                       vectorDouble const & weights3D,
                                       double const B[][6],
                                       arrayDouble const & dPhi,
                                       double  R[][8] ) const;
  #else
  int  gradPhiGradPhi( const int & nPointsPerElement,
                       const int & order,
                       vectorDouble  & weights3D,
                       double const B[][6],
                       arrayDouble  & dPhi,
                       double  R[][8]) const;
  #endif

  // compute dx
  #ifdef USE_RAJA
  LVARRAY_HOST_DEVICE int  computeDs( const int & iFace,
                                      const int & order,
                                      arrayIntView const & faceInfos,
                                      int  numOfBasisFunctionOnFace[],
                                      float  Js[][6],
                                      arrayRealView   const & globalNodesCoords,
                                      arrayDoubleView const & derivativeBasisFunction2DX,
                                      arrayDoubleView const & derivativeBasisFunction2DY,
                                      float  ds[] ) const;
  #elif defined USE_KOKKOS
  KOKKOS_FUNCTION int  computeDs( const int & iFace,
                                  const int & order,
                                  arrayInt  const & faceInfos,
                                  int  numOfBasisFunctionOnFace[],
                                  float  Js[][6],
                                  arrayReal   const & globalNodesCoords,
                                  arrayDouble const & derivativeBasisFunction2DX,
                                  arrayDouble const & derivativeBasisFunction2DY,
                                  float  ds[]  ) const;
  #else
  int  computeDs( const int & iFace,
                  const int & order,
                  arrayInt  & faceInfos,
                  int  numOfBasisFunctionOnFace[],
                  float  Js[][6],
                  arrayReal    & globalNodesCoords,
                  arrayDouble  & derivativeBasisFunction2DX,
                  arrayDouble  & derivativeBasisFunction2DY,
                  float  ds[]  ) const;
  #endif
  };
}
#endif //QKGL_HPP_
