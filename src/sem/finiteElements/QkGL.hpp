  #ifndef QKGL_HPP_
  #define QKGL_HPP_

  #include    "dataType.hpp"
  #include    "commonMacro.hpp"
  #include    <cmath>

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
    PROXY_HOST_DEVICE QkGL(){};
    PROXY_HOST_DEVICE ~QkGL(){};

  // get Gauss Lobatto quadrature points
  void gaussLobattoQuadraturePoints( int order, vectorDouble & quadraturePoints ) const;

  // get Gauss Lobatto quadrature weights
  void  gaussLobattoQuadratureWeights( int order, vectorDouble & weights ) const;

  // compute  1d shape Functions and derivatives
  std::vector<double>  shapeFunction1D( int order, double xi ) const;
  std::vector<double>  derivativeShapeFunction1D( int order, double xi ) const;

  // get  1d shape Functions and derivatives for all quadrature points
  void  getBasisFunction1D( int order, vectorDouble const & quadraturePoints ,arrayDouble & basisFunction1D) const;
  void getDerivativeBasisFunction1D( int order, vectorDouble const & quadraturePoints, 
                                                arrayDouble & derivativeBasisFunction1D ) const;

  // compute 2D gauss-lobatto weights
  void getGaussLobattoWeights2D( vectorDouble const & quadraturePoints,
                                           vectorDouble const & weights,
                                           vectorDouble & W )const;

  // compute 3D gauss-lobatto weights
  void getGaussLobattoWeights3D( vectorDouble const & quadraturePoints,
                                           vectorDouble const & weights,
                                           vectorDouble & W )const;

  // get  2d shape Functions  for all quadrature points
  void getBasisFunction2D( vectorDouble const & quadraturePoints,
                           arrayDouble & a,
                           arrayDouble & b,
                           arrayDouble & c) const;

  //2D
  // compute B and M
  PROXY_HOST_DEVICE int computeB( const int & elementNumber,
                                const int & order,
  #if defined(USE_RAJA) || defined(USE_KOKKOS)
                                arrayIntView     const & nodesList,
                                arrayRealView    const & nodesCoords,
                                vectorDoubleView const & weights2D,
                                arrayDoubleView  const & dPhi,
  #else
                                arrayIntView     & nodesList,
                                arrayRealView    & nodesCoords,
                                vectorDoubleView & weights2D,
                                arrayDoubleView  & dPhi,
  #endif
                                float massMatrixLocal[],
                                float   B[][4]) const;

  //3D
  PROXY_HOST_DEVICE int computeB( const int & elementNumber,
                                const int & order,
  #if defined(USE_RAJA) || defined(USE_KOKKOS)
                                arrayIntView     const & nodesList,
                                arrayRealView    const & nodesCoords,
                                vectorDoubleView const & weights3D,
                                arrayDoubleView  const & dPhi,
  #else
                                arrayIntView     & nodesList,
                                arrayRealView    & nodesCoords,
                                vectorDoubleView & weights3D,
                                arrayDoubleView  & dPhi,
  #endif
                                float massMatrixLocal[],
                                float   B[][6]) const;

  // 2D
  // compute the matrix $R_{i,j}=\int_{K}{\nabla{\phi_i}.\nabla{\phi_j}dx}$
  // Marc Durufle Formulae
  PROXY_HOST_DEVICE int  gradPhiGradPhi( const int & nPointsPerElement,
                                       const int & order,
  #if defined(USE_RAJA) || defined(USE_KOKKOS)
                                       vectorDoubleView const & weights2D,
                                       arrayDoubleView const & dPhi,
  #else
                                       vectorDoubleView & weights2D,
                                       arrayDoubleView & dPhi,
  #endif
                                       float const B[][4],
                                       float const pnLocal[],
                                       float R[],
                                       float Y[]) const;

  // 3D version
  // compute the matrix $R_{i,j}=\int_{K}{\nabla{\phi_i}.\nabla{\phi_j}dx}$
  // Marc Durufle Formulae
  PROXY_HOST_DEVICE int gradPhiGradPhi( const int & nPointsPerElement,
                                      const int & order,
  #if defined(USE_RAJA) || defined(USE_KOKKOS)
                                      vectorDoubleView const & weights3D,
                                      arrayDoubleView const & dPhi,
  #else
                                      vectorDoubleView & weights3D,
                                      arrayDoubleView & dPhi,
  #endif
                                      float const B[][6],
                                      float const pnLocal[],
                                      float R[],
                                      float Y[]) const;

  // compute dx
  PROXY_HOST_DEVICE int  computeDs( const int & iFace,
                                  const int & order,
                                  arrayIntView  & faceInfos,
                                  int  numOfBasisFunctionOnFace[],
                                  float  Js[][6],
                                  arrayRealView   & globalNodesCoords,
                                  arrayDoubleView & derivativeBasisFunction2DX,
                                  arrayDoubleView & derivativeBasisFunction2DY,
                                  float  ds[]  ) const;
  };
}
#endif //QKGL_HPP_
