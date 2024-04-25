#ifndef QKGL_HPP_
#define QKGL_HPP_

#include "dataType.hpp"
using namespace std;

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
  void gaussLobattoQuadraturePoints( int order, vectorDoubleView const & quadraturePoints ) const;

  // get Gauss Lobatto quadrature weights
  void  gaussLobattoQuadratureWeights( int order, vectorDoubleView const & weights ) const;

  // compute  1d shape Functions and derivatives
  vector<double>  shapeFunction1D( int order, double xi ) const;
  vector<double>  derivativeShapeFunction1D( int order, double xi ) const;

  // get  1d shape Functions and derivatives for all quadrature points
  void  getBasisFunction1D( int order, vectorDoubleView const & quadraturePoints ,arrayDoubleView const & basisFunction1D) const;
  void getDerivativeBasisFunction1D( int order, vectorDoubleView const & quadraturePoints, 
                                                arrayDoubleView const & derivativeBasisFunction1D ) const;

  // compute 2D gauss-lobatto weights
  void getGaussLobattoWeights2D( const int & order,
                                 vectorDoubleView const & weights,
                                 vectorDoubleView const & W )const;

  // compute 3D gauss-lobatto weights
  void getGaussLobattoWeights3D( const int & order,
                                 vectorDoubleView const & weights,
                                 vectorDoubleView const & W )const;

  // get  2d shape Functions  for all quadrature points
  void getBasisFunction2D( const int & order,
                           arrayDoubleView const & a,
                           arrayDoubleView const & b,
                           arrayDoubleView const & c) const;

  //2D
  // compute B and M
  PROXY_HOST_DEVICE int computeB( const int & elementNumber,
                                    const int & order,
                                    arrayIntView     const & nodesList,
                                    arrayRealView    const & nodesCoords,
                                    vectorDoubleView const & weights2D,
                                    arrayDoubleView  const & dPhi,
                                    float massMatrixLocal[],
                                    float   B[][4] ) const;
  //3D
  PROXY_HOST_DEVICE int computeB( const int & elementNumber,
                                    const int & order,
                                    arrayIntView     const & nodesList,
                                    arrayRealView    const & nodesCoords,
                                    vectorDoubleView const & weights3D,
                                    arrayDoubleView  const & dPhi,
                                    float massMatrixLocal[],
                                    float   B[][6] ) const;
  // 2D
  // compute the matrix $R_{i,j}=\int_{K}{\nabla{\phi_i}.\nabla{\phi_j}dx}$
  // Marc Durufle Formulae
  PROXY_HOST_DEVICE int  gradPhiGradPhi( const int & nPointsPerElement,
                                           const int & order,
                                           vectorDoubleView const & weights2D,
                                           arrayDoubleView const & dPhi,
                                           float const B[][4],
                                           float const pnLocal[],
                                           float R[],
                                           float Y[]) const;
  // 3D version
  // compute the matrix $R_{i,j}=\int_{K}{\nabla{\phi_i}.\nabla{\phi_j}dx}$
  // Marc Durufle Formulae
  PROXY_HOST_DEVICE int gradPhiGradPhi( const int & nPointsPerElement,
                                      const int & order,
                                      vectorDoubleView const & weights3D,
                                      arrayDoubleView const & dPhi,
                                      float const B[][6],
                                      float const pnLocal[],
                                      float R[],
                                      float Y[]) const;
  // compute dx
  PROXY_HOST_DEVICE int  computeDs( const int & iFace,
                                      const int & order,
                                      arrayIntView const & faceInfos,
                                      int  numOfBasisFunctionOnFace[],
                                      float  Js[][6],
                                      arrayRealView   const & globalNodesCoords,
                                      arrayDoubleView const & derivativeBasisFunction2DX,
                                      arrayDoubleView const & derivativeBasisFunction2DY,
                                      float  ds[] ) const;
  //computeDs
  PROXY_HOST_DEVICE int computeDs(  const int & iFace,
                                      const int & order,
                                      arrayIntView const & faceInfos,
                                      int  numOfBasisFunctionOnFace[],
                                      float  Js[][6],
                                      arrayRealView   const & globalNodesCoords,
                                      arrayDoubleView const & dPhi,
                                      float  ds[]  ) const;
  };
}
#endif //QKGL_HPP_
