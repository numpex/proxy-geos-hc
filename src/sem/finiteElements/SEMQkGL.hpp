#ifndef SEMQKGL_HPP_
#define SEMQKGL_HPP_

#include "dataType.hpp"
#include "SEMmacros.hpp"
#include "SEMdata.hpp"
using namespace std;

/**
 * This class is the basis class for the hexahedron finite element cells with shape functions defined on Gauss-Lobatto quadrature points.
 */
class SEMQkGL
{
private:
  int order;

public:
  PROXY_HOST_DEVICE SEMQkGL(){};
  PROXY_HOST_DEVICE ~SEMQkGL(){};

  // get Gauss Lobatto quadrature points
  void gaussLobattoQuadraturePoints( int order, vectorDouble const & quadraturePoints ) const;

  // get Gauss Lobatto quadrature weights
  void  gaussLobattoQuadratureWeights( int order, vectorDouble const & weights ) const; 
  // compute  1d shape Functions and derivatives
  vector< double >  shapeFunction1D( int order, double xi ) const;
  vector< double >  derivativeShapeFunction1D( int order, double xi ) const;

  // get  1d shape Functions and derivatives for all quadrature points
  void  getBasisFunction1D( int order, vectorDouble const & quadraturePoints, arrayDouble const & basisFunction1D ) const;
  void getDerivativeBasisFunction1D( int order, vectorDouble const & quadraturePoints,
                                     arrayDouble const & derivativeBasisFunction1D ) const;
  // compute B and M
  PROXY_HOST_DEVICE void computeB( const int & elementNumber,
                                   const int & order,
                                   VECTOR_DOUBLE_VIEW const & weights,
                                   ARRAY_INT_VIEW const & nodesList,
                                   ARRAY_REAL_VIEW const & nodesCoords,
                                   ARRAY_DOUBLE_VIEW const & dPhi,
                                   float massMatrixLocal[],
                                   float B[][COL] ) const;

  // compute the matrix $R_{i,j}=\int_{K}{\nabla{\phi_i}.\nabla{\phi_j}dx}$
  // Marc Durufle Formulae
  PROXY_HOST_DEVICE void gradPhiGradPhi( const int & nPointsPerElement,
                                         const int & order,
                                         VECTOR_DOUBLE_VIEW const & weights,
                                         ARRAY_DOUBLE_VIEW const & dPhi,
                                         float const B[][COL],
                                         float const pnLocal[],
                                         float R[],
                                         float Y[] ) const;
  //computeDs
  PROXY_HOST_DEVICE void computeDs( const int & iFace,
                                   const int & order,
                                   ARRAY_INT_VIEW const & faceInfos,
                                   int numOfBasisFunctionOnFace,
                                   float Js[][6],
                                   ARRAY_REAL_VIEW const & globalNodesCoords,
                                   ARRAY_DOUBLE_VIEW const & dPhi,
                                   float ds[] ) const;

  // compute mass matrix and stiffness matrices and then
  // compute global mass Matrix and global stiffness vector
  void computeMassAndStiffnessMatrices(
                                          const int & i1,
                                          const int & i2,
                                          SEMinfo & myInfo,
                                          ARRAY_INT_VIEW const & listOfElementsByColor,
                                          VECTOR_DOUBLE_VIEW const & weights,
                                          ARRAY_DOUBLE_VIEW const & derivativeBasisFunction1D,
                                          ARRAY_INT_VIEW const & globalNodesList,
                                          ARRAY_REAL_VIEW const & globalNodesCoords,
                                          VECTOR_REAL_VIEW const & massMatrixGlobal,
                                          VECTOR_REAL_VIEW const & yGlobal,
                                          VECTOR_REAL_VIEW const & model,
                                          ARRAY_REAL_VIEW const & pnGlobal) const;
};
#endif //SEMQKGL_HPP_
