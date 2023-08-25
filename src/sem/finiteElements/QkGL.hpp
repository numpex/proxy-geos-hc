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

class QkGL
{
private:
  int order;
public:

  QkGL();
  ~QkGL();

  // get Gauss Lobatto quadrature points
  vectorDouble gaussLobattoQuadraturePoints( int order ) const;

  // get Gauss Lobatto quadrature weights
  vectorDouble  gaussLobattoQuadratureWeights( int order ) const;

  // compute  1d shape Functions and derivatives
  vectorDouble  shapeFunction1D( int order, double xi ) const;
  vectorDouble  derivativeShapeFunction1D( int order, double xi ) const;

  // get  1d shape Functions and derivatives for all quadrature points
  arrayDouble  getBasisFunction1D( int order, vectorDouble & quadraturePoints ) const;
  arrayDouble  getDerivativeBasisFunction1D( int order, vectorDouble & quadraturePoints ) const;

  // compute 2D gauss-lobatto weights
  vectorDouble  getGaussLobattoWeights( vectorDouble & quadraturePoints,
                                        vectorDouble & weights ) const;

  // get  2d shape Functions  for all quadrature points
  arrayDouble  getBasisFunction2D( vectorDouble & quadraturePoints,
                                   arrayDouble & a,
                                   arrayDouble & b ) const;
  // compute Jacobian Matrix
  //arrayDouble  computeJacobianMatrix( const int & nPointsPerElement,
  //                                    arrayDouble & Xi,
  //                                    arrayDouble & dxPhi,
  //                                    arrayDouble & dyPhi ) const;

  void computeJacobianMatrix( const int & nPointsPerElement,
                              arrayDouble & Xi,
                              arrayDouble & dxPhi,
                              arrayDouble & dyPhi,
                              arrayDouble & jacobianMatrix ) const;

  // compute determinant of Jacobian Matrix
  //vectorDouble  computeDeterminantOfJacobianMatrix( const int & nPointsPerElement,
  //                                                  arrayDouble & jacobianMatrix ) const;
  void computeDeterminantOfJacobianMatrix( const int & nPointsPerElement,
                                           arrayDouble & jacobianMatrix,
                                           vectorDouble & detJ ) const;

  // compute inverse of Jacobian Matrix
  //arrayDouble  computeInvJacobianMatrix( const int & nPointsPerElement,
  //                                       arrayDouble & jacobianMatrix,
  //                                       vectorDouble & detJ ) const;
  void computeInvJacobianMatrix( const int & nPointsPerElement,
                                         arrayDouble & jacobianMatrix,
                                         vectorDouble & detJ,
                                         arrayDouble & invJacobianMatrix ) const;

  // compute tranposed inverse of Jacobian Matrix
  //arrayDouble  computeTranspInvJacobianMatrix( const int & nPointsPerElement,
  //                                             arrayDouble & jacobianMatrix,
  //                                             vectorDouble & detJ ) const;
  void computeTranspInvJacobianMatrix( const int & nPointsPerElement,
                                       arrayDouble & jacobianMatrix,
                                       vectorDouble & detJ,
                                       arrayDouble  &transpInvJacobianMatrix ) const;

  // compute 𝐵 the matrix containing the geometrical informations
  //arrayDouble  computeB( const int & nPointsPerElement,
  //                       arrayDouble & invJacobianMatrix,
  //                       arrayDouble & transpInvJacobianMatrix,
  //                       vectorDouble & detJ ) const;
  void computeB( const int & nPointsPerElement,
                 arrayDouble & invJacobianMatrix,
                 arrayDouble & transpInvJacobianMatrix,
                 vectorDouble & detJ,
                 arrayDouble & B ) const;

                
  // compute the matrix $R_{i,j}=\int_{K}{\nabla{\phi_i}.\nabla{\phi_j}dx}$
  // Marc Durufle Formulae
  //arrayDouble  gradPhiGradPhi( const int & nPointsPerElement,
  //                             const int & order,
  //                             vectorDouble & weights,
  //                             arrayDouble & B,
  //                             arrayDouble & dPhi ) const;
  void gradPhiGradPhi( const int & nPointsPerElement,
                       const int & order,
                       vectorDouble & weights2D,
                       arrayDouble & B,
                       arrayDouble & dPhi,
                       arrayDouble & R ) const;
  ///**
  // compute the matrix $R_{i,j}=\int_{K}{\nabla{\phi_i}.\nabla{\phi_j}dx}$
  //arrayDouble  gradPhiGradPhi( const int & nPointsPerElement,
  //                             vectorDouble & weights,
  //                             arrayDouble & B,
  //                             arrayDouble & dxPhi,
  //                             arrayDouble & dyPhi ) const;
  void gradPhiGradPhi( const int & nPointsPerElement,
                       vectorDouble & weights,
                       arrayDouble & B,
                       arrayDouble & dxPhi,
                       arrayDouble & dyPhi,
                       arrayDouble & R ) const;
  //**/

  // compute the matrix $M_{i,j}=\int_{K}{{\phi_i}.{\phi_j}dx}$ (optimized formulation)
  //vectorDouble  phiIphiJ( const int & nPointsPerElement,
  //                        vectorDouble & weights2D,
  //                        vectorDouble & detJ )const;
  void phiIphiJ( const int & nPointsPerElement,
                 vectorDouble & weights2D,
                 vectorDouble & detJ,
                 vectorDouble & massMatrixLocal ) const;
  
  // compute the matrix $M_{i,j}= \int_{K}{{\phi_i}.{\phi_j}dx}$ ( non optimized formulation)
  //arrayDouble  phiIphiJ( const int & nPointsPerElement,
  //                       vectorDouble & weights2D,
  //                       arrayDouble & phi,
  //                       vectorDouble & detJ )const;
  void phiIphiJ( const int & nPointsPerElement,
                 vectorDouble & weights2D,
                 arrayDouble & phi,
                 vectorDouble & detJ,
                 arrayDouble  & massMatrixLocal )const;
  
  // compute dx
  //vectorReal computeDs( const int & iFace,
  //                      const int & order,
  //                      arrayInt & faceInfos,
  //                      arrayReal & globalNodesCoords,
  //                      arrayDouble & derivativeBasisFunction2DX,
  //                      arrayDouble & derivativeBasisFunction2DY ) const;
  void computeDs( const int & iFace,
                      const int & order,
                      arrayInt  & faceInfos,
                      vectorInt & numOfBasisFunctionOnFace,
                      arrayReal & Js,
                      arrayReal & globalNodesCoords,
                      arrayDouble & derivativeBasisFunction2DX,
                      arrayDouble & derivativeBasisFunction2DY,
                      vectorReal  & ds ) const;

};

#endif //QKGL_HPP_
