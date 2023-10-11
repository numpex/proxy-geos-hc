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

  QkGL();

#ifdef SEM_USE_RAJA
  ~QkGL();
#elif defined SEM_USE_KOKKOS
KOKKOS_FUNCTION ~QkGL();
#else
  ~QkGL();
#endif


  // get Gauss Lobatto quadrature points
  //vectorDouble gaussLobattoQuadraturePoints( int order ) const;
#ifdef SEM_USE_RAJA
  void gaussLobattoQuadraturePoints( int order, vectorDouble const & quadraturePoints ) const;
#elif defined SEM_USE_KOKKOS
  void gaussLobattoQuadraturePoints( int order, vectorDouble const & quadraturePoints ) const;
#else
  void gaussLobattoQuadraturePoints( int order, vectorDouble & quadraturePoints ) const;
#endif

  // get Gauss Lobatto quadrature weights
  //vectorDouble  gaussLobattoQuadratureWeights( int order ) const;
#ifdef SEM_USE_RAJA
  void  gaussLobattoQuadratureWeights( int order, vectorDouble const & weights ) const;
#elif defined SEM_USE_KOKKOS
  void  gaussLobattoQuadratureWeights( int order, vectorDouble const & weights ) const;
#else
  void  gaussLobattoQuadratureWeights( int order, vectorDouble & weights ) const;
#endif

  // compute  1d shape Functions and derivatives
  std::vector<double>  shapeFunction1D( int order, double xi ) const;
  std::vector<double>  derivativeShapeFunction1D( int order, double xi ) const;

  // get  1d shape Functions and derivatives for all quadrature points
  //arrayDouble  getBasisFunction1D( int order, vectorDouble const & quadraturePoints ) const;
#ifdef SEM_USE_RAJA
  void  getBasisFunction1D( int order, vectorDouble const & quadraturePoints ,arrayDouble const & basisFunction1D) const;
#elif defined SEM_USE_KOKKOS
  void  getBasisFunction1D( int order, vectorDouble const & quadraturePoints ,arrayDouble const & basisFunction1D) const;
#else
  void  getBasisFunction1D( int order, vectorDouble & quadraturePoints ,arrayDouble & basisFunction1D) const;
#endif

  //arrayDouble  getDerivativeBasisFunction1D( int order, vectorDouble const & quadraturePoints ) const;
#ifdef SEM_USE_RAJA
  void getDerivativeBasisFunction1D( int order, vectorDouble const & quadraturePoints, 
                                                arrayDouble const & derivativeBasisFunction1D ) const;
#elif defined SEM_USE_KOKKOS
  void getDerivativeBasisFunction1D( int order, vectorDouble const & quadraturePoints, 
                                                arrayDouble const & derivativeBasisFunction1D ) const;
#else
  void getDerivativeBasisFunction1D( int order, vectorDouble & quadraturePoints, 
                                                arrayDouble & derivativeBasisFunction1D ) const;
#endif

  // compute 2D gauss-lobatto weights
  //vectorDouble  getGaussLobattoWeights( vectorDouble & quadraturePoints,
  //                                      vectorDouble & weights ) const;
#ifdef SEM_USE_RAJA
  void getGaussLobattoWeights( vectorDouble const & quadraturePoints,
                                           vectorDouble const & weights,
                                           vectorDouble const & W )const;
#elif defined SEM_USE_KOKKOS
  void getGaussLobattoWeights( vectorDouble const & quadraturePoints,
                                           vectorDouble const & weights,
                                           vectorDouble const & W )const;
#else
  void getGaussLobattoWeights( vectorDouble & quadraturePoints,
                                           vectorDouble & weights,
                                           vectorDouble & W )const;
#endif

  // get  2d shape Functions  for all quadrature points
  //arrayDouble  getBasisFunction2D( vectorDouble & quadraturePoints,
  //                                 arrayDouble & a,
  //                                 arrayDouble & b ) const;

#ifdef SEM_USE_RAJA
  void getBasisFunction2D( vectorDouble const & quadraturePoints,
                           arrayDouble const & a,
                           arrayDouble const & b,
                           arrayDouble const & c) const;
#elif defined SEM_USE_KOKKOS
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
  // compute Jacobian Matrix
  //arrayDouble  computeJacobianMatrix( const int & nPointsPerElement,
  //                                    arrayDouble & Xi,
  //                                    arrayDouble & dxPhi,
  //                                    arrayDouble & dyPhi ) const;

#ifdef SEM_USE_RAJA
  int  computeJacobianMatrix( const int & nPointsPerElement,
                              double const Xi[][2],
                              arrayDouble const & dxPhi,
                              arrayDouble const & dyPhi,
                              double jacobianMatrix[][4] ) const;
#elif defined SEM_USE_KOKKOS
  KOKKOS_FUNCTION int  computeJacobianMatrix( const int & nPointsPerElement,
                                              double const Xi[][2],
                                              arrayDouble const & dxPhi,
                                              arrayDouble const & dyPhi,
                                              double jacobianMatrix[][4] ) const;
#else
  int  computeJacobianMatrix( const int & nPointsPerElement,
                              double const Xi[][2],
                              arrayDouble  & dxPhi,
                              arrayDouble  & dyPhi,
                              double jacobianMatrix[][4] ) const;

#endif

  // compute determinant of Jacobian Matrix
  //vectorDouble  computeDeterminantOfJacobianMatrix( const int & nPointsPerElement,
  //                                                  arrayDouble & jacobianMatrix ) const;
#ifdef SEM_USE_RAJA
  int  computeDeterminantOfJacobianMatrix( const int & nPointsPerElement,
                                           double  const  jacobianMatrix[][4],
                                           double  detJ[] ) const;
#elif defined SEM_USE_KOKKOS
  KOKKOS_FUNCTION int  computeDeterminantOfJacobianMatrix( const int & nPointsPerElement,
                                                           double const jacobianMatrix[][4],
                                                           double  detJ[] ) const;
#else
  int  computeDeterminantOfJacobianMatrix( const int & nPointsPerElement,
                                           double const  jacobianMatrix[][4],
                                           double detJ[] ) const;
#endif

  // compute inverse of Jacobian Matrix
  //arrayDouble  computeInvJacobianMatrix( const int & nPointsPerElement,
  //                                       arrayDouble & jacobianMatrix,
  //                                       vectorDouble & detJ ) const;
#ifdef SEM_USE_RAJA
  int  computeInvJacobianMatrix( const int & nPointsPerElement,
                                 double const  jacobianMatrix[][4],
                                 double const  detJ[],
                                 double  invJacobianMatrix[][4] ) const;
#elif defined SEM_USE_KOKKOS
  KOKKOS_FUNCTION int  computeInvJacobianMatrix( const int & nPointsPerElement,
                                                 double const  jacobianMatrix[][4],
                                                 double const  detJ[],
                                                 double  invJacobianMatrix[][4] ) const;
#else
  int  computeInvJacobianMatrix( const int & nPointsPerElement,
                                 double const  jacobianMatrix[][4],
                                 double const  detJ[],
                                 double  invJacobianMatrix[][4] ) const;
#endif

  // compute tranposed inverse of Jacobian Matrix
  //arrayDouble  computeTranspInvJacobianMatrix( const int & nPointsPerElement,
  //                                             arrayDouble & jacobianMatrix,
  //                                             vectorDouble & detJ ) const;
#ifdef SEM_USE_RAJA
  int  computeTranspInvJacobianMatrix( const int & nPointsPerElement,
				                               double const  jacobianMatrix[][4],
                                       double const  detJ[],
                                       double  transpInvJacobianMatrix[][4] ) const;
#elif defined SEM_USE_KOKKOS
  KOKKOS_FUNCTION int  computeTranspInvJacobianMatrix( const int & nPointsPerElement,
				                                               double const  jacobianMatrix[][4],
                                                       double const  detJ[],
                                                       double  invJacobianMatrix[][4] ) const;
#else
  int  computeTranspInvJacobianMatrix( const int & nPointsPerElement,
				                               double const  jacobianMatrix[][4],
                                       double const  detJ[],
                                       double  invJacobianMatrix[][4] ) const;
#endif

  // compute 𝐵 the matrix containing the geometrical informations
  //arrayDouble  computeB( const int & nPointsPerElement,
  //                       arrayDouble & invJacobianMatrix,
  //                       arrayDouble & transpInvJacobianMatrix,
  //                       vectorDouble & detJ ) const;
#ifdef SEM_USE_RAJA
  int  computeB( const int & nPointsPerElement,
                 double const invJacobianMatrix[][4],
                 double const transpInvJacobianMatrix[][4],
                 double const detJ[],
                 double B[][4] ) const;
#elif defined SEM_USE_KOKKOS
  KOKKOS_FUNCTION int  computeB( const int & nPointsPerElement,
                                 double const invJacobianMatrix[][4],
                                 double const transpInvJacobianMatrix[][4],
                                 double const detJ[],
                                 double B[][4] ) const;
#else
  int  computeB( const int & nPointsPerElement,
                 double const invJacobianMatrix[][4],
                 double const transpInvJacobianMatrix[][4],
                 double const detJ[],
                 double B[][4] ) const;
#endif

                
  // compute the matrix $R_{i,j}=\int_{K}{\nabla{\phi_i}.\nabla{\phi_j}dx}$
  // Marc Durufle Formulae
  //arrayDouble  gradPhiGradPhi( const int & nPointsPerElement,
  //                             const int & order,
  //                             vectorDouble & weights,
  //                             arrayDouble & B,
  //                             arrayDouble & dPhi ) const;
#ifdef SEM_USE_RAJA
  int  gradPhiGradPhi( const int & nPointsPerElement,
                       const int & order,
                       vectorDouble const & weights2D,
                       double const B[][4],
                       arrayDdouble const & dPhi,
                       double  R[][36] ) const;
#elif defined SEM_USE_KOKKOS
  KOKKOS_FUNCTION int  gradPhiGradPhi( const int & nPointsPerElement,
                                       const int & order,
                                       vectorDouble const & weights2D,
                                       double const B[][4],
                                       arrayDdouble const & dPhi,
                                       double  R[][36] ) const;
#else
  int  gradPhiGradPhi( const int & nPointsPerElement,
                       const int & order,
                       vectorDouble  & weights2D,
                       double const B[][4],
                       arrayDouble  & dPhi,
                       double  R[][36]) const;
#endif
  ///**
  // compute the matrix $R_{i,j}=\int_{K}{\nabla{\phi_i}.\nabla{\phi_j}dx}$
  //arrayDouble  gradPhiGradPhi( const int & nPointsPerElement,
  //                             vectorDouble & weights,
  //                             arrayDouble & B,
  //                             arrayDouble & dxPhi,
  //                             arrayDouble & dyPhi ) const;
#ifdef SEM_USE_RAJA
  int  gradPhiGradPhi( const int & nPointsPerElement,
                       vectorDouble const & weights,
                       double const B[][4],
                       arrayDouble const & dxPhi,
                       arrayDouble const & dyPhi,
                       double  R[][36] ) const;
#elif defined SEM_USE_KOKKOS
  KOKKOS_FUNCTION int  gradPhiGradPhi( const int & nPointsPerElement,
                                       vectorDouble const & weights,
                                       double const B[][4],
                                       arrayDouble const & dxPhi,
                                       arrayDouble const & dyPhi,
                                       double  R[][36] ) const;
#else
  int  gradPhiGradPhi( const int & nPointsPerElement,
                       vectorDouble const & weights,
                       double const B[][4],
                       arrayDouble  & dxPhi,
                       arrayDouble  & dyPhi,
                       double  R[][36] ) const;
#endif
  //**/

  // compute the matrix $M_{i,j}=\int_{K}{{\phi_i}.{\phi_j}dx}$ (optimized formulation)
  //vectorDouble  phiIphiJ( const int & nPointsPerElement,
  //                        vectorDouble & weights2D,
  //                        vectorDouble & detJ )const;
#ifdef SEM_USE_RAJA
  int  phiIphiJ( const int & nPointsPerElement,
                 vectorDouble const & weights2D,
                 double const  detJ[],
                 double  massMatrixLocal[] ) const;
#elif defined SEM_USE_KOKKOS
  KOKKOS_FUNCTION int  phiIphiJ( const int & nPointsPerElement,
                                 vectorDouble const & weights2D,
                                 double const  detJ[],
                                 double  massMatrixLocal[] ) const;
#else
  int  phiIphiJ( const int & nPointsPerElement,
                 vectorDouble  & weights2D,
                 double const  detJ[],
                 double  massMatrixLocal[] ) const;
#endif
  
  // compute dx
  //vectorReal computeDs( const int & iFace,
  //                      const int & order,
  //                      arrayInt & faceInfos,
  //                      arrayReal & globalNodesCoords,
  //                      arrayDouble & derivativeBasisFunction2DX,
  //                      arrayDouble & derivativeBasisFunction2DY ) const;
#ifdef SEM_USE_RAJA
  int  computeDs( const int & iFace,
                  const int & order,
                  arrayInt & faceInfos,
                  int  numOfBasisFunctionOnFace[],
                  float  Js[][6],
                  arrayReal   const & globalNodesCoords,
                  arrayDouble const & derivativeBasisFunction2DX,
                  arrayDouble const & derivativeBasisFunction2DY,
                  float  ds[] ) const;

#elif defined SEM_USE_KOKKOS
  KOKKOS_FUNCTION int  computeDs( const int & iFace,
                                  const int & order,
                                  arrayInt &  faceInfos,
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
