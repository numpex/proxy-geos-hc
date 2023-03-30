#ifndef QKGL_HPP_
#define QKGL_HPP_

#include    <iostream>
#include    <vector>
#include    <cmath>

using    namespace std;

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
       vector<double> gaussLobattoQuadraturePoints(int order);

       // get Gauss Lobatto quadrature weights
       vector<double> gaussLobattoQuadratureWeights(int order);

       // compute  1d shape Functions and derivatives
       vector<double> shapeFunction1D(int order, double xi);
       vector<double> derivativeShapeFunction1D(int order, double xi);

       // get  1d shape Functions and derivatives for all quadrature points
       vector<vector<double>> getBasisFunction1D(int order,const vector<double> & quadraturePoints);
       vector<vector<double>> getDerivativeBasisFunction1D( int order,const vector<double> &quadraturePoints);

       // compute 2D gauss-lobatto weights
       vector<double> getGaussLobattoWeights(const vector<double> &quadraturePoints,
                                             const vector<double> &weights);

       // get  2d shape Functions  for all quadrature points
       vector<vector<double>> getBasisFunction2D(const vector<double> quadraturePoints,
                                                 const vector<vector<double>> &a,
                                                 const vector<vector<double>> &b);
       // compute Jacobian Matrix
       vector<vector<double>> computeJacobianMatrix(const int &nPointsPerElement,
                                                    const vector<vector<double>> &Xi ,
                                                    const  vector<vector<double>> &dxPhi ,
                                                    const  vector<vector<double>> &dyPhi );

       // compute determinant of Jacobian Matrix
       vector<double> computeDeterminantOfJacobianMatrix(const int & nPointsPerElement,
		                                         const vector<vector<double>> &jacobianMatrix);

       // compute inverse of Jacobian Matrix
       vector<vector<double>> computeInvJacobianMatrix(const int & nPointsPerElement,
                                                       const vector<vector<double>> &jacobianMatrix,
                                                       const vector<double> &detJ);

       // compute tranposed inverse of Jacobian Matrix
       vector<vector<double>> computeTranspInvJacobianMatrix(const int & nPointsPerElement,
                                                             const vector<vector<double>> &jacobianMatrix,
                                                             const vector<double> &detJ);

       // compute 𝐵 the matrix containing the geometrical informations
       vector<vector<double>> computeB(const int & nPoinsPerElement,
                                       const vector<vector<double>> &invJacobianMatrix,
                                       const vector<vector<double>> &transpInvJacobianMatrix,
                                       const vector<double> &detJ);
     
       // compute the matrix $R_{i,j}=\int_{K}{\nabla{\phi_i}.\nabla{\phi_j}dx}$
       vector<vector<double>> gradPhiGradPhi(const int & nPointsPerElement,
                                             const vector<double> &weights,
                                             const vector<vector<double>> &B,
                                             const vector<vector<double>> &dxPhi,
                                             const vector<vector<double>> &dyPhi);

       // compute the matrix $M_{i,j}=:w
       // \int_{K}{{\phi_i}.{\phi_j}dx}$
       vector<vector<double>> phiIphiJ(const int & nPointsPerElement,
                                             const vector<double> &weights,
                                             const vector<vector<double>> &phi,
                                             const vector<double> &detJ);

};

#endif //QKGL_HPP_

