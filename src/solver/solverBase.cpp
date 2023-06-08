//************************************************************************
//  SEM proxy application v.0.0.1
//
//  solverBase.cpp: simple 2D acoustive wave equation solver
//
//  the solverBase class servers as a base class for the SEM solver
//
//************************************************************************

#include "solverBase.hpp"

void solverBase::computeFEInit( const int & order,
                                simpleMesh mesh,
                                QkGL Qk )
{

  // get infos from mesh
  numberOfNodes=mesh.getNumberOfNodes();
  numberOfElements=mesh.getNumberOfElements();
  numberOfInteriorNodes=mesh.getNumberOfInteriorNodes();
  globalNodesList=mesh.globalNodesList( numberOfElements );
  listOfInteriorNodes=mesh.getListOfInteriorNodes( numberOfInteriorNodes );
  globalNodesCoords=mesh.nodesCoordinates( numberOfNodes );

  // get model
  model=mesh.getModel( numberOfElements );

  //get infos about finite element order of approximation
  numberOfPointsPerElement = ( order + 1 ) * ( order + 1 );

  // get quadrature points and weights
  quadraturePoints=Qk.gaussLobattoQuadraturePoints( order );
  weights=Qk.gaussLobattoQuadratureWeights( order );
  weights2D=Qk.getGaussLobattoWeights( quadraturePoints, weights );

  // get basis function and corresponding derivatives
  basisFunction1D=Qk.getBasisFunction1D( order, quadraturePoints );
  derivativeBasisFunction1D=Qk.getDerivativeBasisFunction1D( order, quadraturePoints );
  basisFunction2D=Qk.getBasisFunction2D( quadraturePoints, basisFunction1D, basisFunction1D );

  derivativeBasisFunction2DX=Qk.getBasisFunction2D( quadraturePoints, derivativeBasisFunction1D, basisFunction1D );
  derivativeBasisFunction2DY=Qk.getBasisFunction2D( quadraturePoints, basisFunction1D, derivativeBasisFunction1D );

}

// add right and side
void solverBase::addRightAndSides( const int & timeStep,
                                   const int & numberOfRHS,
                                   const int & i2,
                                   const float & timeSample,
                                   arrayReal & pnGlobal,
                                   arrayReal & rhsTerm,
                                   arrayReal & rhsLocation,
                                   simpleMesh mesh )
{
  static int numberOfNodes=mesh.getNumberOfNodes();
  static int numberOfElements=mesh.getNumberOfElements();
  static vectorReal model=mesh.getModel( numberOfElements );
  static arrayInt nodeList=mesh.globalNodesList( numberOfElements );
  int i, rhsElement;
  for( int i=0; i<numberOfRHS; i++ )
  {
    //extract element number for current rhs
    float x=rhsLocation[i][0];
    float y=rhsLocation[i][1];
    int rhsElement=mesh.getElementNumberFromPoints( x, y );
    // compute global node numbe to add source term to
    int nodeRHS=nodeList[rhsElement][0];
    pnGlobal[nodeRHS][i2]+=timeSample*timeSample*model[rhsElement]*model[rhsElement]*rhsTerm[i][timeStep];
  }
}
