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

  //interior elements
  numberOfNodes=mesh.getNumberOfNodes();
  numberOfElements=mesh.getNumberOfElements();
  numberOfInteriorNodes=mesh.getNumberOfInteriorNodes();
  globalNodesList=mesh.globalNodesList( numberOfElements );
  listOfInteriorNodes=mesh.getListOfInteriorNodes( numberOfInteriorNodes );
  globalNodesCoords=mesh.nodesCoordinates( numberOfNodes );

  // boundary elements
  numberOfBoundaryNodes=mesh.getNumberOfBoundaryNodes();
  numberOfBoundaryFaces=mesh.getNumberOfBoundaryFaces();
  listOfBoundaryNodes=mesh.getListOfBoundaryNodes( numberOfBoundaryNodes );
  faceInfos=mesh.getBoundaryFacesInfos();
  localFaceNodeToGlobalFaceNode=mesh.getLocalFaceNodeToGlobalFaceNode();

  // get model
  model=mesh.getModel( numberOfElements );

  //allocate mesh arrays used in kernel
  numberOfPointsPerElement = ( order + 1 ) * ( order + 1 );
  localToGlobal=allocateVector<vectorInt>(numberOfPointsPerElement);
  Xi=allocateArray2D<arrayDouble>( numberOfPointsPerElement, 2 );

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

  jacobianMatrix=allocateArray2D<arrayDouble>( 4, numberOfPointsPerElement );

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
