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
  numberOfPointsPerElement = mesh.getNumberOfPointsPerElement();
  numberOfNodes=mesh.getNumberOfNodes();
  numberOfElements=mesh.getNumberOfElements();
  numberOfInteriorNodes=mesh.getNumberOfInteriorNodes();
  numberOfInteriorElements=mesh.getNumberOfInteriorElements();

  globalNodesList=allocateArray2D<arrayInt>(numberOfElements,numberOfPointsPerElement);
  mesh.globalNodesList( numberOfElements, globalNodesList );

  listOfInteriorNodes=allocateVector<vectorInt>(numberOfInteriorNodes);
  mesh.getListOfInteriorNodes( numberOfInteriorNodes, listOfInteriorNodes );

  listOfInteriorElements=allocateVector<vectorInt>(numberOfInteriorElements);
  mesh.getListOfInteriorElements( numberOfInteriorElements, listOfInteriorElements );

  globalNodesCoords=allocateArray2D<arrayReal>(numberOfNodes,2);
  mesh.nodesCoordinates( numberOfNodes, globalNodesCoords );

  // boundary elements
  numberOfBoundaryNodes=mesh.getNumberOfBoundaryNodes();
  numberOfBoundaryFaces=mesh.getNumberOfBoundaryFaces();

  listOfBoundaryNodes=allocateVector<vectorInt>(numberOfBoundaryNodes);
  mesh.getListOfBoundaryNodes( numberOfBoundaryNodes, listOfBoundaryNodes );
  
  faceInfos=allocateArray2D<arrayInt>(numberOfBoundaryFaces,2+(order+1));
  mesh.getBoundaryFacesInfos(numberOfBoundaryFaces, faceInfos);

  localFaceNodeToGlobalFaceNode=allocateArray2D<arrayInt>(numberOfBoundaryFaces,order+1);
  mesh.getLocalFaceNodeToGlobalFaceNode(numberOfBoundaryFaces, localFaceNodeToGlobalFaceNode);

  // get model
  model=mesh.getModel( numberOfElements );

  //allocate mesh arrays used in kernel
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

  jacobianMatrix=allocateArray2D<arrayDouble>(4, numberOfPointsPerElement);
  detJ=allocateVector<vectorDouble>(numberOfPointsPerElement);
  invJacobianMatrix=allocateArray2D<arrayDouble>(4, numberOfPointsPerElement);
  transpInvJacobianMatrix=allocateArray2D<arrayDouble>(4, numberOfPointsPerElement);
  B=allocateArray2D<arrayDouble>(4, numberOfPointsPerElement);
  R=allocateArray2D<arrayDouble>(numberOfPointsPerElement, numberOfPointsPerElement);
  massMatrixLocal=allocateVector<vectorDouble>(numberOfPointsPerElement);
  massMatrixGlobal=allocateVector<vectorReal>( numberOfNodes );
  yGlobal=allocateVector<vectorReal>( numberOfNodes );
  pnLocal=allocateVector<vectorReal>( numberOfPointsPerElement );
  Y=allocateVector<vectorReal>( numberOfPointsPerElement );

  ShGlobal=allocateVector<vectorReal>( numberOfBoundaryNodes );

  ds=allocateVector<vectorReal>( order+1 );
  Sh=allocateVector<vectorReal>( order+1 );
  numOfBasisFunctionOnFace=allocateVector<vectorInt>( order+1 );
  Js=allocateArray2D<arrayReal>( 2, order+1 );
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
  //static int numberOfNodes=mesh.getNumberOfNodes();
  //static int numberOfElements=mesh.getNumberOfElements();
  static vectorReal model=mesh.getModel( numberOfElements );
  //static arrayInt nodeList=mesh.globalNodesList( numberOfElements );
  int i, rhsElement;
  for( int i=0; i<numberOfRHS; i++ )
  {
    //extract element number for current rhs
    float x=rhsLocation(i,0);
    float y=rhsLocation(i,1);
    int rhsElement=mesh.getElementNumberFromPoints( x, y );
    // compute global node numbe to add source term to
    int nodeRHS=globalNodesList(rhsElement,0);
    pnGlobal(nodeRHS,i2)+=timeSample*timeSample*model[rhsElement]*model[rhsElement]*rhsTerm(i,timeStep);
  }
}
