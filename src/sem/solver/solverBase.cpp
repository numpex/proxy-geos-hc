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

  globalNodesList=allocateArray2D<arrayInt>(numberOfElements,(order+1)*(order+1));
  mesh.globalNodesList( numberOfElements, globalNodesList );

  listOfInteriorNodes=allocateVector<vectorInt>(numberOfInteriorNodes);
  mesh.getListOfInteriorNodes( numberOfInteriorNodes, listOfInteriorNodes );

  globalNodesCoords=allocateArray2D<arrayReal>(numberOfNodes,2);
  mesh.nodesCoordinates( numberOfNodes, globalNodesCoords );

  // boundary elements
  numberOfBoundaryNodes=mesh.getNumberOfBoundaryNodes();
  numberOfBoundaryFaces=mesh.getNumberOfBoundaryFaces();

  listOfBoundaryNodes=allocateVector<vectorInt>(numberOfBoundaryNodes);
  mesh.getListOfBoundaryNodes( numberOfBoundaryNodes, listOfBoundaryNodes );

  faceInfos=allocateArray2D<arrayInt>(numberOfBoundaryFaces,2+(order+1));
  mesh.getBoundaryFacesInfos(faceInfos);

  localFaceNodeToGlobalFaceNode=allocateArray2D<arrayInt>(numberOfBoundaryFaces, order+1 );
  mesh.getLocalFaceNodeToGlobalFaceNode(localFaceNodeToGlobalFaceNode);

  // get model
  model=allocateVector<vectorReal>(numberOfElements);
  mesh.getModel( numberOfElements, model );


  // get quadrature points and weights
  numberOfPointsPerElement = ( order + 1 ) * ( order + 1 );
  quadraturePoints=allocateVector<vectorDouble>(order+1);
  Qk.gaussLobattoQuadraturePoints( order, quadraturePoints );

  weights=allocateVector<vectorDouble>(order+1);
  Qk.gaussLobattoQuadratureWeights( order, weights );
  weights2D=allocateVector<vectorDouble>(numberOfPointsPerElement);
  Qk.getGaussLobattoWeights( quadraturePoints, weights, weights2D );

  // get basis function and corresponding derivatives
  basisFunction1D=allocateArray2D<arrayDouble>(order+1,order+1);
  Qk.getBasisFunction1D( order, quadraturePoints,basisFunction1D );

  derivativeBasisFunction1D=allocateArray2D<arrayDouble>(order+1,order+1);
  Qk.getDerivativeBasisFunction1D( order, quadraturePoints, derivativeBasisFunction1D );

  int nBasisFunctions=(order+1)*(order+1); 
  basisFunction2D=allocateArray2D<arrayDouble>(nBasisFunctions,nBasisFunctions);
  Qk.getBasisFunction2D( quadraturePoints, basisFunction1D, basisFunction1D, basisFunction2D );

  derivativeBasisFunction2DX=allocateArray2D<arrayDouble>(nBasisFunctions,nBasisFunctions);
  Qk.getBasisFunction2D( quadraturePoints, derivativeBasisFunction1D, basisFunction1D, derivativeBasisFunction2DX );
  
  derivativeBasisFunction2DY=allocateArray2D<arrayDouble>(nBasisFunctions,nBasisFunctions);
  Qk.getBasisFunction2D( quadraturePoints, basisFunction1D, derivativeBasisFunction1D, derivativeBasisFunction2DY );

  cout<<"fin allocation  shared:"<<endl;


  //private arrays
#ifndef SEM_USE_OMP
  localToGlobal=allocateArray2D<arrayInt>(numberOfThreads,numberOfPointsPerElement);
  Xi=allocateArray3D<array3DDouble>( numberOfThreads, numberOfPointsPerElement, 2 );

  jacobianMatrix=allocateArray3D<array3DDouble>(numberOfThreads,4, numberOfPointsPerElement);
  detJ=allocateArray2D<arrayDouble>(numberOfThreads,numberOfPointsPerElement);
  invJacobianMatrix=allocateArray3D<array3DDouble>(numberOfThreads,4, numberOfPointsPerElement);
  transpInvJacobianMatrix=allocateArray3D<array3DDouble>(numberOfThreads,4, numberOfPointsPerElement);

  B=allocateArray3D<array3DDouble>(numberOfThreads,4, numberOfPointsPerElement);
  R=allocateArray3D<array3DDouble>(numberOfThreads,numberOfPointsPerElement, numberOfPointsPerElement);

  massMatrixLocal=allocateArray2D<arrayDouble>(numberOfThreads,numberOfPointsPerElement);
  

  pnLocal=allocateArray2D<arrayReal>(numberOfThreads, numberOfPointsPerElement );
  Y=allocateArray2D<arrayReal>(numberOfThreads, numberOfPointsPerElement );

  ds=allocateArray2D<arrayReal>(numberOfThreads, order+1 );
  Sh=allocateArray2D<arrayReal>(numberOfThreads, order+1 );
  numOfBasisFunctionOnFace=allocateArray2D<arrayInt>(numberOfThreads, order+1 );
  Js=allocateArray3D<array3DReal>(numberOfThreads, 2, order+1 );
  cout<<"end of private arrays allocation\n";
#endif


  //shared arrays
  massMatrixGlobal=allocateVector<vectorReal>( numberOfNodes );
  yGlobal=allocateVector<vectorReal>( numberOfNodes );
  ShGlobal=allocateVector<vectorReal>( numberOfBoundaryNodes );
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
  mesh.getModel( numberOfElements, model );
  mesh.globalNodesList( numberOfElements, globalNodesList );
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
