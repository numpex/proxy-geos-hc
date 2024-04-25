//************************************************************************
//   proxy application v.0.0.1
//
//  solverBase.cpp: simple 2D acoustive wave equation solver
//
//  the solverBase class servers as a base class for the SEM solver
//
//************************************************************************

#include "solverBase.hpp"

void solverBase::computeFEInit( const int & order, simpleMesh mesh, QkGL Qk)
{
/*

  // get infos from mesh
  //interior elements
  numberOfNodes=mesh.getNumberOfNodes();
  numberOfElements=mesh.getNumberOfElements();
  numberOfInteriorNodes=mesh.getNumberOfInteriorNodes();
  printf("numberOfNodes %d numberOfElements %d \n",numberOfNodes,numberOfElements);

  // number Of elements by color
  // sort element by color 
  int numberMaxOfElementsByColor=mesh.getNumberOfElementsByColor();
  listOfElementsByColor=allocateArray2D<arrayInt>(numberOfColors,numberMaxOfElementsByColor);
  mesh.sortElementsByColor(numberOfElementsByColor,listOfElementsByColor); 
  printf("number of elements color red %d\n",numberOfElementsByColor[0]);
  printf("number of elements color green %d\n",numberOfElementsByColor[1]);
  printf("number of elements color blue %d\n",numberOfElementsByColor[2]);
  printf("number of elements color yellow %d\n",numberOfElementsByColor[3]);
  
  numberOfPointsPerElement=mesh.getNumberOfPointsPerElement();
  globalNodesList=allocateArray2D<arrayInt>(numberOfElements,numberOfPointsPerElement);
  mesh.globalNodesList( numberOfElements, globalNodesList );

  listOfInteriorNodes=allocateVector<vectorInt>(numberOfInteriorNodes);
  mesh.getListOfInteriorNodes( numberOfInteriorNodes, listOfInteriorNodes );

  globalNodesCoords=allocateArray2D<arrayReal>(numberOfNodes,3);
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
  quadraturePoints=allocateVector<vectorDouble>(order+1);
  Qk.gaussLobattoQuadraturePoints( order, quadraturePoints );

  weights=allocateVector<vectorDouble>(order+1);
  Qk.gaussLobattoQuadratureWeights( order, weights );
  weights2D=allocateVector<vectorDouble>((order+1)*(order+1));
  Qk.getGaussLobattoWeights2D( quadraturePoints, weights, weights2D );
  weights3D=allocateVector<vectorDouble>((order+1)*(order+1)*(order+1));
  Qk.getGaussLobattoWeights3D( quadraturePoints, weights, weights3D );

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

  //shared arrays
  massMatrixGlobal=allocateVector<vectorReal>( numberOfNodes );
  yGlobal=allocateVector<vectorReal>( numberOfNodes );
  ShGlobal=allocateVector<vectorReal>( numberOfBoundaryNodes );
  std::cout<<"end of shared arrays initialization\n";
*/
}

