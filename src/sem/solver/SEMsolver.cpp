//************************************************************************
//   proxy application v.0.0.1
//
//  SEMsolver.cpp: simple 2D acoustive wave equation solver
//
//  the SEMsolver class servers as a base class for the SEM solver
//
//************************************************************************

#include "SEMsolver.hpp"

void SEMsolver::computeFEInit( const int & order, SEMmesh mesh )
{

  // get infos from mesh
  //interior elements
  numberOfNodes=mesh.getNumberOfNodes();
  numberOfElements=mesh.getNumberOfElements();
  numberOfInteriorNodes=mesh.getNumberOfInteriorNodes();
  printf("numberOfNodes %d numberOfElements %d \n",numberOfNodes,numberOfElements);

  // number Of elements by color
  // sort element by color 
  int numberMaxOfElementsByColor=mesh.getNumberOfElementsByColor();
  listOfElementsByColor=allocateArray2D<arrayIntView>(numberOfColors,numberMaxOfElementsByColor);
  mesh.sortElementsByColor(numberOfElementsByColor,listOfElementsByColor); 
  printf("number of elements color red %d\n",numberOfElementsByColor[0]);
  printf("number of elements color green %d\n",numberOfElementsByColor[1]);
  printf("number of elements color blue %d\n",numberOfElementsByColor[2]);
  printf("number of elements color yellow %d\n",numberOfElementsByColor[3]);
  
  numberOfPointsPerElement=mesh.getNumberOfPointsPerElement();

  globalNodesList=allocateArray2D<arrayIntView>(numberOfElements,numberOfPointsPerElement);
  mesh.globalNodesList( numberOfElements, globalNodesList );

  listOfInteriorNodes=allocateVector<vectorIntView>(numberOfInteriorNodes);
  mesh.getListOfInteriorNodes( numberOfInteriorNodes, listOfInteriorNodes );

  globalNodesCoords=allocateArray2D<arrayRealView>(numberOfNodes,3);
  mesh.nodesCoordinates( numberOfNodes, globalNodesCoords );

  // boundary elements
  numberOfBoundaryNodes=mesh.getNumberOfBoundaryNodes();
  numberOfBoundaryFaces=mesh.getNumberOfBoundaryFaces();

  listOfBoundaryNodes=allocateVector<vectorIntView>(numberOfBoundaryNodes);
  mesh.getListOfBoundaryNodes( numberOfBoundaryNodes, listOfBoundaryNodes );

  faceInfos=allocateArray2D<arrayIntView>(numberOfBoundaryFaces,2+(order+1));
  mesh.getBoundaryFacesInfos(faceInfos);

  localFaceNodeToGlobalFaceNode=allocateArray2D<arrayIntView>(numberOfBoundaryFaces, order+1 );
  mesh.getLocalFaceNodeToGlobalFaceNode(localFaceNodeToGlobalFaceNode);

  // get model
  model=allocateVector<vectorRealView>(numberOfElements);
  mesh.getModel( numberOfElements, model );


  // get quadrature points and weights
  quadraturePoints=allocateVector<vectorDoubleView>(order+1);
  myQk.gaussLobattoQuadraturePoints( order, quadraturePoints );

  weights=allocateVector<vectorDoubleView>(order+1);
  myQk.gaussLobattoQuadratureWeights( order, weights );
  weights2D=allocateVector<vectorDoubleView>((order+1)*(order+1));
  myQk.getGaussLobattoWeights2D( order, weights, weights2D );
  weights3D=allocateVector<vectorDoubleView>((order+1)*(order+1)*(order+1));
  myQk.getGaussLobattoWeights3D( order, weights, weights3D );

  // get basis function and corresponding derivatives
  basisFunction1D=allocateArray2D<arrayDoubleView>(order+1,order+1);
  myQk.getBasisFunction1D( order, quadraturePoints,basisFunction1D );

  derivativeBasisFunction1D=allocateArray2D<arrayDoubleView>(order+1,order+1);
  myQk.getDerivativeBasisFunction1D( order, quadraturePoints, derivativeBasisFunction1D );

  int nBasisFunctions=(order+1)*(order+1); 
  basisFunction2D=allocateArray2D<arrayDoubleView>(nBasisFunctions,nBasisFunctions);
  myQk.getBasisFunction2D( order, basisFunction1D, basisFunction1D, basisFunction2D );

  derivativeBasisFunction2DX=allocateArray2D<arrayDoubleView>(nBasisFunctions,nBasisFunctions);
  myQk.getBasisFunction2D( order, derivativeBasisFunction1D, basisFunction1D, derivativeBasisFunction2DX );
  
  derivativeBasisFunction2DY=allocateArray2D<arrayDoubleView>(nBasisFunctions,nBasisFunctions);
  myQk.getBasisFunction2D( order, basisFunction1D, derivativeBasisFunction1D, derivativeBasisFunction2DY );

  cout<<"fin allocation  shared:"<<endl;

  //shared arrays
  massMatrixGlobal=allocateVector<vectorRealView>( numberOfNodes );
  yGlobal=allocateVector<vectorRealView>( numberOfNodes );
  ShGlobal=allocateVector<vectorRealView>( numberOfBoundaryNodes );
  std::cout<<"end of shared arrays initialization\n";

}

// compute one step of the time dynamic wave equation solver
void SEMsolver::computeOneStep(  const int & timeStep,
                                  const float & timeSample,
                                  const int & order,
                                  int & i1,
                                  int & i2,
                                  const int & numberOfRHS,
                                  vectorIntView & rhsElement,
                                  arrayRealView & rhsTerm,
                                  arrayRealView & pnGlobal)
{


  Kokkos::parallel_for( numberOfNodes, KOKKOS_CLASS_LAMBDA ( const int i )
  {
    massMatrixGlobal[i]=0;
    yGlobal[i]=0;
  } );

  // update pnGLobal with right hade side
  Kokkos::parallel_for( numberOfRHS,KOKKOS_CLASS_LAMBDA (const int i)
  {
    int nodeRHS=globalNodesList(rhsElement[i],0);
    pnGlobal(nodeRHS,i2)+=timeSample*timeSample*model[rhsElement[i]]*model[rhsElement[i]]*rhsTerm(i,timeStep);
  });
 
  int numberOfPointsPerElement=(order+1)*(order+1);
  Kokkos::parallel_for( numberOfElements, KOKKOS_CLASS_LAMBDA ( const int e ) 
  {
    // start parallel section
    float B[36][4];
    float R[36];
    float massMatrixLocal[36];
    float pnLocal[36];
    float Y[36];

    // get pnGlobal to pnLocal
    for( int i=0; i<numberOfPointsPerElement; i++ )
    {
      int localToGlobal=globalNodesList(e,i);
      pnLocal[i]=pnGlobal(localToGlobal,i2);
    }

    // compute Jacobian, massMatrix and B
    int o=myQk.computeB( e,order,globalNodesList,globalNodesCoords,weights2D,
                       derivativeBasisFunction1D,massMatrixLocal,B );
    // compute stifness  matrix ( durufle's optimization)
    int p=myQk.gradPhiGradPhi( numberOfPointsPerElement, order, weights2D, derivativeBasisFunction1D, B, pnLocal, R, Y );
    // get pnGlobal to pnLocal
    for( int i=0; i<numberOfPointsPerElement; i++ )
    {
      int localToGlobal=globalNodesList(e,i);
      pnLocal[i]=pnGlobal(localToGlobal,i2);
    }
    //compute gloval mass Matrix and global stiffness vector
    for( int i=0; i<numberOfPointsPerElement; i++ )
    {
      int gIndex=globalNodesList(e,i);
      massMatrixLocal[i]/=(model[e]*model[e]);
      Kokkos::atomic_add(&massMatrixGlobal[gIndex],massMatrixLocal[i]);
      Kokkos::atomic_add(&yGlobal[gIndex],Y[i]);
    } 
  });

  // update pressure
  Kokkos::parallel_for( range_policy(0,numberOfInteriorNodes), KOKKOS_CLASS_LAMBDA ( const int i )
  {
    int I=listOfInteriorNodes[i];
    float tmp=timeSample*timeSample;
    pnGlobal(I,i1)=2*pnGlobal(I,i2)-pnGlobal(I,i1)-tmp*yGlobal[I]/massMatrixGlobal[I];
  } );

  Kokkos::parallel_for( range_policy(0,numberOfBoundaryNodes), KOKKOS_CLASS_LAMBDA ( const int i )
  {
    ShGlobal[i]=0;
  } );
  
  Kokkos::parallel_for (numberOfBoundaryFaces, KOKKOS_CLASS_LAMBDA (const int iFace)
  {
    //get ds
    float ds[6];
    float Sh[6];
    int numOfBasisFunctionOnFace[6];
    float Js[2][6];

    int i=myQk.computeDs( iFace, order, faceInfos,numOfBasisFunctionOnFace,
                  Js, globalNodesCoords, derivativeBasisFunction2DX,
                  derivativeBasisFunction2DY,
                  ds );
    //
    //compute Sh and ShGlobal
    for( int i=0; i<order+1; i++ )
    {
      int gIndexFaceNode=localFaceNodeToGlobalFaceNode(iFace,i);
      Sh[i]=weights[i]*ds[i]/(model[faceInfos(iFace,0)]);
      Kokkos::atomic_add(&ShGlobal[gIndexFaceNode],Sh[i]);
    }
  } );

  // update pressure @ boundaries;
  float tmp=timeSample*timeSample;
  Kokkos::parallel_for( range_policy(0,numberOfBoundaryNodes), KOKKOS_CLASS_LAMBDA  ( const int i )
  {
    int I=listOfBoundaryNodes[i];
    float invMpSh=1/(massMatrixGlobal[I]+timeSample*ShGlobal[i]*0.5);
    float MmSh=massMatrixGlobal[I]-timeSample*ShGlobal[i]*0.5;
    pnGlobal(I,i1)=invMpSh*(2*massMatrixGlobal[I]*pnGlobal(I,i2)-MmSh*pnGlobal(I,i1)-tmp*yGlobal[I]);
  } );

  Kokkos::fence();
}
