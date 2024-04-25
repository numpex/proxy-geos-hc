//************************************************************************
//   proxy application v.0.0.1
//
//  SEMsolver.cpp: simple 2D acoustive wave equation solver
//
//  the SEMsolver class servers as a base class for the SEM solver
//
//************************************************************************

#include "SEMsolver.hpp"

void SEMsolver::computeFEInit( SEMmeshinfo &myMeshinfo, SEMmesh mesh )
{

  order = myMeshinfo.myOrderNumber;

  //interior elements
  globalNodesList=allocateArray2D<arrayIntView>(myMeshinfo.numberOfElements,myMeshinfo.numberOfPointsPerElement);
  mesh.globalNodesList( myMeshinfo.numberOfElements, globalNodesList );

  listOfInteriorNodes=allocateVector<vectorIntView>(myMeshinfo.numberOfInteriorNodes);
  mesh.getListOfInteriorNodes( myMeshinfo.numberOfInteriorNodes, listOfInteriorNodes );

  globalNodesCoords=allocateArray2D<arrayRealView>(myMeshinfo.numberOfNodes,3);
  mesh.nodesCoordinates( myMeshinfo.numberOfNodes, globalNodesCoords );

  // boundary elements
  listOfBoundaryNodes=allocateVector<vectorIntView>(myMeshinfo.numberOfBoundaryNodes);
  mesh.getListOfBoundaryNodes( myMeshinfo.numberOfBoundaryNodes, listOfBoundaryNodes );

  faceInfos=allocateArray2D<arrayIntView>(myMeshinfo.numberOfBoundaryFaces,2+(order+1));
  mesh.getBoundaryFacesInfos(faceInfos);

  localFaceNodeToGlobalFaceNode=allocateArray2D<arrayIntView>(myMeshinfo.numberOfBoundaryFaces, order+1 );
  mesh.getLocalFaceNodeToGlobalFaceNode(localFaceNodeToGlobalFaceNode);

  // get model
  model=allocateVector<vectorRealView>(myMeshinfo.numberOfElements);
  mesh.getModel( myMeshinfo.numberOfElements, model );

  // get quadrature points 
  quadraturePoints=allocateVector<vectorDoubleView>(order+1);
  myQk.gaussLobattoQuadraturePoints( order, quadraturePoints );

  // get gauss-lobatto weights
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

  basisFunction2D=allocateArray2D<arrayDoubleView>(myMeshinfo.nBasisFunctions,myMeshinfo.nBasisFunctions);
  myQk.getBasisFunction2D( order, basisFunction1D, basisFunction1D, basisFunction2D );

  derivativeBasisFunction2DX=allocateArray2D<arrayDoubleView>(myMeshinfo.nBasisFunctions,myMeshinfo.nBasisFunctions);
  myQk.getBasisFunction2D( order, derivativeBasisFunction1D, basisFunction1D, derivativeBasisFunction2DX );
  
  derivativeBasisFunction2DY=allocateArray2D<arrayDoubleView>(myMeshinfo.nBasisFunctions,myMeshinfo.nBasisFunctions);
  myQk.getBasisFunction2D( order, basisFunction1D, derivativeBasisFunction1D, derivativeBasisFunction2DY );

  cout<<"fin allocation  shared:"<<endl;

  //shared arrays
  massMatrixGlobal=allocateVector<vectorRealView>( myMeshinfo.numberOfNodes );
  yGlobal=allocateVector<vectorRealView>( myMeshinfo.numberOfNodes );
  ShGlobal=allocateVector<vectorRealView>( myMeshinfo.numberOfBoundaryNodes );
  std::cout<<"end of shared arrays initialization\n";

}

// compute one step of the time dynamic wave equation solver
void SEMsolver::computeOneStep(  const int & timeStep,
                                 SEMmeshinfo &myMeshinfo,
                                 int & i1,
                                 int & i2,
                                 vectorIntView & rhsElement,
                                 arrayRealView & rhsTerm,
                                 arrayRealView & pnGlobal)
{

  // update pressure @ boundaries;
  float tmp=myMeshinfo.myTimeStep * myMeshinfo.myTimeStep;

  Kokkos::parallel_for( myMeshinfo.numberOfNodes, KOKKOS_CLASS_LAMBDA ( const int i )
  {
    massMatrixGlobal[i]=0;
    yGlobal[i]=0;
  } );

  // update pnGLobal with right hade side
  Kokkos::parallel_for( myMeshinfo.myNumberOfRHS,KOKKOS_CLASS_LAMBDA (const int i)
  {
    int nodeRHS=globalNodesList(rhsElement[i],0);
    pnGlobal(nodeRHS,i2)+=tmp*model[rhsElement[i]]*model[rhsElement[i]]*rhsTerm(i,timeStep);
  });
 
  int numberOfPointsPerElement=(order+1)*(order+1);
  Kokkos::parallel_for( myMeshinfo.numberOfElements, KOKKOS_CLASS_LAMBDA ( const int e ) 
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
  Kokkos::parallel_for( range_policy(0,myMeshinfo.numberOfInteriorNodes), KOKKOS_CLASS_LAMBDA ( const int i )
  {
    int I=listOfInteriorNodes[i];
    pnGlobal(I,i1)=2*pnGlobal(I,i2)-pnGlobal(I,i1)-tmp*yGlobal[I]/massMatrixGlobal[I];
  } );

  Kokkos::parallel_for( range_policy(0,myMeshinfo.numberOfBoundaryNodes), KOKKOS_CLASS_LAMBDA ( const int i )
  {
    ShGlobal[i]=0;
  } );
  
  Kokkos::parallel_for (myMeshinfo.numberOfBoundaryFaces, KOKKOS_CLASS_LAMBDA (const int iFace)
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

  Kokkos::parallel_for( range_policy(0,myMeshinfo.numberOfBoundaryNodes), KOKKOS_CLASS_LAMBDA  ( const int i )
  {
    int I=listOfBoundaryNodes[i];
    float invMpSh=1/(massMatrixGlobal[I]+myMeshinfo.myTimeStep*ShGlobal[i]*0.5);
    float MmSh=massMatrixGlobal[I]-myMeshinfo.myTimeStep*ShGlobal[i]*0.5;
    pnGlobal(I,i1)=invMpSh*(2*massMatrixGlobal[I]*pnGlobal(I,i2)-MmSh*pnGlobal(I,i1)-tmp*yGlobal[I]);
  } );

  Kokkos::fence();
}


void SEMsolver::outputPnValues(  const int & indexTimeStep,
                                 int & i1,
                                 int & myElementSource, 
                                 arrayIntView & nodeList,
                                 arrayRealView & pnGlobal)
{
    //writes debugging ascii file.
    if( indexTimeStep%50==0 )
    {   
      cout<<"TimeStep="<<indexTimeStep<<endl;
    }   
    if( indexTimeStep%100==0 )
    {   
      cout<<" pnGlobal @ elementSource location "<<myElementSource
          <<" after computeOneStep = "<< pnGlobal(nodeList(myElementSource,0),i1)<<endl;
      //myUtils.saveSnapShot( indexTimeStep, i1, pnGlobal, myMesh );
    }  
}

