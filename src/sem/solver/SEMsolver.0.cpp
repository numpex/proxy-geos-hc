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
  numberOfPointsPerElement=pow((order+1),DIMENSION);

  allocateFEarrays( myMeshinfo );
  initFEarrays( myMeshinfo, mesh);

}

// compute one step of the time dynamic wave equation solver
void SEMsolver::computeOneStep(  const int & timeStep,
                                 SEMmeshinfo &myMeshinfo,
                                 int & i1,
                                 int & i2,
                                 arrayReal & myRHSTerm,
                                 arrayReal & myPnGlobal,
                                 vectorInt & myRhsElement)
{
  CREATEVIEWS
   
  printf("myPnGlobal(0,0)=%f; i1=%d;  ",myPnGlobal(globalNodesList(2178,0),i1), i1);
  printf("1. pnGlobal(0,0)=%f  ;",pnGlobal(globalNodesList(2178,0),i1));

  // update pressure @ boundaries;
  LOOPHEAD( myMeshinfo.numberOfNodes, i)
    massMatrixGlobal[i]=0;
    yGlobal[i]=0;
  LOOPEND

  // update pnGLobal with right hade side
  LOOPHEAD( myMeshinfo.myNumberOfRHS, i)
    int nodeRHS=globalNodesList(rhsElement[i],0);
    pnGlobal(nodeRHS,i2)+=myMeshinfo.myTimeStep*myMeshinfo.myTimeStep*model[rhsElement[i]]*model[rhsElement[i]]*rhsTerm(i,timeStep);
  LOOPEND
  printf("2. pnGlobal(0,0)=%f;  ",pnGlobal(globalNodesList(2178,0),i1));
 
  LOOPHEAD( myMeshinfo.numberOfElements, e)
    // start parallel section
    float B[ROW][COL];
    float R[ROW];
    float massMatrixLocal[ROW];
    float pnLocal[ROW];
    float Y[ROW];

    // get pnGlobal to pnLocal
    for( int i=0; i<pow((myMeshinfo.myOrderNumber+1),DIMENSION); i++ )
    {
      int localToGlobal=globalNodesList(e,i);
      pnLocal[i]=pnGlobal(localToGlobal,i2);
    }

    // compute Jacobian, massMatrix and B
    myQk.computeB( e,myMeshinfo.myOrderNumber,DIMENSION, weights, globalNodesList,globalNodesCoords,derivativeBasisFunction1D,massMatrixLocal,B );

    // compute stifness  matrix ( durufle's optimization)
    myQk.gradPhiGradPhi( pow((myMeshinfo.myOrderNumber+1),DIMENSION), myMeshinfo.myOrderNumber, DIMENSION,weights,derivativeBasisFunction1D, B, pnLocal, R, Y );

    //compute gloval mass Matrix and global stiffness vector
    for( int i=0; i<pow((myMeshinfo.myOrderNumber+1),DIMENSION); i++ )
    {
      int gIndex=globalNodesList(e,i);
      massMatrixLocal[i]/=(model[e]*model[e]);
      ATOMICADD( massMatrixGlobal[gIndex], massMatrixLocal[i] );
      ATOMICADD( yGlobal[gIndex], Y[i]);
    } 

  LOOPEND

  // update pressure
  LOOPHEAD( myMeshinfo.numberOfInteriorNodes, i)
    int I=listOfInteriorNodes[i];
    pnGlobal(I,i1)=2*pnGlobal(I,i2)-pnGlobal(I,i1)-myMeshinfo.myTimeStep*myMeshinfo.myTimeStep*yGlobal[I]/massMatrixGlobal[I];
  LOOPEND
  printf("3. pnGlobal(0,0)=%f;  ",pnGlobal(globalNodesList(2178,0),i1));

  if (DIMENSION==2) {
  LOOPHEAD( myMeshinfo.numberOfBoundaryNodes, i)
    ShGlobal[i]=0;
  LOOPEND

  LOOPHEAD( myMeshinfo.numberOfBoundaryFaces, iFace)
    //get ds
    float ds[6];
    float Sh[6];
    int numOfBasisFunctionOnFace[6];
    float Js[2][6];

    int i=myQk.computeDs( iFace, myMeshinfo.myOrderNumber, faceInfos,numOfBasisFunctionOnFace,
                  Js, globalNodesCoords, derivativeBasisFunction1D,
                  ds );

    //compute Sh and ShGlobal
    for( int i=0; i<myMeshinfo.myOrderNumber+1; i++ )
    {
      int gIndexFaceNode=localFaceNodeToGlobalFaceNode(iFace,i);
      Sh[i]=weights[i]*ds[i]/(model[faceInfos(iFace,0)]);
      ATOMICADD(ShGlobal[gIndexFaceNode], Sh[i]);
    }
  LOOPEND

  LOOPHEAD( myMeshinfo.numberOfBoundaryNodes, i)
    int I=listOfBoundaryNodes[i];
    float invMpSh=1/(massMatrixGlobal[I]+myMeshinfo.myTimeStep*ShGlobal[i]*0.5);
    float MmSh=massMatrixGlobal[I]-myMeshinfo.myTimeStep*ShGlobal[i]*0.5;
    pnGlobal(I,i1)=invMpSh*(2*massMatrixGlobal[I]*pnGlobal(I,i2)-MmSh*pnGlobal(I,i1)-myMeshinfo.myTimeStep*myMeshinfo.myTimeStep*yGlobal[I]);
  LOOPEND
  }
  //printf("4. pnGlobal(0,0)=%f\n",pnGlobal(globalNodesList(2178,0),i1));
  myPnGlobal=pnGlobal;
  printf("4. myPnGlobal(0,0)=%f\n",myPnGlobal(globalNodesList(2178,0),i1));
  FENCE
}


void SEMsolver::outputPnValues(  SEMmesh mesh,
		                 const int & indexTimeStep,
                                 int & i1,
                                 int & myElementSource, 
                                 arrayReal & pnGlobal)
{
    //writes debugging ascii file.
    if( indexTimeStep%50==0 )
    {   
      cout<<"TimeStep="<<indexTimeStep<<endl;
    }   
    if( indexTimeStep%100==0 )
    {   
      cout<<" pnGlobal @ elementSource location "<<myElementSource
          <<" after computeOneStep = "<< pnGlobal(globalNodesList(myElementSource,0),i1)<<endl;
      //mesh.saveSnapShot( indexTimeStep, i1, pnGlobal );
    }  
}


void SEMsolver::initFEarrays( SEMmeshinfo &myMeshinfo, SEMmesh mesh )
{
  //interior elements
  mesh.globalNodesList( myMeshinfo.numberOfElements, globalNodesList );
  mesh.getListOfInteriorNodes( myMeshinfo.numberOfInteriorNodes, listOfInteriorNodes );
  mesh.nodesCoordinates( myMeshinfo.numberOfNodes, globalNodesCoords );
  // boundary elements
  mesh.getListOfBoundaryNodes( myMeshinfo.numberOfBoundaryNodes, listOfBoundaryNodes );
  mesh.getBoundaryFacesInfos(faceInfos);
  mesh.getLocalFaceNodeToGlobalFaceNode(localFaceNodeToGlobalFaceNode);
  // get model
  mesh.getModel( myMeshinfo.numberOfElements, model );
  // get quadrature points 
  myQk.gaussLobattoQuadraturePoints( order, quadraturePoints );
  // get gauss-lobatto weights
  myQk.gaussLobattoQuadratureWeights( order, weights );
  // get basis function and corresponding derivatives
  myQk.getBasisFunction1D( order, quadraturePoints,basisFunction1D );
  myQk.getDerivativeBasisFunction1D( order, quadraturePoints, derivativeBasisFunction1D );

}

void SEMsolver::allocateFEarrays( SEMmeshinfo &myMeshinfo )
{
  //interior elements
  globalNodesList=allocateArray2D<arrayInt>(myMeshinfo.numberOfElements,myMeshinfo.numberOfPointsPerElement);
  listOfInteriorNodes=allocateVector<vectorInt>(myMeshinfo.numberOfInteriorNodes);
  globalNodesCoords=allocateArray2D<arrayReal>(myMeshinfo.numberOfNodes,3);
  listOfBoundaryNodes=allocateVector<vectorInt>(myMeshinfo.numberOfBoundaryNodes);
  faceInfos=allocateArray2D<arrayInt>(myMeshinfo.numberOfBoundaryFaces,2+(order+1));
  localFaceNodeToGlobalFaceNode=allocateArray2D<arrayInt>(myMeshinfo.numberOfBoundaryFaces, order+1 );
  model=allocateVector<vectorReal>(myMeshinfo.numberOfElements);
  quadraturePoints=allocateVector<vectorDouble>(order+1);
  weights=allocateVector<vectorDouble>(order+1);
  basisFunction1D=allocateArray2D<arrayDouble>(order+1,order+1);
  derivativeBasisFunction1D=allocateArray2D<arrayDouble>(order+1,order+1);
  //shared arrays
  massMatrixGlobal=allocateVector<vectorReal>( myMeshinfo.numberOfNodes );
  yGlobal=allocateVector<vectorReal>( myMeshinfo.numberOfNodes );
  ShGlobal=allocateVector<vectorReal>( myMeshinfo.numberOfBoundaryNodes );
}

