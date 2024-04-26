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
  tmp=myMeshinfo.myTimeStep * myMeshinfo.myTimeStep;
  numberOfPointsPerElement=pow((order+1),DIMENSION);

  allocateFEarrays( myMeshinfo );
  initFEarrays( myMeshinfo, mesh);

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
  LOOPHEAD( myMeshinfo.numberOfNodes, i)
    massMatrixGlobal[i]=0;
    yGlobal[i]=0;
  LOOPEND

  // update pnGLobal with right hade side
  LOOPHEAD( myMeshinfo.myNumberOfRHS, i)
    int nodeRHS=globalNodesList(rhsElement[i],0);
    pnGlobal(nodeRHS,i2)+=tmp*model[rhsElement[i]]*model[rhsElement[i]]*rhsTerm(i,timeStep);
  LOOPEND
 
  LOOPHEAD( myMeshinfo.numberOfElements, e)
    // start parallel section
    float B[ROW][COL];
    float R[ROW];
    float massMatrixLocal[ROW];
    float pnLocal[ROW];
    float Y[ROW];

    // get pnGlobal to pnLocal
    for( int i=0; i<numberOfPointsPerElement; i++ )
    {
      int localToGlobal=globalNodesList(e,i);
      pnLocal[i]=pnGlobal(localToGlobal,i2);
    }

    // compute Jacobian, massMatrix and B
    myQk.computeB( e,order,DIMENSION, weights, globalNodesList,globalNodesCoords,derivativeBasisFunction1D,massMatrixLocal,B );
    // compute stifness  matrix ( durufle's optimization)
    myQk.gradPhiGradPhi( numberOfPointsPerElement, order, DIMENSION,weights,derivativeBasisFunction1D, B, pnLocal, R, Y );
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
      ATOMICADD( massMatrixGlobal[gIndex], massMatrixLocal[i] );
      ATOMICADD( yGlobal[gIndex], Y[i]);
    } 
  LOOPEND

  // update pressure
  LOOPHEAD( myMeshinfo.numberOfInteriorNodes, i)
    int I=listOfInteriorNodes[i];
    pnGlobal(I,i1)=2*pnGlobal(I,i2)-pnGlobal(I,i1)-tmp*yGlobal[I]/massMatrixGlobal[I];
  LOOPEND

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
      ATOMICADD(ShGlobal[gIndexFaceNode], Sh[i]);
    }
  LOOPEND

  LOOPHEAD( myMeshinfo.numberOfBoundaryNodes, i)
    int I=listOfBoundaryNodes[i];
    float invMpSh=1/(massMatrixGlobal[I]+myMeshinfo.myTimeStep*ShGlobal[i]*0.5);
    float MmSh=massMatrixGlobal[I]-myMeshinfo.myTimeStep*ShGlobal[i]*0.5;
    pnGlobal(I,i1)=invMpSh*(2*massMatrixGlobal[I]*pnGlobal(I,i2)-MmSh*pnGlobal(I,i1)-tmp*yGlobal[I]);
  LOOPEND
  }
  FENCE
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
  myQk.getBasisFunction2D( order, basisFunction1D, basisFunction1D, basisFunction2D );
  myQk.getBasisFunction2D( order, derivativeBasisFunction1D, basisFunction1D, derivativeBasisFunction2DX );
  myQk.getBasisFunction2D( order, basisFunction1D, derivativeBasisFunction1D, derivativeBasisFunction2DY );

}

void SEMsolver::allocateFEarrays( SEMmeshinfo &myMeshinfo )
{
  #ifdef USE_RAJA
  h_globalNodesList=allocateArray2D<arrayInt>(myMeshinfo.numberOfElements,myMeshinfo.numberOfPointsPerElement);
  h_listOfInteriorNodes=allocateVector<vectorInt>(myMeshinfo.numberOfInteriorNodes);
  h_globalNodesCoords=allocateArray2D<arrayReal>(myMeshinfo.numberOfNodes,3);
  h_listOfBoundaryNodes=allocateVector<vectorInt>(myMeshinfo.numberOfBoundaryNodes);
  h_faceInfos=allocateArray2D<arrayInt>(myMeshinfo.numberOfBoundaryFaces,2+(order+1));
  h_localFaceNodeToGlobalFaceNode=allocateArray2D<arrayInt>(myMeshinfo.numberOfBoundaryFaces, order+1 );
  h_model=allocateVector<vectorReal>(myMeshinfo.numberOfElements);
  h_quadraturePoints=allocateVector<vectorDouble>(order+1);
  h_weights=allocateVector<vectorDouble>(order+1);
  h_basisFunction1D=allocateArray2D<arrayDouble>(order+1,order+1);
  h_derivativeBasisFunction1D=allocateArray2D<arrayDouble>(order+1,order+1);
  h_basisFunction2D=allocateArray2D<arrayDouble>(myMeshinfo.nBasisFunctions,myMeshinfo.nBasisFunctions);
  h_derivativeBasisFunction2DX=allocateArray2D<arrayDouble>(myMeshinfo.nBasisFunctions,myMeshinfo.nBasisFunctions);
  h_derivativeBasisFunction2DY=allocateArray2D<arrayDouble>(myMeshinfo.nBasisFunctions,myMeshinfo.nBasisFunctions);
  h_massMatrixGlobal=allocateVector<vectorReal>( myMeshinfo.numberOfNodes );
  h_yGlobal=allocateVector<vectorReal>( myMeshinfo.numberOfNodes );
  h_ShGlobal=allocateVector<vectorReal>( myMeshinfo.numberOfBoundaryNodes );

  //get Views
  globalNodesList=h_globalNodesList.toView();
  listOfInteriorNodes=h_listOfInteriorNodes.toView();
  globalNodesCoords=h_globalNodesCoords.toView();
  listOfBoundaryNodes=h_listOfBoundaryNodes.toView();
  faceInfos=h_faceInfos.toView();
  localFaceNodeToGlobalFaceNode=h_localFaceNodeToGlobalFaceNode.toView();
  model=h_model.toView();
  quadraturePoints=h_quadraturePoints.toView();
  weights=h_weights.toView();
  basisFunction1D=h_basisFunction1D.toView();
  derivativeBasisFunction1D=h_derivativeBasisFunction1D.toView();
  basisFunction2D=h_basisFunction2D.toView();
  derivativeBasisFunction2DX=h_derivativeBasisFunction2DX.toView();
  derivativeBasisFunction2DY=h_derivativeBasisFunction2DY.toView();
  //sharedarrays//sharedarrays
  massMatrixGlobal=h_massMatrixGlobal.toView();
  yGlobal=h_yGlobal.toView();
  ShGlobal=h_ShGlobal.toView();

  #else

  //interior elements
  globalNodesList=allocateArray2D<arrayIntView>(myMeshinfo.numberOfElements,myMeshinfo.numberOfPointsPerElement);
  listOfInteriorNodes=allocateVector<vectorIntView>(myMeshinfo.numberOfInteriorNodes);
  globalNodesCoords=allocateArray2D<arrayRealView>(myMeshinfo.numberOfNodes,3);
  listOfBoundaryNodes=allocateVector<vectorIntView>(myMeshinfo.numberOfBoundaryNodes);
  faceInfos=allocateArray2D<arrayIntView>(myMeshinfo.numberOfBoundaryFaces,2+(order+1));
  localFaceNodeToGlobalFaceNode=allocateArray2D<arrayIntView>(myMeshinfo.numberOfBoundaryFaces, order+1 );
  model=allocateVector<vectorRealView>(myMeshinfo.numberOfElements);
  quadraturePoints=allocateVector<vectorDoubleView>(order+1);
  weights=allocateVector<vectorDoubleView>(order+1);
  basisFunction1D=allocateArray2D<arrayDoubleView>(order+1,order+1);
  derivativeBasisFunction1D=allocateArray2D<arrayDoubleView>(order+1,order+1);
  basisFunction2D=allocateArray2D<arrayDoubleView>(myMeshinfo.nBasisFunctions,myMeshinfo.nBasisFunctions);
  derivativeBasisFunction2DX=allocateArray2D<arrayDoubleView>(myMeshinfo.nBasisFunctions,myMeshinfo.nBasisFunctions);
  derivativeBasisFunction2DY=allocateArray2D<arrayDoubleView>(myMeshinfo.nBasisFunctions,myMeshinfo.nBasisFunctions);
  //shared arrays
  massMatrixGlobal=allocateVector<vectorRealView>( myMeshinfo.numberOfNodes );
  yGlobal=allocateVector<vectorRealView>( myMeshinfo.numberOfNodes );
  ShGlobal=allocateVector<vectorRealView>( myMeshinfo.numberOfBoundaryNodes );

  #endif
}


