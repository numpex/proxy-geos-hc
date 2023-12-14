//************************************************************************
//   proxy application v.0.0.1
//
//  solverRaja.cpp: simple 2D acoustive wave equation solver
//
//  the solverRaja class is derived from the solverBase class
//  with the RAJA implementation of the solver
//
//************************************************************************

#include "solverRaja.hpp"

// compute one step of the time dynamic wave equation solver
void solverRaja::computeOneStep(  const int & timeStep,
                                  const float & timeSample,
                                  const int & order,
                                  int & i1,
                                  int & i2,
                                  const int & numberOfRHS,
                                  vectorInt & rhsElement,
                                  arrayReal & rhsTerm,
                                  arrayReal & pnGlobal)
{

  vectorIntView d_rhsElement=rhsElement.toView();
  arrayRealView d_rhsTerm=rhsTerm.toView();
  arrayRealView d_pnGlobal=pnGlobal.toView();

  //shared arrays
  arrayIntView d_globalNodesList=globalNodesList.toView();
  arrayRealView d_globalNodesCoords=globalNodesCoords.toView();
  vectorIntView d_listOfInteriorNodes=listOfInteriorNodes.toView();
  vectorIntView d_listOfBoundaryNodes=listOfBoundaryNodes.toView();
  arrayIntView d_faceInfos=faceInfos.toView();
  arrayIntView d_localFaceNodeToGlobalFaceNode=localFaceNodeToGlobalFaceNode.toView();
  
  // get model
  vectorRealView d_model=model.toView();

  // get quadrature points and weights
  vectorDoubleView d_quadraturePoints=quadraturePoints.toView();
  vectorDoubleView d_weights=weights.toView();
  vectorDoubleView d_weights2D=weights2D.toView();

  // get basis function and corresponding derivatives
  arrayDoubleView d_basisFunction1D=basisFunction1D.toView();
  arrayDoubleView d_derivativeBasisFunction1D=derivativeBasisFunction1D.toView();
  arrayDoubleView d_basisFunction2D=basisFunction2D.toView();
  arrayDoubleView d_derivativeBasisFunction2DX=derivativeBasisFunction2DX.toView();
  arrayDoubleView d_derivativeBasisFunction2DY=derivativeBasisFunction2DY.toView();

  //shared arrays
  vectorDoubleView d_massMatrixGlobal=massMatrixGlobal.toView();
  vectorDoubleView d_yGlobal=yGlobal.toView();
  vectorRealView d_ShGlobal=ShGlobal.toView();

  RAJA::forall< deviceExecPolicy >( RAJA::RangeSegment( 0, numberOfNodes ),  [=] LVARRAY_HOST_DEVICE  ( int i ) {
    d_massMatrixGlobal[i]=0;
    d_yGlobal[i]=0;
  } );

  // update pnGLobal with right hade side
  RAJA::forall< deviceExecPolicy >( RAJA::RangeSegment( 0, numberOfRHS ), [=] LVARRAY_HOST_DEVICE  ( int i ) 
  {
    int nodeRHS=d_globalNodesList(d_rhsElement[i],0);
    d_pnGlobal(nodeRHS,i2)+=timeSample*timeSample*d_model[d_rhsElement[i]]*d_model[d_rhsElement[i]]*d_rhsTerm(i,timeStep);
  });
 
  int numberOfPointsPerElement=(order+1)*(order+1);


  // loop over mesh elements
  RAJA::forall< deviceExecPolicy >( RAJA::RangeSegment( 0, numberOfElements ), [=] LVARRAY_HOST_DEVICE ( int e )
  {
    int nPointsPerElement=(order+1)*(order+1);
    // start parallel section
    int  localToGlobal[36];
    double Xi[36][2];

    double jacobianMatrix[36][4];
    double detJ[36];
    double invJacobianMatrix[36][4];
    double transpInvJacobianMatrix[36][4];

    double B[36][4];
    double R[36][36];

    double massMatrixLocal[36];

    double pnLocal[36];
    double Y[36];

    // extract global coordinates of element e
    // get local to global indexes of nodes of element e
    int i=mesh.localToGlobalNodes( e, numberOfPointsPerElement, d_globalNodesList, localToGlobal );
    //get global coordinates Xi of element e
    int j=mesh.getXi( nPointsPerElement, d_globalNodesCoords, localToGlobal, Xi );
    // compute jacobian Matrix
    int k=Qk.computeJacobianMatrix( numberOfPointsPerElement, Xi,
                              d_derivativeBasisFunction2DX,
                              d_derivativeBasisFunction2DY,
                              jacobianMatrix );
    // compute determinant of jacobian Matrix
    int l=Qk.computeDeterminantOfJacobianMatrix( numberOfPointsPerElement,
                                           jacobianMatrix,
                                           detJ );
    // compute inverse of Jacobian Matrix
    int m=Qk.computeInvJacobianMatrix( numberOfPointsPerElement,
                                 jacobianMatrix,
                                 detJ,
                                 invJacobianMatrix );
    // compute transposed inverse of Jacobian Matrix
    int n=Qk.computeTranspInvJacobianMatrix( numberOfPointsPerElement,
                                       jacobianMatrix,
                                       detJ,
                                       transpInvJacobianMatrix );
    // compute  geometrical transformation matrix
    int o=Qk.computeB( numberOfPointsPerElement, invJacobianMatrix, transpInvJacobianMatrix, detJ,B );
    // compute stifness and mass matrix ( durufle's optimization)
    int p=Qk.gradPhiGradPhi( numberOfPointsPerElement, order, d_weights2D, B, d_derivativeBasisFunction1D, R );
    // compute local mass matrix ( used optimez version)
    int q=Qk.phiIphiJ( numberOfPointsPerElement, d_weights2D, detJ, massMatrixLocal );
    // get pnGlobal to pnLocal
    for( int i=0; i<nPointsPerElement; i++ )
    {
      massMatrixLocal[i]/=(d_model[e]*d_model[e]);
      pnLocal[i]=d_pnGlobal(localToGlobal[i],i2);
    }
    // compute Y=R*pnLocal
    for( int i=0; i<nPointsPerElement; i++ )
    {
      Y[i]=0;
      for( int j=0; j<nPointsPerElement; j++ )
      {
        Y[i]+=R[i][j]*pnLocal[j];
      }
    }
    //compute gloval mass Matrix and global stiffness vector
    for( int i=0; i<nPointsPerElement; i++ )
    {
      int gIndex=localToGlobal[i];
      //massMatrixGlobal[gIndex]+=massMatrixLocal(threadId,i)
      //yGlobal[gIndex]+=Y(threadId,i);
      RAJA::atomicAdd< deviceAtomicPolicy >(&d_massMatrixGlobal[gIndex],massMatrixLocal[i]);
      RAJA::atomicAdd< deviceAtomicPolicy>(&d_yGlobal[gIndex],Y[i]);
    } 
  } );
  // update pressure
  RAJA::forall< deviceExecPolicy>( RAJA::RangeSegment( 0, numberOfInteriorNodes ), [=] LVARRAY_HOST_DEVICE ( int i ) {
    int I=d_listOfInteriorNodes[i];
    float tmp=timeSample*timeSample;
    d_pnGlobal[I][i1]=2*d_pnGlobal[I][i2]-d_pnGlobal[I][i1]-tmp*d_yGlobal[I]/d_massMatrixGlobal[I];
  } );
  RAJA::forall< deviceExecPolicy>( RAJA::RangeSegment( 0, numberOfBoundaryNodes ), [=] LVARRAY_HOST_DEVICE ( int i ) {
    d_ShGlobal[i]=0;
  } );
  RAJA::forall< deviceExecPolicy >( RAJA::RangeSegment( 0, numberOfBoundaryFaces ), [=] LVARRAY_HOST_DEVICE ( int iFace ){
    //get ds
    float ds[6];
    float Sh[6];
    int numOfBasisFunctionOnFace[6];
    float Js[2][6];

    
    int i=Qk.computeDs( iFace, order, d_faceInfos,numOfBasisFunctionOnFace,
                  Js, d_globalNodesCoords, d_derivativeBasisFunction2DX,
                  d_derivativeBasisFunction2DY,
                  ds );
    //compute Sh and ShGlobal
    for( int i=0; i<order+1; i++ )
    {
      int gIndexFaceNode=d_localFaceNodeToGlobalFaceNode(iFace,i);
      Sh[i]=d_weights[i]*ds[i]/(d_model[d_faceInfos(iFace,0)]);
      RAJA::atomicAdd< deviceAtomicPolicy >(&d_ShGlobal[gIndexFaceNode],Sh[i]);
    }
  } );
  // update pressure @ boundaries;
  float tmp=timeSample*timeSample;
  RAJA::forall< deviceExecPolicy >( RAJA::RangeSegment( 0, numberOfBoundaryNodes ), [=] LVARRAY_HOST_DEVICE ( int i ) {
    int I=d_listOfBoundaryNodes[i];
    float invMpSh=1/(d_massMatrixGlobal[I]+timeSample*d_ShGlobal[i]*0.5);
    float MmSh=d_massMatrixGlobal[I]-timeSample*d_ShGlobal[i]*0.5;
    d_pnGlobal[I][i1]=invMpSh*(2*d_massMatrixGlobal[I]*d_pnGlobal[I][i2]-MmSh*d_pnGlobal[I][i1]-tmp*d_yGlobal[I]);
  } );
  if(timeStep%100==0)
  {
     int nodeRHS=d_globalNodesList(d_rhsElement[0],0);
     RAJA::forall< RAJA::seq_exec>( RAJA::RangeSegment( 0, numberOfNodes ), [=] LVARRAY_HOST_DEVICE ( int i )
     {
	  pnGlobal(nodeRHS,i1)=d_pnGlobal(nodeRHS,i1);
     });
  }
}
