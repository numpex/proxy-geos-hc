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
  vectorDoubleView d_weights3D=weights3D.toView();

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
 
  int numberOfPointsPerElement=(order+1)*(order+1)*(order+1);

  RAJA::forall< deviceExecPolicy >( RAJA::RangeSegment( 0, numberOfElements ), [=] LVARRAY_HOST_DEVICE ( int e )
  {
    int nPointsPerElement=(order+1)*(order+1)*(order+1);
    // start parallel section
    double B[125][6];
    double R[125][125];
    double massMatrixLocal[125];
    double pnLocal[125];
    double Y[125];
    // compute Jacobian, massMatrix and B
    int o=Qk.computeB( e,order,d_globalNodesList,d_globalNodesCoords,d_weights3D,
                       d_derivativeBasisFunction1D,massMatrixLocal,B );
    // compute stifness  matrix ( durufle's optimization)
    int p=Qk.gradPhiGradPhi( numberOfPointsPerElement, order, d_weights3D, B, d_derivativeBasisFunction1D, R );
    // get pnGlobal to pnLocal
    for( int i=0; i<nPointsPerElement; i++ )
    {
      int localToGlobal=d_globalNodesList(e,i);
      massMatrixLocal[i]/=(d_model[e]*d_model[e]);
      pnLocal[i]=d_pnGlobal(localToGlobal,i2);
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
      int gIndex=d_globalNodesList(e,i);
      RAJA::atomicAdd< deviceAtomicPolicy >(&d_massMatrixGlobal[gIndex],massMatrixLocal[i]);
      RAJA::atomicAdd< deviceAtomicPolicy>(&d_yGlobal[gIndex],Y[i]);
    } 
  });
  // update pressure
  RAJA::forall< deviceExecPolicy>( RAJA::RangeSegment( 0, numberOfInteriorNodes ), [=] LVARRAY_HOST_DEVICE ( int i ) {
    int I=d_listOfInteriorNodes[i];
    float tmp=timeSample*timeSample;
    d_pnGlobal[I][i1]=2*d_pnGlobal[I][i2]-d_pnGlobal[I][i1]-tmp*d_yGlobal[I]/d_massMatrixGlobal[I];
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
