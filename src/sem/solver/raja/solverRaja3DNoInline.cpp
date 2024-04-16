//************************************************************************
//   proxy application v.0.0.1
//
//  solverRaja.cpp: 3D acoustive wave equation solver
//
//  the solverRaja class is derived from the solverBase class
//  with the RAJA implementation of the solver
//
//************************************************************************

#include "solver.hpp"

// compute one step of the time dynamic wave equation solver
void solverRaja::computeOneStep(  const int & timeStep,
                                  const float & timeSample,
                                  const int & order,
                                  int & it1,
                                  int & it2,
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
  //shared arrays
  vectorRealView d_massMatrixGlobal=massMatrixGlobal.toView();
  vectorRealView d_yGlobal=yGlobal.toView();

  RAJA::forall< deviceExecPolicy >( RAJA::RangeSegment( 0, numberOfNodes ),  [=] LVARRAY_HOST_DEVICE  ( int i ) {
    d_massMatrixGlobal[i]=0;
    d_yGlobal[i]=0;
  } );

  // update pnGLobal with right hade side
  RAJA::forall< deviceExecPolicy >( RAJA::RangeSegment( 0, numberOfRHS ), [=] LVARRAY_HOST_DEVICE  ( int i ) 
  {
    int nodeRHS=d_globalNodesList(d_rhsElement[i],0);
    d_pnGlobal(nodeRHS,it2)+=timeSample*timeSample*d_model[d_rhsElement[i]]*d_model[d_rhsElement[i]]*d_rhsTerm(i,timeStep);
  });
 
  int numberOfPointsPerElement=(order+1)*(order+1)*(order+1);

  RAJA::forall< deviceExecPolicy >( RAJA::RangeSegment( 0, numberOfElements ), [=] LVARRAY_HOST_DEVICE ( int e )
  {
    int nPointsPerElement=(order+1)*(order+1)*(order+1);
    // start parallel section
    float B[64][6];
    float R[64];
    float massMatrixLocal[64];
    float pnLocal[64];
    float Y[64];

    // get pnGlobal to pnLocal
    for( int i=0; i<numberOfPointsPerElement; i++ )
    {
         int localToGlobal=d_globalNodesList(e,i);
         pnLocal[i]=d_pnGlobal(localToGlobal,it2);
    }
    // compute Jacobian, massMatrix and B
    int o=Qk.computeB( e,order,d_globalNodesList,d_globalNodesCoords,d_weights3D,
                       d_derivativeBasisFunction1D,massMatrixLocal,B );
    // compute stifness  matrix ( durufle's optimization)
    int p=Qk.gradPhiGradPhi( numberOfPointsPerElement, order, d_weights3D, d_derivativeBasisFunction1D, B, pnLocal, R, Y );

    //compute global mass Matrix and global stiffness vector
    for( int i=0; i<nPointsPerElement; i++ )
    {
      int gIndex=d_globalNodesList(e,i);
      massMatrixLocal[i]/=(d_model[e]*d_model[e]);
      RAJA::atomicAdd< deviceAtomicPolicy >(&d_massMatrixGlobal[gIndex],massMatrixLocal[i]);
      RAJA::atomicAdd< deviceAtomicPolicy>(&d_yGlobal[gIndex],Y[i]);
    } 
  });
  // update pressure
  RAJA::forall< deviceExecPolicy>( RAJA::RangeSegment( 0, numberOfInteriorNodes ), [=] LVARRAY_HOST_DEVICE ( int i ) {
    int I=d_listOfInteriorNodes[i];
    float tmp=timeSample*timeSample;
    d_pnGlobal[I][it1]=2*d_pnGlobal[I][it2]-d_pnGlobal[I][it1]-tmp*d_yGlobal[I]/d_massMatrixGlobal[I];
  } );

  if(timeStep%100==0)
  {
     int nodeRHS=d_globalNodesList(d_rhsElement[0],0);
     RAJA::forall< RAJA::seq_exec>( RAJA::RangeSegment( 0, numberOfNodes ), [=] LVARRAY_HOST_DEVICE ( int i )
     {
	  pnGlobal(nodeRHS,it1)=d_pnGlobal(nodeRHS,it1);
     });
  }
}
