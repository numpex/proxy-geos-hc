//************************************************************************
//   proxy application v.0.0.1
//
//  solverRaja.cpp: simple 2D acoustive wave equation solver
//
//  the solverRaja class is derived from the solverBase class
//  with the RAJA implementation of the solver
//
//************************************************************************

#include "solver.hpp"

// compute one step of the time dynamic wave equation solver
void solverKokkos::computeOneStep(  const int & timeStep,
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
    int o=Qk.computeB( e,order,globalNodesList,globalNodesCoords,weights2D,
                       derivativeBasisFunction1D,massMatrixLocal,B );
    // compute stifness  matrix ( durufle's optimization)
    int p=Qk.gradPhiGradPhi( numberOfPointsPerElement, order, weights2D, derivativeBasisFunction1D, B, pnLocal, R, Y );
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

    int i=Qk.computeDs( iFace, order, faceInfos,numOfBasisFunctionOnFace,
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
