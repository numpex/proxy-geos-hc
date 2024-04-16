//************************************************************************
//   proxy application v.0.0.1
//
//  solverKokkos.hpp: simple 2D acoustive wave equation solver
//
//  the solverKokkos class is derived from the solverBase class
//  with the KOKKOS implementation of the solver
//
//************************************************************************

#include "solverKokkos.hpp"
#include <cstdio>
 
// compute one step of the time dynamic wave equation solver
void solverKokkos::computeOneStep( const int & timeStep,
                                   const float & timeSample,
                                   const int & order,
                                   int & i1,
                                   int & i2,
                                   const int & numberOfRHS,
                                   vectorInt & rhsElement,
                                   arrayReal & rhsTerm,
                                   arrayReal & pnGlobal)
{
  
  Kokkos::parallel_for( numberOfNodes, KOKKOS_CLASS_LAMBDA ( const int i )
  {
    massMatrixGlobal[i]=0;
    yGlobal[i]=0;
  } );
  // update pnGLobal with right hade side
  Kokkos::parallel_for(numberOfRHS,KOKKOS_CLASS_LAMBDA (const int i)
  {
    int nodeRHS=globalNodesList(rhsElement[i],0);
    pnGlobal(nodeRHS,i2)+=timeSample*timeSample*model[rhsElement[i]]*model[rhsElement[i]]*rhsTerm(i,timeStep);
  });
  for (int color=0; color<numberOfColors;color++)
  {
    Kokkos::parallel_for( numberOfElementsByColor[color], KOKKOS_CLASS_LAMBDA ( const int eColor )
    {
      // start parallel section
      int e=listOfElementsByColor(color,eColor);
      float B[36][4];
      float R[36];
      float massMatrixLocal[36];
      float pnLocal[36];
      float Y[36];
      // get pnGlobal to pnLocal
      for( int i=0; i<numberOfPointsPerElement; i++ )
      {
        int   localToGlobal=globalNodesList(e,i);
        pnLocal[i]=pnGlobal(localToGlobal,i2);
      }
      // compute Jacobian, massMatrix and B
      int o=Qk.computeB( e,order,globalNodesList,globalNodesCoords,weights2D,
                         derivativeBasisFunction1D,massMatrixLocal,B );
      // compute stifness  matrix ( durufle's optimization)
      int p=Qk.gradPhiGradPhi( numberOfPointsPerElement, order, weights2D, derivativeBasisFunction1D, B, pnLocal, R, Y );
      //compute gloval mass Matrix and global stiffness vector
      for( int i=0; i<numberOfPointsPerElement; i++ )
      {
        int gIndex=globalNodesList(e,i);
        massMatrixLocal[i]/=(model[e]*model[e]);
        massMatrixGlobal[gIndex]+=massMatrixLocal[i];
        yGlobal[gIndex]+=Y[i];
      } 
    });
    //Kokkos::fence();
  } // end color loop
  // update pressure
  Kokkos::parallel_for( range_policy(0,numberOfInteriorNodes), KOKKOS_CLASS_LAMBDA ( const int i )
  {
    int I=listOfInteriorNodes[i];
    float tmp=timeSample*timeSample;
    pnGlobal(I,i1)=2*pnGlobal(I,i2)-pnGlobal(I,i1)-tmp*yGlobal[I]/massMatrixGlobal[I];
  });
  //Kokkos::fence();
  // damping terms
  Kokkos::parallel_for( range_policy(0,numberOfBoundaryNodes), KOKKOS_CLASS_LAMBDA ( const int i )
  {
    ShGlobal[i]=0;
  });
  Kokkos::parallel_for (numberOfBoundaryFaces, KOKKOS_CLASS_LAMBDA (const int iFace)
  {
    float ds[6];
    float Sh[6];
    int numOfBasisFunctionOnFace[6];
    float Js[2][6];
    //get ds
    int i=Qk.computeDs( iFace, order, faceInfos,numOfBasisFunctionOnFace,
                        Js, globalNodesCoords, derivativeBasisFunction2DX,
                        derivativeBasisFunction2DY,
                        ds );

    //compute Sh and ShGlobal
    for( int i=0; i<order+1; i++ )
    {
      int gIndexFaceNode=localFaceNodeToGlobalFaceNode(iFace,i);
      Sh[i]=weights[i]*ds[i]/(model[faceInfos(iFace,0)]);
      //ShGlobal[gIndexFaceNode]+=Sh(threadId,i);
      Kokkos::atomic_add(&ShGlobal[gIndexFaceNode],Sh[i]);
    }
  });
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
