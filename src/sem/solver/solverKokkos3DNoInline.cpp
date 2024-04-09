//************************************************************************
//   proxy application v.0.0.1
//
//  solverRaja.cpp: simple 2D acoustive wave equation solver
//
//  the solverRaja class is derived from the solverBase class
//  with the RAJA implementation of the solver
//
//************************************************************************

#include "solverKokkos.hpp"

// compute one step of the time dynamic wave equation solver
void solverKokkos::computeOneStep(  const int & timeStep,
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
  Kokkos::parallel_for( numberOfRHS,KOKKOS_CLASS_LAMBDA (const int i)
  {
    int nodeRHS=globalNodesList(rhsElement[i],0);
    pnGlobal(nodeRHS,i2)+=timeSample*timeSample*model[rhsElement[i]]*model[rhsElement[i]]*rhsTerm(i,timeStep);
  });
 
  int numberOfPointsPerElement=(order+1)*(order+1);
  Kokkos::parallel_for( numberOfElements, KOKKOS_CLASS_LAMBDA ( const int e ) 
  {
    // start parallel section
    double B[125][6];
    double R[125][125];
    double massMatrixLocal[125];
    double pnLocal[125];
    double Y[125];

    // compute Jacobian, massMatrix and B
    int o=Qk.computeB( e,order,globalNodesList,globalNodesCoords,weights2D,
                       derivativeBasisFunction1D,massMatrixLocal,B );
    // compute stifness  matrix ( durufle's optimization)
    int p=Qk.gradPhiGradPhi( numberOfPointsPerElement, order, weights3D, B, derivativeBasisFunction1D, R );
    // get pnGlobal to pnLocal
    for( int i=0; i<numberOfPointsPerElement; i++ )
    {
      int localToGlobal=globalNodesList(e,i);
      massMatrixLocal[i]/=(model[e]*model[e]);
      pnLocal[i]=pnGlobal(localToGlobal,i2);
    }
    // compute Y=R*pnLocal
    for( int i=0; i<numberOfPointsPerElement; i++ )
    {
      Y[i]=0;
      for( int j=0; j<numberOfPointsPerElement; j++ )
      {
        Y[i]+=R[i][j]*pnLocal[j];
      }
    }
    //compute gloval mass Matrix and global stiffness vector
    for( int i=0; i<numberOfPointsPerElement; i++ )
    {
      int gIndex=globalNodesList(e,i);
      //massMatrixGlobal[gIndex]+=massMatrixLocal(threadId,i)
      //yGlobal[gIndex]+=Y(threadId,i);
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

  Kokkos::fence();
}
