//************************************************************************
//   proxy application v.0.0.1
//
//  solverKokkos.cpp: 3D acoustive wave equation solver
//
//  the solverKokkos class is derived from the solverBase class
//  with the KOKKOS implementation of the solver
//
//************************************************************************

#include "solver.hpp"
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
      float B[64][4];
      float R[64];
      float massMatrixLocal[64];
      float pnLocal[64];
      float Y[64];

      // get pnGlobal to pnLocal
      for( int i=0; i<numberOfPointsPerElement; i++ )
      {
        pnLocal[i]=pnGlobal(localToGlobal[i],i2);
      }

      // compute Jacobian, massMatrix and B
      int o=Qk.computeB( e,order,globalNodesList,globalNodesCoords,weights3D,
                         derivativeBasisFunction1D,massMatrixLocal,B );

      // compute stifness  matrix ( durufle's optimization)
      int p=Qk.gradPhiGradPhi( numberOfPointsPerElement, order, weights3D,  derivativeBasisFunction1D,B, pnLocal, R, Y );

      //compute gloval mass Matrix and global stiffness vector
      for( int i=0; i<numberOfPointsPerElement; i++ )
      {
        int gIndex=localToGlobal[i];
        massMatrixLocal[i]/=(model[e]*model[e]);
        massMatrixGlobal[gIndex]+=massMatrixLocal[i];
        yGlobal[gIndex]+=Y[i];
      } 
    });
  } // end color loop

  // update pressure
  Kokkos::parallel_for( range_policy(0,numberOfInteriorNodes), KOKKOS_CLASS_LAMBDA ( const int i )
  {
    int I=listOfInteriorNodes[i];
    float tmp=timeSample*timeSample;
    pnGlobal(I,i1)=2*pnGlobal(I,i2)-pnGlobal(I,i1)-tmp*yGlobal[I]/massMatrixGlobal[I];
  });
 Kokkos::fence();
}
