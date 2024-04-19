//************************************************************************
//   proxy application v.0.0.1
//
//  solverOMP.cpp: 3D acoustive wave equation solver
//
//  the solverOMP class is derived from the solverBase class
//  with the openMP implementation of the solver
//
//************************************************************************

#include "solver.hpp"

// compute one step of the time dynamic wave equation solver
void solverOMP::computeOneStep( const int & timeStep,
                                const float & timeSample,
                                const int & order,
                                int & i1,
                                int & i2,
                                const int & numberOfRHS,
                                vectorInt & rhsElement,
                                arrayReal & rhsTerm,
                                arrayReal & pnGlobal)
{
#pragma omp parallel
  { // start parallel section

    float B[64][6];
    float R[64];
    float massMatrixLocal[64];
    float pnLocal[64];
    float Y[64];


#pragma omp for 
  for( int i=0; i<numberOfNodes; i++ )
  {
    massMatrixGlobal[i]=0;
    yGlobal[i]=0;
  }

  // update pnGLobal with right hade side
#pragma omp for 
  for( int i=0; i<numberOfRHS;i++)
  { 
    int nodeRHS=globalNodesList(rhsElement[i],0);
    pnGlobal(nodeRHS,i2)+=timeSample*timeSample*model[rhsElement[i]]*model[rhsElement[i]]*rhsTerm(i,timeStep);
  }
  // loop over mesh elements
#pragma omp for
  for( int e=0; e<numberOfElements; e++ )
  {

    // get pnGlobal to pnLocal
    for( int i=0; i<numberOfPointsPerElement; i++ )
    {
      int localToGlobal=globalNodesList(e,i);
      pnLocal[i]=pnGlobal(localToGlobal,i2);
    }

    // compute Jacobian, massMatrix and B
    int o=Qk.computeB( e,order,globalNodesList,globalNodesCoords,weights3D,
                       derivativeBasisFunction1D,massMatrixLocal,B );

    // compute stifness  matrix ( durufle's optimization)
    int p=Qk.gradPhiGradPhi( numberOfPointsPerElement, order, weights3D, derivativeBasisFunction1D, B, pnLocal, R, Y );

    //compute global mass Matrix and global stiffness vector
    for( int i=0; i<numberOfPointsPerElement; i++ )
    {
      int gIndex=globalNodesList(e,i);
      massMatrixLocal[i]/=(model[e]*model[e]);
      massMatrixGlobal[gIndex]+=massMatrixLocal[i];
      yGlobal[gIndex]+=Y[i];
    }
  }
  
  // update pressure
#pragma omp  for
  for( int i=0; i<numberOfInteriorNodes; i++ )
  {
    int I=listOfInteriorNodes[i];
    float tmp=timeSample*timeSample;
    pnGlobal[I][i1]=2*pnGlobal[I][i2]-pnGlobal[I][i1]-tmp*yGlobal[I]/massMatrixGlobal[I];
  }
  
}// end parallel region
}
