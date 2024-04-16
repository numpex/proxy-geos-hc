//************************************************************************
//   proxy application v.0.0.1
//
//  solverOMP.cpp: simple 2D acoustive wave equation solver
//
//  the solverOMP class is derived from the solverBase class
//  with the openMP implementation of the solver
//
//************************************************************************

#include "solverSEQUENTIAL.hpp"

// compute one step of the time dynamic wave equation solver

void solverSEQUENTIAL::computeOneStep( const int & timeStep,
                                const float & timeSample,
                                const int & order,
                                int & i1,
                                int & i2,
                                const int & numberOfRHS,
                                vectorInt & rhsElement,
                                arrayReal & rhsTerm,
                                arrayReal & pnGlobal)
{
  float B[36][4];
  float R[36];
  float massMatrixLocal[36];
  float pnLocal[36];
  float Y[36];

  float ds[6];
  float Sh[6];
  int numOfBasisFunctionOnFace[6];
  float Js[2][6];

  for( int i=0; i<numberOfNodes; i++ )
  {
    massMatrixGlobal[i]=0;
    yGlobal[i]=0;
  }

  // update pnGLobal with right hade side
  for( int i=0; i<numberOfRHS;i++)
  { 
    int nodeRHS=globalNodesList(rhsElement[i],0);
    pnGlobal(nodeRHS,i2)+=timeSample*timeSample*model[rhsElement[i]]*model[rhsElement[i]]*rhsTerm(i,timeStep);
  }
  // loop over mesh elements
  for( int e=0; e<numberOfElements; e++ )
  {
    for( int i=0; i<numberOfPointsPerElement; i++ )
    {
      int   localToGlobal=globalNodesList(e,i);
      pnLocal[i]=pnGlobal[localToGlobal[i]][i2];
    }
    // compute Jacobian, massMatrix and B
    int o=Qk.computeB( e,order,globalNodesList,globalNodesCoords,weights2D,
                       derivativeBasisFunction1D,massMatrixLocal,B );
    // compute stifness  matrix ( durufle's optimization)
    int p=Qk.gradPhiGradPhi( numberOfPointsPerElement, order, weights2D, derivativeBasisFunction1D, B, pnLocal, R, Y );
    // get pnGlobal to pnLocal
    //compute global mass Matrix and global stiffness vector
    for( int i=0; i<numberOfPointsPerElement; i++ )
    {
      int gIndex=localToGlobal[i];
      massMatrixLocal[i]/=(model[e]*model[e]);
      massMatrixGlobal[gIndex]+=massMatrixLocal[i];
      yGlobal[gIndex]+=Y[i];
    }
  }
  
  // update pressure
  for( int i=0; i<numberOfInteriorNodes; i++ )
  {
    int I=listOfInteriorNodes[i];
    float tmp=timeSample*timeSample;
    pnGlobal[I][i1]=2*pnGlobal[I][i2]-pnGlobal[I][i1]-tmp*yGlobal[I]/massMatrixGlobal[I];
  }
  //cout<<"pressure="<<pnGlobal[5][i1]<<endl;
  
  // damping terms
  for( int i=0; i<numberOfBoundaryNodes; i++ )
  {
    ShGlobal[i]=0;
  }

  // Note: this loop is data parallel.
  for( int iFace=0; iFace<numberOfBoundaryFaces; iFace++ )
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
      int gIndexFaceNode=localFaceNodeToGlobalFaceNode[iFace][i];
      Sh[i]=weights[i]*ds[i]/(model[faceInfos[iFace][0]]);
      ShGlobal[gIndexFaceNode]+=Sh[i];
    }
  }

  // update pressure @ boundaries;
  float tmp=timeSample*timeSample;
  for( int i=0; i<numberOfBoundaryNodes; i++ )
  {
    int I=listOfBoundaryNodes[i];
    float invMpSh=1/(massMatrixGlobal[I]+timeSample*ShGlobal[i]*0.5);
    float MmSh=massMatrixGlobal[I]-timeSample*ShGlobal[i]*0.5;
    pnGlobal[I][i1]=invMpSh*(2*massMatrixGlobal[I]*pnGlobal[I][i2]-MmSh*pnGlobal[I][i1]-tmp*yGlobal[I]);
  }
}
