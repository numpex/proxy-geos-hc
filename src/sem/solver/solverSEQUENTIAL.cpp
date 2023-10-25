//************************************************************************
//  SEM proxy application v.0.0.1
//
//  solverSEQUENTIAL.cpp: simple 2D acoustive wave equation solver
//
//  the solverSEQUENTIAL class is derived from the solverBase class
//  with the openSEQUENTIALimplementation of the solver
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
  int localToGlobal[36];
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
    // extract global coordinates of element e
    // get local to global indexes of nodes of element e
    int a=mesh.localToGlobalNodes( e, numberOfPointsPerElement, globalNodesList, localToGlobal );
    
    //get global coordinates Xi of element e
    int b=mesh.getXi( numberOfPointsPerElement, globalNodesCoords, localToGlobal, Xi );
    //if(e<10)printf("%lf %lf %lf %lf %lf %lf %lf %lf\n",Xi[0][0],Xi[0][1],Xi[1][0],Xi[1][1],Xi[2][0],Xi[2][1],Xi[3][0],Xi[3][1]);

    // compute jacobian Matrix
    int c=Qk.computeJacobianMatrix( numberOfPointsPerElement, Xi,
                              derivativeBasisFunction2DX,
                              derivativeBasisFunction2DY,
                              jacobianMatrix );

    // compute determinant of jacobian Matrix
    int d=Qk.computeDeterminantOfJacobianMatrix( numberOfPointsPerElement,
                                           jacobianMatrix,
                                           detJ );
    // compute inverse of Jacobian Matrix
    int f=Qk.computeInvJacobianMatrix( numberOfPointsPerElement,
                                 jacobianMatrix,
                                 detJ,
                                 invJacobianMatrix );
                                 
    // compute transposed inverse of Jacobian Matrix
    int g=Qk.computeTranspInvJacobianMatrix( numberOfPointsPerElement,
                                       jacobianMatrix,
                                       detJ,
                                       transpInvJacobianMatrix );
                        
    // compute  geometrical transformation matrix
    int h=Qk.computeB( numberOfPointsPerElement, invJacobianMatrix, transpInvJacobianMatrix, detJ,B );

    // compute stifness and mass matrix ( durufle's optimization)
    int i=Qk.gradPhiGradPhi( numberOfPointsPerElement, order, weights2D, B, derivativeBasisFunction1D, R );

    // compute local mass matrix ( used optimez version)
    int j=Qk.phiIphiJ( numberOfPointsPerElement, weights2D, detJ, massMatrixLocal );

    // get pnGlobal to pnLocal
    for( int i=0; i<numberOfPointsPerElement; i++ )
    {
      massMatrixLocal[i]/=(model[e]*model[e]);
      pnLocal[i]=pnGlobal[localToGlobal[i]][i2];
      //cout<<"element "<<e<<" massM "<<i<<" "<<massMatrixLocal[i];
      //cout<<"locToGLob "<<i<<" "<<localToGlobal[i]<<" "<<i2;
      //cout<<" pnLocal "<<i<<" "<<pnLocal[i]<<endl;
    }

    // compute Y=R*pnLocal
    for( int i=0; i<numberOfPointsPerElement; i++ )
    {
      Y[i]=0;
      for( int j=0; j<numberOfPointsPerElement; j++ )
      {
        Y[i]+=R[i][j]*pnLocal[j];
      }
           //cout<<" element "<<e<<" Y "<<i<<" "<<Y[i];
    }
    //cout<<endl;

    //compute gloval mass Matrix and global stiffness vector
    for( int i=0; i<numberOfPointsPerElement; i++ )
    {
      int gIndex=localToGlobal[i];
      massMatrixGlobal[gIndex]+=massMatrixLocal[i];
      yGlobal[gIndex]+=Y[i];
      //cout<<"element "<<e<<" massG "<<gIndex<<" "<<massMatrixGlobal[gIndex];
      //cout<<" yGLobal "<<gIndex<<" "<<yGlobal[gIndex]<<endl;
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
