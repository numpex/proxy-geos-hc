//************************************************************************
//  SEM proxy application v.0.0.1
//
//  solverOMP.cpp: simple 2D acoustive wave equation solver
//
//  the solverOMP class is derived from the solverBase class
//  with the openMP implementation of the solver
//
//************************************************************************

#include "solverOMP.hpp"

// compute one step of the time dynamic wave equation solver

void solverOMP::computeOneStep( const float & timeSample,
                                const int & order,
                                int & i1,
                                int & i2,
                                arrayReal & pnGlobal,
                                simpleMesh mesh,
                                QkGL Qk )
{
  

#pragma omp parallel
  { // start parallel section
    vectorInt localToGlobal(numberOfPointsPerElement);
    arrayDouble Xi(numberOfPointsPerElement, 2 );

    arrayDouble jacobianMatrix(4, numberOfPointsPerElement);
    vectorDouble detJ(numberOfPointsPerElement);
    arrayDouble invJacobianMatrix(4, numberOfPointsPerElement);
    arrayDouble transpInvJacobianMatrix(4, numberOfPointsPerElement);

    arrayDouble B(4, numberOfPointsPerElement);
    arrayDouble R(numberOfPointsPerElement, numberOfPointsPerElement);

    vectorDouble massMatrixLocal(numberOfPointsPerElement);

    vectorReal pnLocal( numberOfPointsPerElement );
    vectorReal Y( numberOfPointsPerElement );

    vectorReal ds( order+1 );
    vectorReal Sh( order+1 );
    vectorInt numOfBasisFunctionOnFace( order+1 );
    arrayReal Js( 2, order+1 );

  #pragma omp for 
  for( int i=0; i<numberOfNodes; i++ )
  {
    massMatrixGlobal[i]=0;
    yGlobal[i]=0;
  }

  // loop over mesh elements
#pragma omp for
  for( int e=0; e<numberOfElements; e++ )
  {

    // extract global coordinates of element e
    // get local to global indexes of nodes of element e
    int a=mesh.localToGlobalNodes( e, numberOfPointsPerElement, globalNodesList, localToGlobal );

    //get global coordinates Xi of element e
    int b=mesh.getXi( numberOfPointsPerElement, globalNodesCoords, localToGlobal, Xi );

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
      int gIndex=localToGlobal[i];
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
  //cout<<"pressure="<<pnGlobal[5][i1]<<endl;
  
  // damping terms
#pragma omp  for
  for( int i=0; i<numberOfBoundaryNodes; i++ )
  {
    ShGlobal[i]=0;
  }

  // Note: this loop is data parallel.
#pragma omp  for
  for( int iFace=0; iFace<numberOfBoundaryFaces; iFace++ )
  {
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
#pragma omp  for
  for( int i=0; i<numberOfBoundaryNodes; i++ )
  {
    int I=listOfBoundaryNodes[i];
    float invMpSh=1/(massMatrixGlobal[I]+timeSample*ShGlobal[i]*0.5);
    float MmSh=massMatrixGlobal[I]-timeSample*ShGlobal[i]*0.5;
    pnGlobal[I][i1]=invMpSh*(2*massMatrixGlobal[I]*pnGlobal[I][i2]-MmSh*pnGlobal[I][i1]-tmp*yGlobal[I]);
  }
  }// end parallel region
}
