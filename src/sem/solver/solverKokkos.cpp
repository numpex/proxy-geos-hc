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
      int e=listOfElementsByColor(color,eColor);
      // extract global coordinates of element e
      // get local to global indexes of nodes of element e
      int i=mesh.localToGlobalNodes( e, numberOfPointsPerElement, globalNodesList, localToGlobal );
      //get global coordinates Xi of element e
      int j=mesh.getXi( numberOfPointsPerElement, globalNodesCoords, localToGlobal, Xi );
      // compute jacobian Matrix
      int k=Qk.computeJacobianMatrix( numberOfPointsPerElement, Xi,
                                      derivativeBasisFunction2DX,
                                      derivativeBasisFunction2DY,
                                      jacobianMatrix );
      // compute determinant of jacobian Matrix
      int l=Qk.computeDeterminantOfJacobianMatrix( numberOfPointsPerElement,
                                                   jacobianMatrix,
						   detJ);
      // compute inverse of Jacobian Matrix
      int m=Qk.computeInvJacobianMatrix( numberOfPointsPerElement,
                                         jacobianMatrix,
                                         detJ,
                                         invJacobianMatrix );
      // compute transposed inverse of Jacobian Matrix
      int n=Qk.computeTranspInvJacobianMatrix( numberOfPointsPerElement,
                                               jacobianMatrix,
                                               detJ,
                                               transpInvJacobianMatrix );
      // compute  geometrical transformation matrix
      int p=Qk.computeB( numberOfPointsPerElement, invJacobianMatrix, transpInvJacobianMatrix, detJ,B );
      // compute stifness and mass matrix ( durufle's optimization)
      int q=Qk.gradPhiGradPhi( numberOfPointsPerElement, order, weights2D, B, derivativeBasisFunction1D, R );
      // compute local mass matrix ( used optimez version)
      int r=Qk.phiIphiJ( numberOfPointsPerElement, weights2D, detJ, massMatrixLocal );
      // get pnGlobal to pnLocal
      for( int i=0; i<numberOfPointsPerElement; i++ )
      {
        massMatrixLocal[i]/=(model[e]*model[e]);
        pnLocal[i]=pnGlobal(localToGlobal[i],i2);
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
        //Kokkos::atomic_add(&massMatrixGlobal[gIndex],massMatrixLocal[i]);
        //Kokkos::atomic_add(&yGlobal[gIndex],Y[i]);
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
  } );
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
