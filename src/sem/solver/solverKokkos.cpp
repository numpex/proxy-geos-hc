//************************************************************************
//  SEM proxy application v.0.0.1
//
//  solverRaja.cpp: simple 2D acoustive wave equation solver
//
//  the solverRaja class is derived from the solverBase class
//  with the RAJA implementation of the solver
//
//************************************************************************

#include "solverKokkos.hpp"
#include <cstdio>

// compute one step of the time dynamic wave equation solver
void solverKokkos::computeOneStep(  const int & timeStep,
                                  const float & timeSample,
                                  const int & order,
                                  int & i1,
                                  int & i2,
                                  const int & numberOfRHS,
                                  vectorInt & rhsElement,
                                  arrayReal & rhsTerm,
                                  arrayReal & pnGlobal,
                                  simpleMesh mesh,
                                  QkGL Qk )
{


  Kokkos::parallel_for(numberOfNodes,KOKKOS_CLASS_LAMBDA (const int i)
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
 

  // loop over mesh elements
  Kokkos::parallel_for(numberOfElements,KOKKOS_CLASS_LAMBDA (const int e)
  {
    int nPointsPerElement=(order+1)*(order+1);
    // start parallel section
    int  localToGlobal[36];
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
    // extract global coordinates of element e
    // get local to global indexes of nodes of element e
    //int i=mesh.localToGlobalNodes( e, numberOfPointsPerElement, globalNodesList, localToGlobal );
    //get global coordinates Xi of element e
    //mesh.getXi( nPointsPerElement, globalNodesCoords, localToGlobal, Xi );
    for( int i=0; i<nPointsPerElement; i++ )
    {
      localToGlobal[i]=globalNodesList(e,i);
      Xi[i][0]=globalNodesCoords(localToGlobal[i],0);
      Xi[i][1]=globalNodesCoords(localToGlobal[i],1);
    }

    // compute jacobian Matrix
    //Qk.computeJacobianMatrix( numberOfPointsPerElement, Xi,
    //                          derivativeBasisFunction2DX,
    //                          derivativeBasisFunction2DY,
    //                          jacobianMatrix );
    for( int i=0; i<nPointsPerElement; i++ )
    {
      jacobianMatrix[i][0]=0;
      jacobianMatrix[i][1]=0;
      jacobianMatrix[i][2]=0;
      jacobianMatrix[i][3]=0;
      for( int j=0; j<nPointsPerElement; j++ )
      {
        jacobianMatrix[i][0]+=Xi[j][0]*derivativeBasisFunction2DX(j,i);
        jacobianMatrix[i][1]+=Xi[j][0]*derivativeBasisFunction2DY(j,i);
        jacobianMatrix[i][2]+=Xi[j][1]*derivativeBasisFunction2DX(j,i);
        jacobianMatrix[i][3]+=Xi[j][1]*derivativeBasisFunction2DY(j,i);
      }
    }

    // compute determinant of jacobian Matrix
    //Qk.computeDeterminantOfJacobianMatrix( numberOfPointsPerElement,
    //                                       jacobianMatrix,
    //                                       detJ );
    for( int i=0; i<nPointsPerElement; i++ )
    {
      detJ[i]=(jacobianMatrix[i][0]*jacobianMatrix[i][3]-jacobianMatrix[i][2]*jacobianMatrix[i][1]);
    }

    // compute inverse of Jacobian Matrix
    //Qk.computeInvJacobianMatrix( numberOfPointsPerElement,
    //                             jacobianMatrix,
    //                             detJ,
    //                             invJacobianMatrix );
    for( int i=0; i<nPointsPerElement; i++ )
    {
      invJacobianMatrix[i][0]=(jacobianMatrix[i][3]/detJ[i]);
      invJacobianMatrix[i][1]=(-jacobianMatrix[i][1]/detJ[i]);
      invJacobianMatrix[i][2]=(-jacobianMatrix[i][2]/detJ[i]);
      invJacobianMatrix[i][3]=(jacobianMatrix[i][0]/detJ[i]);
    }
                                 
    // compute transposed inverse of Jacobian Matrix
    //Qk.computeTranspInvJacobianMatrix( numberOfPointsPerElement,
    //                                   jacobianMatrix,
    //                                   detJ,
    //                                   transpInvJacobianMatrix );
    for( int i=0; i<nPointsPerElement; i++ )
    {
      transpInvJacobianMatrix[i][0]=(jacobianMatrix[i][3]/detJ[i]);
      transpInvJacobianMatrix[i][1]=(-jacobianMatrix[i][2]/detJ[i]);
      transpInvJacobianMatrix[i][2]=(-jacobianMatrix[i][1]/detJ[i]);
      transpInvJacobianMatrix[i][3]=(jacobianMatrix[i][0]/detJ[i]);
    }
                        
    // compute  geometrical transformation matrix
    //Qk.computeB( numberOfPointsPerElement, invJacobianMatrix, transpInvJacobianMatrix, detJ,B );
    for( int i=0; i<nPointsPerElement; i++ )
    {
      B[i][0]=(abs( detJ[i] )*(invJacobianMatrix[i][0]*transpInvJacobianMatrix[i][0]+
                               invJacobianMatrix[i][1]*transpInvJacobianMatrix[i][2]));
      B[i][1]=(abs( detJ[i] )*(invJacobianMatrix[i][0]*transpInvJacobianMatrix[i][1]+
                               invJacobianMatrix[i][1]*transpInvJacobianMatrix[i][3]));
      B[i][2]=(abs( detJ[i] )*(invJacobianMatrix[i][2]*transpInvJacobianMatrix[i][0]+
                               invJacobianMatrix[i][3]*transpInvJacobianMatrix[i][2]));
      B[i][3]=(abs( detJ[i] )*(invJacobianMatrix[i][2]*transpInvJacobianMatrix[i][1]+
                               invJacobianMatrix[i][3]*transpInvJacobianMatrix[i][3]));
  }

    // compute stifness and mass matrix ( durufle's optimization)
    //Qk.gradPhiGradPhi( numberOfPointsPerElement, order, weights2D, B, derivativeBasisFunction1D, R );
    for( int i=0;i<nPointsPerElement;i++)
    {
       for(int j=0; j<nPointsPerElement;j++)
       {
         R[j][i]=0;
       }
    }
    // B11
    for( int i1=0; i1<order+1; i1++ )
    {
      for( int i2=0; i2<order+1; i2++ )
      {
        int i=i1+i2*(order+1);
        for( int j1=0; j1<order+1; j1++ )
        {
          int j=j1+i2*(order+1);
          for( int m=0; m<order+1; m++ )
          {
            R[j][i]+=weights2D[m+i2*(order+1)]*(B[m+i2*(order+1)][0]*derivativeBasisFunction1D(i1,m)*derivativeBasisFunction1D(j1,m));
          }
        }
      }
    }
    // B21
    for( int i1=0; i1<order+1; i1++ )
    {
      for( int i2=0; i2<order+1; i2++ )
      {
        int i=i1+i2*(order+1);
        for( int j1=0; j1<order+1; j1++ )
        {
          for( int j2=0; j2<order+1; j2++ )
          {
            int j=j1+j2*(order+1);
            R[j][i]+=weights2D[i1+j2*(order+1)]*(B[i1+j2*(order+1)][1]*derivativeBasisFunction1D(i2,j2)*derivativeBasisFunction1D(j1,i1));
          }
        }
      }
    }
    // B12
    for( int i1=0; i1<order+1; i1++ )
    {
      for( int i2=0; i2<order+1; i2++ )
      {
        int i=i1+i2*(order+1);
        for( int j1=0; j1<order+1; j1++ )
        {
          for( int j2=0; j2<order+1; j2++ )
          {
            int j=j1+j2*(order+1);
            R[j][i]+=weights2D[i2+j1*(order+1)]*(B[i2+j1*(order+1)][2]*derivativeBasisFunction1D(i1,j1)*derivativeBasisFunction1D(j2,i2));
          }
        }
      }
    }
    // B22
    for( int i1=0; i1<order+1; i1++ )
    {
      for( int i2=0; i2<order+1; i2++ )
      {
        int i=i1+i2*(order+1);
        for( int j2=0; j2<order+1; j2++ )
        {
          int j=i1+j2*(order+1);
          for( int n=0; n<order+1; n++ )
          {
          R[j][i]+=weights2D[i1+n*(order+1)]*(B[i1+n*(order+1)][3]*derivativeBasisFunction1D(i2,n)*derivativeBasisFunction1D(j2,n));
          }
        }
      }
    }
    // compute local mass matrix ( used optimez version)
    //Qk.phiIphiJ( numberOfPointsPerElement, weights2D, detJ, massMatrixLocal );
    for( int i=0; i<nPointsPerElement; i++ )
    {
      massMatrixLocal[i]=weights2D[i]*abs( detJ[i] );
    }

    // get pnGlobal to pnLocal
    for( int i=0; i<nPointsPerElement; i++ )
    {
      massMatrixLocal[i]/=(model[e]*model[e]);
      pnLocal[i]=pnGlobal(localToGlobal[i],i2);
    }

    // compute Y=R*pnLocal
    for( int i=0; i<nPointsPerElement; i++ )
    {
      Y[i]=0;
      for( int j=0; j<nPointsPerElement; j++ )
      {
        Y[i]+=R[i][j]*pnLocal[j];
      }
    }

    //compute gloval mass Matrix and global stiffness vector
    for( int i=0; i<nPointsPerElement; i++ )
    {
      int gIndex=localToGlobal[i];
      //massMatrixGlobal[gIndex]+=massMatrixLocal(threadId,i)
      //yGlobal[gIndex]+=Y(threadId,i);
      Kokkos::atomic_add(&massMatrixGlobal[gIndex],massMatrixLocal[i]);
      Kokkos::atomic_add(&yGlobal[gIndex],Y[i]);
    } 
  } );

  // update pressure
  Kokkos::parallel_for(numberOfInteriorNodes,KOKKOS_CLASS_LAMBDA (const int i)
  {
    int I=listOfInteriorNodes[i];
    float tmp=timeSample*timeSample;
    pnGlobal(I,i1)=2*pnGlobal(I,i2)-pnGlobal(I,i1)-tmp*yGlobal[I]/massMatrixGlobal[I];
  } );

  Kokkos::parallel_for(numberOfBoundaryNodes,KOKKOS_CLASS_LAMBDA (const int i)
  {
    ShGlobal[i]=0;
  } );
  // Note: this loop is data parallel.
  Kokkos::parallel_for(numberOfBoundaryFaces,KOKKOS_CLASS_LAMBDA (const int iFace)
  {
    //get ds
    float ds[6];
    float Sh[6];
    int numOfBasisFunctionOnFace[6];
    float Js[2][6];

    //Qk.computeDs( iFace, order, faceInfos,numOfBasisFunctionOnFace,
    //              Js, globalNodesCoords, derivativeBasisFunction2DX,
    //              derivativeBasisFunction2DY,
    //              ds );
    int face=faceInfos(iFace,1);
    // get basis functions on Boundary faces
    switch( face )
    {
      case 0:     // left
        for( int i=0; i<order+1; i++ )
        {
          numOfBasisFunctionOnFace[i]=i*(order+1);
        }
        break;
      case 1:     // bottom
        for( int i=0; i<order+1; i++ )
        {
          numOfBasisFunctionOnFace[i]=i;
        }
        break;
      case 2:         //right
        for( int i=0; i<order+1; i++ )
        {
          numOfBasisFunctionOnFace[i]=order+i*(order+1);
        }
        break;
      case 3:         //top
        for( int i=0; i<order+1; i++ )
        {
          numOfBasisFunctionOnFace[i]=i+order*(order+1);
        }
        break;
      default:
        //cout<<"error in element flag, should be set to: 0, 1, 2, 3"<<endl;
        break;
    }
    // compute ds
    for( int j=0; j<order+1; j++ )
    {
      Js[0][j]=0;    // x
      Js[1][j]=0;    // y
      for( int i=0; i<order+1; i++ )
      {
        float xi=globalNodesCoords(faceInfos(iFace,2+i),0);
        float yi=globalNodesCoords(faceInfos(iFace,2+i),1);
        if( face==0 || face==2 )
        {
          Js[0][j]+=derivativeBasisFunction2DY(numOfBasisFunctionOnFace[i],numOfBasisFunctionOnFace[j])*xi;
          Js[1][j]+=derivativeBasisFunction2DY(numOfBasisFunctionOnFace[i],numOfBasisFunctionOnFace[j])*yi;
        }
        if( face==1 || face==3 )
        {
          Js[0][j]+=derivativeBasisFunction2DX(numOfBasisFunctionOnFace[i],numOfBasisFunctionOnFace[j])*xi;
          Js[1][j]+=derivativeBasisFunction2DX(numOfBasisFunctionOnFace[i],numOfBasisFunctionOnFace[j])*yi;
        }
      }
      ds[j]=sqrt( Js[0][j]*Js[0][j]+Js[1][j]*Js[1][j] );
    }
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
  Kokkos::parallel_for(numberOfBoundaryNodes,KOKKOS_CLASS_LAMBDA (const int i)
  {
    int I=listOfBoundaryNodes[i];
    float invMpSh=1/(massMatrixGlobal[I]+timeSample*ShGlobal[i]*0.5);
    float MmSh=massMatrixGlobal[I]-timeSample*ShGlobal[i]*0.5;
    pnGlobal(I,i1)=invMpSh*(2*massMatrixGlobal[I]*pnGlobal(I,i2)-MmSh*pnGlobal(I,i1)-tmp*yGlobal[I]);
  } );
  Kokkos::fence();
}
