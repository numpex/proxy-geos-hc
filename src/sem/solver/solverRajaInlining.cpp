//************************************************************************
//  SEM proxy application v.0.0.1
//
//  solverRaja.cpp: simple 2D acoustive wave equation solver
//
//  the solverRaja class is derived from the solverBase class
//  with the RAJA implementation of the solver
//
//************************************************************************

#include "solverRaja.hpp"

// compute one step of the time dynamic wave equation solver
void solverRaja::computeOneStep(  const int & timeStep,
                                  const float & timeSample,
                                  const int & order,
                                  int & i1,
                                  int & i2,
                                  const int & numberOfRHS,
                                  vectorInt & rhsElement,
                                  arrayReal & rhsTerm,
                                  arrayReal & pnGlobal)
{

  vectorIntView d_rhsElement=rhsElement.toView();
  arrayRealView d_rhsTerm=rhsTerm.toView();
  arrayRealView d_pnGlobal=pnGlobal.toView();

  //shared arrays
  arrayIntView d_globalNodesList=globalNodesList.toView();
  arrayRealView d_globalNodesCoords=globalNodesCoords.toView();
  vectorIntView d_listOfInteriorNodes=listOfInteriorNodes.toView();
  vectorIntView d_listOfBoundaryNodes=listOfBoundaryNodes.toView();
  arrayIntView d_faceInfos=faceInfos.toView();
  arrayIntView d_localFaceNodeToGlobalFaceNode=localFaceNodeToGlobalFaceNode.toView();
  
  // get model
  vectorRealView d_model=model.toView();

  // get quadrature points and weights
  vectorDoubleView d_quadraturePoints=quadraturePoints.toView();
  vectorDoubleView d_weights=weights.toView();
  vectorDoubleView d_weights2D=weights2D.toView();

  // get basis function and corresponding derivatives
  arrayDoubleView d_basisFunction1D=basisFunction1D.toView();
  arrayDoubleView d_derivativeBasisFunction1D=derivativeBasisFunction1D.toView();
  arrayDoubleView d_basisFunction2D=basisFunction2D.toView();
  arrayDoubleView d_derivativeBasisFunction2DX=derivativeBasisFunction2DX.toView();
  arrayDoubleView d_derivativeBasisFunction2DY=derivativeBasisFunction2DY.toView();

  //shared arrays
  vectorDoubleView d_massMatrixGlobal=massMatrixGlobal.toView();
  vectorDoubleView d_yGlobal=yGlobal.toView();
  vectorRealView d_ShGlobal=ShGlobal.toView();

  RAJA::forall< exec_policy >( RAJA::RangeSegment( 0, numberOfNodes ),  [=] LVARRAY_HOST_DEVICE  ( int i ) {
    d_massMatrixGlobal[i]=0;
    d_yGlobal[i]=0;
  } );
  // update pnGLobal with right hade side
  RAJA::forall< exec_policy >( RAJA::RangeSegment( 0, numberOfRHS ), [=] LVARRAY_HOST_DEVICE  ( int i ) 
  {
    int nodeRHS=d_globalNodesList(d_rhsElement[i],0);
    d_pnGlobal(nodeRHS,i2)+=timeSample*timeSample*d_model[d_rhsElement[i]]*d_model[d_rhsElement[i]]*d_rhsTerm(i,timeStep);
  });
  int numberOfPointsPerElement=(order+1)*(order+1);
  // loop over mesh elements
  RAJA::forall< exec_policy >( RAJA::RangeSegment( 0, numberOfElements ), [=] LVARRAY_HOST_DEVICE ( int e )
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
    //get global coordinates Xi of element e
    for( int i=0; i<nPointsPerElement; i++ )
    {
      localToGlobal[i]=d_globalNodesList(e,i);
      Xi[i][0]=d_globalNodesCoords(localToGlobal[i],0);
      Xi[i][1]=d_globalNodesCoords(localToGlobal[i],1);
    }
    // compute jacobian Matrix
    for( int i=0; i<nPointsPerElement; i++ )
    {
      jacobianMatrix[i][0]=0;
      jacobianMatrix[i][1]=0;
      jacobianMatrix[i][2]=0;
      jacobianMatrix[i][3]=0;
      for( int j=0; j<nPointsPerElement; j++ )
      {
        jacobianMatrix[i][0]+=Xi[j][0]*d_derivativeBasisFunction2DX(j,i);
        jacobianMatrix[i][1]+=Xi[j][0]*d_derivativeBasisFunction2DY(j,i);
        jacobianMatrix[i][2]+=Xi[j][1]*d_derivativeBasisFunction2DX(j,i);
        jacobianMatrix[i][3]+=Xi[j][1]*d_derivativeBasisFunction2DY(j,i);
      }
    }
    // compute determinant of jacobian Matrix
    for( int i=0; i<nPointsPerElement; i++ )
    {
      detJ[i]=(jacobianMatrix[i][0]*jacobianMatrix[i][3]-jacobianMatrix[i][2]*jacobianMatrix[i][1]);
    }
    // compute inverse of Jacobian Matrix
    for( int i=0; i<nPointsPerElement; i++ )
    {
      invJacobianMatrix[i][0]=(jacobianMatrix[i][3]/detJ[i]);
      invJacobianMatrix[i][1]=(-jacobianMatrix[i][1]/detJ[i]);
      invJacobianMatrix[i][2]=(-jacobianMatrix[i][2]/detJ[i]);
      invJacobianMatrix[i][3]=(jacobianMatrix[i][0]/detJ[i]);
    }
    // compute transposed inverse of Jacobian Matrix
    for( int i=0; i<nPointsPerElement; i++ )
    {
      transpInvJacobianMatrix[i][0]=(jacobianMatrix[i][3]/detJ[i]);
      transpInvJacobianMatrix[i][1]=(-jacobianMatrix[i][2]/detJ[i]);
      transpInvJacobianMatrix[i][2]=(-jacobianMatrix[i][1]/detJ[i]);
      transpInvJacobianMatrix[i][3]=(jacobianMatrix[i][0]/detJ[i]);
    }
    // compute  geometrical transformation matrix
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
            R[j][i]+=d_weights2D[m+i2*(order+1)]*(B[m+i2*(order+1)][0]*d_derivativeBasisFunction1D(i1,m)*d_derivativeBasisFunction1D(j1,m));
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
            R[j][i]+=d_weights2D[i1+j2*(order+1)]*(B[i1+j2*(order+1)][1]*d_derivativeBasisFunction1D(i2,j2)*d_derivativeBasisFunction1D(j1,i1));
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
            R[j][i]+=d_weights2D[i2+j1*(order+1)]*(B[i2+j1*(order+1)][2]*d_derivativeBasisFunction1D(i1,j1)*d_derivativeBasisFunction1D(j2,i2));
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
          R[j][i]+=d_weights2D[i1+n*(order+1)]*(B[i1+n*(order+1)][3]*d_derivativeBasisFunction1D(i2,n)*d_derivativeBasisFunction1D(j2,n));
          }
        }
      }
    }
    // compute local mass matrix ( used optimez version)
    for( int i=0; i<nPointsPerElement; i++ )
    {
      massMatrixLocal[i]=d_weights2D[i]*abs( detJ[i] );
    }
    // get pnGlobal to pnLocal
    for( int i=0; i<nPointsPerElement; i++ )
    {
      massMatrixLocal[i]/=(d_model[e]*d_model[e]);
      pnLocal[i]=d_pnGlobal(localToGlobal[i],i2);
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
      RAJA::atomicAdd< atomic_policy >(&d_massMatrixGlobal[gIndex],massMatrixLocal[i]);
      RAJA::atomicAdd< atomic_policy>(&d_yGlobal[gIndex],Y[i]);
    } 
  } );
  // update pressure
  RAJA::forall< exec_policy>( RAJA::RangeSegment( 0, numberOfInteriorNodes ), [=] LVARRAY_HOST_DEVICE ( int i ) {
    int I=d_listOfInteriorNodes[i];
    float tmp=timeSample*timeSample;
    d_pnGlobal[I][i1]=2*d_pnGlobal[I][i2]-d_pnGlobal[I][i1]-tmp*d_yGlobal[I]/d_massMatrixGlobal[I];
  } );

  RAJA::forall< exec_policy>( RAJA::RangeSegment( 0, numberOfBoundaryNodes ), [=] LVARRAY_HOST_DEVICE ( int i ) {
    d_ShGlobal[i]=0;
  } );
  RAJA::forall< exec_policy >( RAJA::RangeSegment( 0, numberOfBoundaryFaces ), [=] LVARRAY_HOST_DEVICE ( int iFace ){
    //get ds
    float ds[6];
    float Sh[6];
    int numOfBasisFunctionOnFace[6];
    float Js[2][6];
    int face=d_faceInfos(iFace,1);
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
        float xi=d_globalNodesCoords(d_faceInfos(iFace,2+i),0);
        float yi=d_globalNodesCoords(d_faceInfos(iFace,2+i),1);
        if( face==0 || face==2 )
        {
          Js[0][j]+=d_derivativeBasisFunction2DY(numOfBasisFunctionOnFace[i],numOfBasisFunctionOnFace[j])*xi;
          Js[1][j]+=d_derivativeBasisFunction2DY(numOfBasisFunctionOnFace[i],numOfBasisFunctionOnFace[j])*yi;
        }
        if( face==1 || face==3 )
        {
          Js[0][j]+=d_derivativeBasisFunction2DX(numOfBasisFunctionOnFace[i],numOfBasisFunctionOnFace[j])*xi;
          Js[1][j]+=d_derivativeBasisFunction2DX(numOfBasisFunctionOnFace[i],numOfBasisFunctionOnFace[j])*yi;
        }
      }
      ds[j]=sqrt( Js[0][j]*Js[0][j]+Js[1][j]*Js[1][j] );
    }
    //compute Sh and ShGlobal
    for( int i=0; i<order+1; i++ )
    {
      int gIndexFaceNode=d_localFaceNodeToGlobalFaceNode(iFace,i);
      Sh[i]=d_weights[i]*ds[i]/(d_model[d_faceInfos(iFace,0)]);
      RAJA::atomicAdd< atomic_policy >(&d_ShGlobal[gIndexFaceNode],Sh[i]);
    }
  } );
  // update pressure @ boundaries;
  float tmp=timeSample*timeSample;
  RAJA::forall< exec_policy >( RAJA::RangeSegment( 0, numberOfBoundaryNodes ), [=] LVARRAY_HOST_DEVICE ( int i ) {
    int I=d_listOfBoundaryNodes[i];
    float invMpSh=1/(d_massMatrixGlobal[I]+timeSample*d_ShGlobal[i]*0.5);
    float MmSh=d_massMatrixGlobal[I]-timeSample*d_ShGlobal[i]*0.5;
    d_pnGlobal[I][i1]=invMpSh*(2*d_massMatrixGlobal[I]*d_pnGlobal[I][i2]-MmSh*d_pnGlobal[I][i1]-tmp*d_yGlobal[I]);
  } );
  if(timeStep%900==0)
  {
     int nodeRHS=d_globalNodesList(d_rhsElement[0],0);
     RAJA::forall< RAJA::seq_exec>( RAJA::RangeSegment( 0, numberOfNodes ), [=] LVARRAY_HOST_DEVICE ( int i )
     {
	  pnGlobal(nodeRHS,i1)=d_pnGlobal(nodeRHS,i1);
     });
  }
}
