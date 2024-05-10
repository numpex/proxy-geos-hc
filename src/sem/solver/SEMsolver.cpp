//************************************************************************
//   proxy application v.0.0.1
//
//  SEMsolver.cpp: simple 2D acoustive wave equation solver
//
//  the SEMsolver class servers as a base class for the SEM solver
//
//************************************************************************

#include "SEMsolver.hpp"

void SEMsolver::computeFEInit( SEMmeshinfo & myMeshinfo, SEMmesh mesh )
{
  order = myMeshinfo.myOrderNumber;
  allocateFEarrays( myMeshinfo );
  initFEarrays( myMeshinfo, mesh );

}


// compute one step of the time dynamic wave equation solver
#ifndef USE_SEM_INLINE
void SEMsolver::computeOneStepNoInline( 
                                const int & timeSample,
                                const int & order,
                                const int & nPointsPerElement,
                                const int & i1,
                                const int & i2,
                                SEMmeshinfo & myMeshinfo,
                                const arrayReal & RHS_Term,
                                arrayReal const & PN_Global,
                                const vectorInt & RHS_Element )
{
  CREATEVIEWS

  // update pressure @ boundaries;
  LOOPHEAD( myMeshinfo.numberOfNodes, i )
  massMatrixGlobal[i]=0;
  yGlobal[i]=0;
  LOOPEND

  // update pnGLobal with right hade side
  LOOPHEAD( myMeshinfo.myNumberOfRHS, i )
  int nodeRHS=globalNodesList( rhsElement[i], 0 );
  pnGlobal( nodeRHS, i2 )+=myMeshinfo.myTimeStep*myMeshinfo.myTimeStep*model[rhsElement[i]]*model[rhsElement[i]]*rhsTerm( i, timeSample );
  LOOPEND

  LOOPHEAD( myMeshinfo.numberOfElements, elementNumber )
  // start parallel section
  float B[ROW][COL];
  float R[ROW];
  float massMatrixLocal[ROW];
  float pnLocal[ROW];
  float Y[ROW];

  // get pnGlobal to pnLocal
  for( int i=0; i<nPointsPerElement; i++ )
  {
    int localToGlobal=globalNodesList( elementNumber, i );
    pnLocal[i]=pnGlobal( localToGlobal, i2 );
  }

  // compute Jacobian, massMatrix and B
  myQk.computeB( elementNumber, order, weights, globalNodesList, globalNodesCoords, derivativeBasisFunction1D, massMatrixLocal, B );

  // compute stifness  matrix ( durufle's optimization)
  myQk.gradPhiGradPhi( nPointsPerElement, order, weights, derivativeBasisFunction1D, B, pnLocal, R, Y );


  //compute gloval mass Matrix and global stiffness vector
  for( int i=0; i<nPointsPerElement; i++ )
  {
    int gIndex=globalNodesList( elementNumber, i );
    massMatrixLocal[i]/=(model[elementNumber]*model[elementNumber]);
    ATOMICADD( massMatrixGlobal[gIndex], massMatrixLocal[i] );
    ATOMICADD( yGlobal[gIndex], Y[i] );
  }

  LOOPEND

  // update pressure
  LOOPHEAD( myMeshinfo.numberOfInteriorNodes, i )
  int I=listOfInteriorNodes[i];
  pnGlobal( I, i1 )=2*pnGlobal( I, i2 )-pnGlobal( I, i1 )-myMeshinfo.myTimeStep*myMeshinfo.myTimeStep*yGlobal[I]/massMatrixGlobal[I];
  LOOPEND

  #ifdef SEM2D
  {
    LOOPHEAD( myMeshinfo.numberOfBoundaryNodes, i )
    ShGlobal[i]=0;
    LOOPEND

    LOOPHEAD( myMeshinfo.numberOfBoundaryFaces, iFace )
    //get ds
    float ds[6];
    float Sh[6];
    float Js[2][6];

    // compute ds
    myQk.computeDs( iFace, order, faceInfos, (order+1)*(order+1), Js, globalNodesCoords, derivativeBasisFunction1D, ds );

    //compute Sh and ShGlobal
    for( int i=0; i<order+1; i++ )
    {
      int gIndexFaceNode=localFaceNodeToGlobalFaceNode( iFace, i );
      Sh[i]=weights[i]*ds[i]/(model[faceInfos( iFace, 0 )]);
      ATOMICADD( ShGlobal[gIndexFaceNode], Sh[i] );
    }
    LOOPEND

    LOOPHEAD( myMeshinfo.numberOfBoundaryNodes, i )
    int I=listOfBoundaryNodes[i];
    float invMpSh=1/(massMatrixGlobal[I]+myMeshinfo.myTimeStep*ShGlobal[i]*0.5);
    float MmSh=massMatrixGlobal[I]-myMeshinfo.myTimeStep*ShGlobal[i]*0.5;
    pnGlobal( I, i1 )=invMpSh*(2*massMatrixGlobal[I]*pnGlobal( I, i2 )-MmSh*pnGlobal( I, i1 )-myMeshinfo.myTimeStep*myMeshinfo.myTimeStep*yGlobal[I]);
    LOOPEND
  }
  #endif
  FENCE
}

#else

// compute one step of the time dynamic wave equation solver
void SEMsolver::computeOneStepInline( 
                                const int & timeSample,
                                const int & order,
                                const int & nPointsPerElement,
                                const int & i1,
                                const int & i2,
                                SEMmeshinfo & myMeshinfo,
                                const arrayReal & RHS_Term,
                                arrayReal const & PN_Global,
                                const vectorInt & RHS_Element )
{
  CREATEVIEWS

  // update pressure @ boundaries;
  LOOPHEAD( myMeshinfo.numberOfNodes, i )
  massMatrixGlobal[i]=0;
  yGlobal[i]=0;
  LOOPEND

  // update pnGLobal with right hade side
  LOOPHEAD( myMeshinfo.myNumberOfRHS, i )
  int nodeRHS=globalNodesList( rhsElement[i], 0 );
  pnGlobal( nodeRHS, i2 )+=myMeshinfo.myTimeStep*myMeshinfo.myTimeStep*model[rhsElement[i]]*model[rhsElement[i]]*rhsTerm( i, timeSample );
  LOOPEND

  LOOPHEAD( myMeshinfo.numberOfElements, elementNumber )
  // start parallel section
  float B[ROW][COL];
  float R[ROW];
  float massMatrixLocal[ROW];
  float pnLocal[ROW];
  float Y[ROW];

  // get pnGlobal to pnLocal
  for( int i=0; i<nPointsPerElement; i++ )
  {
    int localToGlobal=globalNodesList( elementNumber, i );
    pnLocal[i]=pnGlobal( localToGlobal, i2 );
  }

  #ifdef SEM2D
  // compute Jacobian, massMatrix and B
  for( int i2=0; i2<order+1; i2++ )
  {
    for( int i1=0; i1<order+1; i1++ )
    {
      // compute jacobian matrix
      double jac0=0;
      double jac1=0;
      double jac2=0;
      double jac3=0;
      int i=i1+i2*(order+1);
      for( int j1=0; j1<order+1; j1++ )
      {
        int j=j1+i2*(order+1);
        int localToGlobal=globalNodesList( elementNumber, j );
        double X=globalNodesCoords( localToGlobal, 0 );
        double Y=globalNodesCoords( localToGlobal, 1 );
        jac0+=X*derivativeBasisFunction1D( j1, i1 );
        jac2+=Y*derivativeBasisFunction1D( j1, i1 );
      }
      for( int j2=0; j2<order+1; j2++ )
      {
        int j=i1+j2*(order+1);
        int localToGlobal=globalNodesList( elementNumber, j );
        double X=globalNodesCoords( localToGlobal, 0 );
        double Y=globalNodesCoords( localToGlobal, 1 );
        jac1+=X*derivativeBasisFunction1D( j2, i2 );
        jac3+=Y*derivativeBasisFunction1D( j2, i2 );
      }
      // detJ
      double detJ=abs( jac0*jac3-jac2*jac1 );
      double invJac0=jac3;
      double invJac1=-jac1;
      double invJac2=-jac2;
      double invJac3=jac0;
      double transpInvJac0=jac3;
      double transpInvJac1=-jac2;
      double transpInvJac2=-jac1;
      double transpInvJac3=jac0;
      double detJM1=1./detJ;
      // B
      B[i][0]=(invJac0*transpInvJac0+invJac1*transpInvJac2)*detJM1;
      B[i][1]=(invJac0*transpInvJac1+invJac1*transpInvJac3)*detJM1;
      B[i][2]=(invJac2*transpInvJac0+invJac3*transpInvJac2)*detJM1;
      B[i][3]=(invJac2*transpInvJac1+invJac3*transpInvJac3)*detJM1;
      //M
      massMatrixLocal[i]=weights[i1]*weights[i2]*detJ;
    }
  }
  #else
  {
    for( int i3=0; i3<order+1; i3++ )
    {
      for( int i2=0; i2<order+1; i2++ )
      {
        for( int i1=0; i1<order+1; i1++ )
        {
          int i=i1+i2*(order+1)+i3*(order+1)*(order+1);
          // compute jacobian matrix
          double jac00=0;
          double jac01=0;
          double jac02=0;
          double jac10=0;
          double jac11=0;
          double jac12=0;
          double jac20=0;
          double jac21=0;
          double jac22=0;

          for( int j1=0; j1<order+1; j1++ )
          {
            int j=j1+i2*(order+1)+i3*(order+1)*(order+1);
            int localToGlobal=globalNodesList( elementNumber, j );
            double X=globalNodesCoords( localToGlobal, 0 );
            double Y=globalNodesCoords( localToGlobal, 2 );
            double Z=globalNodesCoords( localToGlobal, 1 );
            jac00+=X*derivativeBasisFunction1D( j1, i1 );
            jac10+=Z*derivativeBasisFunction1D( j1, i1 );
            jac20+=Y*derivativeBasisFunction1D( j1, i1 );
          }
          for( int j2=0; j2<order+1; j2++ )
          {
            int j=i1+j2*(order+1)+i3*(order+1)*(order+1);
            int localToGlobal=globalNodesList( elementNumber, j );
            double X=globalNodesCoords( localToGlobal, 0 );
            double Y=globalNodesCoords( localToGlobal, 2 );
            double Z=globalNodesCoords( localToGlobal, 1 );
            jac01+=X*derivativeBasisFunction1D( j2, i2 );
            jac11+=Z*derivativeBasisFunction1D( j2, i2 );
            jac21+=Y*derivativeBasisFunction1D( j2, i2 );
          }
          for( int j3=0; j3<order+1; j3++ )
          {
            int j=i1+i2*(order+1)+j3*(order+1)*(order+1);
            int localToGlobal=globalNodesList( elementNumber, j );
            double X=globalNodesCoords( localToGlobal, 0 );
            double Y=globalNodesCoords( localToGlobal, 2 );
            double Z=globalNodesCoords( localToGlobal, 1 );
            jac02+=X*derivativeBasisFunction1D( j3, i3 );
            jac12+=Z*derivativeBasisFunction1D( j3, i3 );
            jac22+=Y*derivativeBasisFunction1D( j3, i3 );
          }
          // detJ
          double detJ=abs( jac00*(jac11*jac22-jac21*jac12)
                           -jac01*(jac10*jac22-jac20*jac12)
                           +jac02*(jac10*jac21-jac20*jac11));

          // inv of jac is equal of the minors of the transposed of jac
          double invJac00=jac11*jac22-jac12*jac21;
          double invJac01=jac02*jac21-jac01*jac22;
          double invJac02=jac01*jac12-jac02*jac11;
          double invJac10=jac12*jac20-jac10*jac22;
          double invJac11=jac00*jac22-jac02*jac20;
          double invJac12=jac02*jac10-jac00*jac12;
          double invJac20=jac10*jac21-jac11*jac20;
          double invJac21=jac01*jac20-jac00*jac21;
          double invJac22=jac00*jac11-jac01*jac10;

          double transpInvJac00=invJac00;
          double transpInvJac01=invJac10;
          double transpInvJac02=invJac20;
          double transpInvJac10=invJac01;
          double transpInvJac11=invJac11;
          double transpInvJac12=invJac21;
          double transpInvJac20=invJac02;
          double transpInvJac21=invJac12;
          double transpInvJac22=invJac22;

          double detJM1=1./detJ;

          // B
          B[i][0]=(invJac00*transpInvJac00+invJac01*transpInvJac10+invJac02*transpInvJac20)*detJM1;    //B11
          B[i][1]=(invJac10*transpInvJac01+invJac11*transpInvJac11+invJac12*transpInvJac21)*detJM1;    //B22
          B[i][2]=(invJac20*transpInvJac02+invJac21*transpInvJac12+invJac22*transpInvJac22)*detJM1;    //B33
          B[i][3]=(invJac00*transpInvJac01+invJac01*transpInvJac11+invJac02*transpInvJac21)*detJM1;    //B12,B21
          B[i][4]=(invJac00*transpInvJac02+invJac01*transpInvJac12+invJac02*transpInvJac22)*detJM1;    //B13,B31
          B[i][5]=(invJac10*transpInvJac02+invJac11*transpInvJac12+invJac12*transpInvJac22)*detJM1;    //B23,B32

          //M
          massMatrixLocal[i]=weights[i1]*weights[i2]*weights[i3]*detJ;
        }
      }
    }
  }
  #endif

  // compute stifness  matrix ( durufle's optimization)

  #ifdef SEM2D
  for( int i2=0; i2<order+1; i2++ )
  {
    for( int i1=0; i1<order+1; i1++ )
    {
      for( int j=0; j<nPointsPerElement; j++ )
      {
        R[j]=0;
      }
      for( int j1=0; j1<order+1; j1++ )
      {
        int j=j1+i2*(order+1);
        for( int m=0; m<order+1; m++ )
        {
          R[j]+=weights[m]*weights[i2]*(B[m+i2*(order+1)][0]*derivativeBasisFunction1D( i1, m )*derivativeBasisFunction1D( j1, m ));
        }
      }
      // B21
      for( int j1=0; j1<order+1; j1++ )
      {
        for( int j2=0; j2<order+1; j2++ )
        {
          int j=j1+j2*(order+1);
          R[j]+=weights[i1]*weights[j2]*(B[i1+j2*(order+1)][1]*derivativeBasisFunction1D( i2, j2 )*derivativeBasisFunction1D( j1, i1 ));
        }
      }
      // B12
      for( int j1=0; j1<order+1; j1++ )
      {
        for( int j2=0; j2<order+1; j2++ )
        {
          int j=j1+j2*(order+1);
          R[j]+=weights[i2]*weights[j1]*(B[i2+j1*(order+1)][2]*derivativeBasisFunction1D( i1, j1 )*derivativeBasisFunction1D( j2, i2 ));
        }
      }
      // B22
      for( int j2=0; j2<order+1; j2++ )
      {
        int j=i1+j2*(order+1);
        for( int n=0; n<order+1; n++ )
        {
          R[j]+=weights[i1]*weights[n]*(B[i1+n*(order+1)][3]*derivativeBasisFunction1D( i2, n )*derivativeBasisFunction1D( j2, n ));
        }
      }
      int i=i1+i2*(order+1);
      Y[i]=0;
      for( int j=0; j<nPointsPerElement; j++ )
      {
        Y[i]+=R[j]*pnLocal[j];
      }
    }
  }
  #else
  {
    int orderPow2=(order+1)*(order+1);
    for( int i3=0; i3<order+1; i3++ )
    {
      for( int i2=0; i2<order+1; i2++ )
      {
        for( int i1=0; i1<order+1; i1++ )
        {
          for( int j=0; j<nPointsPerElement; j++ )
          {
            R[j]=0;
          }

          //B11
          for( int j1=0; j1<order+1; j1++ )
          {
            int j=j1+i2*(order+1)+i3*orderPow2;
            for( int l=0; l<order+1; l++ )
            {
              int ll=l+i2*(order+1)+i3*orderPow2;
              R[j]+=weights[l]*weights[i2]*weights[i3]*(B[ll][0]*derivativeBasisFunction1D( i1, l )*derivativeBasisFunction1D( j1, l ));
            }
          }
          //B22
          for( int j2=0; j2<order+1; j2++ )
          {
            int j=i1+j2*(order+1)+i3*orderPow2;
            for( int m=0; m<order+1; m++ )
            {
              int mm=i1+m*(order+1)+i3*orderPow2;
              R[j]+=weights[i1]*weights[m]*weights[i3]*(B[mm][1]*derivativeBasisFunction1D( i2, m )*derivativeBasisFunction1D( j2, m ));
            }
          }
          //B33
          for( int j3=0; j3<order+1; j3++ )
          {
            int j=i1+i2*(order+1)+j3*orderPow2;
            for( int n=0; n<order+1; n++ )
            {
              int nn=i1+i2*(order+1)+n*orderPow2;
              R[j]+=weights[i1]*weights[i2]*weights[n]*(B[nn][2]*derivativeBasisFunction1D( i3, n )*derivativeBasisFunction1D( j3, n ));
            }
          }
          // B12,B21 (B[][3])
          for( int j2=0; j2<order+1; j2++ )
          {
            for( int j1=0; j1<order+1; j1++ )
            {
              int j=j1+j2*(order+1)+i3*orderPow2;
              int k=j1+i2*(order+1)+i3*orderPow2;
              int l=i1+j2*(order+1)+i3*orderPow2;
              R[j]+=weights[j1]*weights[i2]*weights[i3]*(B[k][3]*derivativeBasisFunction1D( i1, j1 )*derivativeBasisFunction1D( j2, i2 ))+
                     weights[i1]*weights[j2]*weights[i3]*(B[l][3]*derivativeBasisFunction1D( j1, i1 )*derivativeBasisFunction1D( i2, j2 ));
            }
          }
          // B13,B31 (B[][4])
          for( int j3=0; j3<order+1; j3++ )
          {
            for( int j1=0; j1<order+1; j1++ )
            {
              int j=j1+i2*(order+1)+i3*orderPow2;
              int k=j1+i2*(order+1)+i3*orderPow2;
              int l=j1+i2*(order+1)+j3*orderPow2;
              R[j]+=weights[j1]*weights[i2]*weights[i3]*(B[k][4]*derivativeBasisFunction1D( j1, i1 )*derivativeBasisFunction1D( j3, i3 ))+
                     weights[j1]*weights[i2]*weights[j3]*(B[l][4]*derivativeBasisFunction1D( j1, i1 )*derivativeBasisFunction1D( i3, j3 ));
            }
          }
          // B23,B32 (B[][5])
          for( int j3=0; j3<order+1; j3++ )
          {
            for( int j2=0; j2<order+1; j2++ )
            {
              int j=i1+j2*(order+1)+j3*orderPow2;
              int k=i1+j2*(order+1)+i3*orderPow2;
              int l=i1+i2*(order+1)+j3*orderPow2;
              R[j]+=weights[i1]*weights[j2]*weights[i3]*(B[k][5]*derivativeBasisFunction1D( i2, i2 )*derivativeBasisFunction1D( j3, i3 ))+
                     weights[i1]*weights[i2]*weights[j3]*(B[l][5]*derivativeBasisFunction1D( j2, i2 )*derivativeBasisFunction1D( i3, j3 ));
            }
          }

          int i=i1+i2*(order+1)+i3*orderPow2;
          Y[i]=0;
          for( int j=0; j<nPointsPerElement; j++ )
          {
            Y[i]+=R[j]*pnLocal[j];
          }

        }
      }
    }
  }
  #endif

  //compute gloval mass Matrix and global stiffness vector
  for( int i=0; i<nPointsPerElement; i++ )
  {
    int gIndex=globalNodesList( elementNumber, i );
    massMatrixLocal[i]/=(model[elementNumber]*model[elementNumber]);
    ATOMICADD( massMatrixGlobal[gIndex], massMatrixLocal[i] );
    ATOMICADD( yGlobal[gIndex], Y[i] );
  }

  LOOPEND

  // update pressure
  LOOPHEAD( myMeshinfo.numberOfInteriorNodes, i )
  int I=listOfInteriorNodes[i];
  pnGlobal( I, i1 )=2*pnGlobal( I, i2 )-pnGlobal( I, i1 )-myMeshinfo.myTimeStep*myMeshinfo.myTimeStep*yGlobal[I]/massMatrixGlobal[I];
  LOOPEND

  #ifdef SEM2D
  {
    LOOPHEAD( myMeshinfo.numberOfBoundaryNodes, i )
    ShGlobal[i]=0;
    LOOPEND

    LOOPHEAD( myMeshinfo.numberOfBoundaryFaces, iFace )
    //get ds
    float ds[6];
    float Sh[6];
    float Js[2][6];

    // compute ds
    int face=faceInfos( iFace, 1 );
    for( int j=0; j<order+1; j++ )
    {
      Js[0][j]=0;  // x
      Js[1][j]=0;  // y
      for( int i=0; i<order+1; i++ )
      {
        float xi=globalNodesCoords( faceInfos( iFace, 2+i ), 0 );
        float yi=globalNodesCoords( faceInfos( iFace, 2+i ), 1 );
        if( face==0 || face==2 )
        {
          Js[0][j]+=derivativeBasisFunction1D( i, j )*xi;
          Js[1][j]+=derivativeBasisFunction1D( i, j )*yi;
        }
        if( face==1 || face==3 )
        {
          Js[0][j]+=derivativeBasisFunction1D( i, j )*xi;
          Js[1][j]+=derivativeBasisFunction1D( i, j )*yi;
        }
      }
      ds[j]=sqrt( Js[0][j]*Js[0][j]+Js[1][j]*Js[1][j] );
    }

    //compute Sh and ShGlobal
    for( int i=0; i<order+1; i++ )
    {
      int gIndexFaceNode=localFaceNodeToGlobalFaceNode( iFace, i );
      Sh[i]=weights[i]*ds[i]/(model[faceInfos( iFace, 0 )]);
      ATOMICADD( ShGlobal[gIndexFaceNode], Sh[i] );
    }
    LOOPEND

    LOOPHEAD( myMeshinfo.numberOfBoundaryNodes, i )
    int I=listOfBoundaryNodes[i];
    float invMpSh=1/(massMatrixGlobal[I]+myMeshinfo.myTimeStep*ShGlobal[i]*0.5);
    float MmSh=massMatrixGlobal[I]-myMeshinfo.myTimeStep*ShGlobal[i]*0.5;
    pnGlobal( I, i1 )=invMpSh*(2*massMatrixGlobal[I]*pnGlobal( I, i2 )-MmSh*pnGlobal( I, i1 )-myMeshinfo.myTimeStep*myMeshinfo.myTimeStep*yGlobal[I]);
    LOOPEND
  }
  #endif
  FENCE
}

#endif

void SEMsolver::outputPnValues( SEMmesh mesh,
                                const int & indexTimeStep,
                                int & i1,
                                int & myElementSource,
                                const arrayReal & pnGlobal )
{
  //writes debugging ascii file.
  if( indexTimeStep%100==0 )
  {
    cout<<"TimeStep="<<indexTimeStep<<";  pnGlobal @ elementSource location "<<myElementSource
        <<" after computeOneStep = "<< pnGlobal( globalNodesList( myElementSource, 0 ), i1 )<<endl;
    #ifdef SEM_SAVE_SNAPSHOTS
    mesh.saveSnapShot( indexTimeStep, i1, pnGlobal );
    #endif
  }
}


void SEMsolver::initFEarrays( SEMmeshinfo & myMeshinfo, SEMmesh mesh )
{
  //interior elements
  mesh.globalNodesList( myMeshinfo.numberOfElements, globalNodesList );
  mesh.getListOfInteriorNodes( myMeshinfo.numberOfInteriorNodes, listOfInteriorNodes );
  mesh.nodesCoordinates( myMeshinfo.numberOfNodes, globalNodesCoords );
  // boundary elements
  mesh.getListOfBoundaryNodes( myMeshinfo.numberOfBoundaryNodes, listOfBoundaryNodes );
  mesh.getBoundaryFacesInfos( faceInfos );
  mesh.getLocalFaceNodeToGlobalFaceNode( localFaceNodeToGlobalFaceNode );
  // get model
  mesh.getModel( myMeshinfo.numberOfElements, model );
  // get quadrature points
  myQk.gaussLobattoQuadraturePoints( order, quadraturePoints );
  // get gauss-lobatto weights
  myQk.gaussLobattoQuadratureWeights( order, weights );
  // get basis function and corresponding derivatives
  myQk.getBasisFunction1D( order, quadraturePoints, basisFunction1D );
  myQk.getDerivativeBasisFunction1D( order, quadraturePoints, derivativeBasisFunction1D );

}

void SEMsolver::allocateFEarrays( SEMmeshinfo & myMeshinfo )
{
  //interior elements
  cout<<"Allocate host memory for arrays in the solver ..."<<endl;
  globalNodesList=allocateArray2D< arrayInt >( myMeshinfo.numberOfElements, myMeshinfo.numberOfPointsPerElement, "globalNodesList" );
  listOfInteriorNodes=allocateVector< vectorInt >( myMeshinfo.numberOfInteriorNodes, "listOfInteriorNodes" );
  globalNodesCoords=allocateArray2D< arrayReal >( myMeshinfo.numberOfNodes, 3, "globalNodesCoords" );
  listOfBoundaryNodes=allocateVector< vectorInt >( myMeshinfo.numberOfBoundaryNodes, "listOfBoundaryNodes" );
  faceInfos=allocateArray2D< arrayInt >( myMeshinfo.numberOfBoundaryFaces, 2+(order+1), "faceInfos" );
  localFaceNodeToGlobalFaceNode=allocateArray2D< arrayInt >( myMeshinfo.numberOfBoundaryFaces, order+1, "localFaceNodeToGlobalFaceNode" );
  model=allocateVector< vectorReal >( myMeshinfo.numberOfElements, "model" );
  quadraturePoints=allocateVector< vectorDouble >( order+1, "quadraturePoints" );
  weights=allocateVector< vectorDouble >( order+1, "weights" );
  basisFunction1D=allocateArray2D< arrayDouble >( order+1, order+1, "basisFunction1D" );
  derivativeBasisFunction1D=allocateArray2D< arrayDouble >( order+1, order+1, "derivativeBasisFunction1D" );
  //shared arrays
  massMatrixGlobal=allocateVector< vectorReal >( myMeshinfo.numberOfNodes, "massMatrixGlobal" );
  yGlobal=allocateVector< vectorReal >( myMeshinfo.numberOfNodes, "yGlobal" );
  ShGlobal=allocateVector< vectorReal >( myMeshinfo.numberOfBoundaryNodes, "ShGlobal" );
}
