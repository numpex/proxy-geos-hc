//************************************************************************
//   proxy application v.0.0.1
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

  RAJA::forall< deviceExecPolicy >( RAJA::RangeSegment( 0, numberOfNodes ),  [=] LVARRAY_HOST_DEVICE  ( int i ) {
    d_massMatrixGlobal[i]=0;
    d_yGlobal[i]=0;
  } );

  // update pnGLobal with right hade side
  RAJA::forall< deviceExecPolicy >( RAJA::RangeSegment( 0, numberOfRHS ), [=] LVARRAY_HOST_DEVICE  ( int i ) 
  {
    int nodeRHS=d_globalNodesList(d_rhsElement[i],0);
    d_pnGlobal(nodeRHS,i2)+=timeSample*timeSample*d_model[d_rhsElement[i]]*d_model[d_rhsElement[i]]*d_rhsTerm(i,timeStep);
  });
 
  int numberOfPointsPerElement=(order+1)*(order+1);

  const int TEAM_SZ = 256;
  const int GRID_SZ = RAJA_DIVIDE_CEILING_INT(numberOfElements,TEAM_SZ);
  RAJA::TypedRangeSegment<int> e_Range(0, numberOfElements);
  RAJA::TypedRangeSegment<int> qx_Range(0, order+1);
  RAJA::TypedRangeSegment<int> qy_Range(0, order+1);
  /*
  #ifdef ENABLE_CUDA
    using cuda_thread_e = RAJA::LoopPolicy<RAJA::cuda_thread_x_loop>;
    using cuda_thread_y = RAJA::LoopPolicy<RAJA::cuda_thread_y_loop>;
    using cuda_thread_z = RAJA::LoopPolicy<RAJA::cuda_thread_z_loop>;
    using launch_policy_cuda = RAJA::LaunchPolicy<RAJA::cuda_launch_t<async>>;
  #elif defined ENABLE_HIP
    using cuda_thread_e = RAJA::LoopPolicy<RAJA::hip_thread_x_loop>;
    using cuda_thread_y = RAJA::LoopPolicy<RAJA::hip_thread_y_loop>;
    using cuda_thread_z = RAJA::LoopPolicy<RAJA::hip_thread_z_loop>;
    using launch_policy_cuda = RAJA::LaunchPolicy<RAJA::hip_launch_t<async>>;
  #endif
  */

  const bool async = false; //execute asynchronously

  //RAJA::launch<launch_policy_cuda>
  //  (RAJA::LaunchParams(RAJA::Teams(2), RAJA::Threads(32,8)),
  //  [=] RAJA_HOST_DEVICE (RAJA::LaunchContext ctx) {

  //    RAJA::loop<cuda_thread_e>(ctx, e_Range, [&] (int e) {
  RAJA::forall< deviceExecPolicy >( RAJA::RangeSegment( 0, numberOfElements ), [=] LVARRAY_HOST_DEVICE ( int e )
  {
        // start parallel section
        double Xi[36][2];
	double B[36][4];
        double R[36][36];
        double massMatrixLocal[36];
	double pnLocal[36];
	double Y[36];
	//RAJA::loop<cuda_thread_y>(ctx, RAJA::RangeSegment(0, order+1), [&](int x) 
	for( int x=0; x<order+1;x++)
	{
           //RAJA::loop<cuda_thread_z>(ctx, RAJA::RangeSegment(0, order+1), [&](int y) 
	   for( int y=0; y<order+1;y++)
	   {
	     int i=x+y*(order+1);
	     int localToGlobal=d_globalNodesList(e,i);
             Xi[i][0]=d_globalNodesCoords(localToGlobal,0);
             Xi[i][1]=d_globalNodesCoords(localToGlobal,1);
	   }//);
        }//);
	
	//RAJA::loop<cuda_thread_y>(ctx, RAJA::RangeSegment(0, numberOfPointsPerElement), [&](int i) 
	for (int i=0;i<numberOfPointsPerElement;i++)
	{
	   // compute jacobian matrix
           double jac0=0;
           double jac1=0;
           double jac2=0;
           double jac3=0;

           //RAJA::loop<cuda_thread_z>(ctx, RAJA::RangeSegment(0, numberOfPointsPerElement), [&](int j) 
	   for (int j=0;j<numberOfPointsPerElement;j++)
	   {
	     jac0+=Xi[j][0]*d_derivativeBasisFunction2DX(j,i);
             jac1+=Xi[j][0]*d_derivativeBasisFunction2DY(j,i);
             jac2+=Xi[j][1]*d_derivativeBasisFunction2DX(j,i);
             jac3+=Xi[j][1]*d_derivativeBasisFunction2DY(j,i);
           }//);
	   // detJ
           double detJ=abs(jac0*jac3-jac2*jac1);

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
           massMatrixLocal[i]=d_weights2D[i]*detJ;

        }//);
        
	// compute stiffness
	//RAJA::loop<cuda_thread_y>(ctx, RAJA::RangeSegment(0, numberOfPointsPerElement), [&](int i) 
	for (int i=0;i<numberOfPointsPerElement;i++)
	{
           //RAJA::loop<cuda_thread_z>(ctx, RAJA::RangeSegment(0, numberOfPointsPerElement), [&](int j) 
	   for (int j=0;j<numberOfPointsPerElement;j++)
	   {
	     R[i][j]=0;
	   }//);
	}//);
	
	//RAJA::loop<cuda_thread_y>(ctx, RAJA::RangeSegment(0, order+1), [&](int i1) 
	for (int i1=0;i1<order+1;i1++)
	{
           //RAJA::loop<cuda_thread_z>(ctx, RAJA::RangeSegment(0, order+1), [&](int i2) 
	   for (int i2=0;i2<order+1;i2++)
	   {
	        int i=i1+i2*(order+1);
                for( int j1=0; j1<order+1; j1++ )
                {
                  int j=j1+i2*(order+1);
                  for( int m=0; m<order+1; m++ )
                  {
                    R[i][j]+=d_weights2D[m+i2*(order+1)]*(B[m+i2*(order+1)][0]*
		             d_derivativeBasisFunction1D(i1,m)*d_derivativeBasisFunction1D(j1,m));
                  }
                }
           }//);
        }//);
        // B21
        //RAJA::loop<cuda_thread_y>(ctx, RAJA::RangeSegment(0, order+1), [&](int i1) 
	for (int i1=0;i1<order+1;i1++)
        {
          //RAJA::loop<cuda_thread_z>(ctx, RAJA::RangeSegment(0, order+1), [&](int i2) 
	  for (int i2=0;i2<order+1;i2++)
          {
           int i=i1+i2*(order+1);
           for( int j1=0; j1<order+1; j1++ )
           {
             for( int j2=0; j2<order+1; j2++ )
             {
               int j=j1+j2*(order+1);
               R[i][j]+=d_weights2D[i1+j2*(order+1)]*(B[i1+j2*(order+1)][1]*
			d_derivativeBasisFunction1D(i2,j2)*d_derivativeBasisFunction1D(j1,i1));
             }
          }
         }//);
        }//);
        // B12
        //RAJA::loop<cuda_thread_y>(ctx, RAJA::RangeSegment(0, order+1), [&](int i1) 
	for (int i1=0;i1<order+1;i1++)
        {
          //RAJA::loop<cuda_thread_z>(ctx, RAJA::RangeSegment(0, order+1), [&](int i2) 
	  for (int i2=0;i2<order+1;i2++)
          {
           int i=i1+i2*(order+1);
           for( int j1=0; j1<order+1; j1++ )
           {
             for( int j2=0; j2<order+1; j2++ )
             {
               int j=j1+j2*(order+1);
               R[i][j]+=d_weights2D[i2+j1*(order+1)]*(B[i2+j1*(order+1)][2]*
                        d_derivativeBasisFunction1D(i1,j1)*d_derivativeBasisFunction1D(j2,i2));
             }
           }
         }//);
        }//);
        // B22
        //RAJA::loop<cuda_thread_y>(ctx, RAJA::RangeSegment(0, order+1), [&](int i1) 
	for (int i1=0;i1<order+1;i1++)
        {
          //RAJA::loop<cuda_thread_z>(ctx, RAJA::RangeSegment(0, order+1), [&](int i2) 
	  for (int i2=0;i2<order+1;i2++)
          {
           int i=i1+i2*(order+1);
           for( int j2=0; j2<order+1; j2++ )
           {
              int j=i1+j2*(order+1);
              for( int n=0; n<order+1; n++ )
              {
              R[i][j]+=d_weights2D[i1+n*(order+1)]*(B[i1+n*(order+1)][3]*
	               d_derivativeBasisFunction1D(i2,n)*d_derivativeBasisFunction1D(j2,n));
              }
           }
          }//);
        }//);

       // get pnGlobal to pnLocal
      for( int i=0; i<numberOfPointsPerElement; i++ )
      {
	int localToGlobal=d_globalNodesList(e,i);
        massMatrixLocal[i]/=(d_model[e]*d_model[e]);
        pnLocal[i]=d_pnGlobal(localToGlobal,i2);
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

      //compute global mass Matrix and global stiffness vector
      for( int i=0; i<numberOfPointsPerElement; i++ )
      {
        int gIndex=d_globalNodesList(e,i);
        //massMatrixGlobal[gIndex]+=massMatrixLocal(threadId,i)
        //yGlobal[gIndex]+=Y(threadId,i);
        RAJA::atomicAdd< deviceAtomicPolicy >(&d_massMatrixGlobal[gIndex],massMatrixLocal[i]);
        RAJA::atomicAdd< deviceAtomicPolicy>(&d_yGlobal[gIndex],Y[i]);
      }

    });
  //});
  // update pressure
  RAJA::forall< deviceExecPolicy>( RAJA::RangeSegment( 0, numberOfInteriorNodes ), [=] LVARRAY_HOST_DEVICE ( int i ) {
    int I=d_listOfInteriorNodes[i];
    float tmp=timeSample*timeSample;
    d_pnGlobal[I][i1]=2*d_pnGlobal[I][i2]-d_pnGlobal[I][i1]-tmp*d_yGlobal[I]/d_massMatrixGlobal[I];
  } );

  RAJA::forall< deviceExecPolicy>( RAJA::RangeSegment( 0, numberOfBoundaryNodes ), [=] LVARRAY_HOST_DEVICE ( int i ) {
    d_ShGlobal[i]=0;
  } );
  
  /*
  RAJA::forall< deviceExecPolicy >( RAJA::RangeSegment( 0, numberOfBoundaryFaces ), [=] LVARRAY_HOST_DEVICE ( int iFace ){
    //get ds
    float ds[6];
    float Sh[6];
    int numOfBasisFunctionOnFace[6];
    float Js[2][6];

    int i=Qk.computeDs( iFace, order, d_faceInfos,numOfBasisFunctionOnFace,
                  Js, d_globalNodesCoords, d_derivativeBasisFunction2DX,
                  d_derivativeBasisFunction2DY,
                  ds );
    //
    //compute Sh and ShGlobal
    for( int i=0; i<order+1; i++ )
    {
      int gIndexFaceNode=d_localFaceNodeToGlobalFaceNode(iFace,i);
      Sh[i]=d_weights[i]*ds[i]/(d_model[d_faceInfos(iFace,0)]);
      RAJA::atomicAdd< deviceAtomicPolicy >(&d_ShGlobal[gIndexFaceNode],Sh[i]);
    }
  } );

  // update pressure @ boundaries;
  float tmp=timeSample*timeSample;
  RAJA::forall< deviceExecPolicy >( RAJA::RangeSegment( 0, numberOfBoundaryNodes ), [=] LVARRAY_HOST_DEVICE ( int i ) {
    int I=d_listOfBoundaryNodes[i];
    float invMpSh=1/(d_massMatrixGlobal[I]+timeSample*d_ShGlobal[i]*0.5);
    float MmSh=d_massMatrixGlobal[I]-timeSample*d_ShGlobal[i]*0.5;
    d_pnGlobal[I][i1]=invMpSh*(2*d_massMatrixGlobal[I]*d_pnGlobal[I][i2]-MmSh*d_pnGlobal[I][i1]-tmp*d_yGlobal[I]);
  } );

  if(timeStep%100==0)
  {
     int nodeRHS=d_globalNodesList(d_rhsElement[0],0);
     RAJA::forall< RAJA::seq_exec>( RAJA::RangeSegment( 0, numberOfNodes ), [=] LVARRAY_HOST_DEVICE ( int i )
     {
          pnGlobal(nodeRHS,i1)=d_pnGlobal(nodeRHS,i1);
     });
  }
  */
}
