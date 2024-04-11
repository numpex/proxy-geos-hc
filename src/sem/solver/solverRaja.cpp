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
                                  int & it1,
                                  int & it2,
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
  vectorDoubleView d_weights3D=weights3D.toView();

  // get basis function and corresponding derivatives
  arrayDoubleView d_basisFunction1D=basisFunction1D.toView();
  arrayDoubleView d_derivativeBasisFunction1D=derivativeBasisFunction1D.toView();

  //shared arrays
  vectorDoubleView d_massMatrixGlobal=massMatrixGlobal.toView();
  vectorDoubleView d_yGlobal=yGlobal.toView();

  RAJA::forall< deviceExecPolicy >( RAJA::RangeSegment( 0, numberOfNodes ),  [=] LVARRAY_HOST_DEVICE  ( int i ) {
    d_massMatrixGlobal[i]=0;
    d_yGlobal[i]=0;
  } );

  // update pnGLobal with right hade side
  RAJA::forall< deviceExecPolicy >( RAJA::RangeSegment( 0, numberOfRHS ), [=] LVARRAY_HOST_DEVICE  ( int i ) 
  {
    int nodeRHS=d_globalNodesList(d_rhsElement[i],0);
    d_pnGlobal(nodeRHS,it2)+=timeSample*timeSample*d_model[d_rhsElement[i]]*d_model[d_rhsElement[i]]*d_rhsTerm(i,timeStep);
  });
 
  int numberOfPointsPerElement=(order+1)*(order+1)*(order+1);
  RAJA::forall< deviceExecPolicy >( RAJA::RangeSegment( 0, numberOfElements ), [=] LVARRAY_HOST_DEVICE ( int e )
  {
        // start parallel section
	double B[NumPointsPerElem][6];
        double R[NumPointsPerElem];
        double massMatrixLocal[NumPointsPerElem];
	double pnLocal[NumPointsPerElem];
	double Y[NumPointsPerElem];

	int    orderPow2=(order+1)*(order+1);
        // get pnGlobal to pnLocal
        for( int i=0; i<numberOfPointsPerElement; i++ )
        {
           int localToGlobal=d_globalNodesList(e,i);
           pnLocal[i]=d_pnGlobal(localToGlobal,it2);
        }

	//compute Jac, B  and Mass Matrix
	for (int i3=0;i3<order+1;i3++)
	{
	   for (int i2=0;i2<order+1;i2++)
	   {
	       for (int i1=0;i1<order+1;i1++)
	       {
	          int i=i1+i2*(order+1)+i3*orderPow2;
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

	          for (int j1=0;j1<order+1;j1++)
	          {
	              int j=j1+i2*(order+1)+i3*orderPow2;
	              int localToGlobal=d_globalNodesList(e,j);
                      double X=d_globalNodesCoords(localToGlobal,0);
                      double Y=d_globalNodesCoords(localToGlobal,2);
                      double Z=d_globalNodesCoords(localToGlobal,1);
	              jac00+=X*d_derivativeBasisFunction1D(j1,i1);
	              jac10+=Z*d_derivativeBasisFunction1D(j1,i1);
	              jac20+=Y*d_derivativeBasisFunction1D(j1,i1);
		  }
	          for (int j2=0;j2<order+1;j2++)
	          {
	              int j=i1+j2*(order+1)+i3*orderPow2;
	              int localToGlobal=d_globalNodesList(e,j);
                      double X=d_globalNodesCoords(localToGlobal,0);
                      double Y=d_globalNodesCoords(localToGlobal,2);
                      double Z=d_globalNodesCoords(localToGlobal,1);
                      jac01+=X*d_derivativeBasisFunction1D(j2,i2);
                      jac11+=Z*d_derivativeBasisFunction1D(j2,i2);
                      jac21+=Y*d_derivativeBasisFunction1D(j2,i2);
	          }
	          for (int j3=0;j3<order+1;j3++)
	          {
	              int j=i1+i2*(order+1)+j3*orderPow2;
	              int localToGlobal=d_globalNodesList(e,j);
                      double X=d_globalNodesCoords(localToGlobal,0);
                      double Y=d_globalNodesCoords(localToGlobal,2);
                      double Z=d_globalNodesCoords(localToGlobal,1);
                      jac02+=X*d_derivativeBasisFunction1D(j3,i3);
                      jac12+=Z*d_derivativeBasisFunction1D(j3,i3);
                      jac22+=Y*d_derivativeBasisFunction1D(j3,i3);
                  }
	          // detJ
                  double detJ=abs(jac00*(jac11*jac22-jac21*jac12)
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
                  B[i][0]=(invJac00*transpInvJac00+invJac01*transpInvJac10+invJac02*transpInvJac20)*detJM1;//B11
                  B[i][1]=(invJac10*transpInvJac01+invJac11*transpInvJac11+invJac12*transpInvJac21)*detJM1;//B22
                  B[i][2]=(invJac20*transpInvJac02+invJac21*transpInvJac12+invJac22*transpInvJac22)*detJM1;//B33
                  B[i][3]=(invJac00*transpInvJac01+invJac01*transpInvJac11+invJac02*transpInvJac21)*detJM1;//B12,B21
                  B[i][4]=(invJac00*transpInvJac02+invJac01*transpInvJac12+invJac02*transpInvJac22)*detJM1;//B13,B31
                  B[i][5]=(invJac10*transpInvJac02+invJac11*transpInvJac12+invJac12*transpInvJac22)*detJM1;//B23,B32

	          //M
                  massMatrixLocal[i]=d_weights3D[i]*detJ;
               }
	   }
        }

	// compute Stiffness matrixAND STIFFNESS VECTOR
	for (int i3=0;i3<order+1;i3++)
	{
	   for (int i2=0;i2<order+1;i2++)
	   {
	       for (int i1=0;i1<order+1;i1++)
	       {
	          int i=i1+i2*(order+1)+i3*orderPow2;

                  for( int j3=0; j3<order+1; j3++ )
                  {
                     for( int j2=0; j2<order+1; j2++ )
                     {
                        for( int j1=0; j1<order+1; j1++ )
                        {
	                   int j=j1+j2*(order+1)+j3*orderPow2;
			   R[j]=0;
			}
	             }
		  }

	          //B11
                  for( int j1=0; j1<order+1; j1++ )
                  {
                      int j=j1+i2*(order+1)+i3*orderPow2;
                      for( int l=0; l<order+1; l++ )
                      {
			  int ll=l+i2*(order+1)+i3*orderPow2;
                          R[j]+=d_weights3D[ll]*(B[ll][0]*d_derivativeBasisFunction1D(i1,l)*d_derivativeBasisFunction1D(j1,l));
                      }
                  }
		  //B22
                  for( int j2=0; j2<order+1; j2++ )
                  {
                      int j=i1+j2*(order+1)+i3*orderPow2;
                      for( int m=0; m<order+1; m++ )
                      {
			  int mm=i1+m*(order+1)+i3*orderPow2;
                          R[j]+=d_weights3D[mm]*(B[mm][1]*d_derivativeBasisFunction1D(i2,m)*d_derivativeBasisFunction1D(j2,m));
                      }
                  }
                  //B33 
                  for( int j3=0; j3<order+1; j3++ )
                  {
                      int j=i1+i2*(order+1)+j3*orderPow2;
                      for( int n=0; n<order+1; n++ )
                      {
			  int nn=i1+i2*(order+1)+n*orderPow2;
                          R[j]+=d_weights3D[nn]*(B[nn][2]*d_derivativeBasisFunction1D(i3,n)*d_derivativeBasisFunction1D(j3,n));
                      }
                  }
                  // B12,B21 (B[][3])
                  for( int j1=0; j1<order+1; j1++ )
                  {
                    for( int j2=0; j2<order+1; j2++ )
                    {
                      int j=j1+j2*(order+1)+i3*orderPow2;
		      int k=j1+i2*(order+1)+i3*orderPow2;
		      int l=i1+j2*(order+1)+i3*orderPow2;
                      R[j]+=d_weights3D[k]*(B[k][3]*d_derivativeBasisFunction1D(i1,j1)*d_derivativeBasisFunction1D(j2,i2))+
                            d_weights3D[l]*(B[l][3]*d_derivativeBasisFunction1D(j1,i1)*d_derivativeBasisFunction1D(i2,j2));
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
                      R[j]+=d_weights3D[k]*(B[k][4]*d_derivativeBasisFunction1D(j1,i1)*d_derivativeBasisFunction1D(j3,i3))+
                            d_weights3D[l]*(B[l][4]*d_derivativeBasisFunction1D(j1,i1)*d_derivativeBasisFunction1D(i3,j3));
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
                      R[j]+=d_weights3D[k]*(B[k][5]*d_derivativeBasisFunction1D(i2,i2)*d_derivativeBasisFunction1D(j3,i3))+
                            d_weights3D[l]*(B[l][5]*d_derivativeBasisFunction1D(j2,i2)*d_derivativeBasisFunction1D(i3,j3));
                    }
                  }

                  Y[i]=0;
                  for( int j=0; j<numberOfPointsPerElement; j++ )
                  {
                    Y[i]+=R[j]*pnLocal[j];
                  }

               }
	   }
        }


      //compute global mass Matrix and global stiffness vector
      for( int i=0; i<numberOfPointsPerElement; i++ )
      {
        int gIndex=d_globalNodesList(e,i);
        massMatrixLocal[i]/=(d_model[e]*d_model[e]);
        RAJA::atomicAdd< deviceAtomicPolicy >(&d_massMatrixGlobal[gIndex],massMatrixLocal[i]);
        RAJA::atomicAdd< deviceAtomicPolicy>(&d_yGlobal[gIndex],Y[i]);
      }

  });
  // update pressure
  RAJA::forall< deviceExecPolicy>( RAJA::RangeSegment( 0, numberOfInteriorNodes ), [=] LVARRAY_HOST_DEVICE ( int i ) {
    int I=d_listOfInteriorNodes[i];
    float tmp=timeSample*timeSample;
    d_pnGlobal[I][it1]=2*d_pnGlobal[I][it2]-d_pnGlobal[I][it1]-tmp*d_yGlobal[I]/d_massMatrixGlobal[I];
  });

  
  if(timeStep%100==0)
  {
     int nodeRHS=d_globalNodesList(d_rhsElement[0],0);
     RAJA::forall< RAJA::seq_exec>( RAJA::RangeSegment( 0, numberOfNodes ), [=] LVARRAY_HOST_DEVICE ( int i )
     {
          pnGlobal(nodeRHS,it1)=d_pnGlobal(nodeRHS,it1);
     });
  }
}
