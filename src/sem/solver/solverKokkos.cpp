//************************************************************************
//   proxy application v.0.0.1
//
//  solverRaja.cpp: simple 2D acoustive wave equation solver
//
//  the solverRaja class is derived from the solverBase class
//  with the RAJA implementation of the solver
//
//************************************************************************

#include "solverKokkos.hpp"

// compute one step of the time dynamic wave equation solver
void solverKokkos::computeOneStep(  const int & timeStep,
                                   const float & timeSample,
                                   const int & order,
                                   int & i1,
                                   int & i2,
                                   const int & numberOfRHS,
                                   vectorInt & rhsElement,
                                   arrayReal & rhsTerm,
                                   arrayReal & pnGlobal)
{

  auto d_rhsElement = Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, rhsElement);
  auto d_rhsTerm = Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, rhsTerm);
  auto d_pnGlobal = Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, pnGlobal);


  auto d_globalNodesList = Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, globalNodesList);
  auto d_globalNodesCoords = Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, globalNodesCoords);
  auto d_listOfInteriorNodes = Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, listOfInteriorNodes);

  auto d_model = Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, model);
  auto d_weights2D = Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, weights2D);


  auto d_derivativeBasisFunction1D = Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, derivativeBasisFunction1D);
  auto d_derivativeBasisFunction2DX = Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, derivativeBasisFunction2DX);
  auto d_derivativeBasisFunction2DY = Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, derivativeBasisFunction2DY);

  auto d_massMatrixGlobal = Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, massMatrixGlobal);
  auto d_yGlobal = Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, yGlobal);
  auto d_ShGlobal = Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, ShGlobal);

  Kokkos::parallel_for( numberOfNodes, KOKKOS_CLASS_LAMBDA ( const int i )
  {
    d_massMatrixGlobal[i]=0;
    d_yGlobal[i]=0;
  } );

  // update pnGLobal with right hade side
  Kokkos::parallel_for(numberOfRHS,KOKKOS_CLASS_LAMBDA (const int i)
  {
    int nodeRHS=d_globalNodesList(d_rhsElement[i],0);
    d_pnGlobal(nodeRHS,i2)+=timeSample*timeSample*d_model[d_rhsElement[i]]*d_model[d_rhsElement[i]]*d_rhsTerm(i,timeStep);
  });
  Kokkos::fence();
  
  // Set up a policy that launches NumberOfNodes teams, with the maximum number
  // of threads per team.
  using team_policy=Kokkos::TeamPolicy<>;
  using Kokkos::TeamThreadRange;
  int nthreads=32;
  int numberOfPointsPerElement=(order+1)*(order+1);
  int leagueSize=(numberOfElements+nthreads-1)/nthreads;
  const team_policy policyElements(leagueSize,nthreads);// Kokkos::AUTO);
  Kokkos::parallel_for( policyElements, KOKKOS_CLASS_LAMBDA ( const Kokkos::TeamPolicy<>::member_type & team )
  {
     int teamSize=team.league_size();
     int teamId=team.league_rank();
     int threadId=team.team_rank();
     int blockThreadSize=team.team_size();
     //if(teamId<=2)printf("Hello World: %d %d // %d %d \n", teamId, threadId,teamSize , blockThreadSize);
     //printf("Hello World: %d %d // %d %d \n", teamId, threadId,teamSize , blockThreadSize);

     Kokkos::parallel_for(TeamThreadRange(team,blockThreadSize), [=] (const int index)
     {
        int e=teamId*blockThreadSize+index;
      if(e<numberOfElements-1){
        // start parallel section
        double Xi[36][2];
        double   B[36][4];
        double R[36][36];
        double massMatrixLocal[36];
        double   pnLocal[36];
        double   Y[36];
      
        for(   int x=0; x<order+1;x++)
        {
           for( int y=0; y<order+1;y++)
           {
             int i=x+y*(order+1);
             int localToGlobal=d_globalNodesList(e,i);
               Xi[i][0]=d_globalNodesCoords(localToGlobal,0);
               Xi[i][1]=d_globalNodesCoords(localToGlobal,1);
           }//);
        }//);
        for   (int i=0;i<numberOfPointsPerElement;i++)
        {
           // compute jacobian matrix
             double jac0=0;
             double jac1=0;
             double jac2=0;
             double jac3=0;

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
          
        //   compute stiffness
        for   (int i=0;i<numberOfPointsPerElement;i++)
        {
           for (int j=0;j<numberOfPointsPerElement;j++)
           {
             R[i][j]=0;
           }//);
        }//);
      
        for   (int i1=0;i1<order+1;i1++)
        {
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
        for   (int i1=0;i1<order+1;i1++)
        {
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
        for   (int i1=0;i1<order+1;i1++)
        {
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
        for   (int i1=0;i1<order+1;i1++)
        {
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
          int   localToGlobal=d_globalNodesList(e,i);
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
          Kokkos::atomic_add(&d_massMatrixGlobal[gIndex],massMatrixLocal[i]);
          Kokkos::atomic_add(&d_yGlobal[gIndex],Y[i]);
        }
      } // end if ( e <numberOfElements - 1 )
    });
    team.team_barrier();
  });
  Kokkos::fence();
  // update pressure
  Kokkos::parallel_for( range_policy(0,numberOfInteriorNodes), KOKKOS_CLASS_LAMBDA ( const int i )
  {
    int I=d_listOfInteriorNodes[i];
    float tmp=timeSample*timeSample;
    d_pnGlobal(I,i1)=2*d_pnGlobal(I,i2)-d_pnGlobal(I,i1)-tmp*d_yGlobal[I]/d_massMatrixGlobal[I];
  } );
  Kokkos::fence();

  Kokkos::parallel_for( range_policy(0,numberOfBoundaryNodes), KOKKOS_CLASS_LAMBDA ( const int i )
  {
    d_ShGlobal[i]=0;
  } );

  Kokkos::deep_copy( pnGlobal, d_pnGlobal);
  /* 
  Kokkos::parallel_for (numberOfBoundaryFaces, KOKKOS_CLASS_LAMBDA (const int iFace)
  {
    //get ds
    float ds[6];
    float Sh[6];
    int numOfBasisFunctionOnFace[6];
    float Js[2][6];

    int i=Qk.computeDs( iFace, order, faceInfos,numOfBasisFunctionOnFace,
                  Js, globalNodesCoords, derivativeBasisFunction2DX,
                  derivativeBasisFunction2DY,
                  ds );
    //
    //compute Sh and ShGlobal
    for( int i=0; i<order+1; i++ )
    {
      int gIndexFaceNode=localFaceNodeToGlobalFaceNode(iFace,i);
      Sh[i]=weights[i]*ds[i]/(model[faceInfos(iFace,0)]);
      Kokkos::atomic_add(&ShGlobal[gIndexFaceNode],Sh[i]);
    }
  } );

  Kokkos::fence();

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
  */
}
