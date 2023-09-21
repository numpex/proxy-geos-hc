//************************************************************************
//  SEM proxy application v.0.0.1
//
//  solverKokkos.hpp: simple 2D acoustive wave equation solver
//
//  the solverKokkos class is derived from the solverBase class
//  with the KOKKOS implementation of the solver
//
//************************************************************************

#include "solverKokkos.hpp"
 
// compute one step of the time dynamic wave equation solver
void solverKokkos::computeOneStep( const int & timeStep,
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
  
  //static vectorReal massMatrixGlobal( numberOfNodes );
  //static vectorReal yGlobal( numberOfNodes );

  Kokkos::parallel_for( numberOfNodes, KOKKOS_LAMBDA ( const int i )
  {
    massMatrixGlobal[i]=0;
    yGlobal[i]=0;
    pnGlobal(i,2)=pnGlobal(i,1); // pnm1=pn
    pnGlobal(i,1)=pnGlobal(i,0);//pn=pnp1
  } );

  // update pnGLobal with right hade side
  Kokkos::parallel_for(numberOfRHS,[=] (const int i)
  {
    int nodeRHS=globalNodesList(rhsElement[i],0);
    pnGlobal(nodeRHS,1)+=timeSample*timeSample*model[rhsElement[i]]*model[rhsElement[i]]*rhsTerm(i,timeStep);
  });

for( int irange=0;irange<numberOfElements;irange+=numberOfThreads)
{
  int rangeMin=irange;
  int rangeMax=irange+numberOfThreads;
  if(rangeMax>numberOfElements)rangeMax=numberOfElements;
  //cout<<rangeMin<<" "<<rangeMax<<endl;
  Kokkos::parallel_for(range_policy(rangeMin,rangeMax), KOKKOS_LAMBDA ( const int e )
  {
    int threadId=e-rangeMin;
    //cd bui  cout<<threadId<<" "<<irange<<" "<<e<<endl;
  
     // extract global coordinates of element e
    // get local to global indexes of nodes of element e
    int i=mesh.localToGlobalNodes(threadId, e, numberOfPointsPerElement, globalNodesList, localToGlobal );

    //get global coordinates Xi of element e
    int j=mesh.getXi(threadId, numberOfPointsPerElement, globalNodesCoords, localToGlobal, Xi );
    
    // compute jacobian Matrix
    int k=Qk.computeJacobianMatrix( threadId, numberOfPointsPerElement, Xi,
                                    derivativeBasisFunction2DX,
                                    derivativeBasisFunction2DY,
                                    jacobianMatrix );

    // compute determinant of jacobian Matrix
    int l=Qk.computeDeterminantOfJacobianMatrix( threadId, numberOfPointsPerElement,
                                                 jacobianMatrix,
                                                 detJ );
    // compute inverse of Jacobian Matrix
    int m=Qk.computeInvJacobianMatrix( threadId, numberOfPointsPerElement,
                                       jacobianMatrix,
                                       detJ,
                                       invJacobianMatrix );
                                 
    // compute transposed inverse of Jacobian Matrix
    int n=Qk.computeTranspInvJacobianMatrix( threadId,numberOfPointsPerElement,
                                             jacobianMatrix,
                                             detJ,
                                             transpInvJacobianMatrix );
                        
    // compute  geometrical transformation matrix
    int p=Qk.computeB(threadId, numberOfPointsPerElement, invJacobianMatrix, transpInvJacobianMatrix, detJ,B );

    // compute stifness and mass matrix ( durufle's optimization)
    int q=Qk.gradPhiGradPhi(threadId, numberOfPointsPerElement, order, weights2D, B, derivativeBasisFunction1D, R );

    // compute local mass matrix ( used optimez version)
    int r=Qk.phiIphiJ(threadId, numberOfPointsPerElement, weights2D, detJ, massMatrixLocal );

    // get pnGlobal to pnLocal
    for( int i=0; i<numberOfPointsPerElement; i++ )
    {
      massMatrixLocal(threadId,i)/=(model[e]*model[e]);
      //pnLocal(threadId,i)=pnGlobal(localToGlobal(threadId,i),i2);
      pnLocal(threadId,i)=pnGlobal(localToGlobal(threadId,i),1);

      //cout<<"element "<<e<<" massM "<<i<<" "<<massMatrixLocal(threadId,i);
      //cout<<"locToGLob "<<i<<" "<<localToGlobal(threadId,i)<<" "<<i2;
      //cout<<" pnLocal "<<i<<" "<<pnLocal(threadId,i)<<endl;
    }

    // compute Y=R*pnLocal
    for( int i=0; i<numberOfPointsPerElement; i++ )
    {
      Y(threadId,i)=0;
      for( int j=0; j<numberOfPointsPerElement; j++ )
      {
        Y(threadId,i)+=R(threadId,i,j)*pnLocal(threadId,j);
      }
      //cout<<"element "<<e<<" Y "<<i<<" "<<Y(threadId,i);
    }
    //cout<<endl;

    //compute gloval mass Matrix and global stiffness vector
    for( int i=0; i<numberOfPointsPerElement; i++ )
    {
      int gIndex=localToGlobal(threadId,i);
      //massMatrixGlobal[gIndex]+=massMatrixLocal(threadId,i)
      //yGlobal[gIndex]+=Y(threadId,i);
      Kokkos::atomic_add(&massMatrixGlobal[gIndex],massMatrixLocal(threadId,i));
      Kokkos::atomic_add(&yGlobal[gIndex],Y(threadId,i));
      //cout<<" element "<<e<<" massG "<<gIndex<<" "<<massMatrixGlobal[gIndex];
      //cout<<" yGLobal "<<gIndex<<" "<<yGlobal[gIndex]<<endl;
    } 
  
  } );
}

  // update pressure
  Kokkos::parallel_for( range_policy(0,numberOfInteriorNodes), KOKKOS_LAMBDA ( const int i )
  {
    int I=listOfInteriorNodes[i];
    float tmp=timeSample*timeSample;
    //pnGlobal(I,i1)=2*pnGlobal(I,i2)-pnGlobal(I,i1)-tmp*yGlobal[I]/massMatrixGlobal[I];
    pnGlobal(I,0)=2*pnGlobal(I,1)-pnGlobal(I,2)-tmp*yGlobal[I]/massMatrixGlobal[I];
  } );
  //cout<<"pressure="<<pnGlobal(5,i1)<<endl;

  // damping terms
  //static vectorReal ShGlobal( numberOfBoundaryNodes );

  for( int i=0; i<numberOfBoundaryNodes; i++ )
  {
    ShGlobal[i]=0;
  }
for( int irange=0;irange<numberOfBoundaryFaces;irange+=numberOfThreads)
{
  int rangeMin=irange;
  int rangeMax=irange+numberOfThreads;
  if(rangeMax>numberOfBoundaryFaces)rangeMax=numberOfBoundaryFaces;
  Kokkos::parallel_for(range_policy(rangeMin,rangeMax), KOKKOS_LAMBDA ( const int iFace )
  {
    int threadId=iFace-rangeMin;
  
    //get ds
    int i=Qk.computeDs( threadId, iFace, order, faceInfos,numOfBasisFunctionOnFace,
                        Js, globalNodesCoords, derivativeBasisFunction2DX,
                        derivativeBasisFunction2DY,
                        ds );
  
    //compute Sh and ShGlobal
    for( int i=0; i<order+1; i++ )
    {
      int gIndexFaceNode=localFaceNodeToGlobalFaceNode(iFace,i);
      Sh(threadId,i)=weights[i]*ds(threadId,i)/(model[faceInfos(iFace,0)]);
      //ShGlobal[gIndexFaceNode]+=Sh(threadId,i);
      Kokkos::atomic_add(&ShGlobal[gIndexFaceNode],Sh(threadId,i));
    }
  } );
}
  // update pressure @ boundaries;
  float tmp=timeSample*timeSample;
  Kokkos::parallel_for( range_policy(0,numberOfBoundaryNodes), KOKKOS_LAMBDA  ( const int i )
  {
    int I=listOfBoundaryNodes[i];
    float invMpSh=1/(massMatrixGlobal[I]+timeSample*ShGlobal[i]*0.5);
    float MmSh=massMatrixGlobal[I]-timeSample*ShGlobal[i]*0.5;
    //pnGlobal(I,i1)=invMpSh*(2*massMatrixGlobal[I]*pnGlobal(I,i2)-MmSh*pnGlobal(I,i1)-tmp*yGlobal[I]);
    pnGlobal(I,0)=invMpSh*(2*massMatrixGlobal[I]*pnGlobal(I,1)-MmSh*pnGlobal(I,2)-tmp*yGlobal[I]);
  } );

  

  /**
     for ( int i=0 ; i< numberOfBoundaryNodes; i++)
     {
     cout<<"ShGlobal["<<i<<"]="<<ShGlobal[i]<<endl;
     }
     for ( int i=0 ; i< numberOfBoundaryNodes; i++)
     {
     int I=listOfBoundaryNodes[i];
     cout<<"i="<<i<<", BoundaryNode="<<I<<endl;
     }
   **/
}
