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
void solverRaja::computeOneStep( const float & timeSample,
                                 const int & order,
                                 int & i1,
                                 int & i2,
                                 arrayReal & pnGlobal,
                                 simpleMesh mesh,
                                 QkGL Qk )
{

  RAJA::forall< RAJA::omp_parallel_for_exec >( RAJA::RangeSegment( 0, numberOfNodes ), [=] ( int i ) {
    massMatrixGlobal[i]=0;
    yGlobal[i]=0;
  } );


  // loop over mesh elements
  RAJA::forall< RAJA::omp_parallel_for_exec >( RAJA::RangeSegment( 0, numberOfElements ), [=, &pnGlobal] ( int e )
  {
    // extract global coordinates of element e
    // get local to global indexes of nodes of element e
    mesh.localToGlobalNodes( e, numberOfPointsPerElement, globalNodesList, localToGlobal );

    //get global coordinates Xi of element e
    mesh.getXi( numberOfPointsPerElement, globalNodesCoords, localToGlobal, Xi );

    // compute jacobian Matrix
    Qk.computeJacobianMatrix( numberOfPointsPerElement, Xi,
                              derivativeBasisFunction2DX,
                              derivativeBasisFunction2DY,
                              jacobianMatrix );

    // compute determinant of jacobian Matrix
    Qk.computeDeterminantOfJacobianMatrix( numberOfPointsPerElement,
                                           jacobianMatrix,
                                           detJ );
    // compute inverse of Jacobian Matrix
    Qk.computeInvJacobianMatrix( numberOfPointsPerElement,
                                 jacobianMatrix,
                                 detJ,
                                 invJacobianMatrix );
                                 
    // compute transposed inverse of Jacobian Matrix
    Qk.computeTranspInvJacobianMatrix( numberOfPointsPerElement,
                                       jacobianMatrix,
                                       detJ,
                                       transpInvJacobianMatrix );
                        
    // compute  geometrical transformation matrix
    Qk.computeB( numberOfPointsPerElement, invJacobianMatrix, transpInvJacobianMatrix, detJ,B );

    // compute stifness and mass matrix ( durufle's optimization)
    Qk.gradPhiGradPhi( numberOfPointsPerElement, order, weights2D, B, derivativeBasisFunction1D, R );

    // compute local mass matrix ( used optimez version)
    Qk.phiIphiJ( numberOfPointsPerElement, weights2D, detJ, massMatrixLocal );

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
        Y[i]+=R(i,j)*pnLocal[j];
      }
    }
  } );

  // update pressure
  RAJA::forall< RAJA::omp_parallel_for_exec >( RAJA::RangeSegment( 0, numberOfInteriorNodes ), [=, &pnGlobal] ( int i ) {
    int I=listOfInteriorNodes[i];
    float tmp=timeSample*timeSample;
    pnGlobal[I][i1]=2*pnGlobal[I][i2]-pnGlobal[I][i1]-tmp*yGlobal[I]/massMatrixGlobal[I];
  } );
  //cout<<"pressure="<<pnGlobal[5][i1]<<endl;

  RAJA::forall< RAJA::omp_parallel_for_exec >( RAJA::RangeSegment( 0, numberOfBoundaryNodes ), [=] ( int i ) {
    ShGlobal[i]=0;
  } );
  // Note: this loop is data parallel.
  RAJA::forall< RAJA::omp_parallel_for_exec >( RAJA::RangeSegment( 0, numberOfBoundaryFaces ), [=] ( int iFace ){
   //get ds
    Qk.computeDs( iFace, order, faceInfos,numOfBasisFunctionOnFace,
                  Js, globalNodesCoords, derivativeBasisFunction2DX,
                  derivativeBasisFunction2DY,
                  ds );
    //compute Sh and ShGlobal
    for( int i=0; i<order+1; i++ )
    {
      int gIndexFaceNode=localFaceNodeToGlobalFaceNode(iFace,i);
      Sh[i]=weights[i]*ds[i]/(model[faceInfos(iFace,0)]);
      ShGlobal[gIndexFaceNode]+=Sh[i];
    }
  } );

  // update pressure @ boundaries;
  float tmp=timeSample*timeSample;
  RAJA::forall< RAJA::omp_parallel_for_exec >( RAJA::RangeSegment( 0, numberOfBoundaryNodes ), [=, &pnGlobal] ( int i ) {
    int I=listOfBoundaryNodes[i];
    float invMpSh=1/(massMatrixGlobal[I]+timeSample*ShGlobal[i]*0.5);
    float MmSh=massMatrixGlobal[I]-timeSample*ShGlobal[i]*0.5;
    pnGlobal[I][i1]=invMpSh*(2*massMatrixGlobal[I]*pnGlobal[I][i2]-MmSh*pnGlobal[I][i1]-tmp*yGlobal[I]);
  } );

}
