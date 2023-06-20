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
void solverKokkos::computeOneStep( const float & timeSample,
                                   const int & order,
                                   int & i1,
                                   int & i2,
                                   arrayReal & pnGlobal,
                                   simpleMesh mesh,
                                   QkGL Qk )
{
  // get infos from mesh
  static int numberOfNodes=mesh.getNumberOfNodes();
  static int numberOfElements=mesh.getNumberOfElements();
  static int numberOfInteriorNodes=mesh.getNumberOfInteriorNodes();
  static arrayInt globalNodesList=mesh.globalNodesList( numberOfElements );
  static vectorInt listOfInteriorNodes=mesh.getListOfInteriorNodes( numberOfInteriorNodes );
  static arrayReal globalNodesCoords=mesh.nodesCoordinates( numberOfNodes );

  // get model
  static vectorReal model=mesh.getModel( numberOfElements );

  //get infos about finite element order of approximation
  int numberOfPointsPerElement;
  if( order==1 )
    numberOfPointsPerElement=4;
  if( order==2 )
    numberOfPointsPerElement=9;
  if( order==3 )
    numberOfPointsPerElement=16;
  if( order==4 )
    numberOfPointsPerElement=25;
  if( order==5 )
    numberOfPointsPerElement=36;

  // get quadrature points and weights
  static vectorDouble quadraturePoints=Qk.gaussLobattoQuadraturePoints( order );
  static vectorDouble weights=Qk.gaussLobattoQuadratureWeights( order );
  static vectorDouble weights2D=Qk.getGaussLobattoWeights( quadraturePoints, weights );
  // get basis function and corresponding derivatives
  static arrayDouble basisFunction1D=Qk.getBasisFunction1D( order, quadraturePoints );
  static arrayDouble derivativeBasisFunction1D=Qk.getDerivativeBasisFunction1D( order, quadraturePoints );
  static arrayDouble basisFunction2D=Qk.getBasisFunction2D( quadraturePoints,
                                                            basisFunction1D,
                                                            basisFunction1D );
  static arrayDouble derivativeBasisFunction2DX=Qk.getBasisFunction2D( quadraturePoints,
                                                                       derivativeBasisFunction1D,
                                                                       basisFunction1D );
  static arrayDouble derivativeBasisFunction2DY=Qk.getBasisFunction2D( quadraturePoints,
                                                                       basisFunction1D,
                                                                       derivativeBasisFunction1D );

  static vectorReal massMatrixGlobal( numberOfNodes );
  static vectorReal yGlobal( numberOfNodes );

  Kokkos::parallel_for( numberOfNodes, [=] ( const int i )
  {
    massMatrixGlobal[i]=0;
    yGlobal[i]=0;
  } );

  // loop over mesh elements
  Kokkos::parallel_for( numberOfElements, [=, &pnGlobal] ( const int e )
  {
    // extract global coordinates of element e
    // get local to global indexes of nodes of element e
    vectorInt localToGlobal=mesh.localToGlobalNodes( e, numberOfPointsPerElement, globalNodesList );

    //get global coordinates Xi of element e
    arrayDouble Xi=mesh.getXi( numberOfPointsPerElement, globalNodesCoords, localToGlobal );

    // compute jacobian Matrix
    arrayDouble jacobianMatrix= Qk.computeJacobianMatrix( numberOfPointsPerElement, Xi,
                                                          derivativeBasisFunction2DX,
                                                          derivativeBasisFunction2DY );
    // compute determinant of jacobian Matrix
    vectorDouble detJ= Qk.computeDeterminantOfJacobianMatrix( numberOfPointsPerElement,
                                                              jacobianMatrix );
    // compute inverse of Jacobian Matrix
    arrayDouble invJacobianMatrix= Qk.computeInvJacobianMatrix( numberOfPointsPerElement,
                                                                jacobianMatrix,
                                                                detJ );
    // compute transposed inverse of Jacobian Matrix
    arrayDouble transpInvJacobianMatrix= Qk.computeTranspInvJacobianMatrix( numberOfPointsPerElement,
                                                                            jacobianMatrix,
                                                                            detJ );
    // compute  geometrical transformation matrix
    arrayDouble B=Qk.computeB( numberOfPointsPerElement, invJacobianMatrix, transpInvJacobianMatrix, detJ );

    /**
       // compute stifness and mass matrix
       vector<vector<double>> const R=Qk.gradPhiGradPhi(numberOfPointsPerElement, weights2D, B, derivativeBasisFunction2DX,
                          derivativeBasisFunction2DY);
     **/

    // compute stifness and mass matrix
    arrayDouble R=Qk.gradPhiGradPhi( numberOfPointsPerElement, order, weights2D, B, derivativeBasisFunction1D );

    // compute local mass matrix
    vectorDouble massMatrixLocal=Qk.phiIphiJ( numberOfPointsPerElement, weights2D, detJ );
    // get pnGlobal to pnLocal
    vectorReal pnLocal( numberOfPointsPerElement );
    vectorReal Y( numberOfPointsPerElement );
    for( int i=0; i<numberOfPointsPerElement; i++ )
    {
      massMatrixLocal[i]/=(model[e]*model[e]);
      pnLocal[i]=pnGlobal[localToGlobal[i]][i2];
    }
    for( int i=0; i<numberOfPointsPerElement; i++ )
    {
      Y[i]=0;
      for( int j=0; j<numberOfPointsPerElement; j++ )
      {
        Y[i]+=R[i][j]*pnLocal[j];

      }
    }
    for( int i=0; i<numberOfPointsPerElement; i++ )
    {
      int gIndex=localToGlobal[i];
      massMatrixGlobal[gIndex]+=massMatrixLocal[i];
      yGlobal[gIndex]+=Y[i];
    }
  } );

  // update pressure
  Kokkos::parallel_for( numberOfInteriorNodes, [=, &pnGlobal] ( const int i )
  {
    int I=listOfInteriorNodes[i];
    float tmp=timeSample*timeSample;
    pnGlobal[I][i1]=2*pnGlobal[I][i2]-pnGlobal[I][i1]-tmp*yGlobal[I]/massMatrixGlobal[I];
  } );
  //cout<<"pressure="<<pnGlobal[5][i1]<<endl;

  // damping terms
  static vectorReal ShGlobal( numberOfBoundaryNodes );

  for( int i=0; i<numberOfBoundaryNodes; i++ )
  {
    ShGlobal[i]=0;
  }
  // Note: this loop is data parallel.
  Kokkos::parallel_for( numberOfBoundaryFaces, [=] ( const int iFace )
  {
    vectorReal ds( order+1 );
    vectorReal Sh( order+1 );
    //get ds
    ds=Qk.computeDs( iFace, order, faceInfos, globalNodesCoords,
                     derivativeBasisFunction2DX,
                     derivativeBasisFunction2DY );
    //compute Sh and ShGlobal
    for( int i=0; i<order+1; i++ )
    {
      int gIndexFaceNode=localFaceNodeToGlobalFaceNode[iFace][i];
      Sh[i]=weights[i]*ds[i]/(model[faceInfos[iFace][0]]);
      ShGlobal[gIndexFaceNode]+=Sh[i];
    }
    /**
       cout<<"iFace="<<iFace<<endl;
       for (int i=0; i<order+1;i++)
       {
       int gIndexFaceNode=localFaceNodeToGlobalFaceNode[iFace][i];
       cout<<"   gIndex="<<gIndexFaceNode<<endl;
       cout<<"   Sh["<<i<<"]="<<Sh[i]<<endl;
       cout<<"   ShGlobal["<<gIndexFaceNode<<"]="<<ShGlobal[gIndexFaceNode]<<endl;
       }
     **/
  } );

  // update pressure @ boundaries;
  float tmp=timeSample*timeSample;
  Kokkos::parallel_for( numberOfBoundaryNodes, [=, &pnGlobal] ( const int i )
  {
    int I=listOfBoundaryNodes[i];
    float invMpSh=1/(massMatrixGlobal[I]+timeSample*ShGlobal[i]*0.5);
    float MmSh=massMatrixGlobal[I]-timeSample*ShGlobal[i]*0.5;
    pnGlobal[I][i1]=invMpSh*(2*massMatrixGlobal[I]*pnGlobal[I][i2]-MmSh*pnGlobal[I][i1]-tmp*yGlobal[I]);
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
