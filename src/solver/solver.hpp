#ifndef SOLVER_HPP_
#define SOLVER_HPP_

#include    <iostream>
#include    <vector>
#include    <cmath>
#include    "dataType.hpp"
#include    "QkGL.hpp"
#include    "simpleMesh.hpp"

using    namespace std;


/*
 *  simple 2D acoustive wave equation solver
 */

class solver
{
private:
  bool preCompute=true;
public:

  solver();
  ~solver();

  // init all FE components for computing mass and stiffness matrices
  // method defines here because shared by all solverXXX.cpp
  void computeFEInit(const int order,
                    simpleMesh mesh,
                    QkGL Qk){
  // get infos from mesh
  numberOfNodes=mesh.getNumberOfNodes();
  numberOfElements=mesh.getNumberOfElements();
  numberOfInteriorNodes=mesh.getNumberOfInteriorNodes();
  globalNodesList=mesh.globalNodesList( numberOfElements );
  listOfInteriorNodes=mesh.getListOfInteriorNodes( numberOfInteriorNodes );
  globalNodesCoords=mesh.nodesCoordinates( numberOfNodes );

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
  // get model
  model=mesh.getModel( numberOfElements );

  // get quadrature points and weights
  quadraturePoints=Qk.gaussLobattoQuadraturePoints( order );
  weights=Qk.gaussLobattoQuadratureWeights( order );
  weights2D=Qk.getGaussLobattoWeights( quadraturePoints,weights );
  // get basis function and corresponding derivatives
  basisFunction1D=Qk.getBasisFunction1D( order, quadraturePoints );
  derivativeBasisFunction1D=Qk.getDerivativeBasisFunction1D( order, quadraturePoints );
  basisFunction2D=Qk.getBasisFunction2D( quadraturePoints,
                                                            basisFunction1D,
                                                            basisFunction1D );
  derivativeBasisFunction2DX=Qk.getBasisFunction2D( quadraturePoints,
                                                                       derivativeBasisFunction1D,
                                                                       basisFunction1D );
  derivativeBasisFunction2DY=Qk.getBasisFunction2D( quadraturePoints,
                                                                       basisFunction1D,
                                                                       derivativeBasisFunction1D );

 
};

  // compute one time step of wave propagation
  //vector<vector<float>>
  void computeOneStep( const float & timeSample,
                       const int & order,
                       int & i1,
                       int & i2,
                       arrayReal & pnGlobal,
                       simpleMesh mesh,
                       QkGL Qk );
//vector< vector< float > > & pnGlobal,
  // add right and side
  void addRightAndSides( const int & timeStep,
                         const int & numberOfRHS,
                         const int & i2,
                         const float & timeSample,
                         arrayReal & pnGlobal,
                         arrayReal & rhsTerm,
                         arrayReal & rhsLocation,
                         simpleMesh mesh );

  int i1=0, i2=1;

private:
  // get infos from mesh
  int numberOfNodes;
  int numberOfElements;
  int numberOfInteriorNodes;
  arrayInt globalNodesList;
  vectorInt listOfInteriorNodes;
  arrayReal globalNodesCoords;

  // get model
  vectorReal   model;

  //get infos about finite element order of approximation
  int numberOfPointsPerElement;
  

  // get quadrature points and weights
  vectorDouble quadraturePoints;
  vectorDouble weights;
  vectorDouble weights2D;
  // get basis function and corresponding derivatives
  arrayDouble basisFunction1D;
  arrayDouble derivativeBasisFunction1D;
  arrayDouble basisFunction2D;
  arrayDouble derivativeBasisFunction2DX;
  arrayDouble derivativeBasisFunction2DY;

  // end init (shoulb be
};
#endif //SOLVER_HPP_
