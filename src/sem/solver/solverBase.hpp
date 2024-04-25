//************************************************************************
//   proxy application v.0.0.1
//
//  solverBase.hpp: simple 2D acoustive wave equation solver
//
//  the solverBase class servers as a base class for the SEM solver
//
//************************************************************************

#ifndef SOLVER_BASE_HPP_
#define SOLVER_BASE_HPP_

#include    "QkGL.hpp"
#include    "simpleMesh.hpp"

using namespace grid;
using namespace FE;

class solverBase
{
public:

PROXY_HOST_DEVICE solverBase(){};
PROXY_HOST_DEVICE ~solverBase(){};

  /**
   * @brief computeFEInit function:
   * init all FE components for computing mass and stiffness matrices
   * method defines here because shared by all derived classes
   */
  void computeFEInit ( const int & order,simpleMesh mesh, QkGL Qk);

  virtual void computeOneStep( const int & indexTimeStep,
                               const float & timeSample,
                               const int & order,
                               int & i1,
                               int & i2,
                               const int & numberOfRHS,
                               vectorIntView & rhsElement,
                               arrayRealView & rhsTerm,
                               arrayRealView & pnGlobal) = 0;


public:

  int i1=0, i2=1;

  // get infos from mesh
  int numberOfNodes;
  int numberOfElements;
  int numberOfPointsPerElement;
  int numberOfInteriorNodes;
  int numberOfBoundaryNodes;
  int numberOfBoundaryFaces;
  const int numberOfColors=4;

  //shared arrays
  int numberOfElementsByColor[4];
  arrayIntView listOfElementsByColor;
  arrayIntView globalNodesList;
  arrayRealView globalNodesCoords;
  vectorIntView listOfIntVieweriorNodes;
  vectorIntView listOfBoundaryNodes;
  arrayIntView faceInfos;
  arrayIntView localFaceNodeToGlobalFaceNode;
  
  // get model
  vectorRealView model;

  // get quadrature points and weights
  vectorDoubleView quadraturePoints;
  vectorDoubleView weights;
  vectorDoubleView weights2D;
  vectorDoubleView weights3D;

  // get basis function and corresponding derivatives
  arrayDoubleView basisFunction1D;
  arrayDoubleView derivativeBasisFunction1D;
  arrayDoubleView basisFunction2D;
  arrayDoubleView derivativeBasisFunction2DX;
  arrayDoubleView derivativeBasisFunction2DY;

  //shared arrays
  vectorRealView massMatrixGlobal;
  vectorRealView yGlobal;
  vectorRealView ShGlobal;
  
  simpleMesh mesh;
  QkGL Qk;
  
};
#endif //SOLVER_BASE_HPP_
