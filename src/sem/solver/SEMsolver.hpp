//************************************************************************
//   proxy application v.0.0.1
//
//  SEMsolver.hpp: simple 2D acoustive wave equation solver
//
//  the SEMsolver class servers as a base class for the SEM solver
//
//************************************************************************

#ifndef SEM_SOLVER_HPP_
#define SEM_SOLVER_HPP_

#include "SEMQkGL.hpp"
#include "SEMmesh.hpp"
#include "utils.hpp"

class SEMsolver
{
public:

PROXY_HOST_DEVICE SEMsolver(){};
PROXY_HOST_DEVICE ~SEMsolver(){};

  /**
   * @brief computeFEInit function:
   * init all FE components for computing mass and stiffness matrices
   */
  void initFiniteElem(const int & order);

   /**
   * @brief computeOneStep function:
   * init all FE components for computing mass and stiffness matrices
   */
  void computeOneStep( const int & indexTimeStep,
                       const float & timeSample,
                       const int & order,
                       int & i1,
                       int & i2,
                       const int & numberOfRHS,
                       vectorIntView & rhsElement,
                       arrayRealView & rhsTerm,
                       arrayRealView & pnGlobal);

private:

  int i1=0, i2=1;

  // get infos from mesh
  int numberOfNodes;
  int numberOfElements;
  int numberOfPointsPerElement;
  int numberOfInteriorNodes;
  int numberOfBoundaryNodes;
  int numberOfBoundaryFaces;

  //shared arrays
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
  
};
#endif //SEM_SOLVER_HPP_
