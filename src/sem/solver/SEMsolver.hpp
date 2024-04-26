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
#include "SEMdata.hpp"



class SEMsolver
{
public:

PROXY_HOST_DEVICE SEMsolver(){};
PROXY_HOST_DEVICE ~SEMsolver(){};


  /**
   * @brief computeFEInit function:
   * init all FE components for computing mass and stiffness matrices
   */
  void computeFEInit ( SEMmeshinfo &myMeshinfo, SEMmesh mesh);

   /**
   * @brief computeOneStep function:
   * init all FE components for computing mass and stiffness matrices
   */
  void computeOneStep( const int & indexTimeStep,
                       SEMmeshinfo &myMeshinfo,
                       int & i1,
                       int & i2,
                       vectorIntView & rhsElement,
                       arrayRealView & rhsTerm,
                       arrayRealView & pnGlobal);

  void outputPnValues (  const int & indexTimeStep,
                         int & i1, 
                         int & myElementSource, 
                         arrayIntView & nodeList,
                         arrayRealView & pnGlobal);

  void initFEarrays( SEMmeshinfo &myMeshinfo, SEMmesh mesh );

  void allocateFEarrays( SEMmeshinfo &myMeshinfo );

private:

  int order; 
  float tmp;
  int numberOfPointsPerElement;

  SEMQkGL myQk;
  
  //shared arrays
  arrayIntView globalNodesList;
  arrayRealView globalNodesCoords;
  vectorIntView listOfInteriorNodes;
  vectorIntView listOfIntVieweriorNodes;
  vectorIntView listOfBoundaryNodes;
  arrayIntView faceInfos;
  arrayIntView localFaceNodeToGlobalFaceNode;
  
  // get model
  vectorRealView model;

  // get quadrature points and weights
  vectorDoubleView quadraturePoints;
  vectorDoubleView weights;

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
