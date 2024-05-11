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
  void computeFEInit ( SEMmeshinfo & myMeshinfo,
                       SEMmesh mesh );

  /**
   * @brief computeOneStep function:
   * init all FE components for computing mass and stiffness matrices
   */

  void computeOneStepNoInline ( 
                       const int & timeSample,
                        const int & order,
                        const int & nPointsPerElement,
                        const int & i1,
                        const int & i2,
                        SEMmeshinfo & myMeshinfo,
                        const arrayReal & myRHSTerm,
                        arrayReal const & myPnGlobal,
                        const vectorInt & myRhsElement );

  void computeOneStepInline ( 
                        const int & timeSample,
                        const int & order,
                        const int & nPointsPerElement,
                        const int & i1,
                        const int & i2,
                        SEMmeshinfo & myMeshinfo,
                        const arrayReal & myRHSTerm,
                        arrayReal const & myPnGlobal,
                        const vectorInt & myRhsElement );

 
  void outputPnValues ( SEMmesh mesh,
                        const int & indexTimeStep,
                        int & i1,
                        int & myElementSource,
                        const arrayReal & pnGlobal );

  void initFEarrays( SEMmeshinfo & myMeshinfo, SEMmesh mesh );

  void allocateFEarrays( SEMmeshinfo & myMeshinfo );

private:

  int order;
  SEMQkGL myQk;

  //shared arrays
  arrayInt globalNodesList;
  arrayReal globalNodesCoords;
  vectorInt listOfInteriorNodes;
  vectorInt listOfIntVieweriorNodes;
  vectorInt listOfBoundaryNodes;
  arrayInt faceInfos;
  arrayInt localFaceNodeToGlobalFaceNode;

  // get model
  vectorReal model;

  // get quadrature points and weights
  vectorDouble quadraturePoints;
  vectorDouble weights;

  // get basis function and corresponding derivatives
  arrayDouble basisFunction1D;
  arrayDouble derivativeBasisFunction1D;

  //shared arrays
  vectorReal massMatrixGlobal;
  vectorReal yGlobal;
  vectorReal ShGlobal;

  arrayInt listOfElementsByColor;
};
#endif //SEM_SOLVER_HPP_
