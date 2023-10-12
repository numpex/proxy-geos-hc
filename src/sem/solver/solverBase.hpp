//************************************************************************
//  SEM proxy application v.0.0.1
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
#include    "dataType.hpp"
#include    "omp.h"

using namespace grid;
using namespace FE;

class solverBase
{
public:

  solverBase(){};
  ~solverBase(){};

  /**
   * @brief computeFEInit function:
   * init all FE components for computing mass and stiffness matrices
   * method defines here because shared by all derived classes
   */
  void computeFEInit ( const int & order,
                       simpleMesh mesh,
                       QkGL Qk );

  #ifdef SEM_USE_OMP
  /**
   * @brief addRightAndSides function:
   * add right and side
   */
  void addRightAndSides( const int & timeStep,
                         const int & numberOfRHS,
                         const int & i2,
                         const float & timeSample,
                         arrayReal & pnGlobal,
                         arrayReal & rhsTerm,
                         arrayReal & rhsLocation,
                         simpleMesh mesh );
  virtual void computeOneStep( const float & timeSample,
                               const int & order,
                               int & i1,
                               int & i2,
                               arrayReal & pnGlobal,
                               simpleMesh mesh,
                               QkGL Qk ) = 0;
  #else
  virtual void computeOneStep( const int & indexTimeStep,
                               const float & timeSample,
                               const int & order,
                               int & i1,
                               int & i2,
                               const int & numberOfRHS,
                               vectorInt & rhsElement,
                               arrayReal & rhsTerm,
                               arrayReal & pnGlobal,
                               simpleMesh mesh,
                               QkGL Qk ) = 0;
  #endif


protected:

  int i1=0, i2=1;
  int numberOfThreads=omp_get_max_threads();

  // get infos from mesh
  int numberOfNodes;
  int numberOfElements;
  int numberOfInteriorNodes;
  int numberOfBoundaryNodes;
  int numberOfBoundaryFaces;
  int numberOfPointsPerElement;

  //shared arrays
  arrayInt globalNodesList;
  arrayReal globalNodesCoords;
  vectorInt listOfInteriorNodes;
  vectorInt listOfBoundaryNodes;
  arrayInt faceInfos;
  arrayInt localFaceNodeToGlobalFaceNode;
  
  // get model
  vectorReal model;

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

  //shared arrays
  vectorDouble massMatrixGlobal;
  vectorDouble yGlobal;
  vectorReal ShGlobal;
//*/
  
  
  
  // end init
};
#endif //SOLVER_BASE_HPP_
