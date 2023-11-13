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
#include    "dataType.hpp"
#include    "omp.h"

using namespace grid;
using namespace FE;

class solverBase
{
public:

  solverBase(){};
#ifdef USE_RAJA
  LVARRAY_HOST_DEVICE ~solverBase(){};
#elif defined USE_KOKKOS
  KOKKOS_FUNCTION ~solverBase(){};
#else
  ~solverBase(){};
#endif

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
                               vectorInt & rhsElement,
                               arrayReal & rhsTerm,
                               arrayReal & pnGlobal) = 0;


public:

  int i1=0, i2=1;

  // get infos from mesh
  int numberOfNodes;
  int numberOfElements;
  int numberOfPointsPerElement;
  int numberOfInteriorNodes;
  int numberOfBoundaryNodes;
  int numberOfBoundaryFaces;

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
  
  simpleMesh mesh;
  QkGL Qk;
  
};
#endif //SOLVER_BASE_HPP_
