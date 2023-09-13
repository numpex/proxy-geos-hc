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

  /**
   * @brief computeFEInit function:
   * compute one time step of wave propagation
   */
  virtual void computeOneStep( const float & timeSample,
                               const int & order,
                               int & i1,
                               int & i2,
                               arrayReal & pnGlobal,
                               simpleMesh mesh,
                               QkGL Qk ) = 0;


protected:

  int i1=0, i2=1;
  int numberOfThreads=32;

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


#ifndef SEM_USE_OMP
  // private arrays (for Kokkos And LvArray)
  arrayInt localToGlobal;
  array3DDouble Xi;

  array3DDouble jacobianMatrix;
  arrayDouble detJ;
  array3DDouble invJacobianMatrix;
  array3DDouble transpInvJacobianMatrix;

  array3DDouble B;
  array3DDouble R;

  arrayDouble massMatrixLocal;
  arrayReal pnLocal;
  arrayReal Y;

  arrayReal ds;
  arrayReal Sh;
  arrayInt numOfBasisFunctionOnFace;
  array3DReal Js;
#endif

  //shared arrays
  vectorReal massMatrixGlobal;
  vectorReal yGlobal;
  vectorReal ShGlobal;
  
  
  
  // end init
};
#endif //SOLVER_BASE_HPP_
