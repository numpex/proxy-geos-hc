//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.hpp: the main interface of SEM proxy application
//
//************************************************************************

#ifndef SEMPROXY_HPP_
#define SEMPROXY_HPP_

#include "solver.hpp"
#include "utils.hpp"

/**
 * @class SEMProxy
 */
class SEMProxy
{
public:

  /**
   * @brief Constructor of the SEMProxy class
   */
  SEMProxy(){};

  /**
   * @brief Destructor of the SEMProxy class
   */
  ~SEMProxy(){};

  /**
   * @brief Initialize the simulation.
   * @post run()
   */
  void init();

  /**
   * @brief Run the simulation.
   * @pre This must be called after init()
   * @post Optional printout performance resutls
   */
  void run();

  void saveSnapShot( const int indexTimeStep, const int i1, arrayRealView pnGlobal, simpleMesh mesh );

protected:

  int i1=0;
  int i2=1;
  int myElementSource;

  int numberOfNodes;
  int numberOfElements;

  const int myNumberOfRHS=1;
  int   sourceOrder=1;
  float f0=10.;
  float myTimeMax=1.0;
  float myTimeStep=0.001;

  const int myNumSamples=myTimeMax/myTimeStep;
  const int myOrderNumber=1;
  
  // initialize mesh
  simpleMesh  myMesh {200, 200,200, 2000, 2000, 2000, myOrderNumber};

  arrayReal myRHSLocation;
  arrayReal myRHSTerm;
  arrayInt nodeList;
  arrayReal pnGlobal;
  vectorInt rhsElement;

  QkGL myQk;
  SolverUtils myUtils;

  SOLVER mySolver;

};

#endif /* SEMPROXY_HPP_ */
