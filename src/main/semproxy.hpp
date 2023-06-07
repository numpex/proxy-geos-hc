//************************************************************************
//  SEM proxy application v.0.0.1
//
//  semproxy.hpp: the main interface of SEM proxy application
//
//************************************************************************

#ifndef SEMPROXY_HPP_
#define SEMPROXY_HPP_

#include <iostream>
#include <vector>
#include <cmath>

#include <fstream>
#include "QkGL.hpp"
#include "simpleMesh.hpp"
#include "utils.hpp"
#include "commonMacro.hpp"

/**
 * @class SEMProxy
 */
class SEMProxy
{
public:

  /**
   * @brief Constructor of the SEMProxy class
   */
  //SEMProxy();
  SEMProxy():myRHSLocation(myNumberOfRHS, 2), myRHSTerm(myNumberOfRHS, myNumSamples){};
  /**
   * @brief Destructor of the SEMProxy class
   */
  ~SEMProxy();

  /**
   * @brief Initialize mesh for the simulation.
   * @post run()
   */
  void initMesh();

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

private:

  float myTimeMax=1.0;
  float myTimeStep=0.001;
  const int myNumSamples=myTimeMax/myTimeStep;
  const int myOrderNumber=2;

  // element number of source term
  int myElementSource;
  const int myNumberOfRHS=1;

  arrayReal myRHSLocation;
  arrayReal myRHSTerm;

  simpleMesh const myMesh {50, 50, 1000, 1000, myOrderNumber};
  //simpleMesh const myMesh {100, 100, 2000, 2000, myOrderNumber};
  int i1=0;
  int i2=1;
  int numberOfNodes;
  int numberOfElements;

  QkGL myQk;
  solverUtils myUtils;

#if defined(SEM_USE_RAJA)
  solverRaja mySolver;
#elif defined(SEM_USE_OMP)
  solverOMP mySolver;
#else
  solverSEQ mySolver;
#endif

};

#endif /* SEMPROXY_HPP_ */
