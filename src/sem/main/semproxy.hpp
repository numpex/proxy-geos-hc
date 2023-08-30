//************************************************************************
//  SEM proxy application v.0.0.1
//
//  semproxy.hpp: the main interface of SEM proxy application
//
//************************************************************************

#ifndef SEMPROXY_HPP_
#define SEMPROXY_HPP_

#include <iostream>
#include <fstream>

#include "commonMacro.hpp"
#include "QkGL.hpp"
#include "simpleMesh.hpp"
#include "utils.hpp"
#include "dataType.hpp"

/**
 * @class SEMProxy
 */
class SEMProxy
{
public:

  /**
   * @brief Constructor of the SEMProxy class
   */
  //SEMProxy(): myRHSLocation( myNumberOfRHS, 2 ), myRHSTerm( myNumberOfRHS, myNumSamples ){};
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

protected:

  float myTimeMax=1.0;
  float myTimeStep=0.001;
  const int myNumSamples=myTimeMax/myTimeStep;
  const int myOrderNumber=2;

  // element number of source term
  int myElementSource;
  const int myNumberOfRHS=1;
  float f0=15.;
  int sourceOrder=1;

  arrayReal myRHSLocation;
  arrayReal myRHSTerm;


  simpleMesh  myMesh {50, 50, 1000, 1000, myOrderNumber};
  //simpleMesh const myMesh {2, 2, 4000, 4000, myOrderNumber};
  int i1=0;
  int i2=1;
  int numberOfNodes;
  int numberOfElements;

  arrayInt nodeList;
  //nodeList=allocateArray2D<arrayInt>(numberOfElements,(myOrderNumber+1)*(myOrderNumber+1));
  //myMesh.globalNodesList( numberOfElements, nodeList );

  arrayReal pnGlobal;
  //pnGLobal=allocateArray2D<arrayReal>( numberOfNodes, 2 );

  QkGL myQk;
  solverUtils myUtils;

#if defined(SEM_USE_RAJA)
  solverRaja mySolver;
#elif defined(SEM_USE_OMP)
  solverOMP mySolver;
#elif defined(SEM_USE_KOKKOS)
  solverKokkos mySolver;
#else
  solverSEQ mySolver;
#endif

};

#endif /* SEMPROXY_HPP_ */
