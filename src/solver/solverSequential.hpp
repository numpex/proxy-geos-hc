//************************************************************************
//  SEM proxy application v.0.0.1
//
//  solverSequential.hpp: simple 2D acoustive wave equation solver
//
//  the solverSEQ class is derived from the solverBase class
//  with the RAJA implementation of the solver
//
//************************************************************************

#ifndef SOLVER_SEQ_HPP_
#define SOLVER_SEQ_HPP_

#include    <iostream>
#include    <vector>
#include    <cmath>
#include    "QkGL.hpp"
#include    "simpleMesh.hpp"
#include    "solverBase.hpp"

using    namespace std;

class solverSEQ : public solverBase
{

public:

  solverSEQ(){};
  ~solverSEQ(){};

  void computeOneStep( const float & timeSample,
                       const int & order,
                       int & i1,
                       int & i2,
                       arrayReal & pnGlobal,
                       simpleMesh mesh,
                       QkGL Qk );
};
#endif //SOLVER_SEQ_HPP_
