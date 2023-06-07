//************************************************************************
//  SEM proxy application v.0.0.1
//
//  solverRaja.hpp: simple 2D acoustive wave equation solver
//
//  the solverRaja class is derived from the solverBase class
//  with the RAJA implementation of the solver
//
//************************************************************************

#ifndef SOLVER_Raja_HPP_
#define SOLVER_Raja_HPP_

#include    <iostream>
#include    <vector>
#include    <cmath>
#include    "QkGL.hpp"
#include    "simpleMesh.hpp"
#include    "solverBase.hpp"

using    namespace std;

class solverRaja : public solverBase
{

public:

  solverRaja(){};
  ~solverRaja(){};

  void computeOneStep( const float & timeSample,
                       const int & order,
                       int & i1,
                       int & i2,
                       arrayReal & pnGlobal,
                       simpleMesh mesh,
                       QkGL Qk );
};
#endif //SOLVER_Raja_HPP_
