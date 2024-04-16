//************************************************************************
//  proxy application v.0.0.1
//
//  solverRaja.hpp: simple 2D acoustive wave equation solver
//
//  the solverRaja class is derived from the solverBase class
//  with the RAJA implementation of the solver
//
//************************************************************************

#ifndef SOLVER_Raja_HPP_
#define SOLVER_Raja_HPP_

#include    "solverBase.hpp"

class solverRaja : public solverBase
{

public:

LVARRAY_HOST_DEVICE  solverRaja(){};
LVARRAY_HOST_DEVICE  ~solverRaja(){};

  void computeOneStep( const int & timeStep,
                       const float & timeSample,
                       const int & order,
                       int & i1,
                       int & i2,
                      const int & numberOfRHS,
                      vectorInt & rhsElement,
                      arrayReal & rhsTerm,
                      arrayReal & pnGlobal);
};
#endif //SOLVER_Raja_HPP_
