//************************************************************************
//  proxy application v.0.0.1
//
//  solver.hpp: 2D or 3D acoustive wave equation solver
//
//  the solver class is derived from the solverBase class
//  with the RAJA, Kokkos, OMP, or Sequential implementation of the solver
//
//************************************************************************

#ifndef SOLVER_HPP_
#define SOLVER_HPP_

#include    "solverBase.hpp"


class SOLVER: public solverBase
{

public:

  PROXY_HOST_DEVICE SOLVER(){};
  PROXY_HOST_DEVICE ~SOLVER(){};

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


#endif //SOLVER_HPP_
