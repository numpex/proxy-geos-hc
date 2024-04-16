//************************************************************************
//  proxy application v.0.0.1
//
//  solverKokkos.hpp: simple 2D acoustive wave equation solver
//
//  the solverKokkos class is derived from the solverBase class
//  with the KOKKOS implementation of the solver
//
//************************************************************************

#ifndef SOLVER_Kokkos_HPP_
#define SOLVER_Kokkos_HPP_

#include    "solverBase.hpp"

class solverKokkos : public solverBase
{

public:

KOKKOS_FUNCTION  solverKokkos(){};
KOKKOS_FUNCTION  ~solverKokkos(){};

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
#endif //SOLVER_Kokkos_HPP_
