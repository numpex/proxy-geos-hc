//************************************************************************
//  SEM proxy application v.0.0.1
//
//  solverOMP.hpp: simple 2D acoustive wave equation solver
//
//  the solverOMP class is derived from the solverBase class
//  with the openMP implementation of the solver
//
//************************************************************************

#ifndef SOLVER_OMP_HPP_
#define SOLVER_OMP_HPP_

#include    <iostream>
#include    <vector>
#include    <cmath>
#include    "QkGL.hpp"
#include    "simpleMesh.hpp"
#include    "solverBase.hpp"

using    namespace std;

class solverOMP : public solverBase
{

public:

  solverOMP(){};
  ~solverOMP(){};

  void computeOneStep( const float & timeSample,
                       const int & order,
                       int & i1,
                       int & i2,
                       arrayReal & pnGlobal,
                       simpleMesh mesh,
                       QkGL Qk );
};
#endif //SOLVER_OMP_HPP_
