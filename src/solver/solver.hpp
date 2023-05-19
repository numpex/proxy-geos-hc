#ifndef SOLVER_HPP_
#define SOLVER_HPP_

#include    <iostream>
#include    <vector>
#include    <cmath>
#include    "QkGL.hpp"
#include    "simpleMesh.hpp"

using    namespace std;

/*
 *  simple 2D acoustive wave equation solver
 */

class solver
{
private:
  bool preCompute=true;
public:

  solver();
  ~solver();

  // compute one time step of wave propagation
  //vector<vector<float>>
  void computeOneStep( const float & timeSample,
                       const int & order,
                       int & i1,
                       int & i2,
                       vector< vector< float > > & pnGlobal,
                       simpleMesh mesh,
                       QkGL Qk );

  // add right and side
  void addRightAndSides( const int & timeStep,
                         const int & numberOfRHS,
                         const int & i2,
                         const float & timeSample,
                         vector< vector< float > > & pnGlobal,
                         const vector< vector< float > > & rhsTerm,
                         const vector< vector< float > > & rhsLocation,
                         simpleMesh mesh );

  int i1=0, i2=1;
};
#endif //SOLVER_HPP_
