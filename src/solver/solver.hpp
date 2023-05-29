#ifndef SOLVER_HPP_
#define SOLVER_HPP_

#include    <iostream>
#include    <vector>
#include    <cmath>
#include    "QkGL.hpp"
#include    "simpleMesh.hpp"
#include    "dataType.hpp"

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
                       arrayReal & pnGlobal,
                       simpleMesh mesh,
                       QkGL Qk );
//vector< vector< float > > & pnGlobal,
  // add right and side
  void addRightAndSides( const int & timeStep,
                         const int & numberOfRHS,
                         const int & i2,
                         const float & timeSample,
                         arrayReal & pnGlobal,
                         arrayReal & rhsTerm,
                         arrayReal & rhsLocation,
                         simpleMesh mesh );

  int i1=0, i2=1;
};
#endif //SOLVER_HPP_
