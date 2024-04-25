//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.hpp: the main interface of SEM proxy application
//
//************************************************************************

#ifndef SEMPROXY_HPP_
#define SEMPROXY_HPP_

#include "solver.hpp"
#include "SEMdata.hpp"

/**
 * @class SEMProxy
 */
class SEMProxy
{
public:

  /**
   * @brief Constructor of the SEMProxy class
   */
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

  SEMdata myData;

  //SEMQkGL myQk;
  QkGL myQk;

  //SEMsolver mySolver;
  
};

#endif /* SEMPROXY_HPP_ */
