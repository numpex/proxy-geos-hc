//************************************************************************
// FD proxy application v.0.0.1
//
// main.cpp: this main file is simply a driver
//************************************************************************

#include"FDTDinit.hpp"

int main( int argc, char *argv[] )
{
  #ifdef USE_KOKKOS
  Kokkos::initialize(argc, argv);
  { 
  #endif
    struct FDTDGRIDS myGrids;
    struct FDTDMODELS myModels;

    // imports utility and kernel modules
    FDTDInit myInit; 
    FDTDKernel myKernel;
    FDTDUtils myFDTDUtils;

    // initialize geometry 
    myInit.init_geometry( argc, argv, myGrids );

    // initialize coefficients
    myInit.init_coefficients( myGrids, myModels );

    // initialize source 
    myInit.init_source( myModels );

    // initialize velocity and pressure models, etc
    myInit.init_models( myGrids, myModels );

    // start timer
    time_point< system_clock > startTime = system_clock::now();

    // main loop for wave propagation on each time step
    for (int itSample=0; itSample<myInit.nSamples; itSample++)
    {
       // add RHS term
       myKernel.addRHS(myGrids, itSample, myModels);

       //compute one step
       myKernel.computeOneStep( myGrids, myModels ); 

       // swap wavefields
       myKernel.swapWavefields( myGrids, myModels );

       // print infos and save wavefields
       myFDTDUtils.output(myGrids, myModels.pn, itSample);

   }

    // print timing information
    cout << "Elapsed Time : "<<( system_clock::now()-startTime ).count()/1E9 <<" seconds.\n"<<endl;

  #ifdef USE_KOKKOS
  }
  Kokkos::finalize();
  #endif
  return (0);
}

