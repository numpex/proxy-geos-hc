
//************************************************************************
// SEM proxy application v.0.0.1
//
// main.cpp: this main file is simply a driver
//************************************************************************

#include"FDTDinit.hpp"
#include"FDTDkernel.hpp"

int main( int argc, char *argv[] )
{

  #ifdef USE_KOKKOS
  Kokkos::initialize(argc, argv);
  { 
  #endif

   struct FDTDGRIDS myGrids;
   struct FDTDMODELS myModels;

   // imports utility and kernel modules
   FDTDInit myFDTDInit; 
   FDTDKernel myKernel;
   FDTDUtils myFDTDUtils;

   // initialize geometry 
   myFDTDInit.init_geometry( argc, argv, myGrids );
   // initialize coefficients
   myFDTDInit.init_coefficients( myGrids, myModels );

   // initialize source 
   myFDTDInit.init_source( myModels );
   // initialize velocity and pressure models, etc
   myFDTDInit.init_models( myGrids, myModels );

   // start timer
   chrono::time_point< chrono::system_clock > startTime = chrono::system_clock::now();

   // main loop for wave propagation on each time step
   for (int itSample=0; itSample<myFDTDInit.nSamples;itSample++)
   {
      // add RHS term
      myKernel.addRHS(myGrids, itSample, myModels.RHSTerm, myModels.vp, myModels.pn);                     

      //compute one step
      myKernel.computeOneStep(myGrids, myFDTDInit.coef0, 
                              myModels.coefx, myModels.coefy, myModels.coefz,
                              myModels.vp, myModels.phi, myModels.eta,
                              myModels.pnp1, myModels.pn, myModels.pnm1);


      // swap wavefields
      myKernel.swapWavefields(myGrids, myModels.pnp1, myModels.pn, myModels.pnm1);

      // print infos and save wavefields
      myFDTDUtils.output(myGrids, myModels.pn, itSample);

   }
   // print timing information
   cout << "Elapsed Time : "<<chrono::duration_cast< chrono::milliseconds >( chrono::system_clock::now()-startTime ).count()/1000.0 <<" seconds.\n"<<endl;

  #ifdef USE_KOKKOS
  }
  Kokkos::finalize();
  #endif
  return (0);
}

