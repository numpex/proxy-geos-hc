//************************************************************************
//  SEM proxy application v.0.0.1
//
//  semproxy.cpp: the main interface of SEM proxy application
//
//************************************************************************

#include "semproxy.hpp"
// Initialize the simulation.
void SEMProxy::init()
{
  SEM_CALIPER_MARK_BEGIN( "InitTime" );

  // get information from mesh
  numberOfNodes=myMesh.getNumberOfNodes();
  numberOfElements=myMesh.getNumberOfElements();

  // allocate arrays and vectors
  myRHSLocation=allocateArray2D<arrayReal>( myNumberOfRHS, 2 );
  myRHSTerm=allocateArray2D<arrayReal>( myNumberOfRHS, myNumSamples );
  nodeList=allocateArray2D<arrayInt>(numberOfElements,(myOrderNumber+1)*(myOrderNumber+1));
  pnGlobal=allocateArray2D<arrayReal>( numberOfNodes, 2 );

  // set number of rhs and location
  myRHSLocation(0,0)=501;
  myRHSLocation(0,1)=501;
  cout << "Source location: "<<myRHSLocation(0,0)<<", "<<myRHSLocation(0,1)<<endl;

  // get element number of source term
  myElementSource=myMesh.getElementNumberFromPoints( myRHSLocation(0,0), myRHSLocation(0,1) );
  cout <<"Element number for the source location: "<<myElementSource<<endl;

  // initialize source term
  vector< float > sourceTerm=myUtils.computeSourceTerm( myNumSamples, myTimeStep, f0, sourceOrder );
  for( int j=0; j<myNumSamples; j++ )
  {
    myRHSTerm(0,j)=sourceTerm[j];
  }

  for( int i=0; i<myNumSamples; i++ )
  {
    if( i%100==0 )
      cout<<"Sample "<<i<<"\t: sourceTerm = "<<sourceTerm[i]<<endl;
  }

  // get nodelist 
  myMesh.globalNodesList( numberOfElements, nodeList );

#ifdef SEM_USE_KOKKOS 
  rhsElement=allocateVector<vectorInt>(myNumberOfRHS);
  for( int i=0; i<myNumberOfRHS; i++ )
  {
    //extract element number for current rhs
    float x=myRHSLocation(i,0);
    float y=myRHSLocation(i,1);
    int rhsE=myMesh.getElementNumberFromPoints( x, y );
    rhsElement[i]=rhsE;
  }
#endif

#ifdef SEM_USE_RAJA
  rhsElement=allocateVector<vectorInt>(myNumberOfRHS);
  for( int i=0; i<myNumberOfRHS; i++ )
  {
    //extract element number for current rhs
    float x=myRHSLocation(i,0);
    float y=myRHSLocation(i,1);
    int rhsE=myMesh.getElementNumberFromPoints( x, y );
    rhsElement[i]=rhsE;
  }
#endif

  SEM_CALIPER_MARK_END( "InitTime" );
}


// Run the simulation.
void SEMProxy::run()
{
  SEM_CALIPER_MARK_BEGIN( "RunTime" );


  mySolver.computeFEInit( myOrderNumber, myMesh, myQk );

  for( int indexTimeStep=0; indexTimeStep<myNumSamples; indexTimeStep++ )
  {
    #ifdef SEM_USE_OMP
      mySolver.addRightAndSides( indexTimeStep, myNumberOfRHS, i2, myTimeStep, pnGlobal, myRHSTerm, myRHSLocation, myMesh );
      mySolver.computeOneStep( myTimeStep, myOrderNumber, i1, i2, pnGlobal, myMesh, myQk );
    #elif defined SEM_USE_KOKKOS
      mySolver.computeOneStep( indexTimeStep, myTimeStep, myOrderNumber, i1, i2, myNumberOfRHS, rhsElement, myRHSTerm, pnGlobal, myMesh, myQk  );
    #else 
      mySolver.computeOneStep( indexTimeStep, myTimeStep, myOrderNumber, i1, i2, myNumberOfRHS, rhsElement, myRHSTerm, pnGlobal, myMesh, myQk  );
    #endif

    //writes debugging ascii file.
    if( indexTimeStep%50==0 )
    {
      cout<<"TimeStep="<<indexTimeStep<<"\t: pnGlobal @ elementSource location "<<myElementSource
          <<" after computeOneStep = "<< pnGlobal(nodeList(myElementSource,0),0)<<endl;
      //myUtils.saveSnapShot( indexTimeStep, i1, pnGlobal, myMesh );
    }
    swap( i1, i2 );
  }

  SEM_CALIPER_MARK_END( "RunTime" );


}
