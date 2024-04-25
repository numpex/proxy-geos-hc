//************************************************************************
//   proxy application v.0.0.1
//
//  semproxy.cpp: the main interface of  proxy application
//
//************************************************************************

#include "semproxy.hpp"

// Initialize the simulation.
void SEMProxy::init()
{
  // get information from mesh
  myData.getMeshInfo();

  // allocate arrays and vectors
  myData.init_arrays();

  // initialize source and RHS
  myData.init_source();
}


// Run the simulation.
void SEMProxy::run()
{
/*
  _CALIPER_MARK_BEGIN( "RunTime" );


  mySolver.computeFEInit( myOrderNumber,myMesh, myQk);

  for( int indexTimeStep=0; indexTimeStep<myNumSamples; indexTimeStep++ )
  {
      mySolver.computeOneStep( indexTimeStep, myTimeStep, myOrderNumber, i1, i2, myNumberOfRHS, rhsElement, myRHSTerm, pnGlobal);

    //writes debugging ascii file.
    if( indexTimeStep%50==0 )
    {
      cout<<"TimeStep="<<indexTimeStep<<endl;
    }
    if( indexTimeStep%100==0 )
    {
      cout<<" pnGlobal @ elementSource location "<<myElementSource
          <<" after computeOneStep = "<< pnGlobal(nodeList(myElementSource,0),i1)<<endl;
      saveSnapShot( indexTimeStep, i1, pnGlobal, myMesh );
    }
    swap( i1, i2 );
  }

  _CALIPER_MARK_END( "RunTime" );

*/
}

/*
void SEMProxy::saveSnapShot( const int indexTimeStep, const int i1, arrayRealView pnGlobal, simpleMesh mesh )
{

    int nx=mesh.getNx();
    int ny=mesh.getNy();
    int nz=mesh.getNz();
    float dx=mesh.getDx();
    float dy=mesh.getDy();
    float dz=mesh.getDz();
    int numberOfNodes=nx*nz;
    std::vector<float> inputVector( numberOfNodes );
    int offset=(ny==1?offset=0:offset=nx*nz*(ny/2-1));
    printf(" nx, ny/2-1, nz %d %d %d offset=%d\n",nx,ny/2-1,nz, offset);
    for( int i = offset; i< offset+numberOfNodes; i++ )
    {
      inputVector[i-offset]=pnGlobal(i,i1);
    }
    std::vector<std::vector<float>> grid=mesh.projectToGrid( numberOfNodes, inputVector );
    fstream snapFile;
    string snapNumber = "snapshot"+to_string( indexTimeStep );
    snapFile.open( snapNumber, ios::out| ios::trunc );
    std::cout<<"nx="<<nx<<" ny="<<ny<<std::endl;
    for( int i=0; i<nx; i++ )
    {
      for( int j=0; j<nz; j++ )
      {
        snapFile<<i*dx<<" "<<j*dz<<" " <<grid[i][j]<<endl;
      }
    }
    snapFile.close();
}
*/
