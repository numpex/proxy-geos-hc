#ifndef SEMDATA_HPP_
#define SEMDATA_HPP_

#include "utils.hpp"
#include "simpleMesh.hpp"

struct SEMdata
{
  const int i1=0;
  const int i2=1;

  const float f0=10.;
  const float myTimeMax=1.;
  const float myTimeStep=0.001;
  const int myNumSamples=myTimeMax/myTimeStep;

  const int sourceOrder=1;
  const int myNumberOfRHS=1;
  const int myOrderNumber=3;


  int myElementSource;
  int numberOfNodes;
  int numberOfElements;
  int numberOfPointsPerElement;
 
  SolverUtils myUtils;

  // initialize mesh
  simpleMesh myMesh {200, 200,200, 2000, 2000, 2000, myOrderNumber};

  #ifdef USE_RAJA
  arrayReal h_myRHSLocation;
  arrayReal h_myRHSTerm;
  vectorInt h_rhsElement;
  arrayReal h_pnGlobal; 
  arrayInt  h_nodeList;     
  #endif
  arrayRealView myRHSLocation;
  arrayRealView myRHSTerm;
  vectorIntView rhsElement;
  arrayRealView pnGlobal;     
  arrayIntView  nodeList;

  void init_source() 
  {
    // set number of rhs and location
    myRHSLocation(0,0)=1001;
    myRHSLocation(0,1)=0;//1001;
    myRHSLocation(0,2)=1001;
    cout << "Source location: "<<myRHSLocation(0,0)<<", "<<myRHSLocation(0,1)<<", "<<myRHSLocation(0,2)<<endl;

    // initialize source term
    vector< float > sourceTerm=myUtils.computeSourceTerm( myNumSamples, myTimeStep, f0, sourceOrder );
    for( int j=0; j<myNumSamples; j++ )
    {
      myRHSTerm(0,j)=sourceTerm[j];
      if( j%100==0 ) cout<<"Sample "<<j<<"\t: sourceTerm = "<<sourceTerm[j]<<endl;
    }
  
    for( int i=0; i<myNumberOfRHS; i++ )
    {
      //extract element number for current rhs
      float x=myRHSLocation(i,0);
      float y=myRHSLocation(i,1);
      float z=myRHSLocation(i,2);
      int rhsE=myMesh.getElementNumberFromPoints( x, y,z );
      rhsElement[i]=rhsE;
      printf(" rhsElement=%d\n",rhsElement[i]);
    }
  }
  
  void getMeshInfo()
  {
    numberOfNodes=myMesh.getNumberOfNodes();
    numberOfElements=myMesh.getNumberOfElements();
    numberOfPointsPerElement=myMesh.getNumberOfPointsPerElement();
    myElementSource=myMesh.getElementNumberFromPoints( myRHSLocation(0,0), myRHSLocation(0,1),myRHSLocation(0,2) );
  }
  
  void init_arrays()
  {
    #ifdef USE_RAJA
    h_myRHSLocation=allocateArray2D<arrayReal>( myNumberOfRHS, 3 );
    h_myRHSTerm=allocateArray2D<arrayReal>( myNumberOfRHS, myNumSamples );
    h_rhsElement=allocateVector<vectorInt>(myNumberOfRHS);
    h_nodeList=allocateArray2D<arrayInt>(numberOfElements,numberOfPointsPerElement);
    h_pnGlobal=allocateArray2D<arrayReal>( numberOfNodes, 2 );
  
    myRHSLocation =  h_myRHSLocation.toView();  
    myRHSTerm = h_myRHSTerm.toView();
    rhsElement = h_rhsElement.toView();
    pnGlobal = h_pnGlobal.toView(); 
    nodeList = h_nodeList.toView();     
    #else
    myRHSLocation=allocateArray2D<arrayRealView>( myNumberOfRHS, 3 );
    myRHSTerm=allocateArray2D<arrayRealView>( myNumberOfRHS, myNumSamples );
    rhsElement=allocateVector<vectorIntView>(myNumberOfRHS);
    nodeList=allocateArray2D<arrayIntView>(numberOfElements,numberOfPointsPerElement);
    pnGlobal=allocateArray2D<arrayRealView>( numberOfNodes, 2 );
    #endif
  
    myMesh.globalNodesList( numberOfElements, nodeList );
  }
};

#endif //SEMDATA_HPP_
