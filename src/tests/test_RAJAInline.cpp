//************************************************************************
//  SEM proxy application v.0.0.1
//
//  main.cpp: this main file is simply a driver
//************************************************************************

#include <iostream>
#include <cstdio>
#include <fstream>
#include <cmath>
#include <chrono>
#include <omp.h>

#include "test_RAJAInline.hpp"


void globalNodesList( const int & order, const int & ex, const int & ey, const int & nx, arrayInt const & nodesList )
{
    //int nDof=(order+1)*(order+1);
    //arrayInt nodesList( numberOfElements, nDof );

    for( int j=0; j<ey; j++ )
    {
      for( int i=0; i<ex; i++ )
      {
        int n0=i+j*ex;
        int offset=i*order+j*order*nx;
        //cout<<"element "<<n0<<endl;
        for( int k=0; k<order+1; k++ )
        {
          for( int l=0; l<order+1; l++ )
          {
            int dofLocal=l+k*(order+1);
            int dofGlobal=offset+l+k*nx;
            nodesList(n0,dofLocal)=dofGlobal;
            //cout<<nodesList[n0][dofLocal]<<", ";
          }
        }
        //cout<<endl;
      }
    }
    //return nodesList;
}


void nodesCoordinates( const int & nx, const int & ny, const int & order, const int & ex, const int & ey,
                       const float & hx , const float & hy, const int & numberOfNodes, arrayReal const & nodeCoords )
{
  //arrayReal nodeCoords( numberOfNodes, 2 );
  std::vector<float> coordX( nx );
  std::vector<float> coordY( ny );
  std::vector<float> xi( order+1 );

  switch( order )
  {
    case 1:
      xi[0]=-1.;
      xi[1]=1.;
      break;
    case 2:
      xi[0]=-1.;
      xi[1]=0.;
      xi[2]=1.;
      break;
    default:
      break;
  }
  for( int i=0; i<ex; i++ )
  {
    //cout<<"elementx "<<i<<endl;
    float x0=i*hx;
    float x1=(i+1)*hx;
    float b=(x1+x0)/2.;
    float a=b-x0;

    for( int j=0; j<order+1; j++ )
    {
      coordX[j+i*order]=a*xi[j]+b;
      //cout<<coordX[j+i*order]<<" nx= "<<nx<< " j+i*order="<<j+i*order<<", ";
    }
    //cout<<endl;
  }
  for( int i=0; i<ey; i++ )
  {
    //cout<<"elementx "<<i<<endl;
    float y0=i*hy;
    float y1=(i+1)*hy;
    float b=(y1+y0)/2.;
    float a=b-y0;
    for( int j=0; j<order+1; j++ )
    {
      coordY[j+i*order]=a*xi[j]+b;
      //cout<<coordY[j+i*order]<<" ny= "<<ny<< " j+i*order="<<j+i*order<<", ";
    }
    //cout<<endl;
  }
  for( int j=0; j<ny; j++ )
  {
    for( int i=0; i<nx; i++ )
    {
      nodeCoords(i+nx*j,0)=coordX[i];
      nodeCoords(i+nx*j,1)=coordY[j];
      //cout<<"Xi["<<i+nx*j<<"][0]="<<nodeCoords[i+nx*j][0]<<", ";
      //cout<<"Xi["<<i+nx*j<<"][1]="<<nodeCoords[i+nx*j][1]<<endl;
    }
  }
  //return nodeCoords;
}


int main( int argc, char *argv[] )
{
  const int nIter=100;

  // first omp loop
  int order=1;
  int ex=100;
  int ey=100;
  int nx=ex*order+1;
  int ny=ey*order+1;
  float hx=10;
  float hy=10;
  int numberOfElements=ex*ey;
  int numberOfPointsPerElement=(order+1)*(order+1);
  int numberOfNodes=(ex*order+1)*(ey*order+1);

  mesh m;

  arrayInt nodesList( numberOfElements, numberOfPointsPerElement );
  arrayReal nodeCoords(numberOfNodes,2);
  globalNodesList( order, ex, ey, nx, nodesList );
  nodesCoordinates( nx, ny, order, ex, ey, hx, hy, numberOfNodes, nodeCoords );

// we expect :
// 0 0 0 10 0 0 10 10 10
// 1 10 0 20 0 10 10 20 10
// 2 20 0 30 0 20 10 30 10
// 3 30 0 40 0 30 10 40 10
// 4 40 0 50 0 40 10 50 10
// 5 50 0 60 0 50 10 60 10
// 6 60 0 70 0 60 10 70 10
// 7 70 0 80 0 70 10 80 10
// 8 80 0 90 0 80 10 90 10
// 9 90 0 100 0 90 10 100 10
//
 
  arrayIntView d_nodesList=nodesList.toView();
  arrayRealView d_nodeCoords=nodeCoords.toView();
  using execPolicy=RAJA::cuda_exec<32>;
  RAJA::forall< execPolicy >( RAJA::RangeSegment( 0, nIter ), [=] LVARRAY_HOST_DEVICE ( int e )
  {
     int localToGlobal[9];
     double Xi[9][2];
     int i=m.localToGlobalNodes( e, numberOfPointsPerElement, d_nodesList, localToGlobal );
     //get global coordinates Xi of element e
     int j=m.getXi(numberOfPointsPerElement, d_nodeCoords, localToGlobal, Xi );

     if(e<10)printf("e=%d %f %f %f  %f %f %f %f %f  \n ",e,Xi[0][0],Xi[0][1],Xi[1][0],Xi[1][1],Xi[2][0],Xi[2][1],Xi[3][0],Xi[3][1]);
  });

  return (0);
}
