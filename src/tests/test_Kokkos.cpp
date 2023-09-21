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
#include <Kokkos_Macros.hpp>
#include <Kokkos_Core.hpp>
#include <omp.h>

//#define MemSpace Kokkos::CudaUVMSpace
#define MemSpace Kokkos::SharedSpace
//#define MemSpace Kokkos::HostSpace


using ExecSpace = MemSpace::execution_space;
using range_policy = Kokkos::RangePolicy<ExecSpace>;

/*
typedef Kokkos::View<int*, Kokkos::LayoutRight, MemSpace>   vectorInt;
typedef Kokkos::View<int**, Kokkos::LayoutRight, MemSpace>   arrayInt;
typedef Kokkos::View<float**, Kokkos::LayoutRight, MemSpace>   arrayReal;
typedef Kokkos::View<double**, Kokkos::LayoutRight, MemSpace>   arrayDouble;
*/
typedef Kokkos::View<int*,MemSpace>   vectorInt;
typedef Kokkos::View<int**,MemSpace>   arrayInt;
typedef Kokkos::View<float**,MemSpace>   arrayReal;
typedef Kokkos::View<double**,MemSpace>   arrayDouble;
typedef Kokkos::View<double***,MemSpace>   array3DDouble;

template<class T>
  T allocateVector(int n1)
  {
     std::cout<<"allocate vector of size "<<n1<<std::endl;
     T vect("v",n1);
    return vect;
  }
template<class T>
  T allocateArray2D(int n1, int n2)
  {
    std::cout<<"allocate array of size "<<n1<<", "<<n2<<std::endl;
    T array("a",n1, n2);
    return array;
  }


KOKKOS_FUNCTION  int  compute(int const n1, int const n2, arrayReal const & arrayViewIn, arrayReal const & arrayViewOut)
{
  for (int i=0; i<n1;i++)
    {
      for ( int j=0; j<n2;j++)
      {
        arrayViewOut(i,j)=sqrt(2*i*j)*arrayViewIn(i,j);
      }
    }
  return(0);
}

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
    case 3:
      static constexpr double sqrt5 = 2.2360679774997897;
      xi[0] = -1.0;
      xi[1] = -1./sqrt5;
      xi[2] = 1./sqrt5;
      xi[3] = 1.;
      break;
    case 4:
      static constexpr double sqrt3_7 = 0.6546536707079771;
      xi[0] = -1.0;
      xi[1] = -sqrt3_7;
      xi[2] = 0.0;
      xi[3] = sqrt3_7;
      xi[4] = 1.0;
      break;
    case 5:
      static constexpr double sqrt__7_plus_2sqrt7__ = 3.50592393273573196;
      static constexpr double sqrt__7_mins_2sqrt7__ = 1.30709501485960033;
      static constexpr double sqrt_inv21 = 0.218217890235992381;
      xi[0] = -1.0;
      xi[1] = -sqrt_inv21*sqrt__7_plus_2sqrt7__;
      xi[2] = -sqrt_inv21*sqrt__7_mins_2sqrt7__;
      xi[3] = sqrt_inv21*sqrt__7_mins_2sqrt7__;
      xi[4] = sqrt_inv21*sqrt__7_plus_2sqrt7__;
      xi[5] = 1.0;
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

KOKKOS_FUNCTION int localToGlobalNodes( const int & threadId,const int & e,
                                        const int & nPointsPerElement, arrayInt const & nodesList,
                                        arrayInt const & localToGlobal)
{
  //vectorInt localToGlobal( nPointsPerElement );
  for( int i=0; i<nPointsPerElement; i++ )
  {
    localToGlobal(threadId,i)=nodesList(e,i);
    //if(e==0)printf("i=%d %d\n",i,localToGlobal(e,i));
  }
  //return localToGlobal;
  return 0;
}

KOKKOS_FUNCTION int getXi( const int & e, const int & numberOfPointsPerElement, arrayReal const & globalNodesCoords,
                           arrayInt const & localToGlobal , array3DDouble const & Xi)
{
  //arrayDouble Xi( numberOfPointsPerElement, 2 );
  for( int i=0; i<numberOfPointsPerElement; i++ )
  {
    Xi(e,i,0)=globalNodesCoords(localToGlobal(e,i),0);
    Xi(e,i,1)=globalNodesCoords(localToGlobal(e,i),1);
  }
  //if(e==0)printf("Xi %f %f %f %f \n",Xi(e,0,0),Xi(e,0,1),Xi(e,1,0),Xi(e,1,1));
  //if(e==1)printf("Xi %f %f %f %f \n",Xi(e,0,0),Xi(e,0,1),Xi(e,1,0),Xi(e,1,1));
  //return Xi;
  return 0;
}


int main( int argc, char *argv[] )
{
  Kokkos::initialize(argc,argv);
  {
  const int nIter=100;
  const int n1=1000;
  const int n2=1000;

  arrayReal arrayViewIn;
  arrayViewIn=allocateArray2D<arrayReal>(n1,n2);
  arrayReal arrayViewOut;
  arrayViewOut=allocateArray2D<arrayReal>(n1,n2);

  //arrayReal::HostMirror h_arrayViewIn=Kokkos::create_mirror_view(arrayViewIn);
  //arrayReal::HostMirror h_arrayViewOut=Kokkos::create_mirror_view(arrayViewOut);

  for (int i=0; i<n1;i++)
  {
    for ( int j=0; j<n2;j++)
    {
      //h_arrayViewIn(i,j)=sqrt(2*i*j);
      //h_arrayViewOut(i,j)=0;
      arrayViewIn(i,j)=sqrt(2*i*j);
      arrayViewOut(i,j)=0;
    }
  }

  //Kokkos::deep_copy(arrayViewIn, h_arrayViewIn);
  //Kokkos::deep_copy(arrayViewOut, h_arrayViewOut);
   // Timer products.
  Kokkos::Timer timer;
  std::cout<<"start loop1 \n";
  Kokkos::parallel_for("loop", nIter , KOKKOS_LAMBDA  ( int e )
    {   for (int i=0; i<n1;i++)
        {
           for ( int j=0; j<n2;j++)
           {
              arrayViewOut(i,j)=sqrt(2*i*j)*arrayViewIn(i,j);
           }
        }
   });
  double time=timer.seconds();
  std::cout << "Elapsed Time Kokkos Cuda loop : "<<time <<" seconds.\n"<<std::endl;
  //Kokkos::deep_copy(h_arrayViewOut,arrayViewOut);
  //std::cout << h_arrayViewOut(500,500)<<std::endl;
  Kokkos::fence();
  std::cout << arrayViewOut(500,500)<<std::endl;

  // test call subroutine
  // -------------------------------
  // Timer products.
  for (int i=0; i<n1;i++)
  {
    for ( int j=0; j<n2;j++)
    {
      //h_arrayViewIn(i,j)=sqrt(2*i*j)*i;
      //h_arrayViewOut(i,j)=0;
      arrayViewIn(i,j)=sqrt(2*i*j)*i;
      arrayViewOut(i,j)=0;
    }
  }
  Kokkos::Timer timer1;
  //Kokkos::deep_copy(arrayViewIn, h_arrayViewIn);
  //Kokkos::deep_copy(arrayViewOut, h_arrayViewOut);
  Kokkos::parallel_for("routine", nIter , KOKKOS_LAMBDA ( int e)
  {
     int i=compute(n1,n2,arrayViewIn,arrayViewOut) ;
  });
  double time1=timer1.seconds();
  std::cout << "Elapsed Time Cuda function call  loop : "<<time1 <<" seconds.\n"<<std::endl;
  // the deep copy is mandatory if not we do not move back to host the values of arrayViewOut
  //Kokkos::deep_copy(h_arrayViewOut,arrayViewOut);
  //std::cout << h_arrayViewOut(500,500)<<std::endl;
  // SYNCHRONIZE 
  Kokkos::fence();
  std::cout << arrayViewOut(500,500)<<std::endl;

  // test OMP
  // -------------------------------
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

  printf("ici \n");
  int nthreads=omp_get_max_threads();
  printf("nthreads %d\n",nthreads);
  arrayInt nodesList("nodeList", numberOfElements, numberOfPointsPerElement );
  arrayReal nodeCoords("nodeCoords",numberOfNodes,2);
  globalNodesList( order, ex, ey, nx, nodesList );
  nodesCoordinates( nx, ny, order, ex, ey, hx, hy, numberOfNodes, nodeCoords );
//
// test 1 omp parallel loop. Xi and localToGlobal are writable by all threads concurrently
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
// in this case results are incorrect. 

  std::cout<<"test1 omp parallel for\n";
  arrayInt localToGlobal("localToGLobal",nthreads,numberOfPointsPerElement);
  array3DDouble Xi("Xi", nthreads,numberOfPointsPerElement, 2 );
  #pragma omp parallel for
  for ( int e=0; e<numberOfElements;e++)
  {
     int threadId=omp_get_thread_num();
     int i=localToGlobalNodes( threadId,e, numberOfPointsPerElement, nodesList, localToGlobal );

    //get global coordinates Xi of element e
    int j=getXi(threadId, numberOfPointsPerElement, nodeCoords, localToGlobal, Xi );
    if(e<10)
    {
            std::cout<<threadId<<" "<<Xi(threadId,0,0)<<" "<<Xi(threadId,0,1)<<" ";
            std::cout<<Xi(threadId,1,0)<<" "<<Xi(threadId,1,1)<<" ";
            std::cout<<Xi(threadId,2,0)<<" "<<Xi(threadId,2,1)<<" ";
            std::cout<<Xi(threadId,3,0)<<" "<<Xi(threadId,3,1)<<std::endl;
    }
  }

  std::cout<<"test2 omp parallel followed with omp for\n";
//
// test 2 omp parallel block and omp for loop. Xi and localToGlobal are writable by all threads concurrently
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
// in this case results are correct. 
#pragma omp parallel 
  {
    arrayInt localToGlobalP("localToGLobalP",nthreads,numberOfPointsPerElement);
    array3DDouble XiP("XiP", nthreads,numberOfPointsPerElement, 2 );

# pragma omp for
    for ( int e=0; e<numberOfElements;e++)
    {
       int threadId=0;
       int i=localToGlobalNodes( threadId,e, numberOfPointsPerElement, nodesList, localToGlobalP );
      //get global coordinates Xi of element e
      int j=getXi(threadId, numberOfPointsPerElement, nodeCoords, localToGlobalP, XiP );
      if(e<10)
      {
            std::cout<<e<<" "<<XiP(threadId,0,0)<<" "<<XiP(threadId,0,1)<<" ";
            std::cout<<XiP(threadId,1,0)<<" "<<XiP(threadId,1,1)<<" ";
            std::cout<<XiP(threadId,2,0)<<" "<<XiP(threadId,2,1)<<" ";
            std::cout<<XiP(threadId,3,0)<<" "<<XiP(threadId,3,1)<<std::endl;
      }
    }
  } // end omp parrallel
  std::cout<<"test3 Kokkos loop rangePolicy\n";
//
// test 3
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
// in this case results are correct. 
  nthreads=2;
  arrayInt localToGlobalK("localToGLobalK",nthreads,numberOfPointsPerElement);
  array3DDouble XiK("XiK", nthreads,numberOfPointsPerElement, 2 );
  array3DDouble XiKOut("XiKOut", nthreads,numberOfPointsPerElement, 2 );


  // Kokkos::parallel_for( numberOfElements , KOKKOS_LAMBDA ( int e) 
  typedef Kokkos::TeamPolicy<ExecSpace>::member_type member_type;
  // Create an instance of the policy
  Kokkos::TeamPolicy<ExecSpace> policy ((numberOfElements+nthreads-1)/nthreads, nthreads );
  Kokkos::parallel_for (policy, KOKKOS_LAMBDA (member_type team_member)
  {
     int e=team_member.league_rank () * team_member.team_size () + team_member.team_rank ();
     int threadId= team_member.team_rank ();
     int i=localToGlobalNodes( threadId,e, numberOfPointsPerElement, nodesList, localToGlobalK );
     //get global coordinates Xi of element e
     int j=getXi(threadId, numberOfPointsPerElement, nodeCoords, localToGlobalK, XiK );

     team_member.team_barrier();
     printf("%d %f %f %f  %f %f %f %f %f  \n ",threadId,XiK(threadId,0,0),XiK(threadId,0,1),XiK(threadId,1,0),XiK(threadId,1,1),XiK(threadId,2,0),XiK(threadId,2,1),XiK(threadId,3,0),XiK(threadId,3,1));
  });

  Kokkos::fence();
}
  Kokkos::finalize();
  return (0);
}
