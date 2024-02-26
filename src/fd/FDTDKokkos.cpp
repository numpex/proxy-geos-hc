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
#include <vector>
#include "utils.hpp"
#include "dataType.hpp"
#include "FDTDutils.hpp"
#include "FDTDKernel.hpp"


int main( int argc, char *argv[] )
{
  Kokkos::initialize(argc,argv);
  {
    int nx= (argc > 1)? std::stoi(argv[1]) : 150;
    int ny=nx;
    int nz=nx;

    int xs=nx/2;
    int ys=ny/2;
    int zs=nz/2;

    constexpr int lx=4;
    constexpr int ly=4;
    constexpr int lz=4;
    constexpr float dx=10;
    constexpr float dy=10;
    constexpr float dz=10;

    constexpr int   sourceOrder=1;
    constexpr float f0=15.;
    constexpr float fmax=2.5*f0;
    constexpr float timeMax=1.0;

    constexpr int ncoefs=5;
    constexpr float vmin=1500;
    constexpr float vmax=4500;
    
    // init pml limits
    constexpr int ntaperx=3;
    constexpr int ntapery=3;
    constexpr int ntaperz=3;
    constexpr float hdx_2=1./(4.*dx*dx);
    constexpr float hdy_2=1./(4.*dy*dy);
    constexpr float hdz_2=1./(4.*dz*dz);
    constexpr float lambdamax=vmax/fmax;
    constexpr int ndampx=ntaperx*lambdamax/dx;
    constexpr int ndampy=ntapery*lambdamax/dy;
    constexpr int ndampz=ntaperz*lambdamax/dz;
    printf("nx=%d ny=%d nz=%d\n",nx, ny,nz);
    printf("ndampx=%d ndampy=%d ndampz=%d\n",ndampx, ndampy,ndampz);

    constexpr int x1=0;
    constexpr int x2=ndampx;
    constexpr int x3=ndampx;
    int x4=nx-ndampx;
    int x5=nx-ndampx;
    int x6=nx;
    constexpr int y1=0;
    constexpr int y2=ndampy;
    constexpr int y3=ndampy;
    int y4=ny-ndampy;
    int y5=ny-ndampy;
    int y6=ny;
    constexpr int z1=0;
    constexpr int z2=ndampz;
    constexpr int z3=ndampz;
    int z4=nz-ndampz;
    int z5=nz-ndampz;
    int z6=nz;

    // allocate vector and arrays 
    // FD coefs
    vectorReal coefx=allocateVector<vectorReal>(ncoefs);
    vectorReal coefy=allocateVector<vectorReal>(ncoefs);
    vectorReal coefz=allocateVector<vectorReal>(ncoefs);

    vectorReal::HostMirror h_coefx = Kokkos::create_mirror_view( coefx );
    vectorReal::HostMirror h_coefy = Kokkos::create_mirror_view( coefy );
    vectorReal::HostMirror h_coefz = Kokkos::create_mirror_view( coefz );

    // imports utils
    solverUtils myUtils;
    FDTDUtils myFDTDUtils;
    FDTDKernel myKernel;

    // extract FD coefs
    myFDTDUtils.init_coef(dx, h_coefx);
    myFDTDUtils.init_coef(dy, h_coefy);
    myFDTDUtils.init_coef(dz, h_coefz);

    float h_coef0=-2.*(h_coefx[1]+h_coefx[2]+h_coefx[3]+h_coefx[4]);
    h_coef0+=-2.*(h_coefy[1]+h_coefy[2]+h_coefy[3]+h_coefy[4]);
    h_coef0+=-2.*(h_coefz[1]+h_coefz[2]+h_coefz[3]+h_coefz[4]);
    printf("h_coef0=%f\n",h_coef0);

    // model
    vectorReal vp=allocateVector<vectorReal>(nx*ny*nz);
    vectorReal::HostMirror h_vp = Kokkos::create_mirror_view( vp );

    // pressure fields
    vectorReal pnp1=allocateVector<vectorReal>((nx+2*lx)*(ny+2*ly)*(nz+2*lz));
    vectorReal pn=allocateVector<vectorReal>((nx+2*lx)*(ny+2*ly)*(nz+2*lz));
    vectorReal::HostMirror h_pnp1 = Kokkos::create_mirror_view( pnp1 );
    vectorReal::HostMirror h_pn = Kokkos::create_mirror_view( pn );

    // PML arrays
    vectorReal phi=allocateVector<vectorReal>(nx*ny*nz);
    vectorReal eta=allocateVector<vectorReal>((nx+2)*(ny+2)*(nz+2));
    vectorReal::HostMirror h_phi = Kokkos::create_mirror_view( phi );
    vectorReal::HostMirror h_eta = Kokkos::create_mirror_view( eta );

    // compute time step
    float timeStep=myFDTDUtils.compute_dt_sch(vmax,h_coefx,h_coefy,h_coefz);
    float timeStep2=timeStep*timeStep;
    const int nSamples=timeMax/timeStep;
    printf("timeStep=%f\n",timeStep);

    // compute source term
    // source term
    vectorReal RHSTerm=allocateVector<vectorReal>(nSamples);
    vectorReal::HostMirror h_RHSTerm = Kokkos::create_mirror_view( RHSTerm );
    std::vector<float> sourceTerm=myUtils.computeSourceTerm(nSamples,timeStep,f0,sourceOrder);

    for(int i=0;i<nSamples;i++)
    {
      h_RHSTerm[i]=sourceTerm[i];
    }

    printf("memory used for vectra and arrays %d bytes\n",
           (nx*ny*nz+3*(nx+2*lx)*(ny+2*ly)*(nz+2*lz)+nSamples+ncoefs)*4);

    // initialize vp and pressure field
    #pragma omp parallel for collapse(3)
    for( int i=0; i<nx;i++)
    {
       for( int j=0; j<ny;j++)
       {
          for( int k=0; k<nz;k++)
          {
            h_vp[IDX3(i,j,k)]=vmin*vmin*timeStep2;
            h_phi[IDX3(i,j,k)]=0.;
          }
       }
    }
    #pragma omp parallel for collapse(3)
    for( int i=-lx; i<nx+lx;i++)
    {
       for( int j=-ly; j<ny+ly;j++)
       {
          for( int k=-lz; k<nz+lz;k++)
          {
            h_pnp1[IDX3_l(i,j,k)]=0.000001;
            h_pn[IDX3_l(i,j,k)]  =0.000001;
          }
       }
    }

    //init pml eta array
    myFDTDUtils.init_eta( nx,  ny,  nz,
		          ndampx,  ndampy, ndampz,
                          x1,  x2,  x3,  x4,  x5,  x6,
                          y1,  y2,  y3,  y4,  y5,  y6,
                          z1,  z2,  z3,  z4,  z5,  z6,
                          dx,  dy,  dz,  timeStep,
                          vmax, h_eta);

    Kokkos::deep_copy( vp, h_vp);
    Kokkos::deep_copy( pn, h_pn);
    Kokkos::deep_copy( pnp1, h_pnp1);
    Kokkos::deep_copy( RHSTerm, h_RHSTerm);

    Kokkos::deep_copy( coefx, h_coefx);
    Kokkos::deep_copy( coefy, h_coefy);
    Kokkos::deep_copy( coefz, h_coefz);
    Kokkos::deep_copy( eta, h_eta);

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    for (int itSample=0; itSample<nSamples;itSample++)
    {

      // add RHS term
      myKernel.addRHS(nx,ny,nz,lx,ly,lz,xs,ys,zs,itSample,RHSTerm,vp,pn);

      //compute one step
      myKernel.computeOneStep(nx,ny,nz,
                              lx,ly,lz,
                              x1,  x2,  x3,
                              x4,  x5,  x6,
                              y1,  y2,  y3,
                              y4,  y5,  y6,
                              z1,  z2,  z3,
                              z4,  z5,  z6,
                              h_coef0,
                              hdx_2,hdy_2,hdz_2,
                              coefx,coefy,coefz,
                              vp,phi,eta,
                              pnp1,pn);

      // swap wavefields
      myKernel.swapWavefields(nx,ny,nz,lx,ly,lz,pnp1,pn);

      // print infos and save wavefields
      if(itSample%50==0)
      {
        Kokkos::fence();
        Kokkos::deep_copy( h_pn, pn);
	printf("result1 %f\n",h_pn[IDX3_l(xs,ys,zs)]);
        //myFDTDUtils.write_io(nx,ny,nz,lx,ly,lz,0,nx,ny/2,ny/2,0,nz,pn,itSample);
      }

    }

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::time_t start_time = std::chrono::system_clock::to_time_t(start);
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "started computation at " << std::ctime(&start_time)<<std::endl
              << "finished computation at " << std::ctime(&end_time)<<std::endl
              << "elapsed time: " << elapsed_seconds.count() << "s\n";

  } 
  Kokkos::finalize();
  return (0);
}
