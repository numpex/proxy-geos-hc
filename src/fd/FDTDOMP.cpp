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
    int nx=std::stoi(argv[1]);
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
    constexpr int   xs=nx/2;
    constexpr int   ys=ny/2;
    constexpr int   zs=nz/2;
    constexpr float f0=15.;
    constexpr float fmax=2.5*f0;
    constexpr float timeMax=0.8;

    constexpr int ncoefs=5;
    constexpr float vmax=1500;
    
    // imports utils
    solverUtils myUtils;
    FDTDUtils myFDTDUtils;
    FDTDKernel myKernel;

    // init pml limits
    constexpr int ntaperx=3;
    constexpr int ntapery=3;
    constexpr int ntaperz=3;
    constexpr float hdx_2=1./(4.*dx*dx);
    constexpr float hdy_2=1./(4.*dy*dy);
    constexpr float hdz_2=1./(4.*dz*dz);
    constexpr float lambdamax=vmax/fmax;
    constexpr int ndampx=ntaperx*lambdamax/dx;
    constexpr int ndampy=ntaperx*lambdamax/dx;
    constexpr int ndampz=ntaperx*lambdamax/dx;
    printf("ndampx=%d ndampy=%d ndampz=%d\n",ndampx, ndampy,ndampz);
    constexpr int x1=0;
    constexpr int x2=ndampx;
    constexpr int x3=ndampx;
    constexpr int x4=nx-ndampx;
    constexpr int x5=nx-ndampx;
    constexpr int x6=nx;
    constexpr int y1=0;
    constexpr int y2=ndampy;
    constexpr int y3=ndampy;
    constexpr int y4=ny-ndampy;
    constexpr int y5=ny-ndampy;
    constexpr int y6=ny;
    constexpr int z1=0;
    constexpr int z2=ndampz;
    constexpr int z3=ndampz;
    constexpr int z4=nz-ndampz;
    constexpr int z5=nz-ndampz;
    constexpr int z6=nz;

    // allocate vector and arrays 
    // FD coefs
    vectorReal coefx=allocateVector<vectorReal>(ncoefs);
    vectorReal coefy=allocateVector<vectorReal>(ncoefs);
    vectorReal coefz=allocateVector<vectorReal>(ncoefs);
    // model
    vectorReal vp=allocateVector<vectorReal>(nx*ny*nz);
    // pressure fields
    vectorReal pnp1=allocateVector<vectorReal>((nx+2*lx)*(ny+2*ly)*(nz+2*lz));
    vectorReal pn=allocateVector<vectorReal>((nx+2*lx)*(ny+2*ly)*(nz+2*lz));
    vectorReal pnm1=allocateVector<vectorReal>((nx+2*lx)*(ny+2*ly)*(nz+2*lz));
    // PML arrays
    vectorReal phi=allocateVector<vectorReal>(nx*ny*nz);
    vectorReal eta=allocateVector<vectorReal>((nx+2)*(ny+2)*(nz+2));

    // extract FD coefs
    myFDTDUtils.init_coef(dx, coefx);
    myFDTDUtils.init_coef(dy, coefy);
    myFDTDUtils.init_coef(dz, coefz);

    double coef0=-2.*(coefx[1]+coefx[2]+coefx[3]+coefx[4]);
    coef0+=-2.*(coefy[1]+coefy[2]+coefy[3]+coefy[4]);
    coef0+=-2.*(coefz[1]+coefz[2]+coefz[3]+coefz[4]);

    // compute time step
    float timeStep=myFDTDUtils.compute_dt_sch(vmax,coefx,coefy,coefz);
    float timeStep2=timeStep*timeStep;
    const int nSamples=timeMax/timeStep;
    printf("timeStep=%f\n",timeStep);
    // compute source term
    // source term
    vectorReal RHSTerm=allocateVector<vectorReal>(nSamples);
    std::vector<float>sourceTerm=myUtils.computeSourceTerm(nSamples,timeStep,f0,sourceOrder);
    for(int i=0;i<nSamples;i++)
    {
      RHSTerm[i]=sourceTerm[i];
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
            vp[IDX3(i,j,k)]=1500.*1500.*timeStep2;
            phi[IDX3(i,j,k)]=0.;
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
            pnp1[IDX3_l(i,j,k)]=0.000001;
            pn[IDX3_l(i,j,k)]  =0.000001;
            pnm1[IDX3_l(i,j,k)]=0.000001;
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
                          vmax, eta);
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
                              coef0,
                              hdx_2,hdy_2,hdz_2,
                              coefx,coefy,coefz,
                              vp,phi,eta,
                              pnp1,pn,pnm1);
      // swap wavefields
      myKernel.swapWavefields(nx,ny,nz,lx,ly,lz,pnp1,pn,pnm1);
      // print infos and save wavefields
      if(itSample%50==0)
      {
        printf("result1 %f\n",pn[IDX3_l(xs,ys,zs)]);
        myFDTDUtils.write_io(nx,ny,nz,lx,ly,lz,0,nx,ny/2,ny/2,0,nz,pn,itSample);
      }

    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::time_t start_time = std::chrono::system_clock::to_time_t(start);
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "started computation at " << std::ctime(&start_time)<<std::endl
              << "finished computation at " << std::ctime(&end_time)<<std::endl
              << "elapsed time: " << elapsed_seconds.count() << "s\n";

    return (0);
}
