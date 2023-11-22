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
#include <Kokkos_ScatterView.hpp>
#include <omp.h>
#include <vector>
#include "utils.hpp"
#include "dataType.hpp"
#include "FDTDutils.hpp"


//inner points
void inner3D(const int nx, const int ny, const int nz,
        const int x3, const int x4,
        const int y3, const int y4,
        const int z3, const int z4,
        const int lx, const int ly, const int lz,
        const float coef0,
        vectorReal const & coefx,
        vectorReal const & coefy,
        vectorReal const & coefz,
        vectorReal const & vp,
        vectorReal const & pnp1,
        vectorReal const & pn,
        vectorReal const & pnm1)
{
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({x3,y3,z3},{x4,y4,z4}),KOKKOS_LAMBDA(int i,int j,int k)
    {
     float lapx=(coefx[1]*(pn[IDX3_l(i+1,j,k)]+pn[IDX3_l(i-1,j,k)])
                +coefx[2]*(pn[IDX3_l(i+2,j,k)]+pn[IDX3_l(i-2,j,k)])
                +coefx[3]*(pn[IDX3_l(i+3,j,k)]+pn[IDX3_l(i-3,j,k)])
                +coefx[4]*(pn[IDX3_l(i+4,j,k)]+pn[IDX3_l(i-4,j,k)]));
     float lapy=(coefy[1]*(pn[IDX3_l(i,j+1,k)]+pn[IDX3_l(i,j-1,k)])
                +coefy[2]*(pn[IDX3_l(i,j+2,k)]+pn[IDX3_l(i,j-2,k)])
                +coefy[3]*(pn[IDX3_l(i,j+3,k)]+pn[IDX3_l(i,j-3,k)])
                +coefy[4]*(pn[IDX3_l(i,j+4,k)]+pn[IDX3_l(i,j-4,k)]));
     float lapz=(coefz[1]*(pn[IDX3_l(i,j,k+1)]+pn[IDX3_l(i,j,k-1)])
                +coefz[2]*(pn[IDX3_l(i,j,k+2)]+pn[IDX3_l(i,j,k-2)])
                +coefz[3]*(pn[IDX3_l(i,j,k+3)]+pn[IDX3_l(i,j,k-3)])
                +coefz[4]*(pn[IDX3_l(i,j,k+4)]+pn[IDX3_l(i,j,k-4)]));
     pnp1[IDX3_l(i,j,k)]=2.*pn[IDX3_l(i,j,k)]-pnm1[IDX3_l(i,j,k)]
                        +vp[IDX3(i,j,k)]*(coef0*pn[IDX3_l(i,j,k)]+lapx+lapy+lapz);
     //if(i==nx/2 && j==ny/2 && k==nz/2)printf("%f %f %f\n",coef0*pn[IDX3_l(i,j,k)],
     //			                          lapx+lapy+lapz,pn[IDX3_l(i,j,k)]);
    });
}

void pml3D(const int nx, const int ny, const int nz,
            const int x3, const int x4, 
            const int y3, const int y4, 
            const int z3, const int z4,
            const int lx, const int ly, const int lz,
            const float hdx_2, const float hdy_2, const float hdz_2,
            vectorReal const & coefx,
            vectorReal const & coefy,
            vectorReal const & coefz,
            vectorReal const & vp,
            vectorReal const & phi,
            vectorReal const & eta,
            vectorReal const & pnp1,
            vectorReal const & pn,
            vectorReal const & pnm1)
{
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({x3,y3,z3},{x4,y4,z4}),KOKKOS_LAMBDA(int i,int j,int k)
    {
      float coef0 = coefx[0] + coefy[0] + coefz[0];
      float lap;
      float lapx=(coefx[1]*(pn[IDX3_l(i+1,j,k)]+pn[IDX3_l(i-1,j,k)])
                 +coefx[2]*(pn[IDX3_l(i+2,j,k)]+pn[IDX3_l(i-2,j,k)])
                 +coefx[3]*(pn[IDX3_l(i+3,j,k)]+pn[IDX3_l(i-3,j,k)])
                 +coefx[4]*(pn[IDX3_l(i+4,j,k)]+pn[IDX3_l(i-4,j,k)]));
      float lapy=(coefy[1]*(pn[IDX3_l(i,j+1,k)]+pn[IDX3_l(i,j-1,k)])
                 +coefy[2]*(pn[IDX3_l(i,j+2,k)]+pn[IDX3_l(i,j-2,k)])
                 +coefy[3]*(pn[IDX3_l(i,j+3,k)]+pn[IDX3_l(i,j-3,k)])
                 +coefy[4]*(pn[IDX3_l(i,j+4,k)]+pn[IDX3_l(i,j-4,k)]));
      float lapz=(coefz[1]*(pn[IDX3_l(i,j,k+1)]+pn[IDX3_l(i,j,k-1)])
                 +coefz[2]*(pn[IDX3_l(i,j,k+2)]+pn[IDX3_l(i,j,k-2)])
                 +coefz[3]*(pn[IDX3_l(i,j,k+3)]+pn[IDX3_l(i,j,k-3)])
                 +coefz[4]*(pn[IDX3_l(i,j,k+4)]+pn[IDX3_l(i,j,k-4)]));

      lap=coef0*pn[IDX3_l(i,j,k)]+lapx+lapy+lapz;

      pnp1[IDX3_l(i,j,k)]=((2.-eta[IDX3_eta1(i,j,k)]*eta[IDX3_eta1(i,j,k)]
                 +2.*eta[IDX3_eta1(i,j,k)])*pn[IDX3_l(i,j,k)]
                 +vp[IDX3(i,j,k)]*(lap+phi[IDX3(i,j,k)]))/(1.+2.*eta[IDX3_eta1(i,j,k)])
                 -pnm1[IDX3_l(i,j,k)];

      phi[IDX3(i,j,k)]=(phi[IDX3(i,j,k)]-((eta[IDX3_eta1(i+1,j,k)]-eta[IDX3_eta1(i-1,j,k)])
                 *(pn[IDX3_l(i+1,j,k)]-pn[IDX3_l(i-1,j,k)])*hdx_2
                 +(eta[IDX3_eta1(i,j+1,k)]-eta[IDX3_eta1(i,j-1,k)])
                 *(pn[IDX3_l(i,j+1,k)]-pn[IDX3_l(i,j-1,k)])*hdy_2
                 +(eta[IDX3_eta1(i,j,k+1)]-eta[IDX3_eta1(i,j,k-1)])
                 *(pn[IDX3_l(i,j,k+1)]-pn[IDX3_l(i,j,k-1)])*hdz_2))
                 /(1.+eta[IDX3_eta1(i,j,k)]);
    });
}

int main( int argc, char *argv[] )
{
  Kokkos::initialize(argc,argv);
  {
    constexpr int nx=150;
    constexpr int ny=150;
    constexpr int nz=150;
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
    constexpr float f0=5.;
    constexpr float fmax=2.5*f0;
    constexpr float timeMax=0.8;

    constexpr int ncoefs=5;
    constexpr float vmax=1500;
    
    // imports utils
    solverUtils myUtils;
    FDTDUtils myFDTDUtils;

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
      Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({xs,xs,zs},{xs+1,ys+1,zs+1}),KOKKOS_LAMBDA(int i,int j,int k)
      {
        pn[IDX3_l(xs,ys,zs)]+=vp[IDX3(xs,ys,zs)]*RHSTerm[itSample];
      });
      //up
      pml3D(nx,ny,nz,0,nx,0,ny,z1,z2,lx,ly,lz,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
      //front
      pml3D(nx,ny,nz,0,nx,y1,y2,z3,z4,lx,ly,lz,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
        Kokkos::fence();
      //left
      pml3D(nx,ny,nz,x1,x2,y3,y4,z3,z4,lx,ly,lz,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
      //inner points
      inner3D(nx,ny,nz,x3,x4,y3,y4,z3,z4,lx,ly,lz,coef0,coefx,coefy,coefz,vp,pnp1,pn,pnm1);
      //right
      pml3D(nx,ny,nz,x5,x6,y3,y4,z3,z4,lx,ly,lz,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
      //back
      pml3D(nx,ny,nz,0,nx,y5,y6,z3,z4,lx,ly,lz,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
      // bottom
      pml3D(nx,ny,nz,0,nx,0,ny,z5,z6,lx,ly,lz,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);

      if(itSample%50==0)
      {
        Kokkos::fence();
	printf("result 1 %f\n",pnp1[IDX3_l(xs,ys,zs)]);
      }
      Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0,0,0},{nx,ny,nz}),KOKKOS_LAMBDA(int i,int j,int k)
      {
         pnm1[IDX3_l(i,j,k)]=pn[IDX3_l(i,j,k)];
         pn[IDX3_l(i,j,k)]=pnp1[IDX3_l(i,j,k)];
      });
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
