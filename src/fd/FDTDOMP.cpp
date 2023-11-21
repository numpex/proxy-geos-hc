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


//inner points
void inner3D(const int nx, const int ny, const int nz,
             const int x3, const int x4,
             const int y3, const int y4,
             const int z3, const int z4,
	     const float timeStep2,
	     const float coef0,
             const vectorReal & coefx,
             const vectorReal & coefy,
             const vectorReal & coefz,
             vectorReal & vp,
             vectorReal & pnp1,
             vectorReal & pn,
             vectorReal & pnm1)
{
  #pragma omp parallel for collapse(3)
  for( int i=x3; i<x4;i++)
  {
     for( int j=y3; j<y4;j++)
     {
        for( int k=z3; k<z4;k++)
        {
          float lapx=(coefx[1]*(pn[IDX3_l(i+1,j,k)]+pn[IDX3_l(i-1,j,k)])
                    +coefx[2]*(pn[IDX3_l(i+2,j,k)]+pn[IDX3_l(i-2,j,k)])
                    +coefx[3]*(pn[IDX3_l(i+3,j,k)]+pn[IDX3_l(i-3,j,k)])
                    +coefx[4]*(pn[IDX3_l(i+4,j,k)]+pn[IDX3_l(i-4,j,k)]));
          float lapy=(coefy[1]*(pn[IDX3_l(i,j+1,k)]+pn[IDX3_l(i,j-1,k)])
                    +coefy[2]*(pn[IDX3_l(i,j+2,k)]+pn[IDX3_l(i,j-2,k)])
                    +coefy[3]*(pn[IDX3_l(i,j+3,k)]+pn[IDX3_l(i,j-3,k)])
                    +coefy[4]*(pn[IDX3_l(i,j+4,k)]+pn[IDX3_l(i,j-4,k)]));
          float lapz=(coefz[1]*(pn[IDX3_l(i,j,k+1)]+pn[IDX3_l(i,j,k-1)]])
                    +coefz[2]*(pn[IDX3_l(i,j,k+2)]+pn[IDX3_l(i,j,k-2)])
                    +coefz[3]*(pn[IDX3_l(i,j,k+3)]+pn[IDX3_l(i,j,k-3)])
                    +coefz[4]*(pn[IDX3_l(i,j,k+4)]+pn[IDX3_l(i,j,k-4)]));

          pnp1[IDX3_l(i,j,k)]=2.*pn[IDX3_l(i,j,k)]-pnm1[IDX3_l(i,j,k)]
		             +timeStep2*vp[IDX3(i,j,k)]*(coef0*pn[IDX3_l(i,j,k)]+lapx+lapy+lapz);
          //if(i==xs && j==ys && k==zs)printf("%f %f %f\n",coef0*pn(i,j,k),lapx+lapy+lapz,pn(i,j,k));
        }
     }
  }
}

void pml3d(const int nx, const int ny, const int nz,
           const int x3, const int x4, 
	   const int y3, const int y4, 
	   const int z3, const int z4,
           const int lx, const int ly, const int lz,
           const float hdx_2, const float hdy_2, const float hdz_2,
           vectorReal & vp,
           vectorReal & phi,
	   vectorReal & eta,
           vectorReal & pnp1,
           vectorReal & pn,
           vectorReal & pnm1)
{
    float coef0 = coefx[0] + coefy[0] + coefz[0];
    float lap;
    for (int i = x3; i < x4; ++i) 
    {
        for (int j = y3; j < y4; ++j) 
	{
            for (int k = z3; k < z4; ++k) 
	    {

		float lapx=(coefx[1]*(pn[IDX3_l(i+1,j,k)]+pn[IDX3_l(i-1,j,k)])
                          +coefx[2]*(pn[IDX3_l(i+2,j,k)]+pn[IDX3_l(i-2,j,k)])
                          +coefx[3]*(pn[IDX3_l(i+3,j,k)]+pn[IDX3_l(i-3,j,k)])
                          +coefx[4]*(pn[IDX3_l(i+4,j,k)]+pn[IDX3_l(i-4,j,k)]));
                float lapy=(coefy[1]*(pn[IDX3_l(i,j+1,k)]+pn[IDX3_l(i,j-1,k)])
                          +coefy[2]*(pn[IDX3_l(i,j+2,k)]+pn[IDX3_l(i,j-2,k)])
                          +coefy[3]*(pn[IDX3_l(i,j+3,k)]+pn[IDX3_l(i,j-3,k)])
                          +coefy[4]*(pn[IDX3_l(i,j+4,k)]+pn[IDX3_l(i,j-4,k)]));
                float lapz=(coefz[1]*(pn[IDX3_l(i,j,k+1)]+pn[IDX3_l(i,j,k-1)]])
                          +coefz[2]*(pn[IDX3_l(i,j,k+2)]+pn[IDX3_l(i,j,k-2)])
                          +coefz[3]*(pn[IDX3_l(i,j,k+3)]+pn[IDX3_l(i,j,k-3)])
                          +coefz[4]*(pn[IDX3_l(i,j,k+4)]+pn[IDX3_l(i,j,k-4)]));

		lap=coef0*pn(i,j,k)+lapx+lapy+lapz;

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
            }
        }
    }
}

int main( int argc, char *argv[] )
{
    constexpr int nx=20;
    constexpr int ny=20;
    constexpr int nz=20;
    constexpr int lx=4;
    constexpr int ly=4;
    constexpr int ly=4;

    
    constexpr int   sourceOrder=1;
    constexpr int   xs=nx/2;
    constexpr int   ys=ny/2;
    constexpr int   zs=nz/2;
    constexpr float f0=10.;
    constexpr float fmax=2.5*f0;
    constexpr float timeMax=2.0;

    constexpr int ncoefs=5;
    vectorReal coefx;
    vectorReal coefy;
    vectorReal coefz;
    coefx=allocateVector<vectorReal>(ncoefs);
    coefy=allocateVector<vectorReal>(ncoefs);
    coefz=allocateVector<vectorReal>(ncoefs);

    constexpr float dx=10;
    constexpr float dy=10;
    constexpr float dz=10;
    FDTDUtils myUtils;
    myUtils.init_coef(dx, coefx);
    myUtils.init_coef(dy, coefy);
    myUtils.init_coef(dz, coefz);

    double coef0=-2.*(coefx[1]+coefx[2]+coefx[3]+coefx[4]);
    coef0+=-2.*(coefy[1]+coefy[2]+coefy[3]+coefy[4]);
    coef0+=-2.*(coefz[1]+coefz[2]+coefz[3]+coefz[4]);

    constexpr float vmax=1500;

    float timeStep=myUtils.compute_dt_sch(vmax,coefx,coefy,coefz);

    float timeStep2=timeStep*timeStep;
    const int nSamples=timeMax/timeStep;
    printf("timeStep=%f\n",timeStep);

    vectorReal RHSTerm=allocateVector<vectorReal>(nSamples);
    // compute source term
    std::vector<float>sourceTerm=myUtils.computeSourceTerm(nSamples,timeStep,f0,sourceOrder);
    for(int i=0;i<nSamples;i++)
    {
      RHSTerm[i]=sourceTerm[i];
    }
    
    // init pml limits
    constexpr int ntaperx=3;
    constexpr int ntapery=3;
    constexpr int ntaperz=3;
    constexpr float hdx_2=1./(4.*dx*dx);
    constexpr float hdy_2=1./(4.*dy*dy);
    constexpr float hdz_2=1./(4.*dz*dz);
    constexpr float lambdamax=vmax/fmax;
    constexpr float ndampx=ntaperx*lambdamax/dx;
    constexpr float ndampy=ntaperx*lambdamax/dx;
    constexpr float ndampz=ntaperx*lambdamax/dx;
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
    vectorReal vp=allocateVector<VectorReal>(nx*ny*nz);
    VectorReal pnp1=allocateVector<VectorReal>((nx+2*lx)*(ny+2*ly)*(nz+2*lz));
    VectorReal pn=allocateVector<VectorReal>((nx+2*lx)*(ny+2*ly)*(nz+2*lz));
    VectorReal pnm1=allocateVector<VectorReal>((nx+2*lx)*(ny+2*ly)*(nz+2*lz));
    // PML arrays
    vectorReal phi=allocateVector<vectorReal>(nx*ny*nz);
    vectorReal eta=allocateVector<vectorReal>((nx+2)*(ny+2)*(nz+2));

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
            vp[IDX3(i,j,k)]=1500.*1500.;
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
            pnp1[IDX3_l(i,j,k)]=0.;
            pn[IDX3_l(i,j,k)]=0.;
            pnm1[IDX3_l(i,j,k)]=0.;
          }
       }
    }

    //init pml eta array
    myUtils.init_eta(int nx, int ny, int nz,
		     int ndampx, int ndampy,int ndampz,
                     int x1, int x2, int x3, int x4, int x5, int x6,
                     int y1, int y2, int y3, int y4, int y5, int y6,
                     int z1, int z2, int z3, int z4, int z5, int z6
                     float dx, float dy, float dz, float dt_sch,
                     float vmax, vectorReal eta)

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    for (int itSample=0; itSample<nSamples;itSample++)
    {
      pn(xs,ys,zs)+=vp(xs,ys,zs)*timeStep*timeStep*RHSTerm[itSample];

      //up
      pml3D(nx,ny,nz,0,nx,0,ny,z1,z2,lx,ly,lz,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
      //front
      pml3D(nx,ny,nz,0,nx,y1,y2,z3,z4,lx,ly,lz,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
      //left
      pml3D(nx,ny,nz,x1,x2,y3,y4,z3,z4,lx,ly,lz,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
      //inner points
      inner3D(nx,ny,nz,x3,x4,y3,y4,z3,z4,timeStep2,coef0,coefx,coefy,coefz,vp,pnp1,pn,pnm1);
      //right
      pml3D(nx,ny,nz,x5,x6,y3,y4,z3,z4,lx,ly,lz,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
      //back
      pml3D(nx,ny,nz,0,nx,y5,y6,z3,z4,lx,ly,lz,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
      // bottom
      pml3D(nx,ny,nz,0,nx,0,ny,z5,z6,lx,ly,lz,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);

      if(itSample%50==0){
      printf("result 1 %f\n",pnp1[IDX3_l(xs,ys,zs)]);}
      #pragma omp parallel for collapse(3)
      for( int i=0; i<nx;i++)
      {
         for( int j=0; j<ny;j++)
         {
            for( int k=0; k<nz;k++)
            {
               pnm1[IDX3_l(i,j,k)]=pn[IDX3_l(i,j,k)];
               pn[IDX3_l(i,j,k)]=pnp1[IDX3_l(i,j,k)];
            }
         }
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
