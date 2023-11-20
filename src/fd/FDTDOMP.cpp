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

#define POW2(x) ((x)*(x))
#define IDX3(i,j,k) (n3*n2*(i) + n3*(j) + (k))
#define IDX3_l(i,j,k) ((n3+2*lz)*(n2+2*ly)*((i)+lx) + (n3+2*lz)*((j)+ly) + ((k)+lz))
#define IDX3_eta1(i,j,k) ((n3+2)*(n2+2)*((i)+1) + (n3+2)*((j)+1) + ((k)+1))
#define IDX3_eta0(i,j,k) ((n3+2)*(n2+2)*(i) + (n3+2)*(j) + (k))

//inner points
void inner3D(const int n1, const int n2, const int n3,
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

void pml3d(const int n1, const int n2, const int n3,
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
    constexpr int n1=20;
    constexpr int n2=20;
    constexpr int n3=20;
    constexpr int lx=4;
    constexpr int ly=4;
    constexpr int ly=4;

    constexpr int x1=0;
    constexpr int x2=lx-1;
    constexpr int x3=lx;
    constexpr int x4=n1-lx;
    constexpr int x5=n1=lx+1;
    constexpr int x6=n1;
    constexpr int y1=0;
    constexpr int y2=ly-1;
    constexpr int y3=ly;
    constexpr int y4=n2-ly;
    constexpr int y5=n2=ly+1;
    constexpr int y6=n2;
    constexpr int z1=0;
    constexpr int z2=lz-1;
    constexpr int z3=lz;
    constexpr int z4=n3-lz;
    constexpr int z5=n3=lz+1;
    constexpr int z6=n3;
    
    constexpr int   sourceOrder=1;
    constexpr int   xs=n1/2;
    constexpr int   ys=n2/2;
    constexpr int   zs=n3/2;
    constexpr float f0=10.;
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
    solverUtils myUtils;
    myUtils.init_coef(dx, coefx);
    myUtils.init_coef(dy, coefy);
    myUtils.init_coef(dz, coefz);

    double coef0=-2.*(coefx[1]+coefx[2]+coefx[3]+coefx[4]);
    coef0+=-2.*(coefy[1]+coefy[2]+coefy[3]+coefy[4]);
    coef0+=-2.*(coefz[1]+coefz[2]+coefz[3]+coefz[4]);

    float vmax=1500;

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
    
    // allocate vector and arrays 
    vectorReal vp=allocateVector<VectorReal>(n1*n2*n3);
    VectorReal pnp1=allocateVector<VectorReal>((n1+2*lx)*(n2+2*ly)*(n3+2*lz));
    VectorReal pn=allocateVector<VectorReal>((n1+2*lx)*(n2+2*ly)*(n3+2*lz));
    VectorReal pnm1=allocateVector<VectorReal>((n1+2*lx)*(n2+2*ly)*(n3+2*lz));
    // PML arrays
    vectorReal phi=allocateVector<vectorReal>(n1*n2*n3);
    vectorReal eta=allocateVector<vectorReal>((n1+2)*(n2+2)*(n3+2));

    printf("memory used for vectra and arrays %d bytes\n",
           (n1*n2*n3+3*(n1+2*lx)*(n2+2*ly)*(n3+2*lz)+nSamples+ncoefs)*4);

    // initialize vp and pressure field
    #pragma omp parallel for collapse(3)
    for( int i=0; i<n1;i++)
    {
       for( int j=0; j<n2;j++)
       {
          for( int k=0; k<n3;k++)
          {
            vp[IDX3(i,j,k)]=1500.*1500.;
          }
       }
    }
    #pragma omp parallel for collapse(3)
    for( int i=-lx; i<n1+lx;i++)
    {
       for( int j=-ly; j<n2+ly;j++)
       {
          for( int k=-lz; k<n3+lz;k++)
          {
            pnp1[IDX3_l(i,j,k)]=0.;
            pn[IDX3_l(i,j,k)]=0.;
            pnm1[IDX3_l(i,j,k)]=0.;
          }
       }
    }

    // init pml
    constexpr float hdx_2=1./(4.*dx*dx);
    constexpr float hdy_2=1./(4.*dy*dy);
    constexpr float hdz_2=1./(4.*dz*dz);

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    for (int itSample=0; itSample<nSamples;itSample++)
    {
      pn(xs,ys,zs)+=vp(xs,ys,zs)*timeStep*timeStep*RHSTerm[itSample];

      //up
      pml3D(n1,n2,n3,0,n1,0,n2,z1,z2,lx,ly,lz,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
      //front
      pml3D(n1,n2,n3,0,n1,y1,y2,z3,z4,lx,ly,lz,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
      //left
      pml3D(n1,n2,n3,x1,x2,y3,y4,z3,z4,lx,ly,lz,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
      //inner points
      inner3D(n1,n2,n3,x3,x4,y3,y4,z3,z4,timeStep2,coef0,coefx,coefy,coefz,vp,pnp1,pn,pnm1);
      //right
      pml3D(n1,n2,n3,x5,x6,y3,y4,z3,z4,lx,ly,lz,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
      //back
      pml3D(n1,n2,n3,0,n1,y5,y6,z3,z4,lx,ly,lz,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
      // bottom
      pml3D(n1,n2,n3,0,n1,0,n2,z5,z6,lx,ly,lz,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);

      if(itSample%50==0){
      printf("result 1 %f\n",pnp1[IDX3_l(xs,ys,zs)]);}
      #pragma omp parallel for collapse(3)
      for( int i=0; i<n1;i++)
      {
         for( int j=0; j<n2;j++)
         {
            for( int k=0; k<n3;k++)
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
