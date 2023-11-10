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


int main( int argc, char *argv[] )
{
    constexpr int n1=208;
    constexpr int n2=208;
    constexpr int n3=208;
    constexpr float dx=10;
    constexpr float dy=10;
    constexpr float dz=10;

    constexpr float invDx2=1./(dx*dx);
    constexpr float invDy2=1./(dy*dy);
    constexpr float invDz2=1./(dz*dz);
    
    int   sourceOrder=1;
    constexpr int   xs=n1/2;
    constexpr int   ys=n2/2;
    constexpr int   zs=n3/2;
    constexpr float f0=10.;
    constexpr float timeMax=2.0;

    const int ncoefs=5;
    vectorReal coefx;
    vectorReal coefy;
    vectorReal coefz;
    coefx=allocateVector<vectorReal>(ncoefs);
    coefy=allocateVector<vectorReal>(ncoefs);
    coefz=allocateVector<vectorReal>(ncoefs);

    solverUtils myUtils;
    myUtils.init_coef(dx, coefx);
    myUtils.init_coef(dy, coefy);
    myUtils.init_coef(dz, coefz);

    double coef0=-2.*invDx2*(coefx[1]+coefx[2]+coefx[3]+coefx[4]);
    coef0+=-2.*invDy2*(coefy[1]+coefy[2]+coefy[3]+coefy[4]);
    coef0+=-2.*invDz2*(coefz[1]+coefz[2]+coefz[3]+coefz[4]);

    float vmax=1500;

    float timeStep=myUtils.compute_dt_sch(vmax,coefx,coefy,coefz);

    float timeStep2=timeStep*timeStep;
    const int nSamples=timeMax/timeStep;
    printf("timeStep=%f\n",timeStep);

    vectorReal RHSTerm;
    RHSTerm=allocateVector<vectorReal>(nSamples);
    // compute source term
    std::vector<float>sourceTerm=myUtils.computeSourceTerm(nSamples,timeStep,f0,sourceOrder);
    for(int i=0;i<nSamples;i++)
    {
      RHSTerm[i]=sourceTerm[i];
    }
    
    // allocate vector and arrays 
    array3DReal vp;
    vp=allocateArray3D<array3DReal>(n1,n2,n3);
    array3DReal pnp1;
    pnp1=allocateArray3D<array3DReal>(n1,n2,n3);
    array3DReal pn;
    pn=allocateArray3D<array3DReal>(n1,n2,n3);
    array3DReal pnm1;
    pnm1=allocateArray3D<array3DReal>(n1,n2,n3);

    printf("memory used for vectra and arrays %d bytes\n",(4*n1*n2*n3+nSamples+ncoefs)*4);

    // initialize vp and pressure field
#pragma omp parallel for collapse(3)
    for( int i=0; i<n1;i++)
    {
       for( int j=0; j<n2;j++)
       {
          for( int k=0; k<n3;k++)
          {
            vp(i,j,k)=1500.*1500.;
            pnp1(i,j,k)=0.;
            pn(i,j,k)=0.;
            pnm1(i,j,k)=0.;
          }
       }
    }
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    for (int itSample=0; itSample<nSamples;itSample++)
    {
      pn(xs,ys,zs)+=vp(xs,ys,zs)*timeStep*timeStep*RHSTerm[itSample];
#pragma omp parallel for collapse(3)
      for( int i=4; i<n1-4;i++)
      {
         for( int j=4; j<n2-4;j++)
         {
            for( int k=4; k<n3-4;k++)
            {
              float lapx=(coefx[1]*(pn(i+1,j,k)+pn(i-1,j,k))
                        +coefx[2]*(pn(i+2,j,k)+pn(i-2,j,k))
                        +coefx[3]*(pn(i+3,j,k)+pn(i-3,j,k))
                        +coefx[4]*(pn(i+4,j,k)+pn(i-4,j,k)))*invDx2;
              float lapy=(coefy[1]*(pn(i,j+1,k)+pn(i,j-1,k))
                        +coefy[2]*(pn(i,j+2,k)+pn(i,j-2,k))
                        +coefy[3]*(pn(i,j+3,k)+pn(i,j-3,k)) 
                        +coefy[4]*(pn(i,j+4,k)+pn(i,j-4,k)))*invDy2;
              float lapz=(coefz[1]*(pn(i,j,k+1)+pn(i,j,k-1))
                        +coefz[2]*(pn(i,j,k+2)+pn(i,j,k-2))
                        +coefz[3]*(pn(i,j,k+3)+pn(i,j,k-3))
                        +coefz[4]*(pn(i,j,k+4)+pn(i,j,k-4)))*invDz2;
              pnp1(i,j,k)=2.*pn(i,j,k)-pnm1(i,j,k)+timeStep2*vp(i,j,k)*(coef0*pn(i,j,k)+lapx+lapy+lapz);
              //if(i==xs && j==ys && k==zs)printf("%f %f %f\n",coef0*pn(i,j,k),lapx+lapy+lapz,pn(i,j,k));
            }
         }
      }
      if(itSample%50==0){
      printf("result 1 %f\n",pnp1(xs,ys,zs));}
#pragma omp parallel for collapse(3)
      for( int i=0; i<n1;i++)
      {
         for( int j=0; j<n2;j++)
         {
            for( int k=0; k<n3;k++)
            {
               pnm1(i,j,k)=pn(i,j,k);
               pn(i,j,k)=pnp1(i,j,k);
            }
         }
      }
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "finished computation at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s\n";

    return (0);
}
