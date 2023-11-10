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

#define MemSpace Kokkos::SharedSpace


int main( int argc, char *argv[] )
{
  Kokkos::initialize(argc,argv);
  {
    const int n1=208;
    const int n2=208;
    const int n3=208;
    constexpr float dx=10;
    constexpr float dy=10;
    constexpr float dz=10;

    constexpr float invDx2=1./(dx*dx);
    constexpr float invDy2=1./(dy*dy);
    constexpr float invDz2=1./(dz*dz);
    
    int   sourceOrder=1;
    int   xs=n1/2;
    int   ys=n2/2;
    int   zs=n3/2;
    float f0=10.;
    float timeMax=2.0;

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
    vp=allocateArray3D<array3DReal>(n1,n1,n3);
    array3DReal pnp1;
    pnp1=allocateArray3D<array3DReal>(n1,n2,n3);
    array3DReal pn;
    pn=allocateArray3D<array3DReal>(n1,n2,n3);
    array3DReal pnm1;
    pnm1=allocateArray3D<array3DReal>(n1,n2,n3);

    printf("memory used for vectra and arrays %d bytes\n",(4*n1*n2*n3+nSamples+ncoefs)*4);

    Kokkos::Timer timer1;
    // initialize vp and pressure field
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>> ({0,0,0},{n1,n2,n3}) , KOKKOS_LAMBDA ( int i,int j, int k)
    {
      vp(i,j,k)=1500.*1500.;
      pnp1(i,j,k)=0.;
      pn(i,j,k)=0.;
      pnm1(i,j,k)=0.;
    });

    for (int itSample=0; itSample<nSamples;itSample++)
    {
      Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>> ({xs,xs,zs},{xs+1,ys+1,zs+1}) , KOKKOS_LAMBDA ( int i,int j, int k)
      {
        pn(i,j,k)+=vp(i,j,k)*timeStep*timeStep*RHSTerm[itSample];
      });
      Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>> ({4,4,4},{n1-4,n2-4,n3-4}) , KOKKOS_LAMBDA ( int i,int j, int k)
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
	 if(i==xs && j==ys && k==zs)printf("%f %f %f\n",coef0*pn(i,j,k),lapx+lapy+lapz,pn(i,j,k));
      });
      Kokkos::fence();
      if(itSample%50==0){
      printf("result 1 %f\n",pnp1(xs,ys,zs));}
      Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>> ({0,0,0},{n1,n2,n3}) , KOKKOS_LAMBDA ( int i,int j, int k)
      {
        pnm1(i,j,k)=pn(i,j,k);
        pn(i,j,k)=pnp1(i,j,k);
      });
    }
    Kokkos::fence();
    double time1=timer1.seconds();
    std::cout << "Elapsed Time parallel_for MDRANGE kokkos  loop : "<<time1 <<" seconds.\n\n";

  }
  Kokkos::finalize();
  return (0);
}
