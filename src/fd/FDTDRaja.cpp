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
#include "Array.hpp"
#include "RAJA/RAJA.hpp"
#include "Macros.hpp"
#include "ChaiBuffer.hpp"
#include "utils.hpp"
#include "dataType.hpp"

#define MemSpace Kokkos::SharedSpace


int main( int argc, char *argv[] )
{
    constexpr int n1=208;
    constexpr int n2=208;
    constexpr int n3=208;
    constexpr float dx=10;
    constexpr float dy=10;
    constexpr float dz=10;

    constexpr int   sourceOrder=1;
    constexpr int   xs=n1/2;
    constexpr int   ys=n2/2;
    constexpr int   zs=n3/2;
    float f0=10.;
    float timeMax=2.0;

    const int ncoefs=5;
    vectorReal h_coefx;
    vectorReal h_coefy;
    vectorReal h_coefz;
    h_coefx=allocateVector<vectorReal>(ncoefs);
    h_coefy=allocateVector<vectorReal>(ncoefs);
    h_coefz=allocateVector<vectorReal>(ncoefs);

    solverUtils myUtils;
    myUtils.init_coef(dx, h_coefx);
    myUtils.init_coef(dy, h_coefy);
    myUtils.init_coef(dz, h_coefz);

    double coef0=-2.*(h_coefx[1]+h_coefx[2]+h_coefx[3]+h_coefx[4]);
    coef0+=-2.*(h_coefy[1]+h_coefy[2]+h_coefy[3]+h_coefy[4]);
    coef0+=-2.*(h_coefz[1]+h_coefz[2]+h_coefz[3]+h_coefz[4]);

    float vmax=1500;

    float timeStep=myUtils.compute_dt_sch(vmax,h_coefx,h_coefy,h_coefz);

    float timeStep2=timeStep*timeStep;
    const int nSamples=timeMax/timeStep;
    printf("timeStep=%f\n",timeStep);

    vectorReal h_RHSTerm;
    h_RHSTerm=allocateVector<vectorReal>(nSamples);
    // compute source term
    std::vector<float>sourceTerm=myUtils.computeSourceTerm(nSamples,timeStep,f0,sourceOrder);
    for(int i=0;i<nSamples;i++)
    {
      h_RHSTerm[i]=sourceTerm[i];
    }
    
    // allocate vector and arrays 
    array3DReal h_vp;
    h_vp=allocateArray3D<array3DReal>(n1,n1,n3);
    array3DReal h_pnp1;
    h_pnp1=allocateArray3D<array3DReal>(n1,n2,n3);
    array3DReal h_pn;
    h_pn=allocateArray3D<array3DReal>(n1,n2,n3);
    array3DReal h_pnm1;
    h_pnm1=allocateArray3D<array3DReal>(n1,n2,n3);

    printf("memory used for vectra and arrays %d bytes\n",(4*n1*n2*n3+nSamples+ncoefs)*4);

#pragma omp parallel for collapse(3)
    for(int i=0;i<n1;i++)
    {
       for(int j=0;j<n2;j++)
       {
          for(int k=0;k<n3;k++)
          {
            h_vp(i,j,k)=1500.*1500.;
            h_pnp1(i,j,k)=0.;
            h_pn(i,j,k)=0.;
            h_pnm1(i,j,k)=0.;
          }
       }
    }

    vectorRealView   const coefx=h_coefx.toView();
    vectorRealView   const coefy=h_coefy.toView();
    vectorRealView   const coefz=h_coefz.toView();
    vectorRealView   const RHSTerm=h_RHSTerm.toView();
    array3DRealView  const vp=h_vp.toView();
    array3DRealView  const pnp1=h_pnp1.toView();
    array3DRealView  const pn=h_pn.toView();
    array3DRealView  const pnm1=h_pnm1.toView();

    //RAJA_INDEX_VALUE_T(KIDX, int, "KIDX");
    //RAJA_INDEX_VALUE_T(JIDX, int, "JIDX");
    //RAJA_INDEX_VALUE_T(IIDX, int, "IIDX");
    constexpr int imins = xs;
    constexpr int imaxs = xs+1;
    constexpr int jmins = ys;
    constexpr int jmaxs = ys+1;
    constexpr int kmins = zs;
    constexpr int kmaxs = zs+1;
    RAJA::TypedRangeSegment<int> KRanges(kmins, kmaxs);
    RAJA::TypedRangeSegment<int> JRanges(jmins, jmaxs);
    RAJA::TypedRangeSegment<int> IRanges(imins, imaxs);
    constexpr int imini = 4;
    constexpr int imaxi = n1-4;
    constexpr int jmini = 4;
    constexpr int jmaxi = n2-4;
    constexpr int kmini = 4;
    constexpr int kmaxi = n3-4;
    RAJA::TypedRangeSegment<int> KRangei(kmini, kmaxi);
    RAJA::TypedRangeSegment<int> JRangei(jmini, jmaxi);
    RAJA::TypedRangeSegment<int> IRangei(imini, imaxi);
    constexpr int imin = 0;
    constexpr int imax = n1;
    constexpr int jmin = 0;
    constexpr int jmax = n2;
    constexpr int kmin = 0;
    constexpr int kmax = n3;
    RAJA::TypedRangeSegment<int> KRange(kmin, kmax);
    RAJA::TypedRangeSegment<int> JRange(jmin, jmax);
    RAJA::TypedRangeSegment<int> IRange(imin, imax);

    using EXEC_POL5 =
    RAJA::KernelPolicy<
      RAJA::statement::CudaKernel<
        RAJA::statement::For<2, RAJA::cuda_thread_x_loop,      // k
          RAJA::statement::For<1, RAJA::cuda_thread_y_loop,    // j
            RAJA::statement::For<0, RAJA::cuda_thread_z_loop,  // i
              RAJA::statement::Lambda<0>
            >
          >
        >
      >
    >;

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    for (int itSample=0; itSample<nSamples;itSample++)
    {
      RAJA::kernel<EXEC_POL5>( RAJA::make_tuple(IRanges, JRanges, KRanges), [=] __device__ (int i, int j, int k)
      {
        pn(i,j,k)+=vp(i,j,k)*timeStep*timeStep*RHSTerm[itSample];
      });

      RAJA::kernel<EXEC_POL5>( RAJA::make_tuple(IRangei, JRangei, KRangei), [=] __device__ (int i, int j, int k)
      {
         float lapx=(coefx[1]*(pn(i+1,j,k)+pn(i-1,j,k))
                   +coefx[2]*(pn(i+2,j,k)+pn(i-2,j,k))
                   +coefx[3]*(pn(i+3,j,k)+pn(i-3,j,k))
                   +coefx[4]*(pn(i+4,j,k)+pn(i-4,j,k)));
         float lapy=(coefy[1]*(pn(i,j+1,k)+pn(i,j-1,k))
                   +coefy[2]*(pn(i,j+2,k)+pn(i,j-2,k))
                   +coefy[3]*(pn(i,j+3,k)+pn(i,j-3,k))
                   +coefy[4]*(pn(i,j+4,k)+pn(i,j-4,k)));
         float lapz=(coefz[1]*(pn(i,j,k+1)+pn(i,j,k-1))
                   +coefz[2]*(pn(i,j,k+2)+pn(i,j,k-2))
                   +coefz[3]*(pn(i,j,k+3)+pn(i,j,k-3))
                   +coefz[4]*(pn(i,j,k+4)+pn(i,j,k-4)));
         pnp1(i,j,k)=2.*pn(i,j,k)-pnm1(i,j,k)+timeStep2*vp(i,j,k)*(coef0*pn(i,j,k)+lapx+lapy+lapz);
         //if(i==xs && j==ys && k==zs)printf("%f %f %f\n",coef0*pn(i,j,k),lapx+lapy+lapz,pn(i,j,k));
      });
      if(itSample%50==0){
        RAJA::forall<RAJA::loop_exec>(RAJA::RangeSegment(0,n1), [pnp1] ( int i)
                         {});
      printf("result 1 %f\n",pnp1(xs,ys,zs));}

      RAJA::kernel<EXEC_POL5>( RAJA::make_tuple(IRange, JRange, KRange), [=] __device__ (int i, int j, int k)
      {
        pnm1(i,j,k)=pn(i,j,k);
        pn(i,j,k)=pnp1(i,j,k);
      });
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
