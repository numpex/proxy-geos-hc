#ifndef FDTDKERNEL_HPP
#define FDTDKERNEL_HPP

#include "dataType.hpp"
#ifdef USE_RAJA
#include"Array.hpp"
#include"RAJA/RAJA.hpp"
#include"Macros.hpp"
#include"ChaiBuffer.hpp"
#elif defined USE_KOKKOS
#include <Kokkos_Macros.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>
#endif


struct FDTDKernel
{
  //innerpoints
  int inner3D(const int nx, const int ny, const int nz,
         const int x3, const int x4,
         const int y3, const int y4,
         const int z3, const int z4,
         const int lx, const int ly, const int lz,
         const float coef0,
  #ifdef USE_RAJA
         vectorRealView const & coefx,
         vectorRealView const & coefy,
         vectorRealView const & coefz,
         vectorRealView const & vp,
         vectorRealView const & pnp1,
         vectorRealView const & pn,
         vectorRealView const & pnm1)const
  #elif defined USE_KOKKOS
         vectorReal const & coefx,
         vectorReal const & coefy,
         vectorReal const & coefz,
         vectorReal const & vp,
         vectorReal const & pnp1,
         vectorReal const & pn,
         vectorReal const & pnm1)const
  #else
         vectorReal  & coefx,
         vectorReal  & coefy,
         vectorReal  & coefz,
         vectorReal  & vp,
         vectorReal  & pnp1,
         vectorReal  & pn,
         vectorReal  & pnm1)
#endif
  {
#ifdef USE_RAJA
     RAJA::TypedRangeSegment<int> KRange(z3, z4);
     RAJA::TypedRangeSegment<int> JRange(y3, y4);
     RAJA::TypedRangeSegment<int> IRange(x3, x4);

     using EXEC_POL =
     RAJA::KernelPolicy<
       RAJA::statement::CudaKernel<
         RAJA::statement::For<2, RAJA::cuda_thread_x_loop,      // k
           RAJA::statement::For<1, RAJA::cuda_thread_y_loop,    // j
             RAJA::statement::For<0, RAJA::cuda_thread_z_loop,  // i
               RAJA::statement::Lambda<0,RAJA::Segs<0,1,2>>
             >
           >
         >
       >
     >;
     RAJA::kernel<EXEC_POL>( RAJA::make_tuple(IRange, JRange, KRange), [=] __device__ (int i, int j, int k) 
     {
#elif defined USE_KOKKOS
     Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({z3,x3,y3},{z4,x4,y4}),KOKKOS_LAMBDA(int k,int i,int j)
     {
#else
#ifdef USE_OMP
      #pragma omp parallel for collapse(3)
#endif
      for (int i = x3; i < x4; ++i)
      {
          for (int j = y3; j < y4; ++j)
          {
              for (int k = z3; k < z4; ++k)
              {

#endif
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
#ifdef USE_RAJA
     });
#elif defined USE_KOKKOS
     });
#else
             }
          }
      }
#endif
     return (0);
  }

  int pml3D(const int nx, const int ny, const int nz,
             const int x3, const int x4, 
             const int y3, const int y4, 
             const int z3, const int z4,
             const int lx, const int ly, const int lz,
             const float coef0,
             const float hdx_2, const float hdy_2, const float hdz_2,
#ifdef USE_RAJA
             vectorRealView const & coefx,
             vectorRealView const & coefy,
             vectorRealView const & coefz,
             vectorRealView const & vp,
             vectorRealView const & phi,
             vectorRealView const & eta,
             vectorRealView const & pnp1,
             vectorRealView const & pn,
             vectorRealView const & pnm1)const
#elif defined USE_KOKKOS
             vectorReal const & coefx,
             vectorReal const & coefy,
             vectorReal const & coefz,
             vectorReal const & vp,
             vectorReal const & phi,
             vectorReal const & eta,
             vectorReal const & pnp1,
             vectorReal const & pn,
             vectorReal const & pnm1)const
#else
             vectorReal  & coefx,
             vectorReal  & coefy,
             vectorReal  & coefz,
             vectorReal  & vp,
             vectorReal  & phi,
             vectorReal  & eta,
             vectorReal  & pnp1,
             vectorReal  & pn,
             vectorReal  & pnm1)
#endif
  {
#ifdef USE_RAJA
     RAJA::TypedRangeSegment<int> KRange(z3, z4);
     RAJA::TypedRangeSegment<int> JRange(y3, y4);
     RAJA::TypedRangeSegment<int> IRange(x3, x4);

     using EXEC_POL =
     RAJA::KernelPolicy<
       RAJA::statement::CudaKernel<
         RAJA::statement::For<2, RAJA::cuda_thread_x_loop,      // k
           RAJA::statement::For<1, RAJA::cuda_thread_y_loop,    // j
             RAJA::statement::For<0, RAJA::cuda_thread_z_loop,  // i
               RAJA::statement::Lambda<0,RAJA::Segs<0,1,2>>
             >
           >
         >
       >
     >;
     RAJA::kernel<EXEC_POL>( RAJA::make_tuple(IRange, JRange, KRange), [=] __device__ (int i, int j, int k) 
     {
#elif defined USE_KOKKOS
     Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({z3,x3,y3},{z4,x4,y4}),KOKKOS_LAMBDA(int k,int i,int j)
     {
#else
#ifdef USE_OMP
      #pragma omp parallel for collapse(3)
#endif
      for (int i = x3; i < x4; ++i)
      {
          for (int j = y3; j < y4; ++j)
          {
              for (int k = z3; k < z4; ++k)
              {
#endif
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

       float lap=coef0*pn[IDX3_l(i,j,k)]+lapx+lapy+lapz;

       pnp1[IDX3_l(i,j,k)]=((2.-eta[IDX3_eta1(i,j,k)]*eta[IDX3_eta1(i,j,k)]
                  +2.*eta[IDX3_eta1(i,j,k)])*pn[IDX3_l(i,j,k)]
                  -pnm1[IDX3_l(i,j,k)]
                  +vp[IDX3(i,j,k)]*(lap+phi[IDX3(i,j,k)]))/(1.+2.*eta[IDX3_eta1(i,j,k)]);

       phi[IDX3(i,j,k)]=(phi[IDX3(i,j,k)]-((eta[IDX3_eta1(i+1,j,k)]-eta[IDX3_eta1(i-1,j,k)])
                  *(pn[IDX3_l(i+1,j,k)]-pn[IDX3_l(i-1,j,k)])*hdx_2
                  +(eta[IDX3_eta1(i,j+1,k)]-eta[IDX3_eta1(i,j-1,k)])
                  *(pn[IDX3_l(i,j+1,k)]-pn[IDX3_l(i,j-1,k)])*hdy_2
                  +(eta[IDX3_eta1(i,j,k+1)]-eta[IDX3_eta1(i,j,k-1)])
                  *(pn[IDX3_l(i,j,k+1)]-pn[IDX3_l(i,j,k-1)])*hdz_2))
                  /(1.+eta[IDX3_eta1(i,j,k)]);
#ifdef USE_RAJA
     });
#elif defined USE_KOKKOS
     });
#else
             }
          }
      }
#endif
     return(0);
  }
};
#endif //FDTDKERNE_HPP
