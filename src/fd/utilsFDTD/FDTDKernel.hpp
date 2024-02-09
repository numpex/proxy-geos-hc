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


#ifdef USE_RAJA
// implementation of sqrt via binary search
// copied from https://stackoverflow.com/questions/8622256/in-c11-is-sqrt-defined-as-constexpr
constexpr size_t sqrt_helper(size_t n, size_t lo, size_t hi)
{
  return (lo == hi)
           ? lo // search complete
           : ((n / ((lo + hi + 1) / 2) < ((lo + hi + 1) / 2))
                ? sqrt_helper(n, lo, ((lo + hi + 1) / 2)-1) // search lower half
                : sqrt_helper(n, ((lo + hi + 1) / 2), hi)); // search upper half
}
// constexpr integer sqrt
constexpr size_t sqrt(size_t n)
{               
  return sqrt_helper(n, 0, n/2 + 1); 
}
// implementation of lesser_of_squarest_factor_pair via linear search
constexpr size_t lesser_of_squarest_factor_pair_helper(size_t n, size_t guess)
{
  return ((n / guess) * guess == n)
           ? guess // search complete, guess is a factor
           : lesser_of_squarest_factor_pair_helper(n, guess - 1); // continue searching
}
// constexpr return the lesser of the most square pair of factors of n
// ex. 12 has pairs of factors (1, 12) (2, 6) *(3, 4)* and returns 3
constexpr size_t lesser_of_squarest_factor_pair(size_t n)
{
  return (n == 0)
      ? 0 // return 0 in the 0 case
      : lesser_of_squarest_factor_pair_helper(n, sqrt(n));
}
// constexpr return the greater of the most square pair of factors of n
// ex. 12 has pairs of factors (1, 12) (2, 6) *(3, 4)* and returns 4
constexpr size_t greater_of_squarest_factor_pair(size_t n)
{
  return (n == 0)
      ? 0 // return 0 in the 0 case
      : n / lesser_of_squarest_factor_pair_helper(n, sqrt(n));
}

#define block_size (256)
#define x_block_sz (32)
#define y_block_sz (greater_of_squarest_factor_pair(block_size/x_block_sz))
#define z_block_sz (lesser_of_squarest_factor_pair(block_size/x_block_sz))
//const int x_block_sz=32;
//const int y_block_sz=4;
//const int z_block_sz=2;
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
         vectorRealView const & pn)const
  #elif defined USE_KOKKOS
         vectorReal const & coefx,
         vectorReal const & coefy,
         vectorReal const & coefz,
         vectorReal const & vp,
         vectorReal const & pnp1,
         vectorReal const & pn)const
  #else
         vectorReal  & coefx,
         vectorReal  & coefy,
         vectorReal  & coefz,
         vectorReal  & vp,
         vectorReal  & pnp1,
         vectorReal  & pn)
#endif
  {
#ifdef USE_RAJA
     RAJA::TypedRangeSegment<int> KRange(z3, z4);
     RAJA::TypedRangeSegment<int> JRange(y3, y4);
     RAJA::TypedRangeSegment<int> IRange(x3, x4);


     using EXEC_POL =
     RAJA::KernelPolicy<
        RAJA::statement::CudaKernelFixedAsync<x_block_sz*y_block_sz*z_block_sz,
          RAJA::statement::For<0, RAJA::cuda_global_size_z_direct<x_block_sz>,     //z
            RAJA::statement::For<1, RAJA::cuda_global_size_y_direct<y_block_sz>,   //g
              RAJA::statement::For<2, RAJA::cuda_global_size_x_direct<z_block_sz>,
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
      pnp1[IDX3_l(i,j,k)]=2.*pn[IDX3_l(i,j,k)]-pnp1[IDX3_l(i,j,k)]
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
             vectorRealView const & pn)const
#elif defined USE_KOKKOS
             vectorReal const & coefx,
             vectorReal const & coefy,
             vectorReal const & coefz,
             vectorReal const & vp,
             vectorReal const & phi,
             vectorReal const & eta,
             vectorReal const & pnp1,
             vectorReal const & pn)const
#else
             vectorReal  & coefx,
             vectorReal  & coefy,
             vectorReal  & coefz,
             vectorReal  & vp,
             vectorReal  & phi,
             vectorReal  & eta,
             vectorReal  & pnp1,
             vectorReal  & pn)
#endif
  {
#ifdef USE_RAJA
     RAJA::TypedRangeSegment<int> KRange(z3, z4);
     RAJA::TypedRangeSegment<int> JRange(y3, y4);
     RAJA::TypedRangeSegment<int> IRange(x3, x4);

     using EXEC_POL =
     RAJA::KernelPolicy<
        RAJA::statement::CudaKernelFixedAsync<x_block_sz*y_block_sz*z_block_sz,
          RAJA::statement::For<0, RAJA::cuda_global_size_z_direct<x_block_sz>,     //z
            RAJA::statement::For<1, RAJA::cuda_global_size_y_direct<y_block_sz>,   //g
              RAJA::statement::For<2, RAJA::cuda_global_size_x_direct<z_block_sz>,
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
                  -pnp1[IDX3_l(i,j,k)]
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

  // add RHS term
  int addRHS(const int nx,const int ny,const int nz,
	     const int lx,const int ly,const int lz,
	     const int xs,const int ys,const int zs,const int itSample,
#ifdef USE_RAJA
	     vectorRealView const & RHSTerm,
	     vectorRealView const & vp,
	     vectorRealView const & pn) const
#elif defined USE_KOKKOS
	     vectorReal const & RHSTerm,
	     vectorReal const & vp,
	     vectorReal const & pn) const
#else
	     vectorReal & RHSTerm,
	     vectorReal & vp,
	     vectorReal & pn)
#endif
  {
#ifdef USE_RAJA
     // define policy for source term
     RAJA::TypedRangeSegment<int> KRanges(xs, xs+1);
     RAJA::TypedRangeSegment<int> JRanges(ys, ys+1);
     RAJA::TypedRangeSegment<int> IRanges(zs, zs+1);

     using EXEC_POL =
     RAJA::KernelPolicy<
        RAJA::statement::CudaKernelFixedAsync<x_block_sz*y_block_sz*z_block_sz,
          RAJA::statement::For<0, RAJA::cuda_global_size_z_direct<x_block_sz>,     //z
            RAJA::statement::For<1, RAJA::cuda_global_size_y_direct<y_block_sz>,   //g
              RAJA::statement::For<2, RAJA::cuda_global_size_x_direct<z_block_sz>,
               RAJA::statement::Lambda<0,RAJA::Segs<0,1,2>>
             >
           >
         >
       >
     >;

     RAJA::kernel<EXEC_POL>( RAJA::make_tuple(IRanges, JRanges, KRanges), [=] __device__ (int i, int j, int k)
     {
       pn[IDX3_l(i,j,k)]+=vp[IDX3(i,j,k)]*RHSTerm[itSample];
     });
#elif defined USE_KOKKOS
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({xs,xs,zs},{xs+1,ys+1,zs+1}),KOKKOS_LAMBDA(int i,int j,int k)
    {
    pn[IDX3_l(i,j,k)]+=vp[IDX3(i,j,k)]*RHSTerm[itSample];
    });
#else
    pn[IDX3_l(xs,ys,zs)]+=vp[IDX3(xs,ys,zs)]*RHSTerm[itSample];
#endif
    return(0);
  }

  // swap wavefields
  int swapWavefields(const int nx,const int ny,const int nz,
		     const int lx,const int ly,const int lz,
#ifdef USE_RAJA
		     vectorRealView const & pnp1,
		     vectorRealView const & pn) const
#elif defined USE_KOKKOS
		     vectorReal const & pnp1,
		     vectorReal const & pn) const
#else
		     vectorReal & pnp1,
		     vectorReal & pn)
#endif
  {
#ifdef USE_RAJA
     // define policy for swapping
     RAJA::TypedRangeSegment<int> KRange(0, nx);
     RAJA::TypedRangeSegment<int> JRange(0, ny);
     RAJA::TypedRangeSegment<int> IRange(0, ny);
 
     using EXEC_POL =
     RAJA::KernelPolicy<
        RAJA::statement::CudaKernelFixedAsync<x_block_sz*y_block_sz*z_block_sz,
          RAJA::statement::For<0, RAJA::cuda_global_size_z_direct<x_block_sz>,     //z
            RAJA::statement::For<1, RAJA::cuda_global_size_y_direct<y_block_sz>,   //g
              RAJA::statement::For<2, RAJA::cuda_global_size_x_direct<z_block_sz>,
               RAJA::statement::Lambda<0,RAJA::Segs<0,1,2>>
             >
           >
         >
       >
     >;

     RAJA::kernel<EXEC_POL>( RAJA::make_tuple(IRange, JRange, KRange), [=] __device__ (int i, int j, int k)
     {
         pn[IDX3_l(i,j,k)]=pnp1[IDX3_l(i,j,k)];
     });
#elif defined USE_KOKKOS
     Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0,0,0},{nz+2*lz,nx+2*lx,ny+2*ly}),KOKKOS_LAMBDA(int K,int I,int J)
     {
        int i=I-lx;
        int j=J-ly;
        int k=K-lz;
        pn[IDX3_l(i,j,k)]=pnp1[IDX3_l(i,j,k)];
     });
#else
#ifdef USE_OMP
      #pragma omp parallel for collapse(3)
#endif
      for( int i=0; i<nx;i++)
      {
         for( int j=0; j<ny;j++)
         {
            for( int k=0; k<nz;k++)
            {
               pn[IDX3_l(i,j,k)]=pnp1[IDX3_l(i,j,k)];
            }
         }
      }
#endif
      return(0);
  }

  // compute one step
  int computeOneStep(const int nx, const int ny, const int nz,
                     const int lx, const int ly, const int lz,
	             const int x1, const int x2, const int x3, 
	             const int x4, const int x5, const int x6,
	             const int y1, const int y2, const int y3, 
	             const int y4, const int y5, const int y6,
	             const int z1, const int z2, const int z3, 
	             const int z4, const int z5, const int z6,
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
                     vectorRealView const & pn)const
#elif defined USE_KOKKOS
                     vectorReal const & coefx,
                     vectorReal const & coefy,
                     vectorReal const & coefz,
                     vectorReal const & vp,
                     vectorReal const & phi,
                     vectorReal const & eta,
                     vectorReal const & pnp1,
                     vectorReal const & pn)const
#else
                     vectorReal  & coefx,
                     vectorReal  & coefy,
                     vectorReal  & coefz,
                     vectorReal  & vp,
                     vectorReal  & phi,
                     vectorReal  & eta,
                     vectorReal  & pnp1,
                     vectorReal  & pn)
#endif
  {
    //up
    pml3D(nx,ny,nz,0,nx,0,ny,z1,z2,lx,ly,lz,coef0,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn);
    //front
    pml3D(nx,ny,nz,0,nx,y1,y2,z3,z4,lx,ly,lz,coef0,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn);
    //left
    pml3D(nx,ny,nz,x1,x2,y3,y4,z3,z4,lx,ly,lz,coef0,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn);
    //inner points
    inner3D(nx,ny,nz,x3,x4,y3,y4,z3,z4,lx,ly,lz,coef0,coefx,coefy,coefz,vp,pnp1,pn);
    //right
    pml3D(nx,ny,nz,x5,x6,y3,y4,z3,z4,lx,ly,lz,coef0,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn);
    //back
    pml3D(nx,ny,nz,0,nx,y5,y6,z3,z4,lx,ly,lz,coef0,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn);
    // bottom
    pml3D(nx,ny,nz,0,nx,0,ny,z5,z6,lx,ly,lz,coef0,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn);
    return(0);
  }

};
#endif //FDTDKERNE_HPP
