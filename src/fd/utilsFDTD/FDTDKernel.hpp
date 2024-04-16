#ifndef FDTDKERNEL_HPP
#define FDTDKERNEL_HPP

#include "dataType.hpp"
#include "FDTDMacros.hpp"

struct FDTDKernel
{
  //innerpoints
  int inner3D(const int nx, const int ny, const int nz,
         const int x3, const int x4,
         const int y3, const int y4,
         const int z3, const int z4,
         const int lx, const int ly, const int lz,
         const float coef0,
         vectorRealView const & coefx,
         vectorRealView const & coefy,
         vectorRealView const & coefz,
         vectorRealView const & vp,
         vectorRealView const & pnp1,
         vectorRealView const & pn  ,
         vectorRealView const & pnm1)const
  {
      LOOP3DHEAD
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
      LOOP3DEND
      return 0;
  }

  int pml3D(const int nx, const int ny, const int nz,
             const int x3, const int x4, 
             const int y3, const int y4, 
             const int z3, const int z4,
             const int lx, const int ly, const int lz,
             const float coef0,
             const float hdx_2, const float hdy_2, const float hdz_2,
             vectorRealView const & coefx,
             vectorRealView const & coefy,
             vectorRealView const & coefz,
             vectorRealView const & vp,
             vectorRealView const & phi,
             vectorRealView const & eta,
             vectorRealView const & pnp1,
             vectorRealView const & pn ,
             vectorRealView const & pnm1)const
  {
       LOOP3DHEAD
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
     LOOP3DEND
     return(0);
  }

  // add RHS term
  int addRHS(const int nx,const int ny,const int nz,
	     const int lx,const int ly,const int lz,
	     const int xs,const int ys,const int zs,const int itSample,
	     vectorRealView const & RHSTerm,
	     vectorRealView const & vp,
	     vectorRealView const & pn) const
  {
#ifdef USE_RAJA
     RAJANestedLoop(xs,ys,zs,xs+1,ys+1,zs+1)
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
		     vectorRealView const & pnp1,
		     vectorRealView const & pn  ,
		     vectorRealView const & pnm1) const
  {
#ifdef USE_RAJA
     RAJANestedLoop(0,0,0,nx,ny,nz)
     {
        pnm1[IDX3_l(i,j,k)]=pn[IDX3_l(i,j,k)];
        pn[IDX3_l(i,j,k)]=pnp1[IDX3_l(i,j,k)];
     });
#elif defined USE_KOKKOS
     Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0,0,0},{nz+2*lz,nx+2*lx,ny+2*ly}),KOKKOS_LAMBDA(int K,int I,int J)
     {
        int i=I-lx;
        int j=J-ly;
        int k=K-lz;
        pnm1[IDX3_l(i,j,k)]=pn[IDX3_l(i,j,k)];
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
	       pnm1[IDX3_l(i,j,k)]=pn[IDX3_l(i,j,k)];
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
                     vectorRealView const & coefx,
                     vectorRealView const & coefy,
                     vectorRealView const & coefz,
                     vectorRealView const & vp,
                     vectorRealView const & phi,
                     vectorRealView const & eta,
                     vectorRealView const & pnp1,
                     vectorRealView const & pn  ,
                     vectorRealView const & pnm1)const
  {
    //up
    pml3D(nx,ny,nz,0,nx,0,ny,z1,z2,lx,ly,lz,coef0,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
    //front
    pml3D(nx,ny,nz,0,nx,y1,y2,z3,z4,lx,ly,lz,coef0,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
    //left
    pml3D(nx,ny,nz,x1,x2,y3,y4,z3,z4,lx,ly,lz,coef0,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
    //inner points
    inner3D(nx,ny,nz,x3,x4,y3,y4,z3,z4,lx,ly,lz,coef0,coefx,coefy,coefz,vp,pnp1,pn,pnm1);
    //right
    pml3D(nx,ny,nz,x5,x6,y3,y4,z3,z4,lx,ly,lz,coef0,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
    //back
    pml3D(nx,ny,nz,0,nx,y5,y6,z3,z4,lx,ly,lz,coef0,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
    // bottom
    pml3D(nx,ny,nz,0,nx,0,ny,z5,z6,lx,ly,lz,coef0,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
    return(0);
  }

};
#endif //FDTDKERNE_HPP

