#ifndef FDTDKERNEL_HPP
#define FDTDKERNEL_HPP

struct FDTDKernel
{

  //innerpoints
  int inner3D( FDTDGRIDS &myGrids,
         const int x3, const int x4,
         const int y3, const int y4,
         const int z3, const int z4,
         const float coef0,
         vectorRealView const &coefx, vectorRealView const &coefy, vectorRealView const &coefz, vectorRealView const &vp,
         vectorRealView const &pnp1, vectorRealView const &pn, vectorRealView const &pnm1) const
  {
#ifdef USE_OMP
      #pragma omp parallel for collapse(3)
#endif
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
      pnp1[IDX3_l(i,j,k)]=2.*pn[IDX3_l(i,j,k)]-pnm1[IDX3_l(i,j,k)]+vp[IDX3(i,j,k)]*(coef0*pn[IDX3_l(i,j,k)]+lapx+lapy+lapz);
      LOOP3DEND
      return 0;
  }

  int pml3D( FDTDGRIDS &myGrids,
             const int x3, const int x4, 
             const int y3, const int y4, 
             const int z3, const int z4,
             const float coef0,
             vectorRealView const &coefx, vectorRealView const &coefy, vectorRealView const &coefz, 
             vectorRealView const &vp, vectorRealView const &phi, vectorRealView const &eta,
             vectorRealView const &pnp1, vectorRealView const &pn, vectorRealView const &pnm1) const
  {
#ifdef USE_OMP
      #pragma omp parallel for collapse(3)
#endif
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

       pnp1[IDX3_l(i,j,k)]=((2.-eta[IDX3_eta1(i,j,k)]*eta[IDX3_eta1(i,j,k)]+2.*eta[IDX3_eta1(i,j,k)])*pn[IDX3_l(i,j,k)]
                  -pnm1[IDX3_l(i,j,k)]+vp[IDX3(i,j,k)]*(lap+phi[IDX3(i,j,k)]))/(1.+2.*eta[IDX3_eta1(i,j,k)]);

       phi[IDX3(i,j,k)]=(phi[IDX3(i,j,k)]-((eta[IDX3_eta1(i+1,j,k)]-eta[IDX3_eta1(i-1,j,k)])*(pn[IDX3_l(i+1,j,k)]-pn[IDX3_l(i-1,j,k)])*myGrids.hdx_2
                         +(eta[IDX3_eta1(i,j+1,k)]-eta[IDX3_eta1(i,j-1,k)])*(pn[IDX3_l(i,j+1,k)]-pn[IDX3_l(i,j-1,k)])*myGrids.hdy_2
                         +(eta[IDX3_eta1(i,j,k+1)]-eta[IDX3_eta1(i,j,k-1)])*(pn[IDX3_l(i,j,k+1)]-pn[IDX3_l(i,j,k-1)])*myGrids.hdz_2))
                        /(1.+eta[IDX3_eta1(i,j,k)]);

     LOOP3DEND
     return(0);
  }

  // swap wavefields
  int swapWavefields( FDTDGRIDS &myGrids, 
		     vectorRealView const & pnp1,
		     vectorRealView const & pn  ,
		     vectorRealView const & pnm1) const
  {
    int nx=myGrids.nx;
    int ny=myGrids.ny;
    int nz=myGrids.nz;

#ifdef USE_RAJA
     RAJANestedLoop(0,0,0,nx,ny,nz)
     {
        pnm1[IDX3_l(i,j,k)]=pn[IDX3_l(i,j,k)];
        pn[IDX3_l(i,j,k)]=pnp1[IDX3_l(i,j,k)];
     });
#elif defined USE_KOKKOS
    int lx=myGrids.lx;
    int ly=myGrids.ly;
    int lz=myGrids.lz;

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
  int computeOneStep( FDTDGRIDS &myGrids, const float coef0,
                      vectorRealView const &coefx, vectorRealView const &coefy, vectorRealView const &coefz, 
                      vectorRealView const &vp, vectorRealView const &phi, vectorRealView const &eta, 
                      vectorRealView const &pnp1, vectorRealView const &pn, vectorRealView const &pnm1) const
  {
    //up

    pml3D(myGrids, 0, myGrids.nx, 0, myGrids.ny, myGrids.z1, myGrids.z2, coef0, coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
    //front
    pml3D(myGrids, 0, myGrids.nx, myGrids.y1, myGrids.y2, myGrids.z3, myGrids.z4, coef0, coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
    //left
    pml3D(myGrids, myGrids.x1, myGrids.x2, myGrids.y3, myGrids.y4, myGrids.z3, myGrids.z4, coef0, coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
    //inner points
    inner3D(myGrids, myGrids.x3, myGrids.x4, myGrids.y3, myGrids.y4, myGrids.z3, myGrids.z4, coef0, coefx,coefy,coefz,vp,pnp1,pn,pnm1);
    //right
    pml3D(myGrids, myGrids.x5, myGrids.x6, myGrids.y3, myGrids.y4, myGrids.z3, myGrids.z4, coef0, coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
    //back
    pml3D(myGrids, 0, myGrids.nx, myGrids.y5, myGrids.y6, myGrids.z3, myGrids.z4, coef0, coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
    // bottom
    pml3D(myGrids, 0, myGrids.nx, 0, myGrids.ny, myGrids.z5, myGrids.z6, coef0, coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);

    return(0);
  }

  // add RHS term
  int addRHS( FDTDGRIDS &myGrids,
             const int itSample,
             vectorRealView const & RHSTerm,
             vectorRealView const & vp, 
             vectorRealView const & pn) const
  {

     int xs=myGrids.xs;
     int ys=myGrids.ys;
     int zs=myGrids.zs;

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


};
#endif //FDTDKERNE_HPP

