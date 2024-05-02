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
         VECTOR_REAL_VIEW const &coefx, VECTOR_REAL_VIEW const &coefy, VECTOR_REAL_VIEW const &coefz, VECTOR_REAL_VIEW const &vp,
         VECTOR_REAL_VIEW const &pnp1, VECTOR_REAL_VIEW const &pn, VECTOR_REAL_VIEW const &pnm1) const
  {
      LOOP3DHEAD (x3,y3,z3,x4,y4,z4)
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
             VECTOR_REAL_VIEW const &coefx, VECTOR_REAL_VIEW const &coefy, VECTOR_REAL_VIEW const &coefz, 
             VECTOR_REAL_VIEW const &vp, VECTOR_REAL_VIEW const &phi, VECTOR_REAL_VIEW const &eta,
             VECTOR_REAL_VIEW const &pnp1, VECTOR_REAL_VIEW const &pn, VECTOR_REAL_VIEW const &pnm1) const
  {
       LOOP3DHEAD (x3,y3,z3,x4,y4,z4)
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
		     VECTOR_REAL_VIEW const & pnp1,
		     VECTOR_REAL_VIEW const & pn  ,
		     VECTOR_REAL_VIEW const & pnm1) const
  {
      LOOP3DHEAD(0,0,0,myGrids.nx,myGrids.ny,myGrids.nz)
        pnm1[IDX3_l(i,j,k)]=pn[IDX3_l(i,j,k)];
        pn[IDX3_l(i,j,k)]=pnp1[IDX3_l(i,j,k)];
      LOOP3DEND
      return(0);
  }

  // compute one step
  int computeOneStep( FDTDGRIDS &myGrids, const float coef0,
                      VECTOR_REAL_VIEW const &coefx, VECTOR_REAL_VIEW const &coefy, VECTOR_REAL_VIEW const &coefz, 
                      VECTOR_REAL_VIEW const &vp, VECTOR_REAL_VIEW const &phi, VECTOR_REAL_VIEW const &eta, 
                      VECTOR_REAL_VIEW const &pnp1, VECTOR_REAL_VIEW const &pn, VECTOR_REAL_VIEW const &pnm1) const
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
             VECTOR_REAL_VIEW const & RHSTerm,
             VECTOR_REAL_VIEW const & vp, 
             VECTOR_REAL_VIEW const & pn) const
  {
    LOOP3DHEAD(myGrids.xs,myGrids.ys,myGrids.zs,myGrids.xs+1,myGrids.ys+1,myGrids.zs+1)
       pn[IDX3_l(i,j,k)]+=vp[IDX3(i,j,k)]*RHSTerm[itSample];
    LOOP3DEND
    return(0);
  }


};
#endif //FDTDKERNE_HPP

