#ifndef FDTDKERNEL_HPP
#define FDTDKERNEL_HPP

struct FDTDKernel
{

 //innerpoints
  int inner3D( FDTDGRIDS &myGrids, FDTDMODELS &myModels,
         const int x3, const int x4,
         const int y3, const int y4,
         const int z3, const int z4 ) const
  {
      CREATEVIEWINNER
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
      pnp1[IDX3_l(i,j,k)]=2.*pn[IDX3_l(i,j,k)]-pnGlobal(IDX3_l(i,j,k),0)+vp[IDX3(i,j,k)]*(coef0*pn[IDX3_l(i,j,k)]+lapx+lapy+lapz);
      LOOP3DEND
      return 0;
  }

  int pml3D( FDTDGRIDS &myGrids, FDTDMODELS &myModels,
             const int x3, const int x4, 
             const int y3, const int y4, 
             const int z3, const int z4 ) const
  {
       CREATEVIEWPML
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
                  -pnGlobal(IDX3_l(i,j,k),0)+vp[IDX3(i,j,k)]*(lap+phi[IDX3(i,j,k)]))/(1.+2.*eta[IDX3_eta1(i,j,k)]);

       phi[IDX3(i,j,k)]=(phi[IDX3(i,j,k)]-((eta[IDX3_eta1(i+1,j,k)]-eta[IDX3_eta1(i-1,j,k)])*(pn[IDX3_l(i+1,j,k)]-pn[IDX3_l(i-1,j,k)])*myGrids.hdx_2
                         +(eta[IDX3_eta1(i,j+1,k)]-eta[IDX3_eta1(i,j-1,k)])*(pn[IDX3_l(i,j+1,k)]-pn[IDX3_l(i,j-1,k)])*myGrids.hdy_2
                         +(eta[IDX3_eta1(i,j,k+1)]-eta[IDX3_eta1(i,j,k-1)])*(pn[IDX3_l(i,j,k+1)]-pn[IDX3_l(i,j,k-1)])*myGrids.hdz_2))
                        /(1.+eta[IDX3_eta1(i,j,k)]);

     LOOP3DEND
     return(0);
  }

  // swap wavefields
  int swapWavefields( FDTDGRIDS &myGrids, FDTDMODELS &myModels ) const
  {
      CREATEVIEWAVEFD
      LOOP3DHEAD(0,0,0,myGrids.nx,myGrids.ny,myGrids.nz)
        pnGlobal(IDX3_l(i,j,k),0)=pn[IDX3_l(i,j,k)];
        pn[IDX3_l(i,j,k)]=pnp1[IDX3_l(i,j,k)];
      LOOP3DEND
      return(0);
  }

  // compute one step
  int computeOneStep( FDTDGRIDS &myGrids, FDTDMODELS &myModels ) const
  {
    //up
    pml3D(myGrids, myModels, 0, myGrids.nx, 0, myGrids.ny, myGrids.z1, myGrids.z2);
    //front
    pml3D(myGrids, myModels, 0, myGrids.nx, myGrids.y1, myGrids.y2, myGrids.z3, myGrids.z4);
    //left
    pml3D(myGrids, myModels, myGrids.x1, myGrids.x2, myGrids.y3, myGrids.y4, myGrids.z3, myGrids.z4);
    //inner points
    inner3D(myGrids, myModels, myGrids.x3, myGrids.x4, myGrids.y3, myGrids.y4, myGrids.z3, myGrids.z4);
    //right
    pml3D(myGrids, myModels, myGrids.x5, myGrids.x6, myGrids.y3, myGrids.y4, myGrids.z3, myGrids.z4);
    //back
    pml3D(myGrids, myModels, 0, myGrids.nx, myGrids.y5, myGrids.y6, myGrids.z3, myGrids.z4);
    // bottom
    pml3D(myGrids, myModels, 0, myGrids.nx, 0, myGrids.ny, myGrids.z5, myGrids.z6);

    return(0);
  }

  // add RHS term
  int addRHS( FDTDGRIDS &myGrids, const int itSample, FDTDMODELS &myModels ) 
  {
    CREATEVIEWRHS
    LOOP3DHEAD(myGrids.xs,myGrids.ys,myGrids.zs,myGrids.xs+1,myGrids.ys+1,myGrids.zs+1)
       pn[IDX3_l(i,j,k)]+=vp[IDX3(i,j,k)]*RHSTerm[itSample];
    LOOP3DEND
    return(0);
  }


};
#endif //FDTDKERNE_HPP

