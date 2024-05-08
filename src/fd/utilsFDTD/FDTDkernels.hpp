#ifndef FDTDKERNEL_HPP
#define FDTDKERNEL_HPP

struct FDTDKernel
{

 //innerpoints
  int inner3D( FDTDGRIDS &myGrids, FDTDMODELS &myModels, int & ca, int & cb, 
         const int x3, const int x4,
         const int y3, const int y4,
         const int z3, const int z4, arrayReal const &PN_Global) const
  {
      CREATEVIEWINNER
      LOOP3DHEAD (x3,y3,z3,x4,y4,z4)
      float lapx=(coefx[1]*(pnGlobal(IDX3_l(i+1,j,k),cb)+pnGlobal(IDX3_l(i-1,j,k),cb))
                 +coefx[2]*(pnGlobal(IDX3_l(i+2,j,k),cb)+pnGlobal(IDX3_l(i-2,j,k),cb))
                 +coefx[3]*(pnGlobal(IDX3_l(i+3,j,k),cb)+pnGlobal(IDX3_l(i-3,j,k),cb))
                 +coefx[4]*(pnGlobal(IDX3_l(i+4,j,k),cb)+pnGlobal(IDX3_l(i-4,j,k),cb)));
      float lapy=(coefy[1]*(pnGlobal(IDX3_l(i,j+1,k),cb)+pnGlobal(IDX3_l(i,j-1,k),cb))
                 +coefy[2]*(pnGlobal(IDX3_l(i,j+2,k),cb)+pnGlobal(IDX3_l(i,j-2,k),cb))
                 +coefy[3]*(pnGlobal(IDX3_l(i,j+3,k),cb)+pnGlobal(IDX3_l(i,j-3,k),cb))
                 +coefy[4]*(pnGlobal(IDX3_l(i,j+4,k),cb)+pnGlobal(IDX3_l(i,j-4,k),cb)));
      float lapz=(coefz[1]*(pnGlobal(IDX3_l(i,j,k+1),cb)+pnGlobal(IDX3_l(i,j,k-1),cb))
                 +coefz[2]*(pnGlobal(IDX3_l(i,j,k+2),cb)+pnGlobal(IDX3_l(i,j,k-2),cb))
                 +coefz[3]*(pnGlobal(IDX3_l(i,j,k+3),cb)+pnGlobal(IDX3_l(i,j,k-3),cb))
                 +coefz[4]*(pnGlobal(IDX3_l(i,j,k+4),cb)+pnGlobal(IDX3_l(i,j,k-4),cb)));
      pnGlobal(IDX3_l(i,j,k),ca)=2.*pnGlobal(IDX3_l(i,j,k),cb)-pnGlobal(IDX3_l(i,j,k),ca)
                         +vp[IDX3(i,j,k)]*(coef0*pnGlobal(IDX3_l(i,j,k),cb)+lapx+lapy+lapz);
      LOOP3DEND
      return 0;
  }

  int pml3D( FDTDGRIDS &myGrids, FDTDMODELS &myModels, int & ca, int & cb,
             const int x3, const int x4, 
             const int y3, const int y4, 
             const int z3, const int z4, arrayReal const& PN_Global ) const
  {
       CREATEVIEWPML
       LOOP3DHEAD (x3,y3,z3,x4,y4,z4)
       float lapx=(coefx[1]*(pnGlobal(IDX3_l(i+1,j,k),cb)+pnGlobal(IDX3_l(i-1,j,k),cb))
                 +coefx[2]*(pnGlobal(IDX3_l(i+2,j,k),cb)+pnGlobal(IDX3_l(i-2,j,k),cb))
                 +coefx[3]*(pnGlobal(IDX3_l(i+3,j,k),cb)+pnGlobal(IDX3_l(i-3,j,k),cb))
                 +coefx[4]*(pnGlobal(IDX3_l(i+4,j,k),cb)+pnGlobal(IDX3_l(i-4,j,k),cb)));
      float lapy=(coefy[1]*(pnGlobal(IDX3_l(i,j+1,k),cb)+pnGlobal(IDX3_l(i,j-1,k),cb))
                 +coefy[2]*(pnGlobal(IDX3_l(i,j+2,k),cb)+pnGlobal(IDX3_l(i,j-2,k),cb))
                 +coefy[3]*(pnGlobal(IDX3_l(i,j+3,k),cb)+pnGlobal(IDX3_l(i,j-3,k),cb))
                 +coefy[4]*(pnGlobal(IDX3_l(i,j+4,k),cb)+pnGlobal(IDX3_l(i,j-4,k),cb)));
      float lapz=(coefz[1]*(pnGlobal(IDX3_l(i,j,k+1),cb)+pnGlobal(IDX3_l(i,j,k-1),cb))
                 +coefz[2]*(pnGlobal(IDX3_l(i,j,k+2),cb)+pnGlobal(IDX3_l(i,j,k-2),cb))
                 +coefz[3]*(pnGlobal(IDX3_l(i,j,k+3),cb)+pnGlobal(IDX3_l(i,j,k-3),cb))
                 +coefz[4]*(pnGlobal(IDX3_l(i,j,k+4),cb)+pnGlobal(IDX3_l(i,j,k-4),cb)));

       float lap=coef0*pnGlobal(IDX3_l(i,j,k),cb)+lapx+lapy+lapz;

       pnGlobal(IDX3_l(i,j,k),ca)=((2.-eta[IDX3_eta1(i,j,k)]*eta[IDX3_eta1(i,j,k)]+2.*eta[IDX3_eta1(i,j,k)])
                                     *pnGlobal(IDX3_l(i,j,k),cb)-pnGlobal(IDX3_l(i,j,k),ca)
                                     +vp[IDX3(i,j,k)]*(lap+phi[IDX3(i,j,k)]))/(1.+2.*eta[IDX3_eta1(i,j,k)]);

       phi[IDX3(i,j,k)]=(phi[IDX3(i,j,k)]
         -((eta[IDX3_eta1(i+1,j,k)]-eta[IDX3_eta1(i-1,j,k)])*(pnGlobal(IDX3_l(i+1,j,k),cb)-pnGlobal(IDX3_l(i-1,j,k),cb))*myGrids.hdx_2
          +(eta[IDX3_eta1(i,j+1,k)]-eta[IDX3_eta1(i,j-1,k)])*(pnGlobal(IDX3_l(i,j+1,k),cb)-pnGlobal(IDX3_l(i,j-1,k),cb))*myGrids.hdy_2
          +(eta[IDX3_eta1(i,j,k+1)]-eta[IDX3_eta1(i,j,k-1)])*(pnGlobal(IDX3_l(i,j,k+1),cb)-pnGlobal(IDX3_l(i,j,k-1),cb))*myGrids.hdz_2))
         /(1.+eta[IDX3_eta1(i,j,k)]);

     LOOP3DEND
     return(0);
  }

  // compute one step
  int computeOneStep( FDTDGRIDS &myGrids, FDTDMODELS &myModels, int & ca, int & cb ) const
  {
    //up
    pml3D(myGrids, myModels, ca, cb, 0, myGrids.nx, 0, myGrids.ny, myGrids.z1, myGrids.z2, myModels.pnGlobal);
    //front
    pml3D(myGrids, myModels, ca, cb, 0, myGrids.nx, myGrids.y1, myGrids.y2, myGrids.z3, myGrids.z4, myModels.pnGlobal);
    //left
    pml3D(myGrids, myModels, ca, cb, myGrids.x1, myGrids.x2, myGrids.y3, myGrids.y4, myGrids.z3, myGrids.z4, myModels.pnGlobal);
    //inner points
    inner3D(myGrids, myModels, ca, cb, myGrids.x3, myGrids.x4, myGrids.y3, myGrids.y4, myGrids.z3, myGrids.z4, myModels.pnGlobal);
    //right
    pml3D(myGrids, myModels, ca, cb, myGrids.x5, myGrids.x6, myGrids.y3, myGrids.y4, myGrids.z3, myGrids.z4, myModels.pnGlobal);
    //back
    pml3D(myGrids, myModels, ca, cb, 0, myGrids.nx, myGrids.y5, myGrids.y6, myGrids.z3, myGrids.z4, myModels.pnGlobal);
    // bottom
    pml3D(myGrids, myModels, ca, cb, 0, myGrids.nx, 0, myGrids.ny, myGrids.z5, myGrids.z6, myModels.pnGlobal);

    return(0);
  }

  // add RHS term
  int addRHS( FDTDGRIDS &myGrids, const int itSample, FDTDMODELS &myModels,  int & cb, arrayReal const &PN_Global) 
  {
    CREATEVIEWRHS
    LOOP3DHEAD(myGrids.xs,myGrids.ys,myGrids.zs,myGrids.xs+1,myGrids.ys+1,myGrids.zs+1)
       pnGlobal(IDX3_l(i,j,k), cb)+=vp[IDX3(i,j,k)]*RHSTerm[itSample];
    LOOP3DEND
    return(0);
  }


};
#endif //FDTDKERNE_HPP

