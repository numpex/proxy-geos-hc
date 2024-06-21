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
      float lapx=0;
      for (int l=1; l<coefx.size();l++)
      {
         lapx+=coefx[l]*(pnGlobal(IDX3_l(i+l,j,k),cb)+pnGlobal(IDX3_l(i-l,j,k),cb));
      }
      float lapy=0;
      for (int l=1; l<coefy.size();l++)
      {
         lapy+=coefy[l]*(pnGlobal(IDX3_l(i,j+l,k),cb)+pnGlobal(IDX3_l(i,j-l,k),cb));
      }
      float lapz=0;
      for (int l=1; l<coefz.size();l++)
      {
         lapz+=coefz[l]*(pnGlobal(IDX3_l(i,j,k+l),cb)+pnGlobal(IDX3_l(i,j,k-l),cb));
      }
      pnGlobal(IDX3_l(i,j,k),ca)=2.*pnGlobal(IDX3_l(i,j,k),cb)-pnGlobal(IDX3_l(i,j,k),ca)
                         +vp[IDX3(i,j,k)]*(coef0*pnGlobal(IDX3_l(i,j,k),cb)+lapx+lapy+lapz);
      LOOP3DEND
      return 0;
  }


  // compute one step
  int computeOneStep( FDTDGRIDS &myGrids, FDTDMODELS &myModels, int & ca, int & cb ) const
  {
    //inner points
    inner3D(myGrids, myModels, ca, cb, myGrids.x3, myGrids.x4, myGrids.y3, myGrids.y4, myGrids.z3, myGrids.z4, myModels.pnGlobal);
    // apply sponge boundary to wavefield
    applySponge(myGrids,myModels,ca,cb,myGrids.x3, myGrids.x4, myGrids.y3, myGrids.y4, myGrids.z3, myGrids.z4,
                myModels.pnGlobal);

    return(0);
  }

  //innerpoints
  int applySponge( FDTDGRIDS &myGrids, FDTDMODELS &myModels,int & ca, int & cb, 
                    const int x3, const int x4, const int y3, const int y4, const int z3, const int z4, 
                    arrayReal const &PN_Global) const
  {
      CREATEVIEWSPONGE
      LOOP3DHEAD (x3,y3,z3,x4,y4,z4)
      pnGlobal(IDX3_l(i,j,k),ca)=pnGlobal(IDX3_l(i,j,k),ca)*spongeArray(IDX3(i,j,k));
      pnGlobal(IDX3_l(i,j,k),cb)=pnGlobal(IDX3_l(i,j,k),cb)*spongeArray(IDX3(i,j,k));
      LOOP3DEND
      return 0;
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

