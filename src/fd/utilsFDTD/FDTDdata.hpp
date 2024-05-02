#ifndef FDTDDATA_HPP_
#define FDTDDATA_HPP_

#include "FDTDmacros.hpp"

struct FDTDGRIDS {

   int nx, ny, nz;    // number of grid on the x,y,z direction
   int xs, ys, zs;    // source location on the x,y,z direction
   int lx, ly, lz;    
   float dx, dy, dz;    

   int ntaperx, ntapery, ntaperz;    
   int ndampx, ndampy, ndampz;    
   float hdx_2, hdy_2, hdz_2;    

   int x1, x2, x3, x4, x5, x6;
   int y1, y2, y3, y4, y5, y6;
   int z1, z2, z3, z4, z5, z6;

};

struct FDTDMODELS {

  VECTOR_REAL_VIEW vp; 
  VECTOR_REAL_VIEW phi;
  VECTOR_REAL_VIEW eta;
  VECTOR_REAL_VIEW pnp1;
  VECTOR_REAL_VIEW pn; 
  VECTOR_REAL_VIEW pnm1;

  VECTOR_REAL_VIEW RHSTerm;
  VECTOR_REAL_VIEW coefx;
  VECTOR_REAL_VIEW coefy;
  VECTOR_REAL_VIEW coefz;

};

#endif //FDTDDATA_HPP_
