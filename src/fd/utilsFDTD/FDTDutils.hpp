#ifndef FDTDUTILS_HPP
#define FDTDUTILS_HPP

#include "dataType.hpp"
#include "FDTDdata.hpp"

using namespace std;

using namespace std;

struct FDTDUtils
{
  void init_coef(int L,float dx, vectorReal &coef)
  {
      float dx2 = dx*dx;
      switch(L)
      {
        case 1:
          coef[0]=-2.f/dx2;
          coef[1]=1.f/dx2;
          break;
        case 2:
          coef[0]=-5.f/2.f/dx2;
          coef[1]=4.f/3.f/dx2;
          coef[2]=-1.f/12.f/dx2;
          break;
        case 3:
          coef[0]=-49.f/18.f/dx2;
          coef[1]=3.f/2.f/dx2;
          coef[2]=-3./20./dx2;
          coef[3]=1./90./dx2;
          break;
        case 4:
          coef[0] = -205.f/72.f/dx2;
          coef[1] = 8.f/5.f/dx2;
          coef[2] = -1.f/5.f/dx2;
          coef[3] = 8.f/315.f/dx2;
          coef[4] = -1.f/560.f/dx2;
          break;
        case 5:
          coef[0]=-5269.f/1800.f/dx2;
          coef[1]=5.f/3.f/dx2;
          coef[2]=-5.f/21.f/dx2;
          coef[3]=5.f/126.f/dx2;
          coef[4]=-5.f/1008.f/dx2;
          coef[5]=1.f/3150.f/dx2;
          break;
        case 6:
          coef[0]=-5369.f/1800.f/dx2;
          coef[1]=12.f/7.f/dx2;
          coef[2]=-15.f/56.f/dx2;
          coef[3]=10.f/189.f/dx2; 
          coef[4]=-1.f/112.f/dx2;
          coef[5]=2.f/1925.f/dx2;
          coef[6]=1.f/1663.f/dx2;
          break;
      }
  }
  float compute_dt_sch(const float vmax, vectorReal const &coefx, vectorReal const &coefy, vectorReal const &coefz) 
  {

      float ftmp = 0.;
      float cfl=0.8;
      ftmp += fabsf(coefx[0]) + fabsf(coefy[0]) + fabsf(coefz[0]);
      for (int i = 1; i < coefx.size(); i++) {
          ftmp += 2.f*fabsf(coefx[i]);
      }
      for (int i = 1; i < coefy.size(); i++) {
          ftmp += 2.f*fabsf(coefy[i]);
      }
      for (int i = 1; i < coefz.size(); i++) {
          ftmp += 2.f*fabsf(coefz[i]);
      }
      printf( " coefs computed\n" );
      return 2*cfl/(sqrtf(ftmp)*vmax);
  }

      
  void output(FDTDGRIDS &myGrids, arrayReal const&pnGlobal, int itSample, const int &i1,const bool saveSnapShots)
  {

      if(itSample%50==0)
      {   

        FDFENCE

        printf("TimeStep=%d\t; Pressure value at source [%d %d %d] = %f\n", itSample,
               myGrids.xs, myGrids.ys, myGrids.zs, pnGlobal(IDX3_l(myGrids.xs,myGrids.ys,myGrids.zs),i1));
        //#ifdef FD_SAVE_SNAPSHOTS
        if(saveSnapShots)
            write_io( myGrids, 0, myGrids.nx, myGrids.ny/2, myGrids.ny/2, 0, myGrids.nz, pnGlobal, itSample, i1);
        //#endif

      } 

  }

  void write_io( FDTDGRIDS &myGrids,
		 int x0, int x1, 
		 int y0, int y1, 
		 int z0, int z1, 
                 arrayReal const&pnGlobal, int istep, const int &i1)
  {
      char filename_buf[32];
      snprintf(filename_buf, sizeof(filename_buf), "snapshot_it_%d.H@", istep);
      FILE *snapshot_file = fopen(filename_buf, "wb");
      //printf("write snapshot for: x=[%d %d], y=[%d %d], z=[%d %d]\n",x0,x1,y0,y1,z0,z1);
     
      for (int k = z0; k < z1; ++k) {
          for (int j = y0; j < y1+1; ++j) {
              for (int i = x0; i < x1; ++i) {
		  fwrite(&pnGlobal(IDX3_l(i,j,k),i1), sizeof(float),1, snapshot_file);
              }
          }
      }
      // Clean up 
      fclose(snapshot_file);
  }

};

#endif  //FDTDUTILS_HPP
