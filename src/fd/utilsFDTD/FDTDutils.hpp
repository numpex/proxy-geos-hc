#ifndef FDTDUTILS_HPP
#define FDTDUTILS_HPP

#include "dataType.hpp"
#include "FDTDdata.hpp"
#include "FDTDmacros.hpp"

struct FDTDUtils
{
  void init_coef(float dx, vectorReal &coef)
  {
      float dx2 = dx*dx;
      coef[0] = -205.f/72.f/dx2;
      coef[1] = 8.f/5.f/dx2;
      coef[2] = -1.f/5.f/dx2;
      coef[3] = 8.f/315.f/dx2;
      coef[4] = -1.f/560.f/dx2;
  }
  float compute_dt_sch(const float vmax,const vectorReal &coefx,const vectorReal &coefy,const vectorReal &coefz) 
  {

      float ftmp = 0.;
      float cfl=0.8;
      ftmp += fabsf(coefx[0]) + fabsf(coefy[0]) + fabsf(coefz[0]);
      for (int i = 1; i < 5; i++) {
          ftmp += 2.f*fabsf(coefx[i]);
          ftmp += 2.f*fabsf(coefy[i]);
          ftmp += 2.f*fabsf(coefz[i]);
      }
      return 2*cfl/(sqrtf(ftmp)*vmax);
  }

  void write_io(int nx, int ny, int nz,
		int lx, int ly, int lz,
		int x0, int x1, 
		int y0, int y1, 
		int z0, int z1, 
                vectorReal const &u, int istep)
  {
      char filename_buf[32];
      snprintf(filename_buf, sizeof(filename_buf), "snapshot_it_%d.H@", istep);
      FILE *snapshot_file = fopen(filename_buf, "wb");
      printf(" %d %d %d %d %d %d\n",x0,x1,y0,y1,z0,z1);
      for (int k = z0; k < z1; ++k) {
          for (int j = y0; j < y1+1; ++j) {
              for (int i = x0; i < x1; ++i) {
		  fwrite(&u[IDX3_l(i,j,k)], sizeof(float),1, snapshot_file);
              }
          }
      }
      /* Clean up */
      fclose(snapshot_file);
  }

  void pml_profile_init(vectorReal  & profile, int i_min, int i_max, int n_first, int n_last, float scale)
  {
    int n = i_max-i_min+1;
    int shift = i_min-1;

    int first_beg = 1 + shift;
    int first_end = n_first + shift;
    int last_beg  = n - n_last+1 + shift;
    int last_end  = n + shift;

    for (int i = i_min; i <= i_max; ++i) {
        profile[i] = 0.f;
    }

    float tmp = scale / POW2(first_end-first_beg+1);
    for (int i = 1; i <= first_end-first_beg+1; ++i) {
        profile[first_end-i+1] = POW2(i)*tmp;
    }

    for (int i = 1; i <= last_end-last_beg+1; ++i) {
        profile[last_beg+i-1] = POW2(i)*tmp;
    }
  }

  void pml_profile_extend(int nx, int ny, int nz,
                          vectorReal &  eta, vectorReal & etax, vectorReal & etay, vectorReal & etaz,
                          int xbeg, int xend, int ybeg, int yend, int zbeg, int zend)
  {
    const int n_ghost = 1;
    for (int ix = xbeg-n_ghost; ix <= xend+n_ghost; ++ix) {
        for (int iy = ybeg-n_ghost; iy <= yend+n_ghost; ++iy) {
            for (int iz = zbeg-n_ghost; iz <= zend+n_ghost; ++iz) {
                eta[IDX3_eta0(ix,iy,iz)] = etax[ix] + etay[iy] + etaz[iz];
            }
        }
    }
  }

  void pml_profile_extend_all(int nx, int ny, int nz,
                              vectorReal & eta, vectorReal & etax, vectorReal & etay, vectorReal & etaz,
                              int xmin, int xmax, int ymin, int ymax,
                              int x1, int x2, int x5, int x6,
                              int y1, int y2, int y3, int y4, int y5, int y6,
                              int z1, int z2, int z3, int z4, int z5, int z6)
  {
    // Top.
    if (z1 != -1)
    pml_profile_extend(nx,ny,nz,eta,etax,etay,etaz,xmin,xmax,ymin,ymax,z1,z2);
    // Bottom.
    if (z5 != -5)
    pml_profile_extend(nx,ny,nz,eta,etax,etay,etaz,xmin,xmax,ymin,ymax,z5,z6);
    // Front.
    if ((y1!=-1) && (z3!=-3))
    pml_profile_extend(nx,ny,nz,eta,etax,etay,etaz,xmin,xmax,y1,y2,z3,z4);
    // Back.
    if ((y6!=-6) && (z3!=-3))
    pml_profile_extend(nx,ny,nz,eta,etax,etay,etaz,xmin,xmax,y5,y6,z3,z4);
    // Left.
    if ((x1!=-1) && (y3!=-3) && (z3!=-3))
    pml_profile_extend(nx,ny,nz,eta,etax,etay,etaz,x1,x2,y3,y4,z3,z4);
    // Right.
    if ((x6!=-6) && (y3!=-3) && (z3!=-3))
    pml_profile_extend(nx,ny,nz,eta,etax,etay,etaz,x5,x6,y3,y4,z3,z4);
  }

  void init_eta(int nx, int ny, int nz,
                int ndampx, int ndampy,int ndampz,
                int x1, int x2, int x3, int x4, int x5, int x6,
                int y1, int y2, int y3, int y4, int y5, int y6,
                int z1, int z2, int z3, int z4, int z5, int z6,
                float dx, float dy, float dz, float dt_sch,
                float vmax, vectorReal & eta)
  {
    for (int i = -1; i < nx+1; ++i) {
        for (int j = -1; j < ny+1; ++j) {
            for (int k = -1; k < nz+1; ++k) {
                eta[IDX3_eta1(i,j,k)] = 0.f;
            }
        }
    }

    /* etax */
    float param = dt_sch * 3.f * vmax * logf(1000.f)/(2.f*ndampx*dx);
    printf("param=%f\n",param);
    vectorReal etax=allocateVector<vectorReal>(nx+2);
    pml_profile_init(etax, 0, nx+1, ndampx, ndampx, param);

    /* etay */
    param = dt_sch*3.f*vmax*logf(1000.f)/(2.f*ndampy*dy);
    printf("param=%f\n",param);
    vectorReal etay=allocateVector<vectorReal>(ny+2);
    pml_profile_init(etay, 0, ny+1, ndampy, ndampy, param);

    /* etaz */
    param = dt_sch*3.f*vmax*logf(1000.f)/(2.f*ndampz*dz);
    printf("param=%f\n",param);
    vectorReal etaz=allocateVector<vectorReal>(nz+2);
    pml_profile_init(etaz, 0, nz+1, ndampz, ndampz, param);

    (void)pml_profile_extend_all(nx, ny, nz,
                eta, etax, etay, etaz,
                1, nx, 1, ny,
                x1+1, x2, x5+1, x6,
                y1+1, y2, y3+1, y4, y5+1, y6,
                z1+1, z2, z3+1, z4, z5+1, z6);

  }
      
  //void output(FDTDGRIDS &myGrids, vectorRealView const & pn, int itSample)
  void output(FDTDGRIDS &myGrids, vectorRealView &pn, int itSample)
  {

      if(itSample%50==0)
      {   

        #ifdef USE_RAJA
        RAJA::forall<RAJA::omp_parallel_for_exec>(RAJA::RangeSegment(0,nx), [pn] ( int i){});
        #elif defined USE_KOKKOS
        Kokkos::fence();
        #endif

        //printf("result1 %f\n",pn[IDX3_l(myGrids.xs,myGrids.ys,myGrids.zs)]);

        printf("result1 %f\n", 
        pn[(myGrids.nz+2*myGrids.lz)*(myGrids.ny+2*myGrids.ly)*(myGrids.xs+myGrids.lx) 
                     +(myGrids.nz+2*myGrids.lz)*(myGrids.ys+myGrids.ly)+(myGrids.zs+myGrids.lz)]);

        write_io( myGrids.nx ,myGrids.ny ,myGrids.nz,
                              myGrids.lx ,myGrids.ly ,myGrids.lz,
                              0 ,myGrids.nx ,myGrids.ny/2 ,myGrids.ny/2, 0 ,myGrids.nz,
                              pn, itSample);

      } 

  }



};

#endif  //FDTDUTILS_HPP
