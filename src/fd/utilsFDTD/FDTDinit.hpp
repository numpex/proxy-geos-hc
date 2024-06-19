#ifndef FDTDINIT_HPP
#define FDTDINIT_HPP

#include "utils.hpp"
#include "FDTDutils.hpp"
#include "FDTDkernels.hpp"


struct FDTDInit
{
  FDTDUtils myFDTDUtils;
  SolverUtils myUtils;

  int sourceOrder=1;
  int ncoefsX;
  int ncoefsY;
  int ncoefsZ;

  float f0=15.;
  float fmax=2.5*f0;
  float timeMax=1.0;
  float timeStep;
  int nSamples;

  float vmin=2500;
  float vmax=4500;

  int i1=0;
  int i2=1;

  vectorReal RHSTerm;

  void init_geometry( int argc, char *argv[], FDTDGRIDS & myGrids )
  {
    myGrids.nx = 150;
    myGrids.ny = myGrids.nx;
    myGrids.nz = myGrids.nx;

    myGrids.xs = myGrids.nx/2;
    myGrids.ys = myGrids.ny/2;
    myGrids.zs = myGrids.nz/2;
    myGrids.lx=4;
    myGrids.ly=4;
    myGrids.lz=4;
    ncoefsX=5;
    ncoefsY=5;
    ncoefsZ=5;

    myGrids.dx=10;
    myGrids.dy=10;
    myGrids.dz=10;
    if(argc>1)
    {
        for (int i=1; i<argc; i=i+2)
        {
            // dimensions nx,ny,nz
            printf("arg %s %s  \n",argv[i],argv[i+1]);
            std::string arg=argv[i];
            if (arg=="-nx") 
            {
            printf("arg in %s %s  \n",argv[i],argv[i+1]);
                myGrids.nx=atoi(argv[i+1]);
                myGrids.ny=myGrids.nx;
                myGrids.nz=myGrids.nx; 
            }
            if (arg=="-ny") myGrids.ny=atoi(argv[i+1]);
            if (arg=="-nz") myGrids.nz=atoi(argv[i+1]);
           
            // spatial sampling dx dy dz
            if (arg=="-dx") 
            {
                myGrids.dx=atoi(argv[i+1]) ;
                myGrids.dy=myGrids.dx;
                myGrids.dz=myGrids.dx;
            }
            if (arg=="-dy") myGrids.dy=atoi(argv[i+1]) ;
            if (arg=="-dz") myGrids.dz=atoi(argv[i+1]) ;

            // half stencil length lx ly lz
            if (arg=="-lx") 
            {
                myGrids.lx=atoi(argv[i+1]);
                myGrids.ly=myGrids.lx;
                myGrids.lz=myGrids.lx; 
            }
            if (arg=="-ly") myGrids.ly=atoi(argv[i+1]);
            if (arg=="-lz") myGrids.lz=atoi(argv[i+1]);
            ncoefsX=myGrids.lx+1;
            ncoefsY=myGrids.ly+1;
            ncoefsZ=myGrids.lz+1;

            // source location xs ys zs
            if (arg=="-xs") 
            {
                myGrids.xs=atoi(argv[i+1]) ;
                myGrids.ys=myGrids.xs;
                myGrids.zs=myGrids.zs;
            }
            if (arg=="-ys") myGrids.ys=atoi(argv[i+1]) ;
            if (arg=="-zs") myGrids.ys=atoi(argv[i+1]) ;

        }
    }

    printf( "Number of grids: nx=%d, ny=%d, nz=%d\n", myGrids.nx, myGrids.ny, myGrids.nz);
    printf( "Source location: xs=%d, xy=%d, xz=%d\n", myGrids.xs, myGrids.ys, myGrids.zs);
    float lambdamax=vmin/fmax;

    // init pml limits
    myGrids.ntaperx=4;
    myGrids.ntapery=4;
    myGrids.ntaperz=4;

    myGrids.hdx_2=1./(4. * myGrids.dx * myGrids.dx);
    myGrids.hdy_2=1./(4. * myGrids.dy * myGrids.dy);
    myGrids.hdz_2=1./(4. * myGrids.dz * myGrids.dz);

    myGrids.ndampx=myGrids.ntaperx * lambdamax / myGrids.dx;
    myGrids.ndampy=myGrids.ntapery * lambdamax / myGrids.dy;
    myGrids.ndampz=myGrids.ntaperz * lambdamax / myGrids.dz;

    myGrids.x1=0;
    myGrids.x2=myGrids.ndampx;
    myGrids.x3=myGrids.ndampx;
    myGrids.x4=myGrids.nx-myGrids.ndampx;
    myGrids.x5=myGrids.nx-myGrids.ndampx;
    myGrids.x6=myGrids.nx;

    myGrids.y1=0;
    myGrids.y2=myGrids.ndampy;
    myGrids.y3=myGrids.ndampy;
    myGrids.y4=myGrids.ny-myGrids.ndampy;
    myGrids.y5=myGrids.ny-myGrids.ndampy;
    myGrids.y6=myGrids.ny;

    myGrids.z1=0;
    myGrids.z2=myGrids.ndampz;
    myGrids.z3=myGrids.ndampz;
    myGrids.z4=myGrids.nz-myGrids.ndampz;
    myGrids.z5=myGrids.nz-myGrids.ndampz;
    myGrids.z6=myGrids.nz;

  }

  void init_coefficients( FDTDGRIDS & myGrids, FDTDMODELS & myModels )
  {
    myModels.coefx = allocateVector< vectorReal >( ncoefsX, "coefx" );
    myModels.coefy = allocateVector< vectorReal >( ncoefsY, "coefy" );
    myModels.coefz = allocateVector< vectorReal >( ncoefsZ, "coefz" );

    myFDTDUtils.init_coef( myGrids.lx, myGrids.dx, myModels.coefx );
    myFDTDUtils.init_coef( myGrids.ly, myGrids.dy, myModels.coefy );
    myFDTDUtils.init_coef( myGrids.lz, myGrids.dz, myModels.coefz );

    float tmpX=0;
    for( int i=1;i<ncoefsX;i++ )
    {
       tmpX+=myModels.coefx[i];
    }
    float tmpY=0;
    for( int i=1;i<ncoefsY;i++ )
    {
       tmpY+=myModels.coefy[i];
    }
    float tmpZ=0;
    for( int i=1;i<ncoefsZ;i++ )
    {
       tmpZ+=myModels.coefz[i];
    }
    myModels.coef0 = -2.*(tmpX+tmpY+tmpZ);
    /*
    myModels.coef0 = -2.*(myModels.coefx[1]+myModels.coefx[2]+myModels.coefx[3]+myModels.coefx[4]);
    myModels.coef0+= -2.*(myModels.coefy[1]+myModels.coefy[2]+myModels.coefy[3]+myModels.coefy[4]);
    myModels.coef0+= -2.*(myModels.coefz[1]+myModels.coefz[2]+myModels.coefz[3]+myModels.coefz[4]);
    */

    timeStep=myFDTDUtils.compute_dt_sch(vmax, myModels.coefx, myModels.coefy, myModels.coefz );
    nSamples=timeMax/timeStep;
    printf("init coefs done\n");

  }

  void init_source( FDTDMODELS & myModels )
  {
    // compute source term
    myModels.RHSTerm = allocateVector< vectorReal >( nSamples, "RHSTerm" );

    std::vector< float > sourceTerm=myUtils.computeSourceTerm( nSamples, timeStep, f0, sourceOrder );
    for( int i=0; i<nSamples; i++ )
    {
      myModels.RHSTerm[i]=sourceTerm[i];
      //cout<<"sample "<<i<<"\t: sourceTerm = "<<sourceTerm[i]<< endl;
    }

  }

  void init_models( FDTDGRIDS & myGrids, FDTDMODELS & myModels )
  {
    int modelVolume = myGrids.nx * myGrids.ny * myGrids.nz;
    int extModelVolume = ( myGrids.nx + 2 * myGrids.lx ) *
                         ( myGrids.ny + 2 * myGrids.ly ) *
                         ( myGrids.nz + 2 * myGrids.lz );
    int etaModelVolume = ( myGrids.nx + 2 ) * ( myGrids.ny + 2 ) * ( myGrids.nz + 2 );

    myModels.vp   = allocateVector< vectorReal >( modelVolume, "vp" );
    myModels.pnp1 = allocateVector< vectorReal >( extModelVolume, "pnp1" );
    myModels.pn   = allocateVector< vectorReal >( extModelVolume, "pn" );
    myModels.phi  = allocateVector< vectorReal >( modelVolume, "phi" );
    myModels.eta  = allocateVector< vectorReal >( etaModelVolume, "eta" );
    myModels.pnGlobal = allocateArray2D< arrayReal >( extModelVolume, 2, "pnGlobal" );

    myFDTDUtils.init_eta( myGrids.nx, myGrids.ny, myGrids.nz,
                          myGrids.ndampx, myGrids.ndampy, myGrids.ndampz,
                          myGrids.x1, myGrids.x2, myGrids.x3, myGrids.x4, myGrids.x5, myGrids.x6,
                          myGrids.y1, myGrids.y2, myGrids.y3, myGrids.y4, myGrids.y5, myGrids.y6,
                          myGrids.z1, myGrids.z2, myGrids.z3, myGrids.z4, myGrids.z5, myGrids.z6,
                          myGrids.dx, myGrids.dy, myGrids.dz, timeStep,
                          vmax, myModels.eta );

    float init_vp_value = vmin*vmin*timeStep*timeStep;

    // initialize vp and pressure field
   #pragma omp parallel for collapse(3)
    for( int i=0; i<myGrids.nx; i++ )
    {
      for( int j=0; j<myGrids.ny; j++ )
      {
        for( int k=0; k<myGrids.nz; k++ )
        {
          myModels.vp[IDX3( i, j, k )]=init_vp_value;
        }
      }
    }

   #pragma omp parallel for collapse(3)
    for( int i=-myGrids.lx; i<myGrids.nx+myGrids.lx; i++ )
    {
      for( int j=-myGrids.ly; j<myGrids.ny+myGrids.ly; j++ )
      {
        for( int k=-myGrids.lz; k<myGrids.nz+myGrids.lz; k++ )
        {
          myModels.pnp1[IDX3_l( i, j, k )]=0.000001;
          myModels.pn[IDX3_l( i, j, k )]=0.000001;
          myModels.pnGlobal( IDX3_l( i, j, k ), 0 )=0.000001;
          myModels.pnGlobal( IDX3_l( i, j, k ), 1 )=0.000001;
        }
      }
    }
  }

};

#endif  //FDTDINIT_HPP
