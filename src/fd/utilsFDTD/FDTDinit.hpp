#ifndef FDTDINIT_HPP
#define FDTDINIT_HPP

#include "utils.hpp"
#include "FDTDutils.hpp"


struct FDTDInit
{
  FDTDUtils myFDTDUtils;
  SolverUtils myUtils;

  int   sourceOrder=1;
  int   ncoefs=5;

  float f0=15.;
  float fmax=2.5*f0;   
  float timeMax=1.0;
  float timeStep;
  int   nSamples;

  float vmin=1500;
  float vmax=4500;
  double coef0;

  vectorReal RHSTerm;
  vectorReal coefx;
  vectorReal coefy;
  vectorReal coefz;

  vectorReal vp;
  vectorReal phi;
  vectorReal eta;
  vectorReal pnp1;
  vectorReal pn;
  vectorReal pnm1;

  void init_geometry(int argc, char *argv[], FDTDGRIDS &myGrids)
  {
   myGrids.nx = (argc > 1)? std::stoi(argv[1]) : 150;
   myGrids.ny = myGrids.nx;
   myGrids.nz = myGrids.nx;

   myGrids.xs = myGrids.nx/2;
   myGrids.ys = myGrids.ny/2;
   myGrids.zs = myGrids.nz/2;

   myGrids.lx=4;
   myGrids.ly=4;
   myGrids.lz=4;

   myGrids.dx=10;
   myGrids.dy=10;
   myGrids.dz=10;

   float lambdamax=vmin/fmax;

   // init pml limits
   myGrids.ntaperx=3;
   myGrids.ntapery=3;
   myGrids.ntaperz=3;

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

  void init_coefficients(FDTDGRIDS &myGrids) 
  {
   coefx = allocateVector<vectorReal>(ncoefs, "coefx");
   coefy = allocateVector<vectorReal>(ncoefs, "coefy");
   coefz = allocateVector<vectorReal>(ncoefs, "coefz");

   myFDTDUtils.init_coef(myGrids.dx, coefx);
   myFDTDUtils.init_coef(myGrids.dy, coefy);
   myFDTDUtils.init_coef(myGrids.dz, coefz);

   coef0 = -2.*(coefx[1]+coefx[2]+coefx[3]+coefx[4]);
   coef0+= -2.*(coefy[1]+coefy[2]+coefy[3]+coefy[4]);
   coef0+= -2.*(coefz[1]+coefz[2]+coefz[3]+coefz[4]);

   timeStep=myFDTDUtils.compute_dt_sch(vmax,coefx,coefy,coefz);
   nSamples=timeMax/timeStep;

  }

  void init_source()
  {
   // compute source term
   RHSTerm = allocateVector<vectorReal>(nSamples, "RHSTerm");

   std::vector<float> sourceTerm=myUtils.computeSourceTerm(nSamples,timeStep,f0,sourceOrder);
   for(int i=0;i<nSamples;i++)
   {
     RHSTerm[i]=sourceTerm[i];
   }

  }

  void init_models(FDTDGRIDS &myGrids, FDTDMODELS &myModels)
  {
   int modelVolume = myGrids.nx * myGrids.ny * myGrids.nz;
   int extModelVolume = ( myGrids.nx + 2 * myGrids.lx ) * 
                        ( myGrids.ny + 2 * myGrids.ly ) * 
                        ( myGrids.nz + 2 * myGrids.lz ) ;
   int etaModelVolume = ( myGrids.nx + 2 ) * ( myGrids.ny + 2 ) * ( myGrids.nz + 2 ) ;

   vp   = allocateVector<vectorReal>( modelVolume, "vp");
   pnp1 = allocateVector<vectorReal>( extModelVolume, "pnp1");
   pn   = allocateVector<vectorReal>( extModelVolume, "pn");
   pnm1 = allocateVector<vectorReal>( extModelVolume, "pnm1" );
   phi  = allocateVector<vectorReal>( modelVolume, "phi");
   eta  = allocateVector<vectorReal>( etaModelVolume, "eta");

   myFDTDUtils.init_eta( myGrids.nx,  myGrids.ny,  myGrids.nz,
		         myGrids.ndampx,  myGrids.ndampy, myGrids.ndampz,
                         myGrids.x1,  myGrids.x2,  myGrids.x3,  myGrids.x4,  myGrids.x5,  myGrids.x6,
                         myGrids.y1,  myGrids.y2,  myGrids.y3,  myGrids.y4,  myGrids.y5,  myGrids.y6,
                         myGrids.z1,  myGrids.z2,  myGrids.z3,  myGrids.z4,  myGrids.z5,  myGrids.z6,
                         myGrids.dx,  myGrids.dy,  myGrids.dz,  timeStep,
                         vmax, eta);

   float init_vp_value = vmin*vmin*timeStep*timeStep;

   // initialize vp and pressure field
   #pragma omp parallel for collapse(3)
   for( int i=0; i<myGrids.nx;i++)
   {
      for( int j=0; j<myGrids.ny;j++)
      {
         for( int k=0; k<myGrids.nz;k++)
         {
           vp[IDX3(i,j,k)]=init_vp_value;
         }
      }
   }

   #pragma omp parallel for collapse(3)
   for( int i=-myGrids.lx; i<myGrids.nx+myGrids.lx;i++)
   {
      for( int j=-myGrids.ly; j<myGrids.ny+myGrids.ly;j++)
      {
         for( int k=-myGrids.lz; k<myGrids.nz+myGrids.lz;k++)
         {
           pnp1[IDX3_l(i,j,k)]=0.000001;
           pnm1[IDX3_l(i,j,k)]=0.000001;
           pn[IDX3_l(i,j,k)]=0.000001;
         }
      }
   }

   #if defined (USE_RAJA)
   myModels.vp   = vp.toView();
   myModels.pnp1 = pnp1.toView();
   myModels.pn   = pn.toView();
   myModels.pnm1 = pnm1.toView();
   myModels.phi  = phi.toView();
   myModels.eta  = eta.toView();
   myModels.RHSTerm  = RHSTerm.toView();
   myModels.coefx  = coefx.toView();
   myModels.coefy  = coefy.toView();
   myModels.coefz  = coefz.toView();
   #else
   myModels.vp   = vp;
   myModels.pnp1 = pnp1;
   myModels.pn   = pn;
   myModels.pnm1 = pnm1;
   myModels.phi  = phi;
   myModels.eta  = eta;
   myModels.RHSTerm  = RHSTerm;
   myModels.coefx  = coefx;
   myModels.coefy  = coefy;
   myModels.coefz  = coefz;
   #endif

 }

};

#endif  //FDTDINIT_HPP
