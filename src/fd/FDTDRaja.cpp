//************************************************************************
// SEM proxy application v.0.0.1
//
// main.cpp: this main file is simply a driver
//************************************************************************

#include<iostream>
#include<cstdio>
#include<fstream>
#include<cmath>
#include<chrono>
#include"Array.hpp"
#include"RAJA/RAJA.hpp"
#include"Macros.hpp"
#include"ChaiBuffer.hpp"
#include<omp.h>
#include<vector>
#include"utils.hpp"
#include"dataType.hpp"
#include"FDTDutils.hpp"


//innerpoints
int inner3D(const int nx, const int ny, const int nz,
       const int x3, const int x4,
       const int y3, const int y4,
       const int z3, const int z4,
       const int lx, const int ly, const int lz,
       const float coef0,
       vectorRealView const & coefx,
       vectorRealView const & coefy,
       vectorRealView const & coefz,
       vectorRealView const & vp,
       vectorRealView const & pnp1,
       vectorRealView const & pn,
       vectorRealView const & pnm1)
{
   RAJA::TypedRangeSegment<int> KRange(z3, z4);
   RAJA::TypedRangeSegment<int> JRange(y3, y4);
   RAJA::TypedRangeSegment<int> IRange(x3, x4);

   using EXEC_POL =
   RAJA::KernelPolicy<
     RAJA::statement::CudaKernel<
       RAJA::statement::For<2, RAJA::cuda_thread_x_loop,      // k
         RAJA::statement::For<1, RAJA::cuda_thread_y_loop,    // j
           RAJA::statement::For<0, RAJA::cuda_thread_z_loop,  // i
             RAJA::statement::Lambda<0>
           >
         >
       >
     >
   >;

   RAJA::kernel<EXEC_POL>( RAJA::make_tuple(IRange, JRange, KRange), [=] __device__ (int i, int j, int k) {
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
    pnp1[IDX3_l(i,j,k)]=2.*pn[IDX3_l(i,j,k)]-pnm1[IDX3_l(i,j,k)]
                       +vp[IDX3(i,j,k)]*(coef0*pn[IDX3_l(i,j,k)]+lapx+lapy+lapz);
   });
   return (0);
}

int pml3D(const int nx, const int ny, const int nz,
           const int x3, const int x4, 
           const int y3, const int y4, 
           const int z3, const int z4,
           const int lx, const int ly, const int lz,
           const float coef0,
           const float hdx_2, const float hdy_2, const float hdz_2,
           vectorRealView const & coefx,
           vectorRealView const & coefy,
           vectorRealView const & coefz,
           vectorRealView const & vp,
           vectorRealView const & phi,
           vectorRealView const & eta,
           vectorRealView const & pnp1,
           vectorRealView const & pn,
           vectorRealView const & pnm1)
{
   RAJA::TypedRangeSegment<int> KRange(z3, z4);
   RAJA::TypedRangeSegment<int> JRange(y3, y4);
   RAJA::TypedRangeSegment<int> IRange(x3, x4);

   using EXEC_POL =
   RAJA::KernelPolicy<
     RAJA::statement::CudaKernel<
       RAJA::statement::For<2, RAJA::cuda_thread_x_loop,      // k
         RAJA::statement::For<1, RAJA::cuda_thread_y_loop,    // j
           RAJA::statement::For<0, RAJA::cuda_thread_z_loop,  // i
             RAJA::statement::Lambda<0>
           >
         >
       >
     >
   >;

   RAJA::kernel<EXEC_POL>( RAJA::make_tuple(IRange, JRange, KRange), [=] __device__ (int i, int j, int k)
   {
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

     pnp1[IDX3_l(i,j,k)]=((2.-eta[IDX3_eta1(i,j,k)]*eta[IDX3_eta1(i,j,k)]
                +2.*eta[IDX3_eta1(i,j,k)])*pn[IDX3_l(i,j,k)]
                -pnm1[IDX3_l(i,j,k)]
                +vp[IDX3(i,j,k)]*(lap+phi[IDX3(i,j,k)]))/(1.+2.*eta[IDX3_eta1(i,j,k)]);

     phi[IDX3(i,j,k)]=(phi[IDX3(i,j,k)]-((eta[IDX3_eta1(i+1,j,k)]-eta[IDX3_eta1(i-1,j,k)])
                *(pn[IDX3_l(i+1,j,k)]-pn[IDX3_l(i-1,j,k)])*hdx_2
                +(eta[IDX3_eta1(i,j+1,k)]-eta[IDX3_eta1(i,j-1,k)])
                *(pn[IDX3_l(i,j+1,k)]-pn[IDX3_l(i,j-1,k)])*hdy_2
                +(eta[IDX3_eta1(i,j,k+1)]-eta[IDX3_eta1(i,j,k-1)])
                *(pn[IDX3_l(i,j,k+1)]-pn[IDX3_l(i,j,k-1)])*hdz_2))
                /(1.+eta[IDX3_eta1(i,j,k)]);
   });
   return(0);
}

int main( int argc, char *argv[] )
{
   constexpr int nx=150;
   constexpr int ny=150;
   constexpr int nz=150;
   constexpr int lx=4;
   constexpr int ly=4;
   constexpr int lz=4;
   constexpr float dx=10;
   constexpr float dy=10;
   constexpr float dz=10;

   constexpr int   sourceOrder=1;
   constexpr int   xs=nx/2;
   constexpr int   ys=ny/2;
   constexpr int   zs=nz/2;
   constexpr float f0=15.;
   constexpr float fmax=2.5*f0;
   constexpr float timeMax=0.8;

   constexpr int ncoefs=5;
   constexpr float vmax=1500;
   
   // imports utils
   solverUtils myUtils;
   FDTDUtils myFDTDUtils;

   // init pml limits
   constexpr int ntaperx=3;
   constexpr int ntapery=3;
   constexpr int ntaperz=3;
   constexpr float hdx_2=1./(4.*dx*dx);
   constexpr float hdy_2=1./(4.*dy*dy);
   constexpr float hdz_2=1./(4.*dz*dz);
   constexpr float lambdamax=vmax/fmax;
   constexpr int ndampx=ntaperx*lambdamax/dx;
   constexpr int ndampy=ntapery*lambdamax/dy;
   constexpr int ndampz=ntaperz*lambdamax/dz;
   printf("ndampx=%d ndampy=%d ndampz=%d\n",ndampx, ndampy,ndampz);
   constexpr int x1=0;
   constexpr int x2=ndampx;
   constexpr int x3=ndampx;
   constexpr int x4=nx-ndampx;
   constexpr int x5=nx-ndampx;
   constexpr int x6=nx;
   constexpr int y1=0;
   constexpr int y2=ndampy;
   constexpr int y3=ndampy;
   constexpr int y4=ny-ndampy;
   constexpr int y5=ny-ndampy;
   constexpr int y6=ny;
   constexpr int z1=0;
   constexpr int z2=ndampz;
   constexpr int z3=ndampz;
   constexpr int z4=nz-ndampz;
   constexpr int z5=nz-ndampz;
   constexpr int z6=nz;

   // allocate vector and arrays 
   // FD coefs
   vectorReal h_coefx=allocateVector<vectorReal>(ncoefs);
   vectorReal h_coefy=allocateVector<vectorReal>(ncoefs);
   vectorReal h_coefz=allocateVector<vectorReal>(ncoefs);
   // model
   vectorReal h_vp=allocateVector<vectorReal>(nx*ny*nz);
   // pressure fields
   vectorReal h_pnp1=allocateVector<vectorReal>((nx+2*lx)*(ny+2*ly)*(nz+2*lz));
   vectorReal h_pn=allocateVector<vectorReal>((nx+2*lx)*(ny+2*ly)*(nz+2*lz));
   vectorReal h_pnm1=allocateVector<vectorReal>((nx+2*lx)*(ny+2*ly)*(nz+2*lz));
   // PML arrays
   vectorReal h_phi=allocateVector<vectorReal>(nx*ny*nz);
   vectorReal h_eta=allocateVector<vectorReal>((nx+2)*(ny+2)*(nz+2));

   // extract FD coefs
   myFDTDUtils.init_coef(dx, h_coefx);
   myFDTDUtils.init_coef(dy, h_coefy);
   myFDTDUtils.init_coef(dz, h_coefz);

   double coef0=-2.*(h_coefx[1]+h_coefx[2]+h_coefx[3]+h_coefx[4]);
   coef0+=-2.*(h_coefy[1]+h_coefy[2]+h_coefy[3]+h_coefy[4]);
   coef0+=-2.*(h_coefz[1]+h_coefz[2]+h_coefz[3]+h_coefz[4]);

   // compute time step
   float timeStep=myFDTDUtils.compute_dt_sch(vmax,h_coefx,h_coefy,h_coefz);
   float timeStep2=timeStep*timeStep;
   const int nSamples=timeMax/timeStep;
   printf("timeStep=%f\n",timeStep);

   // compute source term
   // source term
   vectorReal h_RHSTerm=allocateVector<vectorReal>(nSamples);
   std::vector<float>sourceTerm=myUtils.computeSourceTerm(nSamples,timeStep,f0,sourceOrder);
   for(int i=0;i<nSamples;i++)
   {
     h_RHSTerm[i]=sourceTerm[i];
   }

   printf("memory used for vectra and arrays %d bytes\n",
          (nx*ny*nz+3*(nx+2*lx)*(ny+2*ly)*(nz+2*lz)+nSamples+ncoefs)*4);

   // initialize vp and pressure field
   #pragma omp parallel for collapse(3)
   for( int i=0; i<nx;i++)
   {
      for( int j=0; j<ny;j++)
      {
         for( int k=0; k<nz;k++)
         {
           h_vp[IDX3(i,j,k)]=1500.*1500.*timeStep2;
           h_phi[IDX3(i,j,k)]=0.;
         }
      }
   }
   #pragma omp parallel for collapse(3)
   for( int i=-lx; i<nx+lx;i++)
   {
      for( int j=-ly; j<ny+ly;j++)
      {
         for( int k=-lz; k<nz+lz;k++)
         {
           h_pnp1[IDX3_l(i,j,k)]=0.000001;
           h_pn[IDX3_l(i,j,k)]  =0.000001;
           h_pnm1[IDX3_l(i,j,k)]=0.000001;
         }
      }
   }

   //init pml eta array
   myFDTDUtils.init_eta( nx,  ny,  nz,
		         ndampx,  ndampy, ndampz,
                         x1,  x2,  x3,  x4,  x5,  x6,
                         y1,  y2,  y3,  y4,  y5,  y6,
                         z1,  z2,  z3,  z4,  z5,  z6,
                         dx,  dy,  dz,  timeStep,
                         vmax, h_eta);

   char filename_buf[32];
   snprintf(filename_buf, sizeof(filename_buf), "debug_eta.H@");
   printf("\n");
   printf("eta file n1=%d n2=%d\n",nx,nz);
   printf("\n");
   FILE *dbg = fopen(filename_buf, "wb");
   for (int k = 0; k < nz; ++k) {
       for (int j = ny/2; j < ny/2+1; ++j) {
           for (int i = 0; i < nx; ++i) {
               fwrite(&h_eta[IDX3_eta1(i,j,k)], sizeof(float),1, dbg);
           }
        }
   }
   /* Clean up */
   fclose(dbg);


   // define policy for source term
   RAJA::TypedRangeSegment<int> KRanges(xs, xs+1);
   RAJA::TypedRangeSegment<int> JRanges(ys, ys+1);
   RAJA::TypedRangeSegment<int> IRanges(zs, zs+1);
   
   // define policy for swapping 
   RAJA::TypedRangeSegment<int> KRange(0, nx);
   RAJA::TypedRangeSegment<int> JRange(0, ny);
   RAJA::TypedRangeSegment<int> IRange(0, ny);

   using EXEC_POL =
   RAJA::KernelPolicy<
     RAJA::statement::CudaKernel<
       RAJA::statement::For<2, RAJA::cuda_thread_x_loop,      // i
         RAJA::statement::For<1, RAJA::cuda_thread_y_loop,    // j
           RAJA::statement::For<0, RAJA::cuda_thread_z_loop,  // k
             RAJA::statement::Lambda<0>
           >
         >
       >
     >
   >;

   vectorRealView const RHSTerm=h_RHSTerm.toView();
   vectorRealView const coefx=h_coefx.toView();
   vectorRealView const coefy=h_coefy.toView();
   vectorRealView const coefz=h_coefz.toView();
   vectorRealView const vp=h_vp.toView();
   vectorRealView const phi=h_phi.toView();
   vectorRealView const eta=h_eta.toView();
   vectorRealView  pnp1=h_pnp1.toView();
   vectorRealView  pn=h_pn.toView();
   vectorRealView  pnm1=h_pnm1.toView();
   std::chrono::time_point<std::chrono::system_clock> start, end;
   start = std::chrono::system_clock::now();
   for (int itSample=0; itSample<nSamples;itSample++)
   {

     RAJA::kernel<EXEC_POL>( RAJA::make_tuple(IRanges, JRanges, KRanges), [=] __device__ (int i, int j, int k)
     {
       pn[IDX3_l(i,j,k)]+=vp[IDX3(i,j,k)]*RHSTerm[itSample];
     });
     //
     //up
     pml3D(nx,ny,nz,0,nx,0,ny,z1,z2,lx,ly,lz,coef0,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
     //front
     pml3D(nx,ny,nz,0,nx,y1,y2,z3,z4,lx,ly,lz,coef0,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
     //left
     pml3D(nx,ny,nz,x1,x2,y3,y4,z3,z4,lx,ly,lz,coef0,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
     //inner points
     inner3D(nx,ny,nz,x3,x4,y3,y4,z3,z4,lx,ly,lz,coef0,coefx,coefy,coefz,vp,pnp1,pn,pnm1);
     //right
     pml3D(nx,ny,nz,x5,x6,y3,y4,z3,z4,lx,ly,lz,coef0,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
     //back
     pml3D(nx,ny,nz,0,nx,y5,y6,z3,z4,lx,ly,lz,coef0,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);
     // bottom
     pml3D(nx,ny,nz,0,nx,0,ny,z5,z6,lx,ly,lz,coef0,hdx_2,hdy_2,hdz_2,coefx,coefy,coefz,vp,phi,eta,pnp1,pn,pnm1);

     if(itSample%50==0)
     {
	RAJA::forall<RAJA::loop_exec>(RAJA::RangeSegment(0,nx), [pnp1,phi] ( int i)
                         {});
	printf("result1 %f\n",pnp1[IDX3_l(xs,ys,zs)]);
	char filename_buf[32];
        snprintf(filename_buf, sizeof(filename_buf), "snapshot.it%d.H@", itSample);
        printf("snapshot file nx=%d nz=%d\n",nx,nz);
        FILE *snapshot_file = fopen(filename_buf, "wb");
        for (int k = 0; k < nz; ++k) {
            for (int j = ny/2; j < ny/2+1; ++j) {
                for (int i = 0; i < nx; ++i) {
                    fwrite(&pnp1[IDX3_l(i,j,k)], sizeof(float),1, snapshot_file);
                }
            }
        }
        /* Clean up */
        fclose(snapshot_file);

	snprintf(filename_buf, sizeof(filename_buf), "snapshotPhi.it%d.H@", itSample);
        printf("snapshotPhi file size n1=%d n2=%d\n",nx,nz);
        FILE *snapshotPhi_file = fopen(filename_buf, "wb");
        for (int k = 0; k < nz; ++k) {
            for (int j = ny/2; j < ny/2+1; ++j) {
                for (int i = 0; i < nx; ++i) {
                    fwrite(&phi[IDX3(i,j,k)], sizeof(float),1, snapshotPhi_file);
                }
            }
        }
        /* Clean up */
        fclose(snapshot_file);
     }
     // swap
     RAJA::kernel<EXEC_POL>( RAJA::make_tuple(IRange, JRange, KRange), [=] __device__ (int i, int j, int k)
     {
        pnm1[IDX3_l(i,j,k)]=pn[IDX3_l(i,j,k)];
        pn[IDX3_l(i,j,k)]=pnp1[IDX3_l(i,j,k)];
     });
   }
   end = std::chrono::system_clock::now();
   std::chrono::duration<double> elapsed_seconds = end - start;
   std::time_t start_time = std::chrono::system_clock::to_time_t(start);
   std::time_t end_time = std::chrono::system_clock::to_time_t(end);

   std::cout << "started computation at " << std::ctime(&start_time)<<std::endl
             << "finished computation at " << std::ctime(&end_time)<<std::endl
             << "elapsed time: " << elapsed_seconds.count() << "s\n";

   return (0);
}
