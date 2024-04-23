#ifndef UTILS_HPP_
#define UTILS_HPP_

#include "dataType.hpp"
#include "../sem/mesh/simpleMesh.hpp"

using namespace grid;
using namespace std::chrono;

struct SolverUtils
{

  float evaluateRicker( float const & time_n, float const & f0, int order )
  {
    float const o_tpeak = 1.0/f0;
    float pulse = 0.0;
    if((time_n <= -0.9*o_tpeak) || (time_n >= 2.9*o_tpeak))
    {
      return pulse;
    }

    constexpr float pi = M_PI;
    float const lam = (f0*pi)*(f0*pi);

    switch( order )
    {
      case 2:
      {
        pulse = 2.0*lam*(2.0*lam*(time_n-o_tpeak)*(time_n-o_tpeak)-1.0)*exp( -lam*(time_n-o_tpeak)*(time_n-o_tpeak));
      }
      break;
      case 1:
      {
        pulse = -2.0*lam*(time_n-o_tpeak)*exp( -lam*(time_n-o_tpeak)*(time_n-o_tpeak));
      }
      break;
      case 0:
      {
        pulse = -(time_n-o_tpeak)*exp( -2*lam*(time_n-o_tpeak)*(time_n-o_tpeak) );
      }
      break;
      default:
      std::cout<<"This option is not supported yet, rickerOrder must be 0, 1 or 2"<<std::endl;
        break;
    }

    return pulse;
  }

  std::vector< float > computeSourceTerm( const int nSamples,const float timeSample,const float f0,const int order )
  {
    std::vector< float > sourceTerm( nSamples );
    for( int i=0; i<nSamples; i++ )
    {
      float time_n=i*timeSample;
      sourceTerm[i]=evaluateRicker( time_n, f0, order );
    }
    return sourceTerm;
  }


  void saveSnapShot( const int indexTimeStep, const int i1, arrayRealView pnGlobal, simpleMesh mesh )
  {

    int nx=mesh.getNx();
    int ny=mesh.getNy();
    int nz=mesh.getNz();
    float dx=mesh.getDx();
    float dy=mesh.getDy();
    float dz=mesh.getDz();
    int numberOfNodes=nx*nz;
    std::vector<float> inputVector( numberOfNodes );
    int offset=(ny==1?offset=0:offset=nx*nz*(ny/2-1));
    printf(" nx, ny/2-1, nz %d %d %d offset=%d\n",nx,ny/2-1,nz, offset);
    for( int i = offset; i< offset+numberOfNodes; i++ )
    {
      inputVector[i-offset]=pnGlobal(i,i1);
    }
    std::vector<std::vector<float>> grid=mesh.projectToGrid( numberOfNodes, inputVector );
    fstream snapFile;
    string snapNumber = "snapshot"+to_string( indexTimeStep );
    snapFile.open( snapNumber, ios::out| ios::trunc );
    std::cout<<"nx="<<nx<<" ny="<<ny<<std::endl;
    for( int i=0; i<nx; i++ )
    {
      for( int j=0; j<nz; j++ )
      {
        snapFile<<i*dx<<" "<<j*dz<<" " <<grid[i][j]<<endl;
      }
    }
    snapFile.close();
  }


};
#endif //UTILS_HPP_
