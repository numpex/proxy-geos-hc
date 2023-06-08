#ifndef UTILS_HPP_
#define UTILS_HPP_

#include  <cmath>
#include "dataType.hpp"

using namespace std;

struct solverUtils
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
        cout<<"This option is not supported yet, rickerOrder must be 0, 1 or 2"<<endl;
        break;
    }

    return pulse;
  }

  vector< float > computeSourceTerm( const int nSamples, const float timeSample, const float f0, const int order )
  {
    vector< float > sourceTerm( nSamples );
    for( int i=0; i<nSamples; i++ )
    {
      float time_n=i*timeSample;
      sourceTerm[i]=evaluateRicker( time_n, f0, order );
    }
    return sourceTerm;
  }

  void saveSnapShot( const int indexTimeStep, const int i1, arrayReal pnGlobal, simpleMesh mesh )
  {

    int numberOfNodes=mesh.getNumberOfNodes();
    vectorReal inputVector( numberOfNodes );
    int nx=mesh.getNx();
    int ny=mesh.getNy();
    float dx=mesh.getDx();
    float dy=mesh.getDy();
    for( int i = 0; i< numberOfNodes; i++ )
    {
      inputVector[i]=pnGlobal[i][i1];
    }
    arrayReal grid=mesh.projectToGrid( numberOfNodes, inputVector );
    fstream snapFile;
    string snapNumber = "snapshot"+to_string( indexTimeStep );
    snapFile.open( snapNumber, ios::out| ios::trunc );
    cout<<"nx="<<nx<<" ny="<<ny<<endl;
    for( int i=0; i<nx; i++ )
    {
      for( int j=0; j<ny; j++ )
      {
        snapFile<<i*dx<<" "<<j*dy<<" " <<grid[i][j]<<endl;
      }
    }
    snapFile.close();
  }
};
#endif //UTILS_HPP_
