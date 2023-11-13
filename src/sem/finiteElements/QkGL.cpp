// C++ Code generated from Python Code:
#include <iostream>
#include <vector>
#include <cmath>
#include "QkGL.hpp"
#include "commonMacro.hpp"


namespace FE
{

QkGL::QkGL(){};

#ifdef USE_RAJA
 QkGL::~QkGL(){};
#elif defined USE_KOKKOS
KOKKOS_FUNCTION QkGL::~QkGL(){};
#else
  QkGL::~QkGL(){};
#endif

#ifdef USE_RAJA
void QkGL::gaussLobattoQuadraturePoints( int order, vectorDouble const & quadraturePoints ) const
#elif defined USE_KOKKOS
void QkGL::gaussLobattoQuadraturePoints( int order, vectorDouble const & quadraturePoints ) const
#else
void QkGL::gaussLobattoQuadraturePoints( int order, vectorDouble & quadraturePoints ) const
#endif
{
  if( order == 1 )
  {
    quadraturePoints[0]=-1.0;
    quadraturePoints[1]=1.0;
  }
  if( order == 2 )
  {
    quadraturePoints[0]=-1.0;
    quadraturePoints[1]=0.0;
    quadraturePoints[2]=1.0;
  }
  if( order == 3 )
  {
    quadraturePoints[0]=-1.0;
    quadraturePoints[1]=-0.4472136;
    quadraturePoints[2]=0.4472136;
    quadraturePoints[3]=1.0;
  }
  if( order == 4 )
  {
    quadraturePoints[0]=-1.0;
    quadraturePoints[1]=-0.65465367;
    quadraturePoints[2]=0.0;
    quadraturePoints[3]=0.65465367;
    quadraturePoints[4]=1.0;
  }
  if( order == 5 )
  {
    quadraturePoints[0]=-1.0;
    quadraturePoints[1]=-0.76505532;
    quadraturePoints[2]=-0.28523152;
    quadraturePoints[3]=0.28523152;
    quadraturePoints[4]=0.76505532;
    quadraturePoints[5]=1.0;
  }
}
#ifdef USE_RAJA
void QkGL::gaussLobattoQuadratureWeights( int order, vectorDouble const & weights ) const
#elif defined USE_KOKKOS
void QkGL::gaussLobattoQuadratureWeights( int order, vectorDouble const & weights ) const
#else
void QkGL::gaussLobattoQuadratureWeights( int order, vectorDouble & weights ) const
#endif
{
  if( order == 1 )
  {
    weights[0]=1.0;
    weights[1]=1.0;
  }
  if( order == 2 )
  {
    weights[0]=0.33333333;
    weights[1]=1.33333333;
    weights[2]= 0.33333333;
  }
  if( order == 3 )
  {
    weights[0]=0.16666667;
    weights[1]=0.83333333;
    weights[2]=0.83333333;
    weights[3]=0.16666667;
  }
  if( order == 4 )
  {
    weights[0]=0.1;
    weights[1]=0.54444444;
    weights[2]=0.71111111;
    weights[3]=0.54444444;
    weights[4]=0.1;
  }
  if( order == 5 )
  {
    weights[0]=0.06666667;
    weights[1]=0.37847496;
    weights[2]=0.55485838;
    weights[3]=0.55485838;
    weights[4]=0.37847496;
    weights[5]=0.06666667;
  }
}
std::vector<double> QkGL::shapeFunction1D( int order, double xi ) const
{
  std::vector<double> shapeFunction( order+1 );
  if( order==1 )
  {
    shapeFunction[0]=0.5*(1.0-xi);

    shapeFunction[1]=0.5*(1.0+xi);
  }
  if( order==2 )
  {
    shapeFunction[0]= -1.0*xi*(0.5 - 0.5*xi);

    shapeFunction[1]=(1.0 - 1.0*xi)*(1.0*xi + 1.0);

    shapeFunction[2]= 1.0*xi*(0.5*xi + 0.5);
  }
  if( order==3 )
  {
    shapeFunction[0]=(0.309016994374947 - 0.690983005625053*xi)*(0.5 - 0.5*xi)
                      *(-1.80901699437495*xi - 0.809016994374947);

    shapeFunction[1]=(0.5 - 1.11803398874989*xi)*(0.690983005625053 - 0.690983005625053*xi)
                      *(1.80901699437495*xi + 1.80901699437495);

    shapeFunction[2]=(1.80901699437495 - 1.80901699437495*xi)
                      *(0.690983005625053*xi + 0.690983005625053)*(1.11803398874989*xi + 0.5);

    shapeFunction[3]=(0.5*xi + 0.5)*(0.690983005625053*xi + 0.309016994374947)
                      *(1.80901699437495*xi - 0.809016994374947);
  }
  if( order==4 )
  {
    shapeFunction[0]=1.0*xi*(0.39564392373896 - 0.60435607626104*xi)*(0.5 - 0.5*xi)
                      *(-2.89564392373896*xi - 1.89564392373896);

    shapeFunction[1]=-1.52752523165195*xi*(0.5 - 0.763762615825973*xi)*(0.60435607626104 - 0.60435607626104*xi)
                      *(2.89564392373896*xi + 2.89564392373896);

    shapeFunction[2]=(1.0 - 1.52752523165195*xi)*(1.0 - 1.0*xi)*(1.0*xi + 1.0)*(1.52752523165195*xi + 1.0);

    shapeFunction[3]= 1.52752523165195*xi*(2.89564392373896 - 2.89564392373896*xi)
                      *(0.60435607626104*xi + 0.60435607626104)*(0.763762615825973*xi + 0.5);

    shapeFunction[4]= 1.0*xi*(0.5*xi + 0.5)*(0.60435607626104*xi + 0.39564392373896)
                      *(2.89564392373896*xi - 1.89564392373896);
  }
  if( order==5 )
  {
    shapeFunction[0]=(0.221930066935875 - 0.778069933064125*xi)*(0.433445520691247 - 0.566554479308753*xi)
                      *(0.5 - 0.5*xi)*(-4.25632117622354*xi - 3.25632117622354)
                      *(-1.39905441140358*xi - 0.399054411403579);

    shapeFunction[1]=(0.271574874126072 - 0.952120850728289*xi)*(0.5 - 0.6535475074298*xi)
                      *(0.566554479308753 - 0.566554479308753*xi)*(-2.0840983387567*xi - 0.594450529658367)
                      *(4.25632117622354*xi + 4.25632117622354);

    shapeFunction[2]=(0.5 - 1.75296196636787*xi)*(0.728425125873928 - 0.952120850728289*xi)
                      *(0.778069933064125 - 0.778069933064125*xi)*(1.39905441140358*xi + 1.39905441140358)
                      *(2.0840983387567*xi + 1.59445052965837);

    shapeFunction[3]=(1.39905441140358 - 1.39905441140358*xi)*(1.59445052965837 - 2.0840983387567*xi)
                      *(0.778069933064125*xi + 0.778069933064125)*(0.952120850728289*xi + 0.728425125873928)
                      *(1.75296196636787*xi + 0.5);

    shapeFunction[4]=(4.25632117622354 - 4.25632117622354*xi)*(0.566554479308753*xi + 0.566554479308753)
                      *(0.6535475074298*xi + 0.5)*(0.952120850728289*xi + 0.271574874126072)
                      *(2.0840983387567*xi - 0.594450529658367);

    shapeFunction[5]=(0.5*xi + 0.5)*(0.566554479308753*xi + 0.433445520691247)
                      *(0.778069933064125*xi + 0.221930066935875)*(1.39905441140358*xi - 0.399054411403579)
                      *(4.25632117622354*xi - 3.25632117622354);
  }
  return shapeFunction;
}
std::vector<double> QkGL::derivativeShapeFunction1D( int order, double xi ) const
{
  std::vector<double> derivativeShapeFunction( order+1 );

  if( order == 1 )
  {
    derivativeShapeFunction[0]=-0.5;
    derivativeShapeFunction[1]=0.5;
  }
  if( order == 2 )
  {
    derivativeShapeFunction[0]=1.0*xi - 0.5;
    derivativeShapeFunction[1]=-2.0*xi;
    derivativeShapeFunction[2]=1.0*xi + 0.5;
  }
  if( order == 3 )
  {
    derivativeShapeFunction[0]=-1.80901699437495*(0.309016994374947 - 0.690983005625053*xi)*(0.5 - 0.5*xi)
                                + (-1.80901699437495*xi - 0.809016994374947)*(0.345491502812526*xi - 0.345491502812526)
                                + (-1.80901699437495*xi - 0.809016994374947)*(0.345491502812526*xi - 0.154508497187474);

    derivativeShapeFunction[1]=1.80901699437495*(0.5 - 1.11803398874989*xi)*(0.690983005625053 - 0.690983005625053*xi)
                                + (0.772542485937369*xi - 0.772542485937369)*(1.80901699437495*xi + 1.80901699437495)
                                + (0.772542485937369*xi - 0.345491502812526)*(1.80901699437495*xi + 1.80901699437495);

    derivativeShapeFunction[2]=(1.80901699437495 - 1.80901699437495*xi)*(0.772542485937369*xi + 0.345491502812526) +
                                (1.80901699437495 - 1.80901699437495*xi)*(0.772542485937369*xi + 0.772542485937369) -
                                1.80901699437495*(0.690983005625053*xi + 0.690983005625053)*(1.11803398874989*xi + 0.5);

    derivativeShapeFunction[3]=(0.345491502812526*xi + 0.154508497187474)*(1.80901699437495*xi - 0.809016994374947) +
                                (0.345491502812526*xi + 0.345491502812526)*(1.80901699437495*xi - 0.809016994374947) +
                                1.80901699437495*(0.5*xi + 0.5)*(0.690983005625053*xi + 0.309016994374947);
  }
  if( order == 4 )
  {
    derivativeShapeFunction[0]=2.89564392373896*xi*(0.39564392373896 - 0.60435607626104*xi)*(0.5 - 0.5*xi) +
                                0.5*xi*(0.39564392373896 - 0.60435607626104*xi)*(-2.89564392373896*xi - 1.89564392373896)
                                + 0.60435607626104*xi*(0.5 - 0.5*xi)*(-2.89564392373896*xi - 1.89564392373896) +
                                (0.39564392373896 - 0.60435607626104*xi)*(-2.89564392373896*xi - 1.89564392373896)*(0.5*xi - 0.5);

    derivativeShapeFunction[1]=-4.42316915539091*xi*(0.5 - 0.763762615825973*xi)*(0.60435607626104 - 0.60435607626104*xi)
                                + 0.923169155390906*xi*(0.5 - 0.763762615825973*xi)*(2.89564392373896*xi + 2.89564392373896)
                                + 1.16666666666667*xi*(0.60435607626104 - 0.60435607626104*xi)*(2.89564392373896*xi + 2.89564392373896)
                                + (0.60435607626104 - 0.60435607626104*xi)*(1.16666666666667*xi - 0.763762615825973)
                                *(2.89564392373896*xi + 2.89564392373896);

    derivativeShapeFunction[2]=(1.0 - 1.52752523165195*xi)*(1.0 - 1.0*xi)*(1.52752523165195*xi + 1.0) +
                                (1.0 - 1.52752523165195*xi)*(1.0 - 1.0*xi)*(1.52752523165195*xi + 1.52752523165195)
                                - 1.0*(1.0 - 1.52752523165195*xi)*(1.0*xi + 1.0)*(1.52752523165195*xi + 1.0)
                                - 1.52752523165195*(1.0 - 1.0*xi)*(1.0*xi + 1.0)*(1.52752523165195*xi + 1.0);

    derivativeShapeFunction[3]=1.16666666666667*xi*(2.89564392373896 - 2.89564392373896*xi)*(0.60435607626104*xi + 0.60435607626104)
                                + 0.923169155390906*xi*(2.89564392373896 - 2.89564392373896*xi)*(0.763762615825973*xi + 0.5)
                                - 4.42316915539091*xi*(0.60435607626104*xi + 0.60435607626104)*(0.763762615825973*xi + 0.5)
                                + (2.89564392373896 - 2.89564392373896*xi)*(0.60435607626104*xi + 0.60435607626104)
                                *(1.16666666666667*xi + 0.763762615825973);

    derivativeShapeFunction[4]=2.89564392373896*xi*(0.5*xi + 0.5)*(0.60435607626104*xi + 0.39564392373896)
                                + 0.60435607626104*xi*(0.5*xi + 0.5)*(2.89564392373896*xi - 1.89564392373896)
                                + 0.5*xi*(0.60435607626104*xi + 0.39564392373896)*(2.89564392373896*xi - 1.89564392373896)
                                + (0.5*xi + 0.5)*(0.60435607626104*xi + 0.39564392373896)*(2.89564392373896*xi - 1.89564392373896);
  }
  if( order == 5 )
  {
    derivativeShapeFunction[0]=-1.39905441140358*(0.221930066935875 - 0.778069933064125*xi)*(0.433445520691247 - 0.566554479308753*xi)
                                *(0.5 - 0.5*xi)*(-4.25632117622354*xi - 3.25632117622354)
                                - 4.25632117622354*(0.221930066935875 - 0.778069933064125*xi)*(0.433445520691247 - 0.566554479308753*xi)
                                *(0.5 - 0.5*xi)*(-1.39905441140358*xi - 0.399054411403579) + (0.221930066935875 - 0.778069933064125*xi)
                                *(-4.25632117622354*xi - 3.25632117622354)*(-1.39905441140358*xi - 0.399054411403579)
                                *(0.283277239654376*xi - 0.283277239654376) + (0.221930066935875 - 0.778069933064125*xi)
                                *(-4.25632117622354*xi - 3.25632117622354)*(-1.39905441140358*xi - 0.399054411403579)
                                *(0.283277239654376*xi - 0.216722760345624) - 0.778069933064125*(0.433445520691247 - 0.566554479308753*xi)
                                *(0.5 - 0.5*xi)*(-4.25632117622354*xi - 3.25632117622354)*(-1.39905441140358*xi - 0.399054411403579);

    derivativeShapeFunction[1]= -2.0840983387567*(0.271574874126072 - 0.952120850728289*xi)*(0.5 - 0.6535475074298*xi)
                                *(0.566554479308753 - 0.566554479308753*xi)*(4.25632117622354*xi + 4.25632117622354)
                                - 0.566554479308753*(0.271574874126072 - 0.952120850728289*xi)*(0.5 - 0.6535475074298*xi)
                                *(-2.0840983387567*xi - 0.594450529658367)*(4.25632117622354*xi + 4.25632117622354)
                                +(0.271574874126072 - 0.952120850728289*xi)*(0.566554479308753 - 0.566554479308753*xi)
                                *(2.12816058811177 - 2.78170809554157*xi)*(-2.0840983387567*xi - 0.594450529658367)
                                +(0.271574874126072 - 0.952120850728289*xi)*(0.566554479308753 - 0.566554479308753*xi)
                                *(-2.78170809554157*xi - 2.78170809554157)*(-2.0840983387567*xi - 0.594450529658367)
                                - 0.952120850728289*(0.5 - 0.6535475074298*xi)*(0.566554479308753 - 0.566554479308753*xi)
                                *(-2.0840983387567*xi - 0.594450529658367)*(4.25632117622354*xi + 4.25632117622354);

    derivativeShapeFunction[2]= 2.0840983387567*(0.5 - 1.75296196636787*xi)*(0.728425125873928 - 0.952120850728289*xi)
                                *(0.778069933064125 - 0.778069933064125*xi)*(1.39905441140358*xi + 1.39905441140358)
                                + 1.39905441140358*(0.5 - 1.75296196636787*xi)*(0.728425125873928 - 0.952120850728289*xi)
                                *(0.778069933064125 - 0.778069933064125*xi)*(2.0840983387567*xi + 1.59445052965837)
                                - 0.952120850728289*(0.5 - 1.75296196636787*xi)*(0.778069933064125 - 0.778069933064125*xi)
                                *(1.39905441140358*xi + 1.39905441140358)*(2.0840983387567*xi + 1.59445052965837)
                                + (0.728425125873928 - 0.952120850728289*xi)*(1.3639269998358*xi - 1.3639269998358)
                                *(1.39905441140358*xi + 1.39905441140358)*(2.0840983387567*xi + 1.59445052965837)
                                + (0.728425125873928 - 0.952120850728289*xi)*(1.3639269998358*xi - 0.389034966532063)
                                *(1.39905441140358*xi + 1.39905441140358)*(2.0840983387567*xi + 1.59445052965837);

    derivativeShapeFunction[3]=0.952120850728289*(1.39905441140358 - 1.39905441140358*xi)*(1.59445052965837 - 2.0840983387567*xi)
                                *(0.778069933064125*xi + 0.778069933064125)*(1.75296196636787*xi + 0.5)
                                + (1.39905441140358 - 1.39905441140358*xi)*(1.59445052965837 - 2.0840983387567*xi)
                                *(0.952120850728289*xi + 0.728425125873928)*(1.3639269998358*xi + 0.389034966532063)
                                + (1.39905441140358 - 1.39905441140358*xi)*(1.59445052965837 - 2.0840983387567*xi)
                                *(0.952120850728289*xi + 0.728425125873928)*(1.3639269998358*xi + 1.3639269998358)
                                - 2.0840983387567*(1.39905441140358 - 1.39905441140358*xi)*(0.778069933064125*xi + 0.778069933064125)
                                *(0.952120850728289*xi + 0.728425125873928)*(1.75296196636787*xi + 0.5)
                                - 1.39905441140358*(1.59445052965837 - 2.0840983387567*xi)*(0.778069933064125*xi + 0.778069933064125)
                                *(0.952120850728289*xi + 0.728425125873928)*(1.75296196636787*xi + 0.5);

    derivativeShapeFunction[4]=(2.78170809554157 - 2.78170809554157*xi)*(0.566554479308753*xi + 0.566554479308753)
                                *(0.952120850728289*xi + 0.271574874126072)*(2.0840983387567*xi - 0.594450529658367)
                                + 2.0840983387567*(4.25632117622354 - 4.25632117622354*xi)*(0.566554479308753*xi + 0.566554479308753)
                                *(0.6535475074298*xi + 0.5)*(0.952120850728289*xi + 0.271574874126072)
                                + 0.952120850728289*(4.25632117622354 - 4.25632117622354*xi)*(0.566554479308753*xi
                                                                                              + 0.566554479308753)*(0.6535475074298*xi + 0.5)*(2.0840983387567*xi - 0.594450529658367)
                                + 0.566554479308753*(4.25632117622354 - 4.25632117622354*xi)*(0.6535475074298*xi + 0.5)
                                *(0.952120850728289*xi + 0.271574874126072)*(2.0840983387567*xi - 0.594450529658367)
                                + (-2.78170809554157*xi - 2.12816058811177)*(0.566554479308753*xi + 0.566554479308753)
                                *(0.952120850728289*xi + 0.271574874126072)*(2.0840983387567*xi - 0.594450529658367);

    derivativeShapeFunction[5]=(0.283277239654376*xi + 0.216722760345624)*(0.778069933064125*xi + 0.221930066935875)
                                *(1.39905441140358*xi - 0.399054411403579)*(4.25632117622354*xi - 3.25632117622354)
                                + (0.283277239654376*xi + 0.283277239654376)*(0.778069933064125*xi + 0.221930066935875)
                                *(1.39905441140358*xi - 0.399054411403579)*(4.25632117622354*xi - 3.25632117622354)
                                + 4.25632117622354*(0.5*xi + 0.5)*(0.566554479308753*xi + 0.433445520691247)
                                *(0.778069933064125*xi + 0.221930066935875)*(1.39905441140358*xi - 0.399054411403579)
                                + 1.39905441140358*(0.5*xi + 0.5)*(0.566554479308753*xi + 0.433445520691247)
                                *(0.778069933064125*xi + 0.221930066935875)*(4.25632117622354*xi - 3.25632117622354)
                                + 0.778069933064125*(0.5*xi + 0.5)*(0.566554479308753*xi + 0.433445520691247)
                                *(1.39905441140358*xi - 0.399054411403579)*(4.25632117622354*xi - 3.25632117622354);
  }
  return derivativeShapeFunction;
}
// get 1D basis functions @ quadrature points
// returns 2D vector basisDunction1D of dimensions nBasisFunction1D,nQuadraturePoints
#ifdef USE_RAJA
void QkGL::getBasisFunction1D( int order, vectorDouble const & quadraturePoints, arrayDouble const & basisFunction1D ) const
#elif defined USE_KOKKOS
void QkGL::getBasisFunction1D( int order, vectorDouble const & quadraturePoints, arrayDouble const & basisFunction1D ) const
#else
void QkGL::getBasisFunction1D( int order, vectorDouble & quadraturePoints, arrayDouble & basisFunction1D ) const
#endif
{
  // loop over quadrature points
  for( int i = 0; i < order+1; i++ )
  {
    std::vector<double> tmp( order+1 );
    //extract all basis functions  for current quadrature point
    tmp=shapeFunction1D( order, quadraturePoints[i] );
    for( int j=0; j<order+1; j++ )
    {
      basisFunction1D(j,i)=tmp[j];
    }
  }
}

// get derivative of 1D basis functions @ quadrature points
// returns 2D vector derivativeBasisDunction1D of dimensions nBasisFunction1D,nQuadraturePoints
#ifdef USE_RAJA
void QkGL::getDerivativeBasisFunction1D( int order, vectorDouble const & quadraturePoints, 
                                          arrayDouble const & derivativeBasisFunction1D ) const
#elif defined USE_KOKKOS
void QkGL::getDerivativeBasisFunction1D( int order, vectorDouble const & quadraturePoints, 
                                          arrayDouble const & derivativeBasisFunction1D ) const
#else
void QkGL::getDerivativeBasisFunction1D( int order, vectorDouble & quadraturePoints, 
                                          arrayDouble & derivativeBasisFunction1D ) const
#endif
{
  // loop over quadrature points
  for( int i = 0; i < order+1; i++ )
  {
    std::vector<double> tmp( order+1 );
    //extract all basis functions  for current quadrature point
    tmp=derivativeShapeFunction1D( order, quadraturePoints[i] );
    for( int j=0; j<order+1; j++ )
    {
      derivativeBasisFunction1D(j,i)=tmp[j];
    }
  }
}
// compute 2D gauss-lobatto weights
#ifdef USE_RAJA
void QkGL::getGaussLobattoWeights( vectorDouble const & quadraturePoints,
                                           vectorDouble const & weights,
                                           vectorDouble const & W )const
#elif defined USE_KOKKOS
void QkGL::getGaussLobattoWeights( vectorDouble const & quadraturePoints,
                                           vectorDouble const & weights,
                                           vectorDouble const & W )const
#else
void QkGL::getGaussLobattoWeights( vectorDouble & quadraturePoints,
                                           vectorDouble & weights,
                                           vectorDouble & W )const
#endif
{
  for( int j=0; j<quadraturePoints.size(); j++ )
  {
    for( int i=0; i<quadraturePoints.size(); i++ )
    {
      W[i+j*quadraturePoints.size()]= weights[i]*weights[j];
    }
  }
}

// returns 2D vector basisFunction2D of dimensions nBasisFunctions,nQuadraturePoints
#ifdef USE_RAJA
void QkGL::getBasisFunction2D( vectorDouble const & quadraturePoints,
                               arrayDouble const & a,
                               arrayDouble const & b,
                               arrayDouble const & c )const                                      
#elif defined USE_KOKKOS
void QkGL::getBasisFunction2D( vectorDouble const & quadraturePoints,
                               arrayDouble const & a,
                               arrayDouble const & b,
                               arrayDouble const & c )const                                      
#else
void QkGL::getBasisFunction2D( vectorDouble & quadraturePoints,
                               arrayDouble & a,
                               arrayDouble & b,
                               arrayDouble & c )const                                      
#endif
{
  for( int j = 0; j < quadraturePoints.size(); j++ )
  {
    for( int i = 0; i<quadraturePoints.size(); i++ )
    {
      for( int k = 0; k< quadraturePoints.size(); k++ )
      {
        for( int l=0; l<quadraturePoints.size(); l++ )
        {
          c(i+quadraturePoints.size()*j,l+quadraturePoints.size()*k)=a(i,l)*b(j,k);
        }
      }
    }
  }
}
// compute jacobian matrix for element to ref element coordinates
// Xi[0,..,1][0,..,nPointsPerElement], global coordinate of element
// dxPhi,dyPhi 2D derivative of basis Functions
#ifdef USE_RAJA
LVARRAY_HOST_DEVICE int QkGL::computeJacobianMatrix(  const int & nPointsPerElement,
                                                      double const  Xi[][2],
                                                      arrayDoubleView const & dxPhi,
                                                      arrayDoubleView const & dyPhi,
                                                      double  jacobianMatrix[][4] )const
#elif defined USE_KOKKOS
KOKKOS_FUNCTION int QkGL::computeJacobianMatrix( const int & nPointsPerElement,
                                                 double const Xi[][2],
                                                 arrayDouble const & dxPhi,
                                                 arrayDouble const & dyPhi,
                                                 double  jacobianMatrix[][4] )const
#else
int QkGL::computeJacobianMatrix(  const int & nPointsPerElement,
                                  double const  Xi[][2],
                                  arrayDouble & dxPhi,
                                  arrayDouble & dyPhi,
                                  double  jacobianMatrix[][4] )const
#endif
{
  for( int i=0; i<nPointsPerElement; i++ )
  {
    jacobianMatrix[i][0]=0;
    jacobianMatrix[i][1]=0;
    jacobianMatrix[i][2]=0;
    jacobianMatrix[i][3]=0;
    for( int j=0; j<nPointsPerElement; j++ )
    {
      jacobianMatrix[i][0]+=Xi[j][0]*dxPhi(j,i);
      jacobianMatrix[i][1]+=Xi[j][0]*dyPhi(j,i);
      jacobianMatrix[i][2]+=Xi[j][1]*dxPhi(j,i);
      jacobianMatrix[i][3]+=Xi[j][1]*dyPhi(j,i);
    }
  }
  return 0;
}

// compute jacobian matrix determinant
#ifdef USE_RAJA
LVARRAY_HOST_DEVICE int QkGL::computeDeterminantOfJacobianMatrix(  const int & nPointsPerElement,
                                                                   double const jacobianMatrix[][4],
                                                                   double detJ[] ) const
#elif defined USE_KOKKOS
KOKKOS_FUNCTION int QkGL::computeDeterminantOfJacobianMatrix(  const int & nPointsPerElement,
                                                               double const jacobianMatrix[][4],
                                                               double detJ[] ) const
#else
int QkGL::computeDeterminantOfJacobianMatrix( const int & nPointsPerElement,
                                              double const jacobianMatrix[][4],
                                              double detJ[]  ) const
#endif
{
  for( int i=0; i<nPointsPerElement; i++ )
  {
    detJ[i]=(jacobianMatrix[i][0]*jacobianMatrix[i][3]-jacobianMatrix[i][2]*jacobianMatrix[i][1]);
  }
  return 0;
}

// compute inverse of Jacobian Matrix
#ifdef USE_RAJA
LVARRAY_HOST_DEVICE int QkGL::computeInvJacobianMatrix( const int & nPointsPerElement,
                                                        double const  jacobianMatrix[][4],
                                                        double const  detJ[],
                                                        double invJacobianMatrix [][4]) const
#elif defined USE_KOKKOS
KOKKOS_FUNCTION int QkGL::computeInvJacobianMatrix( const int & nPointsPerElement,
                                                    double const  jacobianMatrix[][4],
                                                    double const  detJ[],
                                                    double invJacobianMatrix [][4] ) const
#else
int QkGL::computeInvJacobianMatrix( const int & nPointsPerElement,
                                    double const  jacobianMatrix[][4],
                                    double const  detJ[],
                                    double invJacobianMatrix [][4]  ) const
#endif
{
  for( int i=0; i<nPointsPerElement; i++ )
  {
    invJacobianMatrix[i][0]=(jacobianMatrix[i][3]/detJ[i]);
    invJacobianMatrix[i][1]=(-jacobianMatrix[i][1]/detJ[i]);
    invJacobianMatrix[i][2]=(-jacobianMatrix[i][2]/detJ[i]);
    invJacobianMatrix[i][3]=(jacobianMatrix[i][0]/detJ[i]);
  }
  return 0;
}

// compute inverse of Jacobian Matrix
#ifdef USE_RAJA
LVARRAY_HOST_DEVICE int QkGL::computeTranspInvJacobianMatrix( const int & nPointsPerElement,
                                                              double const  jacobianMatrix[][4],
                                                              double const detJ[],
                                                              double  transpInvJacobianMatrix[][4] ) const
#elif defined USE_KOKKOS
KOKKOS_FUNCTION int QkGL::computeTranspInvJacobianMatrix( const int & nPointsPerElement,
                                                          double const  jacobianMatrix[][4],
                                                          double const detJ[],
                                                          double transpInvJacobianMatrix[][4] ) const
#else
int QkGL::computeTranspInvJacobianMatrix( const int & nPointsPerElement,
                                          double const  jacobianMatrix[][4],
                                          double const detJ[],
                                          double  transpInvJacobianMatrix[][4] ) const
#endif
{
  for( int i=0; i<nPointsPerElement; i++ )
  {
    transpInvJacobianMatrix[i][0]=(jacobianMatrix[i][3]/detJ[i]);
    transpInvJacobianMatrix[i][1]=(-jacobianMatrix[i][2]/detJ[i]);
    transpInvJacobianMatrix[i][2]=(-jacobianMatrix[i][1]/detJ[i]);
    transpInvJacobianMatrix[i][3]=(jacobianMatrix[i][0]/detJ[i]);
  }
  return 0;
}

// compute B the matrix containing the geometrical informations
#ifdef USE_RAJA
LVARRAY_HOST_DEVICE int QkGL::computeB( const int & nPointsPerElement,
                                        double const invJacobianMatrix[][4],
                                        double const transpInvJacobianMatrix[][4],
                                        double const detJ[],
                                        double   B[][4] ) const
#elif defined USE_KOKKOS
KOKKOS_FUNCTION int QkGL::computeB( const int & nPointsPerElement,
                                    double const invJacobianMatrix[][4],
                                    double const transpInvJacobianMatrix[][4],
                                    double const detJ[],
                                    double   B[][4]) const
#else
int QkGL::computeB( const int & nPointsPerElement,
                    double const invJacobianMatrix[][4],
                    double const transpInvJacobianMatrix[][4],
                    double const detJ[],
                    double   B[][4]) const
#endif
{
  for( int i=0; i<nPointsPerElement; i++ )
  {
    B[i][0]=(abs( detJ[i] )*(invJacobianMatrix[i][0]*transpInvJacobianMatrix[i][0]+
                             invJacobianMatrix[i][1]*transpInvJacobianMatrix[i][2]));
    B[i][1]=(abs( detJ[i] )*(invJacobianMatrix[i][0]*transpInvJacobianMatrix[i][1]+
                             invJacobianMatrix[i][1]*transpInvJacobianMatrix[i][3]));
    B[i][2]=(abs( detJ[i] )*(invJacobianMatrix[i][2]*transpInvJacobianMatrix[i][0]+
                             invJacobianMatrix[i][3]*transpInvJacobianMatrix[i][2]));
    B[i][3]=(abs( detJ[i] )*(invJacobianMatrix[i][2]*transpInvJacobianMatrix[i][1]+
                             invJacobianMatrix[i][3]*transpInvJacobianMatrix[i][3]));
  }
  return 0;
}

// compute the matrix $R_{i,j}=\int_{K}{\nabla{\phi_i}.\nabla{\phi_j}dx}$
// Marc Durufle Formulae
#ifdef USE_RAJA
LVARRAY_HOST_DEVICE int QkGL::gradPhiGradPhi( const int & nPointsPerElement,
                                              const int & order,
                                              vectorDoubleView const & weights2D,
                                              double const  B[][4],
                                              arrayDoubleView const & dPhi,
                                              double   R[][36] ) const
#elif defined USE_KOKKOS
KOKKOS_FUNCTION int QkGL::gradPhiGradPhi( const int & nPointsPerElement,
                                          const int & order,
                                          vectorDouble const & weights2D,
                                          double const  B[][4],
                                          arrayDouble const & dPhi,
                                          double   R[][36] ) const
#else
int QkGL::gradPhiGradPhi( const int & nPointsPerElement,
                          const int & order,
                          vectorDouble  & weights2D,
                          double const  B[][4],
                          arrayDouble & dPhi,
                          double   R[][36] ) const
#endif
{
  for (int i=0;i<nPointsPerElement;i++)
  {
      for (int j=0; j<nPointsPerElement;j++)
      {
        R[j][i]=0;
      }
  }
  // B11
  for( int i1=0; i1<order+1; i1++ )
  {
    for( int i2=0; i2<order+1; i2++ )
    {
      int i=i1+i2*(order+1);
      for( int j1=0; j1<order+1; j1++ )
      {
        int j=j1+i2*(order+1);
        for( int m=0; m<order+1; m++ )
        {
          R[j][i]+=weights2D[m+i2*(order+1)]*(B[m+i2*(order+1)][0]*dPhi(i1,m)*dPhi(j1,m));
        }
      }
    }
  }
  // B21
  for( int i1=0; i1<order+1; i1++ )
  {
    for( int i2=0; i2<order+1; i2++ )
    {
      int i=i1+i2*(order+1);
      for( int j1=0; j1<order+1; j1++ )
      {
        for( int j2=0; j2<order+1; j2++ )
        {
          int j=j1+j2*(order+1);
          R[j][i]+=weights2D[i1+j2*(order+1)]*(B[i1+j2*(order+1)][1]*dPhi(i2,j2)*dPhi(j1,i1));
        }
      }
    }
  }
  // B12
  for( int i1=0; i1<order+1; i1++ )
  {
    for( int i2=0; i2<order+1; i2++ )
    {
      int i=i1+i2*(order+1);
      for( int j1=0; j1<order+1; j1++ )
      {
        for( int j2=0; j2<order+1; j2++ )
        {
          int j=j1+j2*(order+1);
          R[j][i]+=weights2D[i2+j1*(order+1)]*(B[i2+j1*(order+1)][2]*dPhi(i1,j1)*dPhi(j2,i2));
        }
      }
    }
  }
  // B22
  for( int i1=0; i1<order+1; i1++ )
  {
    for( int i2=0; i2<order+1; i2++ )
    {
      int i=i1+i2*(order+1);
      for( int j2=0; j2<order+1; j2++ )
      {
        int j=i1+j2*(order+1);
        for( int n=0; n<order+1; n++ )
        {
          R[j][i]+=weights2D[i1+n*(order+1)]*(B[i1+n*(order+1)][3]*dPhi(i2,n)*dPhi(j2,n));
        }
      }
    }
  }
  return 0;
}
// compute the matrix $R_{i,j}=\int_{K}{\nabla{\phi_i}.\nabla{\phi_j}dx}$
#ifdef USE_RAJA
LVARRAY_HOST_DEVICE int QkGL::gradPhiGradPhi( const int & nPointsPerElement,
                                              vectorDoubleView const & weights2D,
                                              double const  B[][4],
                                              arrayDoubleView const & dxPhi,
                                              arrayDoubleView const & dyPhi,
                                              double  R[][36] ) const
#elif defined USE_KOKKOS
KOKKOS_FUNCTION int QkGL::gradPhiGradPhi( const int & nPointsPerElement,
                                          vectorDouble const & weights2D,
                                          double const  B[][4],
                                          arrayDouble const & dxPhi,
                                          arrayDouble const & dyPhi,
                                          double  R[][36] ) const
#else
int QkGL::gradPhiGradPhi( const int & nPointsPerElement,
                          vectorDouble const & weights2D,
                          double const  B[][4],
                          arrayDouble  & dxPhi,
                          arrayDouble  & dyPhi,
                          double  R[][36]) const
#endif
{
  for( int i=0; i<nPointsPerElement; i++ )
  {
    for( int j=0; j<nPointsPerElement; j++ )
    {
      double tmp=0;
      for( int r=0; r<nPointsPerElement; r++ )
      {
        tmp+=weights2D[r]*(B[r][0]*dxPhi(i,r)*dxPhi(j,r)+
                           B[r][1]*dxPhi(i,r)*dyPhi(j,r)+
                           B[r][2]*dyPhi(i,r)*dxPhi(j,r)+
                           B[r][3]*dyPhi(i,r)*dyPhi(j,r));

      }
      R[j][i]=tmp;
    }
  }
  return 0;
}
// compute the matrix $M_{i,j}=\int_{K}{{\phi_i}.{\phi_j}dx}$ (optimized formulation)
#ifdef USE_RAJA
LVARRAY_HOST_DEVICE int QkGL::phiIphiJ( const int & nPointsPerElement,
                                        vectorDoubleView const & weights2D,
                                        double const  detJ[],
                                        double massMatrixLocal[] ) const
#elif defined USE_KOKKOS
KOKKOS_FUNCTION int QkGL::phiIphiJ( const int & nPointsPerElement,
                                    vectorDouble const & weights2D,
                                    double const  detJ[],
                                    double massMatrixLocal[] ) const
#else
int QkGL::phiIphiJ( const int & nPointsPerElement,
                    vectorDouble & weights2D,
                    double const  detJ[],
                    double massMatrixLocal[]) const
#endif
{
  for( int i=0; i<nPointsPerElement; i++ )
  {
    massMatrixLocal[i]=weights2D[i]*abs( detJ[i] );
  }
  return 0;
}
//computeDs
#ifdef USE_RAJA
LVARRAY_HOST_DEVICE int QkGL::computeDs(  const int & iFace,
                                          const int & order,
                                          arrayIntView const & faceInfos,
                                          int  numOfBasisFunctionOnFace[],
                                          float  Js[][6],
                                          arrayRealView   const & globalNodesCoords,
                                          arrayDoubleView const & derivativeBasisFunction2DX,
                                          arrayDoubleView const & derivativeBasisFunction2DY,
                                          float  ds[]  ) const
#elif defined USE_KOKKOS
KOKKOS_FUNCTION int QkGL::computeDs(  const int & iFace,
                                      const int & order,
                                      arrayInt  const & faceInfos,
                                      int  numOfBasisFunctionOnFace[],
                                      float  Js[][6],
                                      arrayReal   const & globalNodesCoords,
                                      arrayDouble const & derivativeBasisFunction2DX,
                                      arrayDouble const & derivativeBasisFunction2DY,
                                      float  ds[] ) const
#else
int QkGL::computeDs(  const int & iFace,
                      const int & order,
                      arrayInt & faceInfos,
                      int  numOfBasisFunctionOnFace[],
                      float  Js[][6],
                      arrayReal    & globalNodesCoords,
                      arrayDouble  & derivativeBasisFunction2DX,
                      arrayDouble  & derivativeBasisFunction2DY,
                      float  ds[] ) const
#endif
{
  int face=faceInfos(iFace,1);
  // get basis functions on Boundary faces
  switch( face )
  {
    case 0:     // left
      for( int i=0; i<order+1; i++ )
      {
        numOfBasisFunctionOnFace[i]=i*(order+1);
      }
      break;
    case 1:     // bottom
      for( int i=0; i<order+1; i++ )
      {
        numOfBasisFunctionOnFace[i]=i;
      }
      break;
    case 2:         //right
      for( int i=0; i<order+1; i++ )
      {
        numOfBasisFunctionOnFace[i]=order+i*(order+1);
      }
      break;
    case 3:         //top
      for( int i=0; i<order+1; i++ )
      {
        numOfBasisFunctionOnFace[i]=i+order*(order+1);
      }
      break;
    default:
      //cout<<"error in element flag, should be set to: 0, 1, 2, 3"<<endl;
      break;
  }
  // compute ds
  for( int j=0; j<order+1; j++ )
  {
    Js[0][j]=0;    // x
    Js[1][j]=0;    // y
    for( int i=0; i<order+1; i++ )
    {
      float xi=globalNodesCoords(faceInfos(iFace,2+i),0);
      float yi=globalNodesCoords(faceInfos(iFace,2+i),1);
      if( face==0 || face==2 )
      {
        Js[0][j]+=derivativeBasisFunction2DY(numOfBasisFunctionOnFace[i],numOfBasisFunctionOnFace[j])*xi;
        Js[1][j]+=derivativeBasisFunction2DY(numOfBasisFunctionOnFace[i],numOfBasisFunctionOnFace[j])*yi;
      }
      if( face==1 || face==3 )
      {
        Js[0][j]+=derivativeBasisFunction2DX(numOfBasisFunctionOnFace[i],numOfBasisFunctionOnFace[j])*xi;
        Js[1][j]+=derivativeBasisFunction2DX(numOfBasisFunctionOnFace[i],numOfBasisFunctionOnFace[j])*yi;
      }
    }
    ds[j]=sqrt( Js[0][j]*Js[0][j]+Js[1][j]*Js[1][j] );
    //cout<<"j="<<j<<", ds="<<ds[j]<<", "<<Js[0][j]<<", "<<Js[1][j]<<endl;
  }
  return 0;
}
} // namespace FE
