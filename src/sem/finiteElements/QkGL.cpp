// C++ Code generated from Python Code:
#include <iostream>
#include <vector>
#include <cmath>
#include "QkGL.hpp"
#include "commonMacro.hpp"


namespace FE
{

QkGL::QkGL(){};
QkGL::~QkGL(){};

//vectorDouble QkGL::gaussLobattoQuadraturePoints( int order ) const
#ifdef SEM_USE_RAJA
void QkGL::gaussLobattoQuadraturePoints( int order, vectorDouble const & quadraturePoints ) const
#elif defined SEM_USE_KOKKOS
void QkGL::gaussLobattoQuadraturePoints( int order, vectorDouble const & quadraturePoints ) const
#else
void QkGL::gaussLobattoQuadraturePoints( int order, vectorDouble & quadraturePoints ) const
#endif
{
  //vectorDouble quadraturePoints( order+1 );
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
  //return quadraturePoints;
}
//vectorDouble QkGL::gaussLobattoQuadratureWeights( int order ) const
#ifdef SEM_USE_RAJA
void QkGL::gaussLobattoQuadratureWeights( int order, vectorDouble const & weights ) const
#elif defined SEM_USE_KOKKOS
void QkGL::gaussLobattoQuadratureWeights( int order, vectorDouble const & weights ) const
#else
void QkGL::gaussLobattoQuadratureWeights( int order, vectorDouble & weights ) const
#endif
{
  //vectorDouble weights( order+1 );

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
  //return weights;
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
//arrayDouble QkGL::getBasisFunction1D( int order, vectorDouble & quadraturePoints ) const
#ifdef SEM_USE_RAJA
void QkGL::getBasisFunction1D( int order, vectorDouble const & quadraturePoints, arrayDouble const & basisFunction1D ) const
#elif defined SEM_USE_KOKKOS
void QkGL::getBasisFunction1D( int order, vectorDouble const & quadraturePoints, arrayDouble const & basisFunction1D ) const
#else
void QkGL::getBasisFunction1D( int order, vectorDouble & quadraturePoints, arrayDouble & basisFunction1D ) const
#endif
{
  //int nBasisFunction1D=order+1;
  //arrayDouble basisFunction1D( nBasisFunction1D, nBasisFunction1D );
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
  //return basisFunction1D;
}

// get derivative of 1D basis functions @ quadrature points
// returns 2D vector derivativeBasisDunction1D of dimensions nBasisFunction1D,nQuadraturePoints
//arrayDouble QkGL::getDerivativeBasisFunction1D( int order, vectorDouble & quadraturePoints ) const
#ifdef SEM_USE_RAJA
void QkGL::getDerivativeBasisFunction1D( int order, vectorDouble const & quadraturePoints, 
                                          arrayDouble const & derivativeBasisFunction1D ) const
#elif defined SEM_USE_KOKKOS
void QkGL::getDerivativeBasisFunction1D( int order, vectorDouble const & quadraturePoints, 
                                          arrayDouble const & derivativeBasisFunction1D ) const
#else
void QkGL::getDerivativeBasisFunction1D( int order, vectorDouble & quadraturePoints, 
                                          arrayDouble & derivativeBasisFunction1D ) const
#endif
{
  //int nBasisFunction1D=order+1;
  //arrayDouble derivativeBasisFunction1D( nBasisFunction1D, nBasisFunction1D );
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
  //return derivativeBasisFunction1D;
}


// compute 2D gauss-lobatto weights
//vectorDouble QkGL::getGaussLobattoWeights( vectorDouble & quadraturePoints,
//                                           vectorDouble & weights )const
#ifdef SEM_USE_RAJA
void QkGL::getGaussLobattoWeights( vectorDouble const & quadraturePoints,
                                           vectorDouble const & weights,
                                           vectorDouble const & W )const
#elif defined SEM_USE_KOKKOS
void QkGL::getGaussLobattoWeights( vectorDouble const & quadraturePoints,
                                           vectorDouble const & weights,
                                           vectorDouble const & W )const
#else
void QkGL::getGaussLobattoWeights( vectorDouble & quadraturePoints,
                                           vectorDouble & weights,
                                           vectorDouble & W )const
#endif
{
  //vectorDouble W( quadraturePoints.size()*quadraturePoints.size());
  for( int j=0; j<quadraturePoints.size(); j++ )
  {
    for( int i=0; i<quadraturePoints.size(); i++ )
    {
      W[i+j*quadraturePoints.size()]= weights[i]*weights[j];
    }
  }
  //return W;
}

// returns 2D vector basisFunction2D of dimensions nBasisFunctions,nQuadraturePoints
//arrayDouble QkGL::getBasisFunction2D( vectorDouble & quadraturePoints,
//                                      arrayDouble & a,
//                                      arrayDouble & b )const
#ifdef SEM_USE_RAJA
void QkGL::getBasisFunction2D( vectorDouble const & quadraturePoints,
                               arrayDouble const & a,
                               arrayDouble const & b,
                               arrayDouble const & c )const                                      
#elif defined SEM_USE_KOKKOS
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
  int nBasisFunctions=quadraturePoints.size()*quadraturePoints.size();
  //arrayDouble c( nBasisFunctions, nBasisFunctions );
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
  //return c;
}
// compute jacobian matrix for element to ref element coordinates
// Xi[0,..,1][0,..,nPointsPerElement], global coordinate of element
// dxPhi,dyPhi 2D derivative of basis Functions
//arrayDouble QkGL::computeJacobianMatrix( const int & nPointsPerElement,
//                                         arrayDouble & Xi,
//                                         arrayDouble & dxPhi,
//                                         arrayDouble & dyPhi )const

#ifdef SEM_USE_RAJA
int QkGL::computeJacobianMatrix( const int & nPointsPerElement,
                                  arrayDouble const & Xi,
                                  arrayDouble const & dxPhi,
                                  arrayDouble const & dyPhi,
                                  arrayDouble const & jacobianMatrix )const
#elif defined SEM_USE_KOKKOS
KOKKOS_FUNCTION int QkGL::computeJacobianMatrix( const int & nPointsPerElement,
                                  arrayDouble const & Xi,
                                  arrayDouble const & dxPhi,
                                  arrayDouble const & dyPhi,
                                  arrayDouble const & jacobianMatrix )const
#else
int QkGL::computeJacobianMatrix( const int & nPointsPerElement,
                                  arrayDouble & Xi,
                                  arrayDouble & dxPhi,
                                  arrayDouble & dyPhi,
                                  arrayDouble & jacobianMatrix )const
#endif
{
  //arrayDouble jacobianMatrix( 4, nPointsPerElement );

  for( int i=0; i<nPointsPerElement; i++ )
  {
    jacobianMatrix(0,i)=0;
    jacobianMatrix(1,i)=0;
    jacobianMatrix(2,i)=0;
    jacobianMatrix(3,i)=0;
    for( int j=0; j<nPointsPerElement; j++ )
    {
      jacobianMatrix(0,i)+=Xi(j,0)*dxPhi(j,i);
      jacobianMatrix(1,i)+=Xi(j,0)*dyPhi(j,i);
      jacobianMatrix(2,i)+=Xi(j,1)*dxPhi(j,i);
      jacobianMatrix(3,i)+=Xi(j,1)*dyPhi(j,i);
    }
  }
  return 0;
  //return jacobianMatrix;
}

// compute jacobian matrix determinant
//vectorDouble QkGL::computeDeterminantOfJacobianMatrix( const int & nPointsPerElement,
//                                                      arrayDouble & jacobianMatrix ) const
#ifdef SEM_USE_RAJA
int QkGL::computeDeterminantOfJacobianMatrix( const int & nPointsPerElement,
                                               arrayDouble const & jacobianMatrix,
                                               vectorDouble const & detJ ) const
#elif defined SEM_USE_KOKKOS
KOKKOS_FUNCTION int QkGL::computeDeterminantOfJacobianMatrix( const int & nPointsPerElement,
                                               arrayDouble const & jacobianMatrix,
                                               vectorDouble const & detJ ) const
#else
int QkGL::computeDeterminantOfJacobianMatrix( const int & nPointsPerElement,
                                               arrayDouble & jacobianMatrix,
                                               vectorDouble & detJ ) const
#endif
{
  //vectorDouble detJ( nPointsPerElement );
  for( int i=0; i<nPointsPerElement; i++ )
  {
    detJ[i]=(jacobianMatrix(0,i)*jacobianMatrix(3,i)-jacobianMatrix(2,i)*jacobianMatrix(1,i));
  }
  return 0;
  //return detJ;
}

// compute inverse of Jacobian Matrix
//arrayDouble QkGL::computeInvJacobianMatrix( const int & nPointsPerElement,
//                                            arrayDouble & jacobianMatrix,
//                                            vectorDouble & detJ ) const
#ifdef SEM_USE_RAJA
int QkGL::computeInvJacobianMatrix( const int & nPointsPerElement,
                                     arrayDouble const & jacobianMatrix,
                                     vectorDouble const & detJ,
                                     arrayDouble const & invJacobianMatrix ) const
#elif defined SEM_USE_KOKKOS
KOKKOS_FUNCTION int QkGL::computeInvJacobianMatrix( const int & nPointsPerElement,
                                     arrayDouble const & jacobianMatrix,
                                     vectorDouble const & detJ,
                                     arrayDouble const & invJacobianMatrix ) const
#else
int QkGL::computeInvJacobianMatrix( const int & nPointsPerElement,
                                     arrayDouble & jacobianMatrix,
                                     vectorDouble & detJ,
                                     arrayDouble & invJacobianMatrix ) const
#endif
{
  //arrayDouble invJacobianMatrix( 4, nPointsPerElement );
  for( int i=0; i<nPointsPerElement; i++ )
  {
    invJacobianMatrix(0,i)=(jacobianMatrix(3,i)/detJ[i]);
    invJacobianMatrix(1,i)=(-jacobianMatrix(1,i)/detJ[i]);
    invJacobianMatrix(2,i)=(-jacobianMatrix(2,i)/detJ[i]);
    invJacobianMatrix(3,i)=(jacobianMatrix(0,i)/detJ[i]);
  }
  return 0;
  //return invJacobianMatrix;
}

// compute inverse of Jacobian Matrix
//arrayDouble QkGL::computeTranspInvJacobianMatrix( const int & nPointsPerElement,
//                                                  arrayDouble & jacobianMatrix,
//                                                  vectorDouble & detJ ) const
#ifdef SEM_USE_RAJA
int QkGL::computeTranspInvJacobianMatrix( const int & nPointsPerElement,
                                                  arrayDouble const & jacobianMatrix,
                                                  vectorDouble const & detJ,
                                                  arrayDouble  const &transpInvJacobianMatrix ) const
#elif defined SEM_USE_KOKKOS
KOKKOS_FUNCTION int QkGL::computeTranspInvJacobianMatrix( const int & nPointsPerElement,
                                                  arrayDouble const & jacobianMatrix,
                                                  vectorDouble const & detJ,
                                                  arrayDouble  const &transpInvJacobianMatrix ) const
#else
int QkGL::computeTranspInvJacobianMatrix( const int & nPointsPerElement,
                                                  arrayDouble & jacobianMatrix,
                                                  vectorDouble & detJ,
                                                  arrayDouble  &transpInvJacobianMatrix ) const
#endif
{
  //arrayDouble transpInvJacobianMatrix( 4, nPointsPerElement );
  for( int i=0; i<nPointsPerElement; i++ )
  {
    transpInvJacobianMatrix(0,i)=(jacobianMatrix(3,i)/detJ[i]);
    transpInvJacobianMatrix(1,i)=(-jacobianMatrix(2,i)/detJ[i]);
    transpInvJacobianMatrix(2,i)=(-jacobianMatrix(1,i)/detJ[i]);
    transpInvJacobianMatrix(3,i)=(jacobianMatrix(0,i)/detJ[i]);
  }
  return 0;
  //return transpInvJacobianMatrix;
}

// compute B the matrix containing the geometrical informations
//arrayDouble QkGL::computeB( const int & nPointsPerElement,
//                            arrayDouble & invJacobianMatrix,
//                            arrayDouble & transpInvJacobianMatrix,
//                            vectorDouble & detJ ) const
#ifdef SEM_USE_RAJA
int QkGL::computeB( const int & nPointsPerElement,
                 arrayDouble const & invJacobianMatrix,
                 arrayDouble const & transpInvJacobianMatrix,
                 vectorDouble const & detJ,
                 arrayDouble const & B ) const
#elif defined SEM_USE_KOKKOS
KOKKOS_FUNCTION int QkGL::computeB( const int & nPointsPerElement,
                 arrayDouble const & invJacobianMatrix,
                 arrayDouble const & transpInvJacobianMatrix,
                 vectorDouble const & detJ,
                 arrayDouble const & B ) const
#else
int QkGL::computeB( const int & nPointsPerElement,
                 arrayDouble & invJacobianMatrix,
                 arrayDouble & transpInvJacobianMatrix,
                 vectorDouble & detJ,
                 arrayDouble & B ) const
#endif
{
  //arrayDouble B( 4, nPointsPerElement );
  for( int i=0; i<nPointsPerElement; i++ )
  {
    B(0,i)=(abs( detJ[i] )*(invJacobianMatrix(0,i)*transpInvJacobianMatrix(0,i)+
                             invJacobianMatrix(1,i)*transpInvJacobianMatrix(2,i)));
    B(1,i)=(abs( detJ[i] )*(invJacobianMatrix(0,i)*transpInvJacobianMatrix(1,i)+
                             invJacobianMatrix(1,i)*transpInvJacobianMatrix(3,i)));
    B(2,i)=(abs( detJ[i] )*(invJacobianMatrix(2,i)*transpInvJacobianMatrix(0,i)+
                             invJacobianMatrix(3,i)*transpInvJacobianMatrix(2,i)));
    B(3,i)=(abs( detJ[i] )*(invJacobianMatrix(2,i)*transpInvJacobianMatrix(1,i)+
                             invJacobianMatrix(3,i)*transpInvJacobianMatrix(3,i)));
  }
  return 0;
  //return B;

}

// compute the matrix $R_{i,j}=\int_{K}{\nabla{\phi_i}.\nabla{\phi_j}dx}$
// Marc Durufle Formulae
//arrayDouble QkGL::gradPhiGradPhi( const int & nPointsPerElement,
//                                  const int & order,
//                                  vectorDouble & weights2D,
//                                  arrayDouble & B,
//                                  arrayDouble & dPhi ) const
#ifdef SEM_USE_RAJA
int QkGL::gradPhiGradPhi( const int & nPointsPerElement,
                                  const int & order,
                                  vectorDouble const & weights2D,
                                  arrayDouble const & B,
                                  arrayDouble const & dPhi,
                                  arrayDouble const & R ) const
#elif defined SEM_USE_KOKKOS
KOKKOS_FUNCTION int QkGL::gradPhiGradPhi( const int & nPointsPerElement,
                                  const int & order,
                                  vectorDouble const & weights2D,
                                  arrayDouble const & B,
                                  arrayDouble const & dPhi,
                                  arrayDouble const & R ) const
#else
int QkGL::gradPhiGradPhi( const int & nPointsPerElement,
                                  const int & order,
                                  vectorDouble & weights2D,
                                  arrayDouble & B,
                                  arrayDouble & dPhi,
                                  arrayDouble & R ) const
#endif
{
  //arrayDouble R( nPointsPerElement, nPointsPerElement );
  for (int i=0;i<nPointsPerElement;i++)
  {
      for (int j=0; j<nPointsPerElement;j++)
      {
        R(i,j)=0;
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
          R(i,j)+=weights2D[m+i2*(order+1)]*(B(0,m+i2*(order+1))*dPhi(i1,m)*dPhi(j1,m));
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
          R(i,j)+=weights2D[i1+j2*(order+1)]*(B(1,i1+j2*(order+1))*dPhi(i2,j2)*dPhi(j1,i1));
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
          R(i,j)+=weights2D[i2+j1*(order+1)]*(B(2,i2+j1*(order+1))*dPhi(i1,j1)*dPhi(j2,i2));
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
          R(i,j)+=weights2D[i1+n*(order+1)]*(B(3,i1+n*(order+1))*dPhi(i2,n)*dPhi(j2,n));
        }
      }
    }
  }
  return 0;
  //return R;
}
///**
// compute the matrix $R_{i,j}=\int_{K}{\nabla{\phi_i}.\nabla{\phi_j}dx}$
//arrayDouble QkGL::gradPhiGradPhi( const int & nPointsPerElement,
//                                  vectorDouble & weights2D,
//                                  arrayDouble & B,
//                                  arrayDouble & dxPhi,
//                                  arrayDouble & dyPhi ) const
#ifdef SEM_USE_RAJA
int QkGL::gradPhiGradPhi( const int & nPointsPerElement,
                                  vectorDouble const & weights2D,
                                  arrayDouble const & B,
                                  arrayDouble const & dxPhi,
                                  arrayDouble const & dyPhi,
                                  arrayDouble const & R ) const
#elif defined SEM_USE_KOKKOS
KOKKOS_FUNCTION int QkGL::gradPhiGradPhi( const int & nPointsPerElement,
                                  vectorDouble const & weights2D,
                                  arrayDouble const & B,
                                  arrayDouble const & dxPhi,
                                  arrayDouble const & dyPhi,
                                  arrayDouble const & R ) const
#else
int QkGL::gradPhiGradPhi( const int & nPointsPerElement,
                                  vectorDouble & weights2D,
                                  arrayDouble & B,
                                  arrayDouble & dxPhi,
                                  arrayDouble & dyPhi,
                                  arrayDouble & R ) const
#endif
{
  //arrayDouble R( nPointsPerElement, nPointsPerElement );
  for( int i=0; i<nPointsPerElement; i++ )
  {
    for( int j=0; j<nPointsPerElement; j++ )
    {
      double tmp=0;
      for( int r=0; r<nPointsPerElement; r++ )
      {
        tmp+=weights2D[r]*(B(0,r)*dxPhi(i,r)*dxPhi(j,r)+
                           B(1,r)*dxPhi(i,r)*dyPhi(j,r)+
                           B(2,r)*dyPhi(i,r)*dxPhi(j,r)+
                           B(3,r)*dyPhi(i,r)*dyPhi(j,r));
      }
      R(i,j)=tmp;
    }
  }
  return 0;
  //return R;
}

//**/
// compute the matrix $M_{i,j}=\int_{K}{{\phi_i}.{\phi_j}dx}$ (optimized formulation)
//vectorDouble QkGL::phiIphiJ( const int & nPointsPerElement,
//                             vectorDouble & weights2D,
//                             vectorDouble & detJ ) const
#ifdef SEM_USE_RAJA
int QkGL::phiIphiJ( const int & nPointsPerElement,
                             vectorDouble const & weights2D,
                             vectorDouble const & detJ,
                             vectorDouble const & massMatrixLocal ) const
#elif defined SEM_USE_KOKKOS
KOKKOS_FUNCTION int QkGL::phiIphiJ( const int & nPointsPerElement,
                             vectorDouble const & weights2D,
                             vectorDouble const & detJ,
                             vectorDouble const & massMatrixLocal ) const
#else
int QkGL::phiIphiJ( const int & nPointsPerElement,
                             vectorDouble & weights2D,
                             vectorDouble & detJ,
                             vectorDouble & massMatrixLocal ) const
#endif
{
  //vectorDouble massMatrixLocal( nPointsPerElement );
  for( int i=0; i<nPointsPerElement; i++ )
  {
    massMatrixLocal[i]=weights2D[i]*abs( detJ[i] );
  }
  return 0;
  //return MassMatrixLocal;
}

///**
// compute the matrix $M_{i,j}=\int_{K}{{\phi_i}.{\phi_j}dx}$ (non optimized formulation)
//arrayDouble QkGL::phiIphiJ( const int & nPointsPerElement,
//                            vectorDouble & weights2D,
//                            arrayDouble & phi,
//                            vectorDouble & detJ ) const
#ifdef SEM_USE_RAJA
int QkGL::phiIphiJ( const int & nPointsPerElement,
                            vectorDouble const & weights2D,
                            arrayDouble const & phi,
                            vectorDouble const & detJ,
                            arrayDouble  const & massMatrixLocal) const
#elif defined SEM_USE_KOKKOS
KOKKOS_FUNCTION int QkGL::phiIphiJ( const int & nPointsPerElement,
                            vectorDouble const & weights2D,
                            arrayDouble const & phi,
                            vectorDouble const & detJ,
                            arrayDouble  const & massMatrixLocal) const
#else
int QkGL::phiIphiJ( const int & nPointsPerElement,
                            vectorDouble & weights2D,
                            arrayDouble & phi,
                            vectorDouble & detJ,
                            arrayDouble  & massMatrixLocal) const
#endif
{
  //arrayDouble massMatrixLocal( nPointsPerElement, nPointsPerElement );
  for( int i=0; i<nPointsPerElement; i++ )
  {
    for( int j=0; j<nPointsPerElement; j++ )
    {
      double tmp=0;
      for( int r=0; r<nPointsPerElement; r++ )
      {
        tmp+=weights2D[r]*(phi(i,r)*phi(j,r))*abs( detJ[r] );
      }
      massMatrixLocal(i,j)=tmp;
    }
  }
  return 0;
  //return massMatrixLocal;
}
//**/

//vectorReal QkGL::computeDs( const int & iFace,
//                            const int & order,
//                            arrayInt & faceInfos,
//                            arrayReal & globalNodesCoords,
//                            arrayDouble & derivativeBasisFunction2DX,
//                            arrayDouble & derivativeBasisFunction2DY ) const

#ifdef SEM_USE_RAJA
int QkGL::computeDs( const int & iFace,
                      const int & order,
                      arrayInt  const & faceInfos,
                      vectorInt const & numOfBasisFunctionOnFace,
                      arrayReal const & Js,
                      arrayReal const & globalNodesCoords,
                      arrayDouble const & derivativeBasisFunction2DX,
                      arrayDouble const & derivativeBasisFunction2DY,
                      vectorReal  const & ds ) const
#elif defined SEM_USE_KOKKOS
KOKKOS_FUNCTION int QkGL::computeDs( const int & iFace,
                      const int & order,
                      arrayInt  const & faceInfos,
                      vectorInt const & numOfBasisFunctionOnFace,
                      arrayReal const & Js,
                      arrayReal const & globalNodesCoords,
                      arrayDouble const & derivativeBasisFunction2DX,
                      arrayDouble const & derivativeBasisFunction2DY,
                      vectorReal  const & ds ) const
#else
int QkGL::computeDs( const int & iFace,
                      const int & order,
                      arrayInt  & faceInfos,
                      vectorInt & numOfBasisFunctionOnFace,
                      arrayReal & Js,
                      arrayReal & globalNodesCoords,
                      arrayDouble & derivativeBasisFunction2DX,
                      arrayDouble & derivativeBasisFunction2DY,
                      vectorReal  & ds ) const
#endif
{
  //vectorInt numOfBasisFunctionOnFace( order+1 );
  //arrayReal Js( 2, order+1 );
  //vectorReal ds( order+1 );

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
    Js(0,j)=0;    // x
    Js(1,j)=0;    // y
    for( int i=0; i<order+1; i++ )
    {
      float xi=globalNodesCoords(faceInfos(iFace,2+i),0);
      float yi=globalNodesCoords(faceInfos(iFace,2+i),1);
      if( face==0 || face==2 )
      {
        Js(0,j)+=derivativeBasisFunction2DY(numOfBasisFunctionOnFace[i],numOfBasisFunctionOnFace[j])*xi;
        Js(1,j)+=derivativeBasisFunction2DY(numOfBasisFunctionOnFace[i],numOfBasisFunctionOnFace[j])*yi;
      }
      if( face==1 || face==3 )
      {
        Js(0,j)+=derivativeBasisFunction2DX(numOfBasisFunctionOnFace[i],numOfBasisFunctionOnFace[j])*xi;
        Js(1,j)+=derivativeBasisFunction2DX(numOfBasisFunctionOnFace[i],numOfBasisFunctionOnFace[j])*yi;
      }
    }
    ds[j]=sqrt( Js(0,j)*Js(0,j)+Js(1,j)*Js(1,j) );
    //cout<<"j="<<j<<", ds="<<ds[j]<<", "<<Js[0][j]<<", "<<Js[1][j]<<endl;
  }
  return 0;
  //return ds;
}

}
