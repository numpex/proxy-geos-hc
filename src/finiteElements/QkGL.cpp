// C++ Code generated from Python Code: 
#include <iostream>
#include <vector>
#include <cmath>
#include "QkGL.hpp"
#include "commonMacro.hpp"

using namespace std;

QkGL::QkGL(){};
QkGL::~QkGL(){};

vector<double> QkGL::gaussLobattoQuadraturePoints(int order)
{
    vector<double> quadraturePoints;
    quadraturePoints.reserve(order+1);
    if (order == 1)
    {
        quadraturePoints.push_back(-1.0);
        quadraturePoints.push_back(1.0);
    }
    if (order == 2)
    {
        quadraturePoints.push_back(-1.0);
        quadraturePoints.push_back(0.0);
        quadraturePoints.push_back(1.0);
    }
    if (order == 3)
    {
        quadraturePoints.push_back(-1.0);
        quadraturePoints.push_back(-0.4472136);
        quadraturePoints.push_back(0.4472136);
        quadraturePoints.push_back(1.0);
    }
    if (order == 4)
    {
        quadraturePoints.push_back(-1.0);
        quadraturePoints.push_back(-0.65465367);
        quadraturePoints.push_back(0.0);
        quadraturePoints.push_back(0.65465367);
        quadraturePoints.push_back(1.0);
    }
    if (order == 5)
    {
        quadraturePoints.push_back(-1.0);
        quadraturePoints.push_back(-0.76505532);
        quadraturePoints.push_back(-0.28523152);
        quadraturePoints.push_back(0.28523152);
        quadraturePoints.push_back(0.76505532);
        quadraturePoints.push_back(1.0);
    }
    return quadraturePoints;
}
vector<double> QkGL::gaussLobattoQuadratureWeights(int order)
{
    vector<double> weights;
    weights.reserve(order+1);
    if (order == 1)
    {
        weights.push_back(1.0);
        weights.push_back(1.0);
    }
    if (order == 2)
    {
        weights.push_back(0.33333333);
        weights.push_back(1.33333333);
        weights.push_back(0.33333333);
    }
    if (order == 3)
    {
        weights.push_back(0.16666667);
        weights.push_back(0.83333333);
        weights.push_back(0.83333333);
        weights.push_back(0.16666667);
    }
    if (order == 4)
    {
        weights.push_back(0.1);
        weights.push_back(0.54444444);
        weights.push_back(0.71111111);
        weights.push_back(0.54444444);
        weights.push_back(0.1);
    }
    if (order == 5)
    {
        weights.push_back(0.06666667);
        weights.push_back(0.37847496);
        weights.push_back(0.55485838);
        weights.push_back(0.55485838);
        weights.push_back(0.37847496);
        weights.push_back(0.06666667);
    }
    return weights;
}

vector<double> QkGL::shapeFunction1D(int order, double xi)
{
   vector<double> shapeFunction;
   shapeFunction.reserve(order+1);
   if (order==1)
   {
      shapeFunction.push_back(0.5*(1.0-xi));

      shapeFunction.push_back(0.5*(1.0+xi));
   }
   if (order==2)
   {
      shapeFunction.push_back(-1.0*xi*(0.5 - 0.5*xi));

      shapeFunction.push_back((1.0 - 1.0*xi)*(1.0*xi + 1.0));

      shapeFunction.push_back(1.0*xi*(0.5*xi + 0.5));
   }
   if (order==3)
   {
      shapeFunction.push_back((0.309016994374947 - 0.690983005625053*xi)*(0.5 - 0.5*xi)
		             *(-1.80901699437495*xi - 0.809016994374947));

      shapeFunction.push_back((0.5 - 1.11803398874989*xi)*(0.690983005625053 - 0.690983005625053*xi)
		             *(1.80901699437495*xi + 1.80901699437495));

      shapeFunction.push_back((1.80901699437495 - 1.80901699437495*xi)
		             *(0.690983005625053*xi + 0.690983005625053)*(1.11803398874989*xi + 0.5));

      shapeFunction.push_back((0.5*xi + 0.5)*(0.690983005625053*xi + 0.309016994374947)
		             *(1.80901699437495*xi - 0.809016994374947));
   }
   if (order==4)
   {
      shapeFunction.push_back(1.0*xi*(0.39564392373896 - 0.60435607626104*xi)*(0.5 - 0.5*xi)
		            *(-2.89564392373896*xi - 1.89564392373896));

      shapeFunction.push_back(-1.52752523165195*xi*(0.5 - 0.763762615825973*xi)*(0.60435607626104 - 0.60435607626104*xi)
		            *(2.89564392373896*xi + 2.89564392373896));

      shapeFunction.push_back((1.0 - 1.52752523165195*xi)*(1.0 - 1.0*xi)*(1.0*xi + 1.0)*(1.52752523165195*xi + 1.0));

      shapeFunction.push_back(1.52752523165195*xi*(2.89564392373896 - 2.89564392373896*xi) 
		            *(0.60435607626104*xi + 0.60435607626104)*(0.763762615825973*xi + 0.5));

      shapeFunction.push_back(1.0*xi*(0.5*xi + 0.5)*(0.60435607626104*xi + 0.39564392373896) 
		            *(2.89564392373896*xi - 1.89564392373896));
   }
   if (order==5)
   {
      shapeFunction.push_back((0.221930066935875 - 0.778069933064125*xi)*(0.433445520691247 - 0.566554479308753*xi)
                             *(0.5 - 0.5*xi)*(-4.25632117622354*xi - 3.25632117622354)
                             *(-1.39905441140358*xi - 0.399054411403579));

      shapeFunction.push_back((0.271574874126072 - 0.952120850728289*xi)*(0.5 - 0.6535475074298*xi)
                             *(0.566554479308753 - 0.566554479308753*xi)*(-2.0840983387567*xi - 0.594450529658367)
                             *(4.25632117622354*xi + 4.25632117622354));

      shapeFunction.push_back((0.5 - 1.75296196636787*xi)*(0.728425125873928 - 0.952120850728289*xi)
                             *(0.778069933064125 - 0.778069933064125*xi)*(1.39905441140358*xi + 1.39905441140358)
                             *(2.0840983387567*xi + 1.59445052965837));

      shapeFunction.push_back((1.39905441140358 - 1.39905441140358*xi)*(1.59445052965837 - 2.0840983387567*xi)
                             *(0.778069933064125*xi + 0.778069933064125)*(0.952120850728289*xi + 0.728425125873928)
                             *(1.75296196636787*xi + 0.5));

      shapeFunction.push_back((4.25632117622354 - 4.25632117622354*xi)*(0.566554479308753*xi + 0.566554479308753)
                             *(0.6535475074298*xi + 0.5)*(0.952120850728289*xi + 0.271574874126072)
                             *(2.0840983387567*xi - 0.594450529658367));

      shapeFunction.push_back((0.5*xi + 0.5)*(0.566554479308753*xi + 0.433445520691247)
                             *(0.778069933064125*xi + 0.221930066935875)*(1.39905441140358*xi - 0.399054411403579)
                             *(4.25632117622354*xi - 3.25632117622354));
   }
   return shapeFunction;
}

vector<double> QkGL::derivativeShapeFunction1D(int order,double xi)
{
   vector<double> derivativeShapeFunction;
   derivativeShapeFunction.reserve(order+1);
   if (order == 1)
   {
      derivativeShapeFunction.push_back(-0.5);
      derivativeShapeFunction.push_back(0.5);
   }
   if (order == 2)
   {
      derivativeShapeFunction.push_back(1.0*xi - 0.5);
      derivativeShapeFunction.push_back(-2.0*xi);
      derivativeShapeFunction.push_back(1.0*xi + 0.5);
   }
   if (order == 3)
   {
      derivativeShapeFunction.push_back(-1.80901699437495*(0.309016994374947 - 0.690983005625053*xi)*(0.5 - 0.5*xi)
                                 + (-1.80901699437495*xi - 0.809016994374947)*(0.345491502812526*xi - 0.345491502812526)
                                 + (-1.80901699437495*xi - 0.809016994374947)*(0.345491502812526*xi - 0.154508497187474));

      derivativeShapeFunction.push_back(1.80901699437495*(0.5 - 1.11803398874989*xi)*(0.690983005625053 - 0.690983005625053*xi) 
                                 + (0.772542485937369*xi - 0.772542485937369)*(1.80901699437495*xi + 1.80901699437495) 
                                 + (0.772542485937369*xi - 0.345491502812526)*(1.80901699437495*xi + 1.80901699437495));
				 
      derivativeShapeFunction.push_back((1.80901699437495 - 1.80901699437495*xi)*(0.772542485937369*xi + 0.345491502812526) + 
                                 (1.80901699437495 - 1.80901699437495*xi)*(0.772542485937369*xi + 0.772542485937369) -
                                 1.80901699437495*(0.690983005625053*xi + 0.690983005625053)*(1.11803398874989*xi + 0.5));

      derivativeShapeFunction.push_back((0.345491502812526*xi + 0.154508497187474)*(1.80901699437495*xi - 0.809016994374947) + 
                                 (0.345491502812526*xi + 0.345491502812526)*(1.80901699437495*xi - 0.809016994374947) + 
                                 1.80901699437495*(0.5*xi + 0.5)*(0.690983005625053*xi + 0.309016994374947));
                                       
   }
   if (order == 4)
   {
      derivativeShapeFunction.push_back(2.89564392373896*xi*(0.39564392373896 - 0.60435607626104*xi)*(0.5 - 0.5*xi) + 
                                 0.5*xi*(0.39564392373896 - 0.60435607626104*xi)*(-2.89564392373896*xi - 1.89564392373896)
                                 + 0.60435607626104*xi*(0.5 - 0.5*xi)*(-2.89564392373896*xi - 1.89564392373896) + 
                                 (0.39564392373896 - 0.60435607626104*xi)*(-2.89564392373896*xi - 1.89564392373896)*(0.5*xi - 0.5));

      derivativeShapeFunction.push_back(-4.42316915539091*xi*(0.5 - 0.763762615825973*xi)*(0.60435607626104 - 0.60435607626104*xi)
                                 + 0.923169155390906*xi*(0.5 - 0.763762615825973*xi)*(2.89564392373896*xi + 2.89564392373896) 
                                 + 1.16666666666667*xi*(0.60435607626104 - 0.60435607626104*xi)*(2.89564392373896*xi + 2.89564392373896) 
                                 + (0.60435607626104 - 0.60435607626104*xi)*(1.16666666666667*xi - 0.763762615825973)
                                 *(2.89564392373896*xi + 2.89564392373896));

      derivativeShapeFunction.push_back((1.0 - 1.52752523165195*xi)*(1.0 - 1.0*xi)*(1.52752523165195*xi + 1.0) + 
                                 (1.0 - 1.52752523165195*xi)*(1.0 - 1.0*xi)*(1.52752523165195*xi + 1.52752523165195) 
                                 - 1.0*(1.0 - 1.52752523165195*xi)*(1.0*xi + 1.0)*(1.52752523165195*xi + 1.0)
                                 - 1.52752523165195*(1.0 - 1.0*xi)*(1.0*xi + 1.0)*(1.52752523165195*xi + 1.0));

      derivativeShapeFunction.push_back(1.16666666666667*xi*(2.89564392373896 - 2.89564392373896*xi)*(0.60435607626104*xi + 0.60435607626104) 
                                 + 0.923169155390906*xi*(2.89564392373896 - 2.89564392373896*xi)*(0.763762615825973*xi + 0.5) 
                                 - 4.42316915539091*xi*(0.60435607626104*xi + 0.60435607626104)*(0.763762615825973*xi + 0.5) 
                                 + (2.89564392373896 - 2.89564392373896*xi)*(0.60435607626104*xi + 0.60435607626104)
                                 *(1.16666666666667*xi + 0.763762615825973));

      derivativeShapeFunction.push_back(2.89564392373896*xi*(0.5*xi + 0.5)*(0.60435607626104*xi + 0.39564392373896) 
                                 + 0.60435607626104*xi*(0.5*xi + 0.5)*(2.89564392373896*xi - 1.89564392373896) 
                                 + 0.5*xi*(0.60435607626104*xi + 0.39564392373896)*(2.89564392373896*xi - 1.89564392373896)
                                 + (0.5*xi + 0.5)*(0.60435607626104*xi + 0.39564392373896)*(2.89564392373896*xi - 1.89564392373896));
   }
   if (order == 5)
   {
      derivativeShapeFunction.push_back(-1.39905441140358*(0.221930066935875 - 0.778069933064125*xi)*(0.433445520691247 - 0.566554479308753*xi)
                                      *(0.5 - 0.5*xi)*(-4.25632117622354*xi - 3.25632117622354) 
                                      - 4.25632117622354*(0.221930066935875 - 0.778069933064125*xi)*(0.433445520691247 - 0.566554479308753*xi)
                                      *(0.5 - 0.5*xi)*(-1.39905441140358*xi - 0.399054411403579) + (0.221930066935875 - 0.778069933064125*xi)
                                      *(-4.25632117622354*xi - 3.25632117622354)*(-1.39905441140358*xi - 0.399054411403579)
                                      *(0.283277239654376*xi - 0.283277239654376) + (0.221930066935875 - 0.778069933064125*xi)
                                      *(-4.25632117622354*xi - 3.25632117622354)*(-1.39905441140358*xi - 0.399054411403579)
                                      *(0.283277239654376*xi - 0.216722760345624) - 0.778069933064125*(0.433445520691247 - 0.566554479308753*xi)
                                      *(0.5 - 0.5*xi)*(-4.25632117622354*xi - 3.25632117622354)*(-1.39905441140358*xi - 0.399054411403579));

     derivativeShapeFunction.push_back(-2.0840983387567*(0.271574874126072 - 0.952120850728289*xi)*(0.5 - 0.6535475074298*xi)
                                      *(0.566554479308753 - 0.566554479308753*xi)*(4.25632117622354*xi + 4.25632117622354) 
                                      - 0.566554479308753*(0.271574874126072 - 0.952120850728289*xi)*(0.5 - 0.6535475074298*xi)
                                      *(-2.0840983387567*xi - 0.594450529658367)*(4.25632117622354*xi + 4.25632117622354) 
                                      +(0.271574874126072 - 0.952120850728289*xi)*(0.566554479308753 - 0.566554479308753*xi)
                                      *(2.12816058811177 - 2.78170809554157*xi)*(-2.0840983387567*xi - 0.594450529658367) 
                                      +(0.271574874126072 - 0.952120850728289*xi)*(0.566554479308753 - 0.566554479308753*xi)
                                      *(-2.78170809554157*xi - 2.78170809554157)*(-2.0840983387567*xi - 0.594450529658367) 
                                      - 0.952120850728289*(0.5 - 0.6535475074298*xi)*(0.566554479308753 - 0.566554479308753*xi)
                                      *(-2.0840983387567*xi - 0.594450529658367)*(4.25632117622354*xi + 4.25632117622354));

     derivativeShapeFunction.push_back(2.0840983387567*(0.5 - 1.75296196636787*xi)*(0.728425125873928 - 0.952120850728289*xi)
                                      *(0.778069933064125 - 0.778069933064125*xi)*(1.39905441140358*xi + 1.39905441140358)
                                      + 1.39905441140358*(0.5 - 1.75296196636787*xi)*(0.728425125873928 - 0.952120850728289*xi)
                                      *(0.778069933064125 - 0.778069933064125*xi)*(2.0840983387567*xi + 1.59445052965837)
                                      - 0.952120850728289*(0.5 - 1.75296196636787*xi)*(0.778069933064125 - 0.778069933064125*xi)
                                      *(1.39905441140358*xi + 1.39905441140358)*(2.0840983387567*xi + 1.59445052965837) 
                                      + (0.728425125873928 - 0.952120850728289*xi)*(1.3639269998358*xi - 1.3639269998358)
                                      *(1.39905441140358*xi + 1.39905441140358)*(2.0840983387567*xi + 1.59445052965837)
                                      + (0.728425125873928 - 0.952120850728289*xi)*(1.3639269998358*xi - 0.389034966532063)
                                      *(1.39905441140358*xi + 1.39905441140358)*(2.0840983387567*xi + 1.59445052965837));

     derivativeShapeFunction.push_back(0.952120850728289*(1.39905441140358 - 1.39905441140358*xi)*(1.59445052965837 - 2.0840983387567*xi)
                                      *(0.778069933064125*xi + 0.778069933064125)*(1.75296196636787*xi + 0.5) 
                                      + (1.39905441140358 - 1.39905441140358*xi)*(1.59445052965837 - 2.0840983387567*xi)
                                      *(0.952120850728289*xi + 0.728425125873928)*(1.3639269998358*xi + 0.389034966532063) 
                                      + (1.39905441140358 - 1.39905441140358*xi)*(1.59445052965837 - 2.0840983387567*xi)
                                      *(0.952120850728289*xi + 0.728425125873928)*(1.3639269998358*xi + 1.3639269998358) 
                                      - 2.0840983387567*(1.39905441140358 - 1.39905441140358*xi)*(0.778069933064125*xi + 0.778069933064125)
                                      *(0.952120850728289*xi + 0.728425125873928)*(1.75296196636787*xi + 0.5) 
                                      - 1.39905441140358*(1.59445052965837 - 2.0840983387567*xi)*(0.778069933064125*xi + 0.778069933064125)
                                      *(0.952120850728289*xi + 0.728425125873928)*(1.75296196636787*xi + 0.5));

     derivativeShapeFunction.push_back((2.78170809554157 - 2.78170809554157*xi)*(0.566554479308753*xi + 0.566554479308753)
                                      *(0.952120850728289*xi + 0.271574874126072)*(2.0840983387567*xi - 0.594450529658367) 
                                      + 2.0840983387567*(4.25632117622354 - 4.25632117622354*xi)*(0.566554479308753*xi + 0.566554479308753)
                                      *(0.6535475074298*xi + 0.5)*(0.952120850728289*xi + 0.271574874126072)
                                      + 0.952120850728289*(4.25632117622354 - 4.25632117622354*xi)*(0.566554479308753*xi 
                                      + 0.566554479308753)*(0.6535475074298*xi + 0.5)*(2.0840983387567*xi - 0.594450529658367)
                                      + 0.566554479308753*(4.25632117622354 - 4.25632117622354*xi)*(0.6535475074298*xi + 0.5)
                                      *(0.952120850728289*xi + 0.271574874126072)*(2.0840983387567*xi - 0.594450529658367) 
                                      + (-2.78170809554157*xi - 2.12816058811177)*(0.566554479308753*xi + 0.566554479308753)
                                      *(0.952120850728289*xi + 0.271574874126072)*(2.0840983387567*xi - 0.594450529658367));

     derivativeShapeFunction.push_back((0.283277239654376*xi + 0.216722760345624)*(0.778069933064125*xi + 0.221930066935875)
                                      *(1.39905441140358*xi - 0.399054411403579)*(4.25632117622354*xi - 3.25632117622354)
                                      + (0.283277239654376*xi + 0.283277239654376)*(0.778069933064125*xi + 0.221930066935875)
                                      *(1.39905441140358*xi - 0.399054411403579)*(4.25632117622354*xi - 3.25632117622354)
                                      + 4.25632117622354*(0.5*xi + 0.5)*(0.566554479308753*xi + 0.433445520691247)
                                      *(0.778069933064125*xi + 0.221930066935875)*(1.39905441140358*xi - 0.399054411403579)
                                      + 1.39905441140358*(0.5*xi + 0.5)*(0.566554479308753*xi + 0.433445520691247)
                                      *(0.778069933064125*xi + 0.221930066935875)*(4.25632117622354*xi - 3.25632117622354)
                                      + 0.778069933064125*(0.5*xi + 0.5)*(0.566554479308753*xi + 0.433445520691247)
                                      *(1.39905441140358*xi - 0.399054411403579)*(4.25632117622354*xi - 3.25632117622354));
   }
   return derivativeShapeFunction;
}

// get 1D basis functions @ quadrature points
// returns 2D vector basisDunction1D of dimensions nBasisFunction1D,nQuadraturePoints
vector<vector<double>> QkGL::getBasisFunction1D(int order, const vector<double> & quadraturePoints)
{
    int nBasisFunction1D=order+1;
    vector<vector<double>> basisFunction1D(nBasisFunction1D);
    vector<vector<double>> tmp(order+1);
    // loop over quadrature points
    for (int i = 0; i < order+1; i++)
    {
       //extract all basis functions  for current quadrature point
       tmp[i]=shapeFunction1D(order,quadraturePoints[i]);
     }
    // transpose to get  basis Function values at quadrature nodes
    for ( int j=0; j < order+1; j++)
    {
	for ( int i=0; i < nBasisFunction1D; i++)
        {
                basisFunction1D[i].push_back(tmp[j][i]);
	}
    }
    return basisFunction1D;
}

// get derivative of 1D basis functions @ quadrature points
// returns 2D vector derivativeBasisDunction1D of dimensions nBasisFunction1D,nQuadraturePoints
vector<vector<double>> QkGL::getDerivativeBasisFunction1D(int order,const vector<double> &quadraturePoints)
{
    int nBasisFunction1D=order+1;
    vector<vector<double>> derivativeBasisFunction1D(order+1);
    vector<vector<double>> tmp(order+1);
    // loop over quadrature points
    for (int i = 0; i < order+1; i++)
    {
       //extract all basis functions  for current quadrature point
       tmp[i]=derivativeShapeFunction1D(order,quadraturePoints[i]);
     }
    // transpose to get  basis Function values at quadrature nodes
    for ( int j=0; j < order+1; j++)
    {
	for ( int i=0; i < nBasisFunction1D; i++)
        {
                derivativeBasisFunction1D[i].push_back(tmp[j][i]);
	}
    }
    return derivativeBasisFunction1D;
}


// compute 2D gauss-lobatto weights
vector<double> QkGL::getGaussLobattoWeights(const vector<double> &quadraturePoints,
                                            const vector<double> &weights)
{
    vector<double>W;
    for (int j=0; j<quadraturePoints.size(); j++)
    {
        for (int i=0; i<quadraturePoints.size(); i++)
        {
            W.push_back(weights[i]*weights[j]);
        }
    }
    return W;
}
   
// returns 2D vector basisFunction2D of dimensions nBasisFunctions,nQuadraturePoints
vector<vector<double>> QkGL::getBasisFunction2D(const vector<double> quadraturePoints,
		                                        const vector<vector<double>> &a,
					                            const vector<vector<double>> &b)
{
    int nBasisFunctions=quadraturePoints.size()*quadraturePoints.size();
    vector<vector<double>> c(nBasisFunctions);

    for (int j = 0; j < quadraturePoints.size(); j++)
    {	    
        for (int k = 0; k<quadraturePoints.size();k++)
	{	
            for ( int i = 0; i< quadraturePoints.size(); i++)
            {
	        for ( int l=0; l<quadraturePoints.size(); l++)
                { 
	            c[i+quadraturePoints.size()*j].push_back(a[i][l]*b[j][k]);
            	}
            }
        }
    }
    return c;
}

// compute jacobian matrix for element to ref element coordinates
// Xi[0,..,1][0,..,nPointsPerElement], global coordinate of element
// dxPhi,dyPhi 2D derivative of basis Functions
vector<vector<double>> QkGL::computeJacobianMatrix(const int &nPointsPerElement,
	                                               const vector<vector<double>> &Xi ,
	                                               const vector<vector<double>> &dxPhi ,
	                                               const vector<vector<double>> &dyPhi )
{
    vector<vector<double>> jacobianMatrix(4,vector<double>(nPointsPerElement,0));

    for (int i=0; i<nPointsPerElement; i++)
    {
        for (int j=0; j<nPointsPerElement; j++) 
        {
            jacobianMatrix[0][i]+=Xi[j][0]*dxPhi[j][i];
            jacobianMatrix[1][i]+=Xi[j][0]*dyPhi[j][i];
            jacobianMatrix[2][i]+=Xi[j][1]*dxPhi[j][i];
            jacobianMatrix[3][i]+=Xi[j][1]*dyPhi[j][i];
        }
    }
    return jacobianMatrix;
}

// compute jacobian matrix determinant
vector<double>  QkGL::computeDeterminantOfJacobianMatrix(const int &nPointsPerElement,
	                                                     const vector<vector<double>> &jacobianMatrix)
{
    vector<double>  detJ(nPointsPerElement,0);
    for (int i=0; i<nPointsPerElement; i++) 
    {
        detJ[i]=(jacobianMatrix[0][i]*jacobianMatrix[3][i]-jacobianMatrix[2][i]*jacobianMatrix[1][i]);
    }
    return detJ;
}

// compute inverse of Jacobian Matrix
vector<vector<double>> QkGL::computeInvJacobianMatrix(const int &nPointsPerElement,
	                                                  const vector<vector<double>> &jacobianMatrix,
				                                      const vector<double> &detJ)
{
    vector<vector<double>> invJacobianMatrix(4,vector<double>(nPointsPerElement,0));
    for (int i=0; i<nPointsPerElement; i++)
    {
	invJacobianMatrix[0][i]=(jacobianMatrix[3][i]/detJ[i]);
	invJacobianMatrix[1][i]=(-jacobianMatrix[1][i]/detJ[i]);
	invJacobianMatrix[2][i]=(-jacobianMatrix[2][i]/detJ[i]);
	invJacobianMatrix[3][i]=(jacobianMatrix[0][i]/detJ[i]);
    }
    return invJacobianMatrix;
}

// compute inverse of Jacobian Matrix
vector<vector<double>> QkGL::computeTranspInvJacobianMatrix(const int &nPointsPerElement,
		                                                    const vector<vector<double>> &jacobianMatrix,
					                                        const vector<double> &detJ)
{
    vector<vector<double>> transpInvJacobianMatrix(4,vector<double>(nPointsPerElement,0));
    for (int i=0; i<nPointsPerElement; i++)
    {
	transpInvJacobianMatrix[0][i]=(jacobianMatrix[3][i]/detJ[i]);
	transpInvJacobianMatrix[1][i]=(-jacobianMatrix[2][i]/detJ[i]);
	transpInvJacobianMatrix[2][i]=(-jacobianMatrix[1][i]/detJ[i]);
	transpInvJacobianMatrix[3][i]=(jacobianMatrix[0][i]/detJ[i]);
    }
    return transpInvJacobianMatrix;
}

// compute B the matrix containing the geometrical informations
vector<vector<double>> QkGL::computeB(const int & nPointsPerElement,
                                      const vector<vector<double>> &invJacobianMatrix,
                                      const vector<vector<double>> &transpInvJacobianMatrix,
                                      const vector<double> &detJ)
{
    vector<vector<double>> B(4,vector<double>(nPointsPerElement,0));
    for (int i=0; i<nPointsPerElement; i++)
    {
        B[0][i]=(abs(detJ[i])*(invJacobianMatrix[0][i]*transpInvJacobianMatrix[0][i]+
				     invJacobianMatrix[1][i]*transpInvJacobianMatrix[2][i]));
        B[1][i]=(abs(detJ[i])*(invJacobianMatrix[0][i]*transpInvJacobianMatrix[1][i]+
			             invJacobianMatrix[1][i]*transpInvJacobianMatrix[3][i]));
        B[2][i]=(abs(detJ[i])*(invJacobianMatrix[2][i]*transpInvJacobianMatrix[0][i]+
			             invJacobianMatrix[3][i]*transpInvJacobianMatrix[2][i]));
        B[3][i]=(abs(detJ[i])*(invJacobianMatrix[2][i]*transpInvJacobianMatrix[1][i]+
			             invJacobianMatrix[3][i]*transpInvJacobianMatrix[3][i]));
    }
    return B;

}

// compute the matrix $R_{i,j}=\int_{K}{\nabla{\phi_i}.\nabla{\phi_j}dx}$ 
vector<vector<double>> QkGL::gradPhiGradPhi(const int & nPointsPerElement,
                                            const vector<double> &weights2D,
                                            const vector<vector<double>> &B,
                                            const vector<vector<double>> &dxPhi,
					                        const vector<vector<double>> &dyPhi)
{
    vector<vector<double>> R(nPointsPerElement,vector<double>(nPointsPerElement,0));
    for (int i=0; i<nPointsPerElement; i++)
    {
        for (int j=0; j<nPointsPerElement; j++)
	    {
	        double tmp=0;
            for (int r=0; r<nPointsPerElement; r++)
	        {
		        R[i][j]+=weights2D[r]*(B[0][r]*dxPhi[i][r]*dxPhi[j][r]+
                                       B[1][r]*dxPhi[i][r]*dyPhi[j][r]+
                                       B[2][r]*dyPhi[i][r]*dxPhi[j][r]+
                                       B[3][r]*dyPhi[i][r]*dyPhi[j][r]);
	        }
	    }
    }
    return R;
}

// compute the matrix $M_{i,j}=\int_{K}{{\phi_i}.{\phi_j}dx}$ 
vector<vector<double>> QkGL::phiIphiJ(const int & nPointsPerElement,
                                      const vector<double> &weights2D,
					                  const vector<vector<double>> &phi,
					                  const vector<double> &detJ)
{
    vector<vector<double>> M(nPointsPerElement,vector<double>(nPointsPerElement,0));
    for (int i=0; i<nPointsPerElement; i++)
    {
        for (int j=0; j<nPointsPerElement; j++)
	    {
            for (int r=0; r<nPointsPerElement; r++)
	        {
		        M[i][j]+=weights2D[r]*(phi[i][r]*phi[j][r])*abs(detJ[r]);
	        }
	    }
    }
    return M;
}


 
vector<float>QkGL::computeDs(const int & iFace,
                            const int & order,
                            const vector<vector<int>> & faceInfos,
                            const vector<vector<float>> & globalNodesCoords,
                            const vector<vector<double>>& derivativeBasisFunction2DX,
                            const vector<vector<double>>& derivativeBasisFunction2DY)
{
    vector<int>numOfBasisFunctionOnFace(order+1,0);
    vector<vector<float>>Js(2,vector<float>(order+1,0));
    vector<float>ds(order+1,0);

    int face=faceInfos[iFace][1];
    // get basis functions on Boundary faces
    switch (face)
    {
        case 0: // left
            for (int i=0;i<order+1;i++)
            {
                numOfBasisFunctionOnFace[i]=i*(order+1);
            }
            break;
        case 1: // bottom
            for (int i=0;i<order+1;i++)
            {
                numOfBasisFunctionOnFace[i]=i;
            }
            break;
            case 2: //right
                for (int i=0;i<order+1;i++)
                {
                    numOfBasisFunctionOnFace[i]=order+i*(order+1);
                }
                break;
            case 3: //top
                for (int i=0;i<order+1;i++)
                {
                    numOfBasisFunctionOnFace[i]=i+order*(order+1);
                }   
                break;               
            default :
                cout<<"error in element flag, should be set to: 0, 1, 2, 3"<<endl;
                break;
    }
    // compute ds
    for (int j=0; j<order+1; j++)
    {
        Js[0][j]=0;// x 
        Js[1][j]=0;// y
        for (int i=0; i<order+1; i++)
        {
            float xi=globalNodesCoords[faceInfos[iFace][2+i]][0];
            float yi=globalNodesCoords[faceInfos[iFace][2+i]][1];
            if ( face==0 || face==2)
            {
                Js[0][j]+=derivativeBasisFunction2DY[numOfBasisFunctionOnFace[i]][numOfBasisFunctionOnFace[j]]*xi;
                Js[1][j]+=derivativeBasisFunction2DY[numOfBasisFunctionOnFace[i]][numOfBasisFunctionOnFace[j]]*yi;
            }  
            if ( face==1 || face==3)
            {
                Js[0][j]+=derivativeBasisFunction2DX[numOfBasisFunctionOnFace[i]][numOfBasisFunctionOnFace[j]]*xi; 
                Js[1][j]+=derivativeBasisFunction2DX[numOfBasisFunctionOnFace[i]][numOfBasisFunctionOnFace[j]]*yi;
            }               
        }
        ds[j]=sqrt(Js[0][j]*Js[0][j]+Js[1][j]*Js[1][j]);
        //cout<<"j="<<j<<", ds="<<ds[j]<<endl;
    }
    return ds;
}
            

