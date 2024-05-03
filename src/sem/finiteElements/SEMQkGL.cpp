// C++ Code generated from Python Code:
#include "SEMQkGL.hpp"

void SEMQkGL::gaussLobattoQuadraturePoints( int order, vectorDouble const & quadraturePoints ) const
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
void SEMQkGL::gaussLobattoQuadratureWeights( int order, vectorDouble const & weights ) const
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
vector<double> SEMQkGL::shapeFunction1D( int order, double xi ) const
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
vector<double> SEMQkGL::derivativeShapeFunction1D( int order, double xi ) const
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

void SEMQkGL::getBasisFunction1D( int order, vectorDouble const & quadraturePoints, arrayDouble const & basisFunction1D ) const
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

void SEMQkGL::getDerivativeBasisFunction1D( int order, vectorDouble const & quadraturePoints, 
                                          arrayDouble const & derivativeBasisFunction1D ) const
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

// compute B and M  
PROXY_HOST_DEVICE void SEMQkGL::computeB(const int & elementNumber,
		                       const int & order,
		                       const int & dimension,
                                       VECTOR_DOUBLE_VIEW const & weights,
			               ARRAY_INT_VIEW const & nodesList,
			               ARRAY_REAL_VIEW const & nodesCoords,
                                       ARRAY_DOUBLE_VIEW const & dPhi,
				       float massMatrixLocal[],
                                       float B[][COL] ) const
{
   if ( dimension==2 ) {
   for (int i2=0;i2<order+1;i2++)
   {
       for (int i1=0;i1<order+1;i1++)
       {
          // compute jacobian matrix
          double jac0=0;
          double jac1=0;
          double jac2=0;
          double jac3=0;
          int i=i1+i2*(order+1);
          for (int j1=0;j1<order+1;j1++)
          {
              int j=j1+i2*(order+1);
              int localToGlobal=nodesList(elementNumber,j);
              double X=nodesCoords(localToGlobal,0);
              double Y=nodesCoords(localToGlobal,1);
              jac0+=X*dPhi(j1,i1);
              jac2+=Y*dPhi(j1,i1);
          }
          for (int j2=0; j2<order+1;j2++)
          {
              int j=i1+j2*(order+1);
              int localToGlobal=nodesList(elementNumber,j);
              double X=nodesCoords(localToGlobal,0);
              double Y=nodesCoords(localToGlobal,1);
              jac1+=X*dPhi(j2,i2);
              jac3+=Y*dPhi(j2,i2);
          }
          // detJ
          double detJ=abs(jac0*jac3-jac2*jac1);
          double invJac0=jac3;
          double invJac1=-jac1;
          double invJac2=-jac2;
          double invJac3=jac0;
          double transpInvJac0=jac3;
          double transpInvJac1=-jac2;
          double transpInvJac2=-jac1;
          double transpInvJac3=jac0;
          double detJM1=1./detJ;
          // B
          B[i][0]=(invJac0*transpInvJac0+invJac1*transpInvJac2)*detJM1;
          B[i][1]=(invJac0*transpInvJac1+invJac1*transpInvJac3)*detJM1;
          B[i][2]=(invJac2*transpInvJac0+invJac3*transpInvJac2)*detJM1;
          B[i][3]=(invJac2*transpInvJac1+invJac3*transpInvJac3)*detJM1;
          //M
          massMatrixLocal[i]=weights[i1]*weights[i2]*detJ;
       }
  }}
  else { //3D case
  for (int i3=0;i3<order+1;i3++)
  {
      for (int i2=0;i2<order+1;i2++)
      {
          for (int i1=0;i1<order+1;i1++)
          {
              int i=i1+i2*(order+1)+i3*(order+1)*(order+1);
              // compute jacobian matrix
              double jac00=0;
              double jac01=0;
              double jac02=0;
              double jac10=0;
              double jac11=0;
              double jac12=0;
              double jac20=0;
              double jac21=0;
              double jac22=0;

              for (int j1=0;j1<order+1;j1++)
              {
                  int j=j1+i2*(order+1)+i3*(order+1)*(order+1);
                  int localToGlobal=nodesList(elementNumber,j);
                  double X=nodesCoords(localToGlobal,0);
                  double Y=nodesCoords(localToGlobal,2);
                  double Z=nodesCoords(localToGlobal,1);
                  jac00+=X*dPhi(j1,i1);
                  jac10+=Z*dPhi(j1,i1);
                  jac20+=Y*dPhi(j1,i1);
              }
              for (int j2=0;j2<order+1;j2++)
              {
                  int j=i1+j2*(order+1)+i3*(order+1)*(order+1);
                  int localToGlobal=nodesList(elementNumber,j);
                  double X=nodesCoords(localToGlobal,0);
                  double Y=nodesCoords(localToGlobal,2);
                  double Z=nodesCoords(localToGlobal,1);
                  jac01+=X*dPhi(j2,i2);
                  jac11+=Z*dPhi(j2,i2);
                  jac21+=Y*dPhi(j2,i2);
              }
              for (int j3=0;j3<order+1;j3++)
              {
                  int j=i1+i2*(order+1)+j3*(order+1)*(order+1);
                  int localToGlobal=nodesList(elementNumber,j);
                  double X=nodesCoords(localToGlobal,0);
                  double Y=nodesCoords(localToGlobal,2);
                  double Z=nodesCoords(localToGlobal,1);
                  jac02+=X*dPhi(j3,i3);
                  jac12+=Z*dPhi(j3,i3);
                  jac22+=Y*dPhi(j3,i3);
              }
              // detJ
              double detJ=abs(jac00*(jac11*jac22-jac21*jac12)
                             -jac01*(jac10*jac22-jac20*jac12)
                             +jac02*(jac10*jac21-jac20*jac11));

              // inv of jac is equal of the minors of the transposed of jac
              double invJac00=jac11*jac22-jac12*jac21;
              double invJac01=jac02*jac21-jac01*jac22;
              double invJac02=jac01*jac12-jac02*jac11;
              double invJac10=jac12*jac20-jac10*jac22;
              double invJac11=jac00*jac22-jac02*jac20;
              double invJac12=jac02*jac10-jac00*jac12;
              double invJac20=jac10*jac21-jac11*jac20;
              double invJac21=jac01*jac20-jac00*jac21;
              double invJac22=jac00*jac11-jac01*jac10;

              double transpInvJac00=invJac00;
              double transpInvJac01=invJac10;
              double transpInvJac02=invJac20;
              double transpInvJac10=invJac01;
              double transpInvJac11=invJac11;
              double transpInvJac12=invJac21;
              double transpInvJac20=invJac02;
              double transpInvJac21=invJac12;
              double transpInvJac22=invJac22;

              double detJM1=1./detJ;

              // B
              B[i][0]=(invJac00*transpInvJac00+invJac01*transpInvJac10+invJac02*transpInvJac20)*detJM1;//B11
              B[i][1]=(invJac10*transpInvJac01+invJac11*transpInvJac11+invJac12*transpInvJac21)*detJM1;//B22
              B[i][2]=(invJac20*transpInvJac02+invJac21*transpInvJac12+invJac22*transpInvJac22)*detJM1;//B33
              B[i][3]=(invJac00*transpInvJac01+invJac01*transpInvJac11+invJac02*transpInvJac21)*detJM1;//B12,B21
              B[i][4]=(invJac00*transpInvJac02+invJac01*transpInvJac12+invJac02*transpInvJac22)*detJM1;//B13,B31
              B[i][5]=(invJac10*transpInvJac02+invJac11*transpInvJac12+invJac12*transpInvJac22)*detJM1;//B23,B32

              //M
              massMatrixLocal[i]=weights[i1]*weights[i2]*weights[i3]*detJ;
          }
      }
  }}
}

// compute the matrix $R_{i,j}=\int_{K}{\nabla{\phi_i}.\nabla{\phi_j}dx}$
// Marc Durufle Formulae
PROXY_HOST_DEVICE void SEMQkGL::gradPhiGradPhi( const int & nPointsPerElement,
                                              const int & order,
		                              const int & dimension,
                                              VECTOR_DOUBLE_VIEW const & weights,
                                              ARRAY_DOUBLE_VIEW const & dPhi,
                                              float const  B[][COL],
			                      float const pnLocal[],
                                              float R[],
	                                      float Y[]) const
{
  if ( dimension==2 ) {
  // B11
  for( int i2=0; i2<order+1; i2++ )
  {
    for( int i1=0; i1<order+1; i1++ )
    {
      for (int j=0; j<nPointsPerElement;j++)
      {
        R[j]=0;
      }
      for( int j1=0; j1<order+1; j1++ )
      {
        int j=j1+i2*(order+1);
        for( int m=0; m<order+1; m++ )
        {
          R[j]+=weights[m]*weights[i2]*(B[m+i2*(order+1)][0]*dPhi(i1,m)*dPhi(j1,m));
        }
      }
      // B21
      for( int j1=0; j1<order+1; j1++ )
      {
        for( int j2=0; j2<order+1; j2++ )
        {
          int j=j1+j2*(order+1);
          R[j]+=weights[i1]*weights[j2]*(B[i1+j2*(order+1)][1]*dPhi(i2,j2)*dPhi(j1,i1));
        }
      }
      // B12
      for( int j1=0; j1<order+1; j1++ )
      {
        for( int j2=0; j2<order+1; j2++ )
        {
          int j=j1+j2*(order+1);
          R[j]+=weights[i2]*weights[j1]*(B[i2+j1*(order+1)][2]*dPhi(i1,j1)*dPhi(j2,i2));
        }
      }
      // B22
      for( int j2=0; j2<order+1; j2++ )
      {
        int j=i1+j2*(order+1);
        for( int n=0; n<order+1; n++ )
        {
          R[j]+=weights[i1]*weights[n]*(B[i1+n*(order+1)][3]*dPhi(i2,n)*dPhi(j2,n));
        }
      }
      int i=i1+i2*(order+1);
      Y[i]=0;
      for( int j=0; j<nPointsPerElement; j++ )
      {
         Y[i]+=R[j]*pnLocal[j];
      }
    }
  }}
  else{
  int orderPow2=(order+1)*(order+1);
  for (int i3=0;i3<order+1;i3++)
  {
      for (int i2=0;i2<order+1;i2++)
      {
          for (int i1=0;i1<order+1;i1++)
          {
              for( int j=0; j<nPointsPerElement; j++ )
              {
                 R[j]=0;
              }

	      //B11
              for( int j1=0; j1<order+1; j1++ )
              {
                  int j=j1+i2*(order+1)+i3*orderPow2;
                  for( int l=0; l<order+1; l++ )
                  {
                      int ll=l+i2*(order+1)+i3*orderPow2;
                      R[j]+=weights[l]*weights[i2]*weights[i3]*(B[ll][0]*dPhi(i1,l)*dPhi(j1,l));
                  }
              }
              //B22
              for( int j2=0; j2<order+1; j2++ )
              {
                  int j=i1+j2*(order+1)+i3*orderPow2;
                  for( int m=0; m<order+1; m++ )
                  {
                      int mm=i1+m*(order+1)+i3*orderPow2;
                      R[j]+=weights[i1]*weights[m]*weights[i3]*(B[mm][1]*dPhi(i2,m)*dPhi(j2,m));
                  }
              }
              //B33
              for( int j3=0; j3<order+1; j3++ )
              {
                  int j=i1+i2*(order+1)+j3*orderPow2;
                  for( int n=0; n<order+1; n++ )
                  {
                      int nn=i1+i2*(order+1)+n*orderPow2;
                      R[j]+=weights[i1]*weights[i2]*weights[n]*(B[nn][2]*dPhi(i3,n)*dPhi(j3,n));
                  }
              }
              // B12,B21 (B[][3])
              for( int j2=0; j2<order+1; j2++ )
              {
                for( int j1=0; j1<order+1; j1++ )
                {
                  int j=j1+j2*(order+1)+i3*orderPow2;
                  int k=j1+i2*(order+1)+i3*orderPow2;
                  int l=i1+j2*(order+1)+i3*orderPow2;
                  R[j]+=weights[j1]*weights[i2]*weights[i3]*(B[k][3]*dPhi(i1,j1)*dPhi(j2,i2))+
                        weights[i1]*weights[j2]*weights[i3]*(B[l][3]*dPhi(j1,i1)*dPhi(i2,j2));
                }
              }
              // B13,B31 (B[][4])
              for( int j3=0; j3<order+1; j3++ )
              {
                for( int j1=0; j1<order+1; j1++ )
                {
                  int j=j1+i2*(order+1)+i3*orderPow2;
                  int k=j1+i2*(order+1)+i3*orderPow2;
                  int l=j1+i2*(order+1)+j3*orderPow2;
                  R[j]+=weights[j1]*weights[i2]*weights[i3]*(B[k][4]*dPhi(j1,i1)*dPhi(j3,i3))+
                        weights[j1]*weights[i2]*weights[j3]*(B[l][4]*dPhi(j1,i1)*dPhi(i3,j3));
	        }
              }
              // B23,B32 (B[][5])
              for( int j3=0; j3<order+1; j3++ )
              {
                for( int j2=0; j2<order+1; j2++ )
                {
                  int j=i1+j2*(order+1)+j3*orderPow2;
                  int k=i1+j2*(order+1)+i3*orderPow2;
                  int l=i1+i2*(order+1)+j3*orderPow2;
                  R[j]+=weights[i1]*weights[j2]*weights[i3]*(B[k][5]*dPhi(i2,i2)*dPhi(j3,i3))+
                        weights[i1]*weights[i2]*weights[j3]*(B[l][5]*dPhi(j2,i2)*dPhi(i3,j3));
                }
              }

              int i=i1+i2*(order+1)+i3*orderPow2;
              Y[i]=0;
              for( int j=0; j<nPointsPerElement; j++ )
              {
                Y[i]+=R[j]*pnLocal[j];
              }

          }
      }
  }}
}

//computeDs
PROXY_HOST_DEVICE int SEMQkGL::computeDs( const int & iFace,
                                          const int & order,
                                          ARRAY_INT_VIEW const & faceInfos,
                                          int  numOfBasisFunctionOnFace[],
                                          float  Js[][6],
                                          ARRAY_REAL_VIEW const & globalNodesCoords,
                                          ARRAY_DOUBLE_VIEW const & dPhi,
                                          float  ds[]  ) const
{
  // compute ds
  int face=faceInfos(iFace,1);
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
        Js[0][j]+=dPhi(i,j)*xi;
        Js[1][j]+=dPhi(i,j)*yi;
      }
      if( face==1 || face==3 )
      {
        Js[0][j]+=dPhi(i,j)*xi;
        Js[1][j]+=dPhi(i,j)*yi;
      }
    }
    ds[j]=sqrt( Js[0][j]*Js[0][j]+Js[1][j]*Js[1][j] );
    //cout<<"j="<<j<<", ds="<<ds[j]<<", "<<Js[0][j]<<", "<<Js[1][j]<<endl;
  }
  return 0;
}

