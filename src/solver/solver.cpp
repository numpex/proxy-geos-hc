#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>

#include "solver.hpp"

using namespace std;

solver::solver() {}
solver::~solver() {}

// compute one step of the time dynamic wave equation solver
//vector<vector<float>> 
void solver::computeOneStep(const float & timeSample,
		                    const int &order,
					        int &i1,
					        int &i2,
                            vector<vector<float>> & pnGlobal,
					        simpleMesh mesh,
					        QkGL Qk)
{
    static int numberOfNodes=mesh.getNumberOfNodes();
    static int numberOfElements=mesh.getNumberOfElements();
    static vector<vector<int>> globalNodesList=mesh.globalNodesList(numberOfElements);
    static vector<vector<float>> globalNodesCoords=mesh.nodesCoordinates(numberOfNodes);
    static int numberOfPointsPerElement;
    if(order==1) numberOfPointsPerElement=4;
    if(order==2) numberOfPointsPerElement=9;
    if(order==3) numberOfPointsPerElement=16;
    if(order==4) numberOfPointsPerElement=25;
    if(order==5) numberOfPointsPerElement=36;
    static vector<double> quadraturePoints=Qk.gaussLobattoQuadraturePoints(order);
    static vector<double> weights=Qk.gaussLobattoQuadratureWeights(order);
    static vector<double> weights2D=Qk.getGaussLobattoWeights(quadraturePoints,
                                                              weights);
    static vector<vector<double>> basisFunction1D=Qk.getBasisFunction1D(order, quadraturePoints);
    static vector<vector<double>> derivativeBasisFunction1D=Qk.getDerivativeBasisFunction1D(order, quadraturePoints);
    static vector<vector<double>> basisFunction2D=Qk.getBasisFunction2D(quadraturePoints,
                                                                        basisFunction1D,
                                                                        basisFunction1D);
    static vector<vector<double>> derivativeBasisFunction2DX=Qk.getBasisFunction2D(quadraturePoints,
                                                                                   derivativeBasisFunction1D,
                                                                                   basisFunction1D);
    static vector<vector<double>> derivativeBasisFunction2DY=Qk.getBasisFunction2D(quadraturePoints,
                                                                                   basisFunction1D,
                                                                                   derivativeBasisFunction1D);


    static vector<float> model=mesh.getModel(numberOfElements);
    static vector<float>massMatrixGlobal(numberOfNodes);
    static vector<float>yGlobal(numberOfNodes);
    int i;
    #pragma omp  parallel for shared(massMatrixGlobal,yGlobal) private(i)
    for ( int i=0; i<numberOfNodes; i++)
    {
	     massMatrixGlobal[i]=0;
	     yGlobal[i]=0;
     }
    // loop over mesh elements
    #pragma omp  parallel shared(model,massMatrixGlobal,yGlobal,pnGlobal,globalNodesCoords,globalNodesList)
    {
    int e;
    int gIndex;
    vector<int>localToGlobal;
    vector<vector<double>>Xi(numberOfPointsPerElement,vector<double>(2,0));
    vector<vector<double>> jacobianMatrix;
    vector<double> detJ;
    vector<vector<double>> invJacobianMatrix;
    vector<vector<double>> transpInvJacobianMatrix;
    vector<vector<double>> B;
    vector<vector<double>> R;
    vector<vector<double>> massMatrixLocal;
    vector<float> pnLocal(numberOfPointsPerElement); 
    vector<float>Y(numberOfPointsPerElement);
    vector<int>neighbors(5);
    //Xi, pnLocal and Y cannot be private: segmentatin fault ????
    #pragma omp for schedule(dynamic) private(e,i,gIndex,localToGlobal,jacobianMatrix,detJ)\
                                      private(invJacobianMatrix,transpInvJacobianMatrix)\
                                      private(B,R,massMatrixLocal,neighbors)
    for ( int e=0; e<numberOfElements; e++)
    {
        // extract global coordinates of element e
        // get local to global indexes of nodes of element e
        localToGlobal=mesh.localToGlobalNodes(e,numberOfPointsPerElement,globalNodesList);  
        
        //get global coordinates Xi of element e
        mesh.getXi(numberOfPointsPerElement, globalNodesCoords, localToGlobal, Xi);
        
	    // compute jacobian Matrix
        jacobianMatrix= Qk.computeJacobianMatrix(numberOfPointsPerElement,Xi,
                                                                        derivativeBasisFunction2DX,
                                                                        derivativeBasisFunction2DY);
        //cout<<"Jacobian marix done"<<endl;
	    // compute determinant of jacobian Matrix
        detJ= Qk.computeDeterminantOfJacobianMatrix(numberOfPointsPerElement,
                                                                   jacobianMatrix);
        //cout<<"detJ Done"<<endl;
	    // compute inverse of Jacobian Matrix
        invJacobianMatrix= Qk.computeInvJacobianMatrix(numberOfPointsPerElement,
                                                                              jacobianMatrix,
                                                                              detJ);
        //cout<<"InvdetJ Done"<<endl;     
	    // compute transposed inverse of Jacobian Matrix
        transpInvJacobianMatrix= Qk.computeTranspInvJacobianMatrix(numberOfPointsPerElement,
                                                                                          jacobianMatrix,
                                                                                          detJ);
        //cout<<"transpInvdetJ Done"<<endl;
	    // compute  geometrical transformation matrix
        B=Qk.computeB(numberOfPointsPerElement, invJacobianMatrix, transpInvJacobianMatrix, detJ);
        //cout<<"computeB Done"<<endl;
	    // compute stifness and mass matrix
        R=Qk.gradPhiGradPhi(numberOfPointsPerElement, weights2D, B, derivativeBasisFunction2DX,
                                                   derivativeBasisFunction2DY);

        // compute local mass matrix
        massMatrixLocal=Qk.phiIphiJ(numberOfPointsPerElement, weights2D, basisFunction2D, detJ);
        //cout<<"compute R and massmatrix done"<<endl;
	    // get pnGlobal to pnLocal 
        //#pragma omp for private(massMatrixLocal,pnLocal)
	    for ( int i=0; i<numberOfPointsPerElement; i++)
        {
            massMatrixLocal[i][i]/=(model[e]*model[e]);
            pnLocal[i]=pnGlobal[localToGlobal[i]][i2];
        }
        //cout <<"update pnLocal Done  "<<numberOfPointsPerElement<<endl;
        //#pragma omp for private(Y)
        for ( int i=0; i<numberOfPointsPerElement; i++)
        {
	       Y[i]=0;
           for ( int j=0; j<numberOfPointsPerElement; j++)
           {
               Y[i]-=R[i][j]*pnLocal[j]; 
           }
        }
        //cout <<"update Y Done"<<endl
        for ( int i=0; i<numberOfPointsPerElement; i++)
        {
	        gIndex=localToGlobal[i];
	        massMatrixGlobal[gIndex]+=massMatrixLocal[i][i];
	        yGlobal[gIndex]+=Y[i];
        }
    neighbors=mesh.neighbors(e);
	}
    }
    // update pressure
    float tmp;
    #pragma omp for schedule(dynamic) private(tmp)
    for ( int i=0 ; i< numberOfNodes; i++)
    {
        tmp=timeSample*timeSample;
        pnGlobal[i][i1]=2*pnGlobal[i][i2]-pnGlobal[i][i1]+tmp*yGlobal[i]/massMatrixGlobal[i];    
    }
}

/// add right and side
void solver::addRightAndSides(const int & timeStep,
                              const int & numberOfRHS,
                              const int & i2,
                              const float & timeSample,
                              vector<vector<float>> & pnGlobal,
                              const vector<vector<float>> & rhsTerm,
                              const vector<vector<float>> & rhsLocation,
                              simpleMesh mesh)
{
    static int numberOfNodes=mesh.getNumberOfNodes();
    static int numberOfElements=mesh.getNumberOfElements();
    static vector<float> model=mesh.getModel(numberOfElements);
    static vector<vector<int>> nodeList=mesh.globalNodesList(numberOfElements);
    int  i, rhsElement;
    float tmp=timeSample*timeSample;
    for ( int i=0; i<numberOfRHS;i++)
    {
        //extract element number for current rhs
        float x=rhsLocation[i][0];
        float y=rhsLocation[i][1];
        int rhsElement=mesh.getElementNumberFromPoints(x,y);
        // compute global node numbe to add source term to
        int nodeRHS=nodeList[rhsElement][0];
        pnGlobal[nodeRHS][i2]+=tmp*model[rhsElement]*model[rhsElement]*rhsTerm[i][timeStep];
    }
}