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
   // get infos from mesh
   static int numberOfNodes=mesh.getNumberOfNodes();
   static int numberOfElements=mesh.getNumberOfElements();
   static int numberOfInteriorNodes=mesh.getNumberOfInteriorNodes();
   static int numberOfBoundaryNodes=mesh.getNumberOfBoundaryNodes();
   static int numberOfBoundaryFaces=mesh.getNumberOfBoundaryFaces();

   static vector<vector<int>> globalNodesList=mesh.globalNodesList(numberOfElements);
   static vector<int> listOfInteriorNodes=mesh.getListOfInteriorNodes(numberOfInteriorNodes);

   static vector<int> listOfBoundaryNodes=mesh.getListOfBoundaryNodes(numberOfBoundaryNodes);
   static vector<int> listOfBoundaryLocalNodes=mesh.getListOfBoundaryLocalNodes(numberOfBoundaryNodes);
   static vector<vector<float>> globalNodesCoords=mesh.nodesCoordinates(numberOfNodes);
   static vector<vector<int>>faceInfos=mesh.getBoundaryFacesInfos();

   
   

   // get model
   static vector<float> model=mesh.getModel(numberOfElements);

   //get infos about finite element order of approximation
   static int numberOfPointsPerElement;
   if(order==1) numberOfPointsPerElement=4;
   if(order==2) numberOfPointsPerElement=9;
   if(order==3) numberOfPointsPerElement=16;
   if(order==4) numberOfPointsPerElement=25;
   if(order==5) numberOfPointsPerElement=36;

   // get quadrature points and weights
   static vector<double> quadraturePoints=Qk.gaussLobattoQuadraturePoints(order);
   static vector<double> weights=Qk.gaussLobattoQuadratureWeights(order);
   static vector<double> weights2D=Qk.getGaussLobattoWeights(quadraturePoints,
                                                            weights);

   // get basis function and corresponding derivatives
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


    
   static vector<float>massMatrixGlobal(numberOfNodes);
   static vector<float>yGlobal(numberOfNodes);

   static vector<float>ShGlobal(numberOfBoundaryNodes,0);

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
      vector<int>neighboursOfElementE(5);
      vector<vector<int>> nodesFace(4,vector<int>(order+1,0));
      //Xi, pnLocal and Y cannot be private: segmentatin fault ????
      #pragma omp for schedule(dynamic) private(e,i,gIndex,localToGlobal,jacobianMatrix,detJ)\
                                      private(invJacobianMatrix,transpInvJacobianMatrix)\
                                      private(B,R,massMatrixLocal,neighboursOfElementE,nodesFace)
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
	     for ( int i=0; i<numberOfPointsPerElement; i++)
        {
            massMatrixLocal[i][i]/=(model[e]*model[e]);
            pnLocal[i]=pnGlobal[localToGlobal[i]][i2];
        }
        //cout <<"update pnLocal Done  "<<numberOfPointsPerElement<<endl;
        for ( int i=0; i<numberOfPointsPerElement; i++)
        {
	       Y[i]=0;
           for ( int j=0; j<numberOfPointsPerElement; j++)
           {
               Y[i]+=R[i][j]*pnLocal[j]; 
           }
        }
        //cout <<"update Y Done"<<endl
        for ( int i=0; i<numberOfPointsPerElement; i++)
        {
	        gIndex=localToGlobal[i];
	        massMatrixGlobal[gIndex]+=massMatrixLocal[i][i];
	        yGlobal[gIndex]+=Y[i];
        }
      }
      // update pressure
      float tmp;
      //cout<< "number of interior nodes="<<numberOfInteriorNodes<<endl;
      #pragma omp for schedule(dynamic) private(tmp)
      for ( int i=0 ; i< numberOfInteriorNodes; i++)
      {
         tmp=timeSample*timeSample;
         int I=listOfInteriorNodes[i];
         //cout<<"i="<<i<<" interior global node="<<I<<endl;
         pnGlobal[I][i1]=2*pnGlobal[I][i2]-pnGlobal[I][i1]-tmp*yGlobal[I]/massMatrixGlobal[I];    
      }

      // damping terms
      #pragma omp  parallel for shared(ShGlobal) private(i)
      for ( int i=0; i<numberOfBoundaryNodes; i++)
      {
	      ShGlobal[i]=0;
      }
      // Note: this loop is data parallel.
      #pragma omp  parallel shared(ShGlobal)
      {
         int face=-1;
         int iFace;
         int gIndexFaceNode=0;
         float xi=0;
         float yi=0;
         vector<int>numOfBasisFunctionOnFace(order+1,0);
         vector<vector<float>>Js(2,vector<float>(order+1,0));
         vector<float>ds(order+1,0);
         vector<float>Sh(order+1,0);
         #pragma omp  for schedule(static) private( iFace,xi,yi,face,gIndexFaceNode)
         for (iFace=0; iFace< numberOfBoundaryFaces; iFace++)
         {
            face=faceInfos[iFace][1];
            cout<<"iFace="<<iFace<<" face= "<<face<<endl;
            // get basis functions on Boundary faces
            switch (face)
            {
               case 0: // left
                  //cout<<"left face "<<endl;
                  for (int i=0;i<order+1;i++)
                  {
                     numOfBasisFunctionOnFace[i]=i*(order+1);
                     //cout<<"num of basis function="<<numOfBasisFunctionOnFace[i]<<", ";
                  }
                  //cout<<endl;
                  break;
               case 1: // bottom
                  //cout<<"bottom face"<<endl;
                  for (int i=0;i<order+1;i++)
                  {
                     numOfBasisFunctionOnFace[i]=i;
                     //cout<<"num of basis function="<<numOfBasisFunctionOnFace[i]<<", ";
                  }
                  //cout<<endl;
                  break;
               case 2: //right
                  //cout<<"right face "<<endl;
                  for (int i=0;i<order+1;i++)
                  {
                     numOfBasisFunctionOnFace[i]=order+i*(order+1);
                     //cout<<"num of basis function="<<numOfBasisFunctionOnFace[i]<<", ";
                  }
                  //cout<<endl;
                  break;
               case 3: //top
                  //cout<<"top face "<<endl; 
                  for (int i=0;i<order+1;i++)
                  {
                     numOfBasisFunctionOnFace[i]=i+order*(order+1);
                     //cout<<"num of basis function="<<numOfBasisFunctionOnFace[i]<<", ";
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
               //cout<<"iface="<<iFace<<", j="<<j<<endl;
               //cout<<"------------"<<endl;
               for (int i=0; i<order+1; i++)
               {
                  xi=globalNodesCoords[faceInfos[iFace][2+i]][0];
                  yi=globalNodesCoords[faceInfos[iFace][2+i]][1];
                  //cout<<"i="<<i<<", facinfo  element "<<faceInfos[iFace][0]<<endl;
                  //cout<<"i="<<i<<", facinfo  face "<<faceInfos[iFace][0]<<endl;
                  //cout<<"i="<<i<<", facinfo  numpoint "<<faceInfos[iFace][2+i]<<endl;
                  cout<<"i="<<i<<", face="<<face<<" xi,yi="<<xi<<", "<<yi<<endl;
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
               cout<<"iFace="<<iFace<<", ds["<<j<<"]="<<ds[j]<<endl;
               cout<<"iFace="<<iFace<<", Js[0]["<<j<<"]="<<Js[0][j]<<endl;
            }
            // compute damping matrix
            for (int i=0; i<order+1;i++)
            {
               /*
               float tmp=0;
               for (int r=0;r<order+1;r++)
               {
                  tmp+=weights[r]*ds[r];
               }
               */
               Sh[i]=weights[i]*ds[i]/(model[faceInfos[iFace][0]]);
               gIndexFaceNode=iFace*order+i;
               //cout<<"gIndexFaceNode="<<gIndexFaceNode<<endl;
               ShGlobal[gIndexFaceNode]+=Sh[i];
            }
         }
      } 
      // update pressure @ boundaries 
      float invMpSh;
      float MmSh;
      #pragma omp for schedule(dynamic) private(tmp,invMpSh,MmSh)
      for ( int i=0 ; i< numberOfBoundaryNodes; i++)
      {
         //cout<<i<<" "<<listOfBoundaryNodes[i]<<endl;
         tmp=timeSample*timeSample;
         int I=listOfBoundaryNodes[i];
         invMpSh=1/(massMatrixGlobal[I]+timeSample*ShGlobal[i]*0.5);
         MmSh=massMatrixGlobal[I]-timeSample*ShGlobal[i]*0.5;
         pnGlobal[I][i1]=invMpSh*(2*massMatrixGlobal[I]*pnGlobal[I][i2]-MmSh*pnGlobal[I][i1]-tmp*yGlobal[I]);   
      }
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