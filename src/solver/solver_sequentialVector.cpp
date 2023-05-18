#include "solver.hpp"


solver::solver() {}
solver::~solver() {}

// compute one step of the time dynamic wave equation solver
//vector<vector<float>> 
void solver::computeOneStep(const float & timeSample,
		                      const int &order,
					             int &i1,
					             int &i2,
                            vector<vector<float>> & pnGlobal,
					             simpleMesh  mesh,
					             QkGL Qk)
{
   // get infos from mesh
   static int const numberOfNodes=mesh.getNumberOfNodes();
   static int const numberOfElements=mesh.getNumberOfElements();
   static int const numberOfInteriorNodes=mesh.getNumberOfInteriorNodes();
   static vector<vector<int>> const globalNodesList=mesh.globalNodesList(numberOfElements);
   static vector<int>   const listOfInteriorNodes=mesh.getListOfInteriorNodes(numberOfInteriorNodes);
   static vector<vector<float>>  const globalNodesCoords=mesh.nodesCoordinates(numberOfNodes);

   // get model
   static vector<float>  const model=mesh.getModel(numberOfElements);

   //get infos about finite element order of approximation
   int numberOfPointsPerElement;
   if(order==1) numberOfPointsPerElement=4;
   if(order==2) numberOfPointsPerElement=9;
   if(order==3) numberOfPointsPerElement=16;
   if(order==4) numberOfPointsPerElement=25;
   if(order==5) numberOfPointsPerElement=36;

   // get quadrature points and weights
   static vector<double> const quadraturePoints=Qk.gaussLobattoQuadraturePoints(order);
   static vector<double> const weights=Qk.gaussLobattoQuadratureWeights(order);
   static vector<double> const weights2D=Qk.getGaussLobattoWeights(quadraturePoints,
                                                            weights);

   // get basis function and corresponding derivatives
   static vector<vector<double>> const basisFunction1D=Qk.getBasisFunction1D(order, quadraturePoints);
   static vector<vector<double>> const derivativeBasisFunction1D=Qk.getDerivativeBasisFunction1D(order, quadraturePoints);
   static vector<vector<double>> const basisFunction2D=Qk.getBasisFunction2D(quadraturePoints,
                                                                        basisFunction1D,
                                                                        basisFunction1D);
   static vector<vector<double>> const derivativeBasisFunction2DX=Qk.getBasisFunction2D(quadraturePoints,
                                                                                   derivativeBasisFunction1D,
                                                                                   basisFunction1D);
   static vector<vector<double>> const derivativeBasisFunction2DY=Qk.getBasisFunction2D(quadraturePoints,
                                                                                   basisFunction1D,
                                                                                   derivativeBasisFunction1D);

   static vector<float>massMatrixGlobal(numberOfNodes);
   static vector<float>yGlobal(numberOfNodes);

   for ( int i=0; i<numberOfNodes;i++)
   {
	   massMatrixGlobal[i]=0;
	   yGlobal[i]=0;
   }

   // loop over mesh elements
   for (int e=0; e<numberOfElements;e++)
   {
      // extract global coordinates of element e
      // get local to global indexes of nodes of element e
      vector<int> const localToGlobal=mesh.localToGlobalNodes(e,numberOfPointsPerElement,globalNodesList);  
        
      //get global coordinates Xi of element e
      vector<vector<double>> const Xi=mesh.getXi(numberOfPointsPerElement, globalNodesCoords, localToGlobal);
        
	   // compute jacobian Matrix
      vector<vector<double>> const jacobianMatrix= Qk.computeJacobianMatrix(numberOfPointsPerElement,Xi,
                                                derivativeBasisFunction2DX,
                                                derivativeBasisFunction2DY);
	   // compute determinant of jacobian Matrix
      vector<double> const detJ= Qk.computeDeterminantOfJacobianMatrix(numberOfPointsPerElement,
                                                   jacobianMatrix);
	   // compute inverse of Jacobian Matrix
      vector<vector<double>> const invJacobianMatrix= Qk.computeInvJacobianMatrix(numberOfPointsPerElement,
                                                      jacobianMatrix,
                                                      detJ);   
	   // compute transposed inverse of Jacobian Matrix
      vector<vector<double>> const transpInvJacobianMatrix= Qk.computeTranspInvJacobianMatrix(numberOfPointsPerElement,
                                                                  jacobianMatrix,
                                                                  detJ);
	   // compute  geometrical transformation matrix
      vector<vector<double>> const B=Qk.computeB(numberOfPointsPerElement, invJacobianMatrix, transpInvJacobianMatrix, detJ);

	   /**
	   // compute stifness and mass matrix
      vector<vector<double>> const R=Qk.gradPhiGradPhi(numberOfPointsPerElement, weights2D, B, derivativeBasisFunction2DX,
                           derivativeBasisFunction2DY);
      **/
      // compute stifness and mass matrix
      vector<vector<double>> const R=Qk.gradPhiGradPhi(numberOfPointsPerElement,order, weights2D, B, derivativeBasisFunction1D);


      // compute local mass matrix
      vector<double> massMatrixLocal=Qk.phiIphiJ(numberOfPointsPerElement, weights2D, detJ);
	   // get pnGlobal to pnLocal 
      vector<float> pnLocal(numberOfPointsPerElement); 
      vector<float> Y(numberOfPointsPerElement);
	   for ( int i=0; i<numberOfPointsPerElement; i++)
      {
         massMatrixLocal[i]/=(model[e]*model[e]);
         pnLocal[i]=pnGlobal[localToGlobal[i]][i2];
      }
      for ( int i=0; i<numberOfPointsPerElement; i++)
      {
	      Y[i]=0;
        for ( int j=0; j<numberOfPointsPerElement; j++)   
        {
            Y[i]+=R[i][j]*pnLocal[j]; 
         }
      }
      for ( int i=0; i<numberOfPointsPerElement; i++)
     {
        int gIndex=localToGlobal[i];
	     massMatrixGlobal[gIndex]+=massMatrixLocal[i];
	     yGlobal[gIndex]+=Y[i];
     }
   }

   // update pressure
   for (int i=0; i<numberOfInteriorNodes; i++)
  {
         int I=listOfInteriorNodes[i];
         float tmp=timeSample*timeSample;
         pnGlobal[I][i1]=2*pnGlobal[I][i2]-pnGlobal[I][i1]-tmp*yGlobal[I]/massMatrixGlobal[I];    
   }
      
   // damping terms
   // get infos from mesh
   static int const numberOfBoundaryNodes=mesh.getNumberOfBoundaryNodes();
   static int const numberOfBoundaryFaces=mesh.getNumberOfBoundaryFaces();
   static vector<int>  const listOfBoundaryNodes=mesh.getListOfBoundaryNodes(numberOfBoundaryNodes);
   static vector<vector<int>> const faceInfos=mesh.getBoundaryFacesInfos();
   static vector<vector<int>> const localFaceNodeToGlobalFaceNode=mesh.getLocalFaceNodeToGlobalFaceNode();
   static vector<float>ShGlobal(numberOfBoundaryNodes,0);

   for ( int i=0; i<numberOfBoundaryNodes; i++)
   {
	   ShGlobal[i]=0;
   }
   // Note: this loop is data parallel.
   for ( int iFace=0; iFace<numberOfBoundaryFaces; iFace++)
   {
         vector<float>ds(order+1,0);
         vector<float>Sh(order+1,0);
         //get ds
         ds=Qk.computeDs(iFace,order,faceInfos,globalNodesCoords,
                        derivativeBasisFunction2DX,
                        derivativeBasisFunction2DY);
         //compute Sh and ShGlobal
         for (int i=0; i<order+1;i++)
         {
            int gIndexFaceNode=localFaceNodeToGlobalFaceNode[iFace][i];
            Sh[i]=weights[i]*ds[i]/(model[faceInfos[iFace][0]]);
            ShGlobal[gIndexFaceNode]+=Sh[i];
         }
         /**
         cout<<"iFace="<<iFace<<endl;
         for (int i=0; i<order+1;i++)
         {
            int gIndexFaceNode=localFaceNodeToGlobalFaceNode[iFace][i];
            cout<<"   gIndex="<<gIndexFaceNode<<endl;
            cout<<"   Sh["<<i<<"]="<<Sh[i]<<endl;
            cout<<"   ShGlobal["<<gIndexFaceNode<<"]="<<ShGlobal[gIndexFaceNode]<<endl;
         }
         **/    
   }

   // update pressure @ boundaries; 
   for ( int i=0; i<numberOfBoundaryNodes; i++)
   {
      float tmp=timeSample*timeSample;
      int I=listOfBoundaryNodes[i];
      float invMpSh=1/(massMatrixGlobal[I]+timeSample*ShGlobal[i]*0.5);
      float MmSh=massMatrixGlobal[I]-timeSample*ShGlobal[i]*0.5;
      pnGlobal[I][i1]=invMpSh*(2*massMatrixGlobal[I]*pnGlobal[I][i2]-MmSh*pnGlobal[I][i1]-tmp*yGlobal[I]);  
   }
      
   /**
   for ( int i=0 ; i< numberOfBoundaryNodes; i++)
   {
      cout<<"ShGlobal["<<i<<"]="<<ShGlobal[i]<<endl;
   }
   for ( int i=0 ; i< numberOfBoundaryNodes; i++)
   {
      int I=listOfBoundaryNodes[i];
      cout<<"i="<<i<<", BoundaryNode="<<I<<endl;
   }  
   **/
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
