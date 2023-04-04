#ifndef SIMPLEMESH_HPP_
#define SIMPLEMESH_HPP_

#include <iostream>
#include <vector>
#include <cmath>
#include "simpleMesh.hpp"

using namespace std;

simpleMesh::simpleMesh(const int &ex_in, const int &ey_in, const float &lx_in, const float &ly_in, const int &order_in)
{
	order=order_in;
   ex=ex_in;
	ey=ey_in;
	lx=lx_in;
	ly=ly_in;
	nx=ex_in*order+1;
	ny=ey_in*order+1;
   hx=lx/(1.*ex);
   hy=ly/(1.*ey);
}
simpleMesh::~simpleMesh(){};

int simpleMesh::getNumberOfNodes()
{
	 int numberOfNodes=(ex*order+1)*(ey*order+1);
	 return numberOfNodes;
}

int  simpleMesh::getNumberOfElements()
{
   int numberOfElements=ex*ey;
   return numberOfElements;
}

//get nx
int simpleMesh::getNx()
{return nx;}
//get nxy
int simpleMesh::getNy()
{return ny;}

//get dx
int simpleMesh::getDx()
{return hx;}
//get dy
int simpleMesh::getDy()
{return hy;}

// Initialize nodal coordinates.
vector<vector<float>> simpleMesh::nodesCoordinates(const int &numberOfNodes)
{
      vector<vector<float>> nodeCoords(numberOfNodes,vector<float>(2,0));
      vector<float>coordX(nx);
      vector<float>coordY(ny);
      vector<float>xi(order+1);
         
      switch (order)
      {
         case 1 :
            xi[0]=-1.;
            xi[1]=1.;
            break;
         case 2 :
            xi[0]=-1.;
            xi[1]=0.;
            xi[2]=1.;
            break;
         case 3:
          static constexpr double sqrt5 = 2.2360679774997897;
          xi[0] = -1.0;
          xi[1] = -1./sqrt5;
          xi[2] = 1./sqrt5;
          xi[3] = 1.;
          break;
        case 4:
          static constexpr double sqrt3_7 = 0.6546536707079771;
          xi[0] = -1.0;
          xi[1] = -sqrt3_7;
          xi[2] = 0.0;
          xi[3] = sqrt3_7;
          xi[4] = 1.0;
          break;
        case 5:
          static constexpr double sqrt__7_plus_2sqrt7__ = 3.50592393273573196;
          static constexpr double sqrt__7_mins_2sqrt7__ = 1.30709501485960033;
          static constexpr double sqrt_inv21 = 0.218217890235992381;
          xi[0] = -1.0;
          xi[1] = -sqrt_inv21*sqrt__7_plus_2sqrt7__;
          xi[2] = -sqrt_inv21*sqrt__7_mins_2sqrt7__;
          xi[3] = sqrt_inv21*sqrt__7_mins_2sqrt7__;
          xi[4] = sqrt_inv21*sqrt__7_plus_2sqrt7__;
          xi[5] = 1.0;
          break;
         default:
            break;
      }
      for (int i=0;i<ex;i++)
      {
         //cout<<"elementx "<<i<<endl;
         float x0=i*hx;
         float x1=(i+1)*hx;
         float b=(x1+x0)/2.;
         float a=b-x0;

         for ( int j=0;j<order+1;j++)
         {
            coordX[j+i*order]=a*xi[j]+b;
            //cout<<coordX[j+i*order]<<" nx= "<<nx<< " j+i*order="<<j+i*order<<", ";
            }
         //cout<<endl;
      }  
      for (int i=0;i<ey;i++)
      {
         //cout<<"elementx "<<i<<endl;
         float y0=i*hy;
         float y1=(i+1)*hy;
         float b=(y1+y0)/2.;
         float a=b-y0;
         for ( int j=0;j<order+1;j++)
         {
            coordY[j+i*order]=a*xi[j]+b;
            //cout<<coordY[j+i*order]<<" ny= "<<ny<< " j+i*order="<<j+i*order<<", ";
         }
         //cout<<endl;
      }  
      for ( int j=0; j<ny;j++)
      {
         for (int i=0; i<nx;i++)
         {
            nodeCoords[i+nx*j][0]=coordX[i];
            nodeCoords[i+nx*j][1]=coordY[j];
		      //cout<<"Xi["<<i+nx*j<<"][0]="<<nodeCoords[i+nx*j][0]<<", ";
		      //cout<<"Xi["<<i+nx*j<<"][1]="<<nodeCoords[i+nx*j][1]<<endl;
         }
      }
     return nodeCoords;
}

//  list of global nodes ( vertices) for each element
vector<vector<int>> simpleMesh::globalNodesList(const int &numberOfElements)
{
   int nDof=(order+1)*(order+1);  
   vector<vector<int>> nodesList(numberOfElements,vector<int>(nDof));
   
   for( int j=0;j<ey;j++)
   {
      for ( int i=0; i<ex; i++)
      {
         int n0=i+j*ex;
	      int offset=i*order+j*order*nx;
         //cout<<"element "<<n0<<endl;
         for (int k=0;k<order+1;k++)
         {
            for (int l=0; l<order+1;l++)
            {
               int dofLocal=l+k*(order+1);
               int dofGlobal=offset+l+k*nx;
               nodesList[n0][dofLocal]=dofGlobal;
               //cout<<nodesList[n0][dofLocal]<<", ";
            }
         }
	      //cout<<endl;
      }
   } 
   return nodesList;
}

// local to global
vector<int> simpleMesh::localToGlobalNodes(const int &elementNumber, const int & nPointsPerElement, const vector<vector<int>> &nodesList)
{
  vector<int> localToGlobal(nPointsPerElement,0);
  for ( int i=0; i<nPointsPerElement; i++)
  {
      localToGlobal[i]=nodesList[elementNumber][i];
  }
  return localToGlobal;
}

// compute global node to grid  indexes
int simpleMesh::Itoij(const int &I, int &i, int &j)
{
    i=I%nx;
    j=int((I-i)/nx);
    return 0;
}

// project vector node to grid
vector<vector<float>> simpleMesh::projectToGrid(const int numberOfNodes,const vector<float>inputVector)
{
   vector<vector<float>>grid(nx,vector<float>(ny,0));
   int i,j;
   cout<< " In projecToGrid numberOfNodes="<<numberOfNodes<<endl;
   for ( int node=0; node<numberOfNodes;node++)
   {
      Itoij(node,i,j);
      grid[i][j]=inputVector[node];
   }
   return grid;
}

// compute element e where (x,y) belongs to
int simpleMesh::getElementNumberFromPoints(const float &x,const float &y)
{
   int eX,eY;
   for ( int i=0; i<ex; i++)
   {
      if ( x>= i*hx && x<(i+1)*hx)
         eX=i ;
   }
   for ( int i=0; i<ey; i++)
   {
      if ( y>= i*hy && y<(i+1)*hy)
         eY=i ;
   }
   return eX+eY*ex;
}

// set model
vector<float>simpleMesh::getModel(const int & numberOfElements)
{
   vector<float>model(numberOfElements);
   for ( int i=0; i<numberOfElements;i++)
   {
      model[i]=1500;
   }
   return model;
}
// list neighbors for element e
vector<int> simpleMesh::neighbors(const int & e)
{
   vector<int> neigh(5,0);
   int i=e%ex;
   int j=int((e-i)/ex);
   neigh[0]=e-1;
   neigh[1]=e-ex;
   neigh[2]=e+1;
   neigh[3]=e+ex;
   if(i==0)
   {
      neigh[0]=-1;
      neigh[4]=-1;
   }
   if(i==ex-1)
   {  
      neigh[2]=-2;
      neigh[4]=-3;
   }
   if(j==0)
   {
      neigh[1]=-1;
      neigh[4]=-2;
   }
   if(j==ey-1)
   {
      neigh[3]=-2;
      neigh[4]=-4;
   }
   if(i==0 && j==0)
   {
      neigh[4]=-5;
   }
   if(i==ex-1 && j==0)
   {  
      neigh[4]=-6;
   }
   if(i==0 && j==ey-1)
   {
      neigh[4]=-8;
   }
   if(i==ex-1 && j==ey-1)
   {
      neigh[4]=-7;
   }
   return neigh;

<<<<<<< HEAD
=======
}

// get global coordinates of element e
void simpleMesh::getXi(const int & numberOfPointsPerElement,const vector<vector<float>> & globalNodesCoords,
                       const vector<int> & localToGlobal, vector<vector<double>> & Xi)
{
   for ( int i=0; i<numberOfPointsPerElement; i++)
   {
       Xi[i][0]=globalNodesCoords[localToGlobal[i]][0];
       Xi[i][1]=globalNodesCoords[localToGlobal[i]][1];
       //cout<<" node "<<i<<"  "<<Xi[i][0]<<", "<<Xi[i][1]<<endl;
   }
>>>>>>> ompNew
}
#endif //SIMPLEMESH_HPP_
