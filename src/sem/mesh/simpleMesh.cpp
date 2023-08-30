#ifndef SIMPLEMESH_HPP_
#define SIMPLEMESH_HPP_

#include <iostream>
#include <cmath>
#include "dataType.hpp"
#include "simpleMesh.hpp"

simpleMesh::simpleMesh( const int & ex_in, const int & ey_in, const float & lx_in, const float & ly_in, const int & order_in )
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
  cout <<"simpleMesh initiliazed\n";
  cout<<"nx, ny="<<nx<<", "<<ny<<endl;
}
simpleMesh::~simpleMesh(){};

int simpleMesh::getNumberOfNodes() const
{
  int numberOfNodes=(ex*order+1)*(ey*order+1);
  return numberOfNodes;
}

int simpleMesh::getNumberOfElements() const
{
  int numberOfElements=ex*ey;
  return numberOfElements;
}

// get number of points per element
int simpleMesh::getNumberOfPointsPerElement() const
{
  return((order+1)*(order+1));
}

//get nx
int simpleMesh::getNx() const
{return nx;}

//get nxy
int simpleMesh::getNy() const
{return ny;}

//get dx
int simpleMesh::getDx() const
{return hx;}

//get dy
int simpleMesh::getDy() const
{return hy;}

//get number of interior elements
int simpleMesh::getNumberOfInteriorElements() const
{return (ex-2)*(ey-2);}

//get number of interior Nodes
int simpleMesh::getNumberOfInteriorNodes() const
{return (nx-2)*(ny-2);}

//get number of boundary Faces
int simpleMesh::getNumberOfBoundaryFaces() const
{return 2*(ex+ey);}

// get number of Boundary nodes
int simpleMesh::getNumberOfBoundaryNodes() const
{return 2*(nx+ny)-4;}

// Initialize nodal coordinates.
//arrayReal simpleMesh::nodesCoordinates( const int & numberOfNodes ) const
void simpleMesh::nodesCoordinates( const int & numberOfNodes, arrayReal & nodeCoords ) const
{
  //arrayReal nodeCoords( numberOfNodes, 2 );
  std::vector<float> coordX( nx );
  std::vector<float> coordY( ny );
  vectorReal xi( order+1 );

  switch( order )
  {
    case 1:
      xi[0]=-1.;
      xi[1]=1.;
      break;
    case 2:
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
  for( int i=0; i<ex; i++ )
  {
    //cout<<"elementx "<<i<<endl;
    float x0=i*hx;
    float x1=(i+1)*hx;
    float b=(x1+x0)/2.;
    float a=b-x0;

    for( int j=0; j<order+1; j++ )
    {
      coordX[j+i*order]=a*xi[j]+b;
      //cout<<coordX[j+i*order]<<" nx= "<<nx<< " j+i*order="<<j+i*order<<", ";
    }
    //cout<<endl;
  }
  for( int i=0; i<ey; i++ )
  {
    //cout<<"elementx "<<i<<endl;
    float y0=i*hy;
    float y1=(i+1)*hy;
    float b=(y1+y0)/2.;
    float a=b-y0;
    for( int j=0; j<order+1; j++ )
    {
      coordY[j+i*order]=a*xi[j]+b;
      //cout<<coordY[j+i*order]<<" ny= "<<ny<< " j+i*order="<<j+i*order<<", ";
    }
    //cout<<endl;
  }
  for( int j=0; j<ny; j++ )
  {
    for( int i=0; i<nx; i++ )
    {
      nodeCoords(i+nx*j,0)=coordX[i];
      nodeCoords(i+nx*j,1)=coordY[j];
      //cout<<"Xi["<<i+nx*j<<"][0]="<<nodeCoords[i+nx*j][0]<<", ";
      //cout<<"Xi["<<i+nx*j<<"][1]="<<nodeCoords[i+nx*j][1]<<endl;
    }
  }
  //return nodeCoords;
}

//  list of global nodes ( vertices) for each element
//arrayInt simpleMesh::globalNodesList( const int & numberOfElements ) const
void simpleMesh::globalNodesList( const int & numberOfElements, arrayInt & nodesList ) const
{
  int nDof=(order+1)*(order+1);
  //arrayInt nodesList( numberOfElements, nDof );

  for( int j=0; j<ey; j++ )
  {
    for( int i=0; i<ex; i++ )
    {
      int n0=i+j*ex;
      int offset=i*order+j*order*nx;
      //cout<<"element "<<n0<<endl;
      for( int k=0; k<order+1; k++ )
      {
        for( int l=0; l<order+1; l++ )
        {
          int dofLocal=l+k*(order+1);
          int dofGlobal=offset+l+k*nx;
          nodesList(n0,dofLocal)=dofGlobal;
          //cout<<nodesList[n0][dofLocal]<<", ";
        }
      }
      //cout<<endl;
    }
  }
  //return nodesList;
}

// local to global
//const vectorInt simpleMesh::localToGlobalNodes( const int & elementNumber, const int & nPointsPerElement, arrayInt & nodesList )const
void simpleMesh::localToGlobalNodes( const int & elementNumber, const int & nPointsPerElement, arrayInt & nodesList,  vectorInt &localToGlobal)const
{
  //vectorInt localToGlobal( nPointsPerElement );
  for( int i=0; i<nPointsPerElement; i++ )
  {
    localToGlobal[i]=nodesList(elementNumber,i);
  }
  //return localToGlobal;
}

// compute global node to grid  indexes
int simpleMesh::Itoij( const int & I, int & i, int & j ) const
{
  i=I%nx;
  j=int((I-i)/nx);
  return 0;
}

// project vector node to grid
std::vector<std::vector<float>> simpleMesh::projectToGrid( const int numberOfNodes, const std::vector<float> inputVector ) const
{
  std::vector<vector<float>> grid( nx,std::vector<float> (ny) );
  int i, j;
  cout<< " In projecToGrid numberOfNodes="<<numberOfNodes<<endl;
  for( int node=0; node<numberOfNodes; node++ )
  {
    Itoij( node, i, j );
    grid[i][j]=inputVector[node];
  }
  return grid;
}

// compute element e where (x,y) belongs to
int simpleMesh::getElementNumberFromPoints( const float & x, const float & y ) const
{
  int eX, eY;
  for( int i=0; i<ex; i++ )
  {
    if( x>= i*hx && x<(i+1)*hx )
      eX=i;
  }
  for( int i=0; i<ey; i++ )
  {
    if( y>= i*hy && y<(i+1)*hy )
      eY=i;
  }
  return eX+eY*ex;
}

// set model
//vectorReal simpleMesh::getModel( const int & numberOfElements ) const
void simpleMesh::getModel( const int & numberOfElements, vectorReal & model ) const
{
  //vectorReal model( numberOfElements );
  for( int i=0; i<numberOfElements; i++ )
  {
    model[i]=1500;
  }
  //return model;
}
// list neighbors for element e
// neigh[0] to neigh[3] list neigbours elements
// neigh[4] if not zero indicate that element e is a boundary element
// what type of boundary element
vectorInt simpleMesh::neighbors( const int & e ) const
{
  vectorInt neigh( 5 );
  int i=e%ex;
  int j=int((e-i)/ex);
  neigh[0]=e-1;
  neigh[1]=e-ex;
  neigh[2]=e+1;
  neigh[3]=e+ex;
  neigh[4]=0;
  // if left boundary element
  if( i==0 )
  {
    neigh[0]=-1;
    neigh[4]=-1;
  }
  // if right boundary element
  if( i==ex-1 )
  {
    neigh[2]=-2;
    neigh[4]=-3;
  }
  // if bottom boundary element
  if( j==0 )
  {
    neigh[1]=-1;
    neigh[4]=-2;
  }
  // if top boundary element
  if( j==ey-1 )
  {
    neigh[3]=-2;
    neigh[4]=-4;
  }
  return neigh;

}

// get global coordinates of element e
//arrayDouble simpleMesh::getXi( const int & numberOfPointsPerElement, arrayReal & globalNodesCoords,
//                              vectorInt & localToGlobal ) const
void simpleMesh::getXi( const int & numberOfPointsPerElement, arrayReal & globalNodesCoords,
                     vectorInt & localToGlobal , arrayDouble & Xi) const
{
  //arrayDouble Xi( numberOfPointsPerElement, 2 );
  for( int i=0; i<numberOfPointsPerElement; i++ )
  {
    Xi(i,0)=globalNodesCoords(localToGlobal[i],0);
    Xi(i,1)=globalNodesCoords(localToGlobal[i],1);
  }
  //return Xi;
}


// get global DOF belonging to the faces of element e
//    _____3_____
//   |           |
//   |           |
// 0 |           | 2
//   |           |
//   |______1____|
//
void simpleMesh::getGlobalDofOfFace( const int & e,
                                         arrayInt  & globalNodesList,
                                         vectorInt & localToGlobal,
                                         arrayInt  & nodesFace ) const
{
  //arrayInt nodesFace( 4, order+1 );

  //left face
  for( int i=0; i<order+1; i++ )
  {
    int dofLocal=i*(order+1);
    nodesFace(0,i)=globalNodesList(e,dofLocal);
  }
  //bottom face
  for( int i=0; i<order+1; i++ )
  {
    int dofLocal=i;
    nodesFace(1,i)=globalNodesList(e,dofLocal);
  }
  //right face
  for( int i=0; i<order+1; i++ )
  {
    int dofLocal=order+i*(order+1);
    nodesFace(2,i)=globalNodesList(e,dofLocal);
  }
  //top face
  for( int i=0; i<order+1; i++ )
  {
    int dofLocal=i+order*(order+1);
    nodesFace(3,i)=globalNodesList(e,dofLocal);
  }
  //return nodesFace;
}

// provides informations about boundary  faces:
// element number,
// orientation of the face
// list of global indexes
// this method is sequential only for omp !!!
//arrayInt simpleMesh::getBoundaryFacesInfos() const
void  simpleMesh::getBoundaryFacesInfos(arrayInt & faceInfos) const
{
  //int numberOfBoundaryFaces=getNumberOfBoundaryFaces();
  int numFace=0;
  //arrayInt faceInfos( numberOfBoundaryFaces, 2+(order+1));
  // bottom, j=0, l=0
  for( int i=0; i<ex; i++ )
  {
    numFace=i;
    faceInfos(numFace,0)=i;
    faceInfos(numFace,1)=1;
    int offset=i*order;
    for( int j=0; j<order+1; j++ )
    {
      faceInfos(numFace,2+j)=offset+j;
    }
  }
  // right i=ex-1 l=order
  for( int j=0; j<ey; j++ )
  {
    int e=ex-1+j*ex;
    numFace=ex+j;
    int offset=(ex-1)*order+j*order*nx;
    faceInfos(numFace,0)=e;
    faceInfos(numFace,1)=2;
    for( int k=0; k<order+1; k++ )
    {
      faceInfos(numFace,2+k)=offset+order+k*nx;
    }
  }
  // top j=ey-1, k=order
  for( int i=0; i<ex; i++ )
  {
    int e=i+(ey-1)*ex;
    numFace=ey+ex+i;
    faceInfos(numFace,0)=e;
    faceInfos(numFace,1)=3;
    int offset=i*order+(ey-1)*order*nx;
    int k=order;
    for( int l=0; l<order+1; l++ )
    {
      faceInfos(numFace,2+l)=offset+l+k*nx;
    }
  }
  // left, i=0 and l=0
  for( int j=0; j<ey; j++ )
  {
    int e=j*ex;
    numFace=2*ex+ey+j;
    int offset=j*order*nx;
    faceInfos[numFace][0]=e;
    faceInfos[numFace][1]=0;
    for( int k=0; k<order+1; k++ )
    {
      faceInfos[numFace][2+k]=offset+k*nx;
    }
  }
/**
   for ( int iFace=0;iFace<numberOfBoundaryFaces;iFace++)
   {
      for ( int i=0; i<order+1;i++)
      {
         cout<<"iFace="<<iFace<<" node number="<<faceInfos[iFace][2+i]<<endl;
      }
   }
 **/
  //return faceInfos;
}

// get list of interior Elements
void simpleMesh::getListOfInteriorElements(vectorInt & listOfInteriorElements) const
{
  //int numberOfInteriorElements=getNumberOfInteriorElements();
  //vectorInt listOfInteriorElements( numberOfInteriorElements );
  int k=0;
  for( int j=1; j<ey-1; j++ )
  {
    for( int i=1; i<ex-1; i++ )
      listOfInteriorElements[k]=i+j*ex;
    k++;
  }
  //return listOfInteriorElements;
}

//  get list of global interior nodes
//vectorInt simpleMesh::getListOfInteriorNodes( const int & numberOfInteriorNodes ) const
void simpleMesh::getListOfInteriorNodes( const int & numberOfInteriorNodes,
                                         vectorInt & listOfInteriorNodes ) const
{
  //vectorInt listOfInteriorNodes( numberOfInteriorNodes );
  int k=0;
  for( int j=1; j<ny-1; j++ )
  {
    for( int i=1; i<nx-1; i++ )
    {
      listOfInteriorNodes[k]=i+j*nx;
      k++;
    }
  }
  //return listOfInteriorNodes;
}

//  get list of global boundary nodes
//vectorInt simpleMesh::getListOfBoundaryNodes( const int & numberOfBoundaryNodes ) const
void simpleMesh::getListOfBoundaryNodes( const int & numberOfBoundaryNodes, vectorInt & listOfBoundaryNodes ) const
{
  //vectorInt listOfBoundaryNodes( numberOfBoundaryNodes );
  //cout<<"numberOfBoundaryNOdes="<<numberOfBoundaryNodes<<endl;
  int k=0;
  //bottom
  int j=0;
  for( int i=0; i<nx; i++ )
  {
    listOfBoundaryNodes[k]=i+j*nx;
    k++;
    //cout<<"k="<<k<<endl;
  }
  //right
  int i=nx-1;
  for( int j=1; j<ny; j++ )
  {
    listOfBoundaryNodes[k]=i+j*nx;
    k++;
    //cout<<"k="<<k<<endl;
  }
  //top
  j=ny-1;
  for( int i=0; i<nx-1; i++ )
  {
    listOfBoundaryNodes[k]=i+j*nx;
    k++;
    //cout<<"k="<<k<<endl;
  }
  //left
  i=0;
  for( int j=1; j<ny-1; j++ )
  {
    listOfBoundaryNodes[k]=i+j*nx;
    k++;
    //cout<<"k="<<k<<endl;
  }
  //for( int j=0;j<numberOfBoundaryNodes;j++)
  //{
  //  cout<<"j="<<j<<", "<<listOfBoundaryNodes[j]<<endl;
  //}
  //return listOfBoundaryNodes;
}

// provides a mapping between local node of a face and global node Face:
//arrayInt simpleMesh::getLocalFaceNodeToGlobalFaceNode() const
void simpleMesh::getLocalFaceNodeToGlobalFaceNode(arrayInt &localFaceNodeToGlobalFaceNode) const
{
  //int numberOfBoundaryFaces=getNumberOfBoundaryFaces();
  int numFace=0;
  int offset;
  //arrayInt localFaceNodeToGlobalFaceNode( numberOfBoundaryFaces, order+1 );
  // bottom, j=0, l=0
  for( int i=0; i<ex; i++ )
  {
    numFace=i;
    offset=i*order;
    for( int j=0; j<order+1; j++ )
    {
      localFaceNodeToGlobalFaceNode(numFace,j)=offset+j;
    }
  }
  // right i=ex-1 l=order
  for( int j=0; j<ey; j++ )
  {
    int e=ex-1+j*ex;
    numFace=ex+j;
    int offset=localFaceNodeToGlobalFaceNode(numFace-1,order);
    for( int k=0; k<order+1; k++ )
    {
      localFaceNodeToGlobalFaceNode(numFace,k)=offset+k;
    }
  }
  // top j=ey-1, k=order
  for( int i=0; i<ex; i++ )
  {
    int e=i+(ey-1)*ex;
    numFace=ey+ex+i;
    if( i<ex-1 )
    {
      if( i==0 )
      {
        offset=localFaceNodeToGlobalFaceNode(numFace-1,order)+1;
      }
      else
      {
        offset=localFaceNodeToGlobalFaceNode(numFace-1,order);
      }
      for( int l=0; l<order+1; l++ )
      {
        localFaceNodeToGlobalFaceNode(numFace,l)=offset+l;
      }
    }
    else if( i==ex-1 )
    {
      offset=localFaceNodeToGlobalFaceNode(numFace-1,order);
      for( int l=0; l<order; l++ )
      {
        localFaceNodeToGlobalFaceNode(numFace,l)=offset+l;
      }
      localFaceNodeToGlobalFaceNode(numFace,order)=localFaceNodeToGlobalFaceNode(ex+ey-1,order);
    }

  }
  // left, i=0 and l=0
  for( int j=0; j<ey; j++ )
  {
    numFace=2*ex+ey+j;
    if( j==0 )
    {
      localFaceNodeToGlobalFaceNode(numFace,0)=0;
      for( int k=1; k<order+1; k++ )
      {
        offset=localFaceNodeToGlobalFaceNode(numFace-1,order-1);
        localFaceNodeToGlobalFaceNode(numFace,k)=offset+k;
      }
    }
    else if( j>0 && j<ey-1 )
    {
      for( int k=0; k<order+1; k++ )
      {
        offset=localFaceNodeToGlobalFaceNode(numFace-1,order);
        localFaceNodeToGlobalFaceNode(numFace,k)=offset+k;
      }
    }
    else if( j==ey-1 )
    {
      for( int k=0; k<order; k++ )
      {
        offset=localFaceNodeToGlobalFaceNode(numFace-1,order);
        localFaceNodeToGlobalFaceNode(numFace,k)=offset+k;
      }
      offset=localFaceNodeToGlobalFaceNode(ex+ey,0);
      localFaceNodeToGlobalFaceNode(numFace,order)=offset;
    }

  }
  /**for ( int iFace=0;iFace<numberOfBoundaryFaces;iFace++)
     {
     for ( int i=0; i<order+1;i++)
     {
        cout<<"iFace="<<iFace<<" local number="<<i<<" global number="<<localFaceNodeToGlobalFaceNode[iFace][i]<<endl;
     }
     }**/
  //return localFaceNodeToGlobalFaceNode;
}
#endif //SIMPLEMESH_HPP_
