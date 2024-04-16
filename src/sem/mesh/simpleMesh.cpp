#ifndef SIMPLEMESH_HPP_
#define SIMPLEMESH_HPP_

#include "simpleMesh.hpp"

namespace grid
{
simpleMesh::simpleMesh( const int & ex_in, const int & ey_in, const int & ez_in,  
		        const float & lx_in, const float & ly_in,const float & lz_in,
                        const int & order_in)
{
  orderx=order_in;
  ordery=order_in;
  orderz=order_in;
  order=order_in;
  ex=ex_in;
  ey=ey_in;
  ez=ez_in;
  lx=lx_in;
  ly=ly_in;
  lz=lz_in;
  nx=ex_in*orderx+1;
  ny=((ey==0)?1:ey_in*ordery+1);
  nz=ez_in*orderz+1;
  hx=lx/(1.*ex);
  hy=((ey==0)?1:ly/(1.*ey));
  hz=lz/(1.*ez);
  cout <<"simpleMesh initiliazed\n";
  cout<<"ex, ey, ez="<<ex<<", "<<ey<<", "<<ez<<endl;
  cout<<"nx, ny, nz="<<nx<<", "<<ny<<", "<<nz<<endl;
}

int simpleMesh::getNumberOfNodes() const
{
  int numberOfNodes=(ex*orderx+1)*(ey*ordery+1)*(ez*orderz+1);
  printf("number of nodes %d\n",numberOfNodes);
  return numberOfNodes;
}

int simpleMesh::getNumberOfElements() const
{
  int numberOfElements=((ey==0)?ex*ez:ex*ey*ez);
  printf("number of element %d\n",numberOfElements);
  return numberOfElements;
}

// get number of points per element
int simpleMesh::getNumberOfPointsPerElement() const
{
  return(((ny==1)?(orderx+1)*(orderz+1):(orderx+1)*(ordery+1)*(orderz+1)));
}

//get number of interior elements
int simpleMesh::getNumberOfInteriorElements() const
{return ((ny==1)?(ex-2)*(ez-2):(ex-2)*(ey-2)*(ez-2));}

//get number of interior Nodes
int simpleMesh::getNumberOfInteriorNodes() const
{return ((ny==1)?(nx-2)*(nz-2):(nx-2)*(ny-2)*(nz-2));}

//get nx
int simpleMesh::getNx() const
{return nx;}

//get ny
int simpleMesh::getNy() const
{return ny;}

//get nz
int simpleMesh::getNz() const
{return nz;}

//get dx
int simpleMesh::getDx() const
{return hx;}

//get dy
int simpleMesh::getDy() const
{return hy;}

//get dz
int simpleMesh::getDz() const
{return hz;}

// Initialize nodal coordinates.
void simpleMesh::nodesCoordinates( const int & numberOfNodes, arrayReal const & nodeCoords ) const
{
  std::vector<float> coordX( nx );
  std::vector<float> coordY( ny );
  std::vector<float> coordZ( nz );
  std::vector<float> xi( order+1 );
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
    float x0=i*hx;
    float x1=(i+1)*hx;
    float b=(x1+x0)/2.;
    float a=b-x0;
    for( int j=0; j<order+1; j++ )
    {
      coordX[j+i*order]=a*xi[j]+b;
    }
  }
  for( int i=0; i<ey; i++ )
  {
    float y0=i*hy;
    float y1=(i+1)*hy;
    float b=(y1+y0)/2.;
    float a=b-y0;
    for( int j=0; j<order+1; j++ )
    {
      coordY[j+i*order]=a*xi[j]+b;
    }
  }
  for( int i=0; i<ez; i++ )
  {
    float z0=i*hz;
    float z1=(i+1)*hz;
    float b=(z1+z0)/2.;
    float a=b-z0;
    for( int j=0; j<order+1; j++ )
    {
      coordZ[j+i*order]=a*xi[j]+b;
    }
  }
  //printf("nx ny,nz %d %d %d\n",nx,ny,nz);
  for( int k=0; k<ny; k++ )
  {
    for( int j=0; j<nz; j++ )
    {
      for( int i=0; i<nx; i++ )
      {
        nodeCoords(i+nx*j+k*nx*nz,0)=coordX[i];
        nodeCoords(i+nx*j+k*nx*nz,1)=coordZ[j];
        nodeCoords(i+nx*j+k*nx*nz,2)=coordY[k];
	//printf("i=%d,j=%d,k=%d, n=%d, X %lf,Z %lf,Y %lf\n",i,j,k,i+nx*j+k*nx*nz,coordX[i],coordZ[j],coordY[k]);
      }
    }
  }
}

//  list of global nodes ( vertices) for each element
void simpleMesh::globalNodesList( const int & numberOfElements, arrayInt const & nodesList ) const
{
   for( int j=0; j<((ey==0)?1:ey); j++ )
   {
     for( int k=0; k<ez; k++ )
     {
       for( int i=0; i<ex; i++ )
       {
          int n0=i+k*ex+j*ex*ez;
          int offset=i*order+k*order*nx+j*order*nx*nz;
          for( int m=0; m<((ey==0)?1:order+1); m++ )
          {
             for( int n=0; n<order+1; n++ )
             {
               for( int l=0; l<order+1; l++ )
               {
                 int dofLocal=l+n*(order+1)+m*(order+1)*(order+1);
                 int dofGlobal=offset+l+n*nx+m*nx*nz;
                 nodesList(n0,dofLocal)=dofGlobal;
		 //if(n0==7)printf("offset=%d dofLocal=%d dofGlobal=%d\n",offset,dofLocal,dofGlobal);
               }
             }
          }
        }
      }
   }
}

// local to global
PROXY_HOST_DEVICE int simpleMesh::localToGlobalNodes( const int & elementNumber, 
                                                        const int & nPointsPerElement, 
                                                        arrayIntView const & nodesList,  
                                                        int   localToGlobal[])const
{
  for( int i=0; i<nPointsPerElement; i++ )
  {
    localToGlobal[i]=nodesList(elementNumber,i);
  }
  return 0;
}

// compute global node to grid  indexes
int simpleMesh::Itoijk( const int & I, int & i, int & j, int & k ) const
{
  i=I%nx;
  j=int((I%(nx*nz))/nx);
  k=I/(nx*nz);
  return 0;
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
  std::vector<vector<float>> grid( nx,std::vector<float> (nz) );
  int i, j;
  cout<< " In projecToGrid numberOfNodes="<<numberOfNodes<<endl;
  for( int node=0; node<numberOfNodes; node++ )
  {
     Itoij( node, i, j );
     grid[i][j]=inputVector[node];
  }
  return grid;
}

// compute element e where (x,y,z) belongs to
int simpleMesh::getElementNumberFromPoints( const float & x, const float & y, const float & z ) const
{
  int eX=0, eY=0, eZ=0;
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
  for( int i=0; i<ez; i++ )
  {
    if( z>= i*hz && z<(i+1)*hz )
      eZ=i;
  }
  return eX+eZ*ex+eY*ex*ez;
}

// set model
void simpleMesh::getModel( const int & numberOfElements, vectorReal const & model ) const
{
  for( int i=0; i<numberOfElements; i++ )
  {
    model[i]=1500;
  }
}

//  get list of global interior nodes
int simpleMesh::getListOfInteriorNodes( const int & numberOfInteriorNodes,
                                         vectorInt const & listOfInteriorNodes ) const
{
  int m=0;
  if(ny==1)
  {
    for( int j=1; j<nz-1; j++ )
    {
      for( int i=1; i<nx-1; i++ )
      {
        listOfInteriorNodes[m]=i+j*nx;
        m++;
      }
    }
  }
  else
  {
    for( int k=1; k<ny-1; k++ )
    {
      for( int j=1; j<nz-1; j++ )
      {
        for( int i=1; i<nx-1; i++ )
        {
          listOfInteriorNodes[m]=i+j*nx+k*nx*nz;
          m++;
        }
      }
    }
  }
  return 0;
}

// get list of interior Elements
void simpleMesh::getListOfInteriorElements(vectorInt const & listOfInteriorElements) const
{
  int m=0;
  if(ey==0)
  {
    for( int j=1; j<ez-1; j++ )
    {
      for( int i=1; i<ex-1; i++ )
      {
        listOfInteriorElements[m]=i+j*ex;
        m++;
      }
    }
  }
  else
  {
    for( int k=1; k<ey-1; k++ )
    {
      for( int j=1; j<ez-1; j++ )
      {
        for( int i=1; i<ex-1; i++ )
        {
          listOfInteriorElements[m]=i+j*ex+k*ex*ez;
          m++;
        }
      }
    }
  }
}

// sort element by color
// red=0, green=1, blue=2, yellow=3
// get
int simpleMesh::getNumberOfElementsByColor() const
{return ((ey==0)?(ex/2+ex%2)*(ez/2+ez%2):(ex/2+ex%2)*(ey/2+ey%2)*(ez/2+ez%2));}

//sort
void simpleMesh::sortElementsByColor(int  numberOfElementsByColor[] ,arrayInt const & listOfElementsByColor) const
{
  // red
  int k=0;
  for ( int j=0;j<ez; j=j+2)
  {
      for ( int i=0;i<ex; i=i+2)
      {
          listOfElementsByColor(0,k)=i+j*ex;
          k=k+1;
      }
  }
  numberOfElementsByColor[0]=k;
  // green
  k=0;
  for ( int j=0;j<ez; j=j+2)
  {
      for ( int i=1;i<ex; i=i+2)
      {
          listOfElementsByColor(1,k)=i+j*ex;
          k=k+1;
      }
  }
  numberOfElementsByColor[1]=k;
  // blue
  k=0;
  for ( int j=1;j<ez; j=j+2)
  {
      for ( int i=0;i<ex; i=i+2)
      {
          listOfElementsByColor(2,k)=i+j*ex;
          k=k+1;
      }
  }
  numberOfElementsByColor[2]=k;
  // yellow
  k=0;
  for ( int j=1;j<ez; j=j+2)
  {
      for ( int i=1;i<ex; i=i+2)
      {
          listOfElementsByColor(3,k)=i+j*ex;
          k=k+1;
      }
  }
  numberOfElementsByColor[3]=k;
}


//get number of boundary Faces
int simpleMesh::getNumberOfBoundaryFaces() const
{return ((ey==0)?2*(ex+ez):2*(ex*ey+ey*ez+ex*ez));}

// get number of Boundary nodes
int simpleMesh::getNumberOfBoundaryNodes() const
{return (ey==0?2*(nx+nz)-4:2*((nx-1)*(nz-1)+(nx-1)*(ny-1)+(ny-1)*(ny-1)));}

// get global DOF belonging to the faces of element e
//    _____3_____
//   |           |
//   |           |
// 0 |           | 2
//   |           |
//   |______1____|
//
PROXY_HOST_DEVICE int simpleMesh::getGlobalDofOfFace( const int & e,
                                                        arrayIntView  const & globalNodesList,
                                                        int const localToGlobal[],
                                                        int  nodesFace[][6] ) const
{
  //left face
  for( int i=0; i<order+1; i++ )
  {
    int dofLocal=i*(order+1);
    nodesFace[0][i]=globalNodesList(e,dofLocal);
  }
  //bottom face
  for( int i=0; i<order+1; i++ )
  {
    int dofLocal=i;
    nodesFace[1][i]=globalNodesList(e,dofLocal);
  }
  //right face
  for( int i=0; i<order+1; i++ )
  {
    int dofLocal=order+i*(order+1);
    nodesFace[2][i]=globalNodesList(e,dofLocal);
  }
  //top face
  for( int i=0; i<order+1; i++ )
  {
    int dofLocal=i+order*(order+1);
    nodesFace[3][i]=globalNodesList(e,dofLocal);
  }
  return 0;
}

// list of global indexes
// this method is sequential only for omp !!!
void  simpleMesh::getBoundaryFacesInfos(arrayInt const & faceInfos) const
{
  int numFace=0;
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
  for( int j=0; j<ez; j++ )
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
  // top j=ez-1, k=order
  for( int i=0; i<ex; i++ )
  {
    int e=i+(ez-1)*ex;
    numFace=ez+ex+i;
    faceInfos(numFace,0)=e;
    faceInfos(numFace,1)=3;
    int offset=i*order+(ez-1)*order*nx;
    int k=order;
    for( int l=0; l<order+1; l++ )
    {
      faceInfos(numFace,2+l)=offset+l+k*nx;
    }
  }
  // left, i=0 and l=0
  for( int j=0; j<ez; j++ )
  {
    int e=j*ex;
    numFace=2*ex+ez+j;
    int offset=j*order*nx;
    faceInfos(numFace,0)=e;
    faceInfos(numFace,1)=0;
    for( int k=0; k<order+1; k++ )
    {
      faceInfos(numFace,2+k)=offset+k*nx;
    }
  }
}

//  get list of global boundary nodes
int simpleMesh::getListOfBoundaryNodes( const int & numberOfBoundaryNodes, vectorInt const & listOfBoundaryNodes ) const
{
  int k=0;
  //bottom
  int j=0;
  for( int i=0; i<nx; i++ )
  {
    listOfBoundaryNodes[k]=i+j*nx;
    k++;
  }
  //right
  int i=nx-1;
  for( int j=1; j<nz; j++ )
  {
    listOfBoundaryNodes[k]=i+j*nx;
    k++;
  }
  //top
  j=nz-1;
  for( int i=0; i<nx-1; i++ )
  {
    listOfBoundaryNodes[k]=i+j*nx;
    k++;
  }
  //left
  i=0;
  for( int j=1; j<nz-1; j++ )
  {
    listOfBoundaryNodes[k]=i+j*nx;
    k++;
  }
  return 0;
}

// provides a mapping between local node of a face and global node Face:
void simpleMesh::getLocalFaceNodeToGlobalFaceNode(arrayInt const &localFaceNodeToGlobalFaceNode) const
{
  int numFace=0;
  int offset;
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
  for( int j=0; j<ez; j++ )
  {
    //int e=ex-1+j*ex;
    numFace=ex+j;
    int offset=localFaceNodeToGlobalFaceNode(numFace-1,order);
    for( int k=0; k<order+1; k++ )
    {
      localFaceNodeToGlobalFaceNode(numFace,k)=offset+k;
    }
  }
  // top j=ez-1, k=order
  for( int i=0; i<ex; i++ )
  {
    numFace=ez+ex+i;
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
      localFaceNodeToGlobalFaceNode(numFace,order)=localFaceNodeToGlobalFaceNode(ex+ez-1,order);
    }

  }
  // left, i=0 and l=0
  for( int j=0; j<ez; j++ )
  {
    numFace=2*ex+ez+j;
    if( j==0 )
    {
      localFaceNodeToGlobalFaceNode(numFace,0)=0;
      for( int k=1; k<order+1; k++ )
      {
        offset=localFaceNodeToGlobalFaceNode(numFace-1,order-1);
        localFaceNodeToGlobalFaceNode(numFace,k)=offset+k;
      }
    }
    else if( j>0 && j<ez-1 )
    {
      for( int k=0; k<order+1; k++ )
      {
        offset=localFaceNodeToGlobalFaceNode(numFace-1,order);
        localFaceNodeToGlobalFaceNode(numFace,k)=offset+k;
      }
    }
    else if( j==ez-1 )
    {
      for( int k=0; k<order; k++ )
      {
        offset=localFaceNodeToGlobalFaceNode(numFace-1,order);
        localFaceNodeToGlobalFaceNode(numFace,k)=offset+k;
      }
      offset=localFaceNodeToGlobalFaceNode(ex+ez,0);
      localFaceNodeToGlobalFaceNode(numFace,order)=offset;
    }

  }
}

}
#endif //SIMPLEMESH_HPP_
