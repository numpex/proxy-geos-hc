#ifndef SIMPLE_MESH_
#define SIMPLE_MESH_

#include    <iostream>
#include    <vector>
#include    <cmath>
#include    "dataType.hpp"

using    namespace std;

/*
 *  simple 2D quadrangle element mesh
 */

class simpleMesh
{
private:
  int ex, ey;
  int nx, ny;
  float lx, ly;
  float hx, hy;
  int order;
  int nbFaces;
public:

  simpleMesh( const int & ex_in, const int & ey_in, const float & lx_in, const float & ly_in, const int & order_in );
  ~simpleMesh();

  // Returns number of Nodes in the mesh
  int  getNumberOfNodes() const;

  //  Returns the number of elements of the mesh
  int  getNumberOfElements() const;

  //get nx
  int getNx() const;

  //get ny
  int getNy() const;

  //get Dx
  int getDx()const;

  //get Dy
  int getDy() const;

  //get number of interior elements
  int getNumberOfInteriorElements() const;

  //get number of interior elements
  int getNumberOfInteriorNodes() const;

  // get number of Boundary Faces
  int getNumberOfBoundaryFaces() const;

  // get number of Boundary nodes
  int getNumberOfBoundaryNodes() const;

  // Initialize nodal coordinates.
  arrayReal nodesCoordinates( const int & numberOfNodes ) const;

  //  list of global nodes ( vertices)
  arrayInt globalNodesList( const int & numberOfElements ) const;

  // local to global
  const vectorInt localToGlobalNodes( const int & elementNumber, const int & nPointsPerElement, arrayInt & nodesList ) const;
  // compute global to local node indes
  int Itoij( const int & I, int & i, int & j ) const;

  // project vector node to grid
  arrayReal projectToGrid( const int numberOfNodes, vectorReal inputVector ) const;

  // compute element e where (x,y) belongs to
  int getElementNumberFromPoints( const float & x, const float & y ) const;

  // set model
  vectorReal getModel( const int & numberOfNodes ) const;

  // list of neighbours of element e
  vectorInt neighbors( const int & e ) const;

  // get global coordinates of element e
  arrayDouble getXi( const int & numberOfPointsPerElement, arrayReal & globalNodesCoords,
                     vectorInt & localToGlobal ) const;

  // get global DOF belonging to the faces of element e
  arrayInt getGlobalDofOfFace( const int & e,
                               arrayInt & globalNodesList,
                               vectorInt & localToGlobal ) const;

  // provides informations about boundary  faces:
  // element number,
  // orientation of the face
  arrayInt getBoundaryFacesInfos()const;

  // get list of interior Elements
  vectorInt getListOfInteriorElements() const;

  //  get list of global interior nodes
  vectorInt getListOfInteriorNodes( const int & numberOfInteriorNodes ) const;

  //  get list of global boundary nodes
  vectorInt getListOfBoundaryNodes( const int & numberOfBoundaryNodes ) const;

  // provides a mapping between local node of a face and global node Face:
  arrayInt getLocalFaceNodeToGlobalFaceNode() const;


};
#endif //SIMPLE_MESH_
