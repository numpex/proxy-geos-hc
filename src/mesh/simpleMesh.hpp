#ifndef SIMPLE_MESH_
#define SIMPLE_MESH_

#include    <iostream>
#include    <vector>
#include    <cmath>

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
  vector< vector< float > > nodesCoordinates( const int & numberOfNodes ) const;

  //  list of global nodes ( vertices)
  vector< vector< int > > globalNodesList( const int & numberOfElements ) const;

  // local to global
  const vector< int > localToGlobalNodes( const int & elementNumber, const int & nPointsPerElement, vector< vector< int > >  const & nodesList ) const;
  // compute global to local node indes
  int Itoij( const int & I, int & i, int & j ) const;

  // project vector node to grid
  vector< vector< float > > projectToGrid( const int numberOfNodes, const vector< float > inputVector ) const;

  // compute element e where (x,y) belongs to
  int getElementNumberFromPoints( const float & x, const float & y ) const;

  // set model
  vector< float >getModel( const int & numberOfNodes ) const;

  // list of neighbours of element e
  vector< int > neighbors( const int & e ) const;

  // get global coordinates of element e
  vector< vector< double > > getXi( const int & numberOfPointsPerElement, const vector< vector< float > > & globalNodesCoords,
                                    const vector< int > & localToGlobal ) const;

  // get global DOF belonging to the faces of element e
  vector< vector< int > > getGlobalDofOfFace( const int & e,
                                              const vector< vector< int > > & globalNodesList,
                                              const vector< int > & localToGlobal ) const;

  // provides informations about boundary  faces:
  // element number,
  // orientation of the face
  vector< vector< int > > getBoundaryFacesInfos()const;

  // get list of interior Elements
  vector< int > getListOfInteriorElements() const;

  //  get list of global interior nodes
  vector< int > getListOfInteriorNodes( const int & numberOfInteriorNodes ) const;

  //  get list of global boundary nodes
  vector< int > getListOfBoundaryNodes( const int & numberOfBoundaryNodes ) const;

  // provides a mapping between local node of a face and global node Face:
  vector< vector< int > >getLocalFaceNodeToGlobalFaceNode() const;


};
#endif //SIMPLE_MESH_
