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
namespace grid
{

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
  simpleMesh();
  #ifdef USE_RAJA
  LVARRAY_HOST_DEVICE ~simpleMesh();
  #elif defined USE_KOKKOS
  KOKKOS_FUNCTION ~simpleMesh();
  #else
  ~simpleMesh();
  #endif
  // Returns number of Nodes in the mesh
  int  getNumberOfNodes() const;
  //  Returns the number of elements of the mesh
  int  getNumberOfElements() const;
  // get number of points per element
  int getNumberOfPointsPerElement() const;
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
  // sort element by color
  // red=0, green=1, blue=2, yellow=3
  int getNumberOfElementsByColor() const;
  #ifdef USE_RAJA
  void sortElementsByColor(int  numberOfElementsByColor[] ,arrayInt const & listOfElementsByColor) const;
  #elif defined USE_KOKKOS
  void sortElementsByColor(int  numberOfElementsByColor[] ,arrayInt const & listOfElementsByColor) const;
  #else
  void sortElementsByColor(int  numberOfElementsByColor[] ,arrayInt  & listOfElementsByColor) const;
  #endif
  // Initialize nodal coordinates.
  #ifdef USE_RAJA
  void nodesCoordinates( const int & numberOfNodes, arrayReal const & nodeCoords ) const;
  #elif defined USE_KOKKOS
  void nodesCoordinates( const int & numberOfNodes, arrayReal const & nodeCoords ) const;
  #else
  void nodesCoordinates( const int & numberOfNodes, arrayReal & nodeCoords ) const;
  #endif
  //  list of global nodes ( vertices)
  #ifdef USE_RAJA
  void globalNodesList(const int & numberOfElements, arrayInt const & nodesList ) const;
  #elif defined USE_KOKKOS
  void globalNodesList(const int & numberOfElements, arrayInt const & nodesList ) const;
  #else
  void globalNodesList(const int & numberOfElements, arrayInt & nodesList ) const;
  #endif
  // local to global
  #ifdef USE_RAJA
  LVARRAY_HOST_DEVICE int localToGlobalNodes( const int & elementNumber, 
                                              const int & nPointsPerElement, 
                                              arrayIntView const & nodesList, 
                                              int localToGlobal[]) const;
  #elif defined USE_KOKKOS
  KOKKOS_FUNCTION int localToGlobalNodes( const int & elementNumber, 
                                          const int & nPointsPerElement, 
                                          arrayInt const & nodesList, 
                                          int localToGlobal[]) const;
  #else
  int localToGlobalNodes( const int & elementNumber, 
                          const int & nPointsPerElement, 
                          arrayInt  & nodesList, 
                          int localToGlobal[]) const;
  #endif
  // compute element e where (x,y) belongs to
  int getElementNumberFromPoints( const float & x, const float & y ) const;
  // set model
  #ifdef USE_RAJA
  void getModel( const int & numberOfNodes, vectorReal const & model ) const;
  #elif defined USE_KOKKOS
  void getModel( const int & numberOfNodes, vectorReal const & model ) const;
  #else
  void getModel( const int & numberOfNodes, vectorReal & model ) const;
  #endif
  // list of neighbours of element e
  #ifdef USE_RAJA
  void neighbors( const int & e, vectorInt const & neigh ) const;
  #elif defined USE_KOKKOS
  void neighbors( const int & e, vectorInt const & neigh ) const;
  #else
  void neighbors( const int & e, vectorInt & neigh ) const;
  #endif
  // get global coordinates of element e
  #ifdef USE_RAJA
  LVARRAY_HOST_DEVICE int  getXi( const int & elementNumber,
		                  const int & numberOfPointsPerElement,
				  arrayIntView  const & nodesList,
                                  arrayRealView const & globalNodesCoords,
                                  double  Xi[][2]) const;
  LVARRAY_HOST_DEVICE int  getXi( const int & numberOfPointsPerElement,
                                  arrayRealView const & globalNodesCoords,
                                  int const  localToGlobal[] , 
                                  double  Xi[][2]) const;
  #elif defined USE_KOKKOS
  KOKKOS_FUNCTION int  getXi( const int & elementNumber,
		              const int & numberOfPointsPerElement,
			      arrayInt const  & nodesList,
                              arrayReal const & globalNodesCoords,
                              double  Xi[][2]) const;
  KOKKOS_FUNCTION int  getXi( const int & numberOfPointsPerElement,
                              arrayReal const & globalNodesCoords,
                              int const  localToGlobal[] , 
                              double  Xi[][2]) const;
  #else
  int  getXi( const int & elementNumber,
	      const int & numberOfPointsPerElement,
	      arrayInt  & nodesList,
              arrayReal & globalNodesCoords,
              double  Xi[][2]) const;
  int  getXi( const int & numberOfPointsPerElement,
              arrayReal  & globalNodesCoords,
              int const  localToGlobal[] , 
              double  Xi[][2]) const;
  #endif
  // get global DOF belonging to the faces of element e
  #ifdef USE_RAJA
  LVARRAY_HOST_DEVICE int getGlobalDofOfFace(  const int & e,
                                               arrayIntView  const & globalNodesList,
                                               int const localToGlobal[],
                                               int  nodesFace[][6] ) const;
  #elif defined USE_KOKKOS
  KOKKOS_FUNCTION int getGlobalDofOfFace(  const int & e,
                                           arrayInt  const & globalNodesList,
                                           int const localToGlobal[],
                                           int  nodesFace[][6] ) const;
  #else
  int getGlobalDofOfFace( const int & e,
                          arrayInt  & globalNodesList,
                          int const localToGlobal[],
                          int  nodesFace[][6]) const;
  #endif
  // provides informations about boundary  faces:
  // element number,
  // orientation of the face
  #ifdef USE_RAJA
  void getBoundaryFacesInfos(arrayInt const & faceInfos)const;
  #elif defined USE_KOKKOS
  void getBoundaryFacesInfos(arrayInt const & faceInfos)const;
  #else
  void getBoundaryFacesInfos(arrayInt & faceInfos)const;
  #endif
  // get list of interior Elements
  #ifdef USE_RAJA
  void getListOfInteriorElements(vectorInt const & listOfInteriorElements) const;
  #elif defined USE_KOKKOS
  void getListOfInteriorElements(vectorInt const & listOfInteriorElements) const;
  #else
  void getListOfInteriorElements(vectorInt & listOfInteriorElements) const;
  #endif
  //  get list of global interior nodes
  #ifdef USE_RAJA
  int getListOfInteriorNodes( const int & numberOfInteriorNodes, vectorInt const & listOfInteriorNodes ) const;
  #elif defined USE_KOKKOS
  int getListOfInteriorNodes( const int & numberOfInteriorNodes, vectorInt const & listOfInteriorNodes ) const;
  #else
  int getListOfInteriorNodes( const int & numberOfInteriorNodes, vectorInt & listOfInteriorNodes ) const;
  #endif
  //  get list of global boundary nodes
  #ifdef USE_RAJA
  int getListOfBoundaryNodes( const int & numberOfBoundaryNodes, vectorInt const & listOfBoundaryNodes ) const;
  #elif defined USE_KOKKOS
  int getListOfBoundaryNodes( const int & numberOfBoundaryNodes, vectorInt const & listOfBoundaryNodes ) const;
  #else
  int getListOfBoundaryNodes( const int & numberOfBoundaryNodes, vectorInt & listOfBoundaryNodes ) const;
  #endif
  // provides a mapping between local node of a face and global node Face:
  #ifdef USE_RAJA
  void getLocalFaceNodeToGlobalFaceNode(arrayInt const &localFaceNodeToGlobalFaceNode) const;
  #elif defined USE_KOKKOS
  void getLocalFaceNodeToGlobalFaceNode(arrayInt const &localFaceNodeToGlobalFaceNode) const;
  #else
  void getLocalFaceNodeToGlobalFaceNode(arrayInt &localFaceNodeToGlobalFaceNode) const;
  #endif
  // compute global to local node indes
  int Itoij( const int & I, int & i, int & j ) const;
  // project vector node to grid
  std::vector<std::vector<float>> projectToGrid( const int numberOfNodes, const std::vector<float> inputVector ) const;
};
}
#endif //SIMPLE_MESH_
