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
       int ex,ey;
       int nx,ny;
       float lx,ly;
       float hx,hy;
       int order;
     public:

       simpleMesh(const int & ex_in, const int & ey_in, const float & lx_in, const float & ly_in,const int & order_in);
       ~simpleMesh();

       // Returns number of Nodes in the mesh
       int  getNumberOfNodes();

       //  Returns the number of elements of the mesh
       int  getNumberOfElements();

       // Initialize nodal coordinates.
       vector<vector<float>> nodesCoordinates(const int &numberOfNodes);

       //  list of global nodes ( vertices)
       vector<vector<int>> globalNodesList(const int &numberOfElements);

       // local to global
       vector<int> localToGlobalNodes(const int &elementNumber, const int & nPointsPerElement,  const vector<vector<int>> & nodesList);

       // compute global to local node indes
       int Itoij(const int &I, int &i, int &j);

       // compute element e where (x,y) belongs to
       int getElementNumberFromPoints(const float &x,const float &y);

       // set model
       vector<float>getModel(const int & numberOfNodes);

       // list neighbors
       vector<int> neighbors(const int & e);
};
#endif //SIMPLE_MESH_
