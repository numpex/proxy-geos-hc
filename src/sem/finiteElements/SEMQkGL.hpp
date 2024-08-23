#ifndef SEMQKGL_HPP_
#define SEMQKGL_HPP_

#include "dataType.hpp"
#include "SEMmacros.hpp"
#include "SEMdata.hpp"

using namespace std;

/**
 * This class is the basis class for the hexahedron finite element cells with shape functions defined on Gauss-Lobatto quadrature points.
 */
class SEMQkGL
{
private:
  int order;
  struct SEMinfo infos;

  ////////////////////////////////////////////////////////////////////////////////////
  //  from GEOS implementation
  /////////////////////////////////////////////////////////////////////////////////////
  constexpr static double sqrt5 = 2.2360679774997897;
  // order of polynomial approximation
  //constexpr static  int r=2;
  constexpr static  int r=SEMinfo::myOrderNumber;
  // number of support/quadrature/nodes points in one direction
  constexpr static int numSupport1dPoints=r+1;
  constexpr static int num1dNodes=numSupport1dPoints;
  // Half the number of support points, rounded down. Precomputed for efficiency
  constexpr static int halfNodes = ( numSupport1dPoints - 1 )/ 2;
  // the number of nodes/support points per element
  constexpr static int numSupportPoints=(r+1)*(r+1)*(r+1);

public:
  PROXY_HOST_DEVICE SEMQkGL(){};
  PROXY_HOST_DEVICE ~SEMQkGL(){};

  

  // get Gauss Lobatto quadrature points
  void gaussLobattoQuadraturePoints( int order, vectorDouble const & quadraturePoints ) const;

  // get Gauss Lobatto quadrature weights
  void  gaussLobattoQuadratureWeights( int order, vectorDouble const & weights ) const; 
  // compute  1d shape Functions and derivatives
  vector< double >  shapeFunction1D( int order, double xi ) const;
  vector< double >  derivativeShapeFunction1D( int order, double xi ) const;

  // get  1d shape Functions and derivatives for all quadrature points
  void getDerivativeBasisFunction1D( int order, vectorDouble const & quadraturePoints,
                                     arrayDouble const & derivativeBasisFunction1D ) const;

  
  // V1
  // compute B and M
  PROXY_HOST_DEVICE void computeB( const int & elementNumber,
                                   const int & order,
                                   VECTOR_DOUBLE_VIEW const & weights,
                                   ARRAY_INT_VIEW const & nodesList,
                                   ARRAY_REAL_VIEW const & nodesCoords,
                                   ARRAY_DOUBLE_VIEW const & dPhi,
                                   float massMatrixLocal[],
                                   float B[][COL] ) const;
  
  // V2
  // compute B and M
  PROXY_HOST_DEVICE void computeB( const int & elementNumber,
                                   const int & order,
                                   VECTOR_DOUBLE_VIEW const & weights,
                                   ARRAY_REAL_VIEW const & nodesCoordsX,
                                   ARRAY_REAL_VIEW const & nodesCoordsY,
                                   ARRAY_REAL_VIEW const & nodesCoordsZ,
                                   ARRAY_DOUBLE_VIEW const & dPhi,
                                   float massMatrixLocal[],
                                   float B[][COL] ) const;
  
  // compute the matrix $R_{i,j}=\int_{K}{\nabla{\phi_i}.\nabla{\phi_j}dx}$
  // Marc Durufle Formulae
  PROXY_HOST_DEVICE void gradPhiGradPhi( const int & nPointsPerElement,
                                         const int & order,
                                         VECTOR_DOUBLE_VIEW const & weights,
                                         ARRAY_DOUBLE_VIEW const & dPhi,
                                         float const B[][COL],
                                         float const pnLocal[],
                                         float R[],
                                         float Y[] ) const;

  // V1
  // compute stiffnessVector.
  // returns mass matrix and stiffness vector local to an element
  PROXY_HOST_DEVICE void computeMassMatrixAndStiffnessVector(const int & elementNumber,
                                                             const int & order,
                                                             const int & nPointsPerElement,
                                                             ARRAY_INT_VIEW const & nodesList,
                                                             ARRAY_REAL_VIEW const & nodesCoords,
                                                             VECTOR_DOUBLE_VIEW const & weights,
                                                             ARRAY_DOUBLE_VIEW const & dPhi,
                                                             float massMatrixLocal[],
                                                             float const pnLocal[],
                                                             float Y[]) const;

  // V2
  // compute stiffnessVector.
  // returns mass matrix and stiffness vector local to an element
  PROXY_HOST_DEVICE void computeMassMatrixAndStiffnessVector(const int & elementNumber,
                                                             const int & order,
                                                             const int & nPointsPerElement,
                                                             ARRAY_REAL_VIEW const & nodesCoordsX,
                                                             ARRAY_REAL_VIEW const & nodesCoordsY,
                                                             ARRAY_REAL_VIEW const & nodesCoordsZ,
                                                             VECTOR_DOUBLE_VIEW const & weights,
                                                             ARRAY_DOUBLE_VIEW const & dPhi,
                                                             float massMatrixLocal[],
                                                             float const pnLocal[],
                                                             float Y[]) const;

  //computeDs
  PROXY_HOST_DEVICE void computeDs( const int & iFace,
                                    const int & order,
                                    ARRAY_INT_VIEW const & faceInfos,
                                    int numOfBasisFunctionOnFace,
                                    float Js[][6],
                                    ARRAY_REAL_VIEW const & globalNodesCoords,
                                    ARRAY_DOUBLE_VIEW const & dPhi,
                                    float ds[] ) const;


  /////////////////////////////////////////////////////////////////////////////////////
  //  from GEOS implementation
  /////////////////////////////////////////////////////////////////////////////////////
  //constexpr static double parentSupportCoord( const int order, const int supportPointIndex ) 
  template<typename SEMinfo>
  PROXY_HOST_DEVICE
  constexpr static double parentSupportCoord(  const int supportPointIndex ) 
  {
      double result=0.0;
      switch( SEMinfo::myOrderNumber )
      {
        case 1:
          return -1.0 + 2.0 * (supportPointIndex & 1);
        case 2:
          switch( supportPointIndex )
          {
            case 0:
              return -1.0;
              break;
            case 2:
              return 1.0;
            case 1:
            default:
              return 0.0;
          }
        case 3:
          switch( supportPointIndex )
          {
            case 0:
               result = -1.0;
               break;
            case 1:
              result = -1.0/sqrt5;
              break;
            case 2:
              result = 1.0/sqrt5;
              break;
            case 3:
              result = 1.0;
              break;
            default:
              break;
          }
        default:
           return 0;
      }
      return result;
   }


  /**
   * @brief The gradient of the basis function for a support point evaluated at
   *  a given support point. By symmetry, p is assumed to be in 0, ..., (N-1)/2
   * @param q The index of the basis function
   * @param p The index of the support point
   * @return The gradient of basis function.
  */
  template<typename SEMinfo>
  PROXY_HOST_DEVICE  
  constexpr static double gradientAt( const int q, const int p ) 
  {
      switch( SEMinfo::myOrderNumber )
      {
        case 1:
          return q == 0 ? -0.5 : 0.5;
        case 2:
          switch( q )
          {
            case 0:
              return p == 0 ? -1.5 : -0.5;
            case 1:
              return p == 0 ? 2.0 : 0.0;
            case 2:
              return p == 0 ? -0.5 : 0.5;
            default:
              return 0;
          }
        case 3:
          switch( q )
          {
            case 0:
              return p == 0 ? -3.0 : -0.80901699437494742410;
            case 1:
              return p == 0 ? 4.0450849718747371205 : 0.0;
            case 2:
              return p == 0 ? -1.5450849718747371205 : 1.1180339887498948482;
            case 3:
              return p == 0 ? 0.5 : -0.30901699437494742410;
            default:
              return 0;
          }
        default:
           return 0;
      }
   }


  /*
   * @brief Compute the 1st derivative of the q-th 1D basis function at quadrature point p
   * @param q the index of the 1D basis funcion
   * @param p the index of the 1D quadrature point
   * @return The derivative value
  */
  PROXY_HOST_DEVICE 
  constexpr static double basisGradientAt( const int order, const int q, const int p )
  {
     if( p <= halfNodes )
     {
       return gradientAt<SEMinfo>( q, p );
     }
     else
     {
       return -gradientAt<SEMinfo>( numSupport1dPoints - 1 - q, numSupport1dPoints - 1 - p );
     }
  }

  /**
   * @brief The value of the weight for the given support point
   * @param q The index of the support point
   * @return The value of the weight
  */
  template<typename SEMinfo>
  PROXY_HOST_DEVICE
  constexpr static double weight( const int q )
  {
      switch(SEMinfo::myOrderNumber)
      {
        case 1:
           return 1;
        case 2:
           switch( q )
           {
             case 0:
             case 2:
               return 1.0/3.0;
             default:
               return 4.0/3.0;
           }
        case 3:
           switch( q )
           {
             case 1:
             case 2:
               return 5.0/6.0;
             default:
              return 1.0/6.0;
           }
       default:
           return 0;
      }
   }

  /**
   * @brief Calculates the linear index for support/quadrature points from ijk
   *   coordinates.
   * @param r order of polynomial approximation
   * @param i The index in the xi0 direction (0,r)
   * @param j The index in the xi1 direction (0,r)
   * @param k The index in the xi2 direction (0,r)
   * @return The linear index of the support/quadrature point (0-(r+1)^3)
  */
  PROXY_HOST_DEVICE 
  constexpr static int linearIndex( const int r,
                                    const int i,
                                    const int j,
                                    const int k ) 
  {
           return i + (r+1) * j + (r+1)*(r+1) * k;
  }

  /**
   * @brief Calculate the Cartesian/TensorProduct index given the linear index
   *   of a support point.
   * @param linearIndex The linear index of support point
   * @param r order of polynomial approximation
   * @param i0 The Cartesian index of the support point in the xi0 direction.
   * @param i1 The Cartesian index of the support point in the xi1 direction.
   * @param i2 The Cartesian index of the support point in the xi2 direction.
  */
  PROXY_HOST_DEVICE 
  constexpr static void multiIndex( int const r,int const linearIndex, int & i0, int & i1, int & i2 ) 
  {
        i2 = linearIndex/((r+1)*(r+1));
        i1 = (linearIndex%((r+1)*(r+1)))/(r+1);
        i0 = (linearIndex%((r+1)*(r+1)))%(r+1);
  }

  /**
   * @brief Compute the interpolation coefficients of the q-th quadrature point in a given direction
   * @param q the index of the quadrature point in 1D
   * @param k the index of the interval endpoint (0 or 1)
   * @return The interpolation coefficient
  */
  PROXY_HOST_DEVICE 
  constexpr static double interpolationCoord( const int order, const int q, const int k ) 
  {
      const double alpha = (parentSupportCoord<SEMinfo>( q ) + 1.0 ) / 2.0;
      return k == 0 ? ( 1.0 - alpha ) : alpha;
   }

  /**
   * @brief Compute the 1D factor of the coefficient of the jacobian on the q-th quadrature point,
   * with respect to the k-th interval endpoint (0 or 1). The computation depends on the position
   * in the basis tensor product of this term (i, equal to 0, 1 or 2) and on the direction in which
   * the gradient is being computed (dir, from 0 to 2)
   * @param q The index of the quadrature point in 1D
   * @param i The index of the position in the tensor product
   * @param k The index of the interval endpoint (0 or 1)
   * @param dir The direction in which the derivatives are being computed
   * @return The value of the jacobian factor
  */
  PROXY_HOST_DEVICE 
  constexpr static double jacobianCoefficient1D( const int order, const int q, const int i, const int k, const int dir )
  {
      if( i == dir )
      {
        return k== 0 ? -1.0/2.0 : 1.0/2.0;
      }
      else
      {
        return interpolationCoord( order, q, k );
      }
  }

  PROXY_HOST_DEVICE static double determinant(double  m[3][3]); //const;

  /**
   * @brief Invert the symmetric matrix @p srcSymMatrix and store the result in @p dstSymMatrix.
   * @param dstSymMatrix The 3x3 symmetric matrix to write the inverse to.
   * @param srcSymMatrix The 3x3 symmetric matrix to take the inverse of.
   * @return The determinant.
   * @note @p srcSymMatrix can contain integers but @p dstMatrix must contain floating point values.
  */
  PROXY_HOST_DEVICE static void symInvert( double  dstSymMatrix[6], double  srcSymMatrix[6] );// const;

  /**
   * @brief Invert the symmetric matrix @p symMatrix overwritting it.
   * @param symMatrix The 3x3 symmetric matrix to take the inverse of and overwrite.
   * @return The determinant.
   * @note @p symMatrix can contain integers but @p dstMatrix must contain floating point values.
  */
  PROXY_HOST_DEVICE static void symInvert0( double  symMatrix[6] ); //const;

  /**
   * @brief Calculates the isoparametric "Jacobian" transformation
   *  matrix/mapping from the parent space to the physical space.
   * @param qa The 1d quadrature point index in xi0 direction (0,1)
   * @param qb The 1d quadrature point index in xi1 direction (0,1)
   * @param qc The 1d quadrature point index in xi2 direction (0,1)
   * @param X Array containing the coordinates of the mesh support points.
   * @param J Array to store the Jacobian transformation.
  */
  PROXY_HOST_DEVICE static void jacobianTransformation(  int e, int const r,
                              int const qa,
                              int const qb,
                              int const qc,
                              double const (&X)[8][3],
                              double ( & J )[3][3] ); //const;

  /**
   * @brief computes the non-zero contributions of the d.o.f. indexd by q to the
   *   mass matrix M, i.e., the superposition matrix of the shape functions.
   * @param q The quadrature point index
   * @param X Array containing the coordinates of the mesh support points.
   * @return The diagonal mass term associated to q
  */
   PROXY_HOST_DEVICE static double computeMassTerm( int e, int const r,int const q, double const (&X)[8][3] ); //const;

  /**
   * @brief Calculates the isoparametric "geometrical" transformation
   *  matrix/mapping from the parent space to the physical space.
   * @param qa The 1d quadrature point index in xi0 direction (0,1)
   * @param qb The 1d quadrature point index in xi1 direction (0,1)
   * @param qc The 1d quadrature point index in xi2 direction (0,1)
   * @param X Array containing the coordinates of the mesh support points.
   * @param J Array to store the Jacobian transformation.
   * @param B Array to store the  the geometrical symetic matrix=detJ*J^{-1}J^{-T}.
  */
  PROXY_HOST_DEVICE static void computeBMatrix( int e, int const r,
                      int const qa,
                      int const qb,
                      int const qc,
                      double const (&X)[8][3],
                      double (& J)[3][3],
                      double (& B)[6] ); //const;

  /**
   * @brief compute  gradPhi*B*gradPhi
   * returns value implemented into func
   */
  template<typename FUNC>
  PROXY_HOST_DEVICE
  static void computeGradPhiBGradPhi( int const e,
                               int const r,
                               int const qa,
                               int const qb,
                               int const qc,
                               double const (&B)[6],
                               FUNC && func ); //const;

  /**
   * @brief compute  stiffnessTerm
   * returns mass matrix and stiffness vector local to an element
   */
  template<typename FUNC>
  PROXY_HOST_DEVICE 
  static void computeStiffnessTerm( int e,
                             int r,
                             int const q,
                             double const (&X)[8][3],
                             FUNC && func ); //const;

  /**
   * @brief compute  mass Matrix stiffnessVector.
   * returns mass matrix and stiffness vector local to an element
   */
  PROXY_HOST_DEVICE static void computeMassMatrixAndStiffnessVector(const int & elementNumber,
                                                          const int & order,
                                                          const int & nPointsPerElement,
                                                          ARRAY_REAL_VIEW const & nodesCoordsX,
                                                          ARRAY_REAL_VIEW const & nodesCoordsY,
                                                          ARRAY_REAL_VIEW const & nodesCoordsZ,
                                                          float massMatrixLocal[],
                                                          float pnLocal[],
                                                          float Y[]); //const;

  ////////////////////////////////////////////////////////////////////////////////////
  // END OF GEOS FUNCTIONS
  //////////////////////////////////////////////////////////////////////////////////
};

#endif //SEMQKGL_HPP_
