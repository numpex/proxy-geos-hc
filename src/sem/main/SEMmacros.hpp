#include "commonConfig.hpp"

// define Macros test SEM proxy 2D case
#if defined (SEM2D)
#define DIMENSION  2
#define ROW 36
#define COL 4
#define ZEROED2D 0
// test SEM proxy 3D case
#else
#define DIMENSION  3
#define ROW 64
#define COL 6
#define ZEROED2D 1
#endif

// create views only for RAJA
#if defined (USE_RAJA)
#define CREATEVIEWS \
        arrayRealView rhsTerm=RHS_Term.toView(); \
        arrayRealView pnGlobal=PN_Global.toView(); \
        vectorIntView rhsElement=RHS_Element.toView(); \
        vectorRealView massMatrixGlobal=this->massMatrixGlobal.toView(); \
        vectorRealView yGlobal=this->yGlobal.toView(); \
        vectorRealView ShGlobal=this->ShGlobal.toView(); \
        arrayIntView globalNodesList=this->globalNodesList.toView(); \
        vectorRealView model=this->model.toView(); \
        vectorIntView listOfInteriorNodes=this->listOfInteriorNodes.toView(); \
        vectorIntView listOfBoundaryNodes=this->listOfBoundaryNodes.toView(); \
        arrayIntView faceInfos=this->faceInfos.toView(); \
        arrayIntView localFaceNodeToGlobalFaceNode=this->localFaceNodeToGlobalFaceNode.toView(); \
        vectorDoubleView weights=this->weights.toView(); \
        arrayRealView globalNodesCoords=this->globalNodesCoords.toView(); \
        arrayDoubleView derivativeBasisFunction1D=this->derivativeBasisFunction1D.toView();
#else
#define CREATEVIEWS
#define RHS_Term rhsTerm
#define PN_Global pnGlobal
#define RHS_Element rhsElement
#endif
