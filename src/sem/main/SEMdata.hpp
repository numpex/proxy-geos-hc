struct SEMinfo
{
  // get infos from mesh
  int numberOfNodes;
  int numberOfElements;
  int numberOfPointsPerElement;
  int numberOfInteriorNodes;
  int numberOfBoundaryNodes;
  int numberOfBoundaryFaces;

  const int myNumberOfRHS=1;
  const int myOrderNumber=2;
  const float myTimeStep=0.001;
  const int nPointsPerElement = pow((myOrderNumber+1), DIMENSION );

  const float f0=10.;
  const float myTimeMax=1.;
  const int sourceOrder=1;

  int myNumSamples=myTimeMax/myTimeStep;
  int myElementSource;

  #ifdef SEM_MESHCOLOR
  int numberMaxOfElementsByColor;
  const int numberOfColors=4;
  int numberOfElementsByColor[4];
  #endif
};
