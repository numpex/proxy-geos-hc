struct SEMmeshinfo
{
  // get infos from mesh
  int numberOfNodes;
  int numberOfElements;
  int numberOfPointsPerElement;
  int numberOfInteriorNodes;
  int numberOfBoundaryNodes;
  int numberOfBoundaryFaces;
  int numberMaxOfElementsByColor;

  const int myNumberOfRHS=1;
  const int myOrderNumber=3;
  const float myTimeStep=0.001;
  const int nPointsPerElement = pow((myOrderNumber+1), DIMENSION );
  const int numberOfColors=4;
  int numberOfElementsByColor[4];
};
