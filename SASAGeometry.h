// Mark Lewis 11-10-11

#ifndef SASAGEOMETRY_H
#define SASAGEOMETRY_H

#include <eigen2/Eigen/Core>
#include <eigen2/Eigen/LU>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <string>
#include <list>

#include "sasa_transformMatrix.h"

class SASAGeometry
{
public:
//methods
  SASAGeometry() ;
  int makeFromFiles(char *, char *, char *) ;
  ~SASAGeometry() ;
//globals
  std::list<Eigen::Vector3d> xyzPositions, abcPositions, ref_xyzPositions ;
  std::list<double> Radii ;
  Eigen::Matrix3d nonortho_to_ortho, ortho_to_nonortho ;
  Eigen::Vector3d basisLengths, basisAngles ;
private:
//methods
  int radiiPopulator() ;
  int xyzRead(char *) ;
  int datRead(char *) ;
  int atomsRead(char *) ;
//locals
  std::list<std::string> atomicTypes, elements ;
  std::list<double> atom_radii ;
};
#endif
