// Mark Lewis 20-10-11 

#include <list>
#include <Eigen/Core>


pi = 3.141569 ;
struct geometry{ 
  int * nTrials ;	    
  double * probeRadius ;    
  list<Vector3d> * xyz ;    
  list<Vector3d> * abc ;    
  list<double> * Radii ;    
  Matrix3d * ortho_to_none ;
  Matrix3d * non_to_ortho ; 
  Vector3d * basisLengths ; 
  Vector3d * basisAngles ;};

class structSetup
{
 private:
  //int geometrySetup() ;
 public:
  structSetup() ;
  struct geometry geometrySetup() ;
};

#include <Random.h>

structSetup::structSetup(){} ;

int structSetup::geometrySetup(struct geometry * system)
{
  Random rndm ;
  int nTrials = 5 ;	    
  double probeRadius = 1.2 ;    
  list<Vector3d> xyzPositions ;    
  list<Vector3d> abcPositions ;    
  list<double> Radii ;    
  Matrix3d ortho_to_non ;
  Matrix3d non_to_ortho ; 
  Vector3d basisLengths ; 
  Vector3d basisAngles ;
  list<Vector3d>::iterator xyz, abc ;
  list<double>::iterator R ;
  
  xyz = xyzPositions.begin() ;
  abc = abcPositions.begin() ;
  R = Radii.begin() ;
  
  for(int i = 0 ; i < 3 ; ++i, xyz++, abc++, R++){
    basisLength[i] = 1 ;
    basisAngle[i] = pi/2 ;
    for(int j = 0 ; j <3 ; ++j){
      *xyz[j] = rndm.Normal() ;    
      *R = rndm.Normal() + 1 ;    
      if(i = j) ortho_to_non[i][j] = 1 ; non_to_ortho[i][j] = 1 ;
      else ortho_to_non[i][j] = 0 ; non_to_ortho[i][j] = 0 ;
    }
    *abc = ortho_to_non*(*xyz) ;
  }
}
