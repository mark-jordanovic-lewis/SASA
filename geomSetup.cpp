//Provides access to the SASAGeometry object for testing.

#include <cstdlib>
#include <iostream>
#include "SASAGeometry.h"

using namespace Eigen ;
using namespace std ;

/*struct geometry*/ int main(int argv, char * argc[])
{
  list<Vector3d>::iterator listIterator ;
  char * xyz_file = argc[1] ;
  char * dat_file = argc[2] ;
  char * atoms_file = argc[3] ;
  SASAGeometry geomMaker ;
  int geomErr ;

  geomErr = geomMaker.makeFromFiles(xyz_file, dat_file, atoms_file) ;
	
  cout << "xyz's" << endl ;
  cout << "====" << endl ; cout << endl ;  
  for(listIterator = geomMaker.xyzPositions.begin() ; listIterator != geomMaker.xyzPositions.begin() ; listIterator++){
    cout << (*listIterator)[0] << " " <<  (*listIterator)[1] << " " <<  (*listIterator)[2] << endl ;
  }
  cout << endl ;

  cout << "abc's" << endl ;
  cout << "====" << endl ; cout << endl ;  
  for(listIterator = geomMaker.abcPositions.begin() ; listIterator != geomMaker.abcPositions.begin() ; listIterator++){
    cout <<(*listIterator)[0] << " " <<  (*listIterator)[1] << " " <<  (*listIterator)[2] << endl ;
  }
  cout << endl ;
	
  cout << "basis Angles" << endl ;
  cout << "==========" << endl ; cout << endl ;
  cout <<  geomMaker.basisAngles[0] << " " <<  geomMaker.basisAngles[1] << " " <<   geomMaker.basisAngles[2] << " " << endl ;
  cout << endl ;
	
  cout << "basis Lengths" << endl ;
  cout << "==========" << endl ; cout << endl ;
  cout <<  geomMaker.basisLengths[0] << " " <<  geomMaker.basisLengths[1] << " " <<   geomMaker.basisLengths[2] << " " << endl ;
  cout << endl ;

  return 0 ;
}
