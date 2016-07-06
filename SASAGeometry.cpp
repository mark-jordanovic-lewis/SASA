//Mark Lewis 11-10-11

#include "SASAGeometry.h"

using namespace std ;
using namespace Eigen ;

SASAGeometry::SASAGeometry(){}

// === Workhorse methods === \\

/*
  calls file readers
  builds transform matrices from file data
  transforms xyz atomic euclidean coords into abc atomic crystal coords
  calls radiiPopulator method
  checks for system coherence
  return 0 if success, return 1 if failure
*/
int SASAGeometry::makeFromFiles(char * xyz_file, char * dat_file, char * atoms_file)
{
  sasa_transformMatrix basisMaker ;
  list<Vector3d>::iterator pos ;
  
  xyzRead(xyz_file) ;
  datRead(dat_file) ;
  atomsRead(atoms_file) ;
  
  nonortho_to_ortho = basisMaker.makeTransformationMatrix(basisAngles, basisLengths) ;			
  ortho_to_nonortho = nonortho_to_ortho.inverse() ;
  
  //omp ||-ise and cuda-ise matrix ops. in the push_back
  for(pos = ref_xyzPositions.begin() ; pos != ref_xyzPositions.end() ; pos++) abcPositions.push_back(ortho_to_nonortho*(*pos)) ;
  
  //omp & cuda ||-ise
  for(pos = abcPositions.begin() ; pos != abcPositions.end() ; pos++) {
    if((*pos)(0) < 0.0) (*pos)(0) += basisLengths(0) ;
    if((*pos)(0) >= basisLengths(0)) (*pos)(0) -= basisLengths(0) ;
    if((*pos)(1) < 0.0) (*pos)(1) += basisLengths(1) ;
    if((*pos)(1) >= basisLengths(1)) (*pos)(1) -= basisLengths(1) ;
    if((*pos)(2) < 0.0) (*pos)(2) += basisLengths(2) ;
    if((*pos)(2) >= basisLengths(2)) (*pos)(2) -= basisLengths(2) ;
  }
  
  radiiPopulator() ;
  
}

/*
  compares element names with element symbols and matches up with their radii.
  return 0 if sucess, return 1 and quits attempting to find any more matches if a single element is unmatched.
*/
int SASAGeometry::radiiPopulator()
{
  if(Radii.size() > 0) Radii.clear() ;
  int atom_radiusCount = 0, max_radiusCount = atom_radii.size() ;
  bool typeProcessed = false ;
  list<double>::iterator rads ;
  list<string>::iterator types, e, e_ ;

  //omp ||-ise 
  for(types = atomicTypes.begin() ; types != atomicTypes.end() ; types++, atom_radiusCount++){
    typeProcessed = false ; 
    atom_radiusCount = 0 ;
    e = elements.begin() ; 
    rads = atom_radii.begin() ;
    
    while(!typeProcessed){
      if(!strcmp((*types).c_str(), (*e).c_str())){
	Radii.push_back(*rads) ;
	typeProcessed = true ;
      }
      if(atom_radiusCount >= max_radiusCount){
	cerr << "no match of elemental type: " << *types << " in list: " << endl ;
	for(e_ = elements.begin() ; e_ != elements.end() ; e_++){
	  if(e_ == elements.begin()) cerr << *e_ << endl ;
	  else cerr << "                                        " << *e_ << endl ;
	}
	cerr << "check .atoms file or geometrySetup::radiiPopulator()" << endl ;
	return 1 ;
      }
      e++ ; rads ++ ;
    }
  }
  return 0 ;
}


// === File Readers === \\

/*
  reads positions of atoms in ensemble from xyz coord file
  return 0 if success, return 1 if file wrong format, return -1 if file not found
*/
int SASAGeometry::xyzRead(char * xyz_file)
{
  if(ref_xyzPositions.size() > 0) ref_xyzPositions.clear() ;
  if(atomicTypes.size() > 0) atomicTypes.clear() ;

  FILE * xyzFile ;
  char currentTypeBuff[5] ;
  Vector3d currentXYZ ;
  char * temp ;
  int readErr ;
                                          
  xyzFile = fopen(xyz_file, "r") ;
  fgets(temp, 0, xyzFile) ;
  fgets(temp, 0, xyzFile) ;
  fgets(temp, 0, xyzFile) ;

  if(xyz_file == NULL){
    perror("xyz_file open error") ;
    return -1 ;
  }else{//not sure about this bit, think it may be a bit buggy
    do{
      readErr = fscanf(xyzFile, "%s%f%f%f", currentTypeBuff, &currentXYZ[0], &currentXYZ[1], &currentXYZ[2]) ;
      if(readErr < 4){
	cout << "xyz file read error! Read " << readErr << " tokens: " << currentTypeBuff << ", instead of 4." << endl ;
	return 1 ;
      }
      atomicTypes.push_back(currentTypeBuff) ;
      ref_xyzPositions.push_back(currentXYZ) ;
    }while(!feof(xyzFile)) ;
  }
  ref_xyzPositions.pop_back() ;
  atomicTypes.pop_back() ;
  fclose(xyzFile) ;
  cout << xyz_file << " processed" << endl ;
  return 0 ;
}

/*
  reads crystal data from dat file
  return 0 if sucessful, return -1 if lengths or angles not correctly read, return 1 if file opening error
*/
int SASAGeometry::datRead(char * dat_file)
{
  FILE * datFile ;			
  string currentLine ;
  int readErr ;
  char * temp ;
  
  datFile = fopen(dat_file, "r") ;
  
  if(datFile == NULL){
    cout << "dat file read error" << endl ;
    return 1 ;
  }else{
    fgets(temp, 0, datFile) ;
    fgets(temp, 0, datFile) ;
    fgets(temp, 0, datFile) ;
    fgets(temp, 0, datFile) ;
    
    if(fscanf(datFile, "%f%f%f", &basisLengths[0], &basisLengths[1], &basisLengths[2]) < 3){
      cout << "dat file format error! Lengths." << endl ;
      return -1 ;
    }
    if(fscanf(datFile, "%f%f%f", &basisLengths[0], &basisLengths[1], &basisLengths[2]) < 3){
      cout << "dat file format error! Angles." << endl ;
      return -1 ;
    }
    fclose(datFile) ;
    cout << dat_file << " processed" << endl ;
    return 0 ;
  }
}
/*
  reads atomic radius data from atoms file
  return 0 if sucessful, return -1 file format error, return 1 if file read error
*/
int SASAGeometry::atomsRead(char * atoms_file)
{
  if(elements.size() > 0 ) elements.clear() ;
  if(atom_radii.size() > 0) atom_radii.clear() ;
  
  FILE * atomsFile ;			
  string currentType ;
  int readErr ;
  float currentRadius ;
  
  atomsFile = fopen(atoms_file, "r") ;
  if(atomsFile == NULL){
    cout << "atoms file read error!" << endl ;
    return 1 ;
  }else{
    do{
      if(fscanf(atomsFile, "%s%f", currentType.c_str(), &currentRadius) < 2){
	cout << "Atoms file format error!" << endl ;
	return 1 ;
      }else{
	elements.push_back(currentType) ;
	atom_radii.push_back(currentRadius/2) ;
      }
    }while(!feof(atomsFile)) ;
    elements.pop_back() ;
    atom_radii.pop_back() ;
  }
  
  list<string>::iterator e ;
  
  for(e = elements.begin() ; e != elements.end() ; e++){
    if(!strcmp((*e).c_str(), "Hydrogen")) (*e) = "H" ;
    else if(!strcmp((*e).c_str(), "Carbon")) (*e) = "C" ;
    else if(!strcmp((*e).c_str(), "Oxygen")) (*e) = "O" ;
    else if(!strcmp((*e).c_str(), "Nitrogen")) (*e) = "N" ;
    
    else if(!strcmp((*e).c_str(), "Cadmium")) (*e) = "Cd" ;
    else if(!strcmp((*e).c_str(), "Chromium")) (*e) = "Cr" ;
    else if(!strcmp((*e).c_str(), "Zinc")) (*e) = "Zn" ;
    /*
      else if(!strcmp((*e).c_str(), "Helium")) (*e) = "He" ;
      else if(!strcmp((*e).c_str(), "Lithium")) (*e) = "Li" ;
      else if(!strcmp((*e).c_str(), "Beryllium")) (*e) = "Be" ;
      else if(!strcmp((*e).c_str(), "Boron")) (*e) = "B" ;
      else if(!strcmp((*e).c_str(), "Flourine")) (*e) = "F" ;
      else if(!strcmp((*e).c_str(), "Neon")) (*e) = "Ne" ;
      else if(!strcmp((*e).c_str(), "Sodium")) (*e) = "Na" ;
      else if(!strcmp((*e).c_str(), "Magnesium")) (*e) = "Mg" ;
      else if(!strcmp((*e).c_str(), "Aluminium")) (*e) = "Al" ;
      else if(!strcmp((*e).c_str(), "Silicon")) (*e) = "Si" ;
      else if(!strcmp((*e).c_str(), "Phosphorous")) (*e) = "P" ;
      else if(!strcmp((*e).c_str(), "Sulphur") | !strcmp((*e).c_str(), "Sulfur")) (*e) = "S" ;
      else if(!strcmp((*e).c_str(), "Chlorine")) (*e) = "Cl" ;
      else if(!strcmp((*e).c_str(), "Argon")) (*e) = "Ar" ;*/
  }
  
  fclose(atomsFile) ;
  cout << atoms_file << " processed" << endl ;
  return 0 ;
}

// === Destructor === \\ -- not neccessary as yet

SASAGeometry::~SASAGeometry(){}

