#include "sasa_geometrySetup.h"

using namespace std ;
using namespace Eigen ;

sasa_geometrySetup::sasa_geometrySetup(){}

// === FILE READERS === \\
	
bool sasa_geometrySetup::Read(string xyz_file, string dat_file, string atoms_file)							//setup geometry of system by reading & processing data files
{						
	if(abcPositions.size() > 0) abcPositions.clear() ; 
	
	sasa_transformMatrix transformMaker ;
	list<Vector3d>::iterator pos ;
	
	const char* xyz_FILE = xyz_file.c_str() ; 														//store data filenames
	const char* dat_FILE = dat_file.c_str() ;
	const char* atoms_FILE = atoms_file.c_str() ;
	
	xyz(xyz_FILE) ;																			//process data file contents
	dat(dat_FILE) ;
	atom(atoms_FILE) ;

	nonortho_to_ortho = transformMaker.makeTransformationMatrix(basisAngles, basisLengths) ;			//make transform matrices
	ortho_to_nonortho = nonortho_to_ortho.inverse() ;

	for(pos = ref_xyzPositions.begin() ; pos != ref_xyzPositions.end() ; pos++) abcPositions.push_back(ortho_to_nonortho*(*pos)) ; 	//transform : x y z positions => a b c positions 
	
	for(pos = abcPositions.begin() ; pos != abcPositions.end() ; pos++)
	{
		if((*pos)(0) < 0.0) (*pos)(0) += basisLengths(0) ;  											//checking kth Atom center sits in unit cell, if not ammend. Would be nice if we could do this to the whole boundary
		if((*pos)(0) >= basisLengths(0)) (*pos)(0) -= basisLengths(0) ;
		if((*pos)(1) < 0.0) (*pos)(1) += basisLengths(1) ;
		if((*pos)(1) >= basisLengths(1)) (*pos)(1) -= basisLengths(1) ;
		if((*pos)(2) < 0.0) (*pos)(2) += basisLengths(2) ;
		if((*pos)(2) >= basisLengths(2)) (*pos)(2) -= basisLengths(2) ;
	}
	if(xyzPositions.size() > 0) xyzPositions.clear() ;
	Vector3d POS ;
	for(pos = abcPositions.begin() ; pos != abcPositions.end() ; pos++)
	{
		POS(0) = (*pos)(0)/basisLengths(0) ;											
		POS(0) = (*pos)(1)/basisLengths(1) ;
		POS(0) = (*pos)(2)/basisLengths(2) ;
		xyzPositions.push_back(nonortho_to_ortho*POS) ;
	}
	
	
	radiiPopulator() ;
																							//check for list correspondance
	if(Radii.size() < xyzPositions.size()) cerr << "Radii < xyzPositions : mismatch " << Radii.size() << " < " << xyzPositions.size() << " check geometrysetup::xyz() or ::radiiPopulator()" << endl ;
	else if(Radii.size() > xyzPositions.size()) cerr << "Radii > xyzPositions : mismatch " << Radii.size() << " > " << xyzPositions.size() << " check geometrysetup::xyz() or ::radiiPopulator()" << endl ;
	
	cout << "Read achieved" << endl ;	
	return true ;
}

void sasa_geometrySetup::radiiPopulator()													
{				
	if(Radii.size() > 0) Radii.clear() ;
	int atom_radiusCount = 0 ;
	bool typeProcessed = false ;
	list<double>::iterator rads ;
	list<string>::iterator types, e, e_ ;
	
	for(types = atomicTypes.begin() ; types != atomicTypes.end() ; types++, atom_radiusCount++) 							//compare atomicTypes & elements and store corresponding radii in list
	{
		typeProcessed = false ; atom_radiusCount = 0 ;
		e = elements.begin() ; rads = atom_radii.begin() ;
		while(!typeProcessed) 
		{
			if(!strcmp((*types).c_str(), (*e).c_str()))
			{
				Radii.push_back(*rads) ;
				typeProcessed = true ;
			}
			if(atom_radiusCount >= atom_radii.size())
			{
				cerr << "no match of elemental type: " << *types << " in list: " << endl ;
				for(e_ = elements.begin() ; e_ != elements.end() ; e_++)
				{
					if(e_ == elements.begin()) cerr << *e_ << endl ;
					else cerr << "                                        " << *e_ << endl ;
				}
				cerr << "check .atoms file or geometrySetup::radiiPopulator()" << endl ;
				typeProcessed = true ;
			}
			e++ ; rads ++ ;
		}
	}
}

bool sasa_geometrySetup::xyz(const char* fileName)												//read .xyz and store data in atomicTypes and xyzPositions lists
{	
	if(ref_xyzPositions.size() > 0) ref_xyzPositions.clear() ;
	if(atomicTypes.size() > 0) atomicTypes.clear() ;
	
	string currentLine ;
	Vector3d atomicPosition ;
	vector<string> ref_xyzPositionstr ;																											//LINE CONTENTS
	fstream in_file(fileName) ;	
	
	getline(in_file, currentLine, ' ') ; 																										//number of atoms in the molecule
	getline(in_file, currentLine, ' ') ;
	getline(in_file, currentLine, '\n') ; 																										//first two lines omitted 
	
	do{ 																					
		getline(in_file, currentLine, ' ') ; 																									//element symbol
		if(strcmp(currentLine.c_str(), " ") == -1) getline(in_file, currentLine, ' ' ) ;
		atomicTypes.push_back(currentLine.c_str()) ; 						
		
		getline(in_file, currentLine) ; 																										//xyz coords
	  	dataTokens.tokenize(ref_xyzPositionstr, currentLine.c_str(), " ") ; 							 
	  
		if(ref_xyzPositionstr.size() == 0) 																//if no coords, push 0 0 0 onto string list
		{
			ref_xyzPositionstr.push_back("0") ;
			ref_xyzPositionstr.push_back("0") ;
			ref_xyzPositionstr.push_back("0") ;
		}
		
		from_string<double>(atomicPosition(0), ref_xyzPositionstr[0], dec) ; 								//string list => Vector3d
		from_string<double>(atomicPosition(1), ref_xyzPositionstr[1], dec) ;
		from_string<double>(atomicPosition(2), ref_xyzPositionstr[2], dec) ;
		
		ref_xyzPositions.push_back(atomicPosition) ;													
	}while(!in_file.eof()) ;
	ref_xyzPositions.pop_back() ; 																	//remove final list entries: EOF and 0 0 0
	atomicTypes.pop_back() ;
	
	cout << fileName << " processed" << endl ;
	return true ;
}



bool sasa_geometrySetup::dat(const char* fileName) 													//reads .dat file for basis data 
{	
	string currentLine ;
	vector<string> baseStr ;
	ifstream file_in(fileName) ;																											//First four lines are dealt with in GUI :
	
	getline(file_in, currentLine, '\n') ; 																										//force field filename
	getline(file_in, currentLine, '\n') ; 																										//.xyz filename
	getline(file_in, currentLine, '\n') ; 																										//probe size
	getline(file_in, currentLine, '\n') ; 																										//number of trials
	
	getline(file_in, currentLine, '\n') ; 																//Lines of crystal data :						//a b c basis lengths
	dataTokens.tokenize(baseStr, currentLine.c_str(), " ") ; 										
	from_string<double>(basisLengths(0), baseStr[0], dec) ;									
	from_string<double>(basisLengths(1), baseStr[1], dec) ;
	from_string<double>(basisLengths(2), baseStr[2], dec) ;
	
	getline(file_in, currentLine, '\n') ;																										//alpha beta gamma basis angles
	dataTokens.tokenize(baseStr, currentLine.c_str(), " ") ;										
	from_string<double>(basisAngles(0), baseStr[0], dec) ;										
	from_string<double>(basisAngles(1), baseStr[1], dec) ;
	from_string<double>(basisAngles(2), baseStr[2], dec) ;
	
	getline(file_in, currentLine, '\n') ; 																										//crystal density
	cout << fileName << " processed" << endl ;
	return true ;
}



bool sasa_geometrySetup::atom(const char* fileName) 												//reads in the .atom file, put element symbols and VdW atom_radii in lists
{
	if(elements.size() > 0 ) elements.clear() ;
	if(atom_radii.size() > 0) atom_radii.clear() ;
	
	ifstream in_file(fileName) ;
	string elementStr, atom_radiusStr, currentLine ;
	int start, end ;
															
	while(strcmp(elementStr.c_str(), "EOF")) 														//each while iteration processes one line		//ELEMENT_NAME ELEMENT_DIAMETER
	{
		getline(in_file, currentLine, '\n') ;
		start = 0 ; 																	
		end = currentLine.find_first_of(' ', start) ;
		elementStr = currentLine.substr(start, end) ;
		elements.push_back(elementStr) ;													

		start = currentLine.find_first_not_of(' ',  end) ;											
		end = currentLine.length() -1 ;
		atom_radiusStr = currentLine.substr(start, end)  ;

		atom_radii.push_back((double)atof(atom_radiusStr.c_str())/2) ; 								
	}

	elements.pop_back() ;																		//final line removal
	atom_radii.pop_back() ;																

	list<string>::iterator e ;
	
	for(e = elements.begin() ; e != elements.end() ; e++) 											//changing element names to symbols 
	{
		if(!strcmp((*e).c_str(), "Hydrogen")) (*e) = "H" ;												//common expected elements
		else if(!strcmp((*e).c_str(), "Carbon")) (*e) = "C" ;
		else if(!strcmp((*e).c_str(), "Oxygen")) (*e) = "O" ;
		else if(!strcmp((*e).c_str(), "Nitrogen")) (*e) = "N" ;
		
		else if(!strcmp((*e).c_str(), "Cadmium")) (*e) = "Cd" ;										//Elements required by test cases
		else if(!strcmp((*e).c_str(), "Chromium")) (*e) = "Cr" ;
		else if(!strcmp((*e).c_str(), "Zinc")) (*e) = "Zn" ;
		/*
		else if(!strcmp((*e).c_str(), "Helium")) (*e) = "He" ;											//the rest of the first three rows of the periodic table.
		else if(!strcmp((*e).c_str(), "Lithium")) (*e) = "Li" ;											//commented out due to unneccessary processing during testing
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
	
	cout << fileName << " processed" << endl ;
	return true ;
}

// === DESTRUCTOR === \\

sasa_geometrySetup::~sasa_geometrySetup(){}
