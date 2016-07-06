//Mark lewis - sasa calculator extension for Avogadro
// this will handle the neccessary files required by the calculator
// .xyz, .dat and .atom 

#ifndef SASA_GEOMETRYSETUP_H
#define SASA_GEOMETRYSETUP_H

#include <cstdlib>
#include <vector>
#include <sstream>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/LU>
#include <cmath>
#include <cstring>
#include <string>
#include <list>

#include "sasa_transformMatrix.h"
#include "sasa_tokenUtil.h"

class sasa_geometrySetup
{
public:
	bool Read(std::string, std::string, std::string) ;
	sasa_geometrySetup() ;
	~sasa_geometrySetup() ;
	bool consolodate() ;
	
	//globals
	std::list<Eigen::Vector3d> xyzPositions, abcPositions, ref_xyzPositions ;
	std::list<double> Radii ;
	Eigen::Matrix3d nonortho_to_ortho, ortho_to_nonortho ;
	Eigen::Vector3d basisLengths, basisAngles ;
	
private:
	bool dat(const char* fileName) ;
	bool xyz(const char* fileName) ;
	bool atom(const char* fileName) ;
	bool produceGeometry() ;
	void radiiPopulator() ;

	//locals
	tokenUtil dataTokens ;
	std::list<std::string> atomicTypes, elements ;
	std::list<double> atom_radii ;
};
#endif
