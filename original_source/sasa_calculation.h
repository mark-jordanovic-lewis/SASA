#ifndef SASA_CALCULATION_H
#define SASA_CALCULATION_H

#include <ctime>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/LU>
#include <cmath>
#include <cstring>
#include <string>
#include <list>

#include "sasa_RandomGenerator.h"
#include "sasa_tokenUtil.h"

class sasa_calculation
{
public:
	//method instantiators
	sasa_calculation() ;
	void setParams(int*, double*, std::list<Eigen::Vector3d>*, std::list<Eigen::Vector3d>*, std::list<double>*, Eigen::Matrix3d*, Eigen::Matrix3d*, Eigen::Vector3d*, Eigen::Vector3d*) ;
	double calculateSA() ;
	~sasa_calculation() ;
		
private:
	//method instantiators
	void reduceABClist(int) ;
	void probeSetup() ;
	void manipulateProbe() ;
	void probeTest(int) ;
	
	//local general calculation variables
	const static double pi = 3.141592653589793 ;
	int Natoms ;
	double iRadius ;
	int successfulTrialCount ;
	std::list<Eigen::Vector3d> availableSites ;

	//local reduced test atom lists
	std::list<Eigen::Vector3d> potentialOverlapAtoms ;
	std::list<double> potentialOverlapRadii ;

	//local suface o/p file
	fstream surfaceFile('sasa_surface.dat', ios::out|ios::trunc) ; 
	
	//local random probe vars
	Eigen::Vector3d probe, site, nonortho_probe, manipulated_probe, final_probe  ;
	double dx, dy, dz, probeSize, checkDist ;
	
	//local pointer vars, to be converted to geometry type asap
	int *nTrials ;
	double *probeRadius ;
	std::list<Eigen::Vector3d> *xyzPositions ;
	std::list<Eigen::Vector3d> *abcPositions ;
	std::list<double> *Radii ;
	Eigen::Matrix3d *ortho_to_nonortho, *nonortho_to_ortho ;
	Eigen::Vector3d *basisLengths, *basisAngles ;
	
	//local iterators to iterate through the public geometry obj's
	std::list<Eigen::Vector3d>::iterator ithPos ;
	std::list<Eigen::Vector3d>::iterator xyzPos ;
	std::list<double>::iterator iRad ;
	std::list<Eigen::Vector3d>::iterator kthPos ;
	std::list<double>::iterator kRad ;
} ;
#endif
