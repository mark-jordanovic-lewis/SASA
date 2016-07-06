#include "sasa_calculation.h"

using namespace std ;
using namespace Eigen ;

Random rnd(clock()) ;

sasa_calculation::sasa_calculation(){}

void sasa_calculation::setParams(int *nT, double *pR, list<Vector3d> *xyz, list<Vector3d> *abc, list<double> *R, Matrix3d *otn, Matrix3d *nto, Vector3d *bL, Vector3d *bA)
{																				//local variables for calculation
  nTrials = nT ;	
  probeRadius = pR ;
  xyzPositions = xyz ;
  abcPositions = abc ;
  Radii = R ;
  ortho_to_nonortho = otn ;
  nonortho_to_ortho = nto ;
  basisLengths = bL ;
  basisAngles = bA ;
}

double sasa_calculation::calculateSA()													//calls methods to provide calculation
{
  cout << "Crystal Setup, Surface Area Calculation commencing" << endl ;
  cout << "-----------------------------------------------------" << endl ;
	
  Natoms = (*xyzPositions).size() ;													//count limit for loops
  ithPos = (*abcPositions).begin() ;													//list iterators
  xyzPos = (*xyzPositions).begin() ;
  iRad = (*Radii).begin() ;
	
  double totAvailableSurface = 0 ;
	
  for(int i = 0 ; i < Natoms ; ++i, ithPos++, iRad++, xyzPos++)								//lists iterated in parallel, novel : abc, xyz & radius per i-loop
    {
      iRadius = *iRad ;															//store ith atom radius
      successfulTrialCount = 0 ;
		
      reduceABClist(i) ;															//call reduction method, i-loop iterator as arg
		
      for(int t = 0 ; t < (*nTrials) ; ++t)												//no geometry iterated in t-loop header
	{
	  probeSetup() ;															//call probeSetup
	  probeTest(i) ;															//call to probeTest, i-loop iterator as arg 
	}

      double SAfrac = (double)successfulTrialCount/(*nTrials) ;								//calculation of ith atom contribution to total surface
      double iSurfaceRadius = (iRadius + *probeRadius) ;
      double availableSurface = 4*pi*(iSurfaceRadius)*(iSurfaceRadius)*SAfrac ;
      totAvailableSurface += availableSurface ;
    }
  cout << "total Probe Accessable Surface: " << totAvailableSurface << endl ;
	
  return  totAvailableSurface ;
}
	
//called from i-loop
void sasa_calculation::reduceABClist(int i)												//reduces calculations by checking dist between ith_xyz and kth_xyz and maxProbeLength
{
  potentialOverlapAtoms.clear() ;														//each i-loop => new reduced list	
  potentialOverlapRadii.clear() ;
  kthPos = (*xyzPositions).begin() ;
  kRad = (*Radii).begin() ;
  list<Vector3d>::iterator pushPos = (*abcPositions).begin() ;
  //i-loop(i_xyz++, i_abc++, i_rad++){
  for(int k = 0 ; k < Natoms ; ++k, kthPos++,  pushPos++, kRad++)						//iterate lists in parallel, novel : xyz, abc, radius per k-loop				k-loop(k_xyz++, k_abc++, k_rad++)
    {
      if(k == i){++k ; kthPos++ ; pushPos++ ; kRad++ ;}
		
      Vector3d i2kVect = *xyzPos - *kthPos ;											// ithPos -> kthPos xyz vector (xyzPos accessed, not changed)
      double i2kDist = sqrt(i2kVect.dot(i2kVect)) ;
      double maxProbeLength = 2*(*iRad) + 2*(*probeRadius) + 2*(*kRad) ;							//look @ var name....
		
      if(maxProbeLength > i2kDist)  												//if maxProbeLength larger than dist between ith and kth atom
	{
	  potentialOverlapAtoms.push_back(*pushPos) ;									//store abc of kth atom
	  potentialOverlapRadii.push_back(*kRad) ;										//store rad of kth atom
	}
    }																			//return to i-loop
}



void sasa_calculation::probeSetup()														//random creation of vector pointing from ith abcPos
{
  double phi, costheta, theta ;
		
  phi = 2*pi*rnd.Uniform() ;                  												//random angles : theta from x -axis, phi from z-axis
  costheta = 1 - 2*rnd.Uniform() ;
  theta = acos(costheta) ;
	
  probe(0) = sin(theta)*cos(phi) ;													//3d euclidean unit vector origin 0,0,0
  probe(1) = sin(theta)*sin(phi) ;													//probe
  probe(2) = cos(theta) ;

  site = probe*iRadius + *xyzPos ;											//t'th-xyzPos of probe-atom surface contact
  /*site = (*ortho_to_nonortho)*(probe*iRadius) ;
    site = site.cwise()*(*basisLengths) ;
    site += *ithPos ;*/														//t'th abcPos of probe-atom surface contact
	
  probe *= ((*probeRadius) + iRadius) ; 													//ith Atom-specific probe radiusorigin 0, 0, 0
  //transforming the probe to the crystal coord system, origin 0,0,0 : probe -> nonortho_probe
  nonortho_probe = (*ortho_to_nonortho)*probe ; 
  nonortho_probe = nonortho_probe.cwise()*(*basisLengths) ;													//rescaling vector to crys basis lengths
  nonortho_probe+= *ithPos ; 															//ammending probe by ith abcPos, nonortho_probe coord position is probe atom center
	
  if(nonortho_probe(0) < 0.0) nonortho_probe(0) += (*basisLengths)(0) ;  						//checking probe center sits the unit cell, if not ammend
  if(nonortho_probe(0) >= (*basisLengths)(0)) nonortho_probe(0) -= (*basisLengths)(0) ;
  if(nonortho_probe(1) < 0.0) nonortho_probe(1) += (*basisLengths)(1) ;
  if(nonortho_probe(1) >= (*basisLengths)(1)) nonortho_probe(1) -= (*basisLengths)(1) ;
  if(nonortho_probe(2) < 0.0) nonortho_probe(2) += (*basisLengths)(2) ;
  if(nonortho_probe(2) >= (*basisLengths)(2)) nonortho_probe(2) -= (*basisLengths)(2) ;
	
	
}																				//reenter t-loop with fully constructed random vector, nonortho_probe


																				//called from t-loop
void sasa_calculation::probeTest(int i)													//more probe manipulation and then testing for k-atom overlap
{
  bool failedTrial = false ;
  kthPos = potentialOverlapAtoms.begin() ;
  kRad = potentialOverlapRadii.begin() ;
  int kcount = potentialOverlapAtoms.size() ;
	
  for(int k = 0 ; k < kcount ; ++k, kthPos++, kRad++)									//iterates list in parallel novel:  abc, radius per k-loop								
    {
      cout.precision(8) ;
      if(k == i){++k; kthPos++ ; kRad++ ;} 											//testing to see if the test and origin atom are the same

      manipulateProbe() ;															//nonortho_probe is manipulated again.		
		
      probeSize = sqrt(final_probe.dot(final_probe)) ;

      checkDist = *kRad + *probeRadius ;

      if(probeSize < 0.999*checkDist) 												//check for overlap of test atom and probe
	{
	  failedTrial = true ;
	  k = kcount ;
	}
    }
  if(!failedTrial)							//THIS NEEDS TO BE ACCOUNTED FOR WHEN PARALLELISING THE CODE
    {
      ++successfulTrialCount ;
      surfaceFile << site[0] << ' ' << site[1] << ' ' <<site[2] << endl ;		//THIS NEEDS COMPILING TO GET TEH OP TO RUN SASASnake
      availableSites.push_back(site) ;													//storing accesssible positions in a list for later use
    }
}
	


void sasa_calculation::manipulateProbe()   // === THIS NEEDS TO BE CHECKED BIG TIME, IF ITS WRONG THIS IS A HUGE ISSUE === \\
{
  manipulated_probe = nonortho_probe-(*kthPos) ;									//vector pointing between probe center and test atom center
  //nonortho_probe retains value for next call
  Vector3d fudge ; //dodgy FORTRAN int() mimicer --- this is the application of the minimum image convention treatment.

  fudge = 2*(manipulated_probe.cwise()/(*basisLengths)) ; //checking for the entire diameter vector of the probe-atom system, scaled relative to the basis lengths
	
  if(fudge(0) < 0) fudge(0) = ceil(fudge(0)) ;
  else fudge(0) = floor(fudge(0)) ;
  if(fudge(1) < 0) fudge(1) = ceil(fudge(1)) ;
  else fudge(1) = floor(fudge(1)) ;
  if(fudge(2) < 0) fudge(2) = ceil(fudge(2)) ;
  else fudge(2) = floor(fudge(2)) ;   			  //

  manipulated_probe -= (*basisLengths).cwise()*fudge ;
									
  manipulated_probe = manipulated_probe.cwise()/(*basisLengths) ;	
  //manipulated_probe -> final_probe
  final_probe = (*nonortho_to_ortho)*manipulated_probe ;
}
	

sasa_calculation::~sasa_calculation(){}
