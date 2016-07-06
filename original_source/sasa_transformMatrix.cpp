#include "sasa_transformMatrix.h"

using namespace std ;
using namespace Eigen ;

double volume ;

sasa_transformMatrix::sasa_transformMatrix()
{/* EMPTY CONSTRUCTOR*/}

Matrix3d sasa_transformMatrix::makeTransformationMatrix(Vector3d angles, Vector3d lengths)
{
	double alpha, beta, gamma, bi, bj, ci, cj ;
	Matrix3d returnMatrix ;
	double pi = 3.1415926535897932384626433832795 ;

	alpha = pi*angles(0)/180.0 ;
	beta = pi*angles(1)/180.0 ;
	gamma = pi*angles(2)/180.0 ;

	bi = lengths(1)*cos(gamma) ;
	bj = lengths(1)*sin(gamma) ;
	ci = lengths(2)*cos(beta) ;
	cj = (lengths(1)*lengths(2)*cos(alpha) - bi*ci)/bj ;

	returnMatrix(0,0) = lengths(0) ;
	returnMatrix(0,1) = bi ;
	returnMatrix(0,2) = ci ;
	returnMatrix(1,0) = 0 ;
	returnMatrix(1,1) = bj ;
	returnMatrix(1,2) = cj ;
	returnMatrix(2,0) = 0 ;
	returnMatrix(2,1) = 0 ;
	returnMatrix(2,2) = sqrt(lengths(2)*lengths(2)-ci*ci-cj*cj) ;

	volume = abs(lengths(0)*bj*returnMatrix(2,2)) ; //there are zero terms which are omitted here due to the
                                                    //matrix being upper triangular
	return returnMatrix ;
}

double sasa_transformMatrix::getVol()
{
	return volume ;
}
