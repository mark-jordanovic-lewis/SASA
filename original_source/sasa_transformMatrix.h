//Mark Lewis - sasa_transforomation class
//produces a 3x3 Eigen Matrix object from crystal data unit cell parameters

#ifndef SASA_TRANSFORMATIONMATRIX_H
#define SASA_TRANSFORMATIONMATRIX_H

#include <cmath>
#include <iostream>
#include <Eigen/Core>

class sasa_transformMatrix
{
private:

public:
	sasa_transformMatrix() ;
	Eigen::Matrix3d makeTransformationMatrix(Eigen::Vector3d, Eigen::Vector3d) ;
	double getVol() ;
} ;
#endif
