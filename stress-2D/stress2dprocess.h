#include"stressfunction.h"
#include<vector>
#include<memory>
/*-----------------------------------*/
/*          �������̿���             */
/*-----------------------------------*/
struct stress2d
{
	std::vector<double> sxx;
	std::vector<double> syy;
	std::vector<double> sxy;
	int tag; //tag=1������Ч��0Ϊ��Ч
	stress2d() :sxx(), syy(), sxy(), tag(0){}
};
Matrix<double, 3, 1> Stress_Calculate(const Matrix<double, 3, 1> &sumStress, const Matrix<double, 3, 1> &sumStrain, const Matrix<double, 3, 1> &dStrain, const double &Temperature, const double &dt);
std::shared_ptr<stress2d> Stress_Process(const std::vector<double>& exx, const std::vector<double>& eyy, const std::vector<double>& exy, const std::vector<double>& temperature,int length);
