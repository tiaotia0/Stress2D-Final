#include"SUS304function.h"
/*-----------------------------------*/
/*          Ӧ����ؼ��㺯��         */
/*-----------------------------------*/
struct Dep_and_ST
{
	//���ڼ���Dep���¶ȱ仯�����Ӧ������StressInrease_T(ST)ʽ�Ӽ���һ�����ʿ�ͬʱ���㣬��װ�ڴ˽ṹ����
	Matrix<double, 3, 3> Dep;
	Matrix<double, 3, 1> StressIncrease_T;
};
class  Euqal_Strain_Ele_Property
{
	//���ڴ��һ��ĵ�ЧӦ�䣬ǰӦ����ǰӦ�䣬ǰ�¶�
private:
	static double eq_plasticstrain, preTemperature;
	static Matrix<double, 3, 1> preStress, preStrain;
public:
	Euqal_Strain_Ele_Property() {}
	static void Change_eq_plasticstrain(const double& c);
	static void Change_preTemperatrue(const double& c);
	static void Change_preStress(const Matrix<double, 3, 1>& c);
	static void Change_preStrain(const Matrix<double, 3, 1>& c);
	static double Query_eq_plasticstrain();
	static double Query_preTemperatrue();
	static Matrix<double, 3, 1> Query_preStress();
	static Matrix<double, 3, 1> Query_preStrain();
	static void Reset_all();
};
Matrix<double, 3, 1> Stress_Derivative(const Matrix<double, 3, 1> &stress);
double Equal_PlasticStrain(const Matrix<double, 3, 1> &Stress, const Matrix<double, 3, 1> &Strain, const double &Temperature);
double Equal_Stress(const Matrix<double, 3, 1> &Stress);
Dep_and_ST Dep_and_StressInrease_T(const Matrix<double, 3, 1> &Stress, const Matrix<double, 3, 1> &Strain, const double &Temperature, const double &dt);
Matrix<double, 3, 1> StrainIncrease_Temperature(const Matrix<double, 3, 1> &stress, const double &Temperature, const double &dt);