#include"stress2dprocess.h"
/*-----------------------------------*/
/*          �������̿���             */
/*-----------------------------------*/
Matrix<double, 3, 1> Stress_Calculate(const Matrix<double, 3, 1> &sumStress, const Matrix<double, 3, 1> &sumStrain, const Matrix<double, 3, 1> &dStrain, const double &Temperature, const double &dt)
{
	//���Ʊ������ȣ���С���ݹ���ü���dStress��
	Matrix<double, 3, 1> dStress;
	dStress.setZero();
	//����¶���Ӧ���Ƿ����仯��������仯�����ֺ���㲢����
	double threshold;
	if (Temperature > 473)  threshold = 1e-5; else threshold = 1e-5;
	int n = abs(dt) / 30, n0 = abs(dStrain(0, 0) / threshold), n1 = abs(dStrain(1, 0) / threshold), n2 = abs(dStrain(2, 0) / threshold);
	if (n0 > n) n = n0; if (n1 > n) n = n1; if (n2 > n) n = n2;
	if (n != 0)
	{ //����仯
		Matrix<double, 3, 1> strainstep = dStrain / (n + 1);
		double dt_new = dt / (n + 1);
		Matrix<double, 3, 1> sumStress_temp = sumStress, sumStrain_temp = sumStrain;
		double Temperature_temp = Temperature;
		for (int i = 0; i < n + 1; i++)
		{
			Matrix<double, 3, 1> dStress_temp = Stress_Calculate(sumStress_temp, sumStrain_temp, strainstep, Temperature_temp, dt_new);
			sumStress_temp += dStress_temp;
			sumStrain_temp += strainstep;
			Temperature_temp += dt_new;
			dStress += dStress_temp;
		}
		return dStress;
	}
	else
	{ //С���仯
		double e = Expansion_In_T(Temperature), enext = Expansion_In_T(Temperature + dt);
		Matrix<double, 3, 1> dexpansion(enext - e, enext - e, 0);
		//�Լ�������
		Matrix<double, 3, 1> de0 = StrainIncrease_Temperature(sumStress, Temperature, dt);
		Dep_and_ST Dep_ST = Dep_and_StressInrease_T(sumStress, sumStrain, Temperature, dt);
		dStress = Dep_ST.Dep*(dStrain - de0 - dexpansion) + Dep_ST.StressIncrease_T;
		//������ء�ж��״̬���ж�����
		double s1 = sumStress(0, 0), s2 = sumStress(1, 0), s3 = sumStress(2, 0),
			ds1 = dStress(0, 0), ds2 = dStress(1, 0), ds3 = dStress(2, 0);
		double CRI = (2 * s1*ds1 + 2 * s2*ds2 - s1*ds2 - s2*ds1 + 6 * s3*ds3) / 3;  //dJ2
		if (CRI >= 0)
		{//���ع���
			return dStress;
		}
		else
		{//ж�ع��̣�������������뵥��ж�أ�
			Matrix<double, 3, 3> De = De_In_T(Temperature);
			dStress = De*(dStrain - de0 - dexpansion); 
			double ds1_unload = dStress(0, 0), ds2_unload = dStress(1, 0), ds3_unload = dStress(2, 0);
			double CRI_unload = (2 * s1*ds1_unload + 2 * s2*ds2_unload - s1*ds2_unload - s2*ds1_unload + 6 * s3*ds3_unload) / 3;
			if (CRI_unload < 0) 
			{ 
				return dStress;
			}//����ж��
			else
			{ //������ع��̣��ٶ���ж����Ӧ����ȫΪ0״̬���Ҵ�ʱ�¶ȱ仯�������̵�һ��
				Matrix<double, 3, 1> zero, d_elastic_strain, mid_Strain, dStress1, dStress2;
				zero.setZero(); d_elastic_strain.setZero();
				dStress1 = zero - sumStress;
				d_elastic_strain = De.inverse()*dStress1;
				mid_Strain = sumStrain + d_elastic_strain + de0 / 2 + dexpansion / 2;
				Matrix<double, 3, 1> dStrain_now = sumStrain + dStrain - mid_Strain;
				dStress2 = Stress_Calculate(zero, mid_Strain, dStrain_now, Temperature + dt / 2, dt / 2);
				return dStress1 + dStress2;
			}
		}    
	}
}
std::shared_ptr<stress2d> Stress_Process(const std::vector<double>& exx, const std::vector<double>& eyy, const std::vector<double>& exy, const std::vector<double>& temperature,int length)
{
	/*����һ���Ӧ�����¶ȣ�������Ӧ���仯����*/
	using std::shared_ptr;
	using std::make_shared;
	 //���exx,eyy,exy,temperature���ݳ����Ƿ�һ��
	shared_ptr<stress2d> p = make_shared<stress2d>();
	if (!(exx.size()>=length&&eyy.size()>=length&&exy.size()>=length))
	{
		std::cout << "Data input length inconsistencies !"<<std::endl;
		p->tag = 0;
		return 0;
	}
	 //��ʼ����Ӧ���仯
	p->tag = 1;
	int ele_size = length;
	Matrix<double, 3, 1> sumStress,exStrain;
	sumStress.setZero();
	exStrain.setZero();
	double exTemperature = 300;
	for (int i = 0; i < ele_size; i++)
	{
		//�¶�Temperature��������dt
		double Temperature = temperature[i] + 273;
		Temperature = Temperature > 300 ? Temperature : 300;  //С��27��ʱ�¶�����Ϊ27�ȣ�300K
		double dt = Temperature - exTemperature;
		//Ӧ��sumStrain��������dStrain
		Matrix<double, 3, 1> sumStrain(exx[i], eyy[i], exy[i]);
		Matrix<double, 3, 1> dStrain=sumStrain-exStrain;
		//����Ӧ������
		Matrix<double, 3, 1> dStress = Stress_Calculate(sumStress, sumStrain, dStrain, Temperature, dt);
		//������Ӧ������Ӧ��,�¶�
		sumStress += dStress;
		exStrain = sumStrain;
		exTemperature = Temperature;
		//�����������������ֵ
		p->sxx.push_back(sumStress(0, 0));
		p->syy.push_back(sumStress(1, 0));
		p->sxy.push_back(sumStress(2, 0));
	}
	//�����굥��Ӧ���󣬵��ú������������Ϣ��static���������㣬��һ�㷽�ɼ��㣻
	Euqal_Strain_Ele_Property::Reset_all();
	return p;
}
