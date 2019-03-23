#include"SUS304function.h"
/*-----------------------------------*/
/*       SUS304���ܼ���ֵ���溯��    */
/*-----------------------------------*/
const spline1dinterpolant& Interpolant_Yongmodulus()
{
	//��ֵ����ģ��-T������
	static bool Is_caculated = false;
	static spline1dinterpolant s;
	if (Is_caculated == false)  //ֻ��Ҫ��һ�μ��㱣�棬��߲�ѯ����
	{
		extern const std::vector<double> YoungModulus;
		extern const std::vector<double> TemperInPandY;
		auto i = YoungModulus.size(), j = TemperInPandY.size();
		double *p1 = new double[i];
		double *p2 = new double[j];
		memcpy(p1, &YoungModulus[0], i * sizeof(double));
		memcpy(p2, &TemperInPandY[0], j * sizeof(double));
		real_1d_array y, T;
		y.setcontent(i, p1);
		T.setcontent(j, p2);
		spline1dbuildcubic(T, y, s);
		Is_caculated = true;
		delete p1;
		delete p2;
	}
	return s;
}
double Youngmodulus_vs_T(const double &Temperature)
{
	//��������ģ�����¶ȵĵ���
	double r, ds, d2s; //r��T�¶��µ�����ģ����ds��һ�׵�����d2s�Ƕ��׵�
	if (Temperature < 300) return 0; //�ٶ�С��300Kʱ����ģ���ޱ仯
	else
	{
		const spline1dinterpolant s(Interpolant_Yongmodulus());
		spline1ddiff(s, Temperature, r, ds, d2s);
		return ds;
	}
}
double Youngmodulus_In_T(const double &Temperature)
{
	//��ֵ�õ�T�¶��µ�����ģ��
	const spline1dinterpolant s(Interpolant_Yongmodulus());
	if (Temperature <= 300) return spline1dcalc(s, 300); //�ٶ�С��300Kʱ����ģ��Ϊ300Kʱ��ֵ
	else return spline1dcalc(s, Temperature);
}
const spline1dinterpolant& Interpolant_PoissonRatio()
{
	//��ֵ�õ�T�¶��µĲ��ɱ�-T����
	static bool Is_caculated = false;
	static spline1dinterpolant s;
	if (Is_caculated == false)  //ֻ��Ҫ��һ�μ��㱣�棬��߲�ѯ����
	{
		extern const std::vector<double> PoissonRatio;
		extern const std::vector<double> TemperInPandY;
		auto i = PoissonRatio.size(), j = TemperInPandY.size();
		double *p1 = new double[i];
		double *p2 = new double[j];
		memcpy(p1, &PoissonRatio[0], i * sizeof(double));
		memcpy(p2, &TemperInPandY[0], j * sizeof(double));
		real_1d_array y, T;
		y.setcontent(i, p1);
		T.setcontent(j, p2);
		spline1dbuildcubic(T, y, s);
		Is_caculated = true;
		delete p1;
		delete p2;
	}
	return s;
}
double PoissonRatio_vs_T(const double &Temperature)
{
	//���㲴�ɱȶ��¶ȵĵ���
	double r, ds, d2s; //r��T�¶��µĲ��ɱȣ�ds��һ�׵�����d2s�Ƕ��׵�
	if (Temperature < 300) return 0;  //�ٶ�С��300Kʱ���ɱ��ޱ仯
	else
	{
		const spline1dinterpolant s(Interpolant_PoissonRatio());
		spline1ddiff(s, Temperature, r, ds, d2s);
		return ds;
	}
}
double PoissonRatio_In_T(const double &Temperature)
{
	//��ֵ�õ�T�¶��µĲ��ɱ�
	const spline1dinterpolant s(Interpolant_PoissonRatio());
	if (Temperature <= 300) return spline1dcalc(s, 300); //�ٶ�С��300Kʱ���ɱ�Ϊ300Kʱ��ֵ
	else return spline1dcalc(s, Temperature);
}
const spline1dinterpolant& Interpolant_Expantion()
{
	//��ֵ����������
	static bool Is_caculated = false;
	static spline1dinterpolant s;
	if (Is_caculated == false)  //ֻ��Ҫ��һ�μ��㱣�棬��߲�ѯ����
	{
		extern const std::vector<double> Expansion;
		extern const std::vector<double> TemperInExpan;
		auto i = Expansion.size(), j = TemperInExpan.size();
		double *p1 = new double[i];
		double *p2 = new double[j];
		memcpy(p1, &Expansion[0], i * sizeof(double));
		memcpy(p2, &TemperInExpan[0], j * sizeof(double));
		real_1d_array y, T;
		y.setcontent(i, p1);
		T.setcontent(j, p2);
		spline1dbuildcubic(T, y, s);
		Is_caculated = true;
	}
	return s;
}
double Expansion_In_T(const double &Temperature)
{
	//��ֵT�¶��µ�������
	const spline1dinterpolant s(Interpolant_Expantion());
	if (Temperature <= 300) return 0;  //��׼zero=300k
	else return spline1dcalc(s, Temperature)*(Temperature - 300);
}
double T300K_TensileStress_vs_Strain(const double &Strain)
{
	//��Ϊ��ֵ������Ӧ��������T�͵�Ч����Ӧ��ģ�������T<=300Kʱ�޵��������趨��T<=300Kʱ������Ӧ���Ե�Ч����Ӧ��ƫ������300Kʱ������Ӧ�������Ե�ЧӦ���ƫ��������Ӧ����T����Ϊ0
	//����Ч����Ӧ��Ϊ0ʱ���趨����Ӧ���Ե�Ч����Ӧ�䵼��Ϊ����
	//ע���˴��Ƕ�300Kʱ����Ӧ�������Ե�ЧӦ���󵼣�Strain�ǵ�Ч����Ӧ��
	static bool Is_caculated = false;
	static spline1dinterpolant s;
	//��ֵ�õ�300kʱ����Ӧ��-����Ӧ������
	if (Is_caculated == false)  //ֻ��Ҫ��һ�μ��㱣�棬��߲�ѯ����
	{
		extern const std::vector<double> STRESS300K;
		extern const std::vector<double> STRAIN;
		auto i = STRESS300K.size(), j = STRAIN.size();
		double *p1 = new double[i];
		double *p2 = new double[j];
		memcpy(p1, &STRESS300K[0], i * sizeof(double));
		memcpy(p2, &STRAIN[0], j * sizeof(double));
		real_1d_array y, S; //S������Ӧ��
		y.setcontent(i, p1);
		S.setcontent(j, p2);
		spline1dbuildcubic(S, y, s);  //splinecubic��ֵ
		Is_caculated = true;
		delete p1;
		delete p2;
	}
	double r, ds, d2s; //r��Sʱ������Ӧ��ֵ��ds��һ�׵�����d2s�Ƕ��׵�
	if (Strain == 0) return std::numeric_limits<double>::max(); //Strain�ǵ�Ч����Ӧ��
	else
	{
		spline1ddiff(s, Strain, r, ds, d2s);
		return ds;
	}
}
const spline2dinterpolant& Interpolant_TensileStress()
{
	//��ֵ �¶�������Ӧ��Ϊ������ʵ��Ӧ�����棬��Ӧ�����漴Ϊ����Ӧ������
	static bool Is_caculated = false;
	static spline2dinterpolant s;
	if (Is_caculated == false)  //ֻ��Ҫ��һ�μ��㱣�棬��߲�ѯ����
	{
		extern const std::vector<double> TEMPERATURE;
		extern const std::vector<double> STRAIN;
		extern const std::vector<double> STRESS;
		auto i = TEMPERATURE.size(), j = STRAIN.size(), k = STRESS.size();
		double *pi = new double[i];
		double *pj = new double[j];
		double *pk = new double[k];
		memcpy(pi, &TEMPERATURE[0], i * sizeof(double));
		memcpy(pj, &STRAIN[0], j * sizeof(double));
		memcpy(pk, &STRESS[0], k * sizeof(double));
		real_1d_array T, strain, stress;
		T.setcontent(i, pi);
		strain.setcontent(j, pj);
		stress.setcontent(k, pk);
		spline2dbuildbicubicv(strain, j, T, i, stress, 1, s);
		Is_caculated = true;
		delete pi;
		delete pj;
		delete pk;
	}
	return s;
}
double TensileStress(const double &T, const double &Strain)
{
	//����T�¶��µ�Ч����Ӧ��ΪStrainʱ������Ӧ��
	const spline2dinterpolant s(Interpolant_TensileStress());
	//��������Ե�ЧӦ���ܴ��ڵ���0,�Ұ��¶ȵ���300kʱ��״̬����Ϊ��300kʱ״̬
	if (T <= 300) return spline2dcalc(s, Strain, 300);
	else return spline2dcalc(s, Strain, T);
}
double TensileStress_vs_Strain(const double &T, const double &Strain)
{
	//����T�¶��µ�Ч����Ӧ��ΪStrainʱ��ʵ��Ӧ��������Ӧ������������Ӧ���ƫ����
	//��Ϊ��ֵ������Ӧ��������T<=300Kʱ�޵��������趨��T<=300Kʱ������Ӧ���Ե�Ч����Ӧ�����300Kʱ������Ӧ�������Ե�ЧӦ���ƫ��������Ӧ����T����Ϊ0��
	//����Ч����Ӧ��Ϊ0ʱ���趨����Ӧ���Ե�Ч����Ӧ�䵼��Ϊ���
	if (Strain == 0) { return std::numeric_limits<double>::max(); } //δ����ʱ������Ϊ����
	else
	{
		if (T <= 300) return T300K_TensileStress_vs_Strain(Strain);
		else {
			double Dtemper, Dstrain, TensileStress, Dtemper_and_strain;
			const spline2dinterpolant s(Interpolant_TensileStress());
			spline2ddiff(s, Strain, T, TensileStress, Dstrain, Dtemper, Dtemper_and_strain);
			return Dstrain;
		}
	}
}
double TensileStress_vs_Temper(const double &T, const double &Strain)
{
	//����T�¶��� /��Ч����Ӧ��/ ΪStrainʱ��  /ʵ��Ӧ��������Ӧ����/  �����¶ȵ�ƫ��
	//��Ϊ��ֵ������Ӧ����T�͵�Ч����Ӧ���������T<=300Kʱ�޵��������趨��T<=300Kʱ������Ӧ���Ե�Ч����Ӧ�����300Kʱ������Ӧ�������Ե�ЧӦ���ƫ��������Ӧ����T����Ϊ0
	//����Ч����Ӧ��Ϊ0ʱ���趨����Ӧ���Ե�Ч����Ӧ�䵼��Ϊ����
	if (Strain == 0) return 0;  //δ����ʱ������Ӧ�����¶ȵ���Ϊ0
	else
	{
		if (T < 300) return 0;
		else
		{
			double Dtemper, Dstrain, TensileStress, Dtemper_and_strain;
			const spline2dinterpolant s(Interpolant_TensileStress());
			spline2ddiff(s, Strain, T, TensileStress, Dstrain, Dtemper, Dtemper_and_strain);
			return Dtemper;
		}
	}
}
Matrix<double, 3, 3> De_In_T(const double &T)
{
	//����T�¶�ʱ�ĵ��Ծ���De��uΪ���ɱȣ�EΪ����ģ��
	double E, u;
	E = Youngmodulus_In_T(T);
	u = PoissonRatio_In_T(T);
	Matrix<double, 3, 3> De;
	De.setZero();
	double temp = E / (1 - u*u);
	De(0, 0) = 1; De(1, 1) = 1; De(2, 2) = (1 - u) / 2;
	De(0, 1) = u; De(1, 0) = u;
	return De*temp;
}