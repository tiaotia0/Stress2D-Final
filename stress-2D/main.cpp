#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include"stress2dprocess.h"
#include<time.h>

//std::vector<std::vector<double> > query_temperature_vector(int points,int poin_steps)
//{
//	//��ȡ�¶ȳ��ļ��е����ݴ����һ����άvector�У��¶ȳ��ļ�һ������¶�Ϊһ�У�����˳��Ϊ�������ң����ϵ���
//	// pointsΪ������point_stepsΪÿ������¶���������
//	std::vector<std::vector<double> > v_temperature(points);
//	std::ifstream temperature_filein("C:\\Users\\Songyulin\\Desktop\\temperature1.txt");
//	if (temperature_filein.fail()) std::cout << "Temperature load fail!" << std::endl;
//	double tmp;
//	for (int i = 0; i < points; i++)
//	{
//		for (int j = 0; j<poin_steps; j++)
//		{
//			temperature_filein >> tmp;
//			v_temperature[i].push_back(tmp);
//		}
//	}
//	temperature_filein.close();
//	return v_temperature;
//}
//std::vector<int> query_temperature_vector_tag()
//{
//	//�����¶ȳ��ĵ�����Ӧ�䳡�ĵĵ�����һ�����¶ȳ�����<Ӧ�䳡�ĵ㣩���轫�¶ȳ�����Ӧ�䳡�����Ӧ��
//	//��Լ����Ӧ�䳡�������˳�򣨷Ǳ��������ƣ����������ṩһ�����˳���Ӧ���¶ȳ�vector���±�����
//	std::vector<int> v1, v2;
//	for (int i = 0; i < 225; i++)
//	{
//		v1.push_back(i);
//		v1.push_back(i);
//	}
//	auto p = v1.begin();
//	while (p != v1.end())
//	{
//		v2.insert(v2.end(), p, p + 30);
//		v2.insert(v2.end(), p, p + 30);
//		p += 30;
//	}
//	return v2;
//}
//std::string strain_file_name(int row, int col)
//{
//	//����ָ������������Ӧ�����ļ�����
//	char des[256] = "C:\\Users\\Songyulin\\Desktop\\dst\\";
//	char txt[] = ".txt";
//	std::stringstream srow, scol;
//	srow.fill('0');
//	scol.fill('0');
//	srow.width(2);
//	scol.width(2);
//	srow << row;
//	scol << col;
//	char crow[10], ccol[10];
//	srow >> crow;
//	scol >> ccol;
//	strcat_s(des, crow);
//	strcat_s(des, ccol);
//	strcat_s(des, txt);
//	return des;
//}
//void stress_file_out(const stress2d& result,int col)
//{
//	//��stress�����xx,yy,xy�����ļ������ÿ�����col������
//	std::ofstream outxx("C:\\Users\\Songyulin\\Desktop\\sxx.txt"),
//		outyy("C:\\Users\\Songyulin\\Desktop\\syy.txt"),
//		outxy("C:\\Users\\Songyulin\\Desktop\\sxy.txt");
//	int size = result.sxx.size(), cnt = 0;
//	for (int i = 0; i < size; i++)
//	{
//		outxx << result.sxx[i];
//		cnt++;
//		if (cnt == col) {
//			cnt = 0;
//			outxx << '\n';
//		}
//		else outxx << '\t';
//	}
//	cnt = 0;
//	for (int i = 0; i < size; i++)
//	{
//		outyy << result.syy[i];
//		cnt++;
//		if (cnt == col) {
//			cnt = 0;
//			outyy << '\n';
//		}
//		else outyy << '\t';
//	}
//	cnt = 0;
//	for (int i = 0; i < size; i++)
//	{
//		outxy << result.sxy[i];
//		cnt++;
//		if (cnt == col) {
//			cnt = 0;
//			outxy << '\n';
//		}
//		else outxy << '\t';
//	}
//	outxx.close();
//	outyy.close();
//	outxy.close();
//}
int main()
{
	/*std::vector<std::vector<double> > v_temperature = query_temperature_vector(225, 70);
	std::vector<int> v_tag = query_temperature_vector_tag();
	int tag = 0;
	stress2d result;
	for (int row = 0; row <29; row++)
	{
		
		for (int col = 29; col >=0; col--)
		{
			std::vector<double> exx, eyy, exy;
			double exx_tmp, eyy_tmp, exy_tmp;
			std::string name = strain_file_name(row + 1, col + 1);
			std::ifstream strain_filein(name);
			while (strain_filein >> exx_tmp >> eyy_tmp >> exy_tmp )
			{
				exx.push_back(exx_tmp);
				eyy.push_back(eyy_tmp);
				exy.push_back(exy_tmp);
			}
			auto stress_evolution = Stress_Process(exx, eyy, exy, v_temperature[v_tag[tag++]],70);
			if (stress_evolution->tag == 0)
			{
				std::cout << "caculation error!" << std::endl;
				return 0;
			}
			result.sxx.push_back(stress_evolution->sxx.back());
			result.syy.push_back(stress_evolution->syy.back());
			result.sxy.push_back(stress_evolution->sxy.back());
			strain_filein.close();
			std::cout << "POINT " << tag << " COMPLETE!" << std::endl;
		}
	}
	stress_file_out(result,30);*/
	std::vector<double> exx, eyy, exy,T;
	double exx_tmp, eyy_tmp, exy_tmp,T_tmp;
	std::ifstream strain_filein("C:\\Users\\Songyulin\\Desktop\\14.txt");
	auto a = strain_filein.good();
	while (strain_filein >> exx_tmp >> eyy_tmp >> exy_tmp >> T_tmp)
	{
		exx.push_back(exx_tmp);
		eyy.push_back(eyy_tmp);
		exy.push_back(exy_tmp);
		T.push_back(T_tmp);
	}
	clock_t t1 = clock();
	auto stress_evolution = Stress_Process(exx, eyy, exy, T, T.size());
	clock_t t2 = clock();
	std::cout << double(t2 - t1) / CLOCKS_PER_SEC;
	if (stress_evolution->tag == 0)
	{
		std::cout << "caculation error!" << std::endl;
		return 0;
	}
	std::ofstream strain_fileout("C:\\Users\\Songyulin\\Desktop\\14out50.txt");
	int cnt = 0;
	while (cnt<stress_evolution->sxx.size())
	{
		strain_fileout << stress_evolution->sxx[cnt] << '\n';
		cnt++;
	}
	system("pause");
	return 0;
}
