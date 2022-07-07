#include <iostream>
#include <complex>
#include <valarray>
#include <fstream>

#include "RK45.h"
#include "RK45_1D.h"
#include "RK45_1D_Spin.h"
#include "RK45_2D_Spin.h"

int main()
{
	//std::valarray<comp> v(4);
	//for (int i = 0; i < v.size(); i++)
	//{
	//	if (i % 2 == 0)
	//	{
	//		v[i] = 1.0;
	//	}
	//}

	//RK45_1D_Spin rk(0.0001, 0.000000000001, v, 1.0, 1.0);

	//std::cout << rk.getTime() << ' ' << (1 - rk.getNorm()) << " | " << rk.printNorms() << '\n';

	//while (rk.getTime() <= 5.0)
	//{
	//	int n = rk.full_step();
	//
	//	std::cout << rk.getTime() << ' ' << std::abs(1 - rk.getNorm()) << " | " << rk.printNorms() << '\n';
	//}

	//rk.groundState();
	//std::cout << rk.getMu() << ' ' << (1 - rk.getNorm()) << " | " << rk.printNorms() << '\n';

	//std::cout << rk.runTime(10, true);

	int L = 60;

	std::valarray<comp> v(2*L*L);
	v[0] = 1.0;
	for (int i = 0; i < v.size(); i++)
	{
		v[i] = double(std::rand() % 10000) / 10000;
	}
	//for (int i = 0; i < v.size(); i++)
	//{
	//	if (i % 2 == 1)
	//	{
	//		v[i] = 0.001;
	//	}
	//	else
	//	{
	//		v[i] = 1.0;
	//	}
	//}
	
	for (int n = 0; n < 9; n++)
	{
		double args[6] = { 2.0 * L * L, 1.0, (n / 3), 0.1, 0.1, -0.401 + (0.001 * (n % 3)) };

		RK45_2D_Spin rk(0.0001, 0.000000000001, v, L, args);

		rk.groundState();

		std::ofstream file;
		std::string filename = (
			"n(" + std::to_string(args[0] / (L * L)) +
			")_U(" + std::to_string(args[2]) +
			")_tz(" + std::to_string(args[3]) +
			")_tso(" + std::to_string(args[4]) +
			")_Mz(" + std::to_string(args[5]) + ").txt");
		//std::string filename = "Test.txt";

		std::cout << filename << '\n';
		file.open(filename);

		file << rk.getMu() << '\n';

		for (int i = 0; i < rk.getVector().size(); i++)
		{
			file << rk.getVector()[i] << ' ';
		}

		file.close();
	}

	return 0;
}