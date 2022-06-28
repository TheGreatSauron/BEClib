#include <iostream>
#include <complex>
#include <valarray>

#include "RK45.h"
#include "RK45_1D.h"
#include "RK45_1D_Spin.h"
#include "RK45_2D_Spin.h"

typedef std::complex<double> comp;

std::complex<double> dpdt(std::complex<double> psi, std::complex<double> psi_1, double psi_2, double Q[3])
{
	std::complex<double> sum = (Q[0] + Q[2]) * psi + Q[1] * (psi_1 + psi_2) + 2 * Q[2] * psi * std::norm(psi);

	return sum;
}

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

	std::valarray<comp> v(200);
	//v[0] = 1.0;
	for (int i = 0; i < v.size(); i++)
	{
		v[i] = double(std::rand() % 1000) / 1000;
	}

	double args[5] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
	RK45_2D_Spin rk(0.0001, 0.000000000001, v, 10, args);

	std::cout << (1 - rk.getNorm()) << " | " << rk.printNorms() << '\n';

	rk.groundState();

	std::cout << (1 - rk.getNorm()) << " | " << rk.printNorms() << '\n';
	std::cout << (1 - rk.getNorm()) << " | " << rk.printPhases() << '\n';

	return 0;
}