#include <iostream>
#include <complex>
#include <valarray>

#include "RK45.h"
#include "RK45_1D.h"
#include "RK45_1D_Spin.h"

typedef std::complex<double> comp;

std::complex<double> dpdt(std::complex<double> psi, std::complex<double> psi_1, double psi_2, double Q[3])
{
	std::complex<double> sum = (Q[0] + Q[2]) * psi + Q[1] * (psi_1 + psi_2) + 2 * Q[2] * psi * std::norm(psi);

	return sum;
}

int main()
{
	/*RK45 rk(0.1, 0.000001, {1});

	std::cout << rk.getTime() << ' ' << rk.getNorm() << " | " << rk.printVector() << ' ' << std::exp(rk.getTime()) << '\n';

	for (int i = 0; i < 20; i++)
	{
		rk.full_step();

		std::cout << rk.getTime() << ' ' << rk.getNorm() << " | " << rk.printVector() << ' ' << std::exp(rk.getTime()) << '\n';
	}*/

	std::valarray<comp> v(4);
	for (int i = 0; i < v.size(); i++)
	{
		if (i % 2 == 0)
		{
			v[i] = 1.0;
		}
	}

	RK45_1D_Spin rk(0.0001, 0.000000000001, v, 1.0, 1.0);

	std::cout << rk.getTime() << ' ' << (1 - rk.getNorm()) << " | " << rk.printNorms() << '\n';

	while (rk.getTime() <= 5.0)
	{
		int n = rk.full_step();
	
		std::cout << rk.getTime() << ' ' << std::abs(1 - rk.getNorm()) << " | " << rk.printNorms() << '\n';
	}

	//rk.groundState();
	//std::cout << rk.getMu() << ' ' << (1 - rk.getNorm()) << " | " << rk.printNorms() << '\n';

	//std::cout << rk.runTime(10, true);

	return 0;
}