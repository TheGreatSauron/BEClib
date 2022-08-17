#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <complex>
#include <valarray>
#include <fstream>

#include "Bose_Hubbard_RK45.h"

int main()
{
	int L = 60;	
	int len_args = 2;
	double dt = 1e-4;
	double acc = 1e-12;

	std::valarray<comp> initial(2*L*L);
	//initial[0] = 1.0;
	for (int i = 0; i < initial.size(); i++)
	{
		//initial[i] = double(std::rand() % 10000) / 10000;
		//initial[i] = 1.0;

		// Analytical ground-state (Mz < 4*tz)
		if (i % 2 == 1)
		{
			initial[i] = 0.0;
		}
		else
		{
			initial[i] = std::polar(1.0, M_PI * (((i /2)%L) + ((i /2)/L)));
		}

		// Analytical ground-state (Mz > 4*tz)
		//if (i % 2 == 0)
		//{
		//	initial[i] = 0.0;
		//}
		//else
		//{
		//	initial[i] = std::polar(1.0, M_PI * ((((i-1) / 2) % L) + (((i - 1) / 2) / L)));
		//}
	}
	
	double args[8] = { 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	RK45_2D_Spin rk(dt, acc, initial, L, args);

	for (int n = 0; n < 6; n++)
	{
		//double args[8] = { 1.0 * L * L, 1.0, 0.0, 0.1, 0.1, 0.2, -M_PI + (n%5) * (M_PI/2.0), -M_PI + (n/5) * (M_PI / 2.0) };
		//double args[8] = { 1.0, 1.0, 0.0, 0.1, 0.1, 0.2, M_PI, M_PI };
		double newArgs[8] = {1.0 * L * L, 1.0, n * 0.1, 0.1, 0.1, 0.2, M_PI, M_PI};
		//double args[8] = { 1.0 * L * L, 1.0, 0.0, (n * 0.1/5.0), (n * 0.1/5.0), (n * 0.2/5.0), M_PI, M_PI};

		// This step resets the vector for the method, 
		// do not activvate if you desire a recursive feedback between ground-states
		//rk.setVector(initial);
		
		rk.setArgs(newArgs);
		rk.groundState();

		std::ofstream file;
		std::string filename = (
			"n(" + std::to_string(newArgs[0] / (L * L)) +
			")_U(" + std::to_string(newArgs[2]) +
			")_tz(" + std::to_string(newArgs[3]) +
			")_tso(" + std::to_string(newArgs[4]) +
			")_Mz(" + std::to_string(newArgs[5]) +
			")_kx(" + std::to_string(newArgs[6]) +
			")_ky(" + std::to_string(newArgs[7]) +
			").txt");
		//std::string filename = "Test.txt";

		std::cout << filename << '\n';
		file.open(filename);
		file.precision(18);

		file << rk.getMu() << '\n';
		file << std::scientific;

		for (int i = 0; i < rk.getVector().size(); i++)
		{
			comp z = rk.getVector()[i];

			if (z.imag() >= 0)
			{
				file << "(" << z.real() << "+" << z.imag() << "j) ";
			}
			else
			{
				file << "(" << z.real() << "-" << std::abs(z.imag()) << "j) ";
			}
		}

		file.close();
	}

	return 0;
}
