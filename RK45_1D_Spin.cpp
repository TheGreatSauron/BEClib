#include "RK45_1D_Spin.h"

#include <iostream>

RK45_1D_Spin::RK45_1D_Spin(double step_size, double accuracy, std::valarray<comp> initial, double mu, double U) :
    RK45_1D(step_size, accuracy, initial, mu, U)
{
	// Ensure the vector doesn't have an odd number of entries
	if (y.size() % 2 != 0)
	{
		std::cout << "Error: Odd vector for 2-component spin vector.";

		y = {};
	}
}

std::valarray<comp> RK45_1D_Spin::func(std::valarray<comp> y1)
{
	int n = y1.size();
	std::valarray<comp> der = y1;

	for (int i = 0; i < n; i++)
	{
		// Even numbers component 1, odd component 2
		if (i % 2 == 0)
		{
			der[i] = (-mu * y1[i] - J * (y1[((i - 2) % n + n) % n] + y1[((i + 2) % n + n) % n]) + U * y1[i] * (std::norm(y[i]) + std::norm(y[i + 1]))) / i1;
		}
		else
		{
			der[i] = (-mu * y1[i] - J * (y1[((i - 2) % n + n) % n] + y1[((i + 2) % n + n) % n]) + U * y1[i] * (std::norm(y[i]) + std::norm(y[i - 1]))) / i1;
		}
	}

	return der;
}