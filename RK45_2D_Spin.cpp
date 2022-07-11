#include "RK45_2D_Spin.h"

#include <iostream>

RK45_2D_Spin::RK45_2D_Spin(double step_size, double accuracy, std::valarray<comp> initial, unsigned int L, double args[6]) :
	RK45_1D(step_size, accuracy, initial, args[0], args[1], args[2]), width(L), tz(args[3]), tso(args[4]), Mz(args[5])
{
	// Ensure the lattice properly fits the given dimensions
	if (y.size() % width != 0)
	{
		std::cout << "Error: Vector does not match given lattice width.\n";

		y = {};
	}

	// Ensure the vector has an even number of cells
	if (y.size() % 2 != 0)
	{
		std::cout << "Error: Odd vector for 2-component spin vector.\n";

		y = {};
	}
}

std::valarray<comp> RK45_2D_Spin::func(std::valarray<comp> y1)
{
	int n = y1.size();
	std::valarray<comp> der = y1;

	for (int i = 0; i < n; i++)
	{
		// Calculate each component individually
		comp H1, H2, H3, Hz, Hso, Hd;

		// First two terms are spin independant
		H1 = -mu * y1[i];
		H2 = -J * std::polar(1.0, M_PI) * (y1[((i - 2) % n + n) % n] + y1[((i + 2) % n + n) % n]
			+ y1[((i - 2 * width) % n + n) % n] + y1[((i + 2 * width) % n + n) % n]);

		// Even spin up, odd spin down
		// Index +1 means flipped down, index -1 means flipped up
		if (i % 2 == 0)
		{
			H3 = U * y1[i] * (std::norm(y[i]) + std::norm(y[i + 1]));
			Hz = tz * (y1[((i - 2) % n + n) % n] + y1[((i + 2) % n + n) % n]
				+ y1[((i - 2 * width) % n + n) % n] + y1[((i + 2 * width) % n + n) % n]);
			Hd = Mz * y1[i];
			Hso = tso * (y1[((i - 2 + 1) % n + n) % n] + y1[((i + 2 + 1) % n + n) % n]
				- y1[((i - 2 * width + 1) % n + n) % n] - y1[((i + 2 * width + 1) % n + n) % n]
				- i1 * y1[((i - 2 * width - 2 + 1) % n + n) % n] - i1 * y1[((i + 2 * width + 2 + 1) % n + n) % n]
				+ i1 * y1[((i + 2 * width - 2 + 1) % n + n) % n] + i1 * y1[((i - 2 * width + 2 + 1) % n + n) % n]);
		}
		else
		{
			H3 = U * y1[i] * (std::norm(y[i]) + std::norm(y[i - 1]));
			Hz = -tz * (y1[((i - 2) % n + n) % n] + y1[((i + 2) % n + n) % n]
				+ y1[((i - 2 * width) % n + n) % n] + y1[((i + 2 * width) % n + n) % n]);
			Hd = -Mz * y1[i];
			Hso = tso * (y1[((i - 2 - 1) % n + n) % n] + y1[((i + 2 - 1) % n + n) % n]
				- y1[((i - 2 * width - 1) % n + n) % n] - y1[((i + 2 * width - 1) % n + n) % n]
				+ i1 * y1[((i - 2 * width - 2 - 1) % n + n) % n] + i1 * y1[((i + 2 * width + 2 - 1) % n + n) % n]
				- i1 * y1[((i + 2 * width - 2 - 1) % n + n) % n] - i1 * y1[((i - 2 * width + 2 - 1) % n + n) % n]);
		}

		// Return Hamiltonian divided by imaginary unit
		der[i] = (H1 + H2 + H3 + Hz + Hd + Hso) / i1;
	}

	return der;
}
