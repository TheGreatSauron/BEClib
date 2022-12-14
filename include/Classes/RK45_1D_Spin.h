#ifndef RK45_1D_SPIN_H
#define RK45_1D_SPIN_H

#include "RK45_1D.h"

namespace bec
{

// Extension of standard Bose-Hubbard model, with two component condensate
class RK45_1D_Spin :
    public RK45_1D
{
public:

     // Standard constructor
    RK45_1D_Spin(double step_size, double accuracy, std::valarray<comp> initial, double norm, double mu, double U) :
		RK45_1D(step_size, accuracy, initial, norm, mu, U)
	{
		// Ensure the vector doesn't have an odd number of entries
		if (y.size() % 2 != 0)
		{
			std::cout << "Error: Odd vector for 2-component spin vector.";

			y = {};
		}
	}

    // Override of derivative function
    virtual std::valarray<comp> func(std::valarray<comp> y1) override
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
};

// Vector has first component on even indexes and second component on odd indexes

} // namespace bec

#endif // RK45_1D_SPIN_H
