#include "RK45_1D_Spin.h"

RK45_1D_Spin::RK45_1D_Spin(double step_size, double accuracy, std::valarray<comp> initial, double mu, double U) :
    RK45_1D(step_size, accuracy, initial, mu, U)
{
}

std::valarray<comp> RK45_1D_Spin::func(std::valarray<comp> y1)
{
    // Even numbers component 1, odd component 2

	int n = y1.size();
	std::valarray<comp> y2 = y1;

	for (int i = 0; i < n; i++)
	{
		if (i == 0)
		{
			y2[i] = (-mu * y1[i] - J * (y1[n - 1] + y1[(i + 1) % n]) + U * y1[i] * std::norm(y[i])) / (comp(0.0, 1.0));
		}
		else
		{
			y2[i] = (-mu * y1[i] - J * (y1[(i - 1) % n] + y1[(i + 1) % n]) + U * y1[i] * std::norm(y[i])) / (comp(0.0, 1.0));
		}
	}

	return y2;
}
