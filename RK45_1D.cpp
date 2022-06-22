#include "RK45_1D.h"

#include <iostream>

RK45_1D::RK45_1D(double step_size, double accuracy, std::valarray<comp> initial,
	double mu, double U)
	: RK45(step_size, accuracy, initial), mu(mu), U(U), J(1.0)
{
	double norm = 0.0;

	for (int i = 0; i < y.size(); i++)
	{
		norm += std::norm(y[i]);
	}

	if (norm != 0.0)
	{
		for (int i = 0; i < y.size(); i++)
		{
			y[i] /= std::sqrt(norm);
		}
	}
}

std::valarray<comp> RK45_1D::func(std::valarray<comp> y1)
{
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

void RK45_1D::groundState()
{
	double diff;

	do {
		diff = 0.0;
		std::valarray<comp> y1 = y;

		full_step(comp(0.0, -1.0));

		if (std::abs(mu) <= getAcc())
		{
			mu = -mu/std::abs(mu);
		}

		if (mu > 0)
		{
			mu = mu / getNorm();
		}
		else
		{ 
			mu = mu * getNorm();
		}

		//std::cout << func(y1)[0] << ' ' << getMu() << ' ' << (1 - getNorm()) << " | " << printNorms() << '\n';

		y /= std::sqrt(getNorm());

		y1 = y - y1;

		for (int i = 0; i < y1.size(); i++)
		{
			diff += std::norm(y1[i]);
		}

		diff = std::sqrt(diff / y1.size());
	} while (diff > getAcc());

	setTime(0.0);
}

double RK45_1D::getMu() const
{
	return mu;
}