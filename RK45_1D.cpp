#include "RK45_1D.h"

#include <iostream>

RK45_1D::RK45_1D(double step_size, double accuracy, std::valarray<comp> initial, double mu, double U) : 
	RK45(step_size, accuracy, initial), mu(mu), U(U), J(1.0)
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
		y2[i] = (-mu * y1[i] - J * (y1[((i - 1) % n + n) % n] + y1[((i + 1) % n + n) % n]) + U * y1[i] * std::norm(y[i])) / i1;
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

		mu = mu / getNorm();
		y /= std::sqrt(getNorm());

		//std::cout << getMu() << ' ' << (1 - getNorm()) << " | " << printNorms() << '\n';

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

double RK45_1D::getU() const
{
	return U;
}
