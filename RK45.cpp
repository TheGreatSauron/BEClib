#include "RK45.h"

#include <iostream>

// Define static constants
double RK45::B[5][5] = {
	{2.0/9.0, 0.0, 0.0, 0.0, 0.0},
	{1.0/12.0, 1.0/4.0, 0.0, 0.0, 0.0},
	{69.0/128.0, -243.0/128.0, 135.0/64.0, 0.0, 0.0},
	{-17.0/12.0, 27.0/4.0, -27.0/5.0, 16.0/15.0, 0.0},
	{65.0/432.0, -5.0/16.0, 13.0/16.0, 4.0/27.0, 5.0/144.0} };

double RK45::CH[6] = {47.0/450.0, 0.0, 12.0/25.0, 32.0/225.0, 1.0/30.0, 6.0/25.0};
double RK45::CT[6] = {-1.0/150.0, 0.0, 3.0/100.0, -16.0/75.0, -1.0/20.0, 6.0/25.0};

RK45::RK45(double step_size, double accuracy, std::valarray<comp> initial)
	: t(0.0), h(step_size), acc(accuracy), y(initial)
{
}

int RK45::full_step(comp factor)
{
	int n = 0;
	std::valarray<comp> y2;
	double TE = 0.0;

	// Begin stepping forward
	do {
		// Track number of steps attempted
		n++;

		// Initialize derivative weights
		std::valarray<comp> k[6];

		y2 = y;
		TE = 0.0;

		// Marginally faster to calculate this way than through loops,
		// it cuts down on unnecessary calculations of zero
		k[0] = factor * h * func(y2);
		k[1] = factor * h * func(y2 + B[0][0] * k[0]);
		k[2] = factor * h * func(y2 + B[1][0] * k[0] + B[1][1] * k[1]);
		k[3] = factor * h * func(y2 + B[2][0] * k[0] + B[2][1] * k[1] + B[2][2] * k[2]);
		k[4] = factor * h * func(y2 + B[3][0] * k[0] + B[3][1] * k[1] + B[3][2] * k[2] + B[3][3] * k[3]);
		k[5] = factor * h * func(y2 + B[4][0] * k[0] + B[4][1] * k[1] + B[4][2] * k[2] + B[4][3] * k[3] + B[4][4] * k[4]);

		// Calculate error between orders
		std::valarray<comp> err(y2.size());

		for (int n = 0; n < 6; n++)
		{
			err += CT[n] * k[n];
		}

		for (int i = 0; i < err.size(); i++)
		{
			TE += std::norm(err[i]);
		}

		TE = std::sqrt(TE / err.size());

		if (TE <= acc)
		{
			// Increment time
			t += h;

			// Calculate new vector
			for (int n = 0; n < 6; n++)
			{
				y2 += CH[n] * k[n];
			}
		}

		// Update step size
		if (TE != 0.0)
		{
			h = 0.9 * h * std::pow(acc / TE, 1.0 / 5.0);
		}
	} while (TE > acc);

	// Set new vector
	y = y2;

	return n;
}

double RK45::getRMS() const
{
	double x2 = 0.0;
	double x = 0.0;

	for (int i = 0; i < y.size(); i++)
	{
		x2 += i*i * std::norm(y[i]);
		x += i * std::norm(y[i]);
	}

	return std::sqrt(x2 - x * x);
}

std::valarray<comp> RK45::getVector() const
{
	return y;
}

double RK45::getStepSize() const
{
	return h;
}

void RK45::setAcc(double accuracy)
{
	acc = accuracy;

	return;
}

void RK45::setTime(double time)
{
	t = time;

	return;
}

std::string RK45::printNorms() const
{
	std::string str = "";

	for (int i = 0; i < y.size(); i++)
	{
		str += std::to_string(std::norm(y[i]));
		str += ' ';
	}

	return str;
}

std::string RK45::printVector() const
{
	std::string str = "";

	for (int i = 0; i < y.size(); i++)
	{
		str += "(";
		str += std::to_string(std::real(y[i]));
		str += ',';
		str += std::to_string(std::imag(y[i]));
		str += ") ";
	}

	return str;
}

std::string RK45::printPhases() const
{
	std::string str = "";

	for (int i = 0; i < y.size(); i++)
	{
		str += std::to_string(std::arg(y[i]));
		str += " ";
	}

	return str;
}

void RK45::setStepSize(double step_size)
{
	h = step_size;
}

double RK45::getNorm() const
{
	double norm = 0.0;

	for (int i = 0; i < y.size(); i++)
	{
		norm += std::norm(y[i]);
	}
	
	return norm;
}

double RK45::getTime() const
{
	return t;
}

double RK45::getAcc() const
{
	return acc;
}

std::valarray<comp> RK45::func(std::valarray<comp> y1)
{
	// Default derivative function returns itself
	return y1;
}