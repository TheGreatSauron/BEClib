#pragma once

#include "RK45.h"

class RK45_1D : 
	public RK45
{
public:

	RK45_1D(double step_size, double accuracy, std::valarray<comp> initial, double mu, double U);

	//virtual comp func(int i, comp y1) override;
	virtual std::valarray<comp> func(std::valarray<comp> y1) override;

	void groundState();

	double getMu() const;
	double getU() const;

protected:

	double mu, U, J;
};
