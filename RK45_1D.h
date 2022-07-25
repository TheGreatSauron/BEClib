#ifndef RK45_1D_H
#define RK45_1D_H

#include "RK45.h"

// Extension of RK45 for 1 dimensional Bose-Hubbard model
class RK45_1D : 
	public RK45
{
public:

	// Standard constructor
	RK45_1D(double step_size, double accuracy, std::valarray<comp> initial, double norm, double mu, double U);

	// Override of function to be integrated -> [H psi]/i
	virtual std::valarray<comp> func(std::valarray<comp> y1) override;

	// Imaginary time evolution to find the ground-state wave-function,
	// sets current time to zero at completion
	void groundState();

	// Returns mu
	double getMu() const;
	// Returns U
	double getU() const;
	// Returns N
	double getN() const;

protected:

	// Constants
	double mu, U, J, N;
};

#endif // RK45_1D_H