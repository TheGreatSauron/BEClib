#ifndef RK45_1D_H
#define RK45_1D_H

#include "RK45.h"

namespace bec
{

// Extension of RK45 for 1 dimensional Bose-Hubbard model
class RK45_1D : 
	public RK45
{
public:

	// Standard constructor
	RK45_1D(double step_size, double accuracy, std::valarray<comp> initial, double norm, double mu, double U) :
		RK45(step_size, accuracy, initial), mu(mu), U(U), J(1.0), N(norm)
	{
		// Ensure the initial wave-function is normalized
		double sum = 0.0;

		for (int i = 0; i < y.size(); i++)
		{
			sum += std::norm(y[i]);
		}

		if (sum != 0.0)
		{
			for (int i = 0; i < y.size(); i++)
			{
				y[i] *= std::sqrt(N / sum);
			}
		}
	}

	// Override of function to be integrated -> [H psi]/i
	virtual std::valarray<comp> func(std::valarray<comp> y1) override
	{
		int n = y1.size();
		std::valarray<comp> der = y1;

		// Calculate derivative for every term
		for (int i = 0; i < n; i++)
		{
			der[i] = (-mu * y1[i] - J * (y1[((i - 1) % n + n) % n] + y1[((i + 1) % n + n) % n]) + U * y1[i] * std::norm(y[i])) / i1;
		}

		return der;
	}

	// Imaginary time evolution to find the ground-state wave-function,
	// sets current time to zero at completion
	void groundState()
	{
		double mu_diff, norm_diff;
		double param = getAcc() * 10000.0;
		int n = 0;

		do {
			// Ensure that the imaginary time propagation eventually ends,
			// limit is arbitrary
			if (n > 10000)
			{
				std::cout << "Warning: step overflow, no convergence.\n";
				break;
			}

			norm_diff = 0.0;
			mu_diff = 0.0;
			double old_mu = mu;

			n += full_step(comp(0.0, -1.0));

			// Flip the sign on mu and reset it if mu becomes too small
			if (std::abs(mu) <= std::pow(10.0, -6))
			{
				mu = -mu / std::abs(mu);
			}

			// Rescale mu every step
			if (mu > 0)
			{
				mu = mu / (getNorm() / N);
			}
			else
			{
				mu = mu * (getNorm() / N);
			}

			mu_diff = std::abs(mu - old_mu);
			norm_diff = std::abs(1.0 - (getNorm() / N));

			//std::cout << func(y1)[0] << ' ' << getMu() << ' ' << (1 - getNorm()) << " | " << printNorms() << '\n';

			// Ensure normalization
			y *= std::sqrt(N / getNorm());

		} while ((norm_diff > param) || (mu_diff > param));

		setTime(0.0);
	}

	// Returns mu
	double getMu() const
	{
		return mu;
	}
	// Returns U
	double getU() const
	{
		return U;
	}
	// Returns N
	double getN() const
	{
		return N;
	}

protected:

	// Constants
	double mu, U, J, N;
};

} // namespace bec

#endif // RK45_1D_H
