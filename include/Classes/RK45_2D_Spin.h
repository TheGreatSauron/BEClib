#ifndef RK45_2D_SPIN_H
#define RK45_2D_SPIN_H

#include "RK45_1D.h"

namespace bec
{

// Extension of standard Bose-Hubbard model for 2D, two component condensate
class RK45_2D_Spin :
    public RK45_1D
{
public:

    // Standard constructor, x sets width of lattice in x-direction,
    // args array is {N, mu, U, tz, tso, Mz, kx, ky} in that order
    RK45_2D_Spin(double step_size, double accuracy, std::valarray<comp> initial, unsigned int L, double args[8]) :
        RK45_1D(step_size, accuracy, initial, args[0], args[1], args[2]), width(L), tz(args[3]), tso(args[4]), Mz(args[5]), kx(args[6]), ky(args[7])
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

    // Override of the derivative function -> [H psi]/i,
    // includes multiple new terms including spin-orbit coupling
    virtual std::valarray<comp> func(std::valarray<comp> y1) override
	{
		int n = y1.size();
		std::valarray<comp> der = y1;

		for (int i = 0; i < n; i++)
		{
			// Calculate each component individually
			comp H1, H2, H3, Hz, Hso, Hd;

			// First two terms are spin independant
			H1 = -mu * y1[i];
			H2 = -J * (std::polar(1.0, kx) * (y1[((i - 2) % n + n) % n] + y1[((i + 2) % n + n) % n])
				+ std::polar(1.0, ky) * (y1[((i - 2 * width) % n + n) % n] + y1[((i + 2 * width) % n + n) % n]));

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

    // Sets all arguments for the Hamiltonian, and renormalizes the vector
    // args array is {N, mu, U, tz, tso, Mz, kx, ky} in that order
    void setArgs(double args[8])
	{
		N = args[0];
		mu = args[1];
		U = args[2];
		tz = args[3];
		tso = args[4];
		Mz = args[5];
		kx = args[6];
		ky = args[7];

		// Ensure the wave-function is normalized
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

protected:

    // The width of the lattice
    unsigned int width;

    // Constants
    double tz, tso, Mz, kx, ky;
};

// Vector has first component on even indexes and second component on odd indexes, each row is listed sequentially

} // namespace bec

#endif // RK45_2D_SPIN_H
