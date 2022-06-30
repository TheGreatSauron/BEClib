#ifndef RK45_2D_SPIN_H
#define RK45_2D_SPIN_H

#include "RK45_1D.h"

// Extension of standard Bose-Hubbard model for 2D, two component condensate
class RK45_2D_Spin :
    public RK45_1D
{
public:

    // Standard constructor, x sets width of lattice in x-direction,
    // args array is {N, mu, U, tz, tso, Mz} in that order
    RK45_2D_Spin(double step_size, double accuracy, std::valarray<comp> initial, unsigned int x, double args[6]);

    // Override of the derivative function -> [H psi]/i,
    // includes multiple new terms including spin-orbit coupling
    virtual std::valarray<comp> func(std::valarray<comp> y1) override;

private:

    // The width of the lattice
    unsigned int width;

    // Constants
    double tz, tso, Mz;
};

// Vector has first component on even indexes and second component on odd indexes, each row is listed sequentially

#endif // RK45_2D_SPIN_H