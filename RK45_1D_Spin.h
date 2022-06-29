#ifndef RK45_1D_SPIN_H
#define RK45_1D_SPIN_H

#include "RK45_1D.h"

// Extension of standard Bose-Hubbard model, with two component condensate
class RK45_1D_Spin :
    public RK45_1D
{
public:

     // Standard constructor
    RK45_1D_Spin(double step_size, double accuracy, std::valarray<comp> initial, double mu, double U);

    // Override of derivative function
    virtual std::valarray<comp> func(std::valarray<comp> y1) override;
};

// Vector has first component on even indexes and second component on odd indexes

#endif // RK45_1D_SPIN_H