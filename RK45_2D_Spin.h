#pragma once

#include "RK45_1D.h"

class RK45_2D_Spin :
    public RK45_1D
{
public:

    RK45_2D_Spin(double step_size, double accuracy, std::valarray<comp> initial, unsigned int x, double args[5]);

    virtual std::valarray<comp> func(std::valarray<comp> y1) override;

private:

    unsigned int width;

    double tz, tso, Mz;
};