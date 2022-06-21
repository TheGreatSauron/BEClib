#pragma once

#include "RK45_1D.h"

class RK45_1D_Spin :
    public RK45_1D
{
public:

    RK45_1D_Spin(double step_size, double accuracy, std::valarray<comp> initial, double mu, double U);

    virtual std::valarray<comp> func(std::valarray<comp> y1) override;
};