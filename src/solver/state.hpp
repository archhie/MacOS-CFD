#pragma once

#include "field.hpp"

struct State {
    Field2D<double> u, v, p;    // primary fields
    Field2D<double> rhs, tmp;   // work buffers
    Field2D<float> scalar;      // visualization buffer
};

