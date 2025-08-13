#pragma once

#include "grid.hpp"
#include "field.hpp"

// Compute cell-centered divergence of the MAC velocity field.
void divergence(const Grid &g, const Field2D<double> &u,
                const Field2D<double> &v, Field2D<double> &rhs);

// Velocity correction: u -= dt * grad p
void subtract_grad_p(const Grid &g, Field2D<double> &u, Field2D<double> &v,
                     const Field2D<double> &p, double dt);

