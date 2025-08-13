#pragma once

#include "grid.hpp"
#include "field.hpp"

double compute_cfl(const Grid&, const Field2D<double>& u, const Field2D<double>& v, double Re, double CFL_target);

// Compute maximum speed magnitude at cell centers.
double max_velocity(const Grid&, const Field2D<double>& u, const Field2D<double>& v);

// L2 norm of divergence of velocity field (after projection should be small)
double divergence_l2(const Grid&, const Field2D<double>& u,
                     const Field2D<double>& v);

