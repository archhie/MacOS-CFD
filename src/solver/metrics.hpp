#pragma once

#include "grid.hpp"
#include "field.hpp"

double compute_cfl(const Grid&, const Field2D<double>& u, const Field2D<double>& v, double Re, double CFL_target);

