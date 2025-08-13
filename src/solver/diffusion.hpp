#pragma once

#include "grid.hpp"
#include "field.hpp"

void diffuse_u(const Grid&, const Field2D<double>& u, Field2D<double>& du_dt, double invRe);
void diffuse_v(const Grid&, const Field2D<double>& v, Field2D<double>& dv_dt, double invRe);

