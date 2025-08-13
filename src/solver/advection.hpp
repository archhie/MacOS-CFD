#pragma once

#include "grid.hpp"
#include "field.hpp"

void advect_u(const Grid&, const Field2D<double>& u, const Field2D<double>& v, Field2D<double>& du_dt);
void advect_v(const Grid&, const Field2D<double>& u, const Field2D<double>& v, Field2D<double>& dv_dt);

