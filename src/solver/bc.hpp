#pragma once

#include "grid.hpp"
#include "field.hpp"

enum class BCType { Wall, Moving, Inflow, Outflow, Periodic };

struct BC {
    BCType left = BCType::Wall;
    BCType right = BCType::Wall;
    BCType bottom = BCType::Wall;
    BCType top = BCType::Wall;
    double movingU = 1.0, movingV = 0.0;
    double inflowUx = 1.0, inflowUy = 0.0;
};

void apply_bc_u(const Grid&, Field2D<double>&, const BC&);
void apply_bc_v(const Grid&, Field2D<double>&, const BC&);
void apply_bc_p(const Grid&, Field2D<double>&, const BC&);

