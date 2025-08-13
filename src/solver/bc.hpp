#pragma once

#include "grid.hpp"
#include "field.hpp"

enum class BCType { Wall, Moving, Inflow, Outflow, Periodic };

struct BCSide {
    BCType type = BCType::Wall;
    double moving = 0.0;   // tangential speed for Moving
    double inflow_u = 0.0; // specified inflow velocities
    double inflow_v = 0.0;
};

struct BC {
    BCSide left{BCType::Wall, 0.0, 1.0, 0.0};
    BCSide right;
    BCSide bottom;
    BCSide top;
    // Jet inflow parameters (used on left when type==Inflow)
    double jet_center = 0.5; // y0
    double jet_width = 0.2;  // w
    double jet_thickness = 0.02; // smoothing thickness delta
    double jet_eps = 0.03;       // perturbation amplitude
    int jet_k = 8;               // perturbation wavenumber
    double jet_phase = 0.0;      // perturbation phase
};

// Jet velocity profile function
double jet_velocity(const Grid& g, const BC& bc, double y);

void apply_bc_u(const Grid&, Field2D<double>&, const BC&);
void apply_bc_v(const Grid&, Field2D<double>&, const BC&);
void apply_bc_p(const Grid&, Field2D<double>&, const BC&);

