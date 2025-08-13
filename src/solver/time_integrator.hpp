#pragma once

#include "advection.hpp"
#include "bc.hpp"
#include "diffusion.hpp"
#include "grid.hpp"
#include "metrics.hpp"
#include "pressure.hpp"
#include "project.hpp"
#include "state.hpp"

// Helper function to sanitize field values (replace non-finite with 0)
void sanitize_field(Field2D<double>& field);

// RK2 (Heun) time integrator with projection for incompressible flow.
struct TimeIntegrator {
    const Grid &g;
    Field2D<double> du_dt, dv_dt;  // work arrays for derivatives
    Field2D<double> u1, v1;        // intermediate velocities

    explicit TimeIntegrator(const Grid &grid);

    // Zero out all work arrays
    void clear_work();

    // Advance one step. Returns chosen dt. pressure_residual is output.
    double step(State &s, const BC &bc, double Re, double CFL,
                PressureSolver &pressure, const PressureParams &pp,
                double &pressure_residual, double dt_override = 0.0);
};

