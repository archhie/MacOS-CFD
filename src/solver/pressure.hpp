#pragma once

#include "field.hpp"
#include "grid.hpp"

enum class PressureSolverType { Multigrid, PCG };

struct PressureParams {
    PressureSolverType type = PressureSolverType::Multigrid;
    int mg_levels = 3;
    int vcycles = 2;
    int rbgs_sweeps = 2;
    double tol = 1e-6;
    int pcg_max_iters = 200;
};

class PressureSolver {
  public:
    explicit PressureSolver(const Grid &g);

    // Solve Laplace(p) = rhs. Returns final residual norm.
    double solve(Field2D<double> &p, const Field2D<double> &rhs,
                 const PressureParams &params);

  private:
    const Grid &g;
    Field2D<double> r, z, s, Ap; // PCG work buffers

    void apply_laplacian(const Field2D<double> &in,
                         Field2D<double> &out) const;
    double dot(const Field2D<double> &a, const Field2D<double> &b) const;
    double smooth(Field2D<double> &p, const Field2D<double> &rhs,
                  int sweeps);
};

