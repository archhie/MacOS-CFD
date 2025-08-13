#pragma once

#include <vector>

#include "solver/grid.hpp"
#include "solver/state.hpp"

enum class ScalarField { U, V, Speed, Pressure, Vorticity };

// Compute a scalar field for visualization. Returns min/max in vmin/vmax.
void compute_scalar(const Grid &g, const State &s, ScalarField field,
                    std::vector<double> &out, double &vmin, double &vmax);
