#pragma once

#include <vector>

#include "solver/grid.hpp"
#include "solver/state.hpp"

enum class ScalarField { U, V, Speed, Pressure, Vorticity, Streamlines, Temperature };

// Compute a scalar field for visualization. Returns min/max in vmin/vmax.
void compute_scalar(const Grid &g, const State &s, ScalarField field,
                    std::vector<double> &out, double &vmin, double &vmax);

// Compute streamlines (simplified implementation)
void compute_streamlines(const Grid &g, const State &s, std::vector<double> &out, 
                        double &vmin, double &vmax, int num_streamlines = 50);
