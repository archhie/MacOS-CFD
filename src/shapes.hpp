#pragma once
// Placeholder for shape definitions

#include "solver/grid.hpp"
#include "solver/state.hpp"
#include <vector>

// Forward declaration
struct Shape;

// Function to apply shapes to velocity field
void apply_shape_boundary_conditions(const Grid& grid, State& state, const std::vector<Shape>& shapes);
