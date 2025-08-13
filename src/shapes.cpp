#include "shapes.hpp"
#include "solver/grid.hpp"
#include "solver/state.hpp"
#include "gui/gui.hpp"
#include <cmath>

void apply_shape_boundary_conditions(const Grid& grid, State& state, const std::vector<Shape>& shapes) {
    for (const auto& shape : shapes) {
        if (!shape.enabled) continue;
        
        // Convert shape coordinates to grid indices
        int start_i = static_cast<int>(shape.pos.x);
        int start_j = static_cast<int>(shape.pos.y);
        int end_i, end_j;
        
        if (shape.type == 0) { // Rectangle
            end_i = static_cast<int>(shape.pos.x + shape.size.x);
            end_j = static_cast<int>(shape.pos.y + shape.size.y);
        } else { // Circle
            float radius = shape.size.x;
            end_i = static_cast<int>(shape.pos.x + radius);
            end_j = static_cast<int>(shape.pos.y + radius);
        }
        
        // Clamp to grid bounds
        start_i = std::max(0, std::min(start_i, grid.nx - 1));
        start_j = std::max(0, std::min(start_j, grid.ny - 1));
        end_i = std::max(0, std::min(end_i, grid.nx - 1));
        end_j = std::max(0, std::min(end_j, grid.ny - 1));
        
        // Apply boundary conditions to cells inside the shape
        for (int j = start_j; j <= end_j; ++j) {
            for (int i = start_i; i <= end_i; ++i) {
                // Check if point is inside shape
                bool inside = false;
                
                if (shape.type == 0) { // Rectangle
                    inside = (i >= start_i && i <= end_i && j >= start_j && j <= end_j);
                } else { // Circle
                    float center_x = shape.pos.x + shape.size.x / 2.0f;
                    float center_y = shape.pos.y + shape.size.y / 2.0f;
                    float radius = shape.size.x / 2.0f;
                    float dx = (i + 0.5f) - center_x;
                    float dy = (j + 0.5f) - center_y;
                    inside = (dx * dx + dy * dy <= radius * radius);
                }
                
                if (inside) {
                    int ii = i + grid.ngx;
                    int jj = j + grid.ngy;
                    
                    if (shape.bcType == 0) { // Solid wall
                        // Set velocity to zero for solid walls
                        if (i < grid.u_nx()) {
                            state.u.at_raw(ii, jj) = 0.0;
                        }
                        if (j < grid.v_ny()) {
                            state.v.at_raw(ii, jj) = 0.0;
                        }
                    } else { // Porous
                        // Reduce velocity by a factor for porous media
                        const double porous_factor = 0.1; // 90% reduction
                        if (i < grid.u_nx()) {
                            state.u.at_raw(ii, jj) *= porous_factor;
                        }
                        if (j < grid.v_ny()) {
                            state.v.at_raw(ii, jj) *= porous_factor;
                        }
                    }
                }
            }
        }
    }
}
