#include "viz.hpp"

#include <algorithm>
#include <cmath>

void compute_scalar(const Grid &g, const State &s, ScalarField field,
                    std::vector<double> &out, double &vmin, double &vmax) {
    int nx = g.nx;
    int ny = g.ny;
    out.resize(static_cast<size_t>(nx) * ny);
    vmin = 1e30;
    vmax = -1e30;
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            double val = 0.0;
            
            switch (field) {
            case ScalarField::U:
                val = 0.5 * (s.u.at_raw(ii, jj) + s.u.at_raw(ii + 1, jj));
                break;
            case ScalarField::V:
                val = 0.5 * (s.v.at_raw(ii, jj) + s.v.at_raw(ii, jj + 1));
                break;
            case ScalarField::Speed: {
                double uc = 0.5 * (s.u.at_raw(ii, jj) + s.u.at_raw(ii + 1, jj));
                double vc = 0.5 * (s.v.at_raw(ii, jj) + s.v.at_raw(ii, jj + 1));
                val = std::sqrt(uc * uc + vc * vc);
                break;
            }
            case ScalarField::Pressure:
                val = s.p.at_raw(ii, jj);
                break;
            case ScalarField::Vorticity: {
                double dv_dx = (s.v.at_raw(ii + 1, jj) - s.v.at_raw(ii - 1, jj)) /
                               (2.0 * g.dx);
                double du_dy = (s.u.at_raw(ii, jj + 1) - s.u.at_raw(ii, jj - 1)) /
                               (2.0 * g.dy);
                val = dv_dx - du_dy;
                break;
            }
            case ScalarField::Streamlines:
                // For streamlines, we'll use a simplified approach
                // This is a placeholder - actual streamlines would require integration
                val = 0.0; // Placeholder
                break;
            case ScalarField::Temperature:
                // Placeholder for temperature field (could be pressure or other scalar)
                val = s.p.at_raw(ii, jj); // Using pressure as placeholder
                break;
            }
            
            if (!std::isfinite(val))
                val = 0.0;
            out[i + j * nx] = val;
            vmin = std::min(vmin, val);
            vmax = std::max(vmax, val);
        }
    }
    
    double span = std::max(1e-6, vmax - vmin);
    vmax = vmin + span;
}

void compute_streamlines(const Grid &g, const State &s, std::vector<double> &out, 
                        double &vmin, double &vmax, int num_streamlines) {
    int nx = g.nx;
    int ny = g.ny;
    out.resize(static_cast<size_t>(nx) * ny);
    vmin = 1e30;
    vmax = -1e30;
    
    // Simplified streamline computation using velocity magnitude
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            
            double uc = 0.5 * (s.u.at_raw(ii, jj) + s.u.at_raw(ii + 1, jj));
            double vc = 0.5 * (s.v.at_raw(ii, jj) + s.v.at_raw(ii, jj + 1));
            double speed = std::sqrt(uc * uc + vc * vc);
            
            // Simple streamline density based on velocity magnitude
            double val = speed;
            
            if (!std::isfinite(val))
                val = 0.0;
            out[i + j * nx] = val;
            vmin = std::min(vmin, val);
            vmax = std::max(vmax, val);
        }
    }
    
    double span = std::max(1e-6, vmax - vmin);
    vmax = vmin + span;
}
