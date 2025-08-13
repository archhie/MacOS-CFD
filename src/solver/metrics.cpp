#include "metrics.hpp"

#include <algorithm>
#include <cmath>

double compute_cfl(const Grid& g, const Field2D<double>& u, const Field2D<double>& v, double Re, double CFL_target) {
    double umax = 0.0;
#pragma omp parallel for collapse(2) reduction(max : umax)
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.u_nx(); ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            umax = std::max(umax, std::abs(u.at_raw(ii, jj)));
        }
    }
    double vmax = 0.0;
#pragma omp parallel for collapse(2) reduction(max : vmax)
    for (int j = 0; j < g.v_ny(); ++j) {
        for (int i = 0; i < g.nx; ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            vmax = std::max(vmax, std::abs(v.at_raw(ii, jj)));
        }
    }

    double dt_diff = 0.5 * Re * std::min(g.dx * g.dx, g.dy * g.dy);
    double dt_adv = 1e30;
    
    // Only use advective limit if velocities are non-zero and reasonable
    if (umax > 1e-12 && umax < 1e6)
        dt_adv = std::min(dt_adv, g.dx / umax);
    if (vmax > 1e-12 && vmax < 1e6)
        dt_adv = std::min(dt_adv, g.dy / vmax);
    
    // If advective limit is still 1e30 (no valid velocities), use diffusion limit
    if (dt_adv >= 1e30)
        dt_adv = dt_diff;
    
    // Ensure we have a reasonable minimum timestep
    double dt_min = 1e-6;  // Increased minimum timestep
    double dt = std::max(std::min(dt_adv, dt_diff), dt_min);
    
    // Apply CFL target
    dt = CFL_target * dt;
    
    // Final safety check
    dt = std::clamp(dt, dt_min, 1e-3);
    
    return dt;
}

double max_velocity(const Grid& g, const Field2D<double>& u, const Field2D<double>& v) {
    double vmax = 0.0;
#pragma omp parallel for collapse(2) reduction(max : vmax)
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            double uc = 0.5 * (u.at_raw(ii, jj) + u.at_raw(ii + 1, jj));
            double vc = 0.5 * (v.at_raw(ii, jj) + v.at_raw(ii, jj + 1));
            double speed = std::sqrt(uc * uc + vc * vc);
            vmax = std::max(vmax, speed);
        }
    }
    return vmax;
}

double divergence_l2(const Grid& g, const Field2D<double>& u,
                     const Field2D<double>& v) {
    double idx = 1.0 / g.dx;
    double idy = 1.0 / g.dy;
    double sum = 0.0;
#pragma omp parallel for collapse(2) reduction(+ : sum)
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            double divx = u.at_raw(ii + 1, jj) - u.at_raw(ii, jj);
            double divy = v.at_raw(ii, jj + 1) - v.at_raw(ii, jj);
            double d = divx * idx + divy * idy;
            sum += d * d;
        }
    }
    double denom = static_cast<double>(g.nx) * g.ny;
    return std::sqrt(sum / denom);
}
