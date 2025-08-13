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

    double dt_adv = 1e30;
    if (umax > 1e-12)
        dt_adv = std::min(dt_adv, g.dx / umax);
    if (vmax > 1e-12)
        dt_adv = std::min(dt_adv, g.dy / vmax);
    double dt_diff = 0.5 * Re * std::min(g.dx * g.dx, g.dy * g.dy);
    double dt = std::min(dt_adv, dt_diff);
    return CFL_target * dt;
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
