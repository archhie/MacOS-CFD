#include "diffusion.hpp"

#include <algorithm>

void diffuse_u(const Grid& g, const Field2D<double>& u, Field2D<double>& du_dt, double invRe) {
    double idx2 = 1.0 / (g.dx * g.dx);
    double idy2 = 1.0 / (g.dy * g.dy);
    int ngx = g.ngx;
    int ngy = g.ngy;

#pragma omp parallel for collapse(2)
    for (int j = 1; j < g.ny - 1; ++j) {
        for (int i = 1; i < g.u_nx() - 1; ++i) {
            int ii = i + ngx;
            int jj = j + ngy;
            double lap = (u.at_raw(ii + 1, jj) - 2.0 * u.at_raw(ii, jj) + u.at_raw(ii - 1, jj)) * idx2 +
                         (u.at_raw(ii, jj + 1) - 2.0 * u.at_raw(ii, jj) + u.at_raw(ii, jj - 1)) * idy2;
            du_dt.at_raw(ii, jj) += invRe * lap;
        }
    }
}

void diffuse_v(const Grid& g, const Field2D<double>& v, Field2D<double>& dv_dt, double invRe) {
    double idx2 = 1.0 / (g.dx * g.dx);
    double idy2 = 1.0 / (g.dy * g.dy);
    int ngx = g.ngx;
    int ngy = g.ngy;

#pragma omp parallel for collapse(2)
    for (int j = 1; j < g.v_ny() - 1; ++j) {
        for (int i = 1; i < g.nx - 1; ++i) {
            int ii = i + ngx;
            int jj = j + ngy;
            double lap = (v.at_raw(ii + 1, jj) - 2.0 * v.at_raw(ii, jj) + v.at_raw(ii - 1, jj)) * idx2 +
                         (v.at_raw(ii, jj + 1) - 2.0 * v.at_raw(ii, jj) + v.at_raw(ii, jj - 1)) * idy2;
            dv_dt.at_raw(ii, jj) += invRe * lap;
        }
    }
}
