#include "project.hpp"

#include <cstring>

void divergence(const Grid &g, const Field2D<double> &u,
                const Field2D<double> &v, Field2D<double> &rhs) {
    double idx = 1.0 / g.dx;
    double idy = 1.0 / g.dy;

#pragma omp parallel for collapse(2)
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            double divx = u.at_raw(ii + 1, jj) - u.at_raw(ii, jj);
            double divy = v.at_raw(ii, jj + 1) - v.at_raw(ii, jj);
            rhs.at_raw(ii, jj) = (divx * idx + divy * idy);
        }
    }
}

void subtract_grad_p(const Grid &g, Field2D<double> &u, Field2D<double> &v,
                     const Field2D<double> &p, double dt) {
    double idx = 1.0 / g.dx;
    double idy = 1.0 / g.dy;

#pragma omp parallel for collapse(2)
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.u_nx(); ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            double gradp = (p.at_raw(ii, jj) - p.at_raw(ii - 1, jj)) * idx;
            u.at_raw(ii, jj) -= dt * gradp;
        }
    }

#pragma omp parallel for collapse(2)
    for (int j = 0; j < g.v_ny(); ++j) {
        for (int i = 0; i < g.nx; ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            double gradp = (p.at_raw(ii, jj) - p.at_raw(ii, jj - 1)) * idy;
            v.at_raw(ii, jj) -= dt * gradp;
        }
    }
}

