#include "project.hpp"

#include <cmath>
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
            double uR = u.at_raw(ii + 1, jj);
            double uL = u.at_raw(ii, jj);
            double vT = v.at_raw(ii, jj + 1);
            double vB = v.at_raw(ii, jj);
            if (!std::isfinite(uR)) uR = 0.0;
            if (!std::isfinite(uL)) uL = 0.0;
            if (!std::isfinite(vT)) vT = 0.0;
            if (!std::isfinite(vB)) vB = 0.0;
            double divx = uR - uL;
            double divy = vT - vB;
            double val = divx * idx + divy * idy;
            rhs.at_raw(ii, jj) = std::isfinite(val) ? val : 0.0;
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
            if (!std::isfinite(gradp)) gradp = 0.0;
            u.at_raw(ii, jj) -= dt * gradp;
        }
    }

#pragma omp parallel for collapse(2)
    for (int j = 0; j < g.v_ny(); ++j) {
        for (int i = 0; i < g.nx; ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            double gradp = (p.at_raw(ii, jj) - p.at_raw(ii, jj - 1)) * idy;
            if (!std::isfinite(gradp)) gradp = 0.0;
            v.at_raw(ii, jj) -= dt * gradp;
        }
    }
}

