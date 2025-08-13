#include "project.hpp"

#include <cmath>
#include <cstring>
#include <cstdio>

void divergence(const Grid &g, const Field2D<double> &u,
                const Field2D<double> &v, Field2D<double> &rhs) {
    double idx = 1.0 / g.dx;
    double idy = 1.0 / g.dy;

    // MAC divergence at pressure centers: (u(i+1,j)-u(i,j))/dx + (v(i,j+1)-v(i,j))/dy
    // This includes all boundary faces properly
#pragma omp parallel for collapse(2)
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            
            // u-faces: u(i+1,j) - u(i,j) - includes boundary faces
            double uR = u.at_raw(ii + 1, jj);  // right face of pressure cell
            double uL = u.at_raw(ii, jj);      // left face of pressure cell
            
            // v-faces: v(i,j+1) - v(i,j) - includes boundary faces  
            double vT = v.at_raw(ii, jj + 1);  // top face of pressure cell
            double vB = v.at_raw(ii, jj);      // bottom face of pressure cell
            
            // Sanitize values
            if (!std::isfinite(uR)) uR = 0.0;
            if (!std::isfinite(uL)) uL = 0.0;
            if (!std::isfinite(vT)) vT = 0.0;
            if (!std::isfinite(vB)) vB = 0.0;
            
            // MAC divergence: ∇·u = ∂u/∂x + ∂v/∂y
            double divx = (uR - uL) * idx;
            double divy = (vT - vB) * idy;
            double div = divx + divy;
            
            rhs.at_raw(ii, jj) = std::isfinite(div) ? div : 0.0;
        }
    }
}

void subtract_grad_p(const Grid &g, Field2D<double> &u, Field2D<double> &v,
                     const Field2D<double> &p, double dt) {
    double idx = 1.0 / g.dx;
    double idy = 1.0 / g.dy;

    // Velocity correction: subtract grad p from u at faces
    // u-face at (i,j): u -= dt * (p(i,j) - p(i-1,j))/dx
#pragma omp parallel for collapse(2)
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.u_nx(); ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            
            // Pressure gradient at u-face: (p(i,j) - p(i-1,j))/dx
            double pR = p.at_raw(ii, jj);     // pressure at right of u-face
            double pL = p.at_raw(ii - 1, jj); // pressure at left of u-face
            double gradp = (pR - pL) * idx;
            
            if (!std::isfinite(gradp)) gradp = 0.0;
            u.at_raw(ii, jj) -= dt * gradp;
        }
    }

    // Velocity correction: subtract grad p from v at faces
    // v-face at (i,j): v -= dt * (p(i,j) - p(i,j-1))/dy
#pragma omp parallel for collapse(2)
    for (int j = 0; j < g.v_ny(); ++j) {
        for (int i = 0; i < g.nx; ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            
            // Pressure gradient at v-face: (p(i,j) - p(i,j-1))/dy
            double pT = p.at_raw(ii, jj);     // pressure at top of v-face
            double pB = p.at_raw(ii, jj - 1); // pressure at bottom of v-face
            double gradp = (pT - pB) * idy;
            
            if (!std::isfinite(gradp)) gradp = 0.0;
            v.at_raw(ii, jj) -= dt * gradp;
        }
    }
}

