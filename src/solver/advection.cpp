#include "advection.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

// TVD debug toggle
bool tvd_debug = false;

static inline double minmod(double a, double b) {
    return (a * b <= 0.0) ? 0.0 : (std::abs(a) < std::abs(b) ? a : b);
}

// Clip face states to local monotone bounds (no new extrema)
static inline void clip_face_states(double& uL, double& uR, double uim1, double ui, double uip1) {
    double min_val = std::min(uim1, std::min(ui, uip1));
    double max_val = std::max(uim1, std::max(ui, uip1));
    
    if (tvd_debug) {
        double overshoot_before = std::max(0.0, std::max(uL - max_val, min_val - uL));
        overshoot_before = std::max(overshoot_before, std::max(uR - max_val, min_val - uR));
        if (overshoot_before > 1e-10) {
            static double max_overshoot = 0.0;
            max_overshoot = std::max(max_overshoot, overshoot_before);
            std::fprintf(stderr, "TVD max overshoot before clipping: %g\n", max_overshoot);
        }
    }
    
    uL = std::clamp(uL, min_val, max_val);
    uR = std::clamp(uR, min_val, max_val);
}

void advect_u(const Grid& g, const Field2D<double>& u, const Field2D<double>& v, Field2D<double>& du_dt) {
    double idx = 1.0 / g.dx;
    double idy = 1.0 / g.dy;
    int ngx = g.ngx;
    int ngy = g.ngy;

    // Clear the output field first
    size_t sz = static_cast<size_t>(du_dt.pitch) * (du_dt.ny + 2 * du_dt.ngy);
    std::memset(du_dt.data, 0, sz * sizeof(double));

    // Interior pass - no branches in hot loops
#pragma omp parallel for collapse(2)
    for (int j = 1; j < g.ny - 1; ++j) {
        for (int i = 1; i < g.u_nx() - 1; ++i) {
            int ii = i + ngx;
            int jj = j + ngy;

            // X-direction advection (u-faces) - MUSCL with minmod limiter
            // Stencil: [i-2, i-1, i, i+1, i+2] for u at i+1/2 face
            double uim2 = u.at_raw(ii - 2, jj);
            double uim1 = u.at_raw(ii - 1, jj);
            double ui = u.at_raw(ii, jj);
            double uip1 = u.at_raw(ii + 1, jj);
            double uip2 = u.at_raw(ii + 2, jj);

            // MUSCL reconstruction for i+1/2 face
            double s_i = minmod(ui - uim1, uip1 - ui);
            double s_ip1 = minmod(uip1 - ui, uip2 - uip1);
            double uL = ui + 0.5 * s_i;
            double uR = uip1 - 0.5 * s_ip1;
            
            // Clip to local monotone bounds
            clip_face_states(uL, uR, uim1, ui, uip1);
            
            // Upwind flux using face-normal velocity (u at i+1/2 face)
            double u_face = ui;  // u-velocity at i+1/2 face
            double flux_x_p = (u_face > 0.0 ? uL : uR) * u_face;

            // MUSCL reconstruction for i-1/2 face
            double s_im1 = minmod(uim1 - uim2, ui - uim1);
            double uL_m = uim1 + 0.5 * s_im1;
            double uR_m = ui - 0.5 * s_i;
            
            // Clip to local monotone bounds
            clip_face_states(uL_m, uR_m, uim2, uim1, ui);
            
            // Upwind flux using face-normal velocity (u at i-1/2 face)
            double u_face_m = uim1;  // u-velocity at i-1/2 face
            double flux_x_m = (u_face_m > 0.0 ? uL_m : uR_m) * u_face_m;

            // Y-direction advection (v-faces) - MUSCL with minmod limiter
            // Stencil: [j-2, j-1, j, j+1, j+2] for v at j+1/2 face
            double ujm2 = u.at_raw(ii, jj - 2);
            double ujm1 = u.at_raw(ii, jj - 1);
            double uj = ui;
            double ujp1 = u.at_raw(ii, jj + 1);
            double ujp2 = u.at_raw(ii, jj + 2);

            // MUSCL reconstruction for j+1/2 face
            double sy_j = minmod(uj - ujm1, ujp1 - uj);
            double sy_jp1 = minmod(ujp1 - uj, ujp2 - ujp1);
            double uLy = uj + 0.5 * sy_j;
            double uRy = ujp1 - 0.5 * sy_jp1;
            
            // Clip to local monotone bounds
            clip_face_states(uLy, uRy, ujm1, uj, ujp1);
            
            // Upwind flux using face-normal velocity (v at j+1/2 face)
            double v_face = v.at_raw(ii, jj);  // v-velocity at j+1/2 face
            double flux_y_p = (v_face > 0.0 ? uLy : uRy) * v_face;

            // MUSCL reconstruction for j-1/2 face
            double sy_jm1 = minmod(ujm1 - ujm2, uj - ujm1);
            double uLy_m = ujm1 + 0.5 * sy_jm1;
            double uRy_m = uj - 0.5 * sy_j;
            
            // Clip to local monotone bounds
            clip_face_states(uLy_m, uRy_m, ujm2, ujm1, uj);
            
            // Upwind flux using face-normal velocity (v at j-1/2 face)
            double v_face_m = v.at_raw(ii, jj - 1);  // v-velocity at j-1/2 face
            double flux_y_m = (v_face_m > 0.0 ? uLy_m : uRy_m) * v_face_m;

            // Pure FV flux divergence: (F_{i+1/2}-F_{i-1/2})/dx + (F_{j+1/2}-F_{j-1/2})/dy
            double conv = (flux_x_p - flux_x_m) * idx + (flux_y_p - flux_y_m) * idy;
            du_dt.at_raw(ii, jj) = -conv;
        }
    }

    // Boundary pass - handle edges separately to avoid branches in hot loops
    // Left boundary (i=0)
#pragma omp parallel for
    for (int j = 1; j < g.ny - 1; ++j) {
        int ii = ngx;
        int jj = j + ngy;
        
        // X-direction: use first-order upwind at boundary
        double ui = u.at_raw(ii, jj);
        double uip1 = u.at_raw(ii + 1, jj);
        double u_face = ui;
        double flux_x_p = (u_face > 0.0 ? ui : uip1) * u_face;
        double flux_x_m = 0.0; // No flux from left of domain
        
        // Y-direction: same as interior
        double ujm1 = u.at_raw(ii, jj - 1);
        double uj = ui;
        double ujp1 = u.at_raw(ii, jj + 1);
        double v_face = v.at_raw(ii, jj);
        double flux_y_p = (v_face > 0.0 ? uj : ujp1) * v_face;
        double v_face_m = v.at_raw(ii, jj - 1);
        double flux_y_m = (v_face_m > 0.0 ? ujm1 : uj) * v_face_m;
        
        double conv = (flux_x_p - flux_x_m) * idx + (flux_y_p - flux_y_m) * idy;
        du_dt.at_raw(ii, jj) = -conv;
    }
    
    // Right boundary (i=nx-1)
#pragma omp parallel for
    for (int j = 1; j < g.ny - 1; ++j) {
        int ii = ngx + g.u_nx() - 1;
        int jj = j + ngy;
        
        // X-direction: use first-order upwind at boundary
        double uim1 = u.at_raw(ii - 1, jj);
        double ui = u.at_raw(ii, jj);
        double u_face_m = uim1;
        double flux_x_m = (u_face_m > 0.0 ? uim1 : ui) * u_face_m;
        double flux_x_p = 0.0; // No flux to right of domain
        
        // Y-direction: same as interior
        double ujm1 = u.at_raw(ii, jj - 1);
        double uj = ui;
        double ujp1 = u.at_raw(ii, jj + 1);
        double v_face = v.at_raw(ii, jj);
        double flux_y_p = (v_face > 0.0 ? uj : ujp1) * v_face;
        double v_face_m = v.at_raw(ii, jj - 1);
        double flux_y_m = (v_face_m > 0.0 ? ujm1 : uj) * v_face_m;
        
        double conv = (flux_x_p - flux_x_m) * idx + (flux_y_p - flux_y_m) * idy;
        du_dt.at_raw(ii, jj) = -conv;
    }
    
    // Bottom boundary (j=0)
#pragma omp parallel for
    for (int i = 1; i < g.u_nx() - 1; ++i) {
        int ii = i + ngx;
        int jj = ngy;
        
        // X-direction: same as interior
        double uim1 = u.at_raw(ii - 1, jj);
        double ui = u.at_raw(ii, jj);
        double uip1 = u.at_raw(ii + 1, jj);
        double u_face = ui;
        double flux_x_p = (u_face > 0.0 ? ui : uip1) * u_face;
        double u_face_m = uim1;
        double flux_x_m = (u_face_m > 0.0 ? uim1 : ui) * u_face_m;
        
        // Y-direction: use first-order upwind at boundary
        double uj = ui;
        double ujp1 = u.at_raw(ii, jj + 1);
        double v_face = v.at_raw(ii, jj);
        double flux_y_p = (v_face > 0.0 ? uj : ujp1) * v_face;
        double flux_y_m = 0.0; // No flux from bottom of domain
        
        double conv = (flux_x_p - flux_x_m) * idx + (flux_y_p - flux_y_m) * idy;
        du_dt.at_raw(ii, jj) = -conv;
    }
    
    // Top boundary (j=ny-1)
#pragma omp parallel for
    for (int i = 1; i < g.u_nx() - 1; ++i) {
        int ii = i + ngx;
        int jj = ngy + g.ny - 1;
        
        // X-direction: same as interior
        double uim1 = u.at_raw(ii - 1, jj);
        double ui = u.at_raw(ii, jj);
        double uip1 = u.at_raw(ii + 1, jj);
        double u_face = ui;
        double flux_x_p = (u_face > 0.0 ? ui : uip1) * u_face;
        double u_face_m = uim1;
        double flux_x_m = (u_face_m > 0.0 ? uim1 : ui) * u_face_m;
        
        // Y-direction: use first-order upwind at boundary
        double ujm1 = u.at_raw(ii, jj - 1);
        double uj = ui;
        double v_face_m = v.at_raw(ii, jj - 1);
        double flux_y_m = (v_face_m > 0.0 ? ujm1 : uj) * v_face_m;
        double flux_y_p = 0.0; // No flux to top of domain
        
        double conv = (flux_x_p - flux_x_m) * idx + (flux_y_p - flux_y_m) * idy;
        du_dt.at_raw(ii, jj) = -conv;
    }
}

void advect_v(const Grid& g, const Field2D<double>& u, const Field2D<double>& v, Field2D<double>& dv_dt) {
    double idx = 1.0 / g.dx;
    double idy = 1.0 / g.dy;
    int ngx = g.ngx;
    int ngy = g.ngy;

    // Clear the output field first
    size_t sz = static_cast<size_t>(dv_dt.pitch) * (dv_dt.ny + 2 * dv_dt.ngy);
    std::memset(dv_dt.data, 0, sz * sizeof(double));

    // Interior pass - no branches in hot loops
#pragma omp parallel for collapse(2)
    for (int j = 1; j < g.v_ny() - 1; ++j) {
        for (int i = 1; i < g.nx - 1; ++i) {
            int ii = i + ngx;
            int jj = j + ngy;

            // Y-direction advection (v-faces) - MUSCL with minmod limiter
            // Stencil: [j-2, j-1, j, j+1, j+2] for v at j+1/2 face
            double vjm2 = v.at_raw(ii, jj - 2);
            double vjm1 = v.at_raw(ii, jj - 1);
            double vj = v.at_raw(ii, jj);
            double vjp1 = v.at_raw(ii, jj + 1);
            double vjp2 = v.at_raw(ii, jj + 2);

            // MUSCL reconstruction for j+1/2 face
            double sy_j = minmod(vj - vjm1, vjp1 - vj);
            double sy_jp1 = minmod(vjp1 - vj, vjp2 - vjp1);
            double vL = vj + 0.5 * sy_j;
            double vR = vjp1 - 0.5 * sy_jp1;
            
            // Clip to local monotone bounds
            clip_face_states(vL, vR, vjm1, vj, vjp1);
            
            // Upwind flux using face-normal velocity (v at j+1/2 face)
            double v_face = vj;  // v-velocity at j+1/2 face
            double flux_y_p = (v_face > 0.0 ? vL : vR) * v_face;

            // MUSCL reconstruction for j-1/2 face
            double sy_jm1 = minmod(vjm1 - vjm2, vj - vjm1);
            double vL_m = vjm1 + 0.5 * sy_jm1;
            double vR_m = vj - 0.5 * sy_j;
            
            // Clip to local monotone bounds
            clip_face_states(vL_m, vR_m, vjm2, vjm1, vj);
            
            // Upwind flux using face-normal velocity (v at j-1/2 face)
            double v_face_m = vjm1;  // v-velocity at j-1/2 face
            double flux_y_m = (v_face_m > 0.0 ? vL_m : vR_m) * v_face_m;

            // X-direction advection (u-faces) - MUSCL with minmod limiter
            // Stencil: [i-2, i-1, i, i+1, i+2] for u at i+1/2 face
            double vim2 = v.at_raw(ii - 2, jj);
            double vim1 = v.at_raw(ii - 1, jj);
            double vi = vj;
            double vip1 = v.at_raw(ii + 1, jj);
            double vip2 = v.at_raw(ii + 2, jj);

            // MUSCL reconstruction for i+1/2 face
            double sx_i = minmod(vi - vim1, vip1 - vi);
            double sx_ip1 = minmod(vip1 - vi, vip2 - vip1);
            double vLx = vi + 0.5 * sx_i;
            double vRx = vip1 - 0.5 * sx_ip1;
            
            // Clip to local monotone bounds
            clip_face_states(vLx, vRx, vim1, vi, vip1);
            
            // Upwind flux using face-normal velocity (u at i+1/2 face)
            double u_face = u.at_raw(ii, jj);  // u-velocity at i+1/2 face
            double flux_x_p = (u_face > 0.0 ? vLx : vRx) * u_face;

            // MUSCL reconstruction for i-1/2 face
            double sx_im1 = minmod(vim1 - vim2, vi - vim1);
            double vLx_m = vim1 + 0.5 * sx_im1;
            double vRx_m = vi - 0.5 * sx_i;
            
            // Clip to local monotone bounds
            clip_face_states(vLx_m, vRx_m, vim2, vim1, vi);
            
            // Upwind flux using face-normal velocity (u at i-1/2 face)
            double u_face_m = u.at_raw(ii - 1, jj);  // u-velocity at i-1/2 face
            double flux_x_m = (u_face_m > 0.0 ? vLx_m : vRx_m) * u_face_m;

            // Pure FV flux divergence: (F_{i+1/2}-F_{i-1/2})/dx + (F_{j+1/2}-F_{j-1/2})/dy
            double conv = (flux_x_p - flux_x_m) * idx + (flux_y_p - flux_y_m) * idy;
            dv_dt.at_raw(ii, jj) = -conv;
        }
    }

    // Boundary pass - handle edges separately to avoid branches in hot loops
    // Left boundary (i=0)
#pragma omp parallel for
    for (int j = 1; j < g.v_ny() - 1; ++j) {
        int ii = ngx;
        int jj = j + ngy;
        
        // X-direction: use first-order upwind at boundary
        double vi = v.at_raw(ii, jj);
        double vip1 = v.at_raw(ii + 1, jj);
        double u_face = u.at_raw(ii, jj);
        double flux_x_p = (u_face > 0.0 ? vi : vip1) * u_face;
        double flux_x_m = 0.0; // No flux from left of domain
        
        // Y-direction: same as interior
        double vjm1 = v.at_raw(ii, jj - 1);
        double vj = vi;
        double vjp1 = v.at_raw(ii, jj + 1);
        double v_face = vj;
        double flux_y_p = (v_face > 0.0 ? vj : vjp1) * v_face;
        double v_face_m = vjm1;
        double flux_y_m = (v_face_m > 0.0 ? vjm1 : vj) * v_face_m;
        
        double conv = (flux_x_p - flux_x_m) * idx + (flux_y_p - flux_y_m) * idy;
        dv_dt.at_raw(ii, jj) = -conv;
    }
    
    // Right boundary (i=nx-1)
#pragma omp parallel for
    for (int j = 1; j < g.v_ny() - 1; ++j) {
        int ii = ngx + g.nx - 1;
        int jj = j + ngy;
        
        // X-direction: use first-order upwind at boundary
        double vim1 = v.at_raw(ii - 1, jj);
        double vi = v.at_raw(ii, jj);
        double u_face_m = u.at_raw(ii - 1, jj);
        double flux_x_m = (u_face_m > 0.0 ? vim1 : vi) * u_face_m;
        double flux_x_p = 0.0; // No flux to right of domain
        
        // Y-direction: same as interior
        double vjm1 = v.at_raw(ii, jj - 1);
        double vj = vi;
        double vjp1 = v.at_raw(ii, jj + 1);
        double v_face = vj;
        double flux_y_p = (v_face > 0.0 ? vj : vjp1) * v_face;
        double v_face_m = vjm1;
        double flux_y_m = (v_face_m > 0.0 ? vjm1 : vj) * v_face_m;
        
        double conv = (flux_x_p - flux_x_m) * idx + (flux_y_p - flux_y_m) * idy;
        dv_dt.at_raw(ii, jj) = -conv;
    }
    
    // Bottom boundary (j=0)
#pragma omp parallel for
    for (int i = 1; i < g.nx - 1; ++i) {
        int ii = i + ngx;
        int jj = ngy;
        
        // X-direction: same as interior
        double vim1 = v.at_raw(ii - 1, jj);
        double vi = v.at_raw(ii, jj);
        double vip1 = v.at_raw(ii + 1, jj);
        double u_face = u.at_raw(ii, jj);
        double flux_x_p = (u_face > 0.0 ? vi : vip1) * u_face;
        double u_face_m = u.at_raw(ii - 1, jj);
        double flux_x_m = (u_face_m > 0.0 ? vim1 : vi) * u_face_m;
        
        // Y-direction: use first-order upwind at boundary
        double vj = vi;
        double vjp1 = v.at_raw(ii, jj + 1);
        double v_face = vj;
        double flux_y_p = (v_face > 0.0 ? vj : vjp1) * v_face;
        double flux_y_m = 0.0; // No flux from bottom of domain
        
        double conv = (flux_x_p - flux_x_m) * idx + (flux_y_p - flux_y_m) * idy;
        dv_dt.at_raw(ii, jj) = -conv;
    }
    
    // Top boundary (j=ny-1)
#pragma omp parallel for
    for (int i = 1; i < g.nx - 1; ++i) {
        int ii = i + ngx;
        int jj = ngy + g.v_ny() - 1;
        
        // X-direction: same as interior
        double vim1 = v.at_raw(ii - 1, jj);
        double vi = v.at_raw(ii, jj);
        double vip1 = v.at_raw(ii + 1, jj);
        double u_face = u.at_raw(ii, jj);
        double flux_x_p = (u_face > 0.0 ? vi : vip1) * u_face;
        double u_face_m = u.at_raw(ii - 1, jj);
        double flux_x_m = (u_face_m > 0.0 ? vim1 : vi) * u_face_m;
        
        // Y-direction: use first-order upwind at boundary
        double vjm1 = v.at_raw(ii, jj - 1);
        double vj = vi;
        double v_face_m = vjm1;
        double flux_y_m = (v_face_m > 0.0 ? vjm1 : vj) * v_face_m;
        double flux_y_p = 0.0; // No flux to top of domain
        
        double conv = (flux_x_p - flux_x_m) * idx + (flux_y_p - flux_y_m) * idy;
        dv_dt.at_raw(ii, jj) = -conv;
    }
}
