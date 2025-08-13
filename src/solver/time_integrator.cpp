#include "time_integrator.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

void sanitize_field(Field2D<double>& field) {
    size_t total = static_cast<size_t>(field.pitch) * (field.ny + 2 * field.ngy);
    bool has_nan = false;
    
#pragma omp parallel for reduction(|:has_nan)
    for (size_t i = 0; i < total; ++i) {
        if (!std::isfinite(field.data[i])) {
            field.data[i] = 0.0;
            has_nan = true;
        }
    }
    
    if (has_nan) {
        static bool warned = false;
        if (!warned) {
            std::fprintf(stderr, "Non-finite values detected and replaced with 0\n");
            warned = true;
        }
    }
}

TimeIntegrator::TimeIntegrator(const Grid &grid) : g(grid) {
    // allocate work arrays matching velocity layouts
    du_dt.allocate(g.u_nx(), g.ny, g.u_pitch(), g.ngx, g.ngy);
    dv_dt.allocate(g.nx, g.v_ny(), g.v_pitch(), g.ngx, g.ngy);
    u1.allocate(g.u_nx(), g.ny, g.u_pitch(), g.ngx, g.ngy);
    v1.allocate(g.nx, g.v_ny(), g.v_pitch(), g.ngx, g.ngy);
    clear_work();
}

void TimeIntegrator::clear_work() {
    size_t sz_u = static_cast<size_t>(du_dt.pitch) * (du_dt.ny + 2 * du_dt.ngy);
    size_t sz_v = static_cast<size_t>(dv_dt.pitch) * (dv_dt.ny + 2 * dv_dt.ngy);
    std::memset(du_dt.data, 0, sz_u * sizeof(double));
    std::memset(dv_dt.data, 0, sz_v * sizeof(double));
    std::memset(u1.data, 0, sz_u * sizeof(double));
    std::memset(v1.data, 0, sz_v * sizeof(double));
}

double TimeIntegrator::step(State &s, const BC &bc, double Re, double CFL,
                            PressureSolver &pressure, const PressureParams &pp,
                            double &pressure_residual, double dt_override) {
    // Enforce BCs before computing timestep
    apply_bc_u(g, s.u, bc);
    apply_bc_v(g, s.v, bc);
    apply_bc_p(g, s.p, bc);

    double dt_cfl = compute_cfl(g, s.u, s.v, Re, CFL);
    double dt = dt_cfl;
    if (dt_override > 0.0)
        dt = std::min(dt_cfl, dt_override);
    
    // CFL safety: clamp dt to reasonable bounds
    double dt_diff = 0.5 * Re * std::min(g.dx * g.dx, g.dy * g.dy);
    double dt_cap = (dt_override > 0.0) ? dt_override : 1e-3;
    
    if (!std::isfinite(dt) || dt <= 0.0)
        dt = std::min(1e-3, dt_diff);
    dt = std::clamp(dt, 1e-8, dt_cap);

    size_t sz_u = static_cast<size_t>(du_dt.pitch) * (du_dt.ny + 2 * du_dt.ngy);
    size_t sz_v = static_cast<size_t>(dv_dt.pitch) * (dv_dt.ny + 2 * dv_dt.ngy);

    std::memset(du_dt.data, 0, sz_u * sizeof(double));
    std::memset(dv_dt.data, 0, sz_v * sizeof(double));

    // First RK stage
    advect_u(g, s.u, s.v, du_dt);
    diffuse_u(g, s.u, du_dt, 1.0 / Re);
    advect_v(g, s.u, s.v, dv_dt);
    diffuse_v(g, s.v, dv_dt, 1.0 / Re);

    // Sanitize derivatives
    sanitize_field(du_dt);
    sanitize_field(dv_dt);

#pragma omp parallel for collapse(2)
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.u_nx(); ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            u1.at_raw(ii, jj) = s.u.at_raw(ii, jj) + dt * du_dt.at_raw(ii, jj);
        }
    }
#pragma omp parallel for collapse(2)
    for (int j = 0; j < g.v_ny(); ++j) {
        for (int i = 0; i < g.nx; ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            v1.at_raw(ii, jj) = s.v.at_raw(ii, jj) + dt * dv_dt.at_raw(ii, jj);
        }
    }

    // Sanitize intermediate velocities
    sanitize_field(u1);
    sanitize_field(v1);

    // Apply BCs after first stage
    apply_bc_u(g, u1, bc);
    apply_bc_v(g, v1, bc);

    std::memset(du_dt.data, 0, sz_u * sizeof(double));
    std::memset(dv_dt.data, 0, sz_v * sizeof(double));

    // Second RK stage
    advect_u(g, u1, v1, du_dt);
    diffuse_u(g, u1, du_dt, 1.0 / Re);
    advect_v(g, u1, v1, dv_dt);
    diffuse_v(g, v1, dv_dt, 1.0 / Re);

    // Sanitize derivatives
    sanitize_field(du_dt);
    sanitize_field(dv_dt);

#pragma omp parallel for collapse(2)
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.u_nx(); ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            double updated = 0.5 * (s.u.at_raw(ii, jj) +
                                   u1.at_raw(ii, jj) +
                                   dt * du_dt.at_raw(ii, jj));
            s.u.at_raw(ii, jj) = updated;
        }
    }
#pragma omp parallel for collapse(2)
    for (int j = 0; j < g.v_ny(); ++j) {
        for (int i = 0; i < g.nx; ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            double updated = 0.5 * (s.v.at_raw(ii, jj) +
                                   v1.at_raw(ii, jj) +
                                   dt * dv_dt.at_raw(ii, jj));
            s.v.at_raw(ii, jj) = updated;
        }
    }

    // Sanitize updated velocities
    sanitize_field(s.u);
    sanitize_field(s.v);

    // Apply BCs after RK2 update
    apply_bc_u(g, s.u, bc);
    apply_bc_v(g, s.v, bc);

    // Sanitize velocities before projection
    sanitize_field(s.u);
    sanitize_field(s.v);

    // Projection step: ∇²p = (1/dt)*div, with proper logging
    divergence(g, s.u, s.v, s.rhs);
    sanitize_field(s.rhs);
    
    // Log divergence before projection
    double div_before = 0.0;
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            double val = s.rhs.at_raw(ii, jj);
            div_before += val * val;
        }
    }
    div_before = std::sqrt(div_before / (g.nx * g.ny));
    
    // Scale RHS by 1/dt for pressure solve: ∇²p = (1/dt)*div
#pragma omp parallel for collapse(2)
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            double val = s.rhs.at_raw(ii, jj) * (1.0 / dt);
            s.rhs.at_raw(ii, jj) = std::isfinite(val) ? val : 0.0;
        }
    }

    // Pressure solve: ∇²p = (1/dt)*div with Neumann BCs everywhere
    // Pressure is pinned at p(0,0) = 0 to kill nullspace
    s.p.at_raw(g.ngx, g.ngy) = 0.0; // pin pressure at (0,0)
    pressure_residual = pressure.solve(s.p, s.rhs, pp);
    s.p.at_raw(g.ngx, g.ngy) = 0.0; // re-pin after solve
    
    // Velocity correction: subtract grad p from u/v at faces
    subtract_grad_p(g, s.u, s.v, s.p, dt);
    
    // Immediately after projection, re-apply Dirichlet inflow u(i=0,j) and ghost values
    if (bc.left.type == BCType::Inflow) {
        for (int j = 0; j < g.ny; ++j) {
            int jj = j + g.ngy;
            double y = (j + 0.5) * g.dy;
            double val = jet_velocity(g, bc, y);
            s.u.at_raw(g.ngx, jj) = val;      // inflow face
            s.u.at_raw(g.ngx - 1, jj) = val;  // ghost cell
        }
    }
    
    // Apply all boundary conditions after projection
    apply_bc_u(g, s.u, bc);
    apply_bc_v(g, s.v, bc);
    apply_bc_p(g, s.p, bc);

    // Log divergence after projection
    divergence(g, s.u, s.v, s.rhs);
    double div_after = 0.0;
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            double val = s.rhs.at_raw(ii, jj);
            div_after += val * val;
        }
    }
    div_after = std::sqrt(div_after / (g.nx * g.ny));
    
    // Log divergence reduction
    if (div_before > 1e-12) {
        double reduction = div_before / div_after;
        static int log_counter = 0;
        if (log_counter % 60 == 0) {  // Log every 60 steps
            std::fprintf(stderr, "Divergence: before=%.2e, after=%.2e, reduction=%.1fx\n", 
                        div_before, div_after, reduction);
        }
        log_counter++;
    }

    // Final sanitization
    sanitize_field(s.u);
    sanitize_field(s.v);
    sanitize_field(s.p);

    return dt;
}

