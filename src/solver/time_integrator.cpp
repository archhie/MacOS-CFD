#include "time_integrator.hpp"

#include <cstring>

TimeIntegrator::TimeIntegrator(const Grid &grid) : g(grid) {
    // allocate work arrays matching velocity layouts
    du_dt.allocate(g.u_nx(), g.ny, g.u_pitch(), g.ngx, g.ngy);
    dv_dt.allocate(g.nx, g.v_ny(), g.v_pitch(), g.ngx, g.ngy);
    u1.allocate(g.u_nx(), g.ny, g.u_pitch(), g.ngx, g.ngy);
    v1.allocate(g.nx, g.v_ny(), g.v_pitch(), g.ngx, g.ngy);
}

double TimeIntegrator::step(State &s, const BC &bc, double Re, double CFL,
                            PressureSolver &pressure, const PressureParams &pp) {
    double dt = compute_cfl(g, s.u, s.v, Re, CFL);

    size_t sz_u = static_cast<size_t>(du_dt.pitch) * (du_dt.ny + 2 * du_dt.ngy);
    size_t sz_v = static_cast<size_t>(dv_dt.pitch) * (dv_dt.ny + 2 * dv_dt.ngy);

    std::memset(du_dt.data, 0, sz_u * sizeof(double));
    std::memset(dv_dt.data, 0, sz_v * sizeof(double));

    advect_u(g, s.u, s.v, du_dt);
    diffuse_u(g, s.u, du_dt, 1.0 / Re);
    advect_v(g, s.u, s.v, dv_dt);
    diffuse_v(g, s.v, dv_dt, 1.0 / Re);

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

    apply_bc_u(g, u1, bc);
    apply_bc_v(g, v1, bc);

    std::memset(du_dt.data, 0, sz_u * sizeof(double));
    std::memset(dv_dt.data, 0, sz_v * sizeof(double));

    advect_u(g, u1, v1, du_dt);
    diffuse_u(g, u1, du_dt, 1.0 / Re);
    advect_v(g, u1, v1, dv_dt);
    diffuse_v(g, v1, dv_dt, 1.0 / Re);

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

    apply_bc_u(g, s.u, bc);
    apply_bc_v(g, s.v, bc);

    // Projection
    divergence(g, s.u, s.v, s.rhs);
#pragma omp parallel for collapse(2)
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            s.rhs.at_raw(ii, jj) *= 1.0 / dt;
        }
    }

    pressure.solve(s.p, s.rhs, pp);
    apply_bc_p(g, s.p, bc);
    subtract_grad_p(g, s.u, s.v, s.p, dt);
    apply_bc_u(g, s.u, bc);
    apply_bc_v(g, s.v, bc);

    return dt;
}

