#include "bc.hpp"

#include <algorithm>
#include <cmath>

// Smooth top-hat mask for jet inflow (0 outside, ~1 inside)
static inline double jet_mask(const Grid &g, const BC &bc, double y) {
    double delta = bc.jet_thickness;
    double y0 = bc.jet_center;
    double w = bc.jet_width;
    double a = (y - (y0 - 0.5 * w)) / delta;
    double b = (y - (y0 + 0.5 * w)) / delta;
    return 0.5 * (std::tanh(a) - std::tanh(b));
}

static inline double jet_velocity(const Grid &g, const BC &bc, double y) {
    double base = bc.left.inflow_u * jet_mask(g, bc, y);
    const double pi = 3.14159265358979323846;
    double pert = bc.jet_eps * bc.left.inflow_u *
                  std::sin(2.0 * pi * bc.jet_k * (y / g.Ly) + bc.jet_phase);
    return base + pert;
}

void apply_bc_u(const Grid &g, Field2D<double> &u, const BC &bc) {
    int nx = g.nx;
    int ny = g.ny;
    int ngx = g.ngx;
    int ngy = g.ngy;

    // Left/right boundaries for u (normal component)
    for (int j = 0; j < ny; ++j) {
        int jj = j + ngy;
        double y = (j + 0.5) * g.dy;

        // Left
        if (bc.left.type == BCType::Periodic &&
            bc.right.type == BCType::Periodic) {
            // handled after loop
        } else {
            switch (bc.left.type) {
            case BCType::Wall:
            case BCType::Moving:
                u.at_raw(ngx, jj) = 0.0;
                u.at_raw(ngx - 1, jj) = 0.0;
                break;
            case BCType::Inflow: {
                double val = jet_velocity(g, bc, y);
                u.at_raw(ngx, jj) = val;
                u.at_raw(ngx - 1, jj) = val; // mirror ghost
                break;
            }
            case BCType::Outflow: {
                double val = u.at_raw(ngx + 1, jj);
                u.at_raw(ngx, jj) = val;
                u.at_raw(ngx - 1, jj) = val;
                break;
            }
            case BCType::Periodic:
                break;
            }
        }

        // Right
        switch (bc.right.type) {
        case BCType::Wall:
        case BCType::Moving:
            u.at_raw(ngx + nx, jj) = 0.0;
            u.at_raw(ngx + nx + 1, jj) = 0.0;
            break;
        case BCType::Inflow: {
            double val = bc.right.inflow_u;
            u.at_raw(ngx + nx, jj) = val;
            u.at_raw(ngx + nx + 1, jj) = val;
            break;
        }
        case BCType::Outflow: {
            double val = u.at_raw(ngx + nx - 1, jj);
            u.at_raw(ngx + nx, jj) = val;
            u.at_raw(ngx + nx + 1, jj) = val;
            break;
        }
        case BCType::Periodic:
            break;
        }
    }

    // Periodic left-right
    if (bc.left.type == BCType::Periodic && bc.right.type == BCType::Periodic) {
        for (int j = 0; j < ny; ++j) {
            int jj = j + ngy;
            double left_b = u.at_raw(ngx, jj);
            double left_i = u.at_raw(ngx + 1, jj);
            double right_b = u.at_raw(ngx + nx, jj);
            double right_i = u.at_raw(ngx + nx - 1, jj);
            u.at_raw(ngx - 1, jj) = right_i;
            u.at_raw(ngx, jj) = right_b;
            u.at_raw(ngx + nx, jj) = left_b;
            u.at_raw(ngx + nx + 1, jj) = left_i;
        }
    }

    int unx = g.u_nx();
    // Bottom/top boundaries (tangential)
    for (int i = 0; i < unx; ++i) {
        int ii = i + ngx;
        // Bottom j=0
        switch (bc.bottom.type) {
        case BCType::Wall:
            u.at_raw(ii, ngy) = 0.0;
            u.at_raw(ii, ngy - 1) = 0.0;
            break;
        case BCType::Moving: {
            double val = bc.bottom.moving;
            u.at_raw(ii, ngy) = val;
            u.at_raw(ii, ngy - 1) = val;
            break;
        }
        case BCType::Inflow: {
            double val = bc.bottom.inflow_u;
            u.at_raw(ii, ngy) = val;
            u.at_raw(ii, ngy - 1) = val;
            break;
        }
        case BCType::Outflow: {
            double val = u.at_raw(ii, ngy + 1);
            u.at_raw(ii, ngy) = val;
            u.at_raw(ii, ngy - 1) = val;
            break;
        }
        case BCType::Periodic:
            break;
        }

        // Top j=ny-1
        switch (bc.top.type) {
        case BCType::Wall:
            u.at_raw(ii, ngy + ny - 1) = 0.0;
            u.at_raw(ii, ngy + ny) = 0.0;
            break;
        case BCType::Moving: {
            double val = bc.top.moving;
            u.at_raw(ii, ngy + ny - 1) = val;
            u.at_raw(ii, ngy + ny) = val;
            break;
        }
        case BCType::Inflow: {
            double val = bc.top.inflow_u;
            u.at_raw(ii, ngy + ny - 1) = val;
            u.at_raw(ii, ngy + ny) = val;
            break;
        }
        case BCType::Outflow: {
            double val = u.at_raw(ii, ngy + ny - 2);
            u.at_raw(ii, ngy + ny - 1) = val;
            u.at_raw(ii, ngy + ny) = val;
            break;
        }
        case BCType::Periodic:
            break;
        }
    }

    // Periodic top-bottom
    if (bc.bottom.type == BCType::Periodic && bc.top.type == BCType::Periodic) {
        for (int i = 0; i < unx; ++i) {
            int ii = i + ngx;
            double bottom_b = u.at_raw(ii, ngy);
            double bottom_i = u.at_raw(ii, ngy + 1);
            double top_b = u.at_raw(ii, ngy + ny - 1);
            double top_i = u.at_raw(ii, ngy + ny - 2);
            u.at_raw(ii, ngy - 1) = top_i;
            u.at_raw(ii, ngy) = top_b;
            u.at_raw(ii, ngy + ny - 1) = bottom_b;
            u.at_raw(ii, ngy + ny) = bottom_i;
        }
    }
}

void apply_bc_v(const Grid &g, Field2D<double> &v, const BC &bc) {
    int nx = g.nx;
    int ny = g.ny;
    int ngx = g.ngx;
    int ngy = g.ngy;

    // Left/right boundaries (tangential)
    for (int j = 0; j < g.v_ny(); ++j) {
        int jj = j + ngy;
        // Left i=0
        switch (bc.left.type) {
        case BCType::Wall:
            v.at_raw(ngx, jj) = 0.0;
            v.at_raw(ngx - 1, jj) = 0.0;
            break;
        case BCType::Moving: {
            double val = bc.left.moving;
            v.at_raw(ngx, jj) = val;
            v.at_raw(ngx - 1, jj) = val;
            break;
        }
        case BCType::Inflow: {
            double y = j * g.dy;
            double val = bc.left.inflow_v * jet_mask(g, bc, y);
            v.at_raw(ngx, jj) = val;
            v.at_raw(ngx - 1, jj) = val;
            break;
        }
        case BCType::Outflow: {
            double val = v.at_raw(ngx + 1, jj);
            v.at_raw(ngx, jj) = val;
            v.at_raw(ngx - 1, jj) = val;
            break;
        }
        case BCType::Periodic:
            break;
        }

        // Right i=nx-1
        switch (bc.right.type) {
        case BCType::Wall:
            v.at_raw(ngx + nx - 1, jj) = 0.0;
            v.at_raw(ngx + nx, jj) = 0.0;
            break;
        case BCType::Moving: {
            double val = bc.right.moving;
            v.at_raw(ngx + nx - 1, jj) = val;
            v.at_raw(ngx + nx, jj) = val;
            break;
        }
        case BCType::Inflow: {
            double val = bc.right.inflow_v;
            v.at_raw(ngx + nx - 1, jj) = val;
            v.at_raw(ngx + nx, jj) = val;
            break;
        }
        case BCType::Outflow: {
            double val = v.at_raw(ngx + nx - 2, jj);
            v.at_raw(ngx + nx - 1, jj) = val;
            v.at_raw(ngx + nx, jj) = val;
            break;
        }
        case BCType::Periodic:
            break;
        }
    }

    // Periodic left-right for v
    if (bc.left.type == BCType::Periodic && bc.right.type == BCType::Periodic) {
        for (int j = 0; j < g.v_ny(); ++j) {
            int jj = j + ngy;
            double left_b = v.at_raw(ngx, jj);
            double left_i = v.at_raw(ngx + 1, jj);
            double right_b = v.at_raw(ngx + nx - 1, jj);
            double right_i = v.at_raw(ngx + nx - 2, jj);
            v.at_raw(ngx - 1, jj) = right_i;
            v.at_raw(ngx, jj) = right_b;
            v.at_raw(ngx + nx - 1, jj) = left_b;
            v.at_raw(ngx + nx, jj) = left_i;
        }
    }

    int vny = g.v_ny();
    // Bottom/top boundaries (normal component)
    for (int i = 0; i < nx; ++i) {
        int ii = i + ngx;
        // Bottom j=0
        switch (bc.bottom.type) {
        case BCType::Wall:
        case BCType::Moving:
            v.at_raw(ii, ngy) = 0.0;
            v.at_raw(ii, ngy - 1) = 0.0;
            break;
        case BCType::Inflow: {
            double val = bc.bottom.inflow_v;
            v.at_raw(ii, ngy) = val;
            v.at_raw(ii, ngy - 1) = val;
            break;
        }
        case BCType::Outflow: {
            double val = v.at_raw(ii, ngy + 1);
            v.at_raw(ii, ngy) = val;
            v.at_raw(ii, ngy - 1) = val;
            break;
        }
        case BCType::Periodic:
            break;
        }

        // Top j=ny
        switch (bc.top.type) {
        case BCType::Wall:
        case BCType::Moving:
            v.at_raw(ii, ngy + vny - 1) = 0.0;
            v.at_raw(ii, ngy + vny) = 0.0;
            break;
        case BCType::Inflow: {
            double val = bc.top.inflow_v;
            v.at_raw(ii, ngy + vny - 1) = val;
            v.at_raw(ii, ngy + vny) = val;
            break;
        }
        case BCType::Outflow: {
            double val = v.at_raw(ii, ngy + vny - 2);
            v.at_raw(ii, ngy + vny - 1) = val;
            v.at_raw(ii, ngy + vny) = val;
            break;
        }
        case BCType::Periodic:
            break;
        }
    }

    // Periodic top-bottom for v
    if (bc.bottom.type == BCType::Periodic && bc.top.type == BCType::Periodic) {
        for (int i = 0; i < nx; ++i) {
            int ii = i + ngx;
            double bottom_b = v.at_raw(ii, ngy);
            double bottom_i = v.at_raw(ii, ngy + 1);
            double top_b = v.at_raw(ii, ngy + vny - 1);
            double top_i = v.at_raw(ii, ngy + vny - 2);
            v.at_raw(ii, ngy - 1) = top_i;
            v.at_raw(ii, ngy) = top_b;
            v.at_raw(ii, ngy + vny - 1) = bottom_b;
            v.at_raw(ii, ngy + vny) = bottom_i;
        }
    }
}

void apply_bc_p(const Grid &g, Field2D<double> &p, const BC &bc) {
    int nx = g.nx;
    int ny = g.ny;
    int ngx = g.ngx;
    int ngy = g.ngy;

    // Left/right
    for (int j = 0; j < ny; ++j) {
        int jj = j + ngy;
        if (bc.left.type == BCType::Periodic &&
            bc.right.type == BCType::Periodic) {
            // handled later
        } else {
            p.at_raw(ngx - 1, jj) = p.at_raw(ngx, jj);
            p.at_raw(ngx + nx, jj) = p.at_raw(ngx + nx - 1, jj);
        }
    }
    if (bc.left.type == BCType::Periodic && bc.right.type == BCType::Periodic) {
        for (int j = 0; j < ny; ++j) {
            int jj = j + ngy;
            double left = p.at_raw(ngx, jj);
            double right = p.at_raw(ngx + nx - 1, jj);
            p.at_raw(ngx - 1, jj) = right;
            p.at_raw(ngx + nx, jj) = left;
        }
    }

    // Bottom/top
    for (int i = 0; i < nx; ++i) {
        int ii = i + ngx;
        if (bc.bottom.type == BCType::Periodic &&
            bc.top.type == BCType::Periodic) {
            // handled later
        } else {
            p.at_raw(ii, ngy - 1) = p.at_raw(ii, ngy);
            p.at_raw(ii, ngy + ny) = p.at_raw(ii, ngy + ny - 1);
        }
    }
    if (bc.bottom.type == BCType::Periodic && bc.top.type == BCType::Periodic) {
        for (int i = 0; i < nx; ++i) {
            int ii = i + ngx;
            double bottom = p.at_raw(ii, ngy);
            double top = p.at_raw(ii, ngy + ny - 1);
            p.at_raw(ii, ngy - 1) = top;
            p.at_raw(ii, ngy + ny) = bottom;
        }
    }
}

