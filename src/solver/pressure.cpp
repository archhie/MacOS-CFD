#include "pressure.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

PressureSolver::PressureSolver(const Grid &grid) : g(grid) {
    r.allocate(g.nx, g.ny, g.p_pitch(), g.ngx, g.ngy);
    z.allocate(g.nx, g.ny, g.p_pitch(), g.ngx, g.ngy);
    s.allocate(g.nx, g.ny, g.p_pitch(), g.ngx, g.ngy);
    Ap.allocate(g.nx, g.ny, g.p_pitch(), g.ngx, g.ngy);
    size_t total = static_cast<size_t>(r.pitch) * (r.ny + 2 * r.ngy);
    std::memset(r.data, 0, total * sizeof(double));
    std::memset(z.data, 0, total * sizeof(double));
    std::memset(s.data, 0, total * sizeof(double));
    std::memset(Ap.data, 0, total * sizeof(double));
}

void PressureSolver::apply_laplacian(const Field2D<double> &in,
                                     Field2D<double> &out) const {
    double idx2 = 1.0 / (g.dx * g.dx);
    double idy2 = 1.0 / (g.dy * g.dy);

#pragma omp parallel for collapse(2)
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            double c = in.at_raw(ii, jj);
            double l = in.at_raw(ii - 1, jj);
            double rgt = in.at_raw(ii + 1, jj);
            double b = in.at_raw(ii, jj - 1);
            double t = in.at_raw(ii, jj + 1);
            if (!std::isfinite(c)) c = 0.0;
            if (!std::isfinite(l)) l = 0.0;
            if (!std::isfinite(rgt)) rgt = 0.0;
            if (!std::isfinite(b)) b = 0.0;
            if (!std::isfinite(t)) t = 0.0;
            out.at_raw(ii, jj) =
                (rgt - 2.0 * c + l) * idx2 + (t - 2.0 * c + b) * idy2;
        }
    }
}

double PressureSolver::dot(const Field2D<double> &a,
                           const Field2D<double> &b) const {
    double sum = 0.0;
#pragma omp parallel for collapse(2) reduction(+ : sum)
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            double va = a.at_raw(ii, jj);
            double vb = b.at_raw(ii, jj);
            if (!std::isfinite(va)) va = 0.0;
            if (!std::isfinite(vb)) vb = 0.0;
            sum += va * vb;
        }
    }
    return sum;
}

double PressureSolver::smooth(Field2D<double> &p, const Field2D<double> &rhs,
                              int sweeps) {
    double idx2 = 1.0 / (g.dx * g.dx);
    double idy2 = 1.0 / (g.dy * g.dy);
    double denom = 2.0 * (idx2 + idy2);

    for (int sIter = 0; sIter < sweeps; ++sIter) {
        for (int color = 0; color < 2; ++color) {
#pragma omp parallel for collapse(2)
            for (int j = 0; j < g.ny; ++j) {
                for (int i = 0; i < g.nx; ++i) {
                    if ((i + j) % 2 != color)
                        continue;
                    int ii = i + g.ngx;
                    int jj = j + g.ngy;
                    double rhs_val = rhs.at_raw(ii, jj);
                    double pl = p.at_raw(ii - 1, jj);
                    double pr = p.at_raw(ii + 1, jj);
                    double pb = p.at_raw(ii, jj - 1);
                    double pt = p.at_raw(ii, jj + 1);
                    if (!std::isfinite(rhs_val)) rhs_val = 0.0;
                    if (!std::isfinite(pl)) pl = 0.0;
                    if (!std::isfinite(pr)) pr = 0.0;
                    if (!std::isfinite(pb)) pb = 0.0;
                    if (!std::isfinite(pt)) pt = 0.0;
                    double sum = (pr + pl) * idx2 + (pt + pb) * idy2;
                    p.at_raw(ii, jj) = (rhs_val + sum) / denom;
                }
            }
        }
    }

    apply_laplacian(p, Ap);

#pragma omp parallel for collapse(2)
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            double rhs_val = rhs.at_raw(ii, jj);
            double Ap_val = Ap.at_raw(ii, jj);
            if (!std::isfinite(rhs_val)) rhs_val = 0.0;
            if (!std::isfinite(Ap_val)) Ap_val = 0.0;
            r.at_raw(ii, jj) = rhs_val - Ap_val;
        }
    }
    return std::sqrt(dot(r, r));
}

double PressureSolver::solve(Field2D<double> &p, const Field2D<double> &rhs,
                             const PressureParams &params) {
    if (params.type == PressureSolverType::Multigrid) {
        // Multigrid solver
        double res = 1e9;
        for (int vcycle = 0; vcycle < params.vcycles; ++vcycle) {
            res = smooth(p, rhs, 2);
            if (res < params.tol)
                break;
        }
        last_residual_ = std::isfinite(res) ? res : 1e9;
        if (!std::isfinite(res) && !warned_nan_) {
            std::fprintf(stderr, "Non-finite pressure residual\n");
            warned_nan_ = true;
        }
        return last_residual_;
    }

    // PCG solver with proper Neumann BCs and pressure pinning
    // Clear pressure field first
    size_t total = static_cast<size_t>(p.pitch) * (p.ny + 2 * p.ngy);
    std::memset(p.data, 0, total * sizeof(double));
    
    // Pin pressure at (0,0) to kill nullspace
    p.at_raw(g.ngx, g.ngy) = 0.0;
    
    // Apply Neumann boundary conditions: ∂p/∂n = 0
    // This is handled implicitly by the Laplacian operator
    // For walls: ∂p/∂n = 0 means no pressure gradient normal to wall
    // For inflow/outflow: ∂p/∂n = 0 is a reasonable approximation
    
    // Initialize residual: r = rhs - A*p
    apply_laplacian(p, Ap);
#pragma omp parallel for collapse(2)
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            double rhs_val = rhs.at_raw(ii, jj);
            double Ap_val = Ap.at_raw(ii, jj);
            if (!std::isfinite(rhs_val)) rhs_val = 0.0;
            if (!std::isfinite(Ap_val)) Ap_val = 0.0;
            r.at_raw(ii, jj) = rhs_val - Ap_val;
        }
    }

    // Preconditioner: diagonal scaling
    double idx2 = 1.0 / (g.dx * g.dx);
    double idy2 = 1.0 / (g.dy * g.dy);
    double diag = 2.0 * (idx2 + idy2);
    double inv_diag = 1.0 / std::max(diag, 1e-12);

#pragma omp parallel for collapse(2)
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            z.at_raw(ii, jj) = inv_diag * r.at_raw(ii, jj);
            s.at_raw(ii, jj) = z.at_raw(ii, jj);
        }
    }

    double rho = dot(r, z);
    double res_norm = std::sqrt(std::max(0.0, dot(r, r)));
    double initial_res = res_norm;

    for (int iter = 0; iter < params.pcg_max_iters && res_norm > params.tol;
         ++iter) {
        apply_laplacian(s, Ap);
        double denom = dot(s, Ap);
        if (!std::isfinite(denom) || std::abs(denom) < 1e-20) {
            std::fprintf(stderr, "PCG breakdown at iteration %d\n", iter);
            break;
        }
        double alpha = rho / denom;
        
#pragma omp parallel for collapse(2)
        for (int j = 0; j < g.ny; ++j) {
            for (int i = 0; i < g.nx; ++i) {
                int ii = i + g.ngx;
                int jj = j + g.ngy;
                p.at_raw(ii, jj) += alpha * s.at_raw(ii, jj);
                r.at_raw(ii, jj) -= alpha * Ap.at_raw(ii, jj);
            }
        }

        // Re-pin pressure at (0,0) after each iteration to maintain nullspace elimination
        p.at_raw(g.ngx, g.ngy) = 0.0;

        res_norm = std::sqrt(std::max(0.0, dot(r, r)));
        if (!std::isfinite(res_norm) || res_norm < params.tol)
            break;

        // Precondition residual
#pragma omp parallel for collapse(2)
        for (int j = 0; j < g.ny; ++j) {
            for (int i = 0; i < g.nx; ++i) {
                int ii = i + g.ngx;
                int jj = j + g.ngy;
                z.at_raw(ii, jj) = inv_diag * r.at_raw(ii, jj);
            }
        }
        
        double rho_new = dot(r, z);
        if (!std::isfinite(rho_new))
            break;
        double beta = rho_new / rho;
        
#pragma omp parallel for collapse(2)
        for (int j = 0; j < g.ny; ++j) {
            for (int i = 0; i < g.nx; ++i) {
                int ii = i + g.ngx;
                int jj = j + g.ngy;
                s.at_raw(ii, jj) = z.at_raw(ii, jj) + beta * s.at_raw(ii, jj);
            }
        }
        rho = rho_new;
    }
    
    last_residual_ = std::isfinite(res_norm) ? res_norm : 1e9;
    if (!std::isfinite(res_norm) && !warned_nan_) {
        std::fprintf(stderr, "Non-finite pressure residual\n");
        warned_nan_ = true;
    }
    
    // Final pressure pinning
    p.at_raw(g.ngx, g.ngy) = 0.0;
    
    return last_residual_;
}

