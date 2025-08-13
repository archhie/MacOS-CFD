#include "pressure.hpp"

#include <cmath>
#include <cstring>

PressureSolver::PressureSolver(const Grid &grid) : g(grid) {
    r.allocate(g.nx, g.ny, g.p_pitch(), g.ngx, g.ngy);
    z.allocate(g.nx, g.ny, g.p_pitch(), g.ngx, g.ngy);
    s.allocate(g.nx, g.ny, g.p_pitch(), g.ngx, g.ngy);
    Ap.allocate(g.nx, g.ny, g.p_pitch(), g.ngx, g.ngy);
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
            out.at_raw(ii, jj) =
                (in.at_raw(ii + 1, jj) - 2.0 * in.at_raw(ii, jj) +
                 in.at_raw(ii - 1, jj)) * idx2 +
                (in.at_raw(ii, jj + 1) - 2.0 * in.at_raw(ii, jj) +
                 in.at_raw(ii, jj - 1)) * idy2;
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
            sum += a.at_raw(ii, jj) * b.at_raw(ii, jj);
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
                    double sum = (p.at_raw(ii + 1, jj) + p.at_raw(ii - 1, jj)) * idx2 +
                                 (p.at_raw(ii, jj + 1) + p.at_raw(ii, jj - 1)) * idy2;
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
            r.at_raw(ii, jj) = rhs.at_raw(ii, jj) - Ap.at_raw(ii, jj);
        }
    }
    return std::sqrt(dot(r, r));
}

double PressureSolver::solve(Field2D<double> &p, const Field2D<double> &rhs,
                             const PressureParams &params) {
    if (params.type == PressureSolverType::Multigrid) {
        double res = 0.0;
        for (int cycle = 0; cycle < params.vcycles; ++cycle) {
            res = smooth(p, rhs, params.rbgs_sweeps);
            if (res < params.tol)
                break;
        }
        return res;
    }

    // PCG solver
    size_t total = static_cast<size_t>(r.pitch) * (r.ny + 2 * r.ngy);
    std::memset(r.data, 0, total * sizeof(double));
    std::memset(z.data, 0, total * sizeof(double));
    std::memset(s.data, 0, total * sizeof(double));
    std::memset(Ap.data, 0, total * sizeof(double));

    apply_laplacian(p, Ap);
#pragma omp parallel for collapse(2)
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            r.at_raw(ii, jj) = rhs.at_raw(ii, jj) - Ap.at_raw(ii, jj);
        }
    }

    double idx2 = 1.0 / (g.dx * g.dx);
    double idy2 = 1.0 / (g.dy * g.dy);
    double inv_diag = -1.0 / (2.0 * (idx2 + idy2));

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
    double res_norm = std::sqrt(dot(r, r));

    for (int iter = 0; iter < params.pcg_max_iters && res_norm > params.tol;
         ++iter) {
        apply_laplacian(s, Ap);
        double alpha = rho / dot(s, Ap);
#pragma omp parallel for collapse(2)
        for (int j = 0; j < g.ny; ++j) {
            for (int i = 0; i < g.nx; ++i) {
                int ii = i + g.ngx;
                int jj = j + g.ngy;
                p.at_raw(ii, jj) += alpha * s.at_raw(ii, jj);
                r.at_raw(ii, jj) -= alpha * Ap.at_raw(ii, jj);
            }
        }

        res_norm = std::sqrt(dot(r, r));
        if (res_norm < params.tol)
            break;

#pragma omp parallel for collapse(2)
        for (int j = 0; j < g.ny; ++j) {
            for (int i = 0; i < g.nx; ++i) {
                int ii = i + g.ngx;
                int jj = j + g.ngy;
                z.at_raw(ii, jj) = inv_diag * r.at_raw(ii, jj);
            }
        }
        double rho_new = dot(r, z);
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
    return res_norm;
}

