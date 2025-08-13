#include "advection.hpp"

#include <algorithm>
#include <cmath>

static inline double minmod(double a, double b) {
    return (a * b <= 0.0) ? 0.0 : (std::abs(a) < std::abs(b) ? a : b);
}

void advect_u(const Grid& g, const Field2D<double>& u, const Field2D<double>& v, Field2D<double>& du_dt) {
    double idx = 1.0 / g.dx;
    double idy = 1.0 / g.dy;
    int ngx = g.ngx;
    int ngy = g.ngy;

#pragma omp parallel for collapse(2)
    for (int j = 1; j < g.ny - 1; ++j) {
        for (int i = 1; i < g.u_nx() - 1; ++i) {
            int ii = i + ngx;
            int jj = j + ngy;

            double uim2 = u.at_raw(ii - 2, jj);
            double uim1 = u.at_raw(ii - 1, jj);
            double ui = u.at_raw(ii, jj);
            double uip1 = u.at_raw(ii + 1, jj);
            double uip2 = u.at_raw(ii + 2, jj);

            double s_i = minmod(ui - uim1, uip1 - ui);
            double s_ip1 = minmod(uip1 - ui, uip2 - uip1);
            double s_im1 = minmod(uim1 - uim2, ui - uim1);

            double uL = ui + 0.5 * s_i;
            double uR = uip1 - 0.5 * s_ip1;
            double u_face_p = 0.5 * (ui + uip1);
            double flux_x_p = (u_face_p > 0.0 ? uL : uR) * u_face_p;

            double uL_m = uim1 + 0.5 * s_im1;
            double uR_m = ui - 0.5 * s_i;
            double u_face_m = 0.5 * (uim1 + ui);
            double flux_x_m = (u_face_m > 0.0 ? uL_m : uR_m) * u_face_m;

            double ujm2 = u.at_raw(ii, jj - 2);
            double ujm1 = u.at_raw(ii, jj - 1);
            double uj = ui;
            double ujp1 = u.at_raw(ii, jj + 1);
            double ujp2 = u.at_raw(ii, jj + 2);

            double sy_j = minmod(uj - ujm1, ujp1 - uj);
            double sy_jp1 = minmod(ujp1 - uj, ujp2 - ujp1);
            double sy_jm1 = minmod(ujm1 - ujm2, uj - ujm1);

            double uLy = uj + 0.5 * sy_j;
            double uRy = ujp1 - 0.5 * sy_jp1;
            double v_face_p = 0.5 * (v.at_raw(i + ngx, j + 1 + ngy) +
                                      v.at_raw(i - 1 + ngx, j + 1 + ngy));
            double flux_y_p = (v_face_p > 0.0 ? uLy : uRy) * v_face_p;

            double uLy_m = ujm1 + 0.5 * sy_jm1;
            double uRy_m = uj - 0.5 * sy_j;
            double v_face_m = 0.5 * (v.at_raw(i + ngx, j + ngy) +
                                      v.at_raw(i - 1 + ngx, j + ngy));
            double flux_y_m = (v_face_m > 0.0 ? uLy_m : uRy_m) * v_face_m;

            double conv = (flux_x_p - flux_x_m) * idx + (flux_y_p - flux_y_m) * idy;
            du_dt.at_raw(ii, jj) = -conv;
        }
    }
}

void advect_v(const Grid& g, const Field2D<double>& u, const Field2D<double>& v, Field2D<double>& dv_dt) {
    double idx = 1.0 / g.dx;
    double idy = 1.0 / g.dy;
    int ngx = g.ngx;
    int ngy = g.ngy;

#pragma omp parallel for collapse(2)
    for (int j = 1; j < g.v_ny() - 1; ++j) {
        for (int i = 1; i < g.nx - 1; ++i) {
            int ii = i + ngx;
            int jj = j + ngy;

            double vjm2 = v.at_raw(ii, jj - 2);
            double vjm1 = v.at_raw(ii, jj - 1);
            double vj = v.at_raw(ii, jj);
            double vjp1 = v.at_raw(ii, jj + 1);
            double vjp2 = v.at_raw(ii, jj + 2);

            double sy_j = minmod(vj - vjm1, vjp1 - vj);
            double sy_jp1 = minmod(vjp1 - vj, vjp2 - vjp1);
            double sy_jm1 = minmod(vjm1 - vjm2, vj - vjm1);

            double vL = vj + 0.5 * sy_j;
            double vR = vjp1 - 0.5 * sy_jp1;
            double v_face_p = 0.5 * (vj + vjp1);
            double flux_y_p = (v_face_p > 0.0 ? vL : vR) * v_face_p;

            double vL_m = vjm1 + 0.5 * sy_jm1;
            double vR_m = vj - 0.5 * sy_j;
            double v_face_m = 0.5 * (vjm1 + vj);
            double flux_y_m = (v_face_m > 0.0 ? vL_m : vR_m) * v_face_m;

            double vim2 = v.at_raw(ii - 2, jj);
            double vim1 = v.at_raw(ii - 1, jj);
            double vi = vj;
            double vip1 = v.at_raw(ii + 1, jj);
            double vip2 = v.at_raw(ii + 2, jj);

            double sx_i = minmod(vi - vim1, vip1 - vi);
            double sx_ip1 = minmod(vip1 - vi, vip2 - vip1);
            double sx_im1 = minmod(vim1 - vim2, vi - vim1);

            double vHx = vi + 0.5 * sx_i;
            double vHxR = vip1 - 0.5 * sx_ip1;
            double u_face_p = 0.5 * (u.at_raw(i + 1 + ngx, j + ngy) +
                                      u.at_raw(i + 1 + ngx, j - 1 + ngy));
            double flux_x_p = (u_face_p > 0.0 ? vHx : vHxR) * u_face_p;

            double vHxL = vim1 + 0.5 * sx_im1;
            double vHxC = vi - 0.5 * sx_i;
            double u_face_m = 0.5 * (u.at_raw(i + ngx, j + ngy) +
                                      u.at_raw(i + ngx, j - 1 + ngy));
            double flux_x_m = (u_face_m > 0.0 ? vHxL : vHxC) * u_face_m;

            double conv = (flux_x_p - flux_x_m) * idx + (flux_y_p - flux_y_m) * idy;
            dv_dt.at_raw(ii, jj) = -conv;
        }
    }
}
