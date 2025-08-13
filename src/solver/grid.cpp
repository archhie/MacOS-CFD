#include "grid.hpp"

void Grid::init(int Nx, int Ny, double Lx_, double Ly_, int ng) {
    nx = Nx;
    ny = Ny;
    Lx = Lx_;
    Ly = Ly_;
    ngx = ng;
    ngy = ng;
    dx = Lx / static_cast<double>(nx);
    dy = Ly / static_cast<double>(ny);
}

