#pragma once

struct Grid {
    int nx = 0, ny = 0;    // cells in x/y
    int ngx = 1, ngy = 1;  // ghost layers
    double Lx = 1.0, Ly = 1.0;  // domain size
    double dx = 0.0, dy = 0.0;  // cell size

    void init(int Nx, int Ny, double Lx_, double Ly_, int ng = 1);

    // pressure (cell centers)
    inline int p_pitch() const { return nx + 2 * ngx; }
    inline int p_idx(int i, int j) const { return (j + ngy) * p_pitch() + (i + ngx); }

    // u faces (x-directed)
    inline int u_nx() const { return nx + 1; }
    inline int u_pitch() const { return u_nx() + 2 * ngx; }
    inline int u_idx(int i, int j) const { return (j + ngy) * u_pitch() + (i + ngx); }

    // v faces (y-directed)
    inline int v_ny() const { return ny + 1; }
    inline int v_pitch() const { return nx + 2 * ngx; }
    inline int v_idx(int i, int j) const { return (j + ngy) * v_pitch() + (i + ngx); }
};

