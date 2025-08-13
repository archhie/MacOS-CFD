#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>
#include <OpenGL/gl3.h>

#include "gui/gui.hpp"
#include "solver/bc.hpp"
#include "solver/grid.hpp"
#include "solver/metrics.hpp"
#include "solver/pressure.hpp"
#include "solver/state.hpp"
#include "solver/time_integrator.hpp"

static void save_vtk(const Grid &g, const State &s, const std::string &path) {
    std::filesystem::create_directories(std::filesystem::path(path).parent_path());
    std::ofstream out(path);
    out << "# vtk DataFile Version 3.0\nCFD2D snapshot\nASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << g.nx << " " << g.ny << " 1\n";
    out << "ORIGIN 0 0 0\n";
    out << "SPACING " << g.dx << " " << g.dy << " 1\n";
    out << "POINT_DATA " << g.nx * g.ny << "\n";
    out << "SCALARS pressure double\nLOOKUP_TABLE default\n";
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            out << s.p.at_raw(ii, jj) << "\n";
        }
    }
    out << "VECTORS velocity double\n";
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            double uc = 0.5 * (s.u.at_raw(ii, jj) + s.u.at_raw(ii + 1, jj));
            double vc = 0.5 * (s.v.at_raw(ii, jj) + s.v.at_raw(ii, jj + 1));
            out << uc << ' ' << vc << " 0\n";
        }
    }
}

int main(int argc, char **argv) {
    bool no_gui = false;
    double Re = 1000.0;
    double CFL_target = 0.5;
    PressureParams pp;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--no-gui") {
            no_gui = true;
        } else if (arg.rfind("--Re=", 0) == 0) {
            Re = std::stod(arg.substr(5));
        } else if (arg.rfind("--CFL=", 0) == 0) {
            CFL_target = std::stod(arg.substr(6));
        } else if (arg.rfind("--solver=", 0) == 0) {
            std::string val = arg.substr(9);
            if (val == "pcg")
                pp.type = PressureSolverType::PCG;
            else
                pp.type = PressureSolverType::Multigrid;
        } else if (arg.rfind("--mg-levels=", 0) == 0) {
            pp.mg_levels = std::stoi(arg.substr(12));
        } else if (arg.rfind("--vcycles=", 0) == 0) {
            pp.vcycles = std::stoi(arg.substr(10));
        } else if (arg.rfind("--tol=", 0) == 0) {
            pp.tol = std::stod(arg.substr(6));
        } else if (arg.rfind("--pcg-iters=", 0) == 0) {
            pp.pcg_max_iters = std::stoi(arg.substr(12));
        }
    }

    Grid grid;
    grid.init(512, 256, 2.0, 1.0, 1);
    State state;
    state.u.allocate(grid.u_nx(), grid.ny, grid.u_pitch(), grid.ngx, grid.ngy);
    state.v.allocate(grid.nx, grid.v_ny(), grid.v_pitch(), grid.ngx, grid.ngy);
    state.p.allocate(grid.nx, grid.ny, grid.p_pitch(), grid.ngx, grid.ngy);
    state.rhs.allocate(grid.nx, grid.ny, grid.p_pitch(), grid.ngx, grid.ngy);
    state.tmp.allocate(grid.nx, grid.ny, grid.p_pitch(), grid.ngx, grid.ngy);
    state.scalar.allocate(grid.nx, grid.ny, grid.p_pitch(), grid.ngx, grid.ngy);

    PressureSolver pressure(grid);
    TimeIntegrator integrator(grid);
    BC bc;
    bc.left.type = BCType::Inflow;
    bc.right.type = BCType::Outflow;
    bc.top.type = BCType::Wall;
    bc.bottom.type = BCType::Wall;
    bc.left.inflow_u = 1.0;
    bc.left.inflow_v = 0.0;
    bc.jet_center = 0.5;
    bc.jet_width = 0.2;

    if (no_gui) {
        double sim_time = 0.0;
        double pres_res = 0.0;
        for (int step = 0; step < 100; ++step) {
            double dt = integrator.step(state, bc, Re, CFL_target, pressure, pp,
                                        pres_res);
            sim_time += dt;
            double divL2 = divergence_l2(grid, state.u, state.v);
            if (step % 10 == 0) {
                std::printf("step %d time %g dt %g div %g res %g\n", step,
                            sim_time, dt, divL2, pres_res);
            }
        }
        save_vtk(grid, state, "out/final.vtk");
        return 0;
    }

    if (!glfwInit())
        return 1;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GLFW_TRUE);

    GLFWwindow *window =
        glfwCreateWindow(1280, 720, "CFD2D", nullptr, nullptr);
    if (!window) {
        glfwTerminate();
        return 1;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    Gui gui;
    gui.init(window);
    gui.bc = bc;
    gui.Ly = grid.Ly;

    GLuint tex = 0;
    glGenTextures(1, &tex);
    std::vector<unsigned char> texbuf;

    int step = 0;
    double sim_time = 0.0;
    double pres_res = 0.0;
    double last_dt = 0.0;

    while (!glfwWindowShouldClose(window)) {
        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
            glfwSetWindowShouldClose(window, 1);

        glfwPollEvents();
        gui.begin_frame();

        if (gui.running || gui.step) {
            bc = gui.bc;
            last_dt = integrator.step(state, bc, gui.Re, gui.CFL, pressure, pp,
                                      pres_res, gui.dt);
            sim_time += last_dt;
            step++;
            gui.step = false;
        }

        double max_vel = max_velocity(grid, state.u, state.v);
        double divL2 = divergence_l2(grid, state.u, state.v);
        update_field_texture(grid, state, gui.field, tex, texbuf);

        gui.draw(step, sim_time, last_dt, max_vel, divL2, pres_res, tex);

        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        glViewport(0, 0, w, h);
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        gui.end_frame(window);

        if (gui.reset) {
            size_t sz;
            sz = static_cast<size_t>(state.u.pitch) *
                 (state.u.ny + 2 * state.u.ngy);
            std::memset(state.u.data, 0, sz * sizeof(double));
            sz = static_cast<size_t>(state.v.pitch) *
                 (state.v.ny + 2 * state.v.ngy);
            std::memset(state.v.data, 0, sz * sizeof(double));
            sz = static_cast<size_t>(state.p.pitch) *
                 (state.p.ny + 2 * state.p.ngy);
            std::memset(state.p.data, 0, sz * sizeof(double));
            sim_time = 0.0;
            step = 0;
            gui.reset = false;
        }
        gui.bc = bc;
    }

    save_vtk(grid, state, "out/final.vtk");
    glDeleteTextures(1, &tex);
    gui.shutdown();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
