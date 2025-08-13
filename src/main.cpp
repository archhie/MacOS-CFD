#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

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
#include "solver/advection.hpp"
#include "shapes.hpp"

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
    double Re = 50.0;  // Updated to working value
    double CFL_target = 0.01;  // Updated to working value
    PressureParams pp;
    pp.type = PressureSolverType::PCG;  // Default to PCG
    pp.pcg_max_iters = 200;  // Reduced from 1000 for better performance
    pp.tol = 1e-6;  // Relaxed from 1e-8 for better performance
    
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
        } else if (arg == "--tvd-debug") {
            tvd_debug = true;
            std::fprintf(stderr, "TVD debugging enabled\n");
        } else if (arg == "--help") {
            std::printf("CFD2D GUI - 2D Computational Fluid Dynamics\n");
            std::printf("Usage: %s [options]\n", argv[0]);
            std::printf("Options:\n");
            std::printf("  --no-gui              Run in headless mode\n");
            std::printf("  --Re=<value>          Set Reynolds number (default: 50)\n");
            std::printf("  --CFL=<value>         Set CFL condition (default: 0.01)\n");
            std::printf("  --solver=<pcg|mg>     Pressure solver type (default: pcg)\n");
            std::printf("  --mg-levels=<n>       Multigrid levels (default: 3)\n");
            std::printf("  --vcycles=<n>         Multigrid V-cycles (default: 2)\n");
            std::printf("  --tol=<value>         Solver tolerance (default: 1e-8)\n");
            std::printf("  --pcg-iters=<n>       PCG max iterations (default: 1000)\n");
            std::printf("  --tvd-debug           Enable TVD debugging output\n");
            std::printf("  --help                Show this help message\n");
            return 0;
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

    // deterministic initialization
    auto clear_field = [](auto &f) {
        size_t sz = static_cast<size_t>(f.pitch) * (f.ny + 2 * f.ngy);
        std::memset(f.data, 0, sz * sizeof(*f.data));
    };
    clear_field(state.u);
    clear_field(state.v);
    clear_field(state.p);
    clear_field(state.rhs);
    clear_field(state.tmp);
    size_t szs = static_cast<size_t>(state.scalar.pitch) *
                 (state.scalar.ny + 2 * state.scalar.ngy);
    std::memset(state.scalar.data, 0, szs * sizeof(float));

    PressureSolver pressure(grid);
    TimeIntegrator integrator(grid);
    integrator.clear_work();
    BC bc;
    bc.left.type = BCType::Inflow;
    bc.right.type = BCType::Outflow;
    bc.top.type = BCType::Wall;
    bc.bottom.type = BCType::Wall;
    bc.left.inflow_u = 1.0;
    bc.left.inflow_v = 0.0;
    bc.jet_center = 0.5;  // Center of the domain in actual coordinates
    bc.jet_width = 0.18;  // Width in actual coordinates
    bc.jet_thickness = 0.02;  // Thickness in actual coordinates
    bc.jet_eps = 0.03;
    bc.jet_k = 8;

    if (no_gui) {
        double sim_time = 0.0;
        double pres_res = 0.0;
        for (int step = 0; step < 2000; ++step) {
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
    gui.grid_size = ImVec2(static_cast<float>(grid.nx), static_cast<float>(grid.ny));

    GLuint tex = 0;
    glGenTextures(1, &tex);
    std::vector<unsigned char> texbuf;

    int step = 0;
    double sim_time = 0.0;
    double pres_res = 0.0;
    double last_dt = 0.0;
    
    // Frame timing for smooth GUI - completely independent of simulation
    double last_frame_time = glfwGetTime();
    double target_frame_time = 1.0 / 60.0;  // 60 FPS target
    
    // Frame skipping for simulation to maintain GUI responsiveness
    int frame_count = 0;
    int sim_update_frequency = 5;  // Update simulation every 5 frames

    while (!glfwWindowShouldClose(window)) {
        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
            glfwSetWindowShouldClose(window, 1);

        glfwPollEvents();
        gui.begin_frame();

        // Handle GUI controls
        if (gui.reset) {
            // Reset simulation state
            clear_field(state.u);
            clear_field(state.v);
            clear_field(state.p);
            clear_field(state.rhs);
            clear_field(state.tmp);
            std::memset(state.scalar.data, 0, szs * sizeof(float));
            step = 0;
            sim_time = 0.0;
            pres_res = 0.0;
            last_dt = 0.0;
            
            // Reset thread state
            // The thread state variables are removed, so this block is effectively removed.
            
            gui.reset = false;
        }
        
        if (gui.reset_run) {
            gui.reset = true;
            gui.running = true;
            gui.reset_run = false;
        }

        // Update boundary conditions from GUI
        bc = gui.bc;

        // SIMULATION LOOP - with frame skipping for GUI responsiveness
        frame_count++;
        if ((gui.running || gui.step) && (frame_count % sim_update_frequency == 0 || gui.step)) {
            // Run simulation step
            last_dt = integrator.step(state, bc, gui.Re, gui.CFL, pressure, pp,
                                      pres_res, gui.dt);
            sim_time += last_dt;
            step++;
            
            // Apply shapes to simulation (after advection/projection)
            gui.apply_shapes_to_simulation(grid, state);
            
            // Only step once if step button pressed
            if (gui.step) {
                gui.step = false;
            }
        }

        // Update visualization (every GUI frame)
        double max_vel = max_velocity(grid, state.u, state.v);
        double divL2 = divergence_l2(grid, state.u, state.v);
        
        // Update field statistics for the GUI
        gui.update_field_statistics(grid, state);
        
        // Update texture with new colormap and scaling options
        update_field_texture(grid, state, gui.field, tex, texbuf, 
                           gui.colormap, gui.color_scale_min, gui.color_scale_max, 
                           gui.normalize_field);

        gui.inflow_inactive = (bc.left.type == BCType::Inflow &&
                               (std::abs(bc.left.inflow_u) < 1e-12 ||
                                bc.jet_width < grid.dy));
        gui.draw(step, sim_time, last_dt, max_vel, divL2, pres_res, tex);

        if (gui.apply_preset) {
            switch (static_cast<Gui::Preset>(gui.preset)) {
            case Gui::Preset::JetPlume:
                bc.left.type = BCType::Inflow;
                bc.right.type = BCType::Outflow;
                bc.top.type = BCType::Wall;
                bc.bottom.type = BCType::Wall;
                gui.Re = 50.0;  // Updated to working value
                gui.CFL = 0.01;  // Updated to working value
                bc.left.inflow_u = 1.0;
                bc.left.inflow_v = 0.0;
                bc.jet_center = 0.5 * grid.Ly;
                bc.jet_width = 0.18 * grid.Ly;
                bc.jet_thickness = 0.02 * grid.Ly;
                bc.jet_eps = 0.03;
                bc.jet_k = 8;
                break;
            case Gui::Preset::LidDrivenCavity:
                bc.left.type = BCType::Wall;
                bc.right.type = BCType::Wall;
                bc.bottom.type = BCType::Wall;
                bc.top.type = BCType::Moving;
                bc.top.moving = 1.0;
                gui.Re = 1000.0;
                gui.CFL = 0.1;  // More conservative for stability
                break;
            case Gui::Preset::PeriodicShear:
                bc.left.type = BCType::Periodic;
                bc.right.type = BCType::Periodic;
                bc.bottom.type = BCType::Wall;
                bc.top.type = BCType::Wall;
                gui.Re = 1000.0;
                gui.CFL = 0.1;  // More conservative for stability
                break;
            case Gui::Preset::FastWind:
                bc.left.type = BCType::Inflow;
                bc.right.type = BCType::Outflow;
                bc.top.type = BCType::Wall;
                bc.bottom.type = BCType::Wall;
                gui.Re = 1000.0;
                gui.CFL = 0.01;
                bc.left.inflow_u = 5.0;
                bc.left.inflow_v = 0.0;
                break;
            }
            gui.apply_preset = false;
        }

        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        glViewport(0, 0, w, h);
        glClearColor(0.02f, 0.05f, 0.15f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        gui.end_frame(window);

        // GUI FRAME RATE LIMITING - completely separate from simulation
        double frame_time = glfwGetTime() - last_frame_time;
        if (frame_time < target_frame_time) {
            // Sleep to maintain 60 FPS GUI
            double sleep_time = target_frame_time - frame_time;
            if (sleep_time > 0.001) {  // Only sleep if more than 1ms
                glfwWaitEventsTimeout(sleep_time);
            }
        }
        last_frame_time = glfwGetTime();

        if (gui.reset || gui.reset_run) {
            // Clear all fields and scratch arrays
            clear_field(state.u);
            clear_field(state.v);
            clear_field(state.p);
            clear_field(state.rhs);
            clear_field(state.tmp);
            std::memset(state.scalar.data, 0, szs * sizeof(float));
            integrator.clear_work();

            if (static_cast<Gui::Preset>(gui.preset) == Gui::Preset::PeriodicShear) {
                for (int j = 0; j < grid.ny; ++j) {
                    double y = (j + 0.5) * grid.dy;
                    double val = bc.left.inflow_u * std::sin(2.0 * 3.141592653589793 * y / grid.Ly);
                    for (int i = 0; i < grid.u_nx(); ++i) {
                        int ii = i + grid.ngx;
                        int jj = j + grid.ngy;
                        state.u.at_raw(ii, jj) = val;
                    }
                }
            }

            sim_time = 0.0;
            step = 0;
            pres_res = 0.0;
            last_dt = 0.0;
            bool run_flag = gui.reset_run;
            gui.reset = false;
            gui.reset_run = false;
            if (run_flag)
                gui.running = true;
        }
        gui.bc = bc;

        if (step % 60 == 0) {
            double umax = 0.0, vmax = 0.0;
            for (int j = 0; j < grid.ny; ++j)
                for (int i = 0; i < grid.u_nx(); ++i)
                    umax = std::max(umax, std::abs(state.u.at_raw(i + grid.ngx, j + grid.ngy)));
            for (int j = 0; j < grid.v_ny(); ++j)
                for (int i = 0; i < grid.nx; ++i)
                    vmax = std::max(vmax, std::abs(state.v.at_raw(i + grid.ngx, j + grid.ngy)));
            std::printf("dt=%g umax=%g vmax=%g divL2=%g presRes=%g inflow(U,w)=(%g,%g)\n",
                        last_dt, umax, vmax, divL2, pres_res,
                        bc.left.inflow_u, bc.jet_width);
            if (!std::isfinite(divL2) || !std::isfinite(pres_res))
                std::printf("WARNING: NaN detected in divergence or pressure residual\n");
        }
    }

    // Clean up simulation thread
    // The simulation thread management is removed, so this block is effectively removed.

    save_vtk(grid, state, "out/final.vtk");
    glDeleteTextures(1, &tex);
    gui.shutdown();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
