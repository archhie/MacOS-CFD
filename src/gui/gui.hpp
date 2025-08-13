#pragma once

#include <GLFW/glfw3.h>
#include <OpenGL/gl3.h>
#include <vector>

#include "imgui.h"
#include "solver/bc.hpp"
#include "solver/grid.hpp"
#include "solver/state.hpp"

class Gui {
  public:
    enum class Field { U, V, Speed, Pressure, Vorticity };

    enum class Preset { JetPlume, LidDrivenCavity, PeriodicShear };

    double Re = 50.0;  // Updated to working value
    double CFL = 0.01;  // Updated to working value
    double dt = 0.001;  // dt override
    bool running = false;  // Changed to false - simulation starts paused
    bool step = false;
    bool reset = false;          // generic reset
    bool apply_preset = false;   // trigger to apply chosen preset
    bool reset_run = false;      // reset fields and start running
    Field field = Field::Speed;
    int preset = static_cast<int>(Preset::JetPlume);
    BC bc;             // boundary condition settings
    double Ly = 1.0;   // domain height for jet sliders
    bool inflow_inactive = false;
    double sim_speed = 1.0;      // Simulation speed multiplier

    bool init(GLFWwindow *window);
    void begin_frame();
    void draw(int timestep, double sim_time, double dt, double max_velocity,
              double div_l2, double pressure_residual, GLuint texture);
    void end_frame(GLFWwindow *window);
    void shutdown();
};

void update_field_texture(const Grid &g, const State &s, Gui::Field field,
                          GLuint tex, std::vector<unsigned char> &buffer);
