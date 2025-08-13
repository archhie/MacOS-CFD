#pragma once

#include <GLFW/glfw3.h>
#include <OpenGL/gl3.h>
#include <vector>

#include "imgui.h"
#include "solver/grid.hpp"
#include "solver/state.hpp"

class Gui {
  public:
    enum class Field { U, V, Speed, Pressure };

    double Re = 1000.0;
    double CFL = 0.5;
    double dt = 0.001;  // dt override
    bool running = true;
    bool step = false;
    bool reset = false;
    Field field = Field::Speed;

    bool init(GLFWwindow *window);
    void begin_frame();
    void draw(int timestep, double sim_time, double max_velocity,
              double pressure_residual, GLuint texture);
    void end_frame(GLFWwindow *window);
    void shutdown();
};

void update_field_texture(const Grid &g, const State &s, Gui::Field field,
                          GLuint tex, std::vector<unsigned char> &buffer);
