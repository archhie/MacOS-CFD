#pragma once

#include <GLFW/glfw3.h>
#include <OpenGL/gl3.h>
#include <vector>
#include <string>

#include "imgui.h"
#include "solver/bc.hpp"
#include "solver/grid.hpp"
#include "solver/state.hpp"

class Gui {
  public:
    enum class Field { 
        U, V, Speed, Pressure, Vorticity, 
        Streamlines, Temperature  // Added new fields
    };

    enum class Colormap { 
        Viridis, Plasma, Jet, Grayscale, Turbo, RdYlBu 
    };

    enum class Preset { JetPlume, LidDrivenCavity, PeriodicShear };

    // CFD Parameters
    double Re = 50.0;
    double CFL = 0.01;
    double dt = 0.001;
    bool running = false;
    bool step = false;
    bool reset = false;
    bool apply_preset = false;
    bool reset_run = false;
    double sim_speed = 1.0;

    // Visualization Parameters
    Field field = Field::Speed;
    Colormap colormap = Colormap::Viridis;
    bool show_vectors = false;
    bool normalize_field = false;
    float color_scale_min = 0.0f;
    float color_scale_max = 1.0f;
    bool auto_scale = true;
    
    // UI State
    int preset = static_cast<int>(Preset::JetPlume);
    BC bc;
    double Ly = 1.0;
    bool inflow_inactive = false;
    
    // Collapsible sections
    bool show_flow_params = true;
    bool show_timestep_control = true;
    bool show_sim_speed = true;
    bool show_sim_info = true;
    
    // Performance tracking
    double last_frame_time = 0.0;
    double fps = 0.0;
    double iterations_per_sec = 0.0;
    
    // Field statistics
    double field_min = 0.0;
    double field_max = 0.0;

    bool init(GLFWwindow *window);
    void begin_frame();
    void draw(int timestep, double sim_time, double dt, double max_velocity,
              double div_l2, double pressure_residual, GLuint texture);
    void end_frame(GLFWwindow *window);
    void shutdown();
    
    // Helper functions for improved UI
    void setup_modern_style();
    void draw_cfd_controls(int timestep, double sim_time, double dt, 
                          double max_velocity, double div_l2, double pressure_residual);
    void draw_visualization_panel(GLuint texture);
    void draw_boundary_conditions_panel();
    void draw_performance_stats();
    void draw_field_stats();
    void draw_preset_preview();
    
    // Utility functions
    std::string get_field_name(Field field) const;
    std::string get_colormap_name(Colormap cmap) const;
    void update_field_statistics(const Grid &g, const State &s);
};

void update_field_texture(const Grid &g, const State &s, Gui::Field field,
                          GLuint tex, std::vector<unsigned char> &buffer,
                          Gui::Colormap colormap = Gui::Colormap::Viridis,
                          float scale_min = 0.0f, float scale_max = 1.0f,
                          bool normalize = false);
