#pragma once

#include <GLFW/glfw3.h>
#include <OpenGL/gl3.h>
#include <vector>
#include <string>

#include "imgui.h"
#include "solver/bc.hpp"
#include "solver/grid.hpp"
#include "solver/state.hpp"

// Shape structure for object drawer
struct Shape {
    int type; // 0 = rectangle, 1 = circle
    ImVec2 pos;
    ImVec2 size;
    int bcType; // 0 = wall, 1 = porous
    bool visible;
    bool enabled;
    
    Shape() : type(0), pos(0.0f, 0.0f), size(10.0f, 10.0f), bcType(0), visible(true), enabled(true) {}
};

class Gui {
  public:
    enum class Field { 
        U, V, Speed, Pressure, Vorticity, 
        Streamlines, Temperature  // Added new fields
    };

    enum class Colormap { 
        Viridis, Plasma, Jet, Grayscale, Turbo, RdYlBu, RdBu, Coolwarm
    };

    enum class Preset { JetPlume, LidDrivenCavity, PeriodicShear, FastWind };

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
    bool show_contours = false;
    bool show_streamlines = false;
    bool normalize_field = false;
    float color_scale_min = 0.0f;
    float color_scale_max = 1.0f;
    bool auto_scale = true;
    float gamma = 1.0f;
    
    // UI State
    int preset = static_cast<int>(Preset::JetPlume);
    BC bc;
    double Ly = 1.0;
    bool inflow_inactive = false;
    
    // Dockable Interface State
    bool show_controls = true;
    bool show_viz = true;
    bool show_bc = true;
    bool show_object_drawer = true;
    bool layout_initialized = false;
    bool reset_layout = false;
    
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
    double field_mean = 0.0;
    double field_std = 0.0;

    // UI scaling for retina displays
    float ui_scale = 1.0f;
    
    // Font management
    struct {
        ImFont* default_font = nullptr;
        ImFont* mono_font = nullptr;
    } fonts;
    
    // Current tab for legacy support
    int current_tab = 0;

    // Object Drawer functionality
    std::vector<Shape> shapes;
    int selected_shape = -1;
    bool mouse_picking_mode = false;
    ImVec2 grid_size = ImVec2(512.0f, 256.0f); // Will be updated from actual grid
    
    // Current shape creation parameters (for mouse picking)
    int current_shape_type = 0;
    int current_bc_type = 0;
    ImVec2 current_pos = ImVec2(0.0f, 0.0f);
    ImVec2 current_size = ImVec2(10.0f, 10.0f);

    bool init(GLFWwindow *window);
    void begin_frame();
    void draw(int timestep, double sim_time, double dt, double max_velocity,
              double div_l2, double pressure_residual, GLuint texture);
    void end_frame(GLFWwindow *window);
    void shutdown();
    
    // Helper functions for improved UI
    void setup_modern_style();
    void load_fonts();
    void setup_docking_layout();
    void draw_main_menu_bar();
    void draw_simulation_viewport(int timestep, double sim_time, double dt, 
                                 double max_velocity, double div_l2, double pressure_residual,
                                 GLuint texture);
    void draw_cfd_controls_panel(int timestep, double sim_time, double dt, 
                                double max_velocity, double div_l2, double pressure_residual);
    void draw_visualization_panel(GLuint texture);
    void draw_boundary_conditions_panel();
    void draw_object_drawer_panel();
    void draw_performance_stats();
    void draw_field_stats();
    void draw_preset_preview();
    
    // Object drawer functions
    void add_shape(int type, ImVec2 pos, ImVec2 size, int bcType);
    void remove_shape(int index);
    void clear_shapes();
    void apply_shapes_to_simulation(const Grid& grid, State& state);
    ImVec2 screen_to_grid_coords(ImVec2 screen_pos, ImVec2 viewport_pos, ImVec2 viewport_size);
    
    // Utility functions
    std::string get_field_name(Field field) const;
    std::string get_colormap_name(Colormap cmap) const;
    void update_field_statistics(const Grid &g, const State &s);
    std::string generate_ascii_boundary_preview();
    void save_ui_state();
    void load_ui_state();
    void reset_dock_layout(ImGuiID dockspace_id);
};

void update_field_texture(const Grid &g, const State &s, Gui::Field field,
                          GLuint tex, std::vector<unsigned char> &buffer,
                          Gui::Colormap colormap = Gui::Colormap::Viridis,
                          float scale_min = 0.0f, float scale_max = 1.0f,
                          bool normalize = false, float gamma = 1.0f);

// Function to apply shapes to velocity field
void apply_shape_boundary_conditions(const Grid& grid, State& state, const std::vector<Shape>& shapes);
