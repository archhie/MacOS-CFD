#include "gui.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <sstream>
#include <iomanip>

#include "../viz.hpp"
#include "../colormaps.hpp"

#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

bool Gui::init(GLFWwindow *window) {
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    
    // Setup modern dark theme
    setup_modern_style();
    
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");
    return true;
}

void Gui::setup_modern_style() {
    ImGui::StyleColorsDark();
    
    ImGuiStyle& style = ImGui::GetStyle();
    
    // Modern rounded corners and spacing
    style.WindowRounding = 6.0f;
    style.ChildRounding = 4.0f;
    style.FrameRounding = 4.0f;
    style.GrabRounding = 4.0f;
    style.PopupRounding = 4.0f;
    style.ScrollbarRounding = 4.0f;
    style.TabRounding = 4.0f;
    
    // Improved spacing
    style.WindowPadding = ImVec2(12, 12);
    style.FramePadding = ImVec2(8, 4);
    style.ItemSpacing = ImVec2(8, 6);
    style.ItemInnerSpacing = ImVec2(6, 4);
    style.ScrollbarSize = 12.0f;
    style.GrabMinSize = 8.0f;
    
    // Modern colors with better contrast
    ImVec4* colors = style.Colors;
    colors[ImGuiCol_WindowBg] = ImVec4(0.08f, 0.08f, 0.10f, 0.95f);
    colors[ImGuiCol_ChildBg] = ImVec4(0.06f, 0.06f, 0.08f, 0.95f);
    colors[ImGuiCol_PopupBg] = ImVec4(0.10f, 0.10f, 0.12f, 0.95f);
    colors[ImGuiCol_Border] = ImVec4(0.20f, 0.20f, 0.25f, 0.50f);
    colors[ImGuiCol_BorderShadow] = ImVec4(0.00f, 0.00f, 0.00f, 0.00f);
    colors[ImGuiCol_FrameBg] = ImVec4(0.12f, 0.12f, 0.15f, 0.95f);
    colors[ImGuiCol_FrameBgHovered] = ImVec4(0.15f, 0.15f, 0.18f, 0.95f);
    colors[ImGuiCol_FrameBgActive] = ImVec4(0.18f, 0.18f, 0.22f, 0.95f);
    colors[ImGuiCol_TitleBg] = ImVec4(0.10f, 0.10f, 0.12f, 0.95f);
    colors[ImGuiCol_TitleBgActive] = ImVec4(0.12f, 0.12f, 0.15f, 0.95f);
    colors[ImGuiCol_Button] = ImVec4(0.15f, 0.15f, 0.18f, 0.95f);
    colors[ImGuiCol_ButtonHovered] = ImVec4(0.20f, 0.20f, 0.25f, 0.95f);
    colors[ImGuiCol_ButtonActive] = ImVec4(0.25f, 0.25f, 0.30f, 0.95f);
    colors[ImGuiCol_Header] = ImVec4(0.12f, 0.12f, 0.15f, 0.95f);
    colors[ImGuiCol_HeaderHovered] = ImVec4(0.15f, 0.15f, 0.18f, 0.95f);
    colors[ImGuiCol_HeaderActive] = ImVec4(0.18f, 0.18f, 0.22f, 0.95f);
    colors[ImGuiCol_SliderGrab] = ImVec4(0.30f, 0.60f, 0.90f, 0.95f);
    colors[ImGuiCol_SliderGrabActive] = ImVec4(0.40f, 0.70f, 1.00f, 0.95f);
    colors[ImGuiCol_CheckMark] = ImVec4(0.30f, 0.60f, 0.90f, 0.95f);
    colors[ImGuiCol_TextSelectedBg] = ImVec4(0.30f, 0.60f, 0.90f, 0.35f);
}

void Gui::begin_frame() {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
}

void Gui::draw(int timestep, double sim_time, double dt, double max_velocity,
               double div_l2, double pressure_residual, GLuint texture) {
    
    // Update performance stats
    double current_time = ImGui::GetTime();
    if (last_frame_time > 0.0) {
        fps = 1.0 / (current_time - last_frame_time);
    }
    last_frame_time = current_time;
    
    // Draw the three main panels
    draw_cfd_controls(timestep, sim_time, dt, max_velocity, div_l2, pressure_residual);
    draw_visualization_panel(texture);
    draw_boundary_conditions_panel();
}

void Gui::draw_cfd_controls(int timestep, double sim_time, double dt, 
                           double max_velocity, double div_l2, double pressure_residual) {
    ImGui::SetNextWindowSize(ImVec2(320, 500), ImGuiCond_FirstUseEver);
    ImGui::Begin("CFD Controls", nullptr, ImGuiWindowFlags_NoCollapse);
    
    // Flow Parameters Section
    if (ImGui::CollapsingHeader("Flow Parameters", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::PushItemWidth(-1);
        
        float re_float = static_cast<float>(Re);
        if (ImGui::SliderFloat("Reynolds Number (Re)", &re_float, 10.0f, 1000.0f, "%.1f")) {
            Re = static_cast<double>(re_float);
        }
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Reynolds number controls the relative importance of inertial vs viscous forces");
        }
        
        float cfl_float = static_cast<float>(CFL);
        if (ImGui::SliderFloat("CFL Number", &cfl_float, 0.001f, 0.5f, "%.3f")) {
            CFL = static_cast<double>(cfl_float);
        }
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Courant-Friedrichs-Lewy number controls numerical stability");
        }
        
        ImGui::PopItemWidth();
    }
    
    // Timestep Control Section
    if (ImGui::CollapsingHeader("Timestep Control", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::PushItemWidth(-1);
        
        float dt_float = static_cast<float>(this->dt);
        if (ImGui::SliderFloat("dt Override", &dt_float, 0.0001f, 0.01f, "%.4f")) {
            this->dt = static_cast<double>(dt_float);
        }
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Override automatic timestep calculation with fixed value");
        }
        
        ImGui::PopItemWidth();
    }
    
    // Simulation Speed Section
    if (ImGui::CollapsingHeader("Simulation Speed", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::PushItemWidth(-1);
        
        float speed_float = static_cast<float>(sim_speed);
        if (ImGui::SliderFloat("Speed Multiplier", &speed_float, 0.1f, 5.0f, "%.1fx")) {
            sim_speed = static_cast<double>(speed_float);
        }
        
        ImGui::PopItemWidth();
        
        ImGui::Spacing();
        
        // Control buttons
        if (ImGui::Button("Play/Pause", ImVec2(-1, 0))) {
            running = !running;
        }
        
        if (ImGui::Button("Step", ImVec2(-1, 0))) {
            step = true;
        }
        
        if (ImGui::Button("Reset", ImVec2(-1, 0))) {
            reset = true;
        }
        
        if (ImGui::Button("Reset & Run", ImVec2(-1, 0))) {
            reset_run = true;
        }
        
        if (ImGui::Button("Export State", ImVec2(-1, 0))) {
            // TODO: Implement export functionality
        }
    }
    
    // Simulation Info Section
    if (ImGui::CollapsingHeader("Simulation Info", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Text("Timestep: %d", timestep);
        ImGui::Text("Time: %.6f", sim_time);
        ImGui::Text("dt: %.6f", dt);
        ImGui::Text("Max velocity: %.6f", max_velocity);
        ImGui::Text("Divergence L2: %.2e", div_l2);
        ImGui::Text("Pressure residual: %.2e", pressure_residual);
        ImGui::Text("Inflow U=%.3f w=%.3f", bc.left.inflow_u, bc.jet_width);
        
        if (inflow_inactive) {
            ImGui::TextColored(ImVec4(1, 0.3f, 0.3f, 1), "Inflow inactive");
        }
        
        ImGui::Separator();
        
        // Performance stats
        ImGui::Text("FPS: %.1f", fps);
        ImGui::Text("Est. time to next frame: %.3f ms", 1000.0 / std::max(1.0, fps));
    }
    
    ImGui::End();
}

void Gui::draw_visualization_panel(GLuint texture) {
    ImGui::SetNextWindowSize(ImVec2(450, 600), ImGuiCond_FirstUseEver);
    ImGui::Begin("Visualization", nullptr, ImGuiWindowFlags_NoCollapse);
    
    // Field Selection
    ImGui::Text("Field Selection");
    ImGui::PushItemWidth(-1);
    
    const char* field_names[] = {"U (Horizontal Velocity)", "V (Vertical Velocity)", 
                                "Speed (Velocity Magnitude)", "Pressure", "Vorticity", 
                                "Streamlines", "Temperature"};
    int field_int = static_cast<int>(field);
    if (ImGui::Combo("##Field", &field_int, field_names, 7)) {
        field = static_cast<Field>(field_int);
    }
    
    ImGui::PopItemWidth();
    
    ImGui::Spacing();
    
    // Colormap Selection
    ImGui::Text("Colormap");
    ImGui::PushItemWidth(-1);
    
    const char* colormap_names[] = {"Viridis", "Plasma", "Jet", "Grayscale", "Turbo", "RdYlBu"};
    int colormap_int = static_cast<int>(colormap);
    if (ImGui::Combo("##Colormap", &colormap_int, colormap_names, 6)) {
        colormap = static_cast<Colormap>(colormap_int);
    }
    
    ImGui::PopItemWidth();
    
    ImGui::Spacing();
    
    // Visualization Options
    ImGui::Text("Options");
    ImGui::Checkbox("Show Vector Arrows", &show_vectors);
    ImGui::Checkbox("Normalize Field Values", &normalize_field);
    
    ImGui::Spacing();
    
    // Color Scale Control
    ImGui::Text("Color Scale");
    ImGui::Checkbox("Auto Scale", &auto_scale);
    
    if (!auto_scale) {
        ImGui::PushItemWidth(-1);
        ImGui::SliderFloat("Min", &color_scale_min, -10.0f, 10.0f, "%.3f");
        ImGui::SliderFloat("Max", &color_scale_max, -10.0f, 10.0f, "%.3f");
        ImGui::PopItemWidth();
    }
    
    ImGui::Spacing();
    
    // Field Statistics
    draw_field_stats();
    
    ImGui::Spacing();
    
    // Main visualization image
    ImGui::Text("Field Visualization");
    ImGui::Image((void*)(intptr_t)texture, ImVec2(400, 400));
    
    ImGui::End();
}

void Gui::draw_boundary_conditions_panel() {
    ImGui::SetNextWindowSize(ImVec2(350, 500), ImGuiCond_FirstUseEver);
    ImGui::Begin("Boundary Conditions", nullptr, ImGuiWindowFlags_NoCollapse);
    
    if (ImGui::BeginTabBar("BC Tabs")) {
        if (ImGui::BeginTabItem("Presets")) {
            ImGui::Text("Preset Configurations");
            
            const char* preset_names[] = {"Jet Plume", "Lid-Driven Cavity", "Periodic Shear"};
            if (ImGui::Combo("Preset", &preset, preset_names, 3)) {
                // Preset selection changed
            }
            
            ImGui::Spacing();
            
            // Preset preview
            draw_preset_preview();
            
            ImGui::Spacing();
            
            if (ImGui::Button("Apply Preset", ImVec2(-1, 0))) {
                apply_preset = true;
            }
            
            ImGui::Spacing();
            
            ImGui::Text("Custom Presets");
            if (ImGui::Button("Save Current", ImVec2(-1, 0))) {
                // TODO: Implement save functionality
            }
            if (ImGui::Button("Load Custom", ImVec2(-1, 0))) {
                // TODO: Implement load functionality
            }
            
            ImGui::EndTabItem();
        }
        
        if (ImGui::BeginTabItem("Left BC")) {
            ImGui::Text("Left Boundary");
            const char* bc_types[] = {"Wall", "Inflow", "Outflow", "Periodic", "Moving"};
            int left_type = static_cast<int>(bc.left.type);
            if (ImGui::Combo("Type", &left_type, bc_types, 5)) {
                bc.left.type = static_cast<BCType>(left_type);
            }
            
            if (bc.left.type == BCType::Inflow) {
                float uin = static_cast<float>(bc.left.inflow_u);
                if (ImGui::SliderFloat("U_in", &uin, -5.0f, 5.0f)) {
                    bc.left.inflow_u = uin;
                }
            } else if (bc.left.type == BCType::Moving) {
                float moving = static_cast<float>(bc.left.moving);
                if (ImGui::SliderFloat("Speed", &moving, -5.0f, 5.0f)) {
                    bc.left.moving = moving;
                }
            }
            ImGui::EndTabItem();
        }
        
        if (ImGui::BeginTabItem("Right BC")) {
            ImGui::Text("Right Boundary");
            const char* bc_types[] = {"Wall", "Inflow", "Outflow", "Periodic", "Moving"};
            int right_type = static_cast<int>(bc.right.type);
            if (ImGui::Combo("Type", &right_type, bc_types, 5)) {
                bc.right.type = static_cast<BCType>(right_type);
            }
            
            if (bc.right.type == BCType::Inflow) {
                float uin = static_cast<float>(bc.right.inflow_u);
                if (ImGui::SliderFloat("U_in", &uin, -5.0f, 5.0f)) {
                    bc.right.inflow_u = uin;
                }
            } else if (bc.right.type == BCType::Moving) {
                float moving = static_cast<float>(bc.right.moving);
                if (ImGui::SliderFloat("Speed", &moving, -5.0f, 5.0f)) {
                    bc.right.moving = moving;
                }
            }
            ImGui::EndTabItem();
        }
        
        if (ImGui::BeginTabItem("Bottom BC")) {
            ImGui::Text("Bottom Boundary");
            const char* bc_types[] = {"Wall", "Inflow", "Outflow", "Periodic", "Moving"};
            int bottom_type = static_cast<int>(bc.bottom.type);
            if (ImGui::Combo("Type", &bottom_type, bc_types, 5)) {
                bc.bottom.type = static_cast<BCType>(bottom_type);
            }
            
            if (bc.bottom.type == BCType::Inflow) {
                float vin = static_cast<float>(bc.bottom.inflow_v);
                if (ImGui::SliderFloat("V_in", &vin, -5.0f, 5.0f)) {
                    bc.bottom.inflow_v = vin;
                }
            } else if (bc.bottom.type == BCType::Moving) {
                float moving = static_cast<float>(bc.bottom.moving);
                if (ImGui::SliderFloat("Speed", &moving, -5.0f, 5.0f)) {
                    bc.bottom.moving = moving;
                }
            }
            ImGui::EndTabItem();
        }
        
        if (ImGui::BeginTabItem("Top BC")) {
            ImGui::Text("Top Boundary");
            const char* bc_types[] = {"Wall", "Inflow", "Outflow", "Periodic", "Moving"};
            int top_type = static_cast<int>(bc.top.type);
            if (ImGui::Combo("Type", &top_type, bc_types, 5)) {
                bc.top.type = static_cast<BCType>(top_type);
            }
            
            if (bc.top.type == BCType::Inflow) {
                float vin = static_cast<float>(bc.top.inflow_v);
                if (ImGui::SliderFloat("V_in", &vin, -5.0f, 5.0f)) {
                    bc.top.inflow_v = vin;
                }
            } else if (bc.top.type == BCType::Moving) {
                float moving = static_cast<float>(bc.top.moving);
                if (ImGui::SliderFloat("Speed", &moving, -5.0f, 5.0f)) {
                    bc.top.moving = moving;
                }
            }
            ImGui::EndTabItem();
        }
        
        if (ImGui::BeginTabItem("Jet Settings")) {
            ImGui::Text("Jet Configuration");
            
            float uin = static_cast<float>(bc.left.inflow_u);
            float y0 = static_cast<float>(bc.jet_center);
            float w = static_cast<float>(bc.jet_width);
            float delta = static_cast<float>(bc.jet_thickness);
            float eps = static_cast<float>(bc.jet_eps);
            int k = bc.jet_k;
            
            ImGui::SliderFloat("U_in", &uin, -5.0f, 5.0f);
            ImGui::SliderFloat("Center", &y0, 0.0f, static_cast<float>(Ly));
            ImGui::SliderFloat("Width", &w, 0.01f, static_cast<float>(Ly));
            ImGui::SliderFloat("Thickness", &delta, 0.001f, 0.1f);
            ImGui::SliderFloat("Perturbation", &eps, 0.0f, 1.0f);
            ImGui::SliderInt("k", &k, 1, 16);
            
            bc.left.inflow_u = uin;
            bc.jet_center = y0;
            bc.jet_width = w;
            bc.jet_thickness = delta;
            bc.jet_eps = eps;
            bc.jet_k = k;
            
            ImGui::Spacing();
            ImGui::Text("Jet Profile");
            float profile[64];
            for (int i = 0; i < 64; ++i) {
                double y = (i + 0.5) / 64.0 * Ly;
                profile[i] = static_cast<float>(bc.left.inflow_u *
                    0.5 * (std::tanh((y - (bc.jet_center - 0.5 * bc.jet_width)) / bc.jet_thickness) -
                           std::tanh((y - (bc.jet_center + 0.5 * bc.jet_width)) / bc.jet_thickness)));
            }
            ImGui::PlotLines("u_in(y)", profile, 64, 0, nullptr, -5.0f, 5.0f, ImVec2(0, 50));
            ImGui::EndTabItem();
        }
        
        if (ImGui::BeginTabItem("Preset Settings")) {
            ImGui::Text("Preset-Specific Settings");
            
            if (preset == static_cast<int>(Preset::LidDrivenCavity)) {
                float lid = static_cast<float>(bc.top.moving);
                if (ImGui::SliderFloat("Lid speed", &lid, -5.0f, 5.0f)) {
                    bc.top.moving = lid;
                }
            } else if (preset == static_cast<int>(Preset::PeriodicShear)) {
                float u0 = static_cast<float>(bc.left.inflow_u);
                if (ImGui::SliderFloat("U0", &u0, -5.0f, 5.0f)) {
                    bc.left.inflow_u = u0;
                }
            }
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
    }
    ImGui::End();
}

void Gui::draw_field_stats() {
    ImGui::Text("Field Statistics");
    ImGui::Text("Min: %.6f", field_min);
    ImGui::Text("Max: %.6f", field_max);
    ImGui::Text("Range: %.6f", field_max - field_min);
}

void Gui::draw_preset_preview() {
    ImGui::Text("Preview");
    
    // Simple ASCII art preview based on preset
    switch (preset) {
        case static_cast<int>(Preset::JetPlume):
            ImGui::Text("Jet Plume:");
            ImGui::Text("  ┌─────────────┐");
            ImGui::Text("  │    →→→      │");
            ImGui::Text("  │    →→→      │");
            ImGui::Text("  │    →→→      │");
            ImGui::Text("  └─────────────┘");
            break;
        case static_cast<int>(Preset::LidDrivenCavity):
            ImGui::Text("Lid-Driven Cavity:");
            ImGui::Text("  ┌─────────────┐");
            ImGui::Text("  │ →→→→→→→→→→→ │");
            ImGui::Text("  │             │");
            ImGui::Text("  │             │");
            ImGui::Text("  └─────────────┘");
            break;
        case static_cast<int>(Preset::PeriodicShear):
            ImGui::Text("Periodic Shear:");
            ImGui::Text("  ┌─────────────┐");
            ImGui::Text("  │ →→→→→→→→→→→ │");
            ImGui::Text("  │ ←←←←←←←←←←← │");
            ImGui::Text("  │ →→→→→→→→→→→ │");
            ImGui::Text("  └─────────────┘");
            break;
    }
}

std::string Gui::get_field_name(Field field) const {
    switch (field) {
        case Field::U: return "U (Horizontal Velocity)";
        case Field::V: return "V (Vertical Velocity)";
        case Field::Speed: return "Speed (Velocity Magnitude)";
        case Field::Pressure: return "Pressure";
        case Field::Vorticity: return "Vorticity";
        case Field::Streamlines: return "Streamlines";
        case Field::Temperature: return "Temperature";
        default: return "Unknown";
    }
}

std::string Gui::get_colormap_name(Colormap cmap) const {
    switch (cmap) {
        case Colormap::Viridis: return "Viridis";
        case Colormap::Plasma: return "Plasma";
        case Colormap::Jet: return "Jet";
        case Colormap::Grayscale: return "Grayscale";
        case Colormap::Turbo: return "Turbo";
        case Colormap::RdYlBu: return "RdYlBu";
        default: return "Unknown";
    }
}

void Gui::update_field_statistics(const Grid &g, const State &s) {
    // Compute field statistics for the current field selection
    std::vector<double> tmp;
    double vmin = 0.0, vmax = 0.0;
    
    ScalarField sf;
    switch (field) {
    case Field::U: sf = ScalarField::U; break;
    case Field::V: sf = ScalarField::V; break;
    case Field::Speed: sf = ScalarField::Speed; break;
    case Field::Pressure: sf = ScalarField::Pressure; break;
    case Field::Vorticity: sf = ScalarField::Vorticity; break;
    case Field::Streamlines: sf = ScalarField::Streamlines; break;
    case Field::Temperature: sf = ScalarField::Temperature; break;
    }
    
    compute_scalar(g, s, sf, tmp, vmin, vmax);
    
    field_min = vmin;
    field_max = vmax;
    
    // Update auto-scale if enabled
    if (auto_scale) {
        color_scale_min = static_cast<float>(vmin);
        color_scale_max = static_cast<float>(vmax);
    }
}

void Gui::end_frame(GLFWwindow *window) {
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    glfwSwapBuffers(window);
}

void Gui::shutdown() {
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}

void update_field_texture(const Grid &g, const State &s, Gui::Field field,
                          GLuint tex, std::vector<unsigned char> &buffer,
                          Gui::Colormap colormap, float scale_min, float scale_max,
                          bool normalize) {
    int nx = g.nx;
    int ny = g.ny;
    std::vector<double> tmp;
    double vmin = 0.0, vmax = 0.0;
    
    ScalarField sf;
    switch (field) {
    case Gui::Field::U: sf = ScalarField::U; break;
    case Gui::Field::V: sf = ScalarField::V; break;
    case Gui::Field::Speed: sf = ScalarField::Speed; break;
    case Gui::Field::Pressure: sf = ScalarField::Pressure; break;
    case Gui::Field::Vorticity: sf = ScalarField::Vorticity; break;
    case Gui::Field::Streamlines: sf = ScalarField::Streamlines; break;
    case Gui::Field::Temperature: sf = ScalarField::Temperature; break;
    }
    
    compute_scalar(g, s, sf, tmp, vmin, vmax);
    
    // Use provided scale if not auto-scaling
    if (scale_min != 0.0f || scale_max != 1.0f) {
        vmin = scale_min;
        vmax = scale_max;
    }
    
    buffer.resize(static_cast<size_t>(nx) * ny * 3);
    double scale = 1.0 / std::max(1e-6, vmax - vmin);
    
    for (int idx = 0; idx < nx * ny; ++idx) {
        double val = tmp[idx];
        if (normalize) {
            val = std::abs(val);
        }
        
        double norm = (val - vmin) * scale;
        norm = std::clamp(norm, 0.0, 1.0);
        
        unsigned char r, g, b;
        Colormaps::apply_colormap(static_cast<int>(colormap), static_cast<float>(norm), r, g, b);
        
        buffer[3 * idx + 0] = r;
        buffer[3 * idx + 1] = g;
        buffer[3 * idx + 2] = b;
    }
    
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, nx, ny, 0, GL_RGB,
                 GL_UNSIGNED_BYTE, buffer.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
}
