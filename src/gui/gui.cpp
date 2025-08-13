#include "gui.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <sstream>
#include <iomanip>
#include <filesystem>

#include "../viz.hpp"
#include "../colormaps.hpp"

#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"
#include "imgui_internal.h"  // For docking functions

bool Gui::init(GLFWwindow *window) {
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    
    // Enable docking and viewports
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
    io.ConfigFlags |= ImGuiConfigFlags_ViewportsEnable;
    io.FontAllowUserScaling = true;
    io.IniFilename = "imgui.ini"; // persist layout
    
    // Get UI scale for retina displays
    float x_scale, y_scale;
    glfwGetWindowContentScale(window, &x_scale, &y_scale);
    ui_scale = std::max(1.0f, x_scale);
    
    // Setup modern dark theme
    setup_modern_style();
    
    // Load fonts
    load_fonts();
    
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");
    
    // Load saved UI state
    load_ui_state();
    
    return true;
}

void Gui::setup_modern_style() {
    ImGui::StyleColorsDark();
    
    ImGuiStyle& style = ImGui::GetStyle();
    
    // Apply UI scaling
    style.ScaleAllSizes(ui_scale);
    
    // Modern rounded corners and spacing
    style.WindowRounding = 6.0f;
    style.ChildRounding = 4.0f;
    style.FrameRounding = 4.0f;
    style.GrabRounding = 4.0f;
    style.PopupRounding = 4.0f;
    style.ScrollbarRounding = 4.0f;
    style.TabRounding = 4.0f;
    
    // Improved spacing
    style.WindowPadding = ImVec2(12, 10);
    style.FramePadding = ImVec2(8, 6);
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
    
    // Black menu bar for the top toolbar
    colors[ImGuiCol_MenuBarBg] = ImVec4(0.02f, 0.02f, 0.02f, 0.95f);
    
    // When viewports are enabled we tweak WindowRounding/WindowBg so platform windows can look identical to regular ones.
    ImGuiIO& io = ImGui::GetIO();
    if (io.ConfigFlags & ImGuiConfigFlags_ViewportsEnable) {
        style.WindowRounding = 0.0f;
        style.Colors[ImGuiCol_WindowBg].w = 1.0f;
    }
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
    
    // Fill the main viewport with DockSpace
    ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(vp->Pos);
    ImGui::SetNextWindowSize(vp->Size);
    ImGui::SetNextWindowViewport(vp->ID);
    ImGuiWindowFlags dockspace_flags =
        ImGuiWindowFlags_NoDocking |
        ImGuiWindowFlags_NoTitleBar |
        ImGuiWindowFlags_NoCollapse |
        ImGuiWindowFlags_NoResize |
        ImGuiWindowFlags_NoMove |
        ImGuiWindowFlags_NoBringToFrontOnFocus |
        ImGuiWindowFlags_NoNavFocus |
        ImGuiWindowFlags_NoBackground; // we render our own bg (simulation)
    
    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0,0));
    ImGui::Begin("###DockSpaceHost", nullptr, dockspace_flags);
    ImGui::PopStyleVar(2);

    // Central node passthrough so input reaches simulation when nothing overlies it
    ImGuiDockNodeFlags ds_flags = ImGuiDockNodeFlags_PassthruCentralNode;
    ImGuiID dockspace_id = ImGui::GetID("MainDockSpace");
    ImGui::DockSpace(dockspace_id, ImVec2(0,0), ds_flags);

    // Top toolbar (menu bar look, black)
    draw_main_menu_bar();
    
    ImGui::End(); // DockSpace host
    
    // Draw dockable panels
    if (show_controls) {
        draw_cfd_controls_panel(timestep, sim_time, dt, max_velocity, div_l2, pressure_residual);
    }
    
    if (show_viz) {
        draw_visualization_panel(texture);
    }
    
    if (show_bc) {
        draw_boundary_conditions_panel();
    }
    
    // Draw the simulation viewport (background) - must be after DockSpace
    draw_simulation_viewport(timestep, sim_time, dt, max_velocity, div_l2, pressure_residual, texture);
    
    // Handle layout reset
    if (reset_layout) {
        reset_dock_layout(dockspace_id);
        layout_initialized = false;
        reset_layout = false;
    }
    
    // Initialize layout on first run
    if (!layout_initialized) {
        reset_dock_layout(dockspace_id);
        layout_initialized = true;
    }
}

void Gui::draw_main_menu_bar() {
    if (ImGui::BeginMainMenuBar()) {
        // CFD Controls toggle
        if (ImGui::MenuItem("CFD Controls", nullptr, &show_controls)) {
            // Toggle handled by the bool
        }
        
        ImGui::SameLine();
        
        // Visualization toggle
        if (ImGui::MenuItem("Visualization", nullptr, &show_viz)) {
            // Toggle handled by the bool
        }
        
        ImGui::SameLine();
        
        // Boundary Conditions toggle
        if (ImGui::MenuItem("Boundary Conditions", nullptr, &show_bc)) {
            // Toggle handled by the bool
        }
        
        ImGui::Separator();
        
        // Reset Layout button
        if (ImGui::MenuItem("Reset Layout")) {
            reset_layout = true;
        }
        
        ImGui::EndMainMenuBar();
    }
}

void Gui::draw_simulation_viewport(int timestep, double sim_time, double dt, 
                                 double max_velocity, double div_l2, double pressure_residual,
                                 GLuint texture) {
    // Suppress unused parameter warnings
    (void)dt;
    (void)max_velocity;
    (void)div_l2;
    (void)pressure_residual;
    
    // Create the simulation viewport window
    ImGui::SetNextWindowSizeConstraints(ImVec2(400, 300), ImVec2(FLT_MAX, FLT_MAX));
    
    ImGuiWindowFlags window_flags = ImGuiWindowFlags_NoTitleBar | 
                                   ImGuiWindowFlags_NoCollapse |
                                   ImGuiWindowFlags_NoScrollbar | 
                                   ImGuiWindowFlags_NoScrollWithMouse |
                                   ImGuiWindowFlags_NoBringToFrontOnFocus | 
                                   ImGuiWindowFlags_NoNavFocus |
                                   ImGuiWindowFlags_NoMove;
    
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
    
    if (ImGui::Begin("Simulation Viewport", nullptr, window_flags)) {
        // Get the content region for the simulation
        ImVec2 content_size = ImGui::GetContentRegionAvail();
        
        // Set OpenGL viewport to match the content region
        ImVec2 pos = ImGui::GetCursorScreenPos();
        glViewport(pos.x, pos.y, content_size.x, content_size.y);
        
        // Draw the simulation texture
        ImGui::Image((void*)(intptr_t)texture, content_size);
        
        // Small floating overlay in top-left corner
        ImGui::SetCursorPos(ImVec2(10, 10));
        ImGui::BeginChild("SimulationOverlay", ImVec2(180, 100), true, 
                          ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
        
        // Play/Pause button
        if (ImGui::Button(running ? "⏸" : "▶", ImVec2(40, 25))) {
            running = !running;
        }
        
        ImGui::SameLine();
        
        // Step button
        if (ImGui::Button("⏭", ImVec2(40, 25))) {
            step = true;
        }
        
        ImGui::SameLine();
        
        // Reset button
        if (ImGui::Button("⏹", ImVec2(40, 25))) {
            reset = true;
        }
        
        ImGui::Spacing();
        
        // Tiny stats
        ImGui::Text("%s", get_field_name(field).c_str());
        ImGui::Text("FPS: %.0f | Step: %d", fps, timestep);
        ImGui::Text("Time: %.3f", sim_time);
        
        ImGui::EndChild();
    }
    
    ImGui::PopStyleVar();
    ImGui::End();
}

void Gui::draw_cfd_controls_panel(int timestep, double sim_time, double dt, 
                                 double max_velocity, double div_l2, double pressure_residual) {
    ImGui::SetNextWindowSizeConstraints(ImVec2(360, 420), ImVec2(FLT_MAX, FLT_MAX));
    
    if (ImGui::Begin("CFD Controls", &show_controls, ImGuiWindowFlags_NoSavedSettings)) {
        ImGui::BeginChild("##scroll", ImGui::GetContentRegionAvail(), true, 
                          ImGuiWindowFlags_AlwaysVerticalScrollbar);
        
        // Flow Parameters Section
        if (ImGui::CollapsingHeader("Flow Parameters", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::PushItemWidth(-FLT_MIN);
            
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
            ImGui::PushItemWidth(-FLT_MIN);
            
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
            ImGui::PushItemWidth(-FLT_MIN);
            
            float speed_float = static_cast<float>(sim_speed);
            if (ImGui::SliderFloat("Speed Multiplier", &speed_float, 0.1f, 5.0f, "%.1fx")) {
                sim_speed = static_cast<double>(speed_float);
            }
            
            ImGui::PopItemWidth();
            
            ImGui::Spacing();
            
            // Control buttons
            if (ImGui::Button("Play/Pause", ImVec2(-FLT_MIN, 0))) {
                running = !running;
            }
            
            if (ImGui::Button("Step", ImVec2(-FLT_MIN, 0))) {
                step = true;
            }
            
            if (ImGui::Button("Reset", ImVec2(-FLT_MIN, 0))) {
                reset = true;
            }
            
            if (ImGui::Button("Reset & Run", ImVec2(-FLT_MIN, 0))) {
                reset_run = true;
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
        
        ImGui::EndChild();
    }
    ImGui::End();
}

void Gui::draw_visualization_panel(GLuint texture) {
    ImGui::SetNextWindowSizeConstraints(ImVec2(420, 400), ImVec2(FLT_MAX, FLT_MAX));
    
    if (ImGui::Begin("Visualization", &show_viz, ImGuiWindowFlags_NoSavedSettings)) {
        ImGui::BeginChild("##scroll", ImGui::GetContentRegionAvail(), true, 
                          ImGuiWindowFlags_AlwaysVerticalScrollbar);
        
        // Field Selection
        ImGui::Text("Field Selection");
        ImGui::PushItemWidth(-FLT_MIN);
        
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
        ImGui::PushItemWidth(-FLT_MIN);
        
        const char* colormap_names[] = {"Viridis", "Plasma", "Turbo", "Gray", "RdBu", "Coolwarm"};
        int colormap_int = static_cast<int>(colormap);
        if (ImGui::Combo("##Colormap", &colormap_int, colormap_names, 6)) {
            colormap = static_cast<Colormap>(colormap_int);
        }
        
        ImGui::PopItemWidth();
        
        ImGui::Spacing();
        
        // Rendering Options
        ImGui::Text("Rendering Options");
        ImGui::Checkbox("Vector arrows", &show_vectors);
        ImGui::Checkbox("Contour lines", &show_contours);
        ImGui::Checkbox("Streamlines overlay", &show_streamlines);
        ImGui::Checkbox("Normalize values", &normalize_field);
        
        ImGui::Spacing();
        
        // Color Scale Control
        ImGui::Text("Color Scale");
        ImGui::Checkbox("Auto scale (EMA α=0.08)", &auto_scale);
        
        if (!auto_scale) {
            ImGui::PushItemWidth(-FLT_MIN);
            ImGui::SliderFloat("Min", &color_scale_min, -10.0f, 10.0f, "%.3f");
            ImGui::SliderFloat("Max", &color_scale_max, -10.0f, 10.0f, "%.3f");
            ImGui::PopItemWidth();
        }
        
        ImGui::PushItemWidth(-FLT_MIN);
        ImGui::SliderFloat("Gamma", &gamma, 0.8f, 1.2f, "%.2f");
        ImGui::PopItemWidth();
        
        ImGui::Spacing();
        
        // Field Statistics
        draw_field_stats();
        
        ImGui::Spacing();
        
        // Preview of current visualization
        ImGui::Text("Preview");
        ImGui::Image((void*)(intptr_t)texture, ImVec2(200, 200));
        
        ImGui::EndChild();
    }
    ImGui::End();
}

void Gui::draw_boundary_conditions_panel() {
    ImGui::SetNextWindowSizeConstraints(ImVec2(400, 400), ImVec2(FLT_MAX, FLT_MAX));
    
    if (ImGui::Begin("Boundary Conditions", &show_bc, ImGuiWindowFlags_NoSavedSettings)) {
        ImGui::BeginChild("##scroll", ImGui::GetContentRegionAvail(), true, 
                          ImGuiWindowFlags_AlwaysVerticalScrollbar);
        
        // Preset selector
        ImGui::Text("Preset Configurations");
        ImGui::PushItemWidth(-FLT_MIN);
        
        const char* preset_names[] = {"Jet Plume", "Lid-Driven Cavity", "Periodic Shear"};
        if (ImGui::Combo("Preset", &preset, preset_names, 3)) {
            // Preset selection changed
        }
        
        ImGui::PopItemWidth();
        
        ImGui::Spacing();
        
        // ASCII boundary preview
        ImGui::Text("Boundary Preview");
        ImGui::BeginChild("ASCIIPreview", ImVec2(0, 120), true, 
                          ImGuiWindowFlags_HorizontalScrollbar);
        
        if (fonts.mono_font) {
            ImGui::PushFont(fonts.mono_font);
        }
        
        std::string ascii_preview = generate_ascii_boundary_preview();
        ImGui::TextUnformatted(ascii_preview.c_str());
        
        if (fonts.mono_font) {
            ImGui::PopFont();
        }
        
        ImGui::EndChild();
        
        ImGui::Spacing();
        
        if (ImGui::Button("Apply Preset", ImVec2(-FLT_MIN, 0))) {
            apply_preset = true;
        }
        
        ImGui::Spacing();
        
        ImGui::Text("Custom Presets");
        if (ImGui::Button("Save Current", ImVec2(-FLT_MIN, 0))) {
            // TODO: Implement save functionality
        }
        if (ImGui::Button("Load Custom", ImVec2(-FLT_MIN, 0))) {
            // TODO: Implement load functionality
        }
        
        ImGui::Spacing();
        
        // Boundary condition controls
        if (ImGui::CollapsingHeader("Boundary Settings", ImGuiTreeNodeFlags_DefaultOpen)) {
            // Left BC
            ImGui::Text("Left Boundary");
            const char* bc_types[] = {"Wall", "Inflow", "Outflow", "Periodic", "Moving"};
            int left_type = static_cast<int>(bc.left.type);
            if (ImGui::Combo("Left Type", &left_type, bc_types, 5)) {
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
            
            ImGui::Spacing();
            
            // Right BC
            ImGui::Text("Right Boundary");
            int right_type = static_cast<int>(bc.right.type);
            if (ImGui::Combo("Right Type", &right_type, bc_types, 5)) {
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
            
            ImGui::Spacing();
            
            // Bottom BC
            ImGui::Text("Bottom Boundary");
            int bottom_type = static_cast<int>(bc.bottom.type);
            if (ImGui::Combo("Bottom Type", &bottom_type, bc_types, 5)) {
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
            
            ImGui::Spacing();
            
            // Top BC
            ImGui::Text("Top Boundary");
            int top_type = static_cast<int>(bc.top.type);
            if (ImGui::Combo("Top Type", &top_type, bc_types, 5)) {
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
        }
        
        // Jet settings (if applicable)
        if (bc.left.type == BCType::Inflow && ImGui::CollapsingHeader("Jet Settings")) {
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
        }
        
        ImGui::EndChild();
    }
    ImGui::End();
}

void Gui::reset_dock_layout(ImGuiID dockspace_id) {
    // Clear existing layout
    ImGui::DockBuilderRemoveNode(dockspace_id);
    
    // Create the main dock node
    ImGuiDockNodeFlags flags = ImGuiDockNodeFlags_DockSpace;
    flags = static_cast<ImGuiDockNodeFlags>(flags | ImGuiDockNodeFlags_PassthruCentralNode);
    ImGui::DockBuilderAddNode(dockspace_id, flags);
    ImGui::DockBuilderSetNodeSize(dockspace_id, ImGui::GetMainViewport()->Size);

    // Split main dock into left (controls) and right (BC), keep center for simulation
    ImGuiID dock_main_id = dockspace_id;
    ImGuiID dock_left = ImGui::DockBuilderSplitNode(dock_main_id, ImGuiDir_Left, 0.28f, nullptr, &dock_main_id);
    ImGuiID dock_right = ImGui::DockBuilderSplitNode(dock_main_id, ImGuiDir_Right, 0.28f, nullptr, &dock_main_id);

    // Dock our windows
    ImGui::DockBuilderDockWindow("Simulation Viewport", dock_main_id);
    ImGui::DockBuilderDockWindow("CFD Controls", dock_left);
    ImGui::DockBuilderDockWindow("Visualization", dock_left);
    ImGui::DockBuilderDockWindow("Boundary Conditions", dock_right);

    ImGui::DockBuilderFinish(dockspace_id);
}

void Gui::load_fonts() {
    ImGuiIO& io = ImGui::GetIO();
    
    // Load default font
    fonts.default_font = io.Fonts->AddFontDefault();
    
    // Load monospace font for ASCII preview
    std::string font_path = "assets/RobotoMono-Regular.ttf";
    if (std::filesystem::exists(font_path)) {
        fonts.mono_font = io.Fonts->AddFontFromFileTTF(font_path.c_str(), 15.0f * ui_scale);
    } else {
        // Fallback to default font if monospace font not found
        fonts.mono_font = fonts.default_font;
    }
    
    // Build font atlas
    io.Fonts->Build();
}

void Gui::draw_field_stats() {
    ImGui::Text("Field Statistics");
    
    if (ImGui::BeginTable("FieldStats", 2, ImGuiTableFlags_Borders)) {
        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        ImGui::Text("Min");
        ImGui::TableNextColumn();
        ImGui::Text("%.6f", field_min);
        
        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        ImGui::Text("Max");
        ImGui::TableNextColumn();
        ImGui::Text("%.6f", field_max);
        
        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        ImGui::Text("Mean");
        ImGui::TableNextColumn();
        ImGui::Text("%.6f", field_mean);
        
        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        ImGui::Text("Std Dev");
        ImGui::TableNextColumn();
        ImGui::Text("%.6f", field_std);
        
        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        ImGui::Text("Range");
        ImGui::TableNextColumn();
        ImGui::Text("%.6f", field_max - field_min);
        
        ImGui::EndTable();
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
        case Colormap::Turbo: return "Turbo";
        case Colormap::Grayscale: return "Gray";
        case Colormap::RdBu: return "RdBu";
        case Colormap::Coolwarm: return "Coolwarm";
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
    
    // Calculate mean and standard deviation
    double sum = 0.0;
    double sum_sq = 0.0;
    int count = 0;
    
    for (const auto& val : tmp) {
        if (std::isfinite(val)) {
            sum += val;
            sum_sq += val * val;
            count++;
        }
    }
    
    if (count > 0) {
        field_mean = sum / count;
        double variance = (sum_sq / count) - (field_mean * field_mean);
        field_std = std::sqrt(std::max(0.0, variance));
    } else {
        field_mean = 0.0;
        field_std = 0.0;
    }
    
    // Update auto-scale if enabled
    if (auto_scale) {
        color_scale_min = static_cast<float>(vmin);
        color_scale_max = static_cast<float>(vmax);
    }
}

void Gui::end_frame(GLFWwindow *window) {
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    
    // Update and render additional platform windows
    ImGuiIO& io = ImGui::GetIO();
    if (io.ConfigFlags & ImGuiConfigFlags_ViewportsEnable) {
        GLFWwindow* backup_current_context = glfwGetCurrentContext();
        ImGui::UpdatePlatformWindows();
        ImGui::RenderPlatformWindowsDefault();
        glfwMakeContextCurrent(backup_current_context);
    }
    
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
                          bool normalize, float gamma) {
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
        
        // Apply gamma correction
        norm = std::pow(norm, 1.0 / gamma);
        
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

std::string Gui::generate_ascii_boundary_preview() {
    std::string result;
    result.reserve(200);
    
    // 11x7 grid for boundary preview
    const int width = 11;
    const int height = 7;
    
    // Helper function to get character for boundary type
    auto get_bc_char = [](BCType type, bool is_horizontal) -> char {
        (void)is_horizontal; // Suppress unused parameter warning
        switch (type) {
            case BCType::Wall: return '|';
            case BCType::Inflow: return '>';
            case BCType::Outflow: return ' ';
            case BCType::Periodic: return ':';
            case BCType::Moving: return '>';
            default: return '?';
        }
    };
    
    auto get_horizontal_bc_char = [](BCType type) -> char {
        switch (type) {
            case BCType::Wall: return '-';
            case BCType::Inflow: return '>';
            case BCType::Outflow: return ' ';
            case BCType::Periodic: return '=';
            case BCType::Moving: return '>';
            default: return '?';
        }
    };
    
    // Top row
    result += "  ";
    for (int x = 0; x < width; ++x) {
        result += get_horizontal_bc_char(bc.top.type);
    }
    result += "\n";
    
    // Middle rows
    for (int y = 0; y < height; ++y) {
        result += get_bc_char(bc.left.type, false);
        result += " ";
        
        // Check if this row is in the jet slot for left inflow
        bool in_jet_slot = false;
        if (bc.left.type == BCType::Inflow) {
            double y_pos = (y + 0.5) / height * Ly;
            double jet_start = bc.jet_center - 0.5 * bc.jet_width;
            double jet_end = bc.jet_center + 0.5 * bc.jet_width;
            in_jet_slot = (y_pos >= jet_start && y_pos <= jet_end);
        }
        
        for (int x = 0; x < width; ++x) {
            if (in_jet_slot) {
                result += ">";
            } else {
                result += " ";
            }
        }
        result += " ";
        result += get_bc_char(bc.right.type, false);
        result += "\n";
    }
    
    // Bottom row
    result += "  ";
    for (int x = 0; x < width; ++x) {
        result += get_horizontal_bc_char(bc.bottom.type);
    }
    result += "\n";
    
    // Add legend
    result += "\nLegend:\n";
    result += "- | Wall boundary\n";
    result += "> > Inflow/Moving\n";
    result += "  Outflow\n";
    result += ": = Periodic\n";
    
    return result;
}

void Gui::save_ui_state() {
    // TODO: Implement JSON-based UI state saving
    // This would save current tab, panel widths, auto-scale settings, etc.
}

void Gui::load_ui_state() {
    // TODO: Implement JSON-based UI state loading
    // This would restore the saved UI state
} 