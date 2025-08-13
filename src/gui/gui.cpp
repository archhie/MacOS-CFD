#include "gui.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>

#include "../viz.hpp"

#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

bool Gui::init(GLFWwindow *window) {
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");
    return true;
}

void Gui::begin_frame() {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
}

void Gui::draw(int timestep, double sim_time, double dt, double max_velocity,
               double div_l2, double pressure_residual, GLuint texture) {
    ImGui::Begin("CFD Controls");
    
    // Convert doubles to floats for ImGui sliders
    float re_float = static_cast<float>(Re);
    float cfl_float = static_cast<float>(CFL);
    float dt_float = static_cast<float>(this->dt);
    
    ImGui::SliderFloat("Re", &re_float, 1.0f, 10000.0f);
    ImGui::SliderFloat("CFL", &cfl_float, 0.1f, 2.0f);
    ImGui::SliderFloat("dt", &dt_float, 0.0001f, 0.01f);
    
    Re = static_cast<double>(re_float);
    CFL = static_cast<double>(cfl_float);
    this->dt = static_cast<double>(dt_float);
    
    ImGui::Checkbox("Running", &running);
    if (ImGui::Button("Step")) {
        step = true;
    }
    if (ImGui::Button("Reset")) {
        reset = true;
    }
    if (ImGui::Button("Reset & Run")) {
        reset_run = true;
    }
    
    ImGui::Text("Timestep: %d", timestep);
    ImGui::Text("Time: %.6f", sim_time);
    ImGui::Text("Max velocity: %.6f", max_velocity);
    ImGui::Text("Divergence L2: %.2e", div_l2);
    ImGui::Text("Pressure residual: %.2e", pressure_residual);
    ImGui::Text("Inflow U=%.3f w=%.3f", bc.left.inflow_u, bc.jet_width);
    if (inflow_inactive)
        ImGui::TextColored(ImVec4(1, 0.3f, 0.3f, 1), "Inflow inactive");
    
    ImGui::End();
    
    ImGui::Begin("Visualization");
    const char* field_names[] = {"U", "V", "Speed", "Pressure", "Vorticity"};
    int field_int = static_cast<int>(field);
    ImGui::Combo("Field", &field_int, field_names, 5);
    field = static_cast<Field>(field_int);
    
    ImGui::Image((void*)(intptr_t)texture, ImVec2(400, 400));
    ImGui::End();
    
    ImGui::Begin("Boundary Conditions");
    if (ImGui::BeginTabBar("BC Tabs")) {
        if (ImGui::BeginTabItem("Presets")) {
            const char* preset_names[] = {"Jet Plume", "Lid-Driven Cavity", "Periodic Shear"};
            ImGui::Combo("Preset", &preset, preset_names, 3);
            if (ImGui::Button("Apply Preset")) {
                apply_preset = true;
            }
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Left BC")) {
            const char* bc_types[] = {"Wall", "Inflow", "Outflow", "Periodic", "Moving"};
            int left_type = static_cast<int>(bc.left.type);
            ImGui::Combo("Type", &left_type, bc_types, 5);
            bc.left.type = static_cast<BCType>(left_type);
            
            if (bc.left.type == BCType::Inflow) {
                float uin = static_cast<float>(bc.left.inflow_u);
                ImGui::SliderFloat("U_in", &uin, -5.0f, 5.0f);
                bc.left.inflow_u = uin;
            } else if (bc.left.type == BCType::Moving) {
                float moving = static_cast<float>(bc.left.moving);
                ImGui::SliderFloat("Speed", &moving, -5.0f, 5.0f);
                bc.left.moving = moving;
            }
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Right BC")) {
            const char* bc_types[] = {"Wall", "Inflow", "Outflow", "Periodic", "Moving"};
            int right_type = static_cast<int>(bc.right.type);
            ImGui::Combo("Type", &right_type, bc_types, 5);
            bc.right.type = static_cast<BCType>(right_type);
            
            if (bc.right.type == BCType::Inflow) {
                float uin = static_cast<float>(bc.right.inflow_u);
                ImGui::SliderFloat("U_in", &uin, -5.0f, 5.0f);
                bc.right.inflow_u = uin;
            } else if (bc.right.type == BCType::Moving) {
                float moving = static_cast<float>(bc.right.moving);
                ImGui::SliderFloat("Speed", &moving, -5.0f, 5.0f);
                bc.right.moving = moving;
            }
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Bottom BC")) {
            const char* bc_types[] = {"Wall", "Inflow", "Outflow", "Periodic", "Moving"};
            int bottom_type = static_cast<int>(bc.bottom.type);
            ImGui::Combo("Type", &bottom_type, bc_types, 5);
            bc.bottom.type = static_cast<BCType>(bottom_type);
            
            if (bc.bottom.type == BCType::Inflow) {
                float vin = static_cast<float>(bc.bottom.inflow_v);
                ImGui::SliderFloat("V_in", &vin, -5.0f, 5.0f);
                bc.bottom.inflow_v = vin;
            } else if (bc.bottom.type == BCType::Moving) {
                float moving = static_cast<float>(bc.bottom.moving);
                ImGui::SliderFloat("Speed", &moving, -5.0f, 5.0f);
                bc.bottom.moving = moving;
            }
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Top BC")) {
            const char* bc_types[] = {"Wall", "Inflow", "Outflow", "Periodic", "Moving"};
            int top_type = static_cast<int>(bc.top.type);
            ImGui::Combo("Type", &top_type, bc_types, 5);
            bc.top.type = static_cast<BCType>(top_type);
            
            if (bc.top.type == BCType::Inflow) {
                float vin = static_cast<float>(bc.top.inflow_v);
                ImGui::SliderFloat("V_in", &vin, -5.0f, 5.0f);
                bc.top.inflow_v = vin;
            } else if (bc.top.type == BCType::Moving) {
                float moving = static_cast<float>(bc.top.moving);
                ImGui::SliderFloat("Speed", &moving, -5.0f, 5.0f);
                bc.top.moving = moving;
            }
            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem("Jet Settings")) {
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
            if (preset == static_cast<int>(Preset::LidDrivenCavity)) {
                float lid = static_cast<float>(bc.top.moving);
                ImGui::SliderFloat("Lid speed", &lid, -5.0f, 5.0f);
                bc.top.moving = lid;
            } else if (preset == static_cast<int>(Preset::PeriodicShear)) {
                float u0 = static_cast<float>(bc.left.inflow_u);
                ImGui::SliderFloat("U0", &u0, -5.0f, 5.0f);
                bc.left.inflow_u = u0;
            }
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
    }
    ImGui::End();
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
                          GLuint tex, std::vector<unsigned char> &buffer) {
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
    }
    compute_scalar(g, s, sf, tmp, vmin, vmax);
    buffer.resize(static_cast<size_t>(nx) * ny * 3);
    double scale = 1.0 / (vmax - vmin);
    for (int idx = 0; idx < nx * ny; ++idx) {
        double norm = (tmp[idx] - vmin) * scale;
        unsigned char c = static_cast<unsigned char>(std::clamp(norm, 0.0, 1.0) * 255.0);
        buffer[3 * idx + 0] = c;
        buffer[3 * idx + 1] = c;
        buffer[3 * idx + 2] = c;
    }
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, nx, ny, 0, GL_RGB,
                 GL_UNSIGNED_BYTE, buffer.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
}
