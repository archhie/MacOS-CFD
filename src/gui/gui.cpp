#include "gui.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>

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

void Gui::draw(int timestep, double sim_time, double dt_used,
               double max_velocity, double div_l2,
               double pressure_residual, GLuint texture) {
    ImGui::Begin("CFD Controls");
    if (ImGui::BeginTabBar("Tabs")) {
        if (ImGui::BeginTabItem("Controls")) {

    // Convert doubles to floats for ImGui sliders
    float re_float = static_cast<float>(Re);
    float cfl_float = static_cast<float>(CFL);
    float dt_float = static_cast<float>(dt);

    ImGui::SliderFloat("Re", &re_float, 100.0f, 10000.0f);
    ImGui::SliderFloat("CFL", &cfl_float, 0.1f, 1.0f);
    ImGui::SliderFloat("dt", &dt_float, 1e-4f, 1e-1f, "%.5f");

    // Convert back to doubles
    Re = static_cast<double>(re_float);
    CFL = static_cast<double>(cfl_float);
    dt = static_cast<double>(dt_float);

    if (ImGui::Button(running ? "Pause" : "Run"))
        running = !running;
    ImGui::SameLine();
    if (ImGui::Button("Step"))
        step = true;
    ImGui::SameLine();
    if (ImGui::Button("Reset"))
        reset = true;

    ImGui::Separator();
    ImGui::Text("Timestep: %d", timestep);
    ImGui::Text("Sim time: %.3f", sim_time);
    ImGui::Text("dt: %.4e", dt_used);
    ImGui::Text("Max |u|: %.3f", max_velocity);
    ImGui::Text("Div L2: %.3e", div_l2);
    ImGui::Text("Pressure residual: %.3e", pressure_residual);

    ImGui::Separator();
    ImGui::Text("Boundary Conditions");
    const char *bc_items[] = {"Wall", "Moving", "Inflow", "Outflow", "Periodic"};
    int left = static_cast<int>(bc.left.type);
    if (ImGui::Combo("Left", &left, bc_items, IM_ARRAYSIZE(bc_items)))
        bc.left.type = static_cast<BCType>(left);
    if (bc.left.type == BCType::Moving) {
        float sp = static_cast<float>(bc.left.moving);
        ImGui::SliderFloat("Left speed", &sp, -5.0f, 5.0f);
        bc.left.moving = sp;
    } else if (bc.left.type == BCType::Inflow) {
        float u = static_cast<float>(bc.left.inflow_u);
        float v = static_cast<float>(bc.left.inflow_v);
        ImGui::SliderFloat("U_in", &u, -5.0f, 5.0f);
        ImGui::SliderFloat("V_in", &v, -5.0f, 5.0f);
        bc.left.inflow_u = u;
        bc.left.inflow_v = v;
        float y0 = static_cast<float>(bc.jet_center);
        float w = static_cast<float>(bc.jet_width);
        ImGui::SliderFloat("y0", &y0, 0.0f, static_cast<float>(Ly));
        ImGui::SliderFloat("width", &w, 0.0f, static_cast<float>(Ly));
        bc.jet_center = y0;
        bc.jet_width = w;
        float profile[64];
        for (int i = 0; i < 64; ++i) {
            double y = (i + 0.5) / 64.0 * Ly;
            profile[i] = (std::fabs(y - bc.jet_center) <=
                          0.5 * bc.jet_width)
                             ? static_cast<float>(bc.left.inflow_u)
                             : 0.0f;
        }
        ImGui::PlotLines("u_in(y)", profile, 64, 0, nullptr, -5.0f, 5.0f,
                         ImVec2(0, 50));
    }

    int right = static_cast<int>(bc.right.type);
    if (ImGui::Combo("Right", &right, bc_items, IM_ARRAYSIZE(bc_items)))
        bc.right.type = static_cast<BCType>(right);
    if (bc.right.type == BCType::Moving) {
        float sp = static_cast<float>(bc.right.moving);
        ImGui::SliderFloat("Right speed", &sp, -5.0f, 5.0f);
        bc.right.moving = sp;
    } else if (bc.right.type == BCType::Inflow) {
        float u = static_cast<float>(bc.right.inflow_u);
        float v = static_cast<float>(bc.right.inflow_v);
        ImGui::SliderFloat("Right U", &u, -5.0f, 5.0f);
        ImGui::SliderFloat("Right V", &v, -5.0f, 5.0f);
        bc.right.inflow_u = u;
        bc.right.inflow_v = v;
    }

    int bottom = static_cast<int>(bc.bottom.type);
    if (ImGui::Combo("Bottom", &bottom, bc_items, IM_ARRAYSIZE(bc_items)))
        bc.bottom.type = static_cast<BCType>(bottom);
    if (bc.bottom.type == BCType::Moving) {
        float sp = static_cast<float>(bc.bottom.moving);
        ImGui::SliderFloat("Bottom speed", &sp, -5.0f, 5.0f);
        bc.bottom.moving = sp;
    } else if (bc.bottom.type == BCType::Inflow) {
        float u = static_cast<float>(bc.bottom.inflow_u);
        float v = static_cast<float>(bc.bottom.inflow_v);
        ImGui::SliderFloat("Bottom U", &u, -5.0f, 5.0f);
        ImGui::SliderFloat("Bottom V", &v, -5.0f, 5.0f);
        bc.bottom.inflow_u = u;
        bc.bottom.inflow_v = v;
    }

    int top = static_cast<int>(bc.top.type);
    if (ImGui::Combo("Top", &top, bc_items, IM_ARRAYSIZE(bc_items)))
        bc.top.type = static_cast<BCType>(top);
    if (bc.top.type == BCType::Moving) {
        float sp = static_cast<float>(bc.top.moving);
        ImGui::SliderFloat("Top speed", &sp, -5.0f, 5.0f);
        bc.top.moving = sp;
    } else if (bc.top.type == BCType::Inflow) {
        float u = static_cast<float>(bc.top.inflow_u);
        float v = static_cast<float>(bc.top.inflow_v);
        ImGui::SliderFloat("Top U", &u, -5.0f, 5.0f);
        ImGui::SliderFloat("Top V", &v, -5.0f, 5.0f);
        bc.top.inflow_u = u;
        bc.top.inflow_v = v;
    }

    ImGui::Separator();
    const char *items[] = {"u", "v", "speed", "p", "vort"};
    int idx = static_cast<int>(field);
    if (ImGui::Combo("Field", &idx, items, IM_ARRAYSIZE(items)))
        field = static_cast<Field>(idx);

    ImVec2 size(512, 256);
    ImGui::Image(reinterpret_cast<void *>(static_cast<intptr_t>(texture)), size,
                 ImVec2(0, 1), ImVec2(1, 0));
    ImGui::EndTabItem();

        if (ImGui::BeginTabItem("Simulations")) {
            const char *preset_items[] = {"Jet Plume", "Lid-Driven Cavity", "Periodic Shear"};
            ImGui::Combo("Preset", &preset, preset_items, IM_ARRAYSIZE(preset_items));
            if (ImGui::Button("Apply Preset")) apply_preset = true;
            ImGui::SameLine();
            if (ImGui::Button("Reset & Run")) { apply_preset = true; reset_run = true; }

            if (preset == static_cast<int>(Preset::JetPlume)) {
                float uin = static_cast<float>(bc.left.inflow_u);
                float y0 = static_cast<float>(bc.jet_center);
                float w = static_cast<float>(bc.jet_width);
                float delta = static_cast<float>(bc.jet_thickness);
                float eps = static_cast<float>(bc.jet_eps);
                int k = bc.jet_k;
                ImGui::SliderFloat("U_in", &uin, 0.0f, 5.0f);
                ImGui::SliderFloat("y0", &y0, 0.0f, static_cast<float>(Ly));
                ImGui::SliderFloat("w", &w, 0.0f, static_cast<float>(Ly));
                ImGui::SliderFloat("delta", &delta, 0.005f * static_cast<float>(Ly), 0.05f * static_cast<float>(Ly));
                ImGui::SliderFloat("epsilon", &eps, 0.0f, 0.1f);
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
            } else if (preset == static_cast<int>(Preset::LidDrivenCavity)) {
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
    std::vector<double> tmp(static_cast<size_t>(nx) * ny);
    double vmin = 1e30, vmax = -1e30;
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int ii = i + g.ngx;
            int jj = j + g.ngy;
            double val = 0.0;
            switch (field) {
            case Gui::Field::U:
                val = 0.5 * (s.u.at_raw(ii, jj) + s.u.at_raw(ii + 1, jj));
                break;
            case Gui::Field::V:
                val = 0.5 * (s.v.at_raw(ii, jj) + s.v.at_raw(ii, jj + 1));
                break;
            case Gui::Field::Speed: {
                double uc = 0.5 * (s.u.at_raw(ii, jj) + s.u.at_raw(ii + 1, jj));
                double vc = 0.5 * (s.v.at_raw(ii, jj) + s.v.at_raw(ii, jj + 1));
                val = std::sqrt(uc * uc + vc * vc);
                break;
            }
            case Gui::Field::Pressure:
                val = s.p.at_raw(ii, jj);
                break;
            case Gui::Field::Vorticity: {
                double dv_dx = (s.v.at_raw(ii + 1, jj) - s.v.at_raw(ii - 1, jj)) /
                               (2.0 * g.dx);
                double du_dy = (s.u.at_raw(ii, jj + 1) - s.u.at_raw(ii, jj - 1)) /
                               (2.0 * g.dy);
                val = dv_dx - du_dy;
                break;
            }
            }
            tmp[i + j * nx] = val;
            vmin = std::min(vmin, val);
            vmax = std::max(vmax, val);
        }
    }
    buffer.resize(static_cast<size_t>(nx) * ny * 3);
    double scale = (vmax - vmin) < 1e-12 ? 1.0 : 1.0 / (vmax - vmin);
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
