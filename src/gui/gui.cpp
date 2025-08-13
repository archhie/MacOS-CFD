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

void Gui::draw(int timestep, double sim_time, double max_velocity,
               double pressure_residual, GLuint texture) {
    ImGui::Begin("CFD Controls");
    
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
    ImGui::Text("Max velocity: %.3f", max_velocity);
    ImGui::Text("Pressure residual: %.3e", pressure_residual);

    const char *items[] = {"u", "v", "speed", "p"};
    int idx = static_cast<int>(field);
    if (ImGui::Combo("Field", &idx, items, IM_ARRAYSIZE(items)))
        field = static_cast<Field>(idx);

    ImVec2 size(512, 256);
    ImGui::Image(reinterpret_cast<void *>(static_cast<intptr_t>(texture)), size,
                 ImVec2(0, 1), ImVec2(1, 0));
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
