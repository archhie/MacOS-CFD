#include "ui.hpp"
#include <imgui.h>

void Ui::render() {
    ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
    ImGui::Begin("Controls");
    ImGui::BeginDisabled();
    ImGui::Button("Run/Pause");
    ImGui::EndDisabled();
    if (ImGui::Button("Step")) { }
    ImGui::BeginDisabled();
    ImGui::Button("Reset");
    ImGui::EndDisabled();
    ImGui::Text("FPS: %.1f", ImGui::GetIO().Framerate);
    ImGui::End();
}
