#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>
#include <OpenGL/gl3.h>
#include "imgui.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"
#include "solver/bc.hpp"
#include "solver/grid.hpp"
#include "solver/state.hpp"
#include "solver/time_integrator.hpp"
#include "solver/pressure.hpp"

static std::string load_file(const char* path) {
    std::ifstream f(path);
    if (!f.is_open()) {
        return {};
    }
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

static GLuint compile_shader(GLenum type, const char* src) {
    GLuint s = glCreateShader(type);
    glShaderSource(s, 1, &src, nullptr);
    glCompileShader(s);
    return s;
}

static GLuint build_program(const char* vpath, const char* fpath) {
    std::string vsrc = load_file(vpath);
    std::string fsrc = load_file(fpath);
    GLuint vs = compile_shader(GL_VERTEX_SHADER, vsrc.c_str());
    GLuint fs = compile_shader(GL_FRAGMENT_SHADER, fsrc.c_str());
    GLuint prog = glCreateProgram();
    glAttachShader(prog, vs);
    glAttachShader(prog, fs);
    glLinkProgram(prog);
    glDeleteShader(vs);
    glDeleteShader(fs);
    return prog;
}

int main(int argc, char** argv) {
    if (!glfwInit()) return 1;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GLFW_TRUE);

    GLFWwindow* window = glfwCreateWindow(1280, 720, "CFD2D", nullptr, nullptr);
    if (!window) {
        glfwTerminate();
        return 1;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    (void)io;
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");

    GLuint prog = build_program("shaders/quad.vert", "shaders/scalar.frag");

    Grid grid;
    grid.init(512, 256, 2.0, 1.0, 1);
    State state;
    state.u.allocate(grid.u_nx(), grid.ny, grid.u_pitch(), grid.ngx, grid.ngy);
    state.v.allocate(grid.nx, grid.v_ny(), grid.v_pitch(), grid.ngx, grid.ngy);
    state.p.allocate(grid.nx, grid.ny, grid.p_pitch(), grid.ngx, grid.ngy);
    state.rhs.allocate(grid.nx, grid.ny, grid.p_pitch(), grid.ngx, grid.ngy);
    state.tmp.allocate(grid.nx, grid.ny, grid.p_pitch(), grid.ngx, grid.ngy);
    state.scalar.allocate(grid.nx, grid.ny, grid.p_pitch(), grid.ngx, grid.ngy);

    double Re = 1000.0;
    double CFL_target = 0.5;
    PressureParams pp;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.rfind("--Re=", 0) == 0) {
            Re = std::stod(arg.substr(5));
        } else if (arg.rfind("--CFL=", 0) == 0) {
            CFL_target = std::stod(arg.substr(6));
        } else if (arg.rfind("--solver=", 0) == 0) {
            std::string val = arg.substr(9);
            if (val == "pcg")
                pp.type = PressureSolverType::PCG;
            else
                pp.type = PressureSolverType::Multigrid;
        } else if (arg.rfind("--mg-levels=", 0) == 0) {
            pp.mg_levels = std::stoi(arg.substr(12));
        } else if (arg.rfind("--vcycles=", 0) == 0) {
            pp.vcycles = std::stoi(arg.substr(10));
        } else if (arg.rfind("--tol=", 0) == 0) {
            pp.tol = std::stod(arg.substr(6));
        } else if (arg.rfind("--pcg-iters=", 0) == 0) {
            pp.pcg_max_iters = std::stoi(arg.substr(12));
        }
    }

    PressureSolver pressure(grid);
    TimeIntegrator integrator(grid);
    BC bc;
    bc.top = BCType::Moving;
    bc.movingU = 1.0;

    printf("Grid: %d x %d (dx=%g, dy=%g)\n", grid.nx, grid.ny, grid.dx, grid.dy);

    float verts[] = {
        -1.0f, -1.0f, 0.0f, 0.0f,
         1.0f, -1.0f, 1.0f, 0.0f,
        -1.0f,  1.0f, 0.0f, 1.0f,
         1.0f,  1.0f, 1.0f, 1.0f,
    };
    GLuint vao, vbo;
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(verts), verts, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
    glEnableVertexAttribArray(1);

    int frame = 0;
    double sim_time = 0.0;

    while (!glfwWindowShouldClose(window)) {
        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
            glfwSetWindowShouldClose(window, 1);

        double dt = integrator.step(state, bc, Re, CFL_target, pressure, pp);
        sim_time += dt;
        if ((frame++ % 60) == 0) {
            printf("dt: %g time: %g\n", dt, sim_time);
        }

        glfwPollEvents();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        ImGui::Begin("Stats");
        ImGui::Text("FPS: %.1f", io.Framerate);
        ImGui::Text("Time: %.3f", sim_time);
        ImGui::End();
        ImGui::Render();

        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        glViewport(0, 0, w, h);
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        glUseProgram(prog);
        glBindVertexArray(vao);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);
    glDeleteProgram(prog);

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
