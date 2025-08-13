#include "viz.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace {
std::string read_file(const char* path) {
    std::ifstream ifs(path);
    std::stringstream ss;
    ss << ifs.rdbuf();
    return ss.str();
}
}

Viz::Viz() : pixels_(width_ * height_) {
    GLuint vs = compile_shader(GL_VERTEX_SHADER, "shaders/quad.vert");
    GLuint fs = compile_shader(GL_FRAGMENT_SHADER, "shaders/scalar.frag");
    program_ = link_program(vs, fs);
    glDeleteShader(vs);
    glDeleteShader(fs);

    float quad[] = {
        -1.f, -1.f, 0.f, 0.f,
         1.f, -1.f, 1.f, 0.f,
        -1.f,  1.f, 0.f, 1.f,
         1.f,  1.f, 1.f, 1.f,
    };

    glGenVertexArrays(1, &vao_);
    glBindVertexArray(vao_);
    glGenBuffers(1, &vbo_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quad), quad, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
    glBindVertexArray(0);

    glGenTextures(1, &tex_);
    glBindTexture(GL_TEXTURE_2D, tex_);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, width_, height_, 0, GL_RED, GL_UNSIGNED_BYTE, pixels_.data());
    glBindTexture(GL_TEXTURE_2D, 0);
}

Viz::~Viz() {
    if (tex_) glDeleteTextures(1, &tex_);
    if (vbo_) glDeleteBuffers(1, &vbo_);
    if (vao_) glDeleteVertexArrays(1, &vao_);
    if (program_) glDeleteProgram(program_);
}

void Viz::update() {
    static int frame = 0;
    for (int y = 0; y < height_; ++y) {
        for (int x = 0; x < width_; ++x) {
            pixels_[y * width_ + x] = static_cast<unsigned char>((x + frame) % 256);
        }
    }
    ++frame;
    glBindTexture(GL_TEXTURE_2D, tex_);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width_, height_, GL_RED, GL_UNSIGNED_BYTE, pixels_.data());
    glBindTexture(GL_TEXTURE_2D, 0);
}

void Viz::render() {
    glUseProgram(program_);
    glBindVertexArray(vao_);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex_);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
    glBindVertexArray(0);
    glUseProgram(0);
}

GLuint Viz::compile_shader(GLenum type, const char* path) {
    std::string src = read_file(path);
    GLuint shader = glCreateShader(type);
    const char* csrc = src.c_str();
    glShaderSource(shader, 1, &csrc, nullptr);
    glCompileShader(shader);
    GLint ok = 0;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        char log[512];
        glGetShaderInfoLog(shader, 512, nullptr, log);
        throw std::runtime_error(log);
    }
    return shader;
}

GLuint Viz::link_program(GLuint vs, GLuint fs) {
    GLuint prog = glCreateProgram();
    glAttachShader(prog, vs);
    glAttachShader(prog, fs);
    glLinkProgram(prog);
    GLint ok = 0;
    glGetProgramiv(prog, GL_LINK_STATUS, &ok);
    if (!ok) {
        char log[512];
        glGetProgramInfoLog(prog, 512, nullptr, log);
        throw std::runtime_error(log);
    }
    return prog;
}
