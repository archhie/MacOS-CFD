#pragma once
#include <glad/glad.h>
#include <vector>

class Viz {
public:
    Viz();
    ~Viz();
    void update();
    void render();
private:
    GLuint program_{};
    GLuint vao_{};
    GLuint vbo_{};
    GLuint tex_{};
    int width_ = 256;
    int height_ = 256;
    std::vector<unsigned char> pixels_;
    GLuint compile_shader(GLenum type, const char* path);
    GLuint link_program(GLuint vs, GLuint fs);
};
