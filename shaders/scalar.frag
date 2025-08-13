#version 330 core
in vec2 uv;
out vec4 color;
uniform sampler2D tex;
void main() {
    float v = texture(tex, uv).r;
    color = vec4(v, v, v, 1.0);
}
