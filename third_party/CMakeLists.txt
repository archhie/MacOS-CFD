add_library(glad STATIC
    glad/src/glad.c
)
target_include_directories(glad PUBLIC glad/include)

add_library(imgui STATIC
    imgui/imgui.cpp
    imgui/imgui_draw.cpp
    imgui/imgui_tables.cpp
    imgui/imgui_widgets.cpp
    imgui/backends/imgui_impl_glfw.cpp
    imgui/backends/imgui_impl_opengl3.cpp
)

target_include_directories(imgui PUBLIC
    imgui
    imgui/backends
    glad/include
)

target_compile_definitions(imgui PUBLIC IMGUI_IMPL_OPENGL_LOADER_GLAD=1 GLFW_INCLUDE_NONE)
