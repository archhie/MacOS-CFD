cmake_minimum_required(VERSION 3.21)
project(cfd2d_gui LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(imgui
    third_party/imgui/imgui.cpp
    third_party/imgui/imgui_draw.cpp
    third_party/imgui/imgui_widgets.cpp
    third_party/imgui/imgui_tables.cpp
    third_party/imgui/backends/imgui_impl_glfw.cpp
    third_party/imgui/backends/imgui_impl_opengl3.cpp
)

target_include_directories(imgui PUBLIC
    third_party/imgui
    third_party/imgui/backends
)

file(GLOB SRC_FILES
    src/*.cpp
)

set(SOLVER_FILES
    src/solver/grid.cpp
    src/solver/bc.cpp
    src/solver/advection.cpp
    src/solver/diffusion.cpp
    src/solver/metrics.cpp
    src/solver/pressure.cpp
    src/solver/project.cpp
    src/solver/time_integrator.cpp
)

add_executable(cfd2d_gui
    ${SRC_FILES}
    ${SOLVER_FILES}
)

target_include_directories(cfd2d_gui PUBLIC
    src
    third_party/imgui
    third_party/imgui/backends
    third_party
)

find_package(OpenGL REQUIRED)
find_package(glfw3 3.3 REQUIRED CONFIG)

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    target_link_libraries(cfd2d_gui PRIVATE OpenMP::OpenMP_CXX)
endif()

set_target_properties(cfd2d_gui PROPERTIES
    OUTPUT_NAME cfd2d_gui
)

target_link_libraries(cfd2d_gui PRIVATE imgui glfw OpenGL::GL)

target_compile_options(cfd2d_gui PRIVATE -O3 -Wall -Wextra -pedantic)
