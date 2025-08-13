cmake_minimum_required(VERSION 3.16)
project(cfd2d_gui LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(glfw3 CONFIG REQUIRED)
find_package(OpenGL REQUIRED)

add_subdirectory(third_party)
target_link_libraries(imgui PUBLIC glfw)

add_executable(cfd2d_gui
    src/main.cpp
    src/viz.cpp
    src/ui.cpp
)

target_include_directories(cfd2d_gui PRIVATE src)

target_link_libraries(cfd2d_gui PRIVATE
    imgui
    glad
    glfw
    OpenGL::GL
)

target_compile_options(cfd2d_gui PRIVATE -Wall -O3)
target_compile_options(imgui PRIVATE -Wall -O3)
target_compile_options(glad PRIVATE -Wall -O3)
