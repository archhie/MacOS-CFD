# cfd2d_gui

Minimal C++20 OpenGL + Dear ImGui project.

## Build/run
```bash
brew install llvm glfw
cmake -S . -B build \
  -DCMAKE_CXX_COMPILER=$(brew --prefix llvm)/bin/clang++ \
  -DCMAKE_PREFIX_PATH=$(brew --prefix glfw)
cmake --build build -j
./build/cfd2d_gui
```
