# CFD2D Skeleton

A minimal project skeleton for a 2D CFD GUI using GLFW, OpenGL and Dear ImGui.

## Build

```bash
brew install llvm glfw
cmake -S . -B build -DCMAKE_CXX_COMPILER=$(brew --prefix llvm)/bin/clang++ -DCMAKE_PREFIX_PATH=$(brew --prefix glfw) -DOpenMP_CXX_FLAGS="-fopenmp" -DOpenMP_CXX_LIB_NAMES="omp" -DOpenMP_omp_LIBRARY=$(brew --prefix llvm)/lib/libomp.dylib
cmake --build build -j
./build/cfd2d_gui
```

### Troubleshooting

OpenMP link errors: ensure lib from $(brew --prefix llvm)/lib/libomp.dylib and pass OpenMP cmake vars as above.
