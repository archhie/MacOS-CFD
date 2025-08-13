#!/bin/bash

# Build script for CFD2D GUI application

set -e  # Exit on any error

echo "Building CFD2D GUI application..."

# Clean previous build
rm -rf build

# Configure with CMake
cmake -S . -B build \
    -DCMAKE_CXX_COMPILER=$(brew --prefix llvm)/bin/clang++ \
    -DCMAKE_PREFIX_PATH="$(brew --prefix glfw);$(brew --prefix llvm)" \
    -DOpenMP_CXX_FLAGS="-fopenmp" \
    -DOpenMP_CXX_LIB_NAMES="omp" \
    -DOpenMP_omp_LIBRARY=$(brew --prefix llvm)/lib/libomp.dylib

# Build the application
cmake --build build -j

echo "Build completed successfully!"
echo "Run the application with: ./build/cfd2d_gui"
echo "Or run without GUI: ./build/cfd2d_gui --no-gui" 