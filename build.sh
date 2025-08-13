#!/bin/bash

# Build script for CFD2D GUI application with dockable interface
# Updated for modern macOS setup with Homebrew

set -e  # Exit on any error

echo "🚀 Building CFD2D GUI application with dockable interface..."

# Check for required tools
echo "📋 Checking prerequisites..."

if ! command -v cmake &> /dev/null; then
    echo "❌ CMake not found. Please install with: brew install cmake"
    exit 1
fi

if ! command -v brew &> /dev/null; then
    echo "❌ Homebrew not found. Please install Homebrew first."
    exit 1
fi

# Check for GLFW
if ! brew list | grep -q glfw; then
    echo "📦 Installing GLFW..."
    brew install glfw
fi

# Clean previous build
echo "🧹 Cleaning previous build..."
rm -rf build

# Create build directory
mkdir -p build

# Determine compiler setup
COMPILER_FLAGS=""
OPENMP_FLAGS=""

# Check if LLVM is available (for OpenMP support)
if brew list | grep -q llvm; then
    echo "🔧 Using LLVM compiler with OpenMP support..."
    LLVM_PREFIX=$(brew --prefix llvm)
    COMPILER_FLAGS="-DCMAKE_CXX_COMPILER=$LLVM_PREFIX/bin/clang++"
    OPENMP_FLAGS="-DOpenMP_CXX_FLAGS=-fopenmp -DOpenMP_CXX_LIB_NAMES=omp -DOpenMP_omp_LIBRARY=$LLVM_PREFIX/lib/libomp.dylib"
else
    echo "⚠️  LLVM not found. Using system compiler (OpenMP may not be available)..."
    echo "   Install LLVM for OpenMP support: brew install llvm"
fi

# Configure with CMake
echo "⚙️  Configuring with CMake..."
cmake -S . -B build \
    $COMPILER_FLAGS \
    -DCMAKE_PREFIX_PATH="$(brew --prefix glfw)" \
    $OPENMP_FLAGS \
    -DCMAKE_BUILD_TYPE=Release

# Build the application
echo "🔨 Building application..."
cmake --build build -j$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

# Check if build was successful
if [ -f "build/cfd2d_gui" ]; then
    echo ""
    echo "✅ Build completed successfully!"
    echo ""
    echo "🎮 Run the application:"
    echo "   ./build/cfd2d_gui"
    echo ""
    echo "📊 Run without GUI:"
    echo "   ./build/cfd2d_gui --no-gui"
    echo ""
    echo "❓ Show help:"
    echo "   ./build/cfd2d_gui --help"
    echo ""
    echo "🎯 Features available:"
    echo "   • Professional dockable interface"
    echo "   • Real-time CFD simulation"
    echo "   • Multiple visualization options"
    echo "   • Preset configurations"
    echo "   • Layout persistence"
else
    echo "❌ Build failed! Check the error messages above."
    exit 1
fi 