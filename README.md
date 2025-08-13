# CFD2D GUI

A 2D Computational Fluid Dynamics (CFD) application with a graphical user interface using GLFW, OpenGL, and Dear ImGui.

## Features

- 2D incompressible fluid simulation
- Real-time visualization with GUI controls
- Support for different Reynolds numbers and CFL conditions
- Pressure solver with PCG and multigrid options
- VTK output for post-processing
- Command-line mode for batch processing

## Prerequisites

- macOS with Homebrew
- LLVM compiler (for OpenMP support)
- GLFW library

## Installation

1. Install required dependencies:
```bash
brew install llvm glfw
```

2. Clone the repository:
```bash
git clone <repository-url>
cd MacOS-CFD
```

## Building

### Option 1: Using the build script (Recommended)
```bash
./build.sh
```

### Option 2: Manual build
```bash
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
```

## Running

### GUI Mode
```bash
./build/cfd2d_gui
```

### Command-line Mode
```bash
./build/cfd2d_gui --no-gui --Re=100 --CFL=0.5
```

## Configuration

The application supports various command-line options:
- `--Re=<value>`: Set Reynolds number
- `--CFL=<value>`: Set CFL condition
- `--no-gui`: Run in command-line mode
- `--help`: Show all available options

## Output

Simulation results are saved in the `out/` directory in VTK format for post-processing with tools like ParaView.

## Troubleshooting

If you encounter build issues:
1. Ensure all dependencies are installed: `brew install llvm glfw`
2. Clean the build directory: `rm -rf build`
3. Re-run the build script: `./build.sh`

The CMake configuration automatically detects GLFW and other dependencies, making the build process portable across different macOS systems.
