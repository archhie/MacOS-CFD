# CFD2D GUI - 2D Computational Fluid Dynamics with Dockable Interface

A modern, interactive 2D CFD simulation application with a professional dockable interface built using ImGui and OpenGL.

## Features

### 🎯 Professional Dockable Interface
- **Fullscreen DockSpace**: Modern docking system similar to professional IDEs
- **Dockable Panels**: Resizable, movable, and collapsible panels
- **Layout Persistence**: Remembers your layout between sessions
- **Retina Display Support**: Crisp rendering on macOS/Retina displays

### 📊 Simulation Panels
- **CFD Controls**: Flow parameters, timestep control, simulation speed
- **Visualization**: Field selection, colormaps, rendering options, statistics
- **Boundary Conditions**: Presets, ASCII preview, custom BC settings
- **Simulation Viewport**: Real-time simulation display with overlay controls

### 🔧 CFD Features
- **Multiple Solvers**: PCG and Multigrid pressure solvers
- **Preset Configurations**: Jet Plume, Lid-Driven Cavity, Periodic Shear
- **Real-time Visualization**: Multiple field types and colormaps
- **Performance Monitoring**: FPS, convergence metrics, field statistics

## Building

### Prerequisites
- CMake 3.21+
- C++20 compiler
- OpenGL 3.3+
- GLFW 3.3+
- OpenMP (optional, for parallel processing)

### Build Instructions
```bash
mkdir build && cd build
cmake ..
make -j4
```

### Running
```bash
./cfd2d_gui
```

## Usage

### Interface Overview
The application launches with a professional dockable interface:

1. **Top Toolbar**: Black menu bar with panel toggles and Reset Layout
2. **Central Viewport**: Simulation display with floating controls
3. **Left Panels**: CFD Controls and Visualization panels
4. **Right Panel**: Boundary Conditions panel

### Panel Controls

#### CFD Controls Panel
- **Flow Parameters**: Adjust Reynolds number and CFL condition
- **Timestep Control**: Override automatic timestep calculation
- **Simulation Speed**: Control simulation speed multiplier
- **Playback Controls**: Play/Pause, Step, Reset buttons
- **Simulation Info**: Real-time metrics and performance stats

#### Visualization Panel
- **Field Selection**: Choose from U, V, Speed, Pressure, Vorticity, Streamlines, Temperature
- **Colormaps**: Viridis, Plasma, Turbo, Gray, RdBu, Coolwarm
- **Rendering Options**: Vector arrows, contours, streamlines, normalization
- **Color Scale**: Auto-scale or manual min/max with gamma correction
- **Field Statistics**: Min, max, mean, standard deviation

#### Boundary Conditions Panel
- **Preset Configurations**: Quick setup for common scenarios
- **ASCII Preview**: Visual boundary condition representation
- **Custom Settings**: Fine-tune each boundary (left, right, top, bottom)
- **Jet Settings**: Advanced jet plume configuration

### Preset Configurations

#### Jet Plume
- Left: Inflow with jet profile
- Right: Outflow
- Top/Bottom: Walls
- Reynolds: 50, CFL: 0.01

#### Lid-Driven Cavity
- All boundaries: Walls
- Top: Moving lid
- Reynolds: 1000, CFL: 0.1

#### Periodic Shear
- Left/Right: Periodic boundaries
- Top/Bottom: Walls
- Reynolds: 1000, CFL: 0.1

### Keyboard Shortcuts
- **Space**: Play/Pause simulation
- **R**: Reset simulation
- **S**: Single step
- **Escape**: Close application

## Technical Details

### Architecture
- **ImGui Docking Branch**: Modern docking system
- **OpenGL 3.3**: Hardware-accelerated rendering
- **GLFW**: Cross-platform window management
- **C++20**: Modern C++ features

### Performance
- **60 FPS GUI**: Smooth interface responsiveness
- **Separate Simulation Loop**: Independent simulation and GUI frame rates
- **OpenMP Support**: Parallel processing for large grids
- **Memory Efficient**: Optimized data structures and rendering

### File Structure
```
src/
├── main.cpp              # Application entry point
├── gui/
│   ├── gui.hpp          # GUI class declaration
│   └── gui.cpp          # GUI implementation
├── solver/              # CFD solver components
├── viz.hpp              # Visualization utilities
└── colormaps.hpp        # Color mapping functions
```

## Configuration

### Command Line Options
```bash
./cfd2d_gui [options]
  --no-gui              Run in headless mode
  --Re=<value>          Set Reynolds number (default: 50)
  --CFL=<value>         Set CFL condition (default: 0.01)
  --solver=<pcg|mg>     Pressure solver type (default: pcg)
  --mg-levels=<n>       Multigrid levels (default: 3)
  --vcycles=<n>         Multigrid V-cycles (default: 2)
  --tol=<value>         Solver tolerance (default: 1e-8)
  --pcg-iters=<n>       PCG max iterations (default: 1000)
  --tvd-debug           Enable TVD debugging output
  --help                Show help message
```

### Layout Persistence
The application automatically saves your layout preferences to `imgui.ini` in the working directory. This includes:
- Panel positions and sizes
- Docking layout
- Window states (collapsed/expanded)
- UI preferences

## Troubleshooting

### Common Issues
1. **Missing Font**: If monospace font not found, falls back to default
2. **Layout Reset**: Use "Reset Layout" in top toolbar to restore default layout
3. **Performance**: Reduce grid size or disable features for better performance
4. **Visualization**: Try different colormaps or auto-scale for better field visibility

### Debugging
- Enable `--tvd-debug` for detailed solver output
- Check console for convergence and performance metrics
- Monitor FPS in CFD Controls panel

## Development

### Adding New Features
1. **New Fields**: Add to `Field` enum in `gui.hpp`
2. **New Colormaps**: Add to `Colormap` enum and implement in `colormaps.hpp`
3. **New Presets**: Add to `Preset` enum and implement in main loop
4. **New Panels**: Create new panel function and add to main draw loop

### Code Style
- Follow existing C++20 patterns
- Use ImGui best practices for docking
- Maintain 60 FPS GUI responsiveness
- Document new features

## License

This project uses ImGui which is licensed under the MIT License.
