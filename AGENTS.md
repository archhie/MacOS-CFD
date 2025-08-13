TITLE: Modular Agents Plan — Interactive 2D CFD GUI for macOS (C++20 + LLVM/OpenMP)

OVERVIEW
We will build an interactive, real‑time 2D incompressible CFD simulator with a GUI, multiple user‑drawn obstacles, and variable wind inflow. The work is split into sequential, self‑contained AGENT tasks. Each agent outputs compilable code for its scope, with a README snippet and acceptance checks.

GLOBAL TECH & BUILD (FOR ALL AGENTS)

Language: C++20

Window/UI: GLFW + Dear ImGui (vendored in-tree)

Graphics: OpenGL 3.3 core; full-screen textured quad; GLSL for colormap

Parallelism: OpenMP (Homebrew libomp)

Build: CMake; compiler is Homebrew LLVM clang++

Target binary: cfd2d_gui

Minimum macOS deps: brew install llvm glfw

Shared build commands (include verbatim in each agent’s README block):
brew install llvm glfw
cmake -S . -B build -DCMAKE_CXX_COMPILER=$(brew --prefix llvm)/bin/clang++ -DCMAKE_PREFIX_PATH=$(brew --prefix glfw) -DOpenMP_CXX_FLAGS="-fopenmp" -DOpenMP_CXX_LIB_NAMES="omp" -DOpenMP_omp_LIBRARY=$(brew --prefix llvm)/lib/libomp.dylib
cmake --build build -j
./build/cfd2d_gui

PHYSICS & NUMERICS (BASE SPEC SHARED BY AGENTS)

Equations (nondimensional): ∂u/∂t + (u·∇)u = −∇p + (1/Re)∇²u, with ∇·u = 0

Grid: 2D MAC (staggered), 1–2 ghost layers
• u on x-faces (nx+1, ny), v on y-faces (nx, ny+1), p at centers (nx, ny)

Spatial: MUSCL‑TVD (minmod) advection; 2nd‑order central diffusion

Time stepping: RK2 (Heun) on adv+diff, then Chorin projection

Pressure Poisson: ∇²p = (1/Δt)∇·u*; default geometric multigrid (V‑cycle: full‑weight restrict, bilinear prolong‑add, red‑black Gauss‑Seidel smooth); fallback PCG (diagonal precond)

Adaptive Δt via CFL: Δt ≤ CFL * min(dx/|u|max, dy/|v|max, 0.5·Re·min(dx²,dy²))

BCs per side: Wall (no‑slip), Moving wall, Inflow (Dirichlet u,v), Outflow (Neumann u,v; zero‑grad p), Periodic

Obstacles: user shapes → rasterized solid mask; enforce u=v=0 in solids; clamp faces crossing solid boundaries to 0; adjust stencils to skip solid neighbors

Default startup: Lx=2, Ly=1; Nx=512, Ny=256; Re=2000; CFL=0.5; left inflow slot at y∈[0.4,0.6], U_in=1.0; walls top/bottom; outflow right; run immediately

Visualization: scalar = Speed | Vorticity | Pressure; colormaps = Turbo | Viridis | Gray; auto or manual range

Performance: -O3 -march=native; OpenMP; no per‑step allocations; profiled hot loops; multigrid buffers persistent

PROJECT LAYOUT (TARGET STRUCTURE)

CMakeLists.txt

README.md

shaders/  → quad.vert, scalar.frag

third_party/ → Dear ImGui (imgui*.cpp, backends imgui_impl_glfw.cpp/imgui_impl_opengl3.cpp); stb_image_write.h

src/
• main.cpp
• viz.hpp/.cpp, colormaps.hpp
• ui.hpp/.cpp, profiler.hpp/.cpp
• shapes.hpp/.cpp
• solver/
grid.hpp/.cpp, field.hpp, state.hpp, bc.hpp/.cpp,
advection.hpp/.cpp, diffusion.hpp/.cpp, project.hpp/.cpp,
pressure.hpp/.cpp, metrics.hpp/.cpp, stepper.hpp/.cpp

ACCEPTANCE TESTS (GLOBAL)

Builds and runs with the brew/cmake commands above on Apple Silicon or Intel

Startup shows running plume; scalar view toggles; FPS + metrics visible

Divergence L2 small after projection (<1e−6 typical on moderate grids)

Shapes influence the flow; deleting shapes restores baseline

Variable wind parameters clearly affect jet and wake; screenshot saves PNG to out/

──────────────────────────────────────────────────────────────────────────────
AGENT 1 — PROJECT SKELETON & BUILD SYSTEM
SCOPE

Create CMake project producing cfd2d_gui with GLFW, OpenGL, Dear ImGui vendored

Provide shaders/quad.vert & shaders/scalar.frag; basic GL loader; RAII wrappers

Minimal window (1280×720) titled "CFD2D"; ImGui panel with FPS; draw a test texture

DELIVERABLES

CMakeLists.txt with C++20, LLVM/clang++ via brew, OpenMP flags present but optional

third_party/ with Dear ImGui (core + glfw/opengl3 backends), stb_image_write.h

src/main.cpp creating window, GL context, ImGui; render loop draws sample quad

README.md with exact brew/cmake build commands; troubleshooting note for libomp

ACCEPTANCE

App opens a window, shows ImGui panel with FPS; draws colored quad; ESC quits; builds cleanly with -O3 -Wall

──────────────────────────────────────────────────────────────────────────────
AGENT 2 — CORE DATA (GRID, FIELDS, STATE, BC STUBS)
SCOPE

Implement MAC Grid (nx, ny, Lx, Ly, dx, dy, ngx, ngy) and index helpers for p/u/v

Implement Field2D with 64‑byte alignment, pitch padding, ghost layers; RAII

Implement State holding u,v,p, rhs,tmp, advu,advv, Lu,Lv, u1,v1, scalar buffer

Implement BC types enum and stub functions (no physics yet)

Wire allocation in main; print grid info

ACCEPTANCE

Runs; allocates/frees fields; prints grid/dx/dy; no simulation yet

──────────────────────────────────────────────────────────────────────────────
AGENT 3 — ADVECTION, DIFFUSION, CFL (NO PRESSURE)
SCOPE

advection.hpp/.cpp: MUSCL‑TVD (minmod) for u and v on MAC; upwind with local face sign; branch‑free interior loops; separate BC passes

diffusion.hpp/.cpp: 2nd‑order central Laplacian for u and v

metrics.hpp/.cpp: max_velocity(), norms; CFL timestep suggestion

stepper.hpp/.cpp: rk2_advdiff_step(Grid, State, Params, dt_out)

UI: controls for Re, CFL, dt cap; a Step button; print chosen dt

ACCEPTANCE

Clicking Step advances time; dt prints; values finite; app stable

──────────────────────────────────────────────────────────────────────────────
AGENT 4 — PRESSURE PROJECTION (MG + PCG)
SCOPE

project.hpp/.cpp: divergence(u,v) → rhs; subtract_grad_p(p)

pressure.hpp/.cpp: geometric multigrid V‑cycle (full‑weight restrict, bilinear prolong‑add, RBGS smoother); PCG fallback (diag precond)

stepper.hpp/.cpp: project(dt, solver_params) after RK2 adv+diff; add Run/Pause and substep loop to target frame budget

UI: solver choice (MG/PCG), mg_levels, vcycles, tol, pcg_iters; show divergence L2 and MG residual

ACCEPTANCE

With quiescent initial state, projection drives divergence near 0; residual decreases per V‑cycle

──────────────────────────────────────────────────────────────────────────────
AGENT 5 — VISUALIZATION (SCALARS + COLORMAPS)
SCOPE

compute_scalar(mode) where mode ∈ {Speed, Vorticity, Pressure}; speed via face→center interpolation; vorticity ωz ≈ ∂v/∂x − ∂u/∂y

GLSL shader + 1D LUT textures for Turbo, Viridis, Gray (colormaps.hpp)

viz.cpp uploads scalar to GL texture each frame and draws full‑screen quad

UI: scalar radio, colormap combo, Auto range toggle + manual min/max

ACCEPTANCE

Switching scalar/colormap updates colors live; performance sustained

──────────────────────────────────────────────────────────────────────────────
AGENT 6 — BOUNDARY CONDITIONS & DEFAULT JET INFLOW
SCOPE

Implement BCs (Wall, Moving, Inflow, Outflow, Periodic) for u,v,p ghosts; keep hot loops separate from BC passes

Default case: Lx=2, Ly=1; left inflow slot y∈[0.4,0.6] with u=U_in, v=0; walls elsewhere; right outflow; top/bottom walls

UI: controls U_in, slot center/width; per‑side BC selector combos

Start running automatically; show plume within a few seconds on 512×256

ACCEPTANCE

Jet core visible; Vorticity reveals shear rolls; Pressure smooth; stable stepping under CFL

──────────────────────────────────────────────────────────────────────────────
AGENT 7 — SHAPE EDITOR & SOLID MASK COUPLING
SCOPE

Shape toolbox: Select | Rectangle | Circle | Polygon | Eraser | Import PNG/SVG (subset)

Shapes stored with transform/edit gizmos (when paused); draw outlines overlay

Rasterize to solid_mask (uint8 at centers). For polygon, basic scanline fill; circle by radius test; rect by bounds

Solver coupling: inside solids, clamp u=v=0; faces crossing solid boundary set to 0; adjust stencils to skip solid neighbor values

UI: shape list (visibility, delete), Bake Mask button; auto‑rebake on edit if paused

ACCEPTANCE

Placing a rectangle in the jet causes a wake; deleting it restores baseline; stability preserved

──────────────────────────────────────────────────────────────────────────────
AGENT 8 — VARIABLE WIND (SPATIAL + TEMPORAL) & POLISH
SCOPE

Spatial profiles: Top‑hat (slot center/width), Parabolic, Cosine, Custom (8 control points via ImGui curve sampled over y)

Temporal modes: Constant, Sinusoidal (amp, freq), Gusts (randomized amp with seed, update period)

Preview plots: profile(y) and U_in(t)

Keyboard: Space=Run/Pause, N=Step, R=Reset, S=Screenshot (PNG to out/)

Profiler HUD: per‑kernel ms (advect, diffuse, project) and FPS; dt, CFL, max|u|

README polish: features, controls, troubleshooting

ACCEPTANCE

Changing inflow params clearly alters plume; screenshots saved; HUD shows timings; divergence remains small

──────────────────────────────────────────────────────────────────────────────
SHARED CODE & PSEUDOCODE SNIPPETS (FOR REUSE BY AGENTS)

minmod(a,b) = (a*b<=0)?0: (|a|<|b|?a:b)

RK2(Heun): y* = y^n + dt·f(y^n); y^{n+1} = 0.5·( y^n + y* + dt·f(y*) )

Projection: rhs=(1/dt)·div(u**); solve Laplace(p)=rhs; u←u**−dt·grad p

MG V‑cycle: pre‑smooth k times, compute residual, restrict to coarse, solve or smooth more, prolong‑add, post‑smooth

TROUBLESHOOTING (INCLUDE IN README)

OpenMP link errors: ensure lib from $(brew --prefix llvm)/lib/libomp.dylib and pass OpenMP cmake vars as above

Shaders not found: copy shaders/ next to binary or embed GLSL as strings

Performance low: reduce grid (e.g., 256×128), limit max_substeps, tune mg_levels/vcycles

Cursor/Codex token limits: run agents in order; each agent must deliver compilable code and a minimal README update

LICENSE & GIT

Default license: MIT (include LICENSE file)

.gitignore suggestions: /build/, /out/, /cmake-build-*/, *.DS_Store

END OF PLAN

