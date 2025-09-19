# MoFE

**M**odern **F**inite **E**lement library - A C++ project that interfaces with MFEM to solve various engineering problems.

## Overview

MoFE is a C++ library that provides a modern interface to the MFEM (Modular Finite Element Methods) library. It demonstrates how to create clean, object-oriented interfaces for solving finite element problems, starting with computational fluid dynamics (CFD) problems.

## Features

- **CFD Solver Class**: A complete class that demonstrates interfacing with MFEM's core components:
  - `Mesh`: Computational mesh handling
  - `Element`: Finite element operations
  - `FiniteElementSpace`: Finite element space definition
  - `GridFunction`: Solution representation
  - `BilinearForm`: System matrix assembly
  - `LinearForm`: Load vector assembly

## Project Structure

```
MoFE/
├── include/           # Header files
│   └── cfd_solver.hpp # CFD solver class interface
├── src/               # Implementation files  
│   └── cfd_solver.cpp # CFD solver implementation
├── examples/          # Example programs
│   └── cfd_example.cpp # Demonstrates CFD solver usage
├── data/              # Data files (automatically copied to build directory)
│   ├── README.md      # Data directory documentation
│   ├── config.txt     # Sample configuration file
│   └── test.mesh      # Example mesh file
├── mfem/              # MFEM submodule
├── test.mesh          # Example mesh file
└── CMakeLists.txt     # Build configuration
```

## Building

Prerequisites:
- CMake 3.10+
- C++17 compatible compiler
- Git (for submodules)
- Optional: CUDA toolkit (for GPU acceleration)
- Optional: MPI implementation (for parallel computing)
- Optional: AMGX library (for advanced linear solvers)

```bash
# Clone repository with submodules
git clone --recursive https://github.com/Ruochun/MoFE.git
cd MoFE

# Basic build
mkdir build && cd build
cmake ..
make -j4
```

## CMake Configuration Options

MoFE provides several user-configurable options that can be set during CMake configuration:

### Build Type Configuration
The build type can be selected from multiple options and is easily configurable via ccmake GUI:
- `Debug` (default): Debug build with debug symbols
- `Release`: Optimized release build
- `RelWithDebInfo`: Release build with debug information
- `MinSizeRel`: Size-optimized release build

```bash
# Set build type via command line
cmake -DCMAKE_BUILD_TYPE=Release ..

# Or use ccmake for interactive configuration
ccmake .
```

### MFEM Integration Options

**CUDA Support (Automatic)**
- CUDA support is automatically enabled when a CUDA toolkit is detected
- When CUDA toolkit is not available, support is automatically disabled with a warning
- This ensures maximum GPU acceleration when possible while maintaining compatibility

**User-Configurable Options**
- `MOFE_USE_MPI` (default: OFF): Enable MPI parallel computing support
- `MOFE_USE_AMGX` (default: OFF): Enable AMGX advanced linear solver support

```bash
# Enable MPI support
cmake -DMOFE_USE_MPI=ON ..

# Enable both MPI and AMGX
cmake -DMOFE_USE_MPI=ON -DMOFE_USE_AMGX=ON ..

# Interactive configuration with ccmake (recommended)
ccmake .
```

### Data Directory
The `data/` directory is automatically copied to the build directory during compilation, ensuring data files are available to executables. This includes:
- Example mesh files
- Configuration files  
- Test data

## Running Examples

```bash
# Run CFD example with the provided test mesh (from build directory)
./cfd_example -m data/test.mesh -o 1 -out solution.vtk

# Available options:
# -m <mesh_file>    : Mesh file to use
# -o <order>        : Finite element order
# -out <output>     : Output VTK file
```

## Example: CFD Solver Usage

```cpp
#include "cfd_solver.hpp"

int main() {
    // Create CFD solver
    mofe::CFDSolver solver("mesh.mesh", 1);
    
    // Initialize (load mesh, setup FE spaces)
    solver.Initialize();
    
    // Solve the problem
    solver.Solve();
    
    // Save solution
    solver.SaveSolution("output.vtk");
    
    return 0;
}
```

## What the CFD Solver Demonstrates

The `CFDSolver` class shows how to properly interface with MFEM's key classes:

1. **Mesh Management**: Loading and working with computational meshes
2. **Finite Element Spaces**: Setting up H1 finite element collections and spaces
3. **System Assembly**: Using BilinearForm and LinearForm for matrix/vector assembly
4. **Boundary Conditions**: Applying essential (Dirichlet) boundary conditions
5. **Linear Solving**: Using MFEM's built-in linear solvers (PCG with GS smoother)
6. **Solution Output**: Saving results in VTK format for visualization

## Future Extensions

This project serves as a foundation for adding more finite element problem solvers:
- Structural mechanics (elasticity problems)
- Heat transfer
- Electromagnetics
- Multi-physics coupling

Each new solver will demonstrate different aspects of MFEM's capabilities while maintaining clean, modern C++ interfaces.

## License

BSD 3-Clause License - See LICENSE file for details.