#include "cfd_solver.hpp"
#include <iostream>
#include <string>

int main(int argc, char *argv[]) {
    std::cout << "=== MoFE CFD Solver Example ===" << std::endl;
    
    // Initialize MFEM (required for MFEM applications)
    mfem::OptionsParser args(argc, argv);
    
    // Default parameters
    std::string mesh_file = "square.mesh";
    int order = 1;
    std::string output_file = "solution.vtk";
    
    // Parse command line arguments
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&order, "-o", "--order",
                   "Finite element order (polynomial degree).");
    args.AddOption(&output_file, "-out", "--output",
                   "Output file for solution.");
    
    args.Parse();
    if (!args.Good()) {
        args.PrintUsage(std::cout);
        return 1;
    }
    args.PrintOptions(std::cout);
    
    try {
        // Create and initialize CFD solver
        mofe::CFDSolver solver(mesh_file, order);
        
        // Initialize the solver (sets up mesh, finite element spaces, etc.)
        solver.Initialize();
        
        // Solve the CFD problem
        solver.Solve();
        
        // Save solution
        solver.SaveSolution(output_file);
        
        // Print some information about the solution
        std::cout << "\n=== Solution Information ===" << std::endl;
        std::cout << "Mesh dimension: " << solver.GetMesh()->Dimension() << std::endl;
        std::cout << "Number of elements: " << solver.GetMesh()->GetNE() << std::endl;
        std::cout << "Number of DOFs: " << solver.GetFESpace()->GetTrueVSize() << std::endl;
        
        // Compute and print solution statistics
        mfem::GridFunction* sol = solver.GetSolution();
        double min_val = sol->Min();
        double max_val = sol->Max();
        std::cout << "Solution range: [" << min_val << ", " << max_val << "]" << std::endl;
        
        std::cout << "\nCFD simulation completed successfully!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}