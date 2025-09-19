#ifndef CFD_SOLVER_HPP
#define CFD_SOLVER_HPP

#include "mfem.hpp"

namespace mofe {

/**
 * @brief CFD Solver class that interfaces with MFEM classes for solving CFD problems
 * 
 * This class demonstrates how to interface with MFEM's core classes:
 * - Mesh: for computational mesh handling
 * - Element: for finite element operations  
 * - FiniteElementSpace: for finite element space definition
 * - GridFunction: for solution representation
 * - BilinearForm: for assembling system matrices
 * - LinearForm: for assembling load vectors
 */
class CFDSolver {
public:
    /**
     * @brief Constructor
     * @param mesh_file Path to mesh file
     * @param order Finite element order
     */
    CFDSolver(const std::string& mesh_file, int order = 1);

    /**
     * @brief Destructor
     */
    ~CFDSolver();

    /**
     * @brief Initialize the CFD problem
     * Sets up finite element spaces, boundary conditions, etc.
     */
    void Initialize();

    /**
     * @brief Solve the CFD problem
     * Assembles system and solves linear/nonlinear system
     */
    void Solve();

    /**
     * @brief Save solution to file
     * @param filename Output filename
     */
    void SaveSolution(const std::string& filename);

    /**
     * @brief Get the mesh
     * @return Pointer to the mesh
     */
    mfem::Mesh* GetMesh() { return mesh_; }

    /**
     * @brief Get the finite element space
     * @return Pointer to the finite element space
     */
    mfem::FiniteElementSpace* GetFESpace() { return fespace_; }

    /**
     * @brief Get the solution
     * @return Pointer to the solution grid function
     */
    mfem::GridFunction* GetSolution() { return solution_; }

private:
    // MFEM objects - demonstrating interface with core MFEM classes
    mfem::Mesh* mesh_;                    ///< Computational mesh
    mfem::FiniteElementCollection* fec_;  ///< Finite element collection
    mfem::FiniteElementSpace* fespace_;   ///< Finite element space
    mfem::GridFunction* solution_;        ///< Solution grid function
    mfem::BilinearForm* a_;              ///< Bilinear form (system matrix)
    mfem::LinearForm* b_;                ///< Linear form (RHS vector)
    
    int order_;                          ///< Finite element order
    std::string mesh_file_;              ///< Mesh file path

    /**
     * @brief Setup finite element spaces
     */
    void SetupFESpaces();

    /**
     * @brief Setup boundary conditions
     */
    void SetupBoundaryConditions();

    /**
     * @brief Assemble system matrix and RHS
     */
    void AssembleSystem();
};

} // namespace mofe

#endif // CFD_SOLVER_HPP