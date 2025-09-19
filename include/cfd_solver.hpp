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
    // Construct with mesh path and polynomial order (H1)
    CFDSolver(const std::string& mesh_file, int order);
    ~CFDSolver();

    // Build mesh, FE spaces, BCs, and operators (GPU if possible)
    void Initialize();

    // Advance the solution in time (simple semi-implicit viscous step)
    void Solve();

    // Write mesh + velocity field to VTK
    void SaveSolution(const std::string& filename);

    // Accessors
    mfem::Mesh* GetMesh();
    mfem::FiniteElementSpace* GetFESpace();
    mfem::GridFunction* GetSolution();

  private:
    // Internal helpers
    void BuildOperators_(bool is_tensor, bool is_simplex);

    // Mesh / spaces
    mfem::Mesh* mesh_ = nullptr;
    mfem::FiniteElementCollection* fec_ = nullptr;  // H1 FE collection
    mfem::FiniteElementSpace* fes_ = nullptr;       // vector H1 space (vdim=2)

    // Fields
    mfem::GridFunction* u_ = nullptr;      // current velocity
    mfem::GridFunction* u_old_ = nullptr;  // previous step

    // Linear forms / bilinear forms
    mfem::BilinearForm* a_mass_ = nullptr;  // rho * I (vector mass)
    mfem::BilinearForm* a_diff_ = nullptr;  // mu * grad(u):grad(v)
    mfem::LinearForm* b_rhs_ = nullptr;     // (empty container; used by FormLinearSystem)

    // Boundary conditions (Dirichlet on all boundaries for this demo)
    mfem::Array<int> ess_bdr_;   // size = max boundary attribute, 1=essential
    mfem::Array<int> ess_tdof_;  // list of essential true dofs (velocity)

    // Config
    int order_ = 1;
    std::string mesh_file_;
};

}  // namespace mofe

#endif  // CFD_SOLVER_HPP