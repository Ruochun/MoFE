#include "cfd_solver.hpp"
#include <iostream>
#include <memory>

namespace mofe {

CFDSolver::CFDSolver(const std::string& mesh_file, int order)
    : mesh_(nullptr), fec_(nullptr), fespace_(nullptr), 
      solution_(nullptr), a_(nullptr), b_(nullptr),
      order_(order), mesh_file_(mesh_file) {
}

CFDSolver::~CFDSolver() {
    delete b_;
    delete a_;
    delete solution_;
    delete fespace_;
    delete fec_;
    delete mesh_;
}

void CFDSolver::Initialize() {
    std::cout << "Initializing CFD Solver..." << std::endl;
    
    // Load mesh - interfacing with MFEM Mesh class
    mesh_ = new mfem::Mesh(mesh_file_.c_str(), 1, 1);
    
    // Ensure mesh is 2D or 3D
    int dim = mesh_->Dimension();
    std::cout << "Mesh dimension: " << dim << std::endl;
    std::cout << "Number of elements: " << mesh_->GetNE() << std::endl;
    std::cout << "Number of vertices: " << mesh_->GetNV() << std::endl;
    
    // Setup finite element spaces
    SetupFESpaces();
    
    // Setup boundary conditions
    SetupBoundaryConditions();
    
    std::cout << "CFD Solver initialized successfully." << std::endl;
}

void CFDSolver::SetupFESpaces() {
    std::cout << "Setting up finite element spaces..." << std::endl;
    
    // Create finite element collection - interfacing with MFEM FiniteElementCollection
    // Using H1 elements for scalar CFD problems (like heat equation, potential flow)
    fec_ = new mfem::H1_FECollection(order_, mesh_->Dimension());
    
    // Create finite element space - interfacing with MFEM FiniteElementSpace
    fespace_ = new mfem::FiniteElementSpace(mesh_, fec_);
    
    std::cout << "Number of finite element unknowns: " << fespace_->GetTrueVSize() << std::endl;
    
    // Initialize solution grid function - interfacing with MFEM GridFunction
    solution_ = new mfem::GridFunction(fespace_);
    *solution_ = 0.0; // Initialize with zero
}

void CFDSolver::SetupBoundaryConditions() {
    std::cout << "Setting up boundary conditions..." << std::endl;
    
    // For this simple example, we'll set essential boundary conditions
    // on all boundary attributes (Dirichlet BC)
    mfem::Array<int> ess_tdof_list;
    if (mesh_->bdr_attributes.Size()) {
        mfem::Array<int> ess_bdr(mesh_->bdr_attributes.Max());
        ess_bdr = 1; // Set essential BC on all boundaries
        fespace_->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
    }
    
    std::cout << "Number of essential DOFs: " << ess_tdof_list.Size() << std::endl;
}

void CFDSolver::AssembleSystem() {
    std::cout << "Assembling system matrix and RHS..." << std::endl;
    
    // Create bilinear form - interfacing with MFEM BilinearForm
    a_ = new mfem::BilinearForm(fespace_);
    
    // Add integrators for a simple diffusion problem (like heat equation)
    // This demonstrates interfacing with MFEM integrators
    a_->AddDomainIntegrator(new mfem::DiffusionIntegrator());
    
    // Assemble the system matrix
    a_->Assemble();
    
    // Create linear form - interfacing with MFEM LinearForm  
    b_ = new mfem::LinearForm(fespace_);
    
    // Add source term (unit source for demonstration)
    mfem::ConstantCoefficient one(1.0);
    b_->AddDomainIntegrator(new mfem::DomainLFIntegrator(one));
    
    // Assemble the RHS vector
    b_->Assemble();
    
    std::cout << "System assembled successfully." << std::endl;
}

void CFDSolver::Solve() {
    std::cout << "Solving CFD problem..." << std::endl;
    
    // Assemble system
    AssembleSystem();
    
    // Apply boundary conditions
    mfem::Array<int> ess_tdof_list;
    if (mesh_->bdr_attributes.Size()) {
        mfem::Array<int> ess_bdr(mesh_->bdr_attributes.Max());
        ess_bdr = 1; // Essential BC on all boundaries
        fespace_->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
    }
    
    // Get system matrix and RHS vector
    mfem::SparseMatrix A;
    mfem::Vector B, X;
    a_->FormLinearSystem(ess_tdof_list, *solution_, *b_, A, X, B);
    
    std::cout << "System size: " << A.Size() << " x " << A.Size() << std::endl;
    
    // Solve the linear system - interfacing with MFEM linear solvers
    mfem::GSSmoother M(A);
    mfem::PCG(A, M, B, X, 0, 500, 1e-12, 0.0);
    
    // Recover the solution
    a_->RecoverFEMSolution(X, *b_, *solution_);
    
    std::cout << "CFD problem solved successfully." << std::endl;
}

void CFDSolver::SaveSolution(const std::string& filename) {
    std::cout << "Saving solution to " << filename << std::endl;
    
    // Save in VTK format - interfacing with MFEM I/O capabilities
    std::ofstream mesh_ofs(filename.c_str());
    mesh_ofs.precision(8);
    mesh_->PrintVTK(mesh_ofs);
    solution_->SaveVTK(mesh_ofs, "solution", 0);
    mesh_ofs.close();
    
    std::cout << "Solution saved successfully." << std::endl;
}

} // namespace mofe