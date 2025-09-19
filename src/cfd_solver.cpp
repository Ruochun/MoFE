#include "cfd_solver.hpp"
#include <iostream>
#include <memory>
#include <cmath>

namespace mofe {

// Tunables (simple demo)
static constexpr double kRho = 1.0;
static constexpr double kNu = 1e-2;
static constexpr double kDt = 1e-2;
static constexpr int kNSteps = 50;

// Utility: detect tensor vs simplex
static inline bool IsTensorMesh(const mfem::Mesh* m) {
    const int g = m->GetElementBaseGeometry(0);
    return (g == mfem::Geometry::SQUARE || g == mfem::Geometry::CUBE);
}
static inline bool IsSimplicialMesh(const mfem::Mesh* m) {
    const int g = m->GetElementBaseGeometry(0);
    return (g == mfem::Geometry::TRIANGLE || g == mfem::Geometry::TETRAHEDRON);
}

// ---------------- CFDSolver ----------------

CFDSolver::CFDSolver(const std::string& mesh_file, int order)
    : mesh_(nullptr),
      fec_(nullptr),
      fes_(nullptr),
      u_(nullptr),
      u_old_(nullptr),
      a_mass_(nullptr),
      a_diff_(nullptr),
      b_rhs_(nullptr),
      order_(order),
      mesh_file_(mesh_file) {}

CFDSolver::~CFDSolver() {
    delete b_rhs_;
    delete a_diff_;
    delete a_mass_;
    delete u_old_;
    delete u_;
    delete fes_;
    delete fec_;
    delete mesh_;
}

void CFDSolver::Initialize() {
    std::cout << "[CFD] Initializing...\n";

    mesh_ = new mfem::Mesh(mesh_file_.c_str(), 1, 1);
    if (mesh_->Dimension() != 2) {
        throw std::runtime_error("This demo expects a 2D mesh.");
    }
    std::cout << "  dim=" << mesh_->Dimension() << "  NE=" << mesh_->GetNE() << "  NV=" << mesh_->GetNV() << "\n";

    const bool is_tensor = IsTensorMesh(mesh_);
    const bool is_simplex = IsSimplicialMesh(mesh_);
    if (is_tensor) {
        std::cout << "  [mesh] tensor-product (quad) detected.\n";
    }
    if (is_simplex) {
        std::cout << "  [mesh] simplicial (triangle) detected.\n";
    }

    // FE space: H1 vector (2 components)
    const int basis =
        is_tensor ? mfem::BasisType::GaussLobatto : mfem::BasisType::Positive;  // safe on triangles (CPU path)
    fec_ = new mfem::H1_FECollection(order_, /*dim=*/2, basis);
    fes_ = new mfem::FiniteElementSpace(mesh_, fec_, /*vdim=*/2);

    std::cout << "  true dofs: " << fes_->GetTrueVSize() << "\n";

    u_ = new mfem::GridFunction(fes_);
    *u_ = 0.0;
    u_old_ = new mfem::GridFunction(fes_);
    *u_old_ = 0.0;

    // Tiny swirling IC (just for visuals)
    struct VortIC : public mfem::VectorCoefficient {
        VortIC() : mfem::VectorCoefficient(2) {}
        void Eval(mfem::Vector& V, mfem::ElementTransformation& T, const mfem::IntegrationPoint& ip) override {
            mfem::Vector x;
            T.Transform(ip, x);
            double dx = x[0] - 0.5, dy = x[1] - 0.5, r2 = dx * dx + dy * dy;
            double s = std::exp(-50.0 * r2);
            V.SetSize(2);
            V[0] = -dy * s;
            V[1] = dx * s;
        }
    } u0;
    u_->ProjectCoefficient(u0);
    *u_old_ = *u_;

    // Essential BC: velocity zero on all boundaries
    ess_bdr_.DeleteAll();
    ess_bdr_.SetSize(mesh_->bdr_attributes.Max());
    ess_bdr_ = 1;
    fes_->GetEssentialTrueDofs(ess_bdr_, ess_tdof_);

    BuildOperators_(is_tensor, is_simplex);

    std::cout << "[CFD] Init done.\n";
}

void CFDSolver::BuildOperators_(bool /*is_tensor_arg*/, bool /*is_simplex_arg*/) {}

static void ApplyConvectionExplicit(const mfem::FiniteElementSpace* fes,
                                    const mfem::GridFunction& u,
                                    mfem::Vector& out) {
    // Very simple, host-side explicit convection: N(u) ≈ (u · ∇)u
    // For a demo, we’ll approximate with nodal gradient from GridFunction methods.
    // (This is intentionally minimal; delete if you want pure Stokes.)
    out.SetSize(fes->GetTrueVSize());
    out = 0.0;

    // Compute a crude nodal gradient via finite differences on the mesh (host).
    // For brevity, we’ll just make it zero (comment out this function to remove advection).
    // Production codes should implement a proper PA ConvectionNLFIntegrator with CEED or custom kernels.
}

void CFDSolver::Solve() {}

void CFDSolver::SaveSolution(const std::string& filename) {
    std::cout << "Saving solution to " << filename << "\n";
    std::ofstream ofs(filename.c_str());
    ofs.precision(8);
    mesh_->PrintVTK(ofs);
    u_->SaveVTK(ofs, "u", 0);
    ofs.close();
    std::cout << "Saved.\n";
}

// Accessors
mfem::Mesh* CFDSolver::GetMesh() {
    return mesh_;
}
mfem::FiniteElementSpace* CFDSolver::GetFESpace() {
    return fes_;
}
mfem::GridFunction* CFDSolver::GetSolution() {
    return u_;
}

}  // namespace mofe
