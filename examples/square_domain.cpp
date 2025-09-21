#include "mfem.hpp"
#include <iostream>
#include <cmath>

using namespace mfem;
using namespace std;

struct ScopedTimer {
  mfem::StopWatch sw;
  const char* name;
  ScopedTimer(const char* n) : name(n) { sw.Clear(); sw.Start(); }
  ~ScopedTimer() { sw.Stop(); std::cout << name << " = " << sw.RealTime() << " s\n"; }
};

int main(int, char**) {
    bool use_cuda = true; // toggle for A/B
    mfem::Device device(use_cuda ? "cuda" : "cpu");
    device.Print();

    // Params
    int order = 3;       // spline/NURBS order
    int nx = 800, ny = 800;  // knot spans

    // NURBS (IGA) mesh on [0,1]^2
    // last 'true' => create as NURBS
    Mesh mesh = Mesh::MakeCartesian2D(nx, ny, Element::QUADRILATERAL,
                                      /*gen_edges=*/true, 1.0, 1.0, /*want_nurbs=*/true);

    H1_FECollection fec(order, mesh.Dimension());  // IGA with H1 on NURBS
    FiniteElementSpace fes(&mesh, &fec);

    cout << "NE=" << mesh.GetNE() << "  true dofs=" << fes.GetTrueVSize() << endl;

    // Dirichlet BCs everywhere
    Array<int> ess_tdof;
    if (mesh.bdr_attributes.Size()) {
        Array<int> ess_bdr(mesh.bdr_attributes.Max());
        ess_bdr = 1;
        fes.GetEssentialTrueDofs(ess_bdr, ess_tdof);
    }

    // RHS: f = 2*pi^2 * sin(pi x) sin(pi y)
    class Source : public FunctionCoefficient {
      public:
        Source() : FunctionCoefficient(f) {}
        static double f(const Vector& X) { return 2.0 * M_PI * M_PI * sin(M_PI * X[0]) * sin(M_PI * X[1]); }
    } f;

    GridFunction x(&fes);
    x = 0.0;
    LinearForm b(&fes);
    b.AddDomainIntegrator(new DomainLFIntegrator(f));
    b.Assemble();

    // a(u,v) = (grad u, grad v) — Partial Assembly for GPU
    BilinearForm a(&fes);
    a.SetAssemblyLevel(AssemblyLevel::PARTIAL);
    a.AddDomainIntegrator(new DiffusionIntegrator());
    {
        ScopedTimer t("Assemble");
        a.Assemble();
    }

    // Form linear system (matrix-free Operator in A)
    OperatorPtr A;
    Vector X, B;
    a.FormLinearSystem(ess_tdof, x, b, A, X, B);

    // ---- FIX: build Jacobi smoother from the BilinearForm (not from *A) ----
    CGSolver cg;
    cg.SetOperator(*A);
    cg.SetRelTol(1e-8);
    cg.SetMaxIter(500);
    cg.SetPrintLevel(1);

    OperatorJacobiSmoother M(a, ess_tdof /*, damping=1.0*/);
    cg.SetPreconditioner(M);

    {
        ScopedTimer t("Solve");
        cg.Mult(B, X);
    }
    a.RecoverFEMSolution(X, b, x);

    // Error vs exact solution u = sin(pi x) sin(pi y)
    FunctionCoefficient u_exact([](const Vector& X) { return sin(M_PI * X[0]) * sin(M_PI * X[1]); });
    cout << "L2 error = " << x.ComputeL2Error(u_exact) << endl;

    // --- Make sure the mesh has nodal coords and both nodes+solution are HOST-resident

    // make sure data are readable on host
    if (mesh.GetNodes()) {
        mesh.GetNodes()->HostRead();
    }
    x.HostRead();

    mfem::ParaViewDataCollection pvd("iga_poisson", &mesh);
    pvd.SetDataFormat(mfem::VTKFormat::BINARY);

    // Start with low-order output (most robust)
    pvd.SetHighOrderOutput(false);  // set true later if you want HO cells
    pvd.RegisterField("solution", &x);
    pvd.Save();  // writes iga_poisson.pvtu + pieces

    return 0;
}
