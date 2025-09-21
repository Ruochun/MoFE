#include "mfem.hpp"
#include <iostream>
#include <cmath>

using namespace mfem;
using namespace std;

using namespace mfem;
using namespace std;

struct ScopedTimer {
  mfem::StopWatch sw;
  const char* name;
  ScopedTimer(const char* n) : name(n) { sw.Clear(); sw.Start(); }
  ~ScopedTimer() { sw.Stop(); std::cout << name << " = " << sw.RealTime() << " s\n"; }
};

// Helper: fix a single pressure dof to remove nullspace (mean-zero alternative shown below)
static void FixPressureDof(FiniteElementSpace &pfes, Array<int> &ess_p_tdof, int dof_id = -1)
{
    ess_p_tdof.DeleteAll();
    if (dof_id < 0) { dof_id = 0; }
    ess_p_tdof.Append(dof_id);
}

int main(int, char**)
{

}
