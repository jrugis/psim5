#include <fstream>
#include <string>

#include "libsoda/LSODA.h"

#include "cLSODA.hpp"
#include "cLumen.hpp"
#include "global_defs.hpp"
#include "utils.hpp"

/*
 * f routine. Compute function f(t,y).
 */

static void f(double t, double* y, double* ydot, void* data)
{
  cLumen* pt_cLumen = static_cast<cLumen*>(data);
  int ffvars = pt_cLumen->get_nvars();

  MatrixN1d ymat;
  ymat.resize(ffvars, Eigen::NoChange);
  for (int i = 0; i < ffvars; i++) { ymat(i) = y[i]; }

  MatrixN1d ydotmat;
  ydotmat.resize(ffvars, Eigen::NoChange);

  pt_cLumen->fluid_flow_function(t, ymat, ydotmat);

  for (int i = 0; i < ffvars; i++) { ydot[i] = ydotmat(i); }
}

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

cLSODA::cLSODA(cLumen* lumen_, std::ofstream& out_, double abstol_, double reltol_)
  : lumen(lumen_), out(out_), nvars(0), abstol(abstol_), reltol(reltol_)
{
  out << std::scientific;
  out << "<LSODA>: creating LSODA solver" << std::endl;
  out << " tolerances are " << abstol << " (absolute) and ";
  out << reltol << " (relative)" << std::endl;
  out << std::fixed;
}

void cLSODA::init(MatrixN1d& y)
{
  nvars = y.rows();
  yin.resize(nvars);
  yout.resize(nvars);
}

void cLSODA::run(double t, double tout, MatrixN1d& y)
{
  // transfer input values into required format
  for (int i = 0; i < nvars; i++) { yin[i] = y(i); }

  // call the solver
  int istate = 1;
  LSODA lsoda;
  lsoda.lsoda_update(f, nvars, yin, yout, &t, tout, &istate, static_cast<void*>(lumen), reltol, abstol);
  if (istate < 1) { utils::fatal_error("lsoda failed to compute the solution", out); }

  //  extract the results
  for (int i = 0; i < nvars; i++) { y(i) = yout[i + 1]; }
}
