#include <string>
#include <fstream>

#include "libsoda/LSODA.h"

#include "global_defs.hpp"
#include "utils.hpp"
#include "cLumen.hpp"
#include "cLSODA.hpp"


/*
 * f routine. Compute function f(t,y). 
 */

static void f(tCalcs t, tCalcs* y, tCalcs* ydot, void *data) {
  cLumen* pt_cLumen = static_cast<cLumen *>(data);
  int ffvars = pt_cLumen->get_nvars();

  MatrixX1C ymat;
  ymat.resize(ffvars, Eigen::NoChange);
  for (int i = 0; i < ffvars; i++) {
    ymat(i) = y[i];
  }

  MatrixX1C ydotmat;
  ydotmat.resize(ffvars, Eigen::NoChange);
  
  pt_cLumen->fluid_flow_function(t, ymat, ydotmat);
  
  for (int i = 0; i < ffvars; i++) {
    ydot[i] = ydotmat(i);
  }
}

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

cLSODA::cLSODA(cLumen* lumen_, std::ofstream& out_, tCalcs abstol_, tCalcs reltol_) :
    lumen(lumen_), out(out_), nvars(0), abstol(abstol_), reltol(reltol_) {
  out << "<LSODA>: creating LSODA solver" << std::endl;
  out << " tolerances are " << abstol << " (absolute) and ";
  out << reltol << " (relative)" << std::endl;
}

void cLSODA::init(MatrixX1C& y) {
  nvars = y.rows();
  yin.resize(nvars);
  yout.resize(nvars);
}

void cLSODA::run(tCalcs t, tCalcs tout, MatrixX1C& y) {
  // transfer input values into required format
  for (int i = 0; i < nvars; i++) {
    yin[i] = y(i);
  }

  // call the solver
  int istate = 1;
  LSODA lsoda;
  lsoda.lsoda_update(f, nvars, yin, yout, &t, tout, &istate, static_cast<void *>(lumen), reltol, abstol);
  if (istate < 1) {
    utils::fatal_error("lsoda failed to compute the solution", out);
  }

  //  extract the results
  for (int i = 0; i < nvars; i++) {
    y(i) = yout[i + 1];
  }
}
