#ifndef _CLSODA_H_
#define _CLSODA_H_

#include <string>
#include <fstream>

#include "libsoda/LSODA.h"

#include "global_defs.hpp"


class cLumen;

class cLSODA {
  public:
    cLSODA(cLumen* lumen_, std::ofstream& out_, double abstol_, double reltol_);
    void init(MatrixN1d& y);
    void run(double t, double tout, MatrixN1d& y);

  private:
    cLumen *lumen;
    std::ofstream& out;
    int nvars;
    std::vector<double> yin, yout;
    double abstol, reltol;  // tolerances for the solver
};

#endif
