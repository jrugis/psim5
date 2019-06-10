#ifndef _CLSODA_H_
#define _CLSODA_H_

#include <string>
#include <fstream>

#include "libsoda/LSODA.h"

#include "global_defs.hpp"


class cLumen;

class cLSODA {
  public:
    cLSODA(cLumen* lumen_, std::ofstream& out_);
    void init(MatrixX1C& y);
    void run(tCalcs t, tCalcs tout, MatrixX1C& y);

  private:
    cLumen *lumen;
    std::ofstream& out;
    int nvars;
    std::vector<tCalcs> yin, yout;
};

#endif
