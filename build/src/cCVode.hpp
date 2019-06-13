#ifndef _CCVODE_H_
#define _CCVODE_H_

#include <string>
#include <fstream>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */

#include "global_defs.hpp"


class cLumen;

class cCVode {
  public:
    cCVode(cLumen* lumen_, std::ofstream& out_, realtype abstol_, realtype reltol_);
    ~cCVode();
    void init(MatrixX1C& yini);
    void run(realtype t, realtype tend, MatrixX1C& yout);

  private:
    cLumen *lumen;
    std::ofstream& out;
    bool initialised;
    sunindextype nvars;
    N_Vector y;
    SUNMatrix A;
    SUNLinearSolver LS;
    void *cvode_mem;
    realtype abstol, reltol;

    void check_retval(void *returnvalue, std::string funcname, int opt);
    void PrintFinalStatsBrief(void *cvode_mem);
    void PrintFinalStatsDetailed(void *cvode_mem);
};

#endif
