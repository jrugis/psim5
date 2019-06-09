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
    cCVode(cLumen* lumen_, std::ofstream& out_);
    ~cCVode();
    void init(sunindextype nvars_, MatrixX1C& yini);
    void run(realtype tend);
    void get_result(MatrixX1C& yout);

  private:
    cLumen *lumen;
    std::ofstream& out;
    bool initialised;
    sunindextype nvars;
    realtype reltol, t;
    N_Vector y, abstol;
    SUNMatrix A;
    SUNLinearSolver LS;
    void *cvode_mem;
    int retval, retvalr, iout;
        
    void check_retval(void *returnvalue, std::string funcname, int opt);
    void PrintFinalStats(void *cvode_mem);

};

#endif
