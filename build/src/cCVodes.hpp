/*
 * cCVodes.hpp
 *
 *  Created on: 14/12/2018
 *      Author: jrugis
 */

#ifndef CCVODES_H_
#define CCVODES__H_

#include <fstream>

class cCVodes {
public:
  cCVodes(std::ofstream& out);
  ~cCVodes();

private:
  realtype reltol, t, tout;
  N_Vector y, abstol;
  SUNMatrix A;
  SUNLinearSolver LS;
  void *cvode_mem;
  int retval, retvalr, iout;
  int rootsfound[2];

  int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
  int g(realtype t, N_Vector y, realtype *gout, void *user_data);
  int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, 
	               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  int check_retval(void *returnvalue, const char *funcname, int opt);
  int check_ans(N_Vector y, realtype t, realtype rtol, N_Vector atol);
};

#endif /* CCVODES__ */

