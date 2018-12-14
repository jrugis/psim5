/*
 * cCVodes.cpp
 *
 *  Created on: 14/12/2018
 *      Author: jrugis
 *
 *  based on cvodes code example: cvsRoberts_dns.c
 */

#include <cvodes/cvodes.h>             /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */

#include "utils.hpp"
#include "cCVodes.hpp"

/* User-defined vector and matrix accessor macros: Ith, IJth */
/* These macros are defined in order to write code which exactly matches
   the mathematical problem description given above.

   Ith(v,i) references the ith component of the vector v, where i is in
   the range [1..NEQ] and NEQ is defined below. The Ith macro is defined
   using the N_VIth macro in nvector.h. N_VIth numbers the components of
   a vector starting from 0.

   IJth(A,i,j) references the (i,j)th element of the dense matrix A, where
   i and j are in the range [1..NEQ]. The IJth macro is defined using the
   SM_ELEMENT_D macro. SM_ELEMENT_D numbers rows and columns of 
   a dense matrix starting from 0. */
#define Ith(v,i)    NV_Ith_S(v,i-1)         /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) SM_ELEMENT_D(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */

/* Problem Constants */
#define NEQ   3                /* number of equations  */
#define ZERO  RCONST(0.0)

cCVodes::cCVodes(std::ofstream& out) {
  y = abstol = NULL;
  A = NULL;
  LS = NULL;
  cvode_mem = NULL;
  
  /* Create serial vector of length NEQ for I.C. and abstol */
  y = N_VNew_Serial(NEQ);
  if (check_retval((void *)y, "N_VNew_Serial", 0)) utils::fatal_error("<cCVode> ", out);
  abstol = N_VNew_Serial(NEQ); 
  if (check_retval((void *)abstol, "N_VNew_Serial", 0)) utils::fatal_error("<cCVode> ", out);

  /* Initialize y */
  Ith(y,1) = 1.0;
  Ith(y,2) = 0.0;
  Ith(y,3) = 0.0;
  
  /* Set the vector absolute tolerance */
  //Ith(abstol,1) = ATOL1;
  //Ith(abstol,2) = ATOL2;
  //Ith(abstol,3) = ATOL3;

}

cCVodes::~cCVodes() {
}

// f routine. Compute function f(t,y).
int cCVodes::f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype y1, y2, y3, yd1, yd3;

  y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);

  yd1 = Ith(ydot,1) = RCONST(-0.04)*y1 + RCONST(1.0e4)*y2*y3;
  yd3 = Ith(ydot,3) = RCONST(3.0e7)*y2*y2;
        Ith(ydot,2) = -yd1 - yd3;

  return(0);
}

// g routine. Compute functions g_i(t,y) for i = 0,1. 
int cCVodes::g(realtype t, N_Vector y, realtype *gout, void *user_data)
{
  realtype y1, y3;

  y1 = Ith(y,1); y3 = Ith(y,3);
  gout[0] = y1 - RCONST(0.0001);
  gout[1] = y3 - RCONST(0.01);

  return(0);
}

// Jacobian routine. Compute J(t,y) = df/dy. *
int cCVodes::Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, 
               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype y2, y3;

  y2 = Ith(y,2); y3 = Ith(y,3);

  IJth(J,1,1) = RCONST(-0.04);
  IJth(J,1,2) = RCONST(1.0e4)*y3;
  IJth(J,1,3) = RCONST(1.0e4)*y2;

  IJth(J,2,1) = RCONST(0.04); 
  IJth(J,2,2) = RCONST(-1.0e4)*y3-RCONST(6.0e7)*y2;
  IJth(J,2,3) = RCONST(-1.0e4)*y2;

  IJth(J,3,1) = ZERO;
  IJth(J,3,2) = RCONST(6.0e7)*y2;
  IJth(J,3,3) = ZERO;

  return(0);
}

int cCVodes::check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
	      funcname, *retval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}

/* compare the solution at the final time 4e10s to a reference solution computed
   using a relative tolerance of 1e-8 and absoltue tolerance of 1e-14 */
int cCVodes::check_ans(N_Vector y, realtype t, realtype rtol, N_Vector atol)
{
  int      passfail=0;        /* answer pass (0) or fail (1) retval */  
  N_Vector ref;               /* reference solution vector        */
  N_Vector ewt;               /* error weight vector              */
  realtype err;               /* wrms error                       */
  realtype ONE=RCONST(1.0);  

  /* create reference solution and error weight vectors */
  ref = N_VClone(y);
  ewt = N_VClone(y);

  /* set the reference solution data */
  NV_Ith_S(ref,0) = RCONST(5.2083495894337328e-08);
  NV_Ith_S(ref,1) = RCONST(2.0833399429795671e-13);
  NV_Ith_S(ref,2) = RCONST(9.9999994791629776e-01);

  /* compute the error weight vector, loosen atol */
  N_VAbs(ref, ewt);
  N_VLinearSum(rtol, ewt, RCONST(10.0), atol, ewt);
  if (N_VMin(ewt) <= ZERO) {
    fprintf(stderr, "\nSUNDIALS_ERROR: check_ans failed - ewt <= 0\n\n");
    return(-1);
  }
  N_VInv(ewt, ewt);   

  /* compute the solution error */
  N_VLinearSum(ONE, y, -ONE, ref, ref);
  err = N_VWrmsNorm(ref, ewt);

  /* is the solution within the tolerances? */
  passfail = (err < ONE) ? 0 : 1; 

  if (passfail) {
    fprintf(stdout, "\nSUNDIALS_WARNING: check_ans error=%g \n\n", err);
  }

  /* Free vectors */
  N_VDestroy(ref);
  N_VDestroy(ewt);

  return(passfail);
}

