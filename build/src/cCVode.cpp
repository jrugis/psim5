/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Example problem:
 * 
 * The following is a simple example problem, with the coding
 * needed for its solution by CVODE. The problem is from
 * chemical kinetics, and consists of the following three rate
 * equations:         
 *    dy1/dt = -.04*y1 + 1.e4*y2*y3
 *    dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*(y2)^2
 *    dy3/dt = 3.e7*(y2)^2
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions: y1 = 1.0, y2 = y3 = 0. The problem is stiff.
 * While integrating the system, we also use the rootfinding
 * feature to find the points at which y1 = 1e-4 or at which
 * y3 = 0.01. This program solves the problem with the BDF method,
 * Newton iteration with the SUNDENSE dense linear solver, and a
 * user-supplied Jacobian routine.
 * It uses a scalar relative tolerance and a vector absolute
 * tolerance. Output is printed in decades from t = .4 to t = 4.e10.
 * Run statistics (optional outputs) are printed at the end.
 * -----------------------------------------------------------------*/

#include <string>
#include <fstream>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */

#include "global_defs.hpp"
#include "utils.hpp"
#include "cLumen.hpp"
#include "cCVode.hpp"


/*
 * f routine. Compute function f(t,y). 
 */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
  cLumen* pt_cLumen = static_cast<cLumen *>(user_data);
  int nvars = pt_cLumen->get_nvars();

  MatrixX1C ymat;
  ymat.resize(nvars, Eigen::NoChange);
  for (int i = 0; i < nvars; i++) {
    ymat(i) = NV_Ith_S(y, i);
  }

  MatrixX1C ydotmat;
  ydotmat.resize(nvars, Eigen::NoChange);
  
  pt_cLumen->fluid_flow_function(t, ymat, ydotmat);
  
  for (int i = 0; i < nvars; i++) {
    NV_Ith_S(ydot, i) = ydotmat(i);
  }

  return(0);
}


/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

cCVode::cCVode(cLumen* lumen_, std::ofstream& out_) : 
    lumen(lumen_), out(out_), initialised(false), nvars(0) {
  y = NULL;
  A = NULL;
  LS = NULL;
  cvode_mem = NULL;
}

cCVode::~cCVode() {
  if (initialised) {
    /* Free y vector */
    N_VDestroy(y);

    /* Free integrator memory */
    CVodeFree(&cvode_mem);

    /* Free the linear solver memory */
    SUNLinSolFree(LS);

    /* Free the matrix memory */
    SUNMatDestroy(A);
  }
}

void cCVode::init(sunindextype nvars_, MatrixX1C& yini) {
  int retval;
  nvars = nvars_;
  realtype t0 = 0.0;

  /* Create serial vector of length NEQ for I.C. */
  y = N_VNew_Serial(nvars);
  check_retval(static_cast<void *>(y), "N_VNew_Serial", 0);

  /* Initialize y */
  for (sunindextype i = 0; i < nvars; i++) {
    NV_Ith_S(y, i) = yini(i);
  }

  /* Set the scalar relative tolerance */
  realtype reltol = 1e-12;
  /* Set the scalar absolute tolerance */
  realtype abstol = 1e-12;

  /* Call CVodeCreate to create the solver memory and specify the 
   * Backward Differentiation Formula */
  cvode_mem = CVodeCreate(CV_BDF);
  check_retval(static_cast<void *>(cvode_mem), "CVodeCreate", 0);
  
  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  retval = CVodeInit(cvode_mem, f, t0, y);
  check_retval(&retval, "CVodeInit", 1);

  /* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  check_retval(&retval, "CVodeSVtolerances", 1);

  /* Call CVodeSetUserData */
  retval = CVodeSetUserData(cvode_mem, static_cast<void *>(lumen));
  check_retval(&retval, "CVodeSetUserData", 1);

  /* Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(nvars, nvars);
  check_retval(static_cast<void *>(A), "SUNDenseMatrix", 0);

  /* Create dense SUNLinearSolver object for use by CVode */
  LS = SUNLinSol_Dense(y, A);
  check_retval(static_cast<void *>(LS), "SUNLinSol_Dense", 0);

  /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  check_retval(&retval, "CVodeSetLinearSolver", 1);

  initialised = true;
}

void cCVode::run(realtype tend) {
  /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */

//  iout = 0;
//  while(1) {
  int retval = CVode(cvode_mem, tend, y, &t, CV_NORMAL);

  check_retval(&retval, "CVode", 1);
//    if (retval == CV_SUCCESS) {
//      break;
//    }
//  }
  if (retval != CV_SUCCESS) {
    utils::fatal_error("CVode did not succeed", out);
  }

  /* Print some final statistics */
  PrintFinalStats(cvode_mem);
}

void cCVode::get_result(MatrixX1C& yout) {
  for (sunindextype i = 0; i < nvars; i++) {
    yout(i) = NV_Ith_S(y, i);
  }
}

/* 
 * Get and print some final statistics
 */

void cCVode::PrintFinalStats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int retval;

  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1);
  retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_retval(&retval, "CVodeGetNumLinSolvSetups", 1);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_retval(&retval, "CVodeGetNumErrTestFails", 1);
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1);

  retval = CVodeGetNumJacEvals(cvode_mem, &nje);
  check_retval(&retval, "CVodeGetNumJacEvals", 1);
  retval = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
  check_retval(&retval, "CVodeGetNumLinRhsEvals", 1);

  retval = CVodeGetNumGEvals(cvode_mem, &nge);
  check_retval(&retval, "CVodeGetNumGEvals", 1);

  out << "Final Statistics:" << std::endl;
  out << "nst = " << nfe << " nfe = " << nfe << " nsetups = " << nsetups << std::endl;
  out << "nfeLS = " << nfeLS << " nje = " << nje << " nni = " << nni << std::endl;
  out << "ncfn = " << ncfn << " netf = " << netf << " nge = " << nge << std::endl;
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

void cCVode::check_retval(void *returnvalue, std::string funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    utils::fatal_error("SUNDIALS_ERROR: " + funcname + "() failed - returned NULL pointer", out);
  }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      utils::fatal_error("SUNDIALS_ERROR: " + funcname + "() failed with retval = " + std::to_string(*retval), out);
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    utils::fatal_error("MEMORY_ERROR: " + funcname + "() failed - returned NULL pointer", out);
  }
}
