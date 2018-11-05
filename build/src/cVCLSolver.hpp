/*
 * cVCLSolver.hpp
 *
 *  Created on: 23/01/2018
 *      Author: jrugis
 */

#ifndef CVCLSOLVER_H_
#define CVCLSOLVER_H_

#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/linalg/ilu.hpp>

#include "global_defs.hpp"

class cCell_x;

typedef viennacl::compressed_matrix<tCalcs> VclSparse;
typedef viennacl::linalg::ilut_precond<VclSparse> VclPrecond;

class cVCLSolver {
public:
  cVCLSolver(cCell_x* parent);
  virtual ~cVCLSolver();
  void setA(SparseMatrixTCalcs &sparseA);
  void step(MatrixX1C &solvec, MatrixX1C &rhsvec);

private:
	cCell_x* parent;
    int size;                // number of columns
    VclPrecond *vcl_precond; // preconditioner
    VclSparse vcl_sparseA;   // sparse matrix
};

#endif /* CVCLSOLVER_H_ */
