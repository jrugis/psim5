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

#include "viennacl/linalg/gmres.hpp"
//#include "viennacl/linalg/bicgstab.hpp"

#include <viennacl/linalg/ilu.hpp>
//#include <viennacl/linalg/ichol.hpp>

#include "global_defs.hpp"

class cCell_x;

typedef viennacl::compressed_matrix<tCalcs> VclSparse;

typedef viennacl::linalg::ilut_precond<VclSparse> VclPrecond;
//typedef viennacl::linalg::ichol0_precond<VclSparse> VclPrecond;
//typedef viennacl::linalg::jacobi_precond<VclSparse> VclPrecond;
	
typedef viennacl::linalg::gmres_tag VclTag;
typedef viennacl::linalg::gmres_solver<viennacl::vector<tCalcs>> VclSolver;


class cVCLSolver {
public:
  cVCLSolver(cCell_x* parent);
  virtual ~cVCLSolver();
  void setA(SparseMatrixTCalcs &sparseA);
  void step(MatrixX1C &solvec, MatrixX1C &rhsvec);

private:
	cCell_x* parent;
    int size;                // number of columns
    VclSparse vcl_sparseA;   // sparse matrix
    VclPrecond *vcl_precond; // preconditioner
    VclTag *my_tag;          // solver tag
    VclSolver *my_solver;    // solver
};

#endif /* CVCLSOLVER_H_ */
