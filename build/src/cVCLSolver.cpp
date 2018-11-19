/*
 * cVCLSolver.cpp
 *
 *  Created on: 23/01/2018
 *      Author: jrugis
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <Eigen/Sparse>

#define VIENNACL_HAVE_EIGEN // enable Eigen wrappers within ViennaCL

#include "cCell_x.hpp"
#include "cCellMesh.hpp"
#include "cVCLSolver.hpp"

cVCLSolver::cVCLSolver(cCell_x* p){
  parent = p;
  parent->out << "<VCLSolver> initialising a VCL solver..." << std::endl;
  vcl_precond = NULL;
  my_tag = new VclTag(1e-10, 300, 20); // NOTE: select tag typoe by changing typedef in header file
  my_solver = new VclSolver(*my_tag);  // NOTE: select solver by changing typedef in header file
  my_solver->set_monitor(NULL, NULL);
}

cVCLSolver::~cVCLSolver(){
  delete vcl_precond;
}

// set A matrix and preconditioner
void cVCLSolver::setA(SparseMatrixTCalcs &sparseA){
  viennacl::copy(sparseA, vcl_sparseA); // eigen -> vcl    

  // NOTE: also select preconditioner type by changing typedef in header file
  viennacl::linalg::ilut_tag precond_tag(40, 1e-5); // preconditioner
  //viennacl::linalg::ichol0_tag precond_tag(); // preconditioner
  //viennacl::linalg::jacobi_tag precond_tag(); // preconditioner

  if(vcl_precond) delete vcl_precond;
  vcl_precond = new VclPrecond(vcl_sparseA, precond_tag);
}

void cVCLSolver::step(MatrixX1C &solvec, MatrixX1C &rhsvec){
  viennacl::vector<tCalcs> vcl_sol;
  viennacl::vector<tCalcs> vcl_rhs(rhsvec.rows());

  // solve
  viennacl::copy(rhsvec, vcl_rhs); // eigen -> vcl

  // ******************
  // ******************
  //viennacl::linalg::gmres_tag my_gmres_tag(1e-9, 300, 20);
  //viennacl::linalg::gmres_tag my_gmres_tag(1e-10, 300, 20);
  //vcl_sol = viennacl::linalg::solve(vcl_sparseA, vcl_rhs, my_gmres_tag, *vcl_precond);
  // ******************
  my_solver->set_initial_guess(vcl_rhs);
  vcl_sol = (*my_solver)(vcl_sparseA, vcl_rhs, *vcl_precond);
  // ******************
  //viennacl::linalg::bicgstab_tag my_bicgstab_tag(1e-8, 400, 200);
  //vcl_sol = viennacl::linalg::solve(vcl_sparseA, vcl_rhs, my_bicgstab_tag);
  // ******************
  // ******************
  // ******************

  viennacl::copy(vcl_sol, solvec); // vcl -> eigen
 
  //parent->out << "<VCLSolver> iterations: " << my_gmres_tag.iters();
  //parent->out << std::fixed << std::setprecision(8);
  //parent->out << " relative error: " << my_gmres_tag.error() << std::endl;
  parent->out << "<VCLSolver> iterations: " << my_tag->iters();
  parent->out << std::fixed << std::setprecision(8);
  parent->out << " relative error: " << my_tag->error() << std::endl;
}

