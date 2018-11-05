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
#include "viennacl/linalg/gmres.hpp"
#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/linalg/jacobi_precond.hpp"

#include "cCell_x.hpp"
#include "cCellMesh.hpp"
#include "cVCLSolver.hpp"

cVCLSolver::cVCLSolver(cCell_x* p){
  parent = p;
  parent->out << "<VCLSolver> initialising a VCL solver..." << std::endl;
  vcl_precond = NULL;
}

cVCLSolver::~cVCLSolver(){
  delete vcl_precond;
}

// set A matrix and preconditioner
void cVCLSolver::setA(SparseMatrixTCalcs &sparseA){
  viennacl::copy(sparseA, vcl_sparseA); // eigen -> vcl    

  viennacl::linalg::ilut_tag precond_conf(40, 1e-5); // preconditioner
  //viennacl::linalg::jacobi_tag precond_conf(); // preconditioner

  if(vcl_precond) delete vcl_precond;
  vcl_precond = new VclPrecond(vcl_sparseA, precond_conf);
}

void cVCLSolver::step(MatrixX1C &solvec, MatrixX1C &rhsvec){
  viennacl::vector<tCalcs> vcl_sol;
  viennacl::vector<tCalcs> vcl_rhs(rhsvec.rows());

  // solve
  viennacl::copy(rhsvec, vcl_rhs); // eigen -> vcl

  // ******************
  //viennacl::linalg::gmres_tag my_gmres_tag(1e-8, 500, 15);
  //vcl_sol = viennacl::linalg::solve(vcl_sparseA, vcl_rhs, my_gmres_tag, *vcl_precond);
  // dt 0.0001, 1.0, no nan's but ca spikes at ~step 20 
  // ******************
  //viennacl::linalg::gmres_tag my_gmres_tag(1e-8, 300, 40);
  //vcl_sol = viennacl::linalg::solve(vcl_sparseA, vcl_rhs, my_gmres_tag);
  // dt 0.001, 1.0, cell1 nan at step 31
  // ******************
  viennacl::linalg::gmres_tag my_gmres_tag(1e-8, 500, 25);
  vcl_sol = viennacl::linalg::solve(vcl_sparseA, vcl_rhs, my_gmres_tag, *vcl_precond);
  // no PLC, dt 0.0001, 1.0, cell1 spikes at 20
  // no PLC, dt 0.001, 1.0, cell1 nan at step 4
  // ******************
  //viennacl::linalg::bicgstab_tag my_gmres_tag(1e-8, 400, 200);
  //vcl_sol = viennacl::linalg::solve(vcl_sparseA, vcl_rhs, my_gmres_tag);
  // dt 0.001, 1.0, cell1 nan at step 160
  // ******************
  // ******************
  // ******************

  viennacl::copy(vcl_sol, solvec); // vcl -> eigen
 
  parent->out << std::fixed << std::setprecision(8);
  parent->out << "<VCLSolver> " << my_gmres_tag.iters() << " " << my_gmres_tag.error() << std::endl;
}

