/*
 * cCell_x.hpp
 *
 *  Created on: 26/04/2018
 *      Author: jrugis
 */

#ifndef CCELL_X_H_
#define CCELL_X_H_

#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>

class cCellMesh;
class cVCLSolver;

#include "global_defs.hpp"

#define DIFVARS 3 // number of diffusing node variables - c, ip, ce
#define NONDIFVARS 2 // number of non-diffusing variables - g, h
#define VARIABLES (DIFVARS+NONDIFVARS) // total number of node variables
#define REF_MASS_SIZE 4   // reference mass dimension

enum model_element_values{VOL_e, RYR_e, PLC_e, MODELECOUNT};  // element volume and spatial factors
enum model_surface_values{AREA_s, MODELSCOUNT};  // surface triangle area
enum model_node_values{BOOL_apical, MODELNCOUNT}; // apical (boolean)

// some convenience typedefs
typedef Eigen::Array<tCalcs, Eigen::Dynamic, 1> ArrayX1C;
typedef Eigen::Array<tCalcs, 1, DIFVARS> Array1VC;
typedef Eigen::Array<tCalcs, REF_MASS_SIZE, REF_MASS_SIZE> ArrayRefMass;
typedef Eigen::Triplet<tCalcs> Triplet;

struct cfc {int cell; int fcount;}; // other cell, connected face count

class cCell_x {
friend class cCellMesh;
friend class cVCLSolver;
public:
  cCell_x(std::string host_name, int my_rank, int acinus_rank);
  ~cCell_x();
  void run();

private:
  std::string id;
  std::string acinus_id;
  std::ofstream out, ca_file, ip3_file;
  int cell_number, acinus_rank;
  cCellMesh* mesh;
  cVCLSolver* solver;

  tCalcs p[PCOUNT]; // the model parameters array
  std::vector<cfc> cells; // vector of connected cells and face counts

  Eigen::Array<tCalcs, Eigen::Dynamic, MODELECOUNT> element_data;
  Eigen::Array<tCalcs, Eigen::Dynamic, MODELSCOUNT> surface_data;
  Eigen::Array<tCalcs, Eigen::Dynamic, MODELNCOUNT> node_data;

  MatrixX1C solvec, prev_solvec; // solution vector
  SparseMatrixTCalcs sparseA, sparseStiff, sparseMass; // A, stiffness and mass matrices

  void init_solvec();
  void make_matrices();
  void exchange();
  void save_results(std::ofstream &data_file, int var);

  MatrixX1C make_load(tCalcs delta_time, bool plc);
  ArrayRefMass make_ref_mass();
  Array1VC get_body_reactions(tCalcs c, tCalcs ip, tCalcs ce, tCalcs g, tCalcs ryr_f, tCalcs plc_f);
  Array1VC get_apical_reactions(tCalcs c, tCalcs ip, tCalcs ce, tCalcs h);
  tCalcs get_g_reaction(tCalcs c, tCalcs g); // RYR dynamics
  tCalcs get_h_reaction(tCalcs c, tCalcs h);// IPR dynamics (apical)
};

#endif /* CCELL_X_H_ */

