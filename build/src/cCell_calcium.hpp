/*
 * cCell_calcium.hpp
 *
 *  Created on: 26/04/2018
 *      Author: jrugis
 */

#ifndef CCELL_CACIUM_H_
#define CCELL_CACIUM_H_

#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>

class cCellMesh;

#include "global_defs.hpp"

#define DIFVARS 3 // number of diffusing node variables - c, ip, ce
#define NONDIFVARS 2 // number of non-diffusing variables - g, h
#define VARIABLES (DIFVARS+NONDIFVARS) // total number of node variables
#define REF_MASS_SIZE 4   // reference mass dimension

enum model_element_values{VOL_e, RYR_e, PLC_e, MODELECOUNT};  // element volume and spatial factors
enum model_surface_values{AREA_s, MODELSCOUNT};  // surface triangle area
enum model_node_values{BOOL_apical, MODELNCOUNT}; // apical (boolean)
enum model_surface_region_data{AREA_apical, AREA_basal, MODELSRCOUNT}; // apical and basal region surface areas

// some convenience typedefs
typedef Eigen::Array<tCalcs, 1, DIFVARS> Array1VC;
typedef Eigen::Array<tCalcs, REF_MASS_SIZE, REF_MASS_SIZE> ArrayRefMass;
typedef Eigen::Triplet<tCalcs> Triplet;

struct cfc {int cell; int fcount; int sindex;}; // other cell, connected face count, common triangles start index

class cCell_calcium {
friend class cCellMesh;
public:
  cCell_calcium(std::string host_name, int my_rank, int acinus_rank, int lumen_rank);
  ~cCell_calcium();
  void run();

private:
  std::string id;
  std::string acinus_id;
  std::ofstream out, ca_file, ip3_file, cer_file;
  int cell_number, acinus_rank, lumen_rank;
  cCellMesh* mesh;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<tCalcs>> solver;
  tCalcs p[PCOUNT]; // the model parameters array
  tCalcs fp[FPCOUNT]; // fluid flow parameters
  std::vector<cfc> cells; // vector of connected cells and face counts
  std::vector<cfc> cells_apical; // vector of connected cells and face counts for apical region (fluid flow calculations)
  std::vector<tCalcs> cells_apical_area_ratios; // ratio of area of each connected apical region to total apical area
  tCalcs** exchange_send_buffer;  // buffers for exchanging values between connected cells
  tCalcs** exchange_recv_buffer;
  ArrayX1C exchange_load_ip;
  tCalcs cell_volume_term;  // volume term received from lumen

  Eigen::Array<tCalcs, Eigen::Dynamic, MODELECOUNT> element_data;
  Eigen::Array<tCalcs, Eigen::Dynamic, MODELSCOUNT> surface_data;
  Eigen::Array<tCalcs, Eigen::Dynamic, MODELNCOUNT> node_data;
  tCalcs surface_region_data[MODELSRCOUNT];

  MatrixX1C solvec, nd_solvec, prev_solvec, prev_nd_solvec; // solution vectors (for diffusing and non-diffusing)
  SparseMatrixTCalcs sparseA, sparseStiff, sparseMass; // A, stiffness and mass matrices

  void init_solvec();
  void make_matrices();
  void exchange();
  void compute_exchange_values(int cell);
  void compute_exchange_load(int cell);
  void save_results(std::ofstream &data_file, int var);
  void lumen_prep();
  void lumen_exchange();

  MatrixX1C solve_nd(tCalcs delta_time);
  MatrixX1C make_load(tCalcs delta_time, bool plc);
  ArrayRefMass make_ref_mass();
  Array1VC get_body_reactions(tCalcs c, tCalcs ip, tCalcs ce, tCalcs g, tCalcs ryr_f, tCalcs plc_f);
  Array1VC get_apical_reactions(tCalcs c, tCalcs ip, tCalcs ce, tCalcs h);
  tCalcs get_g_reaction(tCalcs c, tCalcs g); // RYR dynamics
  tCalcs get_h_reaction(tCalcs c, tCalcs h);// IPR dynamics (apical)
};

#endif /* CCELL_CACIUM_H_ */

