/*
 * cCell_calcium.hpp
 *
 *	Created on: 26/04/2018
 *	Author: jrugis
 */

#ifndef CCELL_CACIUM_H_
#define CCELL_CACIUM_H_

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <fstream>
#include <string>

#include "global_defs.hpp"

class cCellMesh;

#define DIFVARS 3                        // number of diffusing node variables - c, ip, ce
#define NONDIFVARS 2                     // number of non-diffusing variables - g, h
#define VARIABLES (DIFVARS + NONDIFVARS) // total number of node variables
#define REF_MASS_SIZE 4                  // reference mass dimension

enum model_element_values { VOL_e, RYR_e, PLC_e, MODELECOUNT }; // element volume and spatial factors
enum model_surface_values { AREA_s, MODELSCOUNT };              // surface triangle area
enum model_node_values { BOOL_apical, MODELNCOUNT };            // apical (boolean)

enum cell2cell { // cell to cell mpi message - triangle index, value
  tInd,
  tVal,
  C2CCOUNT
};

// some convenience typedefs
typedef Eigen::Array<double, Eigen::Dynamic, 1> ArrayX1C;
typedef Eigen::Array<double, 1, DIFVARS> Array1VC;
typedef Eigen::Array<double, REF_MASS_SIZE, REF_MASS_SIZE> ArrayRefMass;
typedef Eigen::Triplet<double> Triplet;

struct cfc {
  int cell;
  int fcount;
  int sindex;
}; // other cell, connected face count, common triangles start index

class cCell_calcium {
  friend class cCellMesh;
  friend class cLumenTree;

  public:
  cCell_calcium(std::string host_name, int my_rank, int acinus_rank);
  ~cCell_calcium();
  void run();

  private:
  std::string id;
  std::string acinus_id;
  std::ofstream out, ca_file, ip3_file, cer_file;
  int cell_number, acinus_rank;
  cCellMesh* mesh;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  double p[PCOUNT];              // the model parameters array
  std::vector<cfc> cells;        // vector of connected cells and face counts
  double** exchange_send_buffer; // buffers for exchanging values between connected cells
  double** exchange_recv_buffer;
  ArrayX1C exchange_load_ip;

  Eigen::Array<double, Eigen::Dynamic, MODELECOUNT> element_data;
  Eigen::Array<double, Eigen::Dynamic, MODELSCOUNT> surface_data;
  Eigen::Array<double, Eigen::Dynamic, MODELNCOUNT> node_data;

  MatrixN1d solvec, nd_solvec, prev_solvec, prev_nd_solvec; // solution vectors (for diffusing and non-diffusing)
  SparceMatrixd sparseA, sparseStiff, sparseMass;           // A, stiffness and mass matrices

  void init_solvec();
  void make_matrices();
  void exchange();
  void compute_exchange_values(int cell);
  void compute_exchange_load(int cell);
  void save_results(std::ofstream& data_file, int var);

  MatrixN1d solve_nd(double delta_time);
  MatrixN1d make_load(double delta_time, bool plc);
  ArrayRefMass make_ref_mass();
  Array1VC get_body_reactions(double c, double ip, double ce, double g, double ryr_f, double plc_f);
  Array1VC get_apical_reactions(double c, double ip, double ce, double h);
  double get_g_reaction(double c, double g); // RYR dynamics
  double get_h_reaction(double c, double h); // IPR dynamics (apical)
};

#endif /* CCELL_CACIUM_H_ */
