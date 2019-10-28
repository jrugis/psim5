/*
 * cLumen.hpp
 *
 *	Created on: 06/12/2018
 *			Author: jrugis
 */

#ifndef CLUMEN_H_
#define CLUMEN_H_

#include <fstream>
#include <string>
#include <vector>

#include "global_defs.hpp"

#define LUMENALVARS 3

enum intracellular_variables{Vol, Naplus, Kplus, Clminus, HCO3minus, Hplus, Va, Vb, INTRAVARS};
enum basolateral_fluxes{Qb, JNaK, JNkcc1, JAe4, JNhe1, JBB, JK, Ii, Jwater, BASOFLUXCOUNT};

// some convenience typedefs
typedef Eigen::Array<double, Eigen::Dynamic, 1> ArrayX1C;
typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> ArrayXXC;
typedef Eigen::Array<double, CELLS_COUNT, 1> Array1Cells;
typedef Eigen::Array<double, CELLS_COUNT, CELLS_COUNT> Array2Cells;

class cCVode;
class cLSODA;

class cLumen {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW  // http://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
  cLumen(const std::string host_name, int rank, int c_rank, int c_count);
  ~cLumen();
  void init(int tstride_);
  void iterate(double t, double dt);
  void fluid_flow_function(double t, MatrixN1d &x, MatrixN1d &xdot);
  int get_nvars() { return ffvars; }

private:
  void initx();
  void init_solver();
  void receive_ca_inputs();
  void solve_fluid_flow(double t, double dt);
  void load_adjacency_matrix();
  void prep_cell_calcium();
  void var(MatrixN1d &x);
  void fx_ba();
  void fx_ap();
  void ieq(MatrixN1d &xdot);
  void lum_adj();
  void matrix_add_upper_to_lower(Array2Cells &mat);
  void matrix_move_upper_to_lower(Array2Cells &mat);
  void save_results();
  double compute_flow_rate();

  std::string id;
  std::ofstream out, vars_file;
  int ffvars;  // number of fluid flow variables
  double p[FPCOUNT];  // the fluid flow parameters array
  int my_rank, cell_rank, cell_count;
  int num_compartments;  // number of lumenal compartments (from meshes)
  std::vector<std::vector<int> > neigh;  // connectivity between cells
  std::vector<std::vector<int> > neigh_clust;  // one sided connectivity between cells
  std::vector<std::vector<double> > apical_area_ratios;  // for each cell the ratios of areas of connected apical to total apical
  std::vector<double> basal_areas;  // for each cell the total area of the basal triangles
  std::vector<std::vector<double> > cells_exchange_buffer;  // buffer for receiving Ca input values from cells
  MatrixN1d x_ion;  // solution vector
  Eigen::Array<double, CELLS_COUNT, INTRAVARS> intra;  // intra cellular variables for all cells
  Array2Cells Nal, Kl, Cll;  // lumenal variable arrays
  Eigen::Array<double, CELLS_COUNT, BASOFLUXCOUNT> Jb;  // basolateral fluxes
  Array2Cells JCl, JtNa, JtK, Qa, Qtot;  // apical and tight junctional fluxes
  Array1Cells JCL;
  Array2Cells Sa_p_full;  // array of apical area ratios
  ArrayX1C JtNad_tmp, JtKd_tmp, JCld_tmp, Qtotd_tmp, Nald_tmp, Kld_tmp, Clld_tmp;  // luminal structure equations
  ArrayXXC JtNad, JtKd, JCld, Qtotd, Nald, Kld, Clld;  // luminal structure equations
  ArrayXXC QwNa, QwK, QwCl;  // water/ion influx
  Eigen::ArrayXXi adj;  // adjacency matrix
  bool solver_initialised;
  int tstride;  // how often to store results to file
  int step;  // step number of the simulation
  int solver_flag;  // choose the solver
  cCVode* cvode_solver;
  cLSODA* lsoda_solver;
};

#endif /* CLUMEN_ */
