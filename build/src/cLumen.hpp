/*
 * cLumen.hpp
 *
 *  Created on: 06/12/2018
 *      Author: jrugis
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

class cLumen {
public:
  cLumen(std::string host_name, int rank, int c_rank, int c_count);
  ~cLumen();
  void run();
  void iterate(tCalcs t, tCalcs dt);

private:
  void initx();
  void load_adjacency_matrix();
  void prep_cell_calcium();
  void fluid_flow_function(tCalcs t, MatrixX1C &x);
  void var(MatrixX1C &x);
  void fx_ba();
  void fx_ap();
  void ieq();
  void lum_adj();

  std::string id;
  std::ofstream out;
  tCalcs p[FPCOUNT];  // the fluid flow parameters array
  int my_rank, cell_rank, cell_count;
  int num_compartments;  // number of lumenal compartments (from meshes)
  int ffvars;  // number of fluid flow variables
  tCalcs na; // sodium
  tCalcs k;  // patasium 
  std::vector<std::vector<int> > neigh;  // connectivity between cells
  std::vector<std::vector<int> > neigh_clust;  // one sided connectivity between cells
  std::vector<std::vector<tCalcs> > apical_area_ratios;  // for each cell the ratios of areas of connected apical to total apical
  std::vector<tCalcs> basal_areas;  // for each cell the total area of the basal triangles
  std::vector<std::vector<tCalcs> > cells_exchange_buffer;  // buffer for receiving Ca input values from cells
  MatrixX1C x_ion;  // initial value vector
  MatrixX1C dx_ion;  // solution vector
  MatrixXXC intra;  // intra cellular variables for all cells
  MatrixXXC Nal, Kl, Cll;  // lumenal variable arrays
  Eigen::Matrix<tCalcs, Eigen::Dynamic, BASOFLUXCOUNT> Jb;  // basolateral fluxes
  MatrixXXC JCl, JtNa, JtK, Qa, Qtot;  // apical and tight junctional fluxes
  MatrixX1C JCL;
  Eigen::MatrixXi adj;  // adjacency matrix
};

#endif /* CLUMEN_ */

