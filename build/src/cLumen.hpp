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

#define INTRAVARS 8
#define LUMENALVARS 3

class cLumen {
public:
  cLumen(std::string host_name, int rank, int c_rank, int c_count);
  ~cLumen();
  void run();
  void iterate(tCalcs t, tCalcs dt);

private:
  void initx();
  void prep_cell_calcium();
  void fluid_flow_function(tCalcs t, MatrixX1C &x);
  void var(MatrixX1C &x);

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
  std::vector<std::vector<tCalcs> > cells_exchange_buffer;  // buffer for receiving Ca input values from cells
  MatrixX1C x_ion;  // solution vector
  MatrixXXC intra;  // intra cellular variables for all cells
  MatrixXXC nal, kl, cll;  // lumenal variable arrays
};

#endif /* CLUMEN_ */

