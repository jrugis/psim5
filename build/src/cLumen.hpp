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

class cLumen {
public:
  cLumen(std::string host_name, int rank, int c_rank, int c_count);
  ~cLumen();
  void run();
  void iterate(tCalcs t, tCalcs dt);

private:
  void prep_cell_calcium();
  std::string id;
  std::ofstream out;
  tCalcs p[FPCOUNT];  // the fluid flow parameters array
  int my_rank, cell_rank, cell_count;
  tCalcs na; // sodium
  tCalcs k;  // patasium 
  std::vector<std::vector<int> > neigh;  // connectivity between cells
  std::vector<std::vector<int> > neigh_clust;  // one sided connectivity between cells
};

#endif /* CLUMEN_ */

