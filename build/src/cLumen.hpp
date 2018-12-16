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
#include "cCell_flow.hpp"

class cLumen {
public:
  cLumen(std::string host_name, int my_rank, int cell_count, int acinus_rank);
  ~cLumen();
  void run();

private:
  std::string id;
  std::ofstream out;
  tCalcs p[FPCOUNT];  // the fluid flow parameters array
  cCell_flow *cells[CELLS_COUNT];  // the cells connected to this lumen
  int my_rank, cell_rank, cell_count, acinus_rank;
  tCalcs na; // sodium
  tCalcs k;  // patasium 
};

#endif /* CLUMEN_ */

