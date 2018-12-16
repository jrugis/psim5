/*
 * cLumen.cpp
 *
 *  Created on: 06/12/2018
 *      Author: jrugis
 */

#include <iostream>
#include <string>

#include "global_defs.hpp"
#include "utils.hpp"
#include "cCell_flow.hpp"
#include "cLumen.hpp"

cLumen::cLumen(std::string host_name, int rank, int c_count, int a_rank) {
  my_rank = rank;
  cell_count = c_count;
  acinus_rank = a_rank;
  id = "l1";

  out.open(id + ".out");
  out << "<Lumen> id: " << id << std::endl;
  out << "<Lumen> host_name: " << host_name << std::endl;

  utils::get_parameters(id, flowParms, 1, p, out);
  for(int i = 0; i < CELLS_COUNT; i++) cells[i] = (new cCell_flow(i, p, out));
}

cLumen::~cLumen() {
  for(int i = 0; i < CELLS_COUNT; i++) delete cells[i];
  out.close();
}

void cLumen::run() {
}


