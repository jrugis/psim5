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
#include "cLumen.hpp"

cLumen::cLumen(std::string host_name, int rank, int c_rank, int c_count) {
  my_rank = rank;
  cell_count = c_count;
  cell_rank = c_rank;
  id = "l1";

  out.open(id + ".out");
  out << "<Lumen> id: " << id << std::endl;
  out << "<Lumen> host_name: " << host_name << std::endl;

  utils::get_parameters(id, flowParms, 1, p, out);

  prep_cell_calcium();
}

cLumen::~cLumen() {
  out.close();
}

void cLumen::prep_cell_calcium() {
  // TODO: for each cell we need to get number of neighbours and list of neighbours
  // TODO: we also need to send some parameters to them
  
  for (int i = 0; i < cell_count; i++) {
//    int num
//    MPI_CHECK(MPI_Irecv(msg, mlength, MPI_DOUBLE, dest, CELL_CELL_TAG, MPI_COMM_WORLD, &send_requests[i]));
  }
}

void cLumen::run() {
}


