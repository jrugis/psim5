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
  // receive connectivity information from cells
  neigh.resize(cell_count);
  neigh_clust.resize(cell_count);
  cells_exchange_buffer.resize(cell_count);
  for (int i = 0; i < cell_count; i++) {
    int num_neigh;
    int src = cell_rank + i;
    MPI_Status stat;
    MPI_CHECK(MPI_Recv(&num_neigh, 1, MPI_INT, src, LUMEN_CELL_TAG, MPI_COMM_WORLD, &stat));
    neigh[i].resize(num_neigh);
    MPI_CHECK(MPI_Recv(neigh[i].data(), num_neigh, MPI_INT, src, LUMEN_CELL_TAG, MPI_COMM_WORLD, &stat));
    for (int j = 0; j < num_neigh; j++) {
      int n = neigh[i][j];
      if (n <= i) {
        neigh_clust[i].push_back(n);
      }
    }

    // initialise exchange buffer with correct size
    cells_exchange_buffer[i].resize(num_neigh + 1);
  }
  
  // echo connectivity information
  for (std::vector<std::vector<int> >::size_type i = 0; i < neigh.size(); i++) {
    out << "<Lumen> neigh(" << i + 1 << "):";
    for (std::vector<int>::size_type j = 0; j < neigh[i].size(); j++) {
      out << " " << neigh[i][j] + 1;
    }
    out << std::endl;
  }
  for (std::vector<std::vector<int> >::size_type i = 0; i < neigh_clust.size(); i++) {
    out << "<Lumen> neigh_clust(" << i + 1 << "):";
    for (std::vector<int>::size_type j = 0; j < neigh_clust[i].size(); j++) {
      out << " " << neigh_clust[i][j] + 1;
    }
    out << std::endl;
  }

  // send fluid flow parameters (all for now)
  for (int i = 0; i < cell_count; i++) {
    int dest = cell_rank + i;
    MPI_CHECK(MPI_Send(p, FPCOUNT, MPI_DOUBLE, dest, LUMEN_CELL_TAG, MPI_COMM_WORLD));
  }

  // receive info about common apical region areas / ratios
  apical_area_ratios.resize(cell_count);
  for (int i = 0; i < cell_count; i++) {
    int num_neigh = neigh[i].size();
    apical_area_ratios[i].resize(num_neigh);
    MPI_Status stat;
    MPI_CHECK(MPI_Recv(apical_area_ratios[i].data(), num_neigh, MPI_DOUBLE, cell_rank + i, LUMEN_CELL_TAG, MPI_COMM_WORLD, &stat));
  }

  // TODO: initial conditions for solver variables
  
}

void cLumen::iterate(tCalcs t, tCalcs dt) {
  out << "<Lumen> iteration: t = " << t << " (dt = " << dt << ")" << std::endl;

  // receive Ca inputs
  MPI_Request recv_requests[cell_count];
  for (int i = 0; i < cell_count; i++) {
    MPI_CHECK(MPI_Irecv(cells_exchange_buffer[i].data(), cells_exchange_buffer[i].size(), MPI_DOUBLE,
          cell_rank + i, LUMEN_CELL_TAG, MPI_COMM_WORLD, &recv_requests[i]));
  }

  // wait for data to come in
  for (int i = 0; i < cell_count; i++) {
    MPI_Status status;
    int recv_index;
    MPI_CHECK(MPI_Waitany(cell_count, recv_requests, &recv_index, &status));
    // got input from cell number: recv_index -- do something with it...
  }

  // TODO: solve


  // TODO: send volumes back to cells

}

void cLumen::run() {
}


