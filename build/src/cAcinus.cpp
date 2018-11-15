/*
 * cAcinus.cpp
 *
 *  Created on: 28/04/2018
 *      Author: jrugis
 */

#include <iostream>
#include <string>

#include "global_defs.hpp"
#include "utils.hpp"
#include "cAcinus.hpp"

cAcinus::cAcinus(std::string host_name, int rank, int c_rank, int c_count) {
  my_rank = rank;
  cell_rank = c_rank;
  cell_count = c_count;
  id = "a" + std::to_string(my_rank + 1);

  out.open(id + ".out");
  out << "<Acinus> id: " << id << std::endl;
  out << "<Acinus> host_name: " << host_name << std::endl;

  utils::get_parameters(id, 1, p, out);
  time_file.open(id + "_time.bin", std::ios::binary);
}

cAcinus::~cAcinus() {
  time_file.close();
  out.close();
}

// NOTE: mpi send to all first, then receive from all
tCalcs cAcinus::snd_recv(tCalcs t, tCalcs dt) {
  float msg[ACCOUNT]; 
  MPI_Status stat;

  out << "<Acinus> t: " << t << std::endl;
  msg[dTime] = dt;
  msg[cTime] = t;
  msg[sError] = 0.0;
  for(int n = cell_rank; n < (cell_rank + cell_count); n++){
    MPI_CHECK(MPI_Send(&msg, ACCOUNT, MPI_FLOAT, n, ACINUS_CELL_TAG, MPI_COMM_WORLD));
  }
  for(int n = cell_rank; n < (cell_rank + cell_count); n++){
    MPI_CHECK(MPI_Recv(&msg, ACCOUNT, MPI_FLOAT, n, ACINUS_CELL_TAG, MPI_COMM_WORLD, &stat));
  }
  return(msg[sError]);
}

void cAcinus::run() {
  tCalcs t = 0.0;
  tCalcs solver_dt = p[delT];
  tCalcs prev_dt = solver_dt;
  tCalcs error;

  // simulation time stepping and synchronization
  while(t < p[totalT]) {
    float f = t; // convert to float for reduced file size
    time_file.write(reinterpret_cast<char*>(&f), sizeof(float));
    error = snd_recv(t, solver_dt);
    if(error != 0.0) { // change time step?
      // ...
    }
    t += solver_dt;
  }  

  // instruct cells to finish
  solver_dt = 0.0;
  snd_recv(t, solver_dt);
}


