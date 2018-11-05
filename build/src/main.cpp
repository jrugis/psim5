/*
 * main.cpp
 *
 *  Created on: 26/04/2018
 *      Author: jrugis
 */

#include <iostream>
#include <sys/time.h>
#include <unistd.h>
#include <mpi.h>

#include "global_defs.hpp"
#include "cAcinus.hpp"
#include "cCell_x.hpp"

// one acinus + seven cells (for now, UNTIL WE ADD SOME MORE)
#define ACINUS_RANK 0
#define CELLS_RANK 1
#define CELLS_COUNT 7
#define MPI_NODES (CELLS_COUNT + 1)

#define TEMP_SIZE 40

// the main program function for each mpi node
int main(int argc,char **args){
  int commSize, commRank;
  std::string host_name;
  struct timeval start, end;
  int duration;

  // initialize mpi
  MPI_CHECK(MPI_Init(&argc, &args));
  MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &commSize));
  MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &commRank));
  if(commSize != MPI_NODES) mpi_abort(MPI_NODES_ABORT); // check mpiexec node count

  gettimeofday(&start, NULL); // get the time

  // get the hostname
  char temp[TEMP_SIZE];
  gethostname(temp, TEMP_SIZE);
  host_name = temp;

  //*********************************************************************************
  // This code is running as EITHER an mpi process for the acinus,
  if(commRank == ACINUS_RANK){
    std::cout << "<main> rank " << commRank << " running..." << std::endl;
    cAcinus* acinus = new cAcinus(host_name, commRank, CELLS_RANK, CELLS_COUNT);
    acinus->run();
    delete acinus;
  }
  // OR an mpi process for a cell.
  else{
    cCell_x* cell = new cCell_x(host_name, commRank, ACINUS_RANK);
    cell->run();
    delete cell;
  }
  //*********************************************************************************

  gettimeofday(&end, NULL);
  duration = end.tv_sec - start.tv_sec;
  std::cout << "<main> rank " << commRank << " execution time: " << duration << " sec" << std::endl;

  // shutdown mpi
  MPI_CHECK(MPI_Finalize());

  return 0;
}
