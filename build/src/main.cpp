/*
 * main.cpp
 *
 *  Created on: 26/04/2018
 *      Author: jrugis
 */

#include <iostream>
#include <time.h>
#include <unistd.h>
#include <mpi.h>

#include "global_defs.hpp"
#include "cAcinus.hpp"
//#include "cLumen.hpp"
#include "cCell_calcium.hpp"

// one acinus + seven cells (for now, UNTIL WE ADD SOME MORE)
#define ACINUS_RANK 0
#define CELLS_RANK 1
//#define LUMEN_RANK (CELLS_RANK + CELLS_COUNT)
// Lumen running on same rank as Acinus, can be changed later if needed
#define LUMEN_RANK 0
//#define MPI_NODES (CELLS_COUNT + 2)
#define MPI_NODES (CELLS_COUNT + 1)

#define TEMP_SIZE 40

// the main program function for each mpi node
int main(int argc,char **args){
  int commSize, commRank;
  std::string host_name;
  struct timespec start, end;
  int duration;

  // initialize mpi
  MPI_CHECK(MPI_Init(&argc, &args));
  MPI_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &commSize));
  MPI_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &commRank));
  if(commSize != MPI_NODES) mpi_abort(MPI_NODES_ABORT); // check mpiexec node count

  clock_gettime(CLOCK_REALTIME, &start);

  // get the hostname
  char temp[TEMP_SIZE];
  gethostname(temp, TEMP_SIZE);
  host_name = temp;

  //*********************************************************************************
  // This code is running as EITHER an mpi process for the acinus,
  if(commRank == ACINUS_RANK){
    std::cout << "<main> rank " << commRank << " running..." << std::endl;
    //cAcinus* acinus = new cAcinus(host_name, commRank, CELLS_RANK, CELLS_COUNT, LUMEN_RANK);
    cAcinus* acinus = new cAcinus(host_name, commRank, CELLS_RANK, CELLS_COUNT);
    acinus->run();
    delete acinus;
  }
  // OR an mpi process for the lumenal flow,
//  else if(commRank == LUMEN_RANK){ 
//    cLumen* lumen = new cLumen(host_name, commRank, CELLS_COUNT, ACINUS_RANK);
//    lumen->run();
//    delete lumen;
//  }
  // OR an mpi process for cellular calcium.
  else{
    cCell_calcium* cell = new cCell_calcium(host_name, commRank, ACINUS_RANK, LUMEN_RANK);
    cell->run();
    delete cell;
  }
  //*********************************************************************************

  clock_gettime(CLOCK_REALTIME, &end);
  duration = end.tv_sec - start.tv_sec;
  std::cout << "<main> rank " << commRank << " execution time: " << duration << "s" << std::endl;

  // shutdown mpi
  MPI_CHECK(MPI_Finalize());

  return 0;
}
