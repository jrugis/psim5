/*
 * cLumen.cpp
 *
 *  Created on: 06/12/2018
 *      Author: jrugis
 */

#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

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
  initx();
}

cLumen::~cLumen() {
  out.close();
}

void cLumen::initx() {
  // load initial conditions
  // TODO: currently just loading the txt file, may need to be changed
  out << "<Lumen> intialising variable matrices" << std::endl;

  // number of variables of the system
  ffvars = num_compartments * LUMENALVARS + cell_count * INTRAVARS;
  x_ion.resize(ffvars, Eigen::NoChange);
  dx_ion.resize(ffvars, Eigen::NoChange);
  out << "<Lumen> number of variables: " << ffvars << std::endl;
  
  std::ifstream icfile("IC_new.txt");
  if (not icfile.is_open()) {
    utils::fatal_error("could not open IC_new.txt for fluid flow initial conditions", out);
  }

  std::string line;
  std::vector<tCalcs> ictmp;
  for (int i = 0; i < ffvars; i++) {
    if (not getline(icfile, line)) {
      utils::fatal_error("not enough lines in IC_new.txt", out);
    }
    x_ion(i) = std::stod(line);
    out << "<Lumen> IC value: " << x_ion(i) << std::endl;
  }

  icfile.close();

  // set up variable arrays
  intra.resize(cell_count, INTRAVARS);
  Nal.resize(cell_count, cell_count);
  Kl.resize(cell_count, cell_count);
  Cll.resize(cell_count, cell_count);

  // set up solution arrays
  Jb.resize(cell_count, Eigen::NoChange);
  JCl.resize(cell_count, cell_count);
  JtNa.resize(cell_count, cell_count);
  JtK.resize(cell_count, cell_count);
  Qa.resize(cell_count, cell_count);
  Qtot.resize(cell_count, cell_count);
  JCL.resize(cell_count, Eigen::NoChange);
}


void cLumen::prep_cell_calcium() {
  out << "<Lumen> exchanging info with cells" << std::endl;

  // receive connectivity information from cells
  neigh.resize(cell_count);
  neigh_clust.resize(cell_count);
  cells_exchange_buffer.resize(cell_count);
  num_compartments = 0;
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
        num_compartments++;
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
  out << "<Lumen> number of lumenal compartments: " << num_compartments << std::endl;

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

  // receive info about areas of basal regions of the cells
  basal_areas.resize(cell_count);
  for (int i = 0; i < cell_count; i++) {
    MPI_Status stat;
    MPI_CHECK(MPI_Recv(&basal_areas[i], 1, MPI_DOUBLE, cell_rank + i, LUMEN_CELL_TAG, MPI_COMM_WORLD, &stat));
  }
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
  // just call fluid flow function for testing now
  fluid_flow_function(t, x_ion);


  // TODO: send volumes back to cells

}

void cLumen::run() {
}

void cLumen::fluid_flow_function(tCalcs t, MatrixX1C &x) {
  // arrange inputs (matlab var.m)
  var(x);
  
  // Basolateral fluxes and Membrane Potentials
  fx_ba();

  // Apical fluxes
  fx_ap();

  // Sum all apical chloride fluxes
  JCL = JCl.rowwise().sum();
  Jb.col(Jwater) = Jb.col(Qb) - Qa.rowwise().sum();

  // Set the intracellular concentration differential equations
  ieq();

//  lum_adj();
}

void cLumen::ieq() {
  // This function calculates the differential equations of the intracellular
  // ionic concentrations and the membrane potentials of any given cell
  // using the int matrix of seven columns
  // (that represent each cell) and 8 rows (representing concentrations and PM
  // potentials along with cell volume). 
  // The order is in the intracellular_variables enum.



}

void cLumen::fx_ba() {
  // This function takes as an input the intracellular variables of the cell model and the Ca
  // in order to calculate the membrane ionic fluxes, and the flow rate into the cell, at the
  // basolateral side of the cells.

  // loop over the cells
  for (int cell_no = 0; cell_no < cell_count; cell_no++) {
    // Basolateral flow rate
    Jb(cell_no, Qb) = p[Lb] * (2.0 * (intra(cell_no, 1) + intra(cell_no, 2) + intra(cell_no, 5)) + p[CO20] - p[Ie]);

    // K+ Nernst Potential
    tCalcs VK = RTF * log(p[Ke] / intra(cell_no, 2));

    // Ca2+ Activated K+ Channels open probability
    // Note: PK was precalculated in cCell_calcium and stored in the last element of `cells_exchange_buffer[cell_no]`
    tCalcs PK = cells_exchange_buffer[cell_no].back();

    // Ca2+ Activated K+ Channels Total Flux
    Jb(cell_no, JK) = p[GK] * PK * (intra(cell_no, 7) - VK) / F;

    // 3Na+/2K+ ATPases
    Jb(cell_no, JNaK) = basal_areas[cell_no] * p[aNaK] * (p[r] * pow(p[Ke], 2) * pow(intra(cell_no, 1), 3) / 
        (pow(p[Ke], 2) + p[alpha1] * pow(intra(cell_no, 1), 3)));

    // Na+ K+ 2Cl- Cotransporters
    Jb(cell_no, JNkcc1) = p[aNkcc1] * basal_areas[cell_no] * (p[a1] - p[a2] * intra(cell_no, 1) * intra(cell_no, 2) *
        pow(intra(cell_no, 3), 2)) / (p[a3] + p[a4] * intra(cell_no, 1) * intra(cell_no, 2) * pow(intra(cell_no, 3), 2));

    // Na+-2HCO3-/Cl- Anion Exchanger 4
    Jb(cell_no, JAe4) = basal_areas[cell_no] * p[G4] * ((p[Cle] / (p[Cle] + p[KCl])) * (intra(cell_no, 1) / 
          (intra(cell_no, 1) + p[KNa])) * pow((intra(cell_no, 4) / (intra(cell_no, 4) + p[KB])), 2));

    // Na+/H+ Antiporter 1
    Jb(cell_no, JNhe1) = basal_areas[cell_no] * p[G1] * ((p[Nae] / (p[Nae] + p[KNa])) * (intra(cell_no, 5) /
          (p[KH] + intra(cell_no, 5))) - (intra(cell_no, 1) / (intra(cell_no, 1) + p[KNa])) * (p[He] / (p[KH] + p[He])));

    // CAIV Intracellular Carbonic Anhydrases
    Jb(cell_no, JBB) = intra(cell_no, 0) * p[GB] * (p[kp] * p[CO20] - p[kn] * intra(cell_no, 4) * intra(cell_no, 5));

    // Intracellular Osmolarity (Using electroneutrality principle)
    Jb(cell_no, Ii) = 2.0 * (intra(cell_no, 1) + intra(cell_no, 2) + intra(cell_no, 5)) + p[CO20];
  }
}

void cLumen::fx_ap() {
  // This function takes as an input the intracellular and luminal variables 
  // of the cell model and the Ca in order to calculate the apical 
  // mebrane ionic fluxes, and the flow rate out and into the lumen, 
  // of any particular cell. 

  // loop over cells
  for (int c_no = 0; c_no < cell_count; c_no++) {
    // loop over cell neighbours
    int num_neigh = neigh[c_no].size();
    for (int j = 0; j < num_neigh; j++) {
      int ngh = neigh[c_no][j];
      tCalcs Sa_p = apical_area_ratios[c_no][j];  // ratio of shared apical area with this neighbour to full apical area

      // Tight Junctional Membrane Potential
      tCalcs Vt = intra(c_no, 6) - intra(c_no, 7);

      // Ca2+ Activated Apical Cl- Channels
      // Note: PrCl was computed on cCell_calcium and stored in `cells_exchange_buffer[c_no]`
      tCalcs PrCl = cells_exchange_buffer[c_no][j];
      tCalcs VCl = Sa_p * RTF * log(Cll(c_no, ngh) / intra(c_no, 3));
      JCl(c_no, ngh) = p[GCl] * PrCl * (intra(c_no, 6) + (VCl / Sa_p)) / F;

      // Tight Junctional Fluxes
      tCalcs VtNa = Sa_p * RTF * log(Nal(c_no, ngh) / p[Nae]);
      JtNa(c_no, ngh) = Sa_p * p[GtNa] * p[St] * (Vt - (VtNa / Sa_p)) / F;

      tCalcs VtK = Sa_p * RTF * log(Kl(c_no, ngh) / p[Ke]);
      JtK(c_no, ngh) = Sa_p * p[GtK] * p[St] * (Vt - (VtK / Sa_p)) / F;

      // Luminal Osmolarity (Using electroneutrality principle)
      tCalcs Il = 2.0 * Cll(c_no, ngh) + p[Ul];

      // Flow Rates
      // Apical Flow Rate
      Qa(c_no, ngh) = Sa_p * p[La] * (Il - Jb(c_no, Ii));
      // Tight Junctional Flow Rate
      tCalcs Qt = Sa_p * p[Lt] * (Il - p[Ie]);
      // Total Fluid Flow Rate (into/out of lumen)
      Qtot(c_no, ngh) = Qa(c_no, ngh) + Qt;
    }
  }
}

void cLumen::var(MatrixX1C &x) {
  // ordering of variables in the x vector is like:
  //    `INTRAVARS` intracellular vars for cell 1
  //    `INTRAVARS` intracellular vars for cell 2
  //    ...
  //    `INTRAVARS` intracellular vars for cell `cell_count`
  //    `num_compartments` vars for lumenal sodium
  //    `num_compartments` vars for lumenal potassium
  //    `num_compartments` vars for lumenal chloride

  intra.setZero();
  Nal.setZero();
  Kl.setZero();
  Cll.setZero();
  int n = 0;
  for (int i = 0; i < cell_count; i++) {
    // intracellular variables
    for (int j = 0; j < INTRAVARS; j++) {
      intra(i, j) = x(n++);
    }

    // lumenal variables
    int num_neigh = neigh_clust[i].size();
    for (int j = 0; j < num_neigh; j++) {
      int neigh_index = neigh_clust[i][j];
      Nal(i, neigh_index) = x(n);
      Kl(i, neigh_index) = x(n + num_compartments);
      Cll(i, neigh_index) = x(n + 2 * num_compartments);
      if (i != j) {
        Nal(neigh_index, i) = Nal(i, neigh_index);
        Kl(neigh_index, i) = Kl(i, neigh_index);
        Cll(neigh_index, i) = Cll(i, neigh_index);
      }
      n++;
    }
  }

  // Lumenal concentration matrix
  // NOTE: the lumenal concentration matrix is essentially a lower
  // triangular matrix whose rows represent cell number and its columns
  // represent neighbour. However, for calculation of membrane
  // potentials and tight junctional fluxes some of these must be
  // repeated. However they are not variables of the system. To get
  // around that, I just add the matrix's transpose to the variable
  // matrix and I obtain what I want.

  out << "<Lumen> DEBUGGING: var Nal:" << std::endl;
  out << Nal << std::endl;
}
