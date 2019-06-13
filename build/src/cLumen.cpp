/*
 * cLumen.cpp
 *
 *  Created on: 06/12/2018
 *      Author: jrugis
 */

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <utility>
#include <chrono>
#include <iomanip>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "global_defs.hpp"
#include "utils.hpp"
#include "cCVode.hpp"
#include "cLSODA.hpp"
#include "cLumen.hpp"


cLumen::cLumen(std::string host_name, int rank, int c_rank, int c_count) :
    id("l1"), my_rank(rank), cell_rank(c_rank), cell_count(c_count), solver_initialised(false),
    tstride(1), step(0), solver_flag(-1), cvode_solver(NULL), lsoda_solver(NULL) {
  // file for storing standard output
  out.open(id + ".out");
  out << std::fixed << std::setprecision(6);
  out << "<Lumen> id: " << id << std::endl;
  out << "<Lumen> host_name: " << host_name << std::endl;

  // file for storing results (all variable values)
  vars_file.open(id + "_results.dat");
  vars_file << std::scientific << std::setprecision(15);

  // load fluid flow parameters
  utils::get_parameters(id, flowParms, 1, p, out);
}

cLumen::~cLumen() {
  out.close();
  vars_file.close();
  if (solver_initialised) {
    if (solver_flag == 0) {
      delete cvode_solver;
    }
    else if (solver_flag == 1) {
      delete lsoda_solver;
    }
  }
}

void cLumen::init(int tstride_) {
  tstride = tstride_;
  prep_cell_calcium();
  load_adjacency_matrix();
  initx();
  init_solver();
}

void cLumen::init_solver() {
  // which solver to use
  solver_flag = static_cast<int>(p[odeSolver]);

  if (solver_flag == 0) {
    cvode_solver = new cCVode(this, out, p[odeSolverAbsTol], p[odeSolverRelTol]);
    cvode_solver->init(x_ion);
  }
  else if (solver_flag == 1) {
    lsoda_solver = new cLSODA(this, out, p[odeSolverAbsTol], p[odeSolverRelTol]);
    lsoda_solver->init(x_ion);
  }
  else {
    utils::fatal_error("Unrecognised value for odeSolver fluid flow parameter", out);
  }

  solver_initialised = true;
}

void cLumen::initx() {
  // load initial conditions
  // TODO: currently just loading the txt file, may need to be changed
  out << "<Lumen> intialising variable matrices" << std::endl;

  // number of variables of the system
  ffvars = num_compartments * LUMENALVARS + cell_count * INTRAVARS;
  x_ion.resize(ffvars, Eigen::NoChange);
  out << "<Lumen> number of variables: " << ffvars << std::endl;
  
  // load the initial conditions input file
  std::ifstream icfile("flow_init.dat");
  if (not icfile.is_open()) {
    utils::fatal_error("could not open flow_init.dat", out);
  }

  std::string line;
  int count = 0;
  while (getline(icfile, line)) {
    // remove comments
    if (line.data()[0] == '#') continue;
    int ci = line.find_first_of("#");
	if (ci > 0) line = line.substr(0, ci);

    x_ion(count++) = std::stod(line);
  }

  icfile.close();

  if (count != ffvars) {
    utils::fatal_error("not enough lines in flow_init.dat", out);
  }

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

  JtNad_tmp.resize(num_compartments, Eigen::NoChange);
  JtKd_tmp.resize(num_compartments, Eigen::NoChange);
  JCld_tmp.resize(num_compartments, Eigen::NoChange);
  Qtotd_tmp.resize(num_compartments, Eigen::NoChange);
  Nald_tmp.resize(num_compartments, Eigen::NoChange);
  Kld_tmp.resize(num_compartments, Eigen::NoChange);
  Clld_tmp.resize(num_compartments, Eigen::NoChange);

  JtNad.resize(num_compartments, num_compartments);
  JtKd.resize(num_compartments, num_compartments);
  JCld.resize(num_compartments, num_compartments);
  Qtotd.resize(num_compartments, num_compartments);
  Nald.resize(num_compartments, num_compartments);
  Kld.resize(num_compartments, num_compartments);
  Clld.resize(num_compartments, num_compartments);

  QwNa.resize(num_compartments, num_compartments);
  QwK.resize(num_compartments, num_compartments);
  QwCl.resize(num_compartments, num_compartments);
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

void cLumen::load_adjacency_matrix() {
  out << "<Lumen> loading adjacency matrix" << std::endl;

  adj.resize(num_compartments, num_compartments);

  std::ifstream adj_file("flow_adj.dat");
  if (not adj_file.is_open()) {
    utils::fatal_error("could not open flow_adj.dat", out);
  }

  std::string line;
  for (int i = 0; i < num_compartments; i++) {
    if (getline(adj_file, line)) {
      line = boost::trim_right_copy(line);  // remove trailing whitespace
      std::vector<std::string> tokens; 
      boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
      if (static_cast<int>(tokens.size()) == num_compartments) {
        for (int j = 0; j < num_compartments; j++) {
          adj(i, j) = std::stoi(tokens[j]);
        }
      }
      else {
        utils::fatal_error("wrong number of elements in flow_adj.dat on line " + std::to_string(i + 1), out);
      }
    }
    else {
      utils::fatal_error("not enough lines in flow_adj.dat", out);
    }
  }

  adj_file.close();

  out << "<Lumen> Adjacency matrix:" << std::endl;
  out << adj << std::endl;
}

void cLumen::receive_ca_inputs() {
  // receive Ca inputs (non-blocking)
  MPI_Request recv_requests[cell_count];
  for (int i = 0; i < cell_count; i++) {
    MPI_CHECK(MPI_Irecv(cells_exchange_buffer[i].data(), cells_exchange_buffer[i].size(), MPI_DOUBLE,
          cell_rank + i, LUMEN_CELL_TAG, MPI_COMM_WORLD, &recv_requests[i]));
  }

  // wait for data to come in
  MPI_Status recv_statuses[cell_count];
  MPI_CHECK(MPI_Waitall(cell_count, recv_requests, recv_statuses));
}

void cLumen::iterate(tCalcs t, tCalcs dt) {
  step++;
  out << "<Lumen> step: " << step << " current_time: " << t << "s";
  out << " delta_time: " << dt << "s" << std::endl;

  // receive Ca inputs (non-blocking)
  receive_ca_inputs();

  // solve fluid flow
  solve_fluid_flow(t, dt);

  // send volume terms back to cells
  distribute_volume_terms();

  // save results (TODO: should do while non-blocking distribute)
  if(step % tstride == 0) {
    save_variables();
  }
}

void cLumen::save_variables() {
  for (int i = 0; i < ffvars; i++) {
    vars_file << x_ion(i) << " ";
  }
  vars_file << std::endl;
}

void cLumen::solve_fluid_flow(tCalcs t, tCalcs dt) {
  auto start = std::chrono::system_clock::now();

  // call the solver
  if (solver_flag == 0) {
    lsoda_solver->run(t, t + dt, x_ion);
  }
  else if (solver_flag == 1) {
    cvode_solver->run(t, t + dt, x_ion);
  }
  else {
    utils::fatal_error("Unrecognised value for odeSolver fluid flow parameter", out);
  }

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  out << "<Lumen> solver duration: " << elapsed.count() << "s"<< std::endl;
}

void cLumen::distribute_volume_terms() {
  // compute derivate at solution point too (required for volume term)
  MatrixX1C x_ion_dot;
  x_ion_dot.resize(ffvars, Eigen::NoChange);
  fluid_flow_function(0, x_ion, x_ion_dot);

  // send volume and derivative back to cells
  for (int i = 0; i < cell_count; i++) {
    int dest = cell_rank + i;
    int volume_index = i * INTRAVARS + Vol;
    tCalcs cell_volume_terms[2];
    cell_volume_terms[0] = x_ion(volume_index);
    cell_volume_terms[1] = x_ion_dot(volume_index);
    MPI_CHECK(MPI_Send(cell_volume_terms, 2, MPI_DOUBLE, dest, LUMEN_CELL_TAG, MPI_COMM_WORLD));
  }
}

void cLumen::fluid_flow_function(tCalcs t, MatrixX1C &x, MatrixX1C &xdot) {
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
  ieq(xdot);

  // Set the lumenal concentration equations
  // Use adjacency matrix for connectivity of the lumen.
  lum_adj();

  // compute xdot
  int index = cell_count * INTRAVARS;
  for (int i = 0; i < num_compartments; i++) {
    xdot(index++) = JtNad(i, i) + QwNa.col(i).sum() - Qtotd.col(i).sum() * Nald(i, i);
  }
  for (int i = 0; i < num_compartments; i++) {
    xdot(index++) = JtKd(i, i) + QwK.col(i).sum() - Qtotd.col(i).sum() * Kld(i, i);
  }
  for (int i = 0; i < num_compartments; i++) {
    xdot(index++) = - JCld(i, i) + QwCl.col(i).sum() - Qtotd.col(i).sum() * Clld(i, i);
  }
}

void cLumen::matrix_add_upper_to_lower(MatrixXXC &mat) {
  // mat = (mat - triu(mat)) + triu(mat)';
  int m = mat.rows();
  int n = mat.cols();
  if (m != n) {
    utils::fatal_error("cLumen::matrix_add_upper_to_lower only works for square matrices", out);
  }
  else {
    for (int i = 0; i < m; i++) {
      for (int j = i + 1; j < n; j++) {
        mat(j, i) += mat(i, j);
        mat(i, j) = 0;
      }
    }
  }
}

void cLumen::matrix_move_upper_to_lower(MatrixXXC &mat) {
  // mat = triu(mat)';
  int m = mat.rows();
  int n = mat.cols();
  if (m != n) {
    utils::fatal_error("cLumen::matrix_move_upper_to_lower only works for square matrices", out);
  }
  else {
    for (int i = 0; i < m; i++) {
      for (int j = i + 1; j < n; j++) {
        mat(j, i) = mat(i, j);
      }
    }
  }     
}

void cLumen::lum_adj() {
  // This function uses the adjacency matrix to calculate the lumenal structure
  // equations. 
  // This is a dirty function and should be re-written.

  // add upper to lower triangles, zeroing upper
  matrix_add_upper_to_lower(Qtot);
  matrix_add_upper_to_lower(JtNa);
  matrix_add_upper_to_lower(JtK);
  matrix_add_upper_to_lower(JCl);

  // move upper triangle to lower triangle
  matrix_move_upper_to_lower(Nal);
  matrix_move_upper_to_lower(Kl);
  matrix_move_upper_to_lower(Cll);

  // find locations of non-zero elements in Qtot transposed
  // TODO: why transposed if we switch ordering later??
  std::vector<std::pair<int, int> > nonzeros;
  for (int i = 0; i < Qtot.rows(); i++) {
    for (int j = 0; j < Qtot.cols(); j++) {
      if (Qtot(i, j) != 0) {
        nonzeros.push_back(std::make_pair(i, j));
      }
    }
  }

  if (static_cast<int>(nonzeros.size()) != num_compartments) {
    utils::fatal_error("mismatch between number of nonzeros and number of lumenal compartments", out);
  }

  for (int i = 0; i < num_compartments; i++) {
    int row = nonzeros[i].first;
    int col = nonzeros[i].second;
    Nald_tmp(i) = Nal(row, col);
    Kld_tmp(i) = Kl(row, col);
    Clld_tmp(i) = Cll(row, col);
    Qtotd_tmp(i) = Qtot(row, col);
    JtNad_tmp(i) = JtNa(row, col);
    JtKd_tmp(i) = JtK(row, col);
    JCld_tmp(i) = JCl(row, col);
  }

  // multiply by adjacency matrix
  for (int i = 0; i < num_compartments; i++) {
    Nald.col(i) = adj.col(i).cast<tCalcs>() * Nald_tmp;
    Kld.col(i) = adj.col(i).cast<tCalcs>() * Kld_tmp;
    Clld.col(i) = adj.col(i).cast<tCalcs>() * Clld_tmp;
    Qtotd.col(i) = adj.col(i).cast<tCalcs>() * Qtotd_tmp;
    JtNad.col(i) = adj.col(i).cast<tCalcs>() * JtNad_tmp;
    JtKd.col(i) = adj.col(i).cast<tCalcs>() * JtKd_tmp;
    JCld.col(i) = adj.col(i).cast<tCalcs>() * JCld_tmp;
  }

  // precalculate the water/ion influx  
  QwNa = Qtotd * Nald;
  QwK = Qtotd * Kld;
  QwCl = Qtotd * Clld;
  for (int i = 0; i < num_compartments; i++) {
    QwNa(i, i) = 0.0;
    QwK(i, i) = 0.0;
    QwCl(i, i) = 0.0;
  }
}

void cLumen::ieq(MatrixX1C &xdot) {
  // This function calculates the differential equations of the intracellular
  // ionic concentrations and the membrane potentials of any given cell
  // using the int matrix of seven columns
  // (that represent each cell) and 8 rows (representing concentrations and PM
  // potentials along with cell volume). 
  // The order is in the intracellular_variables enum.

  int index = 0;
  for (int c_no = 0; c_no < cell_count; c_no++) {
    // enum basolateral_fluxes{Qb, JNaK, JNkcc1, JAe4, JNhe1, JBB, JK, Ii, Jwater, BASOFLUXCOUNT};
    // enum intracellular_variables{Vol, Naplus, Kplus, Clminus, HCO3minus, Hplus, Va, Vb, INTRAVARS};

    // dw/dt
    xdot(index++) = Jb(c_no, Jwater);

    // d[Na+]i/dt
    xdot(index++) = (Jb(c_no, JNkcc1) - 3.0 * Jb(c_no, JNaK) + Jb(c_no, JNhe1) - Jb(c_no, JAe4) - 
        Jb(c_no, Jwater) * intra(c_no, Naplus)) / intra(c_no, Vol);

    // d[K+]i/dt
    xdot(index++) = (Jb(c_no, JNkcc1) + 2.0 * Jb(c_no, JNaK) - Jb(c_no, JK) - Jb(c_no, Jwater) * intra(c_no, Kplus)) / intra(c_no, Vol);

    // d[Cl-]i/dt
    xdot(index++) = (2.0 * Jb(c_no, JNkcc1) + Jb(c_no, JAe4) + JCL(c_no) - Jb(c_no, Jwater) * intra(c_no, Clminus)) / intra(c_no, Vol);

    // d[HCO3-]i/dt
    xdot(index++) = (Jb(c_no, JBB) - 2.0 * Jb(c_no, JAe4) - Jb(c_no, Jwater) * intra(c_no, HCO3minus)) / intra(c_no, Vol);

    // d[H+]i/dt
    xdot(index++) = (Jb(c_no, JBB) - Jb(c_no, JNhe1) - Jb(c_no, Jwater) * intra(c_no, Hplus)) / intra(c_no, Vol);

    // dVa/dt
    xdot(index++) = - JCL(c_no) - (JtNa.row(c_no).sum() + JtK.row(c_no).sum());

    // dVb/dt
    xdot(index++) = - Jb(c_no, JNaK) - Jb(c_no, JK) + JtNa.row(c_no).sum() + JtK.row(c_no).sum();
  }
}

void cLumen::fx_ba() {
  // This function takes as an input the intracellular variables of the cell model and the Ca
  // in order to calculate the membrane ionic fluxes, and the flow rate into the cell, at the
  // basolateral side of the cells.

  Jb.setZero();

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
    Jb(cell_no, JK) = p[GK] * PK * (intra(cell_no, 7) - VK) / CONST_F;

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

  JCl.setZero();
  JtNa.setZero();
  JtK.setZero();
  Qa.setZero();
  Qtot.setZero();

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
      JCl(c_no, ngh) = p[GCl] * PrCl * (intra(c_no, 6) + (VCl / Sa_p)) / CONST_F;

      // Tight Junctional Fluxes
      tCalcs VtNa = Sa_p * RTF * log(Nal(c_no, ngh) / p[Nae]);
      JtNa(c_no, ngh) = Sa_p * p[GtNa] * p[St] * (Vt - (VtNa / Sa_p)) / CONST_F;

      tCalcs VtK = Sa_p * RTF * log(Kl(c_no, ngh) / p[Ke]);
      JtK(c_no, ngh) = Sa_p * p[GtK] * p[St] * (Vt - (VtK / Sa_p)) / CONST_F;

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

  // intracellular variables
  intra.setZero();
  int n = 0;
  for (int i = 0; i < cell_count; i++) {
    for (int j = 0; j < INTRAVARS; j++) {
      intra(i, j) = x(n++);
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

  // lumenal variables
  Nal.setZero();
  Kl.setZero();
  Cll.setZero();
  for (int i = 0; i < cell_count; i++) {
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
}
