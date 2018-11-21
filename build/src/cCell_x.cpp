/*
 * cCell_x.cpp
 *
 *  Created on: 26/04/2018
 *      Author: jrugis
 */

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <sys/stat.h>
#include <time.h>
#include <vector>

#include "global_defs.hpp"
#include "utils.hpp"
#include "cCellMesh.hpp"
////#include "cVCLSolver.hpp"
#include "cCell_x.hpp"

cCell_x::cCell_x(std::string host_name, int my_rank, int a_rank) {
  cell_number = my_rank;
  acinus_rank = a_rank;
  acinus_id = "a" + std::to_string(acinus_rank + 1);
  id = acinus_id + "c" + std::to_string(my_rank);
  out.open(id + ".out");
  out << "<Cell_x> id: " << id << std::endl;
  out << "<Cell_x> host_name: " << host_name << std::endl;

  mesh = new cCellMesh(id, this);
  mesh->print_info();
  out << "<Cell_x> common faces with cells:";
  int other_cell = -1;
  int face_count = 0;
  for(int r = 0; r < mesh->common_triangles.rows(); r++) {
    if(mesh->common_triangles(r, oCell) != other_cell) {
      if(other_cell != -1) {
        cells.push_back({other_cell, face_count});
        face_count = 0;
      }
      other_cell = mesh->common_triangles(r, oCell);
      out << " " << other_cell + 1;  // cells are zero indexed
    }
    face_count++;
  }
  cells.push_back({other_cell, face_count}); // one more time
  out << std::endl;

  utils::get_parameters(acinus_id, cell_number, p, out);
  make_matrices();  // create the constant matrices
  init_solvec(); // initialise solution buffer
  ////solver = new cVCLSolver(this);
  solver = new Eigen::SimplicialLDLT<Eigen::SparseMatrix<tCalcs>>;
  ca_file.open(id + "_ca.bin", std::ios::binary);
  ip3_file.open(id + "_ip3.bin", std::ios::binary);
  cer_file.open(id + "_cer.bin", std::ios::binary);
}

cCell_x::~cCell_x() {
  ca_file.close();
  ip3_file.close();
  cer_file.close();
  out.close();
  delete solver;
  delete mesh;
}

void cCell_x::init_solvec(){
  out << "<Cell_x> initialising solution vector..." << std::endl;
  int np = mesh->vertices_count;

  solvec.resize(DIFVARS * np, 1); // NOTE: the variable ordering is c, ip, ce
  prev_solvec.resize(DIFVARS * np, 1);
  prev_solvec.block(0, 0, np, 1) = MatrixX1C().Constant(np, 1, p[c0]);
  prev_solvec.block(np, 0, np, 1) = MatrixX1C().Constant(np, 1, p[ip0]);
  prev_solvec.block(2 * np, 0, np, 1) = MatrixX1C().Constant(np, 1, p[ce0]);
  solvec = prev_solvec;

  nd_solvec.resize(NONDIFVARS * np, 1); // NOTE: the variable ordering is g, h
  prev_nd_solvec.resize(NONDIFVARS * np, 1);
  prev_nd_solvec.block(0, 0, np, 1) = MatrixX1C().Constant(np, 1, p[g0]);
  prev_nd_solvec.block(np, 0, np, 1) = MatrixX1C().Constant(np, 1, 0.0); // default to 0.0
  for(int n = 0; n < np; n++){
    if(node_data(n, BOOL_apical) == 1.0){  // apical nodes only
      prev_nd_solvec(np + n) = p[h0];
    }
  }
  nd_solvec = prev_nd_solvec;
}

ArrayRefMass cCell_x::make_ref_mass(){
  ArrayRefMass ref_mass;
  tCalcs v = (1.0 / 6.0) * 0.25 * 0.25;
  for(int i = 0; i < REF_MASS_SIZE; i++){
    for(int j = 0; j < REF_MASS_SIZE; j++){
      ref_mass(i, j) = v;
    }
  }
  return ref_mass;
}

void cCell_x::make_matrices(){

  out << "<Cell_x> id:" << id << " calculating the spatial factors..." << std::endl;

  int nt = mesh->tetrahedrons_count;
  element_data.resize(nt, Eigen::NoChange);
  int ns = mesh->surface_triangles_count;
  surface_data.resize(ns, Eigen::NoChange);
  int np = mesh->vertices_count;
  node_data.resize(np, Eigen::NoChange);

  // make the reference mass matrix
  ArrayRefMass ref_mass;
  ref_mass = make_ref_mass();

  // make the mass and stiffness matrices and the element constants matrix 
  out << "<Cell_x> calculating the constant matrices..." << std::endl;
  MatrixXXC stiffc, stiffp, small_mass;
  MatrixXXC stiffce;
  stiffc = stiffc.Zero(np, np);
  stiffp = stiffp.Zero(np, np);
  stiffce = stiffce.Zero(np, np);
  small_mass = small_mass.Zero(np, np);

  // --------------------------------
  // for each volume element...
  // --------------------------------
  for(int n = 0; n < (mesh->tetrahedrons_count); n++){ 
    Eigen::Matrix<int,1,4> vi;      // tetrahedron vertex indices
    vi = mesh->tetrahedrons.block<1,4>(n, 0);

    Eigen::Matrix<tCoord,4,3> vert; // tetrahedron vertex coordinates
    for(int i = 0; i < 4; i++)
      vert.block<1,3>(i, 0) = mesh->vertices.block<1,3>(int(vi(i)), 0); // why is typecast needed???

    Eigen::Matrix<tCoord,3,3> J;    // tetrahedron edge vectors
    for(int i = 0; i < 3; i++)
      J.block<1,3>(i, 0) = vert.block<1,3>(i + 1, 0) - vert.block<1,3>(0, 0);
    tCalcs V, Vx6;                  // tetrahedron volume, (6x) volume
    Vx6 = J.determinant();
    V = Vx6 / 6.0;
    element_data(n, VOL_e) = V;   // save the tetrahedron volume

    // RyR and PLC spatial factors per element
    element_data(n, RYR_e) = 
      ((mesh->dfa[n] < p[d_RyR]) ? (1.0 - ((p[d_RyR] - mesh->dfa[n]) / p[d_RyR])) : 1.0);
    element_data(n, PLC_e) = (mesh->dfb[n] < p[d_PLC]) ? 1.0 : 0.0;

    tCalcs Ic = V * p[Dc]; // diffusion coefficients
    tCalcs Ip = V * p[Dp];
    tCalcs Ice = V * p[De];

    Eigen::Matrix<tCoord,4,4> M, C, G;
    M.col(0) << 1, 1, 1, 1;
    M.block<4,3>(0, 1) = vert;
    C = M.inverse();
    G = C.block<3,4>(1, 0).transpose() * C.block<3,4>(1, 0); // gradients of the basis functions

    // construct the mass and stiffness matrix components
    for(int i = 0; i < 4; i++){
      stiffc(vi(i), vi(i)) += G(i, i) * Ic;
      stiffp(vi(i), vi(i)) += G(i, i) * Ip;
      stiffce(vi(i), vi(i)) += G(i, i) * Ice;
      small_mass(vi(i), vi(i)) += ref_mass(i, i) * Vx6;
    }
    for(int i = 0; i < 3; i++){
      stiffc(vi(0), vi(i + 1)) += G(0, i + 1) * Ic;
      stiffp(vi(0), vi(i + 1)) += G(0, i + 1) * Ip;
      stiffce(vi(0), vi(i + 1)) += G(0, i + 1) * Ice;
      small_mass(vi(0), vi(i + 1)) += ref_mass(0, i + 1) * Vx6;
      stiffc(vi(i + 1), vi(0)) = stiffc(vi(0), vi(i + 1));
      stiffp(vi(i + 1), vi(0)) = stiffp(vi(0), vi(i + 1));
      stiffce(vi(i + 1), vi(0)) = stiffce(vi(0), vi(i + 1));
      small_mass(vi(i + 1), vi(0)) = small_mass(vi(0), vi(i + 1));
    }
    for(int i = 0; i < 2; i++){
      stiffc(vi(1), vi(i + 2)) += G(1, i + 2) * Ic;
      stiffp(vi(1), vi(i + 2)) += G(1, i + 2) * Ip;
      stiffce(vi(1), vi(i + 2)) += G(1, i + 2) * Ice;
      small_mass(vi(1), vi(i + 2)) += ref_mass(1, i + 2) * Vx6;
      stiffc(vi(i + 2), vi(1)) = stiffc(vi(1), vi(i + 2));
      stiffp(vi(i + 2), vi(1)) = stiffp(vi(1), vi(i + 2));
      stiffce(vi(i + 2), vi(1)) = stiffce(vi(1), vi(i + 2));
      small_mass(vi(i + 2), vi(1)) = small_mass(vi(1), vi(i + 2));
    }
    stiffc(vi(2), vi(3)) += G(2, 3) * Ic;
    stiffp(vi(2), vi(3)) += G(2, 3) * Ip;
    stiffce(vi(2), vi(3)) += G(2, 3) * Ice;
    small_mass(vi(2), vi(3)) += ref_mass(2, 3) * Vx6;
    stiffc(vi(3), vi(2)) = stiffc(vi(2), vi(3));
    stiffp(vi(3), vi(2)) = stiffp(vi(2), vi(3));
    stiffce(vi(3), vi(2)) = stiffce(vi(2), vi(3));
    small_mass(vi(3), vi(2)) = small_mass(vi(2), vi(3));
  }

  // construct sparse mass matrix from a list of triplets (non zero elements)
  std::vector<Triplet> triplet_list;
  int np2 = np * 2;
  for (int j = 0; j < np; j++) {
    for (int i = 0; i < np; i++) {
      double v_ij = small_mass(i, j);
      if (v_ij != 0) { // add non zeros in first, second and third blocks
        triplet_list.push_back(Triplet(i, j, v_ij));
        triplet_list.push_back(Triplet(np + i, np + j, v_ij));
        triplet_list.push_back(Triplet(np2 + i, np2 + j, v_ij));
      }
    }
  }
  sparseMass.resize(DIFVARS * np, DIFVARS * np);
  sparseMass.setFromTriplets(triplet_list.begin(), triplet_list.end());

  // construct sparse stiffness matrix from list of triplets
  triplet_list.clear();
  for (int j = 0; j < np; j++) {
    for (int i = 0; i < np; i++) {
      double v_stiffc = stiffc(i, j); // first block
      if (v_stiffc != 0) {
        triplet_list.push_back(Triplet(i, j, v_stiffc));
      }
      double v_stiffp = stiffp(i, j); // second block
      if (v_stiffp != 0) {
        triplet_list.push_back(Triplet(np + i, np + j, v_stiffp));
      }
      double v_stiffce = stiffce(i, j); // third block
      if (v_stiffce != 0) {
        triplet_list.push_back(Triplet(np2 + i, np2 + j, v_stiffce));
      }
    }
  }
  sparseStiff.resize(DIFVARS * np, DIFVARS * np);
  sparseStiff.setFromTriplets(triplet_list.begin(), triplet_list.end());

  // set the A matrix size
  sparseA.resize(DIFVARS * np, DIFVARS * np);

  // --------------------------------
  // for each surface triangle...
  // --------------------------------
  for(int n = 0; n < ns; n++){ 
    Eigen::Array<int,1,3> vi; // surface triangle vertex indices
    vi = mesh->surface_triangles.block<1,3>(n, 0);
    Eigen::Array<tCoord,3,3> vert; // triangle vertex coordinates
    for(int i = 0; i < 3; i++)
      vert.block<1,3>(i, 0) = mesh->vertices.block<1,3>(int(vi(i)), 0); // why is typecast needed???
    Eigen::Matrix<tCoord,1,3> side1 = vert.block<1,3>(0, 0) - vert.block<1,3>(1, 0);
    Eigen::Matrix<tCoord,1,3> side2 = vert.block<1,3>(0, 0) - vert.block<1,3>(2, 0);
    tCalcs triangle_area = 0.5 * (side1.cross(side2)).norm();
    surface_data(n, AREA_s) = triangle_area; // save the triangle area
  }
  // --------------------------------
  // for each node...
  // --------------------------------
  for(int n = 0; n < np; n++){ 
    node_data(n, BOOL_apical) = 0.0; // default 
  }
  for(int n = 0; n < mesh->apical_triangles_count; n++){ // for each apical (surface) element...
    Eigen::Array<int,1,3> vi; // apical triangle vertex indices
    vi = mesh->surface_triangles.block<1,3>(mesh->apical_triangles(n), 0);
    for(int i = 0; i < 3; i++){ // for each apical triangle vertex
      node_data(vi(i), BOOL_apical) = 1.0; // flag it as apical
    }
  }
}

Array1VC cCell_x::get_apical_reactions(tCalcs c, tCalcs ip, tCalcs ce, tCalcs h){
  tCalcs phi_c = pow(c, 4) / (pow(p[K_c], 4) + pow(c, 4));
  tCalcs phi_p = pow(ip, 2) / (pow(p[K_p], 2) + pow(ip, 2));
  tCalcs h_alpha = pow(p[K_h], 4) / (pow(p[K_h], 4) + pow(c, 4));
  tCalcs beta  = phi_c * phi_p * h;
  tCalcs alpha = (1 - phi_p) * (1 - (phi_c * h_alpha));
  tCalcs po = beta / (beta + p[k_beta] * (beta + alpha));

  Array1VC reactions;
  reactions(0) = p[kIPR] * po * (ce - c);
  return reactions;
}

tCalcs cCell_x::get_g_reaction(tCalcs c, tCalcs g){ // RYR dynamics
  tCalcs ginf = pow(p[K_hRyR], 2) / (pow(p[K_hRyR], 2) + pow(c, 2));
  return (ginf - g) / p[tau];
}

tCalcs cCell_x::get_h_reaction(tCalcs c, tCalcs h){ // IPR dynamics
  tCalcs hinf = pow(p[K_h], 4) / (pow(p[K_h], 4) + pow(c, 4));
  tCalcs htau = p[tau_max] * pow(p[K_tau], 4) / (pow(p[K_tau], 4) + pow(c, 4));
  return (hinf - h) / htau;
}

Array1VC cCell_x::get_body_reactions(tCalcs c, tCalcs ip, tCalcs ce, tCalcs g, tCalcs ryr_f, tCalcs plc_f){
  tCalcs J_SERCA = p[V_p] * 
                   (pow(c, 2) - p[K_bar] * pow(ce, 2)) / (pow(p[k_p], 2) + pow(c, 2));
  tCalcs J_RYR   = ryr_f * p[V_RyR] * 
                   ( pow(c, p[n_RyR]) / (pow(c, p[n_RyR]) + pow(p[K_RyR], p[n_RyR]))) *
                   ( pow(ce, p[m_RyR]) / (pow(ce, p[m_RyR]) + pow(p[K_RyR2], p[m_RyR]))) * g;
  tCalcs vplc   = plc_f * p[V_PLC] * pow(c, 2) / (pow(p[K_PLC], 2) + pow(c, 2));
  tCalcs vdeg   = (p[V_5K] + (p[V_3K] * pow(c, 2))/(pow(p[K3K], 2) + pow(c, 2))) * ip;

  Array1VC reactions;
  reactions(0) = (J_RYR * (ce - c)) - J_SERCA;
  reactions(1) = vplc - vdeg;
  reactions(2) = -reactions(0) / p[Gamma];

  return reactions;
}

MatrixX1C cCell_x::make_load(tCalcs dt, bool plc){
  int np = mesh->vertices_count;
  ArrayX1C c, ip, g, h;
  ArrayX1C load_c, load_ip;
  ArrayX1C ce;
  ArrayX1C load_ce;
  MatrixX1C load;

  c = prev_solvec.block(0, 0, np, 1);
  ip = prev_solvec.block(np, 0, np, 1);
  ce = prev_solvec.block(2 * np, 0, np, 1);
  g = prev_nd_solvec.block(0, 0, np, 1);
  h = prev_nd_solvec.block(np, 0, np, 1); // note: only apical (surface) nodes used

  load_c = load_c.Zero(np, 1);
  load_ip = load_ip.Zero(np, 1);
  load_ce = load_ce.Zero(np, 1);
  load.resize(DIFVARS * np, Eigen::NoChange);

  // volume reaction terms
  for(int n = 0; n < (mesh->tetrahedrons_count); n++){ // for each volume element...
    Eigen::Array<int,1,4> vi;      // tetrahedron vertex indices
    vi = mesh->tetrahedrons.block<1,4>(n, 0);

    tCalcs cav  = 0.25 * (c(vi(0))  + c(vi(1))  + c(vi(2))  + c(vi(3)));
    tCalcs ipav = 0.25 * (ip(vi(0)) + ip(vi(1)) + ip(vi(2)) + ip(vi(3)));
    tCalcs ceav = 0.25 * (ce(vi(0)) + ce(vi(1)) + ce(vi(2)) + ce(vi(3)));
    tCalcs gav = 0.25 * (g(vi(0)) + g(vi(1)) + g(vi(2)) + g(vi(3)));

    Array1VC reactions = get_body_reactions(cav, ipav, ceav, gav,
        tCalcs(element_data(n, RYR_e)), tCalcs(plc ? element_data(n, PLC_e) : 0.0));

    for(int i = 0; i < 4; i++){ // for each tetrahedron vertex
      load_c(vi(i)) += element_data(n, VOL_e) * 0.25 * reactions(0); // reaction terms, scaled by 1/4 volume
      load_ip(vi(i)) += element_data(n, VOL_e) * 0.25 * reactions(1);
      load_ce(vi(i)) += element_data(n, VOL_e) * 0.25 * reactions(2);
    }
  }

  // apical (surface) reaction terms
  for(int n = 0; n < (mesh->apical_triangles_count); n++){ // for each apical (surface) triangle...
    Eigen::Array<int,1,3> vi; // apical triangle vertex indices
    vi = mesh->surface_triangles.block<1,3>(mesh->apical_triangles(n), 0);

    tCalcs cav  = (c(vi(0)) + c(vi(1)) + c(vi(2))) / 3.0;
    tCalcs ipav  = (ip(vi(0)) + ip(vi(1)) + ip(vi(2))) / 3.0;
    tCalcs ceav  = (ce(vi(0)) + ce(vi(1)) + ce(vi(2))) / 3.0;
    tCalcs hav  = (h(vi(0)) + h(vi(1)) + h(vi(2))) / 3.0;

    Array1VC reactions = get_apical_reactions(cav, ipav, ceav, hav);
    for(int i = 0; i < 3; i++){ // for each apical triangle vertex
      load_c(vi(i)) += // reaction term scaled by 1/3 area
        (surface_data(mesh->apical_triangles(n), AREA_s) / 3.0) * reactions(0); 
    }
  }

  // the diffusing variables
  load.block(0, 0, np, 1) = load_c;
  load.block(np, 0, np, 1) = load_ip;
  load.block(2 * np, 0, np, 1) = load_ce;

  return load;
}

MatrixX1C cCell_x::solve_nd(tCalcs dt){ // the non-diffusing variables
  int np = mesh->vertices_count;
  ArrayX1C c, g, h;
  MatrixX1C svec;

  c = prev_solvec.block(0, 0, np, 1);
  g = prev_nd_solvec.block(0, 0, np, 1);
  h = prev_nd_solvec.block(np, 0, np, 1); // note: only apical (surface) nodes used
  svec.resize(NONDIFVARS * np, Eigen::NoChange);

  for(int n = 0; n < np; n++){ // for each node...
    svec(n) = svec(n) + (dt * get_g_reaction(c(n), g(n))); // g
    if(node_data(n, BOOL_apical) == 1.0){ // only the apical nodes
      svec(np + n) = svec(np + n) + (dt * get_h_reaction(c(n), h(n))); // h
    }
  }
  return svec;
}

void cCell_x::exchange() {
  float* msg;
  MPI_Status stat;
  MPI_Request request;
  // send common face values to other cells (non-blocking)
  for(std::vector<cfc>::iterator it = cells.begin() ; it != cells.end(); ++it) {
    int mlength = C2CCOUNT * it->fcount;
    msg = new float[mlength];
    //...
    MPI_CHECK(MPI_Isend(msg, mlength, MPI_FLOAT, it->cell + 1, CELL_CELL_TAG, MPI_COMM_WORLD, &request));
    delete[] msg;
  }
  // receive common face values back from other cells
  for(std::vector<cfc>::iterator it = cells.begin() ; it != cells.end(); ++it) {
    int mlength = C2CCOUNT * it->fcount;
    msg = new float[mlength];
    MPI_CHECK(MPI_Recv(msg, mlength, MPI_FLOAT, it->cell + 1, CELL_CELL_TAG, MPI_COMM_WORLD, &stat));
    //...
    delete[] msg;
  }
}

void cCell_x::run() {
  float msg[ACCOUNT]; // mpi message from and back to acinus
  tCalcs delta_time, current_time;
  tCalcs prev_delta_time = 0.0;
  MPI_Status stat;
  struct timespec start, end;
  double elapsed;
  int np = mesh->vertices_count;

  MatrixX1C rhs; // the right-hand-side vector
  rhs.resize(DIFVARS * np, Eigen::NoChange);
  bool plc;

  int step = 0;
  while(true) {
    step++;
    prev_solvec = solvec;       // c, ip, ce
	prev_nd_solvec = nd_solvec; // g, h

    // get current time and time step value from acinus
    MPI_CHECK(MPI_Recv(&msg, ACCOUNT, MPI_FLOAT, acinus_rank, ACINUS_CELL_TAG, MPI_COMM_WORLD, &stat));
    delta_time = msg[dTime];
    current_time = msg[cTime];
    if(delta_time == 0.0) { // done?
      MPI_CHECK(MPI_Send(&msg, ACCOUNT, MPI_FLOAT, acinus_rank, ACINUS_CELL_TAG, MPI_COMM_WORLD));
      break;
    }

    // not done
    out << std::fixed << std::setprecision(3);
    out << "<Cell_x> step: " << step << " current_time: " << current_time << "s";
    out << " delta_time: " << delta_time << "s" << std::endl;
    plc = ((current_time >= p[PLCsrt]) and (current_time <= p[PLCfin])); // PLC on or off?

    if(delta_time != prev_delta_time) { // recalculate A matrix if time step changed
      sparseA = sparseMass + (delta_time * sparseStiff);
      ////solver->setA(sparseA);
	  solver->compute(sparseA);
      if(solver->info() != Eigen::Success) { // decomposition failed?
        utils::fatal_error("matrix decomposition failed", out);
      }
      prev_delta_time = delta_time;
    }

    // exchange common face data
    exchange();

    // calculate solution for diffusing variables
    rhs = (sparseMass * prev_solvec) + (delta_time * make_load(delta_time, plc));
    clock_gettime(CLOCK_REALTIME, &start);
    //solver->step(solvec, rhs);                   // VCL solver 
    solvec = solver->solve(rhs);                   // Eigen solver
    if(solver->info() != Eigen::Success) {
      utils::fatal_error("solver failed", out);;
	}
    clock_gettime(CLOCK_REALTIME, &end);
	elapsed = (end.tv_sec - start.tv_sec) + ((end.tv_nsec - start.tv_nsec) / 1000000000.0);
	out << std::fixed << std::setprecision(3);
	out << "<Cell_x> solver duration: " << elapsed << "s"<< std::endl;

    // check solver error and send it to acinus
    // ...
    msg[sError] = 0.0;
    MPI_CHECK(MPI_Send(&msg, ACCOUNT, MPI_FLOAT, acinus_rank, ACINUS_CELL_TAG, MPI_COMM_WORLD));

    // calculate solution for non-diffusing variables
    nd_solvec = solve_nd(delta_time);

    // save results
    if(step % int(p[Tstride]) == 0) {
      save_results(ca_file, 0);   // 0 = calcium
      save_results(ip3_file, 1);  // 1 = ip3
      save_results(cer_file, 2);  // 2 = cer
    }
  }
}

void cCell_x::save_results(std::ofstream &data_file, int var){
  int np = mesh->vertices_count;
  float* fbuf = new float[np];
  for(int n=0; n<np; n++) fbuf[n] = solvec[var*np + n]; // convert to float for reduced file size
  data_file.write(reinterpret_cast<char*>(fbuf), np * sizeof(float));
}

