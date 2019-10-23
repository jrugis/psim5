/*
 * utils.cpp
 *
 *	Created on: 27/04/2018
 *	Author: jrugis
 */

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/tokenizer.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <fstream>
#include <iostream>

#include "global_defs.hpp"
#include "utils.hpp"

// get the center points of triangles
void utils::calc_tri_centers(MatrixN3d& centers, const MatrixN3d& vertices, const MatrixN3i& triangles)
{
  centers.resize(triangles.rows(), Eigen::NoChange);
  for (int n = 0; n < triangles.rows(); n++) {
    Vector3d v1 = Vector3d(vertices.row((triangles)(n, 0)));
    Vector3d v2 = Vector3d(vertices.row((triangles)(n, 1)));
    Vector3d v3 = Vector3d(vertices.row((triangles)(n, 2)));
    centers.row(n) = (v1 + v2 + v3) / 3.0;
  }
}

// NOTE: outputs error message to stderr and the process ".out" file
void utils::fatal_error(const std::string msg, std::ofstream& out)
{
  std::string m = "ERROR: " + msg;
  out << m << std::endl;
  out.close();
  std::cerr << m << std::endl;
  exit(1);
}

// NOTE: used by each of the acinus, lumen and cell objects
void utils::get_parameters(const std::string file_id, int ptype, int cell_num, double* p, std::ofstream& out)
{
  std::string file_name = file_id + ".dat";
  std::ifstream model_file(file_name); // open the model parameters file
  std::string line;                    // file line buffer
  std::string* pnames;                 // pointer to vector of paramemter names
  std::vector<std::string> tokens;     // tokenized line
  int pcount;

  // calcium simulation parameters
  // NOTE: these must match up with the enums in global_defs.hpp !!!
  std::string cpnames[PCOUNT] = {
    "delT",   "totalT", "Tstride", "PLCsrt", "PLCfin", "c0",    "ip0",     "ce0",   "Gamma",  "Dc",    "Dp",
    "De",     "Fc",     "Fip",     "d_RyR",  "V_RyR",  "K_RyR", "K_RyR2",  "m_RyR", "n_RyR",  "IPRdn", "IPRdf",
    "k_beta", "K_p",    "K_c",     "K_h",    "kIPR",   "V_p",   "k_p",     "K_bar", "PLCds",  "PLCdl", "V_3K",
    "V_5K",   "K_PLC",  "K3K",     "V_PLC",  "h0",     "K_tau", "tau_max", "g0",    "K_hRyR", "tau"};

  // fluid flow parameters
  // NOTE: these must match up with the enums in global_defs.hpp !!!
  std::string fpnames[FPCOUNT] = {
    "aNkcc1", "a1",   "a2", "a3",  "a4",  "r",  "alpha1", "aNaK", "GtNa", "GtK",  "GCl", "KCaCC", "eta1", "GK",
    "KCaKC",  "eta2", "G1", "KNa", "KH",  "G4", "KCl",    "KB",   "GB",   "kn",   "kp",  "pHl",   "pHi",  "pHe",
    "HCO3l",  "CO20", "Ul", "Cle", "Nae", "Ke", "HCO3e",  "CO2e", "Hl",   "CO2l", "Hy",  "La",    "Lb",   "Lt"};

  if (ptype == calciumParms) {
    pnames = cpnames;
    pcount = PCOUNT;
  } // calcium simulation?
  if (ptype == flowParms) {
    pnames = fpnames;
    pcount = FPCOUNT;
  } // fluid flow?

  if (not model_file.is_open()) { fatal_error("the model parameters file " + file_name + " could not be opened", out); }

  out << "<utils> reading model parameters..." << std::endl;
  for (int n = 0; n < pcount; n++) p[n] = double(-1.0); // not-hit marker
  while (getline(model_file, line)) {
    if (line.data()[0] == '#') continue;
    int ci = line.find_first_of("#");      // remove comment, if any
    if (ci > 0) line = line.substr(0, ci); //
    line = boost::trim_right_copy(line);   // remove trailing whitespace
    boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
    bool found = false;
    for (int n = 0; n < pcount; n++) {
      if (tokens[0] == pnames[n]) {
        p[n] = atof(tokens[((tokens.size() == 2) ? 1 : cell_num)].c_str());
        found = true;
        break;
      }
    }
    if (!found) fatal_error("invalid parameter: " + tokens[0], out);
  }
  model_file.close();
  for (int n = 0; n < pcount; n++)
    if (p[n] < 0.0) fatal_error("missing parameter: " + pnames[n], out);
}

double utils::get_distance(const Vector3d& p, const Vector3d& v, const Vector3d& w)
{
  // Return minimum distance between line segment vw and point p
  double l2 = (w - v).squaredNorm();      // |w-v|^2		avoid a sqrt
  if (l2 == 0.0) return ((v - p).norm()); // v == w case, return distance(p, v)
  // Consider the line extending the segment, parameterized as v + t (w - v).
  // Find projection of point p onto the line. It falls where t = [(p-v) . (w-v)] / |w-v|^2
  // Clamp t from [0,1] to handle points outside the segment vw.
  double t = std::max(0.0, std::min(1.0, (p - v).dot(w - v) / l2)); // max(0, min(1, dot(p - v, w - v) / l2));
  const Vector3d projection = v + (t * (w - v));                    // Projection falls on the segment
  return ((projection - p).norm());                                 // return distance(p, projection)
}

void utils::save_matrix(std::string file_name, const MatrixNNd& mat)
{
  int rows = mat.rows();
  int cols = mat.cols();
  int rc = rows * cols;
  double* buf = new double[rc]; // temporary buffer
  int n = 0;
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < cols; c++) {
      buf[n++] = mat(r, c); // flatten the data
    }
  }
  std::ofstream data_file;
  data_file.open(file_name, std::ios::binary);
  data_file.write(reinterpret_cast<char*>(buf), rc * sizeof(double));
  data_file.close();
  delete (buf);
}

void utils::save_integer_matrix(std::string file_name, const MatrixNNi& mat)
{
  int rows = mat.rows();
  int cols = mat.cols();
  int rc = rows * cols;
  int* buf = new int[rc]; // temporary buffer
  int n = 0;
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < cols; c++) {
      buf[n++] = mat(r, c); // flatten the data
    }
  }
  std::ofstream data_file;
  data_file.open(file_name, std::ios::binary);
  data_file.write(reinterpret_cast<char*>(buf), rc * sizeof(int));
  data_file.close();
  delete (buf);
}
