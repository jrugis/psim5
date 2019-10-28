/*
 * utils.hpp
 *
 *  Created on: 27/04/2018
 *      Author: jrugis
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <fstream>
#include <iostream>
#include <string>

#include "global_defs.hpp"

namespace utils {
  void calc_tri_centers(MatrixN3d& centers, const MatrixN3d& vertices, const MatrixN3i& triangles);
  void fatal_error(const std::string msg, std::ofstream& out);
  void get_parameters(const std::string file_id, int ptype, int cell_num, double* p, std::ofstream& out);
  double get_distance(const Vector3d& p, const Vector3d& v, const Vector3d& w);
  void read_mesh(const std::string file_name, sMeshVals& mesh_vals, std::ofstream& out);
  void save_matrix(const std::string file_name, const MatrixNNd& mat);
  void save_integer_matrix(const std::string file_name, const MatrixNNi& mat);
} // namespace utils

#endif /* UTILS_H_ */
