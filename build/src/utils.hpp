/*
 * utils.hpp
 *
 *  Created on: 27/04/2018
 *      Author: jrugis
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <iostream>
#include <fstream>
#include <string>

#include "global_defs.hpp"

namespace utils
{
  void fatal_error(const std::string msg, std::ofstream& out);
  void get_parameters(const std::string file_id, int ptype, int cell_num, tCalcs* p, std::ofstream& out);
  tDist get_distance(Eigen::Vector3d p, Eigen::Vector3d v, Eigen::Vector3d w);
  void save_matrix(std::string file_name, MatrixXXC mat);
  void save_integer_matrix(std::string file_name, MatrixXXI mat);
}
#endif /* UTILS_H_ */
