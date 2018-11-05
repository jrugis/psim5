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
  void get_parameters(const std::string file_id, tCalcs* p, std::ofstream& out);
//  void save_matrix(std::string file_name, MatrixXXC mat);
}

#endif /* UTILS_H_ */
