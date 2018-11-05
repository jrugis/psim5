/*
 * utils.cpp
 *
 *  Created on: 27/04/2018
 *      Author: jrugis
 */

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

#include "global_defs.hpp"
#include "utils.hpp"

// NOTE: outputs error message to stderr and the process ".out" file
void utils::fatal_error(const std::string msg, std::ofstream& out){
  std::string m = "ERROR: " + msg;
  out << m << std::endl; out.close();
  std::cerr << m << std::endl;
  exit(1);
}

// NOTE: used by both the acinus and cell objects
void utils::get_parameters(const std::string file_id, tCalcs* p, std::ofstream& out){
  std::string file_name = file_id + ".dat";
  std::ifstream model_file(file_name); // open the model parameters file
  std::string line;                    // file line buffer
  std::vector <std::string> tokens;    // tokenized line

  // check the file is open
  if (not model_file.is_open()) {
    fatal_error("the model parameters file " + file_name + " could not be opened", out);
  }
  out << "<utils> reading model parameters..." << std::endl;
  int n = 0;   // read in the model parameters
  while(getline(model_file, line)){
    if(line.data()[0] == '%') continue;
    boost::split(tokens, line, boost::is_any_of(", "), boost::token_compress_on);
    if((n + tokens.size()) > PCOUNT) fatal_error("too many parameters in the model parameters file", out);
    for(unsigned int m = 0; m < tokens.size(); m++) p[n++] = atof(tokens[m].c_str());
    }
  model_file.close();
  if(n != PCOUNT) fatal_error("too few parameters in the model parameters file", out);
}

/*
void utils::save_matrix(std::string file_name, MatrixXXC mat){
  std::ofstream file(file_name.c_str(), std::ios::binary); // create the file
  tElement rows = mat.rows();
  tElement cols = mat.cols();
  float f;
  int esize = sizeof(f);
  file.write(reinterpret_cast<char*>(&rows), sizeof(rows));
  file.write(reinterpret_cast<char*>(&cols), sizeof(cols));
  for(int j = 0; j < cols; j++){
    for(int i = 0; i < rows; i++){
      f = mat(i, j);  // convert to float for smaller file size
      file.write(reinterpret_cast<char*>(&f), esize); // column order
    }
  }
  file.close();
}
*/

