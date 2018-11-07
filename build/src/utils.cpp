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
  std::string pnames[PCOUNT] = { \
    "delT", "totalT", "Tstride", \
    "PLCsrt", "PLCfin", \
    "c0", "ip0", "ce0", "Gamma", \
    "Dc", "Dp", "De", "Fc", "Fip", \
    "d_RyR", "w_RYR", "V_RyR", "K_RyR", "K_RyR2", "m_RyR", "n_RyR", \
    "k_beta", "K_p", "K_c", "K_h", "kIPR", \
    "V_p", "k_p", "K_bar", \
    "d_PLC", "w_PLC", "V_3K", "V_5K", "K_PLC", "K3K", "VPLC", \
    "h0", "K_tau", "tau_max", \
    "g0", "K_hRyR", "tau"}; // NOTE: these must match up with the enums in global_defs.hpp !!!

  if (not model_file.is_open()) {
    fatal_error("the model parameters file " + file_name + " could not be opened", out);
  }
  out << "<utils> reading model parameters..." << std::endl;
  for(int n = 0; n < PCOUNT; n++) p[n] = tCalcs(-1.0);  // not-hit marker
  while(getline(model_file, line)){
    if(line.data()[0] == '#') continue;
    boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
	bool found = false;
    for(int n = 0; n < PCOUNT; n++) {
      if(tokens[0] == pnames[n]) {
        p[n] = atof(tokens[1].c_str());
		found = true; break;
      }
    }
	if(!found) fatal_error("invalid parameter: " + tokens[0], out);
  }
  model_file.close();
  for(int n = 0; n < PCOUNT; n++)
    if(p[n] < 0.0) fatal_error("missing parameter: " + pnames[n], out);
}

