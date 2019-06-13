/*
 * utils.cpp
 *
 *  Created on: 27/04/2018
 *      Author: jrugis
 */

#include <mpi.h>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>

#include "global_defs.hpp"
#include "utils.hpp"

// NOTE: outputs error message to stderr and the process ".out" file
void utils::fatal_error(const std::string msg, std::ofstream& out){
  std::string m = "ERROR: " + msg;
  out << m << std::endl; out.close();
  std::cerr << m << std::endl;
  MPI_CHECK(MPI_Abort(MPI_COMM_WORLD, 1));  // make sure the whole program fails
  exit(1);
}

// NOTE: used by each of the acinus, lumen and cell objects
void utils::get_parameters(const std::string file_id, int ptype, int cell_num, tCalcs* p, std::ofstream& out){
  std::string file_name = file_id + ".dat";
  std::ifstream model_file(file_name); // open the model parameters file
  std::string line;                    // file line buffer
  std::string *pnames;                 // pointer to vector of paramemter names
  std::vector <std::string> tokens;    // tokenized line
  int pcount;

  // calcium simulation parameters
  std::string cpnames[PCOUNT] = { \
    "delT", "totalT", "Tstride", \
    "PLCsrt", "PLCfin", \
    "c0", "ip0", "ce0", "Gamma", \
    "Dc", "Dp", "De", "Fc", "Fip", \
    "d_RyR", "V_RyR", "K_RyR", "K_RyR2", "m_RyR", "n_RyR", \
    "k_beta", "K_p", "K_c", "K_h", "kIPR", \
    "V_p", "k_p", "K_bar", \
    "d_PLC", "V_3K", "V_5K", "K_PLC", "K3K", "V_PLC", \
    "h0", "K_tau", "tau_max", \
    "g0", "K_hRyR", "tau"}; // NOTE: these must match up with the enums in global_defs.hpp !!!

  // fluid flow parameters
  std::string fpnames[FPCOUNT] = { \
    "odeSolver", "odeSolverAbsTol", "odeSolverRelTol", \
    "aNkcc1", "a1", "a2", "a3", "a4", \
    "r", "alpha1", "aNaK", \
    "GtNa", "GtK", \
    "GCl", "KCaCC", "eta1", \
    "GK", "KCaKC", "eta2", \
    "G1", "KNa", "KH", \
    "G4", "KCl", "KB", \
    "GB", "kn", "kp", \
    "pHl", "pHi", "pHe", "HCO3l", "CO20", "Ul", "Cle", "Nae", \
    "Ke", "HCO3e", "CO2e", "CO2l", "Hy", \
    "La", "Lb", "Lt", \
    "He", "Ie", "Hye", "St", "wl"}; // NOTE: these must match up with the enums in global_defs.hpp !!!

  if(ptype == calciumParms) {pnames = cpnames; pcount = PCOUNT;} // calcium simulation?
  if(ptype == flowParms) {pnames = fpnames; pcount = FPCOUNT;}   // fluid flow?

  if (not model_file.is_open()) {
    fatal_error("the model parameters file " + file_name + " could not be opened", out);
  }
  
  out << "<utils> reading model parameters..." << std::endl;
  for(int n = 0; n < pcount; n++) p[n] = tCalcs(-1.0);  // not-hit marker
  while(getline(model_file, line)){
    if(line.data()[0] == '#') continue;
    int ci = line.find_first_of("#");     // remove comment, if any
	if(ci > 0) line = line.substr(0, ci); //
	line = boost::trim_right_copy(line);  // remove trailing whitespace
    boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
	bool found = false;
    for(int n = 0; n < pcount; n++) {
      if(tokens[0] == pnames[n]) {
        p[n] = atof(tokens[ ((tokens.size() == 2) ? 1 : cell_num) ].c_str());
		found = true; break;
      }
    }
	if(!found) fatal_error("invalid parameter: " + tokens[0], out);
  }
  model_file.close();
  for(int n = 0; n < pcount; n++)
    if(p[n] < 0.0) fatal_error("missing parameter: " + pnames[n], out);
}

