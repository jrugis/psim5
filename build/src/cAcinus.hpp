/*
 * cAcinus.hpp
 *
 *  Created on: 26/04/2018
 *      Author: jrugis
 */

#ifndef CACINUS_H_
#define CACINUS_H_

#include <fstream>
#include <string>

#include "global_defs.hpp"
#include "cLumen.hpp"

class cAcinus {
public:
  //cAcinus(std::string host_name, int my_rank, int cell_rank, int cell_count, int lumen_rank);
  cAcinus(std::string host_name, int my_rank, int cell_rank, int cell_count);
  ~cAcinus();
  void run();

private:
  std::string id;
  std::ofstream out;
  tCalcs p[PCOUNT]; // the calcium model parameters array
  //int my_rank, cell_rank, cell_count, lumen_rank;
  int my_rank, cell_rank, cell_count;
  cLumen* lumen;

  void snd(tCalcs t, tCalcs dt);
  tCalcs recv();
};

#endif /* CACINUS_H_ */

