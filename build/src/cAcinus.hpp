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

class cAcinus {
public:
  cAcinus(std::string host_name, int my_rank, int cell_rank, int cell_count, int lumen_rank);
  ~cAcinus();
  void run();

private:
  std::string id;
  std::ofstream out;
  tCalcs p[PCOUNT]; // the calcium model parameters array
  int my_rank, cell_rank, cell_count, lumen_rank;

  tCalcs snd_recv(tCalcs t, tCalcs dt);
};

#endif /* CACINUS_H_ */

