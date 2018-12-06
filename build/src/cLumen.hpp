/*
 * cLumen.hpp
 *
 *  Created on: 06/12/2018
 *      Author: jrugis
 */

#ifndef CLUMEN_H_
#define CLUMEN_H_

#include <fstream>
#include <string>

#include "global_defs.hpp"

class cLumen {
public:
  cLumen(std::string host_name, int my_rank, int cell_rank, int cell_count, int acinus_rank);
  ~cLumen();
  void run();

private:
  std::string id;
  std::ofstream out;
  tCalcs p[FPCOUNT]; // the fluid flow parameters array
  int my_rank, cell_rank, cell_count, acinus_rank;

  tCalcs snd_recv(tCalcs t, tCalcs dt);
};

#endif /* CLUMEN_ */

