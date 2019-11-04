/*
 * cAcinus.hpp
 *
 *	Created on: 26/04/2018
 *	Author: jrugis
 */

#ifndef CACINUS_H_
#define CACINUS_H_

#include <fstream>
#include <string>

#include "cLumen.hpp"
#include "global_defs.hpp"

class cAcinus {
  public:
  // cAcinus(std::string host_name, int my_rank, int cell_rank, int cell_count, int lumen_rank);
  cAcinus(const std::string host_name, int my_rank, int cell_rank, int cell_count);
  ~cAcinus();
  void run();

  private:
  std::string id;
  std::ofstream out;
  double p[PCOUNT]; // the calcium model parameters array
  // int my_rank, cell_rank, cell_count, lumen_rank;
  int my_rank, cell_rank, cell_count;
  cLumen* lumen;

  void snd(double t, double dt);
  double recv();
};

#endif /* CACINUS_H_ */
