/*
 * cCell_flow.hpp
 *
 *	Created on: 16/12/2018
 *			Author: jrugis
 */

#ifndef CCELL_FLOW_H_
#define CCELL_FLOW_H_

#include <fstream>
#include <string>

#include "global_defs.hpp"

class cCell_flow {
  public:
  cCell_flow(int my_index, double p[], std::ofstream& out);
  ~cCell_flow();
  void run();

  private:
  double* p; // the fluid flow parameters array
  int my_index;
  double omega; // volume
  double na;    // sodium
  double nk;    // potassium
  double cl;    // chloride
  double hc03;  // bicarbonate
  double h;     // hydrogen

  double snd_recv(double t, double dt);
};

#endif /* CCELL_FLOW_H_ */
