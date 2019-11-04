/*
 * cLumenTree.hpp
 *
 *	Created on: 08/10/2019
 *			Author: jrugis
 */

#ifndef CLUMENTREE_H_
#define CLUMENTREE_H_

#include <string>

#include "global_defs.hpp"

class cLumenTree {
  public:
  cLumenTree(std::ofstream& out);
  ~cLumenTree();

  double get_dnl(const Eigen::Vector3d p);

  private:
  std::string id;
  std::ofstream* out;
  int points_count, segments_count; // the number of points and segments in the lumen tree
  MatrixN3d points;                 // 3x coordinate
  MatrixN2i segments;               // line segment indices, 2x points

  void get_segments();
  void print_info();
};

#endif /* CLUMENTREE_ */
