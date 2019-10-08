/*
 * cLumenBase.hpp
 *
 *  Created on: 08/10/2019
 *      Author: jrugis
 */

#ifndef CLUMENBASE_H_
#define CLUMENBASE_H_

#include <string>

#include "global_defs.hpp"
#include "cCell_calcium.hpp"

class cLumenBase {
public:
  cLumenBase(cCell_calcium* p);
  ~cLumenBase();
  
  tDist get_dnl(Eigen::Vector3d p);

private:
  std::string id;
  cCell_calcium* parent;
  int points_count, segments_count; // the number of points and segments in the lumen tree
  Eigen::Array<tCoord, Eigen::Dynamic, 3, Eigen::RowMajorBit> points; // 3x coordinate
  Eigen::Array<int, Eigen::Dynamic, 2, Eigen::RowMajorBit> segments; // line segment indices, 2x points
  
  void get_segments();
  void print_info();
};

#endif /* CLUMENBASE_ */

