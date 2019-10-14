/*
 * cCellMesh.hpp
 *
 *  Created on: 30/03/2018
 *      Author: jrugis
 */

#ifndef CCELLMESH_H_
#define CCELLMESH_H_

#include <string>
#include <fstream>
#include <Eigen/Dense>

#include "global_defs.hpp"

class cCellMesh {
public:
  cCellMesh(std::string mesh_name, std::ofstream& out, int cell_number);
  ~cCellMesh();
  void print_info(std::ofstream& out);

  int vertices_count, tetrahedrons_count;
  int surface_triangles_count, apical_triangles_count, basal_triangles_count;
  int common_triangles_count;
  Eigen::Array<tCoord, Eigen::Dynamic, 3, Eigen::RowMajorBit> vertices; // 3x coordinate
  Eigen::Array<int, Eigen::Dynamic, 3, Eigen::RowMajorBit> surface_triangles; // 3x vertex
  Eigen::Array<int, Eigen::Dynamic, 4, Eigen::RowMajorBit> tetrahedrons; // 4x vertex
  Eigen::Array<int, Eigen::Dynamic, CCONNCOUNT, Eigen::RowMajorBit> common_triangles; // this triangle, other cell, other triangle
  Eigen::Array<int, Eigen::Dynamic, CCONNCOUNT, Eigen::RowMajorBit> common_apical_triangles; // this triangle, other cell, other triangle
  Eigen::Array<int, Eigen::Dynamic, 1> apical_triangles; // surface triangle indicies
  Eigen::Array<int, Eigen::Dynamic, 1> basal_triangles; // surface triangle indicies
  Eigen::Array<tDist, Eigen::Dynamic, 1> dfa; // distance from apical (per element)
  Eigen::Array<tDist, Eigen::Dynamic, 1> dfb; // distance from basal (per element)

private:
  std::string id;
  void get_mesh(std::string file_name, std::ofstream& out, int cell_number);
  void identify_common_apical_triangles(std::ofstream& out, int cell_number);
};

#endif /* CCELLMESH_H_ */
