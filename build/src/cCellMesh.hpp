/*
 * cCellMesh.hpp
 *
 *	Created on: 30/03/2018
 *  Author: jrugis
 */

#ifndef CCELLMESH_H_
#define CCELLMESH_H_

#include <Eigen/Dense>
#include <string>

#include "global_defs.hpp"

class cCell_calcium;

class cCellMesh {
  public:
  cCellMesh(std::string mesh_name, cCell_calcium* parent);
  ~cCellMesh();
  void mesh_calcs();
  void print_info();

  int vertices_count, tetrahedrons_count;
  int surface_triangles_count, apical_triangles_count, basal_triangles_count;
  int common_triangles_count;
  MatrixN3d vertices;          // 3x coordinate
  MatrixN3i surface_triangles; // 3x vertex
  MatrixN4i tetrahedrons;      // 4x vertex
  MatrixNCi common_triangles;  // this triangle, other cell, other triangle
  MatrixN1i apical_triangles;  // surface triangle indicies
  MatrixN1i basal_triangles;   // surface triangle indicies
  MatrixN1d n_dfa;             // distance from apical (per node)
  MatrixN1d e_dfa;             // distance from apical (per element)
  MatrixN1d e_dfb;             // distance from basal (per element)

  private:
  std::string id;
  cCell_calcium* parent;
  void get_mesh(std::string file_name);
  void calc_common();
  void calc_dfa();
  void calc_apical_basal();
  void calc_dfb();
};

#endif /* CCELLMESH_H_ */
