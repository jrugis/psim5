/*
 * cCellMesh.cpp
 *
 *	Created on: 26/04/2018
 *	Author: jrugis
 */

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "cCellMesh.hpp"
#include "cCell_calcium.hpp"
#include "cLumenTree.hpp"
#include "utils.hpp"

cCellMesh::cCellMesh(const std::string mesh_name, cCell_calcium* p)
{
  // initialise member variables
  parent = p;
  id = mesh_name;
  parent->out << "<CellMesh> reading mesh file " + id << std::endl;
  utils::read_mesh(id, vertices_count, vertices, surface_triangles_count, surface_triangles, tetrahedrons_count,
                   tetrahedrons, parent->out);
  mesh_calcs();
}

cCellMesh::~cCellMesh() {}

void cCellMesh::calc_dfa()
{
  // calculate the distance from nodes to the lumen
  cLumenTree lumen(parent->out);
  n_dfa.resize(vertices_count, Eigen::NoChange);
  for (int n = 0; n < vertices_count; n++) { n_dfa(n) = lumen.get_dnl(vertices.row(n)); }
  // calculate the distance from elements to the lumen
  e_dfa.resize(tetrahedrons_count, Eigen::NoChange);
  for (int n = 0; n < tetrahedrons_count; n++) {
    e_dfa(n) =
      (n_dfa(tetrahedrons(n, 0)) + n_dfa(tetrahedrons(n, 1)) + n_dfa(tetrahedrons(n, 2)) + n_dfa(tetrahedrons(n, 3))) / 4.0;
  }
}

// determine the cell-to-cell connetivity data (this_triamgle, other_cell, other_triangle)
//	 do this by comparing surface triangle centers on this cell to those on every other cell
void cCellMesh::calc_common()
{
  MatrixN3d this_centers;  // this cell surface triangle centers
  MatrixN3d other_centers; // other cell surface triangle centers
  MatrixN3d other_vertices;
  int other_surface_triangles_count;
  MatrixN3i other_surface_triangles;
  MatrixN4i other_tetraedrons;                        // required for reading mesh files but not otherwise used
  int other_vertices_count, other_tetrahedrons_count; // ditto

  utils::calc_tri_centers(this_centers, vertices, surface_triangles);
  int ci = id.find_last_of('c') + 1;               // to split the mesh id string, find the position of the last "c"
  int this_cell = atoi(id.substr(ci).c_str()) - 1; // extract this cell index from its mesh id
  common_triangles.resize(surface_triangles_count, Eigen::NoChange); // overkill size for now
  common_triangles_count = 0;
  for (int other_cell = 0; other_cell < CELLS_COUNT; other_cell++) { // check against all the other cell meshes...
    if (other_cell == this_cell) continue;                           // don't check against itself
    std::string other_id = id.substr(0, ci) + std::to_string(other_cell + 1); // the other mesh id
    utils::read_mesh(other_id, other_vertices_count, other_vertices, other_surface_triangles_count,
                     other_surface_triangles, other_tetrahedrons_count, other_tetraedrons, parent->out);
    utils::calc_tri_centers(other_centers, other_vertices, other_surface_triangles);
    for (int this_tri = 0; this_tri < surface_triangles_count; this_tri++) {
      for (int other_tri = 0; other_tri < other_surface_triangles_count; other_tri++) {
        Eigen::Vector3d p1 = this_centers.row(this_tri);
        Eigen::Vector3d p2 = other_centers.row(other_tri);
        if (p1 == p2) {
          common_triangles.row(common_triangles_count++) = Eigen::Vector3i(this_tri, other_cell, other_tri);
        }
      }
    }
  }
  common_triangles.conservativeResize(common_triangles_count, Eigen::NoChange); // actual size
}

void cCellMesh::calc_apical_basal()
{
  // determine the apical and the apical keep-out surface triangle indices
  MatrixN1i apical_keepout;                                        // keep-out triangle indicies
  apical_keepout.resize(surface_triangles_count, Eigen::NoChange); // overkill size for now
  int apical_keepout_count = 0;
  apical_triangles.resize(surface_triangles_count, Eigen::NoChange); // overkill size for now
  apical_triangles_count = 0;
  for (int n = 0; n < surface_triangles_count; n++) {
    double d = // triangle distance from apical
      (n_dfa(surface_triangles(n, 0)) + n_dfa(surface_triangles(n, 1)) + n_dfa(surface_triangles(n, 2))) / 3.0;
    if (d < parent->p[IPRdn]) apical_triangles(apical_triangles_count++) = n;
    if (d < (parent->p[PLCdl])) apical_keepout(apical_keepout_count++) = n;
  }
  apical_triangles.conservativeResize(apical_triangles_count, 1); // actual triangles count
  apical_keepout.conservativeResize(apical_keepout_count, 1);     // actual triangles count

  // determine the basal triangle indices by considering all surface triangles
  //	 then eliminating the common triangles and the triangles that are too close to the lumen
  basal_triangles.resize(surface_triangles_count, Eigen::NoChange); // overkill size for now
  basal_triangles.setOnes(surface_triangles_count);
  for (int n = 0; n < common_triangles_count; n++) { // eliminate the common triangles
    basal_triangles(common_triangles(n, 0)) = 0;
  }
  for (int n = 0; n < apical_keepout_count; n++) { // eliminate the apical keepout triangles
    basal_triangles(apical_keepout(n)) = 0;
  }
  basal_triangles_count = 0;
  for (int n = 0; n < surface_triangles_count; n++) { // convert what's left to a list of indices
    if (basal_triangles(n)) { basal_triangles(basal_triangles_count++) = n; }
  }
  basal_triangles.conservativeResize(basal_triangles_count, Eigen::NoChange); // actual size
}

// calculate the distance from elements to basal surface
//	 use the average of the element-vertex to nearest-basal-triangle-vertex distances
void cCellMesh::calc_dfb()
{
  // get the basal vertices
  MatrixN1i basal_verts;                               // the basal triangle vertices
  basal_verts.resize(vertices_count, Eigen::NoChange); // overkill size for now
  basal_verts.setZero(vertices_count);
  for (int n = 0; n < basal_triangles_count; n++) {
    Eigen::Vector3i vi = Eigen::Vector3i(surface_triangles.row(basal_triangles(n)));
    for (int i = 0; i < 3; i++) basal_verts(vi(i)) = 1; // flag vertex as basal
  }
  int basal_verts_count = 0;
  for (int n = 0; n < vertices_count; n++) { // convert flags to a list of basal indices
    if (basal_verts(n)) basal_verts(basal_verts_count++) = n;
  }
  basal_verts.conservativeResize(basal_verts_count, Eigen::NoChange); // actual size
  // calculate the (per node) node to nearest basal node distance
  MatrixN1d n_dfb;
  n_dfb.resize(vertices_count, Eigen::NoChange);
  for (int n = 0; n < vertices_count; n++) {
    Eigen::Vector3d p1 = Eigen::Vector3d(vertices.row(n));
    n_dfb(n) = 1000.0; // large dummy value
    for (int m = 0; m < basal_verts_count; m++) {
      Eigen::Vector3d p2 = Eigen::Vector3d(vertices.row(basal_verts(m)));
      double d = (p1 - p2).norm();
      if (d < n_dfb(n)) n_dfb(n) = d;
    }
  }
  // for each tet calculate e_dnb as the average it's vertex dnb's
  e_dfb.resize(tetrahedrons_count, Eigen::NoChange);
  for (int n = 0; n < tetrahedrons_count; n++) {
    e_dfb(n) =
      (n_dfb(tetrahedrons(n, 0)) + n_dfb(tetrahedrons(n, 1)) + n_dfb(tetrahedrons(n, 2)) + n_dfb(tetrahedrons(n, 3))) / 4.0;
  }
}

void cCellMesh::mesh_calcs()
{
  parent->out << "<CellMesh> calculating derived properties... " << std::endl;
  calc_common();
  calc_dfa();          // do this first,				distance from apical (per node and per element)
  calc_apical_basal(); // then this,						identify apical and basal triangles
  calc_dfb();          // finally this.				distance from basal (per element)
}

void cCellMesh::print_info()
{
  parent->out << "<CellMesh> vertices_count: " << vertices_count << std::endl;
  parent->out << "<CellMesh> tetrahedrons_count: " << tetrahedrons_count << std::endl;
  parent->out << "<CellMesh> surface_triangles_count: " << surface_triangles_count << std::endl;
  parent->out << "<CellMesh> apical_triangles_count: " << apical_triangles_count << std::endl;
  parent->out << "<CellMesh> basal_triangles_count: " << basal_triangles_count << std::endl;
  parent->out << "<CellMesh> common_triangles_count: " << common_triangles_count << std::endl;

  // ***********************************************************
  utils::save_matrix("vertices_" + id + ".bin", vertices);
  utils::save_integer_matrix("triangles_" + id + ".bin", surface_triangles);
  utils::save_integer_matrix("tetrahedrons_" + id + ".bin", tetrahedrons);
  utils::save_integer_matrix("apical_" + id + ".bin", apical_triangles);
  utils::save_matrix("n_dfa_" + id + ".bin", n_dfa);
  utils::save_matrix("e_dfa_" + id + ".bin", e_dfa);
  utils::save_integer_matrix("common_" + id + ".bin", common_triangles);
  // utils::save_integer_matrix("basal_" + id + ".bin", basal_triangles);
  // utils::save_matrix("e_dfb_" + id + ".bin", e_dfb);
  // ***********************************************************
}
