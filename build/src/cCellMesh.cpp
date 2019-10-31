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
  utils::read_mesh(id, mesh_vals, parent->out);
  mesh_calcs();
}

cCellMesh::~cCellMesh() {}

void cCellMesh::calc_dfa()
{
  // calculate the distance from nodes to the lumen
  cLumenTree lumen(parent->out);
  n_dfa.resize(mesh_vals.vertices_count, Eigen::NoChange);
  for (int n = 0; n < mesh_vals.vertices_count; n++) { n_dfa(n) = lumen.get_dnl(mesh_vals.vertices.row(n)); }
  // calculate the distance from elements to the lumen
  e_dfa.resize(mesh_vals.tetrahedrons_count, Eigen::NoChange);
  for (int n = 0; n < mesh_vals.tetrahedrons_count; n++) {
    e_dfa(n) = (n_dfa(mesh_vals.tetrahedrons(n, 0)) + n_dfa(mesh_vals.tetrahedrons(n, 1)) +
                n_dfa(mesh_vals.tetrahedrons(n, 2)) + n_dfa(mesh_vals.tetrahedrons(n, 3))) /
               4.0;
  }
}

// determine the cell-to-cell connetivity data (this_triamgle, other_cell, other_triangle)
//	 do this by comparing surface triangle centers on this cell to those on every other cell
void cCellMesh::calc_common()
{
  MatrixN3d this_centers;  // this cell surface triangle centers
  MatrixN3d other_centers; // other cell surface triangle centers
  sMeshVals other_mesh;    // other mesh values

  utils::calc_tri_centers(this_centers, mesh_vals.vertices, mesh_vals.surface_triangles);
  int ci = id.find_last_of('c') + 1;               // to split the mesh id string, find the position of the last "c"
  int this_cell = atoi(id.substr(ci).c_str()) - 1; // extract this cell index from its mesh id
  common_triangles.resize(mesh_vals.surface_triangles_count, Eigen::NoChange); // overkill size for now
  common_triangles_count = 0;
  for (int other_cell = 0; other_cell < CELLS_COUNT; other_cell++) { // check against all the other cell meshes...
    if (other_cell == this_cell) continue;                           // don't check against itself
    std::string other_id = id.substr(0, ci) + std::to_string(other_cell + 1); // the other mesh id
    utils::read_mesh(other_id, other_mesh, parent->out);
    utils::calc_tri_centers(other_centers, other_mesh.vertices, other_mesh.surface_triangles);
    for (int this_tri = 0; this_tri < mesh_vals.surface_triangles_count; this_tri++) {
      for (int other_tri = 0; other_tri < other_mesh.surface_triangles_count; other_tri++) {
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
  MatrixN1i apical_keepout;                                                  // keep-out triangle indicies
  apical_keepout.resize(mesh_vals.surface_triangles_count, Eigen::NoChange); // overkill size for now
  int apical_keepout_count = 0;
  apical_triangles.resize(mesh_vals.surface_triangles_count, Eigen::NoChange); // overkill size for now
  apical_triangles_count = 0;
  for (int n = 0; n < mesh_vals.surface_triangles_count; n++) {
    double d = // triangle distance from apical
      (n_dfa(mesh_vals.surface_triangles(n, 0)) + n_dfa(mesh_vals.surface_triangles(n, 1)) +
       n_dfa(mesh_vals.surface_triangles(n, 2))) /
      3.0;
    if (d < parent->p[IPRdn]) apical_triangles(apical_triangles_count++) = n;
    if (d < (parent->p[PLCdl])) apical_keepout(apical_keepout_count++) = n;
  }
  apical_triangles.conservativeResize(apical_triangles_count, 1); // actual triangles count
  apical_keepout.conservativeResize(apical_keepout_count, 1);     // actual triangles count

  // determine the basal triangle indices by considering all surface triangles
  //	 then eliminating the common triangles and the triangles that are too close to the lumen
  basal_triangles.resize(mesh_vals.surface_triangles_count, Eigen::NoChange); // overkill size for now
  basal_triangles.setOnes(mesh_vals.surface_triangles_count);
  for (int n = 0; n < common_triangles_count; n++) { // eliminate the common triangles
    basal_triangles(common_triangles(n, 0)) = 0;
  }
  for (int n = 0; n < apical_keepout_count; n++) { // eliminate the apical keepout triangles
    basal_triangles(apical_keepout(n)) = 0;
  }
  basal_triangles_count = 0;
  for (int n = 0; n < mesh_vals.surface_triangles_count; n++) { // convert what's left to a list of indices
    if (basal_triangles(n)) { basal_triangles(basal_triangles_count++) = n; }
  }
  basal_triangles.conservativeResize(basal_triangles_count, Eigen::NoChange); // actual size
}

// calculate the distance from elements to basal surface
//	 use the average of the element-vertex to nearest-basal-triangle-vertex distances
void cCellMesh::calc_dfb()
{
  // get the basal vertices
  MatrixN1i basal_verts;                                         // the basal triangle vertices
  basal_verts.resize(mesh_vals.vertices_count, Eigen::NoChange); // overkill size for now
  basal_verts.setZero(mesh_vals.vertices_count);
  for (int n = 0; n < basal_triangles_count; n++) {
    Eigen::Vector3i vi = Eigen::Vector3i(mesh_vals.surface_triangles.row(basal_triangles(n)));
    for (int i = 0; i < 3; i++) basal_verts(vi(i)) = 1; // flag vertex as basal
  }
  int basal_verts_count = 0;
  for (int n = 0; n < mesh_vals.vertices_count; n++) { // convert flags to a list of basal indices
    if (basal_verts(n)) basal_verts(basal_verts_count++) = n;
  }
  basal_verts.conservativeResize(basal_verts_count, Eigen::NoChange); // actual size
  // calculate the (per node) node to nearest basal node distance
  MatrixN1d n_dfb;
  n_dfb.resize(mesh_vals.vertices_count, Eigen::NoChange);
  for (int n = 0; n < mesh_vals.vertices_count; n++) {
    Eigen::Vector3d p1 = Eigen::Vector3d(mesh_vals.vertices.row(n));
    n_dfb(n) = 1000.0; // large dummy value
    for (int m = 0; m < basal_verts_count; m++) {
      Eigen::Vector3d p2 = Eigen::Vector3d(mesh_vals.vertices.row(basal_verts(m)));
      double d = (p1 - p2).norm();
      if (d < n_dfb(n)) n_dfb(n) = d;
    }
  }
  // for each tet calculate e_dnb as the average it's vertex dnb's
  e_dfb.resize(mesh_vals.tetrahedrons_count, Eigen::NoChange);
  for (int n = 0; n < mesh_vals.tetrahedrons_count; n++) {
    e_dfb(n) = (n_dfb(mesh_vals.tetrahedrons(n, 0)) + n_dfb(mesh_vals.tetrahedrons(n, 1)) +
                n_dfb(mesh_vals.tetrahedrons(n, 2)) + n_dfb(mesh_vals.tetrahedrons(n, 3))) /
               4.0;
  }
}

void cCellMesh::mesh_calcs()
{
  parent->out << "<CellMesh> calculating derived properties... " << std::endl;
  calc_common();
  calc_dfa();          // do this first,				distance from apical (per node and per element)
  calc_apical_basal(); // then this,						identify apical and basal triangles
  calc_dfb();          // finally this.				  distance from basal (per element)
}

void cCellMesh::print_info()
{
  parent->out << "<CellMesh> vertices_count: " << mesh_vals.vertices_count << std::endl;
  parent->out << "<CellMesh> tetrahedrons_count: " << mesh_vals.tetrahedrons_count << std::endl;
  parent->out << "<CellMesh> surface_triangles_count: " << mesh_vals.surface_triangles_count << std::endl;
  parent->out << "<CellMesh> apical_triangles_count: " << apical_triangles_count << std::endl;
  parent->out << "<CellMesh> basal_triangles_count: " << basal_triangles_count << std::endl;
  parent->out << "<CellMesh> common_triangles_count: " << common_triangles_count << std::endl;

  // ***********************************************************
  //utils::save_matrix("vertices_" + id + ".bin", 3 * mesh_vals.vertices_count * sizeof(double),
  //                   reinterpret_cast<char*>(mesh_vals.vertices.data()));
  //utils::save_matrix("triangles_" + id + ".bin", 3 * mesh_vals.surface_triangles_count * sizeof(int),
  //                   reinterpret_cast<char*>(mesh_vals.surface_triangles.data()));
  //utils::save_matrix("tetrahedrons_" + id + ".bin", 4 * mesh_vals.tetrahedrons_count * sizeof(int),
  //                   reinterpret_cast<char*>(mesh_vals.tetrahedrons.data()));
  utils::save_matrix("apical_" + id + ".bin", apical_triangles_count * sizeof(int),
                     reinterpret_cast<char*>(apical_triangles.data()));
  // utils::save_matrix("n_dfa_" + id + ".bin", mesh_vals.vertices_count * sizeof(double), reinterpret_cast<char*>(n_dfa.data()));
  // utils::save_matrix("e_dfa_" + id + ".bin", mesh_vals.tetrahedrons_count * sizeof(double),
  //                   reinterpret_cast<char*>(e_dfa.data()));
  // utils::save_matrix("common_" + id + ".bin", 3 * common_triangles_count * sizeof(int),
  //                   reinterpret_cast<char*>(common_triangles.data()));
  utils::save_matrix("basal_" + id + ".bin", basal_triangles_count * sizeof(int),
                     reinterpret_cast<char*>(basal_triangles.data()));
  // utils::save_matrix("e_dfb_" + id + ".bin", mesh_vals.tetrahedrons_count * sizeof(double),
  //                   reinterpret_cast<char*>(e_dfb.data()));
  // ***********************************************************
}
