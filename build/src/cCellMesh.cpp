/*
 * cCellMesh.cpp
 *
 *  Created on: 26/04/2018
 *      Author: jrugis
 */

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <cmath>

#include "utils.hpp"
#include "cCell_calcium.hpp"
#include "cCellMesh.hpp"

cCellMesh::cCellMesh(std::string mesh_name, cCell_calcium* p){
  // initialise member variables
  vertices_count = tetrahedrons_count = 0;
  surface_triangles_count = apical_triangles_count = basal_triangles_count = 0;
  common_triangles_count = 0;
  apical_triangles_area = basal_triangles_area = 0.0;
  parent = p;
  id = mesh_name;  
  get_mesh(id + ".bmsh");
  compute_surface_triangle_areas();
  identify_common_apical_triangles();
}

cCellMesh::~cCellMesh(){
}

void cCellMesh::get_mesh(std::string file_name){
  // local variables
  std::ifstream cell_file(file_name.c_str(), std::ios::in | std::ios::binary); // open the mesh file
  uint32_t i32;
  float f32;

  // check the file is open
  if (not cell_file.is_open()) {
    utils::fatal_error("mesh file " + file_name + " could not be opened", parent->out);
  }
  // get the mesh vertices (int32 count, 3x-float32 vertices) 
  cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
  vertices_count = i32;
  vertices.resize(vertices_count, Eigen::NoChange);
  for(int n=0; n<vertices_count; n++){
    for(int m=0; m<3; m++){
      cell_file.read(reinterpret_cast<char *>(&f32), sizeof(f32));
      vertices(n,m) = f32;
    }
  }
  // get the surface triangles (int32 count, 3x-int32 vertex indices)
  cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
  surface_triangles_count = i32;
  surface_triangles.resize(surface_triangles_count, Eigen::NoChange);
  for(int n=0; n<surface_triangles_count; n++){
    for(int m=0; m<3; m++){
      cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
      surface_triangles(n,m) = i32-1; // change to zero indexed
    }
  }
  // get the element tetrahedrons (int32 count, 4x-int32 vertex indices, float32 dfa, float32 dfb)
  cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
  tetrahedrons_count = i32;
  tetrahedrons.resize(tetrahedrons_count, Eigen::NoChange);
  dfa.resize(tetrahedrons_count, Eigen::NoChange);
  dfb.resize(tetrahedrons_count, Eigen::NoChange);
  for(int n=0; n<tetrahedrons_count; n++){
    for(int m=0; m<4; m++){
      cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
      tetrahedrons(n,m) = i32-1; // change to zero indexed
    }
    cell_file.read(reinterpret_cast<char *>(&f32), sizeof(f32));
    dfa(n) = f32;
    cell_file.read(reinterpret_cast<char *>(&f32), sizeof(f32));
    dfb(n) = f32;
  }
  // get the apical triangles (int32 count, int32 triangle indices)
  cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
  apical_triangles_count = i32;
  apical_triangles.resize(apical_triangles_count, Eigen::NoChange);
  for(int n=0; n<apical_triangles_count; n++){
    cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
    apical_triangles(n) = i32-1; // change to zero indexed
  }
  // get the basal triangles (int32 count, int32 triangle indices)
  cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
  basal_triangles_count = i32;
  basal_triangles.resize(basal_triangles_count, Eigen::NoChange);
  for(int n=0; n<basal_triangles_count; n++){
    cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
    basal_triangles(n) = i32-1; // change to zero indexed
  }
  // get the cell-to-cell data (int32 count, 3x-int32 this_triamgle, other_cell, other_triangle)
  cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
  common_triangles_count = i32;
  common_triangles.resize(common_triangles_count, Eigen::NoChange);
  for(int n=0; n<common_triangles_count; n++){
    for(int m=0; m<3; m++){
      cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
      common_triangles(n,m) = i32-1; // change to zero indexed
    }
  }
  cell_file.close();
}

bool cCellMesh::is_triangle_apical(int triangle, int& apical_index) {
  // return True if the given triangle is apical, otherwise False
  // this assumes the list of apical triangle indices is sorted from lowest to highest
  // this has been verified for the current meshes

  // binary search for triangle index in apical triangles array
  int l = 0;
  int r = apical_triangles_count - 1;
  while (l <= r) {
    // cut region in half
    int m = (l + r) / 2;

    // check if index present at midpoint
    if (apical_triangles(m) == triangle) {
      apical_index = m;
      return true;
    }

    // if index greater, search in right half
    if (apical_triangles(m) < triangle) {
      l = m + 1;
    }

    // if index smaller, search in left half
    else {
      r = m - 1;
    }
  }

//  // alternative to binary search
//  for (int i = 0; i < apical_triangles_count; i++) {
//    if (triangle == apical_triangles(i)) {
//      // triangle is apical
//      apical_index = i;
//      return true;
//    }
//  }

  // triangle is not apical
  apical_index = -1;
  return false;
}

void cCellMesh::identify_common_apical_triangles() {
  parent->out << "<CellMesh> building common apical triangles list" << std::endl;

  // NOTE: this could be a feature of the mesh (i.e. included in the mesh file)

//    parent->out << "CHECKING APICAL SORTED...\n";
//    for (int i = 1; i < apical_triangles_count; i++) {
//        if (apical_triangles(i) < apical_triangles(i - 1)) {
//            parent->out << "UNSORTED APICAL TRIANGLES\n";
//        }
//    }

  common_apical_triangles.resize(apical_triangles_count, Eigen::NoChange);

  // record which apical triangles have been included from the common triangles list,
  // so can add any remaining apical triangles to the list
  std::vector<int> apical_mask(apical_triangles_count, 0);

  // start with common triangles that are also apical triangles
  int count = 0;
  for (int i = 0; i < common_triangles_count; i++) {
    int apical_index = -1;
    if (is_triangle_apical(common_triangles(i, tTri), apical_index)) {
      apical_mask[apical_index] = 1;
      for (int j = 0; j < CCONNCOUNT; j++) {
        common_apical_triangles(count, j) = common_triangles(i, j);
      }
      count++;
    }
  }
  parent->out << "<CellMesh> apical triangles common with other cells: " << count << std::endl;

  // now add on the apical triangles not in the common triangles list
  for (int i = 0; i < apical_triangles_count; i++) {
    if (not apical_mask[i]) {
      common_apical_triangles(count, tTri) = apical_triangles(i);
      common_apical_triangles(count, oCell) = parent->cell_number;
      common_apical_triangles(count, oTri) = apical_triangles(i);
      count++;
    }
  }
}

void cCellMesh::compute_surface_triangle_areas() {
  // NOTE: this could be a feature of the mesh (i.e. included in the mesh file)

  // store areas of all surface triangles
  surface_triangle_areas.resize(surface_triangles_count, Eigen::NoChange);
  for (int i = 0; i < surface_triangles_count; i++) {
    // vertices of this triangle
    Eigen::Vector3d a = vertices.row(surface_triangles(i, 0));
    Eigen::Vector3d b = vertices.row(surface_triangles(i, 1));
    Eigen::Vector3d c = vertices.row(surface_triangles(i, 2));

    // two sides of the triangle from a
    Eigen::Vector3d ab = b - a;
    Eigen::Vector3d ac = c - a;

    // compute the area
    surface_triangle_areas(i) = ab.cross(ac).norm() / 2.0;
  }

  // total area of apical triangles
  apical_triangles_area = 0;
  for (int i = 0; i < apical_triangles_count; i++) {
    int index = apical_triangles(i);
    apical_triangles_area += surface_triangle_areas(index);
  }

  // total area of basal triangles
  basal_triangles_area = 0;
  for (int i = 0; i < basal_triangles_count; i++) {
    int index = basal_triangles(i);
    basal_triangles_area += surface_triangle_areas(index);
  }
}

void cCellMesh::print_info(){
  parent->out << "<CellMesh> vertices_count: " << vertices_count << std::endl;
  parent->out << "<CellMesh> tetrahedrons_count: " << tetrahedrons_count << std::endl;
  parent->out << "<CellMesh> surface_triangles_count: " << surface_triangles_count << std::endl;
  parent->out << "<CellMesh> apical_triangles_count: " << apical_triangles_count << std::endl;
  parent->out << "<CellMesh> apical_triangles_area: " << apical_triangles_area << std::endl;
  parent->out << "<CellMesh> basal_triangles_count: " << basal_triangles_count << std::endl;
  parent->out << "<CellMesh> basal_triangles_area: " << basal_triangles_area << std::endl;
  parent->out << "<CellMesh> common_triangles_count: " << common_triangles_count << std::endl;
}
