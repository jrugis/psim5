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
#include "cLumenBase.hpp"
#include "cCell_calcium.hpp"
#include "cCellMesh.hpp"

cCellMesh::cCellMesh(std::string mesh_name, cCell_calcium* p){
  // initialise member variables
  vertices_count = tetrahedrons_count = 0;
  surface_triangles_count = apical_triangles_count = basal_triangles_count = 0;
  common_triangles_count = 0;
  parent = p;
  id = mesh_name;  
  get_mesh(id + ".bmsh");
}

cCellMesh::~cCellMesh(){
}

void cCellMesh::get_mesh(std::string file_name){
  // local variables
  std::ifstream cell_file(file_name.c_str(), std::ios::in | std::ios::binary); // open the mesh file
  uint32_t i32;
  float f32;
  cLumenBase* lumen;  // needed for apical/basal mesh geometry calculations
  lumen = new cLumenBase(parent);

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
    // calculate the distance from node to the lumen
    n_dfa.resize(vertices_count, Eigen::NoChange);
	//n_dfa(n) = lumen->get_dnl(vertices.block<1,3>(n, 0));
	n_dfa(n) = lumen->get_dnl(vertices.row(n));
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
  e_dfa.resize(tetrahedrons_count, Eigen::NoChange);
  e_dfb.resize(tetrahedrons_count, Eigen::NoChange);
  for(int n=0; n<tetrahedrons_count; n++){
    tCalcs d = 0.0; // distance from element to apical
    for(int m=0; m<4; m++){
      cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
      tetrahedrons(n,m) = i32-1;      // change to zero indexed
	  d += n_dfa(tetrahedrons(n,m));  // sum distance from node to apical
    }
    cell_file.read(reinterpret_cast<char *>(&f32), sizeof(f32));
    //e_dfa(n) = f32;
	e_dfa(n) = d / 4; // the average over four vertices
    cell_file.read(reinterpret_cast<char *>(&f32), sizeof(f32));
    e_dfb(n) = f32;
  }
  // get the apical triangles (int32 count, int32 triangle indices)
  cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
  apical_triangles_count = i32;
  //apical_triangles.resize(apical_triangles_count, Eigen::NoChange);
  for(int n=0; n<apical_triangles_count; n++){
	cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
    //apical_triangles(n) = i32-1; // change to zero indexed
  }
  
  apical_triangles.resize(surface_triangles_count, Eigen::NoChange); // overkill size for now
  apical_triangles_count = 0;
  for(int n=0; n<surface_triangles_count; n++){
	double d = (    // triangle distance from apical
		n_dfa(surface_triangles(n,0)) + 
		n_dfa(surface_triangles(n,1)) +
		n_dfa(surface_triangles(n,2))) / 3.0;
    if (d < 0.8) {
	  apical_triangles(apical_triangles_count++) = n;
    }
  }
  apical_triangles.conservativeResize(apical_triangles_count, 1); // actual triangles count
    
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
  delete lumen;
  // **************** DEBUG ************************************
  //utils::save_matrix("vertices_" + id + ".bin", vertices);
  //utils::save_integer_matrix("triangles_" + id + ".bin", surface_triangles);
  //utils::save_integer_matrix("apical_" + id + ".bin", apical_triangles);
  //utils::save_integer_matrix("tetrahedrons_" + id + ".bin", tetrahedrons);
  //utils::save_matrix("n_dfa_" + id + ".bin", n_dfa);
  //utils::save_matrix("e_dfa_" + id + ".bin", e_dfa);
  // ***********************************************************
}

void cCellMesh::print_info(){
  parent->out << "<CellMesh> vertices_count: " << vertices_count << std::endl;
  parent->out << "<CellMesh> tetrahedrons_count: " << tetrahedrons_count << std::endl;
  parent->out << "<CellMesh> surface_triangles_count: " << surface_triangles_count << std::endl;
  parent->out << "<CellMesh> apical_triangles_count: " << apical_triangles_count << std::endl;
  parent->out << "<CellMesh> basal_triangles_count: " << basal_triangles_count << std::endl;
  parent->out << "<CellMesh> common_triangles_count: " << common_triangles_count << std::endl;
}
