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
#include "cLumenTree.hpp"
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

void cCellMesh::mesh_calcs(){
  parent->out << "<CellMesh> calculating derived properties... " << std::endl;

  // calculate the distance from nodes to the lumen
  cLumenTree lumen(parent);
  n_dfa.resize(vertices_count, Eigen::NoChange);
  for(int n=0; n<vertices_count; n++){
    n_dfa(n) = lumen.get_dnl(vertices.row(n));
  }
  // calculate the distance from elements to the lumen  
  e_dfa.resize(tetrahedrons_count, Eigen::NoChange);
  for(int n=0; n<tetrahedrons_count; n++){
    tCalcs d = 0.0; // distance from element to apical
    for(int m=0; m<4; m++){
	  d += n_dfa(tetrahedrons(n,m));  // sum distance from node to apical
    }
	e_dfa(n) = d / 4.0; // the average over four vertices
  }
  // determine the apical surface triangle indices
  apical_triangles.resize(surface_triangles_count, Eigen::NoChange); // overkill size for now
  apical_triangles_count = 0;
  for(int n=0; n<surface_triangles_count; n++){
	double d = (                         // triangle distance from apical
		n_dfa(surface_triangles(n,0)) + 
		n_dfa(surface_triangles(n,1)) +
		n_dfa(surface_triangles(n,2))) / 3.0;
    if (d < parent->p[IPRdn]) {
	  apical_triangles(apical_triangles_count++) = n;
    }
  }
  apical_triangles.conservativeResize(apical_triangles_count, 1); // actual triangles count

  // determine the cell-to-cell connetivity data (this_triamgle, other_cell, other_triangle)
  //   do this by comparing surface triangle centers on this cell to those on the surface of every other cell
  Eigen::Array<tCoord, Eigen::Dynamic, 3, Eigen::RowMajorBit> this_centers; // calculate this cell surface triangle centers
  utils::calc_tri_centers(&this_centers, &vertices, &surface_triangles);
  int ci = id.find_last_of('c') + 1; // to split the mesh id string, find the position of the last "c"
  int this_cell = atoi(id.substr(ci).c_str()) - 1; // extract this cell index from its mesh id
  common_triangles.resize(surface_triangles_count, Eigen::NoChange); // overkill size for now
  common_triangles_count = 0;
  for(int other_cell=0; other_cell<CELLS_COUNT; other_cell++){ // check against all the other cell meshes...
  	if(other_cell == this_cell) continue; // don't check against itself
    std::string other_id = id.substr(0, ci) + std::to_string(other_cell + 1);  // the other mesh id
    cCellMesh* other_mesh = new cCellMesh(other_id, parent); // temp other_mesh
    Eigen::Array<tCoord, Eigen::Dynamic, 3, Eigen::RowMajorBit> other_centers; // calculate other cell surface triangle centers
    utils::calc_tri_centers(&other_centers, &(other_mesh->vertices), &(other_mesh->surface_triangles));
    for(int this_tri = 0; this_tri < surface_triangles_count; this_tri++){
      for(int other_tri = 0; other_tri < other_mesh->surface_triangles_count; other_tri++){
        Eigen::Vector3d p1 = this_centers.row(this_tri);
        Eigen::Vector3d p2 = other_centers.row(other_tri);
        if(p1==p2){
          common_triangles.row(common_triangles_count++) = Eigen::Vector3i(this_tri, other_cell, other_tri);
        }
      }
    }
  	delete other_mesh; // release temp other_mesh
  }
  common_triangles.conservativeResize(common_triangles_count, Eigen::NoChange); // actual size

  // determine the basal triangle indices by considering all surface triangles
  //   then eliminating the common triangles and the triangles that are too close to the lumen
  //
  // eliminate the common triangles
  // ...
  // calculate the apical keep-out traingles and eliminate them
  // ...

  // calculate the distance from elements to basal surface
  //   use the average of the element-vertex to nearest-basal-triangle-vertex distances  
  // 
  // list the basal vertices
  // ...
  // find the distance from every node to the nearest basal vertex (n_dnb)
  // ...
  // for each tet calculate e_dnb as the average it's vertice dnb's
  // ...
}

void cCellMesh::get_mesh(std::string file_name){
  parent->out << "<CellMesh> reading mesh file " + file_name << std::endl;
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
  e_dfb.resize(tetrahedrons_count, Eigen::NoChange);
  for(int n=0; n<tetrahedrons_count; n++){
    for(int m=0; m<4; m++){
      cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
      tetrahedrons(n,m) = i32-1;      // change to zero indexed
    }
   cell_file.read(reinterpret_cast<char *>(&f32), sizeof(f32));  // e_dfa
   cell_file.read(reinterpret_cast<char *>(&f32), sizeof(f32));  // e_dfb
   e_dfb(n) = f32;
  }
  // get the apical triangles (int32 count, int32 triangle indices)
  cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
  apical_triangles_count = i32;
  for(int n=0; n<apical_triangles_count; n++){
	cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
  }
  // get the basal triangles (int32 count, int32 triangle indices)
  cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
  basal_triangles_count = i32;
  basal_triangles.resize(basal_triangles_count, Eigen::NoChange);
  for(int n=0; n<basal_triangles_count; n++){
    cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
    basal_triangles(n) = i32-1; // change to zero indexed
  }
  // get the cell-to-cell data (int32 count, int32 this_triamgle, int32 other_cell, int32 other_triangle)
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

void cCellMesh::print_info(){
  parent->out << "<CellMesh> vertices_count: " << vertices_count << std::endl;
  parent->out << "<CellMesh> tetrahedrons_count: " << tetrahedrons_count << std::endl;
  parent->out << "<CellMesh> surface_triangles_count: " << surface_triangles_count << std::endl;
  parent->out << "<CellMesh> apical_triangles_count: " << apical_triangles_count << std::endl;
  parent->out << "<CellMesh> basal_triangles_count: " << basal_triangles_count << std::endl;
  parent->out << "<CellMesh> common_triangles_count: " << common_triangles_count << std::endl;

  // ***********************************************************
  //utils::save_matrix("vertices_" + id + ".bin", vertices);
  //utils::save_integer_matrix("triangles_" + id + ".bin", surface_triangles);
  //utils::save_integer_matrix("tetrahedrons_" + id + ".bin", tetrahedrons);
  utils::save_integer_matrix("apical_" + id + ".bin", apical_triangles);
  //utils::save_matrix("n_dfa_" + id + ".bin", n_dfa);
  //utils::save_matrix("e_dfa_" + id + ".bin", e_dfa);
  utils::save_integer_matrix("basal_" + id + ".bin", basal_triangles);
  //utils::save_matrix("e_dfb_" + id + ".bin", e_dfb);
  // ***********************************************************
}
