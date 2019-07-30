// library for loading data for the summary plot

#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <Eigen/Dense>

#define APICAL_CUTOFF 0.8
#define BASAL_CUTOFF 0.8

// cell to cell connectivity
// this_triangle, other_cell, other_triangle
enum cell_conn{ \
  tTri, oCell, oTri, \
  CCONNCOUNT};


struct Mesh {
  int vertices_count;
  Eigen::Array<double, Eigen::Dynamic, 3, Eigen::RowMajorBit> vertices; // 3x coordinate
//  int surface_triangles_count;
//  Eigen::Array<int, Eigen::Dynamic, 3, Eigen::RowMajorBit> surface_triangles; // 3x vertex
  int tetrahedrons_count;
  Eigen::Array<int, Eigen::Dynamic, 4, Eigen::RowMajorBit> tetrahedrons; // 4x vertex
//  int common_triangles_count;
//  Eigen::Array<int, Eigen::Dynamic, CCONNCOUNT, Eigen::RowMajorBit> common_triangles; // this triangle, other cell, other triangle
//  int apical_triangles_count;
//  Eigen::Array<int, Eigen::Dynamic, 1> apical_triangles; // surface triangle indicies
//  int basal_triangles_count;
//  Eigen::Array<int, Eigen::Dynamic, 1> basal_triangles; // surface triangle indicies
  Eigen::Array<double, Eigen::Dynamic, 1> dfa; // distance from apical (per element)
  Eigen::Array<double, Eigen::Dynamic, 1> dfb; // distance from basal (per element)
  std::set<int> apical_nodes;  // list of unique nodes belonging to apical tets
  std::set<int> basal_nodes;  // list of unique nodes belonging to basal tets
};



void fatal_error(std::string msg);
void load_mesh(std::string file_name, struct Mesh& mesh);


void fatal_error(std::string msg) {
  std::cerr << msg << std::endl;
  exit(1);
}

void load_mesh(std::string file_name, struct Mesh& mesh){
  std::cout << "Loading mesh file: " << file_name << std::endl;

  // local variables
  std::ifstream cell_file(file_name.c_str(), std::ios::in | std::ios::binary); // open the mesh file
  uint32_t i32;
  float f32;

  // check the file is open
  if (not cell_file.is_open()) {
    fatal_error("mesh file " + file_name + " could not be opened");
  }
  // get the mesh vertices (int32 count, 3x-float32 vertices) 
  cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
  mesh.vertices_count = i32;
  std::cout << "  Number of nodes = " << mesh.vertices_count << std::endl;
  mesh.vertices.resize(mesh.vertices_count, Eigen::NoChange);
  for(int n=0; n<mesh.vertices_count; n++){
    for(int m=0; m<3; m++){
      cell_file.read(reinterpret_cast<char *>(&f32), sizeof(f32));
      mesh.vertices(n,m) = f32;
    }
  }
  // get the surface triangles (int32 count, 3x-int32 vertex indices)
  cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
  int surface_triangles_count = i32;
//  mesh.surface_triangles.resize(mesh.surface_triangles_count, Eigen::NoChange);
  for(int n=0; n<surface_triangles_count; n++){
    for(int m=0; m<3; m++){
      cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
//      mesh.surface_triangles(n,m) = i32-1; // change to zero indexed
    }
  }
  // get the element tetrahedrons (int32 count, 4x-int32 vertex indices, float32 dfa, float32 dfb)
  cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
  mesh.tetrahedrons_count = i32;
  mesh.tetrahedrons.resize(mesh.tetrahedrons_count, Eigen::NoChange);
  mesh.dfa.resize(mesh.tetrahedrons_count, Eigen::NoChange);
  mesh.dfb.resize(mesh.tetrahedrons_count, Eigen::NoChange);
  for(int n=0; n<mesh.tetrahedrons_count; n++){
    for(int m=0; m<4; m++){
      cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
      mesh.tetrahedrons(n,m) = i32-1; // change to zero indexed
    }
    cell_file.read(reinterpret_cast<char *>(&f32), sizeof(f32));
    mesh.dfa(n) = f32;
    cell_file.read(reinterpret_cast<char *>(&f32), sizeof(f32));
    mesh.dfb(n) = f32;
  }
  // get the apical triangles (int32 count, int32 triangle indices)
  cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
  int apical_triangles_count = i32;
//  mesh.apical_triangles.resize(mesh.apical_triangles_count, Eigen::NoChange);
  for(int n=0; n<apical_triangles_count; n++){
    cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
//    mesh.apical_triangles(n) = i32-1; // change to zero indexed
  }
  // get the basal triangles (int32 count, int32 triangle indices)
  cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
  int basal_triangles_count = i32;
//  mesh.basal_triangles.resize(mesh.basal_triangles_count, Eigen::NoChange);
  for(int n=0; n<basal_triangles_count; n++){
    cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
//    mesh.basal_triangles(n) = i32-1; // change to zero indexed
  }
  // get the cell-to-cell data (int32 count, 3x-int32 this_triamgle, other_cell, other_triangle)
  cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
  int common_triangles_count = i32;
//  mesh.common_triangles.resize(mesh.common_triangles_count, Eigen::NoChange);
  for(int n=0; n<common_triangles_count; n++){
    for(int m=0; m<3; m++){
      cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
//     mesh.common_triangles(n,m) = i32-1; // change to zero indexed
    }
  }
  cell_file.close();

  // get the apical/basal nodes
  for (int n = 0; n < mesh.tetrahedrons_count; n++) {
    // volume of tet
    Eigen::Matrix<int,1,4> vi;      // tetrahedron vertex indices
    vi = mesh.tetrahedrons.block<1,4>(n, 0);

    // apical or basal?
    if (mesh.dfa(n) < APICAL_CUTOFF) {
      for (int i = 0; i < 4; i++) {
        mesh.apical_nodes.insert(vi(i));
      }
    }
    if (mesh.dfb(n) < BASAL_CUTOFF) {
      for (int i = 0; i < 4; i++) {
        mesh.basal_nodes.insert(vi(i));
      }
    }

  }
  std::cout << "  Number of apical nodes " << mesh.apical_nodes.size() << std::endl;
  std::cout << "  Number of basal nodes " << mesh.basal_nodes.size() << std::endl;
}

void get_plot_data(std::string file_name, int ntime, struct Mesh& mesh, float* result_apical, float* result_basal) {
  // each row is a time point, each column is a node
  
  std::ifstream result_file(file_name.c_str(), std::ios::in | std::ios::binary); // open the mesh file
  if (result_file) {
    int np = mesh.vertices_count;
    std::vector<float> buffer(np);
    // loop over time steps
    for (int t = 0; t < ntime; t++) {
      // load value at each node for this time point
      result_file.read(reinterpret_cast<char*>(buffer.data()), np * sizeof(float));

      if (result_file) {
        // just average values at nodes
        result_apical[t] = 0.0;
        for (auto index : mesh.apical_nodes) {
          result_apical[t] += buffer[index];
        }
        result_apical[t] /= static_cast<float>(mesh.apical_nodes.size());

        result_basal[t] = 0.0;
        for (auto index : mesh.basal_nodes) {
          result_basal[t] += buffer[index];
        }
        result_basal[t] /= static_cast<float>(mesh.basal_nodes.size());
      }
      else {
        fatal_error("Error reading row");
      }
    }
    result_file.close();
  }
  else {
    fatal_error(std::string("Error loading result file: ") + file_name);
  }
}

// function to be called from Python
extern "C"
int load_summary_plot_data(char* id, int ntime, float* ca_apical, float* ca_basal, float* ip_apical, float* ip_basal) {
  std::string mesh_file = std::string(id) + std::string(".bmsh");
  struct Mesh mesh;
  load_mesh(mesh_file, mesh);

  std::string ca_file = std::string(id) + std::string("_ca.bin");
  get_plot_data(ca_file, ntime, mesh, ca_apical, ca_basal);

  std::string ip_file = std::string(id) + std::string("_ip3.bin");
  get_plot_data(ip_file, ntime, mesh, ip_apical, ip_basal);

  return 0;
}
