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


void fatal_error(std::string msg) {
  std::cerr << msg << std::endl;
  exit(1);
}

struct Mesh {
  int vertices_count;
  Eigen::Array<double, Eigen::Dynamic, 3, Eigen::RowMajorBit> vertices; // 3x coordinate
  int tetrahedrons_count;
  Eigen::Array<int, Eigen::Dynamic, 4, Eigen::RowMajorBit> tetrahedrons; // 4x vertex
  Eigen::Array<double, Eigen::Dynamic, 1> dfa; // distance from apical (per element)
  Eigen::Array<double, Eigen::Dynamic, 1> dfb; // distance from basal (per element)
  std::set<int> apical_nodes;  // list of unique nodes belonging to apical tets
  std::set<int> basal_nodes;  // list of unique nodes belonging to basal tets

  void load(std::string file_name){
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
    vertices_count = i32;
    std::cout << "  Number of nodes = " << vertices_count << std::endl;
    vertices.resize(vertices_count, Eigen::NoChange);
    for(int n=0; n<vertices_count; n++){
      for(int m=0; m<3; m++){
        cell_file.read(reinterpret_cast<char *>(&f32), sizeof(f32));
        vertices(n,m) = f32;
      }
    }
    // get the surface triangles (int32 count, 3x-int32 vertex indices)
    cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
    int surface_triangles_count = i32;
    for(int n=0; n<surface_triangles_count; n++){
      for(int m=0; m<3; m++){
        cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
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
    int apical_triangles_count = i32;
    for(int n=0; n<apical_triangles_count; n++){
      cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
    }
    // get the basal triangles (int32 count, int32 triangle indices)
    cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
    int basal_triangles_count = i32;
    for(int n=0; n<basal_triangles_count; n++){
      cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
    }
    // get the cell-to-cell data (int32 count, 3x-int32 this_triamgle, other_cell, other_triangle)
    cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
    int common_triangles_count = i32;
    for(int n=0; n<common_triangles_count; n++){
      for(int m=0; m<3; m++){
        cell_file.read(reinterpret_cast<char *>(&i32), sizeof(i32));
      }
    }
    cell_file.close();

    // get the apical/basal nodes
    for (int n = 0; n < tetrahedrons_count; n++) {
      // volume of tet
      Eigen::Matrix<int,1,4> vi;      // tetrahedron vertex indices
      vi = tetrahedrons.block<1,4>(n, 0);

      // apical or basal?
      if (dfa(n) < APICAL_CUTOFF) {
        for (int i = 0; i < 4; i++) {
          apical_nodes.insert(vi(i));
        }
      }
      if (dfb(n) < BASAL_CUTOFF) {
        for (int i = 0; i < 4; i++) {
          basal_nodes.insert(vi(i));
        }
      }

    }
    std::cout << "  Number of apical nodes " << apical_nodes.size() << std::endl;
    std::cout << "  Number of basal nodes " << basal_nodes.size() << std::endl;
  }
};

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
        fatal_error(std::string("Error reading row: ") + std::to_string(t + 1));
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
  mesh.load(mesh_file);

  std::string ca_file = std::string(id) + std::string("_ca.bin");
  get_plot_data(ca_file, ntime, mesh, ca_apical, ca_basal);

  std::string ip_file = std::string(id) + std::string("_ip3.bin");
  get_plot_data(ip_file, ntime, mesh, ip_apical, ip_basal);

  return 0;
}
