// library for loading data for the summary plot

#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <Eigen/Dense>
#include "global_defs.hpp"
#include "cCellMesh.hpp"
#include "global_defs.hpp"
#include "utils.hpp"

#define APICAL_CUTOFF 0.8
#define BASAL_CUTOFF 0.8


void my_fatal_error(std::string msg) {
  std::cerr << msg << std::endl;
  exit(1);
}

void get_plot_data(std::string file_name, int ntime, int np, std::set<int>& apical_nodes, std::set<int>& basal_nodes, float* result_apical, float* result_basal) {
  // each row is a time point, each column is a node
  std::cout << "Loading result file: " << file_name << std::endl;
  std::ifstream result_file(file_name.c_str(), std::ios::in | std::ios::binary); // open the mesh file
  if (result_file) {
    std::vector<float> buffer(np);
    // loop over time steps
    for (int t = 0; t < ntime; t++) {
      // load value at each node for this time point
      result_file.read(reinterpret_cast<char*>(buffer.data()), np * sizeof(float));

      if (result_file) {
        // just average values at nodes
        result_apical[t] = 0.0;
        for (auto index : apical_nodes) {
          result_apical[t] += buffer[index];
        }
        result_apical[t] /= static_cast<float>(apical_nodes.size());

        result_basal[t] = 0.0;
        for (auto index : basal_nodes) {
          result_basal[t] += buffer[index];
        }
        result_basal[t] /= static_cast<float>(basal_nodes.size());
      }
      else {
        my_fatal_error(std::string("Error reading row: ") + std::to_string(t + 1));
      }
    }
    result_file.close();
  }
  else {
    my_fatal_error(std::string("Error loading result file: ") + file_name);
  }
}

// function to be called from Python
extern "C"
int load_summary_plot_data(char* id, int cell_number, int ntime, double apical_cutoff,
    float* ca_apical, float* ca_basal, double basal_cutoff, float* ip_apical, float* ip_basal) {
  // load the mesh
  std::ofstream out("summary_plot.out");
  cCellMesh mesh = cCellMesh(id, out, cell_number);

  // get the apical and basal nodes
  std::set<int> apical_nodes;
  std::set<int> basal_nodes;
  for (int n = 0; n < mesh.tetrahedrons_count; n++) {
    // volume of tet
    Eigen::Matrix<int,1,4> vi;      // tetrahedron vertex indices
    vi = mesh.tetrahedrons.block<1,4>(n, 0);

    // apical or basal?
    if (mesh.dfa(n) < apical_cutoff) {
      for (int i = 0; i < 4; i++) {
        apical_nodes.insert(vi(i));
      }
    }
    if (mesh.dfb(n) < basal_cutoff) {
      for (int i = 0; i < 4; i++) {
        basal_nodes.insert(vi(i));
      }
    }
  }
  std::cout << "  Number of apical nodes " << apical_nodes.size() << std::endl;
  std::cout << "  Number of basal nodes " << basal_nodes.size() << std::endl;

  // ca results
  std::string ca_file = std::string(id) + std::string("_ca.bin");
  get_plot_data(ca_file, ntime, mesh.vertices_count, apical_nodes, basal_nodes, ca_apical, ca_basal);

  // ip3 results
  std::string ip_file = std::string(id) + std::string("_ip3.bin");
  get_plot_data(ip_file, ntime, mesh.vertices_count, apical_nodes, basal_nodes, ip_apical, ip_basal);

  return 0;
}
