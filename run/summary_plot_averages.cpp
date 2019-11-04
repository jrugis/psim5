// library for loading data for the summary plot

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>


// function to be called from Python
extern "C"
int load_summary_plot_data(char* result_filename, int ntime, int nnodes, float* avg_apical, float* avg_basal,
    int napical, int* apical, int nbasal, int* basal) {
  // each row is a time point, each column is a node
  std::ifstream result_file(result_filename, std::ios::in | std::ios::binary);
  if (result_file) {
    std::vector<float> buffer(nnodes);
    // loop over time steps
    for (int t = 0; t < ntime; t++) {
      // load value at each node for this time point
      result_file.read(reinterpret_cast<char*>(buffer.data()), nnodes * sizeof(float));

      if (result_file) {
        // just average values at nodes
        avg_apical[t] = 0.0;
        for (int i = 0; i < napical; i++) {
          int index = apical[i];
          avg_apical[t] += buffer[index];
        }
        avg_apical[t] /= static_cast<float>(napical);

        avg_basal[t] = 0.0;
        for (int i = 0; i < nbasal; i++) {
          int index = basal[i];
          avg_basal[t] += buffer[index];
        }
        avg_basal[t] /= static_cast<float>(nbasal);
      }
      else {
        return t + 1; // error reading row
      }
    }
    result_file.close();
  }
  else {
    return -1; // error opening file
  }

  return 0;
}
