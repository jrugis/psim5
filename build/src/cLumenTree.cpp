/*
 * cLumenTree.cpp
 *
 *	Created on: 08/10/19
 *	Author: jrugis
 */

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/tokenizer.hpp>
#include <string>

#include "cCell_calcium.hpp"
#include "cLumenTree.hpp"
#include "utils.hpp"

cLumenTree::cLumenTree(cCell_calcium* p) : parent(p)
{
  id = "l1";
  get_segments();
}

cLumenTree::~cLumenTree() {}

void cLumenTree::get_segments()
{
  std::string file_name = id + ".txt";
  std::ifstream lumen_file(file_name.c_str(), std::ios::in | std::ios::binary); // open the lumen file
  std::string line;                                                             // file line buffer
  std::vector<std::string> tokens;                                              // tokenized line

  // check the file is open
  if (not lumen_file.is_open()) { utils::fatal_error("lumen file " + file_name + " could not be opened", parent->out); }
  // get the lumen tree points
  getline(lumen_file, line);
  boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
  points_count = std::stoi(tokens[0]);
  points.resize(points_count, Eigen::NoChange);
  for (int n = 0; n < points_count; n++) {
    getline(lumen_file, line);
    boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
    for (int i = 0; i < 3; i++) { points(n, i) = std::stof(tokens[i]); }
  }
  // get the lumen tree segments
  getline(lumen_file, line);
  boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
  segments_count = std::stoi(tokens[0]);
  segments.resize(segments_count, Eigen::NoChange);
  for (int n = 0; n < segments_count; n++) {
    getline(lumen_file, line);
    boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
    for (int i = 0; i < 2; i++) {
      segments(n, i) = std::stoi(tokens[i]) - 1; // change to zero Treed indexing
    }
  }
  lumen_file.close();
  print_info();
}

double cLumenTree::get_dnl(Eigen::Vector3d p)
{
  double d = 100.0; // large dummy initial distance
  Eigen::Vector3d w, v;
  for (int n = 0; n < segments_count; n++) {
    v = points.block<1, 3>(segments(n, 0), 0);
    w = points.block<1, 3>(segments(n, 1), 0);
    d = std::min(d, utils::get_distance(p, w, v));
  }
  return (d);
}

void cLumenTree::print_info()
{
  parent->out << "<LumenTree> number of points: " << points_count << std::endl;
  parent->out << "<LumenTree> number of segments: " << segments_count << std::endl;
}