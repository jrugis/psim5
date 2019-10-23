/*
 * cCell_flow.cpp
 *
 *	Created on: 16/12/2018
 *	Author: jrugis
 */

#include "cCell_flow.hpp"

cCell_flow::cCell_flow(int index, double parms[], std::ofstream& out)
{
  out << "<Cell_flow> index: " << index << std::endl;
  my_index = index;
  p = parms;
}

cCell_flow::~cCell_flow() {}
