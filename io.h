#ifndef __io__
#define __io__

#include <vector>
#include <string>

#include "scallop.h"
#include "rational.h"

void write_linear_program(std::vector<std::string>& w, 
                          std::vector<int>& weight, 
                          std::vector<arc>& arc_list, 
                          std::vector<polygon>& polygon_list, 
                          std::string& programFile);

void write_solution_polyhedron(std::vector<std::string>& w, 
                               std::vector<int>& weight, 
                               std::vector<arc>& arc_list, 
                               std::vector<polygon>& polygon_list, 
                               std::string& polyFileName,
                               rational scl);

#endif
