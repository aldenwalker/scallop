#ifndef __io__
#define __io__

#include <vector>
#include <string>

#include "scallop.h"

void write_linear_program(std::vector<std::string>& w, 
                          std::vector<int>& weight, 
                          std::vector<arc>& arc_list, 
                          std::vector<polygon>& polygon_list, 
                          std::string& programFile);

#endif
