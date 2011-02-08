#ifndef __draw__
#define __draw__

#include <vector>
#include <string>

#include "scallop.h"
#include "rational.h"

using namespace std;

void make_graph_and_print(  string drawFile,
                            vector<string> chain, 
                            vector<arc> &arc_list, 
                            vector<polygon>& poly_list, 
                            vector<rational>& weights); 

#endif
