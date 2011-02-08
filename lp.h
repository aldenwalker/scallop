#ifndef __lp__
#define __lp__


#include <vector>

#include "scallop.h"
#include "rational.h"


using namespace std;

enum scallop_lp_solver {GLPK_DOUBLE, GLPK_EXACT, QSOPT_EXACT, EXLP};

void do_linear_program(   vector<string> &w,
                          vector<int>& weight,
                          vector<arc> &arc_list, 
                          vector<polygon> &polygon_list, 
                          rational &scl, 
                          vector<rational> &solutionVector,
                          scallop_lp_solver solver, 
                          int VERBOSE);


#endif
