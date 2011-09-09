#ifndef SCYLLOP_LP_H
#define SCYLLOP_LP_H

#include <vector>

#include "rational.h"
#include "scyllop_classes.h"

enum scyllop_lp_solver {GLPK_DOUBLE, GLPK_EXACT, QSOPT_EXACT, EXLP};


void scyllop_lp(CyclicProduct& G, 
                Chain& C, 
                std::vector<Multiarc>& arcs, 
                std::vector<Polygon>& polys, 
                rational* scl, 
                std::vector<rational>* solution_vector, 
                scyllop_lp_solver solver, 
                int VERBOSE); 

#endif 
