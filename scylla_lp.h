#ifndef SCYLLA_LP_H
#define SCYLLA_LP_H

#include <vector>

#include "rational.h"
#include "scylla_classes.h"

enum scylla_lp_solver {GLPK_SIMPLEX, GLPK_IPT, EXLP};


void scylla_lp(Chain& C, 
                InterfaceEdgeList &IEL ,
                CentralEdgePairList &CEL, 
                std::vector<CentralPolygon> &CP,
                std::vector<GroupTooth> &GT,
                std::vector<GroupRectangle> &GR,
                rational* scl, 
                std::vector<rational>* solution_vector, 
                scylla_lp_solver solver, 
                int VERBOSE,
                int LP_VERBOSE); 

#endif 
