#ifndef SCYLLA_LP_H
#define SCYLLA_LP_H

#include <vector>

#include "rational.h"
#include "scylla_classes.h"

enum scylla_lp_solver {GLPK_DOUBLE, GLPK_EXACT, QSOPT_EXACT, EXLP};


void scylla_lp(Chain& C, 
                std::vector<GroupEdgeList> &GEL, 
                InterfaceEdgeList &IEL ,
                CentralEdgeList &CEL, 
                std::vector<CentralPolygon> &CP,
                std::vector<std::vector<GroupTooth> > &GT,
                std::vector<std::vector<GroupMouth> > &GM,
                std::vector<std::vector<GroupPolygon> > &GP,
                std::vector<std::vector<GroupRectangle> > &GR,
                rational* scl, 
                std::vector<rational>* solution_vector, 
                scylla_lp_solver solver, 
                bool VERBOSE); 

#endif 
