#ifndef scylla_H
#define scylla_H
#include <vector>
#include "scylla_classes.h"

namespace SCYLLA {
  
void compute_group_teeth_and_rectangles(Chain &C,
                                        std::vector<GroupTooth > &GT,
                                        std::vector<GroupRectangle > &GR);
  
  void compute_central_polys(Chain &C,
                             InterfaceEdgeList &IEL, 
                             std::vector<CentralPolygon> &CP,
                             bool limit_central_sides);
  
  void print_central_polys(std::vector<CentralPolygon> &CP, 
                                   std::ostream &os, 
                           int level);
  
  void print_group_teeth_and_rectangles(std::vector<GroupTooth> &GT,
                                                std::vector<GroupRectangle> &GR,
                                                std::ostream &os,
                                        int level);
  
  void scylla_lp(Chain& C, 
                         InterfaceEdgeList &IEL,
                         CentralEdgePairList &CEL, 
                         std::vector<CentralPolygon> &CP,
                         std::vector<GroupTooth> &GT,
                         std::vector<GroupRectangle> &GR,
                         Rational* scl, 
                         std::vector<Rational>* solution_vector, 
                         SparseLPSolver solver, 
                         bool WRITE_LP,
                         std::string LP_filename,
                         int VERBOSE,
                 int LP_VERBOSE);
  
  
void scylla(int argc, char** argv);

}
#endif