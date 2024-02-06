#ifndef scylla_H
#define scylla_H

namespace SCYLLA {
  
// void compute_group_teeth_and_rectangles(Chain &C,
//                                         std::vector<GroupTooth > &GT,
//                                         std::vector<GroupRectangle > &GR);
  
//   void compute_central_polys(Chain &C,
//                              InterfaceEdgeList &IEL, 
//                              std::vector<CentralPolygon> &CP,
//                              bool limit_central_sides);
  
//   void print_central_polys(std::vector<CentralPolygon> &CP, 
//                                    std::ostream &os, 
//                            int level);
  
//   void print_group_teeth_and_rectangles(std::vector<GroupTooth> &GT,
//                                                 std::vector<GroupRectangle> &GR,
//                                                 std::ostream &os,
//                                         int level);
  
//   void write_solution_to_fatgraph(std::string& output_filename,
//                                   Chain& C,
//                                   InterfaceEdgeList& IEL,
//                                   CentralEdgePairList &CEL, 
//                                   std::vector<CentralPolygon> &CP,
//                                   std::vector<GroupTooth> &GT,
//                                   std::vector<GroupRectangle> &GR,
//                                   std::vector<Rational>& solution_vector,
//                                   int verbose );
  
  
//   void scylla_lp(Chain& C, 
//                          InterfaceEdgeList &IEL,
//                          CentralEdgePairList &CEL, 
//                          std::vector<CentralPolygon> &CP,
//                          std::vector<GroupTooth> &GT,
//                          std::vector<GroupRectangle> &GR,
//                          Rational* scl, 
//                          std::vector<Rational>* solution_vector, 
//                          SparseLPSolver solver, 
//                          bool WRITE_LP,
//                          std::string LP_filename,
//                          bool cl_not_scl,
//                          int VERBOSE,
//                          int LP_VERBOSE);
  
  
void scylla(int argc, char** argv);

}
#endif