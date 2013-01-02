#include <vector>
#include <string>
#include "trollop_classes.h"
#include "../rational.h"
#include "../lp.h"

namespace TROLLOP {

  void read_matrix(std::vector<std::vector<int> >& matrix, 
                   std::string matrix_filename, 
                   WordTable &WT,
                   bool rearrange_rows,
                   int VERBOSE);
  
  void read_vector(std::vector<int>& vector, 
                   std::string vector_filename,
                   int VERBOSE);
  
  void compute_triangles(WordTable& WT, 
                         ArcPairList& AL, 
                         int num_copies,
                         std::vector<Triangle>& TR);
  
  void print_triangles(std::vector<Triangle>& TR, std::ostream& os);
  
  void compute_rectangles(WordTable& WT, 
                          ArcPairList& AL, 
                          int num_copies, 
                          std::vector<Rectangle>& RE);
  
  void print_rectangles(std::vector<Rectangle>& RE, std::ostream& os) ;
  
  
  
  
  void collect_dups_and_push(std::vector<int> &temp_ia,
                                    std::vector<int> &temp_ja,
                                    std::vector<int> &temp_ar,
                                    std::vector<int> &ia,
                                    std::vector<int> &ja,
                                    std::vector<int> &ar,
                           int VERBOSE);

  void write_lp(WordTable& WT, 
                WVec& C, 
                ArcPairList& AL, 
                int num_copies,
                std::vector<Triangle>& TR, 
                std::vector<Rectangle>& RE, 
                bool DO_SUP,
                bool MAT_COMP,
                bool MAT_SEPARATE_DOMAIN,
                std::vector<std::vector<int> >& M,
                std::vector<std::vector<int> >& N,
                std::vector<int>& b,
                int VERBOSE,
                int LP_VERBOSE,
                std::string& programFile);
  
  void trollop_lp(WordTable& WT, 
                  WVec& C, 
                  ArcPairList& AL, 
                  int num_copies,
                  std::vector<Triangle>& TR, 
                  std::vector<Rectangle>& RE, 
                  bool DO_SUP,
                  bool MAT_COMP,
                  bool MAT_SEPARATE_DOMAIN,
                  std::vector<std::vector<int> >& M,
                  std::vector<std::vector<int> >& N,
                  std::vector<int>& b,
                  Rational& ans,
                  std::vector<Rational>& solution_vector,
                  SparseLPSolver solver,
                  int VERBOSE,
                  int LP_VERBOSE);
  
  int trollop(int argc, char* argv[]);


}