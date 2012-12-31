#ifndef LP_H
#define LP_H

#include <vector>
#include "rational.h"

enum SparseLPEqualityType {EQ, LE, GE};

enum SparseLPSolver {GLPK, 
                      GLPK_SIMPLEX, 
                      GLPK_IPT, 
                      GUROBI, 
                      GUROBI_SIMPLEX, 
                      GUROBI_IPT, 
                      EXLP};

enum SparseLPSolveCode {LP_OPTIMAL, LP_INFEASIBLE, LP_ERROR};

class SparseLP {

private:
  std::vector<int> ia;
  std::vector<int> ja;
  std::vector<int> ar;
  std::vector<double> double_ar;
  std::vector<int> objective;
  std::vector<double> double_objective;
  std::vector<int> RHS;
  std::vector<double> double_RHS;
  std::vector<Rational> soln_vector;
  std::vector<double> double_soln_vector;
  Rational op_val;
  double double_op_val;
  std::vector<SparseLPEqualityType> eq_type;

  int num_cols;
  int num_rows;
  
  SparseLPSolver solver;

public:
  
  SparseLP(SparseLPSolver s);
  SparseLP(SparseLPSolver s, int nr, int nc);
  void set_num_rows(int nr);
  void set_num_cols(int nc);
  void add_entry(int i, int j, int a);
  void add_entry(int i, int j, double a);
  void add_entry_catch_dups(int i, int j, int a, int how_far_back);
  void add_entry_catch_dups(int i, int j, double a, int how_far_back);
  void set_obj(int i, int v);
  void set_obj(int i, double v);
  void set_RHS(int i, int r);
  void set_RHS(int i, double r);
  void set_equality_type(int i, SparseLPEqualityType et);
  void get_soln_vector(std::vector<double>& sv);
  void get_soln_vector(std::vector<Rational>& sv);
  void get_optimal_value(double& ov);
  void get_optimal_value(Rational& ov);
  
  SparseLPSolveCode solve(int verbose);
  
};

#endif