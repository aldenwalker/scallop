#ifndef LP_H
#define LP_H

#include <vector>
#include <string>
#include "rational.h"

enum SparseLPEqualityType {EQ, LE, GE};

enum SparseLPColumnBoundType {LB, UB, FIX, FREE};

enum SparseLPColumnType {REAL, INT};

enum SparseLPSolver {GLPK, 
                      GLPK_SIMPLEX, 
                      GLPK_IPT, 
                      GUROBI, 
                      GUROBI_SIMPLEX, 
                      GUROBI_IPT, 
                      EXLP};

enum SparseLPSolveCode {LP_OPTIMAL, LP_INFEASIBLE, LP_ERROR, LP_TIME_LIMIT};

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
  std::vector<SparseLPColumnType> col_type;
  int num_ints;
  std::vector<SparseLPColumnBoundType> col_bound_types;
  std::vector<double> col_bounds_double;
  std::vector<int> col_bounds;

  int num_cols;
  int num_rows;
  
  SparseLPSolver solver;

public:
  
  SparseLP(SparseLPSolver s);
  SparseLP(SparseLPSolver s, int nr, int nc);
  void write_to_file(std::string filename);
  void set_num_rows(int nr);
  void set_num_cols(int nc);
  void set_col_type(int c, SparseLPColumnType t);
  void set_col_bound(int c, SparseLPColumnBoundType t);
  void set_col_bound(int c, SparseLPColumnBoundType t, double b);
  void set_col_bound(int c, SparseLPColumnBoundType t, int b);
  void add_entry(int i, int j, Rational& r);
  void add_entry(int i, int j, int a);
  void add_entry(int i, int j, double a);
  void extend_entries_no_dups(std::vector<int>& temp_ia, 
                              std::vector<int>& temp_ja,
                              std::vector<int>& temp_ar);
  void extend_entries_no_dups(std::vector<int>& temp_ia, 
                              std::vector<int>& temp_ja,
                              std::vector<double>& temp_ar);
  void set_obj(int i, int v);
  void set_obj(int i, double v);
  void set_RHS(int i, Rational& r);
  void set_RHS(int i, int r);
  void set_RHS(int i, double r);
  void set_equality_type(int i, SparseLPEqualityType et);
  int get_num_entries();
  void reset_num_entries(int i);
  void get_soln_vector(std::vector<double>& sv);
  void get_soln_vector(std::vector<long double>& sv);
  void get_soln_vector(std::vector<Rational>& sv);
  void get_optimal_value(double& ov);
  void get_optimal_value(Rational& ov);
  
  SparseLPSolveCode solve(int verbose);
  
  void print_LP();
  
};

#endif