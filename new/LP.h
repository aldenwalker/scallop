#include <vector>


enum {EQ, LE, GE} SparseLPEqualityType;

enum {GLPK, GUROBI, EXLP} SparseLPSolver;

enum {LP_OPTIMAL, LP_INFEASIBLE, LP_ERROR} SparseLPSolveCode;

class SparseLP {

private:
  vector<int> ia;
  vector<int> ja;
  vector<int> ar;
  vector<double> double_ar;
  vector<int> objective;
  vector<double> double_objective;
  vector<int> RHS;
  vector<double> double_RHS;
  vector<SparseLPEqualityType> eq_type;


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
  
  SparseLPSolveCode solve(int verbose);
  
  
};