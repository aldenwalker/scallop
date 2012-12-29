#include <vector>
#include <iostream>

SparseLP::SparseLP(SparseLPSolver s) {
  ia.resize(0);
  ja.resize(0);
  ar.resize(0);
  double_ar.resize(0);
  objective.resize(0);
  double_objective.resize(0);
  RHS.resize(0);
  double_RHS.resize(0);
  eq_type.resize(0);
  num_rows = 0;
  num_cols = 0;
  solver = s;
  ia.push_back(-1);
  ja.push_back(-1);
  if (s == EXLP) {
    ar.push_back(-1);
  } else {
    double_ar.push_back(-1);
  }
}

SparseLP::SparseLP(SparseLPSolver s, int nr, int nc) {
  ia.resize(0);
  ja.resize(0);
  ar.resize(0);
  double_ar.resize(0);
  ia.push_back(-1);
  ja.push_back(-1);
  if (s == EXLP) {
    ar.push_back(-1);
    objective.resize(nc);
    double_objective.resize(0);
    RHS.resize(nr);
    double_RHS.resize(0);
  } else {
    double_ar.push_back(-1);
    objective.resize(0);
    double_objective.resize(nc);
    RHS.resize(0);
    double_RHS.resize(nr);
  }
  eq_type.resize(nr);
  num_rows = nr;
  num_cols = nc;
  solver = s;
}

void SparseLP::set_num_rows(int nr) {
  if (solver == EXLP) {
    RHS.resize(nr);
    double_RHS.resize(0);
  } else {
    RHS.resize(0);
    double_RHS.resize(nr);
  }
  equality_type.resize(nr);
  num_rows = nr;
}

void SparseLP::set_num_cols(int nc) {
  if (solver == EXLP) {
    objective.resize(nc);
    double_objective.resize(0);
  } else {
    objective.resize(0);
    double_objective.resize(nc);
  }
  num_cols = nc;
}

void SparseLP::add_entry(int i, int j, int a) {
  if (solver == EXLP) {
    ia.push_back(i);
    ja.push_back(j);
    ar.push_back(a);
  } else {
    ia.push_back(i);
    ja.push_back(j);
    double_ar.push_back((double)a);
  }
}
    
void SparseLP::add_entry(int i, int j, double a) {
  if (solver == EXLP) {
    std::cout << "Can't input a double for rational LP\n";
  } else {
    ia.push_back(i);
    ja.push_back(j);
    double_ar.push_back(a);
  }
}

void SparseLP::add_entry_catch_dups(int i, int j, int a, int how_far_back) {
  if (solver == EXLP) {
    ia.push_back(i);
    ja.push_back(j);
    ar.push_back(a);
  } else {
    ia.push_back(i);
    ja.push_back(j);
    double_ar.push_back((double)a);
  }
}

void SparseLP::add_entry_catch_dups(int i, int j, double a, int how_far_back) {
  if (solver == EXLP) {
    std::cout << "Can't input a double for rational LP\n";
  } else {
    ia.push_back(i);
    ja.push_back(j);
    double_ar.push_back(a);
  }
}  

void SparseLP::set_obj(int i, int v) {
  if (solver == EXLP) {
    objective[i] = v;
  } else {
    double_objective[i] = (double)v;
  }
}

void SparseLP::set_obj(int i, double v) {
  if (solver == EXLP) {
    std::cout << "Can't input a double for rational LP\n";
  } else {
    double_objective[i] = v;
  }
}

void SparseLP::set_RHS(int i, int r) {
  if (solver == EXLP) {
    RHS[i] = r;
  } else {
    double_RHS[i] = (double)r;
  }
}

void SparseLP::set_RHS(int i, double r) {
  if (solver == EXLP) {
    std::cout << "Can't input a double for rational LP\n";
  } else {
    double_RHS[i] = r;
  }
}

void SparseLP::set_equality_type(int i, SparseLPEqualityType et) {
  eq_type[i]= et;
}

SparseLPSolveCode SparseLP::solve(int verbose) {
  return LP_OPTIMAL;
}




