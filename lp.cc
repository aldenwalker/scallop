#include <vector>
#include <string>
#include <iostream>

#include <glpk.h>

#ifdef GUROBI_INSTALLED
extern "C" {
#include <gurobi_c.h>
}
#endif

#include "lp.h"
#include "rational.h"

extern "C" {
#include "exlp-package/lpstruct.h"
}
extern "C" {
#include "exlp-package/solve_lp.h"
}
extern "C" {
#include "exlp-package/mylib.h"
}

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
  soln_vector.resize(0);
  double_soln_vector.resize(0);
  num_rows = 0;
  num_cols = 0;
  solver = s;
  op_val = Rational(-1,1);
  double_op_val = -1;
}

SparseLP::SparseLP(SparseLPSolver s, int nr, int nc) {
  ia.resize(0);
  ja.resize(0);
  ar.resize(0);
  double_ar.resize(0);
  if (s == EXLP) {
    objective.resize(nc);
    double_objective.resize(0);
    RHS.resize(nr);
    double_RHS.resize(0);
    soln_vector.resize(nc);
    double_soln_vector.resize(0);
  } else {
    objective.resize(0);
    double_objective.resize(nc);
    RHS.resize(0);
    double_RHS.resize(nr);
    soln_vector.resize(0);
    double_soln_vector.resize(nc);
  }
  eq_type.resize(nr);
  num_rows = nr;
  num_cols = nc;
  solver = s;
  op_val = Rational(-1,1);
  double_op_val = -1;
  //std::cout << "Made new LP problem with solver: " << solver << "\n";
}

void SparseLP::write_to_file(std::string filename) {
}

void SparseLP::set_num_rows(int nr) {
  if (solver == EXLP) {
    RHS.resize(nr);
    double_RHS.resize(0);
  } else {
    RHS.resize(0);
    double_RHS.resize(nr);
  }
  eq_type.resize(nr);
  num_rows = nr;
}

void SparseLP::set_num_cols(int nc) {
  if (solver == EXLP) {
    objective.resize(nc);
    double_objective.resize(0);
    soln_vector.resize(nc);
    double_soln_vector.resize(0);
  } else {
    objective.resize(0);
    double_objective.resize(nc);
    soln_vector.resize(0);
    double_soln_vector.resize(nc);
  }
  num_cols = nc;
}


void SparseLP::add_entry(int i, int j, Rational& r) {
  if (solver == EXLP) {
    if (r.d() != 1) {
      std::cout << "You can give a rational entry, but it needs to be an integer\n";
      return;
    }
    ia.push_back(i);
    ja.push_back(j);
    ar.push_back(r.n());
  } else {
    ia.push_back(i);
    ja.push_back(j);
    double_ar.push_back(r.get_d());
  }
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


void SparseLP::extend_entries_no_dups(std::vector<int> &temp_ia,
                                       std::vector<int> &temp_ja,
                                       std::vector<int> &temp_ar) {
  int j, k, temp;
  for (j=0; j<(int)temp_ia.size(); j++) {
    if (temp_ar[j] == 0) {
      continue;
    }
    temp = 0;
    for (k=0; k<(int)temp_ia.size(); k++) {
      if (temp_ia[k] == temp_ia[j]) {
        temp += temp_ar[k];
        temp_ar[k] = 0;
      }
    }
    add_entry(temp_ia[j], temp_ja[j], temp);
    //ja.push_back(temp_ja[j]);
    //ia.push_back(temp_ia[j]);
    //ar.push_back(temp);
  }
}

void SparseLP::extend_entries_no_dups(std::vector<int>& temp_ia, 
                                      std::vector<int>& temp_ja,
                                      std::vector<double>& temp_ar) {
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


void SparseLP::set_RHS(int i, Rational& r) {
  if (solver == EXLP) {
    if (r.d() != 1) {
      std::cout << "You can give a rational RHS, but it needs to be an integer\n";
      return;
    }
    RHS[i] = r.n();
  } else {
    double_RHS[i] = r.get_d();
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


int SparseLP::get_num_entries() {
  if (solver == EXLP) {
    return (int)ar.size();
  } else {
    return (int)double_ar.size();
  }
}


void SparseLP::reset_num_entries(int i) {
  ia.resize(i);
  ja.resize(i);
  if (solver == EXLP) {
    ar.resize(i);
  } else {
    double_ar.resize(i);
  }
}


void SparseLP::get_soln_vector(std::vector<double>& sv) {
  sv.resize(num_cols);
  if (solver == EXLP) {
    for (int i; i<num_cols; ++i) {
      sv[i] = soln_vector[i].get_d();
    }
  } else {
    for (int i; i<num_cols; ++i) {
      sv[i] = double_soln_vector[i];
    }
  }
}

void SparseLP::get_soln_vector(std::vector<Rational>& sv) {
  sv.resize(num_cols);
  if (solver == EXLP) {
    for (int i=0; i<num_cols; ++i) {
      sv[i] = soln_vector[i];
    }
  } else {
    for (int i=0; i<num_cols; ++i) {
      sv[i] = approx_rat(double_soln_vector[i]);
    }
  }
}

void SparseLP::get_optimal_value(double& ov) {
  if (solver == EXLP) {
    //std::cout << "Getting double optimal value from EXLP?\n";
    ov = op_val.get_d();
  } else {
    ov = double_op_val;
  }
}

void SparseLP::get_optimal_value(Rational& ov) {
  if (solver != EXLP) {
    //std::cout << "Getting rational optimal value from non-EXLP?\n";
    ov = approx_rat(double_op_val);
  } else {
    ov = op_val;
  }
}

SparseLPSolveCode SparseLP::solve(int verbose) {
  
  /************************************  GLPK *******************************/ 
  
  if (solver == GLPK || solver == GLPK_SIMPLEX || solver == GLPK_IPT) {   
    
    glp_prob *lp;
    glp_smcp parm;
    glp_iptcp ipt_parm;
    
    lp = glp_create_prob();
    
    glp_set_prob_name(lp, "scl");
    glp_set_obj_dir(lp, GLP_MIN);
    
    glp_add_rows(lp, num_rows );
    
    glp_add_cols(lp, num_cols);
    
    for (int i=0; i<num_rows; i++) {
      glp_set_row_bnds(lp, i+1, GLP_FX, double_RHS[i], double_RHS[i]);
    }
    for (int i=0; i<num_cols; i++) {
      glp_set_col_bnds(lp, i+1, GLP_LO, 0.0, 0.0);
      glp_set_obj_coef(lp, i+1, double_objective[i]);
    }
    //rearrange
    ia.push_back(0);
    ja.push_back(0);
    double_ar.push_back(0);
    for (int i=ia.size()-1; i>0; --i) {
      ia[i] = ia[i-1]+1;
      ja[i] = ja[i-1]+1;
      double_ar[i] = double_ar[i-1];
    }
    glp_load_matrix(lp, ia.size()-1, &ia[0], &ja[0], &double_ar[0]);
    //unrearrange
    for (int i=0; i<(int)ia.size()-1; ++i) {
      ia[i] = ia[i+1]-1;
      ja[i] = ja[i+1]-1;
      double_ar[i] = double_ar[i+1];
    }    
    ia.pop_back();
    ja.pop_back();
    double_ar.pop_back();
    
    if (solver == GLPK || solver == GLPK_SIMPLEX) {
      glp_init_smcp(&parm);
      parm.presolve=GLP_ON;
      if (verbose > 1) {
        parm.msg_lev = GLP_MSG_ALL;
      } else {
        parm.msg_lev = GLP_MSG_OFF;
      }
      glp_simplex(lp, &parm);
    } else if (solver == GLPK_IPT) {
      glp_init_iptcp(&ipt_parm);
      if (verbose > 1) {
        ipt_parm.msg_lev = GLP_MSG_ALL;
      } else {
        ipt_parm.msg_lev = GLP_MSG_OFF;
      }
      glp_interior(lp, &ipt_parm);
    }
    
    int stat = glp_get_status(lp);
    if (stat != GLP_OPT) {
      if (stat == GLP_NOFEAS || stat == GLP_UNDEF) {
        glp_delete_prob(lp);
        return LP_INFEASIBLE;
      } else {
        std::cout << "GLPK Linear programming error: " << stat << "\n";
      }
      glp_delete_prob(lp);
      return LP_ERROR;
    }
    
    double_op_val = glp_get_obj_val(lp)/4.0;	
    
    double_soln_vector.resize(num_cols);
    for (int i=0; i<num_cols; i++) {
      double_soln_vector[i] = glp_get_col_prim(lp,i+1);
    }	
    
	  glp_delete_prob(lp);
	  
  /***************************************  EXLP ****************************/  
    
	  
	} else if (solver == EXLP) {
    
    //exlp init
	  mylib_init();
	  
	  LP* lp;
    int  result;
    char buf[100];
    int varNum;
    
    if (verbose>1) 
      std::cout << "About to create a new lp\n";    
    lp = new_lp(NULL);
    
    if (verbose>1) 
      std::cout << "Done\n";
    
    if (verbose>1)
      std::cout << "Init hash\n";
    
    lp_hash_str_init(lp, lp->hash_entries);
    //my_hash_mpq_init(lp->hash_entries);
    
    if (verbose>1)
      std::cout << "Done\n";
    
    
    sprintf(buf, "scl");
    lp_set_name(lp, buf);
    
    //the lp by default (set in lpstruct.c)
    //has all the right stuff, I think
    
    //now we input the rows -- it likes to name them
    for (int i=0; i<num_rows; i++) {
      sprintf(buf, "r%d", i);
      lp_add_row(lp, buf);
      switch (eq_type[i]) {
        case EQ:
          lp_set_row_equality(lp, lp_get_row_num(lp, buf), 'E');
          break;
        case LE:
          lp_set_row_equality(lp, lp_get_row_num(lp, buf), 'L');
          break;
        case GE:
          lp_set_row_equality(lp, lp_get_row_num(lp, buf), 'G');
          break;
      }
    }
    
    if (verbose>1) {
      std::cout << "Entered the row names and equalities\n";
    }
    
    //add the objective function row
    sprintf(buf, "obj");
    lp_set_obj_name(lp, buf);
    
    int* columnIndices = new int[num_cols]; //this is probably useless
    int rowNum;
    int objectiveIndex = lp_get_row_num(lp, (char*)"obj");
    mpq_t entry;
    mpq_init(entry);
    
    for (int i=0; i<num_cols; i++) {
      //since we only enter entries from a column once, we only need to do this once,
      //as opposed to read_columns in mps.c
      sprintf(buf, "col%d", i);
      varNum = lp_add_var_without_A(lp, buf);
      columnIndices[i] = varNum;
      
      //for (j=0; j<nRows; j++) {
      //we're going to enter something in the varNum column, in the right row
      //with coefficient from constraints, but we will know what to do 
      //}
    }
    
    if (verbose>1) {
      std::cout << "Added the columns\n";
    }
    
    //this is from mps.c
    matrix_resize(lp->A, lp->rows, lp->vars);
    vector_resize(lp->b, lp->rows);
    vector_resize(lp->xb, lp->rows);
    vector_resize(lp->cb, lp->rows);
    
    
    //now we actually enter the data
    for (int i=0; i<(int)ia.size(); ++i) {
      varNum = columnIndices[ja[i]];
      sprintf(buf, "r%d", ia[i]);
      rowNum = lp_get_row_num(lp, buf);
      mpq_set_si(entry, ar[i], 1);
      lp_set_coefficient(lp, entry, rowNum, varNum);
    }
    
    //set the objective function
    for (int i=0; i<num_cols; ++i) {
      varNum = columnIndices[i];
      mpq_set_si(entry, objective[i], 1);
      lp_set_coefficient(lp, entry, objectiveIndex, varNum);
    }
    
    //lp->maximize is false, so we reverse the sign
    vector_rev_sgn(lp->c);
    
    //set the right hand sides for the arcs
    for (int i=0; i<num_rows; ++i) {
      sprintf(buf, "r%d", i);
      rowNum = lp_get_row_num(lp, buf);
      mpq_set_si(entry, RHS[i], 1);
      lp_set_rhs(lp, rowNum, entry);
    }
    
    //read the bounds on the columns
    mpq_set_si(entry, 0,1);
    for (int i=0; i<num_cols; i++) {
      varNum = columnIndices[i];
      if (lp->lower.is_valid[varNum] == FALSE) //this will always be the case I think
        mpq_init(lp->lower.bound[varNum]);
      mpq_set(lp->lower.bound[varNum], entry);
      lp->lower.is_valid[varNum] = TRUE;
    }
    
    if (verbose>1) {
      std::cout << "Rows: " << lp->rows << "\n";;
      std::cout << "Vars: " << lp->vars << "\n";
    }
    
    
    result = solve_lp(lp);
    
    if (result != LP_RESULT_OPTIMAL) {
      //std::cout << "got error code " << result << "\n";
      if (result == 2) {
        lp_free(lp);
        return LP_INFEASIBLE;
      } else {
        lp_free(lp);
        return LP_ERROR;
      }
    }
    
    lp_get_object_value(lp, &entry);
    
    op_val = Rational(entry)/Rational(4,1);
    
    for (int i=0; i<num_cols; i++) {
      mpq_set(entry, *vector_get_element_ptr(lp->x, columnIndices[i]));
      soln_vector[i] = Rational(entry);
    }
    
    lp_free(lp);
    
    
    
 /************************  GUROBI *******************************************/   
    
    

  } else if (solver == GUROBI || 
             solver == GUROBI_SIMPLEX || 
             solver == GUROBI_IPT ) {
    
#ifndef GUROBI_INSTALLED

    std::cout << "Not compiled with gurobi support\n";
    return LP_ERROR;

#else
    
    GRBenv   *env   = NULL;
    GRBmodel *model = NULL;
    GRBloadenv( &env, "gurobi.log" );
    
    //create a new model and immediately load in all the columns
    GRBnewmodel( env, &model, "scl", num_cols, &double_objective[0], NULL, NULL, NULL, NULL);
    
    //add the constraints (rows)  here we make them empty equality rows and fix the RHS
    for (int i=0; i<num_rows; i++) {
      switch (eq_type[i]) {
        case EQ: 
          GRBaddconstr( model, 0, NULL, NULL, GRB_EQUAL, double_RHS[i], NULL);
          break;
        case LE:
          GRBaddconstr( model, 0, NULL, NULL, GRB_LESS_EQUAL, double_RHS[i], NULL);
          break;
        case GE:
          GRBaddconstr( model, 0, NULL, NULL, GRB_GREATER_EQUAL, double_RHS[i], NULL);
          break;
      } 
    }
    GRBupdatemodel(model);
    
    //add the matrix:
    
    //first fix the matrix entries, because all the rows and columns are 1-based
    GRBchgcoeffs( model, (int)ia.size(), &ia[0], &ja[0], &double_ar[0] );
    
    //set the correct optimization method
    if (solver == GUROBI_SIMPLEX) {
      GRBsetintparam(GRBgetenv(model), "Method", 1);
    } else if (solver == GUROBI_IPT) {
      GRBsetintparam(GRBgetenv(model), "Method", 2);
    } else { 
      //do nothing -- default
    }
    
    //set the output
    GRBsetintparam(GRBgetenv(model), "OutputFlag", (verbose > 1 ? 1 : 0) );
    GRBsetintparam(GRBgetenv(model), "DisplayInterval", 60 );
    
    
    //optimize
    GRBoptimize( model );
    
    //determine if the problem had a solution
    int problem_status;
    GRBgetintattr( model, GRB_INT_ATTR_STATUS, &problem_status );
    if (problem_status != GRB_OPTIMAL) {
      if (problem_status == GRB_TIME_LIMIT) {
        std::cout << "Time limit\n";
        GRBfreemodel(model);
        GRBfreeenv(env);
        return LP_TIME_LIMIT;
      } else if (problem_status != GRB_INF_OR_UNBD && 
                 problem_status != GRB_INFEASIBLE &&
                 problem_status != GRB_UNBOUNDED) {
        std::cout << "Gurobi Linear programming error\n";
        GRBfreemodel(model);
        GRBfreeenv(env);
        return LP_ERROR;
      } else {
        GRBfreemodel(model);
        GRBfreeenv(env);
        return LP_INFEASIBLE;
      }
    }
    
    
    //get the objective minimum
    GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &double_op_val);
    double_op_val /= 4.0;
    
    //get the solution vector
    double_soln_vector.resize(num_cols);
    GRBgetdblattrarray( model, GRB_DBL_ATTR_X, 0, num_cols, &double_soln_vector[0] );
    
    //free stuff
    GRBfreemodel(model);
    GRBfreeenv(env);
#endif
  }
  
  return LP_OPTIMAL;
}




void SparseLP::print_LP() {
  if (solver != EXLP) {
    std::cout << "Objective: ";
    for (int i=0; i<num_cols; ++i) {
      std::cout << double_objective[i] << " ";
    }
    std::cout << "\n";
    std::cout << "RHS: ";
    for (int i=0; i<num_rows; ++i) {
      std::cout << double_RHS[i] << " ";
    }
    std::cout << "\nEntries: ";
    for (int i=0; i<(int)ia.size(); ++i) {
      std::cout << "(" << ia[i] << "," << ja[i] << "," << double_ar[i] << "), ";
    }
    std::cout << "\n";
  } else {
    std::cout << "Objective: ";
    for (int i=0; i<num_cols; ++i) {
      std::cout << objective[i] << " ";
    }
    std::cout << "\n";
    std::cout << "RHS: ";
    for (int i=0; i<num_rows; ++i) {
      std::cout << RHS[i] << " ";
    }
    std::cout << "\nEntries: ";
    for (int i=0; i<(int)ia.size(); ++i) {
      std::cout << "(" << ia[i] << "," << ja[i] << "," << ar[i] << "), ";
    }
    std::cout << "\n";
  }
}










