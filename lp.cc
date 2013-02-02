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
  col_type.resize(0);
  num_ints = 0;
  col_bounds.resize(0);
  col_bounds_double.resize(0);
  col_bound_types.resize(0);
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
  col_type.resize(0);
  num_ints = 0;
  col_bounds.resize(0);
  col_bounds_double.resize(0);
  col_bound_types.resize(0);
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

void SparseLP::set_col_type(int c, SparseLPColumnType t) {
  if (t == REAL && num_ints == 0) {
    return; //no need to do anything
  }
  if (t == REAL) {
    if (num_ints == 1) {
      num_ints = 0;
      col_type.resize(0);
    } else {
      num_ints += (col_type[c] == INT ? -1 : 0);
      col_type[c] = REAL;
    }
  } else { //t==INT
    if (num_ints == 0) {
      col_type.resize(num_cols, REAL);
      col_type[c] = INT;
      num_ints = 1;
    } else {
      num_ints += (col_type[c] == REAL ? 1 : 0);
      col_type[c] = INT;
    }
  }
}


void SparseLP::set_col_bound(int c, SparseLPColumnBoundType t) {
  if ((int)col_bounds.size() == 0) {
    if (t == LB) return;
    col_bound_types.resize(num_cols, LB);
    if (solver == EXLP) {
      col_bounds.resize(num_cols, 0);
    } else {
      col_bounds_double.resize(num_cols, 0);
    }
  }
  col_bound_types[c] = t;
}

void SparseLP::set_col_bound(int c, SparseLPColumnBoundType t, double b) {
  if ((int)col_bounds.size() == 0) {
    if (t == LB) return;
    col_bound_types.resize(num_cols, LB);
    if (solver == EXLP) {
      col_bounds.resize(num_cols, 0);
    } else {
      col_bounds_double.resize(num_cols, 0);
    }
  }
  col_bound_types[c] = t;
  if (solver == EXLP) {
    std::cout << "Can't input a double column bound for rational LP\n";
  } else {
    col_bounds_double[c] = b;
  }
}

void SparseLP::set_col_bound(int c, SparseLPColumnBoundType t, int b) {
  if ((int)col_bounds.size() == 0) {
    if (t == LB) return;
    col_bound_types.resize(num_cols, LB);
    if (solver == EXLP) {
      col_bounds.resize(num_cols, 0);
    } else {
      col_bounds_double.resize(num_cols, 0);
    }
  }
  col_bound_types[c] = t;
  if (solver == EXLP) {
    col_bounds[c] = b;
  } else {
    col_bounds_double[c] = (double)b;
  }
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
    for (int i=0; i<num_cols; ++i) {
      sv[i] = soln_vector[i].get_d();
    }
  } else {
    for (int i=0; i<num_cols; ++i) {
      sv[i] = double_soln_vector[i];
    }
  }
}

void SparseLP::get_soln_vector(std::vector<long double>& sv) {
  sv.resize(num_cols);
  if (solver == EXLP) {
    for (int i=0; i<num_cols; ++i) {
      sv[i] = soln_vector[i].get_d();
    }
  } else {
    std::cout << "Num cols is " << num_cols << "\n";
    for (int i=0; i<num_cols; ++i) {
      sv[i] = (double)double_soln_vector[i];
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
    glp_iocp int_parm;
    
    lp = glp_create_prob();
    
    glp_set_prob_name(lp, "scl");
    glp_set_obj_dir(lp, GLP_MIN);
    
    glp_add_rows(lp, num_rows );
    
    glp_add_cols(lp, num_cols);
    
    for (int i=0; i<num_rows; i++) {
      glp_set_row_bnds(lp, i+1, GLP_FX, double_RHS[i], double_RHS[i]);
    }
    if ((int)col_bounds.size() == 0) {
      for (int i=0; i<num_cols; i++) {
        glp_set_col_bnds(lp, i+1, GLP_LO, 0.0, 0.0);
        glp_set_obj_coef(lp, i+1, double_objective[i]);
      }
    } else {
      for (int i=0; i<num_cols; i++) {
        switch (col_bound_types[i]) { 
          case LB:
            glp_set_col_bnds(lp, i+1, GLP_LO, col_bounds_double[i], col_bounds_double[i]);
            break;
          case UB:
            glp_set_col_bnds(lp, i+1, GLP_UP, col_bounds_double[i], col_bounds_double[i]);
            break;
          case FREE:
            glp_set_col_bnds(lp, i+1, GLP_FR, 0,0);
            break;
          case FIX:
            glp_set_col_bnds(lp, i+1, GLP_FX, col_bounds_double[i], col_bounds_double[i]);
            break;
        }
        glp_set_obj_coef(lp, i+1, double_objective[i]);
      }
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
    
    if (num_ints==0 && (solver == GLPK || solver == GLPK_SIMPLEX)) {
      glp_init_smcp(&parm);
      parm.presolve=GLP_ON;
      if (verbose > 1) {
        parm.msg_lev = GLP_MSG_ALL;
      } else {
        parm.msg_lev = GLP_MSG_OFF;
      }
      glp_simplex(lp, &parm);
      
    } else if (num_ints == 0 && solver == GLPK_IPT) {
      glp_init_iptcp(&ipt_parm);
      if (verbose > 1) {
        ipt_parm.msg_lev = GLP_MSG_ALL;
      } else {
        ipt_parm.msg_lev = GLP_MSG_OFF;
      }
      glp_interior(lp, &ipt_parm);
      
    } else if (num_ints > 0) {
      for (int i=0; i<num_cols; ++i) {
        glp_set_col_kind(lp, i+1, (col_type[i] == INT ? GLP_IV : GLP_CV));
      }
      glp_init_iocp(&int_parm);
      int_parm.msg_lev = (verbose > 1 ? GLP_MSG_ALL : GLP_MSG_OFF);
      int_parm.presolve = GLP_ON;
      glp_intopt(lp, &int_parm);
    }  
    
    int stat = 0;
    if (num_ints == 0) {
      if (solver == GLPK_IPT) {
        stat = glp_ipt_status(lp);
      } else {
        stat = glp_get_status(lp);
      }
      if (verbose > 1) {
        std::cout << "Got status " << stat << " (optimal = status " << GLP_OPT << ")\n";
      }
    } else {
      stat = glp_mip_status(lp);
      if (verbose > 1) {
        std::cout << "Got status " << stat << " (optimal = status " << GLP_OPT << ")\n";
      }
    }
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
    
    if (num_ints == 0) {
      if (verbose > 1) {
        std::cout << "Retrieving lp value...\n";
      }
      if (solver == GLPK_IPT) {
        double_op_val = glp_ipt_obj_val(lp)/4.0;	
      } else {
        double_op_val = glp_get_obj_val(lp)/4.0;	
      }
    } else {
      if (verbose > 1) {
        std::cout << "Retrieving mip value...\n";
      }
      double_op_val = glp_mip_obj_val(lp)/4.0;
    }
    
    double_soln_vector.resize(num_cols);
    if (num_ints == 0) {
      if (solver == GLPK_IPT) {
        for (int i=0; i<num_cols; i++) {
          double_soln_vector[i] = glp_ipt_col_prim(lp,i+1);
        }	
      } else {
        if (verbose>1) {
          std::cout << "Retriving solution vector from LP\n";
        }
        for (int i=0; i<num_cols; i++) {
          double_soln_vector[i] = glp_get_col_prim(lp,i+1);
        }
      }
    } else {
      for (int i=0; i<num_cols; i++) {
        double_soln_vector[i] = glp_mip_col_val(lp,i+1);
      }	
    }
    
	  glp_delete_prob(lp);
	  
  /***************************************  EXLP ****************************/  
    
	  
	} else if (solver == EXLP) {
    
    if (num_ints > 0) {
      std::cout << "Integer programming not supported with exlp\n";
      return LP_ERROR;
    }
    
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
    
    std::vector<double> lb(0);
    std::vector<double> ub(0);
    if ((int)col_bound_types.size() != 0) {
      lb.resize(num_cols);
      ub.resize(num_cols);
      for (int i=0; i<num_cols; ++i) {
        lb[i] = (col_bound_types[i]  == UB || col_bound_types[i] == FREE ? -GRB_INFINITY : col_bounds_double[i]);
        ub[i] = (col_bound_types[i]  == LB || col_bound_types[i] == FREE ? GRB_INFINITY : col_bounds_double[i]);
      }
    }
    double* lb_pointer = ((int)col_bound_types.size() == 0 ? NULL : &lb[0]);
    double* ub_pointer = ((int)col_bound_types.size() == 0 ? NULL : &ub[0]);
      
    
    //create a new model and immediately load in all the columns
    if (num_ints == 0) {  
      GRBnewmodel( env, &model, "scl", num_cols, &double_objective[0], lb_pointer, ub_double, NULL, NULL);
    } else {
      std::vector<char> var_types(num_cols);
      for (int i=0; i<num_cols; ++i) {
        var_types[i] = (col_type[i] == REAL ? GRB_CONTINUOUS : GRB_INTEGER);
      }
      GRBnewmodel( env, &model, "scl", num_cols, &double_objective[0], lb_pointer, ub_pointer, &var_types[0], NULL);
    }
    
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
    
    if (verbose > 2) {
      std::cout << "got lp solution vector: " << "\n";
      for (int i=0; i<num_cols; ++i) {
        std::cout << double_soln_vector[i] << " ";
      }
      std::cout << "\n";
    }
    
    //free stuff
    GRBfreemodel(model);
    GRBfreeenv(env);
#endif
  }
  
  return LP_OPTIMAL;
}




void SparseLP::print_LP() {
  if (solver != EXLP) {
    std::cout << "Cols: " << num_cols << "\n";
    std::cout << "Rows: " << num_rows << "\n";
    std::cout << "Objective: ";
    for (int i=0; i<num_cols; ++i) {
      std::cout << double_objective[i] << " ";
    }
    std::cout << "\n";
    std::cout << "RHS: ";
    for (int i=0; i<num_rows; ++i) {
      std::cout << double_RHS[i] << " ";
    }
    if (num_ints>0) {
      std::cout << "\nInteger: ";
      for (int i=0; i<num_cols; ++i) {
        std::cout << col_type[i] << " ";
      }
    }
    std::cout << "\nEntries: ";
    for (int i=0; i<(int)ia.size(); ++i) {
      std::cout << "(" << ia[i] << "," << ja[i] << "," << double_ar[i] << "), ";
    }
    std::cout << "\n";
  } else {
    std::cout << "Cols: " << num_cols << "\n";
    std::cout << "Rows: " << num_rows << "\n";
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










