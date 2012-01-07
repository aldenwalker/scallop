#include <vector>

#include <glpk.h>

#include "scylla_lp.h"
#include "scylla_classes.h"




extern "C" {
#include "exlp-package/lpstruct.h"
}
extern "C" {
#include "exlp-package/solve_lp.h"
}
extern "C" {
#include "exlp-package/mylib.h"
}

/***************************************************************************
 a helper function to collect duplicates
 ***************************************************************************/
void collect_dups_and_push(std::vector<int> &temp_ia,
                          std::vector<int> &temp_ja,
                          std::vector<int> &temp_ar,
                          std::vector<int> &ia,
                          std::vector<int> &ja,
                          std::vector<double> &ar,
                          int VERBOSE) {
  int j, k, temp;
  if (VERBOSE < 2) {
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
      ja.push_back(temp_ja[j]);
      ia.push_back(temp_ia[j]);
      ar.push_back(temp);
    }
  } else {
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
      ja.push_back(temp_ja[j]);
      ia.push_back(temp_ia[j]);
      ar.push_back(temp);
      std::cout << "Put " << temp_ia[j] << ", " << temp_ja[j] << ", " << temp << ".\n";
    }
  }    
}


/****************************************************************************
 * The columns are the central polygons, followed by the group teeth, 
 * followed by the group rectangles
 * 
 * the rows are the following:
 * (1) the interface edges
 * (2) the central edge *pairs*
 * (3) the group teeth restrictions
 * (4) the word constraints (N=1)
 * 
 * it's this order for simplicity later
 * ***************************************************************************/
void scylla_lp(Chain& C, 
               InterfaceEdgeList &IEL,
               CentralEdgePairList &CEL, 
               std::vector<CentralPolygon> &CP,
               std::vector<GroupTooth> &GT,
               std::vector<GroupRectangle> &GR,
               rational* scl, 
               std::vector<rational>* solution_vector, 
               scylla_lp_solver solver, 
               int VERBOSE) {
  std::vector<int> ia(0);
	std::vector<int> ja(0);
	std::vector<double> ar(0); 
  int i,j,k,m;
  int ord;
  int num_cols, offset, num_rows;
  int num_equality_rows;
  int i_edge_pairs = IEL.size();
  int c_edge_pairs = CEL.size();
  int num_words = C.num_words();
  
  //ROWS (see above)
  //we need to construct something to tell us the row 
  //number for the group teeth rows
  // group_teeth_rows_reg[i][j][k] gives the row of 
  // group i, letter j base letter k, position 1.
  // all other positions must offset from there.
  std::vector<std::vector<std::vector<int> > > group_teeth_rows_reg;
  std::vector<std::vector<std::vector<int> > > group_teeth_rows_inv;
  num_rows = i_edge_pairs + c_edge_pairs;
  group_teeth_rows_reg.resize((C.G)->num_groups());
  group_teeth_rows_inv.resize((C.G)->num_groups());
  for (i=0; i<(int)(C.G)->num_groups(); i++) {
    ord = (C.G)->index_order(i);
    if (ord == 0) {
      group_teeth_rows_reg[i].resize(0);
      group_teeth_rows_inv[i].resize(0);
      continue;
    }
    group_teeth_rows_reg[i].resize(C.regular_letters[i].size());
    for (j=0; j<(int)C.regular_letters[i].size(); j++) {
      group_teeth_rows_reg[i][j].resize(C.regular_letters[i].size());
      for (k=0; k<(int)C.regular_letters[i].size(); k++) {
        group_teeth_rows_reg[i][j][k] = num_rows;
        num_rows += ord-1;
      }
    }
    group_teeth_rows_inv[i].resize(C.inverse_letters[i].size());
    for (j=0; j<(int)C.inverse_letters[i].size(); j++) {
      group_teeth_rows_inv[i][j].resize(C.inverse_letters[i].size());
      for (k=0; k<(int)C.inverse_letters[i].size(); k++) {
        group_teeth_rows_inv[i][j][k] = num_rows;
        num_rows += ord-1;
      }
    }    
  }   
  
  num_equality_rows = num_rows;
  num_rows += num_words;
  
  num_cols = CP.size() + GT.size() + GR.size();
  
  
  if (VERBOSE>1) {
    std::cout << "Started linear programming setup\n";
  }
  
  if (solver == GLPK_SIMPLEX || solver == GLPK_IPT) {   
    
    //now we actually do the LP setup  
    //note that we need some temporary data
    std::vector<int> temp_ia;
    std::vector<int> temp_ja;
    std::vector<int> temp_ar;
    
	  glp_prob *lp;
    glp_smcp parm;
    glp_iptcp ipt_parm;
	
	  lp = glp_create_prob();
  
	  glp_set_prob_name(lp, "scl");
	  glp_set_obj_dir(lp,GLP_MIN);
    
	  glp_add_rows(lp, num_rows );
	  
    for(i=0; i<(int)num_equality_rows; i++){
		  glp_set_row_bnds(lp, i+1, GLP_FX, 0.0, 0.0);
      if (VERBOSE > 2) {
        std::cout << "Set row " << i+1 << " bounded to " << 0 << "\n";
      }
	  }
	  for(i=0; i<(int)num_words; i++){
		  glp_set_row_bnds(lp, 
                       num_equality_rows+i+1, 
                       GLP_FX, 
                       C.weights[i], 
                       C.weights[i]);	
      if (VERBOSE > 2) {
        std::cout << "Set row " << num_equality_rows+i+1 << " bounded to " << C.weights[i] << "\n";
      }
    }
    
    
    //COLS
	  glp_add_cols(lp, num_cols);
	  for(i=0; i<(int)CP.size(); i++){
		  glp_set_col_bnds(lp, i+1, GLP_LO, 0.0, 0.0);
		  glp_set_obj_coef(lp, i+1, -CP[i].chi_times_2());
      if (VERBOSE>2) {
        std::cout << "Set objective " << i+1 << " to " << -CP[i].chi_times_2() << "\n";
      }
    }
    offset = CP.size();
    for (i=0; i<(int)GT.size(); i++) {
      glp_set_col_bnds(lp, offset+i+1, GLP_LO, 0.0, 0.0);
      glp_set_obj_coef(lp, offset+i+1, -GT[i].chi_times_2(C));
      if (VERBOSE>2) {
        std::cout << "GT Set objective " << offset+i+1 << " to " << -GT[i].chi_times_2(C) << "\n";
      }
    }
    offset = CP.size() + GT.size();
    for (i=0; i<(int)GR.size(); i++) {
      glp_set_col_bnds(lp, offset+i+1, GLP_LO, 0.0, 0.0);
      glp_set_obj_coef(lp, offset+i+1, 0);
      if (VERBOSE>2) {
        std::cout << "GR Set objective " << offset+i+1 << " to " << 0 << "\n";
      }
    }    
    
	  ia.push_back(0);
	  ja.push_back(0);
	  ar.push_back(0);	
	  
    //CENTRAL POLYGONS
    for (i=0; i<(int)CP.size(); i++) {
      temp_ia.resize(0);
      temp_ja.resize(0);
      temp_ar.resize(0);
      CP[i].compute_ia_etc_for_edges(i,
                                     C,
                                     IEL, 
                                     CEL, 
                                     temp_ia,
                                     temp_ja, 
                                     temp_ar);
      if (VERBOSE > 2) {
        std::cout << "CP number " << i << ":\n";
      }
      collect_dups_and_push(temp_ia, temp_ja, temp_ar, ia, ja, ar, VERBOSE);
    }
    
    if (VERBOSE) { 
      std::cout << "Loaded central polygon edge constraints\n";
    }     
    
    //GROUP TEETH and RECTANGLES
    offset = CP.size();
    for (m=0; m<(int)GT.size(); m++) {
      if (GT[m].inverse) {
        GT[m].compute_ia_etc_for_edges(offset + m, 
                                        C, 
                                        IEL, 
                                        group_teeth_rows_inv, 
                                        ia, 
                                        ja, 
                                        ar);
      } else {
        GT[m].compute_ia_etc_for_edges(offset + m, 
                                        C, 
                                        IEL, 
                                        group_teeth_rows_reg, 
                                        ia, 
                                        ja, 
                                        ar);
      }
    }
    offset = CP.size() + GT.size();
    for (m=0; m<(int)GR.size(); m++) {
      GR[m].compute_ia_etc_for_edges(offset + m, IEL, ia, ja, ar);
    }
    
    
    if (VERBOSE>1) {
      std::cout << "Loaded group constraints\n";
    }
    
    //word constraints: for every group rectangle and group polygon, for every edge, put a 1 in the 
    //row corresponding to the word for the first letter
    offset = CP.size();
    for (j=0; j<(int)GT.size(); j++) {
      GT[j].compute_ia_etc_for_words(offset + j, C, num_equality_rows, ia, ja, ar);
    }
    offset += CP.size() + GT.size();
    for (j=0; j<(int)GR.size(); j++) {
      GR[j].compute_ia_etc_for_words(offset + j, 
                                      C, 
                                      num_equality_rows, 
                                      IEL,
                                      ia, 
                                      ja, 
                                      ar);
    }
    
    if (VERBOSE > 1) {
      std::cout << "Loaded word constraints\n";
    }
    
	  if (VERBOSE > 1) {
	    std::cout << "Created " << ia.size() << " nonzeroes on " << num_rows << " rows and " << num_cols << " columns\n";
	  }
  
	  glp_load_matrix(lp, ia.size()-1, &ia[0], &ja[0], &ar[0]);
	
    
    if (solver == GLPK_SIMPLEX) {
      glp_init_smcp(&parm);
      parm.presolve=GLP_ON;
      parm.msg_lev=GLP_MSG_ALL;
      glp_simplex(lp, &parm);
    } else if (solver == GLPK_IPT) {
      glp_init_iptcp(&ipt_parm);
      ipt_parm.msg_lev = GLP_MSG_ALL;
      glp_interior(lp, &ipt_parm);
    }
    
    if (solver == GLPK_SIMPLEX) {
      *scl = approxRat(glp_get_obj_val(lp)/4.0);	
    } else {
      *scl = approxRat_be_nice(glp_ipt_obj_val(lp)/4.0);	
    }
	
    (*solution_vector).resize(num_cols);
    if (solver == GLPK_SIMPLEX) {
      for (i=0; i<num_cols; i++) {
        (*solution_vector)[i] = approxRat(glp_get_col_prim(lp,i+1));
      }	
    } else {
      for (i=0; i<num_cols; i++) {
        (*solution_vector)[i] = approxRat_be_nice(glp_ipt_col_prim(lp,i+1));
      }	
    }
    
    if ( (VERBOSE>1 && (*solution_vector).size() < 1000) 
         || VERBOSE>2) {
      std::cout << "Solution: \n";
      for (i=0; i<(int)CP.size(); i++) {
        if ((*solution_vector)[i] == rational(0,1)) {
          continue;
        }
        std::cout << (*solution_vector)[i] << " * " << CP[i] << "\n";
      }
      offset = CP.size();
      for (j=0; j<(int)GT.size(); j++) {
        if ((*solution_vector)[offset] == rational(0,1)) {
          offset++;
          continue;
        }
        std::cout << (*solution_vector)[offset] << " * " << GT[j] << "\n";
        offset++;
      }
      for (j=0; j<(int)GR.size(); j++) {
        if ((*solution_vector)[offset] == rational(0,1)) {
          offset++;
          continue;
        }
        std::cout << (*solution_vector)[offset] << " * " << GR[j] << "\n";
        offset++;
      }
    }
  
	  glp_delete_prob(lp);
	  
	  
	  
	} else if (solver == EXLP) {
    //not implemented
  }
  
  
}

