#include <vector>

#include <glpk.h>

#include "scyllop_lp.h"
#include "scyllop_classes.h"




extern "C" {
#include "exlp-package/lpstruct.h"
}
extern "C" {
#include "exlp-package/solve_lp.h"
}
extern "C" {
#include "exlp-package/mylib.h"
}



/****************************************************************************
 * The columns are the polygons, followed by the multiarcs 
 * the rows are the edge constraints, followed by the word (chain) constraints,
 * which are accomplished by putting a 1 in the row corresponding to i for every time a polygon
 * (NOT a multiarc, to avoid dups) has the edge (i,j)
 * 
 * The edge constraints are as follows: for every polygon, every time we 
 * see a real edge, we put a 1 down.  Every time we see a blank edge, 
 * we put a 1 or -1 depending on which edge that is in the pair
 * 
 * For the multiarcs, every time we see an edge, we put a -1 down *for its pair*
 * 
 * This conplication is because we glue polys to multiarcs along only real edges
 * and only polys to polys along blank edges
 * 
 * there is one constraint for each real edge (one for each in the pair),
 * and one constraint for each blank edge pair
 * 
 * ***************************************************************************/
void scyllop_lp(CyclicProduct& G, 
                Chain& C, 
                std::vector<std::vector<Multiarc> > &arcs, 
                std::vector<Edge> &edges,
                std::vector<Polygon> &polys, 
                rational* scl, 
                std::vector<rational>* solution_vector, 
                scyllop_lp_solver solver, 
                bool VERBOSE) {
  std::vector<int> ia(0);
	std::vector<int> ja(0);
	std::vector<double> ar(0); 
  int i,j,k,l,row, col, temp;          
  
  int num_edge_pairs;
  int num_polys;
  int num_multiarcs;
  
  std::vector<std::string> words = C.word_list();
  std::vector<int> weights = C.weights_list();
  std::vector<ChainLetter> chain_letters = C.chain_letter_list();
  std::vector<std::vector<int> > real_edges_beginning_with(chain_letters.size());
  std::vector<std::vector<int> > blank_edges_beginning_with(chain_letters.size());
    
  for (i=0; i<(int)chain_letters.size(); i++) {
    real_edges_beginning_with[i].resize(0);
    blank_edges_beginning_with[i].resize(0);
  }
  
  for (i=0; i<(int)edges.size(); i++) {
    if (edges[i].blank) {
      blank_edges_beginning_with[edges[i].first].push_back(i);
    } else {
      real_edges_beginning_with[edges[i].first].push_back(i);
    }
  }
    
  //first, make the stripped arcs
  std::vector<std::vector<int> > stripped_arcs(0);
  std::vector<int> temp_arc;
  int temp1;
  int temp2;
  for (i=0; i<(int)arcs.size(); i++) {
    for (j=0; j<(int)arcs[i].size(); j++) {
      temp_arc.resize(0);
      for (k=0; k<(int)arcs[i][j].letters.size(); k++) {
        temp1 = arcs[i][j].letters[k];
        temp2 = arcs[i][j].letters[(k+1)%(int)arcs[i][j].letters.size()];
        for (l=0; l<(int)real_edges_beginning_with[temp1].size(); l++) {
          if ( edges[
                     real_edges_beginning_with[temp1][l]
                    ].last == temp2 ) {
            break;
          }
        }
        temp_arc.push_back( real_edges_beginning_with[temp1][l] );
      }        
      stripped_arcs.push_back(temp_arc);
      std::cout << "arc: ";
      for (l=0; l<(int)temp_arc.size(); l++) {
        std::cout << temp_arc[l] << ",";
      }
      std::cout << "\n";
    }
  }
  
  
  
  //now, we need to make the edge pairs
  //note that for real edges, the edge (i,j) is paired with (j,i)
  //for fake edges, (i,j) is paired with (j-1, i+1)
  std::vector<std::pair<int, int> > edge_pairs(0);
  std::pair<int, int> temp_pair;
  std::vector<int> which_edge_pair(edges.size(),-1);
  int paired_edge;
  for (i=0; i<(int)edges.size(); i++) {
    if (edges[i].blank && which_edge_pair[i] != -1) {
      continue;
    }
    if (edges[i].blank) {
      temp1 = C.prev_letter(edges[i].last);
      temp2 = C.next_letter(edges[i].first);
      for (j=0; j<(int)blank_edges_beginning_with[temp1].size(); j++) {
        if (edges[ blank_edges_beginning_with[temp1][j] ].last == temp2) {
          break;
        }
      }
      paired_edge = blank_edges_beginning_with[temp1][j];
    
    } else {
      temp1 = edges[i].last;
      temp2 = edges[i].first;
      for (j=0; j<(int)real_edges_beginning_with[temp1].size(); j++) {
        if (edges[ real_edges_beginning_with[temp1][j] ].last == temp2) {
          break;
        }
      }
      paired_edge = real_edges_beginning_with[temp1][j];
    }
    //it will always be the case that the paired edge is larger than i
    temp_pair.first = i;
    temp_pair.second = paired_edge;
    which_edge_pair[i] = edge_pairs.size();
    if (edges[i].blank) {
      which_edge_pair[paired_edge] = edge_pairs.size();
    }
    edge_pairs.push_back(temp_pair);
  }
  
  if (VERBOSE) {
    std::cout << "Edge pairs:\n";
    for (i=0; i<(int)edge_pairs.size(); i++) {
      std::cout << i << ": " << edge_pairs[i].first << ", " << edge_pairs[i].second << "\n";
    }
  }
  
  num_edge_pairs = edge_pairs.size();
  num_polys = polys.size();
  num_multiarcs = stripped_arcs.size();
    
  if (VERBOSE) {
    std::cout << "Started linear programming setup\n";
  }
  
  
  if (solver == GLPK_DOUBLE || solver == GLPK_EXACT) {   
    
    //now we actually do the LP setup  
    //note that we need some temporary data
    std::vector<int> temp_ia;
    std::vector<int> temp_ja;
    std::vector<int> temp_ar;
    
	  glp_prob *lp;
    glp_smcp parm;
	
	  lp = glp_create_prob();
	  glp_init_smcp(&parm);
	  parm.presolve=GLP_OFF;
	  
	  parm.msg_lev=GLP_MSG_ALL;
	  glp_set_prob_name(lp, "scl");
	  glp_set_obj_dir(lp,GLP_MIN);
	
    //ROWS
	  glp_add_rows(lp, num_edge_pairs + words.size() );
	  for(i=1; i<=num_edge_pairs; i++){
		  glp_set_row_bnds(lp,i, GLP_FX, 0.0, 0.0);
      //std::cout << "Set row " << i << " bounded to " << 0 << "\n";
	  }
	  for(i=0; i<(int)words.size(); i++){
		  glp_set_row_bnds(lp, 
                       num_edge_pairs+i+1, 
                       GLP_FX, 
                       words[i].size()*weights[i], 
                       words[i].size()*weights[i]);	
      //std::cout << "Set row " << num_edge_pairs+i+1 << " bounded to " << words[i].size()*weights[i] << "\n";
	  }
    
    //COLS
	  glp_add_cols(lp, num_polys + num_multiarcs);
	  for(i=0; i<num_polys; i++){
		  glp_set_col_bnds(lp, i+1, GLP_LO, 0.0, 0.0);
		  glp_set_obj_coef(lp, i+1, polys[i].edges.size()-2);
	  }
    for (i=0; i<num_multiarcs; i++) {
      glp_set_col_bnds(lp, num_polys+i+1, GLP_LO, 0.0, 0.0);
      glp_set_obj_coef(lp, num_polys+i+1, stripped_arcs[i].size()-2);
    }
    
	  ia.push_back(0);
	  ja.push_back(0);
	  ar.push_back(0);	
	  
    //polygon constraints from edges
    for (i=0; i<(int)polys.size(); i++) {
      col = i;
      temp_ia.resize(0);
      temp_ja.resize(0);
      temp_ar.resize(0);
      for (j=0; j<(int)polys[i].edges.size(); j++) {
        row = which_edge_pair[polys[i].edges[j]];     
        temp_ja.push_back(col+1);
        temp_ia.push_back(row+1);
        if (edges[ polys[i].edges[j] ].blank) {
          //blank edge
          if (edge_pairs[row].first == polys[i].edges[j]) {
            temp_ar.push_back(1);
          } else {
            temp_ar.push_back(-1);
          }
        } else {
          //real edge
          temp_ar.push_back(1);
        }
      }
      //now actually insert, after removing dups
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
        //std::cout << "Put " << temp_ia[j] << ", " << temp_ja[j] << ", " << temp << ".\n";
      }
    }
    
    
    if (VERBOSE) { 
      std::cout << "Loaded polygon edge constraints\n";
    }     
    
    //multiarc constraints from edges
    for (i=0; i<(int)stripped_arcs.size(); i++) {
      col = num_polys + i;
      temp_ia.resize(0);
      temp_ja.resize(0);
      temp_ar.resize(0);
      for (j=0; j<(int)stripped_arcs[i].size(); j++) {
        temp = stripped_arcs[i][j];     //this is the edge on the multiarc
        row = which_edge_pair[edge_pairs[which_edge_pair[temp]].second];    //this is the row (we always put a -1), but we need to row for the *pair*
        temp_ja.push_back(col+1);
        temp_ia.push_back(row+1);
        temp_ar.push_back(-1);
      }
      //remove duplicates
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
        //std::cout << "Put " << temp_ia[j] << ", " << temp_ja[j] << ", " << temp << ".\n";
      }
    }
    
    if (VERBOSE) {
      std::cout << "Loaded multiarc edge constraints\n";
    }
    
    //word constraints: for every multiarc, for every edge, put a 1 in the 
    //row corresponding to the word for the first letter
    for (i=0; i<(int)stripped_arcs.size(); i++) {
      col = num_polys + i;
      temp_ia.resize(0);
      temp_ja.resize(0);
      temp_ar.resize(0);
      for (j=0; j<(int)stripped_arcs[i].size(); j++) {
        temp = stripped_arcs[i][j];
        row = chain_letters[
                            edges[
                                  stripped_arcs[i][j]
                                 ].first
                            ].word;
        row += num_edge_pairs;
        temp_ja.push_back(col+1);
        temp_ia.push_back(row+1);
        temp_ar.push_back(1);
      }
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
        //std::cout << "Put " << temp_ia[j] << ", " << temp_ja[j] << ", " << temp << ".\n";
      }
    }
    
    if (VERBOSE) {
      std::cout << "Loaded word constrains\n";
    }
    
	  if (VERBOSE) {
	    std::cout << "Created " << ia.size() << " constraints\n";
	  }
    
    /* //print the matrix
    for (i=0; i<num_edge_pairs+(int)words.size(); i++) {
      for (j=0; j<num_polys+num_multiarcs; j++) {
        for (k=0; k<(int)ia.size(); k++) {
          if (ia[k] == i+1 && ja[k] == j+1) {
            std::cout << ar[k] << "\t";
            break;
          }
        }
        if (k==(int)ia.size()) {
          std::cout << "0\t";
        }
      }
      std::cout << "\n";
    }
    */
          
  
	  glp_load_matrix(lp, ia.size()-1, &ia[0], &ja[0], &ar[0]);
	
    glp_simplex(lp, &parm);
	
    *scl = approxRat(glp_get_obj_val(lp)/4.0);	
	
    (*solution_vector).resize(num_polys + num_multiarcs);
	  for (i=0; i<num_polys + num_multiarcs; i++) {
	    (*solution_vector)[i] = approxRat(glp_get_col_prim(lp,i+1));
	  }	
    
    if (VERBOSE) {
      std::cout << "Solution: \n";
      for (i=0; i<(int)(*solution_vector).size(); i++) {
        if ((*solution_vector)[i] == rational(0,1)) {
          continue;
        }
        std::cout << (*solution_vector)[i] << " * ";
        std::cout.flush();
        if (i < num_polys) {
          for (j=0; j<(int)polys[i].edges.size(); j++) {
            std::cout << polys[i].edges[j];
            if (edges[polys[i].edges[j]].blank) {
              std::cout << "!";
            }
            std::cout << " ";
          }
          std::cout << "\n";
        } else {
          std::cout << "ma: ";
          for (j=0; j<(int)stripped_arcs[i-num_polys].size(); j++) {
            std::cout << stripped_arcs[i-num_polys][j] << " ";
          }
          std::cout << "\n";
        }
      }
    }
  
	  glp_delete_prob(lp);
	  
	  
	  
	} else if (solver == EXLP) {
    //not implemented
  }
  
  
}
