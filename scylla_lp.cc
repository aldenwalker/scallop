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



/****************************************************************************
 * The columns are the central polygons, followed by the group rectangles, 
 * followed by the group polygons
 * 
 * the rows are the following:
 * (1) the interface edges
 * (2) the central edge *pairs*
 * (3) the group edge *pairs*
 * 
 * Note that some rows won't have any restrictions
 * 
 * ***************************************************************************/
void scyllop_lp(Chain &C, 
                std::vector<GroupEdgeList> &GEL, 
                InterfaceEdgeList &IEL,
                CentralEdgeList &CEL, 
                std::vector<std::vector<GroupPolygon> > &GP,
                std::vector<std::vector<GroupRectangle> > &GR,
                std::vector<CentralPolygon> &CP,
                rational* scl, 
                std::vector<rational>* solution_vector, 
                scyllop_lp_solver solver, 
                bool VERBOSE) {
  std::vector<int> ia(0);
	std::vector<int> ja(0);
	std::vector<double> ar(0); 
  int i,j,k,m,row, col, val;
  int temp_letter_1, temp_letter_2, word_1, word_2, temp;
  int num_cols, offset;
  EdgePair temp_pair;

  //this is the master list of all edge pairs
  std::vector<EdgePair> edge_pairs(IEL.edges.size());
  std::vector<int> central_edge_pairs(CEL.edges.size());
  std::vector<std::vector<int> > group_edge_pairs((C.G)->num_groups());
  
  //add all the interface edges; note these will have the SAME INDICES, THIS IS IMPORTANT
  for (i=0; i<(int)IEL.edges.size(); i++) {
    edge_pairs[i].edge_type = SCYLLA_EDGE_INTERFACE;
    edge_pairs[i].group = -1;
    edge_pairs[i].first = i;
    edge_pairs[i].second = -1;
  }
  
  //pair up the central edges
  for (i=0; i<(int)CEL.edges.size(); i++) {
    central_edge_pairs[i] = -1;
  }
  for (i=0; i<(int)CEL.edges.size(); i++) {
    if (central_edge_pairs[i] == -1) {
      temp_pair.edge_type = SCYLLA_EDGE_CENTRAL;
      temp_pair.group = -1;
      central_edge_pairs[i] = edge_pairs.size();
      temp_pair.first = i;
      temp_letter_1 = C.prev_letter( CEL[i].last );   //this pair is what we're looking for
      temp_letter_2 = C.next_letter( CEL[i].first );
      temp_pair.second = CEL.get_index( temp_letter_1, temp_letter_2 );
      central_edge_pairs[temp_pair.second] = edge_pairs.size();
      edge_pairs.push_back(temp_pair);
    }
  }
  
  //pair up all of the group edges
  for (i=0; i<(C.G)->num_groups(); i++) {
    group_edge_pairs[i].resize(GEL[i].edges.size());
    for (j=0; j<(int)GEL[i].edges.size(); j++) {
      group_edge_pairs[i][j] = -1;
    }
    for (j=0; j<(int)GEL[i].edges.size(); j++) {
      if (group_edge_pairs[i][j] == -1) {
        temp_pair.edge_type = SCYLLA_EDGE_GROUP;
        temp_pair.group = i;
        group_edge_pairs[i][j] = edge_pairs.size();
        temp_pair.first = j;
        temp_letter_1 = GEL[i][j].last;
        temp_letter_2 = GEL[i][j].first;
        temp_pair.second = GEL[i].get_index(temp_letter_1, temp_letter_2);
        group_edge_pairs[i][temp_pair.second] = edge_pairs.size();
        edge_pairs.push_back(temp_pair);
      }
    }
  }
  
  if (VERBOSE) {
    std::cout << "Edge pairs:\n";
    for (i=0; i<(int)edge_pairs.size(); i++) {
      std::cout << i << ": ";
      switch(edge_pairs[i].edge_type) {
        case SCYLLA_EDGE_CENTRAL:
          std::cout << "c ";
          break;
        case SCYLLA_EDGE_INTERFACE:
          std::cout << "i ";
          break;
        case SCYLLA_EDGE_GROUP:
          std::cout << "g" << edge_pairs[i].group << " ";
          break;
      }
      std::cout << edge_pairs[i].first << ", " << edge_pairs[i].second << "\n";
    }
  }
  

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
	  glp_add_rows(lp, edge_pairs.size() + C.num_words() );
	  for(i=0; i<(int)edge_pairs.size(); i++){
		  glp_set_row_bnds(lp, i+1, GLP_FX, 0.0, 0.0);
      //std::cout << "Set row " << i << " bounded to " << 0 << "\n";
	  }
	  for(i=0; i<(int)C.words.size(); i++){
		  glp_set_row_bnds(lp, 
                       edge_pairs.size()+i+1, 
                       GLP_FX, 
                       C.words[i].size()*C.weights[i], 
                       C.words[i].size()*C.weights[i]);	
      //std::cout << "Set row " << num_edge_pairs+i+1 << " bounded to " << words[i].size()*weights[i] << "\n";
	  }
    
    //COLS
    num_cols = CP.size();
    for (i=0; i<(C.G)->num_groups(); i++) {
      num_cols += GR[i].size();
      num_cols += GP[i].size();
    }
	  glp_add_cols(lp, num_cols);
	  for(i=0; i<(int)CP.size(); i++){
		  glp_set_col_bnds(lp, i+1, GLP_LO, 0.0, 0.0);
		  glp_set_obj_coef(lp, i+1, -CP[i].chi_times_2(C, CEL, IEL));
	  }
    offset = CP.size();
    for (i=0; i<(C.G)->num_groups(); i++) {
      for (j=0; j<(int)GR[i].size(); j++) {
        glp_set_col_bnds(lp, offset+1, GLP_LO, 0.0, 0.0);
        glp_set_obj_coef(lp, offset+1, 0);
        offset++;
      }
      for (j=0; j<(int)GP[i].size(); j++) {
        glp_set_col_bnds(lp, offset+1, GLP_LO, 0.0, 0.0);
        glp_set_obj_coef(lp, offset+1, -GP[i][j].chi_times_2(C, GEL[i], IEL));
        offset++;
      }      
    }
    if (offset != num_cols) {
      std::cout << "wrong number of columns\n";
    }
    
	  ia.push_back(0);
	  ja.push_back(0);
	  ar.push_back(0);	
	  
    //CENTRAL POLYGONS
    for (i=0; i<(int)CP.size(); i++) {
      temp_ia.resize(0);
      temp_ja.resize(0);
      temp_ar.resize(0);
      col = i;
      
      for (j=0; j<(int)CP[i].edges.size(); j++) {
        if (CP[i].interface[j]) {
          if (C.next_letter( IEL[CP[i].edges[j]].last ) == IEL[CP[i].edges[j]].first ) {         //don't restrict edges like this
            continue;
          }
          row = CP[i].edges[i];
          val = 1;
        } else {
          if (C.next_letter( CEL[CP[i].edges[j]].first ) == CEL[CP[i].edges[j]].last) {   //nor these edges
            continue;
          }
          row = central_edge_pairs[CP[i].edges[j]];
          if (edge_pairs[row].first == CP[i].edges[j]) { 
            val = 1;
          } else {
            val = -1;
          }
        } 
        temp_ja.push_back(col+1);
        temp_ia.push_back(row+1);
        temp_ar.push_back(val);
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
      std::cout << "Loaded central polygon edge constraints\n";
    }     
    
    offset = CP.size();
    
    //GROUP RECTANGLES AND POLYS
    for (i=0; i<(C.G)->num_groups(); i++) {
      for (j=0; j<(int)GR[i].size(); j++) {
        col = offset;
        row = GR[i][j].first;
        ia.push_back(row+1);
        ja.push_back(col+1);
        ar.push_back(-1);
        
        row = GR[i][j].last;
        ia.push_back(row+1);
        ja.push_back(col+1);
        ar.push_back(-1);
        offset++;
      }
      
      for (m=0; m<(int)GP[i].size(); m++) {
        temp_ia.resize(0);
        temp_ja.resize(0);
        temp_ar.resize(0);
        GP[i][m].get_ia_etc_for_edges(C, IEL, GEL[i], edge_pairs, group_edge_pairs[i], offset, temp_ia, temp_ja, temp_ar);
        offset++;
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
    }
    
    if (VERBOSE) {
      std::cout << "Loaded group constraints\n";
    }
    
    //word constraints: for every group rectangle and group polygon, for every edge, put a 1 in the 
    //row corresponding to the word for the first letter
    offset = CP.size();
    for (i=0; i<(C.G)->num_groups(); i++) {
      for (j=0; j<(int)GR[i].size(); j++) {
        col = offset;
        word_1 = C.chain_letters[IEL.edges[GR[i][j].first].first].word;
        word_2 = C.chain_letters[IEL.edges[GR[i][j].last].first].word;
        if (word_1 == word_2) {
          ia.push_back(word_1 + edge_pairs.size() + 1);
          ja.push_back(offset + 1);
          ar.push_back(2);
        } else {
          ia.push_back(word_1 + edge_pairs.size() + 1);
          ja.push_back(offset + 1);
          ar.push_back(1);
          ia.push_back(word_2 + edge_pairs.size() + 1);
          ja.push_back(offset + 1);
          ar.push_back(1);
        }
        offset++;
      }
      for (m=0; m<(int)GP[i].size(); m++) {
        temp_ia.resize(0);
        temp_ja.resize(0);
        temp_ar.resize(0);
        GP[i][m].get_ia_etc_for_words(C, IEL, GEL[i], edge_pairs, group_edge_pairs[i], offset, temp_ia, temp_ja, temp_ar);
        offset++;
        //remove duplicates
        for (j=0; j<(int)temp_ia.size(); j++) {
          if (temp_ar[j] == 0 ) {
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
    }
    
    if (VERBOSE) {
      std::cout << "Loaded word constrains\n";
    }
    
	  if (VERBOSE) {
	    std::cout << "Created " << ia.size() << " constraints\n";
	  }
  
	  glp_load_matrix(lp, ia.size()-1, &ia[0], &ja[0], &ar[0]);
	
    glp_simplex(lp, &parm);
	
    *scl = approxRat(glp_get_obj_val(lp)/4.0);	
	
    (*solution_vector).resize(num_cols);
	  for (i=0; i<num_cols; i++) {
	    (*solution_vector)[i] = approxRat(glp_get_col_prim(lp,i+1));
	  }	
    
    if (VERBOSE) {
      std::cout << "Solution: \n";
      for (i=0; i<(int)CP.size(); i++) {
        if ((*solution_vector)[i] == rational(0,1)) {
          continue;
        }
        std::cout << (*solution_vector)[i] << " * " << CP[i];
      }
      offset = 0;
      for (i=0; i<(C.G)->num_groups(); i++) {
        for (j=0; j<(int)GR[i].size(); j++) {
          if ((*solution_vector)[i] == rational(0,1)) {
            offset++;
            continue;
          }
          std::cout << (*solution_vector)[offset] << " * " << GR[i][j];
          offset++;
        }
        for (j=0; j<(int)GP[i].size(); j++) {
          if ((*solution_vector)[i] == rational(0,1)) {
            offset++;
            continue;
          } 
          std::cout << (*solution_vector)[offset] << " * " << GP[i][j];
          offset++;
        }
      }
    }
  
	  glp_delete_prob(lp);
	  
	  
	  
	} else if (solver == EXLP) {
    //not implemented
  }
  
  
}

