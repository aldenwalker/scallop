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



void scyllop_lp(CyclicProduct& G, 
                Chain& C, 
                std::vector<std::vector<Multiarc> > &arcs, 
                std::vector<Edge> &edges,
                std::vector<Polygon> &polys, 
                rational* scl, 
                std::vector<rational>* solution_vector, 
                scyllop_lp_solver solver, 
                int VERBOSE) {
  std::vector<int> ia(0);
	std::vector<int> ja(0);
	std::vector<double> ar(0); 
  int i,j,k,l,m,n;          
  
  int num_edge_pairs;
  int num_polys;
  int num_multiarcs;
  
  std::vector<std::string> words = C.word_list();
  std::vector<int> weghts = C.weights_list();
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
                     real_edges_beginning_with[temp1][k]
                    ].last == temp2 ) {
            break;
          }
          temp_arc.push_back( real_edges_beginning_with[temp1][k] );
        }
      }
      stripped_arcs.push_back(temp_arc);
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
    if (which_edge_pair[i] != -1) {
      continue;
    }
    if (edges[i].blank) {
      temp1 = C.prev_letter(edges[i].last);
      temp2 = C.next_letter(edges[i].first);
      for (j=0; j<(int)blank_edges_beginning_with[temp1].size(); j++) {
        if (edges[ blank_edges_beginning_with[temp1][k] ].last == temp2) {
          break;
        }
      }
      paired_edge = blank_edges_beginning_with[temp1][k];
    
    } else {
      temp1 = edges[i].last;
      temp2 = edges[i].first;
      for (j=0; j<(int)real_edges_beginning_with[temp1].size(); j++) {
        if (edges[ real_edges_beginning_with[temp1][k] ].last == temp2) {
          break;
        }
      }
      paired_edge = real_edges_beginning_with[temp1][k];
    }
    //it will always be the case that the paired edge is larger than i
    temp_pair.first = i;
    temp_pair.second = paired_edge;
    which_edge_pair[i] = edge_pairs.size();
    which_edge_pair[paired_edge] = edge_pairs.size();
    edge_pairs.push_back(temp_pair);
  }
  
  
  num_edge_pairs = edge_pairs.size();
  num_polys = polys.size();
  num_multiarcs = stripped_arcs.size();
    
  if (VERBOSE == 1) {
    std::cout << "Started linear programming setup\n";
  }
  
  
  if (solver == GLPK_DOUBLE || solver == GLPK_EXACT) {   
	  glp_prob *lp;
    glp_smcp parm;
	
	  lp = glp_create_prob();
	  glp_init_smcp(&parm);
	  parm.presolve=GLP_OFF;
	  
	  parm.msg_lev=GLP_MSG_OFF;
	  glp_set_prob_name(lp, "scl");
	  glp_set_obj_dir(lp,GLP_MIN);
	
	  glp_add_rows(lp, num_edge_pairs + numWords );
	  for(i=1; i<=num_edge_pairs; i++){
		  glp_set_row_bnds(lp,i, GLP_FX, 0.0, 0.0);
	  }
	
	  // boundary condition
	  for(i=0; i<(int)words.size(); i++){
		  glp_set_row_bnds(lp, 
                       num_edge_pairs+i+1, 
                       GLP_FX, 
                       words[i].size()*weights[i], 
                       words[i].size()*weights[i]);	
	  }
	
	  glp_add_cols(lp, num_polys + num_multiarcs);
    //START HERE
	  for(i=1;i<=polygon_list_length;i++){
		  glp_set_col_bnds(lp,i, GLP_LO, 0.0, 0.0);
		  glp_set_obj_coef(lp,i, (polygon_list[i-1].size-2));
	  };
	  l=0;
	  //start the vectors with the dummy zeroth element
	  ia.push_back(0);
	  ja.push_back(0);
	  ar.push_back(0);	
	  for(j=0;j<polygon_list_length;j++){
		  for(m=0;m<arc_list_length/2;m++){
			  n=0;
			  for(k=0;k<polygon_list[j].size;k++){
				  if(polygon_list[j].arc[k]==2*m){
					  n=n+1;
				  }
				  if(polygon_list[j].arc[k]==2*m+1){
					  n=n-1;
				  }
			  }
			  if(n!=0){
				  l++;
				  //ia[l]=m+1;
				  //ja[l]=j+1;
				  //ar[l]=n;
				  ia.push_back(m+1);
				  ja.push_back(j+1);
				  ar.push_back(n);				
		  //		cout << "entry " << l << " i = " << m+1 << " j = " << j+1 << " entry = " << n << "\n";
			  }
		  }
		
		  for(i=0;i<numWords;i++){
			  s[i]=0;
		  }
		
		  for(k=0;k<polygon_list[j].size;k++){
			  s[arc_list[polygon_list[j].arc[k]].first_word]++;
		  }
		
		  for(i=0;i<numWords;i++){
		    l++;
		    //ia[l]=(arc_list_length/2)+1+i;
		    //ja[l]=j+1;
		    //ar[l]=s[i];
		    ia.push_back((arc_list_length/2)+1+i);
		    ja.push_back(j+1);
		    ar.push_back(s[i]);
		  };

	  };
	
	  if (VERBOSE == 1) {
	    cout << "Created constraints\n";
	  }
	
    //cout << l << " matrix entries \n";
	  //cout <<"vector sizes:" << ia.size() << " " << ja.size() << " " << ar.size() << "\n";
	
	  //the contents of a vector<> are guaranteed to be contiguous in memory, 
	  //so this should be ok
	  glp_load_matrix(lp, ia.size()-1, &(ia[0]), &(ja[0]), &(ar[0]));
	
   
   
    glp_simplex(lp,&parm);
   	if (solver == GLPK_EXACT) {
   //////////////////////  this line might not compile with older glpk ////////
    //glp_exact(lp, &parm);
    ///////////////////////////////////////////////////////////////////////////
	  }
	
	
	
  //	lpx_exact(lp);
    scl = approxRat(glp_get_obj_val(lp)/4.0);	
	
	  for (i=0; i<polygon_list_length; i++) {
	    solutionVector[i] = approxRat(glp_get_col_prim(lp,i+1));
	  }	
	
	  glp_delete_prob(lp);
	  
	  
	  
	} else if (solver == EXLP) {
    //not implemented
  }
  */
  
}
