#include <vector>
#include <iostream>
#include <stdlib.h>

#include <glpk.h>

#include "lp.h"
//#include "sssgmp.h"
#include <gmp.h>
//#include "gmp/gmp-exec/include/gmp.h"




extern "C" {
#include "exlp-package/lpstruct.h"
}
extern "C" {
#include "exlp-package/solve_lp.h"
}
extern "C" {
#include "exlp-package/mylib.h"
}



void do_linear_program(   vector<string> &w,
                          vector<int> &weight,
                          vector<arc> &arc_list, 
                          vector<polygon> &polygon_list, 
                          rational &scl, 
                          vector<rational> &solutionVector,
                          scallop_lp_solver solver,
                          int VERBOSE,
                          int LP_VERBOSE){
     


	vector<int> ia(0);
	vector<int> ja(0);
	vector<double> ar(0); 
	vector<double> s(w.size());  
  int i,j,k,l,m,n;       
  int numWords = w.size();       
  int arc_list_length = arc_list.size();
  int polygon_list_length = polygon_list.size();
    
  if (VERBOSE == 1) {
    cout << "Started linear programming setup\n";
  }
  
  
  if (solver == GLPK_DOUBLE || solver == GLPK_EXACT) {   
	  glp_prob *lp;
    glp_smcp parm;
	
	  //the maximum possible matrix size is polygon_list_length*((arc_list_length/2)+WORD)
	  ia.reserve(polygon_list_length*((arc_list_length/2)+numWords));
	  ja.reserve(polygon_list_length*((arc_list_length/2)+numWords));
	  ar.reserve(polygon_list_length*((arc_list_length/2)+numWords));
	
	  lp = glp_create_prob();
	  glp_init_smcp(&parm);
	  parm.presolve=GLP_OFF;
	  
    if (LP_VERBOSE == 1  || VERBOSE == 1) {
      parm.msg_lev=GLP_MSG_ALL;
    } else { 
      parm.msg_lev=GLP_MSG_OFF;
    }
	  glp_set_prob_name(lp, "scl");
	  glp_set_obj_dir(lp,GLP_MIN);
	
	  glp_add_rows(lp, (arc_list_length/2) + numWords );
	  for(i=1;i<=arc_list_length/2;i++){
		  glp_set_row_bnds(lp,i, GLP_FX, 0.0, 0.0);
	  };
	
	  // boundary condition
	  for(i=0;i<numWords;i++){
		  glp_set_row_bnds(lp, ((arc_list_length)/2)+1+i, GLP_FX, w[i].size()*weight[i], w[i].size()*weight[i]);	
	  };
	
	  glp_add_cols(lp, polygon_list_length);
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
	
    //cout << ia.size() << " nonzeros in " << ((arc_list_length/2) + numWords) << " rows and " << polygon_list_length << " columns\n";
	
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
	  
	  //exlp init
	  mylib_init();
	  
	  LP* lp;
    int  result;
    char buf[100];
    int varNum;

    if (VERBOSE==1) 
      cout << "About to create a new lp\n";    
    lp = new_lp(NULL);
    
    if (VERBOSE==1) 
      cout << "Done\n";
    
    if (VERBOSE)
      cout << "Init hash\n";
        
    lp_hash_str_init(lp, lp->hash_entries);
    //my_hash_mpq_init(lp->hash_entries);
    
    if (VERBOSE)
      cout << "Done\n";

    
    
    sprintf(buf, "scabble");
    lp_set_name(lp, buf);
    
    //the lp by default (set in lpstruct.c)
    //has all the right stuff, I think
    
    //now we input the rows -- it likes to name them
    for (i=0; i<(arc_list_length/2)+numWords; i++) {
      sprintf(buf, "r%d", i);
      lp_add_row(lp, buf);
      lp_set_row_equality(lp, lp_get_row_num(lp, buf), 'E');
    }
    
    if (VERBOSE==1) {
      cout << "Entered the row names and equalities\n";
    }
    
    //add the objective function row
    sprintf(buf, "obj");
    lp_set_obj_name(lp, buf);
    
    int* columnIndices = new int[polygon_list_length]; //this is probably useless
    int rowNum;
    int objectiveIndex = lp_get_row_num(lp, (char*)"obj");
    mpq_t entry;
    mpq_init(entry);
    
    for (i=0; i<polygon_list_length; i++) {
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
    
    if (VERBOSE==1) {
      cout << "Added the columns\n";
    }
    
    //this is from mps.c
    matrix_resize(lp->A, lp->rows, lp->vars);
    vector_resize(lp->b, lp->rows);
    vector_resize(lp->xb, lp->rows);
    vector_resize(lp->cb, lp->rows);
    
    
    //now we actually enter the data
    int numTimesAppears;
    for (i=0; i<polygon_list_length; i++) {
      varNum = columnIndices[i];
      
      //set the arc constraints
      for (j=0; j<(arc_list_length/2); j++) {
        sprintf(buf, "r%d", j);
        rowNum = lp_get_row_num(lp, buf);
        //put the entry in
        numTimesAppears = 0;
        for (k=0; k<polygon_list[i].size; k++) {
          if (polygon_list[i].arc[k]==2*j){
            numTimesAppears++;
          } else if (polygon_list[i].arc[k]==2*j+1){
            numTimesAppears--;
          }
        }
        mpq_set_si(entry, numTimesAppears, 1);
        lp_set_coefficient(lp, entry, rowNum, varNum);
      }
      
      //set the chain constraints
      for (j=0; j<numWords; j++) {
        sprintf(buf, "r%d", (arc_list_length/2)+j);
        rowNum = lp_get_row_num(lp, buf);
        numTimesAppears = 0;
        for (k=0; k<polygon_list[i].size; k++) {
          if (arc_list[polygon_list[i].arc[k]].first_word == j &&
              arc_list[polygon_list[i].arc[k]].first == 0) { //it's the first letter of my word
            numTimesAppears++;
          }
        }
        mpq_set_si(entry, numTimesAppears, 1);
        lp_set_coefficient(lp, entry, rowNum, varNum);
      }
      
      //set the objective function
      mpq_set_si(entry, polygon_list[i].size-2, 1);
      lp_set_coefficient(lp, entry, objectiveIndex, varNum);
    }
    
    //lp->maximize is false, so we reverse the sign
    vector_rev_sgn(lp->c);
    
    //set the right hand sides for the arcs
    mpq_set_si(entry, 0,1);
    for (i=0; i<(arc_list_length/2); i++) {
      sprintf(buf, "r%d", i);
      rowNum = lp_get_row_num(lp, buf);
      lp_set_rhs(lp, rowNum, entry);
    }
    //set it for the chain constraints
    for (i=0; i<numWords; i++) {
      sprintf(buf, "r%d", (arc_list_length/2)+i);
      rowNum = lp_get_row_num(lp, buf);
      mpq_set_si(entry, weight[i], 1);
      lp_set_rhs(lp, rowNum, entry);
    }
    
    
    //read the bounds on the columns
    mpq_set_si(entry, 0,1);
    for (i=0; i<polygon_list_length; i++) {
      varNum = columnIndices[i];
      if (lp->lower.is_valid[varNum] == FALSE) //this will always be the case I think
        mpq_init(lp->lower.bound[varNum]);
      mpq_set(lp->lower.bound[varNum], entry);
      lp->lower.is_valid[varNum] = TRUE;
    }
    
    if (VERBOSE==1) {
      cout << "Rows: " << lp->rows << "\n";;
      cout << "Vars: " << lp->vars << "\n";
    }
    
    
    result = solve_lp(lp);
      
    if (result != LP_RESULT_OPTIMAL) {
      cout << "got error code " << result << "\n";
    }
    
    lp_get_object_value(lp, &entry);
    
    scl = rational(entry)/rational(4,1);
    
    for (i=0; i<polygon_list_length; i++) {
      mpq_set(entry, *vector_get_element_ptr(lp->x, columnIndices[i]));
      solutionVector[i] = rational(entry);
    }
    
    lp_free(lp);
    
  }    
	
	
	
}
	
	
	
	
	
	
	
	
	
	
	
