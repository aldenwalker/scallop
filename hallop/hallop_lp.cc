#include <vector>
#include <iostream>

#include "hallop_lp.h"


void HALLOP::hallop_lp(HALLOP::FreeGroupChain& C,
                       HALLOP::Pieces& P, 
                       SparseLPSolver solver,
                       Rational& scl, 
                       std::vector<Rational>& soln_vec, 
                       int verbose) {
  
  int Nrects = P.rects.size();
  int Ntris = P.tris.size();
  
  //rows and columns in order:
  
  //there is one row for every GluingEdge
  //there is one row for every word in the chain (setting it equal to its weight)
  //there is one row for every relator (setting it equal to the column which records it)
  
  //there is a column for every rectangle
  //there is a column for every triangle
  //there is a column for every relator, recording how often it is used
  
  //it's assumed that the chain C is the usual chain, followed by relators with weight -1
  int relators_start_word = -1;
  int num_relators = 0;
  for (int i=0; i<C.num_words(); ++i) {
    if (C.weights[i] == -1) { 
      relators_start_word = i;
      num_relators = C.num_words() - relators_start_word;
      break;
    }
  }
  
  int relators_start_column = Nrects + Ntris;
  int words_start_row = P.edges.size();
  
  int NR = words_start_row + C.num_words();
  int NC = relators_start_column + num_relators;
  
  SparseLP LP(solver, NR, NC);
  
  /*********** matrix entries ******************/
  
  //add the entries for all the rectangles
  for (int i=0; i<Nrects; ++i) {
    int ind, s;
    
    //these are the gluing edge conditions
    extract_signed_index(P.rects[i].b1, ind, s);
    LP.add_entry(ind, i, s);
    extract_signed_index(P.rects[i].b2, ind, s);
    LP.add_entry(ind, i, s);
    
    //if the rectangle contains first-letters of words, we need to record that
    if (C.chain_letters[P.rects[i].let1].index == 0) {
      LP.add_entry(words_start_row + C.chain_letters[P.rects[i].let1].word, i, 1);
    }
    if (C.chain_letters[P.rects[i].let2].index == 0) {
      LP.add_entry(words_start_row + C.chain_letters[P.rects[i].let2].word, i, 1);
    }
  }
  
  //add entries for the triangles
  for (int i=0; i<Ntris; ++i) {
    int ind, s;
    //these are the gluing edge conditions
    extract_signed_index(P.tris[i].b1, ind, s);
    LP.add_entry(ind, Nrects + i, s);
    extract_signed_index(P.tris[i].b2, ind, s);
    LP.add_entry(ind, Nrects + i, s);
    extract_signed_index(P.tris[i].b3, ind, s);
    LP.add_entry(ind, Nrects + i, s);
  }
  
  //add the entries for the relator columns
  for (int i=0; i<num_relators; ++i) {
    LP.add_entry(words_start_row + relators_start_word + i, relators_start_column + i, -1);
  }
  
  /*************  RHS **************************/
  
  //every GluingEdge row has a RHS of zero
  for (int i=0; i<(int)P.edges.size(); ++i) {
    LP.set_RHS(i, 0);
  }
  //every normal word row has a RHS of the weight
  for (int i=0; i<relators_start_word; ++i) {
    LP.set_RHS(P.edges.size() + i, C.weights[i]);
  }
  //every relator word has a RHS of zero (since it's a free variable >=0)
  for (int i=relators_start_word; i<C.num_words(); ++i) {
    LP.set_RHS(P.edges.size() + i, 0);
  }
  
  /************** objective *********************/
  
  //the objective computes 2(-chi-n), where n is the number of relators
  //so, every triangle contributes 1 
  //and every relator contributes -2
  for (int i=0; i<Nrects; ++i) {
    LP.set_obj(i, 0);
  }
  for (int i=0; i<Ntris; ++i) {
    LP.set_obj(Nrects + i, 1);
  }
  for (int i=0; i<num_relators; ++i) {
    LP.set_obj(Nrects + Ntris + i, -2);
  }

  /************** solving ***********************/

  LP.solve(verbose-1);

  LP.get_optimal_value(scl);
  LP.get_soln_vector(soln_vec);
  
} 

  
  
  
  
  
  
  
  
  
  
  
  
  
