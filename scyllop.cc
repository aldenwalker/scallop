/****************************************************************************
* compute scl in free products of cyclic groups (finite and infinite)
* By Alden Walker
* Implements a generalization of the algorithm described in "scl"
*****************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <ctype.h>

#include "scyllop_classes.h"
#include "rational.h"
#include "scyllop_lp.h"

/*****************************************************************************
 * make the list of multiarcs -- each multiarc is for a different group, and 
 * each one is made up of a collection (possibly including duplicates) of 
 * chunks, where the total number of letters is the order of the group.  
 * only the cyclic order in the multiarc matters.  
 *
 * The current thinking is that the arcs come in three kinds: positive, inverse, 
 * and matched pairs.  Also, wlog we may make all the multiarcs have the maximal
 * number of sides (no fussing with chunks)
 * ***************************************************************************/
void compute_multiarcs(CyclicProduct &G, Chain &C, std::vector<std::vector<Multiarc> > &arcs) {
  std::vector<char> gens = G.gen_list();
  std::vector<int> orders = G.order_list();
  int i,j,k;
  int num_groups = gens.size();
  std::vector<std::vector<int> > group_letters = C.group_letter_list();
  std::vector<ChainLetter> chain_letters = C.chain_letter_list();
  std::vector<int> current_possible_letters;
  Multiset letter_selection;
  Multiarc temp_marc;
  int first_index;
  
  arcs.resize(num_groups);
  
  for (i=0; i<num_groups; i++) {
    temp_marc.group = i;
    arcs[i].resize(0);
    //first, the easy ones -- the pairs of inverses
    for (j=0; j<(int)group_letters[i].size(); j++) {
      for (k=j+1; k<(int)group_letters[i].size(); k++) {
        if (chain_letters[group_letters[i][j]].letter != 
            chain_letters[group_letters[i][k]].letter) {    //they're the same group, so inverse if not the same letter
          temp_marc.letters.resize(2);
          temp_marc.letters[0] = group_letters[i][j];
          temp_marc.letters[1] = group_letters[i][k];
          arcs[i].push_back(temp_marc);
        }
      }
    }
    
    if (orders[i] == 0) {
      continue;
    }
    //now the harder ones -- all groups of exactly the order of the group, 
    //all of the same type (inverse or not)
    //first, the regulars
    int inverses;
    for (inverses=0; inverses<2; inverses++) {
      current_possible_letters.resize(0);
      temp_marc.letters.resize(orders[i]);
      for (j=0; j<(int)group_letters[i].size(); j++) {
        if ( (!inverses && islower(chain_letters[group_letters[i][j]].letter))
            || (inverses && isupper(chain_letters[group_letters[i][j]].letter)) ) {
          current_possible_letters.push_back(group_letters[i][j]);
        }
      }
      if ((int)current_possible_letters.size() == 0) {
        continue;
      }
      for (first_index=0; first_index<(int)group_letters[i].size(); first_index++) {
        letter_selection = Multiset(orders[i]-1, first_index, current_possible_letters.size());
        temp_marc.letters[0] = current_possible_letters[first_index];
        do {
          for (j=0; j<orders[i]-1; j++) {
            temp_marc.letters[j+1] = current_possible_letters[letter_selection[j]];
          }
          arcs[i].push_back(temp_marc);
        } while (1 != letter_selection.next());
      }
    }
    
  }//end of loop through the groups

} 



/*****************************************************************************
 * compute the edges.  there are two kinds -- real edges, which must come from
 * multiarcs, and are given by any pair that can appear there, and blank
 * arcs.  These are given by *all* pairs of letters.
 *****************************************************************************/
void compute_edges(CyclicProduct &G, Chain &C, std::vector<Edge> &edges) {
  int i,j;
  std::vector<ChainLetter> chain_letters = C.chain_letter_list();
  std::vector<int> orders = G.order_list();
  int num_letters = (int)chain_letters.size();
  edges.resize(0);
  Edge temp_edge;
  for (i=0; i<num_letters; i++) {
    temp_edge.first = i;
    for (j=0; j<num_letters; j++) {
      temp_edge.last = j;
      temp_edge.blank = true;
      edges.push_back(temp_edge);
      if (chain_letters[i].group == chain_letters[j].group) {
        temp_edge.blank = false;
        if (orders[chain_letters[i].group] == 0) {
          if (chain_letters[i].letter != chain_letters[j].letter) { //must be inverses
            edges.push_back(temp_edge);
          }
        } else {
          edges.push_back(temp_edge);
        }
      }
    }
  }
}






/*****************************************************************************
 * compute the list of polys.  these come in two types: those with blank
 * edges, and those without
 *****************************************************************************/
void compute_polys(CyclicProduct &G, 
                   Chain &C, 
                   std::vector<Edge> &edges, 
                   std::vector<Polygon> &polys) {
  int i,j,k;
  int num_edges = (int)edges.size();
  int num_real_edges;
  int num_blank_edges;
  int temp1, temp2;
  Polygon temp_poly;
  std::vector<int> word_lens(C.num_words());
  std::vector<int> real_edges(0);
  std::vector<int> blank_edges(0);
  std::vector<ChainLetter> chain_letters = C.chain_letter_list();
  std::vector<std::vector<int> > real_edges_beginning_with(chain_letters.size());
  std::vector<std::vector<int> > blank_edges_beginning_with(chain_letters.size());
  
  for (i=0; i<C.num_words(); i++) {
    word_lens[i] = C[i].size();
  }
  
  for (i=0; i<(int)chain_letters.size(); i++) {
    real_edges_beginning_with[i].resize(0);
    blank_edges_beginning_with[i].resize(0);
  }
  
  for (i=0; i<num_edges; i++) {
    if (edges[i].blank) {
      blank_edges.push_back(i);
      blank_edges_beginning_with[edges[i].first].push_back(i);
    } else {
      real_edges.push_back(i);
      real_edges_beginning_with[edges[i].first].push_back(i);
    }
  }
  num_real_edges = real_edges.size();
  num_blank_edges = blank_edges.size();
  
  polys.resize(0);
  
  //the ones with two blank edges are easy, since we simply pick any two 
  //real edges, and that's it.  
  temp_poly.edges.resize(4);
  for (i=0; i<num_real_edges; i++) {
    temp_poly.edges[0] = real_edges[i];
    for (j=0; j<num_real_edges; j++) {
      temp_poly.edges[2] = real_edges[j];
      //find the blank edges to fill in
      temp1 = edges[temp_poly.edges[0]].last;
      temp2 = edges[temp_poly.edges[2]].first;
      for (k=0; k<(int)blank_edges_beginning_with[temp1].size(); k++) {
        if (edges[ blank_edges_beginning_with[temp1][k] ].last == temp2) {
          temp_poly.edges[1] = blank_edges_beginning_with[temp1][k];
          break;
        }
      }
      temp1 = edges[temp_poly.edges[2]].last;
      temp2 = edges[temp_poly.edges[0]].first;
      for (k=0; k<(int)blank_edges_beginning_with[temp1].size(); k++) {
        if (edges[ blank_edges_beginning_with[temp1][k] ].last == temp2) {
          temp_poly.edges[3] = blank_edges_beginning_with[temp1][k];
          break;
        }
      }
      polys.push_back(temp_poly);
    }
  }
  
  std::cout << "I am done with the 2-blank polys\n";
  for (i=0; i<(int)polys.size(); i++) {
    for (j=0; j<(int)polys[i].edges.size(); j++) {
      std::cout << polys[i].edges[j] << ",";
    }
    std::cout << "\n";
  }
      
  //now the ones with 1 or no blank edges.  Here we just go through all 
  //possibilities, using only real edges.  Note we may assume that 
  //the polygon starts on a minimal letter
  int num_two_blank_polys = polys.size();
  std::vector<int> current_beginning_letters(4);    //this records the first letters of the arc choices
  std::vector<int> current_edges(4);                //this records where we are in the lists real_edges_beginning_with
  int current_len;                                  //this records the current length
  current_len = 1;
  current_beginning_letters[0] = 0;
  current_edges[0] = 0;
  while (true) {
    
    std::cout << "Current attempted polygon:\n";
    for (i=0; i<current_len; i++) {
      std::cout << real_edges_beginning_with[current_beginning_letters[i]][current_edges[i]] << " ";
    }
    std::cout << "\n";
    
    if (current_len > 1) {
      //check if it's allowed to close up
      temp1 = edges[ 
                    real_edges_beginning_with
                                   [current_beginning_letters[current_len-1]]
                                   [current_edges[current_len-1]] 
                   ].last;
      temp2 = edges[ 
                    real_edges_beginning_with
                                   [current_beginning_letters[0]]
                                   [current_edges[0]] 
                    ].first;
      if (temp2 != current_beginning_letters[0]) {
        std::cout << "Badness in polygons creation\n";
      }
      if (C.next_letter(temp1) == temp2) {
        //it does close up, so add it
        temp_poly.edges.resize(current_len);
        for (i=0; i<current_len; i++) {
          temp_poly.edges[i] = real_edges_beginning_with[current_beginning_letters[i]][current_edges[i]];
        }
        polys.push_back(temp_poly);
      }
    }
    //now we advance it: if the list is shorter than maxmal (4), then add
    //one on.  otherwise, step back and leave it short
    //note when we extend it, make sure the index is larger than current_beginning
    if (current_len < 4) { 
      temp1 = edges[
                    real_edges_beginning_with
                         [current_beginning_letters[current_len-1]]
                         [current_edges[current_len-1]]
                    ].last;
      temp2 = C.next_letter(temp1);
      current_beginning_letters[current_len] = temp2;
      std::cout << "Trying to extend with beginning letter " << temp2 << "\n";
      for (i=0; i<(int)real_edges_beginning_with[temp2].size(); i++) {
        if (edges[real_edges_beginning_with[temp2][i]].last > current_beginning_letters[0]) {
          break;
        }
      }
      if (i < (int)real_edges_beginning_with[temp2].size()) {
        current_edges[current_len] = i;
        current_len++;
        continue;
      } else {
        //do nothing, we will have to advance this index
      }
    }
    //if we get here, we need to advance the index current_len-1           
    i = current_len-1;
    while (i >=0 && current_edges[i] == (int)real_edges_beginning_with[current_beginning_letters[i]].size()-1) {
      i--;
    }
    if (i==-1) {
      if (current_beginning_letters[0] == current_len-1) {
        break;
      } else {
        current_beginning_letters[0]++;
        current_edges[0] = 0;
        current_len = 1;
        continue;
      }
    }
    current_edges[i]++;
    current_len = i+1;
  }
      
  
  //OK we have made all the polys with two blanks, and all the polys with no 
  //blanks.  Now we go through, and to all the polys with length 3 or 4, we convert
  //each of the edges in turn into a blank one
  int old_number_of_polys = (int)polys.size();
  int poly_len;
  int edge_backup;
  for (i=num_two_blank_polys; i<old_number_of_polys; i++) {
    if (polys[i].edges.size() > 2) {
      temp_poly.edges = polys[i].edges;
      poly_len = (int)temp_poly.edges.size();
      for (j=0; j<poly_len; j++) {
        temp1 = edges[ temp_poly.edges[j] ].last;
        for (k=0; k<(int)blank_edges_beginning_with[temp1].size(); k++) {
          if (edges[ blank_edges_beginning_with[temp1][k] ].last == 
              edges[ temp_poly.edges[(j+2)%poly_len] ].first) {
            break;
          }
        }
        edge_backup = temp_poly.edges[(j+1)%poly_len];
        temp_poly.edges[(j+1)%poly_len] = blank_edges_beginning_with[temp1][k];
        polys.push_back(temp_poly);
        temp_poly.edges[(j+1)%poly_len] = edge_backup;
      }
    }
  }
  
}



void print_multiarcs(std::ostream &os, std::vector<std::vector<Multiarc> > &arcs) {
  int i,j,k;
  os << "Multiarcs:\n";
  for (i=0; i<(int)arcs.size(); i++) {
    os << "Group " << i << ":\n";
    for (j=0; j<(int)arcs[i].size(); j++) {
      os << "(";
      for (k=0; k<(int)arcs[i][j].letters.size(); k++) {
        os << arcs[i][j].letters[k] << ",";
      }
      os << ")\n";
    }
  }
}

void print_edges(std::ostream &os, std::vector<Edge> &edges) {
  int i;
  os << "Edges: \n";
  for (i=0; i<(int)edges.size(); i++) {
    os << "(" << edges[i].first << "," << edges[i].last << ")";
    if (edges[i].blank) {
      os << "b";
    }
    os << "\n";
  }
}

void print_polys(std::ostream &os, std::vector<Edge> &edges, std::vector<Polygon> &polys) {
  int i,j;
  os << "Polygons: \n";
  for (i=0; i<(int)polys.size(); i++) {
    os << "(";
    for (j=0; j<(int)polys[i].edges.size(); j++) {
      os << polys[i].edges[j];
      if (edges[polys[i].edges[j]].blank) {
        os << "!";
      }
      os << ",";
    }
    os << ")\n";
  }
}

int main(int argc, char* argv[]) {
  int current_arg = 1;
  if (argc < 3 || std::string(argv[1]) == "-h") {
    std::cout << "usage: ./scyllop [-h] <gen string> <chain>\n";
    std::cout << "\twhere <gen string> is of the form <gen1><order1><gen2><order2>...\n";
    std::cout << "\te.g. a5b0 computes in Z/5Z * Z\n";
    std::cout << "\tand <chain> is an integer linear combination of words in the generators\n";
    std::cout << "\te.g. ./scyllop a5b0 aabaaaB\n";
    std::cout << "\t-h: print this message\n";
    exit(0);
  }
  while (argv[current_arg][0] == '-') {
    current_arg++;
    //handle arguments (none yet)
  }
  
  std::string G_in = std::string(argv[current_arg]);
  CyclicProduct G(G_in);                                                         //create the group
  current_arg++;
  
  Chain C(&G, &argv[current_arg], argc-current_arg);                              //process the chain arguments
  
  std::vector<std::vector<Multiarc> > arcs(0);
  std::vector<Edge> edges(0);
  std::vector<Polygon> polys(0);

  std::cout << "Group: " << G << "\n";
  std::cout << "Chain: " << C << "\n";
  C.print_chunks(std::cout);
  
  std::cout << "Letters:\n";
  C.print_letters(std::cout);
  
  std::cout << "Group letters:\n";
  C.print_group_letters(std::cout);
  
  compute_multiarcs(G, C, arcs);                                                 //calls arcs and polys by reference
  print_multiarcs(std::cout, arcs);
  compute_edges(G, C, edges);
  print_edges(std::cout, edges);
  compute_polys(G, C, edges, polys);
  print_polys(std::cout, edges, polys);
  
  rational scl;
  std::vector<rational> solution_vector(polys.size());                           //run the LP
  scyllop_lp(G, C, arcs, polys, &scl, &solution_vector, GLPK_DOUBLE, 0); 
  
  std::cout << "scl( " << C << ") = " << scl << " = " << scl.get_d() << "\n";    //output the answer
  
  return 0;
}
