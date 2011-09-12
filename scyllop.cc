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

#include "scyllop_classes.h"
#include "scyllop_lp.h"
#include "rational.h"

/*****************************************************************************
 * make the list of multiarcs -- each multiarc is for a different group, and 
 * each one is made up of a collection (possibly including duplicates) of 
 * chunks, where the total number of letters is the order of the group.  
 * only the cyclic order in the multiarc matters.  
 * ***************************************************************************/
void compute_multiarcs(CyclicProduct &G, Chain &C, std::vector<Multiarc> &arcs) {
  std::vector<char> gens = G.gen_list();
  std::vector<int> orders = G.order_list();
  int i,j,k;
  int num_groups = gens.size();
  std::vector<ChainChunk> group_chunks(0);
  std::string word;
  ChainChunk temp_chunk;
  Multiset chunk_selection;
  
  arcs.resize(0);
  
  for (i=0; i<num_groups; i++) {
    
    //create a list of chunks in that group
    group_chunks.resize(0);
    for (j=0; j<C.num_words(); j++) {
      word = C[j];
      temp.word = j;
      k=0;
      while (k<word.size()) {
        while (word[k] != gens[i]) {
          k++;
        }
        temp.start_index = k;
        while (word[k] == gens[i]) {
          k++;
        }
        temp.len = k - temp.start_index;
        group_chunks.push_back(temp);
      }
    }
    
    //now go through and create all possible multiarcs; we may assume that the
    //smallest index is first
    for (first_index=0; first_index<group_chunks.size(); first_index++) {
      chunk_selection = Multiset(orders[i], first_index, group_chunks.size());
      do {
        //try this chunk slection
        
      } while (1 != chunk_selection.next());
    }
    
  }

} 

void compute_polys(CyclicProduct &G, 
                   Chain &C, 
                   std::vector<Multiarc> &arcs, 
                   std::vector<Polygon> &polys) {
  
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
  
  CyclicProduct G(std::string(argv[current_arg]));                               //create the group
  current_arg++;
  
  Chain C(&argv[current_arg], argc-current_arg);                                 //process the chain arguments
  
  std::vector<Multiarc> arcs(0);
  std::vector<Polygon> polys(0);

  std::cout << "Group: " << G << "\n";
  std::cout << "Chain: " << C << "\n";
  std::cout.flush();
  
  compute_multiarcs(G, C, arcs);                                                 //calls arcs and polys by reference
  compute_polys(G, C, arcs, polys);
  
  rational scl;
  std::vector<rational> solution_vector(polys.size());                           //run the LP
  scyllop_lp(G, C, arcs, polys, &scl, &solution_vector, GLPK_DOUBLE, 0); 
  
  std::cout << "scl( " << C << ") = " << scl << " = " << scl.get_d() << "\n";    //output the answer
  
  return 0;
}
