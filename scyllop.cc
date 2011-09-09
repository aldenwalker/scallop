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


void compute_multiarcs(CyclicProduct &G, Chain &C, std::vector<Multiarc> &arcs) {
  
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
