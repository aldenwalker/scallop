#include <string>
#include <iostream>
#include <cstdlib>

#include "../lp.h"
#include "../rational.h"

#include "free_group_chain.h"
#include "pieces.h"
#include "hallop_lp.h"
#include "hallop.h"


void HALLOP::hallop(int argc, char** argv) {
  SparseLPSolver solver = GLPK_SIMPLEX;
  int verbose = 1;
  bool lp_verbose = false;
  
  if (argc < 1 || std::string(argv[0]) == "-h") {
    std::cout << "usage: ./scallop -hyp [-m<GLPK,GIPT,EXLP,GUROBI>] [-v[n]] [-R<relator>] <chain>\n";
    std::cout << "\twhere <chain> allows integral weights on the words\n";
    std::cout << "\t-h: print this message\n";
    std::cout << "\t-m: use the LP solver specified (EXLP uses exact arithmetic)\n";
    std::cout << "\t-v[n]: verbosity (if -v isn't used, n=1, if -v but no n, then n=2)\n";
    std::cout << "\t-R relator: add a relator\n";
    std::cout << "\tExample: ./scallop -hyp -RabABcdCD abAB\n";
    exit(0);
  }
  
  std::vector<std::string> relators(0);
  
  int current_arg = 0;
  while (argv[current_arg][0] == '-') {
    if (argv[current_arg][1] == 'R') {
      relators.push_back( std::string(&argv[current_arg][2]) );
    
    } else if (argv[current_arg][1] == 'm') {
      switch (argv[current_arg][3]) {
        case 'L':
          solver = GLPK_SIMPLEX; break;
        case 'I':
          solver = GLPK_IPT; break;
        case 'X':
          solver = EXLP; break;
        case 'U':
          solver = GUROBI; break;
      }
    
    } else if (argv[current_arg][1] == 'v') {
      if (argv[current_arg][2] == '\0') {
        verbose = 2;
      } else {
        verbose = atoi(&argv[current_arg][2]);
      }
    
    } else if (argv[current_arg][1] == 'V') {
      lp_verbose = true;
    }
    current_arg++;
  }
  
  FreeGroupChain C(&argv[current_arg], argc-current_arg);
  
  if (verbose > 1) {
    std::cout << "Got relators:\n";
    for (int i=0; i<(int)relators.size(); ++i) {
      std::cout << i << ": " << relators[i] << "\n";
    }
    std::cout << "Got chain: " << C << "\n";
  }
  
  FreeGroupChain CR = C; 
  
  for (int i=0; i<(int)relators.size(); ++i) {
    CR.add_relator(relators[i]);
  }
  
  if (verbose>1) {
    std::cout << "Chain with relators: " << CR << "\n";
  }
  
  Pieces P(CR);
  
  if (verbose > 2) {
    std::cout << "Pieces from chain:\n";
    P.print(std::cout);
  }
  
  std::vector<Rational> soln_vec(0); //rectangles, then triangles
  Rational scl;
  
  hallop_lp(CR, P, solver, scl, soln_vec, lp_verbose);
  
  if (verbose > 0) {
    std::cout << "Lower bound on scl: " << scl << "\n";
  }
  
  if ((soln_vec.size() < 50 && verbose > 1) || verbose > 2) {
    std::cout << "Solution vector:\n";
    for (int i=0; i<(int)soln_vec.size(); ++i) {
      if (soln_vec[i] == 0) continue;
      if (i < (int)P.rects.size()) {
        std::cout << i << ": " << soln_vec[i] << " * " << P.rects[i] << "\n";
      } else {
        std::cout << i << ": " << soln_vec[i] << " * " << P.tris[i-P.rects.size()] << "\n";
      }
    }
  }
}











