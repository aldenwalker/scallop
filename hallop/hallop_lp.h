#ifndef HALLOP_LP_H
#define HALLOP_LP_H

#include "free_group_chain.h"
#include "pieces.h"
#include "../lp.h"
#include "../rational.h"
  
namespace HALLOP {
  void hallop_lp(HALLOP::FreeGroupChain& C,
                 HALLOP::Pieces& P, 
                 SparseLPSolver solver,
                 Rational& scl, 
                 std::vector<Rational>& soln_vec, 
                 int verbose);
}

#endif


