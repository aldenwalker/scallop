#ifndef SOLVE_LP_H
#define SOLVE_LP_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include "lpstruct.h"

void reinversion(LP* lp);
int  solve_lp(LP* lp);
int  solve_lp_primal(LP* lp);
int  solve_lp_dual(LP* lp);
int solve_lp_core_dual(LP* lp);

#endif
