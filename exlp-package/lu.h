#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include "lpstruct.h"
#include "mylib.h"

void lp_LU_basis(LP* lp);
  /* lp->B をLU分解. 結果は lp->eta に.
     lp->eta->U は上三角行列.
     lp->eta->P[i] は i 番目に i 行と swap する行.
     lp->eta->L[i] は i 番目にかけるエータ行列. */

int  lp_check_redundancy(LP* lp);
void lp_select_dual_feasible_basis(LP* lp);
int lp_select_basis(LP* lp);
