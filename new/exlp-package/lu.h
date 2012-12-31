#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include "lpstruct.h"
#include "mylib.h"

void lp_LU_basis(LP* lp);
  /* lp->B ��LUʬ��. ��̤� lp->eta ��.
     lp->eta->U �Ͼ廰�ѹ���.
     lp->eta->P[i] �� i ���ܤ� i �Ԥ� swap �����.
     lp->eta->L[i] �� i ���ܤˤ����륨��������. */

int  lp_check_redundancy(LP* lp);
void lp_select_dual_feasible_basis(LP* lp);
int lp_select_basis(LP* lp);
