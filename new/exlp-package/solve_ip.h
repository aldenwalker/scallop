#ifndef SOLVE_IP_H
#define SOLVE_IP_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include "lpstruct.h"

#define IP_RESULT_OPTIMAL     10
#define IP_RESULT_UNBOUNDED   11
#define IP_RESULT_INFEASIBLE  12
#define IP_RESULT_UNSUPPORTED 4

int  solve_ip(LP* lp);

#endif
