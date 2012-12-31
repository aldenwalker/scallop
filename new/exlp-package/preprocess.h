#ifndef PREPROCESS_H
#define PREPROCESS_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include "lpstruct.h"

void scaling(LP* lp);
void remove_artificials_from_basis(LP* lp);
int  non_const_vars_in_row(LP* lp, int row, int* var);
int  check_parallel_constraints(LP* lp);
void check_parallel_columns(LP* lp);
void check_trivial_constraints0(LP* lp);
int  check_trivial_constraints1(LP* lp);
int  check_trivial_constraints2(LP* lp);
void check_trivial_constraints3(LP* lp);
void check_tight_constraints1(LP* lp);
void check_tight_constraints2(LP* lp);
int  check_const_vars(LP* lp);

#endif
