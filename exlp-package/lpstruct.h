#ifndef LPSTRUCT_H
#define LPSTRUCT_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include "EXLPvector.h"
#include "matrix.h"
#include "eta_file.h"
#include "hash.h"

#define LP_NAME_LEN_MAX 256

#define LP_EQUALITY_EQ        'E'
#define LP_EQUALITY_LE        'L'
#define LP_EQUALITY_GE        'G'

#define LP_VAR_TYPE_NORMAL     0
#define LP_VAR_TYPE_SLACK      1
#define LP_VAR_TYPE_ARTIFICIAL 2
#define LP_VAR_TYPE_CONST      3

#define LP_PIVOT_BLAND         1
#define LP_PIVOT_BLAND2        2
#define LP_PIVOT_DANTZIG       3
#define LP_PIVOT_DEVEX         8
#define LP_PIVOT_STEEPEST_EDGE 16
#define LP_PIVOT_PROJECTED_STEEPEST_EDGE 32
#define LP_PIVOT_MIX_HEURISTIC 32768

#define LU_MARKOWITZ           0
#define LU_MINIMUM_DEGREE      1

#define LP_RESULT_OPTIMAL      0
#define LP_RESULT_UNBOUNDED    1
#define LP_RESULT_INFEASIBLE   2
#define LP_RESULT_DUAL_INFEASIBLE 3
#define LP_RESULT_UNSUPPORTED  4

typedef struct {
  int*   is_valid;
  mpq_t* bound;
} bounds;

typedef struct {
  FILE*   fp;
  char*   name;
  char*   obj_name;
  char**  row_name;
  char**  var_name;
  int*    row_equality;
  int*    var_type;
  int     rows;
  int     vars;
  int     maximize;
  int     print_sol;
  int     verbose;
  int     perturbation;
  int     scaling;
  int     preprocess;
  int     pivot_rule;
  int     lu_rule;
  int     dual_simplex;
  int     hash_entries;
  int     reinversion_cycle;
  int     phase;
  int     gomory;
  matrix* A;
  EXLPvector* b;
  EXLPvector* b_back;
  EXLPvector* x;
  EXLPvector* xb;
  EXLPvector* c;
  EXLPvector* cb;
  EXLPvector* c_back;
  mpq_t   c_const;
  bounds  upper;
  bounds  lower;
  int*    is_basis;
  int*    basis_column;
  eta_file* eta;
  double* steepest_edge_table;
  double* devex_weight;
  int*    framework;
  EXLPvector_d** A_d;
  EXLPvector_d* c_d;

  int*    is_integer;

  hash_str *hash_str_row_name;
  hash_str *hash_str_var_name;

  int     slacks;
  int*    slackref;

  mpq_t   q_work;

  void    *owner;

} LP;


LP*  new_lp(void *owner);
void lp_free(LP* lp);
void lp_set_name(LP* lp, char* name);
void lp_set_obj_name(LP* lp, char* name);
void lp_resize(LP* lp, int rows, int columns, int usenames);
int  lp_add_row(LP* lp, char* name);
int  lp_add_var(LP* lp, char* name);
void lp_remove_row(LP* lp, int row);
int  lp_get_row_num(LP* lp, char* name);
int  lp_get_var_num(LP* lp, char* name);
void lp_set_coefficient(LP* lp, mpq_t c, int row, int var);
void lp_set_row_equality(LP* lp, int row, int equality);
void lp_set_rhs(LP* lp, int row, mpq_t q);
void lp_set_row_range(LP* lp, int row, mpq_t q);
void lp_get_object_value(LP* lp, mpq_t* q);
int  lp_artificials_are_zeros(LP* lp);
int  is_normal_var(LP* lp, int var);
int  is_slack_var(LP* lp, int var);
int  is_artificial_var(LP* lp, int var);
int  is_const_var(LP* lp, int var);
int  is_free_var(LP* lp, int var);
void swap_basis_columns(LP* lp, int col1, int col2);
void lp_arrange_inequality(LP* lp);
void lp_arrange_inequality2(LP* lp);
void lp_add_slackref(LP* lp, int rowno);
void lp_add_slacks(LP* lp);
void get_valid_rhs(mpq_t* q, LP* lp, int i);
void lp_sol_init(LP* lp);
void lp_sol_init_dual(LP* lp);
void lp_set_rhs_positive(LP* lp);
void lp_add_artificials(LP* lp);
int  lp_normal_vars(LP* lp);
int  lp_slack_vars(LP* lp);
int  lp_artificial_vars(LP* lp);
void lp_print_basis(LP* lp);
void lp_print_nonzero_vars(LP* lp);

int lp_add_var_without_A(LP* lp, char* name);
/* 応急処置です... matrix_add_column(lp->A) しない以外は lp_add_var()
   と同じです. */

void lp_LU_basis(LP* lp);

void lp_hash_str_init(LP* lp, int hash_entries);

#endif
