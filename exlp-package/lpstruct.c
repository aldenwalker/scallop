#include "lpstruct.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mylib.h"
#include "lu.h"
#include "solve_lp.h"
//#include "snprintf.h"


LP* new_lp(void *owner) {
  LP* lp = my_calloc(1, sizeof(LP));
//  LP* lp = my_malloc(sizeof(LP));

  lp->name = my_malloc(LP_NAME_LEN_MAX*sizeof(char));
  lp->obj_name = my_malloc(LP_NAME_LEN_MAX*sizeof(char));
  lp->row_name = NULL;
  lp->hash_str_row_name = NULL;
  lp->var_name = NULL;
  lp->hash_str_var_name = NULL;
  lp->row_equality = NULL;
  lp->var_type = NULL;
  lp->A  = new_matrix(0, 0);
  lp->b  = new_vector(0);
  lp->x  = new_vector(0);
  lp->xb = new_vector(0);
  lp->c  = new_vector(0);
  lp->cb = new_vector(0);
  lp->c_back = new_vector(0);
  mpq_init(lp->c_const);
  lp->upper.is_valid = NULL;
  lp->lower.is_valid = NULL;
  lp->upper.bound = NULL;
  lp->lower.bound = NULL;
  lp->is_basis = NULL;
  lp->basis_column = NULL;

  lp_set_name(lp, "my_lp");
  lp_set_obj_name(lp, "obj");

  lp->rows = 0;
  lp->vars = 0;
  lp->maximize     = FALSE;
  lp->print_sol    = FALSE;
  lp->verbose      = FALSE;
  lp->perturbation = FALSE;
  lp->scaling      = TRUE;
  lp->preprocess   = 1;
  lp->hash_entries = 20011;
  lp->pivot_rule   = LP_PIVOT_MIX_HEURISTIC|LP_PIVOT_STEEPEST_EDGE;
  lp->lu_rule      = LU_MARKOWITZ;
  lp->dual_simplex = FALSE;
  lp->reinversion_cycle = 0;//16;//20;
  lp->gomory       = FALSE;

  lp->steepest_edge_table = NULL;
  lp->devex_weight = NULL;
  lp->framework    = NULL;

  lp->A_d = NULL;

  lp->is_integer = NULL;

  lp->slacks = 0;
  lp->slackref = NULL;

  mpq_init(lp->q_work);

  lp->owner = owner;

  return lp;
}

void lp_free(LP* lp) {
  int  i;

  free(lp->name);
  free(lp->obj_name);
  if (lp->row_name != NULL) {
    for (i = 0; i < lp->rows; i ++)
      free(lp->row_name[i]);
    free(lp->row_name);
    hash_str_free(&lp->hash_str_row_name);
  }
  if (lp->var_name != NULL) {
    for (i = 0; i < lp->vars; i ++)
      free(lp->var_name[i]);
    free(lp->var_name);
    hash_str_free(&lp->hash_str_var_name);
  }
  free(lp->row_equality);
  matrix_free(lp->A);
  vector_free(&lp->b);
  vector_free(&lp->x);
  vector_free(&lp->xb);
  vector_free(&lp->c);
  vector_free(&lp->cb);
  vector_free(&lp->c_back);
  mpq_clear(lp->c_const);
  for (i = 0; i < lp->vars; i ++) {
    if (lp->upper.is_valid[i])
      mpq_clear(lp->upper.bound[i]);
    if (lp->lower.is_valid[i])
      mpq_clear(lp->lower.bound[i]);
  }
  free(lp->upper.is_valid);
  free(lp->lower.is_valid);
  free(lp->upper.bound);
  free(lp->lower.bound);
  free(lp->is_basis);
  free(lp->basis_column);
  free(lp->is_integer);
  if(lp->slackref != NULL)
    free(lp->slackref);
  mpq_clear(lp->q_work);

  free(lp->A_d);
  vector_free(&lp->b_back);
  vector_d_free(&lp->c_d);
  free(lp->var_type);
  if (lp->eta)
    eta_file_free(lp->eta);
  if(lp->hash_str_row_name != NULL)
    hash_str_free(&lp->hash_str_row_name);
  if(lp->hash_str_var_name != NULL)
    hash_str_free(&lp->hash_str_var_name);
  my_hash_mpq_free();

  free(lp);
}

void lp_set_name(LP* lp, char* name) {
  strncpy(lp->name, name, LP_NAME_LEN_MAX);
}

void lp_set_obj_name(LP* lp, char* name) {
  strncpy(lp->obj_name, name, LP_NAME_LEN_MAX);
}

int lp_add_row_without_A(LP* lp, char* name) {
  lp->rows ++;

  if(name != NULL) {
    lp->row_name = my_realloc(lp->row_name, lp->rows*sizeof(char*));
    lp->row_name[lp->rows-1] = my_malloc(LP_NAME_LEN_MAX*sizeof(char));
    if (name == NULL || name[0] == 0) {
      snprintf(lp->row_name[lp->rows-1], LP_NAME_LEN_MAX,
               "r%d", lp->rows);
    } else {
      strncpy(lp->row_name[lp->rows-1], name, LP_NAME_LEN_MAX);
    }
    hash_str_find(lp->hash_str_row_name, name, lp->row_name, lp->rows-1);
  }

  lp->row_equality = my_realloc(lp->row_equality, lp->rows*sizeof(int));
  lp_set_row_equality(lp, lp->rows-1, LP_EQUALITY_EQ);

  lp->basis_column = my_realloc(lp->basis_column, lp->rows*sizeof(int));
  lp->basis_column[lp->rows-1] = 0;
  /* ↑は lp_remove_row() をする際初期化していることが必要なため */
  return(lp->rows - 1);
}

int lp_add_row(LP* lp, char* name) {
  int i = lp_add_row_without_A(lp, name);
  matrix_add_row(lp->A);
  return i;
}

void lp_remove_row(LP* lp, int row) {
  /* phase1 の前に呼ぶのはいいんだけど
     phase1 と phase2 の間に呼んで平気かな? 
     とりあえず lp->eta 以下の面倒みてない. */
  int  i;

  matrix_dec_row_dimension(lp->A, row);
  vector_dec_dimension(lp->b, row);
  vector_dec_dimension(lp->cb, row);
  vector_dec_dimension(lp->xb, row);
  if (0 <= lp->basis_column[row] && lp->basis_column[row] < lp->vars)
    lp->is_basis[lp->basis_column[row]] = FALSE;
  lp->rows --;
  for (i = row; i < lp->rows; i ++) {
    strncpy(lp->row_name[i], lp->row_name[i+1], LP_NAME_LEN_MAX);
    lp->row_equality[i] = lp->row_equality[i+1];
    lp->basis_column[i] = lp->basis_column[i+1];
  }
  free(lp->row_name[lp->rows]);
  lp->row_name = my_realloc(lp->row_name, lp->rows*sizeof(char*));
  lp->row_equality = my_realloc(lp->row_equality, lp->rows*sizeof(int));
  lp->basis_column = my_realloc(lp->basis_column, lp->rows*sizeof(int));
}

int lp_add_var_without_A(LP* lp, char* name) {
  /* 追加されるのは非負変数       */
  /* 戻り値は追加された変数の番号 */
  /* 応急処置です... matrix_add_column(lp->A) しない以外は lp_add_var()
     と同じです. */
  lp->vars ++;

  lp->var_type = my_realloc(lp->var_type, lp->vars*sizeof(int));
  lp->var_type[lp->vars-1] = LP_VAR_TYPE_NORMAL;

  if(name != NULL) {
    lp->var_name = my_realloc(lp->var_name, lp->vars*sizeof(char*));
    lp->var_name[lp->vars-1] = my_malloc(LP_NAME_LEN_MAX*sizeof(char));
    if (name == NULL || name[0] == 0) {
      snprintf(lp->var_name[lp->vars-1], LP_NAME_LEN_MAX,
               "v%d", lp->vars);
    } else {
      strncpy(lp->var_name[lp->vars-1], name, LP_NAME_LEN_MAX);
    }
    hash_str_find(lp->hash_str_var_name, name, lp->var_name, lp->vars-1);
  }

  vector_resize(lp->x, lp->vars);
  vector_resize(lp->c, lp->vars);
  vector_resize(lp->c_back, lp->vars);

  lp->is_basis   = my_realloc(lp->is_basis,   lp->vars*sizeof(int));
  lp->is_integer = my_realloc(lp->is_integer, lp->vars*sizeof(int));

  lp->upper.is_valid = my_realloc(lp->upper.is_valid, lp->vars*sizeof(int));
  lp->lower.is_valid = my_realloc(lp->lower.is_valid, lp->vars*sizeof(int));
  lp->upper.bound = my_realloc(lp->upper.bound, lp->vars*sizeof(mpq_t));
  lp->lower.bound = my_realloc(lp->lower.bound, lp->vars*sizeof(mpq_t));
  lp->upper.is_valid[lp->vars-1] = FALSE;
  lp->lower.is_valid[lp->vars-1] = TRUE;
  mpq_init(lp->lower.bound[lp->vars-1]);

  return lp->vars-1;
}

int lp_add_var(LP* lp, char* name) {
  /* 追加されるのは非負変数       */
  /* 戻り値は追加された変数の番号 */
  int i = lp_add_var_without_A(lp, name);
  matrix_add_column(lp->A);

  return i;
}

void lp_resize(LP* lp, int rows, int columns, int usenames) {
  int  i, oldrows = lp->rows, oldvars = lp->vars;
  char buf[LP_NAME_LEN_MAX];

  /* Row data */
  lp->rows = rows;

  if(usenames) {
    lp->row_name = my_realloc(lp->row_name, lp->rows*sizeof(char*));
    for(i = oldrows; i < rows; i ++) {
      snprintf(buf, LP_NAME_LEN_MAX, "r%d", i + 1);
      lp->row_name[i] = my_malloc(strlen(buf)+1);
      hash_str_find(lp->hash_str_row_name, buf, lp->row_name, i);
    }
  }

  lp->row_equality = my_realloc(lp->row_equality, lp->rows*sizeof(int));
  lp->basis_column = my_realloc(lp->basis_column, lp->rows*sizeof(int));

  for(i = oldrows; i < rows; i ++) {
    lp_set_row_equality(lp, i, LP_EQUALITY_EQ);
    lp->basis_column[i] = 0;
  }

  /* Column data */
  lp->vars = columns;

  if(usenames) {
    lp->var_name = my_realloc(lp->var_name, lp->vars*sizeof(char*));
    for(i = oldvars; i < columns; i ++) {
      snprintf(buf, LP_NAME_LEN_MAX, "v%d", i + 1);
      lp->var_name[i] = my_malloc(strlen(buf)+1);
      hash_str_find(lp->hash_str_var_name, buf, lp->var_name, i);
    }
  }

  lp->var_type = my_realloc(lp->var_type, lp->vars*sizeof(int));

  vector_resize(lp->x, lp->vars);
  vector_resize(lp->c, lp->vars);
  vector_resize(lp->c_back, lp->vars);

  lp->is_basis   = my_realloc(lp->is_basis,   lp->vars*sizeof(int));
  lp->is_integer = my_realloc(lp->is_integer, lp->vars*sizeof(int));

  lp->upper.is_valid = my_realloc(lp->upper.is_valid, lp->vars*sizeof(int));
  lp->lower.is_valid = my_realloc(lp->lower.is_valid, lp->vars*sizeof(int));
  lp->upper.bound = my_realloc(lp->upper.bound, lp->vars*sizeof(mpq_t));
  lp->lower.bound = my_realloc(lp->lower.bound, lp->vars*sizeof(mpq_t));

  for(i = oldvars; i < columns; i ++) {
    lp->var_type[i] = LP_VAR_TYPE_NORMAL;
    lp->upper.is_valid[i] = FALSE;
    lp->lower.is_valid[i] = TRUE;
    mpq_init(lp->lower.bound[i]);
  }

  /* Matrix data */
  matrix_resize(lp->A, lp->rows, lp->vars);
#if 1
  vector_resize(lp->b, lp->rows);
  vector_resize(lp->xb, lp->rows);
  vector_resize(lp->cb, lp->rows);
#endif

}

void lp_hash_str_init(LP *lp, int hash_entries) {
  if(lp->hash_str_row_name != NULL)
    hash_str_free(&lp->hash_str_row_name);
  lp->hash_str_row_name = hash_str_init(hash_entries);
  if(lp->hash_str_var_name != NULL)
    hash_str_free(&lp->hash_str_var_name);
  lp->hash_str_var_name = hash_str_init(hash_entries);
}

int lp_get_row_num(LP* lp, char* name) {
  if (strncmp(lp->obj_name, name, LP_NAME_LEN_MAX) == 0)
    return lp->rows;
  return hash_str_find(lp->hash_str_row_name, name, lp->row_name, -1);
}

int lp_get_var_num(LP* lp, char* name) {
  return hash_str_find(lp->hash_str_var_name, name, lp->var_name, -1);
}

void lp_set_coefficient(LP* lp, mpq_t c, int row, int var) {
  if (row == lp->rows)
    vector_set_element(lp->c, c, var);
  else
    matrix_set_element(lp->A, c, row, var);
}

void lp_set_row_equality(LP* lp, int row, int equality) {
  lp->row_equality[row] = equality;
}

void lp_set_rhs(LP* lp, int row, mpq_t q) {
  if (0 <= row && row < lp->b->dimension)
    vector_set_element(lp->b, q, row);
}

void lp_set_row_range(LP* lp, int row, mpq_t q) {
  mpq_t q2;
  int   var;

  mpq_init(q2);

  if(lp->row_name != NULL) {
    char buf[LP_NAME_LEN_MAX];
    buf[0] = '#'; buf[1] = 'R'; buf[2] = 0;
    strncat(buf, lp->row_name[row], LP_NAME_LEN_MAX-3);
    var = lp_add_var(lp, buf);
  }
  else
    var = lp_add_var(lp, NULL);
  lp_add_slackref(lp, row+1);
  vector_get_element(&q2, lp->b, row);
  lp->upper.is_valid[var] = TRUE;
  mpq_init(lp->upper.bound[var]);
  if (lp->row_equality[row] == LP_EQUALITY_GE) {
    if (mpq_sgn(q) < 0)
      mpq_neg(q, q);
    mpq_set(lp->lower.bound[var], q2);
    mympq_add(q2, q2, q);
    mpq_set(lp->upper.bound[var], q2);
  } else if (lp->row_equality[row] == LP_EQUALITY_LE) {
    if (mpq_sgn(q) < 0)
      mpq_neg(q, q);
    mpq_set(lp->upper.bound[var], q2);
    mympq_sub(q2, q2, q);
    mpq_set(lp->lower.bound[var], q2);
  } else if (mpq_sgn(q) > 0) {
    mpq_set(lp->lower.bound[var], q2);
    mympq_add(q2, q2, q);
    mpq_set(lp->upper.bound[var], q2);        
  } else {
    mpq_set(lp->upper.bound[var], q2);
    mympq_add(q2, q2, q);
    mpq_set(lp->lower.bound[var], q2);
  }
  lp_set_row_equality(lp, row, LP_EQUALITY_EQ);
  lp_set_coefficient(lp, mympq_minus_one, row, var);
  vector_delete_element(lp->b, row);

  mpq_clear(q2);
}

void lp_get_object_value(LP* lp, mpq_t* q) {
  vector_inner_product(q, lp->c_back, lp->x);
  mympq_add(*q, *q, lp->c_const);
  if (lp->maximize == FALSE)
    mpq_neg(*q, *q);
}
/*
void lp_get_dual_object_value(LP* lp, mpq_t* q, EXLPvector* y) {
  EXLPvector* w;
  mpq_t  tmp;
  int  i, j;

  w = new_vector(lp->vars);
  mpq_init(tmp);

  mpq_set_si(q, 0, 1);

  for (i = 0; i < y->nonzeros; i ++) {
    get_valid_rhs(&tmp, lp, y->i[i]);
    mympq_mul(tmp, tmp, y->value[i]);
  }

  for (i = 0; i < lp->vars; i ++) {
  }

  vector_free(&w);
  mpq_clear(tmp);
}
*/
int lp_artificials_are_zeros(LP* lp) {
  int i;

  for (i = 0; i < lp->xb->nonzeros; i ++) {
    if (is_artificial_var(lp, lp->basis_column[lp->xb->i[i]]))
      return FALSE;
  }
  return TRUE;
}

int is_normal_var(LP* lp, int var) {
  if (lp->A->column[var]->nonzeros == 0 &&
      vector_element_is_zero(lp->c, var))
    return FALSE;
  return (lp->var_type[var] == LP_VAR_TYPE_NORMAL);
}

int is_slack_var(LP* lp, int var) {
  return (lp->var_type[var] == LP_VAR_TYPE_SLACK);
}

int is_artificial_var(LP* lp, int var) {
  return (lp->var_type[var] == LP_VAR_TYPE_ARTIFICIAL);
}

int is_const_var(LP* lp, int var) {
  return (lp->var_type[var] == LP_VAR_TYPE_CONST);
}

int is_free_var(LP* lp, int var) {
  return !(lp->upper.is_valid[var] || lp->lower.is_valid[var]);
}

void swap_basis_columns(LP* lp, int col1, int col2) {
  int  tmp;
  double  t;

  tmp = lp->basis_column[col1];
  lp->basis_column[col1] = lp->basis_column[col2];
  lp->basis_column[col2] = tmp;
  vector_swap_elements(lp->xb, col1, col2);
  vector_swap_elements(lp->cb, col1, col2);
  if (lp->dual_simplex && lp->phase == 2) {
    t = lp->steepest_edge_table[col1];
    lp->steepest_edge_table[col1] = lp->steepest_edge_table[col2];
    lp->steepest_edge_table[col2] = t;
  }
}

void lp_arrange_inequality(LP* lp) {
  int  i;

  for (i = 0; i < lp->rows; i ++) {
    if (lp->row_equality[i] == LP_EQUALITY_GE) {
      matrix_rev_row_sgn(lp->A, i);
      vector_rev_element_sgn(lp->b, i);
      lp->row_equality[i] = LP_EQUALITY_LE;
    }
  }
}

void lp_arrange_inequality2(LP* lp) {
  int  i, j, col;

  lp_arrange_inequality(lp);

  for (i = 0; i < lp->rows; i ++) {
    if (lp->row_equality[i] != LP_EQUALITY_LE)
      continue;
    for (j = 0; j < lp->A->row[i]->nonzeros; j ++) {
      col = lp->A->row[i]->i[j];
      if (is_const_var(lp, col))
        continue;
      if (lp->A->column[col]->nonzeros != 1)
        continue;
      if (mpq_sgn(*matrix_get_element_ptr(lp->A, i, col)) < 0)
        continue;
      if (mpq_sgn(*vector_get_element_ptr(lp->c, col)) < 0)
        continue;
      if (!lp->lower.is_valid[col] || mpq_sgn(lp->lower.bound[col]) < 0)
        continue;
      if (lp->upper.is_valid[col])
        continue;
      lp->row_equality[i] = LP_EQUALITY_EQ;
      putchar('_');
      break;
    }
  }
}

void lp_add_slackref(LP* lp, int rowno) {
  int vecsize = 40*(lp->slacks/40);

  if(lp->slacks >= vecsize) {
    vecsize += 40;
    lp->slackref = (int *) realloc(lp->slackref, sizeof(int)*vecsize);
  }
  lp->slackref[lp->slacks++] = rowno;
}

void lp_add_slacks(LP* lp) {
  char buf[LP_NAME_LEN_MAX];
  int   i, j;

  for (i = 0; i < lp->vars; i ++) {
    if (is_const_var(lp, i) ||
        !vector_element_is_zero(lp->c, i) ||
        lp->A->column[i]->nonzeros != 1)
      continue;
    j = lp->A->column[i]->i[0];
    if (lp->row_equality[j] != LP_EQUALITY_LE)
      continue;
    if (mpq_sgn(*matrix_get_element_ptr(lp->A, j, i)) > 0 &&
        lp->upper.is_valid[i])
      continue;
    if (mpq_sgn(*matrix_get_element_ptr(lp->A, j, i)) < 0 &&
        lp->lower.is_valid[i])
      continue;
    lp->row_equality[j] = LP_EQUALITY_EQ;
  }

  j = 0;
  for (i = 0; i < lp->rows; i ++) {
    if (lp->row_equality[i] == LP_EQUALITY_LE)
      j ++;
  }
  matrix_resize(lp->A, lp->rows, lp->vars+j);

  for (i = 0; i < lp->rows; i ++) {
    if (lp->row_equality[i] == LP_EQUALITY_LE) {
      if(lp->row_name != NULL) {
        strncpy(buf, "#S", LP_NAME_LEN_MAX);
        strncat(buf, lp->row_name[i], LP_NAME_LEN_MAX);
        j = lp_add_var_without_A(lp, buf);
      }
      else
        j = lp_add_var_without_A(lp, NULL);
      lp_add_slackref(lp, i+1);
      lp->var_type[j] = LP_VAR_TYPE_SLACK;
      lp_set_coefficient(lp, mympq_one, i, j);
      lp->row_equality[i] = LP_EQUALITY_EQ;
    }
  }
}

void get_valid_rhs(mpq_t* q, LP* lp, int i) {
  mpq_t  r;
  int  j, var;

  mpq_init(r);

  vector_get_element(q, lp->b, i);

  for (j = 0; j < lp->A->row[i]->nonzeros; j ++) {
    var = lp->A->row[i]->i[j];
    if (is_const_var(lp, var))
      continue;
    vector_get_element(&r, lp->x, var);
    mympq_mul(r, r, *matrix_get_element_ptr(lp->A, i, var));
    mympq_sub(*q, *q, r);
  }

  mpq_clear(r);
}

void lp_sol_init(LP* lp) {
  mpq_t* q;
  int  i;

  for (i = 0; i < lp->vars; i ++) {
    lp->is_basis[i] = FALSE;
    if (is_const_var(lp, i))
      continue;
    if (lp->lower.is_valid[i] && lp->upper.is_valid[i]) {
      q = vector_get_element_ptr(lp->c, i);
      if (mpq_sgn(*q) > 0)
        vector_set_element(lp->x, lp->upper.bound[i], i);
      else
        vector_set_element(lp->x, lp->lower.bound[i], i);
    } if (lp->lower.is_valid[i]) {
      vector_set_element(lp->x, lp->lower.bound[i], i);
    } else if (lp->upper.is_valid[i]) {
      vector_set_element(lp->x, lp->upper.bound[i], i);
    }
  }
}

void lp_sol_init_dual(LP* lp) {
  mpq_t  q1, q2;
  mpq_t* q3;
  EXLPvector* v1;
  EXLPvector* v2;
  int  i, j, k;

  if (lp->preprocess < 2) {
    mpq_init(q1);
    //mpq_set_si(q1, 1, 10000);
    for (i = 0; i < lp->vars; i ++) {
      if (is_const_var(lp, i))
        continue;
      mpq_set_si(q1, i, lp->vars*lp->vars*lp->vars*lp->vars);
      if (lp->lower.is_valid[i] && lp->upper.is_valid[i]) {
	q3 = vector_get_element_ptr(lp->c, i);
	if (mpq_sgn(*q3) > 0) {
	  vector_set_element(lp->x, lp->upper.bound[i], i);
          if (lp->is_basis[i]) {
            vector_sub_element(lp->x, q1, i);
            if (mpq_cmp(*vector_get_element_ptr(lp->x, i),
                        lp->lower.bound[i]) < 0)
              vector_add_element(lp->x, q1, i);
	  }
	} else {
	  vector_set_element(lp->x, lp->lower.bound[i], i);
          if (lp->is_basis[i]) {
            vector_add_element(lp->x, q1, i);
            if (mpq_cmp(*vector_get_element_ptr(lp->x, i),
                        lp->upper.bound[i]) > 0)
              vector_sub_element(lp->x, q1, i);
	  }
	}
      } if (lp->lower.is_valid[i]) {
	vector_set_element(lp->x, lp->lower.bound[i], i);
        if (lp->is_basis[i])
          vector_add_element(lp->x, q1, i);
      } else if (lp->upper.is_valid[i]) {
	vector_set_element(lp->x, lp->upper.bound[i], i);
        if (lp->is_basis[i])
          vector_sub_element(lp->x, q1, i);
      }
    }
    mpq_clear(q1);

    for (i = 0; i < lp->rows; i ++) {
      vector_set_element(lp->xb, *vector_get_element_ptr(lp->x, lp->basis_column[i]), i);
    }
  } else {

    reinversion(lp);

    mpq_init(q1);
    mpq_init(q2);
    v1 = new_vector(lp->rows);
    v2 = new_vector(lp->vars);

    eta_file_btran(lp->eta, lp->cb, v1);
    vector_copy(v2, lp->c);
    for (i = 0; i < lp->vars; i ++) {
      if (lp->is_basis[i] || is_const_var(lp, i))
	continue;
      vector_inner_product(&q1, v1, lp->A->column[i]);
      vector_sub_element(v2, q1, i);
    }

    for (i = 0; i < lp->vars; i ++) {
      if (is_const_var(lp, i) || lp->is_basis[i])
	continue;

      if (lp->lower.is_valid[i] && lp->upper.is_valid[i]) {
	q3 = vector_get_element_ptr(v2, i);
	if (mpq_sgn(*q3) > 0)
	  vector_set_element(lp->x, lp->upper.bound[i], i);
	else
	  vector_set_element(lp->x, lp->lower.bound[i], i);

      } else if (lp->lower.is_valid[i]) {
	vector_set_element(lp->x, lp->lower.bound[i], i);

      } else if (lp->upper.is_valid[i]) {
	vector_set_element(lp->x, lp->upper.bound[i], i);
      }
    }

    vector_resize(v2, lp->rows);

    vector_copy(v1, lp->b);
    for (i = 0; i < lp->rows; i ++) {
      mpq_set_si(q1, 0, 1);
      for (j = 0; j < lp->A->row[i]->nonzeros; j ++) {
	k = lp->A->row[i]->i[j];
	if (is_const_var(lp, k) || lp->is_basis[k])
	  continue;
	mympq_mul(q2, *matrix_get_element_ptr(lp->A, i, k),
		  *vector_get_element_ptr(lp->x, k));
	mympq_add(q1, q1, q2);
      }
      vector_sub_element(v1, q1, i);
    }
    eta_file_ftran(lp->eta, v1, v2, NULL);

//mpq_set_si(q1, 1, 10000);
    mpq_set_si(q1, i, lp->vars*lp->vars*lp->vars*lp->vars);
    for (i = 0; i < lp->rows; i ++) {
      k = lp->basis_column[i];
      q3 = vector_get_element_ptr(v2, i);
      //mpq_set_si(q1, 1, 100000);mympq_mul(q1, q1, *q3);
      if (lp->lower.is_valid[k] &&
          mpq_cmp(*q3, lp->lower.bound[k]) <= 0) {
        vector_set_element(lp->x, lp->lower.bound[k], k);
        vector_set_element(lp->xb, lp->lower.bound[k], i);
        vector_add_element(lp->x, q1, k);
        vector_add_element(lp->xb, q1, i);
        if (lp->upper.is_valid[k] &&
            mpq_cmp(lp->upper.bound[k],
                    *vector_get_element_ptr(lp->x, k)) < 0) {
          vector_sub_element(lp->x, q1, k);
          vector_sub_element(lp->xb, q1, i);
        }
        continue;
      }
      if (lp->upper.is_valid[k] &&
          mpq_cmp(*q3, lp->upper.bound[k]) >= 0) {
        vector_set_element(lp->x, lp->upper.bound[k], k);
        vector_set_element(lp->xb, lp->upper.bound[k], i);
        vector_sub_element(lp->x, q1, k);
        vector_sub_element(lp->xb, q1, i);
        if (lp->lower.is_valid[k] &&
            mpq_cmp(lp->lower.bound[k],
                    *vector_get_element_ptr(lp->x, k)) > 0) {
          vector_add_element(lp->x, q1, k);
          vector_add_element(lp->xb, q1, i);
        }
        continue;
      }
      vector_set_element(lp->x, *q3, k);
      vector_set_element(lp->xb, *q3, i);
    }

    mpq_clear(q1);
    mpq_clear(q2);
    vector_free(&v1);
    vector_free(&v2);
  }
}

void lp_set_rhs_positive(LP* lp) {
  mpq_t q;
  int   i;

  mpq_init(q);

  for (i = 0; i < lp->rows; i ++) {
    get_valid_rhs(&q, lp, i);
    if (mpq_sgn(q) < 0) {
      matrix_rev_row_sgn(lp->A, i);
      vector_rev_element_sgn(lp->b, i);
    }
  }

  mpq_clear(q);
}

void lp_add_artificials(LP* lp) {
  mpq_t  q1, q2, q3, best;
  char  buf[LP_NAME_LEN_MAX];
  int  i, j, k, l, m, x, best_var;

  mpq_init(q1);
  mpq_init(q2);
  mpq_init(q3);
  mpq_init(best);
  vector_zero_clear(lp->c);
  vector_zero_clear(lp->cb);

  x = lp_select_basis(lp);

  for (i = 0; i < lp->rows; i ++) {
    if (lp->basis_column[i] != -1)
      continue;
    get_valid_rhs(&q1, lp, i);
    if (mpq_sgn(q1) < 0) {
      matrix_rev_row_sgn(lp->A, i);
      vector_rev_element_sgn(lp->b, i);
      mpq_neg(q1, q1);
    }
    best_var = -1;
    for (m = 0; m < lp->A->row[i]->nonzeros; m ++) {
      j = lp->A->row[i]->i[m];
      if (!is_normal_var(lp, j) ||
          lp->is_basis[j] ||
          lp->A->column[j]->nonzeros == 0)
        continue;
      if (vector_get_first_nonzero_i(lp->A->column[j]) != i)
        continue;
      if (mpq_sgn(*vector_get_element_ptr(lp->x, j)))
        continue;
      mympq_div(q2, q1, *vector_get_element_ptr(lp->A->column[j], i));
      if ((lp->upper.is_valid[j] && mpq_cmp(lp->upper.bound[j], q2) < 0) ||
          (lp->lower.is_valid[j] && mpq_cmp(lp->lower.bound[j], q2) > 0))
        continue;
      if (mpq_sgn(q2)) {
        l = 0;
        for (k = 0; k < lp->A->column[j]->nonzeros; k ++) {
          if (lp->basis_column[lp->A->column[j]->i[k]] == -1 ||
              lp->A->column[j]->i[k] == i)
            continue;
          l = 1;
          break;
        }
        if (l)
          continue;
      }
      vector_get_element(&q3, lp->c_back, j);
      mympq_mul(q3, q3, q2);
      //if (!lp->maximize)
      //  mpq_neg(q3, q3);
      if (best_var >= 0 && mpq_cmp(best, q3) >= 0)
        continue;
      best_var = j;
      mpq_set(best, q3); 
    }
    if (best_var < 0)
      continue;
    x ++;
    j = best_var;
    mympq_div(q2, q1, *vector_get_element_ptr(lp->A->column[j], vector_get_first_nonzero_i(lp->A->column[j])));
    lp->is_basis[j] = TRUE;
    lp->basis_column[i] = j;
    vector_set_element(lp->x,  q2, j);
    vector_set_element(lp->xb, q2, i);
  }
  matrix_resize(lp->A, lp->rows, lp->vars + lp->rows - x);

  for (i = 0; i < lp->rows; i ++) {
    if (lp->basis_column[i] >= 0)
      continue;
    get_valid_rhs(&q1, lp, i);
    if(lp->var_name == NULL)
      j = lp_add_var_without_A(lp, NULL);
    else {
      strncpy(buf, "#A", LP_NAME_LEN_MAX);
      strncat(buf, lp->row_name[i], LP_NAME_LEN_MAX);
      j = lp_add_var_without_A(lp, buf);
    }
    lp_add_slackref(lp, -(i+1));
    lp->var_type[j] = LP_VAR_TYPE_ARTIFICIAL;
    lp->basis_column[i] = j;
    lp->is_basis[j] = TRUE;
    lp_set_coefficient(lp, mympq_one, i, j);
    vector_set_element(lp->x,  q1, j);
    vector_set_element(lp->xb, q1, i);
    vector_set_element(lp->c,  mympq_minus_one, j);
    vector_set_element(lp->cb, mympq_minus_one, i);
  }
#if 0
  lp->vars = MIN(lp->vars, lp->A->columns);
#endif

  mpq_clear(q1);
  mpq_clear(q2);
  mpq_clear(q3);
  mpq_clear(best);
}

int  lp_normal_vars(LP* lp) {
  int  i, s;

  s = 0;
  for (i = 0; i < lp->vars; i ++)
    s += is_normal_var(lp, i);
  return s;
}

int  lp_slack_vars(LP* lp) {
  int  i, s;

  s = 0;
  for (i = 0; i < lp->vars; i ++)
    s += is_slack_var(lp, i);
  return s;
}

int  lp_artificial_vars(LP* lp) {
  int  i, s;

  s = 0;
  for (i = 0; i < lp->vars; i ++)
    s += is_artificial_var(lp, i);
  return s;
}

void lp_print_basis(LP* lp) {
  int  i, var;
  mpq_t q;

  mpq_init(q);
  putchar('\n');
  for (i = 0; i < lp->rows; i ++) {
    var = lp->basis_column[i];
    printf("%s: ", lp->var_name[var]);
    vector_get_element(&q, lp->x, var);
    print_rational_as_float(q, 12);
    printf("\t(");
    if (lp->lower.is_valid[var])
      print_rational_as_float(lp->lower.bound[var], 12);
    else
      printf("  ---  ");
    printf("\t");
    if (lp->upper.is_valid[var])
      print_rational_as_float(lp->upper.bound[var], 12);
    else
      printf("  ---  ");
    printf(")");
    putchar('\n');
  }
  putchar('\n');
  mpq_clear(q);
}

void lp_print_nonzero_vars(LP* lp) {
  int  var;

  for (var = 0; var < lp->vars; var ++) {
    if (vector_element_is_zero(lp->x, var))
      continue;
    printf("%s: ", lp->var_name[var]);
    print_rational_as_float(*vector_get_element_ptr(lp->x, var), 12);
    //print_rational(lp->x->value[i]);
    printf("\t(");
    if (lp->lower.is_valid[var])
      print_rational_as_float(lp->lower.bound[var], 12);
    else
      printf("  ---  ");
    printf("\t");
    if (lp->upper.is_valid[var])
      print_rational_as_float(lp->upper.bound[var], 12);
    else
      printf("  ---  ");
    printf(")");
    if (lp->is_basis[var])
      printf("  basis");
    putchar('\n');
  }
  putchar('\n');
}
