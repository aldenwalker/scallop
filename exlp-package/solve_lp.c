#include "solve_lp.h"
#include <math.h>
#include <string.h>
#include <time.h>
#include "mylib.h"
#include "lpstruct.h"
#include "eta_file.h"
#include "preprocess.h"
#include "lu.h"


/*
#define rdtsc(ll) __asm__ __volatile__ ("rdtsc" : "=A" (ll))

my_clock_t my_clock(void) {
  my_clock_t  t;

  rdtsc(t);
  return t;
}
*/
#define my_clock() clock()


int set_basis_for_presolve(LP* lp) {
  int  r;

  if (lp->preprocess) {
    do {
      r = lp->rows;
      //check_parallel_columns(lp);
      if (check_const_vars(lp) == LP_RESULT_INFEASIBLE)
        return LP_RESULT_INFEASIBLE;
      /*
      if (lp_check_redundancy(lp) == LP_RESULT_INFEASIBLE)
        return LP_RESULT_INFEASIBLE;
      lp_arrange_inequality2(lp);
      check_trivial_constraints0(lp);
      if (check_trivial_constraints1(lp) == LP_RESULT_INFEASIBLE)
        return LP_RESULT_INFEASIBLE;
      if (check_trivial_constraints2(lp) == LP_RESULT_INFEASIBLE)
        return LP_RESULT_INFEASIBLE;
      check_trivial_constraints3(lp);
      check_tight_constraints1(lp);
      check_tight_constraints2(lp);
      if (check_parallel_constraints(lp) == LP_RESULT_INFEASIBLE)
        return LP_RESULT_INFEASIBLE;
      //fprintf(stderr,"  -_-;\n");
      */
    } while (r != lp->rows && 0);
  }
  lp_arrange_inequality(lp);
  lp_add_slacks(lp);
  lp_sol_init(lp);
  lp_set_rhs_positive(lp);
  lp_add_artificials(lp);
  // ↓new_eta_file()を後から呼ぶようにしたんでまずい. 
  //if (lp->preprocess)
  //  remove_artificials_from_basis(lp);
  //fprintf(stderr,"  ^_^;\n");
  if (lp->verbose) {
    printf("variables(after preprocessing): %d + %d(slack) + %d(artificial)\n",
           lp_normal_vars(lp), lp_slack_vars(lp), lp_artificial_vars(lp));
    printf("rows(after preprocessing): %d\n", lp->rows);
  }

  return 0;
}

void set_object_function(LP* lp) {
  mpq_t  q;
  int  i;

  mpq_init(q);
  vector_copy(lp->c, lp->c_back);

  for (i = 0; i < lp->rows; i ++) {
    vector_get_element(&q, lp->c, lp->basis_column[i]);
    vector_set_element(lp->cb, q, i);
  }

  for (i = 0; i < lp->vars; i ++) {
    if (is_artificial_var(lp, i)) {
      mpq_init(lp->upper.bound[i]);
      lp->upper.is_valid[i] = TRUE;
    }
  }

  /*
  remove_artificials_from_basis(lp);
  for (i = 0; i < lp->rows; i ++) {
    if (is_artificial_var(lp, lp->basis_column[i])) {
      lp_remove_row(lp, i);
      i --;
      printf("row %d is removed.\n", i);
    }
  }
  */

  mpq_clear(q);
}

void reinversion(LP* lp) {
  int  i;

  //printf("(%d %d)", lp->eta->k, lp->eta->s);fflush(stdout);
  lp->eta->k = 0;
  lp->eta->init_time = my_clock();
  lp_LU_basis(lp);
  if ((lp->pivot_rule & LP_PIVOT_STEEPEST_EDGE) ||
      (lp->pivot_rule & LP_PIVOT_PROJECTED_STEEPEST_EDGE) ||
      (lp->pivot_rule & LP_PIVOT_DEVEX)) {
    for (i = 0; i < lp->rows; i ++) {
      vector_get_d(lp->eta->U->column[i], lp->eta->U_d[i]);
      vector_get_d(lp->eta->L[i]->eta_vector, lp->eta->L_d[i]->eta_vector);
      lp->eta->L_d[i]->eta_column = lp->eta->L[i]->eta_column;
    }
  }
  for (i = 0; i < lp->eta->Q->dimension; i ++) {
    lp->eta->Q->row[i] = i;
    lp->eta->Q->column[i] = i;
    lp->eta->R->row[i] = i;
    lp->eta->R->column[i] = i;
  }
  lp->eta->s = 0;
  lp->eta->last_time = my_clock();
}

void basis_permutate(matrix* U, int* s, int* t,
                     permutation_matrix* Q, permutation_matrix* R) {
  /* U を Forrest & Tomlin のやり方で置換. *s, *t はそれぞれこぶの始めと終り.
     で, *s, *t には変換後のそれを返す. Q, R には置換行列を返す.
     *s == *t となっていれば上三角になったってことね. */
  int  i;

  for (i = 0; i < *s; i ++)
    Q->column[i] = i;
  Q->column[i++] = *t;
  for ( ; i <= *t; i ++)
    Q->column[i] = i-1;
  for ( ; i < Q->dimension; i ++)
    Q->column[i] = i;
  for (i = 0; i < Q->dimension; i ++)
    Q->row[Q->column[i]] = i;
  Q->s = *s;
  Q->t = *t;

  permutation_matrix_copy(R, Q);
  permutation_matrix_inv(R);

  matrix_permutate(Q, U, R);
}

void basis_permutate_reid(matrix* U, int* s, int* t,
                          permutation_matrix* Q, permutation_matrix* R) {
  /* U を Reid のやり方で置換. *s, *t はそれぞれこぶの始めと終り.
     で, *s, *t には変換後のそれを返す. Q, R には置換行列を返す.
     *s == *t となっていれば上三角になったってことね. */
  int  i, j;
  int  ss, tt;
  EXLPvector* v;
  EXLPvector_d* v_d;

  if (*s == *t) {
    Q->s = Q->t = R->s = R->t = *s;
    return;
  }

  ss = *s;
  tt = *t;

  v_d = U->row[ss];
  if (vector_d_get_nonzeros(v_d, 0, tt) == 1 &&
      !vector_d_element_is_zero(v_d, ss)) {
    Q->column[ss] = tt;
    for (j = ss+1; j <= tt; j ++)
      Q->column[j] = j-1;
    for (j = ss; j <= tt; j ++)
      Q->row[Q->column[j]] = j;
    Q->s = ss;
    Q->t = tt;
    permutation_matrix_copy(R, Q);
    permutation_matrix_inv(R);
    matrix_permutate(Q, U, R);
    *s = *t;
    return;
  }

  for (i = ss+1; i <= tt; i ++) {
    v = U->column[i];
    if (vector_get_nonzeros(v, ss, tt) == 1 &&
        !vector_element_is_zero(v, i)) {
      Q->column[i] = (*s)++;
      Q->row[i] = -1;
      continue;
    }
    v_d = U->row[i];
    if (vector_d_get_nonzeros(v_d, ss, tt) == 1 &&
        !vector_d_element_is_zero(v_d, i)) {
      Q->column[i] = (*t)--;
      Q->row[i] = -1;
      continue;
    }
  }
  j = *s;
  for (i = ss+1; i <= tt; i ++) {
    if (Q->row[i] < 0)
      continue;
    Q->column[i] = j ++;
  }
  Q->column[ss] = j;
  Q->s = ss;
  Q->t = tt;

  for (i = ss; i <= tt; i ++)
    Q->row[Q->column[i]] = i;

  permutation_matrix_copy(R, Q);
  permutation_matrix_inv(R);
  matrix_permutate(Q, U, R);
}

void eta_file_update(LP* lp, EXLPvector* w, int leaving_row) {
  permutation_matrix* Q;
  permutation_matrix* R;
  EXLPvector* v;
  EXLPvector_d** v_d;
  eta_singleton* es;
  int  i, j, k, s, t;
  mpq_t  q, r;
  mpq_t* p;
  //int* tmp;
  //char  tmp[sizeof(mpq_t)];
  my_clock_t  cl;

  cl = my_clock();

  //masashi

  if ((lp->eta->k >= ETA_MAX) ||
      (lp->reinversion_cycle == 0 &&
       cl-lp->eta->init_time <
       (cl-lp->eta->last_time)*(lp->eta->k+1)) ||
      (lp->reinversion_cycle != 0 &&
       lp->eta->k >= lp->reinversion_cycle)) {

    reinversion(lp);
    return;
  }

  leaving_row = lp->eta->R->column[leaving_row];
  t = vector_get_last_nonzero_i(w);

  Q = new_permutation_matrix(lp->eta->Q->dimension);
  R = new_permutation_matrix(lp->eta->R->dimension);

  matrix_set_column(lp->eta->U, w, leaving_row);
  //printf("%d  ", w->nonzeros);fflush(stdout);
  //printf("(%d,%d)", leaving_row, t);fflush(stdout);
  basis_permutate_reid(lp->eta->U, &leaving_row, &t, Q, R);
  //basis_permutate(lp->eta->U, &leaving_row, &t, Q, R);
  //printf("%d,%d  ", leaving_row, t);fflush(stdout);
  //printf("%d  ", t-leaving_row);fflush(stdout);

  mpq_init(q);
  mpq_init(r);

  s = 0;
  for (i = leaving_row; i < t; i ++) {
    v = lp->eta->U->column[i];
    if (vector_element_is_zero(v, t))
      continue;
    matrix_get_element(&q, lp->eta->U, t, i);
    mpq_neg(q, q);
    eta_singleton_set(lp->eta->Ls[lp->eta->s+s], q, t, i);

    for (j = 0; j < lp->eta->U->row[i]->nonzeros; j ++) {
      k = lp->eta->U->row[i]->i[j];
      matrix_get_element(&r, lp->eta->U, i, k);
      mympq_mul(r, r, q);
      matrix_add_element(lp->eta->U, r, t, k);
    }
    s ++;
  }
  matrix_get_element(&r, lp->eta->U, i, i);
  if (mpq_sgn(r)) {
    if (!mpq_equal(r, mympq_one)) {
      mpq_inv(q, r);
      eta_singleton_set(lp->eta->Ls[lp->eta->s+s], q, i, i);
      matrix_row_scalar_product(lp->eta->U, q, i);
      s ++;
    }
  } else {
    printf("something's wrong...\n");
  }

  permutation_matrix_mul(lp->eta->Q, R);
  permutation_matrix_copy(lp->eta->R, lp->eta->Q);
  permutation_matrix_inv(lp->eta->R);

  for (i = 0; i < s; i ++) {
    es = lp->eta->Ls[lp->eta->s+i];
    if (Q->s <= es->row && es->row <= t)
      es->row = lp->eta->Q->column[es->row];
    es->column = lp->eta->Q->column[es->column];
  }

  if ((lp->pivot_rule & LP_PIVOT_STEEPEST_EDGE) ||
      (lp->pivot_rule & LP_PIVOT_PROJECTED_STEEPEST_EDGE) ||
      (lp->pivot_rule & LP_PIVOT_DEVEX)) {
    for (i = 0; i < s; i ++) {
      eta_singleton_get_d(lp->eta->Ls  [lp->eta->s+i],
                          lp->eta->Ls_d[lp->eta->s+i]);
    }
    v_d = my_malloc(Q->dimension*sizeof(EXLPvector_d*));
    for (i = 0; i < Q->dimension; i ++) {
      vector_d_permutate(Q, lp->eta->U_d[i]);
      v_d[i] = lp->eta->U_d[i];
    }
    for (i = 0; i < Q->dimension; i ++)
      lp->eta->U_d[i] = v_d[Q->row[i]];
    vector_get_d(lp->eta->U->column[t], lp->eta->U_d[t]);
    for (i = leaving_row; i < Q->dimension; i ++) {
      p = matrix_get_element_ptr(lp->eta->U, t, i);
      vector_d_set_element(lp->eta->U_d[i], mpq_get_d(*p), t);
    }
    free(v_d);
  }

  lp->eta->s += s;

  permutation_matrix_free(Q);
  permutation_matrix_free(R);
  mpq_clear(q);
  mpq_clear(r);

  lp->eta->last_time = my_clock();
  lp->eta->k ++;
}

int select_entering_column_bland(LP* lp, EXLPvector* y, int* ya_l_c) {
  mpq_t  q1;
  mpq_t* q2;
  int  var, s;

  mpq_init(q1);

  for (var = 0; var < lp->vars; var ++) {
    if (lp->is_basis[var] ||
	//lp->A->column[var]->nonzeros <= 0 ||
        is_artificial_var(lp, var) || is_const_var(lp, var) ||
        is_free_var(lp, var))
      continue;
    vector_inner_product(&q1, y, lp->A->column[var]);
    q2 = vector_get_element_ptr(lp->c, var);
    s = mpq_cmp(q1, *q2);

    if (s < 0) {
      if (lp->upper.is_valid[var]) {
        q2 = vector_get_element_ptr(lp->x, var);
        if (mpq_cmp(*q2, lp->upper.bound[var]) >= 0)
          continue;
      }
      *ya_l_c = TRUE;
      break;
    } else if (s > 0) {
      if (lp->lower.is_valid[var]) {
        q2 = vector_get_element_ptr(lp->x, var);
        if (mpq_cmp(*q2, lp->lower.bound[var]) <= 0)
          continue;
      }
      *ya_l_c = FALSE;
      break;
    }
  }

  mpq_clear(q1);

  if (var == lp->vars)
    return -1;
  return var;
}

int select_entering_column_bland2(LP* lp, EXLPvector* y, int* ya_l_c) {
  mpq_t  q1;
  mpq_t* q2;
  int  var, s;
  static int start = 0;

  mpq_init(q1);

  var = start;
  do {
    var ++;
    if (var >= lp->vars)
      var = 0;

    if (lp->is_basis[var] ||
	//lp->A->column[var]->nonzeros <= 0 ||
        is_artificial_var(lp, var) || is_const_var(lp, var) ||
        is_free_var(lp, var))
      continue;
    vector_inner_product(&q1, y, lp->A->column[var]);
    q2 = vector_get_element_ptr(lp->c, var);
    s = mpq_cmp(q1, *q2);

    if (s < 0) {
      if (lp->upper.is_valid[var]) {
        q2 = vector_get_element_ptr(lp->x, var);
        if (mpq_cmp(*q2, lp->upper.bound[var]) >= 0)
          continue;
      }
      *ya_l_c = TRUE;
      break;
    } else if (s > 0) {
      if (lp->lower.is_valid[var]) {
        q2 = vector_get_element_ptr(lp->x, var);
        if (mpq_cmp(*q2, lp->lower.bound[var]) <= 0)
          continue;
      }
      *ya_l_c = FALSE;
      break;
    }
  } while (var != start);

  mpq_clear(q1);

  if (var == start)
    return -1;

  start = var + 1;
  if (start >= lp->vars)
    start = 0;

  return var;
}

int select_entering_column_dantzig(LP* lp, EXLPvector* y, int* ya_l_c) {
  mpq_t  q1, best;
  mpq_t* q2;
  mpq_t* q3;
  int  i, var;
  int  best_var;
  EXLPvector_d* y_d;
  double  q_d;
  int* t;
  double* d;

  mpq_init(q1);
  mpq_init(best);
  best_var = -1;
  y_d = new_vector_d(y->dimension);
  vector_get_d(y, y_d);

  t = my_malloc(lp->vars*sizeof(int));
  d = my_malloc(lp->vars*sizeof(double));

  for (var = 0; var < lp->vars; var ++) {
    if (lp->is_basis[var] ||
	//lp->A->column[var]->nonzeros <= 0 ||
        is_artificial_var(lp, var) || is_const_var(lp, var) ||
        is_free_var(lp, var)) {
      t[var] = -1;
      d[var] = -2;
      continue;
    }
    t[var] = var;
    q_d = vector_d_inner_product(y_d, lp->A_d[var]);
    q2 = vector_get_element_ptr(lp->c, var);
    q_d -= mpq_get_d(*q2);

    if (q_d < 0) {
      if (lp->upper.is_valid[var]) {
        q2 = vector_get_element_ptr(lp->x, var);
        if (mpq_cmp(*q2, lp->upper.bound[var]) >= 0) {
          d[var] = -1;
          continue;
	}
      }
      d[var] = -q_d;
    } else {
      if (lp->lower.is_valid[var]) {
        q2 = vector_get_element_ptr(lp->x, var);
        if (mpq_cmp(*q2, lp->lower.bound[var]) <= 0) {
          d[var] = -1;
          continue;
	}
      }
      d[var] = q_d;
    }
  }
  my_sort_d2(t, d, 0, lp->vars-1);

  for (i = 0; i < lp->vars; i ++) {
    var = t[i];
    if (var < 0)
      continue;
    if (best_var >= 0 && mpq_get_d(best)-.00001 > d[i])
      continue;

    vector_inner_product(&q1, y, lp->A->column[var]);
    q2 = vector_get_element_ptr(lp->c, var);
    mympq_sub(q1, q1, *q2);

    if (mpq_sgn(q1) < 0) {
      if (lp->upper.is_valid[var]) {
        q2 = vector_get_element_ptr(lp->x, var);
        if (mpq_cmp(*q2, lp->upper.bound[var]) >= 0)
          continue;
      }
      mpq_neg(q1, q1);
      if (best_var >= 0 && mpq_cmp(best, q1) > 0)
        continue;
      if (best_var >= 0 && mpq_equal(best, q1)) {
        q2 = vector_get_element_ptr(lp->c_back, best_var);
        q3 = vector_get_element_ptr(lp->c_back, var);
        if (mpq_cmp(*q2, *q3) >= 0)
          continue;
      }
      *ya_l_c = TRUE;
      mpq_set(best, q1);
      best_var = var;
    } else if (mpq_sgn(q1) > 0) {
      if (lp->lower.is_valid[var]) {
        q2 = vector_get_element_ptr(lp->x, var);
        if (mpq_cmp(*q2, lp->lower.bound[var]) <= 0)
          continue;
      }
      if (best_var >= 0 && mpq_cmp(best, q1) > 0)
        continue;
      if (best_var >= 0 && mpq_equal(best, q1)) {
        q2 = vector_get_element_ptr(lp->c_back, best_var);
        q3 = vector_get_element_ptr(lp->c_back, var);
        if (mpq_cmp(*q2, *q3) >= 0)
          continue;
      }
      *ya_l_c = FALSE;
      mpq_set(best, q1);
      best_var = var;
    }
  }

  mpq_clear(q1);
  mpq_clear(best);
  vector_d_free(&y_d);
  free(t);
  free(d);

  return best_var;
}

void devex_rectify(LP* lp) {
  int  i;

  for (i = 0; i < lp->vars; i ++) {
    lp->devex_weight[i] = 1.0;
    lp->framework[i] = !lp->is_basis[i];
  }
}

void devex_initialize(LP* lp) {
  lp->devex_weight
    = (double*)my_realloc(lp->devex_weight, lp->vars*sizeof(double));
  lp->framework
    = (int*)   my_realloc(lp->framework, lp->vars*sizeof(int));
  devex_rectify(lp);
}

void projected_rectify(LP* lp) {
  int  i;

  for (i = 0; i < lp->vars; i ++) {
    lp->steepest_edge_table[i] = 1.0;
    lp->framework[i] = !lp->is_basis[i];
  }
}

void projected_initialize(LP* lp) {
  lp->steepest_edge_table
    = (double*)my_realloc(lp->steepest_edge_table, lp->vars*sizeof(double));
  lp->framework
    = (int*)   my_realloc(lp->framework, lp->vars*sizeof(int));
  projected_rectify(lp);
}

int devex_test(LP* lp, EXLPvector* d, int e_column) {
  mpq_t* q;
  double  s, t;
  int  i;
  static int  c = 0;

  if (c ++ >= 20) {
    c = 0;
    return TRUE;
  }

  s = 0;
  for (i = 0; i < lp->rows; i ++) {
    if (!lp->framework[lp->basis_column[i]])
      continue;
    q = vector_get_element_ptr(d, i);
    t = mpq_get_d(*q);
    s += t*t;
  }

  t = lp->devex_weight[e_column]-1;

  if (4*t*t < s) {
    c = 0;
    return TRUE;
  }
  return FALSE;
}

void devex_update_weights(LP* lp, EXLPvector* d, int e_col, int l_row) {
  EXLPvector_d* w_d;
  EXLPvector_d* s_d;
  double r, t;
  double a_q, a_j;
  mpq_t* q;
  int  i, var;

  w_d = new_vector_d(lp->rows);
  s_d = new_vector_d(lp->rows);

  r = 1.0;
  for (i = 0; i < d->nonzeros; i ++) {
    var = lp->basis_column[d->i[i]];
    if (!lp->framework[var])
      continue;
    t = mpq_get_d(d->value[i]);
    r += t*t;
  }
  r = sqrt(r);
  vector_d_set_element(w_d, 1.0, l_row);
  eta_file_btran_d(lp->eta, w_d, s_d);
  q = vector_get_element_ptr(d, l_row);
  a_q = mpq_get_d(*q);

  lp->devex_weight[lp->basis_column[l_row]] = MAX(1.0, r/a_q);
  for (var = 0; var < lp->vars; var ++) {
    if (lp->is_basis[var] || var == e_col)
      continue;
    a_j = vector_d_inner_product(s_d, lp->A_d[var]) / a_q;
    a_j /= a_q;
    a_j = ABS(a_j);
    lp->devex_weight[var] = MAX(lp->devex_weight[var], a_j * r);
  }

  vector_d_free(&w_d);
  vector_d_free(&s_d);
}

int select_entering_column_devex(LP* lp, EXLPvector* y, int* ya_l_c) {
  int  var;
  mpq_t  q1;
  mpq_t* q2;
  double  tmp, best;
  int  best_var;
  EXLPvector_d* y_d;
  int  i;
  int* t;
  double* d;

  mpq_init(q1);
  best_var = -1;
  best = 0;
  y_d = new_vector_d(y->dimension);
  vector_get_d(y, y_d);

  t = my_malloc(lp->vars*sizeof(int));
  d = my_malloc(lp->vars*sizeof(double));

  for (var = 0; var < lp->vars; var ++) {
    if (lp->is_basis[var] ||
	//lp->A->column[var]->nonzeros <= 0 ||
        is_artificial_var(lp, var) || is_const_var(lp, var) ||
        is_free_var(lp, var)) {
      t[var] = -1;
      d[var] = -2;
      continue;
    }
    t[var] = var;
    tmp = vector_d_inner_product(y_d, lp->A_d[var]);
    q2 = vector_get_element_ptr(lp->c, var);
    tmp -= mpq_get_d(*q2);

    if (tmp < 0)
      d[var] = -tmp / lp->devex_weight[var];
    else
      d[var] = tmp / lp->devex_weight[var];
  }
  my_sort_d2(t, d, 0, lp->vars-1);

  for (i = 0; i < lp->vars; i ++) {
    var = t[i];
    if (var < 0)
      continue;
    if (best_var >= 0 && best-.00001 > d[i])
      continue;

    vector_inner_product(&q1, y, lp->A->column[var]);
    q2 = vector_get_element_ptr(lp->c, var);
    mympq_sub(q1, q1, *q2);

    if (mpq_sgn(q1) < 0) {
      if (lp->upper.is_valid[var]) {
        q2 = vector_get_element_ptr(lp->x, var);
        if (mpq_cmp(*q2, lp->upper.bound[var]) >= 0)
          continue;
      }
      tmp = - mpq_get_d(q1) / lp->devex_weight[var];
      if (best_var >= 0 && best >= tmp)
        continue;

      *ya_l_c = TRUE;
      best = tmp;
      best_var = var;

    } else if (mpq_sgn(q1) > 0) {
      if (lp->lower.is_valid[var]) {
        q2 = vector_get_element_ptr(lp->x, var);
        if (mpq_cmp(*q2, lp->lower.bound[var]) <= 0)
          continue;
      }
      tmp = mpq_get_d(q1) / lp->devex_weight[var];
      if (best_var >= 0 && best >= tmp)
        continue;
      *ya_l_c = FALSE;
      best = tmp;
      best_var = var;
    }
  }

  mpq_clear(q1);
  vector_d_free(&y_d);
  free(t);
  free(d);

  return best_var;
}

void steepest_edge_initialize(LP* lp) {
  int  var;
  int  i;
  EXLPvector_d* d;
  double  a, b;

  d = new_vector_d(lp->rows);
  lp->steepest_edge_table =
    my_realloc(lp->steepest_edge_table, lp->vars*sizeof(double));

  for (var = 0; var < lp->vars; var ++) {
    if (lp->is_basis[var] ||
        is_artificial_var(lp, var) || is_const_var(lp, var) ||
        is_free_var(lp, var))
      continue;

    eta_file_ftran_d(lp->eta, lp->A_d[var], d);
    a = 1;
    for (i = 0; i < d->nonzeros; i ++) {
      b = d->value[i];
      a += b*b;
    }
    lp->steepest_edge_table[var] = a;
  }

  vector_d_free(&d);
}

int steepest_edge_update_table(LP* lp, EXLPvector* d, int e_col, int l_row) {
/* 誤差がたまって初期化が必要と判断したら非ゼロの値を返す */
  EXLPvector_d* w_d;
  EXLPvector_d* s_d;
  EXLPvector_d* v_d;
  double  a, a_q, g_q, q_d;
  int  var, ret;

  ret = 0;

  w_d = new_vector_d(lp->rows);
  s_d = new_vector_d(lp->rows);
  v_d = new_vector_d(lp->rows);

  vector_get_d(d, w_d);
  a_q = vector_d_get_element(w_d, l_row);
  g_q = vector_d_inner_product(w_d, w_d)+1;
  lp->steepest_edge_table[lp->basis_column[l_row]] = g_q/a_q/a_q;

  vector_d_set_element(v_d, 1.0, l_row);
  eta_file_btran_d(lp->eta, v_d, s_d);
  eta_file_btran_d(lp->eta, w_d, v_d);

  for (var = 0; var < lp->vars; var ++) {
    if (var == e_col || lp->is_basis[var] ||
        is_const_var(lp, var) || is_artificial_var(lp, var) ||
        is_free_var(lp, var))
      continue;
    q_d = vector_d_inner_product(s_d, lp->A_d[var]);
    if (q_d == 0)
      continue;
    a = q_d / a_q;
    q_d = vector_d_inner_product(lp->A_d[var], v_d);
    lp->steepest_edge_table[var] += a*(a*g_q - 2*q_d);
    if (lp->steepest_edge_table[var] <= 0.000001) {
      lp->steepest_edge_table[var] = 0.000001;
      ret = 1;
    }
  }

  vector_d_free(&w_d);
  vector_d_free(&v_d);
  vector_d_free(&s_d);

  return ret;
}

int projected_update_table(LP* lp, EXLPvector* d, int e_col, int l_row) {
/* 誤差がたまって初期化が必要と判断したら非ゼロの値を返す */
  EXLPvector_d* w_d;
  EXLPvector_d* s_d;
  EXLPvector_d* v_d;
  double  a, a_q, g_q, q_d;
  int  i, var, ret;

  ret = 0;

  w_d = new_vector_d(lp->rows);
  s_d = new_vector_d(lp->rows);
  v_d = new_vector_d(lp->rows);

  //vector_get_d(d, w_d);
  for (i = 0; i < d->nonzeros; i ++)
    if (lp->framework[lp->basis_column[d->i[i]]])
      vector_d_set_element(w_d, mpq_get_d(d->value[i]), d->i[i]);

  //a_q = vector_d_get_element(w_d, l_row);
  a_q = mpq_get_d(*vector_get_element_ptr(d, l_row));
  g_q = vector_d_inner_product(w_d, w_d)+lp->framework[e_col];
  lp->steepest_edge_table[lp->basis_column[l_row]]
    = (g_q /*+ lp->framework[lp->basis_column[l_row]] - lp->framework[e_col]*/)
      /a_q /a_q;

  vector_d_set_element(v_d, 1.0, l_row);
  eta_file_btran_d(lp->eta, v_d, s_d);
  eta_file_btran_d(lp->eta, w_d, v_d);

  for (var = 0; var < lp->vars; var ++) {
    if (var == e_col || lp->is_basis[var] ||
        is_const_var(lp, var) || is_artificial_var(lp, var) ||
        is_free_var(lp, var))
      continue;
    q_d = vector_d_inner_product(s_d, lp->A_d[var]);
    if (q_d == 0)
      continue;
    a = q_d / a_q;
    q_d = vector_d_inner_product(lp->A_d[var], v_d);
    lp->steepest_edge_table[var] += a*(a*(g_q /*+ lp->framework[lp->basis_column[l_row]] - lp->framework[e_col]*/) - 2*q_d);
    /*
    if (lp->steepest_edge_table[var] <= 0.000001) {
      lp->steepest_edge_table[var] = 0.000001;
      ret = 1;
    }
    */
  }

  vector_d_free(&w_d);
  vector_d_free(&v_d);
  vector_d_free(&s_d);

  return ret;
}

void steepest_edge_initialize_dual(LP* lp) {
  int  i;
  EXLPvector_d* e;
  EXLPvector_d* r;

  e = new_vector_d(lp->rows);
  r = new_vector_d(lp->rows);
  lp->steepest_edge_table =
    my_realloc(lp->steepest_edge_table, lp->rows*sizeof(double));

  for (i = 0; i < lp->rows; i ++) {
    vector_d_zero_clear(e);
    vector_d_set_element(e, 1, i);
    eta_file_btran_d(lp->eta, e, r);
    lp->steepest_edge_table[i] = vector_d_inner_product(r, r);
  }

  vector_d_free(&e);
  vector_d_free(&r);
}

int steepest_edge_update_table_dual(LP* lp, EXLPvector* d, int e_col, int l_row) {
/* 誤差がたまって初期化が必要と判断したら非ゼロの値を返す */
  EXLPvector_d* w_d;
  EXLPvector_d* s_d;
  EXLPvector_d* t_d;
  double  w_p, w_i, b_p, a;
  int  i, ret;

  ret = 0;

  w_d = new_vector_d(lp->rows);
  s_d = new_vector_d(lp->rows);
  t_d = new_vector_d(lp->rows);

  vector_get_d(d, w_d);
  w_p = vector_d_get_element(w_d, l_row);

  vector_d_set_element(t_d, 1.0, l_row);
  eta_file_btran_d(lp->eta, t_d, s_d);
  eta_file_ftran_d(lp->eta, s_d, t_d);

  b_p = vector_d_inner_product(s_d, s_d);
  if (lp->steepest_edge_table[l_row] < b_p-0.1 ||
      lp->steepest_edge_table[l_row] > b_p+0.1)
    ret = 1;
  //b_p = lp->steepest_edge_table[l_row];
  //printf("%lf %lf\n", vector_d_inner_product(s_d,s_d),b_p);
  lp->steepest_edge_table[l_row] = b_p/w_p/w_p;

  for (i = 0; i < lp->rows; i ++) {
    if (i == l_row)
      continue;
    w_i = vector_d_get_element(w_d, i);
    a = w_i/w_p;
    if (a == 0)
      continue;
    lp->steepest_edge_table[i] += a*(a*b_p-2*vector_d_get_element(t_d, i));
    if (lp->steepest_edge_table[i] <= 0.000001) {
      lp->steepest_edge_table[i] = 0.000001;
      ret = 1;
    }
  }

  vector_d_free(&w_d);
  vector_d_free(&s_d);
  vector_d_free(&t_d);

  return ret;
}

int select_entering_column_steepest_edge(LP* lp, EXLPvector* y, int* ya_l_c) {
  mpq_t  q1;
  mpq_t* q2;
  int  var;
  int  best_var;
  int  i, s;
  EXLPvector_d* y_d;
  double  a;
  int* T;
  double* D;

  mpq_init(q1);
  y_d = new_vector_d(lp->rows);
  vector_get_d(y, y_d);

  best_var = -1;

  T = my_malloc(lp->vars*sizeof(int));
  D = my_malloc(lp->vars*sizeof(double));

  for (var = 0; var < lp->vars; var ++) {
    if (lp->is_basis[var] ||
	//lp->A->column[var]->nonzeros <= 0 ||
        is_artificial_var(lp, var) || is_const_var(lp, var) ||
        is_free_var(lp, var)) {
      T[var] = -1;
      D[var] = -1000000000000.0;
      continue;
    }
    T[var] = var;
    a = vector_d_inner_product(y_d, lp->A_d[var]);
    //vector_inner_product(&q1,y,lp->A->column[var]);
    //a=mpq_get_d(q1);
    a -= vector_d_get_element(lp->c_d, var);
    D[var] = a*a/lp->steepest_edge_table[var];
    if (a < -0.00001) {
      if (lp->upper.is_valid[var]) {
        q2 = vector_get_element_ptr(lp->x, var);
        if (mpq_cmp(*q2, lp->upper.bound[var]) >= 0)
          D[var] -= 100000000;
      }
    }
    if (a > 0.00001) {
      if (lp->lower.is_valid[var]) {
        q2 = vector_get_element_ptr(lp->x, var);
        if (mpq_cmp(*q2, lp->lower.bound[var]) <= 0)
          D[var] -= 100000000;
      }
    }
  }
  my_sort_d2(T, D, 0, lp->vars-1);

  for (i = 0; i < lp->vars; i ++) {
    var = T[i];
    if (var < 0)
      continue;

    vector_inner_product(&q1, y, lp->A->column[var]);
    q2 = vector_get_element_ptr(lp->c, var);
    s = mpq_cmp(q1, *q2);

    if (s < 0) {
      if (lp->upper.is_valid[var]) {
        q2 = vector_get_element_ptr(lp->x, var);
        if (mpq_cmp(*q2, lp->upper.bound[var]) >= 0)
          continue;
      }
      *ya_l_c = TRUE;
      best_var = var;
      break;
    } else if (s > 0) {
      if (lp->lower.is_valid[var]) {
        q2 = vector_get_element_ptr(lp->x, var);
        if (mpq_cmp(*q2, lp->lower.bound[var]) <= 0)
          continue;
      }
      *ya_l_c = FALSE;
      best_var = var;
      break;
    }
  }

  mpq_clear(q1);
  vector_d_free(&y_d);
  free(T);
  free(D);

  return best_var;
}

int select_entering_column_part(LP* lp, EXLPvector* y, int* ya_l_c,
                                int start, int end) {
  /*
  mpq_t  q1;
  mpq_t* q2;
  int  var;

  mpq_init(q1);

  for (var = start; var < end; var ++) {
    if (lp->is_basis[var] ||
	//lp->A->column[var]->nonzeros <= 0 ||
        is_artificial_var(lp, var) || is_const_var(lp, var) ||
        is_free_var(lp, var))
      continue;
    vector_inner_product(&q1, y, lp->A->column[var]);
    q2 = vector_get_element_ptr(lp->c, var);
    mympq_sub(q1, q1, *q2);

    if (mpq_sgn(q1) < 0) {
      if (lp->upper.is_valid[var]) {
        q2 = vector_get_element_ptr(lp->x, var);
        if (mpq_cmp(*q2, lp->upper.bound[var]) >= 0)
          continue;
      }
      *ya_l_c = TRUE;
      break;
    } else if (mpq_sgn(q1) > 0) {
      if (lp->lower.is_valid[var]) {
        q2 = vector_get_element_ptr(lp->x, var);
        if (mpq_cmp(*q2, lp->lower.bound[var]) <= 0)
          continue;
      }
      *ya_l_c = FALSE;
      break;
    }
  }

  mpq_clear(q1);

  if (var == lp->vars)
    return -1;
  return var;
  */

  mpq_t  q1;
  mpq_t* q2;
  int  var;
  int  best_var;
  int  i, s;
  EXLPvector_d* y_d;
  double  a;
  int* T;
  double* D;

  mpq_init(q1);
  y_d = new_vector_d(lp->rows);
  vector_get_d(y, y_d);

  best_var = -1;

  T = my_malloc(lp->vars*sizeof(int));
  D = my_malloc(lp->vars*sizeof(double));

  for (var = start; var < end; var ++) {
    if (lp->is_basis[var] ||
	//lp->A->column[var]->nonzeros <= 0 ||
        is_artificial_var(lp, var) || is_const_var(lp, var) ||
        is_free_var(lp, var)) {
      T[var] = -1;
      D[var] = -1000000000000.0;
      continue;
    }
    T[var] = var;
    a = vector_d_inner_product(y_d, lp->A_d[var]);
    q2 = vector_get_element_ptr(lp->c, var);
    a -= mpq_get_d(*q2);
    D[var] = a*a/lp->steepest_edge_table[var];
    if (a < -0.00001) {
      if (lp->upper.is_valid[var]) {
        q2 = vector_get_element_ptr(lp->x, var);
        if (mpq_cmp(*q2, lp->upper.bound[var]) >= 0)
          D[var] -= 100000000;
      }
    }
    if (a > 0.00001) {
      if (lp->lower.is_valid[var]) {
        q2 = vector_get_element_ptr(lp->x, var);
        if (mpq_cmp(*q2, lp->lower.bound[var]) <= 0)
          D[var] -= 100000000;
      }
    }
  }
  my_sort_d2(T, D, start, end-1);

  for (i = start; i < end; i ++) {
    var = T[i];
    if (var < 0)
      continue;

    vector_inner_product(&q1, y, lp->A->column[var]);
    q2 = vector_get_element_ptr(lp->c, var);
    s = mpq_cmp(q1, *q2);

    if (s < 0) {
      if (lp->upper.is_valid[var]) {
        q2 = vector_get_element_ptr(lp->x, var);
        if (mpq_cmp(*q2, lp->upper.bound[var]) >= 0)
          continue;
      }
      *ya_l_c = TRUE;
      best_var = var;
      break;
    } else if (s > 0) {
      if (lp->lower.is_valid[var]) {
        q2 = vector_get_element_ptr(lp->x, var);
        if (mpq_cmp(*q2, lp->lower.bound[var]) <= 0)
          continue;
      }
      *ya_l_c = FALSE;
      best_var = var;
      break;
    }
  }

  mpq_clear(q1);
  vector_d_free(&y_d);
  free(T);
  free(D);

  return best_var;

}

int check_movability(LP* lp, int ya_l_c, EXLPvector_d* d, int entering) {
  /* 退化してるときに無理矢理入る変数を探すときに使う。*/
  /* 非ゼロを返せば動ける可能性が大きい。              */
  int  i, var, leaving;

  leaving = -1;

  if (ya_l_c) { /* step 正 */

    for (i = 0; i < d->nonzeros; i ++) {
      var = lp->basis_column[d->i[i]];
      if (is_free_var(lp, var))
        continue;
      if (d->value[i] > 0) {
        if (lp->lower.is_valid[var] == FALSE)
          continue;
        if (mpq_equal(*vector_get_element_ptr(lp->xb, d->i[i]),
                      lp->lower.bound[var]))
          return 0;
      } else {
        if (lp->upper.is_valid[var] == FALSE)
          continue;
        if (mpq_equal(*vector_get_element_ptr(lp->xb, d->i[i]),
                      lp->upper.bound[var]))
          return 0;
      }
    }
  } else { /* step 負 */

    for (i = 0; i < d->nonzeros; i ++) {
      var = lp->basis_column[d->i[i]];
      if (is_free_var(lp, var))
        continue;
      if (d->value[i] > 0) {
        if (lp->upper.is_valid[var] == FALSE)
          continue;
        if (mpq_equal(*vector_get_element_ptr(lp->xb, d->i[i]),
                      lp->upper.bound[var]))
          return 0;
      } else {
        if (lp->lower.is_valid[var] == FALSE)
          continue;
        if (mpq_equal(*vector_get_element_ptr(lp->xb, d->i[i]),
                      lp->lower.bound[var]))
          return 0;
      }
    }
  }
  return 1;
}

int select_entering_column_dgn(LP* lp, EXLPvector* y, int* ya_l_c) {
  mpq_t  q1;
  mpq_t* q2;
  int  var;
  int  best_var;
  int  i, s;
  EXLPvector_d* d;
  EXLPvector_d* y_d;
  double  a;
  int* T;
  double* D;

  mpq_init(q1);
  d = new_vector_d(lp->rows);
  y_d = new_vector_d(lp->rows);
  vector_get_d(y, y_d);

  best_var = -1;

  T = my_malloc(lp->vars*sizeof(int));
  D = my_malloc(lp->vars*sizeof(double));

  for (var = 0; var < lp->vars; var ++) {
    if (lp->is_basis[var] ||
	//lp->A->column[var]->nonzeros <= 0 ||
        is_artificial_var(lp, var) || is_const_var(lp, var) ||
        is_free_var(lp, var)) {
      T[var] = -1;
      D[var] = -1000000000000.0;
      continue;
    }
    T[var] = var;
    a = vector_d_inner_product(y_d, lp->A_d[var]);
    q2 = vector_get_element_ptr(lp->c, var);
    a -= mpq_get_d(*q2);
    D[var] = a*a/lp->steepest_edge_table[var];
    if (a < -0.00001) {
      if (lp->upper.is_valid[var]) {
        q2 = vector_get_element_ptr(lp->x, var);
        if (mpq_cmp(*q2, lp->upper.bound[var]) >= 0)
          D[var] -= 100000000;
      }
    }
    if (a > 0.00001) {
      if (lp->lower.is_valid[var]) {
        q2 = vector_get_element_ptr(lp->x, var);
        if (mpq_cmp(*q2, lp->lower.bound[var]) <= 0)
          D[var] -= 100000000;
      }
    }
  }
  my_sort_d2(T, D, 0, lp->vars-1);

  for (i = 0; i < lp->vars; i ++) {
    var = T[i];
    if (var < 0)
      continue;

    vector_inner_product(&q1, y, lp->A->column[var]);
    q2 = vector_get_element_ptr(lp->c, var);
    s = mpq_cmp(q1, *q2);

    if (s < 0) {
      if (lp->upper.is_valid[var]) {
        q2 = vector_get_element_ptr(lp->x, var);
        if (mpq_cmp(*q2, lp->upper.bound[var]) >= 0)
          continue;
      }

      eta_file_ftran_d(lp->eta, lp->A_d[var], d);
      if (!check_movability(lp, TRUE, d, var))
        continue;

      *ya_l_c = TRUE;
      best_var = var;
      break;
    } else if (s > 0) {
      if (lp->lower.is_valid[var]) {
        q2 = vector_get_element_ptr(lp->x, var);
        if (mpq_cmp(*q2, lp->lower.bound[var]) <= 0)
          continue;
      }

      eta_file_ftran_d(lp->eta, lp->A_d[var], d);
      if (!check_movability(lp, FALSE, d, var))
        continue;

      *ya_l_c = FALSE;
      best_var = var;
      break;
    }
  }

  mpq_clear(q1);
  vector_d_free(&d);
  vector_d_free(&y_d);
  free(T);
  free(D);

  //if (best_var < 0)
  //  return select_entering_column_steepest_edge(lp, y, ya_l_c);
  return best_var;
}

int select_entering_column_free_variables(LP* lp, EXLPvector* y, int* ya_l_c) {
  /* 自由変数は一度基底に入ると二度と基底から出ない. そこで自由変数は
     最優先で基底に入れてしまう.
     負の値が返れば基底に入るような自由変数が無かったということ.      */
  mpq_t  q1;
  mpq_t* q2;
  int  var, s;

  mpq_init(q1);

  for (var = 0; var < lp->vars; var ++) {
    if (lp->is_basis[var] ||
	//lp->A->column[var]->nonzeros <= 0 ||
        is_artificial_var(lp, var) || is_const_var(lp, var) ||
        !is_free_var(lp, var))
      continue;
    vector_inner_product(&q1, y, lp->A->column[var]);
    q2 = vector_get_element_ptr(lp->c, var);
    s = mpq_cmp(q1, *q2);

    if (s < 0) {
      *ya_l_c = TRUE;
      mpq_clear(q1);
      return var;
    } else if (s > 0) {
      *ya_l_c = FALSE;
      mpq_clear(q1);
      return var;
    }
  }

  mpq_clear(q1);
  return -1;
}

int select_entering_column(LP* lp, EXLPvector* y, int* ya_l_c) {
  static int  count = 0;
  //static int  part = 0;
  static mpq_t  q1;
  //mpq_t  q2;
  int  c;

  if (count == 0)
    mpq_init(q1);


  /*
  count ++;if (count == 1){putchar('\n');fflush(stdout);}
  mpq_init(q2);
  vector_inner_product(&q2, lp->c, lp->x);
  if (mpq_equal(q1, q2)) {
    printf("o");fflush(stdout);
  } else {
    mpq_set(q1, q2);
  }
  mpq_clear(q2);
  */


  c = select_entering_column_free_variables(lp, y, ya_l_c);
  if (c >= 0)
    return c;

  /*
  if ((lp->pivot_rule & (LP_PIVOT_STEEPEST_EDGE|LP_PIVOT_PROJECTED_STEEPEST_EDGE)) &&
      ((count++) % (lp->rows/16+10) == 0)) {
    mpq_init(q2);
    vector_inner_product(&q2, lp->c, lp->x);
    if (mpq_equal(q1, q2)) {
      mpq_clear(q2);
      c = select_entering_column_dgn(lp, y, ya_l_c);
      if (c >= 0)
        return c;
    }
    mpq_set(q1, q2);
    mpq_clear(q2);
  }
  */
  /*
  if ((lp->pivot_rule & (LP_PIVOT_STEEPEST_EDGE | LP_PIVOT_PROJECTED_STEEPEST_EDGE)) &&
      //((count++) > (lp->rows/16+10))) {
      //((count++) > 10)) {
      ((count++) > 100)) {
    mpq_init(q2);
    vector_inner_product(&q2, lp->c, lp->x);
    if (mpq_equal(q1, q2)) {
      mpq_clear(q2);
      return select_entering_column_bland(lp, y, ya_l_c);
    }
    mpq_set(q1, q2);
    mpq_clear(q2);
    count = 1;
  }
  */
  /*
  if ((lp->pivot_rule & LP_PIVOT_MIX_HUERISTIC) &&
      (count++ > 10)) {
    mpq_init(q2);
    vector_inner_product(&q2, lp->c, lp->x);
    if (mpq_equal(q1, q2)) {
      int  e_c;
      mpq_clear(q2);
      //for ( ; part < 4; ) {
	//printf("part %d\n",part);
        e_c = select_entering_column_part(lp, y, ya_l_c,
                                          part*lp->vars/4,
                                          (part+1)*lp->vars/4);
        if (e_c >= 0)
	  return e_c;
        part ++;
      //}
      if (part == 4)
      part = 0;
      e_c = select_entering_column_dgn(lp, y, ya_l_c);
      if (e_c >= 0)
        return e_c;
      //lp->pivot_rule = LP_PIVOT_MIX_HUERISTIC|LP_PIVOT_BLAND;
    } else {
      //if (lp->pivot_rule & LP_PIVOT_BLAND)
      //  steepest_edge_initialize(lp);
      //lp->pivot_rule = LP_PIVOT_MIX_HUERISTIC|LP_PIVOT_STEEPEST_EDGE;
      count = 1;
      //part = 0;
      mpq_set(q1, q2);
      mpq_clear(q2);
    }
  }
  */

  if (lp->pivot_rule & LP_PIVOT_BLAND)
    return select_entering_column_bland(lp, y, ya_l_c);
  if (lp->pivot_rule & LP_PIVOT_BLAND2)
    return select_entering_column_bland2(lp, y, ya_l_c);
  if (lp->pivot_rule & LP_PIVOT_DANTZIG)
    return select_entering_column_dantzig(lp, y, ya_l_c);
  if (lp->pivot_rule & LP_PIVOT_DEVEX)
    return select_entering_column_devex(lp, y, ya_l_c);
  else // projected でもこれ↓が使えることに注意.
    return select_entering_column_steepest_edge(lp, y, ya_l_c);
}

/*  x[entering] <- x[entering] + step    */
/*  xb          <- xb          - step d  */
int select_leaving_row(LP* lp, mpq_t* step, int ya_l_c,
                       EXLPvector* d, int entering) {
  /* *step には ya_l_c == TRUE なら正, ya_l_c == FALSE なら負の値を返す. */
  /* 戻り値は有界じゃないなら 負, 有界なら 基底の中で出る行.             */
  /* ただし lp->rows がかえったら出る変数無し.                           */
  mpq_t  q;
  int  i, var, leaving;

  mpq_init(q);
  leaving = -1;

  if (ya_l_c) { /* step 正 */
    if (lp->upper.is_valid[entering]) {
      vector_get_element(step, lp->x, entering);
      mympq_sub(*step, lp->upper.bound[entering], *step);
      leaving = lp->rows;
    }
    for (i = 0; i < d->nonzeros; i ++) {
      var = lp->basis_column[d->i[i]];
      if (is_free_var(lp, var))
        continue;
      if (mpq_sgn(d->value[i]) > 0) {
        if (lp->lower.is_valid[var] == FALSE)
          continue;
        vector_get_element(&q, lp->xb, d->i[i]);
        mympq_sub(q, q, lp->lower.bound[var]);
        mympq_div(q, q, d->value[i]);
        if (leaving < 0 || mpq_cmp(q, *step) < 0 ||
	    (mpq_equal(q, *step) && (
             (leaving != lp->rows &&
              lp->basis_column[d->i[i]] < lp->basis_column[leaving]) ||
             (leaving == lp->rows &&
              lp->basis_column[d->i[i]] < entering)))) {
          mpq_set(*step, q);
          leaving = d->i[i];
	}
      } else {
        if (lp->upper.is_valid[var] == FALSE)
          continue;
        vector_get_element(&q, lp->xb, d->i[i]);
        mympq_sub(q, q, lp->upper.bound[var]);
        mympq_div(q, q, d->value[i]);
        if (leaving < 0 || mpq_cmp(q, *step) < 0 ||
	    (mpq_equal(q, *step) && (
             (leaving != lp->rows &&
              lp->basis_column[d->i[i]] < lp->basis_column[leaving]) ||
             (leaving == lp->rows &&
              lp->basis_column[d->i[i]] < entering)))) {
          mpq_set(*step, q);
          leaving = d->i[i];
	}
      }
    }
  } else { /* step 負 */
    if (lp->lower.is_valid[entering]) {
      vector_get_element(step, lp->x, entering);
      mympq_sub(*step, lp->lower.bound[entering], *step);
      leaving = lp->rows;
    }
    for (i = 0; i < d->nonzeros; i ++) {
      var = lp->basis_column[d->i[i]];
      if (is_free_var(lp, var))
        continue;
      if (mpq_sgn(d->value[i]) > 0) {
        if (lp->upper.is_valid[var] == FALSE)
          continue;
        vector_get_element(&q, lp->xb, d->i[i]);
        mympq_sub(q, q, lp->upper.bound[var]);
        mympq_div(q, q, d->value[i]);
        if (leaving < 0 || mpq_cmp(q, *step) > 0 ||
	    (mpq_equal(q, *step) && (
             (leaving != lp->rows &&
              lp->basis_column[d->i[i]] < lp->basis_column[leaving]) ||
             (leaving == lp->rows &&
              lp->basis_column[d->i[i]] < entering)))) {
          mpq_set(*step, q);
          leaving = d->i[i];
	}
      } else {
        if (lp->lower.is_valid[var] == FALSE)
          continue;
        vector_get_element(&q, lp->xb, d->i[i]);
        mympq_sub(q, q, lp->lower.bound[var]);
        mympq_div(q, q, d->value[i]);
        if (leaving < 0 || mpq_cmp(q, *step) > 0 ||
	    (mpq_equal(q, *step) && (
             (leaving != lp->rows &&
              lp->basis_column[d->i[i]] < lp->basis_column[leaving]) ||
             (leaving == lp->rows &&
              lp->basis_column[d->i[i]] < entering)))) {
          mpq_set(*step, q);
          leaving = d->i[i];
	}
      }
    }
  }
  mpq_clear(q);
  return leaving;
}

void get_basis_inverse_column(LP* lp, EXLPvector* v, int col) {
  EXLPvector* e;
  int  i;
  mpq_t  q;

  if (lp->perturbation == 2) {
    vector_copy(v, lp->A->column[lp->basis_column[col]]);
    return;
  }

  e = new_vector(lp->rows);
  mpq_init(q);
  for (i = lp->rows-1; i >= col; i --) {
    //mpq_set_si(q, lp->rows-i, lp->rows);
    mpq_set_si(q, (lp->rows-i+col)*(lp->rows-i+col)+1, lp->rows*lp->rows+1);
    vector_set_element(e, q, i);
  }
  eta_file_ftran(lp->eta, e, v, NULL);
  vector_free(&e);
  mpq_clear(q);

  /*
  vector* e;

  e = new_vector(lp->rows);
  vector_set_element(e, mympq_one, col);
  eta_file_ftran(lp->eta, e, v, NULL);
  vector_free(&e);
  */
}

int sym_perturbation(LP* lp, EXLPvector* d, int ya_l_c, int* dgn_row, int dgns) {
  /* xb = xb - step d なのね */
  EXLPvector* v;
  int* dgn_row2;
  int  dgns2;
  int  i, j, sx, sd, sq, sb;
  mpq_t  q, best;

  dgn_row2 = my_malloc(lp->rows*sizeof(int));
  mpq_init(q);
  mpq_init(best);
  v = new_vector(lp->rows);

  //printf("{");
  for (i = 0; i < lp->rows; i ++) {
    //printf("%d ",dgns);
    get_basis_inverse_column(lp, v, i);
    vector_get_element(&best, v, dgn_row[0]);
    sb = sq = mpq_sgn(best);
    sd = mpq_sgn(*vector_get_element_ptr(d, dgn_row[0]));
    sx = (ya_l_c ? -1: 1);
    sx *= sd;
    // sx は x_b のdgn_row[0]成分の変化の符号
    mympq_div(best, best, *vector_get_element_ptr(d, dgn_row[0]));
    if ((sx < 0 && sq > 0) || (sx > 0 && sq < 0)) {
      mpq_abs(best, best);
    } else
    if ((sx > 0 && sq > 0) || (sx < 0 && sq < 0)) {
      mpq_abs(best, best);
      mpq_neg(best, best);
    }
    dgn_row2[0] = dgn_row[0];
    dgns2 = 1;

    for (j = 1; j < dgns; j ++) {
      vector_get_element(&q, v, dgn_row[j]);
      sq = mpq_sgn(q);
      sd = mpq_sgn(*vector_get_element_ptr(d, dgn_row[j]));
      sx = (ya_l_c ? -1: 1);
      sx *= sd;
      // sx は x_b のdgn_row[j]成分の変化の符号
      mympq_div(q, q, *vector_get_element_ptr(d, dgn_row[j]));
      if ((sx < 0 && sq > 0) || (sx > 0 && sq < 0)) {
	mpq_abs(q, q);
      } else
      if ((sx > 0 && sq > 0) || (sx < 0 && sq < 0)) {
        mpq_abs(q, q);
        mpq_neg(q, q);
      }
      if (mpq_cmp(q, best) < 0) {
	dgn_row2[0] = dgn_row[j];
	dgns2 = 1;
        mpq_set(best, q);
        sb = sq;
      } else if (mpq_equal(q, best)) {
	dgn_row2[dgns2++] = dgn_row[j];
      }
    }
    if (dgns2 == 1)
      break;
    for (j = 0; j < dgns2; j ++)
      dgn_row[j] = dgn_row2[j];
    dgns = dgns2;
  }
  //printf("}");

  free(dgn_row2);
  mpq_clear(q);
  mpq_clear(best);
  vector_free(&v);

  return dgn_row[0];
  /*
  vector* v;
  int* dgn_row2;
  int  dgns2;
  int  i, j;
  mpq_t  q, best;
  int  leaving = 0;

  dgn_row2 = my_malloc(lp->rows*sizeof(int));
  mpq_init(q);
  mpq_init(best);
  v = new_vector(lp->rows);
//putchar('{');
  for (i = 0; i < lp->rows; i ++) {
//printf("%d ",dgns);
    get_basis_inverse_column(lp, v, i);
    vector_get_element(&best, v, dgn_row[0]);
    mympq_div(best, best, *vector_get_element_ptr(d, dgn_row[0]));
    leaving = dgn_row[0];
    dgn_row2[0] = leaving;
    dgns2 = 1;
    for (j = 1; j < dgns; j ++) {
      vector_get_element(&q, v, dgn_row[j]);
      mympq_div(q, q, *vector_get_element_ptr(d, dgn_row[j]));
      if ((!ya_l_c && mpq_cmp(q, best) < 0) ||
	  ( ya_l_c && mpq_cmp(q, best) > 0)) {
      //if (( ya_l_c && mpq_cmp(q, best) < 0) ||
      //  (!ya_l_c && mpq_cmp(q, best) > 0)) {
	leaving = dgn_row[j];
	dgn_row2[0] = leaving;
	dgns2 = 1;
      } else if (mpq_equal(q, best)) {
	dgn_row2[dgns2++] = dgn_row[j];
      }
    }
    if (dgns2 == 1)
      break;
    for (j = 0; j < dgns2; j ++)
      dgn_row[j] = dgn_row2[j];
    dgns = dgns2;
  }
//putchar('}');fflush(stdout);
  free(dgn_row2);
  mpq_clear(q);
  mpq_clear(best);
  vector_free(&v);

  return leaving;
  */
}

void perturbate_aaa(LP* lp, int* dgn_row, int dgns) {
  int  i, var;
  mpq_t  q;

  mpq_init(q);

  for (i = 0; i < dgns; i ++) {
    //mpq_set_si(q, rand()%lp->rows, (rand()%lp->rows)*lp->rows+1);
    //mpq_set_si(q, rand()%lp->rows+1, 1); mpq_div_2exp(q, q, 100);
    mpq_set_si(q, rand()%lp->rows+1, 1); mpq_div_2exp(q, q, 10);
    var = lp->basis_column[i];
    //mpq_set_si(q, var, lp->vars*lp->vars*lp->vars*lp->vars);
    if (lp->upper.is_valid[var] &&
        mpq_equal(*vector_get_element_ptr(lp->xb, dgn_row[i]),
                  lp->upper.bound[var])) {
      vector_sub_element(lp->x, q, var);
      vector_sub_element(lp->xb, q, dgn_row[i]);
      if (lp->lower.is_valid[var] &&
	  mpq_cmp(lp->lower.bound[var],
		  *vector_get_element_ptr(lp->x, var)) > 0) {
	vector_set_element(lp->x, lp->lower.bound[var], var);
	vector_set_element(lp->xb, lp->lower.bound[var], dgn_row[i]);
      }
    } else if (lp->lower.is_valid[var] &&
	       mpq_equal(*vector_get_element_ptr(lp->xb, dgn_row[i]),
			 lp->lower.bound[var])) {
      vector_add_element(lp->x, q, var);
      vector_add_element(lp->xb, q, dgn_row[i]);
      if (lp->upper.is_valid[var] &&
	  mpq_cmp(lp->upper.bound[var],
		  *vector_get_element_ptr(lp->x, var)) > 0) {
	vector_set_element(lp->x, lp->upper.bound[var], var);
	vector_set_element(lp->xb, lp->upper.bound[var], dgn_row[i]);
      }
    }
  }

  mpq_clear(q);
}

/*  x[entering] <- x[entering] + step    */
/*  xb          <- xb          - step d  */
int select_leaving_row_harris(LP* lp, mpq_t* step, int ya_l_c,
                              EXLPvector* d, int entering) {
  /* *step には ya_l_c == TRUE なら正, ya_l_c == FALSE なら負の値を返す. */
  /* 戻り値は有界じゃないなら 負, 有界なら 基底の中で出る行.             */
  /* ただし lp->rows がかえったら出る変数無し.                           */
  /* 後, 記号摂動法導入                                                  */
  mpq_t  q, best;
  int  i, var, leaving;
  int* dgn_row;
  int  dgns;
  //static int  aaa=0;

  dgn_row = my_malloc(lp->rows*sizeof(int));
  dgns = 0;
  mpq_init(q);
  mpq_init(best);
  leaving = -1;

  if (ya_l_c) { /* step 正 */
    if (lp->upper.is_valid[entering]) {
      vector_get_element(step, lp->x, entering);
      mympq_sub(*step, lp->upper.bound[entering], *step);
      leaving = lp->rows;
    }
    for (i = 0; i < d->nonzeros; i ++) {
      var = lp->basis_column[d->i[i]];
      if (is_free_var(lp, var))
        continue;
      if (mpq_sgn(d->value[i]) > 0) {
        if (lp->lower.is_valid[var] == FALSE)
          continue;
        vector_get_element(&q, lp->xb, d->i[i]);
        mympq_sub(q, q, lp->lower.bound[var]);
        mympq_div(q, q, d->value[i]);
        if (leaving < 0 ||
            mpq_cmp(q, *step) < 0 ||
            (mpq_equal(q, *step) &&
             (leaving == lp->rows ||
	     !is_artificial_var(lp, lp->basis_column[leaving])) && (
	      is_artificial_var(lp, lp->basis_column[d->i[i]]) ||
              mpq_cmp(d->value[i], best) > 0 ||
              (mpq_equal(d->value[i], best) &&
                (leaving == lp->rows || var < lp->basis_column[leaving]))))) {
          mpq_set(*step, q);
          mpq_set(best, d->value[i]);
          leaving = d->i[i];
	}
        if (mpq_sgn(q) == 0) {
          dgn_row[dgns++] = d->i[i];
	}
      } else {
        if (lp->upper.is_valid[var] == FALSE)
          continue;
        vector_get_element(&q, lp->xb, d->i[i]);
        mympq_sub(q, q, lp->upper.bound[var]);
        mympq_div(q, q, d->value[i]);
        mpq_neg(d->value[i], d->value[i]);
        if (leaving < 0 ||
            mpq_cmp(q, *step) < 0 ||
            (mpq_equal(q, *step) &&
             (leaving == lp->rows ||
	     !is_artificial_var(lp, lp->basis_column[leaving])) && (
	      is_artificial_var(lp, lp->basis_column[d->i[i]]) ||
              mpq_cmp(d->value[i], best) > 0 ||
              (mpq_equal(d->value[i], best) &&
                (leaving == lp->rows || var < lp->basis_column[leaving]))))) {
          mpq_set(*step, q);
          mpq_set(best, d->value[i]);
          leaving = d->i[i];
	}
        mpq_neg(d->value[i], d->value[i]);
        if (mpq_sgn(q) == 0) {
          dgn_row[dgns++] = d->i[i];
	}
      }
    }
  } else { /* step 負 */
    if (lp->lower.is_valid[entering]) {
      vector_get_element(step, lp->x, entering);
      mympq_sub(*step, lp->lower.bound[entering], *step);
      leaving = lp->rows;
    }
    for (i = 0; i < d->nonzeros; i ++) {
      var = lp->basis_column[d->i[i]];
      if (is_free_var(lp, var))
        continue;
      if (mpq_sgn(d->value[i]) > 0) {
        if (lp->upper.is_valid[var] == FALSE)
          continue;
        vector_get_element(&q, lp->xb, d->i[i]);
        mympq_sub(q, q, lp->upper.bound[var]);
        mympq_div(q, q, d->value[i]);
        if (leaving < 0 ||
            mpq_cmp(q, *step) > 0 ||
            (mpq_equal(q, *step) &&
             (leaving == lp->rows ||
	     !is_artificial_var(lp, lp->basis_column[leaving])) && (
	      is_artificial_var(lp, lp->basis_column[d->i[i]]) ||
              mpq_cmp(d->value[i], best) > 0 ||
              (mpq_equal(d->value[i], best) &&
                (leaving == lp->rows || var < lp->basis_column[leaving]))))) {
          mpq_set(*step, q);
          mpq_set(best, d->value[i]);
          leaving = d->i[i];
	}
        if (mpq_sgn(q) == 0) {
          dgn_row[dgns++] = d->i[i];
	}
      } else {
        if (lp->lower.is_valid[var] == FALSE)
          continue;
        vector_get_element(&q, lp->xb, d->i[i]);
        mympq_sub(q, q, lp->lower.bound[var]);
        mympq_div(q, q, d->value[i]);
        mpq_neg(d->value[i], d->value[i]);
        if (leaving < 0 ||
            mpq_cmp(q, *step) > 0 ||
            (mpq_equal(q, *step) &&
	     (leaving == lp->rows ||
	     !is_artificial_var(lp, lp->basis_column[leaving])) && (
	      is_artificial_var(lp, lp->basis_column[d->i[i]]) ||
              mpq_cmp(d->value[i], best) > 0 ||
              (mpq_equal(d->value[i], best) &&
                (leaving == lp->rows || var < lp->basis_column[leaving]))))) {
          mpq_set(*step, q);
          mpq_set(best, d->value[i]);
          leaving = d->i[i];
	}
        mpq_neg(d->value[i], d->value[i]);
        if (mpq_sgn(q) == 0) {
          dgn_row[dgns++] = d->i[i];
	}
      }
    }
  }

  //if (dgns > 1 && lp->dual_simplex && lp->phase == 1) {
  /*
  if (dgns > 1 && lp->phase == 1) {
    if (aaa++ > 100) {
      perturbate_aaa(lp, dgn_row, dgns);
      leaving = -123;
      aaa = 0;
    }
  } else aaa = 0;
  */

  for (i = 0; i < dgns; i ++) {
    if (is_artificial_var(lp, lp->basis_column[dgn_row[i]])) {
      leaving = dgn_row[i];
      dgns = 1;
      break;
    }
  }
  if (lp->perturbation && dgns > 1) {
    //for (i = 0; i < d->nonzeros; i ++)
    //  vector_neg_element(d, d->i[i]);
    leaving = sym_perturbation(lp, d, ya_l_c, dgn_row, dgns);
    //for (i = 0; i < d->nonzeros; i ++)
    //  vector_neg_element(d, d->i[i]);
  }

  free(dgn_row);
  mpq_clear(q);
  mpq_clear(best);
  return leaving;
}

/*  x[entering] <- x[entering] + step    */
/*  xb          <- xb          - step d  */
void move_vertex(LP* lp, mpq_t step, EXLPvector* d, int e_column) {
  mpq_t  q;
  int  i;

  mpq_init(q);

  vector_add_element(lp->x, step, e_column);

  for (i = 0; i < d->nonzeros; i ++) {
    mympq_mul(q, step, d->value[i]);
    mympq_sub(q, q, *vector_get_element_ptr(lp->xb, d->i[i]));
    mpq_neg(q, q);
    vector_set_element(lp->xb, q, d->i[i]);
    vector_set_element(lp->x,  q, lp->basis_column[d->i[i]]);
  }

  mpq_clear(q);
}

void swap_basis(LP* lp, int e_column, int l_row) {

  vector_set_element(lp->xb, *vector_get_element_ptr(lp->x, e_column), l_row);
  vector_set_element(lp->cb, *vector_get_element_ptr(lp->c, e_column), l_row);

  lp->is_basis[lp->basis_column[l_row]] = FALSE;
  lp->is_basis[e_column] = TRUE;
  lp->basis_column[l_row] = e_column;
}

void update_object_function(LP* lp, EXLPvector* w, int e_column, int l_row) {
  mpq_t  q1, q2;
  int  i, j;

  mpq_init(q1);
  mpq_init(q2);

  vector_get_element(&q1, lp->c, e_column);
  mpq_neg(q1, q1);
  mympq_div(q1, q1, *vector_get_element_ptr(w, e_column));
  vector_set_element(lp->c, q1, lp->basis_column[l_row]);

  for (j = 0; j < w->nonzeros; j ++) {
    i = w->i[j];
    if (lp->is_basis[i] || is_const_var(lp, i))
      continue;
    mympq_mul(q2, w->value[j], q1);
    vector_add_element(lp->c, q2, i);
  }

  mpq_clear(q1);
  mpq_clear(q2);
}

int check_primal_feasibility(LP* lp) {
  int  i, j, k;
  mpq_t  q1, q2;

  mpq_init(q1);
  mpq_init(q2);

  for (i = 0; i < lp->rows; i ++) {
    mpq_set_si(q1, 0, 1);
    for (j = 0; j < lp->A->row[i]->nonzeros; j ++) {
      k = lp->A->row[i]->i[j];
      if (is_const_var(lp, k))
        continue;
      mpq_mul(q2, *vector_get_element_ptr(lp->x, k),
                  *matrix_get_element_ptr(lp->A, i, k));
      mympq_add(q1, q1, q2);
    }
    if (mpq_cmp(q1, *vector_get_element_ptr(lp->b, i)))
      return 0;
  }

  for (i = 0; i < lp->rows; i ++) {
    j = lp->basis_column[i];
    if (lp->upper.is_valid[j] &&
        mpq_cmp(lp->upper.bound[j], *vector_get_element_ptr(lp->xb, i)) < 0)
      return 0;
    if (lp->lower.is_valid[j] &&
        mpq_cmp(lp->lower.bound[j], *vector_get_element_ptr(lp->xb, i)) > 0)
      return 0;
  }

  return 1;
}

int check_dual_feasibility(LP* lp) {
  int  i;
  mpq_t* q;

  for (i = 0; i < lp->vars; i ++) {
    if (lp->is_basis[i] || is_const_var(lp, i))
      continue;
    q = vector_get_element_ptr(lp->c, i);
    if (mpq_sgn(*q) < 0 &&
        (!(lp->lower.is_valid[i]) || mpq_cmp(lp->lower.bound[i], *vector_get_element_ptr(lp->x, i))))
      return 0;
    if (mpq_sgn(*q) > 0 &&
        (!(lp->upper.is_valid[i]) || mpq_cmp(lp->upper.bound[i], *vector_get_element_ptr(lp->x, i))))
      return 0;
  }

  return 1;
}

int solve_lp_core(LP* lp) {
  EXLPvector* y;
  EXLPvector* d;
  EXLPvector* w;
  mpq_t  step;
  int  e_column, l_row;
  int  ya_l_c, result;
  int  i;
  int  co=0;
  int  co2=0;

  y = new_vector(lp->rows);
  d = new_vector(lp->rows);
  w = new_vector(lp->rows);
  mpq_init(step);

  reinversion(lp);

  lp->A_d = my_realloc(lp->A_d, lp->vars*sizeof(EXLPvector_d*));
  for (i = 0; i < lp->A->columns; i ++) {
    lp->A_d[i] = new_vector_d(lp->rows);
    vector_get_d(lp->A->column[i], lp->A_d[i]);
  }

  if (lp->pivot_rule & LP_PIVOT_DEVEX)
    devex_initialize(lp);
  else if (lp->pivot_rule & LP_PIVOT_STEEPEST_EDGE)
    steepest_edge_initialize(lp);
  else if (lp->pivot_rule & LP_PIVOT_PROJECTED_STEEPEST_EDGE)
    projected_initialize(lp);

  for (;;) {

    if (lp->c->nonzeros > 0 && is_artificial_var(lp, lp->c->i[0]) &&
        lp_artificials_are_zeros(lp)) {
      result = LP_RESULT_OPTIMAL;
      break;
    }

    eta_file_btran(lp->eta, lp->cb, y);

   select_entering:

    if (lp->verbose && co % 10 == 0) {
      lp_get_object_value(lp, &step);
      printf("iter: %6d\tobjective value: ", co);
      print_rational_as_float(step, 12);
      printf("\t");
      mpq_set_si(step, 0, 1);
      for (i = 0; i < lp->xb->nonzeros; i ++)
        if (is_artificial_var(lp, lp->basis_column[lp->xb->i[i]]))
          mympq_add(step, step, lp->xb->value[i]);
      print_rational_as_float(step, 12);
      /*
      printf("\t primal: ");
      vector_inner_product(&step, lp->cb, lp->xb);
      mpq_neg(step, step);
      print_rational_as_float(step, 12);
      printf(" dual: ");
      vector_inner_product(&step, y, lp->b);
      print_rational_as_float(step, 12);
      */
      putchar('\n');
    }

    e_column = select_entering_column(lp, y, &ya_l_c);
    if (e_column < 0) {
      result = LP_RESULT_OPTIMAL;
      break;
    }

    eta_file_ftran(lp->eta, lp->A->column[e_column], d, w);
    //if (d->nonzeros < 10) fprintf(stderr, "%d",d->nonzeros); else fprintf(stderr, "X");

    if ((lp->pivot_rule == LP_PIVOT_BLAND) ||
        (lp->pivot_rule == LP_PIVOT_BLAND2))
      l_row = select_leaving_row(lp, &step, ya_l_c, d, e_column);
    else
      l_row = select_leaving_row_harris(lp, &step, ya_l_c, d, e_column);

    if (l_row < 0) {
      result = LP_RESULT_UNBOUNDED;
      break;
    }

    move_vertex(lp, step, d, e_column);

    co++;
    if (mpq_sgn(step) == 0)
      co2 ++;

    if (l_row == lp->rows)
      goto select_entering;

    if (lp->pivot_rule & LP_PIVOT_DEVEX) {
      //if (devex_test(lp, d, e_column))
      if (co % 20 == 19) {
        swap_basis(lp, e_column, l_row);
        eta_file_update(lp, w, l_row);
        devex_rectify(lp);
      }
      else {
        devex_update_weights(lp, d, e_column, l_row);
        swap_basis(lp, e_column, l_row);
        eta_file_update(lp, w, l_row);
      }

    } else if (lp->pivot_rule & LP_PIVOT_STEEPEST_EDGE) {
      i = steepest_edge_update_table(lp, d, e_column, l_row);
      swap_basis(lp, e_column, l_row);
      eta_file_update(lp, w, l_row);
      if (i)
        putchar('!'),steepest_edge_initialize(lp);

    } else if (lp->pivot_rule & LP_PIVOT_PROJECTED_STEEPEST_EDGE) {
      if (co % 1000 == 999) {
        swap_basis(lp, e_column, l_row);
        eta_file_update(lp, w, l_row);
        projected_rectify(lp);
      } else {
        projected_update_table(lp, d, e_column, l_row);
        swap_basis(lp, e_column, l_row);
        eta_file_update(lp, w, l_row);
      }

    } else {
      swap_basis(lp, e_column, l_row);
      eta_file_update(lp, w, l_row);
    }
  }

  if (lp->verbose)
    printf("iteration count(degenerate): %d(%d)\n",co,co2);
  vector_free(&y);
  vector_free(&d);
  vector_free(&w);
  mpq_clear(step);
  for (i = 0; i < lp->A->columns; i ++)
    vector_d_free(&lp->A_d[i]);

  return result;
}

int solve_lp_primal(LP* lp) {

  vector_copy(lp->c_back, lp->c);

  if (lp->scaling)
    scaling(lp);
  if (set_basis_for_presolve(lp) == LP_RESULT_INFEASIBLE)
    return LP_RESULT_INFEASIBLE;
  lp->eta = new_eta_file(lp->rows);
  lp->phase = 1;
  lp->c_d = new_vector_d(lp->c->dimension);
  vector_get_d(lp->c, lp->c_d);
  if (solve_lp_core(lp) != LP_RESULT_OPTIMAL ||
      !lp_artificials_are_zeros(lp)) {
    return LP_RESULT_INFEASIBLE;
  }

  set_object_function(lp);
  //if (lp->scaling)
  //  scaling(lp);
  lp->phase = 2;
  vector_get_d(lp->c, lp->c_d);
  return solve_lp_core(lp);
}

void dual_init_set_fake_rhs(LP* lp) {
  mpq_t  q1, q2;
  int  i, j, k;

  lp->eta = new_eta_file(lp->rows);

  lp_sol_init_dual(lp);

  mpq_init(q1);
  mpq_init(q2);

  for (i = 0; i < lp->rows; i ++) {
    mpq_set_si(q1, 0, 1);
    for (j = 0; j < lp->A->row[i]->nonzeros; j ++) {
      k = lp->A->row[i]->i[j];
      if (is_const_var(lp, k))
        continue;
      vector_get_element(&q2, lp->A->column[k], i);
      mympq_mul(q2, q2, *vector_get_element_ptr(lp->x, k));
      mympq_add(q1, q1, q2);
    }
    vector_set_element(lp->b, q1, i);
  }

  mpq_clear(q1);
  mpq_clear(q2);
}

int set_basis_for_dual_simplex(LP* lp) {
  //int  i, j;
  //mpq_t q1, q2;
  int  r;

  if (lp->preprocess) {
    do {
      r = lp->rows;
      check_parallel_columns(lp);
      if (check_const_vars(lp) == LP_RESULT_INFEASIBLE)
        return LP_RESULT_INFEASIBLE;
      if (lp_check_redundancy(lp) == LP_RESULT_INFEASIBLE)
        return LP_RESULT_INFEASIBLE;
      lp_arrange_inequality2(lp);
      //check_trivial_constraints0(lp);
      if (check_trivial_constraints1(lp) == LP_RESULT_INFEASIBLE)
        return LP_RESULT_INFEASIBLE;
      if (check_trivial_constraints2(lp) == LP_RESULT_INFEASIBLE)
        return LP_RESULT_INFEASIBLE;
      check_trivial_constraints3(lp);
      check_tight_constraints1(lp);
      check_tight_constraints2(lp);
      if (check_parallel_constraints(lp) == LP_RESULT_INFEASIBLE)
        return LP_RESULT_INFEASIBLE;
    } while (r != lp->rows);
  }
  lp_arrange_inequality(lp);
  lp_add_slacks(lp);

  lp_select_dual_feasible_basis(lp);

  if (lp->verbose) {
    printf("variables(after preprocessing): %d + %d(slack)\n",
           lp_normal_vars(lp), lp_slack_vars(lp));
    printf("rows(after preprocessing): %d\n", lp->rows);
  }

  if (lp->scaling) {
    scaling(lp);
  }

  lp->b_back = new_vector(lp->rows);
  vector_copy(lp->b_back, lp->b);
  vector_copy(lp->c_back, lp->c);

  dual_init_set_fake_rhs(lp);

  return 0;
}

void prepare_for_dual_simplex_phase2(LP* lp) {
  EXLPvector* tmp;
  EXLPvector* b;
  int  i, j;
  mpq_t  q1, q2;

  tmp = new_vector(lp->rows);
  b = new_vector(lp->rows);

  vector_copy(b, lp->b_back);
  /*
  for (i = 0; i < lp->b->nonzeros; i ++)
    vector_sub_element(b, lp->b->value[i], lp->b->i[i]);
  */
  mpq_init(q1);
  mpq_init(q2);
  for (i = 0; i < lp->rows; i ++) {
    mpq_set_si(q1, 0, 1);
    for (j = 0; j < lp->vars; j ++) {
      if (is_const_var(lp, j))
        continue;
      mympq_mul(q2, *matrix_get_element_ptr(lp->A, i, j),
                    *vector_get_element_ptr(lp->x, j));
      mympq_add(q1, q1, q2);
    }
    vector_sub_element(b, q1, i);
  }
  mpq_clear(q1);
  mpq_clear(q2);

  eta_file_ftran(lp->eta, b, tmp, NULL);
  for (i = 0; i < tmp->nonzeros; i ++) {
    vector_add_element(lp->xb, tmp->value[i], tmp->i[i]);
    vector_add_element(lp->x, tmp->value[i], lp->basis_column[tmp->i[i]]);
  }

  vector_free(&tmp);
  vector_free(&b);
}

int select_leaving_row_dual_1(LP* lp, int* s) {
  int   i, j;

  for (i = 0; i < lp->rows; i ++) {
    j = lp->basis_column[i];
    if (lp->lower.is_valid[j]) {
      if (mpq_cmp(lp->lower.bound[j], *vector_get_element_ptr(lp->x, j)) > 0) {
        *s = -1;
        return i;
      }
    }
    if (lp->upper.is_valid[j]) {
      if (mpq_cmp(lp->upper.bound[j], *vector_get_element_ptr(lp->x, j)) < 0) {
        *s = 1;
        return i;
      }
    }
  }

  return -1;
}

int select_leaving_row_dual_2(LP* lp, int* s) {
  mpq_t  q, best;
  int  i, j, best_i;

  mpq_init(best);
  mpq_init(q);

  best_i = -1;
  for (i = 0; i < lp->rows; i ++) {
    j = lp->basis_column[i];
    if (lp->lower.is_valid[j]) {
      vector_get_element(&q, lp->x, j);
      mympq_sub(q, lp->lower.bound[j], q);
      if (mpq_sgn(q) > 0 &&
          (best_i < 0 || mpq_cmp(best, q) < 0)) {
        best_i = i;
        *s = -1;
        mpq_set(best, q);
      }
    }
    if (lp->upper.is_valid[j]) {
      vector_get_element(&q, lp->x, j);
      mympq_sub(q, q, lp->upper.bound[j]);
      if (mpq_sgn(q) > 0 &&
          (best_i < 0 || mpq_cmp(best, q) < 0)) {
        best_i = i;
        *s = 1;
        mpq_set(best, q);
      }
    }
  }

  mpq_clear(best);
  mpq_clear(q);
  return best_i;
}

int select_leaving_row_steepest_edge_dual(LP* lp, int* s) {
  mpq_t  q;
  int  i, j, best_i;
  double  d, best;

  mpq_init(q);

  best_i = -1;
  best = 0;

  for (i = 0; i < lp->rows; i ++) {
    j = lp->basis_column[i];
    if (lp->lower.is_valid[j]) {
      vector_get_element(&q, lp->x, j);
      mympq_sub(q, lp->lower.bound[j], q);
      d = mpq_get_d(q);
      d = d*d/lp->steepest_edge_table[i];
      if (mpq_sgn(q) > 0 &&
          (best_i < 0 || best < d)) {
        best_i = i;
        *s = -1;
        best = d;
      }
    }
    if (lp->upper.is_valid[j]) {
      vector_get_element(&q, lp->x, j);
      mympq_sub(q, q, lp->upper.bound[j]);
      d = mpq_get_d(q);
      d = d*d/lp->steepest_edge_table[i];
      if (mpq_sgn(q) > 0 &&
          (best_i < 0 || best < d)) {
        best_i = i;
        *s = 1;
        best = d;
      }
    }
  }

  mpq_clear(q);

  return best_i;
}

int select_leaving_row_dual(LP* lp, int* s) {

  if (lp->pivot_rule & LP_PIVOT_STEEPEST_EDGE)
    return select_leaving_row_steepest_edge_dual(lp, s);
  else
    return select_leaving_row_dual_2(lp, s);
}

int select_entering_column_dual(LP* lp, EXLPvector* w, int s) {
  mpq_t  best, b_w;
  mpq_t* q1;
  mpq_t  q2;
  int  var, best_var, c;

  mpq_init(best);
  mpq_init(b_w);
  mpq_init(q2);
  best_var = -1;

  if (s < 0) {
    for (var = 0; var < lp->vars; var ++) {
      if (lp->is_basis[var] || is_const_var(lp, var))
        continue;

      q1 = vector_get_element_ptr(w, var);
      if (mpq_sgn(*q1) < 0) {
        if (lp->upper.is_valid[var] &&
            mpq_cmp(*vector_get_element_ptr(lp->x, var),
                    lp->upper.bound[var]) >= 0)
          continue;
        vector_get_element(&q2, lp->c, var);
        mympq_div(q2, q2, *q1);
        mpq_abs(q2, q2);
        c = mpq_cmp(q2, best);
        if (best_var >= 0 && c > 0)
          continue;
        if (best_var >= 0 && c == 0) {
          mpq_neg(b_w, b_w);
          if (mpq_cmp(b_w, *q1) < 0) {
            mpq_neg(b_w, b_w);
            continue;
          }
          mpq_neg(b_w, b_w);
        }
        mpq_set(best, q2);
        mpq_neg(b_w, *q1);
        best_var = var;
      } else if (mpq_sgn(*q1) > 0) {
        if (lp->lower.is_valid[var] &&
            mpq_cmp(*vector_get_element_ptr(lp->x, var),
                    lp->lower.bound[var]) <= 0)
          continue;
        vector_get_element(&q2, lp->c, var);
        mympq_div(q2, q2, *q1);
        mpq_abs(q2, q2);
        c = mpq_cmp(q2, best);
        if (best_var >= 0 && c > 0)
          continue;
        if (best_var >= 0 && c == 0) {
          if (mpq_cmp(b_w, *q1) > 0)
            continue;
        }
        mpq_set(best, q2);
        mpq_set(b_w, *q1);
        best_var = var;
      }
    }
  } else {
    for (var = 0; var < lp->vars; var ++) {
      if (lp->is_basis[var] || is_const_var(lp, var))
        continue;

      q1 = vector_get_element_ptr(w, var);
      if (mpq_sgn(*q1) > 0) {
        if (lp->upper.is_valid[var] &&
            mpq_cmp(*vector_get_element_ptr(lp->x, var),
                    lp->upper.bound[var]) >= 0)
          continue;
        vector_get_element(&q2, lp->c, var);
        mympq_div(q2, q2, *q1);
        mpq_abs(q2, q2);
        c = mpq_cmp(q2, best);
        if (best_var >= 0 && c > 0)
          continue;
        if (best_var >= 0 && c == 0) {
          if (mpq_cmp(b_w, *q1) > 0)
            continue;
        }
        mpq_set(best, q2);
        mpq_set(b_w, *q1);
        best_var = var;
      } else if (mpq_sgn(*q1) < 0) {
        if (lp->lower.is_valid[var] &&
            mpq_cmp(*vector_get_element_ptr(lp->x, var),
                    lp->lower.bound[var]) <= 0)
          continue;
        vector_get_element(&q2, lp->c, var);
        mympq_div(q2, q2, *q1);
        mpq_abs(q2, q2);
        c = mpq_cmp(q2, best);
        if (best_var >= 0 && c > 0)
          continue;
        if (best_var >= 0 && c == 0) {
          mpq_neg(b_w, b_w);
          if (mpq_cmp(b_w, *q1) < 0) {
            mpq_neg(b_w, b_w);
            continue;
          }
          mpq_neg(b_w, b_w);
        }
        mpq_set(best, q2);
        mpq_neg(b_w, *q1);
        best_var = var;
      }
    }
  }

  mpq_clear(best);
  mpq_clear(b_w);
  mpq_clear(q2);

  return best_var;
}

int solve_lp_core_dual(LP* lp) {
  EXLPvector* e;
  EXLPvector* v;
  EXLPvector* w;
  EXLPvector* d;
  mpq_t  step;
  mpq_t  q;
  int  e_column, l_row;
  int  result;
  int  i, s;
  int  co = 0;
  int  co2 = 0;

  e = new_vector(lp->rows);
  v = new_vector(lp->rows);
  w = new_vector(lp->vars);
  d = new_vector(lp->rows);
  mpq_init(step);
  mpq_init(q);

  reinversion(lp);

  if (lp->pivot_rule & LP_PIVOT_STEEPEST_EDGE)
    steepest_edge_initialize_dual(lp);

  for ( ; ; ) {

    //lp_print_nonzero_vars(lp);
    //lp_print_basis(lp);
    if ((l_row = select_leaving_row_dual(lp, &s)) < 0) {
      result = LP_RESULT_OPTIMAL;
      break;
    }

    if (lp->verbose && co % 10 == 0) {
      lp_get_object_value(lp, &step);
      printf("iter: %6d\tobjective value: ", co);
      print_rational_as_float(step, 12);
      putchar('\n');
    }
    co ++;

    vector_zero_clear(e);
    vector_set_element(e, mympq_one, l_row);
    eta_file_btran(lp->eta, e, v);
    for (i = 0; i < lp->vars; i ++) {
      if (lp->is_basis[i] || is_const_var(lp, i))
        continue;
      vector_inner_product(&q, v, lp->A->column[i]);
      vector_set_element(w, q, i);
    }

    e_column = select_entering_column_dual(lp, w, s);
    if (e_column < 0) {
      result = LP_RESULT_INFEASIBLE;
      break;
    }

    eta_file_ftran(lp->eta, lp->A->column[e_column], d, v);

    vector_get_element(&step, lp->xb, l_row);
    if (s < 0) {
      mympq_sub(step, step, lp->lower.bound[lp->basis_column[l_row]]);
    } else {
      mympq_sub(step, step, lp->upper.bound[lp->basis_column[l_row]]);
    }
    mympq_div(step, step, *vector_get_element_ptr(w, e_column));
    move_vertex(lp, step, d, e_column);
    if (mpq_sgn(step) == 0)
      co2++;

    if (lp->pivot_rule & LP_PIVOT_STEEPEST_EDGE)
      i = steepest_edge_update_table_dual(lp, d, e_column, l_row);

    update_object_function(lp, w, e_column, l_row);
    swap_basis(lp, e_column, l_row);
    eta_file_update(lp, v, l_row);

    if ((lp->pivot_rule & LP_PIVOT_STEEPEST_EDGE) && i)
      putchar('!'),steepest_edge_initialize_dual(lp);
  }

  if (lp->verbose)
   printf("iteration count(degenerate): %d(%d)\n",co,co2);

  vector_free(&e);
  vector_free(&v);
  vector_free(&w);
  vector_free(&d);
  mpq_clear(step);
  mpq_clear(q);

  return result;
}

int solve_lp_dual(LP* lp) {
  EXLPvector* tmp;
  mpq_t  q;
  int  i;

  lp_arrange_inequality(lp);

  if (set_basis_for_dual_simplex(lp) == LP_RESULT_INFEASIBLE)
    return LP_RESULT_INFEASIBLE;

  lp->phase = 1;
  lp->c_d = new_vector_d(lp->c->dimension);
  vector_get_d(lp->c, lp->c_d);
  if (solve_lp_core(lp) != LP_RESULT_OPTIMAL)
    return LP_RESULT_DUAL_INFEASIBLE;

  tmp = new_vector(lp->rows);

  mpq_init(q);
  eta_file_btran(lp->eta, lp->cb, tmp);
  for (i = 0; i < lp->vars; i ++) {
    if (lp->is_basis[i] || is_const_var(lp, i))
      continue;
    vector_inner_product(&q, tmp, lp->A->column[i]);
    vector_sub_element(lp->c, q, i);
  }
  mpq_clear(q);

  prepare_for_dual_simplex_phase2(lp);

  //vector_copy(lp->b, lp->b_back);

  lp->phase = 2;
  vector_free(&tmp);
  return solve_lp_core_dual(lp);
}

int solve_lp(LP* lp) {
  if (lp->dual_simplex)
    return solve_lp_dual(lp);
  return solve_lp_primal(lp);
}
