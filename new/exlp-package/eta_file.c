#include "eta_file.h"
#include "mylib.h"

eta_file* new_eta_file(int lp_rows) {
  eta_file* ef;
  int  i;

  ef = my_malloc(sizeof(eta_file));

  ef->lp_rows = lp_rows;

  ef->k = 0;

  ef->P = my_malloc(lp_rows*sizeof(int));

  ef->L = my_malloc(lp_rows*sizeof(eta_matrix*));
  for (i = 0; i < lp_rows; i ++) {
    ef->L[i] = new_eta_matrix(lp_rows);
  }
  ef->Ls = my_malloc(lp_rows*ETA_MAX*sizeof(eta_singleton*));
  for (i = 0; i < lp_rows*ETA_MAX; i ++) {
    ef->Ls[i] = new_eta_singleton();
  }

  ef->L_d = my_malloc(lp_rows*sizeof(eta_matrix_d*));
  for (i = 0; i < lp_rows; i ++) {
    ef->L_d[i] = new_eta_matrix_d(lp_rows);
  }
  ef->Ls_d = my_malloc(lp_rows*ETA_MAX*sizeof(eta_singleton_d*));
  for (i = 0; i < lp_rows*ETA_MAX; i ++) {
    ef->Ls_d[i] = new_eta_singleton_d();
  }

  ef->U = new_matrix(lp_rows, lp_rows);
  /* とりあえずこんだけ */

  ef->U_d = (EXLPvector_d**)my_malloc(lp_rows*sizeof(EXLPvector_d*));
  for (i = 0; i < lp_rows; i ++)
    ef->U_d[i] = new_vector_d(lp_rows);

  ef->Q = new_permutation_matrix(lp_rows);
  ef->R = new_permutation_matrix(lp_rows);

  return ef;
}

void eta_file_free(eta_file* ef) {
  int  i;

  free(ef->P);

  for (i = 0; i < ef->lp_rows; i ++) {
    eta_matrix_free(ef->L[i]);
  }
  free(ef->L);
  for (i = 0; i < ef->lp_rows*ETA_MAX; i ++) {
    eta_singleton_free(ef->Ls[i]);
  }
  free(ef->Ls);

  for (i = 0; i < ef->lp_rows; i ++) {
    eta_matrix_d_free(ef->L_d[i]);
  }
  free(ef->L_d);
  for (i = 0; i < ef->lp_rows*ETA_MAX; i ++) {
    eta_singleton_d_free(ef->Ls_d[i]);
  }
  free(ef->Ls_d);

  matrix_free(ef->U);

  for (i = 0; i < ef->lp_rows; i ++)
    vector_d_free(&ef->U_d[i]);
  free(ef->U_d);

  permutation_matrix_free(ef->Q);
  permutation_matrix_free(ef->R);

  free(ef);
}

void eta_file_btran(eta_file* ef, EXLPvector* c, EXLPvector* y) {
  /* yB = c  を解く. */
  eta_matrix* U;
  mpq_t  q;
  int  i;

  U = my_malloc(sizeof(eta_matrix));
  vector_copy(y, c);
  mpq_init(q);

  vector_permutate(ef->R, y);

  for (i = 0; i < ef->U->columns; i ++) {
    U->eta_column = i;
    U->eta_vector = ef->U->column[i];
    eta_matrix_solve_yE_is_equal_to_v(U, y);
  }

  vector_permutate(ef->Q, y);

  for (i = ef->s-1; i >= 0; i --) {
    if (ef->Ls[i]->row == ef->Ls[i]->column) {
      vector_mul_element(y, ef->Ls[i]->value, ef->Ls[i]->row);
    } else {
      vector_get_element(&q, y, ef->Ls[i]->row);
      mympq_mul(q, q, ef->Ls[i]->value);
      vector_add_element(y, q, ef->Ls[i]->column);
    }
  }
  for (i = ef->U->columns-1; i >= 0; i --) {
    if (ef->L[i]->eta_vector->nonzeros == 1) {
      vector_mul_element(y, ef->L[i]->eta_vector->value[0], i);
    } else {
      vector_inner_product(&q, y, ef->L[i]->eta_vector);
      vector_set_element(y, q, i);
    }
    vector_swap_elements(y, i, ef->P[i]);
  }

  mpq_clear(q);
  free(U);
}

void eta_file_btran_d(eta_file* ef, EXLPvector_d* c, EXLPvector_d* y) {
  /* yB = c  を解く. */
  eta_matrix_d* U;
  double  d;
  int  i;

  U = my_malloc(sizeof(eta_matrix_d));
  vector_d_copy(y, c);

  vector_d_permutate(ef->R, y);

  for (i = 0; i < ef->U->columns; i ++) {
    U->eta_column = i;
    U->eta_vector = ef->U_d[i];
    eta_matrix_d_solve_yE_is_equal_to_v(U, y);
  }

  vector_d_permutate(ef->Q, y);

  for (i = ef->s-1; i >= 0; i --) {
    if (ef->Ls_d[i]->row == ef->Ls_d[i]->column) {
      vector_d_mul_element(y, ef->Ls_d[i]->value, ef->Ls_d[i]->row);
    } else {
      d = vector_d_get_element(y, ef->Ls_d[i]->row);
      d *= ef->Ls_d[i]->value;
      vector_d_add_element(y, d, ef->Ls_d[i]->column);
    }
  }

  for (i = ef->U->columns-1; i >= 0; i --) {
    if (ef->L_d[i]->eta_vector->nonzeros == 1) {
      vector_d_mul_element(y, ef->L_d[i]->eta_vector->value[0], i);
    } else {
      d = vector_d_inner_product(y, ef->L_d[i]->eta_vector);
      vector_d_set_element(y, d, i);
    }
    vector_d_swap_elements(y, i, ef->P[i]);
  }

  free(U);
}

void eta_file_ftran(eta_file* ef, EXLPvector* b, EXLPvector* x, EXLPvector* w) {
  /* Bx = b  を解く. */
  eta_matrix* U;
  int  i;
  mpq_t  q;

  mpq_init(q);

  U = my_malloc(sizeof(eta_matrix));
  vector_copy(x, b);

  for (i = 0; i < ef->U->columns; i ++) {
    vector_swap_elements(x, i, ef->P[i]);
    eta_matrix_mul_vector(ef->L[i], x);
  }
  for (i = 0; i < ef->s; i ++) {
    if (ef->Ls[i]->row == ef->Ls[i]->column) {
      vector_mul_element(x, ef->Ls[i]->value, ef->Ls[i]->row);
    } else {
      vector_get_element(&q, x, ef->Ls[i]->column);
      if (mpq_sgn(q)==0)continue;
      mympq_mul(q, q, ef->Ls[i]->value);
      vector_add_element(x, q, ef->Ls[i]->row);
    }
  }

  vector_permutate_inv(ef->Q, x);
  if (w != NULL)
    vector_copy(w, x);
  //if (w!=NULL)printf("(%d ", w->nonzeros),fflush(stdout);

  for (i = ef->U->columns-1; i >= 0; i --) {
    U->eta_column = i;
    U->eta_vector = ef->U->column[i];
    eta_matrix_solve_Ex_is_equal_to_v(U, x);
  }

//fprintf(stderr, "%d ",x->i[0]);
  vector_permutate_inv(ef->R, x);
  //printf(" %d)", x->nonzeros);fflush(stdout);

  mpq_clear(q);
  free(U);
}

void eta_file_ftran_d(eta_file* ef, EXLPvector_d* b, EXLPvector_d* x) {
  /* Bx = b  を解く. */
  eta_matrix_d* U;
  double  d;
  int  i;

  U = my_malloc(sizeof(eta_matrix_d));
  vector_d_copy(x, b);

  for (i = 0; i < ef->U->columns; i ++) {
    vector_d_swap_elements(x, i, ef->P[i]);
    eta_matrix_d_mul_vector(ef->L_d[i], x);
  }
  for (i = 0; i < ef->s; i ++) {
    if (ef->Ls_d[i]->row == ef->Ls_d[i]->column) {
      vector_d_mul_element(x, ef->Ls_d[i]->value, ef->Ls_d[i]->row);
    } else {
      d = vector_d_get_element(x, ef->Ls_d[i]->column);
      d *= ef->Ls_d[i]->value;
      vector_d_add_element(x, d, ef->Ls_d[i]->row);
    }
  }

  vector_d_permutate_inv(ef->Q, x);

  for (i = ef->U->columns-1; i >= 0; i --) {
    U->eta_column = i;
    U->eta_vector = ef->U_d[i];
    eta_matrix_d_solve_Ex_is_equal_to_v(U, x);
  }

  vector_d_permutate_inv(ef->R, x);
  //printf(" %d)", x->nonzeros);fflush(stdout);

  free(U);
}

