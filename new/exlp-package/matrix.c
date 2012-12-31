#include "matrix.h"
#include "EXLPvector.h"
#include "mylib.h"
#include <string.h>

matrix* new_matrix(int rows, int columns) {
  matrix* mat;
  int  i;

  mat = my_malloc(sizeof(matrix));
  mat->rows    = rows;
  mat->columns = columns;

  mat->row =    my_malloc(sizeof(EXLPvector*)*mat->rows);
  mat->column = my_malloc(sizeof(EXLPvector*)*mat->columns);

  for (i = 0; i < mat->rows; i ++)
    mat->row[i] = new_vector_d(columns);
  for (i = 0; i < mat->columns; i ++)
    mat->column[i] = new_vector(rows);

  return mat;
}

void matrix_free(matrix* mat) {
  int  i;

  for (i = 0; i < mat->columns; i ++)
    vector_free(&mat->column[i]);

  for (i = 0; i < mat->rows; i ++)
    vector_d_free(&mat->row[i]);

  free(mat->column);
  free(mat->row);
  free(mat);
}

void matrix_resize(matrix* mat, int rows, int columns) {
  int  i, resize;

  for (i = rows; i < mat->rows; i ++)
    vector_d_free(&mat->row[i]);
  for (i = columns; i < mat->columns; i ++)
    vector_free(&mat->column[i]);

  mat->row    = my_realloc(mat->row, rows*sizeof(EXLPvector*));
  mat->column = my_realloc(mat->column, columns*sizeof(EXLPvector*));


//  for (i = 0; i < mat->rows; i ++)
  resize = MIN(rows, mat->rows);
  for (i = 0; i < resize; i ++)
    vector_d_resize(mat->row[i], columns);
  for ( ; i < rows; i ++)
    mat->row[i] = new_vector_d(columns);

//  for (i = 0; i < mat->columns; i ++)
  resize = MIN(columns, mat->columns);
  for (i = 0; i < resize; i ++)
    vector_resize(mat->column[i], rows);
  for ( ; i < columns; i ++)
    mat->column[i] = new_vector(rows);

  mat->rows    = rows;
  mat->columns = columns;
}

void matrix_zero_clear(matrix* mat) {
  int  i;

  for (i = 0; i < mat->rows; i ++)
    vector_d_zero_clear(mat->row[i]);

  for (i = 0; i < mat->columns; i ++)
    vector_zero_clear(mat->column[i]);
}

void matrix_copy(matrix* to, matrix* from) {
  int  i;

  //matrix_resize(to, from->rows, from->columns);
  if (from->rows != to->rows || from->columns != to->columns) {
    fprintf(stderr,"\n%d %d %d %d\n",from->rows,from->columns,to->rows,to->columns);
    exit(EXIT_FAILURE);
  }

  for (i = 0; i < to->rows; i ++) {
    vector_d_copy(to->row[i], from->row[i]);
  }
  for (i = 0; i < to->columns; i ++) {
    vector_copy(to->column[i], from->column[i]);
  }
}

int matrix_elements(matrix* mat) {
  int  i, n;

  n = 0;
  if (mat->rows < mat->columns) {
    for (i = 0; i < mat->rows; i ++)
      n += mat->row[i]->nonzeros;
  } else {
    for (i = 0; i < mat->columns; i ++)
      n += mat->column[i]->nonzeros;
  }
  return n;
}

void matrix_add_row(matrix* mat) {
  int  i;

  mat->rows ++;

  for (i = 0; i < mat->columns; i ++)
    vector_resize(mat->column[i], mat->rows);

  mat->row = my_realloc(mat->row, mat->rows*sizeof(EXLPvector*));
  mat->row[mat->rows-1] = new_vector_d(mat->columns);
}

void matrix_add_column(matrix* mat) {
  int  i;

  mat->columns ++;

  for (i = 0; i < mat->rows; i ++)
    vector_d_resize(mat->row[i], mat->columns);

  mat->column = my_realloc(mat->column, mat->columns*sizeof(EXLPvector*));
  mat->column[mat->columns-1] = new_vector(mat->rows);
}

void matrix_swap_rows(matrix* mat, int row1, int row2) {
  EXLPvector_d* v;
  int  i, col;
  int* t;

  if (row1 == row2)
    return;

  t = my_malloc(mat->columns*sizeof(int));
  for (i = 0; i < mat->row[row2]->nonzeros; i ++) {
    col = mat->row[row2]->i[i];
    t[col] = 0;
  }
  for (i = 0; i < mat->row[row1]->nonzeros; i ++) {
    col = mat->row[row1]->i[i];
    vector_swap_elements(mat->column[col], row1, row2);
    t[col] = 1;
  }
  for (i = 0; i < mat->row[row2]->nonzeros; i ++) {
    col = mat->row[row2]->i[i];
    if (!t[col])
      vector_swap_elements(mat->column[col], row1, row2);
  }
  free(t);

  v = mat->row[row1];
  mat->row[row1] = mat->row[row2];
  mat->row[row2] = v;
}

void matrix_swap_columns(matrix* mat, int col1, int col2) {
  EXLPvector* v;
  int  i, row;
  int* t;

  if (col1 == col2)
    return;

  t = my_malloc(mat->rows*sizeof(int));
  for (i = 0; i < mat->column[col2]->nonzeros; i ++) {
    row = mat->column[col2]->i[i];
    t[row] = 0;
  }
  for (i = 0; i < mat->column[col1]->nonzeros; i ++) {
    row = mat->column[col1]->i[i];
    vector_d_swap_elements(mat->row[row], col1, col2);
    t[row] = 1;
  }
  for (i = 0; i < mat->column[col2]->nonzeros; i ++) {
    row = mat->column[col2]->i[i];
    if (!t[row])
      vector_d_swap_elements(mat->row[row], col1, col2);
  }
  free(t);

  v = mat->column[col1];
  mat->column[col1] = mat->column[col2];
  mat->column[col2] = v;
}

void matrix_rev_row_sgn(matrix* mat, int row) {
  int  col;
  int  i;

  for (i = 0; i < mat->row[row]->nonzeros; i ++) {
    col = mat->row[row]->i[i];
    //mat//mpq_neg(mat->row[row]->value[i], mat->row[row]->value[i]);
    vector_neg_element(mat->column[col], row);
  }
}

void matrix_set_column(matrix* mat, EXLPvector* vec, int column) {
  int  i;

  for (i = 0; i < mat->column[column]->nonzeros; i ++)
    vector_d_delete_element(mat->row[mat->column[column]->i[i]], column);

  vector_copy(mat->column[column], vec);

  for (i = 0; i < vec->nonzeros; i ++)
    //mat//vector_set_element(mat->row[vec->i[i]], vec->value[i], column);
    vector_d_set_element(mat->row[vec->i[i]], 1, column);
}

void matrix_get_element(mpq_t* val, matrix* mat, int row, int column) {
  vector_get_element(val, mat->column[column], row);
}

mpq_t* matrix_get_element_ptr(matrix* mat, int row, int column) {
  return vector_get_element_ptr(mat->column[column], row);
}

int  matrix_element_is_zero(matrix* mat, int row, int column) {
  return vector_element_is_zero(mat->column[column], row);
}

void matrix_set_element(matrix* mat, mpq_t val, int row, int column) {
  vector_set_element(mat->column[column], val, row);
  if (mpq_sgn(val))
    vector_d_set_element(mat->row[row], /*val*/1, column);
  else
    vector_d_delete_element(mat->row[row], column);
}

void matrix_delete_element(matrix* mat, int row, int column) {
  vector_delete_element(mat->column[column], row);
  vector_d_delete_element(mat->row[row], column);
}

void matrix_add_element(matrix* mat, mpq_t q, int row, int col) {
  vector_add_element(mat->column[col], q, row);
  //mat//vector_add_element(mat->row[row], q, col);
  if (mpq_sgn(*vector_get_element_ptr(mat->column[col], row)))
    vector_d_set_element(mat->row[row], 1, col);
  else
    vector_d_delete_element(mat->row[row], col);
}

void matrix_sub_element(matrix* mat, mpq_t q, int row, int col) {
  vector_sub_element(mat->column[col], q, row);
  //mat//vector_sub_element(mat->row[row], q, col);
  if (mpq_sgn(*vector_get_element_ptr(mat->column[col], row)))
    vector_d_set_element(mat->row[row], 1, col);
  else
    vector_d_delete_element(mat->row[row], col);
}

void matrix_mul_element(matrix* mat, mpq_t q, int row, int col) {
  vector_mul_element(mat->column[col], q, row);
  //mat//vector_add_element(mat->row[row], q, col);
  if (mpq_sgn(*vector_get_element_ptr(mat->column[col], row)))
    vector_d_set_element(mat->row[row], 1, col);
  else
    vector_d_delete_element(mat->row[row], col);
}

void matrix_div_element(matrix* mat, mpq_t q, int row, int col) {
  vector_div_element(mat->column[col], q, row);
  //mat//vector_add_element(mat->row[row], q, col);
}

void matrix_dec_row_dimension(matrix* mat, int row) {
  int  i;

  for (i = 0; i < mat->columns; i ++)
    vector_dec_dimension(mat->column[i], row);

  mat->rows --;
  vector_d_free(&mat->row[row]);
  for (i = row; i < mat->rows; i ++)
    mat->row[i] = mat->row[i+1];
  mat->row = my_realloc(mat->row, mat->rows*sizeof(EXLPvector*));
}

void matrix_dec_column_dimension(matrix* mat, int column) {
  int  i;

  for (i = 0; i < mat->rows; i ++)
    vector_d_dec_dimension(mat->row[i], column);

  mat->columns --;
  vector_free(&mat->column[column]);
  for (i = column; i < mat->columns; i ++)
    mat->column[i] = mat->column[i+1];
  mat->column = my_realloc(mat->column, mat->columns*sizeof(EXLPvector*));
}

void matrix_row_scalar_product(matrix* mat, mpq_t s, int row) {
  int  i;

  if (mpq_equal(s, mympq_one))
    return;

  if (mpq_sgn(s) == 0) {
    for (i = 0; i < mat->row[row]->nonzeros; i ++)
      vector_delete_element(mat->column[mat->row[row]->i[i]], row);
    vector_d_zero_clear(mat->row[row]);
  } else {
    //mat//vector_scalar_product(mat->row[row], s);
    //mat//for (i = 0; i < mat->row[row]->nonzeros; i ++)
    //mat//  vector_set_element(mat->column[mat->row[row]->i[i]],
    //mat//                     mat->row[row]->value[i], row);
    for (i = 0; i < mat->row[row]->nonzeros; i ++)
      vector_mul_element(mat->column[mat->row[row]->i[i]], s, row);
  }
}

void matrix_column_scalar_product(matrix* mat, mpq_t s, int column) {
  int  i;

  if (mpq_equal(s, mympq_one))
    return;

  if (mpq_sgn(s) == 0) {
    for (i = 0; i < mat->column[column]->nonzeros; i ++)
      vector_d_delete_element(mat->row[mat->column[column]->i[i]], column);
    vector_zero_clear(mat->column[column]);
  } else {
    vector_scalar_product(mat->column[column], s);
    //mat//for (i = 0; i < mat->column[column]->nonzeros; i ++)
    //mat//  vector_set_element(mat->row[mat->column[column]->i[i]],
    //mat//                     mat->column[column]->value[i], column);
  }
}

void matrix_print(matrix* mat) {
  int    i, j;
  mpq_t* val;

  for (i = 0; i < mat->rows; i ++) {
    for (j = 0; j < mat->columns; j ++) {
      val = matrix_get_element_ptr(mat, i, j);
      //print_rational_as_float(*val, 3);
      //printf("\t");
      print_rational_as_float(*val, 1);
      //printf(" ");
      printf("\t");
    }
    printf("\n");
  }
}

void matrix_print2(matrix* mat) {
  int    i, j;
  mpq_t* val;

  for (i = 0; i < mat->rows; i ++) {
    for (j = 0; j < mat->columns; j ++) {
      val = matrix_get_element_ptr(mat, i, j);
      if (mpq_sgn(*val))
        putchar('*');
      else
        putchar('.');
    }
    putchar('\n');
  }
}


permutation_matrix* new_permutation_matrix(int dimension) {
  permutation_matrix* p;
  int  i;

  p = my_malloc(sizeof(permutation_matrix));
  p->dimension = dimension;
  p->row    = my_malloc(dimension*sizeof(int));
  p->column = my_malloc(dimension*sizeof(int));

  for (i = 0; i < dimension; i ++) {
    p->row[i]    = i;
    p->column[i] = i;
  }

  p->s = p->dimension;
  p->t = -1;

  return p;
}

void permutation_matrix_free(permutation_matrix* p) {
  free(p->row);
  free(p->column);
  free(p);
}

void permutation_matrix_inv(permutation_matrix* p) {
  int* tmp;

  tmp = p->row;
  p->row = p->column;
  p->column = tmp;
}

void permutation_matrix_copy(permutation_matrix* to,
                             permutation_matrix* from) {
  memmove(to->row,    from->row,    sizeof(int)*to->dimension);
  memmove(to->column, from->column, sizeof(int)*to->dimension);
  to->s = from->s;
  to->t = from->t;
}

void permutation_matrix_mul_core(permutation_matrix* p,
                                 permutation_matrix* q,
                                 int s, int t) {
  /* p = p*q */
  int  i;
  int *tmp;

  tmp = my_malloc(p->dimension*sizeof(int));

  for (i = s; i <= t; i ++)
    tmp[i] = p->column[q->column[i]];

  for (i = s; i <= t; i ++) {
    p->column[i] = tmp[i];
    p->row[tmp[i]] = i;
  }

  free(tmp);
}

void permutation_matrix_mul(permutation_matrix* p, permutation_matrix* q) {
  /* p = p*q */
  p->s = MIN(p->s, q->s);
  p->t = MAX(p->t, q->t);

  permutation_matrix_mul_core(p, q, p->s, p->t);
}

void permutation_matrix_mul2_core(permutation_matrix* p,
                                  permutation_matrix* q,
                                  int s, int t) {
  /* q = p*q */
  int  i;
  int *tmp;

  tmp = my_malloc(p->dimension*sizeof(int));

  for (i = s; i <= t; i ++)
    tmp[i] = p->column[q->column[i]];

  for (i = s; i <= t; i ++) {
    q->column[i] = tmp[i];
    q->row[tmp[i]] = i;
  }

  free(tmp);
}

void permutation_matrix_mul2(permutation_matrix* p, permutation_matrix* q) {
  /* q = p*q */
  q->s = MIN(p->s, q->s);
  q->t = MAX(p->t, q->t);

  permutation_matrix_mul2_core(p, q, q->s, q->t);
}

int cmp(const void* a, const void* b) {
  if (*((int*)a) < *((int*)b))
    return -1;
  return 1;
}

void vector_permutate_core(permutation_matrix* p, EXLPvector* v,
                           int s, int t) {
  int *i;
  int  j;

  if (t <= s)
    return;

  i = my_malloc(v->blocks*VECTOR_BLOCK_SIZE*sizeof(int));

  for (j = 0; j < v->nonzeros; j ++) {
    v->ptr[v->i[j]] = -1;
    i[j] = p->column[v->i[j]];
  }

  free(v->i);
  v->i = i;

  for (j = 0; j < v->nonzeros; j ++)
    v->ptr[v->i[j]] = j;
}

void vector_permutate(permutation_matrix* p, EXLPvector* v) {
  vector_permutate_core(p, v, p->s, p->t);
}

void vector_permutate_inv_core(permutation_matrix* p, EXLPvector* v,
                               int s, int t) {
  // v := p v
  int *i;
  int  j;

  if (t <= s)
    return;

  i = my_malloc(v->blocks*VECTOR_BLOCK_SIZE*sizeof(int));

  for (j = 0; j < v->nonzeros; j ++) {
    v->ptr[v->i[j]] = -1;
    i[j] = p->row[v->i[j]];
  }

  free(v->i);
  v->i = i;

  for (j = 0; j < v->nonzeros; j ++)
    v->ptr[v->i[j]] = j;
}

void vector_permutate_inv(permutation_matrix* p, EXLPvector* v) {
  vector_permutate_inv_core(p, v, 0, p->dimension-1);
}

void vector_d_permutate_core(permutation_matrix* p, EXLPvector_d* v,
                             int s, int t) {
  int *i;
  int  j;

  if (t <= s)
    return;

  i = my_malloc(v->blocks*VECTOR_BLOCK_SIZE*sizeof(int));

  for (j = 0; j < v->nonzeros; j ++) {
    v->ptr[v->i[j]] = -1;
    i[j] = p->column[v->i[j]];
  }

  free(v->i);
  v->i = i;

  for (j = 0; j < v->nonzeros; j ++)
    v->ptr[v->i[j]] = j;
}

void vector_d_permutate(permutation_matrix* p, EXLPvector_d* v) {
  vector_d_permutate_core(p, v, p->s, p->t);
}

void vector_d_permutate_inv_core(permutation_matrix* p, EXLPvector_d* v,
                                 int s, int t) {
  // v := p v
  int *i;
  int  j;

  if (t <= s)
    return;

  i = my_malloc(v->blocks*VECTOR_BLOCK_SIZE*sizeof(int));

  for (j = 0; j < v->nonzeros; j ++) {
    v->ptr[v->i[j]] = -1;
    i[j] = p->row[v->i[j]];
  }

  free(v->i);
  v->i = i;

  for (j = 0; j < v->nonzeros; j ++)
    v->ptr[v->i[j]] = j;
}

void vector_d_permutate_inv(permutation_matrix* p, EXLPvector_d* v) {
  vector_d_permutate_inv_core(p, v, 0, p->dimension-1);
}

void matrix_permutate_core(permutation_matrix* p, matrix* m,
                           permutation_matrix* q,
                           int s, int t) {
  EXLPvector** tmp;
  int  i;

  if (t <= s)
    return;

  for (i = 0; i < m->rows; i ++) {
    vector_permutate_core(p, m->column[i], s, t);
    vector_d_permutate_inv_core(q, m->row[i], s, t);
  }

  tmp = my_malloc(m->rows*sizeof(EXLPvector*));

  for (i = s; i <= t; i ++)
    tmp[i] = (void*)m->row[i];
  for (i = s; i <= t; i ++)
    m->row[i] = (void*)tmp[p->row[i]];

  for (i = s; i <= t; i ++)
    tmp[i] = m->column[i];
  for (i = s; i <= t; i ++)
    m->column[i] = tmp[q->column[i]];

  free(tmp);
}

void matrix_permutate(permutation_matrix* p, matrix* m,
                      permutation_matrix* q) {
  int  s, t;

  s = MIN(p->s, q->s);
  t = MAX(p->t, q->t);
  matrix_permutate_core(p, m, q, s, t);
}

eta_matrix* new_eta_matrix(int dimension) {
  eta_matrix* mat;

  mat = my_malloc(sizeof(eta_matrix));
  mat->dimension = dimension;
  mat->eta_vector = new_vector(dimension);

  return mat;
}

void eta_matrix_free(eta_matrix* E) {
  vector_free(&E->eta_vector);
  free(E);
}

void eta_matrix_mul_vector(eta_matrix* E, EXLPvector* v) {
  /* v := Ev */

  mpq_t  q1, q2;
  int  i;

  if (vector_element_is_zero(v, E->eta_column))
    return;

  mpq_init(q1);
  mpq_init(q2);
  vector_get_element(&q2, v, E->eta_column);

  for (i = 0; E->eta_vector->i[i] != E->eta_column; i ++) {
    mympq_mul(q1, q2, E->eta_vector->value[i]);
    vector_add_element(v, q1, E->eta_vector->i[i]);
  }
  vector_mul_element(v, E->eta_vector->value[i], E->eta_column);
  for (i ++; i < E->eta_vector->nonzeros; i ++) {
    mympq_mul(q1, q2, E->eta_vector->value[i]);
    vector_add_element(v, q1, E->eta_vector->i[i]);
  }

  mpq_clear(q1);
  mpq_clear(q2);
}

void eta_matrix_solve_yE_is_equal_to_v(eta_matrix* E, EXLPvector* v) {
  /* yE = v を解いて解(y)を v に返す. */
  mpq_t  q1, q2;
  mpq_t* q3;

  if (E->eta_vector->nonzeros == 1) {
    vector_div_element(v, E->eta_vector->value[0], E->eta_column);
    return;
  }

  mpq_init(q1);
  mpq_init(q2);

  if (!vector_element_is_zero(v, E->eta_column)) {
    vector_get_element(&q2, v, E->eta_column);
    vector_delete_element(v, E->eta_column);
    vector_inner_product(&q1, v, E->eta_vector);
    mympq_sub(q1, q1, q2);
  } else {
    vector_inner_product(&q1, v, E->eta_vector);
  }

  if (mpq_sgn(q1)) {
    mpq_neg(q1, q1);
    q3 = vector_get_element_ptr(E->eta_vector, E->eta_column);
    mympq_div(q1, q1, *q3);
    vector_set_element(v, q1, E->eta_column);
  }

  mpq_clear(q1);
  mpq_clear(q2);
}

void eta_matrix_solve_Ex_is_equal_to_v(eta_matrix* E, EXLPvector* v) {
  /* Ex = v を解いて解(x)を v に返す. */
  mpq_t q1, q2;
  int  i;

  if (vector_element_is_zero(v, E->eta_column))
    return;

  mpq_init(q1);
  mpq_init(q2);
  vector_get_element(&q1, v, E->eta_column);

  mympq_div(q1, q1, *vector_get_element_ptr(E->eta_vector, E->eta_column));
  vector_set_element(v, q1, E->eta_column);

  for (i = 0; E->eta_vector->i[i] != E->eta_column; i ++) {
    mympq_mul(q2, q1, E->eta_vector->value[i]);
    vector_sub_element(v, q2, E->eta_vector->i[i]);
  }
  for (i ++; i < E->eta_vector->nonzeros; i ++) {
    mympq_mul(q2, q1, E->eta_vector->value[i]);
    vector_sub_element(v, q2, E->eta_vector->i[i]);
  }

  mpq_clear(q1);
  mpq_clear(q2);
}

eta_matrix_d* new_eta_matrix_d(int dimension) {
  eta_matrix_d* mat;

  mat = my_malloc(sizeof(eta_matrix_d));
  mat->eta_vector = new_vector_d(dimension);

  return mat;
}

void eta_matrix_d_free(eta_matrix_d* E) {
  vector_d_free(&E->eta_vector);
  free(E);
}


eta_singleton* new_eta_singleton(void) {
  eta_singleton* es;

  es = my_malloc(sizeof(eta_singleton));
  es->row = es->column = 0;
  mpq_init(es->value);

  return es;
}

void eta_singleton_free(eta_singleton* es) {
  mpq_clear(es->value);
  free(es);
}

void eta_singleton_set(eta_singleton* es, mpq_t q, int row, int column) {
  mpq_set(es->value, q);
  es->row = row;
  es->column = column;
}

eta_singleton_d* new_eta_singleton_d(void) {
  eta_singleton_d* es;

  es = my_malloc(sizeof(eta_singleton_d));
  es->row = es->column = 0;
  es->value = 0;

  return es;
}

void eta_singleton_d_free(eta_singleton_d* es) {
  free(es);
}

void eta_singleton_d_set(eta_singleton_d* es, double d, int row, int column) {
  es->value = d;
  es->row = row;
  es->column = column;
}

void eta_singleton_get_d(eta_singleton* es, eta_singleton_d* es_d) {
  es_d->value = mpq_get_d(es->value);
  es_d->row = es->row;
  es_d->column = es->column;
}


void eta_matrix_d_mul_vector(eta_matrix_d* E, EXLPvector_d* v) {
  /* v := Ev */
  double  q1, q2;
  int  i;

  if (vector_d_element_is_zero(v, E->eta_column))
    return;

  q2 = vector_d_get_element(v, E->eta_column);

  for (i = 0; i < E->eta_vector->nonzeros; i ++) {
    q1 = q2 * E->eta_vector->value[i];
    if (E->eta_vector->i[i] != E->eta_column)
      vector_d_add_element(v, q1, E->eta_vector->i[i]);
    else
      vector_d_set_element(v, q1, E->eta_vector->i[i]);
  }
}

void eta_matrix_d_solve_yE_is_equal_to_v(eta_matrix_d* E, EXLPvector_d* v) {
  /* yE = v を解いて解(y)を v に返す. */
  double  q1, q2;

  if (E->eta_vector->nonzeros == 1) {
    vector_d_div_element(v, E->eta_vector->value[0], E->eta_column);
    return;
  }

  if (!vector_d_element_is_zero(v, E->eta_column)) {
    q2 = vector_d_get_element(v, E->eta_column);
    vector_d_delete_element(v, E->eta_column);
    q1 = vector_d_inner_product(v, E->eta_vector);
    q1 -= q2;
  } else {
    q1 = vector_d_inner_product(v, E->eta_vector);
  }

  if (q1 != 0) {
    q2 = -vector_d_get_element(E->eta_vector, E->eta_column);
    q1 /= q2;
    vector_d_set_element(v, q1, E->eta_column);
  }
}

void eta_matrix_d_solve_Ex_is_equal_to_v(eta_matrix_d* E, EXLPvector_d* v) {
  /* Ex = v を解いて解(x)を v に返す. */
  double  q1, q2;
  int  i;

  if (vector_d_element_is_zero(v, E->eta_column))
    return;

  q1 = vector_d_get_element(v, E->eta_column);

  q2 = vector_d_get_element(E->eta_vector, E->eta_column);
  q1 /= q2;
  vector_d_set_element(v, q1, E->eta_column);

  for (i = 0; i < E->eta_vector->nonzeros; i ++) {
    if (E->eta_vector->i[i] == E->eta_column)
      continue;
    q2 = q1 * E->eta_vector->value[i];
    vector_d_sub_element(v, q2, E->eta_vector->i[i]);
  }
}
