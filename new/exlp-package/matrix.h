#ifndef MATRIX_H
#define MATRIX_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include "EXLPvector.h"

typedef struct {
  int  rows;
  int  columns;
  EXLPvector** column;
  EXLPvector_d** row;
} matrix;

typedef struct {
  int  dimension;
  int  s, t;    // 恒等変換でない範囲
  int* row;     // row[i]:    i行目の 1 の位置
  int* column;  // column[i]: i列目の 1 の位置
} permutation_matrix;

typedef struct {
  int  dimension;
  int  eta_column;
  EXLPvector* eta_vector;
} eta_matrix;

typedef struct {
  int  row;
  int  column;
  mpq_t  value;
} eta_singleton;

typedef struct {
  int  eta_column;
  EXLPvector_d* eta_vector;
} eta_matrix_d;

typedef struct {
  int  row;
  int  column;
  double  value;
} eta_singleton_d;

matrix* new_matrix(int rows, int columns);
void matrix_free(matrix* mat);
void matrix_resize(matrix* mat, int rows, int columns);
void matrix_zero_clear(matrix* mat);
void matrix_copy(matrix* to, matrix* from);
int  matrix_elements(matrix* mat);
void matrix_add_row(matrix* mat);
void matrix_add_column(matrix* mat);
void matrix_swap_rows(matrix* mat, int row1, int row2);
void matrix_swap_columns(matrix* mat, int col1, int col2);
void matrix_rev_row_sgn(matrix* mat, int row);
void matrix_set_column(matrix* mat, EXLPvector* vec, int column);
void matrix_get_element(mpq_t* val, matrix* mat, int row, int column);
mpq_t* matrix_get_element_ptr(matrix* mat, int row, int column);
int  matrix_element_is_zero(matrix* mat, int row, int column);
void matrix_set_element(matrix* mat, mpq_t val, int row, int column);
void matrix_delete_element(matrix* mat, int row, int column);
void matrix_add_element(matrix* mat, mpq_t q, int row, int col);
void matrix_sub_element(matrix* mat, mpq_t q, int row, int col);
void matrix_mul_element(matrix* mat, mpq_t q, int row, int col);
void matrix_div_element(matrix* mat, mpq_t q, int row, int col);
void matrix_dec_row_dimension(matrix* mat, int row);
void matrix_dec_column_dimension(matrix* mat, int column);
void matrix_row_scalar_product(matrix* mat, mpq_t s, int row);
void matrix_column_scalar_product(matrix* mat, mpq_t s, int column);
void matrix_print(matrix* mat);
void matrix_print2(matrix* mat);

permutation_matrix* new_permutation_matrix(int dimension);
void permutation_matrix_free(permutation_matrix* p);
void permutation_matrix_mul(permutation_matrix* p, permutation_matrix* q);
void permutation_matrix_mul2(permutation_matrix* p, permutation_matrix* q);
void permutation_matrix_inv(permutation_matrix* p);
void permutation_matrix_copy(permutation_matrix* to, permutation_matrix* from);
void vector_permutate(permutation_matrix* p, EXLPvector* v);
void vector_permutate_inv_core(permutation_matrix* p, EXLPvector* v, int s, int t);
void vector_permutate_inv(permutation_matrix* p, EXLPvector* v);
void vector_d_permutate(permutation_matrix* p, EXLPvector_d* v);
void vector_d_permutate_inv(permutation_matrix* p, EXLPvector_d* v);
void matrix_permutate(permutation_matrix* p, matrix* m, permutation_matrix* q);

eta_matrix* new_eta_matrix(int dimension);
void eta_matrix_free(eta_matrix* E);
void eta_matrix_mul_vector(eta_matrix* E, EXLPvector* v);
//void matrix_triangular_factorization(matrix* A,
//                                     eta_matrix** L, int* P);
void eta_matrix_solve_yE_is_equal_to_v(eta_matrix* E, EXLPvector* v);
  /* yE = v を解いて解(y)を v に返す. */
void eta_matrix_solve_Ex_is_equal_to_v(eta_matrix* E, EXLPvector* v);
  /* Ex = v を解いて解(x)を v に返す. */

eta_matrix_d* new_eta_matrix_d(int dimension);
void eta_matrix_d_free(eta_matrix_d* E);
void eta_matrix_d_mul_vector(eta_matrix_d* E, EXLPvector_d* v);
void eta_matrix_d_solve_yE_is_equal_to_v(eta_matrix_d* E, EXLPvector_d* v);
void eta_matrix_d_solve_Ex_is_equal_to_v(eta_matrix_d* E, EXLPvector_d* v);

eta_singleton* new_eta_singleton(void);
void eta_singleton_free(eta_singleton* es);
void eta_singleton_set(eta_singleton* es, mpq_t q, int row, int column);
eta_singleton_d* new_eta_singleton_d(void);
void eta_singleton_d_free(eta_singleton_d* es);
void eta_singleton_d_set(eta_singleton_d* es, double d, int row, int column);
void eta_singleton_get_d(eta_singleton* es, eta_singleton_d* es_d);

#endif
