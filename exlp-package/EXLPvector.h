#ifndef VECTOR_H
#define VECTOR_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <stdio.h>
#include <gmp.h>

#define VECTOR_BLOCK_SIZE 16

typedef struct {
  int  dimension;
  int  blocks;
  int  nonzeros;
  int* i;
  mpq_t* value;
  int* ptr;
  mpq_t  tmp;
} EXLPvector;

typedef struct {
  int  dimension;
  int  blocks;
  int  nonzeros;
  int* i;
  double* value;
  int* ptr;
} EXLPvector_d;

void vector_init(EXLPvector* vec, int dimension);
EXLPvector* new_vector(int dimension);
void vector_free(EXLPvector** vec);
void vector_resize(EXLPvector* vec, int dimension);
void vector_copy(EXLPvector* to, EXLPvector* from);
void vector_zero_clear(EXLPvector* vec);
int  vector_element_is_zero(EXLPvector* vec, int i);
int  vector_get_nonzeros(EXLPvector* vec, int s, int t);
int  vector_get_first_nonzero_i(EXLPvector* vec);
int  vector_get_last_nonzero_i(EXLPvector* vec);
int  vector_get_nonzero_pos(EXLPvector* vec, int i);
void vector_rev_sgn(EXLPvector* v);
void vector_rev_element_sgn(EXLPvector* v, int i);
void vector_get_element(mpq_t* val, EXLPvector* vec, int i);
mpq_t* vector_get_element_ptr(EXLPvector* vec, int i);
void vector_set_element(EXLPvector* vec, mpq_t val, int i);
void vector_delete_element(EXLPvector* vec, int i);
void vector_add_element(EXLPvector* vec, mpq_t q, int i);
void vector_sub_element(EXLPvector* vec, mpq_t q, int i);
void vector_mul_element(EXLPvector* vec, mpq_t q, int i);
void vector_div_element(EXLPvector* vec, mpq_t q, int i);
void vector_neg_element(EXLPvector* vec, int i);
void vector_dec_dimension(EXLPvector* vec, int i);
void vector_scalar_product(EXLPvector* v, mpq_t s);
void vector_norm(mpq_t* val, EXLPvector* v);
void vector_inner_product(mpq_t* val, EXLPvector* a, EXLPvector* b);
void vector_swap_elements(EXLPvector* vec, int i1, int i2);
void vector_print(EXLPvector* vec);
int vector_check(EXLPvector* vec);

EXLPvector_d* new_vector_d(int dimension);
void vector_d_free(EXLPvector_d** v_d);
void vector_d_resize(EXLPvector_d* vec, int dimension);
int  vector_d_element_is_zero(EXLPvector_d* vec, int i);
int vector_d_get_nonzeros(EXLPvector_d* vec, int s, int t);
void vector_d_copy(EXLPvector_d* to, EXLPvector_d* from);
void vector_d_zero_clear(EXLPvector_d* vec);
void vector_get_d(EXLPvector* vec, EXLPvector_d* v_d);
double vector_d_get_element(EXLPvector_d* v_d, int i);
void vector_d_set_element(EXLPvector_d* vec, double val, int i);
void vector_d_delete_element(EXLPvector_d* vec, int i);
void vector_d_add_element(EXLPvector_d* vec, double d, int i);
void vector_d_sub_element(EXLPvector_d* vec, double d, int i);
void vector_d_mul_element(EXLPvector_d* v, double d, int i);
void vector_d_div_element(EXLPvector_d* v, double d, int i);
double vector_d_inner_product(EXLPvector_d* a, EXLPvector_d* b);
void vector_d_swap_elements(EXLPvector_d* vec, int i1, int i2);
void vector_d_dec_dimension(EXLPvector_d* v, int i);

#endif
