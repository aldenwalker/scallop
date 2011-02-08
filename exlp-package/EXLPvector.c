#include "EXLPvector.h"
#include "mylib.h"
#include <string.h>


void vector_block_resize(EXLPvector* vec, int blocks) {
  int  i;

  if (blocks == vec->blocks)
    return;

  for (i = blocks*VECTOR_BLOCK_SIZE;
       i < vec->blocks*VECTOR_BLOCK_SIZE; i ++)
    mpq_clear(vec->value[i]);

  vec->i     = my_realloc(vec->i, VECTOR_BLOCK_SIZE*blocks*sizeof(int));
  vec->value = my_realloc(vec->value, VECTOR_BLOCK_SIZE*blocks*
                                      sizeof(mpq_t));

  for (i = vec->blocks*VECTOR_BLOCK_SIZE;
       i < blocks*VECTOR_BLOCK_SIZE; i ++)
    mpq_init(vec->value[i]);

  vec->blocks = blocks;
}

void vector_init(EXLPvector* vec, int dimension) {
  int  i;

  vec->nonzeros = 0;
  vec->i      = NULL;
  vec->value  = NULL;
  vec->blocks = 0;
  vector_block_resize(vec, 1);
  vec->dimension = dimension;
  vec->ptr = my_malloc((dimension/VECTOR_BLOCK_SIZE+1)*VECTOR_BLOCK_SIZE*
                       sizeof(int));
  for (i = 0; i < (dimension/VECTOR_BLOCK_SIZE+1)*VECTOR_BLOCK_SIZE; i ++)
    vec->ptr[i] = -1;
}

EXLPvector* new_vector(int dimension) {
  EXLPvector* vec;

  vec = my_malloc(sizeof(EXLPvector));
  vector_init(vec, dimension);
  mpq_init(vec->tmp);

  return vec;
}

void vector_free(EXLPvector** vec) {
  int  i;

  if (*vec == NULL)
    return;

  for (i = 0; i < (*vec)->blocks*VECTOR_BLOCK_SIZE; i ++) {
    mpq_clear((*vec)->value[i]);
  }

  free((*vec)->i);
  free((*vec)->value);
  free((*vec)->ptr);
  mpq_clear((*vec)->tmp);
  free(*vec);
  (*vec) = NULL;
}

void vector_resize(EXLPvector* vec, int dimension) {
  int  a, b;
  int  i;

  a = vec->dimension/VECTOR_BLOCK_SIZE;
  b =      dimension/VECTOR_BLOCK_SIZE;

  vec->dimension = dimension;

  if (a == b)
    return;

  vec->ptr = my_realloc(vec->ptr, (b+1)*VECTOR_BLOCK_SIZE*sizeof(int));

  for (i = (a+1)*VECTOR_BLOCK_SIZE; i < (b+1)*VECTOR_BLOCK_SIZE; i ++)
    vec->ptr[i] = -1;
}

void vector_copy(EXLPvector* to, EXLPvector* from) {
  int  i;

  vector_block_resize(to, from->blocks);

  to->nonzeros = from->nonzeros;

  for (i = 0; i < from->nonzeros; i ++) {
    to->i[i] = from->i[i];
    mpq_set(to->value[i], from->value[i]);
  }

  to->dimension = from->dimension;
  i = (from->dimension/VECTOR_BLOCK_SIZE+1)*VECTOR_BLOCK_SIZE;
  to->ptr = my_realloc(to->ptr, i*sizeof(int));
  memmove(to->ptr, from->ptr, i*sizeof(int));
  //to->ptr = my_realloc(to->ptr, sizeof(int)*from->dimension);
  //memmove(to->ptr, from->ptr, sizeof(int)*from->dimension);
}

void vector_zero_clear(EXLPvector* vec) {
  int  i;

  for (i = 0; i < vec->nonzeros; i ++)
    vec->ptr[vec->i[i]] = -1;

  vector_block_resize(vec, 1);
  vec->nonzeros = 0;
}

int vector_element_is_zero(EXLPvector* vec, int i) {

  return (((unsigned)(vec->ptr[i])) >= vec->nonzeros);
  // 午後のこーだのページに出てた最適化.
  // 負の数のデータの持ち方を 2 の補数と仮定してる.
  // 条件分岐が減るのがよい.
  //return (vec->ptr[i] < 0 || vec->ptr[i] >= vec->nonzeros);
}

int vector_get_nonzeros(EXLPvector* vec, int s, int t) {
  /* s 以上 t 以下の範囲にいくつの非ゼロ要素があるか */
  int  i, n;

  n = 0;

  if (t-s+1 < vec->nonzeros) {
    for (i = s; i <= t; i ++)
      n += (vector_element_is_zero(vec, i) == 0);

  } else {
    for (i = 0; i < vec->nonzeros; i ++)
      if (s <= vec->i[i] && vec->i[i] <= t)
	n ++;
  }

  return n;
}

int vector_get_first_nonzero_i(EXLPvector* vec) {
  int  i, ii;

  ii = vec->dimension+1;
  for (i = 0; i < vec->nonzeros; i ++)
    if (ii > vec->i[i])
      ii = vec->i[i];

  return ii;
}

int vector_get_last_nonzero_i(EXLPvector* vec) {
  int  i, ii;

  ii = -1;
  for (i = 0; i < vec->nonzeros; i ++)
    if (ii < vec->i[i])
      ii = vec->i[i];

  return ii;
}

int vector_get_nonzero_pos(EXLPvector* vec, int i) {

  if (vector_element_is_zero(vec, i))
    return -1;

  return vec->ptr[i];
}

void vector_rev_sgn(EXLPvector* v) {
  int pos;

  for (pos = 0; pos < v->nonzeros; pos ++)
    mpq_neg(v->value[pos], v->value[pos]);
}

void vector_rev_element_sgn(EXLPvector* v, int i) {

  if (vector_element_is_zero(v, i))
    return;

  mpq_neg(v->value[v->ptr[i]], v->value[v->ptr[i]]);
}

void vector_get_element(mpq_t* val, EXLPvector* vec, int i) {

  if (vector_element_is_zero(vec, i))
    mpq_set_ui(*val, 0, 1);
  else
    mpq_set(*val, vec->value[vec->ptr[i]]);
}

mpq_t* vector_get_element_ptr(EXLPvector* vec, int i) {

  if (vector_element_is_zero(vec, i))
    return &mympq_zero;
  else
    return &(vec->value[vec->ptr[i]]);
}

void vector_set_element(EXLPvector* vec, mpq_t val, int i) {

  if (mpq_sgn(val) == 0) {
    vector_delete_element(vec, i);
    return;
  }

  if (!vector_element_is_zero(vec, i)) {
    mpq_set(vec->value[vec->ptr[i]], val);
    return;
  }

  if (vec->nonzeros >= VECTOR_BLOCK_SIZE*vec->blocks) {
    vector_block_resize(vec, vec->blocks + 1);
  }

  mpq_set(vec->value[vec->nonzeros], val);
  vec->i[vec->nonzeros] = i;
  vec->ptr[i] = vec->nonzeros ++;
}

void vector_delete_element(EXLPvector* vec, int i) {
  mpq_t  tmp;

  if (vector_element_is_zero(vec, i))
    return;

  vec->nonzeros --;

  if (vec->ptr[i] == vec->nonzeros) {
    vec->ptr[i] = -1;
    return;
  }

  memmove(&tmp, &(vec->value[vec->ptr[i]]), sizeof(mpq_t));
  memmove(&(vec->value[vec->ptr[i]]), &(vec->value[vec->nonzeros]),
          sizeof(mpq_t));
  memmove(&(vec->value[vec->nonzeros]), &tmp, sizeof(mpq_t));
  vec->i[vec->ptr[i]] = vec->i[vec->nonzeros];

  vec->ptr[vec->i[vec->nonzeros]] = vec->ptr[i];
  vec->ptr[i] = -1;
}

void vector_add_element(EXLPvector* v, mpq_t q, int i) {

  if (vector_element_is_zero(v, i)) {
    vector_set_element(v, q, i);
    return;
  }

  mympq_add(v->value[v->ptr[i]], v->value[v->ptr[i]], q);

  if (mpq_sgn(v->value[v->ptr[i]]) == 0)
    vector_delete_element(v, i);
}

void vector_sub_element(EXLPvector* v, mpq_t q, int i) {

  if (vector_element_is_zero(v, i)) {
    mpq_neg(q, q);
    vector_set_element(v, q, i);
    mpq_neg(q, q);
    return;
  }

  mympq_sub(v->value[v->ptr[i]], v->value[v->ptr[i]], q);

  if (mpq_sgn(v->value[v->ptr[i]]) == 0)
    vector_delete_element(v, i);
}

void vector_mul_element(EXLPvector* v, mpq_t q, int i) {

  if (mpq_sgn(q) == 0) {
    vector_delete_element(v, i);
    return;
  }

  if (vector_element_is_zero(v, i))
    return;

  mympq_mul(v->value[v->ptr[i]], v->value[v->ptr[i]], q);
}

void vector_div_element(EXLPvector* v, mpq_t q, int i) {

  if (mpq_sgn(q) == 0) {
    fprintf(stderr, "division by 0\n");
    exit(EXIT_FAILURE);
  }

  if (vector_element_is_zero(v, i))
    return;

  mympq_div(v->value[v->ptr[i]], v->value[v->ptr[i]], q);
}

void vector_neg_element(EXLPvector* v, int i) {

  if (vector_element_is_zero(v, i))
    return;

  mpq_neg(v->value[v->ptr[i]], v->value[v->ptr[i]]);
}

void vector_dec_dimension(EXLPvector* v, int i) {
  int  j;

  vector_delete_element(v, i);

  for (j = 0; j < v->nonzeros; j ++) {
    if (v->i[j] > i) {
      v->ptr[v->i[j]] = -1;
      v->i[j] --;
    }
  }

  v->dimension --;

  for (j = 0; j < v->nonzeros; j ++)
    v->ptr[v->i[j]] = j;
}

void vector_scalar_product(EXLPvector* v, mpq_t s) {
  int  i;

  if (mpq_sgn(s) == 0) {
    vector_zero_clear(v);
    return;
  }

  if (mpq_cmp_si(s,  1, 1) == 0)
    return;
  if (mpq_cmp_si(s, -1, 1) == 0) {
    vector_rev_sgn(v);
    return;
  }

  for (i = 0; i < v->nonzeros; i ++)
    mympq_mul(v->value[i], v->value[i], s);
}

void vector_norm(mpq_t* val, EXLPvector* v) {
  int  i;

  mpq_set_ui(*val, 0, 1);
  for (i = 0; i < v->nonzeros; i ++) {
    mympq_mul(v->tmp, v->value[i], v->value[i]);
    mympq_add(*val, *val, v->tmp);
  }
}

/*
void vector_inner_product(mpq_t* val, EXLPvector* a, EXLPvector* b) {
  int  i, j;
  mpq_t** m1;
  mpq_t** m2;
  int* l;

  if (a->nonzeros > b->nonzeros) {
    vector_inner_product(val, b, a);
    return;
  }

  m1 = (mpq_t**)my_malloc((a->nonzeros+1) * sizeof(mpq_t*));
  m2 = (mpq_t**)my_malloc((a->nonzeros+1) * sizeof(mpq_t*));
  l = (int*)my_malloc((a->nonzeros+1) * sizeof(int));

  j = 0;
  for (i = 0; i < a->nonzeros; i ++) {
    if (vector_element_is_zero(b, a->i[i]))
      continue;
    m1[j] = &(a->value[i]);
    m2[j] = &(b->value[b->ptr[a->i[i]]]);
    l[j] = mpz_size(mpq_numref(*(m1[j]))) + mpz_size(mpq_denref(*(m1[j]))) +
           mpz_size(mpq_numref(*(m2[j]))) + mpz_size(mpq_denref(*(m2[j])));
    j ++;
  }

  if (j > 0)
    my_sort2(l, m1, m2, 0, j-1);

  mpq_set_ui(*val, 0, 1);

//fprintf(stderr,"<");
  for (i = 0; i < j; i ++) {

    mympq_mul(a->tmp, *(m1[i]), *(m2[i]));
    mympq_add(*val, *val, a->tmp);
//fprintf(stderr,"%d/%d ", mpz_size(mpq_numref(*val)), mpz_size(mpq_denref(*val)));
  }
//fprintf(stderr,">\n");

  free(m1);
  free(m2);
  free(l);
}
*/

void vector_inner_product(mpq_t* val, EXLPvector* a, EXLPvector* b) {
  int  i;

  if (a->nonzeros > b->nonzeros) {
    vector_inner_product(val, b, a);
    return;
  }

  mpq_set_ui(*val, 0, 1);

  //fprintf(stderr,"<");
  for (i = a->nonzeros-1; i >= 0; i --) {
    if (vector_element_is_zero(b, a->i[i]))
      continue;

    //fprintf(stderr,"%d/%d ", mpz_size(mpq_numref(a->value[i])), mpz_size(mpq_denref(a->value[i])));

    mympq_mul(a->tmp, a->value[i], b->value[b->ptr[a->i[i]]]);
    mympq_add(*val, *val, a->tmp);
  }
  //fprintf(stderr,">\n");
}

void vector_swap_elements(EXLPvector* vec, int i1, int i2) {
  int  tmp;

  if (i1 == i2)
    return;

  if (vector_element_is_zero(vec, i1)) {

    if (vector_element_is_zero(vec, i2))
      return;

    vec->i[vec->ptr[i2]] = i1;
    vec->ptr[i1] = vec->ptr[i2];
    vec->ptr[i2] = -1;
    return;
  }

  if (vector_element_is_zero(vec, i2)) {
    vec->i[vec->ptr[i1]] = i2;
    vec->ptr[i2] = vec->ptr[i1];
    vec->ptr[i1] = -1;
    return;
  }

  tmp = vec->ptr[i1];
  vec->ptr[i1] = vec->ptr[i2];
  vec->ptr[i2] = tmp;
  vec->i[vec->ptr[i1]] = i1;
  vec->i[vec->ptr[i2]] = i2;
}

void vector_print(EXLPvector* vec) {
  int   i;

  for (i = 0; i < vec->dimension; i ++) {
    print_rational_as_float(*vector_get_element_ptr(vec, i), 12);
    putchar('\n');
  }
}

int vector_check(EXLPvector* vec) {
/* for debuging */
  int  i;

  fprintf(stderr, "\nvector check...");
  for (i = 0; i < vec->nonzeros; i ++) {
    if (vec->i[i] < 0 || vec->i[i] >= vec->dimension ||
        vec->ptr[vec->i[i]] != i)
      fprintf(stderr, "aaaaa\n");
  }
  for (i = 0; i < vec->dimension; i ++) {
    if (0 <= vec->ptr[i] && vec->ptr[i] < vec->nonzeros &&
        vec->i[vec->ptr[i]] != i)
      fprintf(stderr, "bbbbb\n");
  }
  fprintf(stderr, "ok\n");
  return 0;
}


void vector_d_block_resize(EXLPvector_d* vec, int blocks) {

  if (blocks == vec->blocks)
    return;

  vec->i     = my_realloc(vec->i, VECTOR_BLOCK_SIZE*blocks*sizeof(int));
  vec->value = my_realloc(vec->value, VECTOR_BLOCK_SIZE*blocks*
                                      sizeof(double));

  vec->blocks = blocks;
}

void vector_d_init(EXLPvector_d* vec, int dimension) {
  int  i;

  vec->nonzeros = 0;
  vec->i      = NULL;
  vec->value  = NULL;
  vec->blocks = 0;
  vector_d_block_resize(vec, 1);
  vec->dimension = dimension;
  vec->ptr = my_malloc((dimension/VECTOR_BLOCK_SIZE+1)*VECTOR_BLOCK_SIZE*
                       sizeof(int));
  for (i = 0; i < (dimension/VECTOR_BLOCK_SIZE+1)*VECTOR_BLOCK_SIZE; i ++)
    vec->ptr[i] = -1;
}

EXLPvector_d* new_vector_d(int dimension) {
  EXLPvector_d* vec;

  vec = my_malloc(sizeof(EXLPvector_d));
  vector_d_init(vec, dimension);

  return vec;
}

void vector_d_free(EXLPvector_d** vec) {
  if ((*vec) == NULL)
    return;
  free((*vec)->i);
  free((*vec)->value);
  free((*vec)->ptr);
  free(*vec);
  (*vec) = NULL;
}

void vector_d_resize(EXLPvector_d* vec, int dimension) {
  int  a, b;
  int  i;

  a = vec->dimension/VECTOR_BLOCK_SIZE;
  b =      dimension/VECTOR_BLOCK_SIZE;

  vec->dimension = dimension;

  if (a == b)
    return;

  vec->ptr = my_realloc(vec->ptr, (b+1)*VECTOR_BLOCK_SIZE*sizeof(int));

  for (i = (a+1)*VECTOR_BLOCK_SIZE; i < (b+1)*VECTOR_BLOCK_SIZE; i ++)
    vec->ptr[i] = -1;
}

void vector_d_copy(EXLPvector_d* to, EXLPvector_d* from) {
  int  i;

  vector_d_block_resize(to, from->blocks);

  to->nonzeros = from->nonzeros;

  for (i = 0; i < from->nonzeros; i ++) {
    to->i[i] = from->i[i];
    to->value[i] = from->value[i];
  }

  to->dimension = from->dimension;
  i = (from->dimension/VECTOR_BLOCK_SIZE+1)*VECTOR_BLOCK_SIZE;
  to->ptr = my_realloc(to->ptr, i*sizeof(int));
  memmove(to->ptr, from->ptr, i*sizeof(int));
}

void vector_d_zero_clear(EXLPvector_d* vec) {
  int  i;

  for (i = 0; i < vec->nonzeros; i ++)
    vec->ptr[vec->i[i]] = -1;

  vector_d_block_resize(vec, 1);
  vec->nonzeros = 0;
}

void vector_get_d(EXLPvector* vec, EXLPvector_d* v_d) {
  int  i;

  vector_d_resize(v_d, vec->dimension);
  vector_d_block_resize(v_d, vec->blocks);
  v_d->nonzeros = vec->nonzeros;

  for (i = 0; i < vec->nonzeros; i ++) {
    v_d->i[i] = vec->i[i];
    v_d->value[i] = mpq_get_d(vec->value[i]);
  }
  memmove(v_d->ptr, vec->ptr,
          (vec->dimension/VECTOR_BLOCK_SIZE+1)*VECTOR_BLOCK_SIZE*sizeof(int));
}

int vector_d_element_is_zero(EXLPvector_d* vec, int i) {

  return (((unsigned)(vec->ptr[i])) >= vec->nonzeros);
  // 午後のこーだのページに出てた最適化.
  // 負の数のデータの持ち方を 2 の補数と仮定してる.
  //return (vec->ptr[i] < 0 || vec->ptr[i] >= vec->nonzeros);
}

int vector_d_get_nonzeros(EXLPvector_d* vec, int s, int t) {
  /* s 以上 t 以下の範囲にいくつの非ゼロ要素があるか */
  int  i, n;

  n = 0;

  if (t-s+1 < vec->nonzeros) {
    for (i = s; i <= t; i ++)
      n += (vector_d_element_is_zero(vec, i) == 0);

  } else {
    for (i = 0; i < vec->nonzeros; i ++)
      if (s <= vec->i[i] && vec->i[i] <= t)
	n ++;
  }

  return n;
}

double vector_d_get_element(EXLPvector_d* vec, int i) {

  if (vector_d_element_is_zero(vec, i))
    return 0;
  else
    return vec->value[vec->ptr[i]];
}

void vector_d_set_element(EXLPvector_d* vec, double val, int i) {

  if (val == 0) {
    vector_d_delete_element(vec, i);
    return;
  }

  if (!vector_d_element_is_zero(vec, i)) {
    vec->value[vec->ptr[i]] = val;
    return;
  }

  if (vec->nonzeros >= VECTOR_BLOCK_SIZE*vec->blocks) {
    vector_d_block_resize(vec, vec->blocks + 1);
  }

  vec->value[vec->nonzeros] = val;
  vec->i[vec->nonzeros] = i;
  vec->ptr[i] = vec->nonzeros ++;
}

void vector_d_delete_element(EXLPvector_d* vec, int i) {
  double  tmp;

  if (vector_d_element_is_zero(vec, i))
    return;

  vec->nonzeros --;

  tmp = vec->value[vec->ptr[i]];
  vec->value[vec->ptr[i]] = vec->value[vec->nonzeros];
  vec->value[vec->nonzeros] = tmp;
  vec->i[vec->ptr[i]] = vec->i[vec->nonzeros];

  vec->ptr[vec->i[vec->nonzeros]] = vec->ptr[i];
  vec->ptr[i] = -1;
}

void vector_d_add_element(EXLPvector_d* v, double d, int i) {

  if (vector_d_element_is_zero(v, i)) {
    vector_d_set_element(v, d, i);
    return;
  }

  v->value[v->ptr[i]] += d;

  if (v->value[v->ptr[i]] == 0)
    vector_d_delete_element(v, i);
}

void vector_d_sub_element(EXLPvector_d* v, double d, int i) {

  if (vector_d_element_is_zero(v, i)) {
    vector_d_set_element(v, -d, i);
    return;
  }

  v->value[v->ptr[i]] -= d;

  if (v->value[v->ptr[i]] == 0)
    vector_d_delete_element(v, i);
}

void vector_d_mul_element(EXLPvector_d* v, double d, int i) {

  if (d == 0) {
    vector_d_delete_element(v, i);
    return;
  }

  if (vector_d_element_is_zero(v, i))
    return;

  v->value[v->ptr[i]] *= d;
}

void vector_d_div_element(EXLPvector_d* v, double d, int i) {

  if (vector_d_element_is_zero(v, i))
    return;

  v->value[v->ptr[i]] /= d;
}

double vector_d_inner_product(EXLPvector_d* a, EXLPvector_d* b) {
  double  val;
  int  i;
  int* bp;
  double* aa;
  double* bb;

  if (a->nonzeros > b->nonzeros)
    return vector_d_inner_product(b, a);

  val = 0;
  bp = b->ptr;
  aa = a->value;
  bb = b->value;

  for (i = a->nonzeros-1; i >= 0; i --) {
    if (!vector_d_element_is_zero(b, a->i[i]))
      val += aa[i] * bb[bp[a->i[i]]];
  }

  return val;
}

void vector_d_swap_elements(EXLPvector_d* vec, int i1, int i2) {
  int  tmp;

  if (i1 == i2)
    return;

  if (vector_d_element_is_zero(vec, i1)) {

    if (vector_d_element_is_zero(vec, i2))
      return;

    vec->i[vec->ptr[i2]] = i1;
    vec->ptr[i1] = vec->ptr[i2];
    vec->ptr[i2] = -1;
    return;
  }

  if (vector_d_element_is_zero(vec, i2)) {
    vec->i[vec->ptr[i1]] = i2;
    vec->ptr[i2] = vec->ptr[i1];
    vec->ptr[i1] = -1;
    return;
  }

  tmp = vec->ptr[i1];
  vec->ptr[i1] = vec->ptr[i2];
  vec->ptr[i2] = tmp;
  vec->i[vec->ptr[i1]] = i1;
  vec->i[vec->ptr[i2]] = i2;
}

void vector_d_dec_dimension(EXLPvector_d* v, int i) {
  int  j;

  vector_d_delete_element(v, i);

  for (j = 0; j < v->nonzeros; j ++) {
    if (v->i[j] > i) {
      v->ptr[v->i[j]] = -1;
      v->i[j] --;
    }
  }

  v->dimension --;

  for (j = 0; j < v->nonzeros; j ++)
    v->ptr[v->i[j]] = j;
}
