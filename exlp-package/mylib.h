#ifndef MY_LIB_H
#define MY_LIB_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <ctype.h>
#include "hash.h"

#ifndef TRUE
#  define TRUE 1
#endif
#ifndef FALSE
#  define FALSE 0
#endif
#ifndef MIN
#  define MIN(a, b) (((a)<(b))?(a):(b))
#endif
#ifndef MAX
#  define MAX(a, b) (((a)>(b))?(a):(b))
#endif
#ifndef ABS
#  define ABS(a) (MAX((a),(-(a))))
#endif

void  mylib_init(void);
int   my_sgn(int i);
void* my_malloc(size_t size);
void* my_realloc(void* ptr, size_t size);
void* my_calloc(size_t n, size_t size);
void  my_sort(int* I, mpq_t* V, int min, int max);
void  my_sort2(int* I, mpq_t** V1, mpq_t** V2, int min, int max);
void  my_sort_d(int* I, double* V, int min, int max);
void  my_sort_d2(int* T, double* D, int min, int max);
void  print_rational(mpq_t q);
void  print_rational_as_float(mpq_t q, int n);
int mympq_cmp_abs(mpq_t a, mpq_t b);
void mympq_floor(mpq_t c, mpq_t q);
  /* c �� q ��ۤ��ʤ�������������֤�. */
void mympq_ceil(mpq_t c, mpq_t q);
  /* c �� q ��꾮�����ʤ��Ǿ����������֤�. */
void mympq_remainder(mpq_t r, mpq_t q);
  /* ��ʬ�� q ����ʬ�� ��ɽ�������ο�ʬ������ q ���֤�. r >= 0 */
int mympq_is_integer(mpq_t q);
void  mympq_simplify(mpq_t to, mpq_t from, int exactness);
  /* exactness ���Ϣʬ��Ÿ���򤹤� ^_^; */
void  mympq_set_float_string(mpq_t q, char* s);
void  mympq_set_string(mpq_t q, char* s);

static char mpq_const_defined = 0;
extern mpq_t mympq_zero;
extern mpq_t mympq_one;
extern mpq_t mympq_minus_one;
extern mpq_t mympq_tmp;

#ifdef NO_GMP_HASH

#define my_hash_mpq_init(n)
#define my_hash_mpq_free()
#define mympq_add mpq_add
#define mympq_sub mpq_sub
#define mympq_mul mpq_mul
#define mympq_div mpq_div

#else

void  my_hash_mpq_init(int hash_entries);
void my_hash_mpq_free(void);
void  mympq_add(mpq_t a, mpq_t b, mpq_t c);
void  mympq_sub(mpq_t a, mpq_t b, mpq_t c);
void  mympq_mul(mpq_t a, mpq_t b, mpq_t c);
void  mympq_div(mpq_t a, mpq_t b, mpq_t c);

#endif

#endif
