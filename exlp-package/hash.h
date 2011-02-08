#ifndef HASH_H
#define HASH_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <stdlib.h>
#include <gmp.h>

#define HASH_OPERATION_NOT_USED 0
#define HASH_OPERATION_ADD      4
#define HASH_OPERATION_MUL      6
#define HASH_OPERATION_SUB      5
#define HASH_OPERATION_DIV      7

typedef struct {
  unsigned long keys;
  int*   operation;
  mpq_t* operand1;
  mpq_t* operand2;
  mpq_t* result;
} hash_mpq;

void   hash_mpq_init(hash_mpq* h, unsigned long keys);
void   hash_mpq_free(hash_mpq* h);
mpq_t* hash_mpq_find(hash_mpq* h,
                     int operation, mpq_t operand1, mpq_t operand2);

typedef struct {
  unsigned long keys;
  int** num;
  int* x;
} hash_str;

hash_str *hash_str_init(unsigned long keys);
int  hash_str_find(hash_str* h, char* str, char** table, int add);
void hash_str_free(hash_str** h);

#endif
