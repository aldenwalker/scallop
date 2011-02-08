#include "hash.h"
#include "mylib.h"
#include <stdio.h>
#include <string.h>

void hash_mpq_init(hash_mpq* h, unsigned long keys) {
  unsigned long  i;

  h->keys = keys;
  h->operation = (int*)  my_malloc(keys*sizeof(int));
  h->operand1  = (mpq_t*)my_malloc(keys*sizeof(mpq_t));
  h->operand2  = (mpq_t*)my_malloc(keys*sizeof(mpq_t));
  h->result    = (mpq_t*)my_malloc(keys*sizeof(mpq_t));

  for (i = 0; i < keys; i ++)
    h->operation[i] = HASH_OPERATION_NOT_USED;
}

void hash_mpq_free(hash_mpq* h) {
  free(h->operation);
  free(h->operand1);
  free(h->operand2);
  free(h->result);
}

mpq_t* hash_mpq_find(hash_mpq* h,
                     int operation, mpq_t operand1, mpq_t operand2) {
  unsigned long  i;

  i  = ((unsigned long)mpz_get_si(mpq_numref(operand1))+101) % h->keys;
  i *= ((unsigned long)mpz_get_si(mpq_denref(operand1))+ 11) % h->keys;
  i *= ((unsigned long)mpz_get_si(mpq_numref(operand2))) % h->keys;
  i *= ((unsigned long)mpz_get_si(mpq_denref(operand2))) % h->keys;
  i %= h->keys;

  if (h->operation[i] == operation        &&
      mpq_equal(h->operand1[i], operand1) &&
      mpq_equal(h->operand2[i], operand2)) {
    //fprintf(stderr,"*");
    return &(h->result[i]);
  }
  //fprintf(stderr,".");

  if (h->operation[i] == HASH_OPERATION_NOT_USED) {
    mpq_init(h->operand1[i]);
    mpq_init(h->operand2[i]);
    mpq_init(h->result  [i]);
  }

  h->operation[i] = operation;
  mpq_set(h->operand1[i], operand1);
  mpq_set(h->operand2[i], operand2);

  if (operation == HASH_OPERATION_ADD)
    mpq_add(h->result[i], operand1, operand2);

  //else if (operation == HASH_OPERATION_SUB)
  //  mpq_sub(h->result[i], operand1, operand2);

  else //if (operation == HASH_OPERATION_MUL)
    mpq_mul(h->result[i], operand1, operand2);

  //else
  //  mpq_div(h->result[i], operand1, operand2);

  return &(h->result[i]);
}

hash_str *hash_str_init(unsigned long keys) {
  hash_str* h = my_malloc(sizeof(*h));

  h->keys = keys;
  h->num  = (int**)my_calloc(keys, sizeof(int*));
  h->x    = (int* )my_calloc(keys, sizeof(int));

  return h;
}

int  hash_str_find(hash_str* h, char* str, char** table, int add) {
  char* s;
  int  i, j;

  //for (i = 1, s = str; *s; s ++)
  //  i = (i**s)%h->keys;
  for (i = 0, s = str; *s; s ++)
    i = i*31 + *s;
  i %= h->keys;

  for (j = 0; j < h->x[i]; j ++) {
    if (strcmp(table[h->num[i][j]], str) == 0)
      return h->num[i][j];
  }

  if (add >= 0) {
    h->num[i]    = (int*)my_realloc(h->num[i], (j+1)*sizeof(int));
    h->num[i][j] = add;
    h->x[i] = j+1;
  }

  return -1;
}

void hash_str_free(hash_str** h) {
  int i, keys;
  if((h == 0) || (*h == 0))
    return;
  keys = (*h)->keys;
  for(i = 0; i < keys; i++) {
    if((*h)->num[i] != 0)
      free((*h)->num[i]);
  }
  free((*h)->num);
  free((*h)->x);
  free((*h));
  *h = NULL;
}
