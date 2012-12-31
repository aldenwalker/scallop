#ifndef ETA_FILE_H
#define ETA_FILE_H

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <time.h>
#include "matrix.h"

#define ETA_MAX 128

//typedef long long my_clock_t;
typedef clock_t my_clock_t;

typedef struct {
  int lp_rows;
  int k;  // update した回数
  int s;  // Ls の数
  int *P;
  eta_matrix** L;
  eta_singleton** Ls;
  eta_matrix_d** L_d;
  eta_singleton_d** Ls_d;
  matrix* U;
  EXLPvector_d** U_d;
  permutation_matrix* Q;
  permutation_matrix* R;
  my_clock_t init_time;
  my_clock_t last_time;
} eta_file;

eta_file* new_eta_file(int lp_rows);
void eta_file_free(eta_file* ef);
void eta_file_btran(eta_file* ef, EXLPvector* c, EXLPvector* y);
void eta_file_ftran(eta_file* ef, EXLPvector* b, EXLPvector* x, EXLPvector* w);

void eta_file_btran_d(eta_file* ef, EXLPvector_d* c, EXLPvector_d* y);
void eta_file_ftran_d(eta_file* ef, EXLPvector_d* b, EXLPvector_d* x);

#endif
