#include "mps.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "mylib.h"
#include "EXLPvector.h"
#include "matrix.h"
#include "lpstruct.h"
#include "solve_lp.h"

int get_word(char* buf, FILE* fp) {
  int  i, c, r;

  r = 0;
  do {
    c = fgetc(fp);
    if (c == '\n')
      r = 1;
  } while (isspace(c));

  for (i = 0; i < LP_NAME_LEN_MAX-1; i ++) {
    if (c == EOF) {
      buf[i] = 0;
      return EOF;
    }
    if (isspace(c)) {
      ungetc(c, fp);
      buf[i] = 0;
      return r;
    }
    buf[i] = c;
    c = fgetc(fp);
  }
  buf[i] = 0;
  return r | 16;
}

void get_until_new_line(char* buf, FILE* fp) {
  int  c, i;

  i = 0;
  while ((c = fgetc(fp)) != EOF && c != '\n'
         && i < LP_NAME_LEN_MAX-1)
    buf[i++] = c;
  buf[i] = 0;
}

void get_first_word(char* buf, FILE* fp) {
NextLine:
  get_word(buf, fp);
  if (buf[0] == '*') {
    get_until_new_line(buf, fp);
    goto NextLine;
  }
}

int trim_word(char* buf) {
  size_t  r = 0, s = strlen(buf);

  while ((r < s) && isspace(buf[r]))
    r++;
  s--;
  while ((s > r) && isspace(buf[s]))
    s--;
  if (s > 0)
    memcpy(buf, buf+r, (s-r+1)*sizeof(char));
  s = s-r+1;
  buf[s] = '\0';
  return((int) s);
}

int reducedigits(char *buf, int kl) {
  int kt, 
      k0 = kl - 1;

  while ((k0 >= 0) && (buf[k0] == '0'))
    k0 --;
  kt = k0;
  while ((kt >= 0) && (buf[kt] != '.'))
    kt --;
  if ((kt > 0) && (buf[kt] == '.')) {
    kl = k0;
    if (kt < k0)
      kl ++;
    buf[kl] = '\0';
  }
  return kl;
}

void read_name(LP* lp, char* buf) {
  get_first_word(buf, lp->fp);
  if (strncmp(buf, "NAME", LP_NAME_LEN_MAX) == 0) {
    get_until_new_line(buf, lp->fp);
    trim_word(buf);
    lp_set_name(lp, buf);
    get_word(buf, lp->fp);
  }
}

void read_rows(LP* lp, char* buf) {
  int  c;

  if (strncmp(buf, "ROWS", LP_NAME_LEN_MAX) != 0) {
    fprintf(stderr, "read error: no row detected.\n");
    exit(1);
  }
  get_word(buf, lp->fp);
  while (strncmp(buf, "N", LP_NAME_LEN_MAX) == 0 ||
         strncmp(buf, "L", LP_NAME_LEN_MAX) == 0 ||
         strncmp(buf, "G", LP_NAME_LEN_MAX) == 0 ||
         strncmp(buf, "E", LP_NAME_LEN_MAX) == 0) {
    c = buf[0];
    get_word(buf, lp->fp);
    if (c == 'N')
      lp_set_obj_name(lp, buf);
    else {
      lp_add_row(lp, buf);
      lp_set_row_equality(lp, lp_get_row_num(lp, buf), c);
    }
    get_word(buf, lp->fp);
  }
}

void read_columns(LP* lp, char* buf) {
  char **buf3 = NULL;
  int  *r = NULL, *v = NULL;
  int  i, kv, kr, 
       kl, nz;

  if (strncmp(buf, "COLUMNS", LP_NAME_LEN_MAX) != 0) {
    fprintf(stderr, "read error: no column detected.\n");
    exit(1);
  }

  nz = 3*lp->rows;
  r = my_realloc(r, (nz+1)*sizeof(int));
  v = my_realloc(v, (nz+1)*sizeof(int));
  buf3 = my_realloc(buf3, (nz+1)*sizeof(char*));

  i = 0;
  get_word(buf, lp->fp);
  while (strncmp(buf, "RHS", LP_NAME_LEN_MAX) != 0) {

    kv = lp_get_var_num(lp, buf);
    if (kv < 0)
      kv = lp_add_var_without_A(lp, buf);

    get_word(buf, lp->fp);
    do {
      kr = lp_get_row_num(lp, buf);
      if (kr < 0) {
        fprintf(stderr, "row '%s' is not defined.\n", buf);
        exit(EXIT_FAILURE);
      }
      get_word(buf, lp->fp);

      if (i == nz) {
        nz += nz / 3;
        r = my_realloc(r, (nz+1)*sizeof(int));
        v = my_realloc(v, (nz+1)*sizeof(int));
        buf3 = my_realloc(buf3, (nz+1)*sizeof(char*));
      }
      /* Simplify (to conserve memory and increase parsing speed) */
      kl = reducedigits(buf, (int) strlen(buf));
      /* Add nz element to triplet arrays (consider using mpq_set_d() instead) */
      buf3[i] = my_malloc(++kl);
      strncpy(buf3[i], buf, kl);
      v[i] = kv;
      r[i] = kr;
      i ++;
    } while ((get_word(buf, lp->fp) & 15) != 1);
  }

  matrix_resize(lp->A, lp->rows, lp->vars);
  vector_resize(lp->b, lp->rows);
  vector_resize(lp->xb, lp->rows);
  vector_resize(lp->cb, lp->rows);

  nz = i;
  for (i = 0; i < nz; i ++) {
    //mympq_set_float_string(lp->q_work, buf3[i]);
    mympq_set_string(lp->q_work, buf3[i]);
    lp_set_coefficient(lp, lp->q_work, r[i], v[i]);
    free(buf3[i]);
  }
  free(buf3);
  free(r);
  free(v);

  if (lp->maximize == FALSE)
    vector_rev_sgn(lp->c);
}

void read_rhs(LP* lp, char* buf) {
  int   row;

  get_word(buf, lp->fp);
  while (strncmp(buf, "RANGES", LP_NAME_LEN_MAX) != 0 &&
         strncmp(buf, "BOUNDS", LP_NAME_LEN_MAX) != 0 &&
         strncmp(buf, "ENDATA", LP_NAME_LEN_MAX) != 0 &&
         buf[0] != 0) {

    /* Processing of the label should be here, 
       ...
       between these comment lines! */

    row = lp_get_row_num(lp, buf);  // This is for blend.mps
    if (row < 0) {                  // which doesn't have label field.
      get_word(buf, lp->fp);
      row = lp_get_row_num(lp, buf);
    }

    for (;;) {
      if (row < 0) {
        fprintf(stderr, "row '%s' is not defined.\n", buf);
        exit(EXIT_FAILURE);
      }
      get_word(buf, lp->fp);
      reducedigits(buf, (int) strlen(buf));
      //mympq_set_float_string(lp->q_work, buf);
      mympq_set_string(lp->q_work, buf);
      lp_set_rhs(lp, row, lp->q_work);
      if ((get_word(buf, lp->fp) & 15) == 1)
	break;
      row = lp_get_row_num(lp, buf);
    }
  }

}

void read_ranges(LP* lp, char* buf) {
  int  row;

  if (strcmp(buf, "RANGES") != 0)
    return;

  get_word(buf, lp->fp);
  while (strncmp(buf, "BOUNDS", LP_NAME_LEN_MAX) != 0 &&
         strncmp(buf, "ENDATA", LP_NAME_LEN_MAX) != 0 &&
         buf[0] != 0) {

    /* Processing of the label should be here, 
       ...
       between these comment lines! */
    get_word(buf, lp->fp);
    do {
      row = lp_get_row_num(lp, buf);
      if (row < 0) {
        fprintf(stderr, "row '%s' is not defined.\n", buf);
        exit(EXIT_FAILURE);
      }
      get_word(buf, lp->fp);
      reducedigits(buf, (int) strlen(buf));
      //mympq_set_float_string(lp->q_work, buf);
      mympq_set_string(lp->q_work, buf);
      lp_set_row_range(lp, row, lp->q_work);
    } while ((get_word(buf, lp->fp) & 15) != 1);
  }

}

void read_bounds(LP* lp, char* buf) {
  char  buf2[LP_NAME_LEN_MAX];
  char  buf3[LP_NAME_LEN_MAX];
  int  var;
  mpq_t  q;

  if (strcmp(buf, "BOUNDS") != 0)
    return;

  mpq_init(q);

  get_word(buf,  lp->fp);
  while (strcmp(buf, "ENDATA") != 0 && buf[0] != 0) {
    get_word(buf2, lp->fp);
    get_word(buf3, lp->fp);
    if (lp_get_var_num(lp, buf3) >= 0) {
      strcpy(buf2, buf3);
      get_word(buf3, lp->fp);
    }
    var = lp_get_var_num(lp, buf2);
    if (var < 0) {
      fprintf(stderr, "variable '%s' not defined.\n", buf2);
      exit(EXIT_FAILURE);
    }
    if (strcmp(buf, "FR") == 0) {
      if (lp->lower.is_valid[var])
        mpq_clear(lp->lower.bound[var]);
      if (lp->upper.is_valid[var])
        mpq_clear(lp->upper.bound[var]);
      lp->lower.is_valid[var] = lp->upper.is_valid[var] = FALSE;
      strcpy(buf, buf3);
      continue;
    } 
    else if (strcmp(buf, "MI") == 0) {
      if (lp->lower.is_valid[var])
        mpq_clear(lp->lower.bound[var]);
      lp->lower.is_valid[var] = FALSE;
      strcpy(buf, buf3);
      continue;
    } 
    else if (strcmp(buf, "PL") == 0) {
      if (lp->upper.is_valid[var])
        mpq_clear(lp->upper.bound[var]);
      lp->upper.is_valid[var] = FALSE;
      strcpy(buf, buf3);
      continue;
    }
    //mympq_set_float_string(q, buf3);
    mympq_set_string(q, buf3);

    if (strcmp(buf, "LO") == 0) {
      if (lp->lower.is_valid[var] == FALSE)
        mpq_init(lp->lower.bound[var]);
      mpq_set(lp->lower.bound[var], q);
      lp->lower.is_valid[var] = TRUE;
    } 
    else if (strcmp(buf, "UP") == 0) {
      if (lp->upper.is_valid[var] == FALSE)
        mpq_init(lp->upper.bound[var]);
      mpq_set(lp->upper.bound[var], q);
      lp->upper.is_valid[var] = TRUE;
    } 
    else if (strcmp(buf, "FX") == 0) {
      if (lp->lower.is_valid[var] == FALSE)
        mpq_init(lp->lower.bound[var]);
      mpq_set(lp->lower.bound[var], q);
      lp->lower.is_valid[var] = TRUE;
      if (lp->upper.is_valid[var] == FALSE)
        mpq_init(lp->upper.bound[var]);
      mpq_set(lp->upper.bound[var], q);
      lp->upper.is_valid[var] = TRUE;
    } 
    else {
      fprintf(stderr, "bound type '%s' is not supported.\n", buf);
      break;
    }
    get_word(buf,  lp->fp);
  }

  mpq_clear(q);
}

void lp_read_mps(LP* lp) {
  char buf[LP_NAME_LEN_MAX];

  lp_hash_str_init(lp, lp->hash_entries);
#ifndef NO_GMP_HASH
  my_hash_mpq_init(lp->hash_entries);
#endif

  read_name(lp, buf);
  read_rows(lp, buf);
  read_columns(lp, buf);
  read_rhs(lp, buf);
  read_ranges(lp, buf);
  read_bounds(lp, buf);

  fclose(lp->fp);
}
