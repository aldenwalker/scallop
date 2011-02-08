#include "mylib.h"
#include <string.h>

void mylib_init(void) {
  if(mpq_const_defined)
    return;
  mpq_init(mympq_zero);
  mpq_init(mympq_one);
  mpq_init(mympq_minus_one);
  mpq_init(mympq_tmp);

  mpq_set_si(mympq_one, 1, 1);
  mpq_set_si(mympq_minus_one, -1, 1);
  mpq_const_defined = 1;
}

int my_sgn(int i) {
  if (i == 0)
    return 0;
  if (i > 0)
    return 1;
  return -1;
}

void* my_malloc(size_t size) {
  void* x;

  if (size <= 0)
    size = 1;

  if (!(x = malloc(size))) {
    fprintf(stderr, "can't alloc the memory.\n");
    exit(EXIT_FAILURE);
  }
  return x;
}

void* my_realloc(void* x, size_t size) {
  void* y;

  if (size <= 0)
    size = 1;

  if (!(y = realloc(x, size))) {
    fprintf(stderr, "can't alloc the memory.\n");
    exit(EXIT_FAILURE);
  }
  return y;
}

void* my_calloc(size_t n, size_t size) {
  void* x;

  if (n <= 0)
    n = 1;

  if (!(x = calloc(n, size))) {
    fprintf(stderr, "can't alloc the memory.\n");
    exit(EXIT_FAILURE);
  }
  return x;
}

void my_sort(int* I, mpq_t* V, int min, int max) {
  /* I[min] から I[max] の値を小さい順に並べ変える.
     その際 V も並べ変えられる. */
  int  i, j, t, x;
  char d[sizeof(mpq_t)];

  x = I[(min + max)/2];
  i = min;
  j = max;
  for ( ; ; ) {
    while (I[i] < x)
      i ++;
    while (x < I[j])
      j --;
    if (i >= j)
      break;
    t = I[i]; I[i] = I[j]; I[j] = t;
    memcpy((void*)d, &(V[i]), sizeof(mpq_t));
    memcpy(&(V[i]), &(V[j]), sizeof(mpq_t));
    memcpy(&(V[j]), (void*)d, sizeof(mpq_t));
    i ++; j --;
  }
  if (min < i-1)
    my_sort(I, V, min, i-1);
  if (j+1 < max)
    my_sort(I, V, j+1, max);
}

void my_sort2(int* I, mpq_t** V1, mpq_t** V2, int min, int max) {
  /* I[min] から I[max] の値を小さい順に並べ変える.
     その際 V1, V2 も並べ変えられる. */
  int  i, j, t, x;
  mpq_t* q;

  x = I[(min + max)/2];
  i = min;
  j = max;
  for ( ; ; ) {
    while (I[i] < x)
      i ++;
    while (x < I[j])
      j --;
    if (i >= j)
      break;
    t =  I[i];  I[i] =  I[j];  I[j] = t;
    //q = V1[i]; V1[i] = V1[j]; V1[j] = q;
    //q = V2[i]; V2[i] = V2[j]; V2[j] = q;
    i ++; j --;
  }
  if (min < i-1)
    my_sort2(I, V1, V2, min, i-1);
  if (j+1 < max)
    my_sort2(I, V1, V2, j+1, max);
}

void my_sort_d(int* I, double* V, int min, int max) {
  /* I[min] から I[max] の値を小さい順に並べ変える.
     その際 V も並べ変えられる. */
  int  i, j, t, x;
  double d;

  x = I[(min + max)/2];
  i = min;
  j = max;
  for ( ; ; ) {
    while (I[i] < x)
      i ++;
    while (x < I[j])
      j --;
    if (i >= j)
      break;
    t = I[i]; I[i] = I[j]; I[j] = t;
    d = V[i]; V[i] = V[j]; V[j] = d;
    i ++; j --;
  }
  if (min < i-1)
    my_sort_d(I, V, min, i-1);
  if (j+1 < max)
    my_sort_d(I, V, j+1, max);
}

void my_sort_d2(int* T, double* D, int min, int max) {
  /* D[min] から D[max] の値を大きい順に並べ変える.
     その際 T も並べ変えられる. */
  int  i, j, t;
  double  x, d;

  x = D[(min + max)/2];
  i = min;
  j = max;
  for ( ; ; ) {
    while (D[i] > x)
      i ++;
    while (x > D[j])
      j --;
    if (i >= j)
      break;
    t = T[i]; T[i] = T[j]; T[j] = t;
    d = D[i]; D[i] = D[j]; D[j] = d;
    i ++; j --;
  }
  if (min < i-1)
    my_sort_d2(T, D, min, i-1);
  if (j+1 < max)
    my_sort_d2(T, D, j+1, max);
}

void print_rational(mpq_t q) {
  mpz_t  z;

  mpz_init(z);

  mpq_get_num(z, q);
  mpz_out_str(stdout, 10, z);
  putchar('/');
  mpq_get_den(z, q);
  mpz_out_str(stdout, 10, z);

  mpz_clear(z);
}

void print_rational_as_float(mpq_t q, int n) {
  /* q を有効数字 n 桁まで表示. ただし整数部は全部表示. */
  /* gmp でやってくれ */
  mpz_t  tmp, num, den;
  int    i;

  mpz_init(tmp);
  mpz_init(num);
  mpz_init(den);

  mpq_get_num(num, q);
  mpq_get_den(den, q);
  if (mpq_sgn(q) < 0) {
    gmp_printf("-");
    mpq_neg(q, q);
    mpq_get_num(num, q);
    mpq_get_den(den, q);
    mpq_neg(q, q);
  }

  mpz_div(tmp, num, den);
  i = mpz_out_str(stdout, 10, tmp);
  if (mpz_cmp(num, den) >= 0)
    n -= i;

  if (n > 0)
    putchar('.');
  while (n -- > 0) {
    mpz_mod(num, num, den);
    mpz_mul_ui(num, num, 10);
    mpz_div(tmp, num, den);
    gmp_printf("%Zd", tmp);
  }

  mpz_clear(tmp);
  mpz_clear(num);
  mpz_clear(den);
}

int mympq_cmp_abs(mpq_t a, mpq_t b) {
  int  ret;

  if (mpq_sgn(a) < 0) {
    if (mpq_sgn(b) < 0)
      return mpq_cmp(b, a);
    mpq_neg(a, a);  // a < 0 で b >= 0 なので &a != &b なの.
    ret = mpq_cmp(a, b);
    mpq_neg(a, a);
    return ret;
  }
  if (mpq_sgn(b) < 0) {
    mpq_neg(b, b);
    ret = mpq_cmp(a, b);
    mpq_neg(b, b);
    return ret;
  }
  return mpq_cmp(a, b);
}

void mympq_floor(mpq_t c, mpq_t q) {
  /* c に q を越えない最大の整数を返す. */
  mpz_t  z;

  mpz_init(z);

  mpz_fdiv_q(z, mpq_numref(q), mpq_denref(q));
  mpq_set_z(c, z);

  mpz_clear(z);
}

void mympq_ceil(mpq_t c, mpq_t q) {
  /* c に q より小さくない最小の整数を返す. */
  mpz_t  z;

  mpz_init(z);

  mpz_cdiv_q(z, mpq_numref(q), mpq_denref(q));
  mpq_set_z(c, z);

  mpz_clear(z);
}

void mympq_remainder(mpq_t r, mpq_t q) {
  /* 仮分数 q を帯分数 で表した時の真分数部を q に返す. r >= 0 */
  mpz_t  z;

  mpz_init(z);

  mpz_fdiv_r(z, mpq_numref(q), mpq_denref(q));
  mpq_set_z(r, z);
  mpq_set_den(r, mpq_denref(q));

  mpz_clear(z);
}

int mympq_is_integer(mpq_t q) {
  // もしかするとまずい. 常に既約になってると仮定してるんで
  return (mpz_cmp_si(mpq_denref(q), 1) == 0);
}

void  mympq_simplify(mpq_t to, mpq_t from, int exactness) {
  /* exactness 桁の連分数展開をする ^_^; */
  mpz_t* b;
  int  i, n;
  int  sgn;
  mpq_t  q1, q2;
  mpz_t  z;
  int  popcount;

  if (exactness <= 0 || mpq_sgn(from) == 0) {
    mpq_set(to, from);
    return;
  }

  mpq_init(q1);
  mpq_init(q2);
  mpz_init(z);
  mpq_set(q1, from);

  sgn = 1;
  if (mpq_sgn(q1) < 0) {
    sgn = -1;
    mpq_neg(q1, q1);
  }
  popcount = mpz_popcount(mpq_numref(q1)) + mpz_popcount(mpq_denref(q1));

  b = my_malloc(sizeof(mpz_t)*exactness);
  for (i = 0; i < exactness; i ++)
    mpz_init(b[i]);

  mpz_div(b[0], mpq_numref(q1), mpq_denref(q1));
  n = exactness;
  for (i = 1; i < n; i ++) {
    mpq_set_z(q2, b[i-1]);
    mympq_sub(q1, q1, q2);
    if (mpq_sgn(q1) == 0) {
      n = i;
      break;
    }
    mpq_inv(q1, q1);
    if(mpq_sgn(q1) < 0)printf("!!!!!\n"),fflush(stdout);
    mpz_div(b[i], mpq_numref(q1), mpq_denref(q1));
  }

  mpq_set_z(to, b[n-1]);
  for (i = n-2; i >= 0; i --) {
    mpz_set(z, b[i]);
    mpz_mul(z, z, mpq_numref(to));
    mpz_add(z, z, mpq_denref(to));
    mpq_set_den(to, mpq_numref(to));
    mpq_set_num(to, z);
  }
  mpq_canonicalize(to);
  if (mpz_popcount(mpq_numref(to)) + mpz_popcount(mpq_denref(to)) > popcount) {
    mpq_set(to, from);
    sgn = 1;
  }

  if (sgn < 0)
    mpq_neg(to, to);

  for (i = 0; i < exactness; i ++)
    mpz_clear(b[i]);
  mpq_clear(q1);
  mpq_clear(q2);
  mpz_clear(z);
}

void mympq_set_float_string(mpq_t q, char* s) {
  mpz_t  num, den;
  int  sgn;

  mpz_init_set_ui(num, 0);
  mpz_init_set_ui(den, 1);

  while (isspace(*s))
    s ++;
  if (*s == '-') {
    sgn = -1;
    s ++;
  } else {
    sgn = 1;
    if (*s == '+')
      s ++;
  }
  while (isspace(*s))
    s ++;

  while (isdigit(*s)) {
    mpz_mul_ui(num, num, 10);
    mpz_add_ui(num, num, *s-'0');
    s ++;
  }

  if (*s == '.') {
    s ++;
    while (isdigit(*s)) {
      mpz_mul_ui(num, num, 10);
      mpz_mul_ui(den, den, 10);
      mpz_add_ui(num, num, *s-'0');
      s ++;
    }
  }

  if (*s == 'e' || *s == 'E') {
    mpz_t  E;
    int  sgn2;
    unsigned int  e;

    mpz_init(E);
    mpz_set_si(E, 10);

    s ++;
    sgn2 = 1;
    if (*s == '+') {
      s ++;
    } else if (*s == '-') {
      s ++;
      sgn2 = -1;
    }
    for (e = 0; isdigit(*s); s ++) {
      e = e*10 + *s-'0';
    }
    mpz_pow_ui(E, E, e);
    if (sgn2 < 0)
      mpz_mul(den, den, E);
    else
      mpz_mul(num, num, E);

    mpz_clear(E);
  }

  mpq_set_num(q, num);
  mpq_set_den(q, den);

  mpq_canonicalize(q);

  if (sgn == -1)
    mpq_neg(q, q);
}

void mympq_set_string(mpq_t q, char* s) {
/* -1/4 とか +3.2/-7.3 とかを認めると. */
  mpq_t  tmp;

  mympq_set_float_string(q, s);

  while (*s != '/' && *s != 0)
    s ++;

  if (*s == '/') {
    mpq_init(tmp);
    mympq_set_float_string(tmp, ++s);
    mympq_div(q, q, tmp);
    mpq_clear(tmp);
  }
}

#ifndef NO_GMP_HASH

static hash_mpq my_hash_add;
static hash_mpq my_hash_mul;

void my_hash_mpq_init(int hash_entries) {
  hash_mpq_init(&my_hash_add, hash_entries);
  hash_mpq_init(&my_hash_mul, hash_entries);
}

void my_hash_mpq_free(void) {
  hash_mpq_free(&my_hash_add);
  hash_mpq_free(&my_hash_mul);
}

void mympq_add(mpq_t a, mpq_t b, mpq_t c) {
  mpq_t* q;

  if (mpq_sgn(b) == 0) {
    mpq_set(a, c);
    return;
  }

  if (mpq_sgn(c) == 0) {
    mpq_set(a, b);
    return;
  }

  if (mpz_size(mpq_numref(b))+mpz_size(mpq_denref(b))+
      mpz_size(mpq_numref(c))+mpz_size(mpq_denref(c)) < 32) {
    mpq_add(a, b, c);
    //fprintf(stderr, "a");
    return;
  }

  //putchar('1');
  q = hash_mpq_find(&my_hash_add, HASH_OPERATION_ADD, b, c);
  mpq_set(a, *q);
}

void mympq_sub(mpq_t a, mpq_t b, mpq_t c) {
  mpq_t* q;

  if (mpq_equal(b, c)) {
    mpq_set_si(a, 0, 1);
    return;
  }
  if (mpq_sgn(b) == 0) {
    mpq_neg(a, c);
    return;
  }
  if (mpq_sgn(c) == 0) {
    mpq_set(a, b);
    return;
  }

  if (mpz_size(mpq_numref(b))+mpz_size(mpq_denref(b))+
      mpz_size(mpq_numref(c))+mpz_size(mpq_denref(c)) < 32) {
    mpq_sub(a, b, c);
    //fprintf(stderr, "s");
    return;
  }

  //putchar('2');
  mpq_neg(mympq_tmp, c);
  q = hash_mpq_find(&my_hash_add, HASH_OPERATION_ADD, b, mympq_tmp);
  mpq_set(a, *q);
}

void mympq_mul(mpq_t a, mpq_t b, mpq_t c) {
  mpq_t* q;

  if (mpq_sgn(b) == 0 || mpq_sgn(c) == 0) {
    mpq_set_ui(a, 0, 1);
    return;
  }
  if (mpq_cmp_si(b,  1, 1) == 0) {mpq_set(a, c); return;}
  if (mpq_cmp_si(c,  1, 1) == 0) {mpq_set(a, b); return;}
  if (mpq_cmp_si(b, -1, 1) == 0) {mpq_neg(a, c); return;}
  if (mpq_cmp_si(c, -1, 1) == 0) {mpq_neg(a, b); return;}

  if (mpz_size(mpq_numref(b))+mpz_size(mpq_denref(b))+
      mpz_size(mpq_numref(c))+mpz_size(mpq_denref(c)) < 32) {
    mpq_mul(a, b, c);
    //fprintf(stderr, "m");
    return;
  }

  //putchar('3');
  q = hash_mpq_find(&my_hash_mul, HASH_OPERATION_MUL, b, c);
  mpq_set(a, *q);
}

void mympq_div(mpq_t a, mpq_t b, mpq_t c) {
  mpq_t* q;

  if (mpq_sgn(b) == 0) {
    mpq_set_ui(a, 0, 1);
    return;
  }
  if (mpq_cmp_si(b, 1, 1) == 0) {mpq_inv(a, c); return;}
  if (mpq_cmp_si(c, 1, 1) == 0) {mpq_set(a, b); return;}

  if (mpz_size(mpq_numref(b))+mpz_size(mpq_denref(b))+
      mpz_size(mpq_numref(c))+mpz_size(mpq_denref(c)) < 32) {
    mpq_div(a, b, c);
    //fprintf(stderr, "d");
    return;
  }

  //if (mpq_cmp_si(b, -1, 1)==0||mpq_cmp_si(c,-1,1)==0)putchar('o');else putchar('x');
  //putchar('4');
  mpq_inv(mympq_tmp, c);
  q = hash_mpq_find(&my_hash_mul, HASH_OPERATION_MUL, b, mympq_tmp);
  mpq_set(a, *q);
}

#endif
