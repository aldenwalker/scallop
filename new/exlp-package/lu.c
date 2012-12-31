#include "lu.h"
#include <stdlib.h>
#include <gmp.h>
#include "matrix.h"
#include "lpstruct.h"

void lu_pivot(matrix* A, int i, EXLPvector* eta_vec) {
  mpq_t  r;
  int  j, k, row, col;

  mpq_init(r);
  vector_zero_clear(eta_vec);

  vector_set_element(eta_vec, *matrix_get_element_ptr(A, i, i), i);

  if (!mpq_equal(eta_vec->value[0], mympq_one)) {
    mpq_inv(eta_vec->value[0], eta_vec->value[0]);
    matrix_row_scalar_product(A, eta_vec->value[0], i);
  }

  for (j = 0; j < A->column[i]->nonzeros; ) {
    row = A->column[i]->i[j];
    if (row <= i) {
      j ++;
      continue;
    }

    matrix_get_element(&r, A, row, i);
    mympq_mul(r, r, *vector_get_element_ptr(eta_vec, i));
    mpq_neg(r, r);
    vector_set_element(eta_vec, r, row);

    for (k = 0; (col = A->row[i]->i[k]) != i; k ++) {
      mympq_mul(r, *matrix_get_element_ptr(A, i, col),
                   *matrix_get_element_ptr(A, row, i));
      matrix_sub_element(A, r, row, col);
    }
    for (k ++; k < A->row[i]->nonzeros; k ++) {
      col = A->row[i]->i[k];
      mympq_mul(r, *matrix_get_element_ptr(A, i, col),
                   *matrix_get_element_ptr(A, row, i));
      matrix_sub_element(A, r, row, col);
    }
    matrix_delete_element(A, row, i);
  }

  mpq_clear(r);
}

void lu_select_pivot_pos_minimum_degree(matrix* A, int i,
                                        int* best_row, int* best_col) {
  int  j, row, col;
  int  best;

  *best_row = -1;

  best = A->columns + 1;
  for (row = i; row < A->rows; row ++) {
    if (A->row[row]->nonzeros < best) {
      *best_row = row;
      best = A->row[row]->nonzeros;
    }
  }

  best = A->rows + 1;
  for (j = A->row[*best_row]->nonzeros-1; j >= 0; j --) {
    col = A->row[*best_row]->i[j];
    if (A->column[col]->nonzeros < best) {
      *best_col = col;
      best = A->column[col]->nonzeros;
    }
  }
}

void lu_select_pivot_pos_minimum_degree_b(matrix* A, int i,
                                        int* best_row, int* best_col) {
  int  j, k, l, n, row, col, row2, col2;
  int  best;

  *best_row = -1;

  best = A->rows * A->columns + 1;
  for (col = i; col < A->columns; col ++) {
    for (j = 0; j < A->column[col]->nonzeros; j ++) {
      row = A->column[col]->i[j];
      if (row < i)
        continue;
      n = 0;
      for (k = 0; k < A->column[col]->nonzeros; k ++) {
        row2 = A->column[col]->i[k];
        if (row2 == row || row2 < i)
          continue;
        for (l = 0; l < A->row[row]->nonzeros; l ++) {
          col2 = A->row[row]->i[l];
          if (col2 == col)
            continue;
          if (matrix_element_is_zero(A, row2, col2))
	    n ++;
	}
      }
      if (n < best ||
          (n == best && A->row[*best_row]->nonzeros > A->row[row]->nonzeros)) {
          //(n == best && A->column[*best_col]->nonzeros > A->column[col]->nonzeros)) {
        best = n;
        *best_row = row;
        *best_col = col;
      }
      if (best == 0)
        return;
    }
  }

  return;
}

void lp_LU_basis_minimum_degree(LP* lp) {
  /* 基底行列をLU分解. 結果は lp->eta に.
     lp->eta->U は上三角行列.
     lp->eta->P[i] は i 番目に i 行と swap する行.
     lp->eta->L[i] は i 番目にかけるエータ行列. */
  /* matrix_swap_rows() を呼んでるのねぇ...
     てなわけで leaving_row の位置が変わっちゃうのよ...
     なもんで l_row が決まってからは reinversion() とか
     呼んじゃダメと言う... */
  int  i, c;

  for (i = 0; i < lp->rows; i ++)
    matrix_set_column(lp->eta->U, lp->A->column[lp->basis_column[i]], i);

  for (i = 0; i < lp->rows; i ++) {

    lu_select_pivot_pos_minimum_degree(lp->eta->U, i, &(lp->eta->P[i]), &c);
    if (lp->eta->P[i] < 0) {
      fprintf(stderr, "rank of the matrix is too small %d\n", i);
      exit(EXIT_FAILURE);
    }

    matrix_swap_columns(lp->eta->U, i, c);
    swap_basis_columns(lp, i, c);

    matrix_swap_rows(lp->eta->U, i, lp->eta->P[i]);

    lu_pivot(lp->eta->U, i, lp->eta->L[i]->eta_vector);
    lp->eta->L[i]->eta_column = i;
  }
}

void lu_select_pivot_pos_markowitz(matrix* A, int i, int* c_table,
                                   int* best_row, int* best_col) {
  int  j, row, col;
  int  best, b;
  int  p, pr, pc;

  best = A->rows * A->columns;
  *best_row = -1;
  b = A->rows;

  for (col = i; col < A->columns; col ++) {
    pc = A->column[col]->nonzeros - c_table[col] - 1;
    for (j = A->column[col]->nonzeros-1; j >= 0; j --) {
      row = A->column[col]->i[j];
      if (row < i)
        continue;
      pr = A->row[row]->nonzeros - 1;
      if (pc == 0 && pr == 0) {
      //if (pc == 0 || pr == 0) {
      //とすると LU分解自体は速くなるんだけど、全体的に遅くなることが多い。
        *best_row = row;
        *best_col = col;
        return;
      }
      p = pc * pr;
      if (p > best)
        continue;
      if (p == best) {
        if (A->row[row]->nonzeros > b)
          continue;
        if (A->row[row]->nonzeros == b) {
          if (A->column[col]->nonzeros >= A->column[*best_col]->nonzeros)
            continue;
	}
      }
      best = p;
      //if (best == 1)putchar('o');else putchar('.');
      b = A->row[row]->nonzeros;
      *best_row = row;
      *best_col = col;
    }
  }
  //if (best == 1)putchar('%');else putchar('X');
  //if (best == 0)printf(" . ");else printf("{%d}", best);
}

void lu_select_pivot_pos_markowitz_b(matrix* A, int i, int* c_table,
                                   int* best_row, int* best_col) {
  int  j, k, row, col, col2;
  double  p, best;
  int  pr, pc, b;

  best = A->rows * A->columns;
  *best_row = -1;
  b = A->rows;

  for (col = i; col < A->columns; col ++) {
    pc = A->column[col]->nonzeros - c_table[col] - 1;
    for (j = A->column[col]->nonzeros-1; j >= 0; j --) {
      row = A->column[col]->i[j];
      if (row < i)
        continue;
      pr = A->row[row]->nonzeros - 1;
      if (pc == 0 && pr == 0) {
      //if (pc == 0 || pr == 0) {
      //とすると LU分解自体は速くなるんだけど、全体的に遅くなることが多い。
        *best_row = row;
        *best_col = col;
        return;
      }

      p = 0;
      for (k = 0; k < A->row[row]->nonzeros; k ++) {
        col2 = A->row[row]->i[k];
        if (col2 == col)
          continue;
        p += ((double)A->rows-i-1 - (A->column[col2]->nonzeros-c_table[col2]-1))/A->rows;
      }
      p *= (double)(A->column[col]->nonzeros-c_table[col]-1);


      if (p > best)
        continue;
      if (p == best) {
        if (A->row[row]->nonzeros > b)
          continue;
        if (A->row[row]->nonzeros == b) {
          if (A->column[col]->nonzeros >= A->column[*best_col]->nonzeros)
            continue;
	}
      }
      best = p;
      //if (best == 1)putchar('o');else putchar('.');
      b = A->row[row]->nonzeros;
      *best_row = row;
      *best_col = col;
    }
  }
  //if (best == 1)putchar('%');else putchar('X');
  //if (best == 0)printf(" . ");else printf("{%d}", best);
}

void lp_LU_basis_markowitz(LP* lp) {
  /* lp_LU_basis_minimum_degree() と同じ。
     ただし、Markowitz の規則でピボットする。
     より疎な LU分解が期待できるけど、結構時間かかる。*/
  int  i, j, c;
  int* c_table;

  c_table = my_calloc(lp->rows, sizeof(int));

  for (i = 0; i < lp->rows; i ++)
    matrix_set_column(lp->eta->U, lp->A->column[lp->basis_column[i]], i);

  for (i = 0; i < lp->rows; i ++) {

    lu_select_pivot_pos_markowitz(lp->eta->U, i, c_table,
                                  &(lp->eta->P[i]), &c);
    if (lp->eta->P[i] < 0) {
      fprintf(stderr, "rank of the matrix is too small %d\n", i);
      exit(EXIT_FAILURE);
    }

    matrix_swap_columns(lp->eta->U, i, c);
    swap_basis_columns(lp, i, c);

    matrix_swap_rows(lp->eta->U, i, lp->eta->P[i]);

    lu_pivot(lp->eta->U, i, lp->eta->L[i]->eta_vector);
    lp->eta->L[i]->eta_column = i;

    c_table[c] = c_table[i];
    for (j = 0; j < lp->eta->U->row[i]->nonzeros; j ++)
      c_table[lp->eta->U->row[i]->i[j]] ++;
  }

  /*
  j = 0;
  for (i = 0; i < lp->rows; i ++) {
    if (lp->eta->U->column[i]->nonzeros != 1 ||
        lp->eta->U->row[i]->nonzeros != 1 ||
        lp->eta->L[i]->eta_vector->nonzeros != 1)
      j = 1;
    else if (lp->eta->U->column[i]->nonzeros == 1 &&
	     lp->eta->U->row[i]->nonzeros == 1 &&
	     lp->eta->L[i]->eta_vector->nonzeros == 1 && j)
      printf("!!!!!!\n");
  }
  */

  free(c_table);
}

int how_dense_b(LP* lp) {
  int  i, j, n;

  n = 0;
  for (j = 0; j < lp->rows; j ++) {
    i = lp->basis_column[j];
    n += lp->A->column[i]->nonzeros;
  }
  return n;
}

int how_dense_a(LP* lp) {
  int  i, n;

  n = 0;
  for (i = 0; i < lp->eta->U->columns; i ++) {
    n += lp->eta->L[i]->eta_vector->nonzeros;
    n += lp->eta->U->column[i]->nonzeros;
    n --;
  }
  return n;
}

void show_how_dense_b(LP* lp) {
  int  i, j, k;

  for (i = 0; i < lp->rows; i ++) {
    for (k = 0; k < lp->rows; k ++) {
      j = lp->basis_column[k];
      if (matrix_element_is_zero(lp->A, i, j))
        putchar('.');
      else if (mpq_equal(*matrix_get_element_ptr(lp->A, i, j), mympq_one) ||
	       mpq_equal(*matrix_get_element_ptr(lp->A, i, j), mympq_minus_one))
        putchar('1');
      else
        putchar('*');
    }
    putchar('\n');
  }
}

void show_how_dense_a(LP* lp) {
  int  i, j;

  for (i = 0; i < lp->rows; i ++) {
    for (j = 0; j < i; j ++) {
      if (vector_element_is_zero(lp->eta->L[j]->eta_vector, i))
        putchar('.');
      else if (mpq_equal(*vector_get_element_ptr(lp->eta->L[j]->eta_vector, i), mympq_one) ||
	       mpq_equal(*vector_get_element_ptr(lp->eta->L[j]->eta_vector, i), mympq_minus_one))
        putchar('1');
      else
        putchar('*');
    }
    for (; j < lp->rows; j ++) {
      if (matrix_element_is_zero(lp->eta->U, i, j))
        putchar('.');
      else if (mpq_equal(*matrix_get_element_ptr(lp->eta->U, i, j), mympq_one) ||
	       mpq_equal(*matrix_get_element_ptr(lp->eta->U, i, j), mympq_minus_one))
        putchar('1');
      else
        putchar('*');
    }
    putchar('\n');
  }
}

void lp_LU_basis(LP* lp) {
  //printf("before\t(%d)\n", how_dense_b(lp));fflush(stdout);//show_how_dense_b(lp);putchar('\n');
  if (lp->lu_rule == LU_MARKOWITZ)
    lp_LU_basis_markowitz(lp);
  else
    lp_LU_basis_minimum_degree(lp);
  //printf("after\t(%d)\n", how_dense_a(lp));fflush(stdout);//show_how_dense_a(lp);putchar('\n');
  //putchar('\n');fflush(stdout);
}

void lp_select_dual_feasible_basis(LP* lp) {
/* dual の phase 1 の前に基底をどれにするか選ぶときに使う.
   とりあえず lp->A を LU 分解する. で, 結果は lp->L と lp->U に.
   と思ったけど lp->U だけでいいや. で始めの反復の前に reinversion()
   してくださいと.
   lp_check_redundancy() と被りすぎな感じ...
   lp->is_basis[] とか lp->basis_column[] の面倒もみてやる.
   で, ランク落ちしてたら lp_remove_row() してやると.
   ところで, ランク落ち部分のせいで infeasible になるときの対処...
*/
  matrix* A;
  EXLPvector* v;
  int  i, j, r, c, t;
  int* rr;
  int* cc;
  double* dummy;
  int* c_table;

  c_table = my_calloc(lp->A->columns, sizeof(int));

  A = new_matrix(lp->A->rows, lp->A->columns);
  for (i = 0; i < lp->vars; i ++) {
    if (is_const_var(lp, i)) {
      matrix_column_scalar_product(A, mympq_zero, i);
      continue;
    }
    matrix_set_column(A, lp->A->column[i], i);
  }

  rr = my_malloc(MAX(lp->rows,lp->vars)*sizeof(int));
  cc = my_malloc(MAX(lp->rows,lp->vars)*sizeof(int));
  for (i = 0; i < lp->rows; i ++)
    rr[i] = i;
  for (i = 0; i < lp->vars; i ++)
    cc[i] = i;

  v = new_vector(lp->rows);

  for (i = 0; i < A->rows; i ++) {

    lu_select_pivot_pos_markowitz(A, i, c_table, &r, &c);
    if (r < 0) {
      dummy = my_malloc(A->rows*sizeof(double));
      my_sort_d(rr, dummy, i, A->rows-1);
      //fprintf(stderr,"here\n");
      for (c = A->rows-1; c >= i; c --)
        lp_remove_row(lp, rr[c]);
      free(dummy);
      break;
    }

    matrix_swap_columns(A, i, c);
    t = cc[i];
    cc[i] = cc[c];
    cc[c] = t;
    matrix_swap_rows(A, i, r);
    t = rr[i];
    rr[i] = rr[r];
    rr[r] = t;

    lu_pivot(A, i, v);

    c_table[c] = c_table[i];
    for (j = 0; j < A->row[i]->nonzeros; j ++)
      c_table[A->row[i]->i[j]] ++;
  }

  for (i = 0; i < lp->vars; i ++)
    lp->is_basis[i] = FALSE;

  for (i = 0; i < lp->rows; i ++) {
    lp->basis_column[i] = cc[i];
    lp->is_basis[cc[i]] = TRUE;
    vector_set_element(lp->cb, *vector_get_element_ptr(lp->c, cc[i]), i);
  }

  matrix_free(A);
  vector_free(&v);
  free(rr);
  free(cc);
  free(c_table);
}

int lp_select_basis(LP* lp) {
/* primal の phase 1 の前に基底をどれにするか選ぶときに使う.
   全部選ぶことはできないので, 
   まず列に一つしか要素がないもので, feasible にできるやつを選び,
   次に, 右辺が 0 の行を集めて LU分解する. ま, そんな感じで.
   なお, 実際には基底が選ばれなかった行については lp->basis_column[.] == -1
   とする.
   返す値は, 実際に基底として選んだ変数の数.
*/
  matrix* A;
  EXLPvector* v;
  mpq_t  q1, q2;
  int  i, j, x;
  int  r, c, t;
  int* rr;
  int* cc;
  int* c_table;

  mpq_init(q1);
  mpq_init(q2);

  x = 0;

  for (i = 0; i < lp->rows; i ++)
    lp->basis_column[i] = -1;

  for (i = 0; i < lp->A->columns; i ++) {
    if (lp->A->column[i]->nonzeros != 1 ||
        is_const_var(lp, i) ||
        lp->basis_column[lp->A->column[i]->i[0]] != -1)
      continue;
    vector_get_element(&q2, lp->x, i);
    mympq_mul(q2, q2, lp->A->column[i]->value[0]);
    j = lp->A->column[i]->i[0];
    get_valid_rhs(&q1, lp, j);
    mympq_add(q1, q1, q2);
    mympq_div(q1, q1, lp->A->column[i]->value[0]);
    if (lp->upper.is_valid[i] &&
        mpq_cmp(q1, lp->upper.bound[i]) > 0)
      continue;
    if (lp->lower.is_valid[i] &&
        mpq_cmp(q1, lp->lower.bound[i]) < 0)
      continue;
    x ++;
    lp->is_basis[i] = TRUE;
    lp->basis_column[j] = i;
    vector_set_element(lp->x,  q1, i);
    vector_set_element(lp->xb, q1, j);
  }
  //return x;

  A = new_matrix(lp->rows, lp->A->columns);
  matrix_copy(A, lp->A);
  for (i = 0; i < lp->rows; i ++) {
    get_valid_rhs(&q1, lp, i);
    if (mpq_sgn(q1) != 0 || lp->basis_column[i] >= 0)
      matrix_row_scalar_product(A, mympq_zero, i);
  }
  for (i = 0; i < lp->vars; i ++) {
    if (lp->is_basis[i] || is_const_var(lp, i))
      matrix_column_scalar_product(A, mympq_zero, i);
  }

  c_table = my_calloc(lp->A->columns, sizeof(int));
  rr = my_malloc(MAX(lp->vars, lp->rows)*sizeof(int));
  cc = my_malloc(MAX(lp->vars, lp->rows)*sizeof(int));
  for (i = 0; i < lp->rows; i ++)
    rr[i] = i;
  for (i = 0; i < lp->vars; i ++)
    cc[i] = i;

  //aaaaa
  v = new_vector(lp->rows);

  for (i = 0; i < A->rows; i ++) {

    lu_select_pivot_pos_markowitz(A, i, c_table, &r, &c);
    if (r < 0)
      break;

    matrix_swap_columns(A, i, c);
    t = cc[i];
    cc[i] = cc[c];
    cc[c] = t;
    matrix_swap_rows(A, i, r);
    t = rr[i];
    rr[i] = rr[r];
    rr[r] = t;

    lu_pivot(A, i, v);

    c_table[c] = c_table[i];
    for (j = 0; j < A->row[i]->nonzeros; j ++)
      c_table[A->row[i]->i[j]] ++;
  }

  for (j = 0; j < i; j ++) {
    if ((lp->lower.is_valid[cc[j]] && mpq_sgn(lp->lower.bound[cc[j]]) > 0) ||
        (lp->upper.is_valid[cc[j]] && mpq_sgn(lp->upper.bound[cc[j]]) < 0))
      continue;
    lp->basis_column[rr[j]] = cc[j];
    lp->is_basis[cc[j]] = TRUE;
    x ++;
    vector_set_element(lp->cb, *vector_get_element_ptr(lp->c, cc[j]), rr[j]);
    vector_delete_element(lp->x, cc[j]);
    vector_delete_element(lp->xb, rr[j]);
  }

  if (lp->preprocess) {
    for (j = 0; j < lp->rows; j ++)
      cc[j] = 0;
    for (j = i; j < lp->A->rows; j ++) {
      get_valid_rhs(&q1, lp, rr[j]);
      if (mpq_sgn(q1) == 0 &&
	  lp->basis_column[rr[j]] < 0)
	//fprintf(stderr,"cccccc %d\n", rr[j]);
	//{lp_remove_row(lp, rr[j]);break;}
	cc[rr[j]] = 1;
    }
    for (j = lp->rows-1; j >= 0; j --)
      if (cc[j])
	lp_remove_row(lp, j);//,fputc('%',stderr);
  }

  matrix_free(A);
  vector_free(&v);
  free(rr);
  free(cc);
  free(c_table);

  return x;
}

int lp_check_redundancy(LP* lp) {
  /* preprocess の一番始めで呼ぶこと. */
  matrix* A;
  EXLPvector* b;
  eta_matrix  E;
  int  i, j, k;
  int  r, c, t;
  int* rr;
  int* ss;
  int* rs;
  int* c_table;

  for (i = lp->rows-1; i >= 0; i --) {
    if (lp->A->row[i]->nonzeros == 0) {
      if (vector_element_is_zero(lp->b, i)) {
	lp_remove_row(lp, i);
	continue;
      }
      if (lp->row_equality[i] == LP_EQUALITY_EQ)
	return LP_RESULT_INFEASIBLE;
      if (lp->row_equality[i] == LP_EQUALITY_LE &&
	  mpq_sgn(*vector_get_element_ptr(lp->b, i)) < 0)
	return LP_RESULT_INFEASIBLE;
      if (lp->row_equality[i] == LP_EQUALITY_GE &&
	  mpq_sgn(*vector_get_element_ptr(lp->b, i)) > 0)
	return LP_RESULT_INFEASIBLE;
      lp_remove_row(lp, i);
      //putchar('X'), fflush(stdout);
    }
  }

  A = new_matrix(lp->rows, lp->A->columns);
  matrix_copy(A, lp->A);
  for (i = 0; i < lp->rows; i ++) {
    if (lp->row_equality[i] != LP_EQUALITY_EQ)
      matrix_row_scalar_product(A, mympq_zero, i);
  }
  b = new_vector(lp->rows);
  vector_copy(b, lp->b);
  E.dimension = lp->rows;
  E.eta_vector = new_vector(lp->rows);

  c_table = my_calloc(lp->A->columns, sizeof(int));
  rr = my_malloc(lp->rows*sizeof(int));
  rs = my_malloc(lp->rows*sizeof(int));
  ss = my_malloc(lp->vars*sizeof(int));
  for (i = 0; i < lp->rows; i ++)
    rr[i] = i;
  for (i = 0; i < lp->vars; i ++)
    ss[i] = i;

  for (i = 0; i < A->rows; i ++) {

    lu_select_pivot_pos_markowitz(A, i, c_table, &r, &c);
    if (r < 0)
      break;

    matrix_swap_columns(A, i, c);
    matrix_swap_rows(A, i, r);
    vector_swap_elements(b, i, r);
    t = rr[i];  rr[i] = rr[r];  rr[r] = t;
    t = ss[i];  ss[i] = ss[c];  ss[c] = t;

    lu_pivot(A, i, E.eta_vector);
    E.eta_column = i;
    eta_matrix_mul_vector(&E, b);

    c_table[c] = c_table[i];
    for (j = 0; j < A->row[i]->nonzeros; j ++)
      c_table[A->row[i]->i[j]] ++;
  }

  for (j = 0; j < i; j ++) {
    if (A->row[j]->nonzeros < lp->A->row[rr[j]]->nonzeros) {
      //continue;
      //putchar('q'),fflush(stdout);
      matrix_row_scalar_product(lp->A, mympq_zero, rr[j]);
      for (k = 0; k < A->row[j]->nonzeros; k ++)
        matrix_set_element(lp->A,
                           *matrix_get_element_ptr(A, j, A->row[j]->i[k]),
                           rr[j], ss[A->row[j]->i[k]]);
      vector_set_element(lp->b, *vector_get_element_ptr(b, j), rr[j]);
    }
  }

  for (j = 0; j < lp->rows; j ++)
    rs[j] = 0;

  for (j = i; j < lp->A->rows; j ++) {
    if (lp->row_equality[rr[j]] == LP_EQUALITY_EQ) {
      if (!vector_element_is_zero(b, j))
        return LP_RESULT_INFEASIBLE;
      rs[rr[j]] = 1;
    }
  }
  for (j = lp->rows-1; j >= 0; j --)
    if (rs[j])
      lp_remove_row(lp, j);//,putchar('R'),fflush(stdout);

  matrix_free(A);
  vector_free(&b);
  vector_free(&(E.eta_vector));
  free(rr);
  free(ss);
  free(rs);
  free(c_table);

  return 0;
}
