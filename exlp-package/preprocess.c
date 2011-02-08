#include "preprocess.h"
#include "mylib.h"
#include "solve_lp.h"
#include <limits.h>


void scaling(LP* lp) {
  int    i, j;
  mpq_t  q1;
  mpq_t* q2;

  mpq_init(q1);

  for (i = 0; i < lp->rows; i ++) {
    mpq_set_si(q1, 0, 1);
    for (j = 0; j < lp->A->row[i]->nonzeros; j ++) {
      q2 = matrix_get_element_ptr(lp->A, i, lp->A->row[i]->i[j]);
      if (mympq_cmp_abs(q1, *q2) < 0)
        mpq_abs(q1, *q2);
    }
    if (mpq_sgn(q1) == 0)
      continue;
    mpq_inv(q1, q1);
    matrix_row_scalar_product(lp->A, q1, i);
    vector_mul_element(lp->b, q1, i);
  }

  for (i = 0; i < lp->vars; i ++) {
    mpq_set_si(q1, 0, 1);
    for (j = 0; j < lp->A->column[i]->nonzeros; j ++) {
      q2 = &(lp->A->column[i]->value[j]);
      if (mympq_cmp_abs(q1, *q2) < 0)
        mpq_abs(q1, *q2);
    }
    if (mpq_sgn(q1) == 0)
      continue;
    mpq_inv(q1, q1);
    matrix_column_scalar_product(lp->A, q1, i);
    vector_mul_element(lp->c, q1, i);
    vector_mul_element(lp->c_back, q1, i);
    if (lp->upper.is_valid[i])
      mympq_div(lp->upper.bound[i], lp->upper.bound[i], q1);
    if (lp->lower.is_valid[i])
      mympq_div(lp->lower.bound[i], lp->lower.bound[i], q1);
    vector_div_element(lp->x, q1, i);
    if (lp->is_basis[i]) {
      for (j = 0; j < lp->rows; j ++) {
        if (lp->basis_column[j] == i) {
          vector_div_element(lp->xb, q1, j);
          vector_mul_element(lp->cb, q1, j);
	}
      }
    }
  }

  mpq_clear(q1);
}

void remove_artificials_from_basis(LP* lp) {
  EXLPvector* e;
  EXLPvector* r;
  mpq_t  q;
  int  i, j;
  int  not_done;
  int  best, best_j;

  e = new_vector(lp->rows);
  r = new_vector(lp->rows);
  mpq_init(q);

  not_done = 1;
  while (not_done) {
    //printf("*");fflush(stdout);
    not_done = 0;
    reinversion(lp);

    for (i = 0; i < lp->rows; i ++) {
      if (!is_artificial_var(lp, lp->basis_column[i]))
        continue;
      vector_get_element(&q, lp->xb, i);
      if (mpq_sgn(q) != 0)
        continue;

      vector_zero_clear(e);
      vector_set_element(e, mympq_one, i);
      eta_file_btran(lp->eta, e, r);

      best   = 2; /* とりあえず密にはしたくないの */ //lp->rows + 1;
      best_j = -1;
      for (j = 0; j < lp->vars; j ++) {
        if (lp->is_basis[j] ||
            is_artificial_var(lp, j) || is_const_var(lp, j) ||
            (best_j >= 0 && best <= lp->A->column[j]->nonzeros))
          continue;
        vector_inner_product(&q, r, lp->A->column[j]);
        if (mpq_sgn(q) == 0)
	  // ほんとはこれを満たす制約って冗長なはず
          continue;
        best   = lp->A->column[j]->nonzeros;
        best_j = j;
      }

//printf("  %d",best);fflush(stdout);
      if (best_j >= 0 && best < 3/* とりあえず密にはしたくないの */) {
        vector_delete_element(lp->cb, i); /* これ, phase1 終了後のとき注意 */
        lp->is_basis[lp->basis_column[i]] = FALSE;
        lp->is_basis[best_j] = TRUE;
        lp->basis_column[i] = best_j;
        vector_get_element(&q, lp->x, best_j);
        vector_set_element(lp->xb, q, i);
        not_done = 1; /* ...わからん */
        break;
      }
    }
  }

  mpq_clear(q);
  vector_free(&e);
  vector_free(&r);
}

int non_const_vars_in_row(LP* lp, int row, int* var) {
  int  i, ret;

  ret = 0;
  for (i = 0; i < lp->A->row[row]->nonzeros; i ++) {
    if (is_const_var(lp, lp->A->row[row]->i[i]))
      continue;
    *var = lp->A->row[row]->i[i];
    ret ++;
  }

  return ret;
}

void check_trivial_constraints0(LP* lp) {
  /* 列に1度しか現れない自由変数. その行は意味無し.
     つーか, 列に1度しか現れない2つの変数x1, x2 >=0 が
       x1 - x2
     みたいになっててもその行意味無し.
     ついでに a1 x1 + ... + an xn < b みたいな制約がいくつかあって
     xk が下に有界じゃない, みたいなときにそれらの制約は無意味. */
  int  i, j, row;
  mpq_t  q1, q2;
  int* rr;
  int  ss;

  if (lp->print_sol)
    return;

  mpq_init(q1);
  mpq_init(q2);

  for (i = 0; i < lp->vars; i ++) {
    if (lp->A->column[i]->nonzeros != 1 || !is_free_var(lp, i))
      continue;
    row = lp->A->column[i]->i[0];
    if (lp->row_equality[row] != LP_EQUALITY_EQ)
      continue;
    vector_get_element(&q1, lp->c_back, i);
    if (mpq_sgn(q1)) {
      for (j = 0; j < lp->A->row[row]->nonzeros; j ++) {
	if (lp->A->row[row]->i[j] == i)
	  continue;
	mympq_div(q2,
                  *matrix_get_element_ptr(lp->A, row, lp->A->row[row]->i[j]),
                  lp->A->column[i]->value[0]);
	mympq_mul(q2, q2, q1);
	vector_sub_element(lp->c_back, q2, lp->A->row[row]->i[j]);
      }
      vector_get_element(&q2, lp->b, row);
      mympq_div(q2, q2, lp->A->column[i]->value[0]);
      mympq_mul(q2, q2, q1);
      //vector_add_element(lp->c_back, q2, INT_MAX);
      mympq_add(lp->c_const, lp->c_const, q2);
      vector_delete_element(lp->c_back, i);
    }
    lp_remove_row(lp, row);
    //fprintf(stderr, "t");
  }

  rr = (int*)malloc(sizeof(int)*lp->rows);
  for (i = 0; i < lp->rows; i ++)
    rr[i] = 0;
  for (i = 0; i < lp->vars; i ++) {
    if (lp->A->column[i]->nonzeros != 1 ||
        !vector_element_is_zero(lp->c_back, i) ||
        lp->upper.is_valid[i] ||
        !lp->lower.is_valid[i] ||
        mpq_sgn(lp->lower.bound[i]))  // もうちっと工夫できるけどめんどくさ
      continue;
    row = lp->A->column[i]->i[0];
    j = mpq_sgn(*matrix_get_element_ptr(lp->A, row, i));
    if (rr[row] == 0) {
      rr[row] = j;
    }
    if (rr[row] == -j) {
      rr[row] = lp->rows;
    }
  }
  for (i = lp->rows-1; i >= 0; i --) {
    if (rr[i] == lp->rows) {
      //fprintf(stderr, "&");
    lp_remove_row(lp, i);
    }
  }
  free(rr);

  //mpq_set_ui(q2, 1, 1);
  //vector_set_element(lp->x, q2, INT_MAX);
  vector_copy(lp->c, lp->c_back);

  for (i = 0; i < lp->vars; i ++) {
    vector_get_element(&q1, lp->c_back, i);
    if (mpq_sgn(q1))
      continue;
    ss = 0;
    for (j = 0; j < lp->A->column[i]->nonzeros; j ++) {
      row = lp->A->column[i]->i[j];
      if (lp->row_equality[row] == LP_EQUALITY_EQ)
        goto next;
      if ((lp->row_equality[row] == LP_EQUALITY_LE &&
	   mpq_sgn(lp->A->column[i]->value[j]) < 0) ||
	  (lp->row_equality[row] == LP_EQUALITY_GE &&
	   mpq_sgn(lp->A->column[i]->value[j]) > 0)) {
	if (ss < 0 || lp->upper.is_valid[i])
	  goto next;
	ss = 1;
      }
      else
      if ((lp->row_equality[row] == LP_EQUALITY_LE &&
	   mpq_sgn(lp->A->column[i]->value[j]) > 0) ||
	  (lp->row_equality[row] == LP_EQUALITY_GE &&
	   mpq_sgn(lp->A->column[i]->value[j]) < 0)) {
	if (ss > 0 || lp->lower.is_valid[i])
	  goto next;
	ss = -1;
      }
    }
    for (j = lp->A->column[i]->nonzeros-1; j >= 0; j --) {
      lp_remove_row(lp, lp->A->column[i]->i[j]);
      //fprintf(stderr, "t");
    }
  next:;
  }

  mpq_clear(q1);
  mpq_clear(q2);
}

int check_trivial_constraints1(LP* lp) {
  /* 3 x <= 5 のような制約式を変数の bound に.
     もし infeasible なことが分かれば LP_RESULT_INFEASIBLE を返す.
     そうでなければ 0 を返す. */
  mpq_t q;
  int  i, var;

  mpq_init(q);

  for (i = 0; i < lp->rows; i ++) {
    if (lp->A->row[i]->nonzeros != 1 ||
        lp->row_equality[i] != LP_EQUALITY_LE)
      continue;
    var = lp->A->row[i]->i[0];
    vector_get_element(&q, lp->b, i);
    mympq_div(q, q, *matrix_get_element_ptr(lp->A, i, var));

    if (mpq_sgn(*matrix_get_element_ptr(lp->A, i, var)) > 0) {
      if (lp->lower.is_valid[var] &&
          mpq_cmp(lp->lower.bound[var], q) > 0) {
        mpq_clear(q);
        return LP_RESULT_INFEASIBLE;
      }
      if (lp->upper.is_valid[var]) {
        if (mpq_cmp(lp->upper.bound[var], q) > 0)
        mpq_set(lp->upper.bound[var], q);
      } else {
        lp->upper.is_valid[var] = TRUE;
        mpq_init(lp->upper.bound[var]);
        mpq_set(lp->upper.bound[var], q);
      }
      lp_remove_row(lp, i);
      //fprintf(stderr, "T");
      i --;
    } else {
      if (lp->upper.is_valid[var] &&
          mpq_cmp(lp->upper.bound[var], q) < 0) {
        mpq_clear(q);
        return LP_RESULT_INFEASIBLE;
      }
      if (lp->lower.is_valid[var]) {
        if (mpq_cmp(lp->lower.bound[var], q) < 0)
        mpq_set(lp->lower.bound[var], q);
      } else {
        lp->lower.is_valid[var] = TRUE;
        mpq_init(lp->lower.bound[var]);
        mpq_set(lp->lower.bound[var], q);
      }
      lp_remove_row(lp, i);
      //fprintf(stderr, "T");
      i --;
    }
  }

  mpq_clear(q);
  return 0;
}

int check_trivial_constraints2(LP* lp) {
  /* a x + b y = 0 から y を削除.
     lp_add_slacks() より後に呼ぶこと. <-- うそ?
     もし infeasible なことが分かれば LP_RESULT_INFEASIBLE を返す.
     そうでなければ 0 を返す. */
  mpq_t q1, q2, q3;
  int  cont, i, j, k, l;

  if (lp->print_sol)  // いや, どーも. だってこれやったら解表示できないよね.
    return 0;

  mpq_init(q1);
  mpq_init(q2);
  mpq_init(q3);

  do {
    cont = 0;

    for (i = 0; i < lp->rows; i ++) {
      //vector_get_element(&q1, lp->b, i);if(mpq_sgn(q1) != 0)continue;
      if (lp->A->row[i]->nonzeros != 2 ||
          lp->row_equality[i] != LP_EQUALITY_EQ)
        continue;

      j = lp->A->row[i]->i[0];
      k = lp->A->row[i]->i[1];
      l = mpq_sgn(*matrix_get_element_ptr(lp->A, i, j)) *
          mpq_sgn(*matrix_get_element_ptr(lp->A, i, k));
      if (lp->lower.is_valid[k]) {
        vector_get_element(&q1, lp->b, i);
        mympq_mul(q2, *matrix_get_element_ptr(lp->A, i, k),
                      lp->lower.bound[k]);
        mympq_sub(q1, q1, q2);
        mympq_div(q1, q1, *matrix_get_element_ptr(lp->A, i, j));
        if (l < 0) {
          if (lp->lower.is_valid[j]) {
            if (mpq_cmp(lp->lower.bound[j], q1) < 0)
              mpq_set(lp->lower.bound[j], q1);
	  } else {
            mpq_init(lp->lower.bound[j]);
            mpq_set(lp->lower.bound[j], q1);
            lp->lower.is_valid[j] = 1;
	  }
	} else {
          if (lp->upper.is_valid[j]) {
            if (mpq_cmp(lp->upper.bound[j], q1) > 0)
              mpq_set(lp->upper.bound[j], q1);
	  } else {
            mpq_init(lp->upper.bound[j]);
            mpq_set(lp->upper.bound[j], q1);
            lp->upper.is_valid[j] = 1;
	  }
	}
      }
      if (lp->upper.is_valid[k]) {
        vector_get_element(&q1, lp->b, i);
        mympq_mul(q2, *matrix_get_element_ptr(lp->A, i, k),
                      lp->upper.bound[k]);
        mympq_sub(q1, q1, q2);
        mympq_div(q1, q1, *matrix_get_element_ptr(lp->A, i, j));
        if (l > 0) {
          if (lp->lower.is_valid[j]) {
            if (mpq_cmp(lp->lower.bound[j], q1) < 0)
              mpq_set(lp->lower.bound[j], q1);
	  } else {
            mpq_init(lp->lower.bound[j]);
            mpq_set(lp->lower.bound[j], q1);
            lp->lower.is_valid[j] = 1;
	  }
	} else {
          if (lp->upper.is_valid[j]) {
            if (mpq_cmp(lp->upper.bound[j], q1) > 0)
              mpq_set(lp->upper.bound[j], q1);
	  } else {
            mpq_init(lp->upper.bound[j]);
            mpq_set(lp->upper.bound[j], q1);
            lp->upper.is_valid[j] = 1;
	  }
	}
      }
      if (lp->lower.is_valid[j] && lp->upper.is_valid[j]) {
        if (mpq_cmp(lp->lower.bound[j], lp->upper.bound[j]) > 0) {
          return LP_RESULT_INFEASIBLE;
	}
      }

      vector_get_element(&q1, lp->b, i);
      mympq_div(q1, q1, *matrix_get_element_ptr(lp->A, i, k));
      mympq_div(q2, *matrix_get_element_ptr(lp->A, i, j),
                    *matrix_get_element_ptr(lp->A, i, k));

      vector_get_element(&q3, lp->c, k);
      mympq_mul(q3, q3, q1);
      mympq_add(lp->c_const, lp->c_const, q3);
      vector_get_element(&q3, lp->c, k);
      mympq_mul(q3, q3, q2);
      vector_sub_element(lp->c, q3, j);
      vector_sub_element(lp->c_back, q3, j);
      vector_delete_element(lp->c, k);
      vector_delete_element(lp->c_back, k);
      for (l = 0; l < lp->A->column[k]->nonzeros; ) {
        if (lp->A->column[k]->i[l] == i) {
          l ++;
          continue;
	}
        mympq_mul(q3, q1, lp->A->column[k]->value[l]);
        vector_sub_element(lp->b, q3, lp->A->column[k]->i[l]);
        mympq_mul(q3, q2, lp->A->column[k]->value[l]);
        matrix_sub_element(lp->A, q3, lp->A->column[k]->i[l], j);
        matrix_delete_element(lp->A, lp->A->column[k]->i[l], k);
      }
      lp_remove_row(lp, i);
      //fprintf(stderr, "u");
      //printf("%s %s\n",lp->var_name[j], lp->var_name[k]);
      i --;

      cont = 1;
    }
  } while (cont);

  mpq_clear(q1);
  mpq_clear(q2);
  mpq_clear(q3);
  return 0;
}

void check_trivial_constraints3(LP* lp) {
  /* a1 x1 + a2 x2 + ... + an xn <= b
     みたいな制約で, 左辺の最大値が b より小さければその行は無意味.
     当然だけどスラック変数をいれる前に呼ぶこと.
  */

  mpq_t q1, q2;
  int  cont, i, j, var;

  mpq_init(q1);
  mpq_init(q2);

  do {
    cont = 0;

    for (i = 0; i < lp->rows; i ++) {
      if (lp->row_equality[i] != LP_EQUALITY_LE)
        continue;

      mpq_set_si(q1, 0, 1);
      for (j = 0; j < lp->A->row[i]->nonzeros; j ++) {
        var = lp->A->row[i]->i[j];
        if (mpq_sgn(*matrix_get_element_ptr(lp->A, i, var)) < 0) {
          if (!lp->lower.is_valid[var])
            goto next;
          mympq_mul(q2, *matrix_get_element_ptr(lp->A, i, var),
                        lp->lower.bound[var]);
          mympq_add(q1, q1, q2);
	} else {
          if (!lp->upper.is_valid[var])
            goto next;
          mympq_mul(q2, *matrix_get_element_ptr(lp->A, i, var),
                        lp->upper.bound[var]);
          mympq_add(q1, q1, q2);
	}
      }
      vector_get_element(&q2, lp->b, i);
      if (mpq_cmp(q1, q2) <= 0) {
	//fprintf(stderr, "U");
        lp_remove_row(lp, i);
        i --;
        cont = 1;
      }
    next:;
    }
  } while (cont);

  mpq_clear(q1);
  mpq_clear(q2);
}

void check_tight_constraints1(LP* lp) {
  /* a1 x1 + a2 x2 + ... + an xn <= b
     みたいな制約で, 左辺の最小値が b ならば x1 から xn の値が決まる. */

  mpq_t q1, q2;
  int  cont, i, j, var;

  mpq_init(q1);
  mpq_init(q2);

  do {
    cont = 0;

    for (i = 0; i < lp->rows; i ++) {
      if (lp->row_equality[i] == LP_EQUALITY_GE)
        continue;

      mpq_set_si(q1, 0, 1);
      for (j = 0; j < lp->A->row[i]->nonzeros; j ++) {
        var = lp->A->row[i]->i[j];
        if (mpq_sgn(*matrix_get_element_ptr(lp->A, i, var)) > 0) {
          if (!lp->lower.is_valid[var])
            goto next;
          mympq_mul(q2, *matrix_get_element_ptr(lp->A, i, var),
                        lp->lower.bound[var]);
          mympq_add(q1, q1, q2);
	} else {
          if (!lp->upper.is_valid[var])
            goto next;
          mympq_mul(q2, *matrix_get_element_ptr(lp->A, i, var),
                        lp->upper.bound[var]);
          mympq_add(q1, q1, q2);
	}
      }

      vector_get_element(&q2, lp->b, i);
      if (mpq_equal(q1, q2)) {
        for (j = 0; j < lp->A->row[i]->nonzeros; j ++) {
          //printf(" %s\n", lp->var_name[lp->A->row[i]->i[j]]);
          var = lp->A->row[i]->i[j];
          if (mpq_sgn(*matrix_get_element_ptr(lp->A, i, var)) > 0) {
            if (lp->upper.is_valid[var]) {
              mpq_set(lp->upper.bound[var], lp->lower.bound[var]);
            } else {
              mpq_init(lp->upper.bound[var]);
              mpq_set(lp->upper.bound[var], lp->lower.bound[var]);
              lp->upper.is_valid[var] = TRUE;
	    }
	  } else {
            if (lp->lower.is_valid[var]) {
              mpq_set(lp->lower.bound[var], lp->upper.bound[var]);
            } else {
              mpq_init(lp->lower.bound[var]);
              mpq_set(lp->lower.bound[var], lp->upper.bound[var]);
              lp->lower.is_valid[var] = TRUE;
	    }
	  }
	}
        lp_remove_row(lp, i);
	//fprintf(stderr, "o");
        i --;
        cont = 1;
      }
    next:;
    }
  } while (cont);

  mpq_clear(q1);
  mpq_clear(q2);
}

void check_tight_constraints2(LP* lp) {
  /* a1 x1 + a2 x2 + ... + an xn >= b
     みたいな制約で, 左辺の最大値が b ならば x1 から xn の値が決まる. */

  mpq_t q1, q2;
  int  cont, i, j, var;

  mpq_init(q1);
  mpq_init(q2);

  do {
    cont = 0;

    for (i = 0; i < lp->rows; i ++) {
      if (lp->row_equality[i] == LP_EQUALITY_LE)
        continue;

      mpq_set_si(q1, 0, 1);
      for (j = 0; j < lp->A->row[i]->nonzeros; j ++) {
        var = lp->A->row[i]->i[j];
        if (mpq_sgn(*matrix_get_element_ptr(lp->A, i, var)) < 0) {
          if (!lp->lower.is_valid[var])
            goto next;
          mympq_mul(q2, *matrix_get_element_ptr(lp->A, i, var),
                        lp->lower.bound[var]);
          mympq_add(q1, q1, q2);
	} else {
          if (!lp->upper.is_valid[var])
            goto next;
          mympq_mul(q2, *matrix_get_element_ptr(lp->A, i, var),
                        lp->upper.bound[var]);
          mympq_add(q1, q1, q2);
	}
      }

      vector_get_element(&q2, lp->b, i);
      if (mpq_equal(q1, q2)) {
        for (j = 0; j < lp->A->row[i]->nonzeros; j ++) {
          //printf(" %s\n", lp->var_name[lp->A->row[i]->i[j]]);
          var = lp->A->row[i]->i[j];
          if (mpq_sgn(*matrix_get_element_ptr(lp->A, i, var)) < 0) {
            if (lp->upper.is_valid[var]) {
              mpq_set(lp->upper.bound[var], lp->lower.bound[var]);
            } else {
              mpq_init(lp->upper.bound[var]);
              mpq_set(lp->upper.bound[var], lp->lower.bound[var]);
              lp->upper.is_valid[var] = TRUE;
	    }
	  } else {
            if (lp->lower.is_valid[var]) {
              mpq_set(lp->lower.bound[var], lp->upper.bound[var]);
            } else {
              mpq_init(lp->lower.bound[var]);
              mpq_set(lp->lower.bound[var], lp->upper.bound[var]);
              lp->lower.is_valid[var] = TRUE;
	    }
	  }
	}
        lp_remove_row(lp, i);
	//fprintf(stderr, "O");
        i --;
        cont = 1;
      }
    next:;
    }
  } while (cont);

  mpq_clear(q1);
  mpq_clear(q2);
}

int check_const_vars(LP* lp) {
  /* 3x = 5 みたいに明らかに定数と分かる変数を定数扱いに.
     もし infeasible なことが分かれば LP_RESULT_INFEASIBLE を返す.
     そうでなければ 0 を返す. */
  mpq_t q1, q2;
  int  cont, i, j, var;

  mpq_init(q1);
  mpq_init(q2);

  for (i = 0; i < lp->vars; i ++) {
    if (!(lp->lower.is_valid[i]) || !(lp->upper.is_valid[i]) ||
        !mpq_equal(lp->lower.bound[i], lp->upper.bound[i]))
      continue;
    mpq_set(q1, lp->lower.bound[i]);
    vector_set_element(lp->x, q1, i);
    lp->var_type[i] = LP_VAR_TYPE_CONST;

    for (j = 0; j < lp->A->column[i]->nonzeros; j ++) {
      mympq_mul(q2, q1, lp->A->column[i]->value[j]);
      vector_sub_element(lp->b, q2, lp->A->column[i]->i[j]);
    }
    matrix_column_scalar_product(lp->A, mympq_zero, i);
  }

  do {
    cont = 0;

    for (i = 0; i < lp->rows; i ++) {
      if (lp->row_equality[i] != LP_EQUALITY_EQ)
        continue;
      j = non_const_vars_in_row(lp, i, &var);
      if (j == 0) {
        vector_get_element(&q1, lp->b, i);
	if (mpq_sgn(q1) != 0)
	  return LP_RESULT_INFEASIBLE;
	lp_remove_row(lp, i);
	//fprintf(stderr, "c");
	i --;
	continue;
      }
      if (j != 1)
	continue;

      matrix_get_element(&q1, lp->A, i, var);
      vector_get_element(&q2, lp->b, i);
      mympq_div(q2, q2, q1);

      if ((lp->upper.is_valid[var] &&
           mpq_cmp(q2, lp->upper.bound[var]) > 0) ||
          (lp->lower.is_valid[var] &&
           mpq_cmp(q2, lp->lower.bound[var]) < 0))
        return LP_RESULT_INFEASIBLE;

      lp_remove_row(lp, i);
      //fprintf(stderr, "c");
      i --;

      cont = 1;
      vector_set_element(lp->x, q2, var);
      lp->var_type[var] = LP_VAR_TYPE_CONST;

      //printf("%s %d  ",lp->var_name[var],mpq_sgn(q2));
      //print_rational_as_float(q2, 10);
      //printf("\n");fflush(stdout);

      for (j = 0; j < lp->A->column[var]->nonzeros; j ++) {
	mympq_mul(q1, q2, lp->A->column[var]->value[j]);
	vector_sub_element(lp->b, q1, lp->A->column[var]->i[j]);
      }
      matrix_column_scalar_product(lp->A, mympq_zero, var);
    }
  } while (cont && 0);

  mpq_clear(q1);
  mpq_clear(q2);
  return 0;
}

int are_parallel_constraints(LP* lp, int row1, int row2) {
  int  i, var;
  mpq_t  q1, q2;

  if (lp->A->row[row1]->nonzeros == 0 ||
      lp->A->row[row1]->nonzeros != lp->A->row[row2]->nonzeros)
    return 0;

  for (i = 0; i < lp->A->row[row1]->nonzeros; i ++) {
    if (vector_d_element_is_zero(lp->A->row[row2], lp->A->row[row1]->i[i]))
      return 0;
  }

  mpq_init(q1);
  mpq_init(q2);
  mympq_div(q1, *matrix_get_element_ptr(lp->A, row2, lp->A->row[row1]->i[0]),
                *matrix_get_element_ptr(lp->A, row1, lp->A->row[row1]->i[0]));

  for (i = 1; i < lp->A->row[row1]->nonzeros; i ++) {
    var = lp->A->row[row1]->i[i];
    mympq_mul(q2, q1, *matrix_get_element_ptr(lp->A, row1, var));
    if (!mpq_equal(q2, *matrix_get_element_ptr(lp->A, row2, var))) {
      mpq_clear(q1);
      mpq_clear(q2);
      return 0;
    }
  }

  mpq_clear(q1);
  mpq_clear(q2);
  return 1;
}

int check_parallel_constraints(LP* lp) {
  int  i, j, k;
  int  v0;
  mpq_t  q1, q2;

  mpq_init(q1);
  mpq_init(q2);

  for (i = 0; i < lp->rows; i ++) {
    if (lp->A->row[i]->nonzeros < 2)
      continue;
    v0 = lp->A->row[i]->i[0];
    for (k = 0; k < lp->A->column[v0]->nonzeros; k ++) {
      j = lp->A->column[v0]->i[k];
      if (j <= i)
        continue;
      if (!are_parallel_constraints(lp, i, j))
        continue;
      if (mpq_sgn(*matrix_get_element_ptr(lp->A, i, v0)) < 0) {
        matrix_rev_row_sgn(lp->A, i);
        vector_neg_element(lp->b, i);
        if (lp->row_equality[i] == LP_EQUALITY_LE)
          lp->row_equality[i] = LP_EQUALITY_GE;
        else if (lp->row_equality[i] == LP_EQUALITY_GE)
          lp->row_equality[i] = LP_EQUALITY_LE;
      }
      if (mpq_sgn(*matrix_get_element_ptr(lp->A, j, v0)) < 0) {
        matrix_rev_row_sgn(lp->A, j);
        vector_neg_element(lp->b, j);
        if (lp->row_equality[j] == LP_EQUALITY_LE)
          lp->row_equality[j] = LP_EQUALITY_GE;
        else if (lp->row_equality[j] == LP_EQUALITY_GE)
          lp->row_equality[j] = LP_EQUALITY_LE;
      }
      vector_get_element(&q1, lp->b, i);
      vector_get_element(&q2, lp->b, j);
      mympq_div(q2, q2, *matrix_get_element_ptr(lp->A, j, v0));
      mympq_mul(q2, q2, *matrix_get_element_ptr(lp->A, i, v0));

      if (lp->row_equality[i] == lp->row_equality[j]) {
        if (lp->row_equality[i] == LP_EQUALITY_EQ) {
          if (!mpq_equal(q1, q2))
            return LP_RESULT_INFEASIBLE;
          lp_remove_row(lp, j);
	  //fprintf(stderr, "p");
          break;
        } else if (lp->row_equality[i] == LP_EQUALITY_LE) {
          if (mpq_cmp(q1, q2) > 0)
            vector_set_element(lp->b, q2, i);
          lp_remove_row(lp, j);
	  //fprintf(stderr, "p");
          break;
        } else {
          if (mpq_cmp(q1, q2) < 0)
            vector_set_element(lp->b, q2, i);
          lp_remove_row(lp, j);
	  //fprintf(stderr, "p");
          break;
        }
      } else {
        if (lp->row_equality[i] == LP_EQUALITY_EQ) {
          if (lp->row_equality[j] == LP_EQUALITY_LE &&
              mpq_cmp(q1, q2) > 0)
            return LP_RESULT_INFEASIBLE;
          if (lp->row_equality[j] == LP_EQUALITY_GE &&
              mpq_cmp(q1, q2) < 0)
            return LP_RESULT_INFEASIBLE;
          lp_remove_row(lp, j);
	  //fprintf(stderr, "p");
          break;
        }
        if (lp->row_equality[j] == LP_EQUALITY_EQ) {
          if (lp->row_equality[i] == LP_EQUALITY_LE &&
              mpq_cmp(q1, q2) < 0)
            return LP_RESULT_INFEASIBLE;
          if (lp->row_equality[i] == LP_EQUALITY_GE &&
              mpq_cmp(q1, q2) > 0)
            return LP_RESULT_INFEASIBLE;
          vector_set_element(lp->b, q2, i);
          lp->row_equality[i] = LP_EQUALITY_EQ;
          lp_remove_row(lp, j);
	  //fprintf(stderr, "p");
          break;
        }
        if (lp->row_equality[i] == LP_EQUALITY_LE)
          mympq_sub(q2, q1, q2);
        else
          mympq_sub(q2, q2, q1);
        lp_set_row_range(lp, i, q2);
        lp_remove_row(lp, j);
	//fprintf(stderr, "p");
        break;
      }
    }
  }

  mpq_clear(q1);
  mpq_clear(q2);
  return 0;
}

int are_parallel_columns(LP* lp, int col1, int col2, int obj) {
/* 平行でないなら 0 を返す.
   平行で向きも同じなら 正の値を返す.
   平行で向きが逆なら負の値を返す.
   obj が非零の場合は目的関数も含めて平行かどうかをチェックする. */
  mpq_t  q1, q2;
  int  i, row;

  if (lp->A->column[col1]->nonzeros == 0 ||
      lp->A->column[col1]->nonzeros != lp->A->column[col2]->nonzeros)
    return 0;

  for (i = 0; i < lp->A->column[col1]->nonzeros; i ++) {
    if (vector_element_is_zero(lp->A->column[col2], lp->A->column[col1]->i[i]))
      return 0;
  }

  mpq_init(q1);
  mpq_init(q2);
  mympq_div(q1,
            *matrix_get_element_ptr(lp->A, lp->A->column[col1]->i[0], col2),
            *matrix_get_element_ptr(lp->A, lp->A->column[col1]->i[0], col1));

  if (obj) {
    mympq_mul(q2, q1, *vector_get_element_ptr(lp->c, col1));
    if (!mpq_equal(q2, *vector_get_element_ptr(lp->c, col2))) {
      mpq_clear(q1);
      mpq_clear(q2);
      return 0;
    }
  }

  for (i = 1; i < lp->A->column[col1]->nonzeros; i ++) {
    row = lp->A->column[col1]->i[i];
    mympq_mul(q2, q1, *matrix_get_element_ptr(lp->A, row, col1));
    if (!mpq_equal(q2, *matrix_get_element_ptr(lp->A, row, col2))) {
      mpq_clear(q1);
      mpq_clear(q2);
      return 0;
    }
  }

  i = mpq_sgn(q1);
  mpq_clear(q1);
  mpq_clear(q2);
  return i;
}

void check_parallel_columns(LP* lp) {
/* 目的関数まで含めて平行な2列はまとめられる. */
  int  i, j, s;
  int  k, l;

  for (i = 0; i < lp->vars-1; i ++) {
    if (is_const_var(lp, i) ||
        lp->A->column[i]->nonzeros == 0 ||
        lp->upper.is_valid[i] ||
        !lp->lower.is_valid[i] ||
        mpq_sgn(lp->lower.bound[i]))  // 考えるのが面倒なので...
      continue;
    k = lp->A->column[i]->i[0];
    for (l = lp->A->row[k]->nonzeros-1; l >= 0; l --) {
      j = lp->A->row[k]->i[l];
      if (i >= j)
        continue;

      if (is_const_var(lp, j) ||
	  lp->upper.is_valid[j] ||
	  !lp->lower.is_valid[j] ||
	  mpq_sgn(lp->lower.bound[j]))  // 考えるのが面倒なので...
	continue;
      s = are_parallel_columns(lp, i, j, 1);
      if (s > 0) {
	//fprintf(stderr, "<%d,%d %d>",i,j,s);
        lp->upper.is_valid[j] = TRUE;
        mpq_init(lp->upper.bound[j]);
        matrix_column_scalar_product(lp->A, mympq_zero, j);
      } else if (s < 0 && !lp->print_sol) {
	//fprintf(stderr, "{%d,%d %d}",i,j,s);
        lp->upper.is_valid[j] = TRUE;
        mpq_init(lp->upper.bound[j]);
        mpq_clear(lp->lower.bound[i]);
        lp->lower.is_valid[i] = FALSE;
        matrix_column_scalar_product(lp->A, mympq_zero, j);
      }
      s = are_parallel_columns(lp, i, j, 0);
      if (s > 0) {
	mpq_t  q;
	mpq_init(q);
	mpq_div(q, *matrix_get_element_ptr(lp->A, lp->A->column[i]->i[0], i),
		   *matrix_get_element_ptr(lp->A, lp->A->column[i]->i[0], j));
        mpq_mul(q, q, *vector_get_element_ptr(lp->c, j));
	if (mpq_cmp(*vector_get_element_ptr(lp->c, i), q) > 0) {
	  //fprintf(stderr, "[%d,%d %d]",i,j,s);
          lp->upper.is_valid[j] = TRUE;
	  mpq_init(lp->upper.bound[j]);
	  matrix_column_scalar_product(lp->A, mympq_zero, j);	  
	}
	mpq_clear(q);
      }
    }
  }
}

int are_almost_parallel_constraints(LP* lp, int row1, int row2) {
/* ほとんど平行かを返す. ほとんどっていうのは, 1変数を除いてってこと.
   ヒューリスティックなので, 全部見付けられるとは限らない.
   ほとんど平行だった場合, 平行でない column を返す.
   それ意外(ぴったり平行だった場合も)では負の値を返す.
   バグってるー. 直せるけどなんか意味無さげなんで放置.
*/
  int  i, j, k, var;
  mpq_t  q1, q2;

  if (lp->A->row[row1]->nonzeros < 2 ||
      (lp->A->row[row1]->nonzeros+1 < lp->A->row[row2]->nonzeros ||
       lp->A->row[row1]->nonzeros-1 > lp->A->row[row2]->nonzeros))
    return -1;

  j = 0;
  k = 0;
  for (i = 0; i < lp->A->row[row1]->nonzeros; i ++) {
    if (vector_d_element_is_zero(lp->A->row[row2], lp->A->row[row1]->i[i])) {
      if (j ++)
	return -1;
      k ++;
    }
  }

  mpq_init(q1);
  mpq_init(q2);
  mympq_div(q1, *matrix_get_element_ptr(lp->A, row2, lp->A->row[row1]->i[k]),
                *matrix_get_element_ptr(lp->A, row1, lp->A->row[row1]->i[k]));
  j = -1;
  for (i = 0; i < lp->A->row[row1]->nonzeros; i ++) {
    if (i == k)
      continue;
    var = lp->A->row[row1]->i[i];
    mympq_mul(q2, q1, *matrix_get_element_ptr(lp->A, row1, var));
    if (!mpq_equal(q2, *matrix_get_element_ptr(lp->A, row2, var))) {
      if (j >= 0) {
	mpq_clear(q1);
	mpq_clear(q2);
	return -1;
      }
      j = i;
    }
  }

  if (j < 0) {
    if (lp->A->row[row1]->nonzeros == lp->A->row[row2]->nonzeros) {
      k = -1; // 平行!
    } else {
      k = -1; // warning 対策
      for (i = 0; i < lp->A->row[row2]->nonzeros; i ++) {
	var = lp->A->row[row2]->i[i];
	if (matrix_element_is_zero(lp->A, row1, var)) {
	  k = var;
	  break;
	}
      }
    }
  } else {
    k = lp->A->row[row1]->i[j];
  }

  mpq_clear(q1);
  mpq_clear(q2);
  if (is_const_var(lp, k))
    return -1;
  return k;
}

int solve_almost_parallel_constraints(LP* lp, int r1, int r2, int c) {
  mpq_t  q1, q2, q3;
  int  i, j, ret;

  return 0;

  mpq_init(q1);
  mpq_init(q2);
  mpq_init(q3);

  for (i = 0; ; i ++) {
    j = lp->A->row[r1]->i[i];
    if (j != c)
      break;
  }

  mympq_div(q1, *matrix_get_element_ptr(lp->A, r1, j),
	        *matrix_get_element_ptr(lp->A, r2, j));

  mympq_mul(q2, *matrix_get_element_ptr(lp->A, r2, c), q1);
  mympq_sub(q2, q2, *matrix_get_element_ptr(lp->A, r1, c));
  mpq_neg(q2, q2);

  mympq_mul(q3, *vector_get_element_ptr(lp->b, r2), q1);
  mympq_sub(q3, q3, *vector_get_element_ptr(lp->b, r1));
  mpq_neg(q3, q3);

  mympq_div(q3, q3, q2);

  ret = 0;

  if (lp->upper.is_valid[c]) {
    if (mpq_cmp(lp->upper.bound[c], q3) < 0)
      fprintf(stderr,"!!!"),ret = LP_RESULT_INFEASIBLE;
    else
      mpq_set(lp->upper.bound[c], q3);
  } else {
    mpq_init(lp->upper.bound[c]);
    mpq_set(lp->upper.bound[c], q3);
    lp->upper.is_valid[c] = TRUE;
  }
  if (lp->lower.is_valid[c]) {
    if (mpq_cmp(lp->lower.bound[c], q3) > 0)
      ret = LP_RESULT_INFEASIBLE;
    else
      mpq_set(lp->lower.bound[c], q3);
  } else {
    mpq_init(lp->lower.bound[c]);
    mpq_set(lp->lower.bound[c], q3);
    lp->lower.is_valid[c] = TRUE;
  }


  mpq_clear(q1);
  mpq_clear(q2);
  mpq_clear(q3);

  return ret;
}

int check_almost_parallel_constraints(LP* lp) {
/* たとえば x+3y+z-w=1, x+3y-z-w=0 とかいう制約があれば,
   z=1/2 として 2本目の制約を除去できる.
   infeasible と分かれば LP_RESULT_INFEASIBLE を返す.
   そうでなければ 0 を返す.
*/
  int  i, r1, r2, c1, c2;

  for (r1 = 0; r1 < lp->rows; r1 ++) {
    if (lp->row_equality[r1] != LP_EQUALITY_EQ ||
	lp->A->row[r1]->nonzeros < 2)
      continue;
    c1 = lp->A->row[r1]->i[0];
    for (i = lp->A->column[c1]->nonzeros-1; i >= 0; i --) {
      r2 = lp->A->column[c1]->i[i];
      if (r1 >= r2)
	continue;
      if (lp->row_equality[r2] != LP_EQUALITY_EQ)
	continue;
      c2 = are_almost_parallel_constraints(lp, r1, r2);
      if (c2 < 0)
	continue;
      fprintf(stderr, "!%d %d,%d %d! ", c2,r1,r2,lp->A->row[r1]->nonzeros);
      if (solve_almost_parallel_constraints(lp, r1, r2, c2) ==
	  LP_RESULT_INFEASIBLE)
	return LP_RESULT_INFEASIBLE;
      //lp_remove_row(lp, r2);
    }
  }
  return 0;
}
