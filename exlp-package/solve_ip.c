#include "mylib.h"
#include "solve_lp.h"
#include "preprocess.h"

void make_all_elements_integer(LP* lp) {
  mpz_t  z;
  mpq_t  q;
  int  i, j, k;

  mpz_init(z);
  mpq_init(q);

  /*
  mpz_set_si(z, 1);
  for (i = 0; i < lp->c->nonzeros; i ++) {
    mpz_lcm(z, z, mpq_denref(lp->c->value[i]));
  }
  mpq_set_z(q, z);
  for (i = 0; i < lp->c->nonzeros; i ++) {
    mympq_mul(lp->c->value[i], lp->c->value[i], q);
  }
  */

  for (i = 0; i < lp->rows; i ++) {
    if (vector_element_is_zero(lp->b, i))
      mpz_set_si(z, 1);
    else
      mpz_set(z, mpq_denref(*vector_get_element_ptr(lp->b, i)));
    for (k = 0; k < lp->A->row[i]->nonzeros; k ++) {
      j = lp->A->row[i]->i[k];
      mpz_lcm(z, z, mpq_denref(*matrix_get_element_ptr(lp->A, i, j)));
    }
    mpq_set_z(q, z);
    vector_mul_element(lp->b, q, i);
    for (k = 0; k < lp->A->row[i]->nonzeros; k ++) {
      j = lp->A->row[i]->i[k];
      matrix_mul_element(lp->A, q, i, j);
    }
  }

  /* ������ɤ���? */
  for (i = 0; i < lp->vars; i ++) {
    if (lp->lower.is_valid[i])
      mympq_ceil(lp->lower.bound[i], lp->lower.bound[i]);
    if (lp->upper.is_valid[i])
      mympq_floor(lp->upper.bound[i], lp->upper.bound[i]);
  }

  mpz_clear(z);
  mpq_clear(q);
}

int is_integer_optimal(LP* lp) {
  /* ���ä�, ��ͳ�ѿ��ä������Ǥ������Ȥϸ¤�ʤ��櫓��.
     �ɤ�����Τ���?
  */
  int  i;

  for (i = 0; i < lp->xb->nonzeros; i ++) {
    if (!lp->is_integer[lp->basis_column[lp->xb->i[i]]])
      continue;
    if (!mympq_is_integer(lp->xb->value[i]))
      return FALSE;
  }

  return TRUE;
}

int gomory_select_source_row(LP* lp, int* c) {
  /*
     source row �ˤǤ���Ԥ��֤�.
     ����Ȳ�����.
     lp->rows �ʤ���Ū�ؿ�. ������ͥ��ˤ�������®����.
     c �ˤ������ˤʤ�ʤ��Ǿ���ź�����֤�. ��������Ū�ؿ��ΤȤ��Ȥ�������.
   */
  int  i, j;
  EXLPvector* v;


  /*
  for (i = 0; i < lp->c->nonzeros; i ++) {
    if (lp->is_basis[i])
      continue;
    if (lp->var_name[lp->basis_column[i]][0] == '#' &&
	lp->var_name[lp->basis_column[i]][1] == 'I')
      continue;
    if (!mympq_is_integer(*vector_get_element_ptr(lp->c, i))) {
      return lp->rows;
    }
  }
  */
  // ���Ȥ��ǤϤ���?
  for (i = 0; i < lp->rows; i ++) {
    if (!mympq_is_integer(*vector_get_element_ptr(lp->xb, i))) {


      EXLPvector* v1;
      EXLPvector* v2;
      mpq_t  q1, q2;
      int s;
      s = i;
      /*
      // �����򥳥��ȥ����Ȥ��Ƥ����� k6.mps ���ߤޤ�.
      // ������...
      if (lp->var_name[lp->basis_column[s]][0] == '#' &&
          lp->var_name[lp->basis_column[s]][1] == 'I')
        continue;
      */
      v1 = new_vector(lp->rows);
      v2 = new_vector(lp->rows);
      mpq_init(q1);
      mpq_init(q2);
      vector_set_element(v2, mympq_one, s);
      eta_file_btran(lp->eta, v2, v1);
      for (i = 0; i < lp->vars-1; i ++) {
	if (lp->is_basis[i])
	  continue;
	vector_inner_product(&q1, v1, lp->A->column[i]);
	mympq_floor(q2, q1);
        if (!mpq_equal(q1, q2)) {
	  mpq_clear(q1);
	  mpq_clear(q2);
	  vector_free(&v1);
	  vector_free(&v2);
          *c = i;
          return s;
	}
      }
      mpq_clear(q1);
      mpq_clear(q2);
      vector_free(&v1);
      vector_free(&v2);
    }
  }
  printf("!!!!!!!!!!!!!!!!!\n");



  v = new_vector(lp->rows);

  /*
  for (i = 0; i < lp->c->nonzeros; i ++) {
    if (!mympq_is_integer(lp->c)) {
      vector_free(&v);
      return lp->rows;
    }
  }
  */

  for (i = 0; i < lp->vars; i ++) {
    if (lp->is_basis[i])
      continue;
    eta_file_ftran(lp->eta, lp->A->column[i], v, NULL);
    for (j = 0; j < v->nonzeros; j ++) {
      // ��������ˤߤĤ��ä��餹��������֤��ΤϤ��äȤΤ�.
      // ����ˤ����ǲ򤤤��������β�̵�̤ˤʤ�Τ��ɤ��ʤ��ΤǤ�?
      // �Ȼפ�����, ̵�̤ˤʤ�ʤ��Τ�. �Ĥޤ�, �����ޤǤ��Ƥ�֤�
      // ���������Ǥ���������. ����Ȥ��Τ�������ɤ��Τ���?
      if (lp->var_name[lp->basis_column[v->i[j]]][0] == '#' &&
          lp->var_name[lp->basis_column[v->i[j]]][1] == 'I')
        continue;
      if (!mympq_is_integer(*vector_get_element_ptr(v, v->i[j]))) {
        *c = i;
        i = v->i[j];
	vector_free(&v);
	return i;
      }
    }
  }

  // �����ˤϤ��ʤ��Ϥ�
  vector_free(&v);
  return -1;
}

void add_gomory_cut(LP* lp) {
  EXLPvector* v1;
  EXLPvector* v2;
  mpq_t  q1;
  mpq_t  q2;
  static int  cuts = 0;
  char  name[LP_NAME_LEN_MAX];
  int  i, s;

  v1 = new_vector(lp->rows);
  v2 = new_vector(lp->rows);
  mpq_init(q1);
  mpq_init(q2);

  s = gomory_select_source_row(lp, &i);

  sprintf(name, "#C%d", cuts);
  lp_add_row(lp, name);
  vector_resize(lp->xb, lp->rows);
  vector_resize(lp->cb, lp->rows);

  sprintf(name, "#I%d", cuts++);
  lp_add_var(lp, name);
  vector_d_resize(lp->c_d, lp->vars);//���ޤ�Τ�...

  if (s == lp->rows-1) {
    for (i = 0; i < lp->vars-1; i ++) {
      if (lp->is_basis[i])
        continue;
      //vector_inner_product(&q1, v1, lp->A->column[i]);
      mpq_set(q1, *vector_get_element_ptr(lp->c, i));
      mympq_floor(q2, q1);
      mympq_sub(q1, q2, q1);
      lp_set_coefficient(lp, q1, lp->rows-1, i);
    }
  } else {
    /*
    // Avis �˶�������Τ����. �Τ��ä�.
    for ( ; i < lp->vars-1; i ++) {
      eta_file_ftran(lp->eta, lp->A->column[i], v1, NULL);
      vector_get_element(&q1, v1, s);
      mympq_floor(q1, q1);
      mympq_sub(q1, q1, *vector_get_element_ptr(v1, s));
      lp_set_coefficient(lp, q1, lp->rows-1, i);
    }
    */

    vector_set_element(v2, mympq_one, s);
    eta_file_btran(lp->eta, v2, v1);
    for ( ; i < lp->vars-1; i ++) {
      if (lp->is_basis[i])
        continue;
      vector_inner_product(&q1, v1, lp->A->column[i]);
      mympq_floor(q2, q1);
      mympq_sub(q1, q2, q1);
      lp_set_coefficient(lp, q1, lp->rows-1, i);
    }

  }
  lp_set_coefficient(lp, mympq_one, lp->rows-1, lp->vars-1);

  vector_get_element(&q1, lp->xb, s);
  mympq_floor(q1, q1);
  mympq_sub(q1, q1, *vector_get_element_ptr(lp->xb, s));
  lp_set_rhs(lp, lp->rows-1, q1);
  vector_set_element(lp->xb, q1, lp->rows-1);
  vector_set_element(lp->x,  q1, lp->vars-1);

  lp->basis_column[lp->rows-1] = lp->vars-1;
  lp->is_basis[lp->vars-1] = TRUE;
  lp->is_integer[lp->vars-1] = FALSE;

  /* �߷פ������Ȥ�����, ��, etafile �˴ؤ��Ƥ� lp_add_row() ������
     �ߤƤ���ʤ���... ��, ���ϼ�ȴ�������ɷ�ɤ��줬�����ɤ��Τ���... */
  eta_file_free(lp->eta);
  lp->eta = new_eta_file(lp->rows);

  //reinversion(lp); �� solve_lp_core_dual() ��Ƭ�ǸƤ�Ǥ����Τ�.

  vector_free(&v1);
  vector_free(&v2);
  mpq_clear(q1);
  mpq_clear(q2);
}

int solve_ip(LP* lp) {
  int  result;
  int  i;
  mpq_t  q;

  mpq_init(q);

  /* Gomory cut ����ľ�ˤ���褦�Ȥ���� �ѿ��β��¤ä� 0 ����ʤ��Ȥ��� */
  /* �ޤ�, ���Фʤ�Ȥ��Ǥ��ʤ��Ϥʤ��������, �Ȥꤢ����              */
  for (i = 0; i < lp->vars; i ++) {
    if (is_const_var(lp, i))
      continue;
    lp->is_integer[i] = TRUE; // -_-;
    if (!lp->lower.is_valid[i] || mpq_sgn(lp->lower.bound[i]) != 0) {
      fprintf(stderr, "not supported type of problem. sorry.\n");
      exit(EXIT_FAILURE);
    }
  }

  make_all_elements_integer(lp);

  lp->c_d = new_vector_d(lp->c->dimension);
  vector_get_d(lp->c, lp->c_d);

  lp->preprocess = FALSE; // �Ȥꤢ������
  lp->scaling = FALSE;    // ����
  result = solve_lp_dual(lp);

  while (result == LP_RESULT_OPTIMAL && !is_integer_optimal(lp)) {
    printf("------ %d %d ------\n", result, lp->rows);
    lp_get_object_value(lp, &q);
    //print_rational_as_float(q, 16);putchar('\n');
    print_rational_as_float(q, 64);putchar('\n');
    //lp_print_nonzero_vars(lp);
    //vector_print(lp->xb);
    add_gomory_cut(lp);
    lp->phase = 2;
    //vector_get_d(lp->c, lp->c_d);//dual �Ǥ� c_d �ȤäƤʤ����餤��ʤ�.
    result = solve_lp_core_dual(lp);
    //lp_print_nonzero_vars(lp);
  }

  mpq_clear(q);
  //return result;
  return LP_RESULT_OPTIMAL;
}
