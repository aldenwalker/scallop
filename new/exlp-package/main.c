#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#include <stdio.h>
#include <string.h>
#include "EXLPvector.h"
#include "matrix.h"
#include "eta_file.h"
#include "mylib.h"
#include "lpstruct.h"
#include "mps.h"
#include "solve_lp.h"
#include "solve_ip.h"


void print_usage(char* com_name) {
  fprintf(stderr, "\n%s  EXact Linear Programming solver (c) 2002-2006 Masashi Kiyomi\n\n\n", com_name);
  fprintf(stderr, "usage: %s [options] [mps_file]\n\n", com_name);
  fprintf(stderr, "  -h, --help          print this help information\n");
  fprintf(stderr, "  -v, --version       print the program version\n\n");
  fprintf(stderr, "  --print-solution    print solution (some preprosessing skipped)\n");
  fprintf(stderr, "  --verbose           verbose output mode\n\n");
  fprintf(stderr, "  -m, --min           minimize objective (defalt)\n");
  fprintf(stderr, "  -M, --max           maximize objective\n\n");
  fprintf(stderr, "  --simplex-pp        use primal simplex for both phases (default)\n");
  fprintf(stderr, "  --simplex-pd        use primal phase I & dual phase II simplex method\n\n");
  fprintf(stderr, "  --perturbation      use lexicographic method when degenerate\n");
  fprintf(stderr, "  --no-scaling        without scaling\n");
  fprintf(stderr, "  --no-preprocess     without preprocessing\n");
  fprintf(stderr, "  --extra-preprocess  do extra preprocessing\n\n");
#ifndef NO_GMP_HASH
  fprintf(stderr, "  --hash-entries <n>  default is 20011; should be a positive prime\n\n");
#endif
  fprintf(stderr, "  --bland             use Bland's minimal index rule\n");
  fprintf(stderr, "  --bland2            use Bland's recursive rule\n");
  fprintf(stderr, "  --dantzig           use Dantzig's rule\n");
  fprintf(stderr, "  --devex             use Harris' Devex rule\n");
  fprintf(stderr, "  --steepest-edge     use Steepest Edge rule (default)\n");
  fprintf(stderr, "  --projected-steepest-edge  use Projected Steepest Edge rule\n\n");
  fprintf(stderr, "  --markowitz         use Markowitz's rule for LU factorization (default)\n");
  fprintf(stderr, "  --minimum-degree    use Minimum Degree Ordering for LU factorization\n");
  fprintf(stderr, "  --reinversion <n>   specify the reinversion cycle (if 0, exlp dinamically decide the timing, defalt = 0)\n\n");
  fprintf(stderr, "  --gomory            assume integer model, and solve using Gomory cuts (not implemented yet)\n\n");
}

void command_line_options(LP* lp, int argc, char* argv[]) {
  int  i;

  for (i = 1; i < argc; i ++) {
    if (argv[i][0] != '-') {
      if (!(lp->fp = fopen(argv[i], "r"))) {
        fprintf(stderr, "can't open '%s'\n", argv[i]);
        exit(EXIT_FAILURE);
      }
      return;
    }
    if (strncmp(&argv[i][1], "h", 32) == 0 ||
        strncmp(&argv[i][1], "-help", 32) == 0) {
      print_usage(argv[0]);
      exit(EXIT_SUCCESS);

    } else if (strncmp(&argv[i][1], "v", 32) == 0 ||
        strncmp(&argv[i][1], "-version", 32) == 0) {
      fprintf(stderr, "%s\n", EXLP_VERSION);
      exit(EXIT_SUCCESS);

    } else if (strncmp(&argv[i][1], "-print-solution", 32) == 0) {
      lp->print_sol = TRUE;
      lp->scaling = FALSE;

    } else if (strncmp(&argv[i][1], "-verbose", 32) == 0) {
      lp->verbose = TRUE;

    } else if (strncmp(&argv[i][1], "m", 32) == 0 ||
        strncmp(&argv[i][1], "-min", 32) == 0) {
      lp->maximize = FALSE;

    } else if (strncmp(&argv[i][1], "M", 32) == 0 ||
        strncmp(&argv[i][1], "-max", 32) == 0) {
      lp->maximize = TRUE;

    } else if (strncmp(&argv[i][1], "-simplex-pp", 32) == 0) {
      lp->dual_simplex = FALSE;

    } else if (strncmp(&argv[i][1], "-simplex-pd", 32) == 0) {
      lp->dual_simplex = TRUE;

    } else if (strncmp(&argv[i][1], "-avis", 32) == 0) {
      lp->perturbation = 2;

    } else if (strncmp(&argv[i][1], "-perturbation", 32) == 0) {
      lp->perturbation = TRUE;

    } else if (strncmp(&argv[i][1], "-no-scaling", 32) == 0) {
      lp->scaling = FALSE;

    } else if (strncmp(&argv[i][1], "-no-preprocess", 32) == 0) {
      lp->preprocess = 0;

    } else if (strncmp(&argv[i][1], "-extra-preprocess", 32) == 0) {
      lp->preprocess = 2;

#ifndef NO_GMP_HASH
    } else if (strncmp(&argv[i][1], "-hash-entries", 32) == 0) {
      lp->hash_entries = atoi(argv[++i]);
      if (lp->hash_entries <= 0)
        lp->hash_entries = 1;
#endif

    } else if (strncmp(&argv[i][1], "-bland", 32) == 0) {
      lp->pivot_rule = LP_PIVOT_BLAND;

    } else if (strncmp(&argv[i][1], "-bland2", 32) == 0) {
      lp->pivot_rule = LP_PIVOT_BLAND2;

    } else if (strncmp(&argv[i][1], "-dantzig", 32) == 0) {
      lp->pivot_rule = LP_PIVOT_DANTZIG;

    } else if (strncmp(&argv[i][1], "-devex", 32) == 0) {
      lp->pivot_rule = LP_PIVOT_DEVEX;

    } else if (strncmp(&argv[i][1], "-steepest-edge", 32) == 0) {
      lp->pivot_rule = LP_PIVOT_STEEPEST_EDGE;

    } else if (strncmp(&argv[i][1], "-projected-steepest-edge", 32) == 0) {
      lp->pivot_rule = LP_PIVOT_PROJECTED_STEEPEST_EDGE;

    } else if (strncmp(&argv[i][1], "-markowitz", 32) == 0) {
      lp->lu_rule = LU_MARKOWITZ;

    } else if (strncmp(&argv[i][1], "-minimum-degree", 32) == 0) {
      lp->lu_rule = LU_MINIMUM_DEGREE;

    } else if (strncmp(&argv[i][1], "-reinversion", 32) == 0) {
      lp->reinversion_cycle = atoi(argv[++i]);
      if (lp->reinversion_cycle < 0 || lp->reinversion_cycle == ETA_MAX) {
	fprintf(stderr, "invalid reinversion cycle (%d)\n", lp->reinversion_cycle);
	exit(EXIT_FAILURE);
      }

    } else if (strncmp(&argv[i][1], "-gomory", 32) == 0) {
      lp->gomory = TRUE;

    } else {
      print_usage(argv[0]);
      exit(EXIT_FAILURE);
    }
  }

  lp->fp = stdin;
}

// DIMACS.MPS
// ex5A1T.mps
// ex5.mps
//int __cdecl main(int argc, char* argv[]) {


/*  commented out for exlp-package


int main(int argc, char* argv[]) {
  LP* lp;
  int  result;

  mylib_init();
  lp = new_lp(NULL);
  command_line_options(lp, argc, argv);
  lp_read_mps(lp);

  printf("model name : '%s'\n", lp->name);
  printf("variables  : %d\n", lp->vars);
  printf("rows       : %d\n", lp->rows);
  if (lp->gomory)
    result = solve_ip(lp);
  else
    result = solve_lp(lp);

  if (result == LP_RESULT_OPTIMAL) {
    lp_get_object_value(lp, &lp->q_work);
    //if (lp->b->i[lp->b->nonzeros-1] == lp->rows)
    //  mympq_sub(lp->q_work, lp->q_work, lp->b->value[lp->b->nonzeros-1]);
    printf("objective  : %s('%s')\n", (lp->maximize ? "MAX" : "MIN"), lp->obj_name); 
    printf("optimum    : ");  print_rational(lp->q_work);
    putchar('\n');
    printf("in decimal : ");  print_rational_as_float(lp->q_work, 12);
    putchar('\n');
    if (lp->verbose)
      printf("size       : %d/%d\n", mpz_size(mpq_numref(lp->q_work)), mpz_size(mpq_denref(lp->q_work)));
    if (lp->print_sol)
      lp_print_nonzero_vars(lp);

  } else if ((result == LP_RESULT_UNBOUNDED) || (result == LP_RESULT_DUAL_INFEASIBLE)) {
    printf("unbounded\n");

  } else if (result == LP_RESULT_INFEASIBLE) {
    printf("infeasible\n");

  } else if (result == LP_RESULT_UNSUPPORTED) {
    printf("not supported to solve this type of problems.\n");
  }

  lp_free(lp);

  return 0;
}

*/



