/*
 * invert_postproc.c
 *
 * Use the R and Q^T b factors computed by invert, apply regularization,
 * compute L-curve, output final coefficients
 *
 * Usage:
 * ./invert_postproc [flags]
 *
 * Flags:
 *   -o coef_output_file
 *   -C config_file
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <math.h>
#include <assert.h>
#include <omp.h>
#include <libconfig.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_oct.h>
#include <mainlib/ml_msynth.h>
#include <mainlib/ml_track.h>
#include <mainlib/ml_magdata.h>
#include <mainlib/ml_matio.h>

#include "invert.h"
#include "invert_residual.h"
#include "invert_multifit.c"
#include "invert_synth.h"
#include "io.h"
#include "pca3d.h"

static int
parse_config_file(const char *filename, invert_parameters *invert_params)
{
  int s;
  config_t cfg;
  double fval;
  int ival;
  const char * sval;

  config_init(&cfg);

  s = config_read_file(&cfg, filename);
  if (s != CONFIG_TRUE)
    {
      fprintf(stderr, "parse_config_file: %s:%d - %s\n",
              config_error_file(&cfg),
              config_error_line(&cfg),
              config_error_text(&cfg));
      config_destroy(&cfg);
      return -1;
    }

  if (config_lookup_float(&cfg, "R", &fval))
    invert_params->R = fval;
  if (config_lookup_int(&cfg, "nfreq", &ival))
    invert_params->nfreq = (size_t) ival;

  if (config_lookup_string(&cfg, "tmode_file", &sval))
    strcpy(invert_params->tmode_file, sval);

  if (config_lookup_int(&cfg, "max_iter", &ival))
    invert_params->max_iter = (size_t) ival;
  if (config_lookup_float(&cfg, "qdlat_fit_cutoff", &fval))
    invert_params->qdlat_fit_cutoff = fval;
  if (config_lookup_int(&cfg, "use_weights", &ival))
    invert_params->use_weights = ival;
  if (config_lookup_int(&cfg, "regularize", &ival))
    invert_params->regularize = ival;

  if (config_lookup_float(&cfg, "weight_X", &fval))
    invert_params->weight_X = fval;
  if (config_lookup_float(&cfg, "weight_Y", &fval))
    invert_params->weight_Y = fval;
  if (config_lookup_float(&cfg, "weight_Z", &fval))
    invert_params->weight_Z = fval;
  if (config_lookup_float(&cfg, "weight_F", &fval))
    invert_params->weight_F = fval;
  if (config_lookup_float(&cfg, "weight_DXDT", &fval))
    invert_params->weight_DXDT = fval;
  if (config_lookup_float(&cfg, "weight_DYDT", &fval))
    invert_params->weight_DYDT = fval;
  if (config_lookup_float(&cfg, "weight_DZDT", &fval))
    invert_params->weight_DZDT = fval;
  if (config_lookup_float(&cfg, "weight_DX", &fval))
    invert_params->weight_DX = fval;
  if (config_lookup_float(&cfg, "weight_DY", &fval))
    invert_params->weight_DY = fval;
  if (config_lookup_float(&cfg, "weight_DZ", &fval))
    invert_params->weight_DZ = fval;

  config_destroy(&cfg);

  return 0;
}

int
build_regularization(gsl_vector * reg, invert_workspace * w)
{
  invert_tmode_workspace * tmode_p = w->tmode_workspace_p;
  invert_smode_workspace * smode_p = w->smode_workspace_p;
  size_t nf, ns, nt;
  char buf[1024];

  for (nf = 0; nf < w->nfreq; ++nf)
    {
      gsl_vector *S;

      /* read singular values for this band */
      sprintf(buf, "%s_%zu", PCA3D_STAGE3B_SVAL_DAT, nf + 1);
      S = pca3d_read_vector(buf);

      for (ns = 0; ns < smode_p->nmodes[nf]; ++ns)
        {
          double Sj = gsl_vector_get(S, ns);
          double Sjinv = 1.0 / Sj;

          for (nt = 0; nt < tmode_p->nmodes[nf]; ++nt)
            {
              size_t idx = invert_coeff_idx(nf, nt, ns, w);
              gsl_vector_set(reg, idx, Sjinv);
            }
        }

      gsl_vector_free(S);
    }

  return 0;
}

/*
solve:

min_x || [ QTb ] - [    R     ] x ||^2
      || [  0  ]   [ lambda L ]   ||
*/

int
compute_reg_solution(const double lambda, const gsl_vector * L, const gsl_matrix * R, const gsl_vector * QTb,
                     gsl_vector * x, double * rnorm, double * snorm)
{
  const size_t p = R->size2;
  gsl_matrix * A = gsl_matrix_calloc(2 * p, p);
  gsl_vector * b = gsl_vector_calloc(2 * p);
  gsl_matrix * T = gsl_matrix_alloc(p, p);
  gsl_vector * sol = gsl_vector_alloc(2 * p);
  gsl_vector * work = gsl_vector_alloc(p);
  gsl_matrix_view m;
  gsl_vector_view v;

  /* A = [ R ; 0 ] */
  m = gsl_matrix_submatrix(A, 0, 0, p, p);
  gsl_matrix_tricpy(CblasUpper, CblasNonUnit, &m.matrix, R);

  /* A = [ R; lambda*L ] */
  m = gsl_matrix_submatrix(A, p, 0, p, p);
  v = gsl_matrix_diagonal(&m.matrix);
  gsl_vector_axpby(lambda, L, 0.0, &v.vector);

  /* b = [ QTb ; 0 ] */
  v = gsl_vector_subvector(b, 0, p);
  gsl_vector_memcpy(&v.vector, QTb);

  /* solve least squares problem */
  gsl_linalg_QR_decomp_r(A, T);
  gsl_linalg_QR_lssolve_r(A, T, b, sol, work);

  v = gsl_vector_subvector(sol, 0, p);
  gsl_vector_memcpy(x, &v.vector);
  *snorm = gsl_blas_dnrm2(x);

  v = gsl_vector_subvector(sol, p, p);
  *rnorm = gsl_blas_dnrm2(&v.vector);

  gsl_matrix_free(A);
  gsl_matrix_free(T);
  gsl_vector_free(b);
  gsl_vector_free(sol);
  gsl_vector_free(work);

  return 0;
}

int
compute_Lcurve(const char * filename, const gsl_vector * L, const gsl_matrix * R, const gsl_vector * QTb)
{
  FILE *fp;
  const size_t N = 100; /* number of points on L-curve */
  double lambda_min = 1.0e-1;
  double lambda_max = 1.0e2;
  double lambda = lambda_min;
  double ratio, rnorm, snorm;
  gsl_vector * x = gsl_vector_alloc(QTb->size);
  size_t i;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "compute_Lcurve: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  ratio = pow(lambda_max / lambda_min, 1.0 / (N - 1.0)); 

  i = 1;
  fprintf(fp, "Field %zu: regularization parameter\n", i++);
  fprintf(fp, "Field %zu: solution norm\n", i++);
  fprintf(fp, "Field %zu: residual norm\n", i++);

  for (i = 0; i < N; ++i)
    {
      compute_reg_solution(lambda, L, R, QTb, x, &rnorm, &snorm);
      lambda *= ratio;

      fprintf(fp, "%.12e %.12e %.12e\n", lambda, snorm, rnorm);
    }

  gsl_vector_free(x);

  fclose(fp);

  return 0;
}

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options]\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --config_file | -C file         - configuration file\n");
  fprintf(stderr, "\t --output_coef_file | -o file    - output coefficient file (binary)\n");
  fprintf(stderr, "\t --matrix | -R file              - R factor file\n");
  fprintf(stderr, "\t --rhs | -b file                 - Q^T b file\n");
}

int
main(int argc, char *argv[])
{
  char *config_file = "INVERT.cfg"; /* default config file */
  char *output_coef_file = "coef_reg.txt";
  char *Lcurve_file = "Lcurve.txt";
  invert_workspace *invert_workspace_p;
  invert_parameters invert_params;
  gsl_vector *reg;    /* regularization matrix */
  gsl_vector *coeffs; /* model coefficients */
  gsl_matrix *R = NULL;   /* R factor */
  gsl_vector *QTb = NULL; /* QTb */
  struct timeval tv0, tv1;
  size_t p;

  invert_init_params(&invert_params);

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "config_file", required_argument, NULL, 'C' },
          { "output_coef_file", required_argument, NULL, 'o' },
          { "matrix", required_argument, NULL, 'R' },
          { "rhs", required_argument, NULL, 'b' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "b:C:o:R:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'C':
            config_file = optarg;
            break;

          case 'i':
            output_coef_file = optarg;
            break;

          case 'b':
            fprintf(stderr, "main: reading %s...", optarg);
            QTb = vecread(optarg);
            fprintf(stderr, "done (vector is length %zu)\n", QTb->size);
            break;

          case 'R':
            fprintf(stderr, "main: reading %s...", optarg);
            R = matread(optarg);
            fprintf(stderr, "done (matrix is %zu-by-%zu)\n", R->size1, R->size2);
            break;

          case 'o':
            output_coef_file = optarg;
            break;

          default:
            print_help(argv);
            exit(1);
            break;
        }
    }

  if (output_coef_file == NULL || R == NULL || QTb == NULL)
    {
      print_help(argv);
      exit(1);
    }

#if 0
  /*XXX*/
  {
    p = R->size1;
    coeffs = gsl_vector_alloc(p);
    gsl_vector_memcpy(coeffs, QTb);
    gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, R, coeffs);
    printv_octave(coeffs, "c1");
    printv_octave(QTb, "QTb");
    printtri_octave(R, "R");
    exit(1);
  }
#endif

  /* parse configuration file */
  fprintf(stderr, "main: parsing configuration file %s...", config_file);
  parse_config_file(config_file, &invert_params);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: nfreq             = %zu\n", invert_params.nfreq);
  fprintf(stderr, "main: output coef file  = %s\n", output_coef_file);

  invert_params.nsat = 0;
  invert_params.invert_data_p = NULL;

  /* allocate invert workspace */
  invert_workspace_p = invert_alloc(&invert_params);
  if (invert_workspace_p == NULL)
    {
      fprintf(stderr, "main: invert_alloc failed\n");
      exit(1);
    }

  p = invert_workspace_p->p;
  fprintf(stderr, "main: number of total parameters:                %zu\n", p);

  coeffs = gsl_vector_alloc(p);
  reg = gsl_vector_alloc(p);

  fprintf(stderr, "main: building regularization matrix...");
  build_regularization(reg, invert_workspace_p);
  fprintf(stderr, "done\n");

  printv_octave(reg, "reg");

  fprintf(stderr, "main: computing L-curve...");
  gettimeofday(&tv0, NULL);
  compute_Lcurve(Lcurve_file, reg, R, QTb);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds, output = %s)\n", time_diff(tv0, tv1), Lcurve_file);

  invert_free(invert_workspace_p);
  gsl_vector_free(coeffs);
  gsl_vector_free(reg);

  return 0;
}
