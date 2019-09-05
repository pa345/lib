/*
 * stage3.c
 *
 * 1. Print cumulative variance curve to determine how
 *    many eigenvectors are needed
 * 2. Solve U*alpha = Q where U is a matrix of principal
 *    eigenvectors
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <getopt.h>
#include <complex.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

#include <mainlib/ml_common.h>
#include <mainlib/ml_green.h>
#include <mainlib/ml_lapack_wrapper.h>

#include "io.h"

static double
norm_fro(const gsl_matrix_complex * m)
{
  const size_t M = m->size1;
  const size_t N = m->size2;
  size_t i, j;
  double sum = 0.0;

  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          gsl_complex z = gsl_matrix_complex_get(m, i, j);
          sum += gsl_complex_abs(z);
        }
    }

  sum = sqrt(sum);

  return sum;
}

int
print_variance(const char *filename, const gsl_vector *eval)
{
  int s = 0;
  const size_t n = eval->size;
  FILE *fp = fopen(filename, "w");
  size_t i;
  double sum_all = 0.0;
  double cumsum = 0.0;

  i = 1;
  fprintf(fp, "# Field %zu: eigenvalue number i\n", i++);
  fprintf(fp, "# Field %zu: ( sum_{j=1}^i lambda_j ) / ( sum_j lambda_j )\n", i++);

  /* compute sum_i lambda_i */
  for (i = 0; i < n; ++i)
    {
      double lambda = gsl_vector_get(eval, i);
      sum_all += lambda;
    }

  for (i = 0; i < n; ++i)
    {
      double lambda = gsl_vector_get(eval, i);
      cumsum += lambda;
      fprintf(fp, "%zu %.12e\n", i + 1, cumsum / sum_all);
    }

  fclose(fp);

  return s;
}

/*
solve_PCA()
  Solve PCA problem by computing alpha, such that

U alpha = Q

in a least squares sense, where U are the principal eigenvectors
of the SDM

Inputs: P        - number of eigenvectors of SDM to use
        Q        - Fourier coefficient matrix, nnm-by-T
        evec     - eigenvector matrix, nnm-by-nnm
        alpha    - (output) alpha matrix, P-by-T
        Qt       - (output) Q~ = U*alpha, nnm-by-T
*/

int
solve_PCA(const size_t P, const gsl_matrix_complex * Q,
          const gsl_matrix_complex * evec, gsl_matrix_complex * alpha,
          gsl_matrix_complex * Qt)
{
  int status = 0;
  const size_t T = Q->size2; /* number of time segments */
  const size_t nnm = evec->size1;
  gsl_matrix_complex *R;     /* R = Q - U*alpha */
  struct timeval tv0, tv1;
  double residual;           /* || Q - U*alpha || */
  int rank;

  /* select largest P eigenvectors of SDM */
  gsl_matrix_complex_const_view U = gsl_matrix_complex_const_submatrix(evec, 0, 0, nnm, P);

  /* solve: U*alpha = Q */
  fprintf(stderr, "solve_PCA: solving PCA problem for alpha...");
  gettimeofday(&tv0, NULL);
  status = lapack_complex_lls(&U.matrix, Q, alpha, &rank);
  gettimeofday(&tv1, NULL);

  /* compute: Q~ = U*alpha */
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, &U.matrix, alpha,
                 GSL_COMPLEX_ZERO, Qt);

  /* compute: R = Q - Q~ */
  R = gsl_matrix_complex_alloc(nnm, T);
  gsl_matrix_complex_memcpy(R, Q);
  gsl_matrix_complex_sub(R, Qt);

  residual = norm_fro(R);

  fprintf(stderr, "done (%g seconds, status = %d, rank = %d, residual = %.12e)\n",
          time_diff(tv0, tv1), status, rank, residual);

  gsl_matrix_complex_free(R);

  return status;
}

int
main(int argc, char *argv[])
{
  const size_t nmax = 60;
  const size_t mmax = GSL_MIN(nmax, 30);
  const double R = R_EARTH_KM;
  green_workspace *green_p = green_alloc(nmax, mmax, R);
  char *knm_file = "data/stage1_knm.dat";
  char *eval_file = "data/stage2_eval.dat";
  char *evec_file = "data/stage2_evec.dat";
  char *Q_file = "data/stage2_Q.dat";

  char *variance_file = "variance.txt";

  /* nnm-by-T matrix storing power at frequency omega for each time segment and each
   * (n,m) channel */
  gsl_matrix_complex *Q;

  gsl_matrix_complex *evec;  /* eigenvectors of SDM */
  gsl_vector *eval;          /* eigenvalues of SDM */
  gsl_matrix_complex *alpha; /* alpha matrix, P-by-T */
  gsl_matrix_complex *Qt;    /* Q~ = U*alpha, nnm-by-T */

  gsl_matrix *knm;           /* knm(t) matrix */

  size_t nnm;
  size_t T;                  /* number of time segments */
  size_t P;                  /* number of principal eigenvectors to use (<= T) */

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          default:
            fprintf(stderr, "Usage: %s <-i stage1_matrix_file> [-t time_segment_length (days)]\n", argv[0]);
            break;
        }
    }

  fprintf(stderr, "main: reading knm matrix from %s...", knm_file);
  knm = pca_read_matrix(knm_file);
  fprintf(stderr, "done (%zu-by-%zu matrix read)\n", knm->size1, knm->size2);

  fprintf(stderr, "main: reading eigenvalues from %s...", eval_file);
  eval = pca_read_vector(eval_file);
  fprintf(stderr, "done (%zu eigenvalues read)\n", eval->size);

  fprintf(stderr, "main: reading eigenvectors from %s...", evec_file);
  evec = pca_read_matrix_complex(evec_file);
  fprintf(stderr, "done (%zu-by-%zu matrix read)\n", evec->size1, evec->size2);

  fprintf(stderr, "main: reading Q matrix from %s...", Q_file);
  Q = pca_read_matrix_complex(Q_file);
  fprintf(stderr, "done (%zu-by-%zu matrix read)\n", Q->size1, Q->size2);

  /* plot a variance curve to help decide how many eigenvectors to keep */
  fprintf(stderr, "main: writing variance curve to %s...", variance_file);
  print_variance(variance_file, eval);
  fprintf(stderr, "done\n");

  T = Q->size2;
  P = GSL_MIN(11, T);
  nnm = evec->size1;

  fprintf(stderr, "main: using %zu largest eigenvectors\n", P);

  alpha = gsl_matrix_complex_alloc(P, T);
  Qt = gsl_matrix_complex_alloc(nnm, T);

  /* find alpha such that || Q - U*alpha || is
   * minimized in a least squares sense */
  solve_PCA(P, Q, evec, alpha, Qt);

  {
    const size_t nt = knm->size2;
    const size_t cidx = green_nmidx(3, 1, green_p);
    const double omega = 2.0*M_PI / 24.0; /* 1 cpd frequency in hr^{-1} */
    size_t i;

    for (i = 0; i < nt; ++i)
      {
        double t = (double) i;
        double knmt = gsl_matrix_get(knm, cidx, i);
        size_t k = (size_t) ((t / nt) * T); /* time segment bin */
        gsl_complex expterm = gsl_complex_rect(cos(omega*t), sin(omega*t));
        gsl_complex Qtnm = gsl_matrix_complex_get(Qt, cidx, k);
        gsl_complex val = gsl_complex_mul(Qtnm, expterm);

        printf("%f %f %f\n",
               t / 24.0,
               knmt,
               GSL_REAL(val));
      }
  }

  gsl_matrix_complex_free(Q);
  gsl_matrix_complex_free(evec);
  gsl_matrix_complex_free(alpha);
  gsl_matrix_complex_free(Qt);
  gsl_matrix_free(knm);
  gsl_vector_free(eval);
  green_free(green_p);

  return 0;
}
