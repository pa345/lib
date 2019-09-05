/*
 * stage2b.c
 *
 * 1. Read in spherical harmonic time series knm from stage1 matrix file (nnm-by-nt)
 * 2. Subtract mean from each knm(t) time series
 * 2. Compute SVD of centered knm matrix
 * 3. Write singular values and singular vectors to disk
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <getopt.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>

#include <fftw3.h>
#include <lapacke/lapacke.h>

#include <mainlib/ml_common.h>
#include <mainlib/ml_oct.h>
#include <mainlib/ml_green.h>
#include <mainlib/ml_lapack_wrapper.h>

#include "io.h"
#include "pca.h"

#include "pca_common.c"

/* define to subtract mean from knm(t) time series prior to SVD */
#define SUBTRACT_MEAN           1

/* define to taper the high degree SH coefficients with a cosine taper (to correct
 * TIEGCM ringing issue) */
#define TAPER_COEFFS            1

/*
subtract_mean()
  Subtract mean values from each row of A and store
means in mu vector

Inputs: A  - N-by-M matrix
        mu - vector, length N
*/

int
subtract_mean(gsl_matrix * A, gsl_vector * mu)
{
  const size_t N = A->size1;
  const size_t M = A->size2;
  size_t i;

  for (i = 0; i < N; ++i)
    {
      gsl_vector_view v = gsl_matrix_row(A, i);
      double mean = gsl_stats_mean(v.vector.data, v.vector.stride, M);
      gsl_vector_add_constant(&v.vector, -mean);
      gsl_vector_set(mu, i, mean);
    }

  return 0;
}

int
correlation_matrix(gsl_matrix * A)
{
  const size_t N = A->size1;
  size_t i, j;
  gsl_vector *d = gsl_vector_alloc(N);

  /* store inverse std deviations D = diag(A)^{-1/2} and set diag(A) = 1 */
  for (i = 0; i < N; ++i)
    {
      double *Aii = gsl_matrix_ptr(A, i, i);

      gsl_vector_set(d, i, 1.0 / sqrt(*Aii));
      *Aii = 1.0;
    }

  /* compute correlation matrix corr = D A D, using lower triangle */
  for (i = 0; i < N; ++i)
    {
      double di = gsl_vector_get(d, i);

      for (j = 0; j < i; ++j)
        {
          double dj = gsl_vector_get(d, j);
          double *Aij = gsl_matrix_ptr(A, i, j);

          *Aij *= di * dj;
        }
    }

  gsl_vector_free(d);

  return 0;
}

/*
calc_pca()

Inputs: ut       - UT hour in [0,23]
        K        - matrix of knm(t), nnm-by-nt
        ut_array - array of UT values, size nt
*/

int
calc_pca(const size_t ut, const gsl_matrix * K, const double *ut_array)
{
  int s = 0;
  const size_t nnm = K->size1;
  const size_t nt_all = K->size2;
  size_t i;
  size_t nt = 0; /* number of timestamps for this UT */
  gsl_matrix *A = gsl_matrix_alloc(nnm, nt_all);
  gsl_vector *mu = gsl_vector_calloc(nnm); /* mean values of knm(t) for this UT */
  gsl_matrix_view m;
  struct timeval tv0, tv1;

  fprintf(stderr, "calc_pca: finding knm coefficients corresponding to %zu UT...", ut);

  for (i = 0; i < nt_all; ++i)
    {
      size_t uti = (size_t) ut_array[i];

      if (ut == uti)
        {
          gsl_vector_const_view v = gsl_matrix_const_column(K, i);
          gsl_vector_view w = gsl_matrix_column(A, nt);

          /* copy set of knm coefficients for this UT hour into the A matrix */
          gsl_vector_memcpy(&w.vector, &v.vector);

          ++nt;
        }
    }

  fprintf(stderr, "done (%zu timestamps found)\n", nt);

  m = gsl_matrix_submatrix(A, 0, 0, nnm, nt);

#if 0
#if SUBTRACT_MEAN
  fprintf(stderr, "calc_pca: subtracting mean from each %zu UT knm time series...", ut);
  subtract_mean(&m.matrix, mu);
  fprintf(stderr, "done\n");
#endif
#endif

  /* compute 1/sqrt(nt) A for SVD computation */
  gsl_matrix_scale(&m.matrix, 1.0 / sqrt((double) nt));

  {
    char sval_txt_file[2048], sval_file[2048], U_file[2048], V_file[2048], mu_file[2048];
    gsl_vector *S = gsl_vector_alloc(GSL_MIN(nnm, nt));
    gsl_matrix *U = gsl_matrix_alloc(nnm, nnm);
    gsl_matrix *V = gsl_matrix_alloc(nt, nt);
    int status;
    FILE *fp;

    sprintf(sval_txt_file, "%s_%02zuUT.txt", PCA_STAGE2B_SVAL_TXT, ut);
    sprintf(sval_file, "%s_%02zuUT.dat", PCA_STAGE2B_SVAL, ut);
    sprintf(U_file, "%s_%02zuUT.dat", PCA_STAGE2B_U, ut);
    sprintf(V_file, "%s_%02zuUT.dat", PCA_STAGE2B_V, ut);
    sprintf(mu_file, "%s_%02zuUT.dat", PCA_STAGE2B_MU, ut);

    fprintf(stderr, "\t performing SVD of matrix for %zu UT (%zu-by-%zu)...", ut, nnm, nt);
    gettimeofday(&tv0, NULL);
    status = lapack_svd(&m.matrix, S, U, V);
    gettimeofday(&tv1, NULL);
    fprintf(stderr, "done (%g seconds, status = %d)\n",
            time_diff(tv0, tv1), status);

    fprintf(stderr, "\t writing singular values in text format to %s...", sval_txt_file);
    fp = fopen(sval_txt_file, "w");
    gsl_vector_fprintf(fp, S, "%.12e");
    fclose(fp);
    fprintf(stderr, "done\n");

    fprintf(stderr, "\t writing singular values in binary format to %s...", sval_file);
    pca_write_vector(sval_file, S);
    fprintf(stderr, "done\n");

    fprintf(stderr, "\t writing left singular vectors in binary format to %s...", U_file);
    pca_write_matrix(U_file, U);
    fprintf(stderr, "done\n");

    fprintf(stderr, "\t writing right singular vectors in binary format to %s...", V_file);
    pca_write_matrix(V_file, V);
    fprintf(stderr, "done\n");

    fprintf(stderr, "\t writing mean values in binary format to %s...", mu_file);
    pca_write_vector(mu_file, mu);
    fprintf(stderr, "done\n");

    gsl_matrix_free(U);
    gsl_matrix_free(V);
    gsl_vector_free(S);
  }

  gsl_matrix_free(A);
  gsl_vector_free(mu);

  return s;
}

int
main(int argc, char *argv[])
{
  size_t nnm;                      /* number of spherical harmonic time series */
  size_t nt;                       /* number of time steps */
  char *infile = PCA_STAGE1_KNM;
  gsl_matrix *K;
  struct timeval tv0, tv1;
  size_t ut;                       /* UT hour */
  double *ut_array;                /* array of UT hours, size nt */
  size_t nmax, mmax;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "i:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          default:
            fprintf(stderr, "Usage: %s <-i stage1_matrix_file>\n", argv[0]);
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i stage1_matrix_file>\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "input file = %s\n", infile);

  fprintf(stderr, "main: reading %s...", infile);
  K = pca_read_matrix(infile);
  fprintf(stderr, "done (%zu-by-%zu matrix)\n", K->size1, K->size2);

  nnm = K->size1;
  nt = K->size2;

  ut_array = malloc(nt * sizeof(double));
  pca_read_data(PCA_STAGE1_DATA, &nmax, &mmax, &nt, ut_array);
  assert(nt == K->size2);

#if TAPER_COEFFS

  /* taper spherical harmonic coefficients to eliminate ringing effect from TIEGCM */
  {
    green_workspace *green_p = green_alloc(60, 30, R_EARTH_KM);
    char *spectrum_file = "spectrum_taper.s";
    gsl_vector_view x = gsl_matrix_column(K, 0);

    fprintf(stderr, "main: tapering knm coefficients...");
    taper_knm(20, K, green_p);
    fprintf(stderr, "done\n");

    fprintf(stderr, "main: writing tapered spectrum to %s...", spectrum_file);
    green_print_spectrum(spectrum_file, &x.vector, green_p);
    fprintf(stderr, "done\n");

    green_free(green_p);
  }

#endif

#if SUBTRACT_MEAN
  {
    gsl_vector *mu = gsl_vector_alloc(nnm);

    fprintf(stderr, "main: subtracting mean from each knm time series...");
    subtract_mean(K, mu);
    fprintf(stderr, "done\n");

    gsl_vector_free(mu);
  }
#endif

#if 1
  /* loop over UT and compute PCA modes for each UT hour */
  for (ut = 0; ut < 24; ++ut)
    {
      calc_pca(ut, K, ut_array);
    }
  exit(1);
#endif

  /* compute 1/sqrt(nt) K for SVD computation */
  gsl_matrix_scale(K, 1.0 / sqrt((double) nt));

  {
    const char *sval_txt_file = "sval_time.txt";
    const char *sval_file = PCA_STAGE2B_SVAL;
    const char *U_file = PCA_STAGE2B_U;
    const char *V_file = PCA_STAGE2B_V;
    gsl_vector *S = gsl_vector_alloc(GSL_MIN(nnm, nt));
    gsl_matrix *U = gsl_matrix_alloc(nnm, nnm);
    gsl_matrix *V = gsl_matrix_alloc(nt, nt);
    int status;
    FILE *fp;

    fprintf(stderr, "main: performing SVD of K matrix (%zu-by-%zu)...", nnm, nt);
    gettimeofday(&tv0, NULL);
    status = lapack_svd(K, S, U, V);
    gettimeofday(&tv1, NULL);
    fprintf(stderr, "done (%g seconds, status = %d)\n",
            time_diff(tv0, tv1), status);

    fprintf(stderr, "main: writing singular values in text format to %s...",
            sval_txt_file);
    fp = fopen(sval_txt_file, "w");
    gsl_vector_fprintf(fp, S, "%.12e");
    fclose(fp);
    fprintf(stderr, "done\n");

    fprintf(stderr, "main: writing singular values in binary format to %s...",
            sval_file);
    pca_write_vector(sval_file, S);
    fprintf(stderr, "done\n");

    fprintf(stderr, "main: writing left singular vectors in binary format to %s...",
            U_file);
    pca_write_matrix(U_file, U);
    fprintf(stderr, "done\n");

    fprintf(stderr, "main: writing right singular vectors in binary format to %s...",
            V_file);
    pca_write_matrix(V_file, V);
    fprintf(stderr, "done\n");

    gsl_matrix_free(U);
    gsl_matrix_free(V);
    gsl_vector_free(S);
  }

  gsl_matrix_free(K);
  free(ut_array);

  return 0;
}
