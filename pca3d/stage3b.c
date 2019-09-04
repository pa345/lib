/*
 * stage3b.c
 *
 * Use results of FFT analysis (stage2b) on SH coefficients to build
 * spectral density matrix and compute SVD
 *
 * ./stage3b [-i fft_data_file]
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <assert.h>
#include <errno.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <mainlib/ml_common.h>
#include <mainlib/ml_bsearch.h>

#include "lapack_wrapper.h"

#include "io.h"
#include "pca3d.h"
#include "tiegcm3d.h"

/*
build_X()

  Build matrix X (3N-by-T, with N = nlm * nr)
for a given frequency

X(w_p) = [ X_1(w_p) X_2(w_p) ... X_T(w_p) ]

with length 3N column vectors

X_t(w_p) = [ Qqt^{t}_l^m(r_i,w_p) ]
           [ Qp^{t}_l^m(r_i,w_p) ]
           [ Qq^{t}_l^m(r_i,w_p) ]

Inputs: ifreq - index of desired frequency
        data  - fft data
        X     - (output) matrix, 3N-by-T

Return: success/error
*/

int
build_X(const size_t ifreq, const pca3d_fft_data * data, gsl_matrix_complex * X)
{
  int s = 0;
  const size_t N = data->nlm * data->nr;
  const size_t T = data->T;               /* number of time window segments */
  const size_t nfreq = data->nfreq;       /* FFT output buffer size (number of frequencies) */
  size_t t;                               /* window index \in [0,T-1] */

  /*
   * The time window index t is the column index of X; the indices
   * (ir,l,m) give the row index of X
   */
  for (t = 0; t < T; ++t)
    {
      size_t ir;

      for (ir = 0; ir < data->nr; ++ir)
        {
          size_t l;

          for (l = data->lmin; l <= data->lmax; ++l)
            {
              int M = (int) GSL_MIN(l, data->mmax);
              int m;

              for (m = 0; m <= M; ++m)
                {
                  size_t lmidx = magfield_lmidx(l, m, data->mmax);
                  size_t row_idx = CIDX2(ir, data->nr, lmidx, data->nlm);
                  size_t grid_idx = CIDX4(t, T, ifreq, nfreq, ir, data->nr, lmidx, data->nlm);
                  gsl_complex Qqt = data->Qqt[grid_idx];
                  gsl_complex Qp = data->Qp[grid_idx];
                  gsl_complex Qq = data->Qq[grid_idx];

                  gsl_matrix_complex_set(X, row_idx, t, Qqt);
                  gsl_matrix_complex_set(X, row_idx + N, t, Qp);
                  gsl_matrix_complex_set(X, row_idx + 2*N, t, Qq);
                }
            }
        }
    }

  return s;
}

int
proc_svd(const size_t ifreq, const pca3d_fft_data * data)
{
  int status;
  const double window_size = data->window_size;
  const double window_shift = data->window_shift;
  const double freq = ifreq / window_size;                   /* desired frequency in cpd */
  const size_t N = data->nr * data->nlm;                     /* spatial grid size */
  const size_t T = data->T;                                  /* number of time window segments */
  gsl_matrix_complex *X = gsl_matrix_complex_alloc(3 * N, T);
  gsl_vector *S = gsl_vector_alloc(T);
  gsl_matrix_complex *U = gsl_matrix_complex_alloc(3 * N, T);
  gsl_matrix_complex *V = gsl_matrix_complex_alloc(T, T);
  struct timeval tv0, tv1;
  char buf[2048];

  fprintf(stderr, "proc_svd: === PROCESSING frequency %.3f [cpd] [index: %zu/%zu] ===\n", freq, ifreq, data->nfreq);

  fprintf(stderr, "proc_svd: building matrix X (%zu-by-%zu) for frequency %.3f [cpd] (index = %zu)...",
          X->size1, X->size2, freq, ifreq);
  gettimeofday(&tv0, NULL);
  build_X(ifreq, data, X);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "proc_svd: performing SVD of X for frequency %g [cpd]...", freq);
  gettimeofday(&tv0, NULL);
  status = lapack_complex_svd_thin(X, S, U, V);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds, condition number = %g, status = %d)\n",
          time_diff(tv0, tv1), gsl_vector_get(S, 0) / gsl_vector_get(S, T - 1), status);

  sprintf(buf, "%s_%zu", PCA3D_STAGE3B_SVAL_TXT, ifreq);
  fprintf(stderr, "proc_svd: writing singular values for frequency %g [cpd] in text format to %s...", freq, buf);
  pca3d_write_S(buf, data->lmax, data->mmax, freq, window_size, window_shift, S);
  fprintf(stderr, "done\n");

  sprintf(buf, "%s_%zu", PCA3D_STAGE3B_SVAL_DAT, ifreq);
  fprintf(stderr, "proc_svd: writing singular values for frequency %g [cpd] in binary format to %s...", freq, buf);
  pca3d_write_vector(buf, S);
  fprintf(stderr, "done\n");

  sprintf(buf, "%s_%zu", PCA3D_STAGE3B_U, ifreq);
  fprintf(stderr, "proc_svd: writing U matrix for frequency %g [cpd] in binary format to %s...", freq, buf);
  pca3d_write_matrix_complex(buf, U);
  fprintf(stderr, "done\n");

  gsl_matrix_complex_free(X);
  gsl_matrix_complex_free(U);
  gsl_matrix_complex_free(V);
  gsl_vector_free(S);

  return 0;
}

int
main(int argc, char *argv[])
{
  /*
   * For the PCA project, Gary provided temporal modes with the frequencies (in cpd):
   * , 0.625, 1, 1.5, 2, 2.5, ..., 5.5, 6
   *
   * These frequency bands correspond to the following indices in the FFT output array, where
   *
   * Freq[idx] = idx * Fs / N
   *
   * where Fs = 24 samples/day is the sampling rate, and N = 192 samples/window (8 days) in the window size
   */
  const size_t indices[] = { 3, 5, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48 };
  const size_t nidx = 13;
  pca3d_fft_data data;
  struct timeval tv0, tv1;
  char *infile = PCA3D_STAGE2B_FFT_DATA;
  size_t i;

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
            fprintf(stderr, "Usage: %s [-i fft_data_file]\n", argv[0]);
            exit(1);
            break;
        }
    }

  fprintf(stderr, "main: reading %s...", infile);
  gettimeofday(&tv0, NULL);
  data = pca3d_read_fft_data2(infile);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "main: nr  = %zu\n", data.nr);
  fprintf(stderr, "main: nlm = %zu\n", data.nlm);

  for (i = 0; i < nidx; ++i)
    proc_svd(indices[i], &data);

  free(data.t);
  free(data.r);
  free(data.Qq);
  free(data.Qqt);
  free(data.Qp);
  gsl_vector_free(data.window);

  return 0;
}
