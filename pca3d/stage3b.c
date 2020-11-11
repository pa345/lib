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
#include <mainlib/ml_lapack_wrapper.h>

#include "io.h"
#include "pca3d.h"
#include "tiegcm3d.h"

#include "common.c"

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
  const size_t N = data->nlm_complex * data->nr;
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
                  size_t lmidx_c, row_idx, grid_idx;
                  gsl_complex Qqt, Qp, Qq;

                  lmidx_c = lmidx_complex(l, m, data->mmax);
                  row_idx = CIDX2(ir, data->nr, lmidx_c, data->nlm_complex);

                  grid_idx = CIDX4(t, T, ifreq, nfreq, ir, data->nr, lmidx, data->nlm);
                  Qqt = data->Qqt[grid_idx];
                  Qp = data->Qp[grid_idx];
                  Qq = data->Qq[grid_idx];

                  gsl_matrix_complex_set(X, row_idx, t, Qqt);
                  gsl_matrix_complex_set(X, row_idx + N, t, Qp);
                  gsl_matrix_complex_set(X, row_idx + 2*N, t, Qq);

                  if (m > 0)
                    {
                      if (ifreq != 0)
                        {
                          /*
                           * compute Q_{n,-m}(omega) = conj(Q_{nm}(-omega))
                           * negative frequencies can be obtained by substituting nfreq - ifreq
                           */
                          grid_idx = CIDX4(t, T, nfreq - ifreq, nfreq, ir, data->nr, lmidx, data->nlm);
                          Qqt = data->Qqt[grid_idx];
                          Qp = data->Qp[grid_idx];
                          Qq = data->Qq[grid_idx];
                        }

                      lmidx_c = lmidx_complex(l, -m, data->mmax);
                      row_idx = CIDX2(ir, data->nr, lmidx_c, data->nlm_complex);

                      gsl_matrix_complex_set(X, row_idx, t, gsl_complex_conjugate(Qqt));
                      gsl_matrix_complex_set(X, row_idx + N, t, gsl_complex_conjugate(Qp));
                      gsl_matrix_complex_set(X, row_idx + 2*N, t, gsl_complex_conjugate(Qq));
                    }
                }
            }
        }
    }

  return s;
}

/*
proc_svd()
  Compute SVD of spectral density matrix, combining a set of frequencies

Inputs iband       - number of this frequency band for output file
       idx_start   - starting index corresponding to first frequency in [0,nfreq-1]
       idx_end     - end index corresponding to last frequency in [0,nfreq-1]
       data        - pca3d data
       nominalFreq - (output) frequency of center of band (cpd)
*/

int
proc_svd(const size_t iband, const size_t idx_start, const size_t idx_end,
         const pca3d_fft_data * data, double * nominalFreq)
{
  int status = 0;
  const double window_size = data->window_size;
  const double window_shift = data->window_shift;
  const size_t nfreq = idx_end - idx_start + 1;              /* number of frequencies to combine */
  const size_t N = data->nr * data->nlm_complex;             /* spatial grid size */
  const size_t T = data->T;                                  /* number of time window segments */
  gsl_matrix_complex *X, *U, *V;
  gsl_vector *S;
  size_t i;
  struct timeval tv0, tv1;
  char buf[2048];

  *nominalFreq = 0.0;

  /* X = [ Q[start] ... Q[end] ] */
  X = gsl_matrix_complex_alloc(3 * N, nfreq * T);

  S = gsl_vector_alloc(GSL_MIN(X->size1, X->size2));
  U = gsl_matrix_complex_alloc(X->size1, X->size2);
  V = gsl_matrix_complex_alloc(X->size2, X->size2);

  for (i = 0; i < nfreq; ++i)
    {
      gsl_matrix_complex_view Xv = gsl_matrix_complex_submatrix(X, 0, i * T, 3 * N, T);
      size_t ifreq = i + idx_start;
      double freq = ifreq / window_size;

      /* build Q[ifreq] and store in X */
      fprintf(stderr, "proc_svd: building matrix X (%zu-by-%zu) for frequency %.3f [cpd] (index = %zu)...",
              3 * N, T, freq, ifreq);
      gettimeofday(&tv0, NULL);
      build_X(ifreq, data, &Xv.matrix);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

      *nominalFreq += freq;
    }

  *nominalFreq /= (double) nfreq;

  fprintf(stderr, "proc_svd: performing SVD of X (%zu-by-%zu) for nominal frequency %g [cpd]...",
          X->size1, X->size2, *nominalFreq);
  gettimeofday(&tv0, NULL);
  status = lapack_complex_svd_thin(X, S, U, V);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds, condition number = %g, status = %d)\n",
          time_diff(tv0, tv1), gsl_vector_get(S, 0) / gsl_vector_get(S, T - 1), status);

  sprintf(buf, "%s_%zu", PCA3D_STAGE3B_SVAL_TXT, iband);
  fprintf(stderr, "proc_svd: writing singular values for nominal frequency %g [cpd] in text format to %s...", *nominalFreq, buf);
  pca3d_write_S(buf, data->lmax, data->mmax, *nominalFreq, window_size, window_shift, S);
  fprintf(stderr, "done\n");

  sprintf(buf, "%s_%zu", PCA3D_STAGE3B_SVAL_DAT, iband);
  fprintf(stderr, "proc_svd: writing singular values for nominal frequency %g [cpd] in binary format to %s...", *nominalFreq, buf);
  pca3d_write_vector(buf, S);
  fprintf(stderr, "done\n");

  sprintf(buf, "%s_%zu", PCA3D_STAGE3B_U, iband);
  fprintf(stderr, "proc_svd: writing U matrix for nominal frequency %g [cpd] in binary format to %s...", *nominalFreq, buf);
  pca3d_write_matrix_complex(buf, U);
  fprintf(stderr, "done\n");

  sprintf(buf, "%s_%zu", PCA3D_STAGE3B_V_TXT, iband);
  fprintf(stderr, "main: writing right singular vectors for nominal frequency %g [cpd] in text format to %s...", *nominalFreq, buf);
  pca3d_write_complex_V(buf, data->lmax, data->mmax, *nominalFreq, window_size, window_shift, V);
  fprintf(stderr, "done\n");

  gsl_matrix_complex_free(X);
  gsl_matrix_complex_free(U);
  gsl_matrix_complex_free(V);
  gsl_vector_free(S);

  return status;
}

int
main(int argc, char *argv[])
{
  pca3d_fft_data data;
  struct timeval tv0, tv1;
  char *infile = PCA3D_STAGE2B_FFT_DATA;
  double bandFreqs[50]; /* band center frequencies (cpd) */
  gsl_vector_view v;

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
  fprintf(stderr, "main: nlm = %zu\n", data.nlm_complex);

  proc_svd(1, 47, 49, &data, &bandFreqs[0]);
  proc_svd(2, 42, 46, &data, &bandFreqs[1]);
  proc_svd(3, 39, 41, &data, &bandFreqs[2]);
  proc_svd(4, 34, 38, &data, &bandFreqs[3]);
  proc_svd(5, 31, 33, &data, &bandFreqs[4]);
  proc_svd(6, 26, 30, &data, &bandFreqs[5]);
  proc_svd(7, 23, 25, &data, &bandFreqs[6]);
  proc_svd(8, 18, 22, &data, &bandFreqs[7]);
  proc_svd(9, 15, 17, &data, &bandFreqs[8]);
  proc_svd(10, 10, 14, &data, &bandFreqs[9]);
  proc_svd(11, 7, 9, &data, &bandFreqs[10]);
  proc_svd(12, 4, 6, &data, &bandFreqs[11]);
  proc_svd(13, 2, 3, &data, &bandFreqs[12]);

  /* write band center frequencies to output */
  fprintf(stderr, "main: writing band center frequencies to %s...", PCA3D_STAGE3B_BANDS);
  v = gsl_vector_view_array(bandFreqs, 13);
  pca3d_write_vector(PCA3D_STAGE3B_BANDS, &v.vector);
  fprintf(stderr, "done\n");

  free(data.t);
  free(data.r);
  free(data.Qq);
  free(data.Qqt);
  free(data.Qp);
  gsl_vector_free(data.window);

  return 0;
}
