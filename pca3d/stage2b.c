/*
 * stage2b.c
 *
 * 1. Read previously calculated SH coefficients for each timestep (stage1b)
 * 2. Compute windowed FFTs of q~_l^m(t,r), p_l^m(t,r) and q_l^m(t,r) time series
 * 3. Store FFT data to disk
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <sys/time.h>
#include <assert.h>
#include <omp.h>
#include <complex.h>

#include <fftw3.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_test.h>

#include <common/common.h>
#include <common/oct.h>

#include <magfield/magfield.h>

#include "io.h"
#include "pca3d.h"
#include "window.h"

#include "common.c"

/*
copy_sh_grids()
  Copy spherical harmonic coefficients from magfield workspace to temporal grids

Inputs: it      - time index in [0,nt-1]
        grid_q  - (output) where to store q_l^m(t,r) coefficients
        grid_qt - (output) where to store qt_l^m(t,r) coefficients
        grid_p  - (output) where to store p_l^m(t,r) coefficients
        w       - workspace
*/

int
copy_sh_grids(const size_t it, gsl_complex * grid_q, gsl_complex * grid_qt, gsl_complex * grid_p,
              const magfield_workspace * w)
{
  size_t ir, l;

  for (ir = 0; ir < w->nr; ++ir)
    {
      for (l = w->lmin; l <= w->lmax; ++l)
        {
          int M = (int) GSL_MIN(l, w->mmax);
          int m;

          for (m = 0; m <= M; ++m)
            {
              size_t lmidx = magfield_lmidx(l, m, w->mmax);

              /* nt is not used in the CIDX3 macro */
              grid_q[CIDX3(it, nt, ir, w->nr, lmidx, w->nlm)] = gsl_vector_complex_get(w->qcoeff, MAG_COEFIDX(ir, lmidx, w));
              grid_qt[CIDX3(it, nt, ir, w->nr, lmidx, w->nlm)] = gsl_vector_complex_get(w->qtcoeff, MAG_COEFIDX(ir, lmidx, w));
              grid_p[CIDX3(it, nt, ir, w->nr, lmidx, w->nlm)] = gsl_vector_complex_get(w->pcoeff, MAG_COEFIDX(ir, lmidx, w));
            }
        }
    }

  return 0;
}

/*
do_transforms()
  Perform windowed FFTs on SH time series stored in grids

Inputs: data         - data including timestamps, size nt
        fs           - sampling frequency (samples/day)
        window_size  - size of windows in days
        window_shift - amount to shift window forward in days
        grid_q       - grid of q_l^m(t,r) values, nt-by-nr-by-nlm
        grid_qt      - grid of q~_l^m(t,r) values, nt-by-nr-by-nlm
        grid_p       - grid of p_l^m(t,r) values, nt-by-nr-by-nlm
        w            - workspace

Return: success/error
*/

int
do_transforms(const pca3d_data * data, const double fs, const double window_size, const double window_shift,
              gsl_complex *grid_q, gsl_complex *grid_qt, gsl_complex *grid_p,
              const magfield_workspace * w)
{
  int s = 0;
  const size_t nt = data->nt;
  const size_t nwindow = (size_t) (window_size * fs);                /* optimal number of samples per window */
  const size_t nforward = (size_t) (window_shift * fs);              /* number of samples to slide forward */
  const size_t nfreq = nwindow;                                      /* number of frequencies computed by FFT */
  const size_t T = count_windows(nt, fs, window_size, window_shift); /* number of time windows */
  const int max_threads = omp_get_max_threads();
  gsl_vector *window = gsl_vector_alloc(nwindow);                    /* window function */
  struct timeval tv0, tv1;
  size_t i;
  gsl_complex *Qq;                                                   /* FT'd q_l^m(t,r), T-by-nfreq-by-nr-by-nlm */
  gsl_complex *Qqt;                                                  /* FT'd q~_l^m(t,r), T-by-nfreq-by-nr-by-nlm */
  gsl_complex *Qp;                                                   /* FT'd p_l^m(t,r), T-by-nfreq-by-nr-by-nlm */
  pca3d_fft_data fft_data;

  /* thread specific variables */
  gsl_vector_complex **workq = malloc(max_threads * sizeof(gsl_vector_complex *));
  gsl_vector_complex **workqt = malloc(max_threads * sizeof(gsl_vector_complex *));
  gsl_vector_complex **workp = malloc(max_threads * sizeof(gsl_vector_complex *));
  fftw_complex **fft_q = malloc(max_threads * sizeof(fftw_complex *));
  fftw_complex **fft_qt = malloc(max_threads * sizeof(fftw_complex *));
  fftw_complex **fft_p = malloc(max_threads * sizeof(fftw_complex *));
  fftw_plan *plan_q = malloc(max_threads * sizeof(fftw_plan));
  fftw_plan *plan_qt = malloc(max_threads * sizeof(fftw_plan));
  fftw_plan *plan_p = malloc(max_threads * sizeof(fftw_plan));

  /* allocate thread variables */
  for (i = 0; i < (size_t) max_threads; ++i)
    {
      workq[i] = gsl_vector_complex_alloc(nwindow);
      workqt[i] = gsl_vector_complex_alloc(nwindow);
      workp[i] = gsl_vector_complex_alloc(nwindow);

      fft_q[i] = fftw_malloc(sizeof(fftw_complex) * nfreq);
      fft_qt[i] = fftw_malloc(sizeof(fftw_complex) * nfreq);
      fft_p[i] = fftw_malloc(sizeof(fftw_complex) * nfreq);

      plan_q[i] = fftw_plan_dft_1d(nwindow, (fftw_complex *) workq[i]->data, fft_q[i], FFTW_FORWARD, FFTW_ESTIMATE);
      plan_qt[i] = fftw_plan_dft_1d(nwindow, (fftw_complex *) workqt[i]->data, fft_qt[i], FFTW_FORWARD, FFTW_ESTIMATE);
      plan_p[i] = fftw_plan_dft_1d(nwindow, (fftw_complex *) workp[i]->data, fft_p[i], FFTW_FORWARD, FFTW_ESTIMATE);
    }

  Qq = malloc(T * nfreq * w->nr * w->nlm * sizeof(gsl_complex));
  Qqt = malloc(T * nfreq * w->nr * w->nlm * sizeof(gsl_complex));
  Qp = malloc(T * nfreq * w->nr * w->nlm * sizeof(gsl_complex));

  fft_data.nt = nt;
  fft_data.nfreq = nfreq;
  fft_data.nr = w->nr;
  fft_data.lmin = w->lmin;
  fft_data.lmax = w->lmax;
  fft_data.mmax = w->mmax;
  fft_data.nlm = w->nlm;
  fft_data.T = T;
  fft_data.fs = fs;
  fft_data.window_size = window_size;
  fft_data.window_shift = window_shift;
  fft_data.nwindow = nwindow;
  fft_data.window = window;
  fft_data.t = data->t;
  fft_data.r = w->r;
  fft_data.grid_q = grid_q;
  fft_data.grid_qt = grid_qt;
  fft_data.grid_p = grid_p;
  fft_data.Qq = Qq;
  fft_data.Qqt = Qqt;
  fft_data.Qp = Qp;

  fprintf(stderr, "do_transforms: samples per window   = %zu\n", nwindow);
  fprintf(stderr, "do_transforms: sample slide forward = %zu\n", nforward);
  fprintf(stderr, "do_transforms: number of freqs      = %zu\n", nfreq);
  fprintf(stderr, "do_transforms: time samples         = %zu [hourly]\n", nt);
  fprintf(stderr, "do_transforms: window segments (T)  = %zu\n", T);

  /* compute window function */
  /*apply_ps1(NULL, window);*/
  apply_modsinsq(NULL, window);

  fprintf(stderr, "do_transforms: computing FFTs of windowed data...\n");
  gettimeofday(&tv0, NULL);

  for (i = 0; i < w->nr; ++i)
    {
      int thread_id = omp_get_thread_num();
      size_t l;

      progress_bar(stdout, (double) i / (double) w->nr, 80);

      for (l = w->lmin; l <= w->lmax; ++l)
        {
          int M = (int) GSL_MIN(l, w->mmax);
          int m;

          for (m = 0; m <= M; ++m)
            {
              size_t lmidx = magfield_lmidx(l, m, w->mmax);
              size_t start_idx = 0;
              size_t t;

              /* loop over window segments for this (r,l,m) coefficient */
              for (t = 0; t < T; ++t)
                {
                  size_t end_idx = GSL_MIN(start_idx + nwindow - 1, nt - 1);
                  size_t n = end_idx - start_idx + 1; /* size of actual window */
                  size_t it, ifreq;

                  assert(start_idx < end_idx);

                  if (n < nwindow)
                    {
                      /* could happen at the end of the time series; zero pad input buffer */
                      gsl_vector_complex_set_zero(workq[thread_id]);
                      gsl_vector_complex_set_zero(workqt[thread_id]);
                      gsl_vector_complex_set_zero(workp[thread_id]);
                    }

                  /* copy current time window into work arrays */
                  for (it = start_idx; it <= end_idx; ++it)
                    {
                      size_t grid_idx = CIDX3(it, nt, i, w->nr, lmidx, w->nlm);

                      gsl_vector_complex_set(workq[thread_id], it - start_idx, grid_q[grid_idx]);
                      gsl_vector_complex_set(workqt[thread_id], it - start_idx, grid_qt[grid_idx]);
                      gsl_vector_complex_set(workp[thread_id], it - start_idx, grid_p[grid_idx]);
                    }

                  /* apply window function */
                  {
                    size_t jj;

                    for (jj = 0; jj < n; ++jj)
                      {
                        gsl_complex *qptr = gsl_vector_complex_ptr(workq[thread_id], jj);
                        gsl_complex *qtptr = gsl_vector_complex_ptr(workqt[thread_id], jj);
                        gsl_complex *pptr = gsl_vector_complex_ptr(workp[thread_id], jj);
                        double wj = gsl_vector_get(window, jj);

                        *qptr = gsl_complex_mul_real(*qptr, wj);
                        *qtptr = gsl_complex_mul_real(*qtptr, wj);
                        *pptr = gsl_complex_mul_real(*pptr, wj);
                      }
                  }

                  /* compute FFT of this windowed data */
                  fftw_execute(plan_q[thread_id]);
                  fftw_execute(plan_qt[thread_id]);
                  fftw_execute(plan_p[thread_id]);

                  /* store FFT result in Q grids */
                  for (ifreq = 0; ifreq < nfreq; ++ifreq)
                    {
                      size_t grid_idx = CIDX4(t, T, ifreq, nfreq, i, w->nr, lmidx, w->nlm);
                      complex double qval = ((complex double) fft_q[thread_id][ifreq]);
                      complex double qtval = ((complex double) fft_qt[thread_id][ifreq]);
                      complex double pval = ((complex double) fft_p[thread_id][ifreq]);

                      Qq[grid_idx] = gsl_complex_rect(creal(qval), cimag(qval));
                      Qqt[grid_idx] = gsl_complex_rect(creal(qtval), cimag(qtval));
                      Qp[grid_idx] = gsl_complex_rect(creal(pval), cimag(pval));
                    }

                  /* advance window */
                  start_idx += nforward;
                }
            }
        }
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "do_transforms: writing FFT grids to %s...", PCA3D_STAGE2B_FFT_DATA);
  pca3d_write_fft_data2(PCA3D_STAGE2B_FFT_DATA, &fft_data, 0);
  fprintf(stderr, "done\n");

  fprintf(stderr, "do_transforms: writing FFT metadata to %s...", PCA3D_STAGE2B_FFT_DATA_LIGHT);
  pca3d_write_fft_data2(PCA3D_STAGE2B_FFT_DATA_LIGHT, &fft_data, 1);
  fprintf(stderr, "done\n");

  for (i = 0; i < (size_t) max_threads; ++i)
    {
      gsl_vector_complex_free(workq[i]);
      gsl_vector_complex_free(workqt[i]);
      gsl_vector_complex_free(workp[i]);

      fftw_free(fft_q[i]);
      fftw_free(fft_qt[i]);
      fftw_free(fft_p[i]);

      fftw_destroy_plan(plan_q[i]);
      fftw_destroy_plan(plan_qt[i]);
      fftw_destroy_plan(plan_p[i]);
    }

  gsl_vector_free(window);
  free(workq);
  free(workqt);
  free(workp);
  free(fft_q);
  free(fft_qt);
  free(fft_p);
  free(plan_q);
  free(plan_qt);
  free(plan_p);
  free(Qq);
  free(Qqt);
  free(Qp);

  return s;
}

/*
main_proc()

Inputs: data         - pca3d data
        prefix       - input file prefix
        fs           - sampling frequency (samples/day)
        window_size  - number of days in each time window
        window_shift - how many days to advance window

Return: success/error
*/

int
main_proc(const pca3d_data * data, const char * prefix, const double fs, const double window_size, const double window_shift)
{
  int s = 0;
  const size_t nt = data->nt;
  char buf[2048];
  size_t i;
  struct timeval tv0, tv1;
  magfield_workspace *w = data->w;
  gsl_complex *grid_q, *grid_qt, *grid_p;

  grid_q = malloc(nt * w->nr * w->nlm * sizeof(gsl_complex));
  grid_p = malloc(nt * w->nr * w->nlm * sizeof(gsl_complex));
  grid_qt = malloc(nt * w->nr * w->nlm * sizeof(gsl_complex));

  for (i = 0; i < nt; ++i)
    {
      sprintf(buf, "%s_%03zu.dat", prefix, i + 1);
      fprintf(stderr, "main_proc: reading %s...", buf);
      gettimeofday(&tv0, NULL);

      /* read previously calculated SH coefficients for this timestep */
      magfield_read_SH(buf, w);

      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

      /* copy SH coefficients for this timestep into (t,r,l,m) grids */
      copy_sh_grids(i, grid_q, grid_qt, grid_p, w);
    }

  /* compute windowed FFTs of q_l^m(t,r), qt_l^m(t,r), p_l^m(t,r) */
  do_transforms(data, fs, window_size, window_shift, grid_q, grid_qt, grid_p, w);

  magfield_free(w);
  free(grid_q);
  free(grid_qt);
  free(grid_p);

  return s;
}

int
main(int argc, char *argv[])
{
  char *input_prefix = PCA3D_STAGE1B_SH_PREFIX;
  const double fs = 24.0;    /* sample frequency (samples/day) */
  double window_size = 8.0;  /* number of days in each time segment */
  double window_shift = 4.0; /* number of days to shift forward in time */
  pca3d_data data;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "t:s:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 't':
            window_size = atof(optarg);
            break;

          case 's':
            window_shift = atof(optarg);
            break;

          default:
            fprintf(stderr, "Usage: %s [-t window_size (days)] [-s window_shift (days)]\n", argv[0]);
            exit(1);
            break;
        }
    }

  fprintf(stderr, "main: allocating magfield workspace from %s...", PCA3D_STAGE1B_DATA);
  data = pca3d_read_data(PCA3D_STAGE1B_DATA);
  fprintf(stderr, "done (%zu timestamps)\n", data.nt);

  main_proc(&data, input_prefix, fs, window_size, window_shift);

  return 0;
}
