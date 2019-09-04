/*
 * stage2a.c
 *
 * Perform PCA (frequency domain) directly on the 3D current grids from TIEGCM
 *
 * 1. Read 3D current grids for some time interval
 * 2. Divide total time interval into T smaller segments; perform windowed FFT
 *    of each time series segment for all grid points
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <sys/time.h>
#include <assert.h>
#include <omp.h>

#include <fftw3.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <mainlib/ml_common.h>

#include "io.h"
#include "pca3d.h"
#include "tiegcm3d.h"
#include "window.h"

#include "common.c"

/*
do_transforms()
  Perform FFTs on each time window segment for all grid points in TIEGCM grid

Inputs: fs           - sampling frequency (samples/day)
        window_size  - number of days in each time window
        window_shift - how many days to advance window
        data         - TIEGCM data
*/

static int
do_transforms(const double fs, const double window_size, const double window_shift,
              const tiegcm3d_data * data)
{
  int s = 0;
  const size_t nt = data->nt;
  const size_t nwindow = (size_t) (window_size * fs);                /* optimal number of samples per window */
  const size_t nforward = (size_t) (window_shift * fs);              /* number of samples to slide forward */
  const size_t nfreq = nwindow / 2 + 1;                              /* number of frequencies computed by FFT */
  const size_t T = count_windows(nt, fs, window_size, window_shift); /* number of time windows */
  const int max_threads = omp_get_max_threads();
  gsl_vector *window = gsl_vector_alloc(nwindow);                    /* window function */

  /* thread specific variables */
  gsl_vector **workr = malloc(max_threads * sizeof(gsl_vector *));
  gsl_vector **workt = malloc(max_threads * sizeof(gsl_vector *));
  gsl_vector **workp = malloc(max_threads * sizeof(gsl_vector *));
  fftw_complex **fft_r = malloc(max_threads * sizeof(fftw_complex *));
  fftw_complex **fft_t = malloc(max_threads * sizeof(fftw_complex *));
  fftw_complex **fft_p = malloc(max_threads * sizeof(fftw_complex *));
  fftw_plan *plan_r = malloc(max_threads * sizeof(fftw_plan));
  fftw_plan *plan_t = malloc(max_threads * sizeof(fftw_plan));
  fftw_plan *plan_p = malloc(max_threads * sizeof(fftw_plan));

  struct timeval tv0, tv1;
  gsl_complex *Qr; /* FT'd Jr grid, nfreq-by-nr-by-nlat-by-nlon */
  gsl_complex *Qt; /* FT'd Jt grid, nfreq-by-nr-by-nlat-by-nlon */
  gsl_complex *Qp; /* FT'd Jp grid, nfreq-by-nr-by-nlat-by-nlon */
  tiegcm3d_fft_data fft_data;
  size_t ir;
  int i;

  /* allocate thread variables */
  for (i = 0; i < max_threads; ++i)
    {
      workr[i] = gsl_vector_alloc(nwindow);
      workt[i] = gsl_vector_alloc(nwindow);
      workp[i] = gsl_vector_alloc(nwindow);

      fft_r[i] = fftw_malloc(sizeof(fftw_complex) * nfreq);
      fft_t[i] = fftw_malloc(sizeof(fftw_complex) * nfreq);
      fft_p[i] = fftw_malloc(sizeof(fftw_complex) * nfreq);

      plan_r[i] = fftw_plan_dft_r2c_1d(nwindow, workr[i]->data, fft_r[i], FFTW_ESTIMATE);
      plan_t[i] = fftw_plan_dft_r2c_1d(nwindow, workt[i]->data, fft_t[i], FFTW_ESTIMATE);
      plan_p[i] = fftw_plan_dft_r2c_1d(nwindow, workp[i]->data, fft_p[i], FFTW_ESTIMATE);
    }

  Qr = malloc(T * nfreq * data->nr * data->nlat * data->nlon * sizeof(gsl_complex));
  Qt = malloc(T * nfreq * data->nr * data->nlat * data->nlon * sizeof(gsl_complex));
  Qp = malloc(T * nfreq * data->nr * data->nlat * data->nlon * sizeof(gsl_complex));

  fprintf(stderr, "do_transforms: samples per window   = %zu\n", nwindow);
  fprintf(stderr, "do_transforms: sample slide forward = %zu\n", nforward);
  fprintf(stderr, "do_transforms: number of freqs      = %zu\n", nfreq);
  fprintf(stderr, "do_transforms: time samples         = %zu [hourly]\n", nt);
  fprintf(stderr, "do_transforms: window segments (T)  = %zu\n", T);

  /* compute window function */
  apply_ps1(NULL, window);
  /*apply_hamming(NULL, window);*/

  fft_data.nt = data->nt;
  fft_data.nfreq = nfreq;
  fft_data.nr = data->nr;
  fft_data.nlat = data->nlat;
  fft_data.nlon = data->nlon;
  fft_data.T = T;
  fft_data.fs = fs;
  fft_data.window_size = window_size;
  fft_data.window_shift = window_shift;
  fft_data.nwindow = nwindow;
  fft_data.t = data->t;
  fft_data.r = data->r;
  fft_data.glat = data->glat;
  fft_data.glon = data->glon;
  fft_data.window = window;
  fft_data.Jr = data->Jr;
  fft_data.Jt = data->Jt;
  fft_data.Jp = data->Jp;
  fft_data.Qr = Qr;
  fft_data.Qt = Qt;
  fft_data.Qp = Qp;

  fprintf(stderr, "do_transforms: computing FFTs of windowed data...");
  gettimeofday(&tv0, NULL);

#pragma omp parallel for private(ir)
  for (ir = 0; ir < data->nr; ++ir)
    {
      int thread_id = omp_get_thread_num();
      size_t ilat, ilon;

      for (ilat = 0; ilat < data->nlat; ++ilat)
        {
          for (ilon = 0; ilon < data->nlon; ++ilon)
            {
              size_t start_idx = 0; /* starting time index */
              size_t t;

              /* loop over window segments for this (r,theta,phi) point */
              for (t = 0; t < T; ++t)
                {
                  size_t end_idx = GSL_MIN(start_idx + nwindow - 1, nt - 1);
                  size_t n = end_idx - start_idx + 1; /* size of actual window */
                  double sqrtn = sqrt((double) n);
                  size_t it, ifreq;

                  assert(start_idx < end_idx);

                  if (n < nwindow)
                    {
                      /* could happen at the end of the time series; zero pad input buffer */
                      gsl_vector_set_zero(workr[thread_id]);
                      gsl_vector_set_zero(workt[thread_id]);
                      gsl_vector_set_zero(workp[thread_id]);
                    }

                  /* copy current time window into work arrays */
                  for (it = start_idx; it <= end_idx; ++it)
                    {
                      size_t idx = TIEGCM3D_IDX(it, ir, ilat, ilon, data);

                      gsl_vector_set(workr[thread_id], it - start_idx, data->Jr[idx]);
                      gsl_vector_set(workt[thread_id], it - start_idx, data->Jt[idx]);
                      gsl_vector_set(workp[thread_id], it - start_idx, data->Jp[idx]);
                    }

                  /* apply window function */
                  gsl_vector_mul(workr[thread_id], window);
                  gsl_vector_mul(workt[thread_id], window);
                  gsl_vector_mul(workp[thread_id], window);

                  /* compute FFT of this windowed data */
                  fftw_execute(plan_r[thread_id]);
                  fftw_execute(plan_t[thread_id]);
                  fftw_execute(plan_p[thread_id]);

                  /* store FFT result in Q grids */
                  for (ifreq = 0; ifreq < nfreq; ++ifreq)
                    {
                      size_t idx = TIEGCM3D_FREQIDX(t, ifreq, ir, ilat, ilon, data, T, nfreq);

                      Qr[idx] = gsl_complex_rect(fft_r[thread_id][ifreq][0] / sqrtn, fft_r[thread_id][ifreq][1] / sqrtn);
                      Qt[idx] = gsl_complex_rect(fft_t[thread_id][ifreq][0] / sqrtn, fft_t[thread_id][ifreq][1] / sqrtn);
                      Qp[idx] = gsl_complex_rect(fft_p[thread_id][ifreq][0] / sqrtn, fft_p[thread_id][ifreq][1] / sqrtn);
                    }

                  start_idx += nforward;
                }
            }
        }
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "do_transforms: writing FFT grids to %s...", PCA3D_STAGE2A_FFT_DATA);
  pca3d_write_fft_data(PCA3D_STAGE2A_FFT_DATA, &fft_data, 0);
  fprintf(stderr, "done\n");

  fprintf(stderr, "do_transforms: writing FFT metadata to %s...", PCA3D_STAGE2A_FFT_DATA_LIGHT);
  pca3d_write_fft_data(PCA3D_STAGE2A_FFT_DATA_LIGHT, &fft_data, 1);
  fprintf(stderr, "done\n");

  for (i = 0; i < max_threads; ++i)
    {
      gsl_vector_free(workr[i]);
      gsl_vector_free(workt[i]);
      gsl_vector_free(workp[i]);

      fftw_free(fft_r[i]);
      fftw_free(fft_t[i]);
      fftw_free(fft_p[i]);

      fftw_destroy_plan(plan_r[i]);
      fftw_destroy_plan(plan_t[i]);
      fftw_destroy_plan(plan_p[i]);
    }

  gsl_vector_free(window);
  free(Qr);
  free(Qt);
  free(Qp);
  free(workr);
  free(workt);
  free(workp);
  free(fft_r);
  free(fft_t);
  free(fft_p);
  free(plan_r);
  free(plan_t);
  free(plan_p);

  fftw_cleanup();

  return s;
}

int
main(int argc, char *argv[])
{
  const double fs = 24.0;    /* sample frequency (samples/day) */
  char *infile = NULL;
  double window_size = 2.0;  /* number of days in each time segment */
  double window_shift = 1.0; /* number of days to shift forward in time */
  struct timeval tv0, tv1;
  tiegcm3d_data *data;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "i:t:s:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          case 't':
            window_size = atof(optarg);
            break;

          case 's':
            window_shift = atof(optarg);
            break;

          default:
            fprintf(stderr, "Usage: %s <-i tiegcm3d_nc_file> [-t window_size (days)] [-s window_shift (days)]\n", argv[0]);
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i tiegcm3d_nc_file> [-t window_size (days)] [-s window_shift (days)]\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "main: input file          = %s\n", infile);
  fprintf(stderr, "main: sample frequency    = %g [samples/day]\n", fs);
  fprintf(stderr, "main: window size         = %g [days]\n", window_size);
  fprintf(stderr, "main: window shift        = %g [days]\n", window_shift);

  fprintf(stderr, "main: reading %s...", infile);
  gettimeofday(&tv0, NULL);

  data = tiegcm3d_read(infile, NULL);
  if (!data)
    {
      fprintf(stderr, "main: error reading %s\n", infile);
      exit(1);
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu records read, %g seconds)\n", data->nt,
          time_diff(tv0, tv1));

  do_transforms(fs, window_size, window_shift, data);

  return 0;
}
