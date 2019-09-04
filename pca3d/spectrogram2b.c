/*
 * spectrogram2b.c
 *
 * Print a spectrogram of q~_l^m(t,r), p_l^m(t,r), q_l^m(t,r) from TIEGCM 3D run
 *
 * 1. Read q~,p,q time series for given (l,m,r)
 * 2. Divide total time interval into T smaller segments; perform windowed FFT
 *    of each time series segment for all grid points
 * 3. Print frequency spectra for each window
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <sys/time.h>
#include <assert.h>
#include <errno.h>
#include <string.h>

#include <fftw3.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <mainlib/ml_bsearch.h>
#include <mainlib/ml_common.h>

#include "io.h"
#include "pca3d.h"
#include "tiegcm3d.h"
#include "window.h"

#include "common.c"

/*
spectrogram()
  Perform FFTs on each time window segment for all grid points in TIEGCM grid

Inputs: filename     - output data file for spectrogram
        fs           - sampling frequency (samples/day)
        window_size  - number of days in each time window
        window_shift - how many days to advance window
        data         - TIEGCM data
*/

static int
spectrogram(const char *filename, const size_t ir, const int l, const int m,
            const double fs, const double window_size, const double window_shift,
            const pca3d_fft_data * data)
{
  int s = 0;
  const size_t lmidx = magfield_lmidx(l, m, data->mmax);
  const size_t nt = data->nt;
  const size_t nwindow = (size_t) (window_size * fs);                /* optimal number of samples per window */
  const size_t nforward = (size_t) (window_shift * fs);              /* number of samples to slide forward */
  const size_t nfreq = nwindow;                                      /* number of frequencies computed by FFT */
  const size_t T = count_windows(nt, fs, window_size, window_shift); /* number of time windows */
  gsl_vector *window = gsl_vector_alloc(nwindow);                    /* window function */
  gsl_vector_complex *workq = gsl_vector_complex_alloc(nwindow);
  gsl_vector_complex *workqt = gsl_vector_complex_alloc(nwindow);
  gsl_vector_complex *workp = gsl_vector_complex_alloc(nwindow);
  fftw_complex *fft_q = fftw_malloc(sizeof(fftw_complex) * nfreq);
  fftw_complex *fft_qt = fftw_malloc(sizeof(fftw_complex) * nfreq);
  fftw_complex *fft_p = fftw_malloc(sizeof(fftw_complex) * nfreq);
  fftw_plan plan_q, plan_qt, plan_p;
  struct timeval tv0, tv1;
  size_t start_idx = 0; /* starting time index */
  size_t t;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "spectrogram: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  plan_q = fftw_plan_dft_1d(nwindow, (fftw_complex *) workq->data, fft_q, FFTW_FORWARD, FFTW_ESTIMATE);
  plan_qt = fftw_plan_dft_1d(nwindow, (fftw_complex *) workqt->data, fft_qt, FFTW_FORWARD, FFTW_ESTIMATE);
  plan_p = fftw_plan_dft_1d(nwindow, (fftw_complex *) workp->data, fft_p, FFTW_FORWARD, FFTW_ESTIMATE);

  fprintf(stderr, "spectrogram: samples per window   = %zu\n", nwindow);
  fprintf(stderr, "spectrogram: sample slide forward = %zu\n", nforward);
  fprintf(stderr, "spectrogram: number of freqs      = %zu\n", nfreq);
  fprintf(stderr, "spectrogram: time samples         = %zu [hourly]\n", nt);
  fprintf(stderr, "spectrogram: window segments (T)  = %zu\n", T);

  /* compute window function */
  apply_ps1(NULL, window);
  /*apply_hamming(NULL, window);*/

  t = 1;
  fprintf(fp, "# SH degree:     %d\n", l);
  fprintf(fp, "# SH order:      %d\n", m);
  fprintf(fp, "# Radius:        %.2f (km) [%.2f km altitude]\n", data->r[ir], data->r[ir] - R_EARTH_KM);
  fprintf(fp, "# window size:   %g [days]\n", window_size);
  fprintf(fp, "# window shift:  %g [days]\n", window_shift);
  fprintf(fp, "# Field %zu: timestamp of window center (UT seconds since 1970-01-01 00:00:00 UTC)\n", t++);
  fprintf(fp, "# Field %zu: frequency (days^{-1})\n", t++);
  fprintf(fp, "# Field %zu: period (days)\n", t++);
  fprintf(fp, "# Field %zu: Power/Frequency of q~_l^m(t,r) (dB/days)\n", t++);
  fprintf(fp, "# Field %zu: Power/Frequency of p_l^m(t,r) (dB/days)\n", t++);
  fprintf(fp, "# Field %zu: Power/Frequency of q_l^m(t,r) (dB/days)\n", t++);
  fprintf(fp, "# Field %zu: Power in q~_l^m(t,r) (nT^2)\n", t++);
  fprintf(fp, "# Field %zu: Power in p_l^m(t,r) (nT^2)\n", t++);
  fprintf(fp, "# Field %zu: Power in q_l^m(t,r) (nT^2)\n", t++);

  fprintf(stderr, "spectrogram: computing FFTs of windowed data...");
  gettimeofday(&tv0, NULL);

  for (t = 0; t < T; ++t)
    {
      size_t end_idx = GSL_MIN(start_idx + nwindow - 1, nt - 1);
      size_t n = end_idx - start_idx + 1; /* size of actual window */
      size_t it, ifreq;

      assert(start_idx < end_idx);

      if (n < nwindow)
        {
          /* could happen at the end of the time series; zero pad input buffer */
          gsl_vector_complex_set_zero(workq);
          gsl_vector_complex_set_zero(workqt);
          gsl_vector_complex_set_zero(workp);
        }

      /* copy current time window into work arrays */
      for (it = start_idx; it <= end_idx; ++it)
        {
          size_t grid_idx = CIDX3(it, nt, ir, data->nr, lmidx, data->nlm);

          gsl_vector_complex_set(workq, it - start_idx, data->grid_q[grid_idx]);
          gsl_vector_complex_set(workqt, it - start_idx, data->grid_qt[grid_idx]);
          gsl_vector_complex_set(workp, it - start_idx, data->grid_p[grid_idx]);
        }

      /* apply window function */
      {
        size_t jj;

        for (jj = 0; jj < n; ++jj)
          {
            gsl_complex *qptr = gsl_vector_complex_ptr(workq, jj);
            gsl_complex *qtptr = gsl_vector_complex_ptr(workqt, jj);
            gsl_complex *pptr = gsl_vector_complex_ptr(workp, jj);
            double wj = gsl_vector_get(window, jj);

            *qptr = gsl_complex_mul_real(*qptr, wj);
            *qtptr = gsl_complex_mul_real(*qtptr, wj);
            *pptr = gsl_complex_mul_real(*pptr, wj);
          }
      }

      /* compute FFT of this windowed data */
      fftw_execute(plan_q);
      fftw_execute(plan_qt);
      fftw_execute(plan_p);

      /* print FFT result for this time window */
      for (ifreq = 0; ifreq < nfreq; ++ifreq)
        {
          gsl_complex Qq = gsl_complex_mul_real(gsl_complex_rect(fft_q[ifreq][0], fft_q[ifreq][1]), 1.0e9);
          gsl_complex Qqt = gsl_complex_mul_real(gsl_complex_rect(fft_qt[ifreq][0], fft_qt[ifreq][1]), 1.0e9);
          gsl_complex Qp = gsl_complex_mul_real(gsl_complex_rect(fft_p[ifreq][0], fft_p[ifreq][1]), 1.0e9);
          double freq = ifreq * (fs / nwindow);
          double period = 1.0 / freq;

          fprintf(fp, "%ld %f %f %.12e %.12e %.12e %.12e %.12e %.12e\n",
                  data->t[start_idx + n/2],
                  freq,
                  period,
                  10.0 * log10(gsl_complex_abs2(Qqt) / (fs * nwindow)),
                  10.0 * log10(gsl_complex_abs2(Qp) / (fs * nwindow)),
                  10.0 * log10(gsl_complex_abs2(Qq) / (fs * nwindow)),
                  gsl_complex_abs2(Qqt),
                  gsl_complex_abs2(Qp),
                  gsl_complex_abs2(Qq));
        }

      fprintf(fp, "\n");

      start_idx += nforward;
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds, output file is %s)\n", time_diff(tv0, tv1), filename);

  fclose(fp);

  gsl_vector_free(window);
  gsl_vector_complex_free(workq);
  gsl_vector_complex_free(workqt);
  gsl_vector_complex_free(workp);
  fftw_free(fft_q);
  fftw_free(fft_qt);
  fftw_free(fft_p);
  fftw_destroy_plan(plan_q);
  fftw_destroy_plan(plan_qt);
  fftw_destroy_plan(plan_p);

  fftw_cleanup();

  return s;
}

/* print time series of q~_l^m(t,r), p_l^m(t,r), q_l^m(t,r) */
int
print_time(const char * filename, const size_t ir, const int l, const int m, const pca3d_fft_data * data)
{
  const size_t lmidx = magfield_lmidx(l, m, data->mmax);
  FILE *fp;
  size_t i;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_time: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# SH degree: %d\n", l);
  fprintf(fp, "# SH order:  %d\n", m);
  fprintf(fp, "# Radius: %.2f (km) [%.2f km altitude]\n", data->r[ir], data->r[ir] - R_EARTH_KM);
  fprintf(fp, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
  fprintf(fp, "# Field %zu: Re q~_l^m(t,r) (nT)\n", i++);
  fprintf(fp, "# Field %zu: Im q~_l^m(t,r) (nT)\n", i++);
  fprintf(fp, "# Field %zu: Re p_l^m(t,r) (nT)\n", i++);
  fprintf(fp, "# Field %zu: Im p_l^m(t,r) (nT)\n", i++);
  fprintf(fp, "# Field %zu: Re q_l^m(t,r) (nT)\n", i++);
  fprintf(fp, "# Field %zu: Im q_l^m(t,r) (nT)\n", i++);

  for (i = 0; i < data->nt; ++i)
    {
      size_t idx = CIDX3(i, data->nt, ir, data->nr, lmidx, data->nlm);

      fprintf(fp, "%ld %16.4e %16.4e %16.4e %16.4e %16.4e %16.4e\n",
              data->t[i],
              GSL_REAL(data->grid_qt[idx]) * 1.0e9,
              GSL_IMAG(data->grid_qt[idx]) * 1.0e9,
              GSL_REAL(data->grid_p[idx]) * 1.0e9,
              GSL_IMAG(data->grid_p[idx]) * 1.0e9,
              GSL_REAL(data->grid_q[idx]) * 1.0e9,
              GSL_IMAG(data->grid_q[idx]) * 1.0e9);
    }

  fclose(fp);

  return 0;
}

int
main(int argc, char *argv[])
{
  const double fs = 24.0;    /* sample frequency (samples/day) */
  char *infile = PCA3D_STAGE2B_FFT_DATA;
  char *outfile_time = "plots/spectrogram2b_time.txt";
  char *outfile = "plots/spectrogram2b.txt";
  double window_size = 2.0;   /* number of days in each time segment */
  double window_shift = 0.25; /* number of days to shift forward in time */
  struct timeval tv0, tv1;
  pca3d_fft_data data;
  int l = 1;                     /* desired SH degree */
  int m = 0;                     /* desired SH order */
  double r = R_EARTH_KM + 110.0; /* desired radius */
  int r_idx;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "altitude", required_argument, NULL, 'a' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "a:i:l:m:o:t:s:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'a':
            r = R_EARTH_KM + atof(optarg);
            break;

          case 'i':
            infile = optarg;
            break;

          case 'o':
            outfile = optarg;
            break;

          case 'l':
            l = atoi(optarg);
            break;

          case 'm':
            m = atoi(optarg);
            break;

          case 't':
            window_size = atof(optarg);
            break;

          case 's':
            window_shift = atof(optarg);
            break;

          default:
            fprintf(stderr, "Usage: %s <-i stage2b_fft_file> [-a altitude (km)] [-l sh_degree] [-m sh_order] [-t window_size (days)] [-s window_shift (days)] [-o output_file]\n", argv[0]);
            exit(1);
            break;
        }
    }

  fprintf(stderr, "main: input file          = %s\n", infile);
  fprintf(stderr, "main: output file         = %s\n", outfile);
  fprintf(stderr, "main: sample frequency    = %g [samples/day]\n", fs);
  fprintf(stderr, "main: window size         = %g [days]\n", window_size);
  fprintf(stderr, "main: window shift        = %g [days]\n", window_shift);

  fprintf(stderr, "main: reading %s...", infile);
  gettimeofday(&tv0, NULL);
  data = pca3d_read_fft_data2(infile);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  /* locate index of desired altitude */
  r_idx = bsearch_double(data.r, r, 0, data.nr - 1);

  fprintf(stderr, "main: r_idx     = %d (%.2f [km])\n", r_idx, data.r[r_idx] - R_EARTH_KM);
  fprintf(stderr, "main: SH degree = %d\n", l);
  fprintf(stderr, "main: SH order  = %d\n", m);

  fprintf(stderr, "main: printing time series to %s...", outfile_time);
  print_time(outfile_time, r_idx, l, m, &data);
  fprintf(stderr, "done\n");

  spectrogram(outfile, r_idx, l, m, fs, window_size, window_shift, &data);

  return 0;
}
