/*
 * print_fft2b.c
 *
 * Print results of FFT analysis (stage2b) on SH time series
 *
 * ./print_fft2b [-i fft_data_file]
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
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <mainlib/ml_common.h>
#include <mainlib/ml_bsearch.h>

#include "io.h"
#include "pca3d.h"

/*
print_windows()
  Print data files (for each time window segment) of:

1. original SH(l,m) time series
2. windowed SH(l,m) time series
3. PSD of SH(l,m)
*/

int
print_windows(const size_t ir, const size_t l, const int m, const pca3d_fft_data *data)
{
  int s = 0;
  const char *file_time = "plots/window_time_SH.txt";
  const char *file_freq = "plots/window_freq_SH.txt";
  const size_t lmidx = magfield_lmidx(l, m, data->mmax);
  const size_t T = data->T;   /* number of time window segments */
  const size_t nt = data->nt; /* number of time steps in SH time series */
  const double fs = data->fs;
  const size_t nwindow = (size_t) (data->window_size * fs);   /* optimal number of samples per window */
  const size_t nforward = (size_t) (data->window_shift * fs); /* number of samples to slide forward */
  const size_t nfreq = nwindow;                               /* FFT output buffer size (number of frequencies) */
  size_t start_idx = 0;       /* starting time index */
  size_t t;
  FILE *fp_t, *fp_f;

  fp_t = fopen(file_time, "w");
  fp_f = fopen(file_freq, "w");

  t = 1;
  fprintf(fp_t, "# SH degree:    %zu\n", l);
  fprintf(fp_t, "# SH order:     %d\n", m);
  fprintf(fp_t, "# Radius:       %.2f (km) [%.2f km altitude]\n", data->r[ir], data->r[ir] - R_EARTH_KM);
  fprintf(fp_t, "# window size:  %g [days]\n", data->window_size);
  fprintf(fp_t, "# window shift: %g [days]\n", data->window_shift);
  fprintf(fp_t, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", t++);
  fprintf(fp_t, "# Field %zu: Re qtlm (nT)\n", t++);
  fprintf(fp_t, "# Field %zu: Im qtlm (nT)\n", t++);
  fprintf(fp_t, "# Field %zu: Re plm (nT)\n", t++);
  fprintf(fp_t, "# Field %zu: Im plm (nT)\n", t++);
  fprintf(fp_t, "# Field %zu: Re qlm (nT)\n", t++);
  fprintf(fp_t, "# Field %zu: Im qlm (nT)\n", t++);
  fprintf(fp_t, "# Field %zu: windowed Re qtlm (nT)\n", t++);
  fprintf(fp_t, "# Field %zu: windowed Im qtlm (nT)\n", t++);
  fprintf(fp_t, "# Field %zu: windowed Re plm (nT)\n", t++);
  fprintf(fp_t, "# Field %zu: windowed Im plm (nT)\n", t++);
  fprintf(fp_t, "# Field %zu: windowed Re qlm (nT)\n", t++);
  fprintf(fp_t, "# Field %zu: windowed Im qlm (nT)\n", t++);

  t = 1;
  fprintf(fp_f, "# SH degree:    %zu\n", l);
  fprintf(fp_f, "# SH order:     %d\n", m);
  fprintf(fp_f, "# Radius:       %.2f (km) [%.2f km altitude]\n", data->r[ir], data->r[ir] - R_EARTH_KM);
  fprintf(fp_f, "# window size:  %g [days]\n", data->window_size);
  fprintf(fp_f, "# window shift: %g [days]\n", data->window_shift);
  fprintf(fp_f, "# Field %zu: frequency (days^{-1})\n", t++);
  fprintf(fp_f, "# Field %zu: period (days)\n", t++);
  fprintf(fp_f, "# Field %zu: Power/Frequency of q~_l^m(t, r) (dB/days)\n", t++);
  fprintf(fp_f, "# Field %zu: Power/Frequency of p_l^m(t, r) (dB/days)\n", t++);
  fprintf(fp_f, "# Field %zu: Power/Frequency of q_l^m(t, r) (dB/days)\n", t++);

  for (t = 0; t < T; ++t)
    {
      size_t end_idx = GSL_MIN(start_idx + nwindow - 1, nt - 1);
      size_t it, ifreq;

      assert(start_idx < end_idx);

      /* print original and windowed data for this time window */
      for (it = start_idx; it <= end_idx; ++it)
        {
          size_t idx = CIDX3(it, nt, ir, data->nr, lmidx, data->nlm);
          double wi = gsl_vector_get(data->window, it - start_idx);

          fprintf(fp_t, "%ld %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n",
                 data->t[it],
                 GSL_REAL(data->grid_qt[idx]) * 1.0e9,
                 GSL_IMAG(data->grid_qt[idx]) * 1.0e9,
                 GSL_REAL(data->grid_p[idx]) * 1.0e9,
                 GSL_IMAG(data->grid_p[idx]) * 1.0e9,
                 GSL_REAL(data->grid_q[idx]) * 1.0e9,
                 GSL_IMAG(data->grid_q[idx]) * 1.0e9,
                 GSL_REAL(data->grid_qt[idx]) * wi * 1.0e9,
                 GSL_IMAG(data->grid_qt[idx]) * wi * 1.0e9,
                 GSL_REAL(data->grid_p[idx]) * wi * 1.0e9,
                 GSL_IMAG(data->grid_p[idx]) * wi * 1.0e9,
                 GSL_REAL(data->grid_q[idx]) * wi * 1.0e9,
                 GSL_IMAG(data->grid_q[idx]) * wi * 1.0e9);
        }

      for (ifreq = 0; ifreq < nfreq; ++ifreq)
        {
          size_t idx = CIDX4(t, T, ifreq, nfreq, ir, data->nr, lmidx, data->nlm);
          gsl_complex Qq = gsl_complex_mul_real(data->Qq[idx], 1.0e9);
          gsl_complex Qqt = gsl_complex_mul_real(data->Qqt[idx], 1.0e9);
          gsl_complex Qp = gsl_complex_mul_real(data->Qp[idx], 1.0e9);
          double freq = ifreq * (fs / nwindow);
          double period = 1.0 / freq;

          fprintf(fp_f, "%f %f %.12e %.12e %.12e\n",
                  freq,
                  period,
                  10.0 * log10(gsl_complex_abs2(Qqt) / (fs * nwindow)),
                  10.0 * log10(gsl_complex_abs2(Qp) / (fs * nwindow)),
                  10.0 * log10(gsl_complex_abs2(Qq) / (fs * nwindow)));
        }

      if (t != T - 1)
        {
          fprintf(fp_t, "\n\n");
          fprintf(fp_f, "\n\n");
        }

      start_idx += nforward;
    }

  fclose(fp_t);
  fclose(fp_f);

  fprintf(stderr, "print_windows: wrote %s\n", file_time);
  fprintf(stderr, "print_windows: wrote %s\n", file_freq);

  return s;
}

int
main(int argc, char *argv[])
{
  pca3d_fft_data data;
  struct timeval tv0, tv1;
  char *infile = PCA3D_STAGE2B_FFT_DATA;
  int l = 1;          /* desired SH degree */
  int m = 0;          /* desired SH order */
  double alt = 110.0; /* desired altitude (km) */
  int r_idx;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "a:i:l:m:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          case 'a':
            alt = atof(optarg);
            break;

          case 'l':
            l = atoi(optarg);
            break;

          case 'm':
            m = atoi(optarg);
            break;

          default:
            fprintf(stderr, "Usage: %s [-i fft_data_file] [-a altitude (km)] [-l degree_l] [-m order_m]\n", argv[0]);
            exit(1);
            break;
        }
    }

  fprintf(stderr, "main: reading %s...", infile);
  gettimeofday(&tv0, NULL);
  data = pca3d_read_fft_data2(infile);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  /* locate index of desired alt/lat/lon */
  r_idx = bsearch_double(data.r, alt + R_EARTH_KM, 0, data.nr - 1);

  fprintf(stderr, "main: r_idx     = %d (%.2f [km])\n", r_idx, data.r[r_idx] - R_EARTH_KM);
  fprintf(stderr, "main: SH degree = %d\n", l);
  fprintf(stderr, "main: SH order  = %d\n", m);

  print_windows(r_idx, (size_t) l, m, &data);

  free(data.t);
  free(data.r);
  free(data.Qq);
  free(data.Qqt);
  free(data.Qp);
  gsl_vector_free(data.window);

  return 0;
}
