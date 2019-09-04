/*
 * spectrogram.c
 *
 * Print a spectrogram of a given grid point time series from 3D TIEGCM run
 *
 * 1. Read 3D current grids for some time interval
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
spectrogram(const char *filename, const size_t ir, const size_t ilat, const size_t ilon,
            const double fs, const double window_size, const double window_shift,
            const tiegcm3d_data * data)
{
  int s = 0;
  const size_t nt = data->nt;
  const size_t nwindow = (size_t) (window_size * fs);                /* optimal number of samples per window */
  const size_t nforward = (size_t) (window_shift * fs);              /* number of samples to slide forward */
  const size_t nfreq = nwindow / 2 + 1;                              /* number of frequencies computed by FFT */
  const size_t T = count_windows(nt, fs, window_size, window_shift); /* number of time windows */
  gsl_vector *window = gsl_vector_alloc(nwindow);                    /* window function */
  gsl_vector *workr = gsl_vector_alloc(nwindow);
  gsl_vector *workt = gsl_vector_alloc(nwindow);
  gsl_vector *workp = gsl_vector_alloc(nwindow);
  fftw_complex *fft_r = fftw_malloc(sizeof(fftw_complex) * nfreq);
  fftw_complex *fft_t = fftw_malloc(sizeof(fftw_complex) * nfreq);
  fftw_complex *fft_p = fftw_malloc(sizeof(fftw_complex) * nfreq);
  fftw_plan plan_r, plan_t, plan_p;
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

  plan_r = fftw_plan_dft_r2c_1d(nwindow, workr->data, fft_r, FFTW_ESTIMATE);
  plan_t = fftw_plan_dft_r2c_1d(nwindow, workt->data, fft_t, FFTW_ESTIMATE);
  plan_p = fftw_plan_dft_r2c_1d(nwindow, workp->data, fft_p, FFTW_ESTIMATE);

  fprintf(stderr, "spectrogram: samples per window   = %zu\n", nwindow);
  fprintf(stderr, "spectrogram: sample slide forward = %zu\n", nforward);
  fprintf(stderr, "spectrogram: number of freqs      = %zu\n", nfreq);
  fprintf(stderr, "spectrogram: time samples         = %zu [hourly]\n", nt);
  fprintf(stderr, "spectrogram: window segments (T)  = %zu\n", T);

  /* compute window function */
  /*apply_ps1(NULL, window);*/
  apply_hamming(NULL, window);

  t = 1;
  fprintf(fp, "# Latitude: %.2f (deg)\n", data->glat[ilat]);
  fprintf(fp, "# Longitude: %.2f (deg)\n", data->glon[ilon]);
  fprintf(fp, "# Radius: %.2f (km) [%.2f km altitude]\n", data->r[ir], data->r[ir] - R_EARTH_KM);
  fprintf(fp, "# Field %zu: timestamp of window center (UT seconds since 1970-01-01 00:00:00 UTC)\n", t++);
  fprintf(fp, "# Field %zu: period (days)\n", t++);
  fprintf(fp, "# Field %zu: frequency (days^{-1})\n", t++);
  fprintf(fp, "# Field %zu: Power in J_r (uA/m^2)\n", t++);
  fprintf(fp, "# Field %zu: Power in J_t (uA/m^2)\n", t++);
  fprintf(fp, "# Field %zu: Power in J_p (uA/m^2)\n", t++);

  fprintf(stderr, "spectrogram: computing FFTs of windowed data...");
  gettimeofday(&tv0, NULL);

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
          gsl_vector_set_zero(workr);
          gsl_vector_set_zero(workt);
          gsl_vector_set_zero(workp);
        }

      /* copy current time window into work arrays */
      for (it = start_idx; it <= end_idx; ++it)
        {
          size_t idx = TIEGCM3D_IDX(it, ir, ilat, ilon, data);

          gsl_vector_set(workr, it - start_idx, data->Jr[idx]);
          gsl_vector_set(workt, it - start_idx, data->Jt[idx]);
          gsl_vector_set(workp, it - start_idx, data->Jp[idx]);
        }

      /* apply window function */
      gsl_vector_mul(workr, window);
      gsl_vector_mul(workt, window);
      gsl_vector_mul(workp, window);

      /* compute FFT of this windowed data */
      fftw_execute(plan_r);
      fftw_execute(plan_t);
      fftw_execute(plan_p);

      /* print FFT result for this time window */
      for (ifreq = 0; ifreq < nfreq; ++ifreq)
        {
          gsl_complex Qr = gsl_complex_rect(fft_r[ifreq][0] / sqrtn, fft_r[ifreq][1] / sqrtn);
          gsl_complex Qt = gsl_complex_rect(fft_t[ifreq][0] / sqrtn, fft_t[ifreq][1] / sqrtn);
          gsl_complex Qp = gsl_complex_rect(fft_p[ifreq][0] / sqrtn, fft_p[ifreq][1] / sqrtn);
          double freq = ifreq * (fs / n);
          double period = 1.0 / freq;
          double fac = (ifreq == 0) ? 1.0 : 2.0;

          fprintf(fp, "%ld %f %f %.12e %.12e %.12e\n",
                  data->t[start_idx + n/2],
                  period,
                  freq,
                  fac * gsl_complex_abs(Qr),
                  fac * gsl_complex_abs(Qt),
                  fac * gsl_complex_abs(Qp));
        }

      fprintf(fp, "\n");

      start_idx += nforward;
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds, output file is %s)\n", time_diff(tv0, tv1), filename);

  fclose(fp);

  gsl_vector_free(window);
  gsl_vector_free(workr);
  gsl_vector_free(workt);
  gsl_vector_free(workp);
  fftw_free(fft_r);
  fftw_free(fft_t);
  fftw_free(fft_p);
  fftw_destroy_plan(plan_r);
  fftw_destroy_plan(plan_t);
  fftw_destroy_plan(plan_p);

  fftw_cleanup();

  return s;
}

int
main(int argc, char *argv[])
{
  const double fs = 24.0;    /* sample frequency (samples/day) */
  char *infile = NULL;
  char *outfile_time = "plots/spectrogram_time.txt";
  char *outfile = "plots/spectrogram.txt";
  double window_size = 2.0;   /* number of days in each time segment */
  double window_shift = 0.25; /* number of days to shift forward in time */
  struct timeval tv0, tv1;
  tiegcm3d_data *data;
  double lon = 150.0;        /* desired longitude */
  double lat = 8.0;          /* desired latitude */
  double r = R_EARTH_KM + 110.0; /* desired radius */
  int r_idx, lat_idx, lon_idx;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "altitude", required_argument, NULL, 'a' },
          { "longitude", required_argument, NULL, 'b' },
          { "latitude", required_argument, NULL, 'c' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "a:b:c:i:o:t:s:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'a':
            r = R_EARTH_KM + atof(optarg);
            break;

          case 'b':
            lon = atof(optarg);
            break;

          case 'c':
            lat = atof(optarg);
            break;

          case 'i':
            infile = optarg;
            break;

          case 'o':
            outfile = optarg;
            break;

          case 't':
            window_size = atof(optarg);
            break;

          case 's':
            window_shift = atof(optarg);
            break;

          default:
            fprintf(stderr, "Usage: %s <-i tiegcm3d_nc_file> [-t window_size (days)] [-s window_shift (days)] [-o output_file]\n", argv[0]);
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i tiegcm3d_nc_file> [-t window_size (days)] [-s window_shift (days)] [-o output_file]\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "main: input file          = %s\n", infile);
  fprintf(stderr, "main: output file         = %s\n", outfile);
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

  /* locate index of desired grid location */
  r_idx = bsearch_double(data->r, r, 0, data->nr - 1);
  lat_idx = bsearch_double(data->glat, lat, 0, data->nlat - 1);
  lon_idx = bsearch_double(data->glon, lon, 0, data->nlon - 1);

  fprintf(stderr, "main: r_idx   = %d (%.2f [km] altitude)\n", r_idx, data->r[r_idx] - R_EARTH_KM);
  fprintf(stderr, "main: lat_idx = %d (%.2f [deg])\n", lat_idx, data->glat[lat_idx]);
  fprintf(stderr, "main: lon_idx = %d (%.2f [deg])\n", lon_idx, data->glon[lon_idx]);

  fprintf(stderr, "main: printing time series to %s...", outfile_time);
  tiegcm3d_print_time(outfile_time, data, r_idx, lat_idx, lon_idx);
  fprintf(stderr, "done\n");

  spectrogram(outfile, r_idx, lat_idx, lon_idx, fs, window_size, window_shift, data);

  return 0;
}
