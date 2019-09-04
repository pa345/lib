/*
 * print_fft.c
 *
 * Print results of FFT analysis (stage2a) on J grids
 *
 * ./print_fft [-i fft_data_file]
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
#include "tiegcm3d.h"

/*
print_windows()
  Print data files (for each time window segment) of:

1. original Jr/Jt/Jp time series
2. windowed Jr/Jt/Jp time series
3. PSD of Jr/Jt/Jp
*/

int
print_windows(const size_t ir, const size_t ilat, const size_t ilon, const tiegcm3d_fft_data *data)
{
  int s = 0;
  const char *file_time = "plots/window_time.txt";
  const char *file_freq = "plots/window_freq.txt";
  const size_t T = data->T;   /* number of time window segments */
  const size_t nt = data->nt; /* number of time steps in J grids */
  const double fs = data->fs;
  const size_t nwindow = (size_t) (data->window_size * fs);   /* optimal number of samples per window */
  const size_t nforward = (size_t) (data->window_shift * fs); /* number of samples to slide forward */
  const size_t nfreq = nwindow / 2 + 1;                       /* FFT output buffer size (number of frequencies) */
  size_t start_idx = 0;       /* starting time index */
  size_t t;
  FILE *fp_t, *fp_f;

  fp_t = fopen(file_time, "w");
  fp_f = fopen(file_freq, "w");

  t = 1;
  fprintf(fp_t, "# Latitude:  %.2f [deg]\n", data->glat[ilat]);
  fprintf(fp_t, "# Longitude: %.2f [deg]\n", data->glon[ilon]);
  fprintf(fp_t, "# Altitude:  %.2f [km]\n", data->r[ir] - R_EARTH_KM);
  fprintf(fp_t, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", t++);
  fprintf(fp_t, "# Field %zu: J_r (uA/m^2)\n", t++);
  fprintf(fp_t, "# Field %zu: J_t (uA/m^2)\n", t++);
  fprintf(fp_t, "# Field %zu: J_p (uA/m^2)\n", t++);
  fprintf(fp_t, "# Field %zu: windowed J_r (uA/m^2)\n", t++);
  fprintf(fp_t, "# Field %zu: windowed J_t (uA/m^2)\n", t++);
  fprintf(fp_t, "# Field %zu: windowed J_p (uA/m^2)\n", t++);

  t = 1;
  fprintf(fp_f, "# Field %zu: frequency (days^{-1})\n", t++);
  fprintf(fp_f, "# Field %zu: period (days)\n", t++);
  fprintf(fp_f, "# Field %zu: Power in J_r (uA/m^2)\n", t++);
  fprintf(fp_f, "# Field %zu: Power in J_t (uA/m^2)\n", t++);
  fprintf(fp_f, "# Field %zu: Power in J_p (uA/m^2)\n", t++);

  for (t = 0; t < T; ++t)
    {
      size_t end_idx = GSL_MIN(start_idx + nwindow - 1, nt - 1);
      size_t n = end_idx - start_idx + 1; /* size of actual window */
      size_t it, ifreq;

      assert(start_idx < end_idx);

      /* print original and windowed J data for this time window */
      for (it = start_idx; it <= end_idx; ++it)
        {
          size_t idx = TIEGCM3D_IDX(it, ir, ilat, ilon, data);
          double wi = gsl_vector_get(data->window, it - start_idx);

          fprintf(fp_t, "%ld %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n",
                 data->t[it],
                 data->Jr[idx],
                 data->Jt[idx],
                 data->Jp[idx],
                 data->Jr[idx] * wi,
                 data->Jt[idx] * wi,
                 data->Jp[idx] * wi);
        }

      for (ifreq = 0; ifreq < nfreq; ++ifreq)
        {
          size_t idx = TIEGCM3D_FREQIDX(t, ifreq, ir, ilat, ilon, data, T, nfreq);
          gsl_complex Qr = data->Qr[idx];
          gsl_complex Qt = data->Qt[idx];
          gsl_complex Qp = data->Qp[idx];
          double freq = ifreq * (fs / n);
          double period = 1.0 / freq;
          double fac = (ifreq == 0) ? 1.0 : 2.0;

          fprintf(fp_f, "%f %f %.12e %.12e %.12e\n",
                  freq,
                  period,
                  fac * gsl_complex_abs(Qr),
                  fac * gsl_complex_abs(Qt),
                  fac * gsl_complex_abs(Qp));
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
  tiegcm3d_fft_data data;
  struct timeval tv0, tv1;
  char *infile = PCA3D_STAGE2A_FFT_DATA;
  double lon = 150.0; /* desired longitude */
  double lat = 8.0;   /* desired latitude */
  double alt = 110.0; /* desired altitude */
  int r_idx, lat_idx, lon_idx;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "a:b:c:i:r:", long_options, &option_index);
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

          case 'b':
            lat = atof(optarg);
            break;

          case 'c':
            lon = atof(optarg);
            break;

          default:
            fprintf(stderr, "Usage: %s [-i fft_data_file] [-a altitude (km)] [-b latitude (deg)] [-c longitude (deg)] [-o output_file]\n", argv[0]);
            exit(1);
            break;
        }
    }

  fprintf(stderr, "main: reading %s...", infile);
  gettimeofday(&tv0, NULL);
  data = pca3d_read_fft_data(infile);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  /* locate index of desired alt/lat/lon */
  r_idx = bsearch_double(data.r, alt + R_EARTH_KM, 0, data.nr - 1);
  lat_idx = bsearch_double(data.glat, lat, 0, data.nlat - 1);
  lon_idx = bsearch_double(data.glon, lon, 0, data.nlon - 1);

  fprintf(stderr, "main: r_idx = %d (%.2f [km])\n", r_idx, data.r[r_idx] - R_EARTH_KM);
  fprintf(stderr, "main: lat_idx = %d (%.2f [deg])\n", lat_idx, data.glat[lat_idx]);
  fprintf(stderr, "main: lon_idx = %d (%.2f [deg])\n", lon_idx, data.glon[lon_idx]);

  print_windows(r_idx, lat_idx, lon_idx, &data);

  free(data.t);
  free(data.r);
  free(data.glat);
  free(data.glon);
  free(data.Jr);
  free(data.Jt);
  free(data.Jp);
  free(data.Qp);
  gsl_vector_free(data.window);

  return 0;
}
