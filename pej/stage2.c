/*
 * stage2.c
 *
 * Fit spherical harmonic model to previously selected dataset
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <complex.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <satdata/satdata.h>
#include <indices/indices.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rstat.h>

#include <common/bin2d.h>
#include <common/common.h>

#include "magdata.h"
#include "magfit.h"
#include "track.h"

int
print_data(const char *filename, const magdata *data)
{
  int s = 0;
  size_t i;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_data: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
  fprintf(fp, "# Field %zu: magnetic local time (degrees)\n", i++);
  fprintf(fp, "# Field %zu: QD latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: X NEC residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: Y NEC residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: Z NEC residual (nT)\n", i++);

  for (i = 0; i < data->n; ++i)
    {
      time_t unix_time = satdata_epoch2timet(data->t[i]);
      double MLT_deg = 15.0 * data->MLT[i] - 180.0;
      double B[4];

      if (!MAGDATA_ExistVector(data->flags[i]))
        continue;

      magdata_residual(i, B, data);

      fprintf(fp, "%ld %7.4f %7.4f %8.2f %8.2f %8.2f\n",
              unix_time,
              MLT_deg,
              data->qdlat[i],
              B[0],
              B[1],
              B[2]);
    }

  fclose(fp);

  return s;
}

int
print_model(const char *filename, const double r, bin2d_workspace *bin[3], magfit_workspace *w)
{
  int s = 0;
  const size_t nMLT = bin[0]->nx;
  const size_t nQD = bin[0]->ny;
  FILE *fp;
  size_t i, j;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_data: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Field %zu: magnetic local time (degrees)\n", i++);
  fprintf(fp, "# Field %zu: QD latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: X grid measurement (nT)\n", i++);
  fprintf(fp, "# Field %zu: Y grid measurement (nT)\n", i++);
  fprintf(fp, "# Field %zu: Z grid measurement (nT)\n", i++);
  fprintf(fp, "# Field %zu: X model (nT)\n", i++);
  fprintf(fp, "# Field %zu: Y model (nT)\n", i++);
  fprintf(fp, "# Field %zu: Z model (nT)\n", i++);

  for (i = 0; i < nMLT; ++i)
    {
      for (j = 1; j < nQD - 1; ++j)
        {
          double MLT, qdlat;
          double theta, phi;
          double B[3], B_model[3];
          size_t k;

          bin2d_xyval(i, j, &MLT, &qdlat, bin[0]);

          phi = MLT * M_PI / 180.0;
          theta = M_PI / 2.0 - qdlat * M_PI / 180.0;

          for (k = 0; k < 3; ++k)
            B[k] = bin2d_median(MLT, qdlat, bin[k]);

          magfit_eval_B(0.0, r, theta, phi, B_model, w);

          fprintf(fp, "%7.4f %7.4f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n",
                  MLT,
                  qdlat,
                  B[0],
                  B[1],
                  B[2],
                  B_model[0],
                  B_model[1],
                  B_model[2]);
        }

      fprintf(fp, "\n");
    }

  fclose(fp);

  return s;
}

int
make_grid(const size_t nMLT, const size_t nQD, const magdata * data, bin2d_workspace *bin[3], double *meanr)
{
  const double MLT_min = -180.0;
  const double MLT_max = 180.0;
  const double QD_min = -90.0;
  const double QD_max = 90.0;
  gsl_rstat_workspace *rstat_p = gsl_rstat_alloc();
  size_t i;

  for (i = 0; i < 3; ++i)
    bin[i] = bin2d_alloc(MLT_min, MLT_max, nMLT, QD_min, QD_max, nQD);

  for (i = 0; i < data->n; ++i)
    {
      double progress = (double)i / (double)data->n;
      double MLT_deg = 15.0 * data->MLT[i] - 180.0;
      double B[4];
      size_t j;

      if (!MAGDATA_ExistVector(data->flags[i]))
        continue;

      gsl_rstat_add(data->r[i], rstat_p);

      magdata_residual(i, B, data);

      for (j = 0; j < 3; ++j)
        bin2d_add_element(MLT_deg, data->qdlat[i], B[j], bin[j]);
    }

  *meanr = gsl_rstat_mean(rstat_p);

  gsl_rstat_free(rstat_p);

  return 0;
}

int
main_proc(const char *data_file, magdata * data)
{
  const char *data_file_surface = "data_surface.txt";
  const size_t nMLT = 360; /* 1 degree bin spacing */
  const size_t nQD = 180;  /* 1 degree bin spacing */
  const magfit_type *T = magfit_sheet;
  magfit_parameters params = magfit_default_parameters();
  magfit_workspace *w;
  bin2d_workspace *bin[3];
  struct timeval tv0, tv1;
  double r;                /* mean radius of data */
  double rnorm, snorm;
  size_t i, j;

  /* make grids of X, Y, Z */
  fprintf(stderr, "main_proc: gridding X,Y,Z data...");
  gettimeofday(&tv0, NULL);
  make_grid(nMLT, nQD, data, bin, &r);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds, mean radius = %g [km])\n", time_diff(tv0, tv1), r);

  /*params.nmax_int = 10;
  params.mmax_int = 10;*/
  params.nmax_ext = 0;
  params.mmax_ext = 0;

  /* don't fit Y component since it is mostly influenced by FAC */
  params.flags = MAGFIT_FLG_FIT_X | MAGFIT_FLG_FIT_Z;
  
  w = magfit_alloc(T, &params);

  fprintf(stderr, "main_proc: internal nmax = %zu\n", params.nmax_int);
  fprintf(stderr, "main_proc: internal mmax = %zu\n", params.mmax_int);

  fprintf(stderr, "main_proc: building LS system...");
  gettimeofday(&tv0, NULL);

#if 1
  for (i = 0; i < nMLT; ++i)
    {
      for (j = 1; j < nQD - 1; ++j)
        {
          double MLT, qdlat;
          double theta, phi;
          double B[3];
          size_t k;

          bin2d_xyval(i, j, &MLT, &qdlat, bin[0]);

          phi = MLT * M_PI / 180.0;
          theta = M_PI / 2.0 - qdlat * M_PI / 180.0;

          for (k = 0; k < 3; ++k)
            B[k] = bin2d_median(MLT, qdlat, bin[k]);

          magfit_add_datum(0.0, r, theta, phi, qdlat, B, w);
        }
    }
#else
  for (i = 0; i < data->n; ++i)
    {
      double progress = (double)i / (double)data->n;
      double MLT_deg = 15.0 * data->MLT[i] - 180.0;
      double phi = MLT_deg * M_PI / 180.0;
      double theta = M_PI / 2.0 - data->qdlat[i] * M_PI / 180.0;
      double B[4];

      if (!MAGDATA_ExistVector(data->flags[i]))
        continue;

      progress_bar(stdout, progress, 80);

      magdata_residual(i, B, data);

      magfit_add_datum(data->t[i], data->r[i], theta, phi, data->qdlat[i], B, w);
    }
#endif
    
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "main_proc: solving LS system...");
  gettimeofday(&tv0, NULL);
  magfit_fit(&rnorm, &snorm, w);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "main_proc: printing model map to %s...", data_file);
  print_model(data_file, r, bin, w);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main_proc: printing model map at Earth surface to %s...", data_file_surface);
  print_model(data_file_surface, R_EARTH_KM, bin, w);
  fprintf(stderr, "done\n");

  magfit_free(w);

  return 0;
}

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options]\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --input_file | -i magdata_input_file        - magdata data file\n");
}

int
main(int argc, char *argv[])
{
  char *outfile = "data.txt";
  magdata *data = NULL;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "input_file", required_argument, NULL, 'i' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "i:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'i':
            fprintf(stderr, "main: reading %s...", optarg);
            data = magdata_read(optarg, NULL);
            fprintf(stderr, "done (%zu data read)\n", data->n);
            break;

          default:
            break;
        }
    }

  if (!data)
    {
      print_help(argv);
      exit(1);
    }

  main_proc(outfile, data);

  magdata_free(data);

  return 0;
}
