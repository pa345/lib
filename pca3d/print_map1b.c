/*
 * print_map1b.c
 *
 * 1. Read TIEGCM 3D current grids
 * 2. Read SH coefficients computed from stage1b
 * 3. Print maps for comparison
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <sys/time.h>
#include <assert.h>
#include <omp.h>
#include <complex.h>
#include <string.h>
#include <errno.h>

#include <fftw3.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_blas.h>

#include <mainlib/ml_common.h>
#include <mainlib/ml_oct.h>
#include <mainlib/ml_bsearch.h>

#include <magfield/magfield.h>
#include <magfield/magfield_eval.h>

#include "io.h"
#include "pca3d.h"
#include "tiegcm3d.h"

#define MAX_T        100000

/*
print_map()
  Print J grid for a fixed time and altitude from TIEGCM and magfield

Inputs: data - tiegcm data
*/

int
print_map(const char *filename, const tiegcm3d_data *data, const int time_idx, const int ir,
          magfield_workspace * w, magfield_eval_workspace *eval_p)
{
  int s = 0;
  size_t i;
  size_t ilat, ilon;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_map: unable to open %s: %s\n",
              filename, strerror(errno));
    }

  i = 1;
  fprintf(fp, "# Time: %ld (%.6f DOY)\n", data->t[time_idx], data->doy[time_idx] + data->ut[time_idx] / 24.0);
  fprintf(fp, "# Radius: %.2f (km) [%.2f km altitude]\n", data->r[ir] * 1.0e-3, data->r[ir] * 1.0e-3 - R_EARTH_KM);
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: TIEGCM J_r (nA/m^2)\n", i++);
  fprintf(fp, "# Field %zu: TIEGCM J_t (nA/m^2)\n", i++);
  fprintf(fp, "# Field %zu: TIEGCM J_p (nA/m^2)\n", i++);
  fprintf(fp, "# Field %zu: magfield J_r (nA/m^2)\n", i++);
  fprintf(fp, "# Field %zu: magfield J_t (nA/m^2)\n", i++);
  fprintf(fp, "# Field %zu: magfield J_p (nA/m^2)\n", i++);

  for (ilon = 0; ilon < data->nlon; ++ilon)
    {
      double phi = w->phi[ilon];

      for (ilat = 0; ilat < data->nlat; ++ilat)
        {
          size_t idx = TIEGCM3D_IDX(time_idx, ir, ilat, ilon, data);
          double theta = w->theta[ilat];
          double J[4];

          magfield_eval_J(data->r[ir], theta, phi, J, eval_p);

          fprintf(fp, "%8.4f %8.4f %16.4e %16.4e %16.4e %16.4e %16.4e %16.4e\n",
                  data->glon[ilon],
                  data->glat[ilat],
                  data->Jr[idx] * 1.0e9,
                  data->Jt[idx] * 1.0e9,
                  data->Jp[idx] * 1.0e9,
                  J[0] * 1.0e9,
                  J[1] * 1.0e9,
                  J[2] * 1.0e9);
        }

      fprintf(fp, "\n");
    }

  fclose(fp);

  return s;
}

/*
print_alt()
  Print J grid for a fixed time and latitude from TIEGCM and magfield

Inputs: data - tiegcm data
*/

int
print_alt(const char *filename, const tiegcm3d_data *data, const int time_idx, const size_t ilat,
          magfield_workspace * w, magfield_eval_workspace *eval_p)
{
  int s = 0;
  const double theta = w->theta[ilat];
  size_t i;
  size_t ir, ilon;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_alt: unable to open %s: %s\n",
              filename, strerror(errno));
    }

  i = 1;
  fprintf(fp, "# Time: %ld (%.6f DOY)\n", data->t[time_idx], data->doy[time_idx] + data->ut[time_idx] / 24.0);
  fprintf(fp, "# Latitude: %.2f (deg)\n", data->glat[ilat]);
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: radius (km)\n", i++);
  fprintf(fp, "# Field %zu: TIEGCM J_r (nA/m^2)\n", i++);
  fprintf(fp, "# Field %zu: TIEGCM J_t (nA/m^2)\n", i++);
  fprintf(fp, "# Field %zu: TIEGCM J_p (nA/m^2)\n", i++);
  fprintf(fp, "# Field %zu: magfield J_r (nA/m^2)\n", i++);
  fprintf(fp, "# Field %zu: magfield J_t (nA/m^2)\n", i++);
  fprintf(fp, "# Field %zu: magfield J_p (nA/m^2)\n", i++);

  for (ilon = 0; ilon < data->nlon; ++ilon)
    {
      double phi = w->phi[ilon];

      for (ir = 0; ir < data->nr; ++ir)
        {
          size_t idx = TIEGCM3D_IDX(time_idx, ir, ilat, ilon, data);
          double J[4];

          magfield_eval_J(data->r[ir], theta, phi, J, eval_p);

          fprintf(fp, "%8.4f %8.4f %16.4e %16.4e %16.4e %16.4e %16.4e %16.4e\n",
                  data->glon[ilon],
                  data->r[ir] * 1.0e-3 - R_EARTH_KM,
                  data->Jr[idx] * 1.0e9,
                  data->Jt[idx] * 1.0e9,
                  data->Jp[idx] * 1.0e9,
                  J[0] * 1.0e9,
                  J[1] * 1.0e9,
                  J[2] * 1.0e9);
        }

      fprintf(fp, "\n");
    }

  fclose(fp);

  return s;
}

/*
main_proc()
  Compute SH decomposition of each radial shell of TIEGCM 3D current
values; store SH coefficients to disk

Inputs: prefix        - output file prefix
        tidx          - on input, starting time index
                        on output, ending time index
        tsteps        - (output) array of timestamps
        data          - TIEGCM data
        residual_file - (optional) write file of residuals between TIEGCM/magfield for first timestep
        w             - magfield workspace
        eval_p        - magfield eval workspace

Return: success/error
*/

int
main_proc(const size_t tidx, const pca3d_data * pdata, const tiegcm3d_data * tdata, magfield_eval_workspace * eval_p)
{
  int s = 0;
  const char *latlon_file = "map_latlon.txt";
  const char *rlon_file = "map_rlon.txt";
  const double alt = 310.0;
  const double lat = 15.0;
  int ir = (int) bsearch_double(tdata->r, R_EARTH_M + alt * 1.0e3, 0, tdata->nr - 1);
  int ilat = (int) bsearch_double(tdata->glat, lat, 0, tdata->nlat - 1);
  magfield_workspace * w = pdata->w;
  char buf[2048];

  sprintf(buf, "%s_%03zu.dat", PCA3D_STAGE1B_SH_PREFIX, tidx + 1);
  fprintf(stderr, "main_proc: reading %s...", buf);
  magfield_read_SH(buf, w);
  fprintf(stderr, "done\n");

  magfield_eval_init(w->qtcoeff, w->qcoeff, w->pcoeff, eval_p);

  fprintf(stderr, "main_proc: writing lat/lon map for time index %zu to %s...", tidx, latlon_file);
  print_map(latlon_file, tdata, tidx, ir, w, eval_p);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main_proc: writing r/lon map for time index %zu to %s...", tidx, rlon_file);
  print_alt(rlon_file, tdata, tidx, ilat, w, eval_p);
  fprintf(stderr, "done\n");

  return s;
}

int
main(int argc, char *argv[])
{
  char *residual_file = NULL;
  struct timeval tv0, tv1;
  magfield_eval_workspace * magfield_eval_p = NULL;
  size_t tidx = 0;
  size_t eval_lmax = 0;
  size_t eval_mmax = 0;
  tiegcm3d_data *data = NULL;
  pca3d_data pdata;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "residual_file", required_argument, NULL, 'r' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "l:m:r:t:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'r':
            residual_file = optarg;
            break;

          case 't':
            tidx = (size_t) atoi(optarg);
            break;

          case 'l':
            eval_lmax = (size_t) atoi(optarg);
            break;

          case 'm':
            eval_mmax = (size_t) atoi(optarg);
            break;

          default:
            fprintf(stderr, "Usage: %s [-l eval_lmax] [-m eval_mmax] [-r residual_file] [-t tidx] file1.nc file2.nc ...\n", argv[0]);
            exit(1);
            break;
        }
    }

  if (optind >= argc)
    {
      fprintf(stderr, "Usage: %s [-l eval_lmax] [-m eval_mmax] [-r residual_file] [-t tidx] file1.nc file2.nc ...\n", argv[0]);
      exit(1);
    }

  while (optind < argc)
    {
      fprintf(stderr, "main: reading %s...", argv[optind]);
      gettimeofday(&tv0, NULL);

      data = tiegcm3d_read(argv[optind], data);
      if (!data)
        {
          fprintf(stderr, "main: error reading %s\n", argv[optind]);
          exit(1);
        }

      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%zu records read, %g seconds)\n", data->nt,
              time_diff(tv0, tv1));

      ++optind;
    }

  fprintf(stderr, "main: allocating magfield workspace from %s...", PCA3D_STAGE1B_DATA);
  pdata = pca3d_read_data(PCA3D_STAGE1B_DATA);
  fprintf(stderr, "done (%zu timestamps)\n", pdata.nt);

  fprintf(stderr, "main: allocating magfield eval workspace...");
  magfield_eval_p = magfield_eval_alloc(&(pdata.w->params));
  fprintf(stderr, "done\n");

  if (eval_lmax > 0 && eval_mmax > 0)
    magfield_eval_set(eval_lmax, eval_mmax, magfield_eval_p);

  fprintf(stderr, "main: eval_lmax = %zu\n", magfield_eval_p->eval_lmax);
  fprintf(stderr, "main: eval_mmax = %zu\n", magfield_eval_p->eval_mmax);

  main_proc(tidx, &pdata, data, magfield_eval_p);

  tiegcm3d_free(data);

  return 0;
}
