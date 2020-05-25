/*
 * print_ezie_grid.c
 *
 * Compute magnetic field from 3D TIEGCM current and print 3D grid of results
 *
 * Grid ranges:
 * lon: [-180:0.1:180]
 * lat: [-89.5:0.1:89.5]
 * alt: [40:1:110]
 *
 * ./print_ezie_grid <-i tiegcm_nc_file>
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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rstat.h>

#include <mainlib/ml_common.h>
#include <mainlib/ml_bsearch.h>
#include <mainlib/ml_msynth.h>

#include <magfield/magfield.h>
#include <magfield/magfield_eval.h>

#include "magfit.h"
#include "tiegcm3d.h"

#define MAX_PTS            1000

/*
print_data()
  Print J grid for a fixed time and altitude

Inputs: data - tiegcm data
*/

int
print_data(const char *filename, const double r, magfield_eval_workspace *w)
{
  int s = 0;
  const double dlon = 0.5;
  const double dlat = 0.1;
  double lat, lon;
  size_t i;
  FILE *fp;
  msynth_workspace * core_p = msynth_shc_read(MSYNTH_CHAOS_FILE);

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_data: unable to open %s: %s\n",
              filename, strerror(errno));
    }

  msynth_set(1, 15, core_p);

  i = 1;
  fprintf(fp, "# Radius: %.2f [km]\n", r);
  fprintf(fp, "# Field %zu: geocentric latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: B_r (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_t (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_p (nT)\n", i++);

  for (lat = -89.5; lat <= 89.5; lat += dlat)
    {
      double theta = M_PI / 2.0 - lat * M_PI / 180.0;

      for (lon = -180.0; lon <= 180.0; lon += dlon)
        {
          double phi = lon * M_PI / 180.0;
          double B[4];

          magfield_eval_B(r * 1.0e3, theta, phi, B, w);

          fprintf(fp, "%.1f %.1f %.4f %.4f %.4f\n",
                  lat,
                  lon,
                  B[0] * 1.0e9,
                  B[1] * 1.0e9,
                  B[2] * 1.0e9);
        }
    }

  fclose(fp);

  return s;
}

int
main(int argc, char *argv[])
{
  tiegcm3d_data *data;
  struct timeval tv0, tv1;
  char *infile = NULL;
  char *output_dir = "EEJ_grid";
  int time_idx = 0;
  size_t lmax = 200;
  size_t mmax = 170;
  double r = R_EARTH_KM + 80.0;
  magfield_workspace *magfield_p;
  magfield_eval_workspace *magfield_eval_p;
  magfield_params params;
  size_t i, j, k;
  char buf[2048];

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "lmax", required_argument, NULL, 'l' },
          { "mmax", required_argument, NULL, 'm' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "i:l:m:o:r:t:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          case 'l':
            lmax = (size_t) atoi(optarg);
            break;

          case 'm':
            mmax = (size_t) atoi(optarg);
            break;

          case 'r':
            r = atof(optarg);
            break;

          case 't':
            time_idx = atol(optarg);
            break;

          case 'o':
            output_dir = optarg;
            break;

          default:
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i tiegcm3d_nc_file> [-t time_idx] [-r radius (km)] [-l lmax] [-m mmax] [-o output_dir]\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "input file = %s\n", infile);
  fprintf(stderr, "lmax       = %zu\n", lmax);
  fprintf(stderr, "mmax       = %zu\n", mmax);

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

  params.lmin = 1;
  params.lmax = lmax;
  params.mmax = mmax;
  params.nr = data->nr;
  params.ntheta = data->nlat;
  params.nphi = data->nlon;
  params.R = R_EARTH_M;
  params.grid_type = MAGFIELD_GAUSS;

  for (i = 0; i < data->nr; ++i)
    params.r[i] = data->r[i];

  fprintf(stderr, "rmin = %.1f [m]\n", params.r[0]);
  fprintf(stderr, "rmax = %.1f [m]\n", params.r[data->nr - 1]);

  fprintf(stderr, "\t allocating magfield workspace...");
  magfield_p = magfield_alloc(&params);
  magfield_eval_p = magfield_eval_alloc(&params);
  fprintf(stderr, "done\n");

  /* fill current grid */
  fprintf(stderr, "main: filling current grid...");

  for (i = 0; i < data->nr; ++i)
    {
      for (j = 0; j < data->nlat; ++j)
        {
          for (k = 0; k < data->nlon; ++k)
            {
              size_t idx = TIEGCM3D_IDX(0, i, j, k, data);

              magfield_grid_set(MAG_IDX_R, i, j, k, data->Jr[idx], magfield_p);
              magfield_grid_set(MAG_IDX_THETA, i, j, k, data->Jt[idx], magfield_p);
              magfield_grid_set(MAG_IDX_PHI, i, j, k, data->Jp[idx], magfield_p);
            }
        }
    }

  fprintf(stderr, "done\n");

  /* perform SH decomposition */
  fprintf(stderr, "main: performing SH decomposition...");
  magfield_decomp(magfield_p);
  fprintf(stderr, "done\n");

  magfield_eval_init(magfield_p->qtcoeff, magfield_p->qcoeff, magfield_p->pcoeff, magfield_eval_p);

  for (r = R_EARTH_KM + 40.0; r <= R_EARTH_KM + 110.0; r += 1.0)
    {
      sprintf(buf, "%s/grid_%g.txt", output_dir, r);
      fprintf(stderr, "main: writing grid data to %s (time idx = %d, r = %.2f [km])...", buf, time_idx, r);
      print_data(buf, r, magfield_eval_p);
      fprintf(stderr, "done\n");
    }

  tiegcm3d_free(data);
  magfield_free(magfield_p);
  magfield_eval_free(magfield_eval_p);

  return 0;
}
