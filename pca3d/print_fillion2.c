/*
 * print_fillion.c
 *
 * 1) Read file of positions where currents are desired
 * 2) Read 3D TIEGCM grid file
 * 3) Output current values along desired path
 *
 * ./print_fillion2 <-i tiegcm_nc_file> <-j current_path_file>
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

#include <common/common.h>
#include <common/bsearch.h>
#include <common/interp.h>

#include <magfield/magfield.h>

#include "tiegcm3d.h"

#define MAX_N            15000

typedef struct
{
  time_t t[MAX_N];
  double fday[MAX_N];
  double r[MAX_N];
  double lat[MAX_N];
  double lon[MAX_N];
} path_data;

size_t
read_path(const char *filename, path_data *data)
{
  FILE *fp;
  char buffer[2048];
  size_t n = 0;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "read_path: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  while (fgets(buffer, sizeof(buffer), fp) != NULL)
    {
      int c;
      double fday, r, colat, lon;

      c = sscanf(buffer, "%lf %lf %lf %lf", &fday, &r, &colat, &lon);
      if (c != 4)
        continue;

      data->t[n] = fday2timet(fday);
      data->fday[n] = fday;
      data->r[n] = r;
      data->lat[n] = 90.0 - colat;
      data->lon[n] = lon;

      if (++n >= MAX_N)
        {
          fprintf(stderr, "MAX_N too small\n");
          break;
        }
    }

  fclose(fp);

  return n;
}

int
fill_grid(const size_t time_idx, tiegcm3d_data *data, magfield_workspace *magfield_p)
{
  size_t i, j, k;

  for (i = 0; i < data->nr; ++i)
    {
      for (j = 0; j < data->nlat; ++j)
        {
          for (k = 0; k < data->nlon; ++k)
            {
              size_t idx = TIEGCM3D_IDX(time_idx, i, j, k, data);

              magfield_current_set(MAG_IDX_R, i, j, k, data->Jr[idx], magfield_p);
              magfield_current_set(MAG_IDX_THETA, i, j, k, data->Jt[idx], magfield_p);
              magfield_current_set(MAG_IDX_PHI, i, j, k, data->Jp[idx], magfield_p);
            }
        }
    }

  return GSL_SUCCESS;
}

int
main(int argc, char *argv[])
{
  tiegcm3d_data *data;
  path_data data_path;
  struct timeval tv0, tv1;
  char *infile = NULL;
  char *outfileJ = "satelliteJ.txt";
  size_t npath = 0;
  time_t t;
  struct tm *tmp;
  size_t tidx = 0;
  magfield_workspace *magfield_p, *magfield_p2;
  magfield_params params;
  size_t i, k;
  FILE *fp;
  int c;

  while ((c = getopt(argc, argv, "i:j:")) != (-1))
    {
      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          case 'j':
            fprintf(stderr, "main: reading %s...", optarg);
            npath = read_path(optarg, &data_path);
            fprintf(stderr, "done (%zu points read)\n", npath);
            break;

          default:
            break;
        }
    }

  if (!infile || npath == 0)
    {
      fprintf(stderr, "Usage: %s <-i tiegcm3d_nc_file> <-j current_path_file>\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "input file = %s\n", infile);

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
  params.lmax = 60;
  params.mmax = 30;
  params.nr = data->nr;
  params.ntheta = data->nlat;
  params.nphi = data->nlon;
  params.rmin = data->r[0] * 1.0e3;
  params.rmax = data->r[data->nr - 1] * 1.0e3;
  params.R = R_EARTH_M;
  params.grid_type = MAGFIELD_GAUSS;

  fprintf(stderr, "\t allocating magfield workspace...");
  magfield_p = magfield_alloc(&params);
  magfield_p2 = magfield_alloc(&params);
  fprintf(stderr, "done\n");

  magfield_set_r(data->r, magfield_p);
  magfield_set_r(data->r, magfield_p2);
  for (i = 0; i < data->nr; ++i)
    {
      magfield_p->r[i] *= 1.0e3; /* convert to m */
      magfield_p2->r[i] *= 1.0e3; /* convert to m */
    }

  fp = fopen(outfileJ, "w");

  i = 1;
  fprintf(fp, "# Field %zu: timestamp (MJD2000)\n", i++);
  fprintf(fp, "# Field %zu: geocentric radius (km)\n", i++);
  fprintf(fp, "# Field %zu: geocentric latitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: geocentric longitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: J_r (A/m2)\n", i++);
  fprintf(fp, "# Field %zu: J_t (A/m2)\n", i++);
  fprintf(fp, "# Field %zu: J_p (A/m2)\n", i++);
  fprintf(fp, "# Field %zu: B_r (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_t (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_p (nT)\n", i++);

  t = data_path.t[0];
  tmp = gmtime(&t);

  /* TIEGCM run is for 2009, change year and search for month/day in tiegcm3d data */
  tmp->tm_year = 2009 - 1900;
  t = mktime(tmp);

  tidx = bsearch_timet(data->t, t, 0, data->nt - 1);

  fprintf(stderr, "main: filling current grid for timestamp %zu...", tidx);
  fill_grid(tidx, data, magfield_p);
  fill_grid(tidx+1, data, magfield_p2);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: performing SH decomposition for timestamp %zu...", tidx);
  magfield_decomp(magfield_p);
  magfield_decomp(magfield_p2);
  fprintf(stderr, "done\n");

  for (i = 0; i < npath; ++i)
    {
      double fday = data_path.fday[i];
      double r = data_path.r[i];
      double theta = M_PI / 2.0 - data_path.lat[i] * M_PI / 180.0;
      double phi = data_path.lon[i] * M_PI / 180.0;
      double B1[4], B2[4], B[4];
      double J1[3], J2[3], J[3];

      t = data_path.t[i];
      tmp = gmtime(&t);
      tmp->tm_year = 2009 - 1900;
      t = mktime(tmp);

      if (t > data->t[tidx + 1])
        {
          tidx = bsearch_timet(data->t, t, 0, data->nt - 1);

          fprintf(stderr, "main: filling current grid for timestamp %zu...", tidx);
          fill_grid(tidx, data, magfield_p);
          fill_grid(tidx+1, data, magfield_p2);
          fprintf(stderr, "done\n");

          fprintf(stderr, "main: performing SH decomposition for timestamp %zu...", tidx);
          magfield_decomp(magfield_p);
          magfield_decomp(magfield_p2);
          fprintf(stderr, "done\n");
        }

      magfield_eval_B(r * 1.0e3, theta, phi, B1, magfield_p);
      magfield_eval_B(r * 1.0e3, theta, phi, B2, magfield_p2);

      magfield_eval_J(r * 1.0e3, theta, phi, J1, magfield_p);
      magfield_eval_J(r * 1.0e3, theta, phi, J2, magfield_p2);

      for (k = 0; k < 3; ++k)
        {
          B[k] = interp1d((double) data->t[tidx], (double) data->t[tidx+1], B1[k], B2[k], (double) t);
          J[k] = interp1d((double) data->t[tidx], (double) data->t[tidx+1], J1[k], J2[k], (double) t);
        }

      fprintf(fp, "%f %.4f %.4f %.4f %.12e %.12e %.12e %.4f %.4f %.4f\n",
              fday,
              r,
              data_path.lat[i],
              data_path.lon[i],
              J[0],
              J[1],
              J[2],
              B[0] * 1.0e9,
              B[1] * 1.0e9,
              B[2] * 1.0e9);
    }

  fclose(fp);

  tiegcm3d_free(data);
  magfield_free(magfield_p);

  return 0;
}