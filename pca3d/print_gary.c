/*
 * print_magfield.c
 *
 * Compute magnetic field from 3D TIEGCM current and print results
 *
 * ./print_magfield <-i tiegcm_nc_file>
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

#include <magfield/magfield.h>

#include "tiegcm3d.h"

#define MAX_PTS            1000

/*
print_chi()
  Print J grid for a fixed time and altitude

Inputs: data - tiegcm data
*/

int
print_chi(const char *filename, const char *filename_matlab, const tiegcm3d_data *data, const int time_idx, const double r, magfield_workspace *w)
{
  int s = 0;
  const double eps = 1.0e-2;
  const size_t nlon = 360;
  const size_t nlat = 180;
  const double dlon = 360.0 / (double) nlon;
  const double dlat = 180.0 / (double) nlat;
  gsl_matrix *C = gsl_matrix_alloc(nlat, nlon);
  size_t i, j;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_data: unable to open %s: %s\n",
              filename, strerror(errno));
    }

  i = 1;
  fprintf(fp, "# Time: %ld (%.6f DOY)\n", data->t[time_idx], data->doy[time_idx] + data->ut[time_idx] / 24.0);
  fprintf(fp, "# Radius: %.2f (km) [%.2f km altitude]\n", r, r - 6371.2);
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: height-integrated chi (A)\n", i++);

  for (i = 0; i < nlon; ++i)
    {
      double lon = -180.0 + i * dlon;
      double phi = lon * M_PI / 180.0;

      for (j = 0; j < nlat; ++j)
        {
          double lat = -90.0 + j * dlat;
          double theta = M_PI / 2.0 - lat * M_PI / 180.0;
          double alt;
          double chi_HI = 0.0, chi;
          double dr = 10.0;

          for (alt = 90.0; alt <= 200.0; alt += dr)
            {
              magfield_eval_chi((R_EARTH_KM + alt) * 1.0e3, theta, phi, &chi, w);
              chi_HI += chi * dr * 1.0e3;
            }

          fprintf(fp, "%8.4f %8.4f %16.4e\n",
                  lon,
                  lat,
                  chi_HI);

          gsl_matrix_set(C, j, i, chi);
        }

      fprintf(fp, "\n");
    }

  print_octave2(C, "output/matlab", filename_matlab);

  fclose(fp);
  gsl_matrix_free(C);

  return s;
}

int
main(int argc, char *argv[])
{
  tiegcm3d_data *data;
  struct timeval tv0, tv1;
  char *infile = NULL;
  char *outfile_map = "data_chi.txt";
  magfield_workspace *magfield_p;
  magfield_params params;
  size_t i, j, k, tidx;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "i:o:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          case 'o':
            outfile_map = optarg;
            break;

          default:
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i tiegcm3d_nc_file> [-o output_file]\n", argv[0]);
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
  fprintf(stderr, "done\n");

  magfield_set_r(data->r, magfield_p);
  for (i = 0; i < data->nr; ++i)
    magfield_p->r[i] *= 1.0e3; /* convert to m */

  for (tidx = 0; tidx < data->nt; ++tidx)
    {
      char buf[1024], buf2[1024];

      sprintf(buf, "output/chi_%03d.txt", tidx);
      sprintf(buf2, "chi%03d", tidx);

      /* fill current grid */
      fprintf(stderr, "main: filling current grid...");

      for (i = 0; i < data->nr; ++i)
        {
          for (j = 0; j < data->nlat; ++j)
            {
              for (k = 0; k < data->nlon; ++k)
                {
                  size_t idx = TIEGCM3D_IDX(tidx, i, j, k, data);

                  magfield_current_set(MAG_IDX_R, i, j, k, data->Jr[idx], magfield_p);
                  magfield_current_set(MAG_IDX_THETA, i, j, k, data->Jt[idx], magfield_p);
                  magfield_current_set(MAG_IDX_PHI, i, j, k, data->Jp[idx], magfield_p);
                }
            }
        }

      fprintf(stderr, "done\n");

      /* perform SH decomposition */
      fprintf(stderr, "main: performing SH decomposition...");
      magfield_decomp(magfield_p);
      fprintf(stderr, "done\n");

      fprintf(stderr, "main: writing chi grid data to %s (time idx = %d)...", buf, tidx);
      print_chi(buf, buf2, data, tidx, R_EARTH_KM + 110.0, magfield_p);
      fprintf(stderr, "done\n");
    }

  tiegcm3d_free(data);
  magfield_free(magfield_p);

  return 0;
}
