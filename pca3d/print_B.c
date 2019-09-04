/*
 * print_B.c
 *
 * Print contents of tiegcm data files with B grids
 *
 * ./print <-i tiegcm_nc_file>
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

#include <mainlib/ml_common.h>
#include <mainlib/ml_bsearch.h>

#include "tiegcm3d.h"

/*
print_data()
  Print B grid for a fixed time and altitude

Inputs: data - tiegcm data
*/

int
print_data(const char *filename, const tiegcm3d_data *data, const int time_idx, const int ir)
{
  int s = 0;
  size_t i;
  size_t ilat, ilon;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_data: unable to open %s: %s\n",
              filename, strerror(errno));
    }

  i = 1;
  fprintf(fp, "# Time: %ld (%.6f DOY)\n", data->t[time_idx], data->doy[time_idx] + data->ut[time_idx] / 24.0);
  fprintf(fp, "# Radius: %.2f (km) [%.2f km altitude]\n", data->r[ir], data->r[ir] - 6371.2);
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: B_r (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_t (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_p (nT)\n", i++);

  for (ilon = 0; ilon < data->nlon; ++ilon)
    {
      for (ilat = 0; ilat < data->nlat; ++ilat)
        {
          size_t idx = TIEGCM3D_IDX(time_idx, ir, ilat, ilon, data);

          fprintf(fp, "%8.4f %8.4f %16.4e %16.4e %16.4e\n",
                  data->glon[ilon],
                  data->glat[ilat],
                  data->Jr[idx],
                  data->Jt[idx],
                  data->Jp[idx]);
        }

      fprintf(fp, "\n");
    }

  fclose(fp);

  return s;
}

/*
print_time()
  Print time series of Br/Bt/Bp for a fixed location

Inputs: data - tiegcm data
*/

int
print_time(const char *filename, const tiegcm3d_data *data, const int ir, const int ilat, const int ilon)
{
  int s = 0;
  size_t it;
  FILE *fp;
  size_t i;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_time: unable to open %s: %s\n",
              filename, strerror(errno));
    }

  i = 1;
  fprintf(fp, "# Latitude: %.2f (deg)\n", data->glat[ilat]);
  fprintf(fp, "# Longitude: %.2f (deg)\n", data->glon[ilon]);
  fprintf(fp, "# Radius: %.2f (km) [%.2f km altitude]\n", data->r[ir], data->r[ir] - 6371.2);
  fprintf(fp, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
  fprintf(fp, "# Field %zu: B_r (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_t (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_p (nT)\n", i++);

  for (it = 0; it < data->nt; ++it)
    {
      size_t idx = TIEGCM3D_IDX(it, ir, ilat, ilon, data);

      fprintf(fp, "%ld %16.4e %16.4e %16.4e\n",
              data->t[it],
              data->Jr[idx],
              data->Jt[idx],
              data->Jp[idx]);
    }

  fclose(fp);

  return s;
}

/*
print_alt()
  Print altitude/latitude map of Br/Bt/Bp for fixed time and longitude

Inputs: data - tiegcm data
*/

int
print_alt(const char *filename, const tiegcm3d_data *data, const int it, const int ilon)
{
  int s = 0;
  FILE *fp;
  size_t i;
  size_t ir, ilat;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_alt: unable to open %s: %s\n",
              filename, strerror(errno));
    }

  i = 1;
  fprintf(fp, "# Time: %ld (%.6f DOY)\n", data->t[it], data->doy[it] + data->ut[it] / 24.0);
  fprintf(fp, "# Longitude: %.2f (deg)\n", data->glon[ilon]);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: altitude (km)\n", i++);
  fprintf(fp, "# Field %zu: B_r (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_t (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_p (nT)\n", i++);

  for (ilat = 0; ilat < data->nlat; ++ilat)
    {
      for (ir = 0; ir < data->nr; ++ir)
        {
          size_t idx = TIEGCM3D_IDX(it, ir, ilat, ilon, data);

          fprintf(fp, "%8.4f %8.4f %16.4e %16.4e %16.4e\n",
                  data->glat[ilat],
                  data->r[ir] - R_EARTH_KM,
                  data->Jr[idx],
                  data->Jt[idx],
                  data->Jp[idx]);
        }

      fprintf(fp, "\n");
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
  char *outfile_map = "data_map.txt";
  char *outfile_time = "data_time.txt";
  char *outfile_alt = "data_alt.txt";
  int time_idx = 0;
  double lon = 150.0; /* desired longitude */
  double lat = 8.0;   /* desired latitude */
  int r_idx = 0;
  int lat_idx, lon_idx;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "i:o:r:t:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          case 'r':
            r_idx = atol(optarg);
            break;

          case 't':
            time_idx = atol(optarg);
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
      fprintf(stderr, "Usage: %s <-i tiegcm3d_nc_file> [-t time_idx] [-o output_file]\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "input file = %s\n", infile);

  fprintf(stderr, "main: reading %s...", infile);
  gettimeofday(&tv0, NULL);

  data = tiegcm3d_read_B(infile, NULL);
  if (!data)
    {
      fprintf(stderr, "main: error reading %s\n", infile);
      exit(1);
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu records read, %g seconds)\n", data->nt,
          time_diff(tv0, tv1));

  /* locate index of desired lat/lon */
  lat_idx = bsearch_double(data->glat, lat, 0, data->nlat - 1);
  lon_idx = bsearch_double(data->glon, lon, 0, data->nlon - 1);

  fprintf(stderr, "main: lat_idx = %d (%.2f [deg])\n", lat_idx, data->glat[lat_idx]);
  fprintf(stderr, "main: lon_idx = %d (%.2f [deg])\n", lon_idx, data->glon[lon_idx]);

  fprintf(stderr, "main: writing grid data to %s (time idx = %d, r idx = %d)...", outfile_map, time_idx, r_idx);
  print_data(outfile_map, data, time_idx, r_idx);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: writing time series data to %s...", outfile_time);
  print_time(outfile_time, data, r_idx, lat_idx, lon_idx);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: writing altitude data to %s...", outfile_alt);
  print_alt(outfile_alt, data, time_idx, lon_idx);
  fprintf(stderr, "done\n");

  tiegcm3d_free(data);

  return 0;
}
