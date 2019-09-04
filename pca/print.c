/*
 * print.c
 *
 * Print contents of tiegcm data files
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

#include "tiegcm.h"

/*
print_data()
  Print Bx/By/Bz grid for a fixed time

Inputs: data - tiegcm data
*/

int
print_data(const char *filename, const tiegcm_data *data, const int time_idx)
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
  fprintf(fp, "# Time: %ld (%.6f)\n", data->t[time_idx], data->doy[time_idx] + data->ut[time_idx] / 24.0);
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: B_x (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_y (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_z (nT)\n", i++);

  for (ilon = 0; ilon < data->nlon; ++ilon)
    {
      for (ilat = 0; ilat < data->nlat; ++ilat)
        {
          size_t idx = TIEGCM_BIDX(time_idx, ilat, ilon, data);

          fprintf(fp, "%8.4f %8.4f %8.2f %8.2f %8.2f\n",
                  data->glon[ilon],
                  data->glat[ilat],
                  data->Bx[idx] * 1.0e9,
                  data->By[idx] * 1.0e9,
                  data->Bz[idx] * 1.0e9);
        }

      fprintf(fp, "\n");
    }

  fclose(fp);

  return s;
}

int
main(int argc, char *argv[])
{
  tiegcm_data *data;
  struct timeval tv0, tv1;
  char *infile = NULL;
  char *outfile = "data.txt";
  int time_idx = 0;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "i:o:t:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          case 't':
            time_idx = atol(optarg);
            break;

          case 'o':
            outfile = optarg;
            break;

          default:
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i tiegcm_nc_file> [-t time_idx] [-o output_file]\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "input file = %s\n", infile);

  fprintf(stderr, "main: reading %s...", infile);
  gettimeofday(&tv0, NULL);

  data = tiegcm_read(infile, NULL);
  if (!data)
    {
      fprintf(stderr, "main: error reading %s\n", infile);
      exit(1);
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu records read, %g seconds)\n", data->nt,
          time_diff(tv0, tv1));

  fprintf(stderr, "main: writing data to %s...", outfile);
  print_data(outfile, data, time_idx);
  fprintf(stderr, "done\n");

  tiegcm_free(data);

  return 0;
}
