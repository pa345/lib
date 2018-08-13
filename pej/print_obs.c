/*
 * print_obs.c
 *
 * Print MLT/QDlat grid of observatory data
 *
 * ./print_obs <-i magdata_file>
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

#include <common/common.h>
#include <common/bin2d.h>

#include "magdata.h"

/*
print_grid()
  Output data points in ASCII format, only printing
time, position and magnetic measurements for points selected
for main field modeling

Inputs: filename - where to store data
        data     - data

Return: success/error
*/

int
print_grid(const char *filename, const magdata *data)
{
  int s = 0;
  const size_t nMLT = 180;
  const size_t nqd = 90;
  size_t i;
  bin2d_workspace *bin2d_X = bin2d_alloc(-180.0, 180.0, nMLT, -90.0, 90.0, nqd);
  bin2d_workspace *bin2d_Y = bin2d_alloc(-180.0, 180.0, nMLT, -90.0, 90.0, nqd);
  bin2d_workspace *bin2d_Z = bin2d_alloc(-180.0, 180.0, nMLT, -90.0, 90.0, nqd);

  /* make grid */
  for (i = 0; i < data->n; ++i)
    {
      double MLT_deg = 15.0 * (data->MLT[i] - 24.0) + 180.0;

      bin2d_add_element(MLT_deg, data->qdlat[i], data->Bx_nec[i], bin2d_X);
      bin2d_add_element(MLT_deg, data->qdlat[i], data->By_nec[i], bin2d_Y);
      bin2d_add_element(MLT_deg, data->qdlat[i], data->Bz_nec[i], bin2d_Z);
    }

  bin2d_print("data_X.txt", bin2d_X);
  bin2d_print("data_Y.txt", bin2d_Y);
  bin2d_print("data_Z.txt", bin2d_Z);

  bin2d_free(bin2d_X);
  bin2d_free(bin2d_Y);
  bin2d_free(bin2d_Z);

  return s;
}

int
main(int argc, char *argv[])
{
  magdata *data = NULL;
  char *outfile = "grid.dat";
  struct timeval tv0, tv1;
  int c;

  while ((c = getopt(argc, argv, "i:o:")) != (-1))
    {
      switch (c)
        {
          case 'i':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);

            data = magdata_read(optarg, NULL);
            if (!data)
              {
                fprintf(stderr, "main: error reading %s\n", optarg);
                exit(1);
              }

            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu records read, %g seconds)\n", data->n,
                    time_diff(tv0, tv1));
            break;

          case 'o':
            outfile = optarg;
            break;

          default:
            break;
        }
    }

  if (!data)
    {
      fprintf(stderr, "Usage: %s <-i magdata_file> [-o output_file]\n",
              argv[0]);
      exit(1);
    }

  fprintf(stderr, "main: writing output to %s...", outfile);
  print_grid(outfile, data);
  fprintf(stderr, "done\n");

  magdata_free(data);

  return 0;
}
