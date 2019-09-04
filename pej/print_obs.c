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

#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>

#include <mainlib/ml_common.h>
#include <mainlib/ml_bin2d.h>
#include <mainlib/ml_magdata.h>

int
my_bin2d_print(const char *filename, const bin2d_workspace *w)
{
  int s = 0;
  size_t i, j;
  FILE *fp;
  double array[MAX_DATA_PER_BIN];
  double work[3 * MAX_DATA_PER_BIN];
  int work_int[5 * MAX_DATA_PER_BIN];

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "bin2d_print: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Field %zu: x\n", i++);
  fprintf(fp, "# Field %zu: y\n", i++);
  fprintf(fp, "# Field %zu: mean\n", i++);
  fprintf(fp, "# Field %zu: median\n", i++);
  fprintf(fp, "# Field %zu: stddev\n", i++);
  fprintf(fp, "# Field %zu: Q_n\n", i++);
  fprintf(fp, "# Field %zu: number of data points in bin\n", i++);

  for (i = 0; i < w->nx; ++i)
    {
      for (j = 0; j < w->ny; ++j)
        {
          size_t ndata = bin2d_fill(i, j, array, w);
          double x, y;
          double mean, sd, median, Qn;
          size_t n;

          bin2d_xyval(i, j, &x, &y, w);

          mean = bin2d_mean(x, y, w);
          sd = bin2d_sd(x, y, w);
          median = bin2d_median(x, y, w);
          n = bin2d_n(x, y, w);

          gsl_sort(array, 1, ndata);
          Qn = gsl_stats_Qn_from_sorted_data(array, 1, ndata, work, work_int);

          fprintf(fp, "%f %f %.12e %.12e %.12e %.12e %zu\n", x, y, mean, median, sd, Qn, n);
        }

      fprintf(fp, "\n");
    }

  fclose(fp);

  return s;
}

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
  const size_t nMLT = 24;
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

  my_bin2d_print("data_X.txt", bin2d_X);
  my_bin2d_print("data_Y.txt", bin2d_Y);
  my_bin2d_print("data_Z.txt", bin2d_Z);

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
