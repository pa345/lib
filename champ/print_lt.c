/*
 * print_lt.c
 *
 * Print time series of CHAMP ascending node local time
 *
 * ./print_lt <-i champ_index_file>
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

#include <satdata/satdata.h>

#include <common/common.h>

int
main(int argc, char *argv[])
{
  char *outfile = "lt.dat";
  satdata_mag *data;
  struct timeval tv0, tv1;
  int c;
  char *infile = NULL;

  while ((c = getopt(argc, argv, "i:o:")) != (-1))
    {
      switch (c)
        {
          case 'i':
            infile = optarg;
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
      fprintf(stderr, "Usage: %s <-i champ_index_file> [-o output_file]\n",
              argv[0]);
      exit(1);
    }

  fprintf(stderr, "input file = %s\n", infile);

  fprintf(stderr, "Reading %s...", infile);
  gettimeofday(&tv0, NULL);

  data = satdata_swarm_read_idx(infile, 1);
  if (!data)
    {
      fprintf(stderr, "main: error reading %s\n", infile);
      exit(1);
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu records read, %g seconds)\n", data->n,
          time_diff(tv0, tv1));

  fprintf(stderr, "main: printing local time data to %s...", outfile);
  satdata_print_lt(outfile, data);
  fprintf(stderr, "done\n");

  satdata_mag_free(data);

  return 0;
} /* main() */
