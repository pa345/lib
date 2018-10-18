/*
 * swarmcdf2fday.c
 *
 * 1. Read Swarm CDF data file
 * 2. Compute fday of first data value and output
 *
 * Usage:
 *
 * ./swarmcdf2fday <-i swarm_cdf_file>
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
#include <gsl/gsl_interp.h>

#include <satdata/satdata.h>
#include <indices/indices.h>

#include <common/common.h>

int
main(int argc, char *argv[])
{
  satdata_mag *data;
  int c;
  char *infile = NULL;
  double fday;

  while ((c = getopt(argc, argv, "i:")) != (-1))
    {
      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          default:
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i swarm_cdf_file>\n",
              argv[0]);
      exit(1);
    }

  fprintf(stderr, "input file = %s\n", infile);

  fprintf(stderr, "Reading %s...", infile);
  data = satdata_swarm_read(infile, NULL);
  if (!data)
    {
      fprintf(stderr, "main: error reading %s\n", infile);
      exit(1);
    }
  fprintf(stderr, "done (%zu records read)\n", data->n);

  fday = satdata_epoch2fday(data->t[0]);
  printf("%f\n", fday);

  satdata_mag_free(data);

  return 0;
} /* main() */
