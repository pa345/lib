/*
 * invert_write_tmode.c
 *
 * Read temporal mode ASCII files provided by Gary and write them
 * to a binary format
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

#include <mainlib/ml_common.h>

#include "invert_tmode.h"

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options]\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --mode_dir        | -d directory              - directory containing ASCII temporal modes\n");
  fprintf(stderr, "\t --output_file     | -o output_file            - binary output data file (magdata format)\n");
}

int
main(int argc, char *argv[])
{
  const size_t nfreq_max = 20;
  const size_t N = 192672; /* length of time series for each mode */
  size_t * nmodes;
  size_t nfreq;
  char *mode_dir = "/data/palken/Gary_Project/Temporal_Modes";
  char *output_file = NULL;
  double freq; /* frequency of mode in cpd */
  invert_tmode_workspace * w;
  char buf[512];

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "mode_dir", required_argument, NULL, 'd' },
          { "output_file", required_argument, NULL, 'o' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "d:o:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'd':
            mode_dir = optarg;
            break;

          case 'o':
            output_file = optarg;
            break;

          default:
            break;
        }
    }

  if (!output_file)
    {
      print_help(argv);
      exit(1);
    }

  nmodes = calloc(nfreq_max, sizeof(size_t));

  for (nfreq = 0; nfreq < nfreq_max; ++nfreq)
    {
      size_t j;

      for (j = 1; ; j++)
        {
          sprintf(buf, "%s/mode_%02zu_%02zu.txt", mode_dir, nfreq + 1, j);
          if (access(buf, F_OK) != -1)
            {
              nmodes[nfreq]++;
            }
          else
            {
              /* file does not exist */
              break;
            }
        }

      if (j == 1)
        break; /* this frequency band does not exist */
    }

  fprintf(stderr, "number of frequency bands = %zu\n", nfreq);
  
  w = invert_tmode_alloc(nfreq, NULL, nmodes, N);

  /* now read all ASCII files and store in w */
  {
    size_t i, j;

    for (i = 0; i < nfreq; ++i)
      {
        double freq;

        for (j = 0; j < nmodes[i]; ++j)
          {
            sprintf(buf, "%s/mode_%02zu_%02zu.txt", mode_dir, i + 1, j + 1);

            fprintf(stderr, "main: processing %s...", buf);
            invert_tmode_read_ascii(buf, i, j, &freq, w);
            fprintf(stderr, "(%f [cpd]) done\n", freq);
          }

        w->freqs[i] = freq;
      }
  }

  /* write binary file */
  fprintf(stderr, "main: writing %s...", output_file);
  invert_tmode_write_binary(output_file, w);
  fprintf(stderr, "done\n");

  free(nmodes);
  invert_tmode_free(w);

  return 0;
}
