/*
 * emm_write_core.c
 *
 * Write EMM core header files, using as input a series
 * of msynth formatted files
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>

#include <gsl/gsl_math.h>

#include <msynth/msynth.h>

#define EMM_NMAX      15

/* external field parameters */
#define START_YEAR    2000
#define END_YEAR      2017
#define NYEAR         (END_YEAR - START_YEAR + 1)

int
main(int argc, char *argv[])
{
  char *outdir = "EMM2017";
  msynth_workspace *msynth_p = NULL;
  size_t nmodels = 0;
  int c;

  while ((c = getopt(argc, argv, "o:")) != (-1))
    {
      switch (c)
        {
          case 'o':
            outdir = optarg;
            break;

          default:
            break;
        }
    }

  if (optind == argc)
    {
      fprintf(stderr, "main: error: no models specified\n");
      fprintf(stderr, "Usage: %s [-o output_dir] msynth1.txt msynth2.txt ...\n",
              argv[0]);
      exit(1);
    }

  nmodels = argc - optind;

  fprintf(stderr, "main: allocating workspace for %zu models...",
          nmodels);
  msynth_p = msynth_alloc(EMM_NMAX, nmodels, NULL);
  fprintf(stderr, "done\n");

  for (c = optind; c < argc; ++c)
    {
      fprintf(stderr, "main: reading %s...", argv[c]);

      msynth_p = msynth_read2(argv[c], msynth_p);

      fprintf(stderr, "done\n");
    }

  /* recompute SV as differences of main field snapshots */
  fprintf(stderr, "main: computing SV coefficients...");
  msynth_calc_sv(msynth_p);
  fprintf(stderr, "done\n");

  {
    size_t i;
    int year;
    char filename[2048];

    for (i = 0; i < msynth_p->n_epochs; ++i)
      {
        year = START_YEAR + i;

        sprintf(filename, "%s/EMM%4d.COF", outdir, year);
        fprintf(stderr, "main: writing %s...", filename);
        msynth_emm_write(filename, (double) year, msynth_p);
        fprintf(stderr, "done\n");

        sprintf(filename, "%s/EMM%4dSV.COF", outdir, year);
        fprintf(stderr, "main: writing %s...", filename);
        msynth_emm_write_sv(filename, (double) year, msynth_p);
        fprintf(stderr, "done\n");
      }
  }

  return 0;
}
