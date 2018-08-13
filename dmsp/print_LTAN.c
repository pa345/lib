/*
 * print_LTAN.c
 *
 * print LT of ascending node
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <gsl/gsl_math.h>

#include <satdata/satdata.h>
#include <indices/indices.h>
#include <common/common.h>

#include "track.h"

int
main(int argc, char *argv[])
{
  satdata_mag *data = NULL;
  track_workspace *track_p;
  struct timeval tv0, tv1;
  size_t i;
  int c;

  while ((c = getopt(argc, argv, "i:d:")) != (-1))
    {
      switch (c)
        {
          case 'i':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);
            data = satdata_dmsp_read_idx(optarg, 1);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu records read, %g seconds)\n", data->n,
                    time_diff(tv0, tv1));
            break;

          default:
            break;
        }
    }

  if (!data)
    {
      fprintf(stderr, "Usage: %s <-i dmsp_index_file>\n", argv[0]);
      exit(1);
    }

  track_p = track_alloc();

  fprintf(stderr, "main: initializing tracks...");
  gettimeofday(&tv0, NULL);
  track_init(data, NULL, track_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  i = 1;
  printf("# Field %zu: timestamp\n", i++);
  printf("# Field %zu: LTAN (hours)\n", i++);

  for (i = 0; i < track_p->n; ++i)
    {
      track_data *tptr = &(track_p->tracks[i]);

      if (tptr->satdir == -1)
        continue;

      printf("%ld %f\n",
             satdata_epoch2timet(tptr->t_eq),
             tptr->lt_eq);
    }

  track_free(track_p);
  satdata_mag_free(data);

  return 0;
}
