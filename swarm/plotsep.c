/*
 * plotsep.c
 *
 * Compute time series of longitude and local time separation
 * between 2 Swarm satellites
 *
 * Usage:
 *
 * ./plotsep <-a swarm_idx_file1> <-b swarm_idx_file2>
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

#include <satdata/satdata.h>

#include <common/common.h>
#include <msynth/msynth.h>
#include <track/track.h>

/* find track closest in time to t */
size_t
find_idx(const double t, track_workspace *track_p)
{
  size_t i;
  size_t idx = 0;
  double dt_min = 1.0e12;

  for (i = 0; i < track_p->n; ++i)
    {
      track_data *tptr = &(track_p->tracks[i]);
      double dt = t - tptr->t_eq;

      if (fabs(dt) < dt_min)
        {
          idx = i;
          dt_min = fabs(dt);
        }
    }

  return idx;
}

int
find_sep(const char *filename, track_workspace *track1, track_workspace *track2)
{
  FILE *fp;
  size_t i;
  double alpha = 1.0e-3;
  double dlon_avg = 0.0;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "find_sep: cannot open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
  fprintf(fp, "# Field %zu: timestamp (decimal year)\n", i++);
  fprintf(fp, "# Field %zu: time difference between satellites crossing the equator (minutes)\n", i++);
  fprintf(fp, "# Field %zu: local time difference (hours)\n", i++);
  fprintf(fp, "# Field %zu: longitude difference (degrees)\n", i++);
  fprintf(fp, "# Field %zu: longitude difference moving average (degrees)\n", i++);

  for (i = 0; i < track1->n; ++i)
    {
      track_data *tptr1 = &(track1->tracks[i]);
      size_t j = find_idx(tptr1->t_eq, track2);
      track_data *tptr2 = &(track2->tracks[j]);
      double dt, dlt, dlon;

      /* find time difference in minutes */
      dt = tptr1->t_eq - tptr2->t_eq;
      dt /= 1000.0 * 60.0;

      if (fabs(dt) > 20.0)
        continue;

      if (tptr1->satdir != tptr2->satdir)
        continue;

      dlt = tptr1->lt_eq - tptr2->lt_eq;
      dlon = wrap180(tptr1->lon_eq - tptr2->lon_eq);

      if (dlon_avg == 0.0)
        dlon_avg = dlon;
      else
        dlon_avg = alpha * dlon + (1.0 - alpha) * dlon_avg;

      fprintf(fp, "%ld %f %f %f %f %f\n",
              satdata_epoch2timet(tptr1->t_eq),
              satdata_epoch2year(tptr1->t_eq),
              dt,
              dlt,
              dlon,
              dlon_avg);
    }

  fclose(fp);

  return 0;
}

int
main(int argc, char *argv[])
{
  satdata_mag *data1 = NULL;
  satdata_mag *data2 = NULL;
  struct timeval tv0, tv1;
  int c;
  char *outfile = "data.dat";
  track_workspace *track1, *track2;

  while ((c = getopt(argc, argv, "a:b:")) != (-1))
    {
      switch (c)
        {
          case 'a':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);
            data1 = satdata_swarm_read_idx(optarg, 0);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu records read, %g seconds)\n", data1->n,
                    time_diff(tv0, tv1));
            break;

          case 'b':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);
            data2 = satdata_swarm_read_idx(optarg, 0);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu records read, %g seconds)\n", data2->n,
                    time_diff(tv0, tv1));

          default:
            break;
        }
    }

  if (!data1 || !data2)
    {
      fprintf(stderr, "Usage: %s <-a swarm_idx_file1> <-b swarm_idx_file2>\n",
              argv[0]);
      exit(1);
    }

  track1 = track_alloc();
  track2 = track_alloc();

  fprintf(stderr, "main: separating tracks for satellite 1...");
  track_init(data1, NULL, track1);
  fprintf(stderr, "done (%zu tracks found)\n", track1->n);

  fprintf(stderr, "main: separating tracks for satellite 2...");
  track_init(data2, NULL, track2);
  fprintf(stderr, "done (%zu tracks found)\n", track2->n);

  find_sep(outfile, track1, track2);

  fprintf(stderr, "main: output written to %s\n", outfile);

  satdata_mag_free(data1);
  satdata_mag_free(data2);
  track_free(track1);
  track_free(track2);

  return 0;
}
