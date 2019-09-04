/*
 * plot_tracks.c
 *
 * Plot information about individual Swarm tracks, including local
 * time of equator crossing
 *
 * ./plot_tracks <-i swarm_index_file> [-o output_file]
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

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_indices.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_track.h>

int
main(int argc, char *argv[])
{
  satdata_mag *data;
  struct timeval tv0, tv1;
  int c;
  char *infile = NULL;
  char *outfile = "swarm_rms.dat";

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
      fprintf(stderr, "Usage: %s <-i swarm_index_file> [-o output_file]\n",
              argv[0]);
      exit(1);
    }

  fprintf(stderr, "input file = %s\n", infile);
  fprintf(stderr, "output file = %s\n", outfile);

  fprintf(stderr, "Reading %s...", infile);
  gettimeofday(&tv0, NULL);

  data = satdata_swarm_read_idx(infile, 0);
  if (!data)
    {
      fprintf(stderr, "main: error reading %s\n", infile);
      exit(1);
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu records read, %g seconds)\n", data->n,
          time_diff(tv0, tv1));

  {
    size_t nrms;
    const double thresh[] = { 1000.0, 1000.0, 1000.0, 1000.0 };
    track_workspace *track_p = track_alloc();

    track_init(data, NULL, track_p);

    nrms = track_filter_rms(outfile, thresh, data, track_p);
    fprintf(stderr, "main: flagged (%zu/%zu) (%.1f%%) points due to high rms\n",
            nrms, data->n, (double) nrms / (double) data->n * 100.0);
    track_free(track_p);
  }

  satdata_mag_free(data);

  return 0;
} /* main() */
