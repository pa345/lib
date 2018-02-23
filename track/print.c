/*
 * print.c
 *
 * Visualize satellite data track by track.
 *
 * Pre-processing steps are:
 * 1. Instrument flags (recommended CHAMP flags except 1 star camera allowed)
 * 2. Track rms test
 * 3. Downsample by factor 15
 * 4. Compute and store along-track gradients
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <complex.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <satdata/satdata.h>
#include <indices/indices.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_statistics.h>

#include <common/bin3d.h>
#include <common/common.h>

#include "track.h"

/* define this to allow only 1 camera for data selection */
#define PRINT_ONE_CAMERA           1

/* number of seconds for computing along-track differences */
#define POLTOR_GRAD_DT             (40.0)

typedef struct
{
  int all;            /* print all tracks */
  double lt_min;      /* local time minimum (hours) */
  double lt_max;      /* local time maximum (hours) */
  double alt_min;     /* altitude minimum (km) */
  double alt_max;     /* altitude maximum (km) */
  double qd_min;      /* minimum required QD latitude (deg) */
  double qd_max;      /* maximum required QD latitude (deg) */
  double lon_min;     /* minimum longitude (deg) */
  double lon_max;     /* maximum longitude (deg) */
  double kp_min;      /* minimum kp */
  double kp_max;      /* maximum kp */
  double IMF_Bz_min;  /* minimum IMF B_z (nT) */
  double IMF_Bz_max;  /* maximum IMF B_z (nT) */
  size_t downsample;  /* downsampling factor */
  double alpha;       /* smoothing factor for high latitudes */
  double thresh[4];   /* rms thresholds */
} preprocess_parameters;

satdata_mag *
read_swarm(const char *filename)
{
  size_t nflag;
  satdata_mag *data;
  struct timeval tv0, tv1;

  fprintf(stderr, "read_swarm: reading %s...", optarg);
  gettimeofday(&tv0, NULL);
  data = satdata_swarm_read_idx(optarg, 0);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu data read, %g seconds)\n",
          data->n, time_diff(tv0, tv1));

  /* check for instrument flags since we use Stage1 data */

  fprintf(stderr, "read_swarm: filtering for instrument flags...");
  nflag = satdata_swarm_filter_instrument(1, data);
  fprintf(stderr, "done (%zu/%zu (%.1f%%) data flagged)\n",
          nflag, data->n, (double)nflag / (double)data->n * 100.0);

  return data;
} /* read_swarm() */

int
print_track_stats(const satdata_mag *data, const track_workspace *track_p)
{
  size_t nflagged = satdata_nflagged(data);
  size_t nleft = data->n - nflagged;
  size_t nflagged_track = track_nflagged(track_p);
  size_t nleft_track = track_p->n - nflagged_track;

  fprintf(stderr, "preprocess_data: total flagged data: %zu/%zu (%.1f%%)\n",
          nflagged, data->n, (double)nflagged / (double)data->n * 100.0);
  fprintf(stderr, "preprocess_data: total remaining data: %zu/%zu (%.1f%%)\n",
          nleft, data->n, (double)nleft / (double)data->n * 100.0);

  fprintf(stderr, "preprocess_data: total flagged tracks: %zu/%zu (%.1f%%)\n",
          nflagged_track, track_p->n, (double)nflagged_track / (double)track_p->n * 100.0);
  fprintf(stderr, "preprocess_data: total remaining tracks: %zu/%zu (%.1f%%)\n",
          nleft_track, track_p->n, (double)nleft_track / (double)track_p->n * 100.0);

  return 0;
} /* print_track_stats() */

/*
preprocess_data()

Inputs: params - preprocess parameters
          lt_min     - minimum local time
          lt_max     - maximum local time
          downsample - downsampling factor
        data   - satellite data

Return: pointer to sorted track workspace (should be freed by caller)
*/

track_workspace *
preprocess_data(const preprocess_parameters *params, satdata_mag *data)
{
  struct timeval tv0, tv1;
  track_workspace *track_p = track_alloc();

  fprintf(stderr, "preprocess_data: initializing tracks...");
  gettimeofday(&tv0, NULL);
  track_init(data, NULL, track_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  if (params->all)
    {
      print_track_stats(data, track_p);
      return track_p; /* no further filtering */
    }

  {
    const char *rmsfile = "satrms.dat";
    size_t nrms;

    nrms = track_flag_rms(rmsfile, params->thresh, NULL, data, track_p);
    fprintf(stderr, "preprocess_data: flagged (%zu/%zu) (%.1f%%) tracks due to high rms\n",
            nrms, track_p->n, (double) nrms / (double) track_p->n * 100.0);
  }

  /* flag local time */
  if (params->lt_min >= 0.0 && params->lt_max >= 0.0)
    {
      size_t nlt_data;
      size_t nlt = track_flag_lt(params->lt_min, params->lt_max, &nlt_data, data, track_p);

      fprintf(stderr, "preprocess_data: flagged data outside LT window [%g,%g]: %zu/%zu (%.1f%%) tracks flagged)\n",
              params->lt_min, params->lt_max,
              nlt, track_p->n, (double)nlt / (double)track_p->n * 100.0);
    }

  /* flag altitude */
  {
    size_t nalt = track_flag_meanalt(params->alt_min, params->alt_max, data, track_p);

    fprintf(stderr, "preprocess_data: flagged data due to altitude: %zu/%zu (%.1f%%) data flagged)\n",
            nalt, data->n, (double)nalt / (double)data->n * 100.0);
  }

  /* flag high kp data */
  {
    size_t nkp = track_flag_kp(params->kp_min, params->kp_max, data, track_p);

    fprintf(stderr, "preprocess_data: flagged data outside kp window [%g,%g]: %zu/%zu (%.1f%%) tracks flagged)\n",
            params->kp_min, params->kp_max,
            nkp, track_p->n, (double)nkp / (double)track_p->n * 100.0);
  }

  /* flag data according to IMF B */
  {
    size_t nIMF = track_flag_IMF(params->IMF_Bz_min, params->IMF_Bz_max, data, track_p);

    fprintf(stderr, "preprocess_data: flagged data outside IMF B_z window [%g,%g]: %zu/%zu (%.1f%%) tracks flagged)\n",
            params->IMF_Bz_min, params->IMF_Bz_max,
            nIMF, track_p->n, (double)nIMF / (double)track_p->n * 100.0);
  }

  /* flag longitude */
  {
    size_t nlon = track_flag_lon(params->lon_min, params->lon_max, NULL, data, track_p);

    fprintf(stderr, "preprocess_data: flagged data due to longitude: %zu/%zu (%.1f%%) tracks flagged)\n",
            nlon, track_p->n, (double)nlon / (double)track_p->n * 100.0);
  }

  /* print track statistics */
  {
    char *jstat_file = "track_jump_stats.dat";

    fprintf(stderr, "preprocess_data: printing jump track statistics to %s...", jstat_file);
    track_print_stats_flag(jstat_file, TRACK_FLG_JUMP, track_p);
    fprintf(stderr, "done\n");
  }

  /* downsample data */
  {
    size_t i;

    fprintf(stderr, "preprocess_data: downsampling data by factor %zu...", params->downsample);

    for (i = 0; i < data->n; ++i)
      {
        if (i % params->downsample != 0)
          data->flags[i] |= SATDATA_FLG_DOWNSAMPLE;
      }

    fprintf(stderr, "done\n");
  }

  if (params->alpha > 0.0)
    {
      fprintf(stderr, "preprocess_data: smoothing high-latitude residuals...");
      track_smooth(params->alpha, data, track_p);
      fprintf(stderr, "done\n");
    }

#if 0
  fprintf(stderr, "preprocess_data: removing magnetospheric offsets...");
  gettimeofday(&tv0, NULL);
  track_fix_offsets(data, track_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));
#endif

  print_track_stats(data, track_p);

  return track_p;
} /* preprocess_data() */

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options]\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --swarm_file | -s swarm_index_file          - Swarm index file\n");
  fprintf(stderr, "\t --champ_file | -c champ_index_file          - CHAMP index file\n");
  fprintf(stderr, "\t --champ_plp_file | -b champ_plp_index_file  - CHAMP PLP index file\n");
  fprintf(stderr, "\t --dmsp_file | -D dmsp_index_file            - DMSP index file\n");
  fprintf(stderr, "\t --all | -a                                  - print all tracks (no filtering)\n");
  fprintf(stderr, "\t --downsample | -d downsample                - downsampling factor\n");
  fprintf(stderr, "\t --output_file | -o output_file              - output file\n");
  fprintf(stderr, "\t --lt_min | -j lt_min                        - local time minimum\n");
  fprintf(stderr, "\t --lt_max | -k lt_max                        - local time maximum\n");
  fprintf(stderr, "\t --alt_min | -l alt_min                      - altitude minimum\n");
  fprintf(stderr, "\t --alt_max | -m alt_max                      - altitude maximum\n");
  fprintf(stderr, "\t --lon_min | -t lon_min                      - longitude minimum\n");
  fprintf(stderr, "\t --lon_max | -u lon_max                      - longitude maximum\n");
  fprintf(stderr, "\t --kp_min | -v kp_min                        - kp minimum\n");
  fprintf(stderr, "\t --kp_max | -w kp_max                        - kp maximum\n");
  fprintf(stderr, "\t --imf_bz_min | -A imf_bz_min                - minimum IMF B_z (nT)\n");
  fprintf(stderr, "\t --imf_bz_max | -B imf_bz_max                - maximum IMF B_z (nT)\n");
  fprintf(stderr, "\t --alpha | -q alpha                          - smoothing factor for high latitudes\n");
}

int
main(int argc, char *argv[])
{
  char *data_file = "track_data.dat";
  char *stats_file = "track_stats.dat";
  satdata_mag *data = NULL;
  satdata_lp *lp_data = NULL;
  struct timeval tv0, tv1;
  track_workspace *track_p;
  preprocess_parameters params;

  /* defaults */
  params.all = 0;
  params.lt_min = -1.0;
  params.lt_max = -1.0;
  params.alt_min = 0.0;
  params.alt_max = 1000.0;
  params.qd_min = -30.0;
  params.qd_max = 30.0;
  params.lon_min = -200.0;
  params.lon_max = 200.0;
  params.kp_min = 0.0;
  params.kp_max = 20.0;
  params.IMF_Bz_min = -1.0e6;
  params.IMF_Bz_max = 1.0e6;
  params.downsample = 15;
  params.alpha = -1.0;
  params.thresh[0] = 80.0;
  params.thresh[1] = 60.0;
  params.thresh[2] = 30.0;
  params.thresh[3] = 50.0;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "all", no_argument, NULL, 'a' },
          { "champ_file", required_argument, NULL, 'c' },
          { "swarm_file", required_argument, NULL, 's' },
          { "champ_plp_file", required_argument, NULL, 'b' },
          { "dmsp_file", required_argument, NULL, 'D' },
          { "downsample", required_argument, NULL, 'd' },
          { "lt_min", required_argument, NULL, 'j' },
          { "lt_max", required_argument, NULL, 'k' },
          { "alt_min", required_argument, NULL, 'l' },
          { "alt_max", required_argument, NULL, 'm' },
          { "lon_min", required_argument, NULL, 't' },
          { "lon_max", required_argument, NULL, 'u' },
          { "kp_min", required_argument, NULL, 'v' },
          { "kp_max", required_argument, NULL, 'w' },
          { "output_file", required_argument, NULL, 'o' },
          { "alpha", required_argument, NULL, 'q' },
          { "imf_bz_min", required_argument, NULL, 'A' },
          { "imf_bz_max", required_argument, NULL, 'B' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "ab:c:d:D:j:k:l:m:o:q:s:t:u:A:B:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'a':
            params.all = 1;
            break;

          case 's':
            data = read_swarm(optarg);
            break;

          case 'c':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);
            data = satdata_champ_read_idx(optarg, 1);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu data read, %g seconds)\n",
                    data->n, time_diff(tv0, tv1));

            /* check for instrument flags since we use Stage1 data */
            {
              size_t nflag;
              size_t champ_flags = 0;

#if PRINT_ONE_CAMERA
              /* allow only 1 camera in data selection */
              champ_flags = SATDATA_FLG_ONESC;
#endif

              fprintf(stderr, "main: filtering for instrument flags...");
              nflag = satdata_champ_filter_instrument(1, champ_flags, data);
              fprintf(stderr, "done (%zu/%zu (%.1f%%) data flagged)\n",
                      nflag, data->n, (double)nflag / (double)data->n * 100.0);
            }

            break;

          case 'D':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);
            data = satdata_dmsp_read_idx(optarg, 0);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu data read, %g seconds)\n",
                    data->n, time_diff(tv0, tv1));
            break;

          case 'b':
            fprintf(stderr, "main: reading CHAMP PLP data from %s...", optarg);
            lp_data = satdata_champ_plp_read_idx(optarg, 0);
            fprintf(stderr, "done\n");
            break;

          case 'd':
            params.downsample = (size_t) atoi(optarg);
            break;

          case 'o':
            data_file = optarg;
            break;

          case 'j':
            params.lt_min = atof(optarg);
            break;

          case 'k':
            params.lt_max = atof(optarg);
            break;

          case 'l':
            params.alt_min = atof(optarg);
            break;

          case 'm':
            params.alt_max = atof(optarg);
            break;

          case 't':
            params.lon_min = atof(optarg);
            break;

          case 'u':
            params.lon_max = atof(optarg);
            break;

          case 'v':
            params.kp_min = atof(optarg);
            break;

          case 'w':
            params.kp_max = atof(optarg);
            break;

          case 'q':
            params.alpha = atof(optarg);
            break;

          case 'A':
            params.IMF_Bz_min = atof(optarg);
            break;

          case 'B':
            params.IMF_Bz_max = atof(optarg);
            break;

          default:
            break;
        }
    }

  if (!data)
    {
      print_help(argv);
      exit(1);
    }

  fprintf(stderr, "main: LT minimum       = %.1f\n", params.lt_min);
  fprintf(stderr, "main: LT maximum       = %.1f\n", params.lt_max);
  fprintf(stderr, "main: altitude minimum = %.1f\n", params.alt_min);
  fprintf(stderr, "main: altitude maximum = %.1f\n", params.alt_max);
  fprintf(stderr, "main: QD minimum       = %.1f\n", params.qd_min);
  fprintf(stderr, "main: QD maximum       = %.1f\n", params.qd_max);
  fprintf(stderr, "main: lon minimum      = %.1f\n", params.lon_min);
  fprintf(stderr, "main: lon maximum      = %.1f\n", params.lon_max);
  fprintf(stderr, "main: kp minimum       = %.1f\n", params.kp_min);
  fprintf(stderr, "main: kp maximum       = %.1f\n", params.kp_max);
  fprintf(stderr, "main: IMF B_z minimum  = %.1f\n", params.IMF_Bz_min);
  fprintf(stderr, "main: IMF B_z maximum  = %.1f\n", params.IMF_Bz_max);
  fprintf(stderr, "main: smoothing alpha  = %f\n", params.alpha);

  if (lp_data)
    {
      fprintf(stderr, "main: adding electron densities to magnetic data...");
      satdata_mag_fill_ne(data, lp_data);
      fprintf(stderr, "done\n");
    }

  /*XXX*/
  {
    fprintf(stderr, "main: computing dipole tilt along track...");
    gettimeofday(&tv0, NULL);
    track_synth_tilt(data);
    gettimeofday(&tv1, NULL);
    fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

    fprintf(stderr, "main: computing MLT along track...");
    gettimeofday(&tv0, NULL);
    track_synth_MLT(data);
    gettimeofday(&tv1, NULL);
    fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));
  }

  track_p = preprocess_data(&params, data);

  fprintf(stderr, "main: printing data to %s...", data_file);
  track_print(data_file, 0, data, track_p);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: printing statistics to %s...", stats_file);
  track_print_stats_flag(stats_file, 0, track_p);
  fprintf(stderr, "done\n");

  satdata_mag_free(data);
  track_free(track_p);

  return 0;
} /* main() */
