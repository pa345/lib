/*
 * print_tracks.c
 *
 * Visualize satellite data track by track.
 *
 * Pre-processing steps are:
 * 1. Instrument flags
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

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_statistics.h>

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_indices.h>
#include <mainlib/ml_bin3d.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_msynth.h>
#include <mainlib/ml_track.h>

/* define this to allow only 1 camera for data selection */
#define POLTOR_ONE_CAMERA          0

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
      size_t nlt = track_flag_lt(params->lt_min, params->lt_max, NULL, data, track_p);

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

    fprintf(stderr, "preprocess_data: flagged data outside kp window [%g,%g]: %zu/%zu (%.1f%%) data flagged)\n",
            params->kp_min, params->kp_max,
            nkp, data->n, (double)nkp / (double)data->n * 100.0);
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

int
print_data(const char *filename, const size_t flags,
           const satdata_mag *data, track_workspace *w)
{
  int s = 0;
  size_t i, j;
  FILE *fp;
  size_t nflagged;
  /*msynth_workspace *msynth_p = msynth_read(MSYNTH_BOUMME_FILE);*/
  msynth_workspace *msynth_p = msynth_swarm_read(MSYNTH_CHAOS_FILE);

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_data: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  msynth_set(1, 15, msynth_p);

  nflagged = track_nflagged(w);
  fprintf(fp, "# Total tracks:    %zu\n", w->n);
  fprintf(fp, "# Total flagged:   %zu\n", nflagged);
  fprintf(fp, "# Total unflagged: %zu\n", w->n - nflagged);
  fprintf(fp, "# Track selection flags: %zu\n", flags);

  i = 1;
  fprintf(fp, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
  fprintf(fp, "# Field %zu: UT (hours)\n", i++);
  fprintf(fp, "# Field %zu: local time (hours)\n", i++);
  fprintf(fp, "# Field %zu: season (doy)\n", i++);
  fprintf(fp, "# Field %zu: altitude (km)\n", i++);
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: geocentric latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: QD latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: NEC X residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: NEC Y residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: NEC Z residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: F residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: NEC X measurement (nT)\n", i++);
  fprintf(fp, "# Field %zu: NEC Y measurement (nT)\n", i++);
  fprintf(fp, "# Field %zu: NEC Z measurement (nT)\n", i++);
  fprintf(fp, "# Field %zu: NEC X main field (nT)\n", i++);
  fprintf(fp, "# Field %zu: NEC Y main field (nT)\n", i++);
  fprintf(fp, "# Field %zu: NEC Z main field (nT)\n", i++);
  fprintf(fp, "# Field %zu: NEC X crustal field (nT)\n", i++);
  fprintf(fp, "# Field %zu: NEC Y crustal field (nT)\n", i++);
  fprintf(fp, "# Field %zu: NEC Z crustal field (nT)\n", i++);
  fprintf(fp, "# Field %zu: NEC X external field (nT)\n", i++);
  fprintf(fp, "# Field %zu: NEC Y external field (nT)\n", i++);
  fprintf(fp, "# Field %zu: NEC Z external field (nT)\n", i++);
  fprintf(fp, "# Field %zu: electron density data (cm^{-3})\n", i++);

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);

      if (flags != 0xffffffff)
        {
          if (flags == 0 && tptr->flags != 0)
            continue;

          if (flags != 0 && !(tptr->flags & flags))
            continue;
        }

      for (j = 0; j < tptr->n; ++j)
        {
          size_t didx = j + tptr->start_idx;
          time_t unix_time = satdata_epoch2timet(data->t[didx]);
          double ut = get_ut(unix_time);
          double lt = get_localtime(unix_time, data->longitude[didx] * M_PI / 180.0);
          double B_main[4];

          if (SATDATA_BadData(data->flags[didx]))
            continue;

#if 0
          B_main[0] = SATDATA_VEC_X(data->B_main, didx),
          B_main[1] = SATDATA_VEC_Y(data->B_main, didx),
          B_main[2] = SATDATA_VEC_Z(data->B_main, didx),
#else
          {
            double r = data->r[didx];
            double theta = M_PI / 2.0 - data->latitude[didx] * M_PI / 180.0;
            double phi = data->longitude[didx] * M_PI / 180.0;

            msynth_eval(satdata_epoch2year(data->t[didx]), r, theta, phi, B_main, msynth_p);
          }
#endif

          fprintf(fp, "%ld %.2f %.2f %.1f %.2f %.3f %.3f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %.5e\n",
                  unix_time,
                  ut,
                  lt,
                  get_season(unix_time),
                  data->r[didx] - R_EARTH_KM,
                  data->longitude[didx],
                  data->latitude[didx],
                  data->qdlat[didx],
                  tptr->Bx[j],
                  tptr->By[j],
                  tptr->Bz[j],
                  tptr->Bf[j],
                  SATDATA_VEC_X(data->B, didx),
                  SATDATA_VEC_Y(data->B, didx),
                  SATDATA_VEC_Z(data->B, didx),
                  B_main[0],
                  B_main[1],
                  B_main[2],
                  SATDATA_VEC_X(data->B_crust, didx),
                  SATDATA_VEC_Y(data->B_crust, didx),
                  SATDATA_VEC_Z(data->B_crust, didx),
                  SATDATA_VEC_X(data->B_ext, didx),
                  SATDATA_VEC_Y(data->B_ext, didx),
                  SATDATA_VEC_Z(data->B_ext, didx),
                  data->ne[didx]);
        }

      fprintf(fp, "\n\n");
    }

  fclose(fp);

  msynth_free(msynth_p);

  return s;
}

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options]\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --swarm_file | -s swarm_index_file          - Swarm index file\n");
  fprintf(stderr, "\t --champ_file | -c champ_index_file          - CHAMP index file\n");
  fprintf(stderr, "\t --champ_plp_file | -b champ_plp_index_file  - CHAMP PLP index file\n");
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
  fprintf(stderr, "\t --alpha | -q alpha                          - smoothing factor for high latitudes\n");
}

int
main(int argc, char *argv[])
{
  char *data_file = "track_data.dat";
  char *stats_file = "track_stats.dat";
  char *output_file = NULL;
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
  params.downsample = 20;
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
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "ab:c:d:j:k:l:m:o:q:s:t:u:", long_options, &option_index);
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
            data = satdata_champ_read_idx(optarg, 0);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu data read, %g seconds)\n",
                    data->n, time_diff(tv0, tv1));

            /* check for instrument flags since we use Stage1 data */
            {
              size_t nflag;
              size_t champ_flags = 0;

#if POLTOR_ONE_CAMERA
              /* allow only 1 camera in data selection */
              champ_flags = SATDATA_FLG_ONESC;
#endif

              fprintf(stderr, "main: filtering for instrument flags...");
              nflag = satdata_champ_filter_instrument(1, champ_flags, data);
              fprintf(stderr, "done (%zu/%zu (%.1f%%) data flagged)\n",
                      nflag, data->n, (double)nflag / (double)data->n * 100.0);
            }

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
            output_file = optarg;
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
  fprintf(stderr, "main: smoothing alpha  = %f\n", params.alpha);

  if (lp_data)
    {
      fprintf(stderr, "main: adding electron densities to magnetic data...");
      satdata_mag_fill_ne(data, lp_data);
      fprintf(stderr, "done\n");
    }

  track_p = preprocess_data(&params, data);

  fprintf(stderr, "main: printing data to %s...", data_file);
  print_data(data_file, 0, data, track_p);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: printing statistics to %s...", stats_file);
  track_print_stats_flag(stats_file, 0, track_p);
  fprintf(stderr, "done\n");

  satdata_mag_free(data);
  track_free(track_p);

  return 0;
} /* main() */
