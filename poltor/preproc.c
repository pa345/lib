/*
 * preproc.c
 *
 * Pre-process satellite data and save in magdata format. This
 * program should be called with Stage1 data, in order to do track
 * rms test, and along-track gradient calculation
 *
 * Pre-processing steps are:
 * 1. Instrument flags
 * 2. Track rms test
 * 3. Downsample by factor 15
 * 4. UT selection (near 11 UT)
 * 5. kp <= 2
 * 6. Compute and store along-track gradients
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <complex.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <flow/flow.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_statistics.h>

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_indices.h>
#include <mainlib/ml_bin3d.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_track.h>
#include <mainlib/ml_euler.h>
#include <mainlib/ml_magdata.h>

#include "poltor.h"

/* define this to allow only 1 camera for data selection */
#define POLTOR_ONE_CAMERA          0

/* number of seconds for computing along-track differences */
#define POLTOR_GRAD_DT             (40.0)

#define MIN_KP                     (0.0)
#define MAX_KP                     (4.0)

typedef struct
{
  int all;            /* print all tracks */
  double lt_min;      /* local time minimum (hours) */
  double lt_max;      /* local time maximum (hours) */
  double alt_min;     /* altitude minimum (km) */
  double alt_max;     /* altitude maximum (km) */
  double qd_min;      /* minimum required QD latitude (deg) */
  double qd_max;      /* maximum required QD latitude (deg) */
  double season_min;  /* season minimum [0,366] */
  double season_max;  /* season maximum [0,366] */
  double season_min2; /* season minimum [0,366] */
  double season_max2; /* season maximum [0,366] */
  double jump_thresh; /* jump threshold (nT) */
  double ut;          /* UT (hours) */
  size_t downsample;  /* downsampling factor */
  double smooth_alpha; /* smooth factor for high latitudes */
} preprocess_parameters;

int
replace_longitude(const double ut, satdata_mag *data)
{
  int s = 0;
  size_t i;

  for (i = 0; i < data->n; ++i)
    {
      time_t unix_time = satdata_epoch2timet(data->t[i]);
      double lt = get_localtime(unix_time, data->longitude[i] * M_PI / 180.0);

      data->longitude[i] = wrap180(15.0 * (lt - ut));
    }

  return s;
}

size_t
copy_track(const size_t track_idx, const satdata_mag *data,
           const track_workspace *track_p, magdata *mdata)
{
  int s = 0;
  size_t i;
  track_data *tptr = &(track_p->tracks[track_idx]);
  const size_t start_idx = tptr->start_idx;
  const size_t end_idx = tptr->end_idx;
  const size_t grad_idx = (size_t) POLTOR_GRAD_DT;
  const double grad_dt_min = POLTOR_GRAD_DT - 2.0;
  const double grad_dt_max = POLTOR_GRAD_DT + 2.0;
  magdata_datum datum;
  size_t ngrad = 0;
  int flagged_start = 0;

  for (i = start_idx; i <= end_idx; ++i)
    {
      size_t flags = MAGDATA_FLG_F;

      /* ignore bad data */
      if (SATDATA_BadData(data->flags[i]))
        continue;

      if (data->flags[i] & SATDATA_FLG_DOWNSAMPLE)
        continue;

      /* initialize to 0 */
      magdata_datum_init(&datum);

      if (!flagged_start)
        {
          flags |= MAGDATA_FLG_TRACK_START;
          flagged_start = 1;
        }

      /*
       * here there is at minimum a scalar measurement available; check
       * for vector measurement
       */
      if (!(data->flags[i] & SATDATA_FLG_NOVEC_NEC))
        {
          flags |= MAGDATA_FLG_X | MAGDATA_FLG_Y | MAGDATA_FLG_Z;

          datum.B_nec[0] = SATDATA_VEC_X(data->B, i);
          datum.B_nec[1] = SATDATA_VEC_Y(data->B, i);
          datum.B_nec[2] = SATDATA_VEC_Z(data->B, i);
          datum.B_vfm[0] = SATDATA_VEC_X(data->B_VFM, i);
          datum.B_vfm[1] = SATDATA_VEC_Y(data->B_VFM, i);
          datum.B_vfm[2] = SATDATA_VEC_Z(data->B_VFM, i);
        }

      datum.F = data->F[i];

      datum.B_model[0] = SATDATA_VEC_X(data->B_main, i) +
                         SATDATA_VEC_X(data->B_crust, i) +
                         SATDATA_VEC_X(data->B_ext, i);
      datum.B_model[1] = SATDATA_VEC_Y(data->B_main, i) +
                         SATDATA_VEC_Y(data->B_crust, i) +
                         SATDATA_VEC_Y(data->B_ext, i);
      datum.B_model[2] = SATDATA_VEC_Z(data->B_main, i) +
                         SATDATA_VEC_Z(data->B_crust, i) +
                         SATDATA_VEC_Z(data->B_ext, i);

      {
        size_t j = GSL_MIN(i + grad_idx, data->n - 1);
        double dt = (data->t[j] - data->t[i]) / 1000.0; /* in s */

        /* check if along-track point should be rejected */
        if (!SATDATA_BadData(data->flags[j]) &&
            (dt >= grad_dt_min && dt <= grad_dt_max))
          {
            /* set flag to indicate scalar gradient information available */
            flags |= MAGDATA_FLG_DF_NS;

            /* check for vector measurement at along-track point */
            if (MAGDATA_ExistVector(flags) && !(data->flags[j] & SATDATA_FLG_NOVEC_NEC))
              {
                assert(data->flags[j] == 0 || data->flags[j] == SATDATA_FLG_DOWNSAMPLE);
                flags |= MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS;

                datum.B_nec_ns[0] = SATDATA_VEC_X(data->B, j);
                datum.B_nec_ns[1] = SATDATA_VEC_Y(data->B, j);
                datum.B_nec_ns[2] = SATDATA_VEC_Z(data->B, j);
              }

            datum.F_ns = data->F[j];

            datum.B_model_ns[0] = SATDATA_VEC_X(data->B_main, j) +
                                  SATDATA_VEC_X(data->B_crust, j) +
                                  SATDATA_VEC_X(data->B_ext, j);
            datum.B_model_ns[1] = SATDATA_VEC_Y(data->B_main, j) +
                                  SATDATA_VEC_Y(data->B_crust, j) +
                                  SATDATA_VEC_Y(data->B_ext, j);
            datum.B_model_ns[2] = SATDATA_VEC_Z(data->B_main, j) +
                                  SATDATA_VEC_Z(data->B_crust, j) +
                                  SATDATA_VEC_Z(data->B_ext, j);

            datum.r_ns = data->altitude[j] + data->R;
            datum.theta_ns = M_PI / 2.0 - data->latitude[j] * M_PI / 180.0;
            datum.phi_ns = data->longitude[j] * M_PI / 180.0;
            datum.qdlat_ns = data->qdlat[j];
          }
      }

      datum.t = data->t[i];
      datum.r = data->altitude[i] + data->R;
      datum.theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
      datum.phi = data->longitude[i] * M_PI / 180.0;
      datum.qdlat = data->qdlat[i];
      datum.ne = data->ne[i];
      datum.satdir = satdata_mag_satdir(i, data);
      datum.flags = flags;

      s = magdata_add(&datum, mdata);
      if (s == 0 && flags & MAGDATA_FLG_DZ_NS)
        ++ngrad;
    }

  return ngrad;
} /* copy_track() */

magdata *
copy_data(const satdata_mag *data, const track_workspace *track_p)
{
  size_t ndata = 0;
  size_t ngrad = 0;
  magdata *mdata;
  size_t i;
  
  /* count number of remaining data */
  for (i = 0; i < data->n; ++i)
    {
      if (SATDATA_BadData(data->flags[i]) ||
          (data->flags[i] & SATDATA_FLG_DOWNSAMPLE))
        continue;

      ++ndata;
    }

  mdata = magdata_alloc(ndata, data->R);
  if (!mdata)
    return 0;

  for (i = 0; i < track_p->n; ++i)
    {
      track_data *tptr = &(track_p->tracks[i]);

      /* discard flagged tracks */
      if (tptr->flags)
        continue;

      ngrad += copy_track(i, data, track_p, mdata);
    }

  assert(ndata == mdata->n);

  fprintf(stderr, "\n");
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) along-track gradient data available\n",
          ngrad, ndata, (double)ngrad / (double)ndata * 100.0);

  return mdata;
}

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

satdata_mag *
read_champ(const char *filename)
{
  size_t nflag;
  satdata_mag *data;
  struct timeval tv0, tv1;
  size_t champ_flags = 0;

  fprintf(stderr, "read_champ: reading %s...", optarg);
  gettimeofday(&tv0, NULL);
  data = satdata_champ_read_idx(optarg, 0);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu data read, %g seconds)\n",
          data->n, time_diff(tv0, tv1));

  /* check for instrument flags since we use Stage1 data */

#if POLTOR_ONE_CAMERA
  /* allow only 1 camera in data selection */
  champ_flags = SATDATA_FLG_ONESC;
#endif

  fprintf(stderr, "read_champ: filtering for instrument flags...");
  nflag = satdata_champ_filter_instrument(1, champ_flags, data);
  fprintf(stderr, "done (%zu/%zu (%.1f%%) data flagged)\n",
          nflag, data->n, (double)nflag / (double)data->n * 100.0);

  return data;
}

int
callback_season(const double doy, const void *params)
{
  preprocess_parameters *p = (preprocess_parameters *) params;
  int s = -1;

  if (p->season_min < 0.0 && p->season_max < 0.0 &&
      p->season_min2 < 0.0 && p->season_max2 < 0.0)
    return 0; /* no seasonal selection, allow everything */

  if (p->season_min >= 0.0 && p->season_max >= 0.0)
    {
      if (doy >= p->season_min && doy <= p->season_max)
        return 0; /* allow */
    }

  if (p->season_min2 >= 0.0 && p->season_max2 >= 0.0)
    {
      if (doy >= p->season_min2 && doy <= p->season_max2)
        return 0; /* allow */
    }

  return s;
} /* callback_season() */

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
          ut         - UT in hours (set to -1 to take all UT)
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

  /* detect bad tracks with rms test */
  {
    const char *rmsfile = "satrms.dat";
    size_t nrms;
    double thresh[] = { 80.0, 60.0, 30.0, 50.0 };

    nrms = track_flag_rms(rmsfile, thresh, NULL, data, track_p);
    fprintf(stderr, "preprocess_data: flagged (%zu/%zu) (%.1f%%) tracks due to high rms\n",
            nrms, track_p->n, (double) nrms / (double) track_p->n * 100.0);
  }

  /*
   * filter data jumps; some CHAMP tracks in 2009 and 2010 have sudden jumps
   * larger than 5 nT (star camera issue possibly?)
   */
  {
    char *jump_file = "jumps.dat";
    size_t njump = track_flag_jumps(params->jump_thresh, data, track_p);

    fprintf(stderr, "preprocess_data: flagged (%zu/%zu) (%.1f%%) points due to jumps (jump threshold = %.1f nT)\n",
            njump, data->n, (double) njump / (double) data->n * 100.0, params->jump_thresh);

#if 0
    fprintf(stderr, "preprocess_data: writing jump tracks to %s...", jump_file);
    track_print(jump_file, TRACK_FLG_JUMP, data, track_p);
    fprintf(stderr, "done\n");
#endif
  }

  /* flag universal time */
  if (params->ut >= 0.0)
    {
      size_t nut;
      const double dut = 0.5;
      const double ut_min = params->ut - dut;
      const double ut_max = params->ut + dut;

      nut = track_flag_ut(ut_min, ut_max, data, track_p);

      fprintf(stderr, "preprocess_data: flagged data outside UT window [%g,%g]: %zu/%zu (%.1f%%) data flagged)\n",
              ut_min, ut_max,
              nut, data->n, (double)nut / (double)data->n * 100.0);
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

  /* flag incomplete tracks */
  {
    size_t ninc = track_flag_incomplete(params->qd_min, params->qd_max, data, track_p);

    fprintf(stderr, "preprocess_data: flagged data due to missing data: %zu/%zu (%.1f%%) data flagged)\n",
            ninc, data->n, (double)ninc / (double)data->n * 100.0);
  }

  /* flag season */
  {
    size_t nseas = track_flag_season(callback_season, (void *) params, data, track_p);

    fprintf(stderr, "preprocess_data: flagged data due to season: %zu/%zu (%.1f%%) data flagged)\n",
            nseas, data->n, (double)nseas / (double)data->n * 100.0);
  }

  /* flag high kp data */
  {
    const double kp_min = MIN_KP;
    const double kp_max = MAX_KP;
    size_t nkp = track_flag_kp(kp_min, kp_max, data, track_p);

    fprintf(stderr, "preprocess_data: flagged data outside kp window [%g,%g]: %zu/%zu (%.1f%%) data flagged)\n",
            kp_min, kp_max,
            nkp, data->n, (double)nkp / (double)data->n * 100.0);
  }

  /* last check: flag tracks with very few good data points left */
  {
    size_t nflag = track_flag_n(2000, data, track_p);

    fprintf(stderr, "preprocess_data: flagged data due to low data points: %zu/%zu (%.1f%%) data flagged)\n",
            nflag, data->n, (double)nflag / (double)data->n * 100.0);
  }

  /* print track statistics */
  {
    char *stat_file_all = "track_stats_all.dat";
    char *stat_file = "track_stats.dat";
    char *jstat_file = "track_jump_stats.dat";

    fprintf(stderr, "preprocess_data: printing all track statistics to %s...", stat_file_all);
    track_print_stats(stat_file_all, track_p);
    fprintf(stderr, "done\n");

    fprintf(stderr, "preprocess_data: printing track statistics to %s...", stat_file);
    track_print_stats_flag(stat_file, 0, track_p);
    fprintf(stderr, "done\n");

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

#if 0
  /* XXX This code is introducting NaN's into the CHAMP B_z vector data
   * on 17 september 2000 ???
   */

  /* smooth high-latitude data */
  if (params->smooth_alpha > 0.0)
    {
      fprintf(stderr, "preprocess_data: smoothing high-latitude data...");
      track_smooth(params->smooth_alpha, data, track_p);
      fprintf(stderr, "done\n");
    }
#endif

#if 0
  fprintf(stderr, "preprocess_data: removing magnetospheric offsets...");
  gettimeofday(&tv0, NULL);
  track_fix_offsets(data, track_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));
#endif

#if 0
  /*
   * XXX: track_filter cannot handle data gaps and sometimes inserts
   * nan's into the data stream
   */
  {
    char *filter_file = "filter.dat";

    fprintf(stderr, "preprocess_data: filtering tracks...");
    gettimeofday(&tv0, NULL);
    track_filter(filter_file, track_p);
    gettimeofday(&tv1, NULL);
    fprintf(stderr, "done (%g seconds, data written to %s)\n",
            time_diff(tv0, tv1), filter_file);
  }
#endif

  print_track_stats(data, track_p);

  return track_p;
} /* preprocess_data() */

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options]\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --swarm_file | -s swarm_index_file  - Swarm index file\n");
  fprintf(stderr, "\t --champ_file | -c champ_index_file  - CHAMP index file\n");
  fprintf(stderr, "\t --all | -a                          - print all tracks (no filtering)\n");
  fprintf(stderr, "\t --downsample | -d downsample        - downsampling factor\n");
  fprintf(stderr, "\t --universal_time | -t decimal UT    - Universal time (hours), or -1\n");
  fprintf(stderr, "\t --euler_file | -e euler_file        - Euler angles file\n");
  fprintf(stderr, "\t --euler_file2 | -f euler_file       - Euler angles file for satellite 2\n");
  fprintf(stderr, "\t --output_file | -o output_file      - output file\n");
  fprintf(stderr, "\t --lt_min | -f lt_min                - local time minimum\n");
  fprintf(stderr, "\t --lt_max | -b lt_max                - local time maximum\n");
  fprintf(stderr, "\t --alt_min | -l alt_min              - altitude minimum\n");
  fprintf(stderr, "\t --alt_max | -m alt_max              - altitude maximum\n");
  fprintf(stderr, "\t --season_min | -g season_min        - season minimum\n");
  fprintf(stderr, "\t --season_max | -h season_max        - season maximum\n");
  fprintf(stderr, "\t --season_min2 | -j season_min2      - season minimum 2\n");
  fprintf(stderr, "\t --season_max2 | -k season_max2      - season maximum 2\n");
  fprintf(stderr, "\t --jump_thresh | -u jump_thresh      - jump threshold (nT)\n");
  fprintf(stderr, "\t --smooth_alpha | -p alpha           - smooth factor for high-latitude filter\n");
}

int
main(int argc, char *argv[])
{
  char *datamap_file = "datamap.dat";
  char *data_prefix = "output";
  char *output_file = NULL;
  satdata_mag *data = NULL;
  magdata *mdata;
  euler_workspace *euler_p = NULL;
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
  params.jump_thresh = 10.0;
  params.ut = -1.0;
  params.downsample = 15;
  params.smooth_alpha = 1.0e-3;

#if 1
  /* solar max */
  params.season_min = 30.0;
  params.season_max = 100.0;
  params.season_min2 = 260.0;
  params.season_max2 = 340.0;
#else
  /* solar min */
  params.season_min = 30.0;
  params.season_max = 150.0;
  params.season_min2 = 244.0;
  params.season_max2 = 315.0;
#endif

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "swarm_file", required_argument, NULL, 's' },
          { "champ_file", required_argument, NULL, 'c' },
          { "all", no_argument, NULL, 'a' },
          { "downsample", required_argument, NULL, 'd' },
          { "output_file", required_argument, NULL, 'o' },
          { "universal_time", required_argument, NULL, 't' },
          { "euler_file", required_argument, NULL, 'e' },
          { "lt_min", required_argument, NULL, 'f' },
          { "lt_max", required_argument, NULL, 'b' },
          { "season_min", required_argument, NULL, 'g' },
          { "season_max", required_argument, NULL, 'h' },
          { "season_min2", required_argument, NULL, 'j' },
          { "season_max2", required_argument, NULL, 'k' },
          { "alt_min", required_argument, NULL, 'l' },
          { "alt_max", required_argument, NULL, 'm' },
          { "jump_thresh", required_argument, NULL, 'u' },
          { "smooth_alpha", required_argument, NULL, 'p' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "ab:c:d:e:f:g:h:j:k:l:m:o:p:t:s:u:", long_options, &option_index);
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
            data = read_champ(optarg);
            break;

            break;

          case 'd':
            params.downsample = (size_t) atoi(optarg);
            break;

          case 'e':
              fprintf(stderr, "main: reading Euler angles from %s...", optarg);
              euler_p = euler_read(optarg);
              if (!euler_p)
                exit(1);
              fprintf(stderr, "done (%zu sets of angles read)\n", euler_p->n);
              break;

          case 'o':
            output_file = optarg;
            break;

          case 't':
            params.ut = atof(optarg);
            break;

          case 'f':
            params.lt_min = atof(optarg);
            break;

          case 'b':
            params.lt_max = atof(optarg);
            break;

          case 'g':
            params.season_min = atof(optarg);
            break;

          case 'h':
            params.season_max = atof(optarg);
            break;

          case 'j':
            params.season_min2 = atof(optarg);
            break;

          case 'k':
            params.season_max2 = atof(optarg);
            break;

          case 'l':
            params.alt_min = atof(optarg);
            break;

          case 'm':
            params.alt_max = atof(optarg);
            break;

          case 'p':
            params.smooth_alpha = atof(optarg);
            break;

          case 'u':
            params.jump_thresh = atof(optarg);
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

  fprintf(stderr, "main: universal time   = %.1f\n", params.ut);
  fprintf(stderr, "main: LT minimum       = %.1f\n", params.lt_min);
  fprintf(stderr, "main: LT maximum       = %.1f\n", params.lt_max);
  fprintf(stderr, "main: altitude minimum = %.1f\n", params.alt_min);
  fprintf(stderr, "main: altitude maximum = %.1f\n", params.alt_max);
  fprintf(stderr, "main: QD minimum       = %.1f\n", params.qd_min);
  fprintf(stderr, "main: QD maximum       = %.1f\n", params.qd_max);
  fprintf(stderr, "main: season min       = %.1f\n", params.season_min);
  fprintf(stderr, "main: season max       = %.1f\n", params.season_max);
  fprintf(stderr, "main: season min 2     = %.1f\n", params.season_min2);
  fprintf(stderr, "main: season max 2     = %.1f\n", params.season_max2);
  fprintf(stderr, "main: smooth alpha     = %f\n", params.smooth_alpha);

  if (euler_p)
    {
      fprintf(stderr, "main: rotating VFM measurements with new Euler angles...");
      euler_apply(data, euler_p);
      fprintf(stderr, "done\n");
    }

  track_p = preprocess_data(&params, data);

#if 1
  fprintf(stderr, "main: recomputing longitude...");
  replace_longitude(12.0, data);
  fprintf(stderr, "done\n");
#endif

  fprintf(stderr, "main: computing vector residuals...");
  gettimeofday(&tv0, NULL);
  mdata = copy_data(data, track_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu residuals copied, %g seconds)\n",
          mdata->n, time_diff(tv0, tv1));

  /* free data after copying arrays to free up memory */
  satdata_mag_free(data);

  magdata_init(mdata);

  fprintf(stderr, "main: computing spatial weighting of data...");
  gettimeofday(&tv0, NULL);
  magdata_calc(mdata);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "main: writing data to %s...", data_prefix);
  magdata_print(data_prefix, mdata);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: writing data map to %s...", datamap_file);
  magdata_map(datamap_file, mdata);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: satellite rmin = %.1f (%.1f) [km]\n",
          mdata->rmin, mdata->rmin - mdata->R);
  fprintf(stderr, "main: satellite rmax = %.1f (%.1f) [km]\n",
          mdata->rmax, mdata->rmax - mdata->R);

  if (output_file)
    {
      fprintf(stderr, "main: writing poltor data to %s...", output_file);
      magdata_write(output_file, mdata);
      fprintf(stderr, "done\n");
    }

  magdata_free(mdata);
  track_free(track_p);

  return 0;
} /* main() */
