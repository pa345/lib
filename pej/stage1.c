/*
 * stage1.c
 *
 * Select satellite data according to IMF critera
 *
 * Steps are:
 * 1. select instrument flags (recommended CHAMP flags except 1 star camera allowed)
 * 2. select track rms test
 * 3. downsample by factor 15
 * 4. convert to magdata format and store to disk
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

#include "magdata.h"
#include "track.h"

/* define this to allow only 1 camera for data selection */
#define PRINT_ONE_CAMERA           1

typedef struct
{
  int all;            /* print all tracks */
  double kp_min;      /* minimum kp */
  double kp_max;      /* maximum kp */
  double IMF_Bz_min;  /* minimum IMF B_z (nT) */
  double IMF_Bz_max;  /* maximum IMF B_z (nT) */
  size_t downsample;  /* downsampling factor */
  double thresh[4];   /* rms thresholds */
} preprocess_parameters;

satdata_mag *
read_swarm(const char *filename)
{
  size_t nflag;
  satdata_mag *data;
  struct timeval tv0, tv1;

  fprintf(stderr, "read_swarm: reading %s...", filename);
  gettimeofday(&tv0, NULL);
  data = satdata_swarm_read_idx(filename, 0);
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

  print_track_stats(data, track_p);

  return track_p;
} /* preprocess_data() */

magdata *
copy_data(const satdata_mag *data, const track_workspace *track_p)
{
  const size_t nflagged = satdata_nflagged(data);
  size_t ndata = data->n - nflagged;
  size_t npts[6] = { 0, 0, 0, 0, 0, 0 };
  magdata_params params;
  magdata *mdata = magdata_alloc(ndata, data->R);
  size_t i;

  if (!mdata)
    return 0;

  params.grad_dt_ns = 10;
  params.grad_dt_ew = 10;
  params.grad_dphi_max = 4.0;
  params.grad_dlat_max = 4.0;
  params.model_main = 1;
  params.model_crust = 1;
  params.model_ext = 1;

  fprintf(stderr, "\n");
  fprintf(stderr, "\t copy_data: copying tracks in magdata format...");

  /* copy tracks into mdata structure */
  for (i = 0; i < track_p->n; ++i)
    {
      track_data *tptr = &(track_p->tracks[i]);

      /* discard flagged tracks */
      if (tptr->flags)
        continue;

      magdata_copy_track(&params, i, data, track_p, mdata, npts);
    }

  fprintf(stderr, "done (ndata = %zu mdata_n = %zu, mdata_ntot = %zu)\n", ndata, mdata->n, mdata->ntot);

  return mdata;
}

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options]\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --swarm_file | -s swarm_index_file          - Swarm index file\n");
  fprintf(stderr, "\t --champ_file | -c champ_index_file          - CHAMP index file\n");
  fprintf(stderr, "\t --dmsp_file | -D dmsp_index_file            - DMSP index file\n");
  fprintf(stderr, "\t --all | -a                                  - print all tracks (no filtering)\n");
  fprintf(stderr, "\t --downsample | -d downsample                - downsampling factor\n");
  fprintf(stderr, "\t --output_file | -o output_file              - output file\n");
  fprintf(stderr, "\t --kp_min | -v kp_min                        - kp minimum\n");
  fprintf(stderr, "\t --kp_max | -w kp_max                        - kp maximum\n");
  fprintf(stderr, "\t --imf_bz_min | -A imf_bz_min                - minimum IMF B_z (nT)\n");
  fprintf(stderr, "\t --imf_bz_max | -B imf_bz_max                - maximum IMF B_z (nT)\n");
}

int
main(int argc, char *argv[])
{
  char *data_file = "data/data.dat";
  satdata_mag *data = NULL;
  struct timeval tv0, tv1;
  track_workspace *track_p;
  preprocess_parameters params;
  magdata *mdata;

  /* defaults */
  params.all = 0;
  params.kp_min = 0.0;
  params.kp_max = 20.0;
  params.IMF_Bz_min = -1.0e6;
  params.IMF_Bz_max = 1.0e6;
  params.downsample = 15;
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
          { "dmsp_file", required_argument, NULL, 'D' },
          { "downsample", required_argument, NULL, 'd' },
          { "kp_min", required_argument, NULL, 'v' },
          { "kp_max", required_argument, NULL, 'w' },
          { "output_file", required_argument, NULL, 'o' },
          { "imf_bz_min", required_argument, NULL, 'A' },
          { "imf_bz_max", required_argument, NULL, 'B' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "ac:d:D:o:s:A:B:", long_options, &option_index);
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

          case 'd':
            params.downsample = (size_t) atoi(optarg);
            break;

          case 'o':
            data_file = optarg;
            break;

          case 'v':
            params.kp_min = atof(optarg);
            break;

          case 'w':
            params.kp_max = atof(optarg);
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

  fprintf(stderr, "main: kp minimum       = %.1f\n", params.kp_min);
  fprintf(stderr, "main: kp maximum       = %.1f\n", params.kp_max);
  fprintf(stderr, "main: IMF B_z minimum  = %.1f\n", params.IMF_Bz_min);
  fprintf(stderr, "main: IMF B_z maximum  = %.1f\n", params.IMF_Bz_max);

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

  fprintf(stderr, "main: converting to magdata format...");
  mdata = copy_data(data, track_p);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: writing data to %s...", data_file);
  magdata_write(data_file, mdata);
  fprintf(stderr, "done\n");

  magdata_free(mdata);
  satdata_mag_free(data);
  track_free(track_p);

  return 0;
}
