/*
 * preproc.c
 *
 * Pre-process satellite data and save in magdata format. This
 * program should be called with Stage1 data, in order to do track
 * rms test, and along-track gradient calculation
 *
 * Pre-processing steps are:
 * 1. Instrument flags (recommended CHAMP flags except 1 star camera allowed)
 * 2. Track rms test
 * 3. filter for WMM criteria
 * 4. Downsample by factor 20
 * 5. Compute and store along-track gradients
 *
 * The result is an output file in magdata format containing all data point to
 * be used in the modeling. All data points will have a MAGDATA_FLG_FIT_xxx flag
 * set, and other flags will vary.
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
#include <flow/flow.h>
#include <indices/indices.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_statistics.h>

#include <common/common.h>
#include <track/track.h>

#include "euler.h"
#include "magdata.h"

static size_t check_data_point(const size_t magdata_flags, const size_t idx,
                               const satdata_mag *data);

/* local-time range for Euler angle fitting */
#define EULER_LT_MIN        (6.0)
#define EULER_LT_MAX        (18.0)

/* maximum QD latitude for fitting Euler angles */
#define EULER_QDLAT         (55.5)

typedef struct
{
  int flag_vec_rms;   /* flag NEC vector rms */
  size_t downsample;  /* downsampling factor */
} preprocess_parameters;

/*
copy_track()
  Copy a single track to magdata structure

Inputs: track_idx - track index
        data      - satellite data
        track_p   - track data
        mdata     - (output) where to store new data
        ndata     - (output) updated to count scalar/vector/euler/gradient
                    data
                    ndata[0] = # scalar measurements for MF modeling
                    ndata[1] = # vector measurements for MF modeling
                    ndata[2] = # vector measurements for Euler angles
                    ndata[3] = # along-track gradient measurements
*/

size_t
copy_track(const size_t track_idx, const satdata_mag *data,
           const track_workspace *track_p, magdata *mdata,
           size_t npts[4])
{
  int s = 0;
  size_t i, k;
  track_data *tptr = &(track_p->tracks[track_idx]);
  const size_t start_idx = tptr->start_idx;
  const size_t end_idx = tptr->end_idx;
  magdata_datum datum;
  size_t ngrad = 0;
  int flagged_start = 0;

  for (i = start_idx; i <= end_idx; ++i)
    {
      size_t flags = 0;     /* final magdata flags */
      size_t fitting_flags;

      if (data->flags[i])
        continue;

      /*
       * determine whether this data point will be used for MF coefficients,
       * Euler angles, or both
       */
      fitting_flags = check_data_point(mdata->global_flags, i, data);
      if (fitting_flags == 0)
        continue; /* this point will not be used in the modeling */

      flags |= fitting_flags;

      /* initialize to 0 */
      magdata_datum_init(&datum);

      if (!flagged_start)
        {
          flags |= MAGDATA_FLG_TRACK_START;
          flagged_start = 1;
        }

      datum.B_nec[0] = SATDATA_VEC_X(data->B, i);
      datum.B_nec[1] = SATDATA_VEC_Y(data->B, i);
      datum.B_nec[2] = SATDATA_VEC_Z(data->B, i);
      datum.B_vfm[0] = SATDATA_VEC_X(data->B_VFM, i);
      datum.B_vfm[1] = SATDATA_VEC_Y(data->B_VFM, i);
      datum.B_vfm[2] = SATDATA_VEC_Z(data->B_VFM, i);
      datum.F = data->F[i];

      for (k = 0; k < 4; ++k)
        datum.q[k] = data->q[4 * i + k];

      /* subtract main field model from measurements */
      datum.B_model[0] = SATDATA_VEC_X(data->B_main, i);
      datum.B_model[1] = SATDATA_VEC_Y(data->B_main, i);
      datum.B_model[2] = SATDATA_VEC_Z(data->B_main, i);

      /* subtract external field model from measurements */
      datum.B_model[0] += SATDATA_VEC_X(data->B_ext, i);
      datum.B_model[1] += SATDATA_VEC_Y(data->B_ext, i);
      datum.B_model[2] += SATDATA_VEC_Z(data->B_ext, i);

      /* subtract crustal field also */
      datum.B_model[0] += SATDATA_VEC_X(data->B_crust, i);
      datum.B_model[1] += SATDATA_VEC_Y(data->B_crust, i);
      datum.B_model[2] += SATDATA_VEC_Z(data->B_crust, i);

      /* scalar measurement available (only for MF fitting) */
      if (fitting_flags & MAGDATA_FLG_FIT_MF)
        flags |= MAGDATA_FLG_F;

      /* vector measurement available for Euler (if below EULER_QDLAT) */
      if (fabs(data->qdlat[i]) <= EULER_QDLAT)
        flags |= MAGDATA_FLG_X | MAGDATA_FLG_Y | MAGDATA_FLG_Z;

      if (!(flags & (MAGDATA_FLG_F|MAGDATA_FLG_X|MAGDATA_FLG_Y|MAGDATA_FLG_Z)))
        {
          fprintf(stderr, "UH OH\n");
        }

      datum.t = data->t[i];
      datum.r = data->altitude[i] + data->R;
      datum.theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
      datum.phi = data->longitude[i] * M_PI / 180.0;
      datum.qdlat = data->qdlat[i];
      datum.satdir = satdata_mag_satdir(i, data);
      datum.flags = flags;

      s = magdata_add(&datum, mdata);
      if (s == 0)
        {
          if (flags & MAGDATA_FLG_FIT_MF)
            {
              if (flags & MAGDATA_FLG_F)
                ++(npts[0]);

              if (flags & (MAGDATA_FLG_X|MAGDATA_FLG_Y|MAGDATA_FLG_Z))
                ++(npts[1]);
            }

          if (flags & MAGDATA_FLG_FIT_EULER)
            {
              if (flags & (MAGDATA_FLG_X|MAGDATA_FLG_Y|MAGDATA_FLG_Z))
                ++(npts[2]);
            }

          if (flags & MAGDATA_FLG_DZ_NS)
            ++npts[3];
        }

      if (s == 0 && flags & MAGDATA_FLG_DZ_NS)
        ++ngrad;
    }

  return ngrad;
} /* copy_track() */

magdata *
copy_data(const size_t magdata_flags, const satdata_mag *data, const track_workspace *track_p)
{
  const size_t nflagged = satdata_nflagged(data);
  size_t ndata = data->n - nflagged;
  magdata *mdata;
  size_t i;
  size_t ngrad = 0;
  size_t npts[4] = { 0, 0, 0, 0 };

  /*
   * for some reason, the loop below adds more than 'ndata'
   * points into mdata, it shouldn't do this but I'm not sure why
   * so add 10000
   */
  mdata = magdata_alloc(ndata + 10000, data->R);
  if (!mdata)
    return 0;

  mdata->global_flags = magdata_flags;

  for (i = 0; i < track_p->n; ++i)
    {
      track_data *tptr = &(track_p->tracks[i]);

      /* discard flagged tracks */
      if (tptr->flags)
        continue;

      ngrad += copy_track(i, data, track_p, mdata, npts);
    }

  ndata = mdata->n;

  fprintf(stderr, "\n");
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) scalar measurements for MF modeling\n",
          npts[0], mdata->n, (double) npts[0] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) vector measurements for MF modeling\n",
          npts[1], mdata->n, (double) npts[1] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) vector measurements for Euler angle modeling\n",
          npts[2], mdata->n, (double) npts[2] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) along-track gradient data available\n",
          ngrad, ndata, (double)ngrad / (double)ndata * 100.0);

  return mdata;
}

satdata_mag *
read_swarm(const char *filename, const int asmv)
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

  /* additional checks for ASM-V data */
  if (asmv)
    {
      size_t i;
      size_t nasmv = 0;

      fprintf(stderr, "read_swarm: checking ASM-V data...");

      for (i = 0; i < data->n; ++i)
        {
          double B_asmv[4];

          B_asmv[0] = SATDATA_VEC_X(data->B_VFM, i);
          B_asmv[1] = SATDATA_VEC_Y(data->B_VFM, i);
          B_asmv[2] = SATDATA_VEC_Z(data->B_VFM, i);
          B_asmv[3] = vec_norm(B_asmv);

          /* check for missing vector data */
          if (B_asmv[3] < 50.0 || data->F[i] < 50.0)
            {
              data->flags[i] |= SATDATA_FLG_OUTLIER;
              nasmv++;
            }

          /* normalize so that |B_asmv| = F */
          SATDATA_VEC_X(data->B_VFM, i) *= data->F[i] / B_asmv[3];
          SATDATA_VEC_Y(data->B_VFM, i) *= data->F[i] / B_asmv[3];
          SATDATA_VEC_Z(data->B_VFM, i) *= data->F[i] / B_asmv[3];
        }

      fprintf(stderr, "done (%zu/%zu (%.1f%%) missing data flagged)\n",
              nasmv, data->n, (double)nasmv / (double)data->n * 100.0);
    }

  return data;
} /* read_swarm() */

satdata_mag *
read_champ(const char *filename)
{
  satdata_mag *data;
  struct timeval tv0, tv1;
  size_t nflag;

  fprintf(stderr, "read_champ: reading %s...", optarg);
  gettimeofday(&tv0, NULL);
  data = satdata_champ_read_idx(optarg, 0);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu data read, %g seconds)\n",
          data->n, time_diff(tv0, tv1));

  /* check for instrument flags since we use Stage1 data */
  fprintf(stderr, "read_champ: filtering for instrument flags...");
  nflag = satdata_champ_filter_instrument(1, 0, data);
  fprintf(stderr, "done (%zu/%zu (%.1f%%) data flagged)\n",
          nflag, data->n, (double)nflag / (double)data->n * 100.0);

  return data;
} /* read_champ() */

/*
check_data_point()
  Check an individual data point to determine if it will be used to
fit various model parameters

Inputs: magdata_flags - MAGDATA_GLOBFLG_xxx
        idx           - index into 'data'
        data          - satellite data

Return: flags indicating fit parameters (MAGDATA_FLG_FIT_xxx)
*/

static size_t
check_data_point(const size_t magdata_flags, const size_t idx,
                 const satdata_mag *data)
{
  size_t flags = 0;
  const time_t unix_time = satdata_epoch2timet(data->t[idx]);
  const double lt = get_localtime(unix_time, data->longitude[idx] * M_PI / 180.0);

  /* check if we should fit Euler angles to this data point */
  if (magdata_flags & MAGDATA_GLOBFLG_EULER &&
      fabs(data->qdlat[idx]) <= EULER_QDLAT &&
      (lt >= EULER_LT_MAX || lt <= EULER_LT_MIN))
    {
      flags |= MAGDATA_FLG_FIT_EULER;
    }

  return flags;
} /* check_data_point() */

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
          downsample - downsampling factor
        data   - satellite data

Return: pointer to sorted track workspace (should be freed by caller)
*/

track_workspace *
preprocess_data(const preprocess_parameters *params, satdata_mag *data)
{
  struct timeval tv0, tv1;
  track_workspace *track_p = track_alloc();
  double rms_thresh[] = { 20.0, 25.0, 15.0, 15.0 };

  fprintf(stderr, "preprocess_data: initializing tracks...");
  gettimeofday(&tv0, NULL);
  track_init(data, NULL, track_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  if (!params->flag_vec_rms)
    {
      rms_thresh[0] = -1.0;
      rms_thresh[1] = -1.0;
      rms_thresh[2] = -1.0;
    }

  /* detect bad tracks with rms test */
  {
    const char *rmsfile = "satrms.dat";
    size_t nrms;

    nrms = track_flag_rms(rmsfile, rms_thresh, NULL, data, track_p);
    fprintf(stderr, "preprocess_data: flagged (%zu/%zu) (%.1f%%) tracks due to high rms\n",
            nrms, track_p->n, (double) nrms / (double) track_p->n * 100.0);
  }

  /* flag according to WMM criteria */
  satdata_filter_wmm(1, data);

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

  /* print track statistics */
  {
    char *stat_file = "track_stats.dat";

    fprintf(stderr, "preprocess_data: printing track statistics to %s...", stat_file);
    track_print_stats(stat_file, track_p);
    fprintf(stderr, "done\n");
  }

  print_track_stats(data, track_p);

  return track_p;
} /* preprocess_data() */

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options]\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --champ_file      | -c champ_index_file       - CHAMP index file\n");
  fprintf(stderr, "\t --swarm_file      | -s swarm_index_file       - Swarm index file\n");
  fprintf(stderr, "\t --swarm_asmv_file | -a swarm_asmv_index_file  - Swarm ASM-V index file\n");
  fprintf(stderr, "\t --downsample      | -d downsample             - downsampling factor\n");
  fprintf(stderr, "\t --euler_file      | -e euler_file             - Euler angles file\n");
  fprintf(stderr, "\t --output_file     | -o output_file            - output file\n");
}

int
main(int argc, char *argv[])
{
  char *datamap_file = "datamap.dat";
  char *data_file = "data.dat";
  char *output_file = NULL;
  satdata_mag *data = NULL;
  magdata *mdata;
  euler_workspace *euler_p = NULL;
  struct timeval tv0, tv1;
  track_workspace *track_p;
  preprocess_parameters params;
  size_t magdata_flags = 0;

  /* defaults */
  params.flag_vec_rms = 1;
  params.downsample = 20;

  magdata_flags = MAGDATA_GLOBFLG_EULER;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "swarm_file", required_argument, NULL, 's' },
          { "swarm_asmv_file", required_argument, NULL, 'a' },
          { "champ_file", required_argument, NULL, 'c' },
          { "downsample", required_argument, NULL, 'd' },
          { "output_file", required_argument, NULL, 'o' },
          { "euler_file", required_argument, NULL, 'e' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "a:c:d:e:o:s:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          /* Swarm ASM-V */
          case 'a':
            data = read_swarm(optarg, 1);
            params.flag_vec_rms = 0; /* no NEC data for rms flagging */
            break;

          /* Swarm official */
          case 's':
            data = read_swarm(optarg, 0);
            break;

          case 'c':
            data = read_champ(optarg);
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

          default:
            break;
        }
    }

  if (!data)
    {
      print_help(argv);
      exit(1);
    }

  fprintf(stderr, "main: LT minimum       = %.1f\n", EULER_LT_MIN);
  fprintf(stderr, "main: LT maximum       = %.1f\n", EULER_LT_MAX);

  if (euler_p)
    {
      fprintf(stderr, "main: rotating VFM measurements with new Euler angles...");
      euler_apply(data, euler_p);
      fprintf(stderr, "done\n");
    }

  track_p = preprocess_data(&params, data);

  fprintf(stderr, "main: computing vector residuals...");
  gettimeofday(&tv0, NULL);
  mdata = copy_data(magdata_flags, data, track_p);
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

  fprintf(stderr, "main: writing data to %s...", data_file);
  magdata_print(data_file, mdata);
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
      fprintf(stderr, "main: writing data to %s...", output_file);
      magdata_write(output_file, mdata);
      fprintf(stderr, "done\n");
    }

  magdata_free(mdata);
  track_free(track_p);

  return 0;
} /* main() */
