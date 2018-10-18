/*
 * preproc.c
 *
 * Routines related to preprocessing data. Calling sequence:
 *
 * 1. magdata_preprocess_default_parameters - build struct of default parameters
 * 2. magdata_preprocess_parse              - parse configuration file
 * 3. magdata_preprocess_check              - check parameters have sane values
 * 4. magdata_preprocess                    - perform data selection and preprocessing
 * 5. magdata_preprocess_fill               - fill magdata structure with preprocessed data
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include <gsl/gsl_math.h>

#include <libconfig.h>
#include <common/common.h>
#include <track/track.h>

#include "magdata.h"

static int print_track_stats(const satdata_mag *data, const track_workspace *track_p);
static int polar_damp_apply(const magdata_preprocess_parameters *params, track_workspace *track_p, satdata_mag *data);

#include "preproc_filter.c"

#define MAGDATA_IDX_X              0
#define MAGDATA_IDX_Y              1
#define MAGDATA_IDX_Z              2
#define MAGDATA_IDX_F              3
#define MAGDATA_IDX_DX_NS          4
#define MAGDATA_IDX_DY_NS          5
#define MAGDATA_IDX_DZ_NS          6
#define MAGDATA_IDX_DF_NS          7
#define MAGDATA_IDX_DX_EW          8
#define MAGDATA_IDX_DY_EW          9
#define MAGDATA_IDX_DZ_EW          10
#define MAGDATA_IDX_DF_EW          11
#define MAGDATA_IDX_B_EULER        12
#define MAGDATA_IDX_END            13

/*
magdata_preprocess_default_parameters()
  Return reasonable default parameters for preprocessing
*/

magdata_preprocess_parameters
magdata_preprocess_default_parameters(void)
{
  magdata_preprocess_parameters params;

  params.downsample = 0;
  params.max_kp = -1.0;
  params.max_dRC = -1.0;
  params.min_LT = -1.0;
  params.max_LT = -1.0;
  params.euler_min_LT = -1.0;
  params.euler_max_LT = -1.0;
  params.qdlat_preproc_cutoff = -1.0;
  params.min_zenith = -1.0;
  params.flag_IMF = 0;
  params.season_min = -1.0;
  params.season_max = -1.0;
  params.season_min2 = -1.0;
  params.season_max2 = -1.0;
  params.rmin = -1.0;
  params.rmax = -1.0;
  params.rms_thresh[0] = -1.0;
  params.rms_thresh[1] = -1.0;
  params.rms_thresh[2] = -1.0;
  params.rms_thresh[3] = -1.0;
  params.gradient_ns = 0;
  params.gradew_dphi_max = -1.0;
  params.gradew_dlat_max = -1.0;
  params.gradew_dt_max = -1.0;
  params.subtract_B_main = -1;
  params.subtract_B_crust = -1;
  params.subtract_B_ext = -1;
  params.polar_damping = 0;
  params.polar_qdlat = -1.0;
  params.pb_flag = -1;
  params.pb_qdmax = -1.0;
  params.pb_thresh[0] = -1.0;
  params.pb_thresh[1] = -1.0;
  params.pb_thresh[2] = -1.0;
  params.pb_thresh[3] = -1.0;

  return params;
}

/*
magdata_parse_config_file()
  Parse configuration file with preprocessing parameters

Inputs: filename - config file name
        params   - (output) configuration parameters

Return: success/error
*/

int
magdata_preprocess_parse(const char *filename, magdata_preprocess_parameters *params)
{
  int s;
  config_t cfg;
  double fval;
  int ival;

  config_init(&cfg);

  s = config_read_file(&cfg, filename);
  if (s != CONFIG_TRUE)
    {
      fprintf(stderr, "magdata_preprocess_parse: %s:%d - %s\n",
              config_error_file(&cfg),
              config_error_line(&cfg),
              config_error_text(&cfg));
      config_destroy(&cfg);
      return -1;
    }

  if (config_lookup_float(&cfg, "max_kp", &fval))
    params->max_kp = fval;
  if (config_lookup_float(&cfg, "max_dRC", &fval))
    params->max_dRC = fval;
  if (config_lookup_float(&cfg, "min_LT", &fval))
    params->min_LT = fval;
  if (config_lookup_float(&cfg, "max_LT", &fval))
    params->max_LT = fval;
  if (config_lookup_float(&cfg, "euler_min_LT", &fval))
    params->euler_min_LT = fval;
  if (config_lookup_float(&cfg, "euler_max_LT", &fval))
    params->euler_max_LT = fval;
  if (config_lookup_float(&cfg, "qdlat_preproc_cutoff", &fval))
    params->qdlat_preproc_cutoff = fval;
  if (config_lookup_float(&cfg, "min_zenith", &fval))
    params->min_zenith = fval;
  if (config_lookup_float(&cfg, "season_min", &fval))
    params->season_min = fval;
  if (config_lookup_float(&cfg, "season_max", &fval))
    params->season_max = fval;
  if (config_lookup_float(&cfg, "season_min2", &fval))
    params->season_min2 = fval;
  if (config_lookup_float(&cfg, "season_max2", &fval))
    params->season_max2 = fval;
  if (config_lookup_float(&cfg, "alt_min", &fval))
    params->rmin = R_EARTH_KM + fval;
  if (config_lookup_float(&cfg, "alt_max", &fval))
    params->rmax = R_EARTH_KM + fval;

  if (config_lookup_int(&cfg, "downsample", &ival))
    params->downsample = (size_t) ival;
  if (config_lookup_int(&cfg, "gradient_ns", &ival))
    params->gradient_ns = (size_t) ival;
  if (config_lookup_int(&cfg, "fit_track_RC", &ival))
    params->fit_track_RC = (size_t) ival;
  if (config_lookup_int(&cfg, "flag_IMF", &ival))
    params->flag_IMF = (size_t) ival;

  if (config_lookup_int(&cfg, "polar_damping", &ival))
    params->polar_damping = (size_t) ival;
  if (config_lookup_float(&cfg, "polar_qdlat", &fval))
    params->polar_qdlat = fval;

  if (config_lookup_int(&cfg, "subtract_B_main", &ival))
    params->subtract_B_main = (size_t) ival;
  if (config_lookup_int(&cfg, "subtract_B_crust", &ival))
    params->subtract_B_crust = (size_t) ival;
  if (config_lookup_int(&cfg, "subtract_B_ext", &ival))
    params->subtract_B_ext = (size_t) ival;

  if (config_lookup_float(&cfg, "gradient_ew_dphi_max", &fval))
    params->gradew_dphi_max = fval;
  if (config_lookup_float(&cfg, "gradient_ew_dlat_max", &fval))
    params->gradew_dlat_max = fval;
  if (config_lookup_float(&cfg, "gradient_ew_dt_max", &fval))
    params->gradew_dt_max = fval;

  if (config_lookup_float(&cfg, "rms_threshold_X", &fval))
    params->rms_thresh[0] = fval;
  if (config_lookup_float(&cfg, "rms_threshold_Y", &fval))
    params->rms_thresh[1] = fval;
  if (config_lookup_float(&cfg, "rms_threshold_Z", &fval))
    params->rms_thresh[2] = fval;
  if (config_lookup_float(&cfg, "rms_threshold_F", &fval))
    params->rms_thresh[3] = fval;

  if (config_lookup_int(&cfg, "pb_flag", &ival))
    params->pb_flag = ival;
  if (config_lookup_float(&cfg, "pb_qdmax", &fval))
    params->pb_qdmax = fval;
  if (config_lookup_float(&cfg, "pb_threshold_dX", &fval))
    params->pb_thresh[0] = fval;
  if (config_lookup_float(&cfg, "pb_threshold_dY", &fval))
    params->pb_thresh[1] = fval;
  if (config_lookup_float(&cfg, "pb_threshold_dZ", &fval))
    params->pb_thresh[2] = fval;
  if (config_lookup_float(&cfg, "pb_threshold_dF", &fval))
    params->pb_thresh[3] = fval;

  config_destroy(&cfg);

  return 0;
}

/*
magdata_preprocess_check()
  Check parameters have reasonable values
*/

int
magdata_preprocess_check(magdata_preprocess_parameters * params)
{
  int s = 0;

  if (params->downsample == 0)
    {
      fprintf(stderr, "magdata_preprocess_check: downsample must be > 0\n");
      ++s;
    }

  if (params->max_kp <= 0.0)
    {
      fprintf(stderr, "magdata_preprocess_check: max_kp must be > 0\n");
      ++s;
    }

  if (params->max_dRC <= 0.0)
    {
      fprintf(stderr, "magdata_preprocess_check: max_dRC must be > 0\n");
      ++s;
    }

  if (params->qdlat_preproc_cutoff < 0.0)
    {
      fprintf(stderr, "magdata_preprocess_check: qdlat_preproc_cutoff must be > 0\n");
      ++s;
    }

  if (params->gradient_ns == 0)
    {
      fprintf(stderr, "magdata_preprocess_check: gradient_ns must be > 0\n");
      ++s;
    }

  if (params->gradew_dphi_max <= 0.0)
    {
      fprintf(stderr, "magdata_preprocess_check: gradient_ew_dphi_max must be > 0\n");
      ++s;
    }

  if (params->gradew_dlat_max <= 0.0)
    {
      fprintf(stderr, "magdata_preprocess_check: gradient_ew_dlat_max must be > 0\n");
      ++s;
    }

  if (params->gradew_dt_max <= 0.0)
    {
      fprintf(stderr, "magdata_preprocess_check: gradient_ew_dt_max must be > 0\n");
      ++s;
    }

  if (params->subtract_B_main < 0)
    {
      fprintf(stderr, "magdata_preprocess_check: subtract_B_main must be 0 or 1\n");
      ++s;
    }

  if (params->subtract_B_crust < 0)
    {
      fprintf(stderr, "magdata_preprocess_check: subtract_B_crust must be 0 or 1\n");
      ++s;
    }

  if (params->subtract_B_ext < 0)
    {
      fprintf(stderr, "magdata_preprocess_check: subtract_B_ext must be 0 or 1\n");
      ++s;
    }

  if (params->pb_flag < 0)
    {
      fprintf(stderr, "magdata_preprocess_check: pb_flag must be 0 or 1\n");
      ++s;
    }

  if (params->pb_qdmax <= 0.0)
    {
      fprintf(stderr, "magdata_preprocess_check: pb_qdmax must be > 0\n");
      ++s;
    }

  if (params->pb_thresh[0] <= 0.0)
    {
      fprintf(stderr, "magdata_preprocess_check: pb_threshold_dX must be > 0\n");
      ++s;
    }

  if (params->pb_thresh[1] <= 0.0)
    {
      fprintf(stderr, "magdata_preprocess_check: pb_threshold_dY must be > 0\n");
      ++s;
    }

  if (params->pb_thresh[2] <= 0.0)
    {
      fprintf(stderr, "magdata_preprocess_check: pb_threshold_dZ must be > 0\n");
      ++s;
    }

  if (params->pb_thresh[3] <= 0.0)
    {
      fprintf(stderr, "magdata_preprocess_check: pb_threshold_dF must be > 0\n");
      ++s;
    }

  if (params->polar_damping && params->polar_qdlat < 0.0)
    {
      fprintf(stderr, "magdata_preprocess_check: polar_damping enabled but polar_qdlat is < 0\n");
      ++s;
    }

  return s;
}

/*
magdata_preprocess()

Inputs: params        - preprocess parameters
        magdata_flags - MAGDATA_GLOBFLG_xxx flags
        data          - satellite data

Return: pointer to sorted track workspace (should be freed by caller)
*/

track_workspace *
magdata_preprocess(const magdata_preprocess_parameters *params, const size_t magdata_flags,
                   satdata_mag *data)
{
  struct timeval tv0, tv1;
  track_workspace *track_p = track_alloc();

  fprintf(stderr, "magdata_preprocess: initializing tracks...");
  gettimeofday(&tv0, NULL);
  track_init(data, NULL, track_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  /* detect bad tracks with rms test */
  {
    const char *rmsfile = "satrms.dat";
    size_t nrms;

    nrms = track_flag_rms(rmsfile, params->rms_thresh, NULL, data, track_p);
    fprintf(stderr, "magdata_preprocess: flagged (%zu/%zu) (%.1f%%) tracks due to high rms\n",
            nrms, track_p->n, (double) nrms / (double) track_p->n * 100.0);
  }

  if (params->polar_damping)
    {
      fprintf(stderr, "magdata_preprocess: applying cosine window to polar data...");
      gettimeofday(&tv0, NULL);
      polar_damp_apply(params, track_p, data);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));
    }

  /* select geomagnetically quiet data */
  fprintf(stderr, "magdata_preprocess: selecting geomagnetically quiet data...");
  magdata_preprocess_filter(magdata_flags, params, track_p, data);
  fprintf(stderr, "done\n");

  /* downsample data */
  {
    size_t i;

    fprintf(stderr, "magdata_preprocess: downsampling data by factor %zu...", params->downsample);

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

    fprintf(stderr, "magdata_preprocess: printing track statistics to %s...", stat_file);
    track_print_stats(stat_file, track_p);
    fprintf(stderr, "done\n");
  }

  print_track_stats(data, track_p);

  return track_p;
}

magdata *
magdata_preprocess_fill(const size_t magdata_flags, const satdata_mag *data, const track_workspace *track_p,
                        const size_t magdata_flags2, const satdata_mag *data2, const track_workspace *track_p2,
                        magdata_preprocess_parameters * preproc_params)
{
  const size_t nflagged = satdata_nflagged(data);
  size_t ndata = data->n - nflagged;
  magdata *mdata;
  magdata_params params;
  size_t i;
  size_t npts[6] = { 0, 0, 0, 0, 0, 0 };
  size_t nmodel[MAGDATA_IDX_END];

  params.grad_dt_ns = (double) preproc_params->gradient_ns;
  params.grad_dt_ew = preproc_params->gradew_dt_max;
  params.grad_dphi_max = preproc_params->gradew_dphi_max;
  params.grad_dlat_max = preproc_params->gradew_dlat_max;

  /* subtract main field from data prior to modeling */
  if (preproc_params->subtract_B_main)
    params.model_main = 1;
  else
    params.model_main = 0;

  /* subtract crustal field model from data prior to modeling */
  if (preproc_params->subtract_B_crust)
    params.model_crust = 1;
  else
    params.model_crust = 0;

  /* subtract external field from data prior to modeling */
  if (preproc_params->subtract_B_ext)
    params.model_ext = 1;
  else
    params.model_ext = 0;

  /* initialize arrays */
  for (i = 0; i < MAGDATA_IDX_END; ++i)
    nmodel[i] = 0;

  mdata = magdata_alloc(ndata, data->R);
  if (!mdata)
    return 0;

  mdata->global_flags = magdata_flags;

  /* copy tracks into mdata structure */
  for (i = 0; i < track_p->n; ++i)
    {
      track_data *tptr = &(track_p->tracks[i]);

      /* discard flagged tracks */
      if (tptr->flags)
        continue;

      if (data2 == NULL)
        magdata_copy_track(&params, i, data, track_p, mdata, npts);
      else
        magdata_copy_track_EW(&params, i, data, track_p, data2, track_p2, mdata, npts);
    }

  fprintf(stderr, "ndata = %zu mdata_n = %zu, mdata_ntot = %zu\n", ndata, mdata->n, mdata->ntot);

  /*
   * now determine which points in mdata will be used for
   * MF modeling, Euler angle fitting, etc
   */
  for (i = 0; i < mdata->n; ++i)
    {
      size_t fitting_flags = 0;

      if (fabs(mdata->qdlat[i]) <= preproc_params->qdlat_preproc_cutoff)
        {
          /*
           * mid-latitude point: check local time of equator crossing to determine whether
           * to fit field model and/or Euler angles
           */

          double LT = mdata->lt_eq[i];
          int fit_MF = check_LT(LT, preproc_params->min_LT, preproc_params->max_LT);
          int fit_euler = check_LT(LT, preproc_params->euler_min_LT, preproc_params->euler_max_LT) &&
                          (mdata->global_flags & MAGDATA_GLOBFLG_EULER);

          if (fit_MF)
            fitting_flags |= MAGDATA_FLG_FIT_MF;

          if (fit_euler)
            fitting_flags |= MAGDATA_FLG_FIT_EULER;
        }
      else
        {
          /* high-latitude point - fit field model but not Euler angles */
          fitting_flags |= MAGDATA_FLG_FIT_MF;
        }

      if (fitting_flags & MAGDATA_FLG_FIT_MF)
        {
          if (MAGDATA_ExistX(mdata->flags[i]))
            ++(nmodel[MAGDATA_IDX_X]);

          if (MAGDATA_ExistY(mdata->flags[i]))
            ++(nmodel[MAGDATA_IDX_Y]);

          if (MAGDATA_ExistZ(mdata->flags[i]))
            ++(nmodel[MAGDATA_IDX_Z]);

          if (MAGDATA_ExistScalar(mdata->flags[i]))
            ++(nmodel[MAGDATA_IDX_F]);

          /* count north-south gradient data */

          if (MAGDATA_ExistDX_NS(mdata->flags[i]))
            ++(nmodel[MAGDATA_IDX_DX_NS]);

          if (MAGDATA_ExistDY_NS(mdata->flags[i]))
            ++(nmodel[MAGDATA_IDX_DY_NS]);

          if (MAGDATA_ExistDZ_NS(mdata->flags[i]))
            ++(nmodel[MAGDATA_IDX_DZ_NS]);

          if (MAGDATA_ExistDF_NS(mdata->flags[i]))
            ++(nmodel[MAGDATA_IDX_DF_NS]);

          /* count east-west gradient data */

          if (MAGDATA_ExistDX_EW(mdata->flags[i]))
            ++(nmodel[MAGDATA_IDX_DX_EW]);

          if (MAGDATA_ExistDY_EW(mdata->flags[i]))
            ++(nmodel[MAGDATA_IDX_DY_EW]);

          if (MAGDATA_ExistDZ_EW(mdata->flags[i]))
            ++(nmodel[MAGDATA_IDX_DZ_EW]);

          if (MAGDATA_ExistDF_EW(mdata->flags[i]))
            ++(nmodel[MAGDATA_IDX_DF_EW]);
        }

      if ((fitting_flags & MAGDATA_FLG_FIT_EULER) &&
          MAGDATA_ExistVector(mdata->flags[i]))
        ++(nmodel[MAGDATA_IDX_B_EULER]);

      mdata->flags[i] |= fitting_flags;
    }

  fprintf(stderr, "\n");
  fprintf(stderr, "\t magdata_preprocess_fill: %zu/%zu (%.1f%%) scalar measurements available\n",
          npts[0], mdata->n, (double) npts[0] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t magdata_preprocess_fill: %zu/%zu (%.1f%%) vector measurements available\n",
          npts[1], mdata->n, (double) npts[1] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t magdata_preprocess_fill: %zu/%zu (%.1f%%) along-track scalar measurements available\n",
          npts[2], mdata->n, (double) npts[2] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t magdata_preprocess_fill: %zu/%zu (%.1f%%) along-track vector measurements available\n",
          npts[3], mdata->n, (double) npts[3] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t magdata_preprocess_fill: %zu/%zu (%.1f%%) east-west scalar measurements available\n",
          npts[4], mdata->n, (double) npts[4] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t magdata_preprocess_fill: %zu/%zu (%.1f%%) east-west vector measurements available\n",
          npts[5], mdata->n, (double) npts[5] / (double) mdata->n * 100.0);

  fprintf(stderr, "\t magdata_preprocess_fill: %zu/%zu (%.1f%%) X vector measurements selected for MF modeling\n",
          nmodel[MAGDATA_IDX_X], mdata->n, (double) nmodel[MAGDATA_IDX_X] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t magdata_preprocess_fill: %zu/%zu (%.1f%%) Y vector measurements selected for MF modeling\n",
          nmodel[MAGDATA_IDX_Y], mdata->n, (double) nmodel[MAGDATA_IDX_Y] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t magdata_preprocess_fill: %zu/%zu (%.1f%%) Z vector measurements selected for MF modeling\n",
          nmodel[MAGDATA_IDX_Z], mdata->n, (double) nmodel[MAGDATA_IDX_Z] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t magdata_preprocess_fill: %zu/%zu (%.1f%%) scalar measurements selected for MF modeling\n",
          nmodel[MAGDATA_IDX_F], mdata->n, (double) nmodel[MAGDATA_IDX_F] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t magdata_preprocess_fill: %zu/%zu (%.1f%%) vector measurements selected for Euler angle modeling\n",
          nmodel[MAGDATA_IDX_B_EULER], mdata->n, (double) nmodel[MAGDATA_IDX_B_EULER] / (double) mdata->n * 100.0);

  fprintf(stderr, "\t magdata_preprocess_fill: %zu/%zu (%.1f%%) North/South dX vector measurements selected for MF modeling\n",
          nmodel[MAGDATA_IDX_DX_NS], mdata->n, (double) nmodel[MAGDATA_IDX_DX_NS] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t magdata_preprocess_fill: %zu/%zu (%.1f%%) North/South dY vector measurements selected for MF modeling\n",
          nmodel[MAGDATA_IDX_DY_NS], mdata->n, (double) nmodel[MAGDATA_IDX_DY_NS] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t magdata_preprocess_fill: %zu/%zu (%.1f%%) North/South dZ vector measurements selected for MF modeling\n",
          nmodel[MAGDATA_IDX_DZ_NS], mdata->n, (double) nmodel[MAGDATA_IDX_DZ_NS] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t magdata_preprocess_fill: %zu/%zu (%.1f%%) North/South dF scalar measurements selected for MF modeling\n",
          nmodel[MAGDATA_IDX_DF_NS], mdata->n, (double) nmodel[MAGDATA_IDX_DF_NS] / (double) mdata->n * 100.0);

  fprintf(stderr, "\t magdata_preprocess_fill: %zu/%zu (%.1f%%) East/West dX vector measurements selected for MF modeling\n",
          nmodel[MAGDATA_IDX_DX_EW], mdata->n, (double) nmodel[MAGDATA_IDX_DX_EW] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t magdata_preprocess_fill: %zu/%zu (%.1f%%) East/West dY vector measurements selected for MF modeling\n",
          nmodel[MAGDATA_IDX_DY_EW], mdata->n, (double) nmodel[MAGDATA_IDX_DY_EW] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t magdata_preprocess_fill: %zu/%zu (%.1f%%) East/West dZ vector measurements selected for MF modeling\n",
          nmodel[MAGDATA_IDX_DZ_EW], mdata->n, (double) nmodel[MAGDATA_IDX_DZ_EW] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t magdata_preprocess_fill: %zu/%zu (%.1f%%) East/West dF scalar measurements selected for MF modeling\n",
          nmodel[MAGDATA_IDX_DF_EW], mdata->n, (double) nmodel[MAGDATA_IDX_DF_EW] / (double) mdata->n * 100.0);

  return mdata;
}

static int
print_track_stats(const satdata_mag *data, const track_workspace *track_p)
{
  size_t nflagged = satdata_nflagged(data);
  size_t nleft = data->n - nflagged;
  size_t nflagged_track = track_nflagged(track_p);
  size_t nleft_track = track_p->n - nflagged_track;

  fprintf(stderr, "print_track_stats: total flagged data: %zu/%zu (%.1f%%)\n",
          nflagged, data->n, (double)nflagged / (double)data->n * 100.0);
  fprintf(stderr, "print_track_stats: total remaining data: %zu/%zu (%.1f%%)\n",
          nleft, data->n, (double)nleft / (double)data->n * 100.0);

  fprintf(stderr, "print_track_stats: total flagged tracks: %zu/%zu (%.1f%%)\n",
          nflagged_track, track_p->n, (double)nflagged_track / (double)track_p->n * 100.0);
  fprintf(stderr, "print_track_stats: total remaining tracks: %zu/%zu (%.1f%%)\n",
          nleft_track, track_p->n, (double)nleft_track / (double)track_p->n * 100.0);

  return 0;
}

/*
polar_damp_apply()
  Apply cosine window to polar data to reduce ionospheric signatures

The window function is defined in Chulliat et al, 2016. I use the square
of the window function defined in that paper

Inputs: params  - magdata parameters
        track_p - track workspace
        data    - satellite data
*/

static int
polar_damp_apply(const magdata_preprocess_parameters *params, track_workspace *track_p,
                 satdata_mag *data)
{
  int s = 0;
  const double thetaq0 = M_PI / 2.0 - params->polar_qdlat * M_PI / 180.0; /* cutoff QD co-latitude */
  size_t i, j, k;

  for (i = 0; i < track_p->n; ++i)
    {
      track_data *tptr = &(track_p->tracks[i]);
      size_t start_idx = tptr->start_idx;
      size_t end_idx = tptr->end_idx;

      for (j = start_idx; j <= end_idx; ++j)
        {
          double thetaq = M_PI / 2.0 - data->qdlat[j] * M_PI / 180.0; /* QD co-latitude */
          double B_res[4], B_model[4];
          double *B_obs = &(data->B[3 * j]);
          double window; /* window function */

          /* ignore bad data */
          if (SATDATA_BadData(data->flags[j]))
            continue;

          /* compute window values */
          if (thetaq < thetaq0)
            {
              /* close to north pole */
              window = 0.5 * (1.0 + cos(M_PI * (thetaq / thetaq0 - 1.0)));
            }
          else if (thetaq > (M_PI - thetaq0))
            {
              /* close to south pole */
              window = 0.5 * (1.0 + cos(M_PI * (thetaq / thetaq0 - M_PI / thetaq0 + 1.0)));
            }
          else
            {
              /* mid/low latitudes, leave signal intact */
              window = 1.0;
            }

          window *= window;

          /* calculate model values for this point */
          satdata_mag_model(j, B_model, data);

          if (SATDATA_ExistVector(data->flags[j]))
            {
              /* calculate vector residual for this point */
              satdata_mag_residual(j, B_res, data);

              /* apply correction to vector satellite data */
              for (k = 0; k < 3; ++k)
                {
                  /* apply window to residual at high latitudes */
                  B_res[k] *= window;

                  /* apply correction to satellite data */
                  B_obs[k] = B_res[k] + B_model[k];
                }

              /* recalculate scalar measurement */
              data->F[j] = gsl_hypot3(B_obs[0], B_obs[1], B_obs[2]);
            }
          else if (SATDATA_ExistScalar(data->flags[j]))
            {
              /* calculate scalar residual for this point */
              B_res[3] = data->F[j] - B_model[3];

              /* apply window function */
              B_res[3] *= window;

              /* recalculate scalar measurement */
              data->F[j] = B_res[3] + B_model[3];
            }
          else
            {
              fprintf(stderr, "polar_damp_apply: error: neither scalar nor vector data exist\n");
              exit(1);
            }
        }
    }

  return s;
}
