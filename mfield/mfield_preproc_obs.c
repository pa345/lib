/*
 * mfield_preproc_obs.c
 *
 * Pre-process observatory data and save in magdata format.
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
#include <omp.h>
#include <libconfig.h>

#include <apex/apex.h>
#include <flow/flow.h>
#include <indices/indices.h>
#include <obsdata/obsdata.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>

#include <common/common.h>
#include <common/solarpos.h>
#include <msynth/msynth.h>

#include "magdata.h"
#include "magfit.h"
#include "track.h"

#include "mfield.h"
#include "mfield_data.h"

typedef struct
{
  size_t downsample;      /* downsampling factor */
  double min_LT;          /* minimum local time for field modeling */
  double max_LT;          /* maximum local time for field modeling */
  double qdlat_preproc_cutoff; /* QD latitude cutoff for high-latitudes */
  double min_zenith;      /* minimum zenith angle for high-latitude data selection */

  int subtract_B_main;    /* subtract a-priori main field from data */
  int subtract_B_crust;   /* subtract a-priori crustal field from data */
  int subtract_B_ext;     /* subtract a-priori external field from data */

  double max_kp;          /* maximum kp */
  double max_dRC;         /* maximum dRC/dt (nT/hour) */
} preprocess_parameters;

static int mfield_check_LT(const double lt, const double lt_min, const double lt_max);
static size_t model_flags(const size_t magdata_flags, const double t,
                          const double theta, const double phi, const double qdlat,
                          const preprocess_parameters * params);

#define MFIELD_IDX_X              0
#define MFIELD_IDX_Y              1
#define MFIELD_IDX_Z              2
#define MFIELD_IDX_F              3
#define MFIELD_IDX_DX_NS          4
#define MFIELD_IDX_DY_NS          5
#define MFIELD_IDX_DZ_NS          6
#define MFIELD_IDX_DF_NS          7
#define MFIELD_IDX_DX_EW          8
#define MFIELD_IDX_DY_EW          9
#define MFIELD_IDX_DZ_EW          10
#define MFIELD_IDX_DF_EW          11
#define MFIELD_IDX_B_EULER        12
#define MFIELD_IDX_END            13

/* Global */
solarpos_workspace *solarpos_workspace_p = NULL;

#if 0
#include "mfield_preproc_filter.c"
#endif

static int
check_parameters(preprocess_parameters * params)
{
  int s = 0;

  if (params->downsample == 0)
    {
      fprintf(stderr, "check_parameters: downsample must be > 0\n");
      ++s;
    }

  if (params->max_kp <= 0.0)
    {
      fprintf(stderr, "check_parameters: max_kp must be > 0\n");
      ++s;
    }

  if (params->max_dRC <= 0.0)
    {
      fprintf(stderr, "check_parameters: max_dRC must be > 0\n");
      ++s;
    }

  if (params->min_LT < 0.0)
    {
      fprintf(stderr, "check_parameters: min_LT must be > 0\n");
      ++s;
    }

  if (params->max_LT < 0.0)
    {
      fprintf(stderr, "check_parameters: max_LT must be > 0\n");
      ++s;
    }

  if (params->qdlat_preproc_cutoff < 0.0)
    {
      fprintf(stderr, "check_parameters: qdlat_preproc_cutoff must be > 0\n");
      ++s;
    }

  if (params->min_zenith < 0.0)
    {
      fprintf(stderr, "check_parameters: min_zenith must be > 0\n");
      ++s;
    }

  if (params->subtract_B_main < 0)
    {
      fprintf(stderr, "check_parameters: subtract_B_main must be 0 or 1\n");
      ++s;
    }

  if (params->subtract_B_crust < 0)
    {
      fprintf(stderr, "check_parameters: subtract_B_crust must be 0 or 1\n");
      ++s;
    }

  if (params->subtract_B_ext < 0)
    {
      fprintf(stderr, "check_parameters: subtract_B_ext must be 0 or 1\n");
      ++s;
    }

  return s;
}

magdata *
copy_data(const size_t magdata_flags, const obsdata_station *station, preprocess_parameters * preproc_params)
{
#if 1

  size_t ndata = station->n_sv;
  magdata *mdata;
  magdata_params params;
  size_t npts[6] = { 0, 0, 0, 0, 0, 0 };
  size_t i;

  params.model_main = 0;
  params.model_crust = 0;
  params.model_ext = 0;

  mdata = magdata_alloc(ndata, R_EARTH_KM);
  if (!mdata)
    return 0;

  mdata->global_flags = magdata_flags;

  magdata_copy_station(&params, station, mdata, npts);

  /* now determine which points in mdata will be used for MF modeling */
  for (i = 0; i < mdata->n; ++i)
    {
      size_t fitting_flags = 0;

      if (fabs(mdata->qdlat[i]) <= preproc_params->qdlat_preproc_cutoff)
        {
          /* mid-latitude point: check local time to determine whether to fit MF model */

          double LT = mdata->lt[i];
          int fit_MF = mfield_check_LT(LT, preproc_params->min_LT, preproc_params->max_LT);

          if (fit_MF)
            fitting_flags |= MAGDATA_FLG_FIT_MF;
        }
      else
        {
          /* high-latitude point - fit field model */
          fitting_flags |= MAGDATA_FLG_FIT_MF;
        }

      mdata->flags[i] |= fitting_flags;
    }

  return mdata;

#elif 0
  const size_t nflagged = obsdata_station_nflagged(station);
  size_t ndata = station->n - nflagged;
  magdata *mdata;
  magdata_params params;
  size_t i;
  size_t npts[6] = { 0, 0, 0, 0, 0, 0 };
  size_t nmodel[MFIELD_IDX_END];

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
  for (i = 0; i < MFIELD_IDX_END; ++i)
    nmodel[i] = 0;

  mdata = magdata_alloc(ndata, R_EARTH_KM);
  if (!mdata)
    return 0;

  mdata->global_flags = magdata_flags;
  
  fprintf(stderr, "\n");
  fprintf(stderr, "\t copy_data: converting station %s to magdata format...", station->name);

#if 0
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
#endif

  fprintf(stderr, "done (ndata = %zu mdata_n = %zu, mdata_ntot = %zu)\n", ndata, mdata->n, mdata->ntot);

  fprintf(stderr, "\t copy_data: flagging data for MF, Euler, etc...");

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
          int fit_MF = mfield_check_LT(LT, preproc_params->min_LT, preproc_params->max_LT);

          if (fit_MF)
            fitting_flags |= MAGDATA_FLG_FIT_MF;
        }
      else
        {
          /* high-latitude point - fit field model */
          fitting_flags |= MAGDATA_FLG_FIT_MF;
        }

      if (fitting_flags & MAGDATA_FLG_FIT_MF)
        {
          if (MAGDATA_ExistX(mdata->flags[i]))
            ++(nmodel[MFIELD_IDX_X]);

          if (MAGDATA_ExistY(mdata->flags[i]))
            ++(nmodel[MFIELD_IDX_Y]);

          if (MAGDATA_ExistZ(mdata->flags[i]))
            ++(nmodel[MFIELD_IDX_Z]);

          if (MAGDATA_ExistScalar(mdata->flags[i]))
            ++(nmodel[MFIELD_IDX_F]);
        }

      mdata->flags[i] |= fitting_flags;
    }

  fprintf(stderr, "done\n");

  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) scalar measurements available\n",
          npts[0], mdata->n, (double) npts[0] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) vector measurements available\n",
          npts[1], mdata->n, (double) npts[1] / (double) mdata->n * 100.0);

  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) X vector measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_X], mdata->n, (double) nmodel[MFIELD_IDX_X] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) Y vector measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_Y], mdata->n, (double) nmodel[MFIELD_IDX_Y] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) Z vector measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_Z], mdata->n, (double) nmodel[MFIELD_IDX_Z] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) scalar measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_F], mdata->n, (double) nmodel[MFIELD_IDX_F] / (double) mdata->n * 100.0);

  return mdata;
#else
  return NULL;
#endif
}

/*
mfield_check_LT()
  Check if a given LT is within [lt_min,lt_max] accounting for mod 24. So it
is possible to have input lt_min < lt_max in order to select data across midnight.

Example: [lt_min,lt_max] = [6,18] will select daytime data between 6am and 6pm
         [lt_min,lt_max] = [18,6] will select nighttime data between 6pm and 6am
         [lt_min,lt_max] = [22,5] will select nighttime data between 10pm and 5am
         [lt_min,lt_max] = [0,5] will select data between midnight and 5am

Return: 1 if LT \in [lt_min,lt_max] (mod 24); 0 otherwise
*/

static int
mfield_check_LT(const double lt, const double lt_min, const double lt_max)
{
  double a, b;

  b = fmod(lt_max - lt_min, 24.0);
  if (b < 0.0)
    b += 24.0;

  a = fmod(lt - lt_min, 24.0);
  if (a < 0.0)
    a += 24.0;

  if (a > b)
    return 0; /* invalid local time */

  /* valid local time */
  return 1;
}

#if 0

/*
model_flags()
  Check an individual data point to determine if it will be used to
fit various model parameters.

Inputs: magdata_flags - MAGDATA_GLOBFLG_xxx
        t             - CDF_EPOCH timestamp
        theta         - colatitude (radians)
        phi           - longitude (radians)
        qdlat         - QD latitude (degrees)
        params        - preprocess parameters

Return: flags indicating fit parameters (MAGDATA_FLG_FIT_xxx)

Notes:
1) If data point is below MFIELD_HIGH_LATITUDE and local time is
within [params->min_LT,params->max_LT], flag is set to MAGDATA_FLG_FIT_MF

2) If data point is higher than MFIELD_HIGH_LATITUDE, zenith angle
is computed and if the point is in darkness, MAGDATA_FLG_FIT_MF is set

3) If we are fitting Euler angles to this satellite, and the data
point satisfies the criteria, MAGDATA_FLG_FIT_EULER is set
*/

static size_t
model_flags(const size_t magdata_flags, const double t,
            const double theta, const double phi, const double qdlat,
            const preprocess_parameters * params)
{
  size_t flags = 0;
  const time_t unix_time = satdata_epoch2timet(t);
  const double lt = get_localtime(unix_time, phi);
  const double lat_deg = 90.0 - theta * 180.0 / M_PI;
  int status;

  /* check if we should fit Euler angles to this data point */
  status = mfield_check_LT(lt, MFIELD_EULER_LT_MIN, MFIELD_EULER_LT_MAX);
  if ((status == 1) &&
      (magdata_flags & MAGDATA_GLOBFLG_EULER) &&
      (fabs(qdlat) <= MFIELD_EULER_QDLAT))
    {
      flags |= MAGDATA_FLG_FIT_EULER;
    }

  /* check if we should fit main field model to this data point */
  if (fabs(lat_deg) <= MFIELD_HIGH_LATITUDE)
    {
      status = mfield_check_LT(lt, params->min_LT, params->max_LT);

      if (status == 1)
        flags |= MAGDATA_FLG_FIT_MF;
    }
  else
    {
      double lat_rad = lat_deg * M_PI / 180.0;
      double zenith;

      solarpos_calc_zenith(unix_time, lat_rad, phi, &zenith, solarpos_workspace_p);
      zenith *= 180.0 / M_PI;
      assert(zenith >= 0.0);

      /* large zenith angle means darkness */
      if (zenith >= MFIELD_MAX_ZENITH)
        flags |= MAGDATA_FLG_FIT_MF;
    }

  return flags;
}

#endif

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
        data   - observatory data

Return: pointer to sorted track workspace (should be freed by caller)
*/

int
preprocess_data(const preprocess_parameters *params, const size_t magdata_flags,
                obsdata *data)
{
  const double interval_SV = 0.5 * 365.25; /* interval for SV differences (days) */
  const size_t min_samples = 1;            /* minimum samples needed for daily means */
  obsdata_select_params obs_params = obsdata_select_default_params();
  struct timeval tv0, tv1;
  size_t i;

  fprintf(stderr, "preprocess_data: selecting quiet time data...");

  for (i = 0; i < data->nstation; ++i)
    {
      obsdata_station * station = data->stations[i];
      size_t nflagged_array[OBSDATA_IDX_END];
      size_t nflagged;

      fprintf(stderr, "\t selecting %s observatory data for quiet conditions...", station->name);
      gettimeofday(&tv0, NULL);
      nflagged = obsdata_station_select(&obs_params, station, data, nflagged_array);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds, %zu/%zu (%.1f%%) data flagged)\n",
              time_diff(tv0, tv1), nflagged, station->n, (double) nflagged / (double) station->n * 100.0);

      fprintf(stderr, "\t computing %s daily means...", station->name);
      gettimeofday(&tv0, NULL);
      obsdata_station_daily(station);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds, %zu means computed)\n",
              time_diff(tv0, tv1), station->n_mean_tot);

      fprintf(stderr, "\t computing %s SV values...", station->name);
      gettimeofday(&tv0, NULL);
      obsdata_station_SV(min_samples, interval_SV, station);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds, %zu SV values computed)\n",
              time_diff(tv0, tv1), station->n_sv);
    }

  fprintf(stderr, "done\n");

#if 0
  /* select geomagnetically quiet data */
  fprintf(stderr, "preprocess_data: selecting geomagnetically quiet data...");
  mfield_preprocess_filter(magdata_flags, params, track_p, data);
  fprintf(stderr, "done\n");
#endif

#if 0
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
#endif

  return 0;
}

int
calc_main(satdata_mag *data)
{
  const int max_threads = omp_get_max_threads();
  msynth_workspace **msynth_p;
  size_t i;

  msynth_p = malloc(max_threads * sizeof(msynth_workspace *));

  for (i = 0; i < (size_t) max_threads; ++i)
    {
      msynth_p[i] = msynth_read(MSYNTH_BOUMME_FILE);
      msynth_set(1, 15, msynth_p[i]);
    }

#pragma omp parallel for private(i)
  for (i = 0; i < data->n; ++i)
    {
      int thread_id = omp_get_thread_num();
      double tyr = satdata_epoch2year(data->t[i]);
      double r = data->r[i];
      double theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
      double phi = data->longitude[i] * M_PI / 180.0;
      double B_core[4];

      /* IMPORTANT: this needs the second check below, since AvailableData() will
       * reject downsampled points, but they could be used for N/S or E/W gradients
       * so the field value must be computed for thos too
       */
#if 0
      if (!SATDATA_AvailableData(data->flags[i]))
        continue;
#else
      if (SATDATA_BadData(data->flags[i]) || (data->flags[i] & SATDATA_FLG_FILTER))
        continue;
#endif

      msynth_eval(tyr, r, theta, phi, B_core, msynth_p[thread_id]);

      SATDATA_VEC_X(data->B_main, i) = B_core[0];
      SATDATA_VEC_Y(data->B_main, i) = B_core[1];
      SATDATA_VEC_Z(data->B_main, i) = B_core[2];
    }

  for (i = 0; i < (size_t) max_threads; ++i)
    msynth_free(msynth_p[i]);

  free(msynth_p);

  return 0;
}

static int
subtract_RC(const char *filename, satdata_mag *data, track_workspace *w)
{
  int s = 0;
  const magfit_type *T = magfit_rc;
  magfit_parameters magfit_params = magfit_default_parameters();
  magfit_workspace *magfit_p;
  size_t i, j;
  FILE *fp;

  magfit_params.rc_p = 1;
  magfit_params.rc_fit_Y = 0;
  magfit_params.rc_subtract_crust = 1;
  magfit_p = magfit_alloc(T, &magfit_params);

  fp = fopen(filename, "w");
  magfit_print_track(1, fp, NULL, data, magfit_p);

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);
      double rnorm, snorm;

      if (tptr->flags != 0)
        continue;

      magfit_reset(magfit_p);

      for (j = 0; j < tptr->n; ++j)
        {
          size_t didx = j + tptr->start_idx;
          double t = data->t[didx];
          double r = data->r[didx];
          double theta = M_PI / 2.0 - data->latitude[didx] * M_PI / 180.0;
          double phi = data->longitude[didx] * M_PI / 180.0;
          double qdlat = data->qdlat[didx];
          double B[3];

          if (SATDATA_BadData(data->flags[didx]) || (data->flags[didx] & SATDATA_FLG_FILTER))
            continue;

          /* only fit RC model to low-latitude data */
          if (fabs(qdlat) > 55.0)
            continue;

          /* start with total measurement */
          B[0] = SATDATA_VEC_X(data->B, didx);
          B[1] = SATDATA_VEC_Y(data->B, didx);
          B[2] = SATDATA_VEC_Z(data->B, didx);

          /* subtract main field */
          B[0] -= SATDATA_VEC_X(data->B_main, didx);
          B[1] -= SATDATA_VEC_Y(data->B_main, didx);
          B[2] -= SATDATA_VEC_Z(data->B_main, didx);

          /* subtract crustal field */
          B[0] -= SATDATA_VEC_X(data->B_crust, didx);
          B[1] -= SATDATA_VEC_Y(data->B_crust, didx);
          B[2] -= SATDATA_VEC_Z(data->B_crust, didx);

          /* subtract external field */
          B[0] -= SATDATA_VEC_X(data->B_ext, didx);
          B[1] -= SATDATA_VEC_Y(data->B_ext, didx);
          B[2] -= SATDATA_VEC_Z(data->B_ext, didx);

          /* add residual to magfit workspace */
          magfit_add_datum(t, r, theta, phi, qdlat, B, magfit_p);
        }

      /* fit RC model */
      s = magfit_fit(&rnorm, &snorm, magfit_p);
      if (s)
        continue;

      magfit_print_track(0, fp, tptr, data, magfit_p);

      /* now add the RC model to the external field model vector */
      for (j = 0; j < tptr->n; ++j)
        {
          size_t didx = j + tptr->start_idx;
          double t = data->t[didx];
          double r = data->r[didx];
          double theta = M_PI / 2.0 - data->latitude[didx] * M_PI / 180.0;
          double phi = data->longitude[didx] * M_PI / 180.0;
          double B[3];

          if (SATDATA_BadData(data->flags[didx]) || (data->flags[didx] & SATDATA_FLG_FILTER))
            continue;

          magfit_eval_B(t, r, theta, phi, B, magfit_p);

          SATDATA_VEC_X(data->B_ext, didx) += B[0];
          SATDATA_VEC_Y(data->B_ext, didx) += B[1];
          SATDATA_VEC_Z(data->B_ext, didx) += B[2];
        }
    }

  magfit_free(magfit_p);
  fclose(fp);

  return s;
}

static int
parse_config_file(const char *filename, preprocess_parameters *params)
{
  int s;
  config_t cfg;
  double fval;
  int ival;

  config_init(&cfg);

  s = config_read_file(&cfg, filename);
  if (s != CONFIG_TRUE)
    {
      fprintf(stderr, "parse_config_file: %s:%d - %s\n",
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
  if (config_lookup_float(&cfg, "qdlat_preproc_cutoff", &fval))
    params->qdlat_preproc_cutoff = fval;
  if (config_lookup_float(&cfg, "min_zenith", &fval))
    params->min_zenith = fval;

  if (config_lookup_int(&cfg, "downsample", &ival))
    params->downsample = (size_t) ival;

  if (config_lookup_int(&cfg, "subtract_B_main", &ival))
    params->subtract_B_main = (size_t) ival;
  if (config_lookup_int(&cfg, "subtract_B_crust", &ival))
    params->subtract_B_crust = (size_t) ival;
  if (config_lookup_int(&cfg, "subtract_B_ext", &ival))
    params->subtract_B_ext = (size_t) ival;

  config_destroy(&cfg);

  return 0;
}

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options]\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --obs_data_file   | -O                        - Observatory data file (%s)\n", OBSDATA_BINARY_FILE);
  fprintf(stderr, "\t --downsample      | -d downsample             - downsampling factor\n");
  fprintf(stderr, "\t --config_file     | -C config_file            - configuration file\n");
}

int
main(int argc, char *argv[])
{
  int status;
  const char *path_dir = "/data/palken/lib/mfield/data/obs";
  char *datamap_file = "datamap.dat";
  char *data_file = "data.dat";
  char *config_file = "OBS.cfg";
  obsdata *data = NULL;
  struct timeval tv0, tv1;
  preprocess_parameters params;
  size_t downsample = 0;     /* downsample factor */
  size_t magdata_flags = 0;  /* MAGDATA_GLOBFLG_xxx */
  int flag_vec_rms = 1;
  char output_file[2048];
  size_t i;

  /* initialize parameters */
  params.downsample = 0;
  params.max_kp = -1.0;
  params.max_dRC = -1.0;
  params.min_LT = -1.0;
  params.max_LT = -1.0;
  params.qdlat_preproc_cutoff = -1.0;
  params.min_zenith = -1.0;
  params.subtract_B_main = -1;
  params.subtract_B_crust = -1;
  params.subtract_B_ext = -1;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "obs_file", no_argument, NULL, 'O' },
          { "downsample", required_argument, NULL, 'd' },
          { "config_file", required_argument, NULL, 'C' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "C:d:O", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'C':
            config_file = optarg;
            break;

          case 'd':
            downsample = (size_t) atoi(optarg);
            break;

          case 'O':
            fprintf(stderr, "main: reading %s...", OBSDATA_BINARY_FILE);
            gettimeofday(&tv0, NULL);
            data = obsdata_read(OBSDATA_BINARY_FILE);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%g seconds, %zu stations read, %zu total measurements)\n",
                    time_diff(tv0, tv1), data->nstation, obsdata_n(data));
            break;

          default:
            break;
        }
    }

  /* parse configuration file */
  parse_config_file(config_file, &params);

  /* replace config values with command-line arguments */
  if (downsample > 0)
    params.downsample = downsample;

  /* check parameters */
  status = check_parameters(&params);
  if (status)
    exit(1);

  if (!data)
    {
      print_help(argv);
      exit(1);
    }

  solarpos_workspace_p = solarpos_alloc();

  fprintf(stderr, "main: downsample       = %zu\n", params.downsample);
  fprintf(stderr, "main: LT minimum       = %.1f\n", params.min_LT);
  fprintf(stderr, "main: LT maximum       = %.1f\n", params.max_LT);

  fprintf(stderr, "main: === PREPROCESSING OBSERVATORY DATA ===\n");
  preprocess_data(&params, magdata_flags, data);

  for (i = 0; i < data->nstation; ++i)
    {
      obsdata_station *station = data->stations[i];
      magdata *mdata = copy_data(magdata_flags, station, &params);

      sprintf(output_file, "%s/%s.dat", path_dir, station->name);
      fprintf(stderr, "main: writing data to %s...", output_file);
      magdata_write(output_file, mdata);
      fprintf(stderr, "done\n");

      magdata_free(mdata);
    }

  /* free data after copying arrays to free up memory */
  obsdata_free(data);

  exit(1);

#if 0
  magdata_init(mdata);

  fprintf(stderr, "main: computing spatial weighting of data...");
  gettimeofday(&tv0, NULL);
  magdata_calc(mdata);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

#if 0
  fprintf(stderr, "main: writing data to %s...", data_file);
  magdata_print(data_file, mdata);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: writing data map to %s...", datamap_file);
  magdata_map(datamap_file, mdata);
  fprintf(stderr, "done\n");
#endif

  fprintf(stderr, "main: observatory rmin = %.1f (%.1f) [km]\n",
          mdata->rmin, mdata->rmin - mdata->R);
  fprintf(stderr, "main: observatory rmax = %.1f (%.1f) [km]\n",
          mdata->rmax, mdata->rmax - mdata->R);

  if (output_file)
    {
      fprintf(stderr, "main: writing data to %s...", output_file);
      magdata_write(output_file, mdata);
      fprintf(stderr, "done\n");
    }

  magdata_free(mdata);
  solarpos_free(solarpos_workspace_p);
#endif

  return 0;
}
