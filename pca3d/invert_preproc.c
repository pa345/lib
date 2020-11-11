/*
 * invert_preproc.c
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
#include <omp.h>
#include <libconfig.h>

#include <flow/flow.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>

#include <mainlib/ml_apex.h>
#include <mainlib/ml_satdata.h>
#include <mainlib/ml_indices.h>
#include <mainlib/ml_obsdata.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_solarpos.h>
#include <mainlib/ml_msynth.h>
#include <mainlib/ml_track.h>
#include <mainlib/ml_magdata.h>

#if 0
#include "mfield.h"
#include "mfield_data.h"
#endif

typedef struct
{
  size_t downsample;      /* downsampling factor */
  double min_LT;          /* minimum local time for field modeling */
  double max_LT;          /* maximum local time for field modeling */
  double qdlat_preproc_cutoff; /* QD latitude cutoff for high-latitudes */
  double min_IMF_Bz;      /* minimum IMF Bz (nT) */
  double max_IMF_Bz;      /* maximum IMF Bz (nT) */
  double min_IMF_By;      /* minimum IMF By (nT) */
  double max_IMF_By;      /* maximum IMF By (nT) */
  size_t gradient_ns;     /* number of seconds between N/S gradient samples */

  char * rms_file;        /* along-track rms output file */
  double rms_thresh[4];   /* rms thresholds (X,Y,Z,F) (nT) */

  double gradew_dphi_max; /* maximum longitude distance for east-west gradients (degrees) */
  double gradew_dlat_max; /* maximum latitude distance for east-west gradients (degrees) */
  double gradew_dt_max;   /* maximum time difference for east-west gradients (seconds) */

  int subtract_B_main;    /* subtract a-priori main field from data */
  int subtract_B_crust;   /* subtract a-priori crustal field from data */
  int subtract_B_ext;     /* subtract a-priori external field from data */

  double max_kp;          /* maximum kp */
  double max_dRC;         /* maximum dRC/dt (nT/hour) */
} preprocess_parameters;

static int invert_check_LT(const double lt, const double lt_min, const double lt_max);
void print_unflagged_data(const char *filename, const satdata_mag *data);

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

#include "invert_preproc_filter.c"

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

  if (params->qdlat_preproc_cutoff < 0.0)
    {
      fprintf(stderr, "check_parameters: qdlat_preproc_cutoff must be > 0\n");
      ++s;
    }

  if (params->gradient_ns == 0)
    {
      fprintf(stderr, "check_parameters: gradient_ns must be > 0\n");
      ++s;
    }

  if (params->gradew_dphi_max <= 0.0)
    {
      fprintf(stderr, "check_parameters: gradient_ew_dphi_max must be > 0\n");
      ++s;
    }

  if (params->gradew_dlat_max <= 0.0)
    {
      fprintf(stderr, "check_parameters: gradient_ew_dlat_max must be > 0\n");
      ++s;
    }

  if (params->gradew_dt_max <= 0.0)
    {
      fprintf(stderr, "check_parameters: gradient_ew_dt_max must be > 0\n");
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
copy_data(const size_t magdata_flags, const satdata_mag *data, const track_workspace *track_p,
          const size_t magdata_flags2, const satdata_mag *data2, const track_workspace *track_p2,
          preprocess_parameters * preproc_params, magdata * mdata)
{
  const size_t nflagged = satdata_nflagged(data);
  size_t ndata = data->n - nflagged;
  magdata_params params;
  size_t i;
  size_t npts[6] = { 0, 0, 0, 0, 0, 0 };
  size_t nmodel[MFIELD_IDX_END];

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
  for (i = 0; i < MFIELD_IDX_END; ++i)
    nmodel[i] = 0;

  if (mdata == NULL)
    mdata = magdata_alloc(ndata, data->R);
  else
    mdata = magdata_realloc(ndata + mdata->n, data->R, mdata);

  if (!mdata)
    return 0;

  mdata->global_flags = magdata_flags;
  
  fprintf(stderr, "\n");
  fprintf(stderr, "\t copy_data: copying tracks in magdata format...");

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
           * to fit field model
           */

          double LT = mdata->lt_eq[i];
          int fit_MF = invert_check_LT(LT, preproc_params->min_LT, preproc_params->max_LT);

          if (fit_MF)
            fitting_flags |= MAGDATA_FLG_FIT_MF;
        }
      else
        {
          /* high-latitude point - fit field model but not Euler angles */
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

          /* count north-south gradient data */

          if (MAGDATA_ExistDX_NS(mdata->flags[i]))
            ++(nmodel[MFIELD_IDX_DX_NS]);

          if (MAGDATA_ExistDY_NS(mdata->flags[i]))
            ++(nmodel[MFIELD_IDX_DY_NS]);

          if (MAGDATA_ExistDZ_NS(mdata->flags[i]))
            ++(nmodel[MFIELD_IDX_DZ_NS]);

          if (MAGDATA_ExistDF_NS(mdata->flags[i]))
            ++(nmodel[MFIELD_IDX_DF_NS]);

          /* count east-west gradient data */

          if (MAGDATA_ExistDX_EW(mdata->flags[i]))
            ++(nmodel[MFIELD_IDX_DX_EW]);

          if (MAGDATA_ExistDY_EW(mdata->flags[i]))
            ++(nmodel[MFIELD_IDX_DY_EW]);

          if (MAGDATA_ExistDZ_EW(mdata->flags[i]))
            ++(nmodel[MFIELD_IDX_DZ_EW]);

          if (MAGDATA_ExistDF_EW(mdata->flags[i]))
            ++(nmodel[MFIELD_IDX_DF_EW]);
        }

      mdata->flags[i] |= fitting_flags;
    }

  fprintf(stderr, "done\n");

  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) scalar measurements available\n",
          npts[0], mdata->n, (double) npts[0] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) vector measurements available\n",
          npts[1], mdata->n, (double) npts[1] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) along-track scalar measurements available\n",
          npts[2], mdata->n, (double) npts[2] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) along-track vector measurements available\n",
          npts[3], mdata->n, (double) npts[3] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) east-west scalar measurements available\n",
          npts[4], mdata->n, (double) npts[4] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) east-west vector measurements available\n",
          npts[5], mdata->n, (double) npts[5] / (double) mdata->n * 100.0);

  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) X vector measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_X], mdata->n, (double) nmodel[MFIELD_IDX_X] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) Y vector measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_Y], mdata->n, (double) nmodel[MFIELD_IDX_Y] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) Z vector measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_Z], mdata->n, (double) nmodel[MFIELD_IDX_Z] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) scalar measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_F], mdata->n, (double) nmodel[MFIELD_IDX_F] / (double) mdata->n * 100.0);

  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) North/South dX vector measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_DX_NS], mdata->n, (double) nmodel[MFIELD_IDX_DX_NS] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) North/South dY vector measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_DY_NS], mdata->n, (double) nmodel[MFIELD_IDX_DY_NS] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) North/South dZ vector measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_DZ_NS], mdata->n, (double) nmodel[MFIELD_IDX_DZ_NS] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) North/South dF scalar measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_DF_NS], mdata->n, (double) nmodel[MFIELD_IDX_DF_NS] / (double) mdata->n * 100.0);

  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) East/West dX vector measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_DX_EW], mdata->n, (double) nmodel[MFIELD_IDX_DX_EW] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) East/West dY vector measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_DY_EW], mdata->n, (double) nmodel[MFIELD_IDX_DY_EW] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) East/West dZ vector measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_DZ_EW], mdata->n, (double) nmodel[MFIELD_IDX_DZ_EW] / (double) mdata->n * 100.0);
  fprintf(stderr, "\t copy_data: %zu/%zu (%.1f%%) East/West dF scalar measurements selected for MF modeling\n",
          nmodel[MFIELD_IDX_DF_EW], mdata->n, (double) nmodel[MFIELD_IDX_DF_EW] / (double) mdata->n * 100.0);

  return mdata;
}

satdata_mag *
read_swarm(const char *filename, const int asmv)
{
  size_t nflag;
  satdata_mag *data;
  struct timeval tv0, tv1;

  fprintf(stderr, "read_swarm: reading %s...", filename);
  gettimeofday(&tv0, NULL);

  if (asmv)
    data = satdata_swarm_asmv_read_idx(filename, 0);
  else
    data = satdata_swarm_read_idx(filename, 1);

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
  satdata_mag *data;
  struct timeval tv0, tv1;
  size_t nflag;

  fprintf(stderr, "read_champ: reading %s...", filename);
  gettimeofday(&tv0, NULL);
  data = satdata_swarm_read_idx(filename, 0);
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
invert_check_LT()
  Check if a given LT is within [lt_min,lt_max] accounting for mod 24. So it
is possible to have input lt_min < lt_max in order to select data across midnight.

Example: [lt_min,lt_max] = [6,18] will select daytime data between 6am and 6pm
         [lt_min,lt_max] = [18,6] will select nighttime data between 6pm and 6am
         [lt_min,lt_max] = [22,5] will select nighttime data between 10pm and 5am
         [lt_min,lt_max] = [0,5] will select data between midnight and 5am

Return: 1 if LT \in [lt_min,lt_max] (mod 24); 0 otherwise
*/

static int
invert_check_LT(const double lt, const double lt_min, const double lt_max)
{
  double a, b;

  if (lt_min < 0.0 && lt_max < 0.0)
    return 1; /* disable LT check */

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
preprocess_data(const preprocess_parameters *params, const size_t magdata_flags,
                satdata_mag *data)
{
  struct timeval tv0, tv1;
  track_workspace *track_p = track_alloc();

  fprintf(stderr, "preprocess_data: initializing tracks...");
  gettimeofday(&tv0, NULL);
  track_init(data, NULL, track_p);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  /* detect bad tracks with rms test */
  {
    size_t nrms;

    nrms = track_flag_rms(params->rms_file, params->rms_thresh, NULL, data, track_p);
    fprintf(stderr, "preprocess_data: flagged (%zu/%zu) (%.1f%%) tracks due to high rms\n",
            nrms, track_p->n, (double) nrms / (double) track_p->n * 100.0);
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
  print_unflagged_data("data.dat.1", data);
#endif

  /* select geomagnetically quiet data */
  fprintf(stderr, "preprocess_data: selecting geomagnetically quiet data...");
  invert_preprocess_filter(magdata_flags, params, track_p, data);
  fprintf(stderr, "done\n");

#if 0
  print_unflagged_data("data.dat.2", data);
#endif

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

/* print unflagged data points for data selection visualization */
void
print_unflagged_data(const char *filename, const satdata_mag *data)
{
  size_t i;
  FILE *fp;

  fp = fopen(filename, "w");

  for (i = 0; i < data->n; ++i)
    {
      double B[3], B_main[3], B_ext[3];

#if 0
      if (!SATDATA_AvailableData(data->flags[i]))
        continue;
#else
      if (SATDATA_BadData(data->flags[i]) || (data->flags[i] & SATDATA_FLG_FILTER))
        continue;
#endif

      B[0] = SATDATA_VEC_X(data->B, i);
      B[1] = SATDATA_VEC_Y(data->B, i);
      B[2] = SATDATA_VEC_Z(data->B, i);

      B_main[0] = SATDATA_VEC_X(data->B_main, i);
      B_main[1] = SATDATA_VEC_Y(data->B_main, i);
      B_main[2] = SATDATA_VEC_Z(data->B_main, i);

      B_ext[0] = SATDATA_VEC_X(data->B_ext, i);
      B_ext[1] = SATDATA_VEC_Y(data->B_ext, i);
      B_ext[2] = SATDATA_VEC_Z(data->B_ext, i);

      fprintf(fp, "%f %f %f %f %f %f %f\n",
              data->qdlat[i],
              B[0],
              B[1],
              B[2],
              B_main[0] + B_ext[0],
              B_main[1] + B_ext[1],
              B_main[2] + B_ext[2]);
    }

  fclose(fp);
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

#if 0 /*XXX*/

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

#endif

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
  if (config_lookup_float(&cfg, "min_IMF_Bz", &fval))
    params->min_IMF_Bz = fval;
  if (config_lookup_float(&cfg, "min_IMF_By", &fval))
    params->min_IMF_By = fval;

  if (config_lookup_int(&cfg, "downsample", &ival))
    params->downsample = (size_t) ival;
  if (config_lookup_int(&cfg, "gradient_ns", &ival))
    params->gradient_ns = (size_t) ival;

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

  config_destroy(&cfg);

  return 0;
}

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options]\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --champ_file      | -c champ_index_file       - CHAMP index file\n");
  fprintf(stderr, "\t --swarm_file      | -s swarm_index_file       - Swarm index file\n");
  fprintf(stderr, "\t --swarm_file2     | -t swarm_index_file2      - Swarm index file 2 (for E/W gradients)\n");
  fprintf(stderr, "\t --swarm_asmv_file | -a swarm_asmv_index_file  - Swarm ASM-V index file\n");
  fprintf(stderr, "\t --oersted_f_file  | -L oersted_f_index_file   - Oersted scalar (f) index file\n");
  fprintf(stderr, "\t --oersted_l_file  | -L oersted_l_index_file   - Oersted vector (l) index file\n");
  fprintf(stderr, "\t --cryosat_file    | -r cryosat_index_file     - Cryosat index file\n");
  fprintf(stderr, "\t --obs_data_file   | -O obs_data_file          - Observatory data file (BGS format)\n");
  fprintf(stderr, "\t --downsample      | -d downsample             - downsampling factor\n");
  fprintf(stderr, "\t --gradient_ns     | -g num_samples            - number of samples between N/S gradient points\n");
  fprintf(stderr, "\t --output_file     | -o output_file            - binary output data file (magdata format)\n");
  fprintf(stderr, "\t --config_file     | -C config_file            - configuration file\n");
  fprintf(stderr, "\t --append_file     | -A append_file            - new data will be appended to this dataset in magdata binary format\n");
  fprintf(stderr, "\t --rms_file        | -R rms_file               - output along-track rms file\n");
}

int
main(int argc, char *argv[])
{
  int status;
  char *datamap_file = "datamap.dat";
  char *data_file = "data.dat";
  char *output_file = NULL;
  char *append_file = NULL;
  char *config_file = "INVERT_preproc.cfg";
  char *rms_file = "satrms.txt";
  char *name = NULL;
  satdata_mag *data = NULL;
  satdata_mag *data2 = NULL;
  obsdata *observatory_data = NULL;
  magdata *mdata;
  struct timeval tv0, tv1;
  track_workspace *track_p = NULL;
  track_workspace *track_p2 = NULL;
  preprocess_parameters params;
  size_t downsample = 0;     /* downsample factor */
  size_t gradient_ns = 0;    /* number of seconds between N/S gradient points */
  size_t magdata_flags = 0;       /* MAGDATA_GLOBFLG_xxx */
  size_t magdata_flags2 = 0;
  int flag_vec_rms = 1;
  size_t nflag;

  /* initialize parameters */
  params.downsample = 0;
  params.max_kp = -1.0;
  params.max_dRC = -1.0;
  params.min_LT = -1.0;
  params.max_LT = -1.0;
  params.qdlat_preproc_cutoff = -1.0;
  params.min_IMF_Bz = -1.0e6;
  params.max_IMF_Bz =  1.0e6;
  params.min_IMF_By = -1.0e6;
  params.max_IMF_By =  1.0e6;
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
  params.rms_file = rms_file;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "swarm_file", required_argument, NULL, 's' },
          { "swarm_file2", required_argument, NULL, 't' },
          { "swarm_asmv_file", required_argument, NULL, 'a' },
          { "oersted_f_file", required_argument, NULL, 'F' },
          { "oersted_l_file", required_argument, NULL, 'L' },
          { "champ_file", required_argument, NULL, 'c' },
          { "cryosat_file", required_argument, NULL, 'r' },
          { "obs_file", required_argument, NULL, 'O' },
          { "downsample", required_argument, NULL, 'd' },
          { "output_file", required_argument, NULL, 'o' },
          { "config_file", required_argument, NULL, 'C' },
          { "gradient_ns", required_argument, NULL, 'g' },
          { "append_file", required_argument, NULL, 'A' },
          { "name", required_argument, NULL, 'N' },
          { "rms_file", required_argument, NULL, 'R' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "a:A:c:C:d:D:F:g:L:N:o:O:r:R:s:t:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          /* Swarm ASM-V */
          case 'a':
            data = read_swarm(optarg, 0);
            magdata_flags = MAGDATA_GLOBFLG_SWARM;
            flag_vec_rms = 0; /* no NEC data for rms flagging */
            break;

          /* Swarm official */
          case 's':
            data = read_swarm(optarg, 0);
            magdata_flags = MAGDATA_GLOBFLG_SWARM;
            break;

          /* For E/W gradients */
          case 't':
            data2 = read_swarm(optarg, 0);
            magdata_flags2 = 0;
            break;

          case 'c':
            magdata_flags = MAGDATA_GLOBFLG_CHAMP;
            data = read_champ(optarg);
            break;

          case 'r':
            magdata_flags = MAGDATA_GLOBFLG_CRYOSAT;
            data = read_swarm(optarg, 0);
            break;

          case 'F':
          case 'L':
            magdata_flags = MAGDATA_GLOBFLG_OERSTED;
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);
            data = satdata_swarm_read_idx(optarg, 1);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu data read, %g seconds)\n",
                    data->n, time_diff(tv0, tv1));

            fprintf(stderr, "main: filtering for instrument flags...");
            nflag = satdata_oersted_filter_instrument(1, data);
            fprintf(stderr, "done (%zu/%zu (%.1f%%) data flagged)\n",
                    nflag, data->n, (double)nflag / (double)data->n * 100.0);

            break;

          case 'C':
            config_file = optarg;
            break;

          case 'd':
            downsample = (size_t) atoi(optarg);
            break;

          case 'g':
            gradient_ns = (size_t) atoi(optarg);
            break;

          case 'o':
            output_file = optarg;
            break;

          case 'O':
            fprintf(stderr, "main: reading %s...", optarg);
            observatory_data = obsdata_swarm_read_idx(optarg, 1);
            fprintf(stderr, "done (%zu stations read, %zu total measurements)\n",
                    observatory_data->nstation, obsdata_n(observatory_data));
            break;

          case 'A':
            append_file = optarg;
            break;

          case 'N':
            name = optarg;
            break;

          case 'R':
            params.rms_file = optarg;
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
  if (gradient_ns > 0)
    params.gradient_ns = gradient_ns;
  if (flag_vec_rms == 0)
    {
      params.rms_thresh[0] = -1.0;
      params.rms_thresh[1] = -1.0;
      params.rms_thresh[2] = -1.0;
    }

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
  fprintf(stderr, "main: gradient_ns      = %zu [s]\n", params.gradient_ns);
  fprintf(stderr, "main: rms X threshold  = %.1f [nT]\n", params.rms_thresh[0]);
  fprintf(stderr, "main: rms Y threshold  = %.1f [nT]\n", params.rms_thresh[1]);
  fprintf(stderr, "main: rms Z threshold  = %.1f [nT]\n", params.rms_thresh[2]);
  fprintf(stderr, "main: rms F threshold  = %.1f [nT]\n", params.rms_thresh[3]);
  fprintf(stderr, "main: LT minimum       = %.1f\n", params.min_LT);
  fprintf(stderr, "main: LT maximum       = %.1f\n", params.max_LT);

  fprintf(stderr, "main: === PREPROCESSING SATELLITE 1 ===\n");
  track_p = preprocess_data(&params, magdata_flags, data);

  if (data2)
    {
      fprintf(stderr, "main: === PREPROCESSING SATELLITE 2 ===\n");
      track_p2 = preprocess_data(&params, magdata_flags2, data2);
    }

  if (append_file)
    {
      fprintf(stderr, "main: reading append file %s...", append_file);
      mdata = magdata_read(append_file, NULL);
      fprintf(stderr, "done (%zu data read)\n", mdata->n);
    }

  fprintf(stderr, "main: computing vector residuals...");
  gettimeofday(&tv0, NULL);
  mdata = copy_data(magdata_flags, data, track_p, magdata_flags2, data2, track_p2, &params, mdata);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu residuals copied, %g seconds)\n",
          mdata->n, time_diff(tv0, tv1));

  /* free data after copying arrays to free up memory */
  satdata_mag_free(data);

  if (data2)
    satdata_mag_free(data2);

  magdata_init(mdata);

  fprintf(stderr, "main: computing spatial weighting of data...");
  gettimeofday(&tv0, NULL);
  magdata_calc(mdata);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  /* set dataset name */
  if (name != NULL)
    strncpy(mdata->name, name, MAGDATA_NAME_LENGTH);

#if 0
  fprintf(stderr, "main: writing data to %s...", data_file);
  magdata_print(data_file, mdata);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: writing data map to %s...", datamap_file);
  magdata_map(datamap_file, mdata);
  fprintf(stderr, "done\n");
#endif

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
  solarpos_free(solarpos_workspace_p);

  if (track_p2)
    track_free(track_p2);

  return 0;
} /* main() */
