/*
 * track.c
 *
 * This module contains routines for:
 *
 * 1) Identifying north/south satellite tracks
 * 2) Computing residuals along the track with a main field model
 * 3) Computing the rms of the track residuals and throwing away tracks
 *    with high rms
 *
 * Calling sequence:
 * 1) track_alloc - allocate track workspace
 * 2) track_init  - parse satdata_mag data and separate into half-orbital tracks;
 *                  compute and store along-track residuals
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <errno.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_statistics.h>

#include <satdata/satdata.h>

#include <common/common.h>
#include <common/ema.h>
#include <common/interp.h>
#include <msynth/msynth.h>

#include "track.h"

static int track_calc_mf(const size_t start_idx, const size_t end_idx,
                         satdata_mag *data, msynth_workspace *w);
static int track_ema(const double alpha, const size_t flags[], double y[], const size_t n);
static int track_ema_reverse(const double alpha, const size_t flags[], double y[], const size_t n);

track_workspace *
track_alloc(void)
{
  track_workspace *w;

  w = calloc(1, sizeof(track_workspace));
  if (!w)
    return 0;

  w->ntot = TRACK_MAX;
  w->tracks = malloc(w->ntot * sizeof(track_data));
  if (!w->tracks)
    {
      track_free(w);
      return 0;
    }

  w->msynth_workspace_p = msynth_alloc(20, 1, NULL);

  w->n = 0;

  return w;
} /* track_alloc() */

void
track_free(track_workspace *w)
{
  size_t i;

  if (w->tracks)
    {
      for (i = 0; i < w->n; ++i)
        {
          track_data *tptr = &(w->tracks[i]);

          if (tptr->Bx)
            free(tptr->Bx);

          if (tptr->By)
            free(tptr->By);

          if (tptr->Bz)
            free(tptr->Bz);

          if (tptr->Bf)
            free(tptr->Bf);
        }

      free(w->tracks);
    }

  if (w->msynth_workspace_p)
    msynth_free(w->msynth_workspace_p);

  free(w);
} /* track_free() */

/*
track_init()
  Go through satellite data, separate data into half-orbital
tracks, store track information in internal structure, and
compute residuals. Points which do not fall into a track are flagged.

Inputs: data     - satellite data
        msynth_p - msynth workspace to compute along-track main
                   field (or NULL if already computed)
        w        - track workspace

Return: number of flagged points

Notes:
1) Tracks with no equator crossing have fields t_eq, lon_eq, lt_eq set to
-1.0e6
*/

size_t
track_init(satdata_mag *data, msynth_workspace *msynth_p, track_workspace *w)
{
  size_t nflagged = 0;
  size_t i = 0, j;

  w->data = data;

  while (i < data->n - 1)
    {
      int s;
      double lon_eq, lat_eq, t_eq, lt_eq;
      size_t sidx, eidx;
      track_data *tptr;

      s = satdata_mag_find_track(i, &eidx, data);
      if (s)
        {
#if TRACK_DEBUG
          fprintf(stderr, "track_init: no more tracks found in [%zu, %zu]\n", sidx, data->n - 1);
#endif
          break;
        }

      /* update index i for next loop */
      sidx = i;
      i = eidx + 1;

      /* check for equator crossing */
      if (data->latitude[sidx] * data->latitude[eidx] < 0.0)
        {
          size_t idx;
          time_t unix_time;

          /* there is an equator crossing, find time and longitude */

          idx = bsearch_double(data->latitude, 0.0, sidx, eidx);

          /* sanity check we found equator crossing */
          assert(data->latitude[idx] * data->latitude[idx + 1] <= 0.0);

          /* interpolate t, but not longitude due to wrapping effects */
          t_eq = interp1d(data->latitude[idx], data->latitude[idx + 1],
                          data->t[idx], data->t[idx + 1], 0.0);
          lon_eq = data->longitude[idx];
          lat_eq = data->latitude[idx];

          /* compute local time */
          unix_time = satdata_epoch2timet(t_eq);
          lt_eq = get_localtime(unix_time, lon_eq * M_PI / 180.0);
        }
      else
        {
          /* no equator crossing */

#if TRACK_DEBUG
          fprintf(stderr, "track_init: flagging partial track [%zu, %zu]\n", sidx, eidx);
#endif
          nflagged += track_flag_data(sidx, eidx, data);
          continue;
        }

      /* compute main field along track in data->F_main */
      if (msynth_p)
        track_calc_mf(sidx, eidx, data, msynth_p);

      /* store this track information */
      tptr = &(w->tracks[w->n]);
      tptr->start_idx = sidx;
      tptr->end_idx = eidx;
      tptr->n = eidx - sidx + 1;
      tptr->t_eq = t_eq;
      tptr->lon_eq = lon_eq;
      tptr->lat_eq = lat_eq;
      tptr->lt_eq = lt_eq;
      tptr->nrms_scal = 0;
      tptr->nrms_vec = 0;
      tptr->flags = 0;
      tptr->k_ext = 0.0;
      tptr->meanalt = gsl_stats_mean(&(data->altitude[sidx]), 1, tptr->n);

      /* use index of track center to compute satellite direction */
      tptr->satdir = satdata_mag_satdir(sidx + tptr->n / 2, data);

      for (j = 0; j < 3; ++j)
        tptr->rms[j] = 0.0;

      assert(tptr->n > 0);

      /* compute and store along-track residuals */
      tptr->Bx = malloc(tptr->n * sizeof(double));
      tptr->By = malloc(tptr->n * sizeof(double));
      tptr->Bz = malloc(tptr->n * sizeof(double));
      tptr->Bf = malloc(tptr->n * sizeof(double));

      track_calc_residuals(tptr, data);

      if (++(w->n) >= w->ntot)
        {
          fprintf(stderr, "track_init: ntot not large enough: %zu\n", w->ntot);
          return nflagged;
        }
    }

#if TRACK_DEBUG
  fprintf(stderr, "track_init: found %zu tracks\n", w->n);
#endif

  return nflagged;
} /* track_init() */

/*
track_calc_mf()
  Calculate main field model along satellite track

Inputs: start_idx - starting index in data
        end_idx   - end index in data
        data      - satellite data
        w         - workspace for main field model

Notes:
1) data->F_main and data->B_main are filled in with main field model
values in nT
*/

static int
track_calc_mf(const size_t start_idx, const size_t end_idx,
              satdata_mag *data, msynth_workspace *w)
{
  size_t i;

  for (i = start_idx; i <= end_idx; ++i)
    {
      double t = satdata_epoch2year(data->t[i]);
      double r = data->altitude[i] + data->R;
      double theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
      double phi = data->longitude[i] * M_PI / 180.0;
      double B_model[4];

      msynth_eval(t, r, theta, phi, B_model, w);

      SATDATA_VEC_X(data->B_main, i) = B_model[0];
      SATDATA_VEC_Y(data->B_main, i) = B_model[1];
      SATDATA_VEC_Z(data->B_main, i) = B_model[2];
      data->F_main[i] = B_model[3];
    }

  return GSL_SUCCESS;
} /* track_calc_mf() */

/*
track_flag_data()
  Flag all satellite data from [sidx,eidx] with SATDATA_FLG_OUTLIER
*/

size_t
track_flag_data(size_t sidx, size_t eidx, satdata_mag *data)
{
  size_t i;
  size_t n = 0;

  for (i = sidx; i <= eidx; ++i)
    {
      /* don't count if already flagged */
      if (data->flags[i] & SATDATA_FLG_OUTLIER)
        continue;

      data->flags[i] |= SATDATA_FLG_OUTLIER;
      ++n;
    }

  return n;
} /* track_flag_data() */

/*
track_flag_track()
  Flag entire track

Inputs: track_idx - track index
        flags     - TRACK_FLG_xxx
        data      - satellite data
        w         - track workspace

Return: number of data flagged

Notes:
1) satellite data is flagged with 'flags'
*/

size_t
track_flag_track(const size_t track_idx, const size_t flags,
                 satdata_mag *data, track_workspace *w)
{
  track_data *tptr = &(w->tracks[track_idx]);
  size_t nflag;

  nflag = track_flag_data(tptr->start_idx, tptr->end_idx, data);
  tptr->flags |= flags;

  return nflag;
} /* track_flag_track() */

int
track_residual(const size_t track_idx, const size_t data_idx, double B[3],
               const satdata_mag *data, track_workspace *w)
{
  int s = 0;
  size_t j;
  double B_obs[3], B_main[3], B_crust[3], B_ext[3];
  double B_ext_corr[3] = { 0.0, 0.0, 0.0 };

  B_obs[0] = SATDATA_VEC_X(data->B, data_idx);
  B_obs[1] = SATDATA_VEC_Y(data->B, data_idx);
  B_obs[2] = SATDATA_VEC_Z(data->B, data_idx);

  B_main[0] = SATDATA_VEC_X(data->B_main, data_idx);
  B_main[1] = SATDATA_VEC_Y(data->B_main, data_idx);
  B_main[2] = SATDATA_VEC_Z(data->B_main, data_idx);

  B_crust[0] = SATDATA_VEC_X(data->B_crust, data_idx);
  B_crust[1] = SATDATA_VEC_Y(data->B_crust, data_idx);
  B_crust[2] = SATDATA_VEC_Z(data->B_crust, data_idx);

  B_ext[0] = SATDATA_VEC_X(data->B_ext, data_idx);
  B_ext[1] = SATDATA_VEC_Y(data->B_ext, data_idx);
  B_ext[2] = SATDATA_VEC_Z(data->B_ext, data_idx);

#if 0
  if (tptr->k_ext != 0.0)
    {
      double t = satdata_epoch2year(data->t[data_idx]);
      double r = data->R + data->altitude[data_idx];
      double theta = M_PI / 2.0 - data->latitude[data_idx] * M_PI / 180.0;
      double phi = data->longitude[data_idx] * M_PI / 180.0;

      dstcorr_calc_model(tptr->k_ext, t, r, theta, phi, B_ext_corr, w->dstcorr_p);
    }
#endif

  for (j = 0; j < 3; ++j)
    B[j] = B_obs[j] - B_main[j] - B_crust[j] - B_ext[j] - B_ext_corr[j];

  return s;
} /* track_residual() */

/*
track_smooth()
  Smooth high-frequency high-latitude features above
55 degrees QD latitude with exponential moving average
filter

Inputs: alpha - filter constant between 0 and 1
                small values = more smoothing
                large values = less smoothing
        w     - w
*/

int
track_smooth(const double alpha, satdata_mag *data, track_workspace *w)
{
  int s = 0;
  size_t i, j;

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);
      size_t sidx = tptr->start_idx;
      size_t eidx = tptr->end_idx;
      size_t idxs, idxn;
      size_t npts;

      if (tptr->satdir == 1)
        {
          idxn = bsearch_double(&data->qdlat[sidx], 55.0, 0, tptr->n - 1);
          idxs = bsearch_double(&data->qdlat[sidx], -55.0, 0, tptr->n - 1);

          npts = tptr->n - idxn;
          track_ema(alpha, &data->flags[sidx + idxn], &tptr->Bf[idxn], npts);
          track_ema(alpha, &data->flags[sidx + idxn], &tptr->Bx[idxn], npts);
          track_ema(alpha, &data->flags[sidx + idxn], &tptr->By[idxn], npts);
          track_ema(alpha, &data->flags[sidx + idxn], &tptr->Bz[idxn], npts);

          npts = idxs + 1;
          track_ema_reverse(alpha, &data->flags[sidx], tptr->Bf, npts);
          track_ema_reverse(alpha, &data->flags[sidx], tptr->Bx, npts);
          track_ema_reverse(alpha, &data->flags[sidx], tptr->By, npts);
          track_ema_reverse(alpha, &data->flags[sidx], tptr->Bz, npts);
        }
      else
        {
          idxn = bsearch_double(&data->qdlat[sidx], 55.0, 0, tptr->n - 1);
          idxs = bsearch_double(&data->qdlat[sidx], -55.0, 0, tptr->n - 1);

          npts = idxn + 1;
          track_ema_reverse(alpha, &data->flags[sidx], tptr->Bf, npts);
          track_ema_reverse(alpha, &data->flags[sidx], tptr->Bx, npts);
          track_ema_reverse(alpha, &data->flags[sidx], tptr->By, npts);
          track_ema_reverse(alpha, &data->flags[sidx], tptr->Bz, npts);

          npts = tptr->n - idxs;
          track_ema(alpha, &data->flags[sidx + idxs], &tptr->Bf[idxs], npts);
          track_ema(alpha, &data->flags[sidx + idxs], &tptr->Bx[idxs], npts);
          track_ema(alpha, &data->flags[sidx + idxs], &tptr->By[idxs], npts);
          track_ema(alpha, &data->flags[sidx + idxs], &tptr->Bz[idxs], npts);
        }

      /* recompute data arrays with smoothed residuals */
      for (j = sidx; j <= eidx; ++j)
        {
          double B_model[4];

          satdata_mag_model(j, B_model, data);

          SATDATA_VEC_X(data->B, j) = tptr->Bx[j - sidx] + B_model[0];
          SATDATA_VEC_Y(data->B, j) = tptr->By[j - sidx] + B_model[1];
          SATDATA_VEC_Z(data->B, j) = tptr->Bz[j - sidx] + B_model[2];
          data->F[j] = tptr->Bf[j - sidx] + B_model[3];
        }

      /* recompute residuals with smoothed data */
      track_calc_residuals(tptr, data);
    }

  return s;
} /* track_smooth() */

/*
track_nflagged()
  Count number of flagged tracks
*/

size_t
track_nflagged(const track_workspace *w)
{
  size_t i;
  size_t nflagged = 0;

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);

      if (tptr->flags != 0)
        ++nflagged;
    }

  return nflagged;
} /* track_nflagged() */

/*
track_data_nflagged()
  Count number of flagged data in a single track
*/

size_t
track_data_nflagged(const track_data *tptr, const satdata_mag *data)
{
  size_t i;
  size_t nflagged = 0;

  for (i = 0; i < tptr->n; ++i)
    {
      size_t didx = i + tptr->start_idx;

      if (!SATDATA_AvailableData(data->flags[didx]))
        ++nflagged;
    }

  return nflagged;
} /* track_data_nflagged() */

/*
track_find()
  Find a track with an equator crossing within dt and dphi
of specified values

Inputs: t_eq   - timestamp of equator crossing (CDF_EPOCH)
        phi_eq - longitude of equator crossing (degrees)
        dt_min - allowed dt difference in minutes
        dphi   - allowed longitude difference in degrees
        idx    - (output) index of track if found
        w      - track workspace

Return: sucess if found, failure if not
*/

int
track_find(const double t_eq, const double phi_eq, const double dt_min,
           const double dphi, size_t *idx, const track_workspace *w)
{
  int s = GSL_FAILURE;
  size_t i;

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);
      double t_diff = fabs(t_eq - tptr->t_eq) / 60000.0; /* convert to minutes */
      double phi_diff = wrap180(phi_eq - tptr->lon_eq);

      if (fabs(phi_diff) > dphi)
        continue;

      if (t_diff > dt_min)
        continue;

      *idx = i;
      s = GSL_SUCCESS;
      break;
    }

  return s;
}

/*
track_find_t()
  Find track whose equator crossing is closest in time to a given timestamp

Inputs: t   - timestamp for comparison (CDF_EPOCH)
        idx - (output) index of track
        w   - track workspace
*/

int
track_find_t(const double t, size_t *idx, const track_workspace *w)
{
  int s = GSL_FAILURE;
  size_t i;
  double dt_min = 1.0e9;

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);
      double dt = fabs(t - tptr->t_eq);

      if (dt < dt_min)
        {
          dt_min = dt;
          *idx = i;
          s = GSL_SUCCESS;
        }
    }

  return s;
}

int
track_print(const char *filename, const size_t flags,
            const satdata_mag *data, track_workspace *w)
{
  int s = 0;
  size_t i, j;
  FILE *fp;
  size_t nflagged;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "track_print: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  nflagged = track_nflagged(w);
  fprintf(fp, "# Total tracks:    %zu\n", w->n);
  fprintf(fp, "# Total flagged:   %zu\n", nflagged);
  fprintf(fp, "# Total unflagged: %zu\n", w->n - nflagged);
  fprintf(fp, "# Track selection flags: %zu\n", flags);

  /* print header */
  track_print_track(1, fp, NULL, data);

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

      /* print this track */
      track_print_track(0, fp, tptr, data);
    }

  fclose(fp);

  return s;
} /* track_print() */

int
track_print_track(const int header, FILE *fp, const track_data *tptr,
                  const satdata_mag *data)
{
  int s = 0;
  size_t j;

  if (header)
    {
      j = 1;
      fprintf(fp, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", j++);
      fprintf(fp, "# Field %zu: local time (hours)\n", j++);
      fprintf(fp, "# Field %zu: local time at equator crossing (hours)\n", j++);
      fprintf(fp, "# Field %zu: season (doy)\n", j++);
      fprintf(fp, "# Field %zu: radius (km)\n", j++);
      fprintf(fp, "# Field %zu: longitude (degrees)\n", j++);
      fprintf(fp, "# Field %zu: geocentric latitude (degrees)\n", j++);
      fprintf(fp, "# Field %zu: Magnetic local time (hours)\n", j++);
      fprintf(fp, "# Field %zu: Dipole tilt angle (degrees)\n", j++);
      fprintf(fp, "# Field %zu: QD latitude (degrees)\n", j++);
      fprintf(fp, "# Field %zu: NEC X residual (nT)\n", j++);
      fprintf(fp, "# Field %zu: NEC Y residual (nT)\n", j++);
      fprintf(fp, "# Field %zu: NEC Z residual (nT)\n", j++);
      fprintf(fp, "# Field %zu: F residual (nT)\n", j++);
      fprintf(fp, "# Field %zu: NEC X measurement (nT)\n", j++);
      fprintf(fp, "# Field %zu: NEC Y measurement (nT)\n", j++);
      fprintf(fp, "# Field %zu: NEC Z measurement (nT)\n", j++);
      fprintf(fp, "# Field %zu: electron density data (cm^{-3})\n", j++);

      return s;
    }

  for (j = 0; j < tptr->n; ++j)
    {
      size_t didx = j + tptr->start_idx;
      time_t unix_time = satdata_epoch2timet(data->t[didx]);
      double ut = get_ut(unix_time);
      double lt = get_localtime(unix_time, data->longitude[didx] * M_PI / 180.0);

      if (data->flags[didx])
        continue;

      fprintf(fp, "%ld %7.4f %7.4f %5.1f %7.2f %7.3f %7.3f %7.4f %7.3f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %.5e\n",
              unix_time,
              lt,
              tptr->lt_eq,
              get_season(unix_time),
              data->r[didx],
              data->longitude[didx],
              data->latitude[didx],
              data->MLT[didx],
              data->tilt[didx] * 180.0 / M_PI,
              data->qdlat[didx],
              tptr->Bx[j],
              tptr->By[j],
              tptr->Bz[j],
              tptr->Bf[j],
              SATDATA_VEC_X(data->B, didx),
              SATDATA_VEC_Y(data->B, didx),
              SATDATA_VEC_Z(data->B, didx),
              data->ne[didx]);
    }

  fprintf(fp, "\n\n");

  return s;
}

/*
track_print_stats_flag()
  Print statistics on satellite tracks which match a given
flag

Inputs: filename - output file
        flag     - only tracks with given flag are printed;
                   if flag = 0, only unflagged tracks are printed
                   if flag = 0xffffffff, all tracks are printed
        w        - workspace
*/

int
track_print_stats_flag(const char *filename, const size_t flag,
                       track_workspace *w)
{
  int s = 0;
  size_t i;
  FILE *fp;
  satdata_mag *data = w->data;
  size_t nflagged;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "track_print_stats: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  nflagged = track_nflagged(w);
  fprintf(fp, "# Total tracks:    %zu\n", w->n);
  fprintf(fp, "# Total flagged:   %zu\n", nflagged);
  fprintf(fp, "# Total unflagged: %zu\n", w->n - nflagged);
  fprintf(fp, "# Track selection flags: %zu\n", flag);

  i = 1;
  fprintf(fp, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
  fprintf(fp, "# Field %zu: longitude of equator crossing (degrees)\n", i++);
  fprintf(fp, "# Field %zu: local time of equator crossing (hours)\n", i++);
  fprintf(fp, "# Field %zu: track mean altitude (km)\n", i++);
  fprintf(fp, "# Field %zu: number of data points in track\n", i++);
  fprintf(fp, "# Field %zu: number of non-flagged data points in track\n", i++);
  fprintf(fp, "# Field %zu: satellite direction (+1 north -1 south)\n", i++);
  fprintf(fp, "# Field %zu: track X rms (nT)\n", i++);
  fprintf(fp, "# Field %zu: track Y rms (nT)\n", i++);
  fprintf(fp, "# Field %zu: track Z rms (nT)\n", i++);
  fprintf(fp, "# Field %zu: track scalar rms (nT)\n", i++);
  fprintf(fp, "# Field %zu: number of along-track data used for scalar rms\n", i++);
  fprintf(fp, "# Field %zu: number of along-track data used for vector rms\n", i++);
  fprintf(fp, "# Field %zu: track flagged due to rms\n", i++);

  for (i = 0; i < w->n; ++i)
    {
      track_data *tptr = &(w->tracks[i]);
      size_t nflag = track_data_nflagged(tptr, data);

      if (flag != 0xffffffff)
        {
          if (flag == 0 && tptr->flags != 0)
            continue;

          if (flag != 0 && !(tptr->flags & flag))
            continue;
        }

      fprintf(fp, "%ld %10.4f %6.2f %6.2f %6zu %6zu %2d %10.4f %10.4f %10.4f %10.4f %6zu %6zu %d\n",
              satdata_epoch2timet(tptr->t_eq),
              wrap180(tptr->lon_eq),
              tptr->lt_eq,
              tptr->meanalt,
              tptr->n,
              tptr->n - nflag,
              tptr->satdir,
              tptr->rms[0],
              tptr->rms[1],
              tptr->rms[2],
              tptr->rms[3],
              tptr->nrms_scal,
              tptr->nrms_vec,
              (tptr->flags & TRACK_FLG_RMS) ? 1 : 0);
    }

  fclose(fp);

  return s;
} /* track_print_stats_flag() */

/*
track_print_stats()
  Print statistics on all satellite tracks

Inputs: filename - output file
        w        - workspace
*/

int
track_print_stats(const char *filename, track_workspace *w)
{
  return track_print_stats_flag(filename, 0xffffffff, w);
}

int
track_calc_residuals(track_data *tptr, const satdata_mag *data)
{
  size_t i;
  double B[4];

  for (i = tptr->start_idx; i <= tptr->end_idx; ++i)
    {
      size_t n = i - tptr->start_idx;

      satdata_mag_residual(i, B, data);

      tptr->Bx[n] = B[0];
      tptr->By[n] = B[1];
      tptr->Bz[n] = B[2];
      tptr->Bf[n] = B[3];
    }

  return 0;
} /* track_calc_residuals() */

static int
track_ema(const double alpha, const size_t flags[], double y[], const size_t n)
{
  size_t i;
  double S;
  int init = 0;
  
  for (i = 1; i < n; ++i)
    {
      /* don't use bad data in filter */
      if (SATDATA_BadData(flags[i]))
        continue;

      if (!init)
        {
          S = y[i];
          init = 1;
        }
      else
        S = alpha * y[i] + (1.0 - alpha) * S;

      y[i] = S;
    }

  return 0;
} /* track_ema() */

/*
track_ema_reverse()
  Compute EMA of array, starting from end of array and working
forward
*/

static int
track_ema_reverse(const double alpha, const size_t flags[], double y[], const size_t n)
{
  size_t i;
  double S;
  int init = 0;

  for (i = 2; i <= n; ++i)
    {
      size_t j = n - i;

      if (SATDATA_BadData(flags[j]))
        continue;

      if (!init)
        {
          S = y[j];
          init = 1;
        }
      else
        S = alpha * y[j] + (1.0 - alpha) * S;

      y[j] = S;
    }

  return 0;
} /* track_ema_reverse() */
