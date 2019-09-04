/*
 * mag_grad.c
 *
 * This module contains routines for inverting a magnetic profile
 * to an EEJ current profile using E/W gradient data
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <cdf.h>

#include <gsl/gsl_math.h>

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_bsearch.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_track.h>

#include "green.h"
#include "inverteef.h"
#include "log.h"
#include "mag.h"
#include "mag_grad.h"
#include "pde.h"

static int mag_track_datagap(const double dlat_max, const size_t sidx,
                             const size_t eidx, const satdata_mag *data);
static int mag_grad_compute_F1(const track_data *tptr, const satdata_mag *data,
                               const track_data *tptr2, const satdata_mag *data2,
                               mag_workspace *w);

/*
mag_grad_proc()
  Main processing loop for satellite data. Read in
tracks which cross the geographic equator and send them
to other modules for processing.

Inputs: data - satellite data
        w    - workspace
*/

int
mag_grad_proc(const mag_params *params, track_workspace *track_p, satdata_mag *data,
              track_workspace *track_p2, satdata_mag *data2, mag_workspace *w)
{
  int status = 0;
  const size_t ndata_min = 1000; /* minimum data points in track */
  const double lat_min = 35.0;   /* minimum latitude coverage needed */
  size_t track_idx, track_idx2;
  size_t ntrack = 0;     /* number of tracks processed */
  size_t nrejgap = 0,    /* number of rejections due to data gap */
         nrejlat = 0,    /* number of rejections due to latitude */
         nrejflag = 0;   /* number of rejections due to various flags */

  if (data->n < ndata_min || data2->n < ndata_min)
    {
      log_proc(w->log_general, "mag_proc: too few data for processing (%zu/%zu)\n",
               data->n, data2->n);
      return -1;
    }

  log_proc(w->log_general, "mag_proc: processing %zu satellite data (%zu tracks)\n",
           data->n, track_p->n);

  /* preprocess tracks using given parameters */
  mag_preproc(params, track_p, data, w);
  mag_preproc(params, track_p2, data2, w);

  for (track_idx = 0; track_idx < track_p->n; ++track_idx)
    {
      int s;
      track_data *tptr = &(track_p->tracks[track_idx]);
      size_t sidx = tptr->start_idx;
      size_t eidx = tptr->end_idx;
      double lon_eq = tptr->lon_eq;
      double t_eq = tptr->t_eq;
      double lt_eq = tptr->lt_eq;
      time_t unix_time = satdata_epoch2timet(tptr->t_eq);

      track_data *tptr2;
      time_t unix_time2;

      double kp;
      char buf[2048];

      /* discard flagged tracks */
      if (tptr->flags)
        {
          ++nrejflag;
          continue;
        }

      /* locate track of second satellite closest in time to the first */
      s = track_find_t(t_eq, &track_idx2, track_p2);
      if (s != GSL_SUCCESS)
        {
          fprintf(stderr, "mag_grad_proc: track for satellite 2 not found, timestamp: %ld\n", unix_time);
          continue;
        }

      tptr2 = &(track_p2->tracks[track_idx2]);
      unix_time2 = satdata_epoch2timet(tptr2->t_eq);

      /* check time difference between tracks */
      {
        double dt = fabs(tptr->t_eq - tptr2->t_eq) * 1.0e-3; /* time difference in seconds */
        if (dt > 60.0)
          {
            /* time difference too large */
            fprintf(stderr, "mag_grad_proc: track for satellite 2 not found due to dt [%.1f sec]\n", dt);
            continue;
          }
      }

      kp_get(unix_time, &kp, w->kp_workspace_p);

      /* reject tracks with insufficient latitude coverage */
      if (fabs(data->qdlat[sidx]) < lat_min || fabs(data->qdlat[eidx]) < lat_min ||
          fabs(data2->qdlat[tptr2->start_idx]) < lat_min || fabs(data2->qdlat[tptr2->end_idx]) < lat_min)
        {
          ++nrejlat;
          continue;
        }

      /* check for data gaps */
      s = mag_track_datagap(params->dlat_max, sidx, eidx, data);
      s += mag_track_datagap(params->dlat_max, tptr2->start_idx, tptr2->end_idx, data2);
      if (s)
        {
          ++nrejgap;
          continue;
        }

      ++ntrack;

      sprintf(buf, "%s", ctime(&unix_time));
      buf[strlen(buf) - 1] = '\0';
      fprintf(stderr, "mag_proc: [sat 1] track %zu, %s, lon = %g, lt = %g, kp = %g, dir = %d\n",
              ntrack, buf, lon_eq, lt_eq, kp, tptr->satdir);

      sprintf(buf, "%s", ctime(&unix_time2));
      buf[strlen(buf) - 1] = '\0';
      fprintf(stderr, "mag_proc: [sat 2] track %zu, %s, lon = %g, lt = %g, kp = %g, dir = %d\n",
              ntrack, buf, tptr2->lon_eq, tptr2->lt_eq, kp, tptr2->satdir);

      /*
       * store track in workspace, compute magnetic coordinates, and
       * F^(1) residuals
       */ 
      s = mag_grad_compute_F1(tptr, data, tptr2, data2, w);
      if (s)
        return s; /* error occurred */

      /* filter out Sq and compute B^(2) or F^(2) residuals */
      s = mag_sqfilt_vector(w, w->sqfilt_scalar_workspace_p);

      if (s)
        return s;

#if 0 /*XXX*/
      /* invert for EEJ height-integrated current density */
      s = mag_eej_proc(&(w->track), w->EEJ, w->eej_workspace_p);
      if (s)
        return s;
#endif

      /* log information for this profile (track number, t_eq, lon_eq) */
      mag_log_profile(0, ntrack, kp, tptr->satdir, w);

      /* print B^(2) residuals to log file */
      mag_log_B2_grad(0, w);

      /* print F^(2) residuals to log file */
      if (!params->use_vector)
        mag_log_F2(0, w);

      /* output Sq SVD, L-curve and corners */
      mag_log_Sq_Lcurve(0, w);
      mag_log_Sq_Lcorner(0, w);
      mag_log_Sq_svd(0, w);

      /* print line current profiles to log file */
      mag_log_LC(0, w);

      /* print EEJ current density to log file */
      mag_log_EEJ(0, w);

      /* output EEJ SVD, L-curve and corners */
      mag_log_EEJ_Lcurve(0, w);
      mag_log_EEJ_Lcorner(0, w);
      mag_log_EEJ_svd(0, w);

      /* stop processing here if we only want profiles */
      if (params->profiles_only)
        continue;

      /* solve PDE */
      s = pde_proc(unix_time, lon_eq * M_PI / 180.0, w->pde_workspace_p);

      /*
       * print PDE solution J(E,u=0) and J(E=0,u) to log file - solution
       * is printed even if pde solver fails to keep indexing
       * consistent in log files
       */
      mag_log_PDE(0, w);

      if (s)
        {
          log_proc(w->log_general, "mag_proc: pde_proc failed on profile %zu\n",
                   ntrack);
        }

      /* invert for EEF */
      {
        gsl_vector_view J_sat = gsl_vector_view_array(w->EEJ, w->ncurr);

        s += inverteef_calc(&J_sat.vector,
                            w->pde_workspace_p->theta_grid,
                            w->pde_workspace_p->J_lat_E,
                            w->pde_workspace_p->J_lat_u,
                            w->inverteef_workspace_p);

        /* compute final EEF value in mV/m */
        w->EEF = w->inverteef_workspace_p->E_scale * EEF_PHI_0 * 1.0e3;

        /* store relative error between modeled and satellite profiles */
        w->RelErr = w->inverteef_workspace_p->RelErr;

        log_proc(w->log_general,
                 "mag_proc: profile %zu: [t,LT,lon,EEF,RelErr,kp,dir] = [%d,%.2f,%7.2f,%6.3f,%.2f,%.1f,%d]\n",
                 ntrack,
                 unix_time,
                 lt_eq,
                 lon_eq,
                 w->EEF,
                 w->RelErr,
                 kp,
                 tptr->satdir);
        fprintf(stderr, "mag_proc: electric field = %f mV/m\n", w->EEF);
        fprintf(stderr, "mag_proc: relative error = %f\n", w->RelErr);

        /* log modeled profile */
        mag_log_model(0, w);

        /* if no errors, log EEF value to output file; since EEF is only printed
         * upon success, this file could be desynced from the other log files
         * which are all synced together for each profile, regardless of failure
         */
        if (s == GSL_SUCCESS)
          mag_log_EEF(0, unix_time, lon_eq * M_PI / 180.0, kp, w);
      }
    }

  log_proc(w->log_general,
           "mag_proc: %zu/%zu tracks rejected due to insufficient latitude coverage\n",
           nrejlat, track_p->n);
  log_proc(w->log_general,
           "mag_proc: %zu/%zu tracks rejected due to latitude gap\n",
           nrejgap, track_p->n);
  log_proc(w->log_general,
           "mag_proc: %zu/%zu tracks rejected due to flags (LT, kp, rms)\n",
           nrejflag, track_p->n);

  return status;
} /* mag_proc() */

/*
mag_track_datagap()
  Check satellite track for data gaps
*/

static int
mag_track_datagap(const double dlat_max,
                  const size_t sidx, const size_t eidx,
                  const satdata_mag *data)
{
  size_t i;

  for (i = sidx + 1; i <= eidx; ++i)
    {
      double prevlat = data->latitude[i - 1];
      double lat = data->latitude[i];

      if (fabs(lat - prevlat) > dlat_max)
        return -1;
    }

  return 0;
} /* mag_track_datagap() */

/*
mag_compute_F1()
  When a track is found for processing, compute magnetic
latitudes and copy the track into mag workspace for the
desired latitudes. Then call a main field model for the internal
and external field values along the track and compute F^(1)
residuals

Inputs: tptr  - track data for satellite 1
        data  - data for satellite 1
        tptr2 - track data for satellite 2
        data2 - data for satellite 2
        w      - workspace

Notes:
1) Upon exit, w->track contains all the information for this
track including the main field model values, QD latitudes and F^(1)
residuals
*/

static int
mag_grad_compute_F1(const track_data *tptr, const satdata_mag *data,
                    const track_data *tptr2, const satdata_mag *data2,
                    mag_workspace *w)
{
  int s = 0;
  const mag_params *params = w->params;
  const double dt_max = 5.0; /* maximum time difference allowed in sec */
  size_t i, j, k;
  size_t idx = 0;
  mag_track *track = &(w->track);

  track->t_eq = tptr->t_eq;
  track->phi_eq = tptr->lon_eq * M_PI / 180.0;

  for (i = tptr->start_idx; i <= tptr->end_idx; ++i)
    {
      double theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
      double phi = data->longitude[i] * M_PI / 180.0;
      double qdlat = data->qdlat[i];
      double B_int[4], B_crust[4], B_ext[4], B_tot[4];
      double B_int_grad[4], B_crust_grad[4], B_ext_grad[4], B_tot_grad[4];

      if (fabs(qdlat) > params->track_qdmax)
        continue;

      /* find measurement for satellite 2 with approximately the same time */
      j = bsearch_double(data2->t, data->t[i], tptr2->start_idx, tptr2->end_idx);
      if (fabs(data->t[i] - data2->t[j]) > dt_max * 1.0e3)
        continue;

      /* store field models for satellite 1 */
      B_int[0] = SATDATA_VEC_X(data->B_main, i);
      B_int[1] = SATDATA_VEC_Y(data->B_main, i);
      B_int[2] = SATDATA_VEC_Z(data->B_main, i);

      B_crust[0] = SATDATA_VEC_X(data->B_crust, i);
      B_crust[1] = SATDATA_VEC_Y(data->B_crust, i);
      B_crust[2] = SATDATA_VEC_Z(data->B_crust, i);

      B_ext[0] = SATDATA_VEC_X(data->B_ext, i);
      B_ext[1] = SATDATA_VEC_Y(data->B_ext, i);
      B_ext[2] = SATDATA_VEC_Z(data->B_ext, i);

      /* store field models for satellite 2 */
      B_int_grad[0] = SATDATA_VEC_X(data2->B_main, j);
      B_int_grad[1] = SATDATA_VEC_Y(data2->B_main, j);
      B_int_grad[2] = SATDATA_VEC_Z(data2->B_main, j);

      B_crust_grad[0] = SATDATA_VEC_X(data2->B_crust, j);
      B_crust_grad[1] = SATDATA_VEC_Y(data2->B_crust, j);
      B_crust_grad[2] = SATDATA_VEC_Z(data2->B_crust, j);

      B_ext_grad[0] = SATDATA_VEC_X(data2->B_ext, j);
      B_ext_grad[1] = SATDATA_VEC_Y(data2->B_ext, j);
      B_ext_grad[2] = SATDATA_VEC_Z(data2->B_ext, j);

      /* add crustal field to main field */
      for (k = 0; k < 3; ++k)
        {
          B_tot[k] = B_int[k] + B_crust[k] + B_ext[k];
          B_tot_grad[k] = B_int_grad[k] + B_crust_grad[k] + B_ext_grad[k];

          B_int[k] += B_crust[k];
          B_int_grad[k] += B_crust_grad[k];
        }

      B_int[3] = gsl_hypot3(B_int[0], B_int[1], B_int[2]);
      B_tot[3] = gsl_hypot3(B_tot[0], B_tot[1], B_tot[2]);

      B_int_grad[3] = gsl_hypot3(B_int_grad[0], B_int_grad[1], B_int_grad[2]);
      B_tot_grad[3] = gsl_hypot3(B_tot_grad[0], B_tot_grad[1], B_tot_grad[2]);

      /* store this track */
      track->t[idx] = data->t[i];
      track->theta[idx] = theta;
      track->lat_deg[idx] = 90.0 - theta*180.0/M_PI;
      track->phi[idx] = phi;
      track->r[idx] = data->r[i];
      track->thetaq[idx] = M_PI / 2.0 - qdlat * M_PI / 180.0;
      track->qdlat[idx] = qdlat;
      track->F[idx] = data->F[i];
      track->X[idx] = SATDATA_VEC_X(data->B, i);
      track->Y[idx] = SATDATA_VEC_Y(data->B, i);
      track->Z[idx] = SATDATA_VEC_Z(data->B, i);
      track->Bx_int[idx] = B_int[0];
      track->By_int[idx] = B_int[1];
      track->Bz_int[idx] = B_int[2];
      track->F_int[idx] = B_int[3];
      track->dF_ext[idx] = vec_dot(B_int, B_ext) / B_int[3];

      track->t_grad[idx] = data2->t[j];
      track->theta_grad[idx] = M_PI / 2.0 - data2->latitude[j] * M_PI / 180.0;
      track->lat_deg_grad[idx] = data2->latitude[j];
      track->phi_grad[idx] = data2->longitude[j] * M_PI / 180.0;
      track->r_grad[idx] = data2->r[j];
      track->thetaq_grad[idx] = M_PI / 2.0 - data2->qdlat[j] * M_PI / 180.0;
      track->qdlat_grad[idx] = data2->qdlat[j];
      track->F_grad[idx] = data2->F[j];
      track->X_grad[idx] = SATDATA_VEC_X(data2->B, j);
      track->Y_grad[idx] = SATDATA_VEC_Y(data2->B, j);
      track->Z_grad[idx] = SATDATA_VEC_Z(data2->B, j);
      track->Bx_int_grad[idx] = B_int_grad[0];
      track->By_int_grad[idx] = B_int_grad[1];
      track->Bz_int_grad[idx] = B_int_grad[2];
      track->F_int_grad[idx] = B_int_grad[3];
      track->dF_ext_grad[idx] = vec_dot(B_int_grad, B_ext_grad) / B_int_grad[3];

      /* compute F^(1) residuals */
      track->F1[idx] = data->F[i] - B_tot[3];
      track->X1[idx] = SATDATA_VEC_X(data->B, i) - B_tot[0];
      track->Y1[idx] = SATDATA_VEC_Y(data->B, i) - B_tot[1];
      track->Z1[idx] = SATDATA_VEC_Z(data->B, i) - B_tot[2];

      track->F1_grad[idx] = data2->F[j] - B_tot_grad[3];
      track->X1_grad[idx] = SATDATA_VEC_X(data2->B, j) - B_tot_grad[0];
      track->Y1_grad[idx] = SATDATA_VEC_Y(data2->B, j) - B_tot_grad[1];
      track->Z1_grad[idx] = SATDATA_VEC_Z(data2->B, j) - B_tot_grad[2];

      if (++idx >= MAG_MAX_TRACK)
        {
          fprintf(stderr, "mag_compute_F1: MAG_MAX_TRACK too small\n");
          return -1;
        }
    }

  track->n = idx;

  return s;
}
