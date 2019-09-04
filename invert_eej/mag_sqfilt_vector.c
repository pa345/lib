/*
 * mag_sqfilt_vector.c
 *
 * Fit spherical harmonic model to satellite scalar or vector residuals
 * to remove Sq and external fields
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_test.h>

#include <mainlib/ml_coord.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_interp.h>

#include "mag.h"

#include "Gdef.h"

#include "gauss.c"

static int sqfilt_calc_B2(mag_workspace *w, gauss_state_t * gauss_p);

/* turn on/off different components in fit */
#define FIT_X                      0
#define FIT_Y                      0
#define FIT_Z                      1

/*
mag_sqfilt_vector()
  Fit a simple spherical harmonic model to a single satellite
track of scalar B^(1) residual data, excluding the EEJ at low-latitudes.

Notes:
1) On output, w->rnorm and w->snorm contain the
residual norm and solution norm

2) On output, mag_p->track.{X,Y,Z}2 contains the B^(2) residuals and
mag_p->track.Sq_model contains the (M + K) Sq model

3) If a datagap is detected, this function returns -1
*/

int
mag_sqfilt_vector(mag_workspace *mag_p, mag_sqfilt_scalar_workspace *w)
{
  int s = 0;
  mag_track *track = &(mag_p->track);
  mag_params *params = mag_p->params;
  const double exclude_qdlat = params->sq_qdmin; /* exclude points below this QD latitude from fit */
  const size_t ntot = track->n;
  size_t i;
  double qd_min, qd_max;
  double B1_min[3], B1_max[3];
  double rnorm, snorm;
  gauss_state_t *gauss_p;
  gauss_parameters gauss_params;

  gauss_params.nmax_int = mag_p->params->sq_nmax_int;
  gauss_params.mmax_int = mag_p->params->sq_mmax_int;
  gauss_params.nmax_ext = mag_p->params->sq_nmax_ext;
  gauss_params.mmax_ext = mag_p->params->sq_mmax_ext;
  gauss_params.flags = 0;
#if FIT_X
  gauss_params.flags |= GAUSS_FLG_FIT_X;
#endif
#if FIT_Y
  gauss_params.flags |= GAUSS_FLG_FIT_Y;
#endif
#if FIT_Z
  gauss_params.flags |= GAUSS_FLG_FIT_Z;
#endif

  gauss_p = gauss_alloc(&gauss_params);

  gauss_reset(gauss_p);

  /* search qdlat array for -exclude_qdlat and store values */
  i = bsearch_double(track->qdlat, -exclude_qdlat, 0, track->n - 1);
  qd_min = track->qdlat[i];
  B1_min[0] = track->X1[i];
  B1_min[1] = track->Y1[i];
  B1_min[2] = track->Z1[i];

  /* search qdlat array for exclude_qdlat and store values */
  i = bsearch_double(track->qdlat, exclude_qdlat, 0, track->n - 1);
  qd_max = track->qdlat[i];
  B1_max[0] = track->X1[i];
  B1_max[1] = track->Y1[i];
  B1_max[2] = track->Z1[i];

  /* check for data gap */
  if ((fabs(qd_min + exclude_qdlat) > 3.0) ||
      (fabs(qd_max - exclude_qdlat) > 3.0))
    {
      fprintf(stderr, "mag_sqfilt_vector: large data gap found, rejecting track\n");
      return -1;
    }

  /* build least squares matrix and RHS */
  for (i = 0; i < ntot; ++i)
    {
      double rhs_B[3];

      if (fabs(track->qdlat[i]) > params->sq_qdmax)
        continue;

      rhs_B[0] = track->X1[i];
      rhs_B[1] = track->Y1[i];
      rhs_B[2] = track->Z1[i];

      /* exclude points in the equatorial electrojet region */
      if (fabs(track->qdlat[i]) < exclude_qdlat)
        {
          /* keep every 10th point in the equatorial region */
          if (i % 10 == 0)
            {
              /* interpolate linearly across EEJ region */
              rhs_B[0] = interp1d(qd_min, qd_max, B1_min[0], B1_max[0], track->qdlat[i]);
              rhs_B[1] = interp1d(qd_min, qd_max, B1_min[1], B1_max[1], track->qdlat[i]);
              rhs_B[2] = interp1d(qd_min, qd_max, B1_min[2], B1_max[2], track->qdlat[i]);
            }
          else
            continue;
        }

      gauss_add_datum(track->t[i], track->r[i], track->theta[i], track->phi[i],
                      track->qdlat[i], rhs_B, gauss_p);
    }

  /* fit the Sq model */
  gauss_fit(&rnorm, &snorm, gauss_p);

  /* calculate B^(2) residuals (eq 9 of paper) */
  sqfilt_calc_B2(mag_p, gauss_p);

  gauss_free(gauss_p);

  return s;
} /* mag_sqfilt_vector() */

/****************************************************
 *       INTERNAL ROUTINES                          *
 ****************************************************/

/*
sqfilt_calc_B2()
  Calculate B^(2) residual (eq. 9 of paper) by subtracting
the model (M + K)

Inputs: c - model coefficients
        w - workspace

Notes:
1) On output, w->track.{X,Y,Z}2 is modified to contain B^(2) residuals
*/

static int
sqfilt_calc_B2(mag_workspace *w, gauss_state_t * gauss_p)
{
  size_t i;
  mag_track *track = &(w->track);

  for (i = 0; i < track->n; ++i)
    {
      double B_int[3], B_ext[3];

      gauss_eval_B_int(track->t[i], track->r[i], track->theta[i], track->phi[i],
                       B_int, gauss_p);
      gauss_eval_B_ext(track->t[i], track->r[i], track->theta[i], track->phi[i],
                       B_ext, gauss_p);

      track->X_Sq_int[i] = B_int[0];
      track->Y_Sq_int[i] = B_int[1];
      track->Z_Sq_int[i] = B_int[2];

      track->X_Sq_ext[i] = B_ext[0];
      track->Y_Sq_ext[i] = B_ext[1];
      track->Z_Sq_ext[i] = B_ext[2];

      track->X2[i] = track->X1[i] - B_int[0] - B_ext[0];
      track->Y2[i] = track->Y1[i] - B_int[1] - B_ext[1];
      track->Z2[i] = track->Z1[i] - B_int[2] - B_ext[2];
    }

  return GSL_SUCCESS;
}
