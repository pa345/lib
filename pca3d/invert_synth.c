/*
 * invert_synth.c
 *
 * Routines for handling synthetic test case
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_test.h>

#include <mainlib/ml_common.h>
#include <mainlib/ml_msynth.h>

#include "invert.h"
#include "invert_synth.h"

static int invert_synth_calc(const int thread_id, const double t, const double r, const double theta, const double phi,
                             const gsl_vector * c, double B[3], invert_workspace *w);

/* fill in coefficient vector of synthetic coefficients */
int
invert_synth_g(gsl_vector * c, invert_workspace * w)
{
  gsl_vector_set_zero(c);

  return 0;
}

/* replace with synthetic data for testing */
int
invert_synth_replace(invert_workspace *w)
{
  gsl_vector * c = gsl_vector_alloc(w->p);
  size_t i, j;

  for (i = 0; i < w->nsat; ++i)
    {
      magdata *mptr = invert_data_ptr(i, w->data_workspace_p);

#pragma omp parallel for private(j)
      for (j = 0; j < mptr->n; ++j)
        {
          int thread_id = omp_get_thread_num();
          double B[3];

          /* synthesize magnetic field vector */
          invert_synth_calc(thread_id, mptr->t[j], mptr->r[j], mptr->theta[j], mptr->phi[j], c, B, w);

          mptr->Bx_nec[j] = B[0];
          mptr->By_nec[j] = B[1];
          mptr->Bz_nec[j] = B[2];
          mptr->F[j] = gsl_hypot3(B[0], B[1], B[2]);

          mptr->Bx_model[j] = 0.0;
          mptr->By_model[j] = 0.0;
          mptr->Bz_model[j] = 0.0;
        }
    }

  gsl_vector_free(c);

  return 0;
}

/* compute B(r,theta,phi) using coefficients 'c' */
static int
invert_synth_calc(const int thread_id, const double t, const double r, const double theta, const double phi, const gsl_vector * c,
                  double B[3], invert_workspace *w)
{
  int s = 0;
  size_t mode_t, mode_s, band;
  gsl_complex alpha;
  gsl_complex Phi[3], B_complex[3];
  size_t i;

  for (i = 0; i < 3; ++i)
    B[i] = 0.0;

  invert_smode_precompute(thread_id, r, theta, phi, w->smode_workspace_p);

#if 1
  /* 1 cpd band */
  mode_t = 1;
  mode_s = 0;
  band = 10;
  invert_smode_get(thread_id, r, theta, phi, band, mode_s, Phi, w->smode_workspace_p);
  alpha = invert_tmode_get(t, band, mode_t, w->tmode_workspace_p);

  for (i = 0; i < 3; ++i)
    {
      B_complex[i] = gsl_complex_mul(alpha, Phi[i]);
      B[i] += GSL_REAL(B_complex[i]);
    }
#endif

#if 1
  /* 2 cpd band */
  mode_t = 5;
  mode_s = 4;
  band = 8;
  invert_smode_get(thread_id, r, theta, phi, band, mode_s, Phi, w->smode_workspace_p);
  alpha = invert_tmode_get(t, band, mode_t, w->tmode_workspace_p);

  for (i = 0; i < 3; ++i)
    {
      B_complex[i] = gsl_complex_mul(alpha, Phi[i]);
      B[i] += GSL_REAL(B_complex[i]);
    }
#endif

  return s;
}
