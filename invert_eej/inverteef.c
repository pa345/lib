/*
 * inverteef.c
 *
 * Steps:
 *
 * 1. Interpolate PDE solutions J_lat_E and J_lat_u to the same grid
 *    used for line currents
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <errno.h>
#include <string.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_spline.h>

#include <mainlib/ml_indices.h>
#include <mainlib/ml_common.h>

#include "inverteef.h"

#include "lse.c"

static inline size_t inverteef_thidx(double theta, inverteef_workspace *w);
static int inverteef_interp_pdesol(const gsl_vector *J_sat, const double *theta_pde,
                                   const gsl_vector *J_lat_E, const gsl_vector *J_lat_u,
                                   inverteef_workspace *w);
static int inverteef_design_matrix(inverteef_workspace *w);
static int inverteef_invert_profile(gsl_vector *J_sat,
                                    inverteef_workspace *w);

inverteef_workspace *
inverteef_alloc(inverteef_parameters *params)
{
  const gsl_interp_type * T = gsl_interp_linear;
  inverteef_workspace *w;
  size_t ncoeffs = 0;

  w = calloc(1, sizeof(inverteef_workspace));
  if (!w)
    {
      fprintf(stderr, "inverteef_alloc: calloc failed: %s\n", strerror(errno));
      return 0;
    }

  w->J_pde = gsl_vector_alloc(params->ncurr);
  w->J_rhs = gsl_vector_alloc(params->ncurr);
  w->J_pde_E = gsl_vector_alloc(params->ncurr);
  w->J_pde_u = gsl_vector_alloc(params->ncurr);
  w->J_diff = gsl_vector_alloc(params->ncurr);
  if (w->J_rhs == 0 || w->J_pde_E == 0 || w->J_pde_u == 0 || w->J_pde == 0 ||
      w->J_diff == 0)
    {
      inverteef_free(w);
      return 0;
    }

#ifdef PDE_INVERT_E_FIELD
  ++ncoeffs;
#endif

#ifdef PDE_ALLOW_DC_SHIFT
  ++ncoeffs;
#endif

  w->X = gsl_matrix_alloc(params->ncurr, ncoeffs);
  w->coeffs = gsl_vector_alloc(ncoeffs);
  w->cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
  w->multifit_p = gsl_multifit_linear_alloc(params->ncurr, ncoeffs);

  w->B = gsl_matrix_alloc(1, 2);
  w->d = gsl_vector_alloc(1);
  w->r = gsl_vector_alloc(params->ncurr);

  w->ntheta = params->ntheta;
  w->theta_min = params->theta_min;
  w->theta_max = params->theta_max;

  w->dtheta = (w->theta_max - w->theta_min) / w->ntheta;

  w->ncurr = params->ncurr;
  w->qdlat_max = params->qdlat_max;
  w->qdlat_step = 2.0 * params->qdlat_max / (params->ncurr - 1.0);

  w->acc = gsl_interp_accel_alloc();
  w->spline_E = gsl_spline_alloc(T, params->ntheta);
  w->spline_u = gsl_spline_alloc(T, params->ntheta);

  return w;
} /* inverteef_alloc() */

void
inverteef_free(inverteef_workspace *w)
{
  if (w->J_pde)
    gsl_vector_free(w->J_pde);

  if (w->J_diff)
    gsl_vector_free(w->J_diff);

  if (w->J_rhs)
    gsl_vector_free(w->J_rhs);

  if (w->J_pde_E)
    gsl_vector_free(w->J_pde_E);

  if (w->J_pde_u)
    gsl_vector_free(w->J_pde_u);

  if (w->X)
    gsl_matrix_free(w->X);

  if (w->cov)
    gsl_matrix_free(w->cov);

  if (w->coeffs)
    gsl_vector_free(w->coeffs);

  if (w->B)
    gsl_matrix_free(w->B);

  if (w->d)
    gsl_vector_free(w->d);

  if (w->r)
    gsl_vector_free(w->r);

  if (w->multifit_p)
    gsl_multifit_linear_free(w->multifit_p);

  if (w->acc)
    gsl_interp_accel_free(w->acc);

  if (w->spline_E)
    gsl_spline_free(w->spline_E);

  if (w->spline_u)
    gsl_spline_free(w->spline_u);

  free(w);
} /* inverteef_free() */

/*
inverteef_thidx()
  Determine index theta in [0, NTHETA - 1] corresponding to given theta

Inputs: theta - colatitude in radians

*/
static inline size_t
inverteef_thidx(double theta, inverteef_workspace *w)
{
  size_t idx;

  idx = (size_t) ((theta - w->theta_min) / w->dtheta);

  return idx;
}

/*
inverteef_interp_pdesol()
  Interpolate PDE solution to the grid used by the line current model

Inputs: J_sat     - satellite line current profile, length ncurr
        theta_pde - theta grid points for PDE solution (radians), length ntheta
        J_lat_E   - PDE solution for E, length ntheta
        J_lat_u   - PDE solution for u, length ntheta

Notes:
1) On output,
w->J_pde_{E,u} contain J_lat_{E,u} interpolated to the line current grid

2) On output,
w->J_rhs contains J_sat - J_pde_u interpolated to line current grid
*/

static int
inverteef_interp_pdesol(const gsl_vector *J_sat, const double *theta_pde,
                        const gsl_vector *J_lat_E, const gsl_vector *J_lat_u,
                        inverteef_workspace *w)
{
  const double qdlat_min = -w->qdlat_max;
  size_t i;

  gsl_spline_init(w->spline_E, theta_pde, J_lat_E->data, w->ntheta);
  gsl_spline_init(w->spline_u, theta_pde, J_lat_u->data, w->ntheta);

  for (i = 0; i < w->ncurr; ++i)
    {
      double qdlat = qdlat_min + i * w->qdlat_step; /* QD latitude on line current grid */
      double theta = M_PI / 2.0 - qdlat * M_PI / 180.0;
      double JE = gsl_spline_eval(w->spline_E, theta, w->acc);
      double Ju = gsl_spline_eval(w->spline_u, theta, w->acc);
      double Js = gsl_vector_get(J_sat, i);

      gsl_vector_set(w->J_pde_E, i, JE);
      gsl_vector_set(w->J_pde_u, i, Ju);

      /* RHS = J_sat_i - J_pde_u_i */
      gsl_vector_set(w->J_rhs, i, Js - Ju);
    }

  return GSL_SUCCESS;
}

/*
inverteef_design_matrix()
  Construct least squares design matrix for the inversion of
the CHAMP profile

The first row of the matrix corresponds to the electric field
scaling factor. The corresponding coefficient is equal to the
factor by which we must scale the electric field.

The second row corresponds to the current DC shift. The
corresponding coefficient is equal to J_DC, where we assume
J_champ = E_scale * J_pde_E + J_pde_u - J_DC

We are minimizing:

min (J_sat - alpha * J_pde_E - J_pde_u + J_DC)

so the least squares equation is:

[ J_E_1 -1.0 ] [ E_scale ] = [ J_champ_1 - J_u_1 ]
[ J_E_2 -1.0 ] [   J_0   ]   [ J_champ_2 - J_u_2 ]
[  ...   ... ]               [        ...        ]
[ J_E_n -1.0 ]               [ J_champ_n - J_u_n ]
*/

static int
inverteef_design_matrix(inverteef_workspace *w)
{
  size_t c; /* column */
  size_t n = w->X->size1;
  size_t i;

  c = 0;

#ifdef PDE_INVERT_E_FIELD

  for (i = 0; i < n; ++i)
    {
      double J = gsl_vector_get(w->J_pde_E, i);

      gsl_matrix_set(w->X, i, c, J);
    }

  ++c;

#endif /* PDE_INVERT_E_FIELD */

#ifdef PDE_ALLOW_DC_SHIFT

  for (i = 0; i < n; ++i)
    {
      gsl_matrix_set(w->X, i, c, -1.0);
    }

  ++c;

#endif /* PDE_ALLOW_DC_SHIFT */

  return GSL_SUCCESS;
} /* inverteef_design_matrix() */

/*
inverteef_invert_profile()
  After the least squares matrix is formed, solve the least
squares problem to find the unknown coefficients

w->coeffs[0] = E field scale factor
w->coeffs[1] = Current DC shift
*/

static int
inverteef_invert_profile(gsl_vector *J_sat, inverteef_workspace *w)
{
  int s;

#ifndef PDE_USE_LSE

  s = gsl_multifit_linear(w->X,
                          w->J_rhs,
                          w->coeffs,
                          w->cov,
                          &(w->chisq),
                          w->multifit_p);

#else

  /* Constrain solution so that J_sol(90 deg) = J_champ(90 deg) */

  {
    double J_E = gsl_vector_get(w->J_pde_E, w->ncurr / 2);
    double J_u = gsl_vector_get(w->J_pde_u, w->ncurr / 2);
    double J_champ = gsl_vector_get(J_sat, w->ncurr / 2);

    gsl_matrix_set(w->B, 0, 0, J_E);
    gsl_matrix_set(w->B, 0, 1, -1.0);
    gsl_vector_set(w->d, 0, J_champ - J_u);

    s = lse(w->X,
            w->B,
            w->J_rhs,
            w->d,
            w->coeffs);

    gsl_multifit_linear_residuals(w->X, w->J_rhs, w->coeffs, w->r);
    gsl_blas_ddot(w->r, w->r, &(w->chisq));

    /* store J_champ to put in datafile */
    w->J_champ = J_champ;
  }

#endif

  if (s)
    return s;

  w->Rsq = 1.0 - w->chisq /
                 gsl_stats_tss(w->J_rhs->data, 1, w->J_rhs->size);

  w->E_scale = gsl_vector_get(w->coeffs, 0);
  w->J_DC = gsl_vector_get(w->coeffs, 1);

  return GSL_SUCCESS;
} /* inverteef_invert_profile() */

/*
inverteef_calc()
  Invert satellite profile

Inputs: J_sat     - satellite current vector profile (length ncurr)
        qdlat_pde - QD latitude grid points for PDE solution
        J_lat_E   - PDE solution (E_0, u = 0) (length ntheta = theta grid points)
        J_lat_u   - PDE solution (E_0 = 0, u) (length ntheta = theta grid points)
        w         - workspace

Notes:
1) On output, w->J_pde contains the modeled profile
2) On output, w->RelErr contains relative error between modeled and
satellite profile
*/

int
inverteef_calc(gsl_vector *J_sat, const double *qdlat_pde,
               gsl_vector *J_lat_E, gsl_vector *J_lat_u, inverteef_workspace *w)
{
  fprintf(stderr, "inverteef_calc: interpolating PDE solution...");

  inverteef_interp_pdesol(J_sat, qdlat_pde, J_lat_E, J_lat_u, w);

  fprintf(stderr, "done\n");
  fflush(stderr);
  fprintf(stderr, "inverteef_calc: constructing least squares matrix...");

  inverteef_design_matrix(w);

  fprintf(stderr, "done\n");
  fflush(stderr);
  fprintf(stderr, "inverteef_calc: inverting profile...");

  inverteef_invert_profile(J_sat, w);

  /* compute final modeled profile and store in w->J_pde */
  {
    size_t j;
    double x[2];
    gsl_vector_view xv = gsl_vector_view_array(x, 2);

    for (j = 0; j < w->ncurr; ++j)
      {
        double J, Jerr;

        x[0] = gsl_vector_get(w->J_pde_E, j);
        x[1] = -1.0;

        gsl_multifit_linear_est(&xv.vector, w->coeffs, w->cov,
                                &J, &Jerr);

        J += gsl_vector_get(w->J_pde_u, j);

        gsl_vector_set(w->J_pde, j, J);
        gsl_vector_set(w->J_diff, j, J - gsl_vector_get(J_sat, j));
      }
  }

  w->R = gsl_stats_correlation(J_sat->data, 1,
                               w->J_pde->data, 1,
                               w->ncurr);

  /* compute relative error between modeled and satellite profiles */
  w->RelErr = gsl_blas_dnrm2(w->J_diff) / gsl_blas_dnrm2(J_sat);

  fprintf(stderr, "done (chisq = %e, Rsq = %f, R = %f)\n",
          w->chisq, w->Rsq, w->R);
  fprintf(stderr, "inverteef_calc: scale factor = %f\n", w->E_scale);
  fflush(stderr);

  return GSL_SUCCESS;
} /* inverteef_calc() */
