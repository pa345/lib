/*
 * sheet.c
 *
 * Fit magnetic field measurements due to a spherical sheet current at a
 * specified altitude.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <satdata/satdata.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>

#include <common/common.h>
#include <common/interp.h>
#include <common/oct.h>

#include "green.h"
#include "lapack_wrapper.h"
#include "track.h"

#include "magfit.h"

#define NULL_VECTOR_VIEW {{0, 0, 0, 0, 0}}

typedef struct
{
  size_t n;         /* total number of measurements in system */
  size_t p;         /* number of coefficients */
  size_t nmax;      /* maximum number of measurements in LS system */
  double b;         /* shell radius (km) */

  size_t lmax_int;  /* maximum spherical harmonic degree for internal field */
  size_t mmax_int;  /* maximum spherical harmonic order for internal field */

  size_t p_int;     /* number of internal coefficients */

  gsl_matrix *X;    /* LS matrix, nmax-by-p */
  gsl_vector *XTy;  /* X^T y */
  gsl_matrix *XTX;  /* X^T X, p-by-p */
  gsl_vector *rhs;  /* rhs vector */
  gsl_vector *wts;  /* weight vector */
  gsl_matrix *cov;  /* covariance matrix */
  gsl_vector *c;    /* solution vector gnm */

  gsl_vector *k;    /* klm coefficients for magnetic field below current shell */

  size_t flags;     /* fitting flags */

  green_workspace *green_p;
} sheet_state_t;

static void *sheet_alloc(const void * params);
static void sheet_free(void * vstate);
static int sheet_reset(void * vstate);
static size_t sheet_ncoeff(void * vstate);
static int sheet_add_datum(const double t, const double r, const double theta, const double phi,
                           const double qdlat, const double B[3], void * vstate);
static int sheet_fit(double * rnorm, double * snorm, void * vstate);
static int sheet_eval_B(const double t, const double r, const double theta, const double phi,
                        double B[3], void * vstate);
static int sheet_eval_J(const double r, const double theta, const double phi,
                        double J[3], void * vstate);
static double sheet_eval_chi(const double theta, const double phi, void * vstate);

static int build_matrix_row(const double r, const double theta, const double phi,
                            gsl_vector *X, gsl_vector *Y, gsl_vector *Z,
                            sheet_state_t *state);

/*
sheet_alloc()
  Allocate sheet workspace

Inputs: params - parameters

Return: pointer to workspace
*/

static void *
sheet_alloc(const void * params)
{
  const magfit_parameters *mparams = (const magfit_parameters *) params;
  const size_t max_observations = 50000;
  sheet_state_t *state;

  state = calloc(1, sizeof(sheet_state_t));
  if (!state)
    return 0;

  state->b = mparams->b;
  state->lmax_int = mparams->nmax_int;
  state->mmax_int = mparams->mmax_int;
  state->nmax = 3 * max_observations;
  state->flags = mparams->flags;
  state->n = 0;

  if (state->lmax_int > 0)
    {
      state->green_p = green_alloc(state->lmax_int, state->mmax_int, R_EARTH_KM);
      state->p_int = green_nnm(state->green_p);
    }
  else
    state->p_int = 0;

  state->p = state->p_int;

  state->X = gsl_matrix_alloc(state->nmax, state->p);
  state->XTX = gsl_matrix_calloc(state->p, state->p);
  state->XTy = gsl_vector_calloc(state->p);
  state->c = gsl_vector_alloc(state->p);
  state->rhs = gsl_vector_alloc(state->nmax);
  state->wts = gsl_vector_alloc(state->nmax);
  state->cov = gsl_matrix_alloc(state->p, state->p);
  state->k = gsl_vector_alloc(state->p);

  return state;
}

static void
sheet_free(void * vstate)
{
  sheet_state_t *state = (sheet_state_t *) vstate;

  if (state->X)
    gsl_matrix_free(state->X);

  if (state->c)
    gsl_vector_free(state->c);

  if (state->rhs)
    gsl_vector_free(state->rhs);

  if (state->wts)
    gsl_vector_free(state->wts);

  if (state->cov)
    gsl_matrix_free(state->cov);

  if (state->k)
    gsl_vector_free(state->k);

  if (state->green_p)
    green_free(state->green_p);

  free(state);
}

/* reset workspace to work on new data set */
static int
sheet_reset(void * vstate)
{
  sheet_state_t *state = (sheet_state_t *) vstate;
  gsl_matrix_set_zero(state->XTX);
  gsl_vector_set_zero(state->XTy);
  state->n = 0;
  return 0;
}

static size_t
sheet_ncoeff(void * vstate)
{
  sheet_state_t *state = (sheet_state_t *) vstate;
  return state->p;
}

/*
sheet_add_datum()
  Add single vector measurement to LS system

Inputs: t      - timestamp (CDF_EPOCH)
        r      - radius (km)
        theta  - colatitude (radians)
        phi    - longitude (radians)
        qdlat  - QD latitude (degrees)
        B      - magnetic field vector NEC (nT)
        vstate - state

Return: success/error

Notes:
1) state->n is updated with the number of total data added
*/

static int
sheet_add_datum(const double t, const double r, const double theta, const double phi,
                   const double qdlat, const double B[3], void * vstate)
{
  sheet_state_t *state = (sheet_state_t *) vstate;
  size_t rowidx = state->n;
  double wi = 1.0;
  gsl_vector_view vx = NULL_VECTOR_VIEW;
  gsl_vector_view vy = NULL_VECTOR_VIEW;
  gsl_vector_view vz = NULL_VECTOR_VIEW;

  (void) t;
  (void) qdlat;

  if (state->flags & MAGFIT_FLG_FIT_X)
    {
      vx = gsl_matrix_row(state->X, rowidx);
      gsl_vector_set(state->rhs, rowidx, B[0]);
      gsl_vector_set(state->wts, rowidx, wi);
      ++rowidx;
    }

  if (state->flags & MAGFIT_FLG_FIT_Y)
    {
      vy = gsl_matrix_row(state->X, rowidx);
      gsl_vector_set(state->rhs, rowidx, B[1]);
      gsl_vector_set(state->wts, rowidx, wi);
      ++rowidx;
    }

  if (state->flags & MAGFIT_FLG_FIT_Z)
    {
      vz = gsl_matrix_row(state->X, rowidx);
      gsl_vector_set(state->rhs, rowidx, B[2]);
      gsl_vector_set(state->wts, rowidx, wi);
      ++rowidx;
    }

  /* build desired rows of the LS matrix */
  build_matrix_row(r, theta, phi, &vx.vector, &vy.vector, &vz.vector, state);

  state->n = rowidx;

  if (state->n >= state->nmax)
    {
      /* fold observations into XTX and XTy */
      gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, state->X, 1.0, state->XTX);
      gsl_blas_dgemv(CblasTrans, 1.0, state->X, state->rhs, 1.0, state->XTy);

      state->n = 0;
    }

  return GSL_SUCCESS;
}

/*
sheet_fit()
  Fit model to previously added tracks

Inputs: rnorm  - residual norm || y - A x ||
        snorm  - solution norm || x ||
        vstate - state

Return: success/error

Notes:
1) Data must be added to workspace via sheet_add_datum()
*/

static int
sheet_fit(double * rnorm, double * snorm, void * vstate)
{
  sheet_state_t *state = (sheet_state_t *) vstate;
  double rcond;

  if (state->n < state->p)
    return -1;

  if (state->n > 0)
    {
      gsl_matrix_view Xv = gsl_matrix_submatrix(state->X, 0, 0, state->n, state->p);
      gsl_vector_view yv = gsl_vector_subvector(state->rhs, 0, state->n);

      /* fold observations into XTX and XTy */
      gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &Xv.matrix, 1.0, state->XTX);
      gsl_blas_dgemv(CblasTrans, 1.0, &Xv.matrix, &yv.vector, 1.0, state->XTy);

      state->n = 0;
    }

  fprintf(stderr, "\n");
  fprintf(stderr, "\t p = %zu\n", state->p);

  lapack_cholesky_solve(state->XTX, state->XTy, state->c, &rcond, NULL);
  *rnorm = 0.0;
  *snorm = gsl_blas_dnrm2(state->c);

  fprintf(stderr, "\t rnorm = %g\n", *rnorm);
  fprintf(stderr, "\t snorm = %g\n", *snorm);
  fprintf(stderr, "\t cond(X) = %g\n", 1.0 / rcond);

  /* store knm coefficients */
  green_g2k(state->b, state->c, state->k, state->green_p);

  return 0;
}

/*
sheet_eval_B()
  Evaluate magnetic field at a given (r,theta,phi) using
previously computed coefficients

Inputs: t      - timestamp (CDF_EPOCH)
        r      - radius (km)
        theta  - colatitude (radians)
        phi    - longitude (radians)
        B      - (output) magnetic field vector (nT)
        vstate - state

Notes:
1) state->c must contain fit coefficients
*/

static int
sheet_eval_B(const double t, const double r, const double theta, const double phi,
             double B[3], void * vstate)
{
  int status = GSL_SUCCESS;
  sheet_state_t *state = (sheet_state_t *) vstate;
  gsl_vector_view vx = gsl_matrix_row(state->X, 0);
  gsl_vector_view vy = gsl_matrix_row(state->X, 1);
  gsl_vector_view vz = gsl_matrix_row(state->X, 2);

  (void) t;

  if (r >= state->b)
    {
      green_calc_int(r, theta, phi, vx.vector.data, vy.vector.data, vz.vector.data, state->green_p);

      gsl_blas_ddot(&vx.vector, state->c, &B[0]);
      gsl_blas_ddot(&vy.vector, state->c, &B[1]);
      gsl_blas_ddot(&vz.vector, state->c, &B[2]);
    }
  else
    {
      green_calc_ext(r, theta, phi, vx.vector.data, vy.vector.data, vz.vector.data, state->green_p);

      gsl_blas_ddot(&vx.vector, state->k, &B[0]);
      gsl_blas_ddot(&vy.vector, state->k, &B[1]);
      gsl_blas_ddot(&vz.vector, state->k, &B[2]);
    }

  return status;
}

/*
sheet_eval_J()
  Evaluate current density at a given (r,theta,phi) using
previously computed coefficients

Inputs: r      - radius (km)
        theta  - colatitude (radians)
        phi    - longitude (radians)
        J      - (output) current density vector [A/km]
        vstate - workspace

Notes:
1) state->c must contain coefficients
*/

static int
sheet_eval_J(const double r, const double theta, const double phi,
                double J[3], void * vstate)
{
  int status;
  sheet_state_t *state = (sheet_state_t *) vstate;

  (void) r; /* unused parameter */

  status = green_eval_sheet_int(R_EARTH_KM + 110.0, theta, phi, state->c, J,
                                state->green_p);

  return status;
}

/*
sheet_eval_chi()
  Evaluate current stream function at a given (theta,phi) using
previously computed coefficients

Inputs: theta  - colatitude (radians)
        phi    - longitude (radians)
        vstate - workspace

Return: current stream function chi in kA/nT

Notes:
1) state->c must contain coefficients
*/

static double
sheet_eval_chi(const double theta, const double phi, void * vstate)
{
  sheet_state_t *state = (sheet_state_t *) vstate;
  double chi;

  chi = green_eval_chi_int(R_EARTH_KM + 110.0, theta, phi, state->c, state->green_p);

  return chi;
}

static int
build_matrix_row(const double r, const double theta, const double phi,
                 gsl_vector *X, gsl_vector *Y, gsl_vector *Z,
                 sheet_state_t *state)
{
  int s = 0;
 
  if (state->p_int > 0)
    {
      gsl_vector_view xv = NULL_VECTOR_VIEW;
      gsl_vector_view yv = NULL_VECTOR_VIEW;
      gsl_vector_view zv = NULL_VECTOR_VIEW;

      if (X->data)
        xv = gsl_vector_subvector(X, 0, state->p_int);

      if (Y->data)
        yv = gsl_vector_subvector(Y, 0, state->p_int);

      if (Z->data)
        zv = gsl_vector_subvector(Z, 0, state->p_int);

      s = green_calc_int(r, theta, phi, xv.vector.data, yv.vector.data, zv.vector.data, state->green_p);
      if (s)
        return s;
    }

  return s;
}

static const magfit_type sheet_type =
{
  "sheet",
  sheet_alloc,
  sheet_reset,
  sheet_ncoeff,
  sheet_add_datum,
  sheet_fit,
  sheet_eval_B,
  sheet_eval_J,
  sheet_eval_chi,
  sheet_free
};

const magfit_type *magfit_sheet = &sheet_type;
