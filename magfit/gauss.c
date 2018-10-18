/*
 * gauss.c
 *
 * Fit internal and external potential field based on Gauss coefficients
 * to magnetic vector data
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
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multilarge.h>
#include <gsl/gsl_blas.h>

#include <common/common.h>
#include <common/interp.h>
#include <common/oct.h>
#include <track/track.h>

#include "green.h"
#include "lapack_wrapper.h"

#include "magfit.h"

#define NULL_VECTOR_VIEW {{0, 0, 0, 0, 0}}

typedef struct
{
  size_t n;         /* total number of measurements in system */
  size_t p;         /* number of coefficients */
  size_t nmax;      /* maximum number of measurements in LS system */

  size_t lmax_int;  /* maximum spherical harmonic degree for internal field */
  size_t mmax_int;  /* maximum spherical harmonic order for internal field */
  size_t lmax_ext;  /* maximum spherical harmonic degree for external field */
  size_t mmax_ext;  /* maximum spherical harmonic order for external field */

  size_t p_int;     /* number of internal coefficients */
  size_t p_ext;     /* number of external coefficients */
  size_t ext_offset; /* offset of external coefficients in solution vector 'c' */

  gsl_matrix *X;    /* LS matrix, nmax-by-p */
  gsl_vector *XTy;  /* X^T y */
  gsl_matrix *XTX;  /* X^T X, p-by-p */
  gsl_vector *rhs;  /* rhs vector */
  gsl_vector *wts;  /* weight vector */
  gsl_matrix *cov;  /* covariance matrix */
  gsl_vector *c;    /* solution vector c = [ c_int ; c_ext ] */

  size_t flags;     /* fitting flags */

  gsl_multifit_linear_workspace *multifit_p;
  gsl_multilarge_linear_workspace *multilarge_p;
  green_workspace *green_int_p;
  green_workspace *green_ext_p;
} gauss_state_t;

static void *gauss_alloc(const void * params);
static void gauss_free(void * vstate);
static int gauss_reset(void * vstate);
static size_t gauss_ncoeff(void * vstate);
static int gauss_add_datum(const double t, const double r, const double theta, const double phi,
                           const double qdlat, const double B[3], void * vstate);
static int gauss_fit(double * rnorm, double * snorm, void * vstate);
static int gauss_eval_B(const double t, const double r, const double theta, const double phi,
                        double B[3], void * vstate);
static int gauss_eval_J(const double r, const double theta, const double phi,
                        double J[3], void * vstate);
static double gauss_eval_chi(const double theta, const double phi, void * vstate);

static int build_matrix_row(const double r, const double theta, const double phi,
                            gsl_vector *X, gsl_vector *Y, gsl_vector *Z,
                            gauss_state_t *state);

/*
gauss_alloc()
  Allocate gauss workspace

Inputs: params - parameters

Return: pointer to workspace
*/

static void *
gauss_alloc(const void * params)
{
  const magfit_parameters *mparams = (const magfit_parameters *) params;
  const size_t max_observations = 50000;
  gauss_state_t *state;

  state = calloc(1, sizeof(gauss_state_t));
  if (!state)
    return 0;

  state->lmax_int = mparams->nmax_int;
  state->mmax_int = mparams->mmax_int;
  state->lmax_ext = mparams->nmax_ext;
  state->mmax_ext = mparams->mmax_ext;
  state->nmax = 3 * max_observations;
  state->flags = mparams->flags;
  state->n = 0;

  if (state->lmax_int > 0)
    {
      state->green_int_p = green_alloc(state->lmax_int, state->mmax_int, R_EARTH_KM);
      state->p_int = green_nnm(state->green_int_p);
    }
  else
    state->p_int = 0;

  if (state->lmax_ext > 0)
    {
      state->green_ext_p = green_alloc(state->lmax_ext, state->mmax_ext, R_EARTH_KM);
      state->p_ext = green_nnm(state->green_ext_p);
    }
  else
    state->p_ext = 0;

  state->ext_offset = state->p_int;

  state->p = state->p_int + state->p_ext;

  state->X = gsl_matrix_alloc(state->nmax, state->p);
  state->XTX = gsl_matrix_calloc(state->p, state->p);
  state->XTy = gsl_vector_calloc(state->p);
  state->c = gsl_vector_alloc(state->p);
  state->rhs = gsl_vector_alloc(state->nmax);
  state->wts = gsl_vector_alloc(state->nmax);
  state->cov = gsl_matrix_alloc(state->p, state->p);
  state->multifit_p = gsl_multifit_linear_alloc(state->nmax, state->p);

  state->multilarge_p = gsl_multilarge_linear_alloc(gsl_multilarge_linear_normal, state->p);

  return state;
}

static void
gauss_free(void * vstate)
{
  gauss_state_t *state = (gauss_state_t *) vstate;

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

  if (state->multifit_p)
    gsl_multifit_linear_free(state->multifit_p);

  if (state->multilarge_p)
    gsl_multilarge_linear_free(state->multilarge_p);

  if (state->green_int_p)
    green_free(state->green_int_p);

  if (state->green_ext_p)
    green_free(state->green_ext_p);

  free(state);
}

/* reset workspace to work on new data set */
static int
gauss_reset(void * vstate)
{
  gauss_state_t *state = (gauss_state_t *) vstate;
  gsl_matrix_set_zero(state->XTX);
  gsl_vector_set_zero(state->XTy);
  state->n = 0;
  return 0;
}

static size_t
gauss_ncoeff(void * vstate)
{
  gauss_state_t *state = (gauss_state_t *) vstate;
  return state->p;
}

/*
gauss_add_datum()
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
gauss_add_datum(const double t, const double r, const double theta, const double phi,
                   const double qdlat, const double B[3], void * vstate)
{
  gauss_state_t *state = (gauss_state_t *) vstate;
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
      gsl_multilarge_linear_accumulate(state->X, state->rhs, state->multilarge_p);

      gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, state->X, 1.0, state->XTX);
      gsl_blas_dgemv(CblasTrans, 1.0, state->X, state->rhs, 1.0, state->XTy);

      state->n = 0;
    }

  return GSL_SUCCESS;
}

/*
gauss_fit()
  Fit model to previously added tracks

Inputs: rnorm  - residual norm || y - A x ||
        snorm  - solution norm || x ||
        vstate - state

Return: success/error

Notes:
1) Data must be added to workspace via gauss_add_datum()
*/

static int
gauss_fit(double * rnorm, double * snorm, void * vstate)
{
  gauss_state_t *state = (gauss_state_t *) vstate;
  double rcond, chisq;

  if (state->n < state->p)
    return -1;

  if (state->n > 0)
    {
      gsl_matrix_view Xv = gsl_matrix_submatrix(state->X, 0, 0, state->n, state->p);
      gsl_vector_view yv = gsl_vector_subvector(state->rhs, 0, state->n);

      /* fold observations into XTX and XTy */
      gsl_multilarge_linear_accumulate(&Xv.matrix, &yv.vector, state->multilarge_p);

      gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, &Xv.matrix, 1.0, state->XTX);
      gsl_blas_dgemv(CblasTrans, 1.0, &Xv.matrix, &yv.vector, 1.0, state->XTy);

      state->n = 0;
    }

  fprintf(stderr, "\n");
  fprintf(stderr, "\t p = %zu\n", state->p);

#if 0
  /* solve system */
  {
    gsl_matrix_view A = gsl_matrix_submatrix(state->X, 0, 0, state->n, state->p);
    gsl_vector_view b = gsl_vector_subvector(state->rhs, 0, state->n);
    gsl_vector_view wts = gsl_vector_subvector(state->wts, 0, state->n);

    gsl_multifit_wlinear(&A.matrix, &wts.vector, &b.vector, state->c, state->cov, &chisq, state->multifit_p);
    rcond = gsl_multifit_linear_rcond(state->multifit_p);
    *rnorm = sqrt(chisq);
    *snorm = gsl_blas_dnrm2(state->c);
  }
#elif 0
  gsl_multilarge_linear_solve(0.0, state->c, rnorm, snorm, state->multilarge_p);
  gsl_multilarge_linear_rcond(&rcond, state->multilarge_p);
#else
  lapack_cholesky_solve(state->XTX, state->XTy, state->c, &rcond, NULL);
  *rnorm = 0.0;
  *snorm = gsl_blas_dnrm2(state->c);
#endif

  fprintf(stderr, "\t rnorm = %g\n", *rnorm);
  fprintf(stderr, "\t snorm = %g\n", *snorm);
  fprintf(stderr, "\t cond(X) = %g\n", 1.0 / rcond);

  return 0;
}

/*
gauss_eval_B()
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
gauss_eval_B(const double t, const double r, const double theta, const double phi,
                double B[3], void * vstate)
{
  int status = GSL_SUCCESS;
  gauss_state_t *state = (gauss_state_t *) vstate;
  gsl_vector_view vx = gsl_matrix_row(state->X, 0);
  gsl_vector_view vy = gsl_matrix_row(state->X, 1);
  gsl_vector_view vz = gsl_matrix_row(state->X, 2);

  (void) t;

  build_matrix_row(r, theta, phi, &vx.vector, &vy.vector, &vz.vector, state);

  gsl_blas_ddot(&vx.vector, state->c, &B[0]);
  gsl_blas_ddot(&vy.vector, state->c, &B[1]);
  gsl_blas_ddot(&vz.vector, state->c, &B[2]);

  return status;
}

/*
gauss_eval_J()
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
gauss_eval_J(const double r, const double theta, const double phi,
                double J[3], void * vstate)
{
  int status;
  gauss_state_t *state = (gauss_state_t *) vstate;

  (void) r; /* unused parameter */

  status = green_eval_sheet_int(R_EARTH_KM + 110.0, theta, phi, state->c, J,
                                state->green_int_p);

  return status;
}

/*
gauss_eval_chi()
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
gauss_eval_chi(const double theta, const double phi, void * vstate)
{
  gauss_state_t *state = (gauss_state_t *) vstate;
  double chi;

  chi = green_eval_chi_int(R_EARTH_KM + 110.0, theta, phi, state->c, state->green_int_p);

  return chi;
}

static int
build_matrix_row(const double r, const double theta, const double phi,
                 gsl_vector *X, gsl_vector *Y, gsl_vector *Z,
                 gauss_state_t *state)
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

      s = green_calc_int(r, theta, phi, xv.vector.data, yv.vector.data, zv.vector.data, state->green_int_p);
      if (s)
        return s;
    }

  if (state->p_ext > 0)
    {
      gsl_vector_view xv = NULL_VECTOR_VIEW;
      gsl_vector_view yv = NULL_VECTOR_VIEW;
      gsl_vector_view zv = NULL_VECTOR_VIEW;

      if (X->data)
        xv = gsl_vector_subvector(X, state->ext_offset, state->p_ext);

      if (Y->data)
        yv = gsl_vector_subvector(Y, state->ext_offset, state->p_ext);

      if (Z->data)
        zv = gsl_vector_subvector(Z, state->ext_offset, state->p_ext);

      s = green_calc_ext(r, theta, phi, xv.vector.data, yv.vector.data, zv.vector.data, state->green_ext_p);
      if (s)
        return s;
    }

  return s;
}

static const magfit_type gauss_type =
{
  "gauss",
  gauss_alloc,
  gauss_reset,
  gauss_ncoeff,
  gauss_add_datum,
  gauss_fit,
  gauss_eval_B,
  gauss_eval_J,
  gauss_eval_chi,
  gauss_free
};

const magfit_type *magfit_gauss = &gauss_type;
