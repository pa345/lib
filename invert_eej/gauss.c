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

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multilarge.h>
#include <gsl/gsl_blas.h>

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_coord.h>

#include "green.h"

#include "Gdef.h"

#define GAUSS_FLG_FIT_X             (1 << 0)
#define GAUSS_FLG_FIT_Y             (1 << 1)
#define GAUSS_FLG_FIT_Z             (1 << 2)

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

  size_t flags;     /* GAUSS_FLG_xxx */

  gsl_matrix *X;    /* LS matrix */
  gsl_vector *rhs;  /* rhs vector */
  gsl_vector *wts;  /* weight vector */
  gsl_matrix *cov;  /* covariance matrix */
  gsl_vector *c;    /* solution vector c = [ c_int ; c_ext ] */

  gsl_vector *workX;
  gsl_vector *workY;
  gsl_vector *workZ;

  gsl_multifit_linear_workspace *multifit_p;
  gsl_multilarge_linear_workspace *multilarge_p;
  green_workspace *green_int_p;
  green_workspace *green_ext_p;
} gauss_state_t;

typedef struct
{
  size_t nmax_int;  /* maximum spherical harmonic degree for internal field */
  size_t mmax_int;  /* maximum spherical harmonic order for internal field */
  size_t nmax_ext;  /* maximum spherical harmonic degree for external field */
  size_t mmax_ext;  /* maximum spherical harmonic order for external field */
  size_t flags;     /* fitting flags */
} gauss_parameters;

static void *gauss_alloc(const gauss_parameters * params);
static void gauss_free(void * vstate);
static int gauss_reset(void * vstate);
static size_t gauss_ncoeff(void * vstate);
static int gauss_add_datum(const double t, const double r, const double theta, const double phi,
                           const double qdlat, const double B[3], void * vstate);
static int gauss_fit(double * rnorm, double * snorm, void * vstate);
static int gauss_eval_B(const double t, const double r, const double theta, const double phi,
                        double B[3], void * vstate);
static int gauss_eval_B_int(const double t, const double r, const double theta, const double phi,
                            double B[3], void * vstate);
static int gauss_eval_B_ext(const double t, const double r, const double theta, const double phi,
                            double B[3], void * vstate);

static int build_matrix_row(const double t, const double r, const double theta, const double phi,
                            gsl_vector *X, gsl_vector *Y, gsl_vector *Z,
                            gauss_state_t *state);

/*
gauss_alloc()
  Allocate gauss workspace

Inputs: params - parameters

Return: pointer to workspace
*/

static void *
gauss_alloc(const gauss_parameters * params)
{
  gauss_state_t *state;

  state = calloc(1, sizeof(gauss_state_t));
  if (!state)
    return 0;

  state->lmax_int = params->nmax_int;
  state->mmax_int = params->mmax_int;
  state->lmax_ext = params->nmax_ext;
  state->mmax_ext = params->mmax_ext;
  state->flags = params->flags;
  state->nmax = 10000;
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
  state->c = gsl_vector_alloc(state->p);
  state->rhs = gsl_vector_alloc(state->nmax);
  state->wts = gsl_vector_alloc(state->nmax);
  state->cov = gsl_matrix_alloc(state->p, state->p);
  state->multifit_p = gsl_multifit_linear_alloc(state->nmax, state->p);

  state->workX = gsl_vector_alloc(state->p);
  state->workY = gsl_vector_alloc(state->p);
  state->workZ = gsl_vector_alloc(state->p);

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

  if (state->workX)
    gsl_vector_free(state->workX);

  if (state->workY)
    gsl_vector_free(state->workY);

  if (state->workZ)
    gsl_vector_free(state->workZ);

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
  gsl_vector_view vx, vy, vz;

  (void) t;
  (void) qdlat;

  if (state->flags & GAUSS_FLG_FIT_X)
    {
      vx = gsl_matrix_row(state->X, rowidx);
      gsl_vector_set(state->rhs, rowidx, B[0]);
      gsl_vector_set(state->wts, rowidx, wi);
      ++rowidx;
    }

  if (state->flags & GAUSS_FLG_FIT_Y)
    {
      vy = gsl_matrix_row(state->X, rowidx);
      gsl_vector_set(state->rhs, rowidx, B[1]);
      gsl_vector_set(state->wts, rowidx, wi);
      ++rowidx;
    }

  if (state->flags & GAUSS_FLG_FIT_Z)
    {
      vz = gsl_matrix_row(state->X, rowidx);
      gsl_vector_set(state->rhs, rowidx, B[2]);
      gsl_vector_set(state->wts, rowidx, wi);
      ++rowidx;
    }

  /* build rows of the LS matrix */
  build_matrix_row(t, r, theta, phi, &vx.vector, &vy.vector, &vz.vector, state);

  state->n = rowidx;

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
  gsl_matrix_view A = gsl_matrix_submatrix(state->X, 0, 0, state->n, state->p);
  gsl_vector_view b = gsl_vector_subvector(state->rhs, 0, state->n);
  gsl_vector_view wts = gsl_vector_subvector(state->wts, 0, state->n);
  double rcond, chisq;

  if (state->n < state->p)
    return -1;

  fprintf(stderr, "\n");
  fprintf(stderr, "\t n = %zu\n", state->n);
  fprintf(stderr, "\t p = %zu\n", state->p);

  /* solve system */
#if 0
  gsl_multifit_wlinear(&A.matrix, &wts.vector, &b.vector, state->c, state->cov, &chisq, state->multifit_p);
#else
  {
    const double tol = 1.0e-2;
    size_t rank;

    gsl_multifit_wlinear_tsvd(&A.matrix, &wts.vector, &b.vector, tol, state->c, state->cov,
                              &chisq, &rank, state->multifit_p);
  }
#endif

  rcond = gsl_multifit_linear_rcond(state->multifit_p);
  *rnorm = sqrt(chisq);
  *snorm = gsl_blas_dnrm2(state->c);

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

  build_matrix_row(t, r, theta, phi, &vx.vector, &vy.vector, &vz.vector, state);

  gsl_blas_ddot(&vx.vector, state->c, &B[0]);
  gsl_blas_ddot(&vy.vector, state->c, &B[1]);
  gsl_blas_ddot(&vz.vector, state->c, &B[2]);

  return status;
}

/*
gauss_eval_B_int()
  Evaluate magnetic field at a given (r,theta,phi) using
previously computed coefficients

Inputs: t      - timestamp (CDF_EPOCH)
        r      - radius (km)
        theta  - colatitude (radians)
        phi    - longitude (radians)
        B      - (output) internal magnetic field vector (nT)
        vstate - state

Notes:
1) state->c must contain fit coefficients
*/

static int
gauss_eval_B_int(const double t, const double r, const double theta, const double phi,
                 double B[3], void * vstate)
{
  int status = GSL_SUCCESS;
  gauss_state_t *state = (gauss_state_t *) vstate;
  gsl_vector_view vx = gsl_matrix_row(state->X, 0);
  gsl_vector_view vy = gsl_matrix_row(state->X, 1);
  gsl_vector_view vz = gsl_matrix_row(state->X, 2);
  gsl_vector_view c = gsl_vector_subvector(state->c, 0, state->p_int);

  build_matrix_row(t, r, theta, phi, &vx.vector, &vy.vector, &vz.vector, state);

  vx = gsl_matrix_subrow(state->X, 0, 0, state->p_int);
  vy = gsl_matrix_subrow(state->X, 1, 0, state->p_int);
  vz = gsl_matrix_subrow(state->X, 2, 0, state->p_int);

  gsl_blas_ddot(&vx.vector, &c.vector, &B[0]);
  gsl_blas_ddot(&vy.vector, &c.vector, &B[1]);
  gsl_blas_ddot(&vz.vector, &c.vector, &B[2]);

  return status;
}

/*
gauss_eval_B_ext()
  Evaluate magnetic field at a given (r,theta,phi) using
previously computed coefficients

Inputs: t      - timestamp (CDF_EPOCH)
        r      - radius (km)
        theta  - colatitude (radians)
        phi    - longitude (radians)
        B      - (output) external magnetic field vector (nT)
        vstate - state

Notes:
1) state->c must contain fit coefficients
*/

static int
gauss_eval_B_ext(const double t, const double r, const double theta, const double phi,
                 double B[3], void * vstate)
{
  int status = GSL_SUCCESS;
  gauss_state_t *state = (gauss_state_t *) vstate;
  gsl_vector_view vx = gsl_matrix_row(state->X, 0);
  gsl_vector_view vy = gsl_matrix_row(state->X, 1);
  gsl_vector_view vz = gsl_matrix_row(state->X, 2);
  gsl_vector_view c;

  if (state->p_ext == 0)
    {
      B[0] = B[1] = B[2] = 0.0;
      return status;
    }
  
  c = gsl_vector_subvector(state->c, state->ext_offset, state->p_ext);

  build_matrix_row(t, r, theta, phi, &vx.vector, &vy.vector, &vz.vector, state);

  vx = gsl_matrix_subrow(state->X, 0, state->ext_offset, state->p_ext);
  vy = gsl_matrix_subrow(state->X, 1, state->ext_offset, state->p_ext);
  vz = gsl_matrix_subrow(state->X, 2, state->ext_offset, state->p_ext);

  gsl_blas_ddot(&vx.vector, &c.vector, &B[0]);
  gsl_blas_ddot(&vy.vector, &c.vector, &B[1]);
  gsl_blas_ddot(&vz.vector, &c.vector, &B[2]);

  return status;
}

static int
build_matrix_row(const double t, const double r, const double theta, const double phi,
                 gsl_vector *X, gsl_vector *Y, gsl_vector *Z, gauss_state_t *state)
{
  int s = 0;
 
  if (state->p_int > 0)
    {
      gsl_vector_view wx = gsl_vector_subvector(state->workX, 0, state->p_int);
      gsl_vector_view wy = gsl_vector_subvector(state->workY, 0, state->p_int);
      gsl_vector_view wz = gsl_vector_subvector(state->workZ, 0, state->p_int);

      /* store green's functions in work arrays */
      s = green_calc_int(r, theta, phi, wx.vector.data, wy.vector.data, wz.vector.data, state->green_int_p);
      if (s)
        return s;

      if (state->flags & GAUSS_FLG_FIT_X)
        {
          gsl_vector_view v = gsl_vector_subvector(X, 0, state->p_int);
          gsl_vector_memcpy(&v.vector, &wx.vector);
        }

      if (state->flags & GAUSS_FLG_FIT_Y)
        {
          gsl_vector_view v = gsl_vector_subvector(Y, 0, state->p_int);
          gsl_vector_memcpy(&v.vector, &wy.vector);
        }

      if (state->flags & GAUSS_FLG_FIT_Z)
        {
          gsl_vector_view v = gsl_vector_subvector(Z, 0, state->p_int);
          gsl_vector_memcpy(&v.vector, &wz.vector);
        }
    }

  if (state->p_ext > 0)
    {
      double fday = satdata_epoch2fday(t);
      double lat = M_PI / 2.0 - theta;
      gsl_vector_view wx = gsl_vector_subvector(state->workX, state->ext_offset, state->p_ext);
      gsl_vector_view wy = gsl_vector_subvector(state->workY, state->ext_offset, state->p_ext);
      gsl_vector_view wz = gsl_vector_subvector(state->workZ, state->ext_offset, state->p_ext);
      double phi_sm, lat_sm, theta_sm;
      size_t i;

      /* compute basis functions for K(r,theta_SM,phi_SM) */
      trans(GEO2SM, fday, phi, lat, &phi_sm, &lat_sm);
      theta_sm = M_PI / 2.0 - lat_sm;

      /* store SM green's functions in work arrays */
      s = green_calc_ext(r, theta_sm, phi_sm, wx.vector.data, wy.vector.data, wz.vector.data, state->green_ext_p);
      if (s)
        return s;

      /* transform back to GEO */
      for (i = 0; i < state->green_ext_p->nnm; ++i)
        {
          double *Xi = gsl_vector_ptr(&wx.vector, i);
          double *Yi = gsl_vector_ptr(&wy.vector, i);
          double *Zi = gsl_vector_ptr(&wz.vector, i);
          double phi_ss, lat_ss;
          double Xt, Yt, Zt;

          trans_vec(SM2GEO, fday, phi_sm, lat_sm, *Xi, *Yi, *Zi, &phi_ss,
                    &lat_ss, &Xt, &Yt, &Zt);

          *Xi = Xt;
          *Yi = Yt;
          *Zi = Zt;
        }

      if (state->flags & GAUSS_FLG_FIT_X)
        {
          gsl_vector_view v = gsl_vector_subvector(X, state->ext_offset, state->p_ext);
          gsl_vector_memcpy(&v.vector, &wx.vector);
        }

      if (state->flags & GAUSS_FLG_FIT_Y)
        {
          gsl_vector_view v = gsl_vector_subvector(Y, state->ext_offset, state->p_ext);
          gsl_vector_memcpy(&v.vector, &wy.vector);
        }

      if (state->flags & GAUSS_FLG_FIT_Z)
        {
          gsl_vector_view v = gsl_vector_subvector(Z, state->ext_offset, state->p_ext);
          gsl_vector_memcpy(&v.vector, &wz.vector);
        }
    }

  return s;
}
