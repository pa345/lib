/*
 * secs2d.c
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
#include <gsl/gsl_blas.h>

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_indices.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_interp.h>
#include <mainlib/ml_oct.h>
#include <mainlib/ml_track.h>

#include "magfit.h"

/* define to use uniform pole spacing in latitude */
#define SECS2D_UNIFORM_THETA         1

/* define to upweight equatorial data */
#define SECS2D_WEIGHT_EEJ            (1.0)

typedef struct
{
  size_t n;         /* total number of measurements in system */
  size_t p;         /* p_df + p_cf */
  size_t nmax;      /* maximum number of measurements in LS system */
  size_t flags;     /* MAGFIT_FLG_xxx */
  size_t ntheta;    /* number of poles in theta direction */
  size_t nphi;      /* number of poles in phi direction */

  size_t df_offset; /* offset in 'c' of DF coefficients */
  size_t cf_offset; /* offset in 'c' of CF coefficients */

  double R;         /* radius of ionosphere (km) */

  double *theta0;   /* theta pole locations, size ntheta */
  double *phi0;     /* phi pole locations, size nphi */

  gsl_matrix *X;    /* LS matrix */
  gsl_vector *rhs;  /* rhs vector */
  gsl_vector *wts;  /* weight vector */
  gsl_matrix *cov;  /* covariance matrix */

  /*
   * solution vector is organized as:
   * c = [ c_df ; c_cf ]
   * where c_df and c_cf are both length npoles = ntheta * nphi
   */
  gsl_vector *c;    /* solution vector */

  gsl_multifit_linear_workspace *multifit_p;
} secs2d_state_t;

static void *secs2d_alloc(const void * params);
static void secs2d_free(void * vstate);
static int secs2d_reset(void * vstate);
static size_t secs2d_ncoeff(void * vstate);
static int secs2d_add_datum(const double t, const double r, const double theta, const double phi,
                            const double qdlat, const double B[4], void * vstate);
static int secs2d_fit(double * rnorm, double * snorm, void * vstate);
static int secs2d_eval_B(const double t, const double r, const double theta, const double phi,
                         double B[3], void * vstate);
static int secs2d_eval_J(const double r, const double theta, const double phi,
                         double J[3], void * vstate);

static int secs2d_green_df(const double r, const double theta, const double phi,
                           const double theta0, const double phi0,
                           double B[3], secs2d_state_t *state);
static int secs2d_green_df_J(const double theta, const double phi, const double theta0, const double phi0,
                             double K[3], secs2d_state_t *state);
static int secs2d_matrix_row_df(const double r, const double theta, const double phi,
                                gsl_vector *X, gsl_vector *Y, gsl_vector *Z, secs2d_state_t *state);
static int secs2d_matrix_row_cf(const double r, const double theta,
                                gsl_vector *Y, secs2d_state_t *state);
static int secs2d_matrix_row_df_J(const double theta, const double phi,
                                  gsl_vector *X, gsl_vector *Y, gsl_vector *Z,
                                  secs2d_state_t *state);
static int secs2d_matrix_row_cf_J(const double theta, gsl_vector *X, gsl_vector *Z, secs2d_state_t *state);
static int secs2d_transform(const double theta0, const double phi0, const double theta, const double phi,
                            double *thetap);
static int secs2d_transform_vec(const double theta0, const double phi0, const double theta, const double phi,
                                const double v_el[3], double v_geo[3]);
static size_t secs2d_idx(const size_t i, const size_t j, const secs2d_state_t *state);

/*
secs2d_alloc()
  Allocate secs 2d workspace

Inputs: params - parameters

Return: pointer to workspace
*/

static void *
secs2d_alloc(const void * params)
{
  const magfit_parameters *mparams = (const magfit_parameters *) params;
  secs2d_state_t *state;
  const size_t flags = mparams->flags;
  const double lat_spacing = mparams->lat_spacing2d;
  const double lat_max = mparams->lat_max;
  const double lat_min = mparams->lat_min;
#if SECS2D_UNIFORM_THETA
  const size_t ntheta = (size_t) ((lat_max - lat_min) / lat_spacing + 1.0);
  const double dtheta = lat_spacing * M_PI / 180.0;
#else
  const size_t ntheta = 31; /* 5 degree spacing outside [-10,10] latitude, 2 degree spacing inside */
#endif
  const double theta_min = M_PI / 2.0 - lat_max * M_PI / 180.0;
  const double lon_spacing = mparams->lon_spacing;
  const double lon_max = mparams->lon_max;
  const double lon_min = mparams->lon_min;
  const double phi_min = lon_min * M_PI / 180.0;
  const double diff = lon_max - lon_min;
  const size_t nphi = GSL_MAX((size_t) (fabs(diff) / lon_spacing + 1.0), 2);
  const double dphi = diff / (nphi - 1.0) * M_PI / 180.0;
  const size_t npoles = ntheta * nphi;
  size_t i;

  state = calloc(1, sizeof(secs2d_state_t));
  if (!state)
    return 0;

  state->nmax = 30000;
  state->n = 0;
  state->p = 0;
  state->ntheta = ntheta;
  state->nphi = nphi;
  state->R = mparams->R;
  state->flags = flags;
  state->df_offset = 0;
  state->cf_offset = 0;

  if (flags & MAGFIT_FLG_SECS_FIT_DF)
    {
      state->df_offset = state->p;
      state->p += npoles;
    }

  if (flags & MAGFIT_FLG_SECS_FIT_CF)
    {
      state->cf_offset = state->p;
      state->p += npoles;
    }

  state->X = gsl_matrix_alloc(state->nmax, state->p);
  state->c = gsl_vector_alloc(state->p);
  state->rhs = gsl_vector_alloc(state->nmax);
  state->wts = gsl_vector_alloc(state->nmax);
  state->cov = gsl_matrix_alloc(state->p, state->p);
  state->multifit_p = gsl_multifit_linear_alloc(state->nmax, state->p);

  state->theta0 = malloc(ntheta * sizeof(double));
  state->phi0 = malloc(nphi * sizeof(double));

  /* fill in theta0 array */
#if SECS2D_UNIFORM_THETA
  for (i = 0; i < ntheta; ++i)
    {
      state->theta0[i] = theta_min + i * dtheta;
    }
#else
  {
    double dtheta = 5.0 * M_PI / 180.0;
    double theta1 = M_PI / 2.0 - 10.0 * M_PI / 180.0;
    double theta2 = M_PI / 2.0 + 10.0 * M_PI / 180.0;

    state->theta0[0] = theta_min;
    for (i = 1; i < ntheta; ++i)
      {
        state->theta0[i] = state->theta0[i - 1] + dtheta;
        if (i >= 10 && i < 20)
          dtheta = 2.0 * M_PI / 180.0;
        else
          dtheta = 5.0 * M_PI / 180.0;
      }
  }
#endif

  /* fill in phi0 array */
  for (i = 0; i < nphi; ++i)
    {
      state->phi0[i] = phi_min + i * dphi;
    }

  fprintf(stderr, "secs2d_alloc: ntheta = %zu\n", state->ntheta);
  fprintf(stderr, "secs2d_alloc: nphi   = %zu\n", state->nphi);
  fprintf(stderr, "secs2d_alloc: ncoeff = %zu\n", state->p);

  return state;
}

static void
secs2d_free(void * vstate)
{
  secs2d_state_t *state = (secs2d_state_t *) vstate;

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

  if (state->theta0)
    free(state->theta0);

  if (state->phi0)
    free(state->phi0);

  if (state->multifit_p)
    gsl_multifit_linear_free(state->multifit_p);

  free(state);
}

/* reset workspace to work on new data set */
static int
secs2d_reset(void * vstate)
{
  secs2d_state_t *state = (secs2d_state_t *) vstate;
  state->n = 0;
  return 0;
}

static size_t
secs2d_ncoeff(void * vstate)
{
  secs2d_state_t *state = (secs2d_state_t *) vstate;
  return state->p;
}

/*
secs2d_add_datum()
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
secs2d_add_datum(const double t, const double r, const double theta, const double phi,
                 const double qdlat, const double B[4], void * vstate)
{
  secs2d_state_t *state = (secs2d_state_t *) vstate;
  size_t rowidx = state->n;
  double wi = 1.0;
  gsl_vector_view vx = gsl_matrix_row(state->X, rowidx);
  gsl_vector_view vy = gsl_matrix_row(state->X, rowidx + 1);
  gsl_vector_view vz = gsl_matrix_row(state->X, rowidx + 2);

  (void) t;

  if (state->flags & MAGFIT_FLG_SECS_FIT_DF)
    {
      /* upweight equatorial data */
      if (fabs(qdlat) < 10.0)
        wi = SECS2D_WEIGHT_EEJ;

      /* set rhs vector */
      gsl_vector_set(state->rhs, rowidx, B[0]);
      gsl_vector_set(state->rhs, rowidx + 1, B[1]);
      gsl_vector_set(state->rhs, rowidx + 2, B[2]);

      /* set weight vector */
      gsl_vector_set(state->wts, rowidx, wi);
      gsl_vector_set(state->wts, rowidx + 1, wi);
      gsl_vector_set(state->wts, rowidx + 2, wi);

      /* build 3 rows of the LS matrix for DF SECS */
      secs2d_matrix_row_df(r, theta, phi, &vx.vector, &vy.vector, &vz.vector, state);

      rowidx += 3;
    }

  if (state->flags & MAGFIT_FLG_SECS_FIT_CF)
    {
      /*XXX*/
    }

  state->n = rowidx;

  return GSL_SUCCESS;
}

/*
secs2d_fit()
  Fit 1D SECS to previously added tracks

Inputs: rnorm  - residual norm || y - A x ||
        snorm  - solution norm || y - A x ||
        vstate - state

Return: success/error

Notes:
1) Data must be added to workspace via
secs2d_add_datum()
*/

static int
secs2d_fit(double * rnorm, double * snorm, void * vstate)
{
  secs2d_state_t *state = (secs2d_state_t *) vstate;
  const size_t npts = 200;
  /* Note: to get a reasonable current map, use tol = 3e-1 */
#if 0
  const double tol = 3.0e-1;
#else
  const double tol = 1.0e-1;
#endif
  gsl_vector *reg_param = gsl_vector_alloc(npts);
  gsl_vector *rho = gsl_vector_alloc(npts);
  gsl_vector *eta = gsl_vector_alloc(npts);
  gsl_vector *G = gsl_vector_alloc(npts);
  gsl_matrix_view A = gsl_matrix_submatrix(state->X, 0, 0, state->n, state->p);
  gsl_vector_view b = gsl_vector_subvector(state->rhs, 0, state->n);
  gsl_vector_view wts = gsl_vector_subvector(state->wts, 0, state->n);
  double lambda_gcv, lambda_l, lambda, G_gcv;
  size_t i;
  const char *lambda_file = "lambda.dat";
  FILE *fp = fopen(lambda_file, "w");
  double s0; /* largest singular value */

  if (state->n < state->p)
    return -1;

  fprintf(stderr, "\n");
  fprintf(stderr, "\t n = %zu\n", state->n);
  fprintf(stderr, "\t p = %zu\n", state->p);

#if 0 /* TSVD */

  {
    double chisq;
    size_t rank;

    gsl_multifit_wlinear_tsvd(&A.matrix, &wts.vector, &b.vector, tol, state->c, state->cov,
                              &chisq, &rank, state->multifit_p);

    *rnorm = sqrt(chisq);
    *snorm = gsl_blas_dnrm2(state->c);

    fprintf(stderr, "secs2d_fit: rank = %zu/%zu\n", rank, state->p);
  }

#else /* Tikhonov / L-curve */

  /* convert to standard form */
  gsl_multifit_linear_applyW(&A.matrix, &wts.vector, &b.vector, &A.matrix, &b.vector);

  fprintf(stderr, "\t computing SVD...");

  /* compute SVD of A */
  gsl_multifit_linear_svd(&A.matrix, state->multifit_p);
  s0 = gsl_vector_get(state->multifit_p->S, 0);

  fprintf(stderr, "done\n");

  /* compute GCV curve */
  gsl_multifit_linear_gcv(&b.vector, reg_param, G, &lambda_gcv, &G_gcv, state->multifit_p);

  /* compute L-curve */
  gsl_multifit_linear_lcurve(&b.vector, reg_param, rho, eta, state->multifit_p);

  fprintf(stderr, "\t secs2d_fit: writing %s...", lambda_file);

  for (i = 0; i < npts; ++i)
    {
      fprintf(fp, "%e %e %e %e\n",
              gsl_vector_get(reg_param, i),
              gsl_vector_get(rho, i),
              gsl_vector_get(eta, i),
              gsl_vector_get(G, i));
    }

  fprintf(stderr, "done\n");

  gsl_multifit_linear_lcorner(rho, eta, &i);
  lambda_l = gsl_vector_get(reg_param, i);

  /* lower bound on lambda */
  lambda = GSL_MAX(lambda_l, tol * s0);

  /* solve regularized system with lambda */
  gsl_multifit_linear_solve(lambda, &A.matrix, &b.vector, state->c, rnorm, snorm, state->multifit_p);

  fprintf(stderr, "\t s0 = %.12e\n", s0);
  fprintf(stderr, "\t lambda_l = %.12e\n", lambda_l);
  fprintf(stderr, "\t lambda_gcv = %.12e\n", lambda_gcv);
  fprintf(stderr, "\t lambda = %.12e\n", lambda);
  fprintf(stderr, "\t rnorm = %.12e\n", *rnorm);
  fprintf(stderr, "\t snorm = %.12e\n", *snorm);
  fprintf(stderr, "\t cond(X) = %.12e\n", 1.0 / gsl_multifit_linear_rcond(state->multifit_p));

#endif

  printv_octave(state->c, "c");

  gsl_vector_free(reg_param);
  gsl_vector_free(rho);
  gsl_vector_free(eta);
  gsl_vector_free(G);

  fclose(fp);

  return 0;
}

/*
secs2d_eval_B()
  Evaluate magnetic field at a given (r,theta) using
previously computed 1D SECS coefficients

Inputs: t      - timestamp (CDF_EPOCH)
        r      - radius (km)
        theta  - colatitude (radians)
        phi    - longitude (radians)
        B      - (output) magnetic field vector (nT)
        vstate - state

Notes:
1) state->c must contain 1D SECS coefficients
*/

static int
secs2d_eval_B(const double t, const double r, const double theta, const double phi,
              double B[3], void * vstate)
{
  secs2d_state_t *state = (secs2d_state_t *) vstate;
  gsl_vector_view vx = gsl_matrix_row(state->X, 0);
  gsl_vector_view vy = gsl_matrix_row(state->X, 1);
  gsl_vector_view vz = gsl_matrix_row(state->X, 2);

  (void) t;   /* unused parameter */

  B[0] = 0.0;
  B[1] = 0.0;
  B[2] = 0.0;

  if (state->flags & MAGFIT_FLG_SECS_FIT_DF)
    {
      secs2d_matrix_row_df(r, theta, phi, &vx.vector, &vy.vector, &vz.vector, state);

      gsl_blas_ddot(&vx.vector, state->c, &B[0]);
      gsl_blas_ddot(&vy.vector, state->c, &B[1]);
      gsl_blas_ddot(&vz.vector, state->c, &B[2]);
    }

  if (state->flags & MAGFIT_FLG_SECS_FIT_CF)
    {
      secs2d_matrix_row_cf(r, theta, &vy.vector, state);

      gsl_blas_ddot(&vy.vector, state->c, &B[1]);
    }

  return 0;
}

/*
secs2d_eval_J()
  Evaluate current density at a given (r,theta) using
previously computed 1D SECS coefficients

Inputs: r      - radius (km)
        theta  - colatitude (radians)
        phi    - longitude (radians)
        J      - (output) current density vector [A/km]
        vstate - workspace

Notes:
1) state->c must contain 1D SECS coefficients
*/

static int
secs2d_eval_J(const double r, const double theta, const double phi,
              double J[3], void * vstate)
{
  secs2d_state_t *state = (secs2d_state_t *) vstate;
  gsl_vector_view vx = gsl_matrix_row(state->X, 0);
  gsl_vector_view vy = gsl_matrix_row(state->X, 1);
  gsl_vector_view vz = gsl_matrix_row(state->X, 2);
  size_t i;

  (void) r;   /* unused parameter */

  J[0] = 0.0;
  J[1] = 0.0;
  J[2] = 0.0;

  if (state->flags & MAGFIT_FLG_SECS_FIT_DF)
    {
      secs2d_matrix_row_df_J(theta, phi, &vx.vector, &vy.vector, &vz.vector, state);

      gsl_blas_ddot(&vx.vector, state->c, &J[0]);
      gsl_blas_ddot(&vy.vector, state->c, &J[1]);
      gsl_blas_ddot(&vz.vector, state->c, &J[2]);
    }

  if (state->flags & MAGFIT_FLG_SECS_FIT_CF)
    {
      secs2d_matrix_row_cf_J(theta, &vx.vector, &vz.vector, state);

      gsl_blas_ddot(&vx.vector, state->c, &J[0]);
      gsl_blas_ddot(&vz.vector, state->c, &J[2]);
    }

  /* convert to units of A/km */
  for (i = 0; i < 3; ++i)
    J[i] *= 1.0e3;

  return 0;
}

/*
secs2d_green_df()
  Compute magnetic field Green's function for a single divergence-free
2D SECS

Inputs: r        - radius (km)
        theta    - colatitude (radians)
        phi      - longitude (radians)
        theta0   - colatitude of SECS pole (radians)
        phi0     - longitude of SECS pole (radians)
        B        - (output) magnetic field Green's function (X,Y,Z)
        w        - workspace

Return: success/error
*/

static int
secs2d_green_df(const double r, const double theta, const double phi,
                const double theta0, const double phi0,
                double B[3], secs2d_state_t *state)
{
  const double s = GSL_MIN(r, state->R) / GSL_MAX(r, state->R);
  double thetap, stp, ctp, d;
  size_t i;

  secs2d_transform(theta0, phi0, theta, phi, &thetap);
  stp = sin(thetap);
  ctp = cos(thetap);
  d = sqrt(1.0 + s*s - 2.0*s*ctp);

  if (r >= state->R)
    {
      /* observation point above ionosphere layer */
      B[0] = ((1.0 - s * ctp) / d - 1.0) / stp;
      B[1] = 0.0;
      B[2] = -s * (1.0 / d - 1.0);
    }
  else
    {
      /* observation point below ionosphere layer */
      B[0] = ((s - ctp) / d + ctp) / stp;
      B[1] = 0.0;
      B[2] = 1.0 - 1.0 / d;
    }

  for (i = 0; i < 3; ++i)
    B[i] *= MAGFIT_MU_0 / (4.0 * M_PI * r);

  /* transform from SECS-centered system to geographic */
  secs2d_transform_vec(theta0, phi0, theta, phi, B, B);

  return GSL_SUCCESS;
}

/*
secs2d_green_df_J()
  Compute current density Green's function for a single divergence-free
1D SECS

Inputs: theta  - colatitude (radians)
        phi    - longitude (radians)
        theta0 - pole colatitude position (radians)
        phi0   - pole longitude position (radians)
        K      - (output) current density Green's function (X,Y,Z)
        w      - workspace

Return: success/error
*/

static int
secs2d_green_df_J(const double theta, const double phi, const double theta0, const double phi0,
                  double K[3], secs2d_state_t *state)
{
  double thetap;

  secs2d_transform(theta0, phi0, theta, phi, &thetap);

  K[0] = 0.0;
  K[1] = 1.0 / (4.0 * M_PI * state->R) * tan(0.5 * (M_PI - thetap));
  K[2] = 0.0;

  /* transform from SECS-centered system to geographic */
  secs2d_transform_vec(theta0, phi0, theta, phi, K, K);

  return GSL_SUCCESS;
}

/*
secs2d_green_cf()
  Compute magnetic field Green's function for a single divergence-free
1D SECS

Inputs: r      - radius (km)
        theta  - colatitude (radians)
        theta0 - pole position (radians)
        B      - (output) magnetic field Green's function (X,Y,Z)
        w      - workspace

Return: success/error
*/

int
secs2d_green_cf(const double r, const double theta, const double theta0,
                double B[3], secs2d_state_t *state)
{
  /*XXX*/
  B[0] = 0.0;
  B[1] = 0.0;
  B[2] = 0.0;

  return GSL_SUCCESS;
}

/*
secs2d_green_cf_J()
  Compute current density Green's function for a single curl-free
1D SECS

Inputs: r      - radius (km)
        theta  - colatitude (radians)
        theta0 - pole position (radians)
        K      - (output) current density Green's function (X,Y,Z)
        w      - workspace

Return: success/error
*/

int
secs2d_green_cf_J(const double r, const double theta, const double theta0, double K[3], secs2d_state_t *state)
{
  /*XXX*/
  K[0] = 0.0;
  K[1] = 0.0;
  K[2] = 0.0;

  if (r >= state->R)
    K[2] = -1.0 / (2.0 * r * r);

  return GSL_SUCCESS;
}

/*
secs2d_matrix_row_df()
  Build matrix rows corresponding to DF SECS

Inputs: r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        X     - (output) X Green's functions for B_df
        Y     - (output) Y Green's functions for B_df
        Z     - (output) Z Green's functions for B_df
        w     - workspace
*/

static int
secs2d_matrix_row_df(const double r, const double theta, const double phi,
                     gsl_vector *X, gsl_vector *Y, gsl_vector *Z, secs2d_state_t *state)
{
  size_t i, j;

  gsl_vector_set_zero(X);
  gsl_vector_set_zero(Y);
  gsl_vector_set_zero(Z);

  for (i = 0; i < state->ntheta; ++i)
    {
      double theta0 = state->theta0[i];

      for (j = 0; j < state->nphi; ++j)
        {
          double phi0 = state->phi0[j];
          size_t pole_idx = secs2d_idx(i, j, state);
          double B[3];

          /*
           * compute divergence-free 2D SECS Green's functions for pole (theta0,phi0)
           * at location (r,theta,phi)
           */
          secs2d_green_df(r, theta, phi, theta0, phi0, B, state);

          gsl_vector_set(X, state->df_offset + pole_idx, B[0]);
          gsl_vector_set(Y, state->df_offset + pole_idx, B[1]);
          gsl_vector_set(Z, state->df_offset + pole_idx, B[2]);
        }
    }

  return 0;
}

/*
secs2d_matrix_row_cf()
  Build matrix rows corresponding to CF SECS

Inputs: r - radius (km)
        theta - colatitude (radians)
        Y     - (output) Y Green's functions for B_cf
        w     - workspace

Notes:
1) B field for CF 1D SECS only has Y component
*/

static int
secs2d_matrix_row_cf(const double r, const double theta,
                    gsl_vector *Y, secs2d_state_t *state)
{
  size_t i;

  gsl_vector_set_zero(Y);

  for (i = 0; i < state->ntheta; ++i)
    {
      double theta0 = state->theta0[i];
      double B[3];

      /* compute curl-free 1D SECS Green's functions */
      secs2d_green_cf(r, theta, theta0, B, state);

      gsl_vector_set(Y, state->cf_offset + i, B[1]);
    }

  return 0;
}

/*
secs2d_matrix_row_df_J()
  Build matrix rows corresponding to DF SECS (current)

Inputs: theta - colatitude (radians)
        phi   - longitude (radians)
        X     - (output) X Green's functions for J_df
        Y     - (output) Y Green's functions for J_df
        Z     - (output) Z Green's functions for J_df
        w     - workspace

Notes:
1) J vector for DF 1D SECS does not have X, Z components
*/

static int
secs2d_matrix_row_df_J(const double theta, const double phi,
                       gsl_vector *X, gsl_vector *Y, gsl_vector *Z,
                       secs2d_state_t *state)
{
  size_t i, j;

  gsl_vector_set_zero(X);
  gsl_vector_set_zero(Y);
  gsl_vector_set_zero(Z);

  for (i = 0; i < state->ntheta; ++i)
    {
      double theta0 = state->theta0[i];

      for (j = 0; j < state->nphi; ++j)
        {
          double phi0 = state->phi0[j];
          size_t pole_idx = secs2d_idx(i, j, state);
          double J[3];

          /* compute divergence-free 1D SECS Green's functions */
          secs2d_green_df_J(theta, phi, theta0, phi0, J, state);

          gsl_vector_set(X, state->df_offset + pole_idx, J[0]);
          gsl_vector_set(Y, state->df_offset + pole_idx, J[1]);
          gsl_vector_set(Z, state->df_offset + pole_idx, J[2]);
        }
    }

  return 0;
}

/*
secs2d_matrix_row_cf_J()
  Build matrix rows corresponding to CF SECS (current)

Inputs: theta - colatitude (radians)
        X     - (output) X Green's functions for J_cf
        Z     - (output) Z Green's functions for J_cf
        w     - workspace

Notes:
1) J vector for CF 1D SECS does not have a Y component
*/

static int
secs2d_matrix_row_cf_J(const double theta, gsl_vector *X, gsl_vector *Z, secs2d_state_t *state)
{
  size_t i;

  gsl_vector_set_zero(X);
  gsl_vector_set_zero(Z);

  for (i = 0; i < state->ntheta; ++i)
    {
      double theta0 = state->theta0[i];
      double J[3];

      /* compute curl-free 1D SECS Green's functions */
      secs2d_green_cf_J(state->R, theta, theta0, J, state);

      gsl_vector_set(X, state->cf_offset + i, J[0]);
      gsl_vector_set(Z, state->cf_offset + i, J[2]);
    }

  return 0;
}

/*
secs2d_transform()
  Given a geographic location (theta,phi), find the point (thetap,phip) in
the SECS-centered system with a pole at geographic location (theta0,phi0)

Inputs: theta0   - geographic colatitude of SECS pole (radians)
        phi0     - geographic longitude of SECS pole (radians)
        theta    - geographic colatitude of desired location
        phi      - geographic longitude of desired location
        thetap   - (output) colatitude of point (theta,phi) in SECS-centered system

Return: success/error
*/

static int
secs2d_transform(const double theta0, const double phi0, const double theta, const double phi,
                 double *thetap)
{
  int s = 0;
  const double ctel = cos(theta0);
  const double stel = sin(theta0);
  const double ct = cos(theta);
  const double st = sin(theta);
  const double dphi = phi0 - phi;
  const double ctp = ct * ctel + st * stel * cos(dphi); /* cos(theta') */

  *thetap = acos(ctp);

  return s;
}

/*
secs2d_transform_vec()
  Transform a vector with components in SECS-centered system into a geographic
system. The point in geographic coordinates is (theta,phi).

Inputs: theta0   - geographic colatitude of SECS pole (radians)
        phi0     - geographic longitude of SECS pole (radians)
        theta    - geographic colatitude of desired location
        phi      - geographic longitude of desired location
        v_el     - 3-vector in SECS-centered system (theta0,phi0) (X',Y',Z')
        v_geo    - (output) 3-vector in geographic system (X,Y,Z)

Return: success/error

Notes:
1) It is allowed for v_geo = v_el for an in-place transform
*/

static int
secs2d_transform_vec(const double theta0, const double phi0, const double theta, const double phi,
                     const double v_el[3], double v_geo[3])
{
  int s = 0;
  const double ctel = cos(theta0);
  const double stel = sin(theta0);
  const double ct = cos(theta);
  const double st = sin(theta);
  const double dphi = phi0 - phi;
  const double ctp = ct * ctel + st * stel * cos(dphi); /* cos(theta') */
  const double stp = sqrt(1.0 - ctp * ctp);             /* sin(theta') */
  const double cosC = (ctel - ct * ctp) / (st * stp);   /* cos(C) */
  const double sinC = stel * sin(dphi) / stp;           /* sin(C) */
  const double Xp = v_el[0];
  const double Yp = v_el[1];

  v_geo[0] = Xp * cosC - Yp * sinC;
  v_geo[1] = Xp * sinC + Yp * cosC;
  v_geo[2] = v_el[2];

  return s;
}

/*
secs2d_idx()
  Return linear index from pole (theta0,phi0) represented
by indices (i,j)

Inputs: i     - theta index in [0,ntheta-1]
        j     - phi index in [0,nphi-1]
        state - state
*/

static size_t
secs2d_idx(const size_t i, const size_t j, const secs2d_state_t *state)
{
  size_t idx = CIDX2(i, state->ntheta, j, state->nphi);
  return idx;
}

static const magfit_type secs2d_type =
{
  "secs2d",
  secs2d_alloc,
  secs2d_reset,
  secs2d_ncoeff,
  secs2d_add_datum,
  secs2d_fit,
  secs2d_eval_B,
  secs2d_eval_J,
  NULL,
  NULL,
  secs2d_free
};

const magfit_type *magfit_secs2d = &secs2d_type;
