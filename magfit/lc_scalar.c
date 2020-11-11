/*
 * lc_scalar.c
 *
 * Fit scalar magnetic field data to line current model
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

#include <mainlib/ml_apex.h>
#include <mainlib/ml_satdata.h>
#include <mainlib/ml_indices.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_geo.h>
#include <mainlib/ml_interp.h>
#include <mainlib/ml_track.h>
#include <mainlib/ml_bsearch.h>

#include "magfit.h"

typedef struct
{
  size_t n;         /* total number of measurements in system */
  size_t p;         /* number of line currents */
  size_t nmax;      /* maximum number of measurements in LS system */
  size_t flags;     /* MAGFIT_FLG_xxx */
  double R;         /* radius of line currents (km) */
  double epoch;     /* epoch (year) of current inversion dataset */

  size_t nlon;
  double dqdlat;    /* spacing between line currents in QD latitude (degrees) */
  double qdlat_max; /* maximum allowed QD latitude (degrees) */

  gsl_matrix *Jx;   /* ncurr-by-nlon */
  gsl_matrix *Jy;   /* ncurr-by-nlon */
  gsl_matrix *Jz;   /* ncurr-by-nlon */
  gsl_matrix *mid_pos_x; /* ncurr-by-nlon */
  gsl_matrix *mid_pos_y; /* ncurr-by-nlon */
  gsl_matrix *mid_pos_z; /* ncurr-by-nlon */
  double * mid_lon;      /* geographic longitude of segment midpoints, length nlon */

  gsl_matrix *X;    /* LS matrix */
  gsl_vector *rhs;  /* rhs vector */
  gsl_vector *wts;  /* weight vector */
  gsl_matrix *cov;  /* covariance matrix */
  gsl_matrix *L;    /* regularization matrix */
  gsl_vector *Ltau; /* regularization matrix Householder scalars */
  gsl_matrix *M;    /* workspace matrix */
  gsl_matrix *Xs;   /* standard form X */
  gsl_vector *bs;   /* standard form b */
  gsl_vector *cs;   /* standard form c */

  gsl_vector *c;    /* solution vector */

  gsl_vector *reg_param; /* vector of regularization parameters, length n_lcurve */
  gsl_vector *rho;       /* vector of residual norms, length n_lcurve */
  gsl_vector *eta;       /* vector of solution norms, length n_lcurve */
  size_t reg_idx;        /* index of chosen regularization parameter */

  double Rsq;            /* coefficient of determination of fit */

  gsl_multifit_linear_workspace *multifit_p;
  apex_workspace * apex_workspace_p;
} lcs_state_t;

static void *lcs_alloc(const void * params);
static void lcs_free(void * vstate);
static int lcs_reset(void * vstate);
static size_t lcs_ncoeff(void * vstate);
static int lcs_add_datum(const double t, const double r, const double theta, const double phi,
                         const double qdlat, const double B[4], void * vstate);
static int lcs_fit(double * rnorm, double * snorm, void * vstate);
static int lcs_dofit(double * rnorm, double * snorm, void * vstate);
static int lcs_eval_B(const double t, const double r, const double theta, const double phi,
                      double B[4], void * vstate);
static int lcs_eval_J(const double r, const double theta, const double phi,
                      double J[3], void * vstate);
static int lcs_matrix_row(const double r, const double theta, const double phi,
                          const double b[3], gsl_vector *v, const lcs_state_t * state);
static int lcs_add_B(const double ri[3], const double rjk[3], const double Jjk[3],
                     double Bij[3]);
static void my_delaz(double lat1, double lon1, double lat2, double lon2,
                     double *delta, double *az);


/*
lcs_alloc()
  Allocate scalar line current workspace

Inputs: params - fit parameters

Return: pointer to workspace
*/

static void *
lcs_alloc(const void * params)
{
  const magfit_parameters *mparams = (const magfit_parameters *) params;
  lcs_state_t *state;
  const double qd_alt = mparams->R - R_EARTH_KM;
  const double qdlon_min = 0.0;
  const double qdlon_max = 359.0;
  const size_t nlon = 360;
  const double dlon = (qdlon_max - qdlon_min) / (nlon - 1.0);
  const size_t ncurr = mparams->lc_ncurr;
  const double qdlat_max = mparams->lc_qdmax;
  const double dqdlat = (2.0 * qdlat_max) / (ncurr - 1.0);
  const size_t flags = mparams->flags;
  size_t j, k;

  state = calloc(1, sizeof(lcs_state_t));
  if (!state)
    return 0;

  state->nmax = 30000;
  state->n = 0;
  state->p = ncurr;
  state->R = mparams->R;
  state->flags = flags;
  state->nlon = nlon;
  state->dqdlat = dqdlat;
  state->qdlat_max = qdlat_max;
  state->epoch = (double) mparams->lc_year;
  state->Rsq = 0.0;

  if (state->p == 0)
    {
      GSL_ERROR_NULL("must specify number of line currents", GSL_EINVAL);
    }

  state->Jx = gsl_matrix_alloc(ncurr, nlon);
  state->Jy = gsl_matrix_alloc(ncurr, nlon);
  state->Jz = gsl_matrix_alloc(ncurr, nlon);
  state->mid_pos_x = gsl_matrix_alloc(ncurr, nlon);
  state->mid_pos_y = gsl_matrix_alloc(ncurr, nlon);
  state->mid_pos_z = gsl_matrix_alloc(ncurr, nlon);
  state->mid_lon = malloc(nlon * sizeof(double));

  state->X = gsl_matrix_alloc(state->nmax, state->p);
  state->c = gsl_vector_alloc(state->p);
  state->rhs = gsl_vector_alloc(state->nmax);
  state->wts = gsl_vector_alloc(state->nmax);
  state->cov = gsl_matrix_alloc(state->p, state->p);
  state->M = gsl_matrix_alloc(state->nmax, state->p);
  state->Xs = gsl_matrix_alloc(state->nmax, state->p);
  state->bs = gsl_vector_alloc(state->nmax);
  state->cs = gsl_vector_alloc(state->p);
  state->multifit_p = gsl_multifit_linear_alloc(state->nmax, state->p);

  state->reg_param = gsl_vector_alloc(mparams->n_lcurve);
  state->rho = gsl_vector_alloc(mparams->n_lcurve);
  state->eta = gsl_vector_alloc(mparams->n_lcurve);

  state->apex_workspace_p = apex_alloc();

  /* regularization matrix L */
  {
    /* single k-derivative norm */
    const size_t k = 2;
    const size_t m = state->p - k; /* L is m-by-p */

    state->L = gsl_matrix_alloc(m, state->p);
    state->Ltau = gsl_vector_alloc(m);

    gsl_multifit_linear_Lk(state->p, k, state->L);
    gsl_multifit_linear_L_decomp(state->L, state->Ltau);
  }

  /*
   * precompute Cartesian positions of each current segment for each
   * line current arc; current arcs follow lines of constant QD
   * latitude and segments are equally spaced in QD longitude. These
   * QD grid points are converted to spherical geocentric points and then
   * to Cartesian.
   *
   * XXX - 24 Jan 2014 - the apex_transform_inv() calls here produce
   * geographic latitudes which aren't as "smooth" as Stefan/GFZ modified
   * apxntrp.f calls - they tend to cluster closely together for
   * certain qdlat inputs as can be seen by plotting
   * (mid_lon(k), lat2-lat1) in the below loop
   */
  for (j = 0; j < ncurr; ++j)
    {
      double qdlat = -qdlat_max + j * dqdlat;

      for (k = 0; k < nlon; ++k)
        {
          double qdlon = qdlon_min + (k + 0.5) * dlon; /* midpoint lon */
          double qdlon1 = qdlon_min + k * dlon;        /* endpoint lons */
          double qdlon2 = qdlon_min + (k + 1.0) * dlon;
          double glat, glon; /* geodetic lat/lon */
          double glat1, glat2, glon1, glon2;
          double r1, r2, lat1, lat2;
          double lat, r;     /* geographic spherical lat/radius */
          double delta, az, delta_km;
          double mid_pos[3]; /* ECEF Cartesian position of segment midpt */
          double J_sph[3];   /* components of J in spherical coords */
          double J_ecef[3];  /* components of J in ECEF Cartesian */

          /* compute geodetic lon/lat of qd point */
          apex_transform_inv_geodetic(state->epoch, qdlat, qdlon, qd_alt,
                                      &glat, &glon, state->apex_workspace_p);

          /* convert geodetic to geocentric spherical */
          geodetic2geo(glat, qd_alt, &lat, &r);

          /* sanity check */
          {
            double tmplat, tmpalt;
            geo2geodetic(lat, glon, r, &tmplat, &tmpalt);
            gsl_test_rel(tmplat, glat, 1.0e-7, "latitude");
            gsl_test_rel(tmpalt, qd_alt, 1.0e-6, "altitude");
          }

          /*
           * now (r,lat,glon) contain the geocentric spherical coordinate
           * of the midpoint of segment k of current j; save ECEF
           * Cartesian position of the midpoints
           */
          sph2ecef(r, M_PI / 2.0 - lat, glon, mid_pos);
          gsl_matrix_set(state->mid_pos_x, j, k, mid_pos[0]);
          gsl_matrix_set(state->mid_pos_y, j, k, mid_pos[1]);
          gsl_matrix_set(state->mid_pos_z, j, k, mid_pos[2]);

          /*
           * save segment longitude to check later if its within 30 deg;
           * this array is overwritten for each current j so its
           * an approximate position
           */
          state->mid_lon[k] = glon * 180.0 / M_PI;

          /* find geocentric coordinates of segment endpoints */
          apex_transform_inv_geodetic(state->epoch, qdlat, qdlon1, qd_alt,
                                      &glat1, &glon1, state->apex_workspace_p);
          apex_transform_inv_geodetic(state->epoch, qdlat, qdlon2, qd_alt,
                                      &glat2, &glon2, state->apex_workspace_p);

          /* convert geodetic to geocentric spherical */
          geodetic2geo(glat1, qd_alt, &lat1, &r1);
          geodetic2geo(glat2, qd_alt, &lat2, &r2);

          lat1 *= 180.0 / M_PI;
          glon1 *= 180.0 / M_PI;
          lat2 *= 180.0 / M_PI;
          glon2 *= 180.0 / M_PI;

          my_delaz(lat1, glon1, lat2, glon2, &delta, &az);

          delta_km = delta * M_PI / 180.0 * r;

          /* compute current components */
          J_sph[0] = 0.0;
          J_sph[1] = -delta_km * cos(az * M_PI / 180.0);
          J_sph[2] = delta_km * sin(az * M_PI / 180.0);

          sph2ecef_vec(M_PI / 2.0 - lat1 * M_PI / 180.0, glon1 * M_PI / 180.0, J_sph, J_ecef);

          /* store current vector "green function" for later use */
          gsl_matrix_set(state->Jx, j, k, J_ecef[0]);
          gsl_matrix_set(state->Jy, j, k, J_ecef[1]);
          gsl_matrix_set(state->Jz, j, k, J_ecef[2]);
        }
    }

  return state;
}

static void
lcs_free(void * vstate)
{
  lcs_state_t *state = (lcs_state_t *) vstate;

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

  if (state->L)
    gsl_matrix_free(state->L);

  if (state->Ltau)
    gsl_vector_free(state->Ltau);

  if (state->M)
    gsl_matrix_free(state->M);

  if (state->Xs)
    gsl_matrix_free(state->Xs);

  if (state->bs)
    gsl_vector_free(state->bs);

  if (state->cs)
    gsl_vector_free(state->cs);

  if (state->reg_param)
    gsl_vector_free(state->reg_param);

  if (state->rho)
    gsl_vector_free(state->rho);

  if (state->eta)
    gsl_vector_free(state->eta);

  if (state->multifit_p)
    gsl_multifit_linear_free(state->multifit_p);

  if (state->apex_workspace_p)
    apex_free(state->apex_workspace_p);

  free(state);
}

/* reset workspace to work on new data set */
static int
lcs_reset(void * vstate)
{
  lcs_state_t *state = (lcs_state_t *) vstate;
  state->n = 0;
  return 0;
}

static size_t
lcs_ncoeff(void * vstate)
{
  lcs_state_t *state = (lcs_state_t *) vstate;
  return state->p;
}

/*
lcs_add_datum()
  Add single vector measurement to LS system

Inputs: t      - timestamp (CDF_EPOCH)
        r      - radius (km)
        theta  - colatitude (radians)
        phi    - longitude (radians)
        qdlat  - QD latitude (degrees)
        B      - unit magnetic field vector NEC (nT) and scalar measurement
                 B = [ x ]
                     [ y ]
                     [ z ]
                     [ F ]
        vstate - state

Return: success/error

Notes:
1) state->n is updated with the number of total data added
*/

static int
lcs_add_datum(const double t, const double r, const double theta, const double phi,
              const double qdlat, const double B[4], void * vstate)
{
  lcs_state_t *state = (lcs_state_t *) vstate;
  size_t rowidx = state->n;
  gsl_vector_view v;
  double wi = 1.0;

  /* ignore data outside [-qdmax,qdmax] */
  if (fabs(qdlat) > state->qdlat_max)
    return 0;

  /* add |F| to rhs vector */
  gsl_vector_set(state->rhs, rowidx, B[3]);

  /* set weight */
  gsl_vector_set(state->wts, rowidx, wi);

  /* construct row of LS matrix */
  v = gsl_matrix_row(state->X, rowidx);
  lcs_matrix_row(r, theta, phi, B, &v.vector, state);

  ++rowidx;
  state->n = rowidx;

  return GSL_SUCCESS;
}

/*
lcs_fit()
  Fit line current model to previously added data

Inputs: rnorm  - residual norm || y - A x ||
        snorm  - solution norm || L x ||
        vstate - state

Return: success/error

Notes:
1) Data must be added to workspace via lcs_add_datum()
2) R^2 is computed and stored in state->Rsq
*/

static int
lcs_fit(double * rnorm, double * snorm, void * vstate)
{
  int status;
  const size_t niter = 1;
  lcs_state_t *state = (lcs_state_t *) vstate;
  gsl_matrix_view A = gsl_matrix_submatrix(state->X, 0, 0, state->n, state->p);
  gsl_vector_view b = gsl_vector_subvector(state->rhs, 0, state->n);
  gsl_vector_view w = gsl_vector_subvector(state->wts, 0, state->n);
  gsl_vector *residual = gsl_vector_alloc(state->n);
  gsl_vector *wts = gsl_vector_alloc(state->n);
  double rms, tss;
  size_t i, iter;

  gsl_vector_memcpy(wts, &w.vector);

  for (iter = 0; iter < niter; ++iter)
    {
      if (iter > 0)
        {
          /* compute residuals from previous iteration */
          gsl_vector_memcpy(residual, &b.vector);
          gsl_blas_dgemv(CblasNoTrans, -1.0, &A.matrix, state->c, 1.0, residual);

          /* residual rms */
          rms = 1.0 / sqrt((double) state->n) * gsl_blas_dnrm2(residual);

          /* update robust weights */
          for (i = 0; i < state->n; ++i)
            {
              double wi = gsl_vector_get(wts, i);
              double ri = gsl_vector_get(residual, i);
              double ei = ri / (1.5 * rms);
              double rwi = GSL_MIN(1.0, 1.0 / fabs(ei));

              gsl_vector_set(state->wts, i, wi * rwi);
            }
        }

      lcs_dofit(rnorm, snorm, vstate);
    }

  /* compute total sum of squares */
  tss = gsl_stats_wtss(wts->data, wts->stride, state->rhs->data, state->rhs->stride, state->n);
  /*tss = gsl_stats_tss(state->rhs->data, state->rhs->stride, state->n);*/

  /* compute R^2 */
  state->Rsq = 1.0 - ((*rnorm) * (*rnorm)) / tss;
  fprintf(stderr, "tss = %e\n", tss);
  fprintf(stderr, "Rsq = %.4f\n", state->Rsq);

  gsl_vector_free(residual);
  gsl_vector_free(wts);

  return status;
}

static int
lcs_dofit(double * rnorm, double * snorm, void * vstate)
{
  lcs_state_t *state = (lcs_state_t *) vstate;
  const size_t npts = 200;
  const double tol = 1.0e-6;
  gsl_vector *G = gsl_vector_alloc(npts);
  gsl_matrix_const_view A = gsl_matrix_const_submatrix(state->X, 0, 0, state->n, state->p);
  gsl_vector_const_view b = gsl_vector_const_subvector(state->rhs, 0, state->n);
  gsl_vector_view wts = gsl_vector_subvector(state->wts, 0, state->n);
  double lambda_l, lambda;
  size_t i;
  const size_t m = state->L->size1;
  gsl_matrix_view M = gsl_matrix_submatrix(state->M, 0, 0, state->n, state->p);
  gsl_matrix_view As = gsl_matrix_submatrix(state->Xs, 0, 0, state->n - state->p + m, m);
  gsl_vector_view bs = gsl_vector_subvector(state->bs, 0, state->n - state->p + m);
  gsl_vector_view cs = gsl_vector_subvector(state->cs, 0, m);
  double s0; /* largest singular value */

  if (state->n < state->p)
    return -1;

#if 0 /* TSVD */

  {
    double chisq;
    size_t rank;

    gsl_multifit_wlinear_tsvd(&A.matrix, &wts.vector, &b.vector, tol, state->c, state->cov,
                              &chisq, &rank, state->multifit_p);

    *rnorm = sqrt(chisq);
    snorm = gsl_blas_dnrm2(state->c);

    fprintf(stderr, "lcs_fit: rank = %zu/%zu\n", rank, state->p);
  }

#else /* Tikhonov / L-curve */

  gsl_multifit_linear_wstdform2(state->L, state->Ltau, &A.matrix, &wts.vector, &b.vector,
                                &As.matrix, &bs.vector, &M.matrix, state->multifit_p);

  /* compute SVD of A */
  gsl_multifit_linear_svd(&As.matrix, state->multifit_p);
  s0 = gsl_vector_get(state->multifit_p->S, 0);

  /* compute L-curve */
  gsl_multifit_linear_lcurve(&bs.vector, state->reg_param, state->rho, state->eta, state->multifit_p);
  gsl_multifit_linear_lcorner(state->rho, state->eta, &i);

  lambda_l = gsl_vector_get(state->reg_param, i);
  lambda_l = GSL_MAX(lambda_l, 1.0e-3 * s0);

  i = bsearch_double(state->reg_param->data, lambda_l, 0, npts - 1);
  lambda = gsl_vector_get(state->reg_param, i);
  state->reg_idx = i;

  /* solve regularized system with lambda */
  gsl_multifit_linear_solve(lambda, &As.matrix, &bs.vector, &cs.vector, rnorm, snorm, state->multifit_p);

  /* convert back to general form */
  gsl_multifit_linear_wgenform2(state->L, state->Ltau, &A.matrix, &wts.vector, &b.vector, &cs.vector, &M.matrix,
                                state->c, state->multifit_p);

  fprintf(stderr, "s0 = %.12e\n", s0);
  fprintf(stderr, "lambda_l = %.12e\n", lambda_l);
  fprintf(stderr, "lambda = %.12e\n", lambda);
  fprintf(stderr, "rnorm = %.12e\n", *rnorm);
  fprintf(stderr, "snorm = %.12e\n", *snorm);
  fprintf(stderr, "cond(X) = %.12e\n", 1.0 / gsl_multifit_linear_rcond(state->multifit_p));

#endif

  gsl_vector_free(G);

  return 0;
}

/*
lcs_eval_B()
  Evaluate magnetic field at a given (r,theta,phi) from
line current model

Inputs: t      - timestamp (CDF_EPOCH)
        r      - radius (km)
        theta  - colatitude (radians)
        phi    - longitude (radians)
        B      - (input/output)
                 On input, first 3 components contain unit main field vector in NEC,
                 B_in = [ bx ]
                        [ by ]
                        [ bz ]
                        [ *  ]
                 On output, scalar magnetic field due to line current model
                 B_out = [ * ]
                         [ * ]
                         [ * ]
                         [ F ]
        vstate - state
*/

static int
lcs_eval_B(const double t, const double r, const double theta, const double phi,
           double B[4], void * vstate)
{
  lcs_state_t *state = (lcs_state_t *) vstate;
  gsl_vector_view v = gsl_matrix_row(state->X, 0);

  (void) t;   /* unused parameter */

  lcs_matrix_row(r, theta, phi, B, &v.vector, state);

  gsl_blas_ddot(&v.vector, state->c, &B[3]);

  return 0;
}

/*
lcs_eval_J()
  Evaluate current density at a given (r,theta,phi) using
line current model

Inputs: r      - radius (km)
        theta  - colatitude (radians)
        phi    - longitude (radians)
        J      - (output) current density vector in NEC
                 J = [ J_x ]
                     [ J_y ]
                     [ J_HI ] (height-integrated eastward current, A/km)
        vstate - workspace
*/

static int
lcs_eval_J(const double r, const double theta, const double phi,
           double J[3], void * vstate)
{
  lcs_state_t *state = (lcs_state_t *) vstate;
  const double curr_dist_km = state->dqdlat * (40008.0) / 360.0;
  const double qdlon_min = 0.0; /* copied from lcs_alloc */
  const double qdlon_max = 359.0;
  const size_t nlon = 360;
  const double dlon = (qdlon_max - qdlon_min) / (nlon - 1.0);
  double alon, alat, qdlat;
  size_t idx_curr, idx_lon;
  double J_ecef[3], J_sph[3], S;
  size_t i;

  J[0] = 0.0;
  J[1] = 0.0;
  J[2] = 0.0;

  apex_transform(state->epoch, theta, phi, r, &alon, &alat, &qdlat,
                 NULL, NULL, NULL, state->apex_workspace_p);

  /* check for quick return */
  if (fabs(qdlat) > state->qdlat_max)
    return 0;

  idx_curr = (size_t) ((qdlat + state->qdlat_max) / state->dqdlat);
  idx_lon = (size_t) ((wrap360(alon) - qdlon_min) / dlon);

  assert(idx_lon < nlon);
  assert(idx_curr < state->p);

  S = gsl_vector_get(state->c, idx_curr);

  J_ecef[0] = gsl_matrix_get(state->Jx, idx_curr, idx_lon);
  J_ecef[1] = gsl_matrix_get(state->Jy, idx_curr, idx_lon);
  J_ecef[2] = gsl_matrix_get(state->Jz, idx_curr, idx_lon);

  /* make unit vector */
  vec_unit(J_ecef, J_ecef);

  /* convert to spherical components */
  ecef2sph_vec(theta, phi, J_ecef, J_sph);

  /* convert to NEC and multiply by factor */
  J[0] = -S * J_sph[1];
  J[1] =  S * J_sph[2];
  J[2] = -S * J_sph[0];

  /* store height-integrated magnetic eastward current in J_z */
  {
    double s0 = gsl_vector_get(state->c, idx_curr);
    double s1 = gsl_vector_get(state->c, idx_curr + 1);
    double qd0 = -state->qdlat_max + idx_curr * state->dqdlat;
    double qd1 = -state->qdlat_max + (idx_curr + 1) * state->dqdlat;
    double val = interp1d(qd0, qd1, s0, s1, qdlat);
    J[2] = val / curr_dist_km;
  }

  /* convert to units of A/km */
  J[2] *= 1.0e3;

  return 0;
}

static int
lcs_postproc(magfit_postproc_params * postproc, void * vstate)
{
  lcs_state_t *state = (lcs_state_t *) vstate;

  postproc->reg_param = state->reg_param;
  postproc->rho = state->rho;
  postproc->eta = state->eta;
  postproc->reg_idx = state->reg_idx;
  postproc->Rsq = state->Rsq;

  return GSL_SUCCESS;
}

/*
lcs_matrix_row()
  Construct row of LS matrix

Inputs: r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        b     - unit magnetic field vector at (r,theta,phi) in NEC
                [ bx ]
                [ by ]
                [ bz ]
        v     - (output) row of LS matrix
        state - state
*/

static int
lcs_matrix_row(const double r, const double theta, const double phi,
               const double b[3], gsl_vector *v, const lcs_state_t * state)
{
  const double max_dlon = 30.0; /* maximum longitude window in deg */
  const double lon_deg = phi * 180.0 / M_PI;
  const double lat_rad = M_PI / 2.0 - theta;
  size_t j, k;
  double ri[3];

  /* compute ECEF cartesian position vector of satellite position */
  sph2ecef(r, theta, phi, ri);

  for (j = 0; j < state->p; ++j)
    {
      double Bij[3]; /* magnetic field at ri due to line j in ECEF */
      double Bij_sph[3]; /* Bij in spherical components */
      double Bij_NEC[3]; /* Bij in NEC components */
      double Fij; /* B_{ij} . b_i */

      Bij[0] = Bij[1] = Bij[2] = 0.0;

      /* compute B_{ij} = sum_k dB_{ijk} (eq 14) */
      for (k = 0; k < state->nlon; ++k)
        {
          double rjk[3], Jjk[3];
          double dlon = wrap180(state->mid_lon[k] - lon_deg);

          /* check if this segment is within max_dlon of obs point */
          if (fabs(dlon) > max_dlon)
            continue;

          rjk[0] = gsl_matrix_get(state->mid_pos_x, j, k);
          rjk[1] = gsl_matrix_get(state->mid_pos_y, j, k);
          rjk[2] = gsl_matrix_get(state->mid_pos_z, j, k);

          Jjk[0] = gsl_matrix_get(state->Jx, j, k);
          Jjk[1] = gsl_matrix_get(state->Jy, j, k);
          Jjk[2] = gsl_matrix_get(state->Jz, j, k);

          lcs_add_B(ri, rjk, Jjk, Bij);
        }

      /* convert B_{ij} to spherical coordinates */
      ecef2sph_vec(theta, phi, Bij, Bij_sph);
      Bij_NEC[0] = -Bij_sph[1];
      Bij_NEC[1] = Bij_sph[2];
      Bij_NEC[2] = -Bij_sph[0];

      /* compute F_{ij} = B_{ij} . b_i (eq 15) */
      Fij = vec_dot(Bij_NEC, b);

      /* add to least squares matrix */
      gsl_vector_set(v, j, Fij);
    } /* for (j = 0; j < state->p; ++j) */

  return 0;
}

/*
lcs_add_B()
  Add contribution of a single line current segment to a running
vector sum using Biot-Savart law. See eq 13 of paper

Inputs: ri  - observer location (ECEF)
        rjk - midpoint of longitudinal segment k of arc current j (ECEF)
        Jjk - unit current vector of flow direction for this segment (ECEF)
        Bij - (output) updated with ECEF X,Y,Z components of magnetic field
              for this line segment at observer location ri (nT)

Notes:
1) on output, Bij = Bij + dBijk
*/

static int
lcs_add_B(const double ri[3], const double rjk[3], const double Jjk[3],
          double Bij[3])
{
  /*
   * the factor of 1e9 gives Bij units of nT if current density Jjk
   * given in km
   */
  const double mu0 = (4.0 * M_PI) * 1.0e-7;
  const double C = 1.0e9 * mu0 / (4.0 * M_PI);

  double dx = ri[0] - rjk[0];
  double dy = ri[1] - rjk[1];
  double dz = ri[2] - rjk[2];
  double invr3 = pow(dx*dx + dy*dy + dz*dz, -1.5);

  Bij[0] += C * (Jjk[1]*dz - Jjk[2]*dy) * invr3;
  Bij[1] += C * (Jjk[2]*dx - Jjk[0]*dz) * invr3;
  Bij[2] += C * (Jjk[0]*dy - Jjk[1]*dx) * invr3;

  return GSL_SUCCESS;
}

static void
my_delaz(double lat1, double lon1, double lat2, double lon2,
         double *delta, double *az)
{

  /* output: 0 <= az < 360.0 */
  /* tested to be correct at high latitudes */

  double theta1,theta2,dlon,dlat,arg;

  theta1 = (90.0-lat1) * M_PI/180.0;
  theta2 = (90.0-lat2) * M_PI/180.0;
  dlon = (lon2-lon1) * M_PI/180.0;
  if (dlon >  M_PI) dlon -= 2*M_PI;
  if (dlon < -M_PI) dlon += 2*M_PI;

  dlat = (lat2-lat1) * M_PI/180.0;

  arg = sin(theta1)*sin(theta2)*cos(dlon) + cos(theta1)*cos(theta2);
  arg = GSL_MAX(-1.0, GSL_MIN(1.0,arg)); 
  *delta = acos(arg) * 180.0/M_PI;

  if (fabs(dlat) > 1.0e-10)
    {
      *az = atan(tan(sin(theta2)*dlon)/sin(dlat)) * 180.0/M_PI;
      if (dlat < 0) *az += 180.0;
    }
  else 
    {
      if (dlon > 0) *az = 90.0; 
      else *az = 270.0; 
    } 
  if (*az < 0) *az += 360.0;
}

static const magfit_type lcs_type =
{
  "line_current_scalar",
  lcs_alloc,
  lcs_reset,
  lcs_ncoeff,
  lcs_add_datum,
  lcs_fit,
  lcs_eval_B,
  lcs_eval_J,
  NULL,
  lcs_postproc,
  lcs_free
};

const magfit_type *magfit_lcs = &lcs_type;
