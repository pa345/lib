/*
 * secs1d.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <satdata/satdata.h>
#include <indices/indices.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>

#include <common/common.h>
#include <common/interp.h>

#include "track.h"
#include "magfit.h"

/* weight factor for low-latitude EEJ region */
#define SECS1D_WEIGHT_EEJ           (10.0)

typedef struct
{
  size_t n;         /* total number of measurements in system */
  size_t p;         /* p_df + p_cf */
  size_t nmax;      /* maximum number of measurements in LS system */
  size_t flags;     /* MAGFIT_FLG_xxx */
  size_t npoles;    /* number of poles */
  double dtheta;    /* spacing of 1D SECS poles in radians */

  size_t df_offset; /* offset in 'c' of DF coefficients */
  size_t cf_offset; /* offset in 'c' of CF coefficients */

  double R_iono;    /* radius of ionosphere (km) */

  size_t lmax;      /* maximum order for legendre functions */
  double *theta0;   /* pole locations */
  double *Pl1theta; /* P_{l,1}(cos(theta)) legendre functions, size lmax + 1 */
  gsl_matrix *Pltheta0; /* P_l(cos(theta0)) functions, npoles-by-(lmax + 1) */
  double *Pltheta;  /* P_l(cos(theta)) functions */

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

  /*
   * solution vector is organized as:
   * c = [ c_df ; c_cf ]
   * where c_df and c_cf are both length 'npoles'
   */
  gsl_vector *c;    /* solution vector */

  gsl_multifit_linear_workspace *multifit_p;
} secs1d_state_t;

static void *secs1d_alloc(const void * params);
static void secs1d_free(void * vstate);
static int secs1d_reset(void * vstate);
static size_t secs1d_ncoeff(void * vstate);
static int secs1d_add_datum(const double t, const double r, const double theta, const double phi,
                            const double qdlat, const double B[3], void * vstate);
static int secs1d_fit(double * rnorm, double * snorm, void * vstate);
static int secs1d_eval_B(const double t, const double r, const double theta, const double phi,
                         double B[3], void * vstate);
static int secs1d_eval_J(const double r, const double theta, const double phi,
                         double J[3], void * vstate);
static double secs1d_eval_chi(const double theta, const double phi, void * vstate);

static int build_matrix_row_df(const double r, const double theta,
                               gsl_vector *X, gsl_vector *Z, secs1d_state_t *state);
static int build_matrix_row_cf(const double r, const double theta,
                               gsl_vector *Y, secs1d_state_t *state);
static int build_matrix_row_df_J(const double theta, gsl_vector *Y, secs1d_state_t *state);
static int build_matrix_row_cf_J(const double theta, gsl_vector *X, gsl_vector *Z, secs1d_state_t *state);
static double secs1d_tancot(const double theta, const double theta0, secs1d_state_t *state);

/*
secs1d_alloc()
  Allocate secs 1d workspace

Inputs: flags        - MAGFIT_FLG_xxx
        lmax         - maximum degree for Legendre functions in expansion
        R_iono       - radius of ionosphere (km)
        pole_spacing - along-orbit latitude spacing of SECS poles (degrees)

Return: pointer to workspace
*/

static void *
secs1d_alloc(const void * params)
{
  const magfit_parameters *mparams = (const magfit_parameters *) params;
  secs1d_state_t *state;
  const size_t flags = mparams->flags;
  const double lat_spacing = mparams->lat_spacing1d;
  const double lat_min = mparams->lat_min;
  const double lat_max = mparams->lat_max;
  const size_t npoles = (size_t) ((lat_max - lat_min) / lat_spacing + 1.0);
  const double dtheta = lat_spacing * M_PI / 180.0;
  const double theta_min = M_PI / 2.0 - lat_max * M_PI / 180.0;
  size_t i;

  state = calloc(1, sizeof(secs1d_state_t));
  if (!state)
    return 0;

  state->nmax = 30000;
  state->n = 0;
  state->p = 0;
  state->npoles = npoles;
  state->R_iono = mparams->R;
  state->lmax = mparams->lmax;
  state->flags = flags;
  state->df_offset = 0;
  state->cf_offset = 0;
  state->dtheta = dtheta;

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

  if (state->p == 0)
    {
      GSL_ERROR_NULL("must specify DF or CF SECS", GSL_EINVAL);
    }

  state->X = gsl_matrix_alloc(state->nmax, state->p);
  state->c = gsl_vector_alloc(state->p);
  state->rhs = gsl_vector_alloc(state->nmax);
  state->wts = gsl_vector_alloc(state->nmax);
  state->cov = gsl_matrix_alloc(state->p, state->p);
  state->multifit_p = gsl_multifit_linear_alloc(state->nmax, state->p);

  /* regularization matrix L */
  {
#if 1
    /* single k-derivative norm */
    const size_t k = 2;
    const size_t m = state->p - k; /* L is m-by-p */
#else
    /* Sobolev norm */
    const size_t kmax = 3;
    const size_t m = state->p;
    const double alpha_data[] = { 0.1, 0.5, 1.0, 1.0 };
    gsl_vector_const_view alpha = gsl_vector_const_view_array(alpha_data, kmax + 1);
#endif

    state->L = gsl_matrix_alloc(m, state->p);
    state->Ltau = gsl_vector_alloc(m);
    state->M = gsl_matrix_alloc(state->nmax, state->p);
    state->Xs = gsl_matrix_alloc(state->nmax, state->p);
    state->bs = gsl_vector_alloc(state->nmax);
    state->cs = gsl_vector_alloc(state->p);

#if 1
    gsl_multifit_linear_Lk(state->p, k, state->L);
#else
    gsl_multifit_linear_Lsobolev(state->p, kmax, &alpha.vector,
                                 state->L, state->multifit_p);
#endif

    gsl_multifit_linear_L_decomp(state->L, state->Ltau);
  }

  state->theta0 = malloc(npoles * sizeof(double));
  state->Pltheta = malloc((state->lmax + 1) * sizeof(double));
  state->Pl1theta = malloc((state->lmax + 1) * sizeof(double));

  state->Pltheta0 = gsl_matrix_alloc(npoles, state->lmax + 1);

  /* fill in theta0 array */
  for (i = 0; i < npoles; ++i)
    {
      gsl_vector_view v = gsl_matrix_row(state->Pltheta0, i);

      state->theta0[i] = theta_min + i * dtheta;

      /* compute P_l(cos(theta0)) */
      gsl_sf_legendre_Pl_array(state->lmax, cos(state->theta0[i]), v.vector.data);
    }

  fprintf(stderr, "secs1d_alloc: npoles = %zu\n", state->npoles);
  fprintf(stderr, "secs1d_alloc: ncoeff = %zu\n", state->p);

  return state;
}

static void
secs1d_free(void * vstate)
{
  secs1d_state_t *state = (secs1d_state_t *) vstate;

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

  if (state->theta0)
    free(state->theta0);

  if (state->Pltheta0)
    gsl_matrix_free(state->Pltheta0);

  if (state->Pltheta)
    free(state->Pltheta);

  if (state->Pl1theta)
    free(state->Pl1theta);

  if (state->multifit_p)
    gsl_multifit_linear_free(state->multifit_p);

  free(state);
}

/* reset workspace to work on new data set */
static int
secs1d_reset(void * vstate)
{
  secs1d_state_t *state = (secs1d_state_t *) vstate;
  state->n = 0;
  return 0;
}

static size_t
secs1d_ncoeff(void * vstate)
{
  secs1d_state_t *state = (secs1d_state_t *) vstate;
  return state->p;
}

/*
secs1d_add_datum()
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
secs1d_add_datum(const double t, const double r, const double theta, const double phi,
                 const double qdlat, const double B[3], void * vstate)
{
  secs1d_state_t *state = (secs1d_state_t *) vstate;
  size_t rowidx = state->n;
  double wi = 1.0;

  (void) t;   /* unused parameter */
  (void) phi; /* unused parameter */

  if (state->flags & MAGFIT_FLG_SECS_FIT_DF)
    {
      gsl_vector_view vx, vz;

      if (fabs(qdlat) < 5.0)
        wi *= SECS1D_WEIGHT_EEJ;

      if (state->flags & MAGFIT_FLG_FIT_X)
        {
          vx = gsl_matrix_row(state->X, rowidx);
          gsl_vector_set(state->rhs, rowidx, B[0]);
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

      /* build rows of the LS matrix for DF SECS */
      build_matrix_row_df(r, theta, &vx.vector, &vz.vector, state);
    }

  if ((state->flags & MAGFIT_FLG_SECS_FIT_CF) &&
      (state->flags & MAGFIT_FLG_FIT_Y))
    {
      gsl_vector_view vy = gsl_matrix_row(state->X, rowidx);

      /* set rhs vector */
      gsl_vector_set(state->rhs, rowidx, B[1]);

      /* set weight vector */
      gsl_vector_set(state->wts, rowidx, wi);

      /* build 1 row of the LS matrix for CF SECS */
      build_matrix_row_cf(r, theta, &vy.vector, state);

      ++rowidx;
    }

  state->n = rowidx;

  return GSL_SUCCESS;
}

/*
secs1d_fit()
  Fit 1D SECS to previously added tracks

Inputs: rnorm  - residual norm || y - A x ||
        snorm  - solution norm || L x ||
        vstate - state

Return: success/error

Notes:
1) Data must be added to workspace via
secs1d_add_datum()
*/

static int
secs1d_fit(double * rnorm, double * snorm, void * vstate)
{
  secs1d_state_t *state = (secs1d_state_t *) vstate;
  const size_t npts = 200;
  const double tol = 1.0e-6;
  gsl_vector *reg_param = gsl_vector_alloc(npts);
  gsl_vector *rho = gsl_vector_alloc(npts);
  gsl_vector *eta = gsl_vector_alloc(npts);
  gsl_vector *G = gsl_vector_alloc(npts);
  gsl_matrix_view A = gsl_matrix_submatrix(state->X, 0, 0, state->n, state->p);
  gsl_vector_view b = gsl_vector_subvector(state->rhs, 0, state->n);
  gsl_vector_view wts = gsl_vector_subvector(state->wts, 0, state->n);
  double lambda_gcv, lambda_l, lambda, G_gcv;
  size_t i;
  const size_t m = state->L->size1;
  gsl_matrix_view M = gsl_matrix_submatrix(state->M, 0, 0, state->n, state->p);
  gsl_matrix_view As = gsl_matrix_submatrix(state->Xs, 0, 0, state->n - state->p + m, m);
  gsl_vector_view bs = gsl_vector_subvector(state->bs, 0, state->n - state->p + m);
  gsl_vector_view cs = gsl_vector_subvector(state->cs, 0, m);
  const char *lambda_file = "lambda.dat";
  FILE *fp = fopen(lambda_file, "a");
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

    fprintf(stderr, "secs1d_fit: rank = %zu/%zu\n", rank, state->p);
  }

#else /* Tikhonov / L-curve */

  gsl_multifit_linear_wstdform2(state->L, state->Ltau, &A.matrix, &wts.vector, &b.vector,
                                &As.matrix, &bs.vector, &M.matrix, state->multifit_p);

  /* compute SVD of A */
  gsl_multifit_linear_svd(&As.matrix, state->multifit_p);
  s0 = gsl_vector_get(state->multifit_p->S, 0);

  /* compute L-curve */
  gsl_multifit_linear_lcurve(&bs.vector, reg_param, rho, eta, state->multifit_p);
  gsl_multifit_linear_lcorner(rho, eta, &i);
  lambda_l = gsl_vector_get(reg_param, i);

  /* compute GCV curve */
  gsl_multifit_linear_gcv(&bs.vector, reg_param, G, &lambda_gcv, &G_gcv, state->multifit_p);

  /* the L-curve method often overdamps the system, not sure why */
  lambda_l = GSL_MAX(lambda_l, 1.0e-5 * s0);
  lambda_l = 1.0e-3*s0;

  lambda = lambda_l;
  lambda = 5.0e-3*s0; /*XXX*/

  /* solve regularized system with lambda */
  gsl_multifit_linear_solve(lambda, &As.matrix, &bs.vector, &cs.vector, rnorm, snorm, state->multifit_p);

  /* convert back to general form */
  gsl_multifit_linear_wgenform2(state->L, state->Ltau, &A.matrix, &wts.vector, &b.vector, &cs.vector, &M.matrix,
                                state->c, state->multifit_p);

  fprintf(stderr, "s0 = %.12e\n", s0);
  fprintf(stderr, "lambda_l = %.12e\n", lambda_l);
  fprintf(stderr, "lambda_gcv = %.12e\n", lambda_gcv);
  fprintf(stderr, "lambda = %.12e\n", lambda);
  fprintf(stderr, "rnorm = %.12e\n", *rnorm);
  fprintf(stderr, "snorm = %.12e\n", *snorm);
  fprintf(stderr, "cond(X) = %.12e\n", 1.0 / gsl_multifit_linear_rcond(state->multifit_p));

  fprintf(stderr, "secs1d_fit: writing %s...", lambda_file);

  for (i = 0; i < npts; ++i)
    {
      fprintf(fp, "%e %e %e %e\n",
              gsl_vector_get(reg_param, i),
              gsl_vector_get(rho, i),
              gsl_vector_get(eta, i),
              gsl_vector_get(G, i));
    }

  fprintf(fp, "\n\n");

  fprintf(stderr, "done\n");

#endif

  gsl_vector_free(reg_param);
  gsl_vector_free(rho);
  gsl_vector_free(eta);
  gsl_vector_free(G);

  fclose(fp);

  return 0;
}

/*
secs1d_eval_B()
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
secs1d_eval_B(const double t, const double r, const double theta, const double phi,
              double B[3], void * vstate)
{
  secs1d_state_t *state = (secs1d_state_t *) vstate;
  gsl_vector_view vx = gsl_matrix_row(state->X, 0);
  gsl_vector_view vy = gsl_matrix_row(state->X, 1);
  gsl_vector_view vz = gsl_matrix_row(state->X, 2);

  (void) t;   /* unused parameter */
  (void) phi; /* unused parameter */

  B[0] = 0.0;
  B[1] = 0.0;
  B[2] = 0.0;

  if (state->flags & MAGFIT_FLG_SECS_FIT_DF)
    {
      build_matrix_row_df(r, theta, &vx.vector, &vz.vector, state);

      gsl_blas_ddot(&vx.vector, state->c, &B[0]);
      gsl_blas_ddot(&vz.vector, state->c, &B[2]);
    }

  if (state->flags & MAGFIT_FLG_SECS_FIT_CF)
    {
      build_matrix_row_cf(r, theta, &vy.vector, state);

      gsl_blas_ddot(&vy.vector, state->c, &B[1]);
    }

  return 0;
}

/*
secs1d_eval_J()
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
secs1d_eval_J(const double r, const double theta, const double phi,
              double J[3], void * vstate)
{
  secs1d_state_t *state = (secs1d_state_t *) vstate;
  gsl_vector_view vx = gsl_matrix_row(state->X, 0);
  gsl_vector_view vy = gsl_matrix_row(state->X, 1);
  gsl_vector_view vz = gsl_matrix_row(state->X, 2);
  size_t i;

  (void) r;   /* unused parameter */
  (void) phi; /* unused parameter */

  J[0] = 0.0;
  J[1] = 0.0;
  J[2] = 0.0;

  if (state->flags & MAGFIT_FLG_SECS_FIT_DF)
    {
      build_matrix_row_df_J(theta, &vy.vector, state);

      gsl_blas_ddot(&vy.vector, state->c, &J[1]);
    }

  if (state->flags & MAGFIT_FLG_SECS_FIT_CF)
    {
      build_matrix_row_cf_J(theta, &vx.vector, &vz.vector, state);

      gsl_blas_ddot(&vx.vector, state->c, &J[0]);
      gsl_blas_ddot(&vz.vector, state->c, &J[2]);
    }

  /* convert to units of A/km */
  for (i = 0; i < 3; ++i)
    J[i] *= 1.0e3;

  return 0;
}

/*
secs1d_eval_chi()
  Evaluate current stream function at a given theta using
previously computed 1D SECS coefficients

Inputs: theta  - colatitude (radians)
        phi    - longitude (radians)
        vstate - workspace

Return: current stream function in kA/nT

Notes:
1) state->c must contain 1D SECS coefficients
*/

static double
secs1d_eval_chi(const double theta, const double phi, void * vstate)
{
  secs1d_state_t *state = (secs1d_state_t *) vstate;
  return 0.0;
}

int
secs1d_green_df_init(const double theta, secs1d_state_t * state)
{
  /* compute P_l(cos(theta)) */
  gsl_sf_legendre_Pl_array(state->lmax, cos(theta), state->Pltheta);

  return 0;
}

/*
secs1d_green_df()
  Compute magnetic field Green's function for a single divergence-free
1D SECS

Inputs: r        - radius (km)
        theta    - colatitude (radians)
        pole_idx - pole position [0,npoles-1]
        B        - (output) magnetic field Green's function (X,Y,Z)
        w        - workspace

Return: success/error

Notes:
1) The array w->Pltheta must be filled with P_l(cos(theta))
   prior to calling this function (see secs1d_green_df_init)

2) DF 1D SECS has no Y component, so B[1] = 0 always
*/

int
secs1d_green_df(const double r, const double theta, const size_t pole_idx,
                double B[3], secs1d_state_t *state)
{
  const double ct = cos(theta);
  const double st = sin(theta);
  gsl_vector_view Ptheta0 = gsl_matrix_row(state->Pltheta0, pole_idx);
  size_t l;

  B[0] = 0.0;
  B[1] = 0.0;
  B[2] = 0.0;

  /* compute P_l^1(cos(theta)) for all l */
  state->Pl1theta[0] = 0.0;
  state->Pl1theta[1] = -st;
  state->Pl1theta[2] = -3.0 * ct * st;

  for (l = 2; l < state->lmax; ++l)
    {
      state->Pl1theta[l + 1] = ((2.0*l + 1.0) * ct * state->Pl1theta[l] -
                            (l + 1.0) * state->Pl1theta[l - 1]) / (double) l;
    }

  if (r > state->R_iono)
    {
      double ratio = state->R_iono / r;
      double rterm = ratio; /* (R / r)^{l+1} */

      for (l = 1; l <= state->lmax; ++l)
        {
          double Pltheta0 = gsl_vector_get(&Ptheta0.vector, l);

          /* (R / r)^{l+1} */
          rterm *= ratio;

          B[0] += rterm / (l + 1.0) * Pltheta0 * state->Pl1theta[l];
          B[2] -= rterm * Pltheta0 * state->Pltheta[l];
        }

      B[0] *= MAGFIT_MU_0 / (2.0 * r);
      B[2] *= MAGFIT_MU_0 / (2.0 * r);
    }
  else
    {
      double ratio = r / state->R_iono;
      double rterm = 1.0; /* (r / R)^l */

      for (l = 1; l <= state->lmax; ++l)
        {
          double Pltheta0 = gsl_vector_get(&Ptheta0.vector, l);

          /* (r / R)^l */
          rterm *= ratio;

          B[0] -= rterm / (double)l * Pltheta0 * state->Pl1theta[l];
          B[2] -= rterm * Pltheta0 * state->Pltheta[l];
        }

      B[0] *= MAGFIT_MU_0 / (2.0 * r);
      B[2] *= MAGFIT_MU_0 / (2.0 * r);
    }

  return GSL_SUCCESS;
}

/*
secs1d_green_df_J()
  Compute current density Green's function for a single divergence-free
1D SECS

Inputs: theta  - colatitude (radians)
        theta0 - pole position (radians)
        K      - (output) current density Green's function (X,Y,Z)
        w      - workspace

Return: success/error
*/

int
secs1d_green_df_J(const double theta, const double theta0, double K[3], secs1d_state_t *state)
{
  K[0] = 0.0;
  K[1] = secs1d_tancot(theta, theta0, state) / (2.0 * state->R_iono);
  K[2] = 0.0;

  return GSL_SUCCESS;
}

/*
secs1d_green_cf()
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
secs1d_green_cf(const double r, const double theta, const double theta0,
                double B[3], secs1d_state_t *state)
{
  B[0] = 0.0;
  B[2] = 0.0;

  if (r > state->R_iono)
    {
      B[1] = -secs1d_tancot(theta, theta0, state);
      B[1] *= MAGFIT_MU_0 / (2.0 * r);
    }
  else
    {
      B[1] = 0.0;
    }

  return GSL_SUCCESS;
}

/*
secs1d_green_cf_J()
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
secs1d_green_cf_J(const double r, const double theta, const double theta0, double K[3], secs1d_state_t *state)
{
  K[0] = -secs1d_tancot(theta, theta0, state) / (2.0 * state->R_iono);
  K[1] = 0.0;
  K[2] = 0.0;

  if (r >= state->R_iono)
    K[2] = -1.0 / (2.0 * r * r);

  return GSL_SUCCESS;
}

/*
build_matrix_row_df()
  Build matrix rows corresponding to DF SECS

Inputs: r     - radius (km)
        theta - colatitude (radians)
        X     - (output) X Green's functions for B_df
        Z     - (output) Z Green's functions for B_df
        w     - workspace

Notes:
1) B field for DF 1D SECS does not have Y component
*/

static int
build_matrix_row_df(const double r, const double theta,
                    gsl_vector *X, gsl_vector *Z, secs1d_state_t *state)
{
  size_t i;

  if (state->flags & MAGFIT_FLG_FIT_X)
    gsl_vector_set_zero(X);

  if (state->flags & MAGFIT_FLG_FIT_Z)
    gsl_vector_set_zero(Z);

  /* initialize needed arrays */
  secs1d_green_df_init(theta, state);

  for (i = 0; i < state->npoles; ++i)
    {
      double B[3];

      /* compute divergence-free 1D SECS Green's functions */
      secs1d_green_df(r, theta, i, B, state);

      if (state->flags & MAGFIT_FLG_FIT_X)
        gsl_vector_set(X, state->df_offset + i, B[0]);

      if (state->flags & MAGFIT_FLG_FIT_Z)
        gsl_vector_set(Z, state->df_offset + i, B[2]);
    }

  return 0;
}

/*
build_matrix_row_cf()
  Build matrix rows corresponding to CF SECS

Inputs: r - radius (km)
        theta - colatitude (radians)
        Y     - (output) Y Green's functions for B_cf
        w     - workspace

Notes:
1) B field for CF 1D SECS only has Y component
*/

static int
build_matrix_row_cf(const double r, const double theta,
                    gsl_vector *Y, secs1d_state_t *state)
{
  size_t i;

  if (!(state->flags & MAGFIT_FLG_FIT_Y))
    return -1;

  gsl_vector_set_zero(Y);

  for (i = 0; i < state->npoles; ++i)
    {
      double theta0 = state->theta0[i];
      double B[3];

      /* compute curl-free 1D SECS Green's functions */
      secs1d_green_cf(r, theta, theta0, B, state);

      gsl_vector_set(Y, state->cf_offset + i, B[1]);
    }

  return 0;
}

/*
build_matrix_row_df_J()
  Build matrix rows corresponding to DF SECS (current)

Inputs: theta - colatitude (radians)
        Y     - (output) Y Green's functions for J_df
        w     - workspace

Notes:
1) J vector for DF 1D SECS does not have X, Z components
*/

static int
build_matrix_row_df_J(const double theta, gsl_vector *Y, secs1d_state_t *state)
{
  size_t i;

  gsl_vector_set_zero(Y);

  for (i = 0; i < state->npoles; ++i)
    {
      double theta0 = state->theta0[i];
      double J[3];

      /* compute divergence-free 1D SECS Green's functions */
      secs1d_green_df_J(theta, theta0, J, state);

      gsl_vector_set(Y, state->df_offset + i, J[1]);
    }

  return 0;
}

/*
build_matrix_row_cf_J()
  Build matrix rows corresponding to CF SECS (current)

Inputs: theta - colatitude (radians)
        X     - (output) X Green's functions for J_cf
        Z     - (output) Z Green's functions for J_cf
        w     - workspace

Notes:
1) J vector for CF 1D SECS does not have a Y component
*/

static int
build_matrix_row_cf_J(const double theta, gsl_vector *X, gsl_vector *Z, secs1d_state_t *state)
{
  size_t i;

  gsl_vector_set_zero(X);
  gsl_vector_set_zero(Z);

  for (i = 0; i < state->npoles; ++i)
    {
      double theta0 = state->theta0[i];
      double J[3];

      /* compute curl-free 1D SECS Green's functions */
      secs1d_green_cf_J(state->R_iono, theta, theta0, J, state);

      gsl_vector_set(X, state->cf_offset + i, J[0]);
      gsl_vector_set(Z, state->cf_offset + i, J[2]);
    }

  return 0;
}

/*
secs1d_tancot()
  Compute a discontinuous function common to several of the 1D SECS basis
functions:

f(theta,theta0) = { -tan(theta/2), theta < theta0
                  {  cot(theta/2), theta > theta0

using linear interpolation across the discontinuous pole region
*/

static double
secs1d_tancot(const double theta, const double theta0, secs1d_state_t *state)
{
  const double dtheta_2 = 0.5 * state->dtheta;
  const double d = theta - theta0;
  double f;

  if (fabs(d) > dtheta_2)
    {
      double cs = cos(theta);
      double sn = sin(theta);

      /*
       * use half angle formulas:
       *
       * -tan(theta/2) = (cos(theta) - 1) / sin(theta)
       *  cot(theta/2) = (cos(theta) + 1) / sin(theta)
       */
      f = (cs + GSL_SIGN(d)) / sn;
    }
  else
    {
      /* use linear interpolation near 1D SECS pole */
      double a = theta0 - dtheta_2;
      double b = theta0 + dtheta_2;
      double fa = -tan(0.5 * a);          /* -tan(a/2) */
      double fb = tan(0.5 * (M_PI - b));  /*  cot(b/2) */

      f = interp1d(a, b, fa, fb, theta);
    }

  return f;
}

static const magfit_type secs1d_type =
{
  "secs1d",
  secs1d_alloc,
  secs1d_reset,
  secs1d_ncoeff,
  secs1d_add_datum,
  secs1d_fit,
  secs1d_eval_B,
  secs1d_eval_J,
#if 0
  secs1d_eval_chi,
#else
  NULL,
#endif
  secs1d_free
};

const magfit_type *magfit_secs1d = &secs1d_type;
