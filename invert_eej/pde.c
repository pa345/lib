/*
 * pde.c
 *
 * Solve the electrodynamic PDE governing the EEJ current system
 * using finite differencing
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <signal.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_spblas.h>
#include <gsl/gsl_splinalg.h>

#include <common/common.h>
#include <common/bsearch.h>
#include <common/interp.h>
#include <common/oct.h>

#include <magfield/magfield.h>

#include "mageq.h"
#include "lisw.h"
#include "superlu.h"

#include "pde.h"
#include "sigma.h"

#include "pde_common.c"
#include "sor.c"

#define PDE_SOLVER_LIS          0
#define PDE_SOLVER_SUPERLU      1
#define PDE_SOLVER_SOR          0

static int pde_initialize(time_t t, double longitude, pde_workspace *w);
static int pde_sigma_tensor(sigma_workspace *sigma_p, pde_workspace *w);
static void pde_compute_wind(pde_workspace *w);
static int pde_scales(sigma_workspace *sigma_p, pde_workspace *w);
static int pde_coefficients(pde_workspace *w);
static int pde_discretize(pde_workspace *w);
static int pde_matrix(pde_workspace *w);
static int pde_compute_psi(gsl_spmatrix *A, gsl_vector *b, pde_workspace *w);
static int pde_check_psi(pde_workspace *w);
static int pde_calc_J(pde_workspace *w);
static int pde_main_field(pde_workspace *w);
static int pde_debug(pde_workspace *w, const char *format, ...);
static int pde_coefs(const double r, const double theta, double coef[GSL_PDE2D_COEF_TOTAL], void *params);
static int pde_bnd(const gsl_pde2d_bnd_t type, const double *v, gsl_vector *alpha, gsl_vector *beta,
                   gsl_vector *gamma, gsl_vector *delta, void *params);
static int pde_magfield(pde_workspace *w);

pde_workspace *
pde_alloc(pde_parameters *params)
{
  pde_workspace *w;
  size_t i, j;
  size_t nrt; /* nr * ntheta */

  w = calloc(1, sizeof(pde_workspace));
  if (!w)
    {
      fprintf(stderr, "pde_alloc: calloc failed: %s\n", strerror(errno));
      return 0;
    }

  w->nr = params->nr;
  w->ntheta = params->ntheta;
  w->rmin = params->rmin;
  w->rmax = params->rmax;
  w->theta_min = params->theta_min;
  w->theta_max = params->theta_max;

  w->dr = (w->rmax - w->rmin) / (w->nr - 1.0);
  w->dtheta = (w->theta_max - w->theta_min) / (w->ntheta - 1.0);

  nrt = w->nr * w->ntheta;

  w->hwm_workspace_p = hwm_alloc(params->f107_file);
  if (!w->hwm_workspace_p)
    {
      pde_free(w);
      return 0;
    }

  /*hwm_set_error_scale(4.0, w->hwm_workspace_p);*/

  w->b = gsl_vector_alloc(nrt);
  w->psi = gsl_vector_alloc(nrt);
  w->J_lat = gsl_vector_alloc(w->ntheta);
  w->J_lat_E = gsl_vector_alloc(w->ntheta);
  w->J_lat_u = gsl_vector_alloc(w->ntheta);
  w->J_r = gsl_matrix_alloc(w->nr, w->ntheta);
  w->J_theta = gsl_matrix_alloc(w->nr, w->ntheta);
  w->J_phi = gsl_matrix_alloc(w->nr, w->ntheta);
  w->E_r = gsl_matrix_alloc(w->nr, w->ntheta);
  w->E_theta = gsl_matrix_alloc(w->nr, w->ntheta);
  w->E_phi = gsl_matrix_alloc(w->nr, w->ntheta);
  if (w->b == 0 || w->psi == 0 || w->J_lat == 0 ||
      w->J_r == 0 || w->J_theta == 0 || w->J_phi == 0 ||
      w->E_r == 0 || w->E_theta == 0 || w->E_phi == 0)
    {
      pde_free(w);
      return 0;
    }

  w->merid = malloc((w->nr + 1) * sizeof(double));
  w->zonal = malloc((w->nr + 1) * sizeof(double));
  if (!w->merid || !w->zonal)
    {
      pde_free(w);
      return 0;
    }

  w->alpha = malloc(nrt * sizeof(double));
  w->beta = malloc(nrt * sizeof(double));
  w->gamma = malloc(nrt * sizeof(double));
  if (!w->alpha || !w->beta || !w->gamma)
    {
      pde_free(w);
      return 0;
    }

  w->f1 = malloc(nrt * sizeof(double));
  w->f2 = malloc(nrt * sizeof(double));
  w->f3 = malloc(nrt * sizeof(double));
  w->f4 = malloc(nrt * sizeof(double));
  w->f5 = malloc(nrt * sizeof(double));
  w->f6 = malloc(nrt * sizeof(double));
  if (!w->f1 || !w->f2 || !w->f3 || !w->f4 || !w->f5 || !w->f6)
    {
      pde_free(w);
      return 0;
    }

  for (i = 0; i < 9; ++i)
    {
      w->DC[i] = malloc(nrt * sizeof(double));
      if (!w->DC[i])
        {
          pde_free(w);
          return 0;
        }
    }

  w->zwind = malloc(nrt * sizeof(double));
  w->mwind = malloc(nrt * sizeof(double));
  w->vwind = malloc(nrt * sizeof(double));
  if (!w->zwind || !w->mwind || !w->vwind)
    {
      pde_free(w);
      return 0;
    }

  w->Br_main = malloc(nrt * sizeof(double));
  w->Bt_main = malloc(nrt * sizeof(double));
  w->Bp_main = malloc(nrt * sizeof(double));
  w->Bf_main = malloc(nrt * sizeof(double));
  if (!w->Br_main || !w->Bt_main || !w->Bp_main || !w->Bf_main)
    {
      pde_free(w);
      return 0;
    }

  w->mageq_workspace_p = mageq_alloc();
  if (!w->mageq_workspace_p)
    {
      pde_free(w);
      return 0;
    }

  w->msynth_workspace_p = msynth_igrf_read(MSYNTH_IGRF_FILE);
  if (!w->msynth_workspace_p)
    {
      pde_free(w);
      return 0;
    }

  w->nprocs = 2;

  w->S = gsl_spmatrix_alloc(nrt, nrt);
  if (!w->S)
    {
      pde_free(w);
      return 0;
    }

  w->sigma = calloc(1, w->nr * w->ntheta * sizeof(gsl_matrix *));
  if (!w->sigma)
    {
      pde_free(w);
      return 0;
    }

  for (i = 0; i < w->nr; ++i)
    {
      for (j = 0; j < w->ntheta; ++j)
        {
          w->sigma[PDE_IDX(i, j, w)] = gsl_matrix_alloc(3, 3);
          if (!w->sigma[PDE_IDX(i, j, w)])
            {
              pde_free(w);
              return 0;
            }
        }
    }

  w->sigma_workspace_p = sigma_alloc(params->f107_file, w->nr, w->ntheta,
                                     w->rmin, w->rmax,
                                     w->theta_min, w->theta_max);

  /* initialize scaling factors to 1 - they are computed later */
  w->r_s = 1.0;
  w->sigma_s = 1.0;
  w->B_s = 1.0;
  w->U_s = 1.0;
  w->E_s = 1.0;
  w->psi_s = 1.0;
  w->J_s = 1.0;

  /* store explicit grid points in array */

  w->r_grid = malloc(w->nr * sizeof(double));
  w->theta_grid = malloc(w->ntheta * sizeof(double));

  for (i = 0; i < w->nr; ++i)
    w->r_grid[i] = pde_r(i, w);

  for (i = 0; i < w->ntheta; ++i)
    w->theta_grid[i] = pde_theta(i, w);

  return w;
} /* pde_alloc() */

void
pde_free(pde_workspace *w)
{
  size_t i;

  if (w->b)
    gsl_vector_free(w->b);

  if (w->psi)
    gsl_vector_free(w->psi);

  if (w->J_lat)
    gsl_vector_free(w->J_lat);

  if (w->J_lat_E)
    gsl_vector_free(w->J_lat_E);

  if (w->J_lat_u)
    gsl_vector_free(w->J_lat_u);

  if (w->J_r)
    gsl_matrix_free(w->J_r);

  if (w->J_theta)
    gsl_matrix_free(w->J_theta);

  if (w->J_phi)
    gsl_matrix_free(w->J_phi);

  if (w->E_r)
    gsl_matrix_free(w->E_r);

  if (w->E_theta)
    gsl_matrix_free(w->E_theta);

  if (w->E_phi)
    gsl_matrix_free(w->E_phi);

  if (w->merid)
    free(w->merid);

  if (w->zonal)
    free(w->zonal);

  if (w->alpha)
    free(w->alpha);

  if (w->beta)
    free(w->beta);

  if (w->gamma)
    free(w->gamma);

  if (w->f1)
    free(w->f1);

  if (w->f2)
    free(w->f2);

  if (w->f3)
    free(w->f3);

  if (w->f4)
    free(w->f4);

  if (w->f5)
    free(w->f5);

  if (w->f6)
    free(w->f6);

  if (w->zwind)
    free(w->zwind);

  if (w->mwind)
    free(w->mwind);

  if (w->vwind)
    free(w->vwind);

  if (w->Br_main)
    free(w->Br_main);

  if (w->Bt_main)
    free(w->Bt_main);

  if (w->Bp_main)
    free(w->Bp_main);

  if (w->Bf_main)
    free(w->Bf_main);

  if (w->S)
    gsl_spmatrix_free(w->S);

  if (w->r_grid)
    free(w->r_grid);

  if (w->theta_grid)
    free(w->theta_grid);

  if (w->hwm_workspace_p)
    hwm_free(w->hwm_workspace_p);

  if (w->mageq_workspace_p)
    mageq_free(w->mageq_workspace_p);

  if (w->msynth_workspace_p)
    msynth_free(w->msynth_workspace_p);

  if (w->sigma_workspace_p)
    sigma_free(w->sigma_workspace_p);

  if (w->sigma)
    {
      size_t i, j;

      for (i = 0; i < w->nr; ++i)
        {
          for (j = 0; j < w->ntheta; ++j)
            {
              if (w->sigma[PDE_IDX(i, j, w)])
                gsl_matrix_free(w->sigma[PDE_IDX(i, j, w)]);
            }
        }

      free(w->sigma);
    }

  for (i = 0; i < 9; ++i)
    {
      if (w->DC[i])
        free(w->DC[i]);
    }

  free(w);
} /* pde_free() */

/*
pde_proc()

Inputs: t         - timestamp of equator crossing
        longitude - longitude of equator crossing
        w         - workspace

Notes:
1) On output, solutions are stored in
   w->J_lat_E = J(E, u = 0)
   w->J_lat_u = J(E = 0, u)
*/

int
pde_proc(const time_t t, const double longitude,
         pde_workspace *w)
{
  int s = 0;

  gsl_vector_set_zero(w->J_lat_E);
  gsl_vector_set_zero(w->J_lat_u);

  pde_debug(w, "pde_proc: computing conductivities on (%zu,%zu) grid...",
            w->nr, w->ntheta);
  s = sigma_calc(t, longitude, w->sigma_workspace_p);
  pde_debug(w, "done (s = %d)\n", s);
  if (s)
    return s;

  pde_debug(w, "pde_proc: initializing PDE parameters...");
  s = pde_initialize(t, longitude, w);
  pde_debug(w, "done (s = %d)\n", s);
  if (s)
    return s;

  /*
   * Record errors in the PDE solution but don't return
   * prematurely - the solution is calculated even if the iterative
   * sparse solver fails, to get an idea of what the solution looks
   * like
   */

#if 0 /*XXX*/

  pde_debug(w, "pde_proc: ---- solving for E_0 (no winds) ----\n");

  s += pde_solve(0, EEF_PHI_0, w);

  /* save J(E_0, u = 0) solution */
  gsl_vector_memcpy(w->J_lat_E, w->J_lat);

  pde_debug(w, "pde_proc: ---- solving for winds (no E_0) ----\n");

  s += pde_solve(1, 0.0, w);

  /* save J(E_0 = 0, u) solution */
  gsl_vector_memcpy(w->J_lat_u, w->J_lat);

#else

  {
    double Ephi0 = 0.1e-3;

    pde_debug(w, "pde_proc: ---- solving for FULL SOLUTION (E0 = %.1f [mV/m]) ----\n", Ephi0 * 1.0e3);
    s += pde_solve(1, Ephi0, w);

#if 0
    pde_debug(w, "pde_proc: ---- computing magnetic field grids ----\n");
    s += pde_magfield(w);
#endif
  }

#endif

  return s;
} /* pde_proc() */

/*
pde_initialize()
  Initialize PDE parameters which don't depend on E_0 or the
winds

Inputs: t         - timestamp of equator crossing
        longitude - longitude of equator crossing (radians)
        w         - pde workspace
*/

static int
pde_initialize(time_t t, double longitude, pde_workspace *w)
{
  int s = 0;
  const double r = 6371.2 + 108.0;
  const double tyr = get_year(t);

  w->t = t;
  w->longitude = longitude;

  pde_debug(w, "pde_init: computing latitude of magnetic equator...");
  w->lat_eq = mageq_calc(longitude, r, tyr, w->mageq_workspace_p);
  pde_debug(w, "done (%f degrees)\n", w->lat_eq * 180.0 / M_PI);

  pde_debug(w, "pde_init: computing EEJ angle...");
  w->eej_angle = mageq_angle(w->longitude, r, tyr, w->mageq_workspace_p);
  pde_debug(w, "done (angle = %f degrees)\n", w->eej_angle * 180.0 / M_PI);

  pde_debug(w, "pde_init: computing magnetic main field...");
  pde_main_field(w);
  pde_debug(w, "done\n");

  pde_debug(w, "pde_init: constructing sigma tensor...");
  s += pde_sigma_tensor(w->sigma_workspace_p, w);
  pde_debug(w, "done\n");
  if (s)
    return s;

  pde_debug(w, "pde_init: computing HWM winds...");
  pde_compute_wind(w);
  pde_debug(w, "done\n");

  pde_debug(w, "pde_init: computing scaling factors...");
  s += pde_scales(w->sigma_workspace_p, w);
  pde_debug(w, "done\n");

  /* now that we have the r scale we can allocate pde2d workspace */
  {
    gsl_pde2d_parameters pde2d_params;

    pde2d_params.nx = w->nr;
    pde2d_params.ny = w->ntheta;
    pde2d_params.xmin = w->rmin / w->r_s;
    pde2d_params.xmax = w->rmax / w->r_s;
    pde2d_params.ymin = w->theta_min;
    pde2d_params.ymax = w->theta_max;
    pde2d_params.coef_function = pde_coefs;
    pde2d_params.coef_params = w;
    pde2d_params.bnd_function = pde_bnd;
    pde2d_params.bnd_params = w;

    w->gsl_pde2d_workspace_p = gsl_pde2d_alloc(&pde2d_params);
  }


  return s;
} /* pde_initialize() */

/*
pde_solve()
  Solve the PDE for the given time and longitude

Inputs: compute_winds - use winds in PDE solution? (0/1)
        E_phi0        - eastward electric field in V/m
        w             - pde workspace

Notes:

1) pde_initialize() must be called prior to this function

2) height-integrated current solution is stored in w->J_lat on output
*/

int
pde_solve(int compute_winds, double E_phi0, pde_workspace *w)
{
  double min, max;
  int s = 0;
  double normb;
  gsl_spmatrix *A;
  gsl_vector *rhs;

  w->compute_winds = compute_winds;
  w->E_phi0 = E_phi0;

  /* compute coefficients of PDE for all grid points */
  pde_debug(w, "pde_solve: computing PDE coefficients...");
  s = pde_coefficients(w);
  pde_debug(w, "done (s = %d)\n", s);
  if (s)
    return s;

  pde_debug(w, "pde_solve: constructing difference equation coefficients...");
  s = pde_discretize(w);
  pde_debug(w, "done (s = %d)\n", s);
  if (s)
    return s;

  pde_debug(w, "pde_solve: constructing PDE matrix...");

  s = pde_matrix(w);

#if 0
  A = w->S;
  rhs = w->b;
#else
  A = w->gsl_pde2d_workspace_p->S;
  rhs = w->gsl_pde2d_workspace_p->rhs;
#endif

  pde_debug(w, "done (non zero elements = %zu/%zu [%.3f%%])\n",
            gsl_spmatrix_nnz(A),
            A->size1 * A->size2,
            (double) gsl_spmatrix_nnz(A) / ((double) (A->size1 * A->size2)) * 100.0);
  if (s)
    return s;

  printsp_octave(w->S, "A1");
  printsp_octave(w->gsl_pde2d_workspace_p->S, "A2");
  printv_octave(w->b, "rhs1");
  printv_octave(w->gsl_pde2d_workspace_p->rhs, "rhs2");

  gsl_spmatrix_minmax(A, &min, &max);
  pde_debug(w, "pde_solve: matrix minimum element = %.4e, maximum element = %.4e\n", min, max);
  gsl_vector_minmax(rhs, &min, &max);
  pde_debug(w, "pde_solve: rhs minimum element = %.4e, maximum element = %.4e\n", min, max);

  normb = gsl_blas_dnrm2(rhs);

  pde_debug(w, "pde_solve: computing psi solution...");
  s += pde_compute_psi(A, rhs, w);
  pde_debug(w, "done (residual norm = %.12e, relative residual norm = %.12e, condition number = %.12e)\n",
            w->residual, w->residual / normb, 1.0 / w->rcond);

#if 0
  pde_debug(w, "pde_solve: checking psi solution...");
  s += pde_check_psi(w);
  pde_debug(w, "done (s = %d)\n", s);
  if (s)
    return s;
#endif

  pde_debug(w, "pde_solve: computing eastward current...");
  s += pde_calc_J(w);
  pde_debug(w, "done (s = %d)\n", s);

  /* J_lat is height-integrated and has units of current * meters */
  pde_debug(w, "pde_solve: scaling output back to physical units...");
  gsl_vector_scale(w->J_lat, w->J_s * w->r_s);
  pde_debug(w, "done\n");

  return s;
} /* pde_solve() */

/*
pde_main_field()
  Compute main magnetic field components at all grid points and
store in w->B{r,t,p,f}
*/

static int
pde_main_field(pde_workspace *w)
{
  size_t i, j;
  double B[4];
  double Bt, Bp;
  double tyear = get_year(w->t);

  for (j = 0; j < w->ntheta; ++j)
    {
      double theta = pde_theta(j, w);

      for (i = 0; i < w->nr; ++i)
        {
          size_t k = PDE_IDX(i, j, w);
          double r = pde_r_km(i, w);

          msynth_eval(tyear, r, theta - w->lat_eq, w->longitude,
                      B, w->msynth_workspace_p);

          /* convert to T */
          B[0] *= 1.0e-9;
          B[1] *= 1.0e-9;
          B[2] *= 1.0e-9;
          B[3] *= 1.0e-9;

          w->Br_main[k] = -B[2];
          w->Bf_main[k] = B[3];

          /*
           * Perform a coordinate rotation on the (B_theta, B_phi)
           * components using the angle between the EEJ and geographic
           * eastward
           */

          Bt = -B[0];
          Bp = B[1];

          w->Bt_main[k] = Bt * cos(w->eej_angle) - Bp * sin(w->eej_angle);
          w->Bp_main[k] = Bt * sin(w->eej_angle) + Bp * cos(w->eej_angle);
        }
    }

  return GSL_SUCCESS;
}

/*
pde_sigma_tensor()
  Compute w->sigma matrices using w->s{0,1,2} previously
computed by 'sigma' program
*/

static int
pde_sigma_tensor(sigma_workspace *sigma_p, pde_workspace *w)
{
  int s = 0;
  size_t i, j;
  double b[3]; /* unit magnetic field vector */
  double s0, s1, s2;

  for (j = 0; j < w->ntheta; ++j)
    {
      double theta = pde_theta(j, w);

      for (i = 0; i < w->nr; ++i)
        {
          size_t k = PDE_IDX(i, j, w);
          gsl_matrix *s = w->sigma[k];
          double dsig;

          /*
           * get the direct, Pedersen and Hall conductivities computed
           * previously by sigma_calc()
           */
          sigma_result(i, j, &s0, &s1, &s2, sigma_p);

          if (s0 < PDE_MIN_CONDUCTIVITY)
            s0 = PDE_MIN_CONDUCTIVITY;
          if (s1 < PDE_MIN_CONDUCTIVITY)
            s1 = PDE_MIN_CONDUCTIVITY;
          if (s2 < PDE_MIN_CONDUCTIVITY)
            s2 = PDE_MIN_CONDUCTIVITY;
          
          dsig = s0 - s1;

          b[IDX_R] = w->Br_main[k] / w->Bf_main[k];
          b[IDX_THETA] = w->Bt_main[k] / w->Bf_main[k];
          b[IDX_PHI] = w->Bp_main[k] / w->Bf_main[k];

          gsl_matrix_set(s, IDX_R, IDX_R,
            dsig * b[IDX_R] * b[IDX_R] + s1);
          gsl_matrix_set(s, IDX_R, IDX_THETA,
            dsig * b[IDX_R] * b[IDX_THETA] - s2 * b[IDX_PHI]);
          gsl_matrix_set(s, IDX_R, IDX_PHI,
            dsig * b[IDX_R] * b[IDX_PHI] + s2 * b[IDX_THETA]);
          gsl_matrix_set(s, IDX_THETA, IDX_R,
            dsig * b[IDX_THETA] * b[IDX_R] + s2 * b[IDX_PHI]);
          gsl_matrix_set(s, IDX_THETA, IDX_THETA,
            dsig * b[IDX_THETA] * b[IDX_THETA] + s1);
          gsl_matrix_set(s, IDX_THETA, IDX_PHI,
            dsig * b[IDX_THETA] * b[IDX_PHI] - s2 * b[IDX_R]);
          gsl_matrix_set(s, IDX_PHI, IDX_R,
            dsig * b[IDX_PHI] * b[IDX_R] - s2 * b[IDX_THETA]);
          gsl_matrix_set(s, IDX_PHI, IDX_THETA,
            dsig * b[IDX_PHI] * b[IDX_THETA] + s2 * b[IDX_R]);
          gsl_matrix_set(s, IDX_PHI, IDX_PHI,
            dsig * b[IDX_PHI] * b[IDX_PHI] + s1);

#ifdef PDE_SIGMA_TAPER

          {
            double theta0 = pde_theta(0, w);
            double thetan = pde_theta(w->ntheta - 1, w);
            double tmp;
            double taperfac;

            if (theta - theta0 < PDE_SIGMA_TAPER_RANGE)
              {
                tmp = sin(M_PI / 2.0 * (theta - theta0) / PDE_SIGMA_TAPER_RANGE);
                taperfac = 0.1 + 0.9 * tmp * tmp;
              }
            else if (thetan - theta < PDE_SIGMA_TAPER_RANGE)
              {
                tmp = sin(M_PI / 2.0 * (thetan - theta) / PDE_SIGMA_TAPER_RANGE);
                taperfac = 0.1 + 0.9 * tmp * tmp;
              }
            else
              taperfac = 1.0;

            gsl_matrix_scale(s, taperfac);
          }

#endif /* PDE_SIGMA_TAPER */

#if 0
          printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
                 90.0 - pde_theta(j,w)*180/M_PI,
                 pde_r_km(i,w) - R_EARTH_KM,
                 w->s0[i][j],
                 w->s1[i][j],
                 w->s2[i][j],
                 b[IDX_R],
                 b[IDX_THETA],
                 b[IDX_PHI],
         /*9*/   gsl_matrix_get(s, IDX_R, IDX_R),
                 gsl_matrix_get(s, IDX_R, IDX_THETA),
                 gsl_matrix_get(s, IDX_R, IDX_PHI),
                 gsl_matrix_get(s, IDX_THETA, IDX_R),
                 gsl_matrix_get(s, IDX_THETA, IDX_THETA),
                 gsl_matrix_get(s, IDX_THETA, IDX_PHI),
                 gsl_matrix_get(s, IDX_PHI, IDX_R),
                 gsl_matrix_get(s, IDX_PHI, IDX_THETA),
                 gsl_matrix_get(s, IDX_PHI, IDX_PHI));
#endif
        } /* for (i = 0; i < w->nr; ++i) */
    } /* for (j = 0; j < w->ntheta; ++j) */
#if 0
  exit(1);
#endif

  return s;
} /* pde_sigma_tensor() */

static void
pde_compute_wind(pde_workspace *w)
{
  size_t i, j, n;
  double latitude;

  for (j = 0; j < w->ntheta; ++j)
    {
      latitude = w->lat_eq + M_PI / 2.0 - pde_theta(j, w);

      n = hwm_call(M_PI / 2.0 - latitude,
                   w->longitude,
                   w->t,
                   w->rmin * 1.0e-3 - R_EARTH_KM,
                   w->dr * 1.0e-3,
                   w->nr,
                   w->merid,
                   w->zonal,
                   w->hwm_workspace_p);
      if (n != w->nr)
        {
          fprintf(stderr, "hwm_call only computed %zu conductivities (nr = %zu)\n",
                  n, w->nr);
          exit(1);
        }

      /*
       * Store wind altitude profiles for this latitude. Since
       * zwind = u_phi and mwind = u_theta, negate merid[i] since
       * HWM outputs meridional winds as positive northward and
       * u_theta points in the southward direction.
       *
       * Also perform a coordinate transformation into wind components
       * parallel and perpendicular to the dip equator at 108km.
       * zwind contains the parallel component and mwind contains
       * the normal component.
       */
      for (i = 0; i < w->nr; ++i)
        {
          size_t k = PDE_IDX(i, j, w);
          double u_phi = w->zonal[i];
          double u_theta = -w->merid[i];

          w->zwind[k] = u_phi * cos(w->eej_angle) +
                        u_theta * sin(w->eej_angle);
          w->mwind[k] = u_theta * cos(w->eej_angle) -
                        u_phi * sin(w->eej_angle);
          w->vwind[k] = 0.0;
        }
    } /* for (j = 0; j < w->ntheta; ++j) */

#ifdef PDE_SYMMETRIC_WINDS

  for (j = 0; j < w->ntheta / 2; ++j)
    {
      for (i = 0; i < w->nr; ++i)
        {
          size_t k1 = PDE_IDX(i, j, w);
          size_t k2 = PDE_IDX(i, w->ntheta - j - 1, w);
          double tmp;
          
          tmp = 0.5 * (w->zwind[k1] + w->zwind[k2]);
          w->zwind[k1] = tmp;
          w->zwind[k2] = tmp;
        }
    }

#endif /* PDE_SYMMETRIC_WINDS */
} /* pde_compute_wind() */

/*
pde_scales()
  Compute scaling factors for non-dimensionalization of PDE equation

r_s = R_max
B_s = max(B_r,B_t,B_p)
U_s = max(U_r,U_t,U_p)
sigma_s = max(sigma)

E_s = U_s * B_s
psi_s = sigma_s * U_s * B_s * r_s^2
*/

static int
pde_scales(sigma_workspace *sigma_p, pde_workspace *w)
{
  int s = 0;
  size_t i, j;
  gsl_matrix_view m;
  double s0_max, s1_max, s2_max;

  w->r_s = pde_r_m(w->nr - 1, w);

  sigma_max(&s0_max, &s1_max, &s2_max, sigma_p);
  w->sigma_s = GSL_MAX(s0_max, GSL_MAX(s1_max, s2_max));

  /* compute B_s = max(|B|) */
  m = gsl_matrix_view_array((double *) w->Bf_main, w->nr, w->ntheta);
  w->B_s = gsl_matrix_max(&m.matrix);

  /* compute U_s = max(U_r,U_t,U_p) */
  m = gsl_matrix_view_array((double *) w->zwind, w->nr, w->ntheta);
  w->U_s = gsl_matrix_max(&m.matrix);

  m = gsl_matrix_view_array((double *) w->mwind, w->nr, w->ntheta);
  w->U_s = GSL_MAX(w->U_s, gsl_matrix_max(&m.matrix));

  m = gsl_matrix_view_array((double *) w->vwind, w->nr, w->ntheta);
  w->U_s = GSL_MAX(w->U_s, gsl_matrix_max(&m.matrix));

  for (i = 0; i < w->nr; ++i)
    {
      for (j = 0; j < w->ntheta; ++j)
        {
          size_t k = PDE_IDX(i, j, w);
          gsl_matrix *sig = w->sigma[k];

          gsl_matrix_scale(sig, 1.0 / w->sigma_s);

          w->Br_main[k] /= w->B_s;
          w->Bt_main[k] /= w->B_s;
          w->Bp_main[k] /= w->B_s;
          w->Bf_main[k] = sqrt(w->Br_main[k] * w->Br_main[k] +
                          w->Bt_main[k] * w->Bt_main[k] +
                          w->Bp_main[k] * w->Bp_main[k]);

          w->zwind[k] /= w->U_s;
          w->mwind[k] /= w->U_s;
          w->vwind[k] /= w->U_s;
        }
    }

  /*
   * these factors are chosen to make the barred equations the same
   * as the unbarred equations
   */
  w->E_s = w->U_s * w->B_s;
  w->psi_s = w->sigma_s * w->U_s * w->B_s * w->r_s * w->r_s;
  w->J_s = w->U_s * w->B_s * w->sigma_s;

  pde_debug(w, "\n\t pde_scales: r_s = %.12e m\n", w->r_s);
  pde_debug(w, "\t pde_scales: sigma_s = %.12e\n", w->sigma_s);
  pde_debug(w, "\t pde_scales: B_s = %.12e T\n", w->B_s);
  pde_debug(w, "\t pde_scales: U_s = %.12e m/s\n", w->U_s);
  pde_debug(w, "\t pde_scales: E_s = %.12e V/m\n", w->E_s);
  pde_debug(w, "\t pde_scales: psi_s = %.12e\n", w->psi_s);
  pde_debug(w, "\t pde_scales: J_s = %.12e A/m2\n", w->J_s);

  return s;
} /* pde_scales() */

/*
pde_coefficients()
  Compute coefficients of PDE: f1-f6
*/

static int
pde_coefficients(pde_workspace *w)
{
  int s = 0;
  size_t i, j;
  double dr = pde_dr(w);
  double dtheta = pde_dtheta(w);

  /* compute parameters alpha, beta, gamma */
  for (i = 0; i < w->nr; ++i)
    {
      double r = pde_r(i, w);

      for (j = 0; j < w->ntheta; ++j)
        {
          size_t k = PDE_IDX(i, j, w);
          double sint = pde_sint(j, w);
          gsl_matrix *s = w->sigma[k];
          double s_rr = gsl_matrix_get(s, IDX_R, IDX_R);
          double s_rt = gsl_matrix_get(s, IDX_R, IDX_THETA);
          double s_tr = gsl_matrix_get(s, IDX_THETA, IDX_R);
          double s_tt = gsl_matrix_get(s, IDX_THETA, IDX_THETA);
          double s_rp = gsl_matrix_get(s, IDX_R, IDX_PHI);
          double s_tp = gsl_matrix_get(s, IDX_THETA, IDX_PHI);
          double sUBr, sUBt; /* [sigma U x B]_r, [sigma U x B]_t */

          w->alpha[k] = r * sint * (s_rr * s_tt - s_tr * s_rt);

          /* E_phi may be 0 here if we are computing wind effects */
          w->beta[k] = r * sint * (s_tr * s_rp - s_rr * s_tp) * E_phi(i, j, w);
          w->gamma[k] = r * sint * (s_tt * s_rp - s_rt * s_tp) * E_phi(i, j, w);

          /* account for wind terms */
          if (w->compute_winds)
            {
              /* [sigma U x B]_r */
              sUBr = s_rr * (w->mwind[k] * w->Bp_main[k] -
                             w->zwind[k] * w->Bt_main[k]) +
                     s_rt * (w->zwind[k] * w->Br_main[k] -
                             w->vwind[k] * w->Bp_main[k]) +
                     s_rp * (w->vwind[k] * w->Bt_main[k] -
                             w->mwind[k] * w->Br_main[k]);
              /* [sigma U x B]_t */
              sUBt = s_tr * (w->mwind[k] * w->Bp_main[k] -
                             w->zwind[k] * w->Bt_main[k]) +
                     s_tt * (w->zwind[k] * w->Br_main[k] -
                             w->vwind[k] * w->Bp_main[k]) +
                     s_tp * (w->vwind[k] * w->Bt_main[k] -
                             w->mwind[k] * w->Br_main[k]);

              w->beta[k] += r * sint * (s_tr * sUBr - s_rr * sUBt);
              w->gamma[k] += r * sint * (s_tt * sUBr - s_rt * sUBt);
            }

          /* for low altitudes < 90km the conductivity could be 0 */
          if (!gsl_finite(w->alpha[k]) || !gsl_finite(w->beta[k]) ||
              !gsl_finite(w->gamma[k]))
            return GSL_FAILURE;

#if 0
          printf("%f %f %e %e %e %e %e\n",
                 pde_theta(j,w)*180/M_PI,
                 pde_r_km(i,w)-R_EARTH_KM,
                 w->alpha[k],
                 w->beta[k],
                 w->gamma[k],
                 w->Br_main[k],
                 w->Bt_main[k]);
#endif
        }
    }
#if 0
  exit(1);
#endif

  /* compute PDE coefficients f1-6 */
  for (i = 0; i < w->nr; ++i)
    {
      double r = pde_r(i, w);

      for (j = 0; j < w->ntheta; ++j)
        {
          size_t k = PDE_IDX(i, j, w);
          gsl_matrix *s = w->sigma[k];
          double s_rr = gsl_matrix_get(s, IDX_R, IDX_R);
          double s_rt = gsl_matrix_get(s, IDX_R, IDX_THETA);
          double s_tr = gsl_matrix_get(s, IDX_THETA, IDX_R);
          double s_tt = gsl_matrix_get(s, IDX_THETA, IDX_THETA);
          double dr4, dr5, dr6; /* d/dr terms in f4, f5, f6 */
          double dt4, dt5, dt6; /* d/dtheta terms in f4, f5, f6 */
          double dadr, dadt; /* d/dr alpha and d/dt alpha */

          w->f1[k] = w->alpha[k] * r * s_rr;
          w->f2[k] = w->alpha[k] * s_tt / r;
          w->f3[k] = w->alpha[k] * (s_tr + s_rt);

          if (i == 0)
            {
              /* use forward differences */

              dadr = 1.0 / dr * (w->alpha[PDE_IDX(i + 1, j, w)] -
                                 w->alpha[k]);

              dr4 = 1.0 / dr *
                    (gsl_matrix_get(w->sigma[PDE_IDX(i + 1, j, w)], IDX_THETA, IDX_R) / (r + dr) -
                     s_tr / r);

              dr5 = 1.0 / dr *
                    (gsl_matrix_get(w->sigma[PDE_IDX(i + 1, j, w)], IDX_R, IDX_R) -
                     s_rr);

              dr6 = 1.0 / dr * (w->beta[PDE_IDX(i + 1, j, w)] - w->beta[k]);
            }
          else if (i == w->nr - 1)
            {
              /* use backward differences */

              dadr = 1.0 / dr * (w->alpha[k] -
                                 w->alpha[PDE_IDX(i - 1, j, w)]);

              dr4 = 1.0 / dr *
                    (s_tr / r -
                     gsl_matrix_get(w->sigma[PDE_IDX(i - 1, j, w)], IDX_THETA, IDX_R) / (r - dr));

              dr5 = 1.0 / dr *
                    (s_rr -
                     gsl_matrix_get(w->sigma[PDE_IDX(i - 1, j, w)], IDX_R, IDX_R));

              dr6 = 1.0 / dr * (w->beta[k] -
                                w->beta[PDE_IDX(i - 1, j, w)]);
            }
          else
            {
              /* use central differences */

              dadr = 0.5 / dr * (w->alpha[PDE_IDX(i + 1, j, w)] -
                                 w->alpha[PDE_IDX(i - 1, j, w)]);

              dr4 = 0.5 / dr *
                    (gsl_matrix_get(w->sigma[PDE_IDX(i + 1, j, w)], IDX_THETA, IDX_R) / (r + dr) -
                     gsl_matrix_get(w->sigma[PDE_IDX(i - 1, j, w)], IDX_THETA, IDX_R) / (r - dr));

              dr5 = 0.5 / dr *
                    (gsl_matrix_get(w->sigma[PDE_IDX(i + 1, j, w)], IDX_R, IDX_R) -
                     gsl_matrix_get(w->sigma[PDE_IDX(i - 1, j, w)], IDX_R, IDX_R));

              dr6 = 0.5 / dr * (w->beta[PDE_IDX(i + 1, j, w)] -
                                w->beta[PDE_IDX(i - 1, j, w)]);
            }

          if (j == 0)
            {
              /* use forward differences */

              dadt = 1.0 / dtheta * (w->alpha[PDE_IDX(i, j + 1, w)] -
                                     w->alpha[k]);

              dt4 = 1.0 / dtheta *
                    (gsl_matrix_get(w->sigma[PDE_IDX(i, j + 1, w)], IDX_THETA, IDX_THETA) / (r + dr) -
                     s_tt / r);

              dt5 = 1.0 / dtheta *
                    (gsl_matrix_get(w->sigma[PDE_IDX(i, j + 1, w)], IDX_R, IDX_THETA) -
                     s_rt);

              dt6 = 1.0 / dtheta * (w->gamma[PDE_IDX(i, j + 1, w)] -
                                    w->gamma[k]);
            }
          else if (j == w->ntheta - 1)
            {
              /* use backward differences */

              dadt = 1.0 / dtheta * (w->alpha[k] -
                                     w->alpha[PDE_IDX(i, j - 1, w)]);

              dt4 = 1.0 / dtheta *
                    (s_tt / r -
                     gsl_matrix_get(w->sigma[PDE_IDX(i, j - 1, w)], IDX_THETA, IDX_THETA) / (r - dr));

              dt5 = 1.0 / dtheta *
                    (s_rt -
                     gsl_matrix_get(w->sigma[PDE_IDX(i, j - 1, w)], IDX_R, IDX_THETA));

              dt6 = 1.0 / dtheta * (w->gamma[k] -
                                    w->gamma[PDE_IDX(i, j - 1, w)]);
            }
          else
            {
              /* use central differences */

              dadt = 0.5 / dtheta * (w->alpha[PDE_IDX(i, j + 1, w)] -
                                     w->alpha[PDE_IDX(i, j - 1, w)]);

              dt4 = 0.5 / dtheta *
                    (gsl_matrix_get(w->sigma[PDE_IDX(i, j + 1, w)], IDX_THETA, IDX_THETA) / (r + dr) -
                     gsl_matrix_get(w->sigma[PDE_IDX(i, j - 1, w)], IDX_THETA, IDX_THETA) / (r - dr));

              dt5 = 0.5 / dtheta *
                    (gsl_matrix_get(w->sigma[PDE_IDX(i, j + 1, w)], IDX_R, IDX_THETA) -
                     gsl_matrix_get(w->sigma[PDE_IDX(i, j - 1, w)], IDX_R, IDX_THETA));

              dt6 = 0.5 / dtheta *
                    (w->gamma[PDE_IDX(i, j + 1, w)] -
                     w->gamma[PDE_IDX(i, j - 1, w)]);
            }

          w->f4[k] = w->alpha[k] * (s_tr / r + r * dr4 + dt4) -
                     s_tr * dadr - s_tt / r * dadt;
          w->f5[k] = w->alpha[k] * (s_rr + r * dr5 + dt5) -
                     r * s_rr * dadr - s_rt * dadt;
          w->f6[k] = w->alpha[k] * (w->beta[k] + r * dr6 + dt6) -
                     r * w->beta[k] * dadr -
                     w->gamma[k] * dadt;

          if (!gsl_finite(w->f1[k]) || !gsl_finite(w->f2[k]) ||
              !gsl_finite(w->f3[k]) || !gsl_finite(w->f4[k]) ||
              !gsl_finite(w->f5[k]) || !gsl_finite(w->f6[k]))
            return GSL_FAILURE;

#if 0
          {
            double sint = pde_sint(j, w);
            double fac = r * r * sint * sint *
                         pow(s_rt*s_tr - s_rr*s_tt, 2.0);

            printf("%f %f %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
                   theta*180/M_PI,
                   pde_r_km(i,w) - R_EARTH_KM,
                   fac,
                   w->f1[k],
                   w->f2[k],
                   w->f3[k],
                   w->f4[k],
                   w->f5[k],
                   w->f6[k],
                   dr4,
                   dr5,
                   dr6,
                   dt4,
                   dt5,
                   dt6);
          }
#endif
        } /* for (j = 0; j < w->ntheta; ++j) */
    } /* for (i = 0; i < w->nr; ++i) */
#if 0
  exit(1);
#endif

  return s;
} /* pde_coefficients() */

/*
pde_discretize()
  Construct difference matrix w->DC

DC[0][k] = coefficient of psi_{i-1,j-1}
DC[1][k] = coefficient of psi_{i,j-1}
DC[2][k] = coefficient of psi_{i+1,j-1}
DC[3][k] = coefficient of psi_{i-1,j}
DC[4][k] = coefficient of psi_{i,j}
DC[5][k] = coefficient of psi_{i+1,j}
DC[6][k] = coefficient of psi_{i-1,j+1}
DC[7][k] = coefficient of psi_{i,j+1}
DC[8][k] = coefficient of psi_{i+1,j+1}
*/

static int
pde_discretize(pde_workspace *w)
{
  int s = 0;
  size_t i, j;
  double dr = pde_dr(w);
  double dr_sq = pde_dr_sq(w);
  double dtheta = pde_dtheta(w);
  double dtheta_sq = pde_dtheta_sq(w);

  for (i = 0; i < w->nr; ++i)
    {
      for (j = 0; j < w->ntheta; ++j)
        {
          size_t k = PDE_IDX(i, j, w);
          double f1 = w->f1[k];
          double f2 = w->f2[k];
          double f3 = w->f3[k];
          double f4 = w->f4[k];
          double f5 = w->f5[k];
          double f6 = w->f6[k];

#if 0
          /* this factor will make the PDE coeffs equal to the old method */
          {
            double fac = w->alpha[k] * w->alpha[k];
            f1 /= fac;
            f2 /= fac;
            f3 /= fac;
            f4 /= fac;
            f5 /= fac;
            f6 /= fac;
            if (!gsl_finite(f1) || !gsl_finite(f2) ||
                !gsl_finite(f3) || !gsl_finite(f4) ||
                !gsl_finite(f5) || !gsl_finite(f6))
              return GSL_FAILURE;
          }
#endif

          w->DC[0][k] = f3 / 4.0 / dr / dtheta;
          w->DC[1][k] = f2 / dtheta_sq - f4 / 2.0 / dtheta;  
          w->DC[2][k] = -f3 / 4.0 / dr / dtheta;  
          w->DC[3][k] = f1 / dr_sq - f5 / 2.0 / dr;  
          w->DC[4][k] = -2.0 * f1 / dr_sq - 2.0 * f2 / dtheta_sq;  
          w->DC[5][k] = f1 / dr_sq + f5 / 2.0 / dr;  
          w->DC[6][k] = -f3 / 4.0 / dr / dtheta;  
          w->DC[7][k] = f2 / dtheta_sq + f4 / 2.0 / dtheta;  
          w->DC[8][k] = f3 / 4.0 / dr / dtheta;  

          gsl_vector_set(w->b, k, -f6);

#if 0
          printf("%f %f %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",
                 90.0 - pde_theta(j,w)*180/M_PI,
                 pde_r_km(i,w) - R_EARTH_KM,
                 w->alpha[k],
                 f1,
                 f2,
                 f3,
                 f4,
                 f5,
                 f6,
                 w->DC[0][k],
                 w->DC[1][k],
                 w->DC[2][k],
                 w->DC[3][k],
                 w->DC[4][k],
                 w->DC[5][k],
                 w->DC[6][k],
                 w->DC[7][k],
                 w->DC[8][k]);
#endif
        }
    }
#if 0
    exit(1);
#endif

  return s;
} /* pde_discretize() */

/*
pde_matrix()
  Construct PDE matrix and RHS and store in w->S and w->b
*/

static int
pde_matrix(pde_workspace *w)
{
  int s = 0;
  size_t i, j, k;

  gsl_pde2d_discretize(w->gsl_pde2d_workspace_p);

  gsl_spmatrix_set_zero(w->S);

  /*
   * add matrix elements for all r grid points, but only interior
   * theta grid points
   */
  for (i = 0; i < w->nr; ++i)
    {
      for (j = 1; j < w->ntheta - 1; ++j)
        {
          /* row index */
          k = PDE_IDX(i, j, w);

          gsl_spmatrix_set(w->S, k, PDE_IDX(i, j - 1, w), w->DC[1][k]);
          gsl_spmatrix_set(w->S, k, PDE_IDX(i, j, w), w->DC[4][k]);
          gsl_spmatrix_set(w->S, k, PDE_IDX(i, j + 1, w), w->DC[7][k]);

          /*
           * For the r boundaries, normally we would add nothing to the
           * matrix and instead add terms to the RHS vector of the form:
           *
           * DC * psi_{rmin} or DC * psi_{rmax}.
           *
           * But psi is 0 on the upper and lower boundaries so this is
           * unnecessary
           */
          if (i > 0)
            {
              gsl_spmatrix_set(w->S, k, PDE_IDX(i - 1, j - 1, w), w->DC[0][k]);
              gsl_spmatrix_set(w->S, k, PDE_IDX(i - 1, j, w), w->DC[3][k]);
              gsl_spmatrix_set(w->S, k, PDE_IDX(i - 1, j + 1, w), w->DC[6][k]);
            }

          if (i < w->nr - 1)
            {
              gsl_spmatrix_set(w->S, k, PDE_IDX(i + 1, j - 1, w), w->DC[2][k]);
              gsl_spmatrix_set(w->S, k, PDE_IDX(i + 1, j, w), w->DC[5][k]);
              gsl_spmatrix_set(w->S, k, PDE_IDX(i + 1, j + 1, w), w->DC[8][k]);
            }
        }
    }

  /* matrix elements for theta boundaries */
  for (i = 0; i < w->nr; ++i)
    {
      /* northern boundary */
      j = 0;
      k = PDE_IDX(i, j, w);

      gsl_spmatrix_set(w->S, k, PDE_IDX(i, j, w), w->DC[1][k] + w->DC[4][k]);
      gsl_spmatrix_set(w->S, k, PDE_IDX(i, j + 1, w), w->DC[7][k]);

      if (i > 0)
        {
          gsl_spmatrix_set(w->S, k, PDE_IDX(i - 1, j, w), w->DC[0][k] + w->DC[3][k]);
          gsl_spmatrix_set(w->S, k, PDE_IDX(i - 1, j + 1, w), w->DC[6][k]);
        }

      if (i < w->nr - 1)
        {
          gsl_spmatrix_set(w->S, k, PDE_IDX(i + 1, j, w), w->DC[2][k] + w->DC[5][k]);
          gsl_spmatrix_set(w->S, k, PDE_IDX(i + 1, j + 1, w), w->DC[8][k]);
        }

      /* southern boundary */
      j = w->ntheta - 1;
      k = PDE_IDX(i, j, w);

      gsl_spmatrix_set(w->S, k, PDE_IDX(i, j - 1, w), w->DC[1][k]);
      gsl_spmatrix_set(w->S, k, PDE_IDX(i, j, w), w->DC[4][k] + w->DC[7][k]);

      if (i > 0)
        {
          gsl_spmatrix_set(w->S, k, PDE_IDX(i - 1, j - 1, w), w->DC[0][k]);
          gsl_spmatrix_set(w->S, k, PDE_IDX(i - 1, j, w), w->DC[3][k] + w->DC[6][k]);
        }

      if (i < w->nr - 1)
        {
          gsl_spmatrix_set(w->S, k, PDE_IDX(i + 1, j - 1, w), w->DC[2][k]);
          gsl_spmatrix_set(w->S, k, PDE_IDX(i + 1, j, w), w->DC[5][k] + w->DC[8][k]);
        }
    }

  return s;
} /* pde_matrix() */

/*
pde_compute_psi()
  Solve the equation A psi = b for psi

Notes: on output, w->residual is set to the residual ||A psi - b||
*/

static int
pde_compute_psi(gsl_spmatrix *A, gsl_vector *b, pde_workspace *w)
{
  int s = 0;
  double bscale;

  /*
   * the RHS tends to have a much smaller order of magnitude than the
   * matrix so scale it
   */
  bscale = gsl_vector_max(b);
  gsl_vector_scale(b, 1.0 / bscale);

#if 1

  /*
   * Sparse solvers which seem to work on this matrix equation:
   * preconditioned iterative: lis
   * direct: pastix, mumps
   *
   * though pastix has some strange jumps in the solution for ntheta=200
   * lis seems to work well with ILU precond and bicgstab/gmres method
   *
   * does NOT work: superlu
   */
#if PDE_SOLVER_LIS

  {
    lis_workspace *lis_p = lis_alloc(A->size1, A->size2);
    const double tol = 1.0e-6;

    s = lis_proc(A, b->data, tol, w->psi->data, lis_p);
    w->residual = lis_p->rnorm;

    mylis_free(lis_p);

    printv_octave(w->psi, "psi_lis");

    if (s)
      fprintf(stderr, "pde_compute_psi: lis_proc failed: s = %d\n", s);
  }

#elif PDE_SOLVER_SUPERLU

  {
    slu_workspace *superlu_p = slu_alloc(A->size1, A->size2, 1);
    gsl_spmatrix *C = gsl_spmatrix_ccs(A);

    s = slu_proc(C, b->data, w->psi->data, superlu_p);
    w->residual = superlu_p->rnorm;
    w->rcond = superlu_p->rcond;

    slu_free(superlu_p);
    gsl_spmatrix_free(C);

    printv_octave(w->psi, "psi_superlu");

    if (s)
      fprintf(stderr, "pde_compute_psi: slu_proc failed: s = %d\n", s);
  }

#elif PDE_SOLVER_SOR

  {
    const size_t max_iter = 500;
    const double tol = 1.0e-6;
    gsl_spmatrix *C = gsl_spmatrix_crs(A);
    gsl_vector *r = gsl_vector_alloc(w->psi->size); /* residual vector */
    size_t iter = 0;
    int status;
    double omega = 1.5;
    double residual;

    gsl_vector_set_zero(w->psi);

    do
      {
        status = pde_sor(omega, tol, C, b, w->psi);

        /* print out residual norm ||A*u - f|| */
        gsl_vector_memcpy(r, b);
        gsl_spblas_dgemv(CblasNoTrans, -1.0, C, w->psi, 1.0, r);
        residual = gsl_blas_dnrm2(r);
        fprintf(stderr, "iter %zu residual = %.12e\n", iter, residual);

        if (status == GSL_SUCCESS)
          fprintf(stderr, "Converged\n");
      }
    while (status == GSL_CONTINUE && ++iter < max_iter);

    {
      FILE *fp = fopen("A.txt", "w");
      gsl_spmatrix_fprintf(fp, A, "%.12e");
      fclose(fp);
      exit(1);
    }

    gsl_spmatrix_free(C);
    gsl_vector_free(r);
  }

#endif

#else

  {
    const double tol = 1.0e-3;
    const size_t max_iter = 100;
    const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;
    gsl_splinalg_itersolve *work =
      gsl_splinalg_itersolve_alloc(T, w->psi->size, 40);
    size_t iter = 0;
    int status;
    double normb = gsl_blas_dnrm2(w->b);

    /* initial guess psi = 0 */
    gsl_vector_set_zero(w->psi);

    do
      {
        status = gsl_splinalg_itersolve_iterate(w->S, w->b, tol,
                                                w->psi, work);

        w->residual = gsl_splinalg_itersolve_normr(work);
        fprintf(stderr, "iter %zu residual = %.12e relerr = %.12e\n",
                iter, w->residual, w->residual / normb);

        if (status == GSL_SUCCESS)
          fprintf(stderr, "Converged in %zu iterations\n", iter + 1);
      }
    while (status == GSL_CONTINUE && ++iter < max_iter);

    printv_octave(w->psi, "psi_gsl");
  }

#endif

  gsl_vector_scale(w->psi, bscale);

  return s;
} /* pde_compute_psi() */

/*
pde_check_psi()
  Check psi solution
*/

static int
pde_check_psi(pde_workspace *w)
{
  size_t i, j;
  double psi_rr,
         psi_tt,
         psi_rt,
         psi_t,
         psi_r;
  double lhs;
  int s = GSL_SUCCESS;
  double dr = pde_dr(w);
  double dr_sq = pde_dr_sq(w);
  double dtheta = pde_dtheta(w);
  double dtheta_sq = pde_dtheta_sq(w);

  for (i = 0; i < w->nr; ++i)
    {
      for (j = 0; j < w->ntheta; ++j)
        {
          size_t k = PDE_IDX(i, j, w);

          /* check boundary conditions */
          if (i == 0 || i == w->nr - 1)
            {
#if 0
              /* psi_{ij} = 0 */
              lhs = PSI_GET(w->psi, i, j);
              if (fabs(lhs) > 1.0e0 * w->nr * w->ntheta * GSL_DBL_EPSILON)
                {
                  fprintf(stderr,
                          "pde_check_psi: BC: (i=%zu,j=%zu) lhs = %.12e\n",
                          i, j, lhs);
                  s = GSL_FAILURE;
                }
#endif

              continue;
            }

          if (j == 0)
            {
#if 0
              lhs = PSI_GET(w->psi, i, j + 1) -
                    PSI_GET(w->psi, i, j);
              if (fabs(lhs) > 1.0e0 * w->nr * w->ntheta * GSL_DBL_EPSILON)
                {
                  fprintf(stderr,
                          "pde_check_psi: BC: (i=%zu,j=%zu) lhs = %.12e\n",
                          i, j, lhs);
                  s = GSL_FAILURE;
                }
#endif

              continue;
            }

          if (j == w->ntheta - 1)
            {
#if 0
              lhs = PSI_GET(w->psi, i, j) -
                    PSI_GET(w->psi, i, j - 1);
              if (fabs(lhs) > 1.0e0 * w->nr * w->ntheta * GSL_DBL_EPSILON)
                {
                  fprintf(stderr,
                          "pde_check_psi: BC: (i=%zu,j=%zu) lhs = %.12e\n",
                          i, j, lhs);
                  s = GSL_FAILURE;
                }
#endif

              continue;
            }

          /* its an interior point */

          psi_rr = 1.0 / dr_sq *
                   (PSI_GET(w->psi, i + 1, j, w) -
                    2.0 * PSI_GET(w->psi, i, j, w) +
                    PSI_GET(w->psi, i - 1, j, w));
          psi_tt = 1.0 / dtheta_sq *
                   (PSI_GET(w->psi, i, j + 1, w) -
                    2.0 * PSI_GET(w->psi, i, j, w) +
                    PSI_GET(w->psi, i, j - 1, w));
          psi_rt = 0.25 / dr / dtheta *
                   (PSI_GET(w->psi, i + 1, j + 1, w) -
                    PSI_GET(w->psi, i - 1, j + 1, w) -
                    PSI_GET(w->psi, i + 1, j - 1, w) +
                    PSI_GET(w->psi, i - 1, j - 1, w));
          psi_t = 0.5 / dtheta *
                  (PSI_GET(w->psi, i, j + 1, w) -
                   PSI_GET(w->psi, i, j - 1, w));
          psi_r = 0.5 / dr *
                  (PSI_GET(w->psi, i + 1, j, w) -
                   PSI_GET(w->psi, i - 1, j, w));

          lhs = w->f1[k] * psi_rr +
                w->f2[k] * psi_tt +
                w->f3[k] * psi_rt +
                w->f4[k] * psi_t +
                w->f5[k] * psi_r +
                w->f6[k];

          if (fabs(lhs) > 1.0e4 * w->nr * w->ntheta * GSL_DBL_EPSILON)
            {
              fprintf(stderr,
                      "pde_check_psi: (i=%zu,j=%zu) lhs = %.12e\n",
                      i, j, lhs);
              s = GSL_FAILURE;
            }
        }
    }

  return s;
} /* pde_check_psi() */

/*
pde_calc_J()
  Compute current density components and store in w->J_{r,theta,phi}. Also
compute height integrated eastward current and store in w->J_lat

J_{r,theta,phi} and J_lat will have dimensionless units
*/

static int
pde_calc_J(pde_workspace *w)
{
  size_t i, j;
  double Er, Et, Ep; /* electric field components */
  double psi_r,  /* @/@r psi */
         psi_t;  /* @/@theta psi */
  double J_r, J_theta, J_phi;
  double dr = pde_dr(w);
  double dtheta = pde_dtheta(w);

  for (i = 0; i < w->nr; ++i)
    {
      double r = pde_r(i, w);

      for (j = 0; j < w->ntheta; ++j)
        {
          size_t k = PDE_IDX(i, j, w);
          double theta = pde_theta(j, w);
          double sint = sin(theta);
          gsl_matrix *s = w->sigma[k];
          double s_rr = gsl_matrix_get(s, IDX_R, IDX_R);
          double s_rt = gsl_matrix_get(s, IDX_R, IDX_THETA);
          double s_tr = gsl_matrix_get(s, IDX_THETA, IDX_R);
          double s_tt = gsl_matrix_get(s, IDX_THETA, IDX_THETA);
          double s_pr = gsl_matrix_get(s, IDX_PHI, IDX_R);
          double s_pt = gsl_matrix_get(s, IDX_PHI, IDX_THETA);
          double s_pp = gsl_matrix_get(s, IDX_PHI, IDX_PHI);

          if (i == 0)
            {
              psi_r = 1.0 / dr *
                      (PSI_GET(w->psi, i + 1, j, w) -
                       PSI_GET(w->psi, i, j, w));
            }
          else if (i == w->nr - 1)
            {
              psi_r = 1.0 / dr *
                      (PSI_GET(w->psi, i, j, w) -
                       PSI_GET(w->psi, i - 1, j, w));
            }
          else
            {
              psi_r = 0.5 / dr *
                      (PSI_GET(w->psi, i + 1, j, w) -
                       PSI_GET(w->psi, i - 1, j, w));
            }

          if (j == 0)
            {
              psi_t = 1.0 / dtheta *
                      (PSI_GET(w->psi, i, j + 1, w) -
                       PSI_GET(w->psi, i, j, w));
            }
          else if (j == w->ntheta - 1)
            {
              psi_t = 1.0 / dtheta *
                      (PSI_GET(w->psi, i, j, w) -
                       PSI_GET(w->psi, i, j - 1, w));
            }
          else
            {
              psi_t = 0.5 / dtheta *
                      (PSI_GET(w->psi, i, j + 1, w) -
                       PSI_GET(w->psi, i, j - 1, w));
            }

          Er = (-s_tt / r * psi_t - s_rt * psi_r - w->gamma[k]) /
               w->alpha[k];

          Et = (s_tr / r * psi_t + s_rr * psi_r + w->beta[k]) /
               w->alpha[k];

          Ep = E_phi(i, j, w);

          J_r = -psi_t / (r * r * sint);
          J_theta = psi_r / (r * sint);
          J_phi = s_pr * Er + s_pt * Et + s_pp * Ep;

          if (w->compute_winds)
            {
              J_phi += s_pr * (w->mwind[k] * w->Bp_main[k] -
                               w->zwind[k] * w->Bt_main[k]) +
                       s_pt * (w->zwind[k] * w->Br_main[k] -
                               w->vwind[k] * w->Bp_main[k]) +
                       s_pp * (w->vwind[k] * w->Bt_main[k] -
                               w->mwind[k] * w->Br_main[k]);
            }

          /* store currents */
          gsl_matrix_set(w->J_r, i, j, J_r);
          gsl_matrix_set(w->J_theta, i, j, J_theta);
          gsl_matrix_set(w->J_phi, i, j, J_phi);

          /* store electric field */
          gsl_matrix_set(w->E_r, i, j, Er);
          gsl_matrix_set(w->E_theta, i, j, Et);
          gsl_matrix_set(w->E_phi, i, j, Ep);

#if 0
          printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
                 pde_theta(j,w)*180/M_PI,
                 pde_r_km(i,w)-R_EARTH_KM,
                 s_pr,
                 s_pt,
                 s_pp,
                 s_rr,
                 s_tt,
                 s_rt,
                 s_tr,
      /*10*/     w->alpha[i][j],
                 w->beta[i][j],
                 w->gamma[i][j],
                 w->delta[i][j],
                 PSI_GET(w->psi, i, j, w),
      /*15*/     J_phi,
                 Er,
                 Et,
                 psi_r,
                 psi_t);
#endif
        }
#if 0
      printf("\n");
#endif
    }
#if 0
  exit(1);
#endif

  /* compute height integrated current density */

  dr = pde_dr(w);

  for (j = 0; j < w->ntheta; ++j)
    {
      double J_sum = 0.0;

      for (i = 0; i < w->nr; ++i)
        J_sum += gsl_matrix_get(w->J_phi, i, j) * dr;

      gsl_vector_set(w->J_lat, j, J_sum);

#if 0
      printf("%e %e\n",
             pde_theta(j,w)*180/M_PI,
             J_sum);
#endif
    }
#if 0
  exit(1);
#endif

  return GSL_SUCCESS;
} /* pde_calc_J() */

static int
pde_debug(pde_workspace *w, const char *format, ...)
{
  int s = 0;
  va_list args;

  if (w->myid != 0)
    return s;

  va_start(args, format);

  vfprintf(stderr, format, args);

  va_end(args);

  fflush(stderr);

  return s;
} /* pde_debug() */

/*
pde_coefs()
  Return coefficients of PDE for a given (r,theta)
*/

static int
pde_coefs(const double r, const double theta, double coef[GSL_PDE2D_COEF_TOTAL], void *params)
{
  pde_workspace *w = (pde_workspace *) params;
  size_t i = bsearch_double(w->r_grid, r * w->r_s, 0, w->nr - 1);
  const size_t j = bsearch_double(w->theta_grid, theta, 0, w->ntheta - 1);
  size_t k;

  assert(theta == w->theta_grid[j]);

  if (fabs(r - w->r_grid[i] / w->r_s) > 1.0e-12)
    ++i;

  assert(fabs(r - w->r_grid[i] / w->r_s) < 1.0e-12);

  k = PDE_IDX(i, j, w);

  coef[GSL_PDE2D_COEF_FXX] = w->f1[k];
  coef[GSL_PDE2D_COEF_FXY] = w->f3[k];
  coef[GSL_PDE2D_COEF_FYY] = w->f2[k];
  coef[GSL_PDE2D_COEF_FX] = w->f5[k];
  coef[GSL_PDE2D_COEF_FY] = w->f4[k];
  coef[GSL_PDE2D_COEF_F] = 0.0;
  coef[GSL_PDE2D_COEF_RHS] = -w->f6[k];

  return GSL_SUCCESS;
}

static int
pde_bnd(const gsl_pde2d_bnd_t type, const double *v, gsl_vector *alpha, gsl_vector *beta,
        gsl_vector *gamma, gsl_vector *delta, void *params)
{
  const size_t n = alpha->size;
  size_t i;

  for (i = 0; i < n; ++i)
    {
      if (type == GSL_PDE2D_BND_X_LOWER)
        {
          /* lower radial boundary: psi = 0 */
          gsl_vector_set(alpha, i, 0.0);
          gsl_vector_set(beta, i, 0.0);
          gsl_vector_set(gamma, i, 1.0);
          gsl_vector_set(delta, i, 0.0);
        }
      else if (type == GSL_PDE2D_BND_X_UPPER)
        {
          /* upper radial boundary: psi = 0 */
          gsl_vector_set(alpha, i, 0.0);
          gsl_vector_set(beta, i, 0.0);
          gsl_vector_set(gamma, i, 1.0);
          gsl_vector_set(delta, i, 0.0);
        }
      else if (type == GSL_PDE2D_BND_Y_LOWER)
        {
          /* lower theta boundary: psi_theta = 0 */
          gsl_vector_set(alpha, i, 0.0);
          gsl_vector_set(beta, i, 1.0);
          gsl_vector_set(gamma, i, 0.0);
          gsl_vector_set(delta, i, 0.0);
        }
      else if (type == GSL_PDE2D_BND_Y_UPPER)
        {
          /* upper theta boundary: psi_theta = 0 */
          gsl_vector_set(alpha, i, 0.0);
          gsl_vector_set(beta, i, 1.0);
          gsl_vector_set(gamma, i, 0.0);
          gsl_vector_set(delta, i, 0.0);
        }
    }

  return GSL_SUCCESS;
}

/*
pde_magfield()
*/

static int
pde_magfield(pde_workspace *w)
{
  int s = 0;
  const size_t nr = w->nr;
  const size_t ntheta = 2 * w->ntheta;
  const size_t nphi = 72;
  magfield_workspace *magfield_p;
  magfield_params params;
  size_t i, j, k;

  params.lmin = 1;
  params.lmax = 45;
  params.mmax = 30;
  params.nr = nr;
  params.ntheta = ntheta;
  params.nphi = nphi;
  params.rmin = w->rmin - 0.001;
  params.rmax = w->rmax + 0.001;
  params.R = R_EARTH_KM * 1.0e3;
  params.grid_type = MAGFIELD_GAUSS;

  fprintf(stderr, "\t allocating magfield workspace...");
  magfield_p = magfield_alloc(&params);
  fprintf(stderr, "done\n");

  /* fill magfield J grid */
  for (i = 0; i < nr; ++i)
    {
      for (j = 0; j < ntheta; ++j)
        {
          double theta = magfield_p->theta[j];
          double Jr, Jt, Jp;

          if (theta < w->theta_min || theta > w->theta_max)
            {
              Jr = 0.0;
              Jt = 0.0;
              Jp = 0.0;
            }
          else
            {
              /*size_t idx = bsearch_double(w->theta_grid, theta, 0, w->ntheta - 1);*/
              size_t idx = pde_thidx(theta, w);

              assert(w->theta_grid[idx] <= theta && theta < w->theta_grid[idx + 1]);

              Jr = interp1d(w->theta_grid[idx], w->theta_grid[idx + 1],
                            gsl_matrix_get(w->J_r, i, idx), gsl_matrix_get(w->J_r, i, idx + 1), theta);
              Jt = interp1d(w->theta_grid[idx], w->theta_grid[idx + 1],
                            gsl_matrix_get(w->J_theta, i, idx), gsl_matrix_get(w->J_theta, i, idx + 1), theta);
              Jp = interp1d(w->theta_grid[idx], w->theta_grid[idx + 1],
                            gsl_matrix_get(w->J_phi, i, idx), gsl_matrix_get(w->J_phi, i, idx + 1), theta);
            }

          for (k = 0; k < nphi; ++k)
            {
              magfield_current_set(MAG_IDX_R, i, j, k, Jr * w->J_s, magfield_p);
              magfield_current_set(MAG_IDX_THETA, i, j, k, Jt * w->J_s, magfield_p);
              magfield_current_set(MAG_IDX_PHI, i, j, k, Jp * w->J_s, magfield_p);
            }
        }
    }

  /* perform SH decomposition */
  fprintf(stderr, "\t performing SH decomposition...");
  magfield_decomp(magfield_p);
  fprintf(stderr, "done\n");

#if 0
  {
    double r;

    for (r = R_EARTH_M; r < (R_EARTH_KM + 750.0) * 1.0e3; r += 0.5*1.0e3)
      {
        complex double p10 = magfield_eval_plmr(r, 1, 0, magfield_p);
        complex double p20 = magfield_eval_plmr(r, 2, 0, magfield_p);
        complex double p30 = magfield_eval_plmr(r, 3, 0, magfield_p);
        complex double qt10 = magfield_eval_qtlmr(r, 1, 0, magfield_p);
        complex double qt20 = magfield_eval_qtlmr(r, 2, 0, magfield_p);
        complex double qt30 = magfield_eval_qtlmr(r, 3, 0, magfield_p);

        printf("%f %.12e %.12e %.12e %.12e %.12e %.12e\n",
               (r - R_EARTH_M) * 1.0e-3,
               creal(qt10),
               creal(qt20),
               creal(qt30),
               creal(p10),
               creal(p20),
               creal(p30));
      }

    exit(1);
  }
#endif

#if 1
  for (i = 0; i < w->nr; ++i)
    {
      double r = w->r_grid[i];

      for (j = 0; j < w->ntheta; ++j)
        {
          double theta = w->theta_grid[j];
          double J_magfield[3], B_magfield[4];

          magfield_eval_J(r, theta, 0.0, J_magfield, magfield_p);
          magfield_eval_B(r, theta, 0.0, B_magfield, magfield_p);

          printf("%f %f %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",
                 r * 1.0e-3 - R_EARTH_KM,
                 90.0 - theta*180.0/M_PI,
                 J_magfield[0],
                 gsl_matrix_get(w->J_r, i, j) * w->J_s,
                 J_magfield[1],
                 gsl_matrix_get(w->J_theta, i, j) * w->J_s,
                 J_magfield[2],
                 gsl_matrix_get(w->J_phi, i, j) * w->J_s,
                 B_magfield[0] * 1.0e9,
                 B_magfield[1] * 1.0e9,
                 B_magfield[2] * 1.0e9,
                 B_magfield[3] * 1.0e9);
        }
      printf("\n");
    }
#endif

#if 0
  {
    double theta;

    for (theta = 0.1; theta <= M_PI - 0.1; theta += 1.0e-3)
      {
        double B[4];

        magfield_eval_B((R_EARTH_KM + 450.0) * 1.0e3, theta, 0.0, B, magfield_p);

        printf("%f %.12e %.12e %.12e %.12e\n",
               90.0 - theta*180.0/M_PI,
               B[0] * 1.0e9,
               B[1] * 1.0e9,
               B[2] * 1.0e9,
               B[3] * 1.0e9);
      }
  }
#endif

  exit(1);

  magfield_free(magfield_p);

  return s;
}
