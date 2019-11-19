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

#include <mainlib/ml_common.h>
#include <mainlib/ml_bsearch.h>
#include <mainlib/ml_interp.h>
#include <mainlib/ml_oct.h>
#include <mainlib/ml_mageq.h>

#include <magfield/magfield.h>

#include "superlu.h"

#include "pde.h"
#include "sigma.h"

#include "pde_common.c"
#include "sor.c"

#define PDE_SOLVER_LIS          0
#define PDE_SOLVER_SUPERLU      1
#define PDE_SOLVER_SOR          0

static int pde_initialize(time_t t, double longitude, pde_workspace *w);
static int pde_solve(const int compute_winds, const double E_phi0, pde_workspace *w);
static int pde_sigma_tensor(sigma_workspace *sigma_p, pde_workspace *w);
static void pde_compute_wind(pde_workspace *w);
static int pde_scales(sigma_workspace *sigma_p, pde_workspace *w);
static int pde_coefficients(pde_workspace *w);
static int pde_rhs_g1(gsl_vector *g1, pde_workspace *w);
static int pde_rhs_g2(const gsl_vector *W_r, const gsl_vector *W_t, gsl_vector *g2, pde_workspace *w);
static int pde_compute_solution(const gsl_spmatrix *A, const gsl_matrix *B, gsl_matrix *X, pde_workspace *w);
static int pde_calc_J(const gsl_vector * psi, pde_workspace *w);
static int pde_calc_W(pde_workspace *w);
static int pde_calc_Jr(const gsl_vector * psi, gsl_vector * Jr, pde_workspace *w);
static int pde_calc_Jt(const gsl_vector * psi, gsl_vector * Jt, pde_workspace *w);
static int pde_calc_Jp(pde_workspace *w);
static int pde_calc_E(pde_workspace *w);
static int pde_main_field(pde_workspace *w);
static int pde_debug(pde_workspace *w, const char *format, ...);
static int pde2d_coefs(const double r, const double theta, double coef[GSL_PDE2D_COEF_TOTAL], void *params);
static int pde2d_rhs(const double r, const double theta, double *rhsval, void *params);
static int pde2d_bnd(const gsl_pde2d_bnd_t type, const double *v, gsl_vector *alpha, gsl_vector *beta,
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
  w->R = params->R;

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

  w->J_lat = gsl_vector_alloc(w->ntheta);
  w->J_lat_E = gsl_vector_alloc(w->ntheta);
  w->J_lat_u = gsl_vector_alloc(w->ntheta);
  w->J_r = gsl_matrix_alloc(w->nr, w->ntheta);
  w->J_theta = gsl_matrix_alloc(w->nr, w->ntheta);
  w->J_phi = gsl_matrix_alloc(w->nr, w->ntheta);
  w->E_r = gsl_matrix_alloc(w->nr, w->ntheta);
  w->E_theta = gsl_matrix_alloc(w->nr, w->ntheta);
  w->E_phi = gsl_matrix_alloc(w->nr, w->ntheta);
  if (w->J_lat == 0 ||
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
  if (!w->f1 || !w->f2 || !w->f3 || !w->f4 || !w->f5)
    {
      pde_free(w);
      return 0;
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
  w->nrhs = 2;

  w->S = gsl_spmatrix_alloc(nrt, nrt);
  w->B = gsl_matrix_alloc(nrt, w->nrhs);
  w->G = gsl_matrix_alloc(nrt, w->nrhs);
  w->PSI = gsl_matrix_alloc(nrt, w->nrhs);
  w->JR = gsl_matrix_alloc(nrt, w->nrhs);
  w->JTHETA = gsl_matrix_alloc(nrt, w->nrhs);
  w->JPHI = gsl_matrix_alloc(nrt, w->nrhs);
  w->JPHI_HI = gsl_matrix_alloc(w->ntheta, w->nrhs);
  w->WR = gsl_matrix_alloc(nrt, w->nrhs);
  w->WTHETA = gsl_matrix_alloc(nrt, w->nrhs);
  w->WPHI = gsl_matrix_alloc(nrt, w->nrhs);
  w->ER = gsl_matrix_alloc(nrt, w->nrhs);
  w->ETHETA = gsl_matrix_alloc(nrt, w->nrhs);
  w->EPHI = gsl_matrix_alloc(nrt, w->nrhs);
  if (!w->S || !w->B || !w->G || !w->PSI || !w->JR || !w->JTHETA || !w->JPHI ||
      !w->JPHI_HI || !w->WR || !w->WTHETA || !w->WPHI || !w->ER || !w->ETHETA || !w->EPHI)
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

  w->superlu_workspace_p = slu_alloc(nrt, nrt, w->nprocs, w->nrhs);

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

  if (w->B)
    gsl_matrix_free(w->B);

  if (w->G)
    gsl_matrix_free(w->G);

  if (w->PSI)
    gsl_matrix_free(w->PSI);

  if (w->JR)
    gsl_matrix_free(w->JR);

  if (w->JTHETA)
    gsl_matrix_free(w->JTHETA);

  if (w->JPHI)
    gsl_matrix_free(w->JPHI);

  if (w->JPHI_HI)
    gsl_matrix_free(w->JPHI_HI);

  if (w->WR)
    gsl_matrix_free(w->WR);

  if (w->WTHETA)
    gsl_matrix_free(w->WTHETA);

  if (w->WPHI)
    gsl_matrix_free(w->WPHI);

  if (w->ER)
    gsl_matrix_free(w->ER);

  if (w->ETHETA)
    gsl_matrix_free(w->ETHETA);

  if (w->EPHI)
    gsl_matrix_free(w->EPHI);

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

  if (w->superlu_workspace_p)
    slu_free(w->superlu_workspace_p);

  free(w);
}

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
    double Ephi0 = 1.0e-3; /* V/m */

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
    pde2d_params.coef_function = pde2d_coefs;
    pde2d_params.coef_params = w;
    pde2d_params.bnd_function = pde2d_bnd;
    pde2d_params.bnd_params = w;
    pde2d_params.rhs_function = pde2d_rhs;
    pde2d_params.rhs_params = w;

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

static int
pde_solve(const int compute_winds, const double E_phi0, pde_workspace *w)
{
  double min, max;
  int s = 0;
  size_t i;

  w->compute_winds = compute_winds;
  w->E_phi0 = E_phi0;

  /* compute coefficients of PDE for all grid points */
  pde_debug(w, "pde_solve: computing PDE coefficients...");
  s = pde_coefficients(w);
  pde_debug(w, "done (s = %d)\n", s);
  if (s)
    return s;

  pde_debug(w, "pde_solve: computing W matrices...");
  s = pde_calc_W(w);
  pde_debug(w, "done (s = %d)\n", s);
  if (s)
    return s;


  {
    gsl_vector_view g1 = gsl_matrix_column(w->G, 0);
    gsl_vector_view g2 = gsl_matrix_column(w->G, 1);
    gsl_vector_const_view Wr = gsl_matrix_const_column(w->WR, 1);
    gsl_vector_const_view Wt = gsl_matrix_const_column(w->WTHETA, 1);

    pde_debug(w, "pde_solve: computing PDE right hand sides...");

    /* compute g_1 */
    s = pde_rhs_g1(&g1.vector, w);

    /* compute g_2(W_0) */
    s += pde_rhs_g2(&Wr.vector, &Wt.vector, &g2.vector, w);

    pde_debug(w, "done (s = %d)\n", s);
    if (s)
      return s;
  }

  pde_debug(w, "pde_solve: initializing PDE...");
  s = gsl_pde2d_init(w->gsl_pde2d_workspace_p);
  pde_debug(w, "done\n");

  pde_debug(w, "pde_solve: discretizing PDE...");

  s = gsl_pde2d_discretize(w->S, w->gsl_pde2d_workspace_p);

  pde_debug(w, "done (non zero elements = %zu/%zu [%.3f%%])\n",
            gsl_spmatrix_nnz(w->S),
            w->S->size1 * w->S->size2,
            (double) gsl_spmatrix_nnz(w->S) / ((double) (w->S->size1 * w->S->size2)) * 100.0);
  if (s)
    return s;

  pde_debug(w, "pde_solve: computing rhs of PDE...");

  {
    gsl_vector_view b1 = gsl_matrix_column(w->B, 0);
    gsl_vector_view b2 = gsl_matrix_column(w->B, 1);

    w->flags = PDE_FLG_RHS_G1;
    s = gsl_pde2d_rhs(&b1.vector, w->gsl_pde2d_workspace_p);
    if (s)
      return s;

    w->flags = PDE_FLG_RHS_G2;
    s = gsl_pde2d_rhs(&b2.vector, w->gsl_pde2d_workspace_p);
    if (s)
      return s;
  }

  pde_debug(w, "done (s = %d)\n", s);

#if 0
  printsp_octave(w->S, "A");
#endif

  gsl_spmatrix_minmax(w->S, &min, &max);
  pde_debug(w, "pde_solve: matrix minimum element = %.4e, maximum element = %.4e\n", min, max);

  gsl_matrix_minmax(w->B, &min, &max);
  pde_debug(w, "pde_solve: rhs minimum element = %.4e, maximum element = %.4e\n", min, max);

  pde_debug(w, "pde_solve: computing PDE solution...\n");

  {
    gsl_vector_const_view b1 = gsl_matrix_const_column(w->B, 0);
    gsl_vector_const_view b2 = gsl_matrix_const_column(w->B, 1);

    s += pde_compute_solution(w->S, w->B, w->PSI, w);

#if 0
    /* psi = E_phi0 * psi1 + psi2 */
    gsl_vector_memcpy(w->psi, &psi1.vector);
    gsl_vector_scale(w->psi, E_phi0 / w->E_s);
    gsl_vector_add(w->psi, &psi2.vector);
#endif

    pde_debug(w, "\t g1 residual = %.12e (relative: %.12e)\n",
              w->superlu_workspace_p->rnorm[0],
              w->superlu_workspace_p->rnorm[0] / gsl_blas_dnrm2(&b1.vector));

    pde_debug(w, "\t g2 residual = %.12e (relative: %.12e)\n",
              w->superlu_workspace_p->rnorm[1],
              w->superlu_workspace_p->rnorm[1] / gsl_blas_dnrm2(&b2.vector));

    pde_debug(w, "\t condition number = %.12e\n", 1.0 / w->rcond);
  }

  pde_debug(w, "done\n");

  pde_debug(w, "pde_solve: computing J_r and J_theta...");

  for (i = 0; i < w->nrhs; ++i)
    {
      gsl_vector_const_view psi = gsl_matrix_const_column(w->PSI, i);
      gsl_vector_view Jr = gsl_matrix_column(w->JR, i);
      gsl_vector_view Jt = gsl_matrix_column(w->JTHETA, i);

      s += pde_calc_Jr(&psi.vector, &Jr.vector, w);
      s += pde_calc_Jt(&psi.vector, &Jt.vector, w);
    }

  pde_debug(w, "done (s = %d)\n", s);

  pde_debug(w, "pde_solve: computing E matrices...");
  s += pde_calc_E(w);
  pde_debug(w, "done (s = %d)\n", s);

  pde_debug(w, "pde_solve: computing J_phi matrix...");
  s += pde_calc_Jp(w);
  pde_debug(w, "done (s = %d)\n", s);

#if 0
  pde_debug(w, "pde_solve: computing eastward current...");
  s += pde_calc_J(w->psi, w);
  pde_debug(w, "done (s = %d)\n", s);

  /* J_lat is height-integrated and has units of current * meters */
  pde_debug(w, "pde_solve: scaling output back to physical units...");
  gsl_vector_scale(w->J_lat, w->J_s * w->r_s);
  pde_debug(w, "done\n");
#endif

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

r_s = R
E_s = 10 mV/m
sigma_s = max(sigma)
B_s = 50,000 nT

psi_s = sigma_s * E_s
U_s = E_s / B_s

Notes:

1) On output, the following parameters are scaled to dimensionless values:

w->sigma
w->B{r,t,p,f}_main
w->{m,z,v}wind
*/

static int
pde_scales(sigma_workspace *sigma_p, pde_workspace *w)
{
  int s = 0;
  size_t i, j;
  double s0_max, s1_max, s2_max;

  w->r_s = w->R;    /* m */
  w->E_s = 10.0e-3; /* V/m */
  w->B_s = 5.0e-5;  /* T */

  sigma_max(&s0_max, &s1_max, &s2_max, sigma_p);

#if 1
  w->sigma_s = GSL_MAX(s1_max, s2_max);
#else
  w->sigma_s = pow(10.0, 0.5*(log10(s0_max) + log10(GSL_MAX(s1_max, s2_max))));
#endif

  w->U_s = w->E_s / w->B_s;

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
          w->Bf_main[k] = gsl_hypot3(w->Br_main[k], w->Bt_main[k], w->Bp_main[k]);

          w->zwind[k] /= w->U_s;
          w->mwind[k] /= w->U_s;
          w->vwind[k] /= w->U_s;
        }
    }

  /*
   * these factors are chosen to make the barred equations the same
   * as the unbarred equations
   */
  w->psi_s = w->sigma_s * w->E_s;
  w->J_s = w->psi_s;

  pde_debug(w, "\n\t pde_scales: r_s     = %.12e [m]\n", w->r_s);
  pde_debug(w, "\t pde_scales: sigma_s = %.12e [S/m]\n", w->sigma_s);
  pde_debug(w, "\t pde_scales: B_s     = %.12e [T]\n", w->B_s);
  pde_debug(w, "\t pde_scales: U_s     = %.12e [m/s]\n", w->U_s);
  pde_debug(w, "\t pde_scales: E_s     = %.12e [V/m]\n", w->E_s);
  pde_debug(w, "\t pde_scales: psi_s   = %.12e [A/m2]\n", w->psi_s);
  pde_debug(w, "\t pde_scales: J_s     = %.12e [A/m2]\n", w->J_s);

  return s;
}

/*
pde_coefficients()
  Compute coefficients of PDE: f1-f5

Notes:
1) On output, the following arrays are initialized
alpha
beta
gamma
W_r
W_t
f1,...,f5
*/

static int
pde_coefficients(pde_workspace *w)
{
  int s = 0;
  const double dr = pde_dr(w);
  const double drinv = 1.0 / dr;
  const double dtheta = pde_dtheta(w);
  const double dtinv = 1.0 / dtheta;
  const double alpha_min = 1.0e-8;
  size_t i, j;

  /* compute parameters alpha, beta, gamma */
  for (i = 0; i < w->nr; ++i)
    {
      for (j = 0; j < w->ntheta; ++j)
        {
          size_t k = PDE_IDX(i, j, w);
          gsl_matrix *s = w->sigma[k];
          double s_rr = gsl_matrix_get(s, IDX_R, IDX_R);
          double s_rt = gsl_matrix_get(s, IDX_R, IDX_THETA);
          double s_tr = gsl_matrix_get(s, IDX_THETA, IDX_R);
          double s_tt = gsl_matrix_get(s, IDX_THETA, IDX_THETA);
          double s_rp = gsl_matrix_get(s, IDX_R, IDX_PHI);
          double s_tp = gsl_matrix_get(s, IDX_THETA, IDX_PHI);

          w->alpha[k] = s_rr * s_tt - s_tr * s_rt;
          w->beta[k] = s_tt * s_rp - s_rt * s_tp;
          w->gamma[k] = s_rp * s_tr - s_rr * s_tp;

          /* for low altitudes < 90km the conductivity could be 0 */
          if (!gsl_finite(w->alpha[k]) || !gsl_finite(w->beta[k]) ||
              !gsl_finite(w->gamma[k]))
            return GSL_FAILURE;

          /* since we divide by alpha to find E_r and E_t, set a lower bound */
          if (w->alpha[k] < alpha_min)
            w->alpha[k] = alpha_min;

#if 0
          /* compute W = sigma (U x B) grids */

          /* [sigma U x B]_r */
          w->W_r[k] = s_rr * (w->mwind[k] * w->Bp_main[k] -
                              w->zwind[k] * w->Bt_main[k]) +
                      s_rt * (w->zwind[k] * w->Br_main[k] -
                              w->vwind[k] * w->Bp_main[k]) +
                      s_rp * (w->vwind[k] * w->Bt_main[k] -
                              w->mwind[k] * w->Br_main[k]);

          /* [sigma U x B]_t */
          w->W_t[k] = s_tr * (w->mwind[k] * w->Bp_main[k] -
                              w->zwind[k] * w->Bt_main[k]) +
                      s_tt * (w->zwind[k] * w->Br_main[k] -
                              w->vwind[k] * w->Bp_main[k]) +
                      s_tp * (w->vwind[k] * w->Bt_main[k] -
                              w->mwind[k] * w->Br_main[k]);
#endif /*XXX*/
        }
    }

  /* compute PDE coefficients f1-5 */
  for (i = 0; i < w->nr; ++i)
    {
      double r = pde_r(i, w);
      double rsq = r * r;

      for (j = 0; j < w->ntheta; ++j)
        {
          size_t k = PDE_IDX(i, j, w);
          double theta = pde_theta(j, w);
          double sint = pde_sint(j, w);
          gsl_matrix *s = w->sigma[k];
          double s_rr = gsl_matrix_get(s, IDX_R, IDX_R);
          double s_rt = gsl_matrix_get(s, IDX_R, IDX_THETA);
          double s_tr = gsl_matrix_get(s, IDX_THETA, IDX_R);
          double s_tt = gsl_matrix_get(s, IDX_THETA, IDX_THETA);
          double dr4, dr5; /* d/dr terms in f4, f5 */
          double dt4, dt5; /* d/dtheta terms in f4, f5 */

          w->f1[k] = rsq * s_rr;
          w->f2[k] = 0.5 * r * (s_tr + s_rt);
          w->f3[k] = s_tt;

          if (i == 0)
            {
              /* use forward differences */

              dr4 = drinv *
                    (gsl_matrix_get(w->sigma[PDE_IDX(i + 1, j, w)], IDX_R, IDX_R) / w->alpha[PDE_IDX(i + 1, j, w)] -
                     s_rr / w->alpha[k]);

              dr5 = drinv *
                    ((gsl_matrix_get(w->sigma[PDE_IDX(i + 1, j, w)], IDX_THETA, IDX_R) / w->alpha[PDE_IDX(i + 1, j, w)]) / (r + dr) -
                     (s_tr / w->alpha[k]) / r);
            }
          else if (i == w->nr - 1)
            {
              /* use backward differences */

              dr4 = drinv *
                    (s_rr / w->alpha[k] -
                     gsl_matrix_get(w->sigma[PDE_IDX(i - 1, j, w)], IDX_R, IDX_R) / w->alpha[PDE_IDX(i - 1, j, w)]);

              dr5 = drinv *
                    ((s_tr / w->alpha[k]) / r -
                     (gsl_matrix_get(w->sigma[PDE_IDX(i - 1, j, w)], IDX_THETA, IDX_R) / w->alpha[PDE_IDX(i - 1, j, w)]) / (r - dr));
            }
          else
            {
              /* use central differences */

              dr4 = 0.5 * drinv *
                    (gsl_matrix_get(w->sigma[PDE_IDX(i + 1, j, w)], IDX_R, IDX_R) / w->alpha[PDE_IDX(i + 1, j, w)] -
                     gsl_matrix_get(w->sigma[PDE_IDX(i - 1, j, w)], IDX_R, IDX_R) / w->alpha[PDE_IDX(i - 1, j, w)]);

              dr5 = 0.5 * drinv *
                    ((gsl_matrix_get(w->sigma[PDE_IDX(i + 1, j, w)], IDX_THETA, IDX_R) / w->alpha[PDE_IDX(i + 1, j, w)]) / (r + dr) -
                     (gsl_matrix_get(w->sigma[PDE_IDX(i - 1, j, w)], IDX_THETA, IDX_R) / w->alpha[PDE_IDX(i - 1, j, w)]) / (r - dr));
            }

          if (j == 0)
            {
              /* use forward differences */

              dt4 = dtinv *
                    ((gsl_matrix_get(w->sigma[PDE_IDX(i, j + 1, w)], IDX_R, IDX_THETA) / w->alpha[PDE_IDX(i, j + 1, w)]) / sin(theta + dtheta) -
                     (s_rt / w->alpha[k]) / sint);

              dt5 = dtinv *
                    ((gsl_matrix_get(w->sigma[PDE_IDX(i, j + 1, w)], IDX_THETA, IDX_THETA) / w->alpha[PDE_IDX(i, j + 1, w)]) / sin(theta + dtheta) -
                     (s_tt / w->alpha[k]) / sint);
            }
          else if (j == w->ntheta - 1)
            {
              /* use backward differences */

              dt4 = dtinv *
                    ((s_rt / w->alpha[k]) / sint -
                     (gsl_matrix_get(w->sigma[PDE_IDX(i, j - 1, w)], IDX_R, IDX_THETA) / w->alpha[PDE_IDX(i, j - 1, w)]) / sin(theta - dtheta));

              dt5 = dtinv *
                    ((s_tt / w->alpha[k]) / sint -
                     (gsl_matrix_get(w->sigma[PDE_IDX(i, j - 1, w)], IDX_THETA, IDX_THETA) / w->alpha[PDE_IDX(i, j - 1, w)]) / sin(theta - dtheta));
            }
          else
            {
              /* use central differences */

              dt4 = 0.5 * dtinv *
                    ((gsl_matrix_get(w->sigma[PDE_IDX(i, j + 1, w)], IDX_R, IDX_THETA) / w->alpha[PDE_IDX(i, j + 1, w)]) / sin(theta + dtheta) -
                     (gsl_matrix_get(w->sigma[PDE_IDX(i, j - 1, w)], IDX_R, IDX_THETA) / w->alpha[PDE_IDX(i, j - 1, w)]) / sin(theta - dtheta));

              dt5 = 0.5 * dtinv *
                    ((gsl_matrix_get(w->sigma[PDE_IDX(i, j + 1, w)], IDX_THETA, IDX_THETA) / w->alpha[PDE_IDX(i, j + 1, w)]) / sin(theta + dtheta) -
                     (gsl_matrix_get(w->sigma[PDE_IDX(i, j - 1, w)], IDX_THETA, IDX_THETA) / w->alpha[PDE_IDX(i, j - 1, w)]) / sin(theta - dtheta));
            }

          w->f4[k] = r * w->alpha[k] * (r * dr4 + sint * dt4);
          w->f5[k] = w->alpha[k] * (rsq * dr5 + sint * dt5);

          if (!gsl_finite(w->f1[k]) || !gsl_finite(w->f2[k]) ||
              !gsl_finite(w->f3[k]) || !gsl_finite(w->f4[k]) ||
              !gsl_finite(w->f5[k]))
            return GSL_FAILURE;
        } /* for (j = 0; j < w->ntheta; ++j) */
    } /* for (i = 0; i < w->nr; ++i) */

  return s;
}

/*
pde_rhs_g1()
  Compute rhs g_1 of PDE

Inputs: g1 - output g_1 vector, length nr*ntheta
        w  - workspace

Notes:
1) pde_coefficients() must be called first to initialize alpha,beta,gamma arrays
*/

static int
pde_rhs_g1(gsl_vector *g1, pde_workspace *w)
{
  int s = 0;
  const double R = w->R / w->r_s;
  const double dr = pde_dr(w);
  const double drinv = 1.0 / dr;
  const double dtheta = pde_dtheta(w);
  const double dtinv = 1.0 / dtheta;
  size_t i, j;

  for (i = 0; i < w->nr; ++i)
    {
      double r = pde_r(i, w);
      double ratio = r / R;

      for (j = 0; j < w->ntheta; ++j)
        {
          size_t k = PDE_IDX(i, j, w);
          double theta = pde_theta(j, w);
          double sint = pde_sint(j, w);
          double dr1, dt1;

          if (i == 0)
            {
              /* use forward differences */

              dr1 = drinv *
                    (w->gamma[PDE_IDX(i + 1, j, w)] / w->alpha[PDE_IDX(i + 1, j, w)] -
                     w->gamma[k] / w->alpha[k]);
            }
          else if (i == w->nr - 1)
            {
              /* use backward differences */

              dr1 = drinv *
                    (w->gamma[k] / w->alpha[k] -
                     w->gamma[PDE_IDX(i - 1, j, w)] / w->alpha[PDE_IDX(i - 1, j, w)]);
            }
          else
            {
              /* use central differences */

              dr1 = 0.5 * drinv *
                    (w->gamma[PDE_IDX(i + 1, j, w)] / w->alpha[PDE_IDX(i + 1, j, w)] -
                     w->gamma[PDE_IDX(i - 1, j, w)] / w->alpha[PDE_IDX(i - 1, j, w)]);
            }

          if (j == 0)
            {
              /* use forward differences */

              dt1 = dtinv *
                    ((w->beta[PDE_IDX(i, j + 1, w)] / w->alpha[PDE_IDX(i, j + 1, w)]) / sin(theta + dtheta) -
                     (w->beta[k] / w->alpha[k]) / sint);
            }
          else if (j == w->ntheta - 1)
            {
              /* use backward differences */

              dt1 = dtinv *
                    ((w->beta[k] / w->alpha[k]) / sint -
                     (w->beta[PDE_IDX(i, j - 1, w)] / w->alpha[PDE_IDX(i, j - 1, w)]) / sin(theta - dtheta));
            }
          else
            {
              /* use central differences */

              dt1 = 0.5 * dtinv *
                    ((w->beta[PDE_IDX(i, j + 1, w)] / w->alpha[PDE_IDX(i, j + 1, w)]) / sin(theta + dtheta) -
                     (w->beta[PDE_IDX(i, j - 1, w)] / w->alpha[PDE_IDX(i, j - 1, w)]) / sin(theta - dtheta));
            }

          /* g_1 rhs vector */
          gsl_vector_set(g1, k, ratio * w->alpha[k] * (r * dr1 + sint * dt1));
        }
    }

  return s;
}

/*
pde_rhs_g2()
  Compute rhs g_2(W)

Inputs: W_r - r component of W = sigma (U x B), length nr*ntheta
        W_t - theta component of W = sigma (U x B), length nr*ntheta
        g2  - (output) output vector g_2(W)
        w   - workspace

Notes:
1) pde_coefficients() must be called first to initialize alpha,beta,gamma arrays
*/

static int
pde_rhs_g2(const gsl_vector *W_r, const gsl_vector *W_t, gsl_vector *g2, pde_workspace *w)
{
  int s = 0;
  const double R = w->R / w->r_s;
  const double dr = pde_dr(w);
  const double drinv = 1.0 / dr;
  const double dtheta = pde_dtheta(w);
  const double dtinv = 1.0 / dtheta;
  size_t i, j;

  for (i = 0; i < w->nr; ++i)
    {
      double r = pde_r(i, w);
      double ratio = r / R;

      for (j = 0; j < w->ntheta; ++j)
        {
          size_t k = PDE_IDX(i, j, w);
          double sint = pde_sint(j, w);
          gsl_matrix *sigma = w->sigma[k];
          double s_rr = gsl_matrix_get(sigma, IDX_R, IDX_R);
          double s_tr = gsl_matrix_get(sigma, IDX_THETA, IDX_R);
          double s_tt = gsl_matrix_get(sigma, IDX_THETA, IDX_THETA);
          double s_rt = gsl_matrix_get(sigma, IDX_R, IDX_THETA);
          double dr2, dt2;

          if (i == 0)
            {
              /* use forward differences */

              dr2 = drinv *
                    ((r + dr) / w->alpha[PDE_IDX(i + 1, j, w)] * (gsl_matrix_get(w->sigma[PDE_IDX(i + 1, j, w)], IDX_THETA, IDX_R) * gsl_vector_get(W_r, PDE_IDX(i + 1, j, w)) -
                                                                  gsl_matrix_get(w->sigma[PDE_IDX(i + 1, j, w)], IDX_R, IDX_R) * gsl_vector_get(W_t, PDE_IDX(i + 1, j, w))) -
                     r / w->alpha[k] * (s_tr * gsl_vector_get(W_r, k) - s_rr * gsl_vector_get(W_t, k)));
            }
          else if (i == w->nr - 1)
            {
              /* use backward differences */

              dr2 = drinv *
                    (r / w->alpha[k] * (s_tr * gsl_vector_get(W_r, k) - s_rr * gsl_vector_get(W_t, k)) -
                     (r - dr) / w->alpha[PDE_IDX(i - 1, j, w)] * (gsl_matrix_get(w->sigma[PDE_IDX(i - 1, j, w)], IDX_THETA, IDX_R) * gsl_vector_get(W_r, PDE_IDX(i - 1, j, w)) -
                                                                  gsl_matrix_get(w->sigma[PDE_IDX(i - 1, j, w)], IDX_R, IDX_R) * gsl_vector_get(W_t, PDE_IDX(i - 1, j, w))));
            }
          else
            {
              /* use central differences */

              dr2 = 0.5 * drinv *
                    ((r + dr) / w->alpha[PDE_IDX(i + 1, j, w)] * (gsl_matrix_get(w->sigma[PDE_IDX(i + 1, j, w)], IDX_THETA, IDX_R) * gsl_vector_get(W_r, PDE_IDX(i + 1, j, w)) -
                                                                  gsl_matrix_get(w->sigma[PDE_IDX(i + 1, j, w)], IDX_R, IDX_R) * gsl_vector_get(W_t, PDE_IDX(i + 1, j, w))) -
                     (r - dr) / w->alpha[PDE_IDX(i - 1, j, w)] * (gsl_matrix_get(w->sigma[PDE_IDX(i - 1, j, w)], IDX_THETA, IDX_R) * gsl_vector_get(W_r, PDE_IDX(i - 1, j, w)) -
                                                                  gsl_matrix_get(w->sigma[PDE_IDX(i - 1, j, w)], IDX_R, IDX_R) * gsl_vector_get(W_t, PDE_IDX(i - 1, j, w))));
            }

          if (j == 0)
            {
              /* use forward differences */

              dt2 = dtinv *
                    (1.0 / w->alpha[PDE_IDX(i, j + 1, w)] * (gsl_matrix_get(w->sigma[PDE_IDX(i, j + 1, w)], IDX_THETA, IDX_THETA) * gsl_vector_get(W_r, PDE_IDX(i, j + 1, w)) -
                                                             gsl_matrix_get(w->sigma[PDE_IDX(i, j + 1, w)], IDX_R, IDX_THETA) * gsl_vector_get(W_t, PDE_IDX(i, j + 1, w))) -
                     1.0 / w->alpha[k] * (s_tt * gsl_vector_get(W_r, k) - s_rt * gsl_vector_get(W_t, k)));
            }
          else if (j == w->ntheta - 1)
            {
              /* use backward differences */

              dt2 = dtinv *
                    (1.0 / w->alpha[k] * (s_tt * gsl_vector_get(W_r, k) - s_rt * gsl_vector_get(W_t, k)) -
                     1.0 / w->alpha[PDE_IDX(i, j - 1, w)] * (gsl_matrix_get(w->sigma[PDE_IDX(i, j - 1, w)], IDX_THETA, IDX_THETA) * gsl_vector_get(W_r, PDE_IDX(i, j - 1, w)) -
                                                             gsl_matrix_get(w->sigma[PDE_IDX(i, j - 1, w)], IDX_R, IDX_THETA) * gsl_vector_get(W_t, PDE_IDX(i, j - 1, w))));
            }
          else
            {
              /* use central differences */

              dt2 = 0.5 * dtinv *
                    (1.0 / w->alpha[PDE_IDX(i, j + 1, w)] * (gsl_matrix_get(w->sigma[PDE_IDX(i, j + 1, w)], IDX_THETA, IDX_THETA) * gsl_vector_get(W_r, PDE_IDX(i, j + 1, w)) -
                                                             gsl_matrix_get(w->sigma[PDE_IDX(i, j + 1, w)], IDX_R, IDX_THETA) * gsl_vector_get(W_t, PDE_IDX(i, j + 1, w))) -
                     1.0 / w->alpha[PDE_IDX(i, j - 1, w)] * (gsl_matrix_get(w->sigma[PDE_IDX(i, j - 1, w)], IDX_THETA, IDX_THETA) * gsl_vector_get(W_r, PDE_IDX(i, j - 1, w)) -
                                                             gsl_matrix_get(w->sigma[PDE_IDX(i, j - 1, w)], IDX_R, IDX_THETA) * gsl_vector_get(W_t, PDE_IDX(i, j - 1, w))));
            }

          gsl_vector_set(g2, k, ratio * ratio * w->alpha[k] * sint * (dr2 + dt2));
        }
    }

  return s;
}

/*
pde_compute_solution()
  Solve the equation A X = B for X

Notes: on output, w->residual is set to the residual ||A psi - b||
*/

static int
pde_compute_solution(const gsl_spmatrix *A, const gsl_matrix *B, gsl_matrix *X, pde_workspace *w)
{
  int s = 0;

#if PDE_SOLVER_SUPERLU

  {
    gsl_spmatrix *C = gsl_spmatrix_ccs(A);

    s = slu_proc(C, B, X, w->superlu_workspace_p);
    w->residual = w->superlu_workspace_p->rnorm[0];
    w->rcond = w->superlu_workspace_p->rcond;

    gsl_spmatrix_free(C);

    if (s)
      fprintf(stderr, "pde_compute_solution: slu_proc failed: s = %d\n", s);
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

  return s;
}

#if 0/*XXX*/
/*
pde_calc_J()
  Compute current density components and store in w->J_{r,theta,phi}. Also
compute height integrated eastward current and store in w->J_lat

J_{r,theta,phi} and J_lat will have dimensionless units

Inputs: psi - PDE solution vector, length nr*ntheta
        w   - workspace
*/

static int
pde_calc_J(const gsl_vector * psi, pde_workspace *w)
{
  size_t i, j;
  const double dr = pde_dr(w);
  const double drinv = 1.0 / dr;
  const double dtheta = pde_dtheta(w);
  const double dtinv = 1.0 / dtheta;
  const double R = w->R / w->r_s;
  double Er, Et, Ep; /* electric field components */
  double psi_r,  /* @/@r psi */
         psi_t;  /* @/@theta psi */
  double J_r, J_t, J_p;

  for (i = 0; i < w->nr; ++i)
    {
      double r = pde_r(i, w);
      double ratio = R / r;
      double ratio_sq = ratio * ratio;

      for (j = 0; j < w->ntheta; ++j)
        {
          size_t k = PDE_IDX(i, j, w);
          double sint = pde_sint(j, w);
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
              psi_r = drinv * (PSI_GET(psi, i + 1, j, w) -
                               PSI_GET(psi, i, j, w));
            }
          else if (i == w->nr - 1)
            {
              psi_r = drinv * (PSI_GET(psi, i, j, w) -
                               PSI_GET(psi, i - 1, j, w));
            }
          else
            {
              psi_r = 0.5 * drinv * (PSI_GET(psi, i + 1, j, w) -
                                     PSI_GET(psi, i - 1, j, w));
            }

          if (j == 0)
            {
              psi_t = dtinv * (PSI_GET(psi, i, j + 1, w) -
                               PSI_GET(psi, i, j, w));
            }
          else if (j == w->ntheta - 1)
            {
              psi_t = dtinv * (PSI_GET(psi, i, j, w) -
                               PSI_GET(psi, i, j - 1, w));
            }
          else
            {
              psi_t = 0.5 * dtinv * (PSI_GET(psi, i, j + 1, w) -
                                     PSI_GET(psi, i, j - 1, w));
            }

          J_r = ratio_sq / sint * psi_t;
          J_t = -ratio * R / sint * psi_r;

          Ep = E_phi(i, j, w);

          Er = (s_tt * J_r - s_rt * J_t - w->beta[k] * Ep) / w->alpha[k];
          Et = -(s_tr * J_r - s_rr * J_t - w->gamma[k] * Ep) / w->alpha[k];

          if (w->compute_winds)
            {
              Er += (-s_tt * w->W_r[k] + s_rt * w->W_t[k]) / w->alpha[k];
              Et -= (-s_tr * w->W_r[k] + s_rr * w->W_t[k]) / w->alpha[k];
            }

          J_p = s_pr * Er + s_pt * Et + s_pp * Ep;

          if (w->compute_winds)
            {
              J_p += s_pr * (w->mwind[k] * w->Bp_main[k] -
                             w->zwind[k] * w->Bt_main[k]) +
                     s_pt * (w->zwind[k] * w->Br_main[k] -
                             w->vwind[k] * w->Bp_main[k]) +
                     s_pp * (w->vwind[k] * w->Bt_main[k] -
                             w->mwind[k] * w->Br_main[k]);
            }

          /* store currents */
          gsl_matrix_set(w->J_r, i, j, J_r);
          gsl_matrix_set(w->J_theta, i, j, J_t);
          gsl_matrix_set(w->J_phi, i, j, J_p);

          /* store electric field */
          gsl_matrix_set(w->E_r, i, j, Er);
          gsl_matrix_set(w->E_theta, i, j, Et);
          gsl_matrix_set(w->E_phi, i, j, Ep);
        }
    }

  /* compute height integrated current density */

  for (j = 0; j < w->ntheta; ++j)
    {
      double J_sum = 0.0;

      for (i = 0; i < w->nr; ++i)
        J_sum += gsl_matrix_get(w->J_phi, i, j) * dr;

      gsl_vector_set(w->J_lat, j, J_sum);
    }

  return GSL_SUCCESS;
}
#endif

/*
pde_calc_W()
  Compute W_r, W_t, and W_p matrices

W_r = [ 0 ; W(u_0)_r ]
W_t = [ 0 ; W(u_0)_t ]
W_p = [ 0 ; W(u_0)_p ]

where u_0 is the HWM wind field and W(u) = sigma*(u x B)

Inputs: w   - workspace

Notes:
1) On output, the w->WR, w->WTHETA and w->WPHI matrices are filled
*/

static int
pde_calc_W(pde_workspace *w)
{
  size_t i, j;

  for (i = 0; i < w->nr; ++i)
    {
      for (j = 0; j < w->ntheta; ++j)
        {
          size_t k = PDE_IDX(i, j, w);
          gsl_matrix *s = w->sigma[k];
          double s_rr = gsl_matrix_get(s, IDX_R, IDX_R);
          double s_rt = gsl_matrix_get(s, IDX_R, IDX_THETA);
          double s_rp = gsl_matrix_get(s, IDX_R, IDX_PHI);
          double s_tr = gsl_matrix_get(s, IDX_THETA, IDX_R);
          double s_tt = gsl_matrix_get(s, IDX_THETA, IDX_THETA);
          double s_tp = gsl_matrix_get(s, IDX_THETA, IDX_PHI);
          double s_pr = gsl_matrix_get(s, IDX_PHI, IDX_R);
          double s_pt = gsl_matrix_get(s, IDX_PHI, IDX_THETA);
          double s_pp = gsl_matrix_get(s, IDX_PHI, IDX_PHI);
          double W_r, W_t, W_p;

          /* column corresponding to E_{phi0} is 0 */
          gsl_matrix_set(w->WR, k, 0, 0.0);
          gsl_matrix_set(w->WTHETA, k, 0, 0.0);
          gsl_matrix_set(w->WPHI, k, 0, 0.0);

          W_r = s_rr * (w->mwind[k] * w->Bp_main[k] - w->zwind[k] * w->Bt_main[k]) +
                s_rt * (w->zwind[k] * w->Br_main[k] - w->vwind[k] * w->Bp_main[k]) +
                s_rp * (w->vwind[k] * w->Bt_main[k] - w->mwind[k] * w->Br_main[k]);

          W_t = s_tr * (w->mwind[k] * w->Bp_main[k] - w->zwind[k] * w->Bt_main[k]) +
                s_tt * (w->zwind[k] * w->Br_main[k] - w->vwind[k] * w->Bp_main[k]) +
                s_tp * (w->vwind[k] * w->Bt_main[k] - w->mwind[k] * w->Br_main[k]);

          W_p = s_pr * (w->mwind[k] * w->Bp_main[k] - w->zwind[k] * w->Bt_main[k]) +
                s_pt * (w->zwind[k] * w->Br_main[k] - w->vwind[k] * w->Bp_main[k]) +
                s_pp * (w->vwind[k] * w->Bt_main[k] - w->mwind[k] * w->Br_main[k]);

          gsl_matrix_set(w->WR, k, 1, W_r);
          gsl_matrix_set(w->WTHETA, k, 1, W_t);
          gsl_matrix_set(w->WPHI, k, 1, W_p);
        }
    }

  return GSL_SUCCESS;
}

/*
pde_calc_Jr()
  Compute J_r component of current density

J_r = (R/r)^2 1/sin(theta) d/dtheta psi

Inputs: psi - PDE solution vector, length nr*ntheta
        Jr  - (output) J_r vector (dimensionless), length nr*ntheta
        w   - workspace
*/

static int
pde_calc_Jr(const gsl_vector * psi, gsl_vector * Jr, pde_workspace *w)
{
  size_t i, j;
  const double dtheta = pde_dtheta(w);
  const double dtinv = 1.0 / dtheta;
  const double R = w->R / w->r_s;
  double psi_t;  /* @/@theta psi */

  for (i = 0; i < w->nr; ++i)
    {
      double r = pde_r(i, w);
      double ratio = R / r;
      double ratio_sq = ratio * ratio;

      for (j = 0; j < w->ntheta; ++j)
        {
          size_t k = PDE_IDX(i, j, w);
          double sint = pde_sint(j, w);

          if (j == 0)
            {
              psi_t = dtinv * (PSI_GET(psi, i, j + 1, w) -
                               PSI_GET(psi, i, j, w));
            }
          else if (j == w->ntheta - 1)
            {
              psi_t = dtinv * (PSI_GET(psi, i, j, w) -
                               PSI_GET(psi, i, j - 1, w));
            }
          else
            {
              psi_t = 0.5 * dtinv * (PSI_GET(psi, i, j + 1, w) -
                                     PSI_GET(psi, i, j - 1, w));
            }

          gsl_vector_set(Jr, k, ratio_sq / sint * psi_t);
        }
    }

  return GSL_SUCCESS;
}

/*
pde_calc_Jt()
  Compute J_t component of current density

J_t = -(R/r) R/sin(theta) d/dr psi

Inputs: psi - PDE solution vector, length nr*ntheta
        Jt  - (output) J_t vector (dimensionless), length nr*ntheta
        w   - workspace
*/

static int
pde_calc_Jt(const gsl_vector * psi, gsl_vector * Jt, pde_workspace *w)
{
  size_t i, j;
  const double dr = pde_dr(w);
  const double drinv = 1.0 / dr;
  const double R = w->R / w->r_s;
  double psi_r;  /* @/@r psi */

  for (i = 0; i < w->nr; ++i)
    {
      double r = pde_r(i, w);
      double ratio = R / r;

      for (j = 0; j < w->ntheta; ++j)
        {
          size_t k = PDE_IDX(i, j, w);
          double sint = pde_sint(j, w);

          if (i == 0)
            {
              psi_r = drinv * (PSI_GET(psi, i + 1, j, w) -
                               PSI_GET(psi, i, j, w));
            }
          else if (i == w->nr - 1)
            {
              psi_r = drinv * (PSI_GET(psi, i, j, w) -
                               PSI_GET(psi, i - 1, j, w));
            }
          else
            {
              psi_r = 0.5 * drinv * (PSI_GET(psi, i + 1, j, w) -
                                     PSI_GET(psi, i - 1, j, w));
            }

          gsl_vector_set(Jt, k, -ratio * R / sint * psi_r);
        }
    }

  return GSL_SUCCESS;
}

/*
pde_calc_Jp()
  Compute J_p component of current density

J_p = s_pr E_r + s_pt E_t + s_pp E_p + W_p

Also compute height-integrated J_p:

J^{HI}_p(theta) = sum_{i=1}^{nr} J_p(r_i,theta)

Inputs: w - workspace

Notes:
1) pde_calc_E() must be called first to compute electric field matrices
2) pde_calc_W() must be called first to compute wind matrices
3) Output is stored in w->JPHI
4) Height-integrated J_phi is stored in w->JPHI_HI
*/

static int
pde_calc_Jp(pde_workspace *w)
{
  const double dr = pde_dr(w);
  size_t h, i, j;

  for (h = 0; h < w->nrhs; ++h)
    {
      for (i = 0; i < w->nr; ++i)
        {
          for (j = 0; j < w->ntheta; ++j)
            {
              size_t k = PDE_IDX(i, j, w);
              gsl_matrix *s = w->sigma[k];
              double s_pr = gsl_matrix_get(s, IDX_PHI, IDX_R);
              double s_pt = gsl_matrix_get(s, IDX_PHI, IDX_THETA);
              double s_pp = gsl_matrix_get(s, IDX_PHI, IDX_PHI);
              double Er = gsl_matrix_get(w->ER, k, h);
              double Et = gsl_matrix_get(w->ETHETA, k, h);
              double Ep = gsl_matrix_get(w->EPHI, k, h);
              double Wp = gsl_matrix_get(w->WPHI, k, h);
              double Jp = s_pr * Er + s_pt * Et + s_pp * Ep + Wp;

              gsl_matrix_set(w->JPHI, k, h, Jp);
            }
        }

      /* compute height integrated J_phi */
      for (j = 0; j < w->ntheta; ++j)
        {
          double J_sum = 0.0;

          for (i = 0; i < w->nr; ++i)
            {
              size_t k = PDE_IDX(i, j, w);
              J_sum += gsl_matrix_get(w->JPHI, k, h) * dr;
            }

          gsl_matrix_set(w->JPHI_HI, j, h, J_sum);
        }
    }

  return GSL_SUCCESS;
}

/*
pde_calc_E()
  Compute electric field matrices E_r and E_t

Inputs: w   - workspace

Notes:
1) w->ER, w->ETHETA and w->EPHI are filled
*/

static int
pde_calc_E(pde_workspace *w)
{
  size_t h, i, j;
  const double R = w->R / w->r_s;

  for (h = 0; h < w->nrhs; ++h)
    {
      for (i = 0; i < w->nr; ++i)
        {
          double r = pde_r(i, w);
          double ratio = R / r;

          for (j = 0; j < w->ntheta; ++j)
            {
              size_t k = PDE_IDX(i, j, w);
              double Jr = gsl_matrix_get(w->JR, k, h);
              double Jt = gsl_matrix_get(w->JTHETA, k, h);
              double Wr = gsl_matrix_get(w->WR, k, h);
              double Wt = gsl_matrix_get(w->WTHETA, k, h);
              double sint = pde_sint(j, w);
              gsl_matrix *s = w->sigma[k];
              double s_rr = gsl_matrix_get(s, IDX_R, IDX_R);
              double s_rt = gsl_matrix_get(s, IDX_R, IDX_THETA);
              double s_tr = gsl_matrix_get(s, IDX_THETA, IDX_R);
              double s_tt = gsl_matrix_get(s, IDX_THETA, IDX_THETA);
              double Er, Et, Ep;

              Er = (s_tt * (Jr - Wr) - s_rt * (Jt - Wt)) / w->alpha[k];
              Et = -(s_tr * (Jr - Wr) - s_rr * (Jt - Wt)) / w->alpha[k];
              Ep = 0.0;

              if (h == 0)
                {
                  Ep = ratio / sint;

                  /* add term multiplying E_{phi0} */
                  Er -= w->beta[k] / w->alpha[k] * Ep;
                  Et += w->gamma[k] / w->alpha[k] * Ep;
                }

              /* store electric field */
              gsl_matrix_set(w->ER, k, h, Er);
              gsl_matrix_set(w->ETHETA, k, h, Et);
              gsl_matrix_set(w->EPHI, k, h, Ep);
            }
        }
    }

  return GSL_SUCCESS;
}

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
pde2d_coefs()
  Return coefficients of PDE for a given (r,theta)
*/

static int
pde2d_coefs(const double r, const double theta, double coef[GSL_PDE2D_COEF_TOTAL], void *params)
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
  coef[GSL_PDE2D_COEF_FXY] = 2.0 * w->f2[k];
  coef[GSL_PDE2D_COEF_FYY] = w->f3[k];
  coef[GSL_PDE2D_COEF_FX] = w->f4[k];
  coef[GSL_PDE2D_COEF_FY] = w->f5[k];
  coef[GSL_PDE2D_COEF_F] = 0.0;

  return GSL_SUCCESS;
}

/*
pde2d_rhs()
  Return coefficients of PDE for a given (r,theta)
*/

static int
pde2d_rhs(const double r, const double theta, double *rhsval, void *params)
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

  if (w->flags & PDE_FLG_RHS_G1)
    *rhsval = gsl_matrix_get(w->G, k, 0);
  else if (w->flags & PDE_FLG_RHS_G2)
    *rhsval = gsl_matrix_get(w->G, k, 1);
  else
    {
      GSL_ERROR("invalid rhs flag", GSL_EINVAL);
    }

  return GSL_SUCCESS;
}

static int
pde2d_bnd(const gsl_pde2d_bnd_t type, const double *v, gsl_vector *alpha, gsl_vector *beta,
          gsl_vector *gamma, gsl_vector *delta, void *params)
{
  const size_t n = alpha->size;
  size_t i;

  (void) v;
  (void) params;

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

#if 0/*XXX*/
          for (k = 0; k < nphi; ++k)
            {
              magfield_current_set(MAG_IDX_R, i, j, k, Jr * w->J_s, magfield_p);
              magfield_current_set(MAG_IDX_THETA, i, j, k, Jt * w->J_s, magfield_p);
              magfield_current_set(MAG_IDX_PHI, i, j, k, Jp * w->J_s, magfield_p);
            }
#endif
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
