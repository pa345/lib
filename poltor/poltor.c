/*
 * poltor.c
 *
 * Decompose satellite magnetic field measurements to toroidal
 * currents mapped onto a thin shell in the ionosphere
 *
 * Calling sequence:
 *
 * 1. magdata_alloc     - allocate data structure for internal storage
 *                        of satellite magnetic measurements
 * 2. magdata_add       - add single datum to internal data structure
 * 3. magdata_calc      - calculate spatial weights for all data
 * 4. poltor_alloc      - allocate a workspace, specify thin current shell
 *                        radii, maximum spherical harmonic degree,
 *                        and previously constructed data structure
 * 5. poltor_calc       - calculate LS matrix, rhs vector, and store both
 *                        to a file for future quick manipulation (LLS.dat)
 * 6. poltor_solve      - regularize and solve LS system for coefficients
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#include <complex.h>
#include <sys/time.h>
#include <assert.h>
#include <omp.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <common/common.h>
#include <common/oct.h>
#include <indices/indices.h>
#include <track/track_weight.h>

#include "green.h"
#include "lls.h"
#include "magdata_list.h"
#include "poltor.h"

static int poltor_row(const double r, const double theta, const double phi,
                      const double r_ns, const double theta_ns, const double phi_ns,
                      gsl_vector_complex *x, gsl_vector_complex *y, gsl_vector_complex *z,
                      gsl_vector_complex *dx, gsl_vector_complex *dy, gsl_vector_complex *dz,
                      poltor_workspace *w);
static int poltor_row_pint(const double r, const double theta,
                           const double r_ns, const double theta_ns,
                           gsl_vector_complex *vx, gsl_vector_complex *vy, gsl_vector_complex *vz,
                           gsl_vector_complex *vdx, gsl_vector_complex *vdy, gsl_vector_complex *vdz,
                           poltor_workspace *w);
static int poltor_row_pext(const double r, const double theta,
                           const double r_ns, const double theta_ns,
                           gsl_vector_complex *vx, gsl_vector_complex *vy, gsl_vector_complex *vz,
                           gsl_vector_complex *vdx, gsl_vector_complex *vdy, gsl_vector_complex *vdz,
                           poltor_workspace *w);
static int poltor_row_psh(const double r, const double theta,
                          const double r_ns, const double theta_ns,
                          gsl_vector_complex *vx, gsl_vector_complex *vy, gsl_vector_complex *vz,
                          gsl_vector_complex *vdx, gsl_vector_complex *vdy, gsl_vector_complex *vdz,
                          poltor_workspace *w);
static int poltor_row_tor(const double r, const double theta,
                          const double r_ns, const double theta_ns,
                          gsl_vector_complex *vx, gsl_vector_complex *vy, gsl_vector_complex *vz,
                          gsl_vector_complex *vdx, gsl_vector_complex *vdy, gsl_vector_complex *vdz,
                          poltor_workspace *w);
static int poltor_eval_chi(const double r, const double theta, const double phi,
                           const size_t idx, double *chi, poltor_workspace *w);
static int poltor_init(poltor_workspace * w);
static int poltor_init_weights(poltor_workspace * w);
static int poltor_init_solar_flux(poltor_workspace * w);

/*
poltor_alloc()
  Allocate a poltor workspace

Inputs: params - parameters
          R        - reference radius (km)
          nmax_int - maximum internal spherical harmonic degree
          mmax_int - maximum internal spherical harmonic order
          nmax_ext - maximum external spherical harmonic degree
          mmax_ext - maximum external spherical harmonic order
          nmax_tor - maximum spherical harmonic degree for B_tor
          mmax_tor - maximum spherical harmonic order for B_tor
          data     - satellite data to invert
*/

poltor_workspace *
poltor_alloc(const poltor_parameters *params)
{
  poltor_workspace *w;

  w = calloc(1, sizeof(poltor_workspace));
  if (!w)
    return 0;

  w->params = *params;

  w->R = params->R;
  w->b = params->b;
  w->d = params->d;
  w->rmin = params->rmin;
  w->rmax = params->rmax;
  w->nmax_int = params->nmax_int;
  w->mmax_int = params->mmax_int;
  w->nmax_ext = params->nmax_ext;
  w->mmax_ext = params->mmax_ext;
  w->nmax_sh = params->nmax_sh;
  w->mmax_sh = params->mmax_sh;
  w->nmax_tor = params->nmax_tor;
  w->mmax_tor = params->mmax_tor;
  w->shell_J = params->shell_J;
  w->data = params->data;
  w->alpha_int = params->alpha_int;
  w->alpha_sh = params->alpha_sh;
  w->alpha_tor = params->alpha_tor;
  w->flags = params->flags;

  w->rmid = 0.5 * (w->rmin + w->rmax);

  w->nmax_max = GSL_MAX(w->nmax_int,
                  GSL_MAX(w->nmax_ext,
                    GSL_MAX(w->nmax_sh, w->nmax_tor)));
  w->mmax_max = GSL_MAX(w->mmax_int,
                  GSL_MAX(w->mmax_ext,
                    GSL_MAX(w->mmax_sh, w->mmax_tor)));

  if (w->data)
    {
      /* count data and find total number of residuals */
      magdata_list_index(w->data, w->res_cnt);
      w->n = w->res_cnt[MAGDATA_LIST_IDX_TOTAL];
    }
  else
    {
      /* minimal workspace for poltor_eval_B() */
      w->n = 6;
    }

  /* count number of spherical harmonic coefficients for each source */
  w->nnm_int = green_calc_nnm(w->nmax_int, w->mmax_int);
  w->nnm_ext = green_calc_nnm(w->nmax_ext, w->mmax_ext);
  w->nnm_sh = green_calc_nnm(w->nmax_sh, w->mmax_sh);
  w->nnm_tor = green_calc_nnm(w->nmax_tor, w->mmax_tor);

  /* compute total poloidal shell coefficients accouting for Taylor expansion */
  w->p_pint = w->nnm_int;
  w->p_pext = w->nnm_ext;
  w->p_psh = w->nnm_sh * (w->shell_J + 1);
  w->p_tor = w->nnm_tor;

  w->p = w->p_pint + w->p_pext + w->p_psh + w->p_tor;

  w->pint_offset = 0;
  w->pext_offset = w->pint_offset + w->p_pint;
  w->psh_offset = w->pext_offset + w->p_pext;
  w->tor_offset = w->psh_offset + w->p_psh;

  if (w->n > 0)
    {
#if 0/*XXX*/
      w->weights = gsl_vector_alloc(w->nblock);
      if (!w->weights)
        {
          fprintf(stderr, "poltor_alloc: cannot allocate weights vector: %s\n",
                  strerror(errno));
          poltor_free(w);
          return 0;
        }
#endif

      w->dof = (int) w->n - (int) w->p;
    }

  w->c = gsl_vector_complex_alloc(w->p);
  if (!w->c)
    {
      fprintf(stderr, "poltor_alloc: cannot allocate coefficient vector: %s\n",
              strerror(errno));
      poltor_free(w);
      return 0;
    }

  w->L = gsl_vector_alloc(w->p);
  if (!w->L)
    {
      fprintf(stderr, "poltor_alloc: cannot allocate Tikhonov vector: %s\n",
              strerror(errno));
      poltor_free(w);
      return 0;
    }

  /* initialize to identity matrix */
  gsl_vector_set_all(w->L, 1.0);

  w->residuals = gsl_vector_complex_alloc(w->n);
  if (!w->residuals)
    {
      fprintf(stderr, "poltor_alloc: cannot allocate residual vector: %s\n",
              strerror(errno));
      poltor_free(w);
      return 0;
    }

  w->solar_flux = gsl_vector_alloc(w->n);
  w->solar_flux_grad = gsl_vector_alloc(w->n);

  w->wts_spatial = gsl_vector_alloc(w->n);
  w->wts_robust = gsl_vector_alloc(w->n);
  w->wts_final = gsl_vector_alloc(w->n);
  w->f = gsl_vector_alloc(w->n);

  w->cquad_workspace_p = gsl_integration_cquad_workspace_alloc(1000);
  if (!w->cquad_workspace_p)
    {
      fprintf(stderr, "poltor_alloc: cannot allocate integration workspace: %s\n",
              strerror(errno));
      poltor_free(w);
      return 0;
    }

  w->nreg = 200;
  w->reg_param = gsl_vector_alloc(w->nreg);
  w->rho = gsl_vector_alloc(w->nreg);
  w->eta = gsl_vector_alloc(w->nreg);

  w->max_threads = (size_t) omp_get_max_threads();

  w->omp_dX = gsl_matrix_complex_alloc(w->max_threads, w->p);
  w->omp_dY = gsl_matrix_complex_alloc(w->max_threads, w->p);
  w->omp_dZ = gsl_matrix_complex_alloc(w->max_threads, w->p);
  w->omp_dX_grad = gsl_matrix_complex_alloc(w->max_threads, w->p);
  w->omp_dY_grad = gsl_matrix_complex_alloc(w->max_threads, w->p);
  w->omp_dZ_grad = gsl_matrix_complex_alloc(w->max_threads, w->p);
  w->omp_rowidx = malloc(w->max_threads * sizeof(size_t));
  w->omp_nrows = malloc(w->max_threads * sizeof(size_t));

  w->JHJ = gsl_matrix_complex_alloc(w->p, w->p);
  w->JHf = gsl_vector_complex_alloc(w->p);

  {
    const size_t nmax_max = GSL_MAX(w->nmax_int, GSL_MAX(w->nmax_ext, GSL_MAX(w->nmax_sh, w->nmax_tor)));
    const size_t mmax_max = GSL_MAX(w->mmax_int, GSL_MAX(w->mmax_ext, GSL_MAX(w->mmax_sh, w->mmax_tor)));
    size_t i;

    /*
     * maximum observations to accumulate at once in LS system, calculated to make
     * each omp_J matrix approximately of size 'POLTOR_MATRIX_SIZE'
     */
    w->nblock = POLTOR_MATRIX_SIZE / (w->p * sizeof(double));

    w->green_p = malloc(w->max_threads * sizeof(green_complex_workspace *));
    w->green_grad_p = malloc(w->max_threads * sizeof(green_complex_workspace *));
    w->omp_J = malloc(w->max_threads * sizeof(gsl_matrix_complex *));
    w->omp_f = malloc(w->max_threads * sizeof(gsl_vector_complex *));
    w->lls_workspace_p = lls_complex_alloc(w->nblock, w->p);

    for (i = 0; i < w->max_threads; ++i)
      {
        w->green_p[i] = green_complex_alloc(nmax_max, mmax_max, w->R);
        w->green_grad_p[i] = green_complex_alloc(nmax_max, mmax_max, w->R);
        w->omp_J[i] = gsl_matrix_complex_alloc(w->nblock, w->p);
        w->omp_f[i] = gsl_vector_complex_alloc(w->nblock);
      }
  }

  w->f107_workspace_p = f107_alloc(F107_IDX_FILE);

  /* determine if system is linear */
  if (w->res_cnt[MAGDATA_LIST_IDX_F] == 0 &&
      w->res_cnt[MAGDATA_LIST_IDX_DF_NS] == 0 &&
      w->res_cnt[MAGDATA_LIST_IDX_DF_EW] == 0)
    {
      w->lls_solution = 1;
    }
  else
    {
      w->lls_solution = 0;
    }

  if (w->data)
    {
      /* compute spatial weights and solar flux array */
      poltor_init(w);
    }

  return w;
} /* poltor_alloc() */

void
poltor_free(poltor_workspace *w)
{
  if (w->weights)
    gsl_vector_free(w->weights);

  if (w->solar_flux)
    gsl_vector_free(w->solar_flux);

  if (w->solar_flux_grad)
    gsl_vector_free(w->solar_flux_grad);

  if (w->f)
    gsl_vector_free(w->f);

  if (w->wts_spatial)
    gsl_vector_free(w->wts_spatial);

  if (w->wts_robust)
    gsl_vector_free(w->wts_robust);

  if (w->wts_final)
    gsl_vector_free(w->wts_final);

  if (w->c)
    gsl_vector_complex_free(w->c);

  if (w->L)
    gsl_vector_free(w->L);

  if (w->residuals)
    gsl_vector_complex_free(w->residuals);

  if (w->lls_workspace_p)
    lls_complex_free(w->lls_workspace_p);

  if (w->cquad_workspace_p)
    gsl_integration_cquad_workspace_free(w->cquad_workspace_p);

  if (w->reg_param)
    gsl_vector_free(w->reg_param);

  if (w->rho)
    gsl_vector_free(w->rho);

  if (w->eta)
    gsl_vector_free(w->eta);

  if (w->omp_dX)
    gsl_matrix_complex_free(w->omp_dX);

  if (w->omp_dY)
    gsl_matrix_complex_free(w->omp_dY);

  if (w->omp_dZ)
    gsl_matrix_complex_free(w->omp_dZ);

  if (w->omp_dX_grad)
    gsl_matrix_complex_free(w->omp_dX_grad);

  if (w->omp_dY_grad)
    gsl_matrix_complex_free(w->omp_dY_grad);

  if (w->omp_dZ_grad)
    gsl_matrix_complex_free(w->omp_dZ_grad);

  if (w->omp_rowidx)
    free(w->omp_rowidx);

  if (w->omp_nrows)
    free(w->omp_nrows);

  if (w->JHJ)
    gsl_matrix_complex_free(w->JHJ);

  if (w->JHf)
    gsl_vector_complex_free(w->JHf);

  {
    size_t i;

    for (i = 0; i < w->max_threads; ++i)
      {
        if (w->green_p[i])
          green_complex_free(w->green_p[i]);

        if (w->green_grad_p[i])
          green_complex_free(w->green_grad_p[i]);

        if (w->omp_J[i])
          gsl_matrix_complex_free(w->omp_J[i]);

        if (w->omp_f[i])
          gsl_vector_complex_free(w->omp_f[i]);
      }

    free(w->green_p);
    free(w->green_grad_p);
    free(w->omp_J);
    free(w->omp_f);
  }

  if (w->f107_workspace_p)
    f107_free(w->f107_workspace_p);

  free(w);
}

/* initialize parameter structure */
int
poltor_init_params(poltor_parameters * params)
{
  params->R = -1.0;
  params->nmax_int = 0;
  params->mmax_int = 0;
  params->nmax_sh = 0;
  params->mmax_sh = 0;
  params->nmax_tor = 0;
  params->mmax_tor = 0;
  params->max_iter = 0;
  params->regularize = 0;
  params->use_weights = 0;
  params->solar_flux_factor = 0;
  params->alpha_int = 0.0;
  params->alpha_sh = 0.0;
  params->alpha_tor = 0.0;
  params->weight_X = 0.0;
  params->weight_Y = 0.0;
  params->weight_Z = 0.0;
  params->weight_F = 0.0;
  params->weight_DX_NS = 0.0;
  params->weight_DY_NS = 0.0;
  params->weight_DZ_NS = 0.0;
  params->weight_DF_NS = 0.0;
  params->fit_X = 1;
  params->fit_Y = 1;
  params->fit_Z = 1;
  params->fit_F = 1;
  params->fit_DX_NS = 1;
  params->fit_DY_NS = 1;
  params->fit_DZ_NS = 1;
  params->fit_DF_NS = 1;
  params->synth_data = 0;

  return 0;
}

/*
poltor_calc()
  Build LS system: A^H W A and A^H W b

Notes:
1) On output, w->L contains diag(L) the regularization matrix
*/

int
poltor_calc(poltor_workspace *w)
{
  int s = 0;
  const char *lls_file = "LLS.dat";

  fprintf(stderr, "poltor_calc: vector residuals      = %zu\n", w->res_cnt[MAGDATA_LIST_IDX_VECTOR]);
  fprintf(stderr, "poltor_calc: scalar residuals      = %zu\n", w->res_cnt[MAGDATA_LIST_IDX_F]);
  fprintf(stderr, "poltor_calc: vector grad residuals = %zu\n", w->res_cnt[MAGDATA_LIST_IDX_VGRAD]);
  fprintf(stderr, "poltor_calc: scalar grad residuals = %zu\n",
    w->res_cnt[MAGDATA_LIST_IDX_DF_NS] + w->res_cnt[MAGDATA_LIST_IDX_DF_EW]);
  fprintf(stderr, "poltor_calc: total residuals       = %zu\n", w->n);
  fprintf(stderr, "poltor_calc: ncoeff                = %zu\n", w->p);
  fprintf(stderr, "poltor_calc: nres / ncoeff         = %.1f\n",
          (double) w->n / (double) w->p);

#if 0
  /* initialize lls module */
  lls_complex_reset(w->lls_workspace_p);

#if !POLTOR_SYNTH_DATA

  if (w->params.flags & POLTOR_FLG_REGULARIZE)
    {
      fprintf(stderr, "poltor_calc: computing regularization matrix...");

      /* construct L = diag(L) */
      s = poltor_regularize(w->L, w);
      if (s)
        return s;

      fprintf(stderr, "done\n");
    }

#endif

  fprintf(stderr, "poltor_calc: building LS system...");

  s += poltor_build_ls(1, w);
  if (s)
    return s;

  fprintf(stderr, "done\n");

  fprintf(stderr, "poltor_calc: writing LS system to %s...",
          lls_file);
  s = lls_complex_save(lls_file, w->lls_workspace_p);
  fprintf(stderr, "done\n");

  printc_octave(w->lls_workspace_p->AHA, "AHA");
  printcv_octave(w->lls_workspace_p->AHb, "AHb");
#endif/*XXX*/

  return s;
} /* poltor_calc() */

/*
poltor_solve()
  Solve LS system previously constructed by poltor_calc()

Inputs: w - workspace

Notes:
1) On output, w->c contains coefficient vector

2) On output, w->chisq contains the residual sum of squares r^T r
*/

int
poltor_solve(poltor_workspace *w)
{
  int s = 0;
  struct timeval tv0, tv1;
  const double dof = (double) w->n - (double) w->p;
  double lambda = 0.0; /* regularization parameter */

#if !POLTOR_SYNTH_DATA

  if (w->params.flags & POLTOR_FLG_REGULARIZE)
    {
      const char *lcurve_file = "lcurve.dat";

      /* construct L = diag(L); this must be done again here
       * in case user loads LLS system from a file
       */
      s = poltor_regularize(w->L, w);
      if (s)
        return s;

      /* calculate L-curve */
      fprintf(stderr, "poltor_solve: calculating L-curve...");
      lls_complex_lcurve(w->reg_param, w->rho, w->eta, w->lls_workspace_p);
      fprintf(stderr, "done\n");

      /* find L-corve corner */
      fprintf(stderr, "poltor_solve: calculating L-curve corner...");
      gsl_multifit_linear_lcorner(w->rho, w->eta, &(w->reg_idx));
      lambda = gsl_vector_get(w->reg_param, w->reg_idx);
      fprintf(stderr, "done (lambda = %f)\n", lambda);

      fprintf(stderr, "poltor_solve: writing L-curve to %s...", lcurve_file);
      poltor_print_lcurve(lcurve_file, w);
      fprintf(stderr, "done\n");
    }

#endif

  fprintf(stderr, "poltor_solve: solving LS system...");
  gettimeofday(&tv0, NULL);

  /* solve A^H W A c = A^H W b */
  s = lls_complex_solve(lambda, w->c, w->lls_workspace_p);
  if (s)
    {
      fprintf(stderr, "poltor_solve: error solving system: %d\n", s);
      return s;
    }

  /* backtransform standard form system to get original solution vector */
  lls_complex_btransform(w->L, w->c, w->lls_workspace_p);

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds, residual = %.12e, chisq/dof = %.12e)\n",
          time_diff(tv0, tv1),
          w->lls_workspace_p->residual,
          w->lls_workspace_p->chisq / dof);

  fprintf(stderr, "poltor_solve: matrix condition number = %f\n",
          w->lls_workspace_p->cond);

  w->chisq = w->lls_workspace_p->chisq;

  return s;
} /* poltor_solve() */

/*
poltor_eval_J_shell()
  Calculate shell toroidal current vector at a given point

Inputs: r     - radius (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        J     - (output) shell current density vector in units of A/m^2
                (J_r,J_theta,J_phi) - J_r = 0 for toroidal currents
        w     - workspace

Notes:
1) J should be initialized prior to calling this function, since
on output, J <- J + J_tor
*/

int
poltor_eval_J_shell(const double r, const double theta, const double phi,
                    double J[3], poltor_workspace *w)
{
  int s = 0;
  size_t n, j;
  int m;
  const double invsint = 1.0 / sin(theta);
  const double ratio = (r - w->rmid) / w->R;
  double rterm = 1.0; /* [(r - r_0) / R]^j */
  complex double Jt = 0.0, Jp = 0.0;
  green_complex_workspace *green_p = w->green_p[0];
  complex double *Ynm = green_p->Ynm;
  complex double *dYnm = green_p->dYnm;
  complex double Kt = 0.0, Kp = 0.0;

  /* compute Ynm and d/dtheta Ynm */
  green_complex_Ynm_deriv(theta, phi, green_p);

  for (j = 0; j <= w->shell_J; ++j)
    {
      for (n = 1; n <= w->nmax_sh; ++n)
        {
          int M = (int) GSL_MIN(n, w->mmax_sh);

          for (m = -M; m <= M; ++m)
            {
              int mabs = abs(m);
              size_t nmidx = poltor_jnmidx(j, n, m, w);
              size_t pidx = gsl_sf_legendre_array_index(n, mabs);
              gsl_complex cnmj = gsl_vector_complex_get(w->c, nmidx);
              complex double qnmj = GSL_REAL(cnmj) + I * GSL_IMAG(cnmj);
              complex double Y_nm, dY_nm;

              if (m >= 0)
                {
                  Y_nm = Ynm[pidx];
                  dY_nm = dYnm[pidx];
                }
              else
                {
                  Y_nm = conj(Ynm[pidx]);
                  dY_nm = conj(dYnm[pidx]);
                }

              /* compute (theta,phi) components of J */
              Jt -= rterm * I * m * invsint * qnmj * Y_nm;
              Jp += rterm * qnmj * dY_nm;
            }
        }

      rterm *= ratio;
    }

  /* convert to A/m^2 */
  J[1] += creal(Jt) / w->R / POLTOR_MU_0 * 1.0e-3;
  J[2] += creal(Jp) / w->R / POLTOR_MU_0 * 1.0e-3;

  return s;
} /* poltor_eval_J_shell */

/*
poltor_eval_K_tor()
  Calculate toroidal sheet current vector at a given point
on the spherical shell r = b

Inputs: b     - radius of spherical shell (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        K     - (output) sheet current vector in units of kA/km
                (K_r,K_theta,K_phi) - K_r = 0 for toroidal currents
        w     - workspace

Notes:
1) K should be initialized prior to calling this function, since
on output, K <- K + K_tor
*/

int
poltor_eval_K_tor(const double b, const double theta, const double phi,
                  double K[3], poltor_workspace *w)
{
  int s = 0;
#if 1
  const double ratio = w->R / b;
  double rterm = ratio * ratio * ratio; /* (R/b)^{n+3} */
  green_complex_workspace *green_p = w->green_p[0];
  complex double *Ynm = green_p->Ynm;
  complex double *dYnm = green_p->dYnm;
  complex double Kt = 0.0, Kp = 0.0;
  size_t n;

  /* compute Ynm and d/dtheta Ynm */
  green_complex_Ynm_deriv(theta, phi, green_p);

  for (n = 1; n <= w->nmax_int; ++n)
    {
      double nfac = -(2.0 * n + 1.0) / (double) n;
      int M = (int) GSL_MIN(n, w->mmax_int);
      int m;

      /* (R/b)^{n+3} */
      rterm *= ratio;

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          size_t nmidx = poltor_nmidx(POLTOR_IDX_PINT, n, m, w);
          gsl_complex gnm = poltor_get(nmidx, w);
          complex double qnm = nfac * rterm * (GSL_REAL(gnm) + I * GSL_IMAG(gnm));
          complex double Y_nm, dY_nm;

          if (m >= 0)
            {
              Y_nm = Ynm[pidx];
              dY_nm = dYnm[pidx];
            }
          else
            {
              Y_nm = conj(Ynm[pidx]);
              dY_nm = conj(dYnm[pidx]);
            }

          /* compute (theta,phi) components of K */
          Kt -= m * qnm * Y_nm;
          Kp += qnm * dY_nm;
        }
    }

  Kt *= I / sin(theta);

  /* sanity checks */
  assert(fabs(cimag(Kt)) < 1.0e-8);
  assert(fabs(cimag(Kp)) < 1.0e-8);

  /* convert to kA/km */
  K[1] += (b / w->R) * creal(Kt) / POLTOR_MU_0;
  K[2] += (b / w->R) * creal(Kp) / POLTOR_MU_0;
#else

  size_t n, m;
  const double invsint = 1.0 / sin(theta);
  complex double Jt = 0.0, Jp = 0.0;

  /* compute legendre functions */
  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, w->nmax_int,
                                  cos(theta), w->Pnm, w->dPnm);

  /* pre-compute Ynm and dYnm */
  for (m = 0; m <= w->mmax_int; ++m)
    {
      complex double expimphi = cos(m * phi) + I * sin(m * phi);

      for (n = GSL_MAX(m, 1); n <= w->nmax_int; ++n)
        {
          size_t pidx = gsl_sf_legendre_array_index(n, m);

          w->Ynm[pidx] = w->Pnm[pidx] * expimphi;
          w->dYnm[pidx] = w->dPnm[pidx] * expimphi;
        }
    }

  for (n = 1; n <= w->nmax_int; ++n)
    {
      int M = (int) GSL_MIN(n, w->mmax_int);
      int m;

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          size_t nmidx = poltor_nmidx(POLTOR_IDX_PINT, n, m, w);
          gsl_complex cnm = poltor_get(nmidx, w);
          double gnm = GSL_REAL(cnm);
          double qnm = -(2.0*n + 1.0) / (double)n * gnm *
                       pow(w->R / b, n + 3.0);
          complex double Ynm, dYnm;

          if (m >= 0)
            {
              Ynm = w->Ynm[pidx];
              dYnm = w->dYnm[pidx];
            }
          else
            {
              Ynm = conj(w->Ynm[pidx]);
              dYnm = conj(w->dYnm[pidx]);
            }

          /* compute (theta,phi) components of J */
          Jt -= I * m * invsint * qnm * Ynm;
          Jp += qnm * dYnm;
        }
    }

  /* convert to kA/km */
  K[1] += (b / w->R) * creal(Jt) / POLTOR_MU_0;
  K[2] += (b / w->R) * creal(Jp) / POLTOR_MU_0;
#endif

  return s;
} /* poltor_eval_K_tor */

/*
poltor_eval_K_ext()
  Calculate toroidal sheet current vector at a given point
on the spherical shell r = b due to external sources (ie:
use external coefficients in calculation)

Inputs: r     - radius of shell b (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        K     - (output) sheet current vector in units of kA/km
                (K_r,K_theta,K_phi) - K_r = 0 for toroidal currents
        w     - workspace

Notes:
1) K should be initialized prior to calling this function, since
on output, K <- K + K_ext
*/

int
poltor_eval_K_ext(const double r, const double theta, const double phi,
                  double J[3], poltor_workspace *w)
{
  int s = 0;
#if 0/*XXX*/
  size_t n, m;
  const size_t nmax = w->nmax_ext;
  const size_t mmax = w->mmax_ext;
  const double invsint = 1.0 / sin(theta);
  complex double Jt = 0.0, Jp = 0.0;
  const double ratio = r / w->R;
  double rterm = pow(ratio, -2.0);

  /* compute legendre functions */
  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, nmax,
                                  cos(theta), w->Pnm, w->dPnm);

  /* pre-compute Ynm and dYnm */
  for (m = 0; m <= mmax; ++m)
    {
      complex double expimphi = cos(m * phi) + I * sin(m * phi);

      for (n = GSL_MAX(m, 1); n <= nmax; ++n)
        {
          size_t pidx = gsl_sf_legendre_array_index(n, m);

          w->Ynm[pidx] = w->Pnm[pidx] * expimphi;
          w->dYnm[pidx] = w->dPnm[pidx] * expimphi;
        }
    }

  for (n = 1; n <= nmax; ++n)
    {
      int M = (int) GSL_MIN(n, mmax);
      int m;
      double nfac = (2.0 * n + 1.0) / (n + 1.0);

      rterm *= ratio;

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t nmidx = poltor_nmidx(POLTOR_IDX_PEXT, n, m, w);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          gsl_complex cnm = gsl_vector_complex_get(w->c, nmidx);
          complex double qnm;
          complex double Ynm, dYnm;

          qnm = nfac * rterm * (GSL_REAL(cnm) + I * GSL_IMAG(cnm));

          if (m >= 0)
            {
              Ynm = w->Ynm[pidx];
              dYnm = w->dYnm[pidx];
            }
          else
            {
              Ynm = conj(w->Ynm[pidx]);
              dYnm = conj(w->dYnm[pidx]);
            }

          /* compute (theta,phi) components of J */
          Jt -= I * m * invsint * qnm * Ynm;
          Jp += qnm * dYnm;
        }
    }

  /* convert to kA/km */
  J[1] += ratio * creal(Jt) / POLTOR_MU_0;
  J[2] += ratio * creal(Jp) / POLTOR_MU_0;
#endif

  return s;
} /* poltor_eval_K_ext */

/*
poltor_eval_chi_int()
  Calculate toroidal current scalar function chi at given point
on the spherical shell r = b

chi = -Q b / mu_0 

Inputs: b     - radius of spherical shell (km)
        theta - colatitude (radians)
        phi   - longitude (radians)
        chi   - (output) current scalar in units of kA
        w     - workspace

Return: success/error
*/

int
poltor_eval_chi_int(const double b, const double theta, const double phi,
                    double *chi, poltor_workspace *w)
{
  int s = 0;
  complex double sum = 0.0;
  const double ratio = w->R / b;
  double rterm = ratio * ratio * ratio; /* (R/b)^{n+3} */
  green_complex_workspace *green_p = w->green_p[0];
  complex double *Ynm = green_p->Ynm;
  size_t n;
  double alpha = poltor_sflux_factor(120.0, w); /*XXX*/

  /* compute Ynm */
  green_complex_Ynm(theta, phi, green_p);

  for (n = 1; n <= w->nmax_int; ++n)
    {
      int M = (int) GSL_MIN(n, w->mmax_int);
      int m;
      double nfac = -(2.0 * n + 1.0) / (double) n;

      /* (R/b)^{n+3} */
      rterm *= ratio;

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t nmidx = poltor_nmidx(POLTOR_IDX_PINT, n, m, w);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          gsl_complex gnm = gsl_vector_complex_get(w->c, nmidx);
          complex double qnm = nfac * rterm * (GSL_REAL(gnm) + I * GSL_IMAG(gnm));
          complex double Y_nm;

          if (m >= 0)
            Y_nm = Ynm[pidx];
          else
            Y_nm = conj(Ynm[pidx]);

          sum += qnm * Y_nm;
        }
    }

  /* current function = -Qr/mu0 (units of kA) */
  *chi = -alpha * (b / POLTOR_MU_0) * (b / w->R) * creal(sum);

  return s;
} /* poltor_eval_chi_int() */

/*
poltor_eval_chi_sh()
  Calculate toroidal current scalar function chi at given point
on the spherical shell at radius r

chi = -1 / mu_0 * (r/R) * sum_{nm} q_{nm}(r) Y_{nm}

Inputs: r     - radius (km) in [rmin,rmax]
        theta - colatitude (radians)
        phi   - longitude (radians)
        chi   - (output) current scalar in units of kA
        w     - workspace

Return: success/error
*/

int
poltor_eval_chi_sh(const double r, const double theta, const double phi,
                   double *chi, poltor_workspace *w)
{
  int s = 0;
  complex double sum = 0.0;
  const double ratio = (r - w->rmid) / w->R;
  double rterm = 1.0;
  green_complex_workspace *green_p = w->green_p[0];
  complex double *Ynm = green_p->Ynm;
  size_t n, j;
  double alpha = poltor_sflux_factor(120.0, w); /*XXX*/

  /* compute Ynm */
  green_complex_Ynm(theta, phi, green_p);

  for (j = 0; j <= w->shell_J; ++j)
    {
      for (n = 1; n <= w->nmax_sh; ++n)
        {
          int M = (int) GSL_MIN(n, w->mmax_sh);
          int m;

          for (m = -M; m <= M; ++m)
            {
              int mabs = abs(m);
              size_t nmjidx = poltor_jnmidx(j, n, m, w);
              size_t pidx = gsl_sf_legendre_array_index(n, mabs);
              gsl_complex cnmj = gsl_vector_complex_get(w->c, nmjidx);
              complex double qnmj = GSL_REAL(cnmj) + I * GSL_IMAG(cnmj);
              complex double Y_nm;

              if (m >= 0)
                Y_nm = Ynm[pidx];
              else
                Y_nm = conj(Ynm[pidx]);

              sum += rterm * qnmj * Y_nm;
            }
        }

      /* ((r-rmid)/R)^j */
      rterm *= ratio;
    }

  /* current function = -Qr/mu0 (units of kA) */
  *chi = -alpha / POLTOR_MU_0 * (r / w->R) * creal(sum);

  return s;
} /* poltor_eval_chi_sh() */

/*
poltor_eval_chi_ext()
  Calculate toroidal current scalar function chi at given point
on the spherical shell r = b

chi = -Q r / mu_0 

Inputs: b          - radius of current shell (km)
        theta      - colatitude (radians)
        phi        - longitude (radians)
        chi        - (output) current scalar in units of kA
        w          - workspace

Return: success/error
*/

int
poltor_eval_chi_ext(const double b, const double theta, const double phi,
                    double *chi, poltor_workspace *w)
{
  int s = 0;
#if 0/*XXX*/
  size_t n, m;
  complex double sum = 0.0;
  const size_t nmax = w->nmax_ext;
  const size_t mmax = w->mmax_ext;
  const double ratio = b / w->R;
  double rterm = pow(ratio, -2.0);

  /* compute legendre functions */
  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, nmax,
                                  cos(theta), w->Pnm, w->dPnm);

  /* pre-compute Ynm and dYnm */
  for (m = 0; m <= mmax; ++m)
    {
      complex double expimphi = cos(m * phi) + I * sin(m * phi);

      for (n = GSL_MAX(m, 1); n <= nmax; ++n)
        {
          size_t pidx = gsl_sf_legendre_array_index(n, m);

          w->Ynm[pidx] = w->Pnm[pidx] * expimphi;
          w->dYnm[pidx] = w->dPnm[pidx] * expimphi;
        }
    }

  for (n = 1; n <= nmax; ++n)
    {
      int M = (int) GSL_MIN(n, mmax);
      int m;
      double nfac = (2.0 * n + 1.0) / (n + 1.0);

      /* (b/R)^(n-1) */
      rterm *= ratio;

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t nmidx = poltor_nmidx(POLTOR_IDX_PEXT, n, m, w);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          gsl_complex cnm = gsl_vector_complex_get(w->c, nmidx);
          complex double qnm;
          complex double Ynm;

          qnm = nfac * rterm * (GSL_REAL(cnm) + I * GSL_IMAG(cnm));

          if (m >= 0)
            Ynm = w->Ynm[pidx];
          else
            Ynm = conj(w->Ynm[pidx]);

          sum += qnm * Ynm;
        }
    }

  /* current function = -Qr/mu0 (units of kA) */
  sum = -sum / POLTOR_MU_0 * b;
  *chi = creal(sum);
#endif

  return s;
} /* poltor_eval_chi_ext() */

/*
poltor_eval_B()
  Calculate total magnetic field model at a given point

Inputs: t     - timestamp (CDF_EPOCH)
        r     - geocentric radius (km)
        theta - geocentric colatitude (radians)
        phi   - longitude (radians)
        B     - (output) NEC magnetic field model in units of nT
                (B_x,B_y,B_z,|B|)
        w     - workspace

Return: success/error
*/

int
poltor_eval_B(const double t, const double r, const double theta, const double phi,
              double B[4], poltor_workspace *w)
{
  int s;
  double B_int[4], B_ext[4], B_sh[4], B_tor[4];
  size_t k;

  s = poltor_eval_B_all(t, r, theta, phi, B_int, B_ext, B_sh, B_tor, w);
  if (s)
    return s;

  for (k = 0; k < 3; ++k)
    B[k] = B_int[k] + B_ext[k] + B_sh[k] + B_tor[k];

  B[3] = gsl_hypot3(B[0], B[1], B[2]);

  return 0;
}

/*
poltor_eval_B_all()
  Calculate magnetic field model at a given point for all sources

Inputs: t     - timestamp (CDF_EPOCH)
        r     - geocentric radius (km)
        theta - geocentric colatitude (radians)
        phi   - longitude (radians)
        B_int - (output) internal poloidal NEC magnetic field model in units of nT
                (B_x,B_y,B_z,|B|)
        B_ext - (output) external poloidal NEC magnetic field model in units of nT
                (B_x,B_y,B_z,|B|)
        B_sh  - (output) shell poloidal NEC magnetic field model in units of nT
                (B_x,B_y,B_z,|B|)
        B_tor - (output) shell toroidal NEC magnetic field model in units of nT
                (B_x,B_y,B_z,|B|)
        w     - workspace

Return: success/error
*/

int
poltor_eval_B_all(const double t, const double r, const double theta, const double phi,
                  double B_int[4], double B_ext[4], double B_sh[4], double B_tor[4],
                  poltor_workspace *w)
{
  int s = 0;
  green_complex_workspace *green_p = w->green_p[0];
  complex double *Ynm = green_p->Ynm;
  complex double *dYnm = green_p->dYnm;
  gsl_vector_complex_view dX = gsl_matrix_complex_row(w->omp_dX, 0);
  gsl_vector_complex_view dY = gsl_matrix_complex_row(w->omp_dY, 0);
  gsl_vector_complex_view dZ = gsl_matrix_complex_row(w->omp_dZ, 0);
  gsl_vector_complex_view X, Y, Z;
  gsl_vector_complex_view c;
  gsl_complex val[3];
  time_t unix_time = satdata_epoch2timet(t);
  double euvac, alpha; /* 1 + N*EUVAC */
  size_t k;

  f107_get_euvac(unix_time, &euvac, w->f107_workspace_p);
  alpha = poltor_sflux_factor(euvac, w);

  /* initialize */
  for (k = 0; k < 4; ++k)
    {
      B_int[k] = 0.0;
      B_ext[k] = 0.0;
      B_sh[k] = 0.0;
      B_tor[k] = 0.0;
    }

  /* compute Ynm and d/dtheta Ynm */
  green_complex_Ynm_deriv(theta, phi, green_p);

  if (w->p_pint > 0)
    {
      c = gsl_vector_complex_subvector(w->c, w->pint_offset, w->p_pint);

      X = gsl_vector_complex_subvector(&dX.vector, w->pint_offset, w->p_pint);
      Y = gsl_vector_complex_subvector(&dY.vector, w->pint_offset, w->p_pint);
      Z = gsl_vector_complex_subvector(&dZ.vector, w->pint_offset, w->p_pint);

      /* compute internal Green's functions */
      green_complex_calc_int(w->nmax_int, w->mmax_int, w->R, r, theta, Ynm, dYnm,
                             (complex double *) X.vector.data,
                             (complex double *) Y.vector.data,
                             (complex double *) Z.vector.data);

      gsl_blas_zdotu(&X.vector, &c.vector, &val[0]);
      gsl_blas_zdotu(&Y.vector, &c.vector, &val[1]);
      gsl_blas_zdotu(&Z.vector, &c.vector, &val[2]);

      B_int[0] = GSL_REAL(val[0]);
      B_int[1] = GSL_REAL(val[1]);
      B_int[2] = GSL_REAL(val[2]);
    }

  if (w->p_pext > 0)
    {
      c = gsl_vector_complex_subvector(w->c, w->pext_offset, w->p_pext);

      X = gsl_vector_complex_subvector(&dX.vector, w->pext_offset, w->p_pext);
      Y = gsl_vector_complex_subvector(&dY.vector, w->pext_offset, w->p_pext);
      Z = gsl_vector_complex_subvector(&dZ.vector, w->pext_offset, w->p_pext);

      /* compute external Green's functions */
      green_complex_calc_ext(w->nmax_ext, w->mmax_ext, w->R, r, theta, Ynm, dYnm,
                             (complex double *) X.vector.data,
                             (complex double *) Y.vector.data,
                             (complex double *) Z.vector.data);

      gsl_blas_zdotu(&X.vector, &c.vector, &val[0]);
      gsl_blas_zdotu(&Y.vector, &c.vector, &val[1]);
      gsl_blas_zdotu(&Z.vector, &c.vector, &val[2]);

      B_ext[0] = GSL_REAL(val[0]);
      B_ext[1] = GSL_REAL(val[1]);
      B_ext[2] = GSL_REAL(val[2]);
    }

  if (w->p_psh > 0)
    {
      c = gsl_vector_complex_subvector(w->c, w->psh_offset, w->p_psh);

      X = gsl_vector_complex_subvector(&dX.vector, w->psh_offset, w->p_psh);
      Y = gsl_vector_complex_subvector(&dY.vector, w->psh_offset, w->p_psh);
      Z = gsl_vector_complex_subvector(&dZ.vector, w->psh_offset, w->p_psh);

      /* compute shell Green's functions */
      poltor_calc_psh(r, theta, Ynm, dYnm,
                      (complex double *) X.vector.data,
                      (complex double *) Y.vector.data,
                      (complex double *) Z.vector.data,
                      w);

      gsl_blas_zdotu(&X.vector, &c.vector, &val[0]);
      gsl_blas_zdotu(&Y.vector, &c.vector, &val[1]);
      gsl_blas_zdotu(&Z.vector, &c.vector, &val[2]);

      B_sh[0] = GSL_REAL(val[0]);
      B_sh[1] = GSL_REAL(val[1]);
      B_sh[2] = GSL_REAL(val[2]);
    }

  if (w->p_tor > 0)
    {
      c = gsl_vector_complex_subvector(w->c, w->tor_offset, w->p_tor);

      X = gsl_vector_complex_subvector(&dX.vector, w->tor_offset, w->p_tor);
      Y = gsl_vector_complex_subvector(&dY.vector, w->tor_offset, w->p_tor);
      Z = gsl_vector_complex_subvector(&dZ.vector, w->tor_offset, w->p_tor);

      /* compute shell Green's functions */
      poltor_calc_tor(r, theta, Ynm, dYnm,
                      (complex double *) X.vector.data,
                      (complex double *) Y.vector.data,
                      (complex double *) Z.vector.data,
                      w);

      gsl_blas_zdotu(&X.vector, &c.vector, &val[0]);
      gsl_blas_zdotu(&Y.vector, &c.vector, &val[1]);
      gsl_blas_zdotu(&Z.vector, &c.vector, &val[2]);

      B_tor[0] = GSL_REAL(val[0]);
      B_tor[1] = GSL_REAL(val[1]);
      B_tor[2] = GSL_REAL(val[2]);
    }

  /* multiply by time dependent factor and compute magnitudes */
  for (k = 0; k < 3; ++k)
    {
      B_int[k] *= alpha;
      B_ext[k] *= alpha;
      B_sh[k] *= alpha;
      B_tor[k] *= alpha;
    }

  B_int[3] = gsl_hypot3(B_int[0], B_int[1], B_int[2]);
  B_ext[3] = gsl_hypot3(B_ext[0], B_ext[1], B_ext[2]);
  B_sh[3] = gsl_hypot3(B_sh[0], B_sh[1], B_sh[2]);
  B_tor[3] = gsl_hypot3(B_tor[0], B_tor[1], B_tor[2]);

  return s;
} /* poltor_eval_B_all() */

int
poltor_write(const char *filename, poltor_workspace *w)
{
  int s = 0;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "poltor_write: unable to open %s: %s\n",
              filename, strerror(errno));
      return GSL_FAILURE;
    }

  fwrite(&(w->params), sizeof(poltor_parameters), 1, fp);
  gsl_vector_complex_fwrite(fp, w->c);

  fclose(fp);

  return s;
} /* poltor_write() */

poltor_workspace *
poltor_read(const char *filename)
{
  poltor_workspace *w;
  poltor_parameters params;
  FILE *fp;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "poltor_read: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  fread(&params, sizeof(poltor_parameters), 1, fp);

  params.data = NULL;

  w = poltor_alloc(&params);

  gsl_vector_complex_fread(fp, w->c);

  fclose(fp);

  return w;
} /* poltor_read() */

/*
poltor_get()
  Return coefficient corresponding to a given index.
The index should be obtained from poltor_nmidx() and so
it will already account for internal/external etc

Inputs: cidx - coefficient index from poltor_nmidx()
        w    - workspace
*/

gsl_complex
poltor_get(const size_t cidx, const poltor_workspace *w)
{
  gsl_complex coeff;
  
  if (cidx >= w->c->size)
    GSL_SET_COMPLEX(&coeff, 0.0, 0.0);
  else
    coeff = gsl_vector_complex_get(w->c, cidx);

  return coeff;
}

double
poltor_spectrum_int(const size_t n, gsl_vector_complex *c, poltor_workspace *w)
{
  int ni = (int) GSL_MIN(n, w->mmax_int);
  int m;
  double sum = 0.0;

  if (n > w->nmax_int)
    return 0.0;

  for (m = -ni; m <= ni; ++m)
    {
      size_t cidx = poltor_nmidx(POLTOR_IDX_PINT, n, m, w);
      gsl_complex coef = gsl_vector_complex_get(c, cidx);
      double gnm = GSL_REAL(coef);

      sum += gnm * gnm;
    }

  /* see Backus (4.4.22) */
  sum *= (n + 1.0);

  return sum;
}

double
poltor_spectrum_ext(const size_t n, gsl_vector_complex *c, poltor_workspace *w)
{
  int ni = (int) GSL_MIN(n, w->mmax_ext);
  int m;
  double sum = 0.0;

  if (n > w->nmax_ext)
    return 0.0;

  for (m = -ni; m <= ni; ++m)
    {
      size_t cidx = poltor_nmidx(POLTOR_IDX_PEXT, n, m, w);
      gsl_complex coef = gsl_vector_complex_get(c, cidx);
      double knm = GSL_REAL(coef);

      sum += knm * knm;
    }

  /* see Backus (4.4.22) */
  sum *= (n + 1.0);

  return sum;
} /* poltor_spectrum_ext() */

double
poltor_spectrum_sh(const size_t n, gsl_vector_complex *c, poltor_workspace *w)
{
  int ni = (int) GSL_MIN(n, w->mmax_sh);
  int m;
  double sum = 0.0;

  if (n > w->nmax_sh)
    return 0.0;

  for (m = -ni; m <= ni; ++m)
    {
      size_t cidx = poltor_nmidx(POLTOR_IDX_PSH, n, m, w);
      gsl_complex coef = gsl_vector_complex_get(c, cidx);
      double qnm = GSL_REAL(coef);

      sum += qnm * qnm;
    }

  /* see Backus (4.4.22) */
  sum *= (n + 1.0);

  return sum;
} /* poltor_spectrum_sh() */

double
poltor_spectrum_tor(const size_t n, gsl_vector_complex *c, poltor_workspace *w)
{
  int ni = (int) GSL_MIN(n, w->mmax_tor);
  int m;
  double sum = 0.0;

  if (n > w->nmax_tor)
    return 0.0;

  for (m = -ni; m <= ni; ++m)
    {
      size_t cidx = poltor_nmidx(POLTOR_IDX_TOR, n, m, w);
      gsl_complex coef = gsl_vector_complex_get(c, cidx);
      double phinm = GSL_REAL(coef);

      sum += phinm * phinm;
    }

  /* see Backus (4.4.22) */
  sum *= (n + 1.0);

  return sum;
} /* poltor_spectrum_tor() */

int
poltor_print_spectrum(const char *filename, gsl_vector_complex *c, poltor_workspace *w)
{
  const size_t nmax = w->nmax_max;
  size_t n;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "poltor_print_spectrum: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  n = 1;
  fprintf(fp, "# Field %zu: spherical harmonic degree\n", n++);
  fprintf(fp, "# Field %zu: internal spectrum g_n (nT^2)\n", n++);
  fprintf(fp, "# Field %zu: external spectrum k_n (nT^2)\n", n++);
  fprintf(fp, "# Field %zu: shell poloidal spectrum q_n (nT^2)\n", n++);
  fprintf(fp, "# Field %zu: shell toroidal spectrum phi_n (nT^2)\n", n++);

  for (n = 1; n <= nmax; ++n)
    {
      double gn = poltor_spectrum_int(n, c, w);
      double kn = poltor_spectrum_ext(n, c, w);
      double qn = poltor_spectrum_sh(n, c, w);
      double phin = poltor_spectrum_tor(n, c, w);

      fprintf(fp, "%zu %.12e %.12e %.12e %.12e\n",
              n,
              gn,
              kn,
              qn,
              phin);
    }

  fclose(fp);

  return 0;
} /* poltor_print_spectrum() */

int
poltor_print_lcurve(const char *filename, const poltor_workspace *w)
{
  FILE *fp;
  size_t i;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "poltor_print_lcurve: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Field %zu: regularization parameter lambda\n", i++);
  fprintf(fp, "# Field %zu: residual norm ||b - A x||\n", i++);
  fprintf(fp, "# Field %zu: solution norm ||L x||\n", i++);

  for (i = 0; i < w->nreg; ++i)
    {
      fprintf(fp, "%.12e %.12e %.12e\n",
              gsl_vector_get(w->reg_param, i),
              gsl_vector_get(w->rho, i),
              gsl_vector_get(w->eta, i));
    }

  fprintf(fp, "\n\n");

  fprintf(fp, "%.12e %.12e %.12e\n",
          gsl_vector_get(w->reg_param, w->reg_idx),
          gsl_vector_get(w->rho, w->reg_idx),
          gsl_vector_get(w->eta, w->reg_idx));

  fclose(fp);

  return 0;
}

/*
poltor_calc_psh()
  Compute Green's functions for X,Y,Z spherical harmonic expansion due to
toroidal current source inside the shell. These are simply the basis functions multiplying
the q_{nm}^{(j)} coefficients. The Y_{nm} and d/dtheta Y_{nm} functions must be
precomputed and passed in

Inputs: r     - radius (km)
        theta - colatitude (radians)
        Ynm   - P_{nm} exp(i m phi), size nnm
        dYnm  - d/dtheta P_{nm} exp(i m phi), size nnm
        X     - (output) array of X Green's functions, size p_psh
        Y     - (output) array of Y Green's functions, size p_psh
        Z     - (output) array of Z Green's functions, size p_psh
        w     - workspace
*/

int
poltor_calc_psh(const double r, const double theta,
                const complex double *Ynm, const complex double *dYnm,
                complex double *X, complex double *Y, complex double *Z,
                poltor_workspace *w)
{
  int s = 0;
  size_t n, j;
  const double ratio = w->R / r;
  const complex double invisint = I / sin(theta);

  /* outer loop over Taylor series */
  for (j = 0; j <= w->shell_J; ++j)
    {
      for (n = 1; n <= w->nmax_sh; ++n)
        {
          int M = (int) GSL_MIN(n, w->mmax_sh);
          int m;
          double An, Bn, ABsum, dABsum;

          /* compute Taylor series integrals */
          s = poltor_shell_An(n, j, r, &An, w);
          s += poltor_shell_Bn(n, j, r, &Bn, w);
          if (s)
            return s;

          ABsum = An + Bn;
          dABsum = -(double)n * An + (n + 1.0) * Bn;

          for (m = -M; m <= M; ++m)
            {
              int mabs = abs(m);
              size_t nmidx = poltor_jnmidx2(j, n, m, w);
              size_t pidx = gsl_sf_legendre_array_index(n, mabs);
              complex double Y_nm, dY_nm;

              if (m >= 0)
                {
                  Y_nm = Ynm[pidx];
                  dY_nm = dYnm[pidx];
                }
              else
                {
                  Y_nm = conj(Ynm[pidx]);
                  dY_nm = conj(dYnm[pidx]);
                }

              /* compute (X,Y,Z) components of B_{lm} */

              X[nmidx] = ratio * dABsum * dY_nm;
              Y[nmidx] = -ratio * m * invisint * dABsum * Y_nm;
              Z[nmidx] = ratio * n * (n + 1.0) * ABsum * Y_nm;
            }
        }
    }

  return s;
}

/*
poltor_calc_tor()
  Compute Green's functions for X,Y,Z spherical harmonic expansion due to
poloidal current source inside the shell (toroidal field). These are simply the
basis functions multiplying the phi{nm} coefficients.
The Y_{nm} and d/dtheta Y_{nm} functions must be precomputed and passed in.

Inputs: r     - radius (km)
        theta - colatitude (radians)
        Ynm   - P_{nm} exp(i m phi), size nnm
        dYnm  - d/dtheta P_{nm} exp(i m phi), size nnm
        X     - (output) array of X Green's functions, size p_tor
        Y     - (output) array of Y Green's functions, size p_tor
        Z     - (output) array of Z Green's functions, size p_tor
        w     - workspace

Notes:
1) This routine assumes phi_{nm}(r) is constant
*/

int
poltor_calc_tor(const double r, const double theta,
                const complex double *Ynm, const complex double *dYnm,
                complex double *X, complex double *Y, complex double *Z,
                poltor_workspace *w)
{
  int s = 0;
  const complex double invisint = I / sin(theta);
  size_t n;

  for (n = 1; n <= w->nmax_tor; ++n)
    {
      int M = (int) GSL_MIN(n, w->mmax_tor);
      int m;

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t nmidx = green_idx(n, m, w->mmax_tor);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          complex double Y_nm, dY_nm;

          if (m >= 0)
            {
              Y_nm = Ynm[pidx];
              dY_nm = dYnm[pidx];
            }
          else
            {
              Y_nm = conj(Ynm[pidx]);
              dY_nm = conj(dYnm[pidx]);
            }

          /* compute (r,theta,phi) components of B_{lm} */

          X[nmidx] = m * invisint * Y_nm;
          Y[nmidx] = dY_nm;
          Z[nmidx] = 0.0;
        }
    }

  return s;
}

double
poltor_sflux_factor(const double flux, const poltor_workspace *w)
{
  if (w->params.synth_data || !w->params.solar_flux_factor)
    {
      return 1.0;
    }
  else
    {
      const double N_wolf = 14.85e-3; /* Wolf ratio in sfu (Chulliat et al 2016) */
      return (1.0 + N_wolf * flux);
    }
}

/***********************************************
 * INTERNAL ROUTINES                           *
 ***********************************************/

#if 0/*XXX*/
/*
poltor_build_ls()
  Build LS design matrix

Inputs: fold - if set to 1, LS matrix and rhs vector are folded
               into A^H W A and A^H W b
               if set to 0, it is assumed we already solved the
               LS system, and now wish to compute residual norm
               ||b - A x|| for L-curve analysis
        w    - workspace

Notes:
1) If fold = 1, on output, w->lls_workspace_p is updated with A^H W A and A^H W b,
ready for solving

2) If fold = 0, on output, w->chisq contains ||sqrt(W) (b - A x)||^2

3) If fold = 0, on output, w->residuals contains sqrt(W) * (b - A x)

4) Flagged outliers are not added to matrix
*/

int
poltor_build_ls(const int fold, poltor_workspace *w)
{
  int s = 0;
  magdata *data = w->data;
  size_t i;
  struct timeval tv0, tv1;
  size_t data_idx = 0; /* index into 'data' */
  size_t res_idx = 0;  /* index into 'residuals' */
  int nleft;           /* number of vector data left to fold in */

  fprintf(stderr, "\n");

  w->chisq = 0.0;

  nleft = (int) data->n;
  while (nleft > 0)
    {
      gsl_matrix_complex_view A;
      gsl_vector_complex_view rhs;
      gsl_vector_view wts;
      size_t ridx = 0; /* index of current row in A */

      /* number of vector measurements to process */
      size_t ndata = (size_t) GSL_MIN(POLTOR_BLOCK_SIZE, nleft);

      fprintf(stderr, "\t adding measurements [%zu,%zu] of %zu into matrix...",
              data_idx, data_idx + ndata - 1, data->n);
      gettimeofday(&tv0, NULL);

      /* initialize matrix to 0 */
      gsl_matrix_complex_set_zero(w->A);

      for (i = 0; i < ndata; ++i, ++data_idx)
        {
          double theta = poltor_theta(data_idx, w);
          double theta_ns = poltor_theta_ns(data_idx, w);
          gsl_vector_complex *x = NULL,
                             *y = NULL,
                             *z = NULL,
                             *dx = NULL,
                             *dy = NULL,
                             *dz = NULL;
          gsl_vector_complex_view vx, vy, vz, vdx, vdy, vdz;
          gsl_complex val;
          double B[4], dB[4];

          if (MAGDATA_Discarded(data->flags[i]))
            continue;

          /* compute residual: B_i - B_main - B_crust - B_ext */
          magdata_residual(data_idx, B, data);

          if (!gsl_finite(B[0]) || !gsl_finite(B[1]) || !gsl_finite(B[2]))
            {
              fprintf(stderr, "ACK1\n");
            }

          if (data->flags[data_idx] & (MAGDATA_FLG_DX_NS|MAGDATA_FLG_DY_NS|MAGDATA_FLG_DZ_NS))
            {
              /* compute dB residual: B_i_ns - B_i */
              magdata_residual_dB_ns(data_idx, dB, data);
              if (!gsl_finite(dB[0]) || !gsl_finite(dB[1]) || !gsl_finite(dB[2]))
                {
                  fprintf(stderr, "ACK2\n");
                }
            }

          if (data->flags[data_idx] & MAGDATA_FLG_X)
            {
              /* weights vector */
              gsl_vector_set(w->weights, ridx, POLTOR_WEIGHT_X * data->weights[data_idx]);

              /* rhs vector */
              GSL_SET_COMPLEX(&val, B[0], 0.0);
              gsl_vector_complex_set(w->rhs, ridx, val);

              vx = gsl_matrix_complex_row(w->A, ridx++);
              x = &vx.vector;
            }

          if (data->flags[data_idx] & MAGDATA_FLG_DX_NS)
            {
              gsl_vector_set(w->weights, ridx, POLTOR_WEIGHT_DX * data->weights[data_idx]);

              GSL_SET_COMPLEX(&val, dB[0], 0.0);
              gsl_vector_complex_set(w->rhs, ridx, val);

              vdx = gsl_matrix_complex_row(w->A, ridx++);
              dx = &vdx.vector;
            }

          if (data->flags[data_idx] & MAGDATA_FLG_Y)
            {
              /* weights vector */
              gsl_vector_set(w->weights, ridx, POLTOR_WEIGHT_Y * data->weights[data_idx]);

              /* rhs vector */
              GSL_SET_COMPLEX(&val, B[1], 0.0);
              gsl_vector_complex_set(w->rhs, ridx, val);

              vy = gsl_matrix_complex_row(w->A, ridx++);
              y = &vy.vector;
            }

          if (data->flags[data_idx] & MAGDATA_FLG_DY_NS)
            {
              gsl_vector_set(w->weights, ridx, POLTOR_WEIGHT_DY * data->weights[data_idx]);

              GSL_SET_COMPLEX(&val, dB[1], 0.0);
              gsl_vector_complex_set(w->rhs, ridx, val);

              vdy = gsl_matrix_complex_row(w->A, ridx++);
              dy = &vdy.vector;
            }

          if (data->flags[data_idx] & MAGDATA_FLG_Z)
            {
              /* weights vector */
              gsl_vector_set(w->weights, ridx, POLTOR_WEIGHT_Z * data->weights[data_idx]);

              /* rhs vector */
              GSL_SET_COMPLEX(&val, B[2], 0.0);
              gsl_vector_complex_set(w->rhs, ridx, val);

              vz = gsl_matrix_complex_row(w->A, ridx++);
              z = &vz.vector;
            }

          if (data->flags[data_idx] & MAGDATA_FLG_DZ_NS)
            {
              gsl_vector_set(w->weights, ridx, POLTOR_WEIGHT_DZ * data->weights[data_idx]);

              GSL_SET_COMPLEX(&val, dB[2], 0.0);
              gsl_vector_complex_set(w->rhs, ridx, val);

              vdz = gsl_matrix_complex_row(w->A, ridx++);
              dz = &vdz.vector;
            }

          /* build these 6 rows of design matrix */
          poltor_row(data->r[data_idx], theta, data->phi[data_idx],
                     data->r_ns[data_idx], theta_ns, data->phi_ns[data_idx],
                     x, y, z, dx, dy, dz, w);

          if ((x && !gsl_finite(gsl_blas_dznrm2(x))) ||
              (y && !gsl_finite(gsl_blas_dznrm2(y))) ||
              (z && !gsl_finite(gsl_blas_dznrm2(z))) ||
              (dx && !gsl_finite(gsl_blas_dznrm2(dx))) ||
              (dy && !gsl_finite(gsl_blas_dznrm2(dy))) ||
              (dz && !gsl_finite(gsl_blas_dznrm2(dz))))
            {
              fprintf(stderr, "UH OH\n");
              exit(1);
            }
        }

      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

      /* add design matrix to least squares system */

      A = gsl_matrix_complex_submatrix(w->A, 0, 0, ridx, w->p);
      rhs = gsl_vector_complex_subvector(w->rhs, 0, ridx);
      wts = gsl_vector_subvector(w->weights, 0, ridx);

      if (fold)
        {
          /* convert LS system to standard form */
          fprintf(stderr, "\t converting (A,b) to standard form...");
          gettimeofday(&tv0, NULL);
          lls_complex_stdform(&A.matrix, &rhs.vector, &wts.vector, w->L, w->lls_workspace_p);
          gettimeofday(&tv1, NULL);
          fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

          /* compute A^H A and A^H b */
          fprintf(stderr, "\t folding into A^H A and A^H b...");
          gettimeofday(&tv0, NULL);
          lls_complex_fold(&A.matrix, &rhs.vector, w->lls_workspace_p);
          gettimeofday(&tv1, NULL);
          fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));
        }
      else
        {
          gsl_vector_complex_view resv = gsl_vector_complex_subvector(w->residuals, res_idx, ridx);
          gsl_complex one, mone;
          double rnorm;

          GSL_SET_COMPLEX(&one, 1.0, 0.0);
          GSL_SET_COMPLEX(&mone, -1.0, 0.0);

          fprintf(stderr, "\t computing ||sqrt(W) (b - Ac)||^2...");

          /* calculate r = b - A c with previously computed coefficients */
          gsl_vector_complex_memcpy(&resv.vector, &rhs.vector);
          gsl_blas_zgemv(CblasNoTrans, mone, &A.matrix, w->c, one, &resv.vector);

          /* compute weighted residual term in chi^2: r_i' = sqrt(w_i) r_i */
          for (i = 0; i < ridx; ++i)
            {
              double wi = gsl_vector_get(&wts.vector, i);
              gsl_complex ri = gsl_vector_complex_get(&resv.vector, i);

              gsl_vector_complex_set(&resv.vector, i, gsl_complex_mul_real(ri, sqrt(wi)));
            }

          /* compute || sqrt(W) (b - A c) || */
          rnorm = gsl_blas_dznrm2(&resv.vector);

          w->chisq += rnorm * rnorm;

          fprintf(stderr, "done (chi^2 = %f)\n", w->chisq);

          res_idx += ridx;
        }

      nleft -= (int) ndata;
    }

  return s;
} /* poltor_build_ls() */

/*
poltor_row()
  Build 6 rows of LS design matrix for X, Y, Z, DX, DY, DZ components

Inputs: r      - geocentric radius (km)
        theta  - geocentric colatitude (radians)
        phi    - longitude (radians)
        r_ns   - geocentric radius for along-track point (km)
        theta_ns - geocentric colatitude for along-track point (radians)
        phi_ns - longitude for along-track point (radians)
        x      - (output) row corresponding to X (or NULL)
        y      - (output) row corresponding to Y (or NULL)
        z      - (output) row corresponding to Z (or NULL)
        dx     - (output) row corresponding to DX (or NULL)
        dy     - (output) row corresponding to DY (or NULL)
        dz     - (output) row corresponding to DZ (or NULL)
        w      - workspace

Return: success/error

Notes: row vectors should be initialized to 0 prior to calling function
*/

static int
poltor_row(const double r, const double theta, const double phi,
           const double r_ns, const double theta_ns, const double phi_ns,
           gsl_vector_complex *x, gsl_vector_complex *y, gsl_vector_complex *z,
           gsl_vector_complex *dx, gsl_vector_complex *dy, gsl_vector_complex *dz,
           poltor_workspace *w)
{
  int s = 0;
  const size_t nmax = w->nmax_max;
  const size_t mmax = w->mmax_max;
  const int grad = dx || dy || dz;
  size_t n, m;

  /* compute legendre functions */
  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, nmax,
                                  cos(theta), w->Pnm, w->dPnm);

  if (grad)
    {
      gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, nmax,
                                      cos(theta_ns), w->Pnm2, w->dPnm2);
    }

  /* pre-compute Ynm and dYnm */
  for (m = 0; m <= mmax; ++m)
    {
      complex double expimphi = cos(m * phi) + I * sin(m * phi);
      complex double expimphi_ns = cos(m * phi_ns) + I * sin(m * phi_ns);

      for (n = GSL_MAX(m, 1); n <= nmax; ++n)
        {
          size_t pidx = gsl_sf_legendre_array_index(n, m);

          w->Ynm[pidx] = w->Pnm[pidx] * expimphi;
          w->dYnm[pidx] = w->dPnm[pidx] * expimphi;

          if (grad)
            {
              w->Ynm2[pidx] = w->Pnm2[pidx] * expimphi_ns;
              w->dYnm2[pidx] = w->dPnm2[pidx] * expimphi_ns;
            }
        }
    }

  /* build row for B_pol^i coefficients */
  if (w->p_pint > 0)
    poltor_row_pint(r, theta, r_ns, theta_ns, x, y, z, dx, dy, dz, w);

  /* build row for B_pol^e coefficients */
  if (w->p_pext > 0)
    poltor_row_pext(r, theta, r_ns, theta_ns, x, y, z, dx, dy, dz, w);

  /* build row for B_pol^{sh} coefficients */
  if (w->p_psh > 0)
    {
#if 0
      double phi_deg = wrap180(phi * 180.0 / M_PI);
      double qdlat = 90.0 - theta * 180.0 / M_PI;

      /* exclude EEJ region from shell poloidal calculation */
      if (!(fabs(qdlat) <= 10.0 && (phi_deg >= -105.0 && phi_deg <= 105.0)))
#endif
        poltor_row_psh(r, theta, r_ns, theta_ns, x, y, z, dx, dy, dz, w);
    }

  /* build row for B_tor coefficients */
  if (w->p_tor > 0)
    poltor_row_tor(r, theta, r_ns, theta_ns, x, y, z, dx, dy, dz, w);

  return s;
} /* poltor_row() */

/*
poltor_row_pint()
  Build a row of the LS design matrix corresponding to poloidal
internal magnetic field B_pol^i

Inputs: r     - radius (km)
        theta - colatitude (radians)
        vx    - (output) where to store matrix row X component (or NULL)
        vy    - (output) where to store matrix row Y component (or NULL)
        vz    - (output) where to store matrix row Z component (or NULL)
        vdx   - (output) where to store matrix row DX component (or NULL)
        vdy   - (output) where to store matrix row DY component (or NULL)
        vdz   - (output) where to store matrix row DZ component (or NULL)
        w     - workspace

Return: success/error

Notes:
1) w->Ynm, w->dYnm must be initialized prior to calling this function
*/

static int
poltor_row_pint(const double r, const double theta,
                const double r_ns, const double theta_ns,
                gsl_vector_complex *vx, gsl_vector_complex *vy, gsl_vector_complex *vz,
                gsl_vector_complex *vdx, gsl_vector_complex *vdy, gsl_vector_complex *vdz,
                poltor_workspace *w)
{
  int s = 0;
  size_t n;
  gsl_complex val;
  const complex double invisint = I / sin(theta);
  const complex double invisint2 = I / sin(theta_ns);
  const double ratio = w->R / r;
  const double ratio2 = w->R / r_ns;
  double rterm = ratio * ratio;
  double rterm2 = ratio2 * ratio2;

  /*
   * construct 3 rows of the design matrix corresponding to the
   * X,Y,Z data observation
   */
  for (n = 1; n <= w->nmax_int; ++n)
    {
      int M = (int) GSL_MIN(n, w->mmax_int);
      int m;

      /* (R/r)^(n+2) */
      rterm *= ratio;

      /* (R/r_ns)^(n+2) */
      rterm2 *= ratio2;

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t nmidx = poltor_nmidx(POLTOR_IDX_PINT, n, m, w);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          complex double X, Y, Z, Ynm, dYnm, Ynm2, dYnm2;

          if (m >= 0)
            {
              Ynm = w->Ynm[pidx];
              dYnm = w->dYnm[pidx];
              Ynm2 = w->Ynm2[pidx];
              dYnm2 = w->dYnm2[pidx];
            }
          else
            {
              Ynm = conj(w->Ynm[pidx]);
              dYnm = conj(w->dYnm[pidx]);
              Ynm2 = conj(w->Ynm2[pidx]);
              dYnm2 = conj(w->dYnm2[pidx]);
            }

          /* compute (X,Y,Z) components of B_{lm} */
          if (vx)
            {
              X = rterm * dYnm;
              GSL_SET_COMPLEX(&val, creal(X), cimag(X));
              gsl_vector_complex_set(vx, nmidx, val);
            }

          if (vdx)
            {
              X = rterm2 * dYnm2 - rterm * dYnm;
              GSL_SET_COMPLEX(&val, creal(X), cimag(X));
              gsl_vector_complex_set(vdx, nmidx, val);
            }

          if (vy)
            {
              Y = -rterm * m * invisint * Ynm;
              GSL_SET_COMPLEX(&val, creal(Y), cimag(Y));
              gsl_vector_complex_set(vy, nmidx, val);
            }

          if (vdy)
            {
              Y = -m * (rterm2 * invisint2 * Ynm2 -
                        rterm * invisint * Ynm);
              GSL_SET_COMPLEX(&val, creal(Y), cimag(Y));
              gsl_vector_complex_set(vdy, nmidx, val);
            }

          if (vz)
            {
              Z = -rterm * (n + 1.0) * Ynm;
              GSL_SET_COMPLEX(&val, creal(Z), cimag(Z));
              gsl_vector_complex_set(vz, nmidx, val);
            }

          if (vdz)
            {
              Z = -(n + 1.0) * (rterm2 * Ynm2 - rterm * Ynm);
              GSL_SET_COMPLEX(&val, creal(Z), cimag(Z));
              gsl_vector_complex_set(vdz, nmidx, val);
            }
        }
    }

  return s;
} /* poltor_row_pint() */

/*
poltor_row_pext()
  Build a row of the LS design matrix corresponding to poloidal
external magnetic field B_pol^e

Inputs: r     - radius (km)
        theta - colatitude (radians)
        vx    - (output) where to store matrix row X component (or NULL)
        vy    - (output) where to store matrix row Y component (or NULL)
        vz    - (output) where to store matrix row Z component (or NULL)
        w     - workspace

Return: success/error

Notes:
1) w->Ynm, and w->dYnm must be initialized prior to calling this function
*/

static int
poltor_row_pext(const double r, const double theta,
                const double r_ns, const double theta_ns,
                gsl_vector_complex *vx, gsl_vector_complex *vy, gsl_vector_complex *vz,
                gsl_vector_complex *vdx, gsl_vector_complex *vdy, gsl_vector_complex *vdz,
                poltor_workspace *w)
{
  int s = 0;
  double ratio = r / w->R;
  double rterm = 1.0; /* (r / R)^(n-1) */
  double ratio2 = r_ns / w->R;
  double rterm2 = 1.0; /* (r_ns / R)^(n-1) */
  const complex double invisint = I / sin(theta);
  const complex double invisint2 = I / sin(theta_ns);
  size_t n;
  gsl_complex val;

  /*
   * construct 3 rows of the design matrix corresponding to the
   * X,Y,Z data observation
   */
  for (n = 1; n <= w->nmax_ext; ++n)
    {
      int M = (int) GSL_MIN(n, w->mmax_ext);
      int m;

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t nmidx = poltor_nmidx(POLTOR_IDX_PEXT, n, m, w);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          complex double X, Y, Z, Ynm, dYnm, Ynm2, dYnm2;

          if (m >= 0)
            {
              Ynm = w->Ynm[pidx];
              dYnm = w->dYnm[pidx];
              Ynm2 = w->Ynm2[pidx];
              dYnm2 = w->dYnm2[pidx];
            }
          else
            {
              Ynm = conj(w->Ynm[pidx]);
              dYnm = conj(w->dYnm[pidx]);
              Ynm2 = conj(w->Ynm2[pidx]);
              dYnm2 = conj(w->dYnm2[pidx]);
            }

          if (vx)
            {
              X = rterm * dYnm;
              GSL_SET_COMPLEX(&val, creal(X), cimag(X));
              gsl_vector_complex_set(vx, nmidx, val);
            }

          if (vdx)
            {
              X = rterm2 * dYnm2 - rterm * dYnm;
              GSL_SET_COMPLEX(&val, creal(X), cimag(X));
              gsl_vector_complex_set(vdx, nmidx, val);
            }

          if (vy)
            {
              Y = -rterm * m * invisint * Ynm;
              GSL_SET_COMPLEX(&val, creal(Y), cimag(Y));
              gsl_vector_complex_set(vy, nmidx, val);
            }

          if (vdy)
            {
              Y = -m * (rterm2 * invisint2 * Ynm2 -
                        rterm * invisint * Ynm);
              GSL_SET_COMPLEX(&val, creal(Y), cimag(Y));
              gsl_vector_complex_set(vdy, nmidx, val);
            }

          if (vz)
            {
              Z = rterm * n * Ynm;
              GSL_SET_COMPLEX(&val, creal(Z), cimag(Z));
              gsl_vector_complex_set(vz, nmidx, val);
            }

          if (vdz)
            {
              Z = n * (rterm2 * Ynm2 - rterm * Ynm);
              GSL_SET_COMPLEX(&val, creal(Z), cimag(Z));
              gsl_vector_complex_set(vdz, nmidx, val);
            }
        }

      /* (r/R)^{n-1} */
      rterm *= ratio;

      /* (r_ns/R)^{n-1} */
      rterm2 *= ratio2;
    }

  return s;
} /* poltor_row_pext() */

/*
poltor_row_psh()
  Build a row of the LS design matrix corresponding to poloidal
magnetic field B_pol^{sh}

Inputs: r     - radius (km)
        theta - colatitude (radians)
        vx    - (output) where to store matrix row X component
        vy    - (output) where to store matrix row Y component
        vz    - (output) where to store matrix row Z component
        w     - workspace

Return: success/error

Notes:
1) w->Ynm, and w->dYnm must be initialized prior to calling this function
*/

static int
poltor_row_psh(const double r, const double theta,
               const double r_ns, const double theta_ns,
               gsl_vector_complex *vx, gsl_vector_complex *vy, gsl_vector_complex *vz,
               gsl_vector_complex *vdx, gsl_vector_complex *vdy, gsl_vector_complex *vdz,
               poltor_workspace *w)
{
  int s = 0;
  size_t n, j;
  const double ratio = w->R / r;
  const double ratio2 = w->R / r_ns;
  const complex double invisint = I / sin(theta);
  const complex double invisint2 = I / sin(theta_ns);
  gsl_complex val;

  /* outer loop over Taylor series */
  for (j = 0; j <= w->shell_J; ++j)
    {
      for (n = 1; n <= w->nmax_sh; ++n)
        {
          int M = (int) GSL_MIN(n, w->mmax_sh);
          int m;
          double An, Bn, ABsum, dABsum, ABsum2, dABsum2;

          /* compute Taylor series integrals */
          s = poltor_shell_An(n, j, r, &An, w);
          s += poltor_shell_Bn(n, j, r, &Bn, w);
          if (s)
            return s;

          ABsum = An + Bn;
          dABsum = -(double)n * An + (n + 1.0) * Bn;

          s = poltor_shell_An(n, j, r_ns, &An, w);
          s += poltor_shell_Bn(n, j, r_ns, &Bn, w);
          if (s)
            return s;

          ABsum2 = An + Bn;
          dABsum2 = -(double)n * An + (n + 1.0) * Bn;

          for (m = -M; m <= M; ++m)
            {
              int mabs = abs(m);
              size_t nmidx = poltor_jnmidx(j, n, m, w);
              size_t pidx = gsl_sf_legendre_array_index(n, mabs);
              complex double X, Y, Z, Ynm, dYnm, Ynm2, dYnm2;

              if (m >= 0)
                {
                  Ynm = w->Ynm[pidx];
                  dYnm = w->dYnm[pidx];
                  Ynm2 = w->Ynm2[pidx];
                  dYnm2 = w->dYnm2[pidx];
                }
              else
                {
                  Ynm = conj(w->Ynm[pidx]);
                  dYnm = conj(w->dYnm[pidx]);
                  Ynm2 = conj(w->Ynm2[pidx]);
                  dYnm2 = conj(w->dYnm2[pidx]);
                }

              /* compute (X,Y,Z) components of B_{lm} */

              if (vx)
                {
                  X = ratio * dABsum * dYnm;
                  GSL_SET_COMPLEX(&val, creal(X), cimag(X));
                  gsl_vector_complex_set(vx, nmidx, val);
                }

              if (vdx)
                {
                  X = ratio2 * dABsum2 * dYnm2 - ratio * dABsum * dYnm;
                  GSL_SET_COMPLEX(&val, creal(X), cimag(X));
                  gsl_vector_complex_set(vdx, nmidx, val);
                }

              if (vy)
                {
                  Y = -ratio * m * invisint * dABsum * Ynm;
                  GSL_SET_COMPLEX(&val, creal(Y), cimag(Y));
                  gsl_vector_complex_set(vy, nmidx, val);
                }

              if (vdy)
                {
                  Y = -m * (ratio2 * invisint2 * dABsum2 * Ynm2 -
                            ratio * invisint * dABsum * Ynm);
                  GSL_SET_COMPLEX(&val, creal(Y), cimag(Y));
                  gsl_vector_complex_set(vdy, nmidx, val);
                }

              if (vz)
                {
                  Z = ratio * n * (n + 1.0) * ABsum * Ynm;
                  GSL_SET_COMPLEX(&val, creal(Z), cimag(Z));
                  gsl_vector_complex_set(vz, nmidx, val);
                }

              if (vdz)
                {
                  Z = n * (n + 1.0) * (ratio2 * ABsum2 * Ynm2 - ratio * ABsum * Ynm);
                  GSL_SET_COMPLEX(&val, creal(Z), cimag(Z));
                  gsl_vector_complex_set(vdz, nmidx, val);
                }
            }
        }
    }

  return s;
} /* poltor_row_psh() */

/*
poltor_row_tor()
  Build a row of the LS design matrix corresponding to toroidal
magnetic field B_tor

Inputs: r     - radius (km)
        theta - colatitude (radians)
        vx    - (output) where to store matrix row X component (or NULL)
        vy    - (output) where to store matrix row Y component (or NULL)
        vz    - (output) where to store matrix row Z component (or NULL)
        w     - workspace

Return: success/error

Notes:
1) w->Ynm, and w->dYnm must be initialized prior to calling this function
*/

static int
poltor_row_tor(const double r, const double theta,
               const double r_ns, const double theta_ns,
               gsl_vector_complex *vx, gsl_vector_complex *vy, gsl_vector_complex *vz,
               gsl_vector_complex *vdx, gsl_vector_complex *vdy, gsl_vector_complex *vdz,
               poltor_workspace *w)
{
  int s = 0;
  const complex double invisint = I / sin(theta);
  const complex double invisint2 = I / sin(theta_ns);
  size_t n;
  gsl_complex val;

  /*
   * construct 3 rows of the design matrix corresponding to the
   * X,Y,Z data observation
   */
  for (n = 1; n <= w->nmax_tor; ++n)
    {
      int M = (int) GSL_MIN(n, w->mmax_tor);
      int m;

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t nmidx = poltor_nmidx(POLTOR_IDX_TOR, n, m, w);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          complex double X, Y, Z, Ynm, dYnm, Ynm2, dYnm2;

          if (m >= 0)
            {
              Ynm = w->Ynm[pidx];
              dYnm = w->dYnm[pidx];
              Ynm2 = w->Ynm2[pidx];
              dYnm2 = w->dYnm2[pidx];
            }
          else
            {
              Ynm = conj(w->Ynm[pidx]);
              dYnm = conj(w->dYnm[pidx]);
              Ynm2 = conj(w->Ynm2[pidx]);
              dYnm2 = conj(w->dYnm2[pidx]);
            }

          /* compute (r,theta,phi) components of B_{lm} */

          if (vx)
            {
              X = m * invisint * Ynm;
              GSL_SET_COMPLEX(&val, creal(X), cimag(X));
              gsl_vector_complex_set(vx, nmidx, val);
            }

          if (vdx)
            {
              X = m * (invisint2 * Ynm2 - invisint * Ynm);
              GSL_SET_COMPLEX(&val, creal(X), cimag(X));
              gsl_vector_complex_set(vdx, nmidx, val);
            }

          if (vy)
            {
              Y = dYnm;
              GSL_SET_COMPLEX(&val, creal(Y), cimag(Y));
              gsl_vector_complex_set(vy, nmidx, val);
            }

          if (vdy)
            {
              Y = dYnm2 - dYnm;
              GSL_SET_COMPLEX(&val, creal(Y), cimag(Y));
              gsl_vector_complex_set(vdy, nmidx, val);
            }

          if (vz)
            {
              Z = 0.0;
              GSL_SET_COMPLEX(&val, creal(Z), cimag(Z));
              gsl_vector_complex_set(vz, nmidx, val);
            }

          if (vdz)
            {
              Z = 0.0;
              GSL_SET_COMPLEX(&val, creal(Z), cimag(Z));
              gsl_vector_complex_set(vdz, nmidx, val);
            }
        }
    }

  return s;
} /* poltor_row_tor() */

#endif /*XXX*/

/*
poltor_nmidx()
  This function returns a unique index in [0,p-1] corresponding
to a given (n,m) pair. The array will look like:

[(1,-1) (1,0) (1,1) (2,-2) (2,-1) (2,0) (2,1) (2,2) ...]

(the (0,0) coefficient is not solved for)

Things get a little more tricky when mmax != nmax, so the
base indices of each 'n' block are precomputed in _alloc and
stored for easy reference here. The offset of (n,m) in a given
'n' block is:

offset = m + min(n,mmax)

which defaults to the standard m + n for the case mmax = nmax

Inputs: type - POLTOR_IDX_xxx
        n    - SH degree (> 0)
        m    - SH order (-n <= m <= n)

Return: index in [0,nnm-1]
*/

size_t
poltor_nmidx(const size_t type, const size_t n, const int m, poltor_workspace *w)
{
  size_t mmax;
  size_t type_offset; /* offset for coefficient type (internal/external) */
  size_t nmidx;

  if (type == POLTOR_IDX_PINT)
    {
      mmax = w->mmax_int;
      type_offset = w->pint_offset;
    }
  else if (type == POLTOR_IDX_PEXT)
    {
      mmax = w->mmax_ext;
      type_offset = w->pext_offset;
    }
  else if (type == POLTOR_IDX_PSH)
    {
      mmax = w->mmax_sh;
      type_offset = w->psh_offset;
    }
  else if (type == POLTOR_IDX_TOR)
    {
      mmax = w->mmax_tor;
      type_offset = w->tor_offset;
    }
  else
    {
      fprintf(stderr, "poltor_nmidx: error: unknown coefficient type\n");
      return 0;
    }

  if (n == 0)
    {
      fprintf(stderr, "poltor_nmidx: error: n = 0\n");
      return 0;
    }
  else if (abs(m) > (int) mmax)
    {
      fprintf(stderr, "poltor_nmidx: error: m = %d [mmax = %zu]\n", m, mmax);
      return 0;
    }

  if (n <= mmax)
    {
      size_t base1 = n * n; /* index of block for this n */
      int offset1 = m + n;  /* offset within block for this m */

      /* subtract 1 to exclude (0,0) coefficient */
      nmidx = base1 + offset1 - 1;
    }
  else
    {
      size_t base1 = (mmax + 1) * (mmax + 1);
      size_t base2 = (n - mmax - 1) * (2 * mmax + 1);
      int offset1 = m + (int)mmax;

      /* subtract 1 to exclude (0,0) coefficient */
      nmidx = base1 + base2 + offset1 - 1;
    }

  nmidx += type_offset;

  return nmidx;
} /* poltor_nmidx() */

/*
poltor_jnmidx()
  Return index of shell poloidal coefficient corresponding
to (j,n,m) where j is the term in the Taylor series expansion
for q_{nm}(r). Index is in [0,p-1]
*/

size_t
poltor_jnmidx(const size_t j, const size_t n, const int m, poltor_workspace *w)
{
  size_t nmidx = poltor_nmidx(POLTOR_IDX_PSH, n, m, w);

  return nmidx + j * w->nnm_sh;
} /* poltor_jnmidx() */

/*
poltor_jnmidx2()
  Return index of shell poloidal coefficient corresponding
to (j,n,m) where j is the term in the Taylor series expansion
for q_{nm}(r). Index is in [0,p_psh - 1]
*/

size_t
poltor_jnmidx2(const size_t j, const size_t n, const int m, poltor_workspace *w)
{
  size_t nmidx = green_idx(n, m, w->mmax_sh);
  return nmidx + j * w->nnm_sh;
}

/*
poltor_regularize()
  Regularize least-squares system with Tikhonov matrix

Inputs: L         - (output) diag(L) damping matrix
        w         - workspace

Notes:
1) Following Olsen et al, 1997, L = diag(alpha * n^3) so that
alpha^2 n^6 is added to the diagonal of the LS matrix
*/

int
poltor_regularize(gsl_vector *L, poltor_workspace *w)
{
  int s = 0;
  const double alpha_int = w->alpha_int;
  const double alpha_sh = w->alpha_sh;
  const double alpha_tor = w->alpha_tor;
  const double beta = 1.0; /* power of n */
  const double gamma_int = 2000.0;
  const double gamma_sh = 100.0;
  size_t n, j;

  gsl_vector_set_zero(L);

  for (n = 1; n <= w->nmax_int; ++n)
    {
#if 1
      double nfac = pow((double) n, beta);
      double val = alpha_int * nfac;
#else
      double val = alpha_int * pow((double) n, beta) + gamma_int;
#endif
      int M = (int) GSL_MIN(n, w->mmax_int);
      int m;

      for (m = -M; m <= M; ++m)
        {
          size_t nmidx = poltor_nmidx(POLTOR_IDX_PINT, n, m, w);
          gsl_vector_set(L, nmidx, val);
        }
    }

  {
    size_t nmax = w->nmax_sh;

    for (j = 0; j <= w->shell_J; ++j)
      {
        for (n = 1; n <= nmax; ++n)
          {
#if 1
            double nfac = pow((double) n, beta);
            double val = alpha_sh * nfac;
#else
            double val = alpha_sh * pow((double) n, beta) + gamma_sh;
#endif
            int M = (int) GSL_MIN(n, w->mmax_sh);
            int m;

            for (m = -M; m <= M; ++m)
              {
                size_t nmidx = poltor_jnmidx(j, n, m, w);
                gsl_vector_set(L, nmidx, val);
              }
          }
      }
  }

  for (n = 1; n <= w->nmax_tor; ++n)
    {
      double nfac = pow((double) n, beta);
      double val = alpha_tor * nfac;
      int M = (int) GSL_MIN(n, w->mmax_tor);
      int m;

      for (m = -M; m <= M; ++m)
        {
          size_t nmidx = poltor_nmidx(POLTOR_IDX_TOR, n, m, w);
          gsl_vector_set(L, nmidx, val);
        }
    }

  return s;
} /* poltor_regularize() */

/*
poltor_eval_chi()
  Calculate toroidal current scalar function chi at given point
on the spherical shell r

chi = -Q r / mu_0 

Inputs: r          - radius of current shell (km)
        theta      - colatitude (radians)
        phi        - longitude (radians)
        idx        - POLTOR_IDX_xxx
        chi        - (output) current scalar in units of kA
        w          - workspace

Return: success/error
*/

static int
poltor_eval_chi(const double r, const double theta, const double phi,
                const size_t idx, double *chi, poltor_workspace *w)
{
  int s = 0;
#if 0
  size_t n, m;
  complex double sum = 0.0;
  size_t nmax, mmax;

  if (idx == POLTOR_IDX_PINT)
    {
      nmax = w->nmax_int;
      mmax = w->mmax_int;
    }
  else if (idx == POLTOR_IDX_PSH)
    {
      nmax = w->nmax_sh;
      mmax = w->mmax_sh;
    }

  /* compute legendre functions */
  gsl_sf_legendre_deriv_alt_array(GSL_SF_LEGENDRE_SCHMIDT, nmax,
                                  cos(theta), w->Pnm, w->dPnm);

  /* pre-compute Ynm and dYnm */
  for (m = 0; m <= mmax; ++m)
    {
      complex double expimphi = cos(m * phi) + I * sin(m * phi);

      for (n = GSL_MAX(m, 1); n <= nmax; ++n)
        {
          size_t pidx = gsl_sf_legendre_array_index(n, m);

          w->Ynm[pidx] = w->Pnm[pidx] * expimphi;
          w->dYnm[pidx] = w->dPnm[pidx] * expimphi;
        }
    }

  for (n = 1; n <= nmax; ++n)
    {
      int M = (int) GSL_MIN(n, mmax);
      int m;

      for (m = -M; m <= M; ++m)
        {
          int mabs = abs(m);
          size_t nmidx = poltor_nmidx(idx, n, m, w);
          size_t pidx = gsl_sf_legendre_array_index(n, mabs);
          gsl_complex cnm = gsl_vector_complex_get(w->c, nmidx);
          complex double qnm = GSL_REAL(cnm) + I * GSL_IMAG(cnm);
          complex double Ynm;

          if (m >= 0)
            Ynm = w->Ynm[pidx];
          else
            Ynm = conj(w->Ynm[pidx]);

          sum += qnm * Ynm;
        }
    }

  /* current function = -Qr/mu0 (units of kA) */
  sum = -sum / POLTOR_MU_0 * r;
  *chi = creal(sum);
#endif

  return s;
} /* poltor_eval_chi() */

/*
poltor_init()
  Initialize poltor workspace (spatial weights, solar flux values)

Inputs: w - workspace
*/

static int
poltor_init(poltor_workspace * w)
{
  int s;
  const poltor_parameters *params = &(w->params);

  /* initialize spatial weights */
  s = poltor_init_weights(w);
  if (s)
    return s;

  /* initialize array of daily solar flux values */
  s = poltor_init_solar_flux(w);
  if (s)
    return s;

  return 0;
}

/*
poltor_init_weights()
  Initialize spatial weights

Inputs: w - workspace

Notes:
1) On output, w->wts_spatial is initialized with spatial weights
*/

static int
poltor_init_weights(poltor_workspace * w)
{
  int s = 0;
  size_t idx = 0;
  size_t i, j;
  const poltor_parameters *params = &(w->params);
  magdata_list *list = w->data;
  const size_t ntheta = 90; /* 2 deg bins */
  const size_t nphi = 24;   /* 15 deg bins */
  track_weight_workspace *weight_p = track_weight_alloc(ntheta, nphi);

  /* initialize spatial weighting histogram */

  fprintf(stderr, "poltor_init_weights: initializing weighting histogram...");

  for (i = 0; i < list->n; ++i)
    {
      magdata *mptr = magdata_list_ptr(i, list);

      for (j = 0; j < mptr->n; ++j)
        {
          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          if (MAGDATA_ExistX(mptr->flags[j]))
            track_weight_add_data(mptr->theta[j], mptr->phi[j], weight_p);

          if (MAGDATA_ExistY(mptr->flags[j]))
            track_weight_add_data(mptr->theta[j], mptr->phi[j], weight_p);

          if (MAGDATA_ExistZ(mptr->flags[j]))
            track_weight_add_data(mptr->theta[j], mptr->phi[j], weight_p);

          if (MAGDATA_ExistScalar(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
            track_weight_add_data(mptr->theta[j], mptr->phi[j], weight_p);

          if (MAGDATA_ExistDX_NS(mptr->flags[j]) || MAGDATA_ExistDX_EW(mptr->flags[j]))
            track_weight_add_data(mptr->theta[j], mptr->phi[j], weight_p);

          if (MAGDATA_ExistDY_NS(mptr->flags[j]) || MAGDATA_ExistDY_EW(mptr->flags[j]))
            track_weight_add_data(mptr->theta[j], mptr->phi[j], weight_p);

          if (MAGDATA_ExistDZ_NS(mptr->flags[j]) || MAGDATA_ExistDZ_EW(mptr->flags[j]))
            track_weight_add_data(mptr->theta[j], mptr->phi[j], weight_p);

          if (MAGDATA_ExistDF_NS(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
            track_weight_add_data(mptr->theta[j], mptr->phi[j], weight_p);

          if (MAGDATA_ExistDF_EW(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
            track_weight_add_data(mptr->theta[j], mptr->phi[j], weight_p);
        }
    }

  fprintf(stderr, "done\n");

  /* compute data weights with histogram */
  track_weight_calc(weight_p);

  /* calculate spatial weights */

  fprintf(stderr, "poltor_init_weights: calculating spatial weights...");

  for (i = 0; i < list->n; ++i)
    {
      magdata *mptr = magdata_list_ptr(i, list);

      for (j = 0; j < mptr->n; ++j)
        {
          double wt; /* spatial weight */

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

#if 1
          {
            double thetaq = M_PI / 2.0 - mptr->qdlat[j] * M_PI / 180.0;
            wt = sin(thetaq);
          }
#else
          track_weight_get(mptr->phi[j], mptr->theta[j], &wt, weight_p);
#endif

          if (MAGDATA_ExistX(mptr->flags[j]))
            gsl_vector_set(w->wts_spatial, idx++, params->weight_X * wt);

          if (MAGDATA_ExistY(mptr->flags[j]))
            gsl_vector_set(w->wts_spatial, idx++, params->weight_Y * wt);

          if (MAGDATA_ExistZ(mptr->flags[j]))
            gsl_vector_set(w->wts_spatial, idx++, params->weight_Z * wt);

          if (MAGDATA_ExistScalar(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
            gsl_vector_set(w->wts_spatial, idx++, params->weight_F * wt);

          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DX_EW))
            gsl_vector_set(w->wts_spatial, idx++, params->weight_DX_NS * wt);

          if (mptr->flags[j] & (MAGDATA_FLG_DY_NS | MAGDATA_FLG_DY_EW))
            gsl_vector_set(w->wts_spatial, idx++, params->weight_DY_NS * wt);

          if (mptr->flags[j] & (MAGDATA_FLG_DZ_NS | MAGDATA_FLG_DZ_EW))
            gsl_vector_set(w->wts_spatial, idx++, params->weight_DZ_NS * wt);

          if (MAGDATA_ExistDF_NS(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
            gsl_vector_set(w->wts_spatial, idx++, params->weight_DF_NS * wt);

          /*XXX*/
          if (MAGDATA_ExistDF_EW(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
            gsl_vector_set(w->wts_spatial, idx++, params->weight_F * wt);
        }
    }

  fprintf(stderr, "done\n");

  assert(idx == w->n);

  track_weight_free(weight_p);

  return s;
}

/*
poltor_init_solar_flux()
  Initialize w->solar_flux and w->solar_flux_grad with:
  
  1 + N * EUVAC
  
where N is the Wolf ratio constant and EUVAC is the daily
solar flux value in sfu
*/

static int
poltor_init_solar_flux(poltor_workspace * w)
{
  int s = 0;
  magdata_list *list = w->data;
  size_t i, j;
  size_t idx = 0;

  gsl_vector_set_zero(w->solar_flux);
  gsl_vector_set_zero(w->solar_flux_grad);

  fprintf(stderr, "poltor_init_solar_flux: computing daily solar flux values...");

  /* store daily solar flux values for all residuals */
  for (i = 0; i < list->n; ++i)
    {
      magdata *mptr = magdata_list_ptr(i, list);

      for (j = 0; j < mptr->n; ++j)
        {
          time_t unix_time;
          double euvac, alpha, alpha_grad;

          if (MAGDATA_Discarded(mptr->flags[j]))
            continue;

          /* calculate solar flux factor */
          unix_time = satdata_epoch2timet(mptr->t[j]);
          f107_get_euvac(unix_time, &euvac, w->f107_workspace_p);
          alpha = poltor_sflux_factor(euvac, w);

          /* calculate solar flux factor for gradient point, which might be different */
          if (mptr->flags[j] & (MAGDATA_FLG_DX_NS | MAGDATA_FLG_DY_NS | MAGDATA_FLG_DZ_NS | MAGDATA_FLG_DF_NS |
                                MAGDATA_FLG_DX_EW | MAGDATA_FLG_DY_EW | MAGDATA_FLG_DZ_EW | MAGDATA_FLG_DF_EW))
            {
              unix_time = satdata_epoch2timet(mptr->t_ns[j]);
              f107_get_euvac(unix_time, &euvac, w->f107_workspace_p);
              alpha_grad = poltor_sflux_factor(euvac, w);
            }
          else
            {
              alpha_grad = 0.0;
            }

          if (MAGDATA_ExistX(mptr->flags[j]))
            {
              gsl_vector_set(w->solar_flux, idx, alpha);
              gsl_vector_set(w->solar_flux_grad, idx++, alpha_grad);
            }

          if (MAGDATA_ExistY(mptr->flags[j]))
            {
              gsl_vector_set(w->solar_flux, idx, alpha);
              gsl_vector_set(w->solar_flux_grad, idx++, alpha_grad);
            }

          if (MAGDATA_ExistZ(mptr->flags[j]))
            {
              gsl_vector_set(w->solar_flux, idx, alpha);
              gsl_vector_set(w->solar_flux_grad, idx++, alpha_grad);
            }

          if (MAGDATA_ExistScalar(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
            {
              gsl_vector_set(w->solar_flux, idx, alpha);
              gsl_vector_set(w->solar_flux_grad, idx++, alpha_grad);
            }

          if (MAGDATA_ExistDX_NS(mptr->flags[j]) || MAGDATA_ExistDX_EW(mptr->flags[j]))
            {
              gsl_vector_set(w->solar_flux, idx, alpha);
              gsl_vector_set(w->solar_flux_grad, idx++, alpha_grad);
            }

          if (MAGDATA_ExistDY_NS(mptr->flags[j]) || MAGDATA_ExistDY_EW(mptr->flags[j]))
            {
              gsl_vector_set(w->solar_flux, idx, alpha);
              gsl_vector_set(w->solar_flux_grad, idx++, alpha_grad);
            }

          if (MAGDATA_ExistDZ_NS(mptr->flags[j]) || MAGDATA_ExistDZ_EW(mptr->flags[j]))
            {
              gsl_vector_set(w->solar_flux, idx, alpha);
              gsl_vector_set(w->solar_flux_grad, idx++, alpha_grad);
            }

          if (MAGDATA_ExistDF_NS(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
            {
              gsl_vector_set(w->solar_flux, idx, alpha);
              gsl_vector_set(w->solar_flux_grad, idx++, alpha_grad);
            }

          if (MAGDATA_ExistDF_EW(mptr->flags[j]) && MAGDATA_FitMF(mptr->flags[j]))
            {
              gsl_vector_set(w->solar_flux, idx, alpha);
              gsl_vector_set(w->solar_flux_grad, idx++, alpha_grad);
            }
        }
    }

  fprintf(stderr, "done\n");

  assert(idx == w->n);

  return s;
}
