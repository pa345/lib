/*
 * euler_calc.c
 *
 * Calculate Euler angles for a given dataset by solving:
 *
 * min sum_i || R_q R_3(alpha,beta,gamma) B_i^{VFM} - B_i^{model} ||^2
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

#include <common/quat.h>

#include "euler.h"
#include "euler_calc.h"

static int euler_calc_f(const gsl_vector * x, void * params, gsl_vector * f);
static int euler_calc_df(const gsl_vector * x, void * params, gsl_matrix * J);
static void euler_calc_callback(const size_t iter, void * params, const gsl_multifit_nlinear_workspace * w);

/*
euler_calc_alloc()
  Allocate euler calc workspace

Inputs: n - maximum number of vector observations to ingest
*/

euler_calc_workspace *
euler_calc_alloc(const size_t n)
{
  euler_calc_workspace *w;

  w = calloc(1, sizeof(euler_calc_workspace));
  if (!w)
    {
      fprintf(stderr, "euler_calc_alloc: calloc failed: %s\n",
              strerror(errno));
      return NULL;
    }

  w->nmax = n;
  w->n = 0;
  w->p = 3;
  w->flags = EULER_FLG_ZYX;

  w->t = malloc(n * sizeof(double));
  w->qdlat = malloc(n * sizeof(double));
  w->q = malloc(4 * n * sizeof(double));
  w->B_VFM = malloc(3 * n * sizeof(double));
  w->B_model = malloc(3 * n * sizeof(double));
  w->c = gsl_vector_alloc(w->p);

  return w;
}

void
euler_calc_free(euler_calc_workspace * w)
{
  if (w->t)
    free(w->t);

  if (w->qdlat)
    free(w->qdlat);

  if (w->q)
    free(w->q);

  if (w->B_VFM)
    free(w->B_VFM);

  if (w->B_model)
    free(w->B_model);

  if (w->c)
    gsl_vector_free(w->c);

  free(w);
}

/*
euler_calc_add()
  Add a measurement to euler workspace

Inputs: t       - timestamp (CDF_EPOCH)
        qdlat   - QD latitude (degrees)
        B_VFM   - vector measurement in VFM frame
        B_model - field model vector in NEC frame
        q       - rotation quaternions from satellite-fixed to NEC
*/

int
euler_calc_add(const double t, const double qdlat, const double B_VFM[3], const double B_model[3],
               const double q[4], euler_calc_workspace * w)
{
  const size_t n = w->n;
  size_t i;

  if (n >= w->nmax)
    {
      fprintf(stderr, "euler_calc_add: error: too many observations added [%zu]\n", w->nmax);
      return GSL_FAILURE;
    }

  for (i = 0; i < 3; ++i)
    {
      w->B_VFM[3 * n + i] = B_VFM[i];
      w->B_model[3 * n + i] = B_model[i];
    }

  for (i = 0; i < 4; ++i)
    w->q[4 * n + i] = q[i];

  w->t[n] = t;
  w->qdlat[n] = qdlat;

  ++(w->n);

  return GSL_SUCCESS;
}

/*
euler_calc_proc()
  Compute Euler angles for previously specified dataset
*/

int
euler_calc_proc(euler_calc_workspace * w)
{
  int s = GSL_SUCCESS;
  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  const double xtol = 1.0e-8;
  const double gtol = 1.0e-8;
  const double ftol = 0.0;
  const size_t max_iter = 200;
  const size_t nres = 3 * w->n;
  gsl_vector * x;
  gsl_vector * weights = gsl_vector_alloc(nres);
  gsl_multifit_nlinear_parameters params =
    gsl_multifit_nlinear_default_parameters();
  gsl_multifit_nlinear_workspace *nlinear_p;
  gsl_multifit_nlinear_fdf fdf;
  int info;

  gsl_vector_set_all(weights, 1.0);
  gsl_vector_set_zero(w->c);

  fdf.f = euler_calc_f;
  fdf.df = euler_calc_df;
  fdf.fvv = NULL;
  fdf.n = nres;
  fdf.p = w->p;
  fdf.params = w;

  nlinear_p = gsl_multifit_nlinear_alloc(T, &params, nres, w->p);

  gsl_multifit_nlinear_winit(w->c, weights, &fdf, nlinear_p);

  s = gsl_multifit_nlinear_driver(max_iter, xtol, gtol, ftol,
                                  euler_calc_callback, NULL, &info, nlinear_p);

  /* save Euler angles */
  x = gsl_multifit_nlinear_position(nlinear_p);
  gsl_vector_memcpy(w->c, x);

  gsl_multifit_nlinear_free(nlinear_p);
  gsl_vector_free(weights);

  return s;
}

int
euler_calc_print_residuals(const char * filename, euler_calc_workspace * w)
{
  const double alpha = gsl_vector_get(w->c, EULER_IDX_ALPHA);
  const double beta = gsl_vector_get(w->c, EULER_IDX_BETA);
  const double gamma = gsl_vector_get(w->c, EULER_IDX_GAMMA);
  FILE *fp;
  size_t i;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "euler_calc_print_residuals: unable to open %s: %s\n",
              filename, strerror(errno));
      return GSL_FAILURE;
    }

  i = 1;
  fprintf(fp, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
  fprintf(fp, "# Field %zu: QD latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: X VFM original (nT)\n", i++);
  fprintf(fp, "# Field %zu: Y VFM original (nT)\n", i++);
  fprintf(fp, "# Field %zu: Z VFM original (nT)\n", i++);
  fprintf(fp, "# Field %zu: X VFM rotated to spacecraft-fixed with Euler angles (nT)\n", i++);
  fprintf(fp, "# Field %zu: Y VFM rotated to spacecraft-fixed with Euler angles (nT)\n", i++);
  fprintf(fp, "# Field %zu: Z VFM rotated to spacecraft-fixed with Euler angles (nT)\n", i++);
  fprintf(fp, "# Field %zu: X model rotated to spacecraft-fixed with quaternions (nT)\n", i++);
  fprintf(fp, "# Field %zu: Y model rotated to spacecraft-fixed with quaternions (nT)\n", i++);
  fprintf(fp, "# Field %zu: Z model rotated to spacecraft-fixed with quaternions (nT)\n", i++);

  for (i = 0; i < w->n; ++i)
    {
      double *B_VFM = &(w->B_VFM[3*i]);
      double *B_model = &(w->B_model[3*i]);
      double *q = &(w->q[4*i]);
      time_t unix_time = satdata_epoch2timet(w->t[i]);
      double B_VFM_sf[3];   /* VFM to spacecraft-fixed */
      double B_model_sf[3]; /* B_model to spacecraft-fixed */

      /* rotate B_VFM by R_3(alpha,beta,gamma) to spacecraft-fixed system */
      euler_apply_R3(w->flags, alpha, beta, gamma, B_VFM, B_VFM_sf);

      /* rotate B_model by R_q^{-1} to spacecraft-fixed system */
      quat_apply_inverse(q, B_model, B_model_sf);

      fprintf(fp, "%ld %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
              unix_time,
              w->qdlat[i],
              B_VFM[0],
              B_VFM[1],
              B_VFM[2],
              B_VFM_sf[0],
              B_VFM_sf[1],
              B_VFM_sf[2],
              B_model_sf[0],
              B_model_sf[1],
              B_model_sf[2]);
    }

  fclose(fp);

  return GSL_SUCCESS;
}

static int
euler_calc_f(const gsl_vector * x, void * params, gsl_vector * f)
{
  euler_calc_workspace * w = (euler_calc_workspace *) params;
  const double alpha = gsl_vector_get(x, EULER_IDX_ALPHA);
  const double beta = gsl_vector_get(x, EULER_IDX_BETA);
  const double gamma = gsl_vector_get(x, EULER_IDX_GAMMA);
  size_t i, j;

  for (i = 0; i < w->n; ++i)
    {
      double *B_VFM = &(w->B_VFM[3 * i]);
      double *B_model = &(w->B_model[3 * i]);
      double *q = &(w->q[4 * i]);
      double B_NEC[3];

      euler_vfm2nec(w->flags, alpha, beta, gamma, q, B_VFM, B_NEC);

      for (j = 0; j < 3; ++j)
        gsl_vector_set(f, 3 * i + j, B_NEC[j] - B_model[j]);
    }

  return GSL_SUCCESS;
}

static int
euler_calc_df(const gsl_vector * x, void * params, gsl_matrix * J)
{
  euler_calc_workspace * w = (euler_calc_workspace *) params;
  const double alpha = gsl_vector_get(x, EULER_IDX_ALPHA);
  const double beta = gsl_vector_get(x, EULER_IDX_BETA);
  const double gamma = gsl_vector_get(x, EULER_IDX_GAMMA);
  size_t i, j;

  for (i = 0; i < w->n; ++i)
    {
      double *B_VFM = &(w->B_VFM[3 * i]);
      double *q = &(w->q[4 * i]);
      double B_out[3];

      euler_vfm2nec(w->flags | EULER_FLG_DERIV_ALPHA, alpha, beta, gamma, q, B_VFM, B_out);

      for (j = 0; j < 3; ++j)
        gsl_matrix_set(J, 3 * i + j, EULER_IDX_ALPHA, B_out[j]);

      euler_vfm2nec(w->flags | EULER_FLG_DERIV_BETA, alpha, beta, gamma, q, B_VFM, B_out);

      for (j = 0; j < 3; ++j)
        gsl_matrix_set(J, 3 * i + j, EULER_IDX_BETA, B_out[j]);

      euler_vfm2nec(w->flags | EULER_FLG_DERIV_GAMMA, alpha, beta, gamma, q, B_VFM, B_out);

      for (j = 0; j < 3; ++j)
        gsl_matrix_set(J, 3 * i + j, EULER_IDX_GAMMA, B_out[j]);
    }

  return GSL_SUCCESS;
}

static void
euler_calc_callback(const size_t iter, void * params, const gsl_multifit_nlinear_workspace * w)
{
  gsl_vector *x = gsl_multifit_nlinear_position(w);
  gsl_vector *f = gsl_multifit_nlinear_residual(w);
  double alpha = gsl_vector_get(x, EULER_IDX_ALPHA);
  double beta = gsl_vector_get(x, EULER_IDX_BETA);
  double gamma = gsl_vector_get(x, EULER_IDX_GAMMA);
  double rcond;

  (void) params;

  gsl_multifit_nlinear_rcond(&rcond, w);

  fprintf(stderr,
          "iter %zu: alpha = %.2f [deg] beta = %.2f [deg] gamma = %.2f [deg] cond(J) = %g |f| = %g\n",
          iter,
          alpha * 180.0 / M_PI,
          beta * 180.0 / M_PI,
          gamma * 180.0 / M_PI,
          1.0 / rcond,
          gsl_blas_dnrm2(f));
}
