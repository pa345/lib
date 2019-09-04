/*
 * magcal.c
 *
 * This module contains routines for fitting the 9 magnetometer
 * calibration parameters (3 scale factors, 3 offsets,
 * 3 non-orthogonality angles).
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_errno.h>

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_oct.h>

#include "magcal.h"

static int magcal_f(const gsl_vector *m, void *params, gsl_vector *f);
static int magcal_df(const gsl_vector *m, void *params, gsl_matrix *J);
static int magcal_B_eval(const gsl_vector *m, const size_t i,
                         double B[4], const magcal_workspace *w);
static double magcal_f_eval(const gsl_vector *m, const size_t i,
                            const magcal_workspace *w);
static void magcal_callback (const size_t iter, void * params, const gsl_multifit_nlinear_workspace * w);

/*
magcal_alloc()
  Allocate a magcal workspace

Inputs: n - maximum number of data points to process
*/

magcal_workspace *
magcal_alloc(const size_t n)
{
  magcal_workspace *w;

  w = calloc(1, sizeof(magcal_workspace));
  if (!w)
    return 0;

  w->ntot = n;
  w->n = 0;
  w->p = MAGCAL_P;

  w->t = malloc(n * sizeof(double));
  w->Ex = malloc(n * sizeof(double));
  w->Ey = malloc(n * sizeof(double));
  w->Ez = malloc(n * sizeof(double));
  w->F = malloc(n * sizeof(double));

  w->covar = gsl_matrix_alloc(w->p, w->p);

  w->spline_f = gsl_spline_alloc(gsl_interp_akima, n);
  w->acc = gsl_interp_accel_alloc();

  w->lambda = 0.0;

  return w;
} /* magcal_alloc() */

void
magcal_free(magcal_workspace *w)
{
  if (w->t)
    free(w->t);

  if (w->Ex)
    free(w->Ex);

  if (w->Ey)
    free(w->Ey);

  if (w->Ez)
    free(w->Ez);

  if (w->F)
    free(w->F);

  if (w->nlinear_workspace_p)
    gsl_multifit_nlinear_free(w->nlinear_workspace_p);

  if (w->covar)
    gsl_matrix_free(w->covar);

  if (w->spline_f)
    gsl_spline_free(w->spline_f);

  if (w->acc)
    gsl_interp_accel_free(w->acc);

  free(w);
}

/*
magcal_add_datum()
  Add VFM vector measurement to workspace

Inputs: t     - timestamp (CDF_EPOCH)
        B_VFM - spacecraft-fixed vector measurement (nT)
        F     - scalar measurement or model (nT)
        w     - workspace
*/

int
magcal_add_datum(const double t, const double B_VFM[3], const double F, magcal_workspace *w)
{
  int s = 0;
  size_t n = w->n;

  w->t[n] = t;
  w->Ex[n] = B_VFM[0];
  w->Ey[n] = B_VFM[1];
  w->Ez[n] = B_VFM[2];
  w->F[n] = F;

  w->n = ++n;

  return s;
}

/*
magcal_proc()
  Compute the 9 magnetometer calibration parameters

Inputs: c - (input/output) vector of calibration parameters
            initialized on input with initial parameters
        w - workspace
*/

int
magcal_proc(gsl_vector *c, magcal_workspace *w)
{
  int s;
  const gsl_multifit_nlinear_type * T = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_parameters fdf_params =
    gsl_multifit_nlinear_default_parameters();
  gsl_vector *f;
  gsl_multifit_nlinear_fdf fdf;
  magcal_params params;
  const double xtol = 1e-6;
  const double gtol = 1e-6;
  const double ftol = 0.0;
  double chisq0, chisq;
  int info;

  params.w = w;
  params.dt = 0.0;

  fdf.f = &magcal_f;
  fdf.df = &magcal_df;
  fdf.fvv = NULL;
  fdf.n = w->n;
  fdf.p = w->p;
  fdf.params = &params;

  if (w->nlinear_workspace_p)
    gsl_multifit_nlinear_free(w->nlinear_workspace_p);

  fdf_params.solver = gsl_multifit_nlinear_solver_cholesky;
  w->nlinear_workspace_p = gsl_multifit_nlinear_alloc(T, &fdf_params, w->n, w->p);

  f = gsl_multifit_nlinear_residual(w->nlinear_workspace_p);

  gsl_multifit_nlinear_init(c, &fdf, w->nlinear_workspace_p);

  gsl_blas_ddot(f, f, &chisq0);

  s = gsl_multifit_nlinear_driver(500, xtol, gtol, ftol, magcal_callback, NULL,
                                  &info, w->nlinear_workspace_p);

  gsl_blas_ddot(f, f, &chisq);

  if (s != GSL_SUCCESS)
    {
      fprintf(stderr, "magcal_proc: error computing parameters: %s\n",
              gsl_strerror(s));
    }
  else
    {
      fprintf(stderr, "magcal_proc: number of data: %zu\n", w->n);
      fprintf(stderr, "magcal_proc: number of parameters: %zu\n", w->p);
      fprintf(stderr, "magcal_proc: number of iterations: %zu\n",
              gsl_multifit_nlinear_niter(w->nlinear_workspace_p));
      fprintf(stderr, "magcal_proc: function evaluations: %zu\n", fdf.nevalf);
      fprintf(stderr, "magcal_proc: jacobian evaluations: %zu\n", fdf.nevaldf);
      fprintf(stderr, "magcal_proc: reason for convergence: %d\n", info);
      fprintf(stderr, "initial |f(x)| = %.2f [nT]\n", sqrt(chisq0));
      fprintf(stderr, "final   |f(x)| = %.2f [nT]\n", sqrt(chisq));
      fprintf(stderr, "final residual rms = %.2f [nT]\n", sqrt(chisq / f->size));

      /* save calibration parameters */
      gsl_vector_memcpy(c, w->nlinear_workspace_p->x);
    }

  return s;
}

double
magcal_rms(const magcal_workspace * w)
{
  gsl_vector *f = gsl_multifit_nlinear_residual(w->nlinear_workspace_p);
  double chisq, rms;

  gsl_blas_ddot(f, f, &chisq);

  rms = sqrt(chisq / f->size);

  return rms;
}

time_t
magcal_mean_time(const magcal_workspace * w)
{
  double t_mean = gsl_stats_mean(w->t, 1, w->n);
  time_t unix_time = satdata_epoch2timet(t_mean);

  return unix_time;
}

/*
magcal_apply()
  Apply calibration parameters to satellite data

Inputs: m    - calibration parameters
        data - satellite data

Notes:
1) data->B_VFM and data->F are updated with new values
2) Only data not flagged with SATDATA_FLG_TIME are processed
*/

int
magcal_apply(const gsl_vector *m, satdata_mag *data)
{
  size_t i;

  for (i = 0; i < data->n; ++i)
    {
      double E[3], B[4];

      /* ignore flagged data */
      if (data->flags[i] & SATDATA_FLG_TIME)
        continue;

      E[0] = SATDATA_VEC_X(data->B_VFM, i);
      E[1] = SATDATA_VEC_Y(data->B_VFM, i);
      E[2] = SATDATA_VEC_Z(data->B_VFM, i);

      /* this function can be called with nT units */
      magcal_apply_cal(m, E, B);

      /* store new values in data */
      SATDATA_VEC_X(data->B_VFM, i) = B[0];
      SATDATA_VEC_Y(data->B_VFM, i) = B[1];
      SATDATA_VEC_Z(data->B_VFM, i) = B[2];
      data->F[i] = B[3];
    }

  return GSL_SUCCESS;
}

int
magcal_print_residuals(const char *filename, const gsl_vector * m, const magcal_workspace *w)
{
  size_t i;
  FILE *fp = fopen(filename, "w");

  i = 1;
  fprintf(fp, "# Field %zu: timestamp (UT seconds since 1970-01-01)\n", i++);
  fprintf(fp, "# Field %zu: original scalar measurement (nT)\n", i++);
  fprintf(fp, "# Field %zu: calibrated scalar measurement (nT)\n", i++);
  fprintf(fp, "# Field %zu: model scalar measurement (nT)\n", i++);

  for (i = 0; i < w->n; ++i)
    {
      double norm_E = gsl_hypot3(w->Ex[i], w->Ey[i], w->Ez[i]); /* uncorrected scalar field */
      double Fi = magcal_f_eval(m, i, w);

      /* compute magnetometer field magnitude with calibration parameters applied */

      fprintf(fp, "%ld %f %f %f\n",
              satdata_epoch2timet(w->t[i]),
              norm_E,
              Fi,
              w->F[i]);
    }

  fclose(fp);

  return GSL_SUCCESS;
}

/*
magcal_f()
  Compute residuals for least-squares fit of magnetometer calibration
parameters
*/

static int
magcal_f(const gsl_vector *m, void *params, gsl_vector *f)
{
  magcal_params *p = (magcal_params *) params;
  magcal_workspace *w = p->w;
  double dt = p->dt;
  size_t i;

  /* initialize output to 0 */
  gsl_vector_set_zero(f);

  for (i = 0; i < w->n; ++i)
    {
      double Fi;
      double ti = w->t[i] + dt;

      if (ti < w->t[0] || ti > w->t[w->n - 1])
        continue;

      /* compute magnetometer field magnitude with calibration parameters applied */
      Fi = magcal_f_eval(m, i, w);

      gsl_vector_set(f, i, Fi - w->F[i]);
    }

  return GSL_SUCCESS;
} /* magcal_f() */

static int
magcal_df(const gsl_vector *m, void *params, gsl_matrix *J)
{
  magcal_workspace *w = ((magcal_params *) params)->w;
  size_t i;
  double AXY = gsl_vector_get(m, MAGCAL_IDX_AXY);
  double AXZ = gsl_vector_get(m, MAGCAL_IDX_AXZ);
  double AYZ = gsl_vector_get(m, MAGCAL_IDX_AYZ);

  gsl_matrix_set_zero(J);

  for (i = 0; i < w->n; ++i)
    {
      gsl_vector_view v = gsl_matrix_row(J, i);
      double E[3], B[4];

      E[0] = w->Ex[i];
      E[1] = w->Ey[i];
      E[2] = w->Ez[i];

      magcal_apply_cal(m, E, B);

      gsl_vector_set(&v.vector, MAGCAL_IDX_SX, B[0] * E[0] / B[3]);
      gsl_vector_set(&v.vector, MAGCAL_IDX_SY, B[1] * E[1] / B[3]);
      gsl_vector_set(&v.vector, MAGCAL_IDX_SZ, B[2] * E[2] / B[3]);

      gsl_vector_set(&v.vector, MAGCAL_IDX_OX, B[0] / B[3]);
      gsl_vector_set(&v.vector, MAGCAL_IDX_OY, B[1] / B[3]);
      gsl_vector_set(&v.vector, MAGCAL_IDX_OZ, B[2] / B[3]);

      gsl_vector_set(&v.vector, MAGCAL_IDX_AXY, -B[0] * E[1] / B[3] * sin(AXY));
      gsl_vector_set(&v.vector, MAGCAL_IDX_AXZ, -B[0] * E[2] / B[3] * sin(AXZ));
      gsl_vector_set(&v.vector, MAGCAL_IDX_AYZ, -B[1] * E[2] / B[3] * sin(AYZ));
    }

  return GSL_SUCCESS;
} /* magcal_df() */

/*
magcal_B_eval()
  Compute corrected magnetometer field with calibration
parameters

Inputs: m - model parameters
        i - measurement point
        B - (output) corrected/calibrated field vector
            B[0] = B_x_calibrated (dimensionless)
            B[1] = B_y_calibrated (dimensionless)
            B[2] = B_z_calibrated (dimensionless)
            B[3] = |B|_calibrated (dimensionless)
        w - workspace

Return: success or error
*/

static int
magcal_B_eval(const gsl_vector *m, const size_t i,
              double B[4], const magcal_workspace *w)
{
  int s = 0;
  double E[3]; /* uncorrected magnetic measurements */

  E[0] = w->Ex[i];
  E[1] = w->Ey[i];
  E[2] = w->Ez[i];

  s += magcal_apply_cal(m, E, B);

  return s;
} /* magcal_B_eval() */

/*
magcal_f_eval()
  Evaluate magnetometer scalar magnitude with calibration parameters
applied

Inputs: m - model parameters
        i - measurement point
        w - workspace

Return: |B| where B is the magnetometer measured field with calibration
        parameters applied (dimensionless)
*/

static double
magcal_f_eval(const gsl_vector *m, const size_t i,
              const magcal_workspace *w)
{
  double B[4]; /* corrected vector field */

  magcal_B_eval(m, i, B, w);

  return B[3];
} /* magcal_f_eval() */

/*
magcal_apply_cal()
  Apply calibration to vector measurements

Inputs: m - model parameters
        E - original vector measurements
            E[0] = B_x_orig (any units)
            E[1] = B_y_orig
            E[2] = B_z_orig
        B - (output) calibrated vector measurements
            B[0] = B_x_calibrated (same units as E)
            B[1] = B_y_calibrated
            B[2] = B_z_calibrated
            B[3] = F_calibrated

Return: success or error
*/

int
magcal_apply_cal(const gsl_vector *m, const double E[3], double B[4])
{
  int s = 0;
  double SX = gsl_vector_get(m, MAGCAL_IDX_SX);
  double SY = gsl_vector_get(m, MAGCAL_IDX_SY);
  double SZ = gsl_vector_get(m, MAGCAL_IDX_SZ);
  double OX = gsl_vector_get(m, MAGCAL_IDX_OX);
  double OY = gsl_vector_get(m, MAGCAL_IDX_OY);
  double OZ = gsl_vector_get(m, MAGCAL_IDX_OZ);
  double AXY = gsl_vector_get(m, MAGCAL_IDX_AXY);
  double AXZ = gsl_vector_get(m, MAGCAL_IDX_AXZ);
  double AYZ = gsl_vector_get(m, MAGCAL_IDX_AYZ);

  /* apply calibration parameters (see Luhr et al, 2013, Eq 1) */
  B[0] = SX * E[0] + OX + cos(AXY) * E[1] + cos(AXZ) * E[2];
  B[1] = SY * E[1] + OY + cos(AYZ) * E[2];
  B[2] = SZ * E[2] + OZ;

  B[3] = gsl_hypot3(B[0], B[1], B[2]);

  return s;
} /* magcal_apply_cal() */

static void
magcal_callback (const size_t iter, void * params, const gsl_multifit_nlinear_workspace * w)
{
  gsl_vector *f = gsl_multifit_nlinear_residual(w);
  gsl_vector *x = gsl_multifit_nlinear_position(w);
  double rcond;

  (void) params;

  gsl_multifit_nlinear_rcond(&rcond, w);

  fprintf(stderr,
          "iter = %zu\n"
          "S = %15.8f %15.8f %15.8f\n"
          "O = %15.8f %15.8f %15.8f [nT]\n"
          "A = %15.8f %15.8f %15.8f [deg]\n"
          "|f(x)| = %g\n"
          "cond(J) = %g\n",
          iter,
          gsl_vector_get (x, MAGCAL_IDX_SX),
          gsl_vector_get (x, MAGCAL_IDX_SY),
          gsl_vector_get (x, MAGCAL_IDX_SZ),
          gsl_vector_get (x, MAGCAL_IDX_OX),
          gsl_vector_get (x, MAGCAL_IDX_OY),
          gsl_vector_get (x, MAGCAL_IDX_OZ),
          gsl_vector_get (x, MAGCAL_IDX_AXY) * 180.0 / M_PI,
          gsl_vector_get (x, MAGCAL_IDX_AXZ) * 180.0 / M_PI,
          gsl_vector_get (x, MAGCAL_IDX_AYZ) * 180.0 / M_PI,
          gsl_blas_dnrm2 (f),
          1.0 / rcond);
}
