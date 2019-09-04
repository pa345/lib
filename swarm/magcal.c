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
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_errno.h>

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_oct.h>

#include "magcal.h"

/* scale to dimensionless units */
#define MAGCAL_SCALE               1

/* Tikhonov damping parameter */
/*#define MAGCAL_LAMBDA              (1.0e-2)*/
#define MAGCAL_LAMBDA              (0.0)

#define MAGCAL_MAX_LATITUDE        (40.0)

static int magcal_init(const satdata_mag *data, magcal_workspace *w);
static int magcal_scale(const int dir, gsl_vector *m, magcal_workspace *w);
static int magcal_f(const gsl_vector *m, void *params, gsl_vector *f);
static int magcal_df(const gsl_vector *m, void *params, gsl_matrix *J);
static int magcal_B_eval(const gsl_vector *m, const size_t i,
                         double B[4], magcal_workspace *w);
static double magcal_f_eval(const gsl_vector *m, const size_t i,
                            magcal_workspace *w);
static void magcal_print_state (magcal_workspace * w,
                                gsl_multifit_fdfsolver * s);

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

  w->n = n;
  w->p = MAGCAL_P + 1; /* add 1 for dt */
  w->n_max = n;

  w->Ex = malloc(n * sizeof(double));
  w->Ey = malloc(n * sizeof(double));
  w->Ez = malloc(n * sizeof(double));
  w->F = malloc(n * sizeof(double));

  w->covar = gsl_matrix_alloc(w->p, w->p);
  w->fdf_type = gsl_multifit_fdfsolver_lmsder;
  w->fdf_s = gsl_multifit_fdfsolver_alloc(w->fdf_type, n, w->p);
  w->fdf_ridge = gsl_multifit_fdfridge_alloc(w->fdf_type, n, w->p);

  /* this is computed in magcal_init */
  w->B_s = 1.0;

  w->lambda = MAGCAL_LAMBDA;

  return w;
} /* magcal_alloc() */

void
magcal_free(magcal_workspace *w)
{
  if (w->Ex)
    free(w->Ex);

  if (w->Ey)
    free(w->Ey);

  if (w->Ez)
    free(w->Ez);

  if (w->F)
    free(w->F);

  if (w->fdf_s)
    gsl_multifit_fdfsolver_free(w->fdf_s);

  if (w->fdf_ridge)
    gsl_multifit_fdfridge_free(w->fdf_ridge);

  if (w->covar)
    gsl_matrix_free(w->covar);

  free(w);
}

/*
magcal_initcond()
  Set initial conditions for scalar calibration parameters
*/

int
magcal_initcond(gsl_vector *m)
{
  gsl_vector_set_zero(m);

  gsl_vector_set(m, MAGCAL_IDX_SX, 1.0);
  gsl_vector_set(m, MAGCAL_IDX_SY, 1.0);
  gsl_vector_set(m, MAGCAL_IDX_SZ, 1.0);
  gsl_vector_set(m, MAGCAL_IDX_OX, 0.0);
  gsl_vector_set(m, MAGCAL_IDX_OY, 0.0);
  gsl_vector_set(m, MAGCAL_IDX_OZ, 0.0);
  gsl_vector_set(m, MAGCAL_IDX_AXY, M_PI / 2.0);
  gsl_vector_set(m, MAGCAL_IDX_AXZ, M_PI / 2.0);
  gsl_vector_set(m, MAGCAL_IDX_AYZ, M_PI / 2.0);

  return GSL_SUCCESS;
}

/*
magcal_proc()
  Compute the 9 magnetometer calibration parameters

Inputs: m    - (input/output) vector of calibration parameters
               initialized on input with initial parameters
        data - satellite data
        w    - workspace
*/

int
magcal_proc(gsl_vector *m, const satdata_mag *data, magcal_workspace *w)
{
  int s;
  gsl_multifit_function_fdf f;
  gsl_multifit_fdfsolver *fdf_s;
  magcal_params params;
  gsl_vector_view v;
  const double xtol = 1e-10;
  const double gtol = 1e-10;
  const double ftol = 0.0;
  int info;

  /* copy data arrays */
  s = magcal_init(data, w);
  if (s)
    return s;

  /* scale parameters to dimensionless units */
  magcal_scale(1, m, w);

  params.w = w;

  f.f = &magcal_f;
  f.df = &magcal_df;
  f.n = w->n;
  f.p = w->p;
  f.params = &params;

  v = gsl_vector_subvector(m, MAGCAL_IDX_DT, w->p);
  gsl_multifit_fdfsolver_set (w->fdf_s, &f, &v.vector);
  gsl_multifit_fdfridge_set (w->fdf_ridge, &f, &v.vector, w->lambda);

#if 0
  s = gsl_multifit_fdfsolver_driver(w->fdf_s, 500, xtol, gtol, ftol, &info);
  fdf_s = w->fdf_s;
#else
  s = gsl_multifit_fdfridge_driver(w->fdf_ridge, 500, xtol, gtol, ftol, &info);
  fdf_s = w->fdf_ridge->s;
#endif
  if (s != GSL_SUCCESS)
    {
      fprintf(stderr, "magcal_proc: error computing parameters: %s\n",
              gsl_strerror(s));
    }
  else
    {
      magcal_print_state (w, fdf_s);
      fprintf(stderr, "magcal_proc: total data processed: %zu\n", w->n);
      fprintf(stderr, "magcal_proc: number of iterations: %zu\n",
              fdf_s->niter);
      fprintf(stderr, "magcal_proc: function evaluations: %zu\n",
              fdf_s->fdf->nevalf);
      fprintf(stderr, "magcal_proc: jacobian evaluations: %zu\n",
              fdf_s->fdf->nevaldf);
      fprintf(stderr, "magcal_proc: reason for convergence: %d\n",
              info);

      /* save calibration parameters */
      gsl_vector_memcpy(&v.vector, fdf_s->x);

      /* restore time shift parameter */
      gsl_vector_set(m, MAGCAL_IDX_DT, 0.0);

      /* scale offsets back to nT */
      magcal_scale(-1, m, w);

#if 0
      /* compute covariance matrix */
      gsl_multifit_covar(fdf_s->J, 0.0, w->covar);
#endif
    }

  return s;
} /* magcal_proc() */

/*
magcal_apply()
  Apply calibration parameters to satellite data

Inputs: m    - calibration parameters
        data - satellite data

Notes: data->B_VFM is updated with new values
*/

int
magcal_apply(const gsl_vector *m, satdata_mag *data)
{
  size_t i;

  for (i = 0; i < data->n; ++i)
    {
      double E[3], B[4];

      E[0] = SATDATA_VEC_X(data->B_VFM, i);
      E[1] = SATDATA_VEC_Y(data->B_VFM, i);
      E[2] = SATDATA_VEC_Z(data->B_VFM, i);

      /* this function can be called with nT units */
      magcal_apply_cal(m, E, B);

      /* store new values in data */
      SATDATA_VEC_X(data->B_VFM, i) = B[0];
      SATDATA_VEC_Y(data->B_VFM, i) = B[1];
      SATDATA_VEC_Z(data->B_VFM, i) = B[2];
    }

  return GSL_SUCCESS;
}

/*
magcal_init()
  Copy satellite data into local arrays. Only non-flagged data
is copied
*/

static int
magcal_init(const satdata_mag *data, magcal_workspace *w)
{
  int s = 0;
  size_t i;
  size_t n = 0;

  for (i = 0; i < data->n; ++i)
    {
      /* don't store flagged data */
      if (data->flags[i])
        continue;

      /* don't process high latitude data */
      if (fabs(data->latitude[i]) > MAGCAL_MAX_LATITUDE)
        continue;

      w->Ex[n] = SATDATA_VEC_X(data->B_VFM, i);
      w->Ey[n] = SATDATA_VEC_Y(data->B_VFM, i);
      w->Ez[n] = SATDATA_VEC_Z(data->B_VFM, i);

      w->F[n] = data->F[i];

      ++n;
    }

  if (n < 200)
    {
      fprintf(stderr, "magcal_init: insufficient data points for calibration: %zu\n",
              n);
      return -1;
    }

  if (n != w->n)
    {
      gsl_multifit_fdfsolver_free(w->fdf_s);
      gsl_multifit_fdfridge_free(w->fdf_ridge);
      w->fdf_s = gsl_multifit_fdfsolver_alloc(w->fdf_type, n, w->p);
      w->fdf_ridge = gsl_multifit_fdfridge_alloc(w->fdf_type, n, w->p);
      w->n = n;
    }

#if MAGCAL_SCALE
  w->B_s = GSL_MAX(gsl_stats_sd(w->Ex, 1, n),
                   GSL_MAX(gsl_stats_sd(w->Ey, 1, n),
                           gsl_stats_sd(w->Ez, 1, n)));
#endif

  /* center and scale data arrays */
  for (i = 0; i < n; ++i)
    {
      w->Ex[i] /= w->B_s;
      w->Ey[i] /= w->B_s;
      w->Ez[i] /= w->B_s;
      w->F[i] /= w->B_s;
    }

  return s;
} /* magcal_init() */

static int
magcal_scale(const int dir, gsl_vector *m, magcal_workspace *w)
{
  int s = 0;
  gsl_vector_view v = gsl_vector_subvector(m, MAGCAL_IDX_OX, 3);

  if (dir == 1) /* scale to dimensionless */
    gsl_vector_scale(&v.vector, 1.0 / w->B_s);
  else          /* scale to nT */
    gsl_vector_scale(&v.vector, w->B_s);

  return s;
} /* magcal_scale() */

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
  size_t i;

  /* initialize output to 0 */
  gsl_vector_set_zero(f);

  for (i = 0; i < w->n; ++i)
    {
      double Fi;

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
              double B[4], magcal_workspace *w)
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
                    magcal_workspace *w)
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
magcal_apply_cal(const gsl_vector *m, const double E[3],
                 double B[4])
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

  B[3] = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);

  return s;
} /* magcal_apply_cal() */

static void
magcal_print_state (magcal_workspace * w, gsl_multifit_fdfsolver * s)
{
  fprintf(stderr,
          "S = %15.8f %15.8f %15.8f\n"
          "O = %15.8f %15.8f %15.8f [nT]\n"
          "A = %15.8f %15.8f %15.8f [deg]\n"
          "|f(x)| = %g\n",
          gsl_vector_get (s->x, MAGCAL_IDX_SX),
          gsl_vector_get (s->x, MAGCAL_IDX_SY),
          gsl_vector_get (s->x, MAGCAL_IDX_SZ),
          gsl_vector_get (s->x, MAGCAL_IDX_OX) * w->B_s,
          gsl_vector_get (s->x, MAGCAL_IDX_OY) * w->B_s,
          gsl_vector_get (s->x, MAGCAL_IDX_OZ) * w->B_s,
          gsl_vector_get (s->x, MAGCAL_IDX_AXY) * 180.0 / M_PI,
          gsl_vector_get (s->x, MAGCAL_IDX_AXZ) * 180.0 / M_PI,
          gsl_vector_get (s->x, MAGCAL_IDX_AYZ) * 180.0 / M_PI,
          gsl_blas_dnrm2 (s->f));
}
