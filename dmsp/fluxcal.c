/*
 * fluxcal.c
 *
 * This module contains routines for fitting the 9 magnetometer
 * calibration parameters (3 scale factors, 3 offsets,
 * 3 non-orthogonality angles) using the Oersted/Olsen approach:
 *
 * B_calibrated = P S (E + O)
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

#include <satdata/satdata.h>
#include <common/common.h>
#include <common/oct.h>

#include "fluxcal.h"

static int fluxcal_f(const gsl_vector *m, void *params, gsl_vector *f);
static int fluxcal_df(const gsl_vector *m, void *params, gsl_matrix *J);
static int fluxcal_B_eval(const gsl_vector *m, const size_t i,
                          double B[4], const fluxcal_workspace *w);
static double fluxcal_f_eval(const gsl_vector *m, const size_t i,
                             const fluxcal_workspace *w);
static void fluxcal_callback (const size_t iter, void * params, const gsl_multifit_nlinear_workspace * w);
static int fluxcal_Pinv(const double u[3], gsl_matrix * Pinv);
static int fluxcal_Pinv_deriv(const double u[3], gsl_matrix * Pinv_1, gsl_matrix * Pinv_2, gsl_matrix * Pinv_3);
static int fluxcal_apply_datum_Pinv(const gsl_matrix * Pinv, const gsl_vector *m, const double E[3], double B[4]);
static int fluxcal_test_convergence(const double epsrel, const double epsabs, const gsl_vector * a, const gsl_vector * b);
static double fluxcal_bisquare(const double x);
static int fluxcal_robweights(const gsl_vector * r, gsl_vector * w);

/*
fluxcal_alloc()
  Allocate a fluxcal workspace

Inputs: n - maximum number of data points to process
*/

fluxcal_workspace *
fluxcal_alloc(const size_t n)
{
  fluxcal_workspace *w;

  w = calloc(1, sizeof(fluxcal_workspace));
  if (!w)
    return 0;

  w->ntot = n;
  w->n = 0;
  w->p = FLUXCAL_P;

  w->t = malloc(n * sizeof(double));
  w->Ex = malloc(n * sizeof(double));
  w->Ey = malloc(n * sizeof(double));
  w->Ez = malloc(n * sizeof(double));
  w->F = malloc(n * sizeof(double));

  w->weights = gsl_vector_alloc(n);
  w->covar = gsl_matrix_alloc(w->p, w->p);

  w->max_iter = 500;
  w->lambda = 0.0;

  return w;
}

void
fluxcal_free(fluxcal_workspace *w)
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

  if (w->weights)
    gsl_vector_free(w->weights);

  if (w->covar)
    gsl_matrix_free(w->covar);

  free(w);
}

/*
fluxcal_add_datum()
  Add VFM vector measurement to workspace

Inputs: t     - timestamp (CDF_EPOCH)
        B_VFM - spacecraft-fixed vector measurement (nT)
        F     - scalar measurement or model (nT)
        w     - workspace
*/

int
fluxcal_add_datum(const double t, const double B_VFM[3], const double F, fluxcal_workspace *w)
{
  int s = 0;
  size_t n = w->n;

  if (n >= w->ntot)
    {
      fprintf(stderr, "fluxcal_add_datum: error: ntot too small [%zu]\n", w->ntot);
      return -1;
    }

  w->t[n] = t;
  w->Ex[n] = B_VFM[0];
  w->Ey[n] = B_VFM[1];
  w->Ez[n] = B_VFM[2];
  w->F[n] = F;

  w->n = ++n;

  return s;
}

/*
fluxcal_proc()
  Compute the 9 magnetometer calibration parameters

Inputs: c - (input/output) vector of calibration parameters
            initialized on input with initial parameters
        w - workspace
*/

int
fluxcal_proc(gsl_vector *c, fluxcal_workspace *w)
{
  int s;
  const size_t max_robit = 100;
  const double epsrel = 1.0e-3;
  gsl_vector * r = gsl_vector_alloc(w->n);
  gsl_vector * cprev = gsl_vector_alloc(w->p);
  gsl_vector_view weights = gsl_vector_subvector(w->weights, 0, w->n);
  fluxcal_params params;
  double chisq, sigma, sigma0;
  size_t iter = 0;

  params.w = w;
  params.dt = 0.0;

  /* compute initial rms with initial guess vector */
  fluxcal_f(c, &params, r);
  gsl_blas_ddot(r, r, &chisq);
  sigma0 = sqrt(chisq / r->size);

  /* initial weights */
  gsl_vector_set_all(&weights.vector, 1.0);

  do
    {
      fprintf(stderr, "fluxcal_proc: === ROBUST ITERATION %zu ===\n", iter + 1);

      /* save current model parameters */
      gsl_vector_memcpy(cprev, c);

      /* solve for model parameters with current weights */
      s = fluxcal_nls(&weights.vector, c, w);
      if (s)
        return s;

      /* r = residuals / sigma */
      fluxcal_f(c, &params, r);
      gsl_blas_ddot(r, r, &chisq);
      sigma = sqrt(chisq / r->size);
      gsl_vector_scale(r, 1.0 / sigma);

      /* compute new robust weights */
      fluxcal_robweights(r, &weights.vector);

      s = fluxcal_test_convergence(epsrel, 0.0, cprev, c);
    }
  while (++iter < max_robit && s == GSL_CONTINUE);

  fprintf(stderr, "fluxcal_proc: number of robust iterations: %zu\n", iter);
  fprintf(stderr, "fluxcal_proc: initial scalar rms:          %g\n", sigma0);
  fprintf(stderr, "fluxcal_proc: final scalar rms:            %g\n", sigma);

  gsl_vector_free(cprev);
  gsl_vector_free(r);

  return s;
}

/*
fluxcal_nls()
  Compute the 9 magnetometer calibration parameters

Inputs: weights - robust weights, size w->n
        c       - (input/output) vector of calibration parameters
                  initialized on input with initial parameters
        w       - workspace

Notes:
1) On output, w->nlinear_workspace_p is allocated for the current size w->n
*/

int
fluxcal_nls(const gsl_vector * weights, gsl_vector * c, fluxcal_workspace *w)
{
  int s;
  const gsl_multifit_nlinear_type * T = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_parameters fdf_params =
    gsl_multifit_nlinear_default_parameters();
  gsl_vector *f;
  gsl_multifit_nlinear_fdf fdf;
  fluxcal_params params;
  const double xtol = 1e-6;
  const double gtol = 1e-6;
  const double ftol = 0.0;
  double chisq0, chisq;
  int info;

  params.w = w;
  params.dt = 0.0;

  fdf.f = &fluxcal_f;
  fdf.df = &fluxcal_df;
  fdf.fvv = NULL;
  fdf.n = w->n;
  fdf.p = w->p;
  fdf.params = &params;

  if (w->nlinear_workspace_p)
    gsl_multifit_nlinear_free(w->nlinear_workspace_p);

  /*fdf_params.solver = gsl_multifit_nlinear_solver_cholesky;*/
  fdf_params.solver = gsl_multifit_nlinear_solver_qr;
  w->nlinear_workspace_p = gsl_multifit_nlinear_alloc(T, &fdf_params, w->n, w->p);

  f = gsl_multifit_nlinear_residual(w->nlinear_workspace_p);

  gsl_multifit_nlinear_winit(c, weights, &fdf, w->nlinear_workspace_p);

  gsl_blas_ddot(f, f, &chisq0);

  s = gsl_multifit_nlinear_driver(w->max_iter, xtol, gtol, ftol, fluxcal_callback, NULL,
                                  &info, w->nlinear_workspace_p);

  gsl_blas_ddot(f, f, &chisq);

  if (s != GSL_SUCCESS)
    {
      fprintf(stderr, "fluxcal_nls: error computing parameters: %s\n",
              gsl_strerror(s));
    }
  else
    {
      fprintf(stderr, "fluxcal_nls: number of data: %zu\n", w->n);
      fprintf(stderr, "fluxcal_nls: number of parameters: %zu\n", w->p);
      fprintf(stderr, "fluxcal_nls: number of iterations: %zu\n",
              gsl_multifit_nlinear_niter(w->nlinear_workspace_p));
      fprintf(stderr, "fluxcal_nls: function evaluations: %zu\n", fdf.nevalf);
      fprintf(stderr, "fluxcal_nls: jacobian evaluations: %zu\n", fdf.nevaldf);
      fprintf(stderr, "fluxcal_nls: reason for convergence: %d\n", info);
      fprintf(stderr, "initial |f(x)| = %.2f [nT]\n", sqrt(chisq0));
      fprintf(stderr, "final   |f(x)| = %.2f [nT]\n", sqrt(chisq));
      fprintf(stderr, "final residual rms = %.2f [nT]\n", sqrt(chisq / f->size));

      /* save calibration parameters */
      gsl_vector_memcpy(c, w->nlinear_workspace_p->x);
    }

  return s;
}

double
fluxcal_rms(const fluxcal_workspace * w)
{
  gsl_vector *f = gsl_multifit_nlinear_residual(w->nlinear_workspace_p);
  double chisq, rms;

  gsl_blas_ddot(f, f, &chisq);

  rms = sqrt(chisq / f->size);

  return rms;
}

time_t
fluxcal_mean_time(const fluxcal_workspace * w)
{
  double t_mean = gsl_stats_mean(w->t, 1, w->n);
  time_t unix_time = satdata_epoch2timet(t_mean);

  return unix_time;
}

/*
fluxcal_apply()
  Apply calibration parameters to satellite data

Inputs: m    - calibration parameters
        data - satellite data

Notes:
1) data->B_VFM and data->F are updated with new values
2) Only data not flagged with SATDATA_FLG_TIME are processed
*/

int
fluxcal_apply(const gsl_vector *m, satdata_mag *data)
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
      fluxcal_apply_datum(m, E, B);

      /* store new values in data */
      SATDATA_VEC_X(data->B_VFM, i) = B[0];
      SATDATA_VEC_Y(data->B_VFM, i) = B[1];
      SATDATA_VEC_Z(data->B_VFM, i) = B[2];
      data->F[i] = B[3];
    }

  return GSL_SUCCESS;
}

int
fluxcal_print_residuals(const char *filename, const gsl_vector * m, const fluxcal_workspace *w)
{
  size_t i;
  FILE *fp = fopen(filename, "w");

  i = 1;
  fprintf(fp, "# Field %zu: timestamp (UT seconds since 1970-01-01)\n", i++);
  fprintf(fp, "# Field %zu: original scalar measurement (nT)\n", i++);
  fprintf(fp, "# Field %zu: calibrated scalar measurement (nT)\n", i++);
  fprintf(fp, "# Field %zu: model scalar measurement (nT)\n", i++);
  fprintf(fp, "# Field %zu: robust weight\n", i++);

  for (i = 0; i < w->n; ++i)
    {
      double norm_E = gsl_hypot3(w->Ex[i], w->Ey[i], w->Ez[i]); /* uncorrected scalar field */
      double Fi = fluxcal_f_eval(m, i, w);
      double wi = gsl_vector_get(w->weights, i);

      /* compute magnetometer field magnitude with calibration parameters applied */

      fprintf(fp, "%ld %12.4f %12.4f %12.4f %.4f\n",
              satdata_epoch2timet(w->t[i]),
              norm_E,
              Fi,
              w->F[i],
              wi);
    }

  fclose(fp);

  return GSL_SUCCESS;
}

/*
fluxcal_f()
  Compute residuals for least-squares fit of magnetometer calibration
parameters
*/

static int
fluxcal_f(const gsl_vector *m, void *params, gsl_vector *f)
{
  fluxcal_params *p = (fluxcal_params *) params;
  fluxcal_workspace *w = p->w;
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
      Fi = fluxcal_f_eval(m, i, w);

      gsl_vector_set(f, i, Fi - w->F[i]);
    }

  return GSL_SUCCESS;
}

static int
fluxcal_df(const gsl_vector *m, void *params, gsl_matrix *J)
{
  fluxcal_workspace *w = ((fluxcal_params *) params)->w;
  size_t i;
  double S[3], O[3], U[3];
  double Pinv_data[9], Pinv_1_data[9], Pinv_2_data[9], Pinv_3_data[9];
  gsl_matrix_view Pinv = gsl_matrix_view_array(Pinv_data, 3, 3);
  gsl_matrix_view Pinv_1 = gsl_matrix_view_array(Pinv_1_data, 3, 3);
  gsl_matrix_view Pinv_2 = gsl_matrix_view_array(Pinv_2_data, 3, 3);
  gsl_matrix_view Pinv_3 = gsl_matrix_view_array(Pinv_3_data, 3, 3);

  gsl_matrix_set_zero(J);

  S[0] = gsl_vector_get(m, FLUXCAL_IDX_SX);
  S[1] = gsl_vector_get(m, FLUXCAL_IDX_SY);
  S[2] = gsl_vector_get(m, FLUXCAL_IDX_SZ);

  O[0] = gsl_vector_get(m, FLUXCAL_IDX_OX);
  O[1] = gsl_vector_get(m, FLUXCAL_IDX_OY);
  O[2] = gsl_vector_get(m, FLUXCAL_IDX_OZ);

  U[0] = gsl_vector_get(m, FLUXCAL_IDX_U1);
  U[1] = gsl_vector_get(m, FLUXCAL_IDX_U2);
  U[2] = gsl_vector_get(m, FLUXCAL_IDX_U3);

  /* compute P^{-1} */
  fluxcal_Pinv(U, &Pinv.matrix);

  /* compute d/du_j P^{-1} matrices */
  fluxcal_Pinv_deriv(U, &Pinv_1.matrix, &Pinv_2.matrix, &Pinv_3.matrix);

  for (i = 0; i < w->n; ++i)
    {
      gsl_vector_view v = gsl_matrix_row(J, i);
      double E[3], B[4];
      double dB1[4], dB2[4], dB3[4];
      double tmp[3], tmp2[3];
      gsl_vector_view Bv = gsl_vector_view_array(B, 3);
      size_t j;

      E[0] = w->Ex[i];
      E[1] = w->Ey[i];
      E[2] = w->Ez[i];

      /* compute B = P^{-1} S (E - O) */
      fluxcal_apply_datum(m, E, B);

      /*
       * compute:
       * dB1 = (d/du_1 P^{-1}) S (E - O)
       * dB2 = (d/du_2 P^{-1}) S (E - O)
       * dB3 = (d/du_3 P^{-1}) S (E - O)
       */
      fluxcal_apply_datum_Pinv(&Pinv_1.matrix, m, E, dB1);
      fluxcal_apply_datum_Pinv(&Pinv_2.matrix, m, E, dB2);
      fluxcal_apply_datum_Pinv(&Pinv_3.matrix, m, E, dB3);

      /* compute tmp[j] = B . P^{-1}(:,j) */
      for (j = 0; j < 3; ++j)
        {
          gsl_vector_view Pinv_j = gsl_matrix_column(&Pinv.matrix, j);
          gsl_blas_ddot(&Bv.vector, &Pinv_j.vector, &tmp[j]);
        }

      /* compute tmp2[j] = B . [ (d/du_j P^{-1}) S (E - O) ] */
      tmp2[0] = vec_dot(B, dB1);
      tmp2[1] = vec_dot(B, dB2);
      tmp2[2] = vec_dot(B, dB3);

      /* df_i/dS_j = 1/|B| * B . dB/dS_j = (E_j - O_j) / |B| * (B . P^{-1}(:,j)) */

      gsl_vector_set(&v.vector, FLUXCAL_IDX_SX, (E[0] - O[0]) / B[3] * tmp[0]);
      gsl_vector_set(&v.vector, FLUXCAL_IDX_SY, (E[1] - O[1]) / B[3] * tmp[1]);
      gsl_vector_set(&v.vector, FLUXCAL_IDX_SZ, (E[2] - O[2]) / B[3] * tmp[2]);

      /* df_i/dO_j = 1/|B| * (B . dB/dO_j) = -S_j/|B| * (B . P^{-1}(:,j)) */

      gsl_vector_set(&v.vector, FLUXCAL_IDX_OX, -S[0] / B[3] * tmp[0]);
      gsl_vector_set(&v.vector, FLUXCAL_IDX_OY, -S[1] / B[3] * tmp[1]);
      gsl_vector_set(&v.vector, FLUXCAL_IDX_OZ, -S[2] / B[3] * tmp[2]);

      /* df_i/dU_j = 1/|B| * (B . dB/dU_j) = 1/|B| * B . [ d/dU_j P^{-1} S (E - O) ] */

      gsl_vector_set(&v.vector, FLUXCAL_IDX_U1, tmp2[0] / B[3]);
      gsl_vector_set(&v.vector, FLUXCAL_IDX_U2, tmp2[1] / B[3]);
      gsl_vector_set(&v.vector, FLUXCAL_IDX_U3, tmp2[2] / B[3]);
    }

  return GSL_SUCCESS;
}

/*
fluxcal_B_eval()
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
fluxcal_B_eval(const gsl_vector *m, const size_t i,
               double B[4], const fluxcal_workspace *w)
{
  int s = 0;
  double E[3]; /* uncorrected magnetic measurements */

  E[0] = w->Ex[i];
  E[1] = w->Ey[i];
  E[2] = w->Ez[i];

  s += fluxcal_apply_datum(m, E, B);

  return s;
}

/*
fluxcal_f_eval()
  Evaluate magnetometer scalar magnitude with calibration parameters
applied

Inputs: m - model parameters
        i - measurement point
        w - workspace

Return: |B| where B is the magnetometer measured field with calibration
        parameters applied (dimensionless)
*/

static double
fluxcal_f_eval(const gsl_vector *m, const size_t i, const fluxcal_workspace *w)
{
  double B[4]; /* corrected vector field */

  fluxcal_B_eval(m, i, B, w);

  return B[3];
}

/*
fluxcal_apply_datum()
  Apply calibration to a single vector measurement

Inputs: m - model parameters
        E - original vector measurements
            E[0] = B_x_orig (any units)
            E[1] = B_y_orig
            E[2] = B_z_orig
        B - (output) calibrated vector measurements
            B[0] = B_x_calibrated
            B[1] = B_y_calibrated
            B[2] = B_z_calibrated
            B[3] = F_calibrated

Return: success or error

Notes:
1) See Eq. 5 of Olsen et al, 2003
*/

int
fluxcal_apply_datum(const gsl_vector *m, const double E[3], double B[4])
{
  int s;
  double U[3], Pinv_data[9];
  gsl_matrix_view Pinv = gsl_matrix_view_array(Pinv_data, 3, 3);

  U[0] = gsl_vector_get(m, FLUXCAL_IDX_U1);
  U[1] = gsl_vector_get(m, FLUXCAL_IDX_U2);
  U[2] = gsl_vector_get(m, FLUXCAL_IDX_U3);

  /* construct P^{-1} */
  fluxcal_Pinv(U, &Pinv.matrix);

  /* compute B = P^{-1} S (E - O) */
  s = fluxcal_apply_datum_Pinv(&Pinv.matrix, m, E, B);

  return s;
}

static void
fluxcal_callback (const size_t iter, void * params, const gsl_multifit_nlinear_workspace * w)
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
          "U = %15.8f %15.8f %15.8f [deg]\n"
          "|f(x)| = %g\n"
          "cond(J) = %g\n",
          iter,
          gsl_vector_get (x, FLUXCAL_IDX_SX),
          gsl_vector_get (x, FLUXCAL_IDX_SY),
          gsl_vector_get (x, FLUXCAL_IDX_SZ),
          gsl_vector_get (x, FLUXCAL_IDX_OX),
          gsl_vector_get (x, FLUXCAL_IDX_OY),
          gsl_vector_get (x, FLUXCAL_IDX_OZ),
          gsl_vector_get (x, FLUXCAL_IDX_U1) * 180.0 / M_PI,
          gsl_vector_get (x, FLUXCAL_IDX_U2) * 180.0 / M_PI,
          gsl_vector_get (x, FLUXCAL_IDX_U3) * 180.0 / M_PI,
          gsl_blas_dnrm2 (f),
          1.0 / rcond);
}

/*
fluxcal_Pinv()
  Construct inverse matrix of non-orthogonality transformation:

P^{-1} = [                   1                                              0                0  ]
         [                 tan(u1)                                        sec(u1)            0  ]
         [ -1/(w*cos(u1)) * (sin(u1)*sin(u3) + cos(u1)*sin(u2))     -sin(u3)/(w*cos(u1))    1/w ]

where:

w = sqrt(1 - sin(u2)^2 - sin(u3)^2)
*/

static int
fluxcal_Pinv(const double u[3], gsl_matrix * Pinv)
{
  const double s1 = sin(u[0]);
  const double c1 = cos(u[0]);
  const double s2 = sin(u[1]);
  const double s3 = sin(u[2]);
  const double w = sqrt(1.0 - s2*s2 - s3*s3);

  gsl_matrix_set(Pinv, 0, 0, 1.0);
  gsl_matrix_set(Pinv, 0, 1, 0.0);
  gsl_matrix_set(Pinv, 0, 2, 0.0);

  gsl_matrix_set(Pinv, 1, 0, s1 / c1);
  gsl_matrix_set(Pinv, 1, 1, 1.0 / c1);
  gsl_matrix_set(Pinv, 1, 2, 0.0);

  gsl_matrix_set(Pinv, 2, 0, -1.0 / (w * c1) * (s1*s3 + c1*s2));
  gsl_matrix_set(Pinv, 2, 1, -s3/(w*c1));
  gsl_matrix_set(Pinv, 2, 2, 1.0 / w);

  return GSL_SUCCESS;
}

/*
fluxcal_Pinv_deriv()
  Construct matrices d/du_i P^{-1}

Inputs: u - non-orthogonality angles (radians)
        Pinv_1 - (output) d/du_1 P^{-1}
        Pinv_2 - (output) d/du_2 P^{-1}
        Pinv_3 - (output) d/du_3 P^{-1}

Return: success/error
*/

static int
fluxcal_Pinv_deriv(const double u[3], gsl_matrix * Pinv_1, gsl_matrix * Pinv_2, gsl_matrix * Pinv_3)
{
  const double s1 = sin(u[0]);
  const double c1 = cos(u[0]);
  const double s2 = sin(u[1]);
  const double c2 = cos(u[1]);
  const double s3 = sin(u[2]);
  const double c3 = cos(u[2]);
  const double w = sqrt(1.0 - s2*s2 - s3*s3);

  /* d/du_1 P^{-1} */

  gsl_matrix_set(Pinv_1, 0, 0, 0.0);
  gsl_matrix_set(Pinv_1, 0, 1, 0.0);
  gsl_matrix_set(Pinv_1, 0, 2, 0.0);

  gsl_matrix_set(Pinv_1, 1, 0, 1.0 / (c1 * c1));
  gsl_matrix_set(Pinv_1, 1, 1, s1 / (c1 * c1));
  gsl_matrix_set(Pinv_1, 1, 2, 0.0);

  gsl_matrix_set(Pinv_1, 2, 0, -s3 / (w * c1 * c1));
  gsl_matrix_set(Pinv_1, 2, 1, -s3 * s1 / (w * c1 * c1));
  gsl_matrix_set(Pinv_1, 2, 2, 0.0);

  /* d/du_2 P^{-1} */

  gsl_matrix_set(Pinv_2, 0, 0, 0.0);
  gsl_matrix_set(Pinv_2, 0, 1, 0.0);
  gsl_matrix_set(Pinv_2, 0, 2, 0.0);

  gsl_matrix_set(Pinv_2, 1, 0, 0.0);
  gsl_matrix_set(Pinv_2, 1, 1, 0.0);
  gsl_matrix_set(Pinv_2, 1, 2, 0.0);

  gsl_matrix_set(Pinv_2, 2, 0, -c2 / (w * w * w) * (c3 * c3 + s2 * s3 * (s1 / c1)));
  gsl_matrix_set(Pinv_2, 2, 1, -s3 * s2 * c2 / (c1 * w * w * w));
  gsl_matrix_set(Pinv_2, 2, 2, s2 * c2 / (w * w * w));

  /* d/du_3 P^{-1} */

  gsl_matrix_set(Pinv_3, 0, 0, 0.0);
  gsl_matrix_set(Pinv_3, 0, 1, 0.0);
  gsl_matrix_set(Pinv_3, 0, 2, 0.0);

  gsl_matrix_set(Pinv_3, 1, 0, 0.0);
  gsl_matrix_set(Pinv_3, 1, 1, 0.0);
  gsl_matrix_set(Pinv_3, 1, 2, 0.0);

  gsl_matrix_set(Pinv_3, 2, 0, -c3 / (w * w * w) * (s2 * s3 + c2 * c2 * (s1 / c1)));
  gsl_matrix_set(Pinv_3, 2, 1, -c3 * c2 * c2 / (c1 * w * w * w));
  gsl_matrix_set(Pinv_3, 2, 2, s3 * c3 / (w * w * w));

  return GSL_SUCCESS;
}

/*
fluxcal_apply_datum_Pinv()
  Apply calibration to a single vector measurement using a given
P^{-1} matrix:

B = P^{-1} S (E - O)

Inputs: Pinv - P^{-1} matrix
        m    - model parameters
        E    - original vector measurements
               E[0] = B_x_orig (any units)
               E[1] = B_y_orig
               E[2] = B_z_orig
        B    - (output) calibrated vector measurements
               B[0] = B_x_calibrated
               B[1] = B_y_calibrated
               B[2] = B_z_calibrated
               B[3] = F_calibrated

Return: success or error

Notes:
1) See Eq. 5 of Olsen et al, 2003
*/

static int
fluxcal_apply_datum_Pinv(const gsl_matrix * Pinv, const gsl_vector *m, const double E[3], double B[4])
{
  int s = 0;
  double S[3], O[3];
  gsl_vector_view Bv = gsl_vector_view_array(B, 3);
  size_t i;

  S[0] = gsl_vector_get(m, FLUXCAL_IDX_SX);
  S[1] = gsl_vector_get(m, FLUXCAL_IDX_SY);
  S[2] = gsl_vector_get(m, FLUXCAL_IDX_SZ);

  O[0] = gsl_vector_get(m, FLUXCAL_IDX_OX);
  O[1] = gsl_vector_get(m, FLUXCAL_IDX_OY);
  O[2] = gsl_vector_get(m, FLUXCAL_IDX_OZ);

  /* tmp = S (E - O) */
  for (i = 0; i < 3; ++i)
    B[i] = S[i] * (E[i] - O[i]);

  /* B = Pinv * S * (E - 0) */
  gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, Pinv, &Bv.vector);

  B[3] = gsl_hypot3(B[0], B[1], B[2]);

  return s;
}

static int
fluxcal_test_convergence(const double epsrel, const double epsabs, const gsl_vector * a, const gsl_vector * b)
{
  const size_t n = a->size;
  size_t i;

  for (i = 0; i < n; ++i)
    {
      double ai = gsl_vector_get(a, i);
      double bi = gsl_vector_get(b, i);
      
      if (fabs(ai - bi) > GSL_MAX(epsrel * fabs(bi), epsabs))
        return GSL_CONTINUE;
    }

  return GSL_SUCCESS;
}

static double
fluxcal_bisquare(const double x)
{
  if (fabs(x) <= 1.0)
    {
      double f = 1.0 - x*x;
      return (f * f);
    }
  else
    return 0.0;
}

static int
fluxcal_robweights(const gsl_vector * r, gsl_vector * w)
{
  const size_t n = r->size;
  size_t i;

  for (i = 0; i < n; ++i)
    {
      double ri = gsl_vector_get(r, i);
      double wi = fluxcal_bisquare(ri);
      gsl_vector_set(w, i, wi);
    }

  return 0;
}
