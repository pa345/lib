/*
 * test_quat2.c
 *
 * min || B_chaos - R_z(theta_i) B_VFM ||
 *
 * calculate theta_i for each point along orbit
 */

#include "euler.h"
#include "oct.h"
#include "quat.h"

#define INVERT_QUAT   1

struct quat_data
{
  size_t n_data;
  satdata_mag *data;
  size_t downsample;
};

static int test_residuals(const char *filename, const gsl_vector *c, struct quat_data *params, satdata_mag *data);

static int
quat_apply_rot(const double q[4], double *B_VFM, double *B_NEC)
{
  /* B_NEC <- Rz(theta) R3(alpha,beta,gamma) B_VFM */
  quat_apply(q, B_VFM, B_NEC);

  return 0;
}

static int
quat_f(const gsl_vector * x, void * params, gsl_vector * f)
{
  struct quat_data *qdata = (struct quat_data *) params;
  satdata_mag *data = qdata->data;
  const size_t downsample = qdata->downsample;
  size_t i, j;
  size_t idx = 0;

  for (i = 0; i < data->n; i += downsample)
    {
      double B_model[4];
      double *B_VFM = &(data->B_VFM[3*i]);
      double *f3 = gsl_vector_ptr(f, 3 * idx);
#if INVERT_QUAT
      double thetai = gsl_vector_get(x, idx);
      double q[4];
#else
      double *q = &(data->q[4*i]);
#endif

#if INVERT_QUAT
      {
        double Rz_data[9];
        gsl_matrix_view Rz = gsl_matrix_view_array(Rz_data, 3, 3);

        gsl_matrix_set_zero(&Rz.matrix);

        gsl_matrix_set(&Rz.matrix, 0, 0, cos(thetai));
        gsl_matrix_set(&Rz.matrix, 0, 1, -sin(thetai));
        gsl_matrix_set(&Rz.matrix, 1, 0, sin(thetai));
        gsl_matrix_set(&Rz.matrix, 1, 1, cos(thetai));
        gsl_matrix_set(&Rz.matrix, 2, 2, 1.0);

        quat_R2q(&Rz.matrix, q);
      }
#endif

      quat_apply_rot(q, B_VFM, f3);

      satdata_mag_model(i, B_model, data);

      /* f <- Rz(theta) R3(alpha,beta,gamma) B_VFM - B_model */
      for (j = 0; j < 3; ++j)
        f3[j] -= B_model[j];

      idx++;
    }

#if 0
  printv_octave(f, "f");
  exit(1);
#endif

  return GSL_SUCCESS;
}

static int
quat_df(const gsl_vector * x, void * params, gsl_matrix * J)
{
  struct quat_data *qdata = (struct quat_data *) params;
  satdata_mag *data = qdata->data;
  const size_t downsample = qdata->downsample;
  size_t i;
  size_t idx = 0;

  gsl_matrix_set_zero(J);

  for (i = 0; i < data->n; i += downsample)
    {
      double *B_VFM = &(data->B_VFM[3*i]);
      double v[3];
#if INVERT_QUAT
      double thetai = gsl_vector_get(x, idx);
      double q[4];
      double Rz_data[9];
      gsl_matrix_view Rz = gsl_matrix_view_array(Rz_data, 3, 3);
#else
      double *q = &(data->q[4*i]);
#endif

#if INVERT_QUAT
      {
        gsl_matrix_set_zero(&Rz.matrix);

        gsl_matrix_set(&Rz.matrix, 0, 0, -sin(thetai));
        gsl_matrix_set(&Rz.matrix, 0, 1, -cos(thetai));
        gsl_matrix_set(&Rz.matrix, 1, 0, cos(thetai));
        gsl_matrix_set(&Rz.matrix, 1, 1, -sin(thetai));

        quat_R2q(&Rz.matrix, q);

        quat_apply(q, B_VFM, v);
        gsl_matrix_set(J, 3 * idx, idx, v[0]);
        gsl_matrix_set(J, 3 * idx + 1, idx, v[1]);
        gsl_matrix_set(J, 3 * idx + 2, idx, v[2]);
      }
#endif

      ++idx;
    }

  return GSL_SUCCESS;
}

static void
callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w)
{
  const char *res_file = "res.dat";
  struct quat_data *qdata = (struct quat_data *) params;
  gsl_vector *x = gsl_multifit_nlinear_position(w);
  gsl_vector *f = gsl_multifit_nlinear_residual(w);
  double rcond;

  gsl_multifit_nlinear_rcond(&rcond, w);

  fprintf(stderr,
          "iter %zu:\n"
          "\t |x|    = %.12e\n"
          "\t |f(x)| = %.12e\n"
          "\t cond(J) = %.12e\n",
          iter,
          gsl_blas_dnrm2(x),
          gsl_blas_dnrm2(f),
          1.0 / rcond);

  if (iter % 5 == 0)
    {
      fprintf(stderr, "printing residuals to %s...", res_file);
      test_residuals(res_file, x, qdata, qdata->data);
      fprintf(stderr, "done\n");
    }
}

static int
test_residuals(const char *filename, const gsl_vector *c, struct quat_data *params, satdata_mag *data)
{
  FILE *fp = fopen(filename, "w");
  size_t idx = 0;
  size_t i;

  for (i = 0; i < data->n; i += params->downsample)
    {
      double B_model[4];
      double *B_VFM = &(data->B_VFM[3*i]);
      double B_NEC[3];
#if INVERT_QUAT
      double thetai = gsl_vector_get(c, idx);
      double q[4];
#else
      double thetai = 0.0;
      double *q = &(data->q[4*i]);
#endif

#if INVERT_QUAT
      {
        double Rz_data[9];
        gsl_matrix_view Rz = gsl_matrix_view_array(Rz_data, 3, 3);

        gsl_matrix_set_zero(&Rz.matrix);

        gsl_matrix_set(&Rz.matrix, 0, 0, cos(thetai));
        gsl_matrix_set(&Rz.matrix, 0, 1, -sin(thetai));
        gsl_matrix_set(&Rz.matrix, 1, 0, sin(thetai));
        gsl_matrix_set(&Rz.matrix, 1, 1, cos(thetai));
        gsl_matrix_set(&Rz.matrix, 2, 2, 1.0);

        quat_R2q(&Rz.matrix, q);
      }
#endif

      quat_apply_rot(q, B_VFM, B_NEC);

      satdata_mag_model(i, B_model, data);

      fprintf(fp, "%f %f %f %f %f %f %f %f %f %f %f %d\n",
              data->qdlat[i],
              thetai * 180.0 / M_PI,
              B_VFM[0],
              B_VFM[1],
              B_VFM[2],
              B_NEC[0],
              B_NEC[1],
              B_NEC[2],
              B_model[0],
              B_model[1],
              B_model[2],
              satdata_satdir(i, data->n, data->latitude));

      idx++;
    }

  fclose(fp);

  return GSL_SUCCESS;
}

static int
test_quat(satdata_mag *data)
{
  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_workspace *w;
  gsl_multifit_nlinear_fdf fdf;
  gsl_multifit_nlinear_parameters fdf_params =
    gsl_multifit_nlinear_default_parameters();
  const size_t downsample = 60;
  size_t n_data = 0;
  size_t n, p;

  gsl_vector *f;
  gsl_vector *c;
  double chisq, chisq0;
  int status, info;
  struct quat_data params;
  size_t i;

  const double xtol = 1e-12;
  const double gtol = 1e-12;
  const double ftol = 0.0;

  for (i = 0; i < data->n; i += downsample)
    {
      ++n_data;
    }

  n = 3 * n_data;

#if INVERT_QUAT
  p = n_data; /* one rotation parameter for each VFM vector */
#else
  p = 3;
#endif

  c = gsl_vector_calloc(p);

  params.n_data = n_data;
  params.data = data;
  params.downsample = downsample;

  /* define the function to be minimized */
  fdf.f = quat_f;
#if 0
  fdf.df = NULL;
#else
  fdf.df = quat_df;
#endif
  fdf.fvv = NULL;     /* not using geodesic acceleration */
  fdf.n = n;
  fdf.p = p;
  fdf.params = &params;

  fprintf(stderr, "test_quat: n = %zu\n", n);
  fprintf(stderr, "test_quat: p = %zu\n", p);

  /* allocate workspace with default parameters */
  fdf_params.trs = gsl_multifit_nlinear_trs_lm;
  fdf_params.solver = gsl_multifit_nlinear_solver_cholesky;
  fdf_params.scale = gsl_multifit_nlinear_scale_levenberg;
  fdf_params.fdtype = GSL_MULTIFIT_NLINEAR_CTRDIFF;
  w = gsl_multifit_nlinear_alloc (T, &fdf_params, n, p);

  /* initialize solver with starting point and weights */
  gsl_multifit_nlinear_init (c, &fdf, w);

  /* compute initial cost function */
  f = gsl_multifit_nlinear_residual(w);
  gsl_blas_ddot(f, f, &chisq0);

  /* solve the system with a maximum of 20 iterations */
  status = gsl_multifit_nlinear_driver(200, xtol, gtol, ftol,
                                       callback, &params, &info, w);

  gsl_blas_ddot(f, f, &chisq);

  gsl_vector_memcpy(c, w->x);

  fprintf(stderr, "summary from method '%s/%s'\n",
          gsl_multifit_nlinear_name(w),
          gsl_multifit_nlinear_trs_name(w));
  fprintf(stderr, "number of iterations: %zu\n",
          gsl_multifit_nlinear_niter(w));
  fprintf(stderr, "function evaluations: %zu\n", fdf.nevalf);
  fprintf(stderr, "Jacobian evaluations: %zu\n", fdf.nevaldf);
  fprintf(stderr, "reason for stopping: %s\n",
          (info == 1) ? "small step size" : "small gradient");
  fprintf(stderr, "initial |f(x)| = %f\n", sqrt(chisq0));
  fprintf(stderr, "final   |f(x)| = %f\n", sqrt(chisq));

  printv_octave(c, "c");

  test_residuals("res.dat", c, &params, data);

  return GSL_SUCCESS;
}
