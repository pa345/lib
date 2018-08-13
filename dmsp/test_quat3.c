/*
 * test_quat2.c
 *
 * min || B_chaos - R_z(alpha_i) B_VFM ||
 *
 * calculate alpha_i for each point along orbit
 */

#include <gsl/gsl_test.h>

#include "euler.h"
#include "oct.h"
#include "quat.h"
#include "ellipsoid.h"

#define USE_QUAT          1

struct quat_data
{
  double B[2];
  double E[2];
  double *q;
};

static int
apply_q(const double *q, double v[2])
{
  double Rq_data[9];
  gsl_matrix_view Rq = gsl_matrix_view_array(Rq_data, 3, 3);
  double R11, R12, R21, R22;
  double tmp1, tmp2;

  quat_Rq(q, &Rq.matrix);

  R11 = gsl_matrix_get(&Rq.matrix, 0, 0);
  R12 = gsl_matrix_get(&Rq.matrix, 0, 1);
  R21 = gsl_matrix_get(&Rq.matrix, 1, 0);
  R22 = gsl_matrix_get(&Rq.matrix, 1, 1);

  tmp1 = R11*v[0] + R12*v[1];
  tmp2 = R21*v[0] + R22*v[1];

  v[0] = tmp1;
  v[1] = tmp2;

  return 0;
}

static int
quat_f(const gsl_vector * x, void * params, gsl_vector * f)
{
  struct quat_data *qdata = (struct quat_data *) params;
  double alpha = gsl_vector_get(x, 0);
  double st = sin(alpha);
  double ct = cos(alpha);
  double *E = qdata->E;
  double *B = qdata->B;
  double v[2];

  v[0] = E[0]*ct - E[1]*st;
  v[1] = E[0]*st + E[1]*ct;

#if USE_QUAT
  apply_q(qdata->q, v);
#endif

  gsl_vector_set(f, 0, B[0] - v[0]);
  gsl_vector_set(f, 1, B[1] - v[1]);

  return GSL_SUCCESS;
}

static int
quat_df(const gsl_vector * x, void * params, gsl_matrix * J)
{
  struct quat_data *qdata = (struct quat_data *) params;
  double alpha = gsl_vector_get(x, 0);
  double st = sin(alpha);
  double ct = cos(alpha);
  double *E = qdata->E;
  double v[2];

  v[0] = -E[0]*st - E[1]*ct;
  v[1] = E[0]*ct - E[1]*st;

#if USE_QUAT
  apply_q(qdata->q, v);
#endif

  gsl_matrix_set(J, 0, 0, -v[0]);
  gsl_matrix_set(J, 1, 0, -v[1]);

  return GSL_SUCCESS;
}

static void
callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w)
{
#if 0
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
#endif
}

/*
find_alpha()

Solve:

min || B - R(alpha) E ||

where

R(alpha) = [ cos(alpha) -sin(alpha) ]
           [ sin(alpha)  cos(alpha) ]
*/

static double
find_alpha(const double *q, const double B[2], const double E[2])
{
  const size_t n = 2;
  const size_t p = 1;
  const double xtol = 1e-12;
  const double gtol = 1e-12;
  const double ftol = 0.0;
  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_workspace *w;
  gsl_multifit_nlinear_fdf fdf;
  gsl_multifit_nlinear_parameters fdf_params =
    gsl_multifit_nlinear_default_parameters();
  gsl_vector *c = gsl_vector_calloc(p);
  double alpha;
  struct quat_data params;
  int info;
  gsl_vector *f;
  double chisq, chisq0;
  size_t i;

  params.q = q;
  for (i = 0; i < 2; ++i)
    {
      params.B[i] = B[i];
      params.E[i] = E[i];
    }

  /* define the function to be minimized */
  fdf.f = quat_f;
  fdf.df = quat_df;
  fdf.fvv = NULL;     /* not using geodesic acceleration */
  fdf.n = n;
  fdf.p = p;
  fdf.params = &params;

  /* allocate workspace with default parameters */
  fdf_params.trs = gsl_multifit_nlinear_trs_lm;
  fdf_params.solver = gsl_multifit_nlinear_solver_svd;
  fdf_params.scale = gsl_multifit_nlinear_scale_levenberg;
  w = gsl_multifit_nlinear_alloc (T, &fdf_params, n, p);

  /* initialize solver with starting point and weights */
  gsl_multifit_nlinear_init (c, &fdf, w);

  f = gsl_multifit_nlinear_residual(w);
  gsl_blas_ddot(f, f, &chisq0);

  /* solve the system with a maximum of 20 iterations */
  gsl_multifit_nlinear_driver(200, xtol, gtol, ftol, callback, &params, &info, w);

  gsl_blas_ddot(f, f, &chisq);
  gsl_vector_memcpy(c, w->x);

  alpha = gsl_vector_get(c, 0);

#if 0
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
#endif

  gsl_vector_free(c);
  gsl_multifit_nlinear_free(w);

  return alpha;
}

static int
calc_spacecraft_basis(const double r_ECEF[3], const double v_ECEF[3],
                      double s1[3], double s2[3], double s3[3])
{
  int s = 0;
  size_t i;
  double norm;

#if 0 /* use s3 = -e_mu */

  /* compute ECEF components of ellipsoid basis vectors, storing e_mu in s3 */
  ellipsoid_basis_mu(r_ECEF, WGS84_MU, s3, s1, s2);

#else /* use s3 = -rhat */

  /* store rhat in s3 */
  ecef2sph_basis(r_ECEF, s3, s1, s2);

#endif

  /* reverse s3 to point downward */
  for (i = 0; i < 3; ++i)
    s3[i] *= -1.0;

  /* s2 = (s3 x v) / | s3 x v | */
  vec_cross(s3, v_ECEF, s2);

  norm = vec_norm(s2);
  for (i = 0; i < 3; ++i)
    s2[i] /= norm;

  /* s1 = s2 x s3 */
  vec_cross(s2, s3, s1);

  return s;
}

static int
test_quat(satdata_mag *data)
{
  const size_t downsample = 60;
  size_t i;
  eph_data *eph_d = eph_data_read_tena("/data/DMSP/EPH/F15/2012DMSP/D.002.2012.DOP");
  eph_workspace *eph = eph_alloc(eph_d);

  for (i = 0; i < data->n; i += downsample)
    {
      double pos[3], vel[3], rhat[3], that[3], phat[3], vhat[3];
      double r = data->r[i];
      double theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
      double phi = data->longitude[i] * M_PI / 180.0;
      double B_model[4], B_VFM[3], E[3], B_rot[3], b_VFM_ECEF[3];
      double alpha, gamma, delta, gamma2, gamma3, gamma4, gamma5, beta;
      double alpha2, dec_model, dec_vfm;
      int dir = satdata_satdir(i, data->n, data->latitude);
      double Rq_data[9];
      gsl_matrix_view Rq = gsl_matrix_view_array(Rq_data, 3, 3);
      double *q = &(data->q[4*i]);
      size_t j;
      double s1[3], s2[3], s3[3];

      /* compute satellite position and velocity at this time */
      eph_interp(data->t[i], pos, vel, eph);

      calc_spacecraft_basis(pos, vel, s1, s2, s3);

      /* compute unit velocity vector in ECEF */
      for (j = 0; j < 3; ++j)
        vhat[j] = vel[j] / vec_norm(vel);

      ecef2sph_basis(pos, rhat, that, phat);

      /* calculate field model */
      satdata_mag_model(i, B_model, data);

      B_VFM[0] = SATDATA_VEC_X(data->B_VFM, i);
      B_VFM[1] = SATDATA_VEC_Y(data->B_VFM, i);
      B_VFM[2] = SATDATA_VEC_Z(data->B_VFM, i);

      alpha = find_alpha(q, B_model, B_VFM);

      dec_model = atan2(B_model[1], B_model[0]);
      dec_vfm = atan2(B_VFM[1], B_VFM[0]);
      alpha2 = dec_model - dec_vfm;

      {
        double ca = cos(alpha);
        double sa = sin(alpha);
        double normb;

        B_rot[0] = B_VFM[0] * ca - B_VFM[1] * sa;
        B_rot[1] = B_VFM[0] * sa + B_VFM[1] * ca;
        B_rot[2] = B_VFM[2];

#if USE_QUAT
        apply_q(q, B_rot);
#endif

#if 0
        /* convert NEC to sph */
        E[0] = -B_rot[2];
        E[1] = -B_rot[0];
        E[2] = B_rot[1];

        sph2ecef_vec(theta, phi, E, b_VFM_ECEF);

        normb = vec_norm(b_VFM_ECEF);
        for (j = 0; j < 3; ++j)
          b_VFM_ECEF[j] /= normb;
#endif
      }

      gamma = acos(-vec_dot(s1, that));
      delta = acos(vec_dot(s1, phat));

      if (dir == 1)
        beta = gamma - fabs(alpha);
      else
        beta = fabs(alpha) - gamma;

      quat_Rq(q, &Rq.matrix);
      gamma2 = acos(gsl_matrix_get(&Rq.matrix, 0, 0));
      gamma3 = asin(gsl_matrix_get(&Rq.matrix, 0, 1));
      gamma4 = asin(gsl_matrix_get(&Rq.matrix, 1, 0));
      gamma5 = acos(gsl_matrix_get(&Rq.matrix, 1, 1));

      /* sanity checks */
      gsl_test_rel(gamma3, gamma2, 1.0e-12, "quat gamma cos/sin");
      gsl_test_rel(gamma5, gamma2, 1.0e-12, "quat gamma cos/cos");
      gsl_test_rel(-gamma4, gamma3, 1.0e-12, "quat gamma sin/sin");
      gsl_test_rel(-gamma4, gamma2, 1.0e-12, "quat gamma sin/cos");

      if (dir == -1)
        print_octave(&Rq.matrix, "Rq");

      printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n",
             satdata_epoch2year(data->t[i]),
             data->qdlat[i],
             data->latitude[i],
             alpha * 180.0 / M_PI,
             wrap180(alpha2 * 180.0 / M_PI),
             gamma * 180.0 / M_PI,
             beta * 180.0 / M_PI,
/* X_VFM */  B_VFM[0],
/* Y_VFM */  B_VFM[1],
/* Z_VFM */  B_VFM[2],
/* X_NEC */  B_rot[0],
/* Y_NEC */  B_rot[1],
/* Z_NEC */  B_rot[2],
             B_model[0],
             B_model[1],
             B_model[2],
             -vec_dot(s1, rhat) * 180.0 / M_PI, /* chat */
             -vec_dot(s1, that) * 180.0 / M_PI, /* nhat, gamma */
             vec_dot(s1, phat) * 180.0 / M_PI,  /* ehat */
             gamma2 * 180.0 / M_PI,
             acos(vec_dot(s1, b_VFM_ECEF)) * 180.0 / M_PI,
             dir);
    }

  return GSL_SUCCESS;
}
