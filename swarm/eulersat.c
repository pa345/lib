/*
 * eulersat.c
 *
 * This module contains routines for computing Euler angles
 * to rotate to a satellite-fixed basis defined by
 *
 * s1 = v / |v|
 * s2 = s3 x s1
 * s3 = -rhat
 *
 * This module has 2 modes, defined by EULER_CALC_CRF2S
 * If undefined, it calculates the Euler angles going from
 * VFM to s1,s2,s3.
 *
 * Alternatively, if defined, Euler angles are computed
 * to go from CRF to s1,s2,s3, by computing angles such that
 * R_s R(alpha,beta,gamma) =~ R_q
 * This bypasses the magnetometer completely, and no magnetic
 * data, nor any magnetic main field model is used in computing
 * the angles.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <errno.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>

#include <mainlib/ml_ellipsoid.h>
#include <mainlib/ml_euler.h>

static int eulersat_proc(const size_t flags, gsl_vector *m, const satdata_mag *data);
static int eulersat_nonlinear_driver (gsl_multifit_fdfsolver * s,
                                      const size_t maxiter,
                                      const double xtol, const double gtol,
                                      const double ftol, int *info);
static int eulersat_f(const gsl_vector *m, void *params, gsl_vector *f);
static int eulersat_df(const gsl_vector *m, void *params, gsl_matrix *J);
static int eulersat_vfm2sat(const size_t flags, const gsl_vector *m,
                            const double B_vfm[3], double B_sat[3],
                            const satdata_mag *data);
static int eulersat_vfm2nec(const size_t flags, const size_t idx, const gsl_vector *m,
                            double B_nec[3], const satdata_mag *data);
static int eulersat_vel(const size_t idx, double V[3], const satdata_mag *data);
static int eulersat_sph_basis(const size_t idx, double rhat[3], double that[3],
                              double phat[3], const satdata_mag *data);
static int eulersat_basis(const size_t idx, double s1_hat[3], double s2_hat[3],
                          double s3_hat[3], const satdata_mag *data);
static int eulersat_Rs(const size_t idx, gsl_matrix *Rs,
                       const satdata_mag *data);
static void eulersat_print_state (size_t iter, gsl_multifit_fdfsolver * s);
static int eulersat_print_residuals(char *filename, const size_t flags,
                                    const gsl_vector *m, const satdata_mag *data);

#define EULER_IDX_ALPHA          0
#define EULER_IDX_BETA           1
#define EULER_IDX_GAMMA          2

/*
 * Define to compute angles alpha,beta,gamma such that:
 *
 * R_s R(alpha,beta,gamma) =~ R_q
 *
 * so that R(alpha,beta,gamma) rotates from CRF to s1,s2,s3
 * system
 */
#define EULER_CALC_CRF2S         0

static double V_CRF[3] = { -1.0, 1.0, 2.0 };

typedef struct
{
  size_t n;
  const satdata_mag *data;
  size_t flags;
} eulersat_params;

/*
eulersat_proc()
  Compute Euler angles alpha,beta,gamma such that

R_s R_3(alpha,beta,gamma) B_VFM - B_main

is minimized in a least-squares sense.

Inputs: flags - EULER_FLG_xxx
        m     - (output) vector of 3 Euler angles
        data  - satellite data

Return: success/error
*/

static int
eulersat_proc(const size_t flags, gsl_vector *m, const satdata_mag *data)
{
  int s = 0;
  size_t n = data->n - satdata_nflagged(data);
  int info;
  gsl_multifit_fdfsolver *fdf_s;
  gsl_multifit_function_fdf f;
  const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
  eulersat_params params;

  if (n < 500)
    {
      fprintf(stderr, "eulersat_proc: insufficient number of data to process: %zu\n",
              n);
      return -1;
    }

  /* initial values */
  gsl_vector_set(m, EULER_IDX_ALPHA, 12.0 * M_PI / 180.0);
  gsl_vector_set(m, EULER_IDX_BETA, -76.0 * M_PI / 180.0);
  gsl_vector_set(m, EULER_IDX_GAMMA, -12.0 * M_PI / 180.0);

  f.f = &eulersat_f;
#if 0
  f.df = NULL;
#else
  f.df = &eulersat_df;
#endif
  f.n = 3 * n; /* 3 field components */
  f.p = 3;
  f.params = &params;

  params.n = f.n;
  params.data = data;
  params.flags = flags;

  fprintf(stderr, "eulersat_proc: number of data = %zu\n", f.n);

  fdf_s = gsl_multifit_fdfsolver_alloc(T, f.n, f.p);

  fprintf(stderr, "eulersat_proc: initializing fdfsolver...");
  gsl_multifit_fdfsolver_set (fdf_s, &f, m);
  fprintf(stderr, "done\n");

  fprintf(stderr, "eulersat_proc: computing Euler angles...");
  s = eulersat_nonlinear_driver (fdf_s, 100, 1.0e-8, 1.0e-8, 0.0, &info);
  if (s != GSL_SUCCESS)
    fprintf(stderr, "eulersat_proc: error in fdfsolver: %d\n", s);
  fprintf(stderr, "done\n");

  eulersat_print_state(gsl_multifit_fdfsolver_niter(fdf_s), fdf_s);

  /* save angles */
  gsl_vector_memcpy(m, fdf_s->x);

  gsl_multifit_fdfsolver_free(fdf_s);

  return s;
} /* eulersat_proc() */

/*
eulersat_nonlinear_driver()
  Iterate the nonlinear least squares solver until completion

Inputs: s       - fdfsolver workspace
        maxiter - maximum iterations to allow
        xtol    - tolerance in step x
        gtol    - tolerance in gradient
        ftol    - tolerance in ||f||
        info    - (output) info flag on why iteration terminated
                  1 = stopped due to small step size ||dx|
                  2 = stopped due to small gradient
                  3 = stopped due to small change in f
                  GSL_ETOLX = ||dx|| has converged to within machine
                              precision (and xtol is too small)
                  GSL_ETOLG = ||g||_inf is smaller than machine
                              precision (gtol is too small)
                  GSL_ETOLF = change in ||f|| is smaller than machine
                              precision (ftol is too small)

Return: GSL_SUCCESS if converged, GSL_MAXITER if maxiter exceeded without
converging
*/

static int
eulersat_nonlinear_driver (gsl_multifit_fdfsolver * s,
                           const size_t maxiter,
                           const double xtol, const double gtol,
                           const double ftol, int *info)
{
  int status;
  size_t iter = 0;

  do
    {
      if (iter % 5 == 0 || iter == 1)
        eulersat_print_state(iter, s);

      status = gsl_multifit_fdfsolver_iterate (s);

      /*
       * if status is GSL_ENOPROG or GSL_SUCCESS, continue iterating,
       * otherwise the method has converged with a GSL_ETOLx flag
       */
      if (status != GSL_SUCCESS && status != GSL_ENOPROG)
        break;

      /* test for convergence */
      status = gsl_multifit_fdfsolver_test(s, xtol, gtol, ftol, info);
    }
  while (status == GSL_CONTINUE && ++iter < maxiter);

  /*
   * the following error codes mean that the solution has converged
   * to within machine precision, so record the error code in info
   * and return success
   */
  if (status == GSL_ETOLF || status == GSL_ETOLX || status == GSL_ETOLG)
    {
      *info = status;
      status = GSL_SUCCESS;
    }

  /* check if max iterations reached */
  if (iter >= maxiter && status != GSL_SUCCESS)
    status = GSL_EMAXITER;

  return status;
} /* eulersat_nonlinear_driver() */

static int
eulersat_f(const gsl_vector *m, void *params, gsl_vector *f)
{
  int s = GSL_SUCCESS;
  eulersat_params *p = ((eulersat_params *) params);
  const satdata_mag *data = p->data;
  size_t i;
  size_t idx = 0;

  gsl_vector_set_zero(f);

  for (i = 0; i < data->n; ++i)
    {
      double rhs[3], B_nec[3];

      /* skip flagged data */
      if (data->flags[i])
        continue;

      s = eulersat_vfm2nec(p->flags, i, m, B_nec, data);
      if (s)
        return s;

#if EULER_CALC_CRF2S
      {
        double *q = &(data->q[4*i]);
        euler_apply_Rq(q, V_CRF, rhs);
      }
#else
      rhs[0] = SATDATA_VEC_X(data->B_main, i);
      rhs[1] = SATDATA_VEC_Y(data->B_main, i);
      rhs[2] = SATDATA_VEC_Z(data->B_main, i);
#endif

      gsl_vector_set(f, idx++, B_nec[0] - rhs[0]);
      gsl_vector_set(f, idx++, B_nec[1] - rhs[1]);
      gsl_vector_set(f, idx++, B_nec[2] - rhs[2]);
    }

  assert(idx == p->n);

  return s;
} /* eulersat_f() */

static int
eulersat_df(const gsl_vector *m, void *params, gsl_matrix *J)
{
  eulersat_params *p = ((eulersat_params *) params);
  const satdata_mag *data = p->data;
  size_t i;
  size_t idx = 0;

  gsl_matrix_set_zero(J);

  for (i = 0; i < data->n; ++i)
    {
      double B_nec[3];

      /* skip flagged data */
      if (data->flags[i])
        continue;

      /* compute d/dalpha f_i */
      p->flags |= EULER_DERIV_ALPHA;
      eulersat_vfm2nec(p->flags, i, m, B_nec, data);
      p->flags &= ~EULER_DERIV_ALPHA;

      gsl_matrix_set(J, idx, EULER_IDX_ALPHA, B_nec[0]);
      gsl_matrix_set(J, idx + 1, EULER_IDX_ALPHA, B_nec[1]);
      gsl_matrix_set(J, idx + 2, EULER_IDX_ALPHA, B_nec[2]);

      /* compute d/dbeta f_i */
      p->flags |= EULER_DERIV_BETA;
      eulersat_vfm2nec(p->flags, i, m, B_nec, data);
      p->flags &= ~EULER_DERIV_BETA;

      gsl_matrix_set(J, idx, EULER_IDX_BETA, B_nec[0]);
      gsl_matrix_set(J, idx + 1, EULER_IDX_BETA, B_nec[1]);
      gsl_matrix_set(J, idx + 2, EULER_IDX_BETA, B_nec[2]);

      /* compute d/dgamma f_i */
      p->flags |= EULER_DERIV_GAMMA;
      eulersat_vfm2nec(p->flags, i, m, B_nec, data);
      p->flags &= ~EULER_DERIV_GAMMA;

      gsl_matrix_set(J, idx, EULER_IDX_GAMMA, B_nec[0]);
      gsl_matrix_set(J, idx + 1, EULER_IDX_GAMMA, B_nec[1]);
      gsl_matrix_set(J, idx + 2, EULER_IDX_GAMMA, B_nec[2]);

      idx += 3;
    }

  assert(idx == p->n);

  return GSL_SUCCESS;
} /* euler_df() */

/*
eulersat_vfm2sat()
  Convert from VFM axes to satellite fixed frame.
These two coordinate axes are both fixed with respect to the satellite,
and are therefore not expected to change over the orbit. Therefore,
a standard Euler rotation with 3 Euler angles should bring them into
alignment.

The transformation applied is:

B_sat = R_x(alpha) R_y(beta) R_z(gamma) B_vfm

where alpha, beta, gamma are the Euler angles

Inputs: m     - euler parameters (Euler angles)
        B_vfm - input vector in VFM frame
        B_sat - (output) rotated vector in satellite frame
        w     - workspace

Notes: w->flags can control the output of this function. If
the EULER_DERIV_xxx flags are set, we take derivatives
of the corresponding rotation matrices when converting from
the magnetometer frame to satellite frame
*/

static int
eulersat_vfm2sat(const size_t flags, const gsl_vector *m,
                 const double B_vfm[3], double B_sat[3],
                 const satdata_mag *data)
{
  int s = 0;
  double alpha = gsl_vector_get(m, EULER_IDX_ALPHA),
         beta = gsl_vector_get(m, EULER_IDX_BETA),
         gamma = gsl_vector_get(m, EULER_IDX_GAMMA);

  /* apply 3D Euler rotation */
  euler_apply_deriv(flags, alpha, beta, gamma, B_vfm, B_sat);

  return s;
} /* eulersat_vfm2sat() */

/*
eulersat_vfm2nec()
  Convert a B vector measurement in the VFM frame to NEC

Inputs: idx   - measurement index to transform
        m     - euler parameters (for Euler angle rotation)
        B_nec - (output) geocentric magnetic field vector in
                north-east-down coordinates
                B_nec[0] = B_x (dimensionless) (north)
                B_nec[1] = B_y (dimensionless) (east)
                B_nec[2] = B_z (dimensionless) (down)
        w     - euler workspace

Return: success or error

Notes: w->flags can control the output of this function. If
the EULER_DERIV_xxx flags are set, we take derivatives
of the corresponding rotation matrices when converting from
the magnetometer frame to satellite frame
*/

static int
eulersat_vfm2nec(const size_t flags, const size_t idx, const gsl_vector *m,
                 double B_nec[3], const satdata_mag *data)
{
  int s = 0;
  double B_vfm[3],  /* measured B in VFM frame */
         B_sat[3],  /* measured B in satellite-fixed frame */
         B_eci[3];  /* measured B in Cartesian ECI frame */
  double s1_hat[3], s2_hat[3], s3_hat[3];
  double rhat[3], that[3], phat[3];
  size_t j;

#if EULER_CALC_CRF2S
  /* replace B_vfm with V_CRF, a constant vector in the CRF frame */
  B_vfm[0] = V_CRF[0];
  B_vfm[1] = V_CRF[1];
  B_vfm[2] = V_CRF[2];
#else
  B_vfm[0] = SATDATA_VEC_X(data->B_VFM, idx);
  B_vfm[1] = SATDATA_VEC_Y(data->B_VFM, idx);
  B_vfm[2] = SATDATA_VEC_Z(data->B_VFM, idx);
#endif

  if (!(flags & EULER_FLG_ALTBASIS))
    {
      /* use star camera quaternions for VFM->NEC transformation */

      double alpha = gsl_vector_get(m, EULER_IDX_ALPHA);
      double beta = gsl_vector_get(m, EULER_IDX_BETA);
      double gamma = gsl_vector_get(m, EULER_IDX_GAMMA);
      double *q = &(data->q[4*idx]);

      s = euler_vfm2nec_deriv(flags, alpha, beta, gamma, q, B_vfm, B_nec);
      return s;
    }

  /* convert from VFM axes to satellite-fixed frame */
  eulersat_vfm2sat(flags, m, B_vfm, B_sat, data);

#if 0

  {
    double Rs_data[9];
    gsl_matrix_view Rs = gsl_matrix_view_array(Rs_data, 3, 3);
    gsl_vector_view sat = gsl_vector_view_array(B_sat, 3);
    gsl_vector_view nec = gsl_vector_view_array(B_nec, 3);

    /* construct R_s matrix explicitly */
    eulersat_Rs(idx, &Rs.matrix, data);

    /* rotate: B_nec = R_s B_sat */
    gsl_blas_dgemv(CblasNoTrans, 1.0, &Rs.matrix, &sat.vector, 0.0, &nec.vector);
  }

#else

  /*
   * compute Cartesian ECI components of satellite frame basis vectors
   * at time t
   */
  eulersat_basis(idx, s1_hat, s2_hat, s3_hat, data);

  /* convert B_sat to ECI frame */
  for (j = 0; j < 3; ++j)
    {
      B_eci[j] = B_sat[0] * s1_hat[j] +
                 B_sat[1] * s2_hat[j] +
                 B_sat[2] * s3_hat[j];
    }

  /* compute ECI components of spherical basis vectors */
  eulersat_sph_basis(idx, rhat, that, phat, data);

  /* convert to local north-east-down coordinates */
  B_nec[0] = -vec_dot(B_eci, that); /* B_x = -B_t */
  B_nec[1] = vec_dot(B_eci, phat);  /* B_y = B_p */
  B_nec[2] = -vec_dot(B_eci, rhat); /* B_z = -B_r */

#endif

  return s;
} /* eulersat_vfm2nec() */

/* compute ECI velocity using finite differences */
static int
eulersat_vel(const size_t idx, double V[3], const satdata_mag *data)
{
  int s = 0;
  double r0 = data->altitude[idx] + data->R;
  double theta0 = M_PI / 2.0 - data->latitude[idx] * M_PI / 180.0;
  double phi0 = data->longitude[idx] * M_PI / 180.0;
  double r1, theta1, phi1;
  double ecef0[3], ecef1[3];
  double eci0[3], eci1[3];
  time_t t0, t1, dt;
  size_t i;

  assert(idx < data->n - 1);

  t0 = satdata_epoch2timet(data->t[idx]);
  t1 = satdata_epoch2timet(data->t[idx + 1]);
  dt = t1 - t0;

  if (dt > 5)
    {
      fprintf(stderr, "eulersat_vel: warning: dt = %ld [s]\n", dt);
    }

  r1 = data->altitude[idx + 1] + data->R;
  theta1 = M_PI / 2.0 - data->latitude[idx + 1] * M_PI / 180.0;
  phi1 = data->longitude[idx + 1] * M_PI / 180.0;

  ecef0[0] = r0 * sin(theta0) * cos(phi0);
  ecef0[1] = r0 * sin(theta0) * sin(phi0);
  ecef0[2] = r0 * cos(theta0);

  ecef1[0] = r1 * sin(theta1) * cos(phi1);
  ecef1[1] = r1 * sin(theta1) * sin(phi1);
  ecef1[2] = r1 * cos(theta1);

  ecef2eci_pos(t0, ecef0, eci0);
  ecef2eci_pos(t1, ecef1, eci1);

  for (i = 0; i < 3; ++i)
    V[i] = (eci1[i] - eci0[i]) / (double)dt;

  return s;
} /* eulersat_vel() */

/*
eulersat_sph_basis()
  Compute ECI components of spherical basis vectors for a given
time t
*/

static int
eulersat_sph_basis(const size_t idx, double rhat[3], double that[3],
                   double phat[3], const satdata_mag *data)
{
  int s = 0;
  double theta = M_PI / 2.0 - data->latitude[idx] * M_PI / 180.0;
  double phi = data->longitude[idx] * M_PI / 180.0;
  double rhat_ecef[3], that_ecef[3], phat_ecef[3]; /* ECEF components */

  /* compute ECEF components of spherical unit vectors */
  sph_basis(theta, phi, rhat_ecef, that_ecef, phat_ecef);

  {
    time_t unix_time = satdata_epoch2timet(data->t[idx]);

    /* transform to ECI components */
    ecef2eci_pos(unix_time, rhat_ecef, rhat);
    ecef2eci_pos(unix_time, that_ecef, that);
    ecef2eci_pos(unix_time, phat_ecef, phat);
  }

  return s;
} /* eulersat_sph_basis() */

/*
eulersat_basis()
  Compute satellite-frame basis vectors in Cartesian (ECI)
coordinates, defined as:

s1 = v_t / |v_t|
s2 = s3 x s1
s3 = -rhat

where

v_t = v - (v . s3) s3

Inputs: idx    - index specifying satellite position
        s1_hat - (output) Cartesian ECI components of s1_hat
        s2_hat - (output) Cartesian ECI components of s2_hat
        s3_hat - (output) Cartesian ECI components of s3_hat
        data   - data
*/

static int
eulersat_basis(const size_t idx, double s1_hat[3], double s2_hat[3],
               double s3_hat[3], const satdata_mag *data)
{
  int s = 0;
  double V[3],     /* velocity vector in ECI */
         Vt[3];    /* tangent velocity vector in ECI */
  double V_3;      /* V . s3_hat */
  size_t j;
  time_t unix_time = satdata_epoch2timet(data->t[idx]);
  double theta = M_PI / 2.0 - data->latitude[idx] * M_PI / 180.0;
  double phi = data->longitude[idx] * M_PI / 180.0;

  /* get ECI velocity */
  eulersat_vel(idx, V, data);

#if 1

  /* s3 = -e_mu */
  {
    double r = data->R + data->altitude[idx];
    double e_mu[3], e_nu[3], e_phi[3]; /* ECEF components of ellipsoidal basis vectors */
    double e_mu_eci[3];                /* ECI components of e_mu */
    double r_ECEF[3];

    /* compute ECEF position */
    r_ECEF[0] = r * sin(theta) * cos(phi);
    r_ECEF[1] = r * sin(theta) * sin(phi);
    r_ECEF[2] = r * cos(theta);

    /* compute ellipsoidal basis vectors at time t */
    ellipsoid_basis_mu(r_ECEF, WGS84_MU, e_mu, e_nu, e_phi);

    /* transform e_mu from ECEF to ECI */
    ecef2eci_pos(unix_time, e_mu, e_mu_eci);

    /* compute s3 = -e_mu */
    for (j = 0; j < 3; ++j)
      s3_hat[j] = -e_mu_eci[j];
  }

#else

  /* s3 = -rhat */
  {
    double rhat_ecef[3], that_ecef[3], phat_ecef[3];
    double rhat_eci[3];

    /* compute spherical basis vectors in ECEF */
    sph_basis(theta, phi, rhat_ecef, that_ecef, phat_ecef);

    /* convert rhat to ECI */
    ecef2eci_pos(unix_time, rhat_ecef, rhat_eci);

    /* compute s3 = -rhat */
    for (j = 0; j < 3; ++j)
      s3_hat[j] = -rhat_eci[j];
  }

#endif

  /* compute V_3 component of V */
  V_3 = vec_dot(V, s3_hat);

  /* compute V_t = V - V_3 s3_hat */
  for (j = 0; j < 3; ++j)
    Vt[j] = V[j] - V_3 * s3_hat[j];

  /* s1 = V_t / |V_t| */
  for (j = 0; j < 3; ++j)
    s1_hat[j] = Vt[j] / vec_norm(Vt);

  /* s2 = s3 x s1 */
  sphcross(s3_hat, s1_hat, s2_hat);

  return s;
} /* eulersat_basis() */

/*
eulersat_Rs()
  Explicitely construct R_s rotation matrix
to take vector from s1,s2,s3 basis to NEC

Inputs: idx  - index specifying satellite position
        Rs   - (output) rotation matrix R_s
        data - satellite data
*/

static int
eulersat_Rs(const size_t idx, gsl_matrix *Rs,
            const satdata_mag *data)
{
  double s1_hat[3], s2_hat[3], s3_hat[3];
  double rhat_ecef[3], that_ecef[3], phat_ecef[3];
  double rhat_eci[3], that_eci[3], phat_eci[3];
  time_t unix_time = satdata_epoch2timet(data->t[idx]);
  double theta = M_PI / 2.0 - data->latitude[idx] * M_PI / 180.0;
  double phi = data->longitude[idx] * M_PI / 180.0;

  /* compute s1,s2,s3 in ECI */
  eulersat_basis(idx, s1_hat, s2_hat, s3_hat, data);

  /* compute spherical basis vectors in ECEF */
  sph_basis(theta, phi, rhat_ecef, that_ecef, phat_ecef);

  /* convert spherical basis vectors to ECI */
  ecef2eci_pos(unix_time, rhat_ecef, rhat_eci);
  ecef2eci_pos(unix_time, that_ecef, that_eci);
  ecef2eci_pos(unix_time, phat_ecef, phat_eci);

  gsl_matrix_set(Rs, 0, 0, -vec_dot(s1_hat, that_eci));
  gsl_matrix_set(Rs, 0, 1, -vec_dot(s2_hat, that_eci));
  gsl_matrix_set(Rs, 0, 2, -vec_dot(s3_hat, that_eci));

  gsl_matrix_set(Rs, 1, 0, vec_dot(s1_hat, phat_eci));
  gsl_matrix_set(Rs, 1, 1, vec_dot(s2_hat, phat_eci));
  gsl_matrix_set(Rs, 1, 2, vec_dot(s3_hat, phat_eci));

  gsl_matrix_set(Rs, 2, 0, -vec_dot(s1_hat, rhat_eci));
  gsl_matrix_set(Rs, 2, 1, -vec_dot(s2_hat, rhat_eci));
  gsl_matrix_set(Rs, 2, 2, -vec_dot(s3_hat, rhat_eci));

  return 0;
} /* eulersat_Rs() */

static void
eulersat_print_state (size_t iter, gsl_multifit_fdfsolver * s)
{
  fprintf(stderr,
          "iter: %3zu\n"
          "E = %15.8f %15.8f %15.8f [deg]\n"
          "|f(x)| = %g\n",
          iter,
          gsl_vector_get (s->x, EULER_IDX_ALPHA) * 180.0 / M_PI,
          gsl_vector_get (s->x, EULER_IDX_BETA) * 180.0 / M_PI,
          gsl_vector_get (s->x, EULER_IDX_GAMMA) * 180.0 / M_PI,
          gsl_blas_dnrm2 (s->f));
}

static int
eulersat_print_residuals(char *filename, const size_t flags,
                         const gsl_vector *m, const satdata_mag *data)
{
  int s = 0;
  FILE *fp;
  size_t i;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "eulersat_print_residuals: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Field %zu: time (years)\n", i++);
  fprintf(fp, "# Field %zu: altitude (km)\n", i++);
  fprintf(fp, "# Field %zu: latitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: longitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: QD latitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: X residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: Y residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: Z residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: satellite direction (+/- 1)\n", i++);

  for (i = 0; i < data->n; ++i)
    {
      double B_nec[3], B_main[3];
      double B_nec_orig[3];

      if (data->flags[i])
        continue;

      B_main[0] = SATDATA_VEC_X(data->B_main, i);
      B_main[1] = SATDATA_VEC_Y(data->B_main, i);
      B_main[2] = SATDATA_VEC_Z(data->B_main, i);

      B_nec_orig[0] = SATDATA_VEC_X(data->B, i);
      B_nec_orig[1] = SATDATA_VEC_Y(data->B, i);
      B_nec_orig[2] = SATDATA_VEC_Z(data->B, i);

      s = eulersat_vfm2nec(flags, i, m, B_nec, data);

      fprintf(fp, "%f %f %f %f %f %.12e %.12e %.12e %.12e %.12e %.12e %d\n",
              satdata_epoch2year(data->t[i]),
              data->altitude[i],
              data->latitude[i],
              data->longitude[i],
              data->qdlat[i],
              B_nec[0] - B_main[0],
              B_nec[1] - B_main[1],
              B_nec[2] - B_main[2],
              B_nec_orig[0] - B_main[0],
              B_nec_orig[1] - B_main[1],
              B_nec_orig[2] - B_main[2],
              satdata_mag_satdir(i, data));
    }

  fclose(fp);

  return s;
}
