/*
 * stage0.c
 *
 * Convert ASCII DMSP data files into CDF format with possible fixing
 * of ephemeris values
 *
 * Usage: ./stage0 <-i input_ascii_gz_file> [-o output_cdf_file]
 *                 [-n nasa_ephemeris_file] [-s sgp4_ephemeris_file]
 *                 [-b bowman_ephemeris_file]
 *
 * The measurements are read from the ascii file, ephemeris values are
 * possibly modified (lat,lon,alt) if needed, and output is written to CDF
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>
#include <sys/time.h>
#include <zlib.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_test.h>

#include <mainlib/ml_indices.h>
#include <mainlib/ml_satdata.h>
#include <mainlib/ml_att.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_ellipsoid.h>
#include <mainlib/ml_eci.h>
#include <mainlib/ml_julian.h>
#include <mainlib/ml_sofa.h>
#include <mainlib/ml_spice.h>
#include <mainlib/ml_eph.h>
#include <mainlib/ml_eph_data.h>
#include <mainlib/ml_hermite.h>

/*
 * 0 = old/original Greenwich only transformation
 * 1 = libspice
 * 2 = libsofa
 */
#define ECI_LIB       0

/*
 * DMSP ellipsoid semi-major and semi-minor axes (km);
 * inverse flattening = 298.25
 *
 * These are the WGS66 ellipsoid parameters
 */
#define DMSP_A        (6378.145)
#define DMSP_B        (6356.75976948868)

typedef struct
{
  double E;  /* 1 / AX^2 */
  double F;  /* 1 / AY^2 */
  double G;  /* 1 / B^2 */
  double X;  /* ECEF X position (km) */
  double Y;  /* ECEF Y position (km) */
  double Z;  /* ECEF Z position (km) */
} ellipsoid_params;

/* read DMSP SSM Fred Rich ASCII data file */
size_t
dmsp_read_MFR(const char *filename, satdata_mag *data)
{
  const double R = R_EARTH_KM;
  char buf[SATDATA_MAX_BUF];
  size_t n = 0;
  gzFile fp;

  fp = gzopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "unable to open %s: %s\n", filename, strerror(errno));
      return 0;
    }

  while (gzgets(fp, buf, SATDATA_MAX_BUF) != 0)
    {
      int c;
      double date, sec, lat, lon, alt, X, Y, Z;
      double caldate, resX, resY, resZ;
      char as[10], istr[10], dstr[10];
      int mid;
      int year, month, day, doy;
      int hh, mm, ss, msec;

      if (*buf == '#')
        continue;

      c = sscanf(buf, "%lf %lf %d %s %lf %lf %lf %s %s %lf %lf %lf %lf %lf %lf %lf",
                 &date,
                 &sec,
                 &mid,
                 as,
                 &lat,
                 &lon,
                 &alt,
                 istr,
                 dstr,
                 &X,
                 &Y,
                 &Z,
                 &caldate,
                 &resX,
                 &resY,
                 &resZ);
      if (c < 16)
        continue;

      year = (int) (date / 1.0e9);
      doy = (int) ((date - year * 1.0e9) / 1.0e6);

      hh = (int) (sec / 3600.0);
      mm = (int) ((sec - hh * 3600.0) / 60.0);
      ss = (int) (sec - hh * 3600.0 - mm * 60.0);
      msec = (int) ((sec - (int) sec) * 1000.0);

      doy2md(year, doy, &month, &day);

      data->t[n] = computeEPOCH(year, month, day, hh, mm, ss, msec);
      data->latitude[n] = lat;
      data->longitude[n] = lon;
      data->r[n] = R + alt;

      /*
       * the (X,Y,Z) in the DMSP ASCII files are not NEC so they are
       * stored differently below
       */
      SATDATA_VEC_X(data->B_VFM, n) = X; /* down */
      SATDATA_VEC_Y(data->B_VFM, n) = Y; /* +/- velocity */
      SATDATA_VEC_Z(data->B_VFM, n) = Z; /* orbit normal, towards sun */

      data->F[n] = gsl_hypot3(X, Y, Z);

      if (++n >= data->ntot)
        {
          fprintf(stderr, "dmsp_read_MFR: file %s contains too many data records\n",
                  filename);
          return n;
        }
    }

  gzclose(fp);

  data->n = n;
  data->R = R;

  return n;
}

int
ellipsoid_f(const gsl_vector * x, void * params, gsl_vector * f)
{
  const ellipsoid_params * p = (const ellipsoid_params *) params;
  const double E = p->E;
  const double F = p->F;
  const double G = p->G;
  const double XG = p->X;                 /* fixed point */
  const double YG = p->Y;
  const double ZG = p->Z;
  const double XE = gsl_vector_get(x, 0); /* projected point on ellipsoid */
  const double YE = gsl_vector_get(x, 1);
  const double ZE = gsl_vector_get(x, 2);

  gsl_vector_set(f, 0, (XE - XG) * G * ZE - (ZE - ZG) * E * XE);
  gsl_vector_set(f, 1, (YE - YG) * G * ZE - (ZE - ZG) * F * YE);
  gsl_vector_set(f, 2, E * XE * XE + F * YE * YE + G * ZE * ZE - 1.0);

  return GSL_SUCCESS;
}

int
ellipsoid_df(const gsl_vector * x, void * params, gsl_matrix * J)
{
  const ellipsoid_params * p = (const ellipsoid_params *) params;
  const double E = p->E;
  const double F = p->F;
  const double G = p->G;
  const double XG = p->X;                 /* fixed point */
  const double YG = p->Y;
  const double ZG = p->Z;
  const double XE = gsl_vector_get(x, 0); /* projected point on ellipsoid */
  const double YE = gsl_vector_get(x, 1);
  const double ZE = gsl_vector_get(x, 2);

  gsl_matrix_set(J, 0, 0, G * ZE - (ZE - ZG) * E);
  gsl_matrix_set(J, 0, 1, 0.0);
  gsl_matrix_set(J, 0, 2, (XE - XG) * G - E * XE);

  gsl_matrix_set(J, 1, 0, 0.0);
  gsl_matrix_set(J, 1, 1, G * ZE - (ZE - ZG) * F);
  gsl_matrix_set(J, 1, 2, (YE - YG) * G - F * YE);

  gsl_matrix_set(J, 2, 0, 2.0 * E * XE);
  gsl_matrix_set(J, 2, 1, 2.0 * F * YE);
  gsl_matrix_set(J, 2, 2, 2.0 * G * ZE);

  return GSL_SUCCESS;
}

int
ellipsoid_fdf(const gsl_vector * x, void * params, gsl_vector * f, gsl_matrix * J)
{
  ellipsoid_f(x, params, f);
  ellipsoid_df(x, params, J);
  return GSL_SUCCESS;
}

void
print_state(const size_t iter, gsl_multiroot_fdfsolver * s)
{
  fprintf(stderr, "iter = %zu: X = %.4f Y = %.4f Z = %.4f\n",
          iter,
          gsl_vector_get(s->x, 0),
          gsl_vector_get(s->x, 1),
          gsl_vector_get(s->x, 2));
}

/*
calc_local_normal()
  The DMSP X spacecraft axis is defined as:
  
"a line through the spacecraft and normal to the Earth ellipsoid, positive toward Earth"

So, for a given spacecraft position, we need to project that position vector onto
the ellipsoid. This is done by solving a set of nonlinear equations, as detailed in:

M. Ligas, "Cartesian to geodetic coordinates conversion on a triaxial ellipsoid",
J. Geod, 86, 2012.

This routine uses "Case 3" from that paper.

Inputs: AX     - semi-major axis of ellipsoid in X (km)
        AY     - semi-minor axis of ellipsoid in Y (km)
        B      - semi-minor axis of ellipsoid in Z (km)
        r_ECEF - spacecraft position vector in ECEF (km)
        normal - (output) outward normal to ellipsoid in ECEF

Return: success/error
*/

int
calc_local_normal(const double AX, const double AY, const double B, const double r_ECEF[3],
                  double normal[3])
{
  int status = 0;
  const gsl_multiroot_fdfsolver_type * T = gsl_multiroot_fdfsolver_gnewton;
  const size_t n = 3;
  const double r = gsl_hypot3(r_ECEF[0], r_ECEF[1], r_ECEF[2]);
  const double tol = 1.0e-10;
  const size_t max_iter = 100;
  size_t iter = 0;
  ellipsoid_params params;
  gsl_multiroot_function_fdf f;
  gsl_multiroot_fdfsolver * solver = gsl_multiroot_fdfsolver_alloc(T, n);
  gsl_vector_view x = gsl_vector_view_array(normal, 3);

  params.E = (1.0 / AX) / AX;;
  params.F = (1.0 / AY) / AY;;
  params.G = (1.0 / B) / B;;
  params.X = r_ECEF[0];
  params.Y = r_ECEF[1];
  params.Z = r_ECEF[2];

  f.f = ellipsoid_f;
  f.df = ellipsoid_df;
  f.fdf = ellipsoid_fdf;
  f.params = &params;
  f.n = n;

  /* initial conditions */
  normal[0] = AX * r_ECEF[0] / r;
  normal[1] = AY * r_ECEF[1] / r;
  normal[2] = B * r_ECEF[2] / r;

  gsl_multiroot_fdfsolver_set(solver, &f, &x.vector);

  /*print_state(iter, solver);*/

  do
    {
      iter++;

      status = gsl_multiroot_fdfsolver_iterate(solver);

      /*print_state(iter, solver);*/

      if (status)
        break;

      status = gsl_multiroot_test_residual(solver->f, tol);
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  /*
   * compute normal vector,
   *
   * n = 2 [ XE / AX^2 ; YE / AY^2 ; ZE / B^2 ]
   */
  normal[0] = 2.0 * params.E * gsl_vector_get(solver->x, 0);
  normal[1] = 2.0 * params.F * gsl_vector_get(solver->x, 1);
  normal[2] = 2.0 * params.G * gsl_vector_get(solver->x, 2);

  /* test that longitude of both points is the same */
  {
    double phi0 = atan2(r_ECEF[1], r_ECEF[0]);
    double phi = atan2(gsl_vector_get(solver->x, 1), gsl_vector_get(solver->x, 0));
    gsl_test_rel(phi, phi0, 1.0e-12, "phi");
  }

  gsl_multiroot_fdfsolver_free(solver);

  return status;
}

/*
calc_spacecraft_basis_ECI()
  Calculate spacecraft basis vectors to rotate a vector from the spacecraft (S/C) frame
to NEC. We define spacecraft-fixed basis vectors as

s1 = -e_mu                       (local geodetic downward)
s2 = s3 x s1                     (+/- velocity direction)
s3 = +/- (s1 x v) / || s1 x v || (orbit normal, positive toward sun)

Inputs: r_ECI  - position (X,Y,Z) in Cartesian ECI (km)
        v_ECI  - velocity (VX,VY,VZ) in Cartesian ECI (km/s)
        s1     - (output) s1 unit vector (ECI)
        s2     - (output) s2 unit vector (ECI)
        s3     - (output) s3 unit vector (ECI)

Notes:
1) The spacecraft-fixed frame is defined as local-vertical local-horizontal (LVLH)
by the APSM DMSP processing problem. This frame is given by the (s1,s2,s3) basis
vectors defined above.

2) Additional rotation matrices are needed to go from NEC to LVLH as defined below.
I don't fully understand these matrices yet (Jan 24 2019)

B_LVLH = A B1 B2 B_NEC

A = [ s1 ]
    [ s2 ]
    [ s3 ]

B1 = [ cos(T) -sin(T) 0 ]
     [ sin(T)  cos(T) 0 ]
     [   0       0    1 ]

B2 = [ -cos(theta) 0 -sin(theta) ]
     [      0      1       0     ]
     [  sin(theta) 0 -cos(theta) ]

where:

theta is geocentric co-latitude
T = GAST + phi
GAST is Greenwich apparent sidereal time
phi is geocentric longitude
*/

static int
calc_spacecraft_basis_ECI(const time_t t, const double r_ECI[3], const double v_ECI[3],
                          double s1[3], double s2[3], double s3[3])
{
  int s = 0;
  double v[3]; /* unit velocity vector */
  size_t i;

  {
    double r_ECEF[3], tmp[3];
    double R_data[9];
    gsl_matrix_view R = gsl_matrix_view_array(R_data, 3, 3);

    /* convert ECI position to ECEF */
#if ECI_LIB == 0
    eci2ecef(t, r_ECI, r_ECEF);
#elif ECI_LIB == 1
    {
      gsl_vector_const_view in = gsl_vector_const_view_array(r_ECI, 3);
      gsl_vector_view out = gsl_vector_view_array(r_ECEF, 3);
      spice_c2t(t, &R.matrix);
      gsl_blas_dgemv(CblasNoTrans, 1.0, &R.matrix, &in.vector, 0.0, &out.vector);
    }
#else
#endif

    /* compute tmp = e_mu in ECEF */
#if 0
    ellipsoid_basis_geoid(DMSP_A, DMSP_B, r_ECEF, tmp, s2, s3);
#else
    calc_local_normal(DMSP_A, DMSP_A, DMSP_B, r_ECEF, tmp);
#endif

    /* compute s1 = e_mu in ECI */
#if ECI_LIB == 0
    ecef2eci(t, tmp, s1);
#elif ECI_LIB == 1
    {
      gsl_vector_const_view in = gsl_vector_const_view_array(tmp, 3);
      gsl_vector_view out = gsl_vector_view_array(s1, 3);
      gsl_blas_dgemv(CblasTrans, 1.0, &R.matrix, &in.vector, 0.0, &out.vector);
    }
#else
#endif

    /* compute s1 = -e_mu in ECI */
    for (i = 0; i < 3; ++i)
      s1[i] *= -1.0;
  }

  vec_unit(s1, s1);

  for (i = 0; i < 3; ++i)
    v[i] = v_ECI[i];

  vec_unit(v, v);

  /* s3 = (s1 x v) / | s1 x v | = orbit normal, positive toward sun */
  vec_cross(s1, v, s3);
  vec_unit(s3, s3);

  /* s2 = s3 x s1 */
  vec_cross(s3, s1, s2);
  vec_unit(s2, s2);

  return s;
}

/*
calc_quaternions_ECI()
  Calculate quaternions which allow a transformation from spacecraft frame
to NEC:

B_NEC = R_q B_LVLH

where LVLH is the local-vertical local-horizontal reference frame, defined
by the (s1,s2,s3) basis.

B_LVLH = A B1 B2 B_NEC

A = [ s1 ]
    [ s2 ]
    [ s3 ]

B1 = [ cos(T) -sin(T) 0 ]
     [ sin(T)  cos(T) 0 ]
     [   0       0    1 ]

B2 = [ -cos(theta) 0 -sin(theta) ]
     [      0      1       0     ]
     [  sin(theta) 0 -cos(theta) ]

where:

theta is geocentric co-latitude
T = GAST + phi
GAST is Greenwich apparent sidereal time
phi is geocentric longitude

Inputs: t     - timestamp
        theta - geocentric co-latitude (radians)
        phi   - longitude (radians)
        r_ECI - ECI position vector (km)
        v_ECI - ECI velocity vector (km/s)
        q     - (output) quaternions for rotation from LVLH to NEC
*/

int
calc_quaternions_ECI(const time_t t, const double theta, const double phi,
                     const double r_ECI[3], const double v_ECI[3], double q[4])
{
  const double jd = (t / 86400.0) + 2440587.5;
  const double GAST = julian2GAST(jd);
  const double T = GAST + phi;
  gsl_vector_view qv = gsl_vector_view_array(q, 4);
  double A_data[9], B1_data[9], B2_data[9], Rq_data[9], tmp_data[9];
  gsl_matrix_view A = gsl_matrix_view_array(A_data, 3, 3);
  gsl_matrix_view B1 = gsl_matrix_view_array(B1_data, 3, 3);
  gsl_matrix_view B2 = gsl_matrix_view_array(B2_data, 3, 3);
  gsl_matrix_view Rq = gsl_matrix_view_array(Rq_data, 3, 3);
  gsl_matrix_view tmp = gsl_matrix_view_array(tmp_data, 3, 3);
  double s1[3], s2[3], s3[3]; /* spacecraft-fixed basis vectors (ECI) */
  size_t i;

  calc_spacecraft_basis_ECI(t, r_ECI, v_ECI, s1, s2, s3);

  gsl_matrix_set_zero(&B1.matrix);
  gsl_matrix_set_zero(&B2.matrix);

  gsl_matrix_set(&B1.matrix, 0, 0, cos(T));
  gsl_matrix_set(&B1.matrix, 0, 1, -sin(T));
  gsl_matrix_set(&B1.matrix, 1, 0, sin(T));
  gsl_matrix_set(&B1.matrix, 1, 1, cos(T));
  gsl_matrix_set(&B1.matrix, 2, 2, 1.0);

  gsl_matrix_set(&B2.matrix, 0, 0, -cos(theta));
  gsl_matrix_set(&B2.matrix, 0, 2, -sin(theta));
  gsl_matrix_set(&B2.matrix, 1, 1, 1.0);
  gsl_matrix_set(&B2.matrix, 2, 0, sin(theta));
  gsl_matrix_set(&B2.matrix, 2, 2, -cos(theta));

  /*
   * A = [ s1 ]
   *     [ s2 ]
   *     [ s3 ]
   */
  for (i = 0; i < 3; ++i)
    {
      gsl_matrix_set(&A.matrix, 0, i, s1[i]);
      gsl_matrix_set(&A.matrix, 1, i, s2[i]);
      gsl_matrix_set(&A.matrix, 2, i, s3[i]);
    }

  /* tmp = B1 * B2 */
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &B1.matrix, &B2.matrix, 0.0, &tmp.matrix);

  /* Rq^T = A * tmp */
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &A.matrix, &tmp.matrix, 0.0, &Rq.matrix);

  /* transpose for final matrix */
  gsl_matrix_transpose(&Rq.matrix);

  /* convert to quaternions */
  att_quat_R2q(&Rq.matrix, &qv.vector);

  return 0;
}

int
calc_quaternions_ECI2(const time_t t, const double r_ECI[3], const double v_ECI[3], double q[4])
{
  double r_ECEF[3];                                /* position vector in ECEF frame */
  double xhat_ECEF[3], yhat_ECEF[3], zhat_ECEF[3]; /* NEC basis vectors in ECEF frame */
  double xhat_ECI[3], yhat_ECI[3], zhat_ECI[3];    /* NEC basis vectors in ECI frame */
  double s1[3], s2[3], s3[3];                      /* spacecraft-fixed basis vectors in ECI frame */
  double DCM_data[9], R_data[9];
  gsl_matrix_view DCM = gsl_matrix_view_array(DCM_data, 3, 3);
  gsl_matrix_view R = gsl_matrix_view_array(R_data, 3, 3);
  gsl_vector_view qv = gsl_vector_view_array(q, 4);

  /* compute position vector in ECEF frame */
#if ECI_LIB == 0
  eci2ecef(t, r_ECI, r_ECEF);
#elif ECI_LIB == 1
  {
    gsl_vector_const_view in = gsl_vector_const_view_array(r_ECI, 3);
    gsl_vector_view out = gsl_vector_view_array(r_ECEF, 3);
    spice_c2t(t, &R.matrix);
    gsl_blas_dgemv(CblasNoTrans, 1.0, &R.matrix, &in.vector, 0.0, &out.vector);
  }
#else
#endif

  /* compute NEC basis vectors in ECEF frame */
  ecef2nec_basis(r_ECEF, xhat_ECEF, yhat_ECEF, zhat_ECEF);

  /* transform NEC basis vectors to ECI components */
#if ECI_LIB == 0
  ecef2eci(t, xhat_ECEF, xhat_ECI);
  ecef2eci(t, yhat_ECEF, yhat_ECI);
  ecef2eci(t, zhat_ECEF, zhat_ECI);
#elif ECI_LIB == 1
  {
    gsl_vector_const_view xin = gsl_vector_const_view_array(xhat_ECEF, 3);
    gsl_vector_const_view yin = gsl_vector_const_view_array(yhat_ECEF, 3);
    gsl_vector_const_view zin = gsl_vector_const_view_array(zhat_ECEF, 3);
    gsl_vector_view xout = gsl_vector_view_array(xhat_ECI, 3);
    gsl_vector_view yout = gsl_vector_view_array(yhat_ECI, 3);
    gsl_vector_view zout = gsl_vector_view_array(zhat_ECI, 3);

    gsl_blas_dgemv(CblasTrans, 1.0, &R.matrix, &xin.vector, 0.0, &xout.vector);
    gsl_blas_dgemv(CblasTrans, 1.0, &R.matrix, &yin.vector, 0.0, &yout.vector);
    gsl_blas_dgemv(CblasTrans, 1.0, &R.matrix, &zin.vector, 0.0, &zout.vector);
  }
#else
#endif

  calc_spacecraft_basis_ECI(t, r_ECI, v_ECI, s1, s2, s3);

  /* build direction cosine matrix */

  gsl_matrix_set(&DCM.matrix, 0, 0, vec_dot(xhat_ECI, s1));
  gsl_matrix_set(&DCM.matrix, 0, 1, vec_dot(xhat_ECI, s2));
  gsl_matrix_set(&DCM.matrix, 0, 2, vec_dot(xhat_ECI, s3));

  gsl_matrix_set(&DCM.matrix, 1, 0, vec_dot(yhat_ECI, s1));
  gsl_matrix_set(&DCM.matrix, 1, 1, vec_dot(yhat_ECI, s2));
  gsl_matrix_set(&DCM.matrix, 1, 2, vec_dot(yhat_ECI, s3));

  gsl_matrix_set(&DCM.matrix, 2, 0, vec_dot(zhat_ECI, s1));
  gsl_matrix_set(&DCM.matrix, 2, 1, vec_dot(zhat_ECI, s2));
  gsl_matrix_set(&DCM.matrix, 2, 2, vec_dot(zhat_ECI, s3));

  /* convert to quaternions */
  att_quat_R2q(&DCM.matrix, &qv.vector);

  return 0;
}

int
interp_eph(satdata_mag *data, eph_data *eph)
{
  int s = 0;
  size_t i, j;
  eph_workspace *w = eph_alloc(eph);
  sofa_workspace * sofa_p = sofa_alloc();

  for (i = 0; i < data->n; ++i)
    {
      double pos[3], vel[3]; /* position and velocity (ECI or ECEF) */
      double r, theta, phi;
      double q[4];           /* quaternions for rotation S/C to NEC */

      /* interpolate ephemeris data to time ti */
      s = eph_interp(data->t[i], pos, vel, w);
      if (s)
        {
          data->flags[i] |= SATDATA_FLG_NOEPH;
          continue;
        }

      if (w->data->flags & EPH_DATA_FLG_ECEF)
        {
          /* compute (r,theta,phi) for this point from ECEF position */
          r = gsl_hypot3(pos[0], pos[1], pos[2]);
          theta = acos(pos[2] / r);
          phi = atan2(pos[1], pos[0]);

          /*XXXcalc_quaternions_ECEF(pos, vel, q);*/
        }
      else if (w->data->flags & EPH_DATA_FLG_ECI)
        {
          time_t unix_time = epoch2timet(data->t[i]);
          double r_sph[3];

          /* compute ECEF (r,theta,phi) from ECI position */
#if ECI_LIB == 0
          eci2sph_pos(unix_time, pos, r_sph);
#elif ECI_LIB == 1
          spice_c2sph(unix_time, pos, r_sph);
#else
          sofa_c2sph(unix_time, pos, r_sph, sofa_p);
#endif

          r = r_sph[0];
          theta = r_sph[1];
          phi = r_sph[2];

#if 0
          calc_quaternions_ECI(unix_time, theta, phi, pos, vel, q);
#elif 1
          calc_quaternions_ECI2(unix_time, pos, vel, q);
#endif
        }
      else
        {
          fprintf(stderr, "interp_eph: unknown interpolation type\n");
          return -1;
        }

      data->r[i] = r;
      data->latitude[i] = 90.0 - theta * 180.0 / M_PI;
      data->longitude[i] = wrap180(phi * 180.0 / M_PI);

      for (j = 0; j < 4; ++j)
        data->q[4*i + j] = q[j];
    }

  eph_free(w);
  sofa_free(sofa_p);

  return 0;
}

int
main(int argc, char *argv[])
{
  char *infile = NULL;
  char *outfile = NULL;
  satdata_mag *data = NULL, *data_final;
  eph_data *eph = NULL;
  int c;
  struct timeval tv0, tv1;
  spice_workspace * spice_p;

  while ((c = getopt(argc, argv, "i:o:b:t:")) != (-1))
    {
      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          case 'b':
            fprintf(stderr, "main: reading Bowman ephemerides from %s...", optarg);
            gettimeofday(&tv0, NULL);
            eph = eph_data_read_bowman(optarg);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu read, %g seconds)\n", eph->n, time_diff(tv0, tv1));
            break;

          case 't':
            fprintf(stderr, "main: reading TENA ephemerides from %s...", optarg);
            gettimeofday(&tv0, NULL);
            eph = eph_data_read_tena(optarg);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu read, %g seconds)\n", eph->n, time_diff(tv0, tv1));
            break;

          case 'o':
            outfile = optarg;
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i DMSP_ascii_gz_file> [-o output_cdf_file] [-b bowman_ephemeris_file] [-t tena_ephemeris_file]\n",
              argv[0]);
      exit(1);
    }

  spice_p = spice_alloc();

  data = satdata_mag_alloc(86400);
  data_final = satdata_mag_alloc(86400);

  fprintf(stderr, "main: reading %s...", infile);
  gettimeofday(&tv0, NULL);
  dmsp_read_MFR(infile, data);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu records read, %g seconds)\n", data->n,
                   time_diff(tv0, tv1));

  if (data->n > 0)
    {
      if (eph)
        {
          fprintf(stderr, "main: interpolating ephemeris of %s...", infile);
          interp_eph(data, eph);
          fprintf(stderr, "done\n");
        }

      fprintf(stderr, "main: copying data...");
      satdata_select_filter(data, data_final);
      fprintf(stderr, "done (%zu data thrown out due to no ephemeris)\n",
              data->n - data_final->n);
 
      if (outfile && data_final->n > 0)
        {
          fprintf(stderr, "main: writing %s...", outfile);
          gettimeofday(&tv0, NULL);
          satdata_dmsp_write(0, outfile, data_final);
          gettimeofday(&tv1, NULL);
          fprintf(stderr, "done (%zu records written, %g seconds)\n", data_final->n,
                  time_diff(tv0, tv1));
        }
    }

  spice_free(spice_p);
  satdata_mag_free(data);
  satdata_mag_free(data_final);
  if (eph)
    eph_data_free(eph);

  return 0;
}
