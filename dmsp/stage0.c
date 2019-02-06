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
#include <gsl/gsl_test.h>

#include <indices/indices.h>
#include <satdata/satdata.h>

#include <common/common.h>
#include <common/ellipsoid.h>
#include <common/quat.h>
#include <common/eci.h>
#include <common/julian.h>

#include "eph.h"
#include "eph_data.h"
#include "hermite.h"

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
      SATDATA_VEC_X(data->B_VFM, n) = Y; /* velocity direction */
      SATDATA_VEC_Y(data->B_VFM, n) = Z; /* orbit normal */
      SATDATA_VEC_Z(data->B_VFM, n) = X; /* down */

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

/*
calc_spacecraft_basis_ECI()
  Calculate spacecraft basis vectors to rotate a vector from the spacecraft (S/C) frame
to NEC. We define spacecraft-fixed basis vectors as

s1 = s2 x s3                   (velocity direction)
s2 = (s3 x v) / | s3 x v |     (negative orbit normal direction)
s3 = -e_mu                     (local geodetic downward)

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

#if 0
  {
    double r_ECEF[3], tmp[3];

    /* convert ECI position to ECEF */
    eci2ecef(t, r_ECI, r_ECEF);

    /* compute tmp = e_mu in ECEF */
    ellipsoid_basis(r_ECEF, tmp, s1, s2);
    /*ellipsoid_basis_mu(r_ECEF, WGS84_MU, tmp, s1, s2);*/

    /* compute s3 = e_mu in ECI */
    ecef2eci(t, tmp, s3);

    /* compute s3 = -e_mu in ECI */
    for (i = 0; i < 3; ++i)
      s3[i] *= -1.0;
  }

#elif 1
  {
    double r_ECEF[3], tmp[3];

    /* convert ECI position to ECEF */
    eci2ecef(t, r_ECI, r_ECEF);

    /* compute tmp = e_mu in ECEF */
    ellipsoid_basis_mu(r_ECEF, WGS84_MU, tmp, s1, s2);

    /* compute s3 = e_mu in ECI */
    ecef2eci(t, tmp, s3);

    /* compute s3 = -e_mu in ECI */
    for (i = 0; i < 3; ++i)
      s3[i] *= -1.0;
  }

#else

  /* store -rhat in s3 */
  for (i = 0; i < 3; ++i)
    s3[i] = -r_ECI[i];

#endif
    
  vec_unit(s3, s3);

  for (i = 0; i < 3; ++i)
    v[i] = v_ECI[i];

  vec_unit(v, v);

  /* s2 = (s3 x v) / | s3 x v | = minus orbit normal */
  vec_cross(s3, v, s2);
  vec_unit(s2, s2);

  /* s1 = s2 x s3 */
  vec_cross(s2, s3, s1);
  vec_unit(s1, s1);

  return s;
}

/*
calc_quaternions_ECI()
  Calculate quaternions which allow a transformation from spacecraft frame
to NEC:

B_NEC = R_q B_LVLH

where LVLH is the local-vertical local-horizontal reference frame, defined
by the (s1,s2,s3) basis:

s1 = s2 x s3                   (velocity direction)
s2 = (s3 x v) / | s3 x v |     (negative orbit normal direction)
s3 = -e_mu                     (local geodetic downward)

Inputs: t - timestamp
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
  quat_R2q(&Rq.matrix, q);

  return 0;
}

int
calc_quaternions_ECI2(const time_t t, const double theta, const double phi,
                     const double r_ECI[3], const double v_ECI[3], double q[4])
{
  const double jd = (t / 86400.0) + 2440587.5;
  const double GAST = julian2GAST(jd);
  const double T = GAST + phi;
  double r_ECEF[3];                                /* position vector in ECEF frame */
  double xhat_ECEF[3], yhat_ECEF[3], zhat_ECEF[3]; /* NEC basis vectors in ECEF frame */
  double xhat_ECI[3], yhat_ECI[3], zhat_ECI[3];    /* NEC basis vectors in ECI frame */
  double s1[3], s2[3], s3[3];                      /* spacecraft-fixed basis vectors in ECI frame */
  double DCM_data[9];
  gsl_matrix_view DCM = gsl_matrix_view_array(DCM_data, 3, 3);
  size_t i;

  /* compute position vector in ECEF frame */
  eci2ecef(t, r_ECI, r_ECEF);

  /* compute NEC basis vectors in ECEF frame */
  ecef2nec_basis(r_ECEF, xhat_ECEF, yhat_ECEF, zhat_ECEF);

  /* transform NEC basis vectors to ECI components */
  ecef2eci(t, xhat_ECEF, xhat_ECI);
  ecef2eci(t, yhat_ECEF, yhat_ECI);
  ecef2eci(t, zhat_ECEF, zhat_ECI);

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
  quat_R2q(&DCM.matrix, q);

  return 0;
}

int
interp_eph(satdata_mag *data, eph_data *eph)
{
  int s = 0;
  size_t i, j;
  eph_workspace *w = eph_alloc(eph);

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
          eci2sph_pos(unix_time, pos, r_sph);

          r = r_sph[0];
          theta = r_sph[1];
          phi = r_sph[2];

#if 0
          calc_quaternions_ECI(unix_time, theta, phi, pos, vel, q);
#else
          calc_quaternions_ECI2(unix_time, theta, phi, pos, vel, q);
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

  satdata_mag_free(data);
  satdata_mag_free(data_final);
  if (eph)
    eph_data_free(eph);

  return 0;
}
