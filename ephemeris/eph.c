/*
 * eph.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp.h>

#include <common/common.h>
#include <common/eci.h>

#include "eph.h"
#include "eph_data.h"
#include "hermite.h"

eph_workspace *
eph_alloc(const eph_data *data)
{
  eph_workspace *w;
  size_t i;

  w = calloc(1, sizeof(eph_workspace));
  if (!w)
    return 0;

  w->data = data;
  w->degree = 3;
  w->n = data->n;

  w->hermite_x = hermite_alloc(w->n, w->degree);
  w->hermite_y = hermite_alloc(w->n, w->degree);
  w->hermite_z = hermite_alloc(w->n, w->degree);
  w->acc = gsl_interp_accel_alloc();
  if (!w->hermite_x || !w->hermite_y || !w->hermite_z || !w->acc)
    {
      fprintf(stderr, "eph_alloc: error allocating hermite interpolation: %s\n", 
              strerror(errno));
      eph_free(w);
      return 0;
    }

  w->t = malloc(w->n * sizeof(double));
  if (!w->t)
    {
      fprintf(stderr, "eph_alloc: error allocating timestamp array: %s\n", 
              strerror(errno));
      eph_free(w);
      return 0;
    }

  /* convert ephemeris time to s since velocities are in km/s */
  for (i = 0; i < w->n; ++i)
    w->t[i] = data->t[i] / 1000.0;

  /* initialize interpolation */
  hermite_init(w->t, data->X, data->VX, data->n, w->hermite_x);
  hermite_init(w->t, data->Y, data->VY, data->n, w->hermite_y);
  hermite_init(w->t, data->Z, data->VZ, data->n, w->hermite_z);

  return w;
} /* eph_alloc() */

void
eph_free(eph_workspace *w)
{
  if (w->hermite_x)
    hermite_free(w->hermite_x);

  if (w->hermite_y)
    hermite_free(w->hermite_y);

  if (w->hermite_z)
    hermite_free(w->hermite_z);

  if (w->acc)
    gsl_interp_accel_free(w->acc);

  if (w->t)
    free(w->t);

  free(w);
} /* eph_free() */

/*
eph_interp()
  Interpolate ephemeris to a given time t

Inputs: t - timestamp (CDF_EPOCH)
        r - (output) X,Y,Z of interpolated time t (ECI or ECEF)
        v - (output) VX,VY,VZ of interpolated time t (ECI or ECEF)
        w - workspace

Return: success/error
*/

int
eph_interp(const double t, double r[3], double v[3], eph_workspace *w)
{
  int s = 0;
  const eph_data *data = w->data;
  const size_t n = data->n;
  const double t_sec = t / 1000.0; /* convert to sec */

  if (t < data->t[0] || t > data->t[n - 1])
    {
      fprintf(stderr, "eph_interp: t outside allowed range\n");
      return -1;
    }

  r[0] = hermite_eval(w->t, data->X, data->VX, t_sec, w->acc, w->hermite_x);
  r[1] = hermite_eval(w->t, data->Y, data->VY, t_sec, w->acc, w->hermite_y);
  r[2] = hermite_eval(w->t, data->Z, data->VZ, t_sec, w->acc, w->hermite_z);

  v[0] = hermite_eval_deriv(w->t, data->X, data->VX, t_sec, w->acc, w->hermite_x);
  v[1] = hermite_eval_deriv(w->t, data->Y, data->VY, t_sec, w->acc, w->hermite_y);
  v[2] = hermite_eval_deriv(w->t, data->Z, data->VZ, t_sec, w->acc, w->hermite_z);

  return s;
}

/*
eph_interp_sph()
  Interpolate ephemeris to a given time t; return result as ECEF
(r,theta,phi)

Inputs: t     - timestamp (CDF_EPOCH)
        r_sph - (output) position of interpolated time t
                r_sph[0] = geocentric radius (km)
                r_sph[1] = geocentric co-latitude (rad)
                r_sph[2] = geocentric longitude (rad)
        w     - workspace

Return: success/error
*/

int
eph_interp_sph(const double t, double r_sph[3], eph_workspace *w)
{
  int s = 0;
  double pos[3]; /* position (X,Y,Z) (ECI or ECEF) */
  double vel[3]; /* velocity (X,Y,Z) (ECI or ECEF) */

  s = eph_interp(t, pos, vel, w);
  if (s)
    return s;

  if (w->data->flags & EPH_DATA_FLG_ECEF)
    {
      /* position is ECEF */
      r_sph[0] = gsl_hypot3(pos[0], pos[1], pos[2]);
      r_sph[1] = acos(pos[2] / r_sph[0]);
      r_sph[2] = atan2(pos[1], pos[0]);
    }
  else if (w->data->flags & EPH_DATA_FLG_ECI)
    {
      /* position is ECI */
      time_t unix_time = epoch2timet(t);
      eci2sph_pos(unix_time, pos, r_sph);
    }
  else
    {
      fprintf(stderr, "eph_interp_sh: error: unknown ephemeris type\n");
      return -1;
    }

  return GSL_SUCCESS;
}
