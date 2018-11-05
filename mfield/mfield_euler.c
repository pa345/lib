/*
 * mfield_euler.c
 *
 * This module contains routines for fitting Euler angles
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>

#include <common/common.h>
#include <bspline2/gsl_bspline2.h>

#include "mfield_euler.h"

static size_t mfield_euler_bin(const size_t sat_idx, const double t,
                               const mfield_workspace *w);

/*
mfield_euler_idx()
  Return index of Euler angles for a given satellite and time period

Inputs: sat_idx - satellite index in [0,nsat-1]
        t       - time (CDF_EPOCH)
        w       - workspace

Return: index so that:
alpha = w->c[index]
beta = w->c[index + 1]
gamma = w->c[index + 2]

Notes:
1) The Euler angles are organized inside the coefficient vector as:

c = [ MF | SV | SA | Euler | Ext ]

and

Euler = [ Euler_0 | Euler_1 | ... | Euler_n ]

where Euler_k are the Euler angles for satellite k. Now

Euler_k = [ Euler_{k0} | Euler_{k1} | ... | Euler_{km} ]

and m = nbins - 1. The appropriate bin is selected according to
the time, as Euler angles are computed for regular time intervals,
specified by the euler_period parameter. Finally, for satellite k and
bin l, the Euler angles are:

Euler_{kl} = [ alpha beta gamma ]
*/

size_t
mfield_euler_idx(const size_t sat_idx, const double t, const mfield_workspace *w)
{
  /*size_t idx = w->euler_offset + 3 * sat_idx;*/
  size_t bin = mfield_euler_bin(sat_idx, t, w);
  size_t idx = w->euler_offset + w->offset_euler[sat_idx] + 3 * bin;

  assert(idx >= w->euler_offset && idx < w->euler_offset + w->p_euler);

  return idx;
} /* mfield_euler_idx() */

int
mfield_euler_print(const char *filename, const size_t sat_idx,
                   const mfield_workspace *w)
{
  FILE *fp;
  gsl_bspline2_workspace *spline_p = w->euler_spline_workspace_p[CIDX2(sat_idx, w->nsat, 0, w->max_threads)];
  const size_t ncontrol = gsl_bspline2_ncontrol(spline_p);
  const size_t euler_idx = w->euler_offset + w->offset_euler[sat_idx];
  const double dt = 5.0 * 86400000; /* 5 days in ms */

  double euler_data[EULER_P];
  gsl_vector_view euler_params = gsl_vector_view_array(euler_data, EULER_P);

  /* Euler angle control points for this dataset */
  gsl_vector_const_view tmp = gsl_vector_const_subvector(w->c, euler_idx, EULER_P * ncontrol);
  gsl_matrix_const_view control_pts = gsl_matrix_const_view_vector(&tmp.vector, EULER_P, ncontrol);

  double t0 = w->data_workspace_p->t0[sat_idx];
  double t1 = w->data_workspace_p->t1[sat_idx];
  double t;
  size_t i;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "mfield_euler_print: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Euler angles for satellite %zu\n", sat_idx);
  fprintf(fp, "# Field %zu: timestamp (CDF_EPOCH)\n", i++);
  fprintf(fp, "# Field %zu: time (decimal years)\n", i++);
  fprintf(fp, "# Field %zu: alpha (degrees)\n", i++);
  fprintf(fp, "# Field %zu: beta (degrees)\n", i++);
  fprintf(fp, "# Field %zu: gamma (degrees)\n", i++);

  for (t = t0; t < t1; t += dt)
    {
      gsl_bspline2_vector_eval(t, &control_pts.matrix, &euler_params.vector, spline_p);

      fprintf(fp, "%f %f %.12e %.12e %.12e\n",
              t,
              satdata_epoch2year(t),
              wrap180(euler_data[0] * 180.0 / M_PI),
              wrap180(euler_data[1] * 180.0 / M_PI),
              wrap180(euler_data[2] * 180.0 / M_PI));
    }

  fclose(fp);

  return 0;
}

/*
mfield_euler_bin()
  Return bin number of Euler angles for a given satellite and time

Inputs: sat_idx - satellite index in [0,nsat-1]
        t       - time (CDF_EPOCH)
        w       - workspace

Return: bin number corresponding to Euler angle triplet
*/

static size_t
mfield_euler_bin(const size_t sat_idx, const double t, const mfield_workspace *w)
{
  double t0 = w->data_workspace_p->t0[sat_idx];
  double t1 = w->data_workspace_p->t1[sat_idx];
  double ratio = (t - t0) / (t1 - t0);
  size_t bin = (size_t) (w->nbins_euler[sat_idx] * ratio);

  if (bin >= w->nbins_euler[sat_idx])
    bin = w->nbins_euler[sat_idx] - 1;

  return bin;
}
