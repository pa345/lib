/*
 * mfield_align.c
 *
 * This module contains routines for fitting alignment parameters
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

#include <mainlib/ml_common.h>
#include <bspline2/gsl_bspline2.h>

#include "mfield_align.h"

int
mfield_align_print(const char *filename, const size_t sat_idx,
                   const mfield_workspace *w)
{
  FILE *fp;
  const magdata *mptr = mfield_data_ptr(sat_idx, w->data_workspace_p);
  const double *t_year = w->data_workspace_p->t_year[sat_idx];
  gsl_bspline2_workspace *spline_p = w->align_spline_workspace_p[CIDX2(sat_idx, w->nsat, 0, w->max_threads)];
  const size_t ncontrol = gsl_bspline2_ncontrol(spline_p);
  const size_t align_idx = w->align_offset + w->offset_align[sat_idx];
  const double dt = 5.0 / 365.25; /* 5 days in years */

  double align_data[ALIGN_P];
  gsl_vector_view align_params = gsl_vector_view_array(align_data, ALIGN_P);

  /* Euler angle control points for this dataset */
  gsl_vector_const_view tmp = gsl_vector_const_subvector(w->c, align_idx, ALIGN_P * ncontrol);
  gsl_matrix_const_view control_pts = gsl_matrix_const_view_vector(&tmp.vector, ALIGN_P, ncontrol);

  double t0 = t_year[0];
  double t1 = t_year[mptr->n - 1];
  double t;
  size_t i;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "mfield_align_print: unable to open %s: %s\n",
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
      gsl_bspline2_vector_eval(t, &control_pts.matrix, &align_params.vector, spline_p);

      fprintf(fp, "%f %f %.12e %.12e %.12e\n",
              t,
              satdata_epoch2year(t),
              wrap180(align_data[0] * 180.0 / M_PI),
              wrap180(align_data[1] * 180.0 / M_PI),
              wrap180(align_data[2] * 180.0 / M_PI));
    }

  fclose(fp);

  return 0;
}

int
mfield_align_nec2vfm(const gsl_vector * align_params, const double * q, const double B_nec[3], double B_vfm[3], const att_workspace * att_p,
                     const mfield_workspace * w)
{
  gsl_vector_const_view in = gsl_vector_const_view_array(B_nec, 3);
  gsl_vector_view out = gsl_vector_view_array(B_vfm, 3);
  gsl_vector_const_view qv = gsl_vector_const_view_array(q, 4);

  /* apply R_q^T */
  att_rotate_inverse(&qv.vector, &in.vector, &out.vector, w->att_quat_p);

  /* apply R(alpha)^T */
  att_rotate_inverse(align_params, &out.vector, &out.vector, att_p);

  return 0;
}

int
mfield_align_vfm2nec(const gsl_vector * align_params, const double * q, const double B_vfm[3], double B_nec[3], const att_workspace * att_p,
                     const mfield_workspace * w)
{
  gsl_vector_const_view in = gsl_vector_const_view_array(B_vfm, 3);
  gsl_vector_view out = gsl_vector_view_array(B_nec, 3);
  gsl_vector_const_view qv = gsl_vector_const_view_array(q, 4);

  /* apply R(alpha) */
  att_rotate(align_params, &in.vector, &out.vector, att_p);

  /* apply R_q */
  att_rotate(&qv.vector, &out.vector, &out.vector, w->att_quat_p);

  return 0;
}

int
mfield_align_matrix_vfm2nec(const gsl_vector * align_params, const double * q, const gsl_matrix * in, gsl_matrix * out, const att_workspace * att_p,
                            const mfield_workspace * w)
{
  gsl_vector_const_view qv = gsl_vector_const_view_array(q, 4);

  /* apply R(alpha) */
  att_rotate_matrix(align_params, in, out, att_p);

  /* apply R_q */
  att_rotate_matrix(&qv.vector, out, out, w->att_quat_p);

  return 0;
}

int
mfield_align_deriv_vfm2nec(const size_t idx, const gsl_vector * align_params, const double * q, const double B_vfm[3], double B_nec[3],
                           const att_workspace * att_p, const mfield_workspace * w)
{
  gsl_vector_const_view in = gsl_vector_const_view_array(B_vfm, 3);
  gsl_vector_view out = gsl_vector_view_array(B_nec, 3);
  gsl_vector_const_view qv = gsl_vector_const_view_array(q, 4);

  /* apply R(alpha) */
  att_rotate_deriv(idx, align_params, &in.vector, &out.vector, att_p);

  /* apply R_q */
  att_rotate(&qv.vector, &out.vector, &out.vector, w->att_quat_p);

  return 0;
}
