/*
 * attitude.c
 *
 * Perform a rough attitude correction for the X/Y components of DMSP data
 *
 * For each track, the following model is fitted:
 *
 * delta = c1 * sin(c2 * theta + c3)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlinear.h>

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_quat.h>
#include <mainlib/ml_track.h>

static double attitude_calc_delta(const double A[3], const double B[3]);
static int attitude_model_calc(const track_data * tptr, const satdata_mag * data, double c[3]);
static double attitude_model_eval(const double theta, const double c[3]);
static int attitude_f(const gsl_vector * x, void * params, gsl_vector * f);
static void attitude_callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w);
static int attitude_rotation(const double *q, const double delta, gsl_matrix * R);

typedef struct
{
  size_t start_idx;
  size_t end_idx;
  const satdata_mag * data;
  double qdlat_max;
  double * rhs;
} attitude_params;

/*
attitude_correct()
  Perform rough attitude correction for satellite data

Inputs: filename - data file to write
        data     - satellite data
        track_p  - track workspace

Return: success/error
*/

static int
attitude_correct(const char *filename, const satdata_mag * data, const track_workspace * track_p)
{
  int s = 0;
  FILE *fp = filename ? fopen(filename, "w") : NULL;
  size_t i;

  if (fp)
    {
      i = 1;
      fprintf(fp, "# Field %zu: timestamp (UT seconds since 1970-01-01 00:00:00 UTC)\n", i++);
      fprintf(fp, "# Field %zu: geographic latitude (degrees)\n", i++);
      fprintf(fp, "# Field %zu: QD latitude (degrees)\n", i++);
      fprintf(fp, "# Field %zu: VFM X residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: VFM Y residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: VFM Z residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: attitude correction angle (degrees)\n", i++);
      fprintf(fp, "# Field %zu: attitude correction model (degrees)\n", i++);
      fprintf(fp, "# Field %zu: corrected VFM X residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: corrected VFM Y residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: corrected VFM Z residual (nT)\n", i++);
      fprintf(fp, "# Field %zu: satellite direction (+1 north, -1 south)\n", i++);
    }

  for (i = 0; i < track_p->n; ++i)
    {
      track_data *tptr = &(track_p->tracks[i]);
      size_t start_idx = tptr->start_idx;
      size_t end_idx = tptr->end_idx;
      size_t j;
      double c[3]; /* attitude correction model coefficients */

      attitude_model_calc(tptr, data, c);

      for (j = start_idx; j <= end_idx; ++j)
        {
          double B_model[4], B_model_SC[3], B_model_VFM[3];
          double R_data[9];
          gsl_matrix_view R = gsl_matrix_view_array(R_data, 3, 3);
          double *B_VFM = &(data->B_VFM[3 * j]);
          double *q = &(data->q[4 * j]);
          double theta = M_PI / 2.0 - data->latitude[j] * M_PI / 180.0;
          double delta_model = attitude_model_eval(theta, c);
          double delta;

          /* calculate field model */
          satdata_mag_model(j, B_model, data);

          /* rotate B_model into SC frame */
          quat_apply_inverse(q, B_model, B_model_SC);

          /* compute attitude correction angle delta for this data point */
          delta = attitude_calc_delta(B_VFM, B_model_SC);

          /* compute new rotation matrix (old quaternions plus delta correction) */
          attitude_rotation(q, delta_model, &R.matrix);

          /* compute new quaternions */
          quat_R2q(&R.matrix, q);

#if 0
          /* rotate B_model_SC into VFM frame using correction model */
          B_model_VFM[0] = cos(delta_model) * B_model_SC[0] - sin(delta_model) * B_model_SC[1];
          B_model_VFM[1] = sin(delta_model) * B_model_SC[0] + cos(delta_model) * B_model_SC[1];
          B_model_VFM[2] = B_model_SC[2];
#else
          /* rotate B_model into VFM frame with new quaternions */
          quat_apply_inverse(q, B_model, B_model_VFM);
#endif

          if (fp)
            {
              fprintf(fp, "%ld %8.4f %8.4f %10.4f %10.4f %10.4f %8.4f %8.4f %10.4f %10.4f %10.4f %2d\n",
                      satdata_epoch2timet(data->t[j]),
                      data->latitude[j],
                      data->qdlat[j],
                      B_VFM[0] - B_model_SC[0],
                      B_VFM[1] - B_model_SC[1],
                      B_VFM[2] - B_model_SC[2],
                      delta * 180.0 / M_PI,
                      delta_model * 180.0 / M_PI,
                      B_VFM[0] - B_model_VFM[0],
                      B_VFM[1] - B_model_VFM[1],
                      B_VFM[2] - B_model_VFM[2],
                      satdata_satdir(j, data->n, data->latitude));
            }
        }

      if (fp)
        fprintf(fp, "\n\n");
    }

  if (fp)
    fclose(fp);

  return s;
}

/*
attitude_calc_delta()
  Compute a rotation angle delta such that:

delta = argmin || A - R_z(delta) B ||

and

R_z(delta) = [ cos(delta) -sin(delta) 0 ]
             [ sin(delta)  cos(delta) 0 ]
             [     0           0      1 ]

Inputs: A - vector A
        B - vector B

Return: angle delta in radians
*/

static double
attitude_calc_delta(const double A[3], const double B[3])
{
  const double fac = 1.0 / (B[0]*B[0] + B[1]*B[1]);
  const double cosd = fac * (B[0]*A[0] + B[1]*A[1]);
  const double sind = fac * (-B[1]*A[0] + B[0]*A[1]);
  const double delta = atan2(sind, cosd);

  return delta;
}

/*
attitude_model_calc()
  Calculate a model to correct attitude quaternions for a single
satellite track. The model is:

delta = c1 * sin(c2 * theta + c3)

where theta is co-latitude, and delta is the angle needed to rotate
about the geodetic vertical to align the VFM X and Y with an a-priori
main field model

Inputs: tptr - track data
        data - satellite data
        c    - (output) model coefficients
*/

static int
attitude_model_calc(const track_data * tptr, const satdata_mag * data, double c[3])
{
  int s = 0;
  const double qdlat_max = 40.0; /* only fit data in [-qdlat_max,qdlat_max] */
  const double xtol = 1e-12;
  const double gtol = 1e-12;
  const double ftol = 0.0;
  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  const size_t start_idx = tptr->start_idx;
  const size_t end_idx = tptr->end_idx;
  const size_t p = 3;
  size_t n = 0;
  gsl_multifit_nlinear_workspace *w;
  gsl_multifit_nlinear_fdf fdf;
  gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
  gsl_vector_view cv = gsl_vector_view_array(c, 3);
  double * delta = malloc(tptr->n * sizeof(double));
  attitude_params params;
  size_t i;
  int info;

  /* count number of residuals and compute array of correction angles */
  for (i = start_idx; i <= end_idx; ++i)
    {
      if (fabs(data->qdlat[i]) <= qdlat_max)
        {
          double B_model[4], B_model_SC[3];
          double *B_VFM = &(data->B_VFM[3 * i]);
          double *q = &(data->q[4 * i]);

          /* calculate field model */
          satdata_mag_model(i, B_model, data);

          /* rotate B_model into VFM frame */
          quat_apply_inverse(q, B_model, B_model_SC);

          /* compute attitude correction angle delta for this data point */
          delta[n++] = attitude_calc_delta(B_VFM, B_model_SC);
        }
    }

  params.start_idx = start_idx;
  params.end_idx = end_idx;
  params.data = data;
  params.qdlat_max = qdlat_max;
  params.rhs = delta;

  fdf.f = attitude_f;
  /*fdf.df = attitude_df;*/
  fdf.df = NULL;
  fdf.fvv = NULL;
  fdf.n = n;
  fdf.p = p;
  fdf.params = &params;

  /* allocate multifit workspace */
  fdf_params.trs = gsl_multifit_nlinear_trs_lm;
  fdf_params.solver = gsl_multifit_nlinear_solver_svd;
  fdf_params.scale = gsl_multifit_nlinear_scale_levenberg;
  w = gsl_multifit_nlinear_alloc (T, &fdf_params, n, p);

  /* initialize solver */
  c[0] = 1.0;
  c[1] = 1.0;
  c[2] = 0.0;
  gsl_multifit_nlinear_init(&cv.vector, &fdf, w);

  /* solve system */
  gsl_multifit_nlinear_driver(200, xtol, gtol, ftol, attitude_callback, &params, &info, w);

#if 0
  for (i = start_idx; i <= end_idx; ++i)
    {
      double B_model[4], B_model_SC[3], B_model_VFM[3];
      double *B_VFM = &(data->B_VFM[3 * i]);
      double *q = &(data->q[4 * i]);
      double delta;

      /* calculate field model */
      satdata_mag_model(i, B_model, data);

      /* rotate B_model into VFM frame */
      quat_apply_inverse(q, B_model, B_model_SC);

      /* compute attitude correction angle delta for this data point */
      delta = attitude_calc_delta(B_VFM, B_model_SC);
    }
#endif

  /* save model parameters */
  {
    gsl_vector * x = gsl_multifit_nlinear_position(w);

    for (i = 0; i < p; ++i)
      c[i] = gsl_vector_get(x, i);
  }

  free(delta);
  gsl_multifit_nlinear_free(w);

  return s;
}

static double
attitude_model_eval(const double theta, const double c[3])
{
  double f = c[0] * sin(c[1] * theta + c[2]);
  return f;
}

static int
attitude_f(const gsl_vector * x, void * params, gsl_vector * f)
{
  attitude_params * p = (attitude_params *) params;
  size_t i;
  size_t idx = 0;

  for (i = p->start_idx; i <= p->end_idx; ++i)
    {
      double theta, model;

      if (fabs(p->data->qdlat[i]) > p->qdlat_max)
        continue;

      theta = M_PI / 2.0 - p->data->latitude[i] * M_PI / 180.0;
      model = attitude_model_eval(theta, x->data);

      gsl_vector_set(f, idx, model - p->rhs[idx]);
      ++idx;
    }

  return GSL_SUCCESS;
}

static void
attitude_callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w)
{
  (void) iter;
  (void) params;
  (void) w;
}

/*
attitude_rotation()
  Compute new rotation matrix from VFM to NEC, including
correction model delta:

R_new = R_q R_z(delta)^T

Inputs: q     - old quaternions
        delta - correction angle (radians)
        R     - (output) rotation matrix

Return: success/error
*/

static int
attitude_rotation(const double *q, const double delta, gsl_matrix * R)
{
  const double sd = sin(delta);
  const double cd = cos(delta);
  double Rz_data[9], Rq_data[9];
  gsl_matrix_view Rq = gsl_matrix_view_array(Rq_data, 3, 3);
  gsl_matrix_view Rz = gsl_matrix_view_array(Rz_data, 3, 3);

  /* compute R_q */
  quat_q2R(q, &Rq.matrix);

  /* compute R_z(delta) */
  gsl_matrix_set(&Rz.matrix, 0, 0, cd);
  gsl_matrix_set(&Rz.matrix, 0, 1, -sd);
  gsl_matrix_set(&Rz.matrix, 0, 2, 0.0);

  gsl_matrix_set(&Rz.matrix, 1, 0, sd);
  gsl_matrix_set(&Rz.matrix, 1, 1, cd);
  gsl_matrix_set(&Rz.matrix, 1, 2, 0.0);

  gsl_matrix_set(&Rz.matrix, 2, 0, 0.0);
  gsl_matrix_set(&Rz.matrix, 2, 1, 0.0);
  gsl_matrix_set(&Rz.matrix, 2, 2, 1.0);

  /* R = R_q R_z(delta)^T */
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &Rq.matrix, &Rz.matrix, 0.0, R);

  return GSL_SUCCESS;
}
