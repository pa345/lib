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
#include <string.h>
#include <errno.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>

#include <mainlib/ml_common.h>
#include <bspline2/gsl_bspline2.h>

#include "mfield_fluxcal.h"

static int fluxcal_P(const double u[3], gsl_matrix * P);
static int fluxcal_Pinv(const double u[3], gsl_matrix * Pinv);
static int fluxcal_Pinv_deriv(const double u[3], gsl_matrix * Pinv_1, gsl_matrix * Pinv_2, gsl_matrix * Pinv_3);
static int fluxcal_apply_datum_Pinv(const gsl_matrix * Pinv, const gsl_vector *m, const double E[3], double B[4]);

int
mfield_fluxcal_print(const char *filename, const size_t sat_idx,
                     const mfield_workspace *w)
{
  FILE *fp;
  gsl_bspline2_workspace *spline_p = w->fluxcal_spline_workspace_p[CIDX2(sat_idx, w->nsat, 0, w->max_threads)];
  const size_t ncontrol = gsl_bspline2_ncontrol(spline_p);
  const size_t fluxcal_idx = w->fluxcal_offset + w->offset_fluxcal[sat_idx];

  double cal_data[FLUXCAL_P];
  gsl_vector_view cal_params = gsl_vector_view_array(cal_data, FLUXCAL_P);

  /* fluxgate calibration control points for this dataset */
  gsl_vector_const_view tmp = gsl_vector_const_subvector(w->c, fluxcal_idx, FLUXCAL_P * ncontrol);
  gsl_matrix_const_view control_pts = gsl_matrix_const_view_vector(&tmp.vector, FLUXCAL_P, ncontrol);

  const double dt = 5.0 * 8.64e7; /* 5 days in ms */
  const double t0 = w->data_workspace_p->t0[sat_idx];
  const double t1 = w->data_workspace_p->t1[sat_idx];
  double t;
  size_t i;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "mfield_fluxcal_print: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Fluxgate calibration parameters for satellite %zu\n", sat_idx);
  fprintf(fp, "# Field %zu: timestamp (CDF_EPOCH)\n", i++);
  fprintf(fp, "# Field %zu: time (decimal years)\n", i++);
  fprintf(fp, "# Field %zu: SX (dimensionless)\n", i++);
  fprintf(fp, "# Field %zu: SY (dimensionless)\n", i++);
  fprintf(fp, "# Field %zu: SZ (dimensionless)\n", i++);
  fprintf(fp, "# Field %zu: OX (nT)\n", i++);
  fprintf(fp, "# Field %zu: OY (nT)\n", i++);
  fprintf(fp, "# Field %zu: OZ (nT)\n", i++);
  fprintf(fp, "# Field %zu: U1 (degrees)\n", i++);
  fprintf(fp, "# Field %zu: U2 (degrees)\n", i++);
  fprintf(fp, "# Field %zu: U3 (degrees)\n", i++);

  for (t = t0; t < t1; t += dt)
    {
      double t_year = epoch2year(t);

      gsl_bspline2_vector_eval(t_year, &control_pts.matrix, &cal_params.vector, spline_p);

      fprintf(fp, "%f %f %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\n",
              t,
              t_year,
              cal_data[FLUXCAL_IDX_SX],
              cal_data[FLUXCAL_IDX_SY],
              cal_data[FLUXCAL_IDX_SZ],
              cal_data[FLUXCAL_IDX_OX],
              cal_data[FLUXCAL_IDX_OY],
              cal_data[FLUXCAL_IDX_OZ],
              cal_data[FLUXCAL_IDX_U1] * 180.0 / M_PI,
              cal_data[FLUXCAL_IDX_U2] * 180.0 / M_PI,
              cal_data[FLUXCAL_IDX_U3] * 180.0 / M_PI);
    }

  fclose(fp);

  return 0;
}

/*
fluxcal_apply_datum()
  Apply calibration to a single vector measurement

B = P^{-1} S (E - O)

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
2) It is allowed for E == B for in-place transform
*/

int
mfield_fluxcal_apply_datum(const gsl_vector *m, const double E[3], double B[4])
{
  int s;
  double U[3], Pinv_data[9];
  gsl_matrix_view Pinv = gsl_matrix_view_array(Pinv_data, 3, 3);

  U[0] = gsl_vector_get(m, FLUXCAL_IDX_U1);
  U[1] = gsl_vector_get(m, FLUXCAL_IDX_U2);
  U[2] = gsl_vector_get(m, FLUXCAL_IDX_U3);

  /* construct P^{-1} */
  s = fluxcal_Pinv(U, &Pinv.matrix);
  if (s)
    return s;

  /* compute B = P^{-1} S (E - O) */
  s = fluxcal_apply_datum_Pinv(&Pinv.matrix, m, E, B);

  return s;
}

/*
fluxcal_invapply_datum()
  Apply inverse calibration to a single vector measurement

E = S^{-1} P B + O

Inputs: m - model parameters
        B - calibrated vector measurements
            B[0] = B_x_calibrated
            B[1] = B_y_calibrated
            B[2] = B_z_calibrated
            B[3] = F_calibrated
        E - (output) uncalibrated vector measurements
            E[0] = B_x_uncalibrated
            E[1] = B_y_uncalibrated
            E[2] = B_z_uncalibrated

Return: success or error

Notes:
1) See Eq. 1 of Olsen et al, 2003
2) It is allowed for E == B for in-place transform
*/

int
mfield_fluxcal_invapply_datum(const gsl_vector *m, const double B[3], double E[3])
{
  int s = 0;
  double S[3], O[3], U[3], P_data[9];
  gsl_matrix_view P = gsl_matrix_view_array(P_data, 3, 3);
  gsl_vector_view Ev = gsl_vector_view_array(E, 3);
  size_t i;

  S[0] = gsl_vector_get(m, FLUXCAL_IDX_SX);
  S[1] = gsl_vector_get(m, FLUXCAL_IDX_SY);
  S[2] = gsl_vector_get(m, FLUXCAL_IDX_SZ);

  O[0] = gsl_vector_get(m, FLUXCAL_IDX_OX);
  O[1] = gsl_vector_get(m, FLUXCAL_IDX_OY);
  O[2] = gsl_vector_get(m, FLUXCAL_IDX_OZ);

  U[0] = gsl_vector_get(m, FLUXCAL_IDX_U1);
  U[1] = gsl_vector_get(m, FLUXCAL_IDX_U2);
  U[2] = gsl_vector_get(m, FLUXCAL_IDX_U3);

  /* construct P */
  fluxcal_P(U, &P.matrix);

  for (i = 0; i < 3; ++i)
    E[i] = B[i];

  /* compute E = P B */
  gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, &P.matrix, &Ev.vector);

  /* compute E = S^{-1} P B + O */
  for (i = 0; i < 3; ++i)
    E[i] = E[i] / S[i] + O[i];

  return s;
}

/*
fluxcal_P()
  Construct inverse matrix of non-orthogonality transformation:

P = [    1         0        0  ]
    [ -sin(u1)   cos(u1)    0  ]
    [  sin(u2)   sin(u3)    w  ]

where:

w = sqrt(1 - sin(u2)^2 - sin(u3)^2)

Notes:
1) Only the lower triangle of P is filled in
*/

static int
fluxcal_P(const double u[3], gsl_matrix * P)
{
  const double s1 = sin(u[0]);
  const double c1 = cos(u[0]);
  const double s2 = sin(u[1]);
  const double s3 = sin(u[2]);
  const double w = sqrt(1.0 - s2*s2 - s3*s3);

  gsl_matrix_set(P, 0, 0, 1.0);

  gsl_matrix_set(P, 1, 0, -s1);
  gsl_matrix_set(P, 1, 1,  c1);

  gsl_matrix_set(P, 2, 0, s2);
  gsl_matrix_set(P, 2, 1, s3);
  gsl_matrix_set(P, 2, 2, w);

  return GSL_SUCCESS;
}

/*
fluxcal_Pinv()
  Construct inverse matrix of non-orthogonality transformation:

P^{-1} = [                   1                                              0                0  ]
         [                 tan(u1)                                        sec(u1)            0  ]
         [ -1/(w*cos(u1)) * (sin(u1)*sin(u3) + cos(u1)*sin(u2))     -sin(u3)/(w*cos(u1))    1/w ]

where:

w = sqrt(1 - sin(u2)^2 - sin(u3)^2)

Return: 0 on success, -1 if angles u cause a singularity in the P matrix
*/

static int
fluxcal_Pinv(const double u[3], gsl_matrix * Pinv)
{
  const double s1 = sin(u[0]);
  const double c1 = cos(u[0]);
  const double s2 = sin(u[1]);
  const double s3 = sin(u[2]);
  const double sterm = s2*s2 + s3*s3;
  double w;

  if (sterm > 1.0)
    return -1; /* invalid (u2,u3) parameters */

  w = sqrt(1.0 - sterm);

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

Notes:
1) Only lower triangular portion of matrices are filled
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

  gsl_matrix_set(Pinv_1, 1, 0, 1.0 / (c1 * c1));
  gsl_matrix_set(Pinv_1, 1, 1, s1 / (c1 * c1));

  gsl_matrix_set(Pinv_1, 2, 0, -s3 / (w * c1 * c1));
  gsl_matrix_set(Pinv_1, 2, 1, -s3 * s1 / (w * c1 * c1));
  gsl_matrix_set(Pinv_1, 2, 2, 0.0);

  /* d/du_2 P^{-1} */

  gsl_matrix_set(Pinv_2, 0, 0, 0.0);

  gsl_matrix_set(Pinv_2, 1, 0, 0.0);
  gsl_matrix_set(Pinv_2, 1, 1, 0.0);

  gsl_matrix_set(Pinv_2, 2, 0, -c2 / (w * w * w) * (c3 * c3 + s2 * s3 * (s1 / c1)));
  gsl_matrix_set(Pinv_2, 2, 1, -s3 * s2 * c2 / (c1 * w * w * w));
  gsl_matrix_set(Pinv_2, 2, 2, s2 * c2 / (w * w * w));

  /* d/du_3 P^{-1} */

  gsl_matrix_set(Pinv_3, 0, 0, 0.0);

  gsl_matrix_set(Pinv_3, 1, 0, 0.0);
  gsl_matrix_set(Pinv_3, 1, 1, 0.0);

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
2) It is allowed for E == B for in-place transform
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

/*
fluxcal_jac()
  Compute elements of the Jacobian matrix corresponding to B(m)

with B(m) = P^{-1} S (E - O)

Inputs: m   - calibration parameters
        E   - uncalibrated vector
        jac - (output) d/dm B(m), 3-by-FLUXCAL_P
              jac(1,j) = d/dm_j B_x(m)
              jac(2,j) = d/dm_j B_y(m)
              jac(3,j) = d/dm_j B_z(m)
              j = 1,2,...,FLUXCAL_P
*/

int
mfield_fluxcal_jac(const gsl_vector *m, const double E[3], gsl_matrix *jac)
{
  double S[3], O[3], U[3];
  double Pinv_data[9], Pinv_1_data[9], Pinv_2_data[9], Pinv_3_data[9];
  gsl_matrix_view Pinv = gsl_matrix_view_array(Pinv_data, 3, 3);
  gsl_matrix_view Pinv_1 = gsl_matrix_view_array(Pinv_1_data, 3, 3);
  gsl_matrix_view Pinv_2 = gsl_matrix_view_array(Pinv_2_data, 3, 3);
  gsl_matrix_view Pinv_3 = gsl_matrix_view_array(Pinv_3_data, 3, 3);
  double tmp1[4], tmp2[4], tmp3[4];
  size_t j;

  S[0] = gsl_vector_get(m, FLUXCAL_IDX_SX);
  S[1] = gsl_vector_get(m, FLUXCAL_IDX_SY);
  S[2] = gsl_vector_get(m, FLUXCAL_IDX_SZ);

  O[0] = gsl_vector_get(m, FLUXCAL_IDX_OX);
  O[1] = gsl_vector_get(m, FLUXCAL_IDX_OY);
  O[2] = gsl_vector_get(m, FLUXCAL_IDX_OZ);

  U[0] = gsl_vector_get(m, FLUXCAL_IDX_U1);
  U[1] = gsl_vector_get(m, FLUXCAL_IDX_U2);
  U[2] = gsl_vector_get(m, FLUXCAL_IDX_U3);

  /* construct P^{-1} */
  fluxcal_Pinv(U, &Pinv.matrix);

  for (j = 0; j < 3; ++j)
    {
      size_t k;

      for (k = 0; k < 3; ++k)
        {
          double Pinv_kj = gsl_matrix_get(&Pinv.matrix, k, j);

          /* compute jac_S = d/dS eps(m) = (E_j - O_j) R_q R_3 P^{-1}_j */
          gsl_matrix_set(jac, k, FLUXCAL_IDX_SX + j, (E[j] - O[j]) * Pinv_kj);

          /* compute jac_O = d/dO eps(m) = -S_j R_q R_3 P^{-1}_j */
          gsl_matrix_set(jac, k, FLUXCAL_IDX_OX + j, -S[j] * Pinv_kj);
        }
    }

  /* compute d/du_j P^{-1} matrices */
  fluxcal_Pinv_deriv(U, &Pinv_1.matrix, &Pinv_2.matrix, &Pinv_3.matrix);

  /* compute tmp_j = d/dU_j P^{-1} * S * (E - O) */
  fluxcal_apply_datum_Pinv(&Pinv_1.matrix, m, E, tmp1);
  fluxcal_apply_datum_Pinv(&Pinv_2.matrix, m, E, tmp2);
  fluxcal_apply_datum_Pinv(&Pinv_3.matrix, m, E, tmp3);

  /* store in output matrices */
  for (j = 0; j < 3; ++j)
    {
      gsl_matrix_set(jac, j, FLUXCAL_IDX_U1, tmp1[j]);
      gsl_matrix_set(jac, j, FLUXCAL_IDX_U2, tmp2[j]);
      gsl_matrix_set(jac, j, FLUXCAL_IDX_U3, tmp3[j]);
    }

  return GSL_SUCCESS;
}
