/*
 * stage2b.c
 *
 * 1. read DMSP CDF data file
 * 2. apply previously computed calibration parameters to data
 * 3. write new CDF with calibrated data
 *
 * Usage:
 *
 * ./stage2b <-i dmsp_cdf_file> [-o output_cdf_file]
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <assert.h>
#include <errno.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include <satdata/satdata.h>

#include <common/common.h>
#include <common/bsearch.h>
#include <common/interp.h>
#include <track/track.h>

#include "magcal.h"

#define MAX_PARAM             1000
#define MAX_BUFFER            2048

typedef struct
{
  time_t t[MAX_PARAM];   /* timestamp */
  double SX[MAX_PARAM];  /* scale factor X */
  double SY[MAX_PARAM];  /* scale factor Y */
  double SZ[MAX_PARAM];  /* scale factor Z */
  double OX[MAX_PARAM];  /* offset X (nT) */
  double OY[MAX_PARAM];  /* offset Y (nT) */
  double OZ[MAX_PARAM];  /* offset Z (nT) */
  double AXY[MAX_PARAM]; /* angle AXY (deg) */
  double AXZ[MAX_PARAM]; /* angle AXZ (deg) */
  double AYZ[MAX_PARAM]; /* angle AYZ (deg) */
  size_t n;              /* number of parameters */
} param_data;

int
param_read(const char * filename, param_data * param)
{
  FILE *fp;
  char buf[MAX_BUFFER];
  size_t n = 0;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "param_read: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  while (fgets(buf, MAX_BUFFER, fp) != NULL)
    {
      int c;
      time_t t;
      double rms, SX, SY, SZ, OX, OY, OZ, AXY, AXZ, AYZ;

      if (*buf == '#')
        continue;

      c = sscanf(buf, "%ld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                 &t,
                 &rms,
                 &SX,
                 &SY,
                 &SZ,
                 &OX,
                 &OY,
                 &OZ,
                 &AXY,
                 &AXZ,
                 &AYZ);
      if (c < 11)
        continue;

      param->t[n] = t;
      param->SX[n] = SX;
      param->SY[n] = SY;
      param->SZ[n] = SZ;
      param->OX[n] = OX;
      param->OY[n] = OY;
      param->OZ[n] = OZ;
      param->AXY[n] = AXY;
      param->AXZ[n] = AXZ;
      param->AYZ[n] = AYZ;

      if (++n > MAX_PARAM)
        {
          fprintf(stderr, "param_read: MAX_PARAM too small\n");
          break;
        }
    }

  param->n = n;

  fclose(fp);

  return 0;
}


int
euler_Rq(const double *q, gsl_matrix *Rq)
{
  const double q1 = q[0];
  const double q2 = q[1];
  const double q3 = q[2];
  const double q4 = q[3];

  gsl_matrix_set(Rq, 0, 0, 1.0 - 2.0*q2*q2 - 2.0*q3*q3);
  gsl_matrix_set(Rq, 0, 1, 2.0*(q1*q2 + q3*q4));
  gsl_matrix_set(Rq, 0, 2, 2.0*(q1*q3 - q2*q4));

  gsl_matrix_set(Rq, 1, 0, 2.0*(q1*q2 - q3*q4));
  gsl_matrix_set(Rq, 1, 1, 1.0 - 2.0*q1*q1 - 2.0*q3*q3);
  gsl_matrix_set(Rq, 1, 2, 2.0*(q2*q3 + q1*q4));

  gsl_matrix_set(Rq, 2, 0, 2.0*(q1*q3 + q2*q4));
  gsl_matrix_set(Rq, 2, 1, 2.0*(q2*q3 - q1*q4));
  gsl_matrix_set(Rq, 2, 2, 1.0 - 2.0*q1*q1 - 2.0*q2*q2);

  return GSL_SUCCESS;
}

int
main(int argc, char *argv[])
{
  satdata_mag *data = NULL;
  param_data param;
  struct timeval tv0, tv1;
  int c;
  char *outfile = NULL;
  gsl_vector *coef;
  size_t i;

  param.n = 0;

  while ((c = getopt(argc, argv, "i:o:p:")) != (-1))
    {
      switch (c)
        {
          case 'i':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);

            data = satdata_dmsp_read(optarg, NULL);
            if (!data)
              {
                fprintf(stderr, "main: error reading %s\n", optarg);
                exit(1);
              }

            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu records read, %g seconds)\n", data->n,
                    time_diff(tv0, tv1));
            break;

          case 'o':
            outfile = optarg;
            break;

          case 'p':
            fprintf(stderr, "main: reading parameters from %s...", optarg);
            param_read(optarg, &param);
            fprintf(stderr, "done (%zu parameters read)\n", param.n);
            break;

          default:
            break;
        }
    }

  if (!data || !outfile || param.n == 0)
    {
      fprintf(stderr, "main: usage: %s <-i dmsp_cdf_file> <-p parameter_file> <-o output_cdf_file>\n",
              argv[0]);
      exit(1);
    }

  fprintf(stderr, "main: output file = %s\n", outfile);

  if (data->n == 0)
    {
      fprintf(stderr, "main: no data to process\n");
      exit(1);
    }

  coef = gsl_vector_alloc(MAGCAL_P);

  for (i = 0; i < data->n; ++i)
    {
      time_t unix_time = satdata_epoch2timet(data->t[i]);
      double E[3], B[4];
      double SX, SY, SZ, OX, OY, OZ, AXY, AXZ, AYZ;

      E[0] = SATDATA_VEC_X(data->B_VFM, i);
      E[1] = SATDATA_VEC_Y(data->B_VFM, i);
      E[2] = SATDATA_VEC_Z(data->B_VFM, i);

      if (unix_time < param.t[0])
        {
          SX = param.SX[0];
          SY = param.SY[0];
          SZ = param.SZ[0];
          OX = param.OX[0];
          OY = param.OY[0];
          OZ = param.OZ[0];
          AXY = param.AXY[0];
          AXZ = param.AXZ[0];
          AYZ = param.AYZ[0];
        }
      else if (unix_time > param.t[param.n - 1])
        {
          SX = param.SX[param.n - 1];
          SY = param.SY[param.n - 1];
          SZ = param.SZ[param.n - 1];
          OX = param.OX[param.n - 1];
          OY = param.OY[param.n - 1];
          OZ = param.OZ[param.n - 1];
          AXY = param.AXY[param.n - 1];
          AXZ = param.AXZ[param.n - 1];
          AYZ = param.AYZ[param.n - 1];
        }
      else
        {
          size_t idx = bsearch_timet(param.t, unix_time, 0, param.n - 1);

          /* interpolate calibration parameters to this timestamp */
          SX = interp1d((double) param.t[idx], (double) param.t[idx + 1],
                        param.SX[idx], param.SX[idx + 1], (double) unix_time);
          SY = interp1d((double) param.t[idx], (double) param.t[idx + 1],
                        param.SY[idx], param.SY[idx + 1], (double) unix_time);
          SZ = interp1d((double) param.t[idx], (double) param.t[idx + 1],
                        param.SZ[idx], param.SZ[idx + 1], (double) unix_time);
          OX = interp1d((double) param.t[idx], (double) param.t[idx + 1],
                        param.OX[idx], param.OX[idx + 1], (double) unix_time);
          OY = interp1d((double) param.t[idx], (double) param.t[idx + 1],
                        param.OY[idx], param.OY[idx + 1], (double) unix_time);
          OZ = interp1d((double) param.t[idx], (double) param.t[idx + 1],
                        param.OZ[idx], param.OZ[idx + 1], (double) unix_time);
          AXY = interp1d((double) param.t[idx], (double) param.t[idx + 1],
                         param.AXY[idx], param.AXY[idx + 1], (double) unix_time);
          AXZ = interp1d((double) param.t[idx], (double) param.t[idx + 1],
                         param.AXZ[idx], param.AXZ[idx + 1], (double) unix_time);
          AYZ = interp1d((double) param.t[idx], (double) param.t[idx + 1],
                         param.AYZ[idx], param.AYZ[idx + 1], (double) unix_time);
        }

      /* set initial values of calibration parameters */
      gsl_vector_set(coef, MAGCAL_IDX_SX, SX);
      gsl_vector_set(coef, MAGCAL_IDX_SY, SY);
      gsl_vector_set(coef, MAGCAL_IDX_SZ, SZ);
      gsl_vector_set(coef, MAGCAL_IDX_OX, OX);
      gsl_vector_set(coef, MAGCAL_IDX_OY, OY);
      gsl_vector_set(coef, MAGCAL_IDX_OZ, OZ);
      gsl_vector_set(coef, MAGCAL_IDX_AXY, AXY * M_PI / 180.0);
      gsl_vector_set(coef, MAGCAL_IDX_AXZ, AXZ * M_PI / 180.0);
      gsl_vector_set(coef, MAGCAL_IDX_AYZ, AYZ * M_PI / 180.0);

      /* this function can be called with nT units */
      magcal_apply_cal(coef, E, B);

      /* store new values in data */
      SATDATA_VEC_X(data->B_VFM, i) = B[0];
      SATDATA_VEC_Y(data->B_VFM, i) = B[1];
      SATDATA_VEC_Z(data->B_VFM, i) = B[2];
      data->F[i] = B[3];
    }

  fprintf(stderr, "main: applying scalar calibration to data...");
  gettimeofday(&tv0, NULL);
  magcal_apply(coef, data);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  /* XXX: copy Z component of B_VFM to B_NEC - ignore X and Y */
  {
    size_t i;
    double Rq_data[9];
    gsl_matrix_view Rq = gsl_matrix_view_array(Rq_data, 3, 3);

    for (i = 0; i < data->n; ++i)
      {
        gsl_vector_view B_VFM = gsl_vector_view_array(&(data->B_VFM[3 * i]), 3);
        gsl_vector_view B_NEC = gsl_vector_view_array(&(data->B[3 * i]), 3);
        double *q = &(data->q[4 * i]);

        euler_Rq(q, &Rq.matrix);
        gsl_blas_dgemv(CblasNoTrans, 1.0, &Rq.matrix, &B_VFM.vector, 0.0, &B_NEC.vector);
      }
  }

  fprintf(stderr, "main: writing %s...", outfile);
  gettimeofday(&tv0, NULL);
  satdata_dmsp_write(1, outfile, data);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu records written, %g seconds)\n", data->n,
          time_diff(tv0, tv1));

  gsl_vector_free(coef);
  satdata_mag_free(data);

  return 0;
}
