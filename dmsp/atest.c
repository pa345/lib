/*
 * atest.c
 *
 * Compute Euler angles for DMSP dataset
 *
 * 1. Read DMSP file(s)
 * 2. Compute residuals in spacecraft-fixed frame:
 *
 * eps_i = R_3(alpha,beta,gamma) B_i^{VFM} - R_q^T B^{model}(t_i,r_i)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_filter.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_rstat.h>
#include <gsl/gsl_statistics.h>

#include <satdata/satdata.h>
#include <indices/indices.h>
#include <common/common.h>
#include <common/quat.h>
#include <msynth/msynth.h>
#include <track/track.h>

#include "euler.h"
#include "euler_calc.h"
#include "eph.h"

int
proc(satdata_mag * data, double euler[3])
{
  int s = 0;
  euler_calc_workspace * w = euler_calc_alloc(data->n);
  size_t i;

  for (i = 0; i < data->n; ++i)
    {
      double *q = &(data->q[4 * i]);
      double *B_VFM = &(data->B_VFM[3 * i]);
      double B_model[4];

      satdata_mag_model(i, B_model, data);

      if (fabs(data->qdlat[i]) > 50.0)
        continue;

      euler_calc_add(data->t[i], data->qdlat[i], B_VFM, B_model, q, w);
    }

  euler_calc_proc(w);

  for (i = 0; i < 3; ++i)
    euler[i] = gsl_vector_get(w->c, i);

  euler_calc_free(w);

  return s;
}

/* print model residuals in spacecraft frame (s1,s2,s3) */
int
print_data(satdata_mag * data, double euler[3])
{
  int s = 0;
  size_t i;

  i = 1;
  printf("# Field %zu: timestamp\n", i++);
  printf("# Field %zu: geocentric latitude (degrees)\n", i++);
  printf("# Field %zu: QD latitude (degrees)\n", i++);
  printf("# Field %zu: X in spacecraft frame\n", i++);
  printf("# Field %zu: Y in spacecraft frame\n", i++);
  printf("# Field %zu: Z in spacecraft frame\n", i++);
  printf("# Field %zu: model X in spacecraft frame\n", i++);
  printf("# Field %zu: model Y in spacecraft frame\n", i++);
  printf("# Field %zu: model Z in spacecraft frame\n", i++);

  for (i = 0; i < data->n; ++i)
    {
      time_t unix_time = satdata_epoch2timet(data->t[i]);
      double *q = &(data->q[4 * i]);
      double *B_VFM = &(data->B_VFM[3 * i]);
      double B_model[4], B_model_s[3], B_VFM_s[3];

      satdata_mag_model(i, B_model, data);

      /* rotate B_model into s frame */
      quat_apply_inverse(q, B_model, B_model_s);

      /* apply Euler angle rotation to B_VFM */
      euler_apply_R3(EULER_FLG_ZYX, euler[0], euler[1], euler[2], B_VFM, B_VFM_s);

      printf("%ld %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",
             unix_time,
             data->latitude[i],
             data->qdlat[i],
             B_VFM_s[0],
             B_VFM_s[1],
             B_VFM_s[2],
             B_model_s[0],
             B_model_s[1],
             B_model_s[2]);
    }

  return s;
}

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s <-i dmsp_index_file>\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -i dmsp_index_file     - Input DMSP index file\n");
}

int
main(int argc, char *argv[])
{
  satdata_mag *data = NULL;
  double euler[3];
  struct timeval tv0, tv1;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "i:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'i':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);
            data = satdata_dmsp_read_idx(optarg, 1);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu records read, %g seconds)\n", data->n,
                    time_diff(tv0, tv1));
            break;

          default:
            break;
        }
    }

  if (!data)
    {
      print_help(argv);
      exit(1);
    }

  proc(data, euler);

  print_data(data, euler);

  return 0;
}
