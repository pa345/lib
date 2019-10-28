/*
 * calc_euler.c
 *
 * Compute Euler angles for CHAMP dataset
 *
 * 1. Read CHAMP file(s)
 * 2. Compute residuals in NEC frame:
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
#include <gsl/gsl_multifit.h>

#include <mainlib/ml_att_calc.h>
#include <mainlib/ml_satdata.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_quat.h>

int
proc_att(satdata_mag * data, att_calc_workspace * w, gsl_vector * attitude)
{
  int s = 0;
  att_workspace * att_quat_p = att_alloc(att_quaternion);
  size_t i;

  for (i = 0; i < data->n; ++i)
    {
      gsl_vector_view q = gsl_vector_view_array(&(data->q[4 * i]), 4);
      double * B_VFM = &(data->B_VFM[3 * i]);
      double B_model[4];
      gsl_vector_view v = gsl_vector_view_array(B_model, 3);

      if (!SATDATA_ExistVector(data->flags[i]))
        continue;

      if (fabs(data->qdlat[i]) > 50.0)
        continue;

      satdata_mag_model(i, B_model, data);

      /* rotate B_model to star camera (S) frame */
      att_rotate_inverse(&q.vector, &v.vector, &v.vector, att_quat_p);

      att_calc_add(B_VFM, B_model, w);
    }

  att_calc_proc(w);

  gsl_vector_memcpy(attitude, w->c);

  att_free(att_quat_p);

  return s;
}

/* print model residuals in NEC frame */
int
print_att_data(satdata_mag * data, const gsl_vector * attitude, att_calc_workspace * w)
{
  int s = 0;
  size_t i;

  i = 1;
  printf("# Field %zu: timestamp\n", i++);
  printf("# Field %zu: geocentric latitude (degrees)\n", i++);
  printf("# Field %zu: QD latitude (degrees)\n", i++);
  printf("# Field %zu: X rotated to NEC frame\n", i++);
  printf("# Field %zu: Y rotated to NEC frame\n", i++);
  printf("# Field %zu: Z rotated to NEC frame\n", i++);
  printf("# Field %zu: model X in NEC frame\n", i++);
  printf("# Field %zu: model Y in NEC frame\n", i++);
  printf("# Field %zu: model Z in NEC frame\n", i++);
  printf("# Field %zu: X in NEC frame\n", i++);
  printf("# Field %zu: Y in NEC frame\n", i++);
  printf("# Field %zu: Z in NEC frame\n", i++);

  for (i = 0; i < data->n; ++i)
    {
      time_t unix_time = satdata_epoch2timet(data->t[i]);
      double *q = &(data->q[4 * i]);
      gsl_vector_view B_VFM = gsl_vector_view_array(&(data->B_VFM[3 * i]), 3);
      double B_NEC_data[3];
      gsl_vector_view B_NEC = gsl_vector_view_array(B_NEC_data, 3);
      double * B_NEC2 = &(data->B[3 * i]);
      double B_model[4];

      if (!SATDATA_ExistVector(data->flags[i]))
        continue;

      satdata_mag_model(i, B_model, data);

      /* rotate B_VFM into spacecraft frame */
      att_rotate(attitude, &B_VFM.vector, &B_NEC.vector, w->att_workspace_p);

      /* rotate B_VFM into NEC frame */
      quat_apply(q, B_NEC_data, B_NEC_data);

      printf("%ld %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",
             unix_time,
             data->latitude[i],
             data->qdlat[i],
             B_NEC_data[0],
             B_NEC_data[1],
             B_NEC_data[2],
             B_model[0],
             B_model[1],
             B_model[2],
             B_NEC2[0],
             B_NEC2[1],
             B_NEC2[2]);
    }

  return s;
}

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s <-i champp_index_file>\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -i champp_index_file     - Input CHAMP index file\n");
}

int
main(int argc, char *argv[])
{
  satdata_mag *data = NULL;
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
            data = satdata_champ_read_idx(optarg, 1);
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

  satdata_champ_filter_instrument(1, 0, data);
  satdata_champ_fix_q(0, data);

  {
#if 1
    att_calc_workspace * att_calc_p = att_calc_alloc(data->n, att_euler_zyx);
#else
    att_calc_workspace * att_calc_p = att_calc_alloc(data->n, att_mrp);
#endif
    gsl_vector * attitude = gsl_vector_alloc(att_calc_p->p);

    proc_att(data, att_calc_p, attitude);
    print_att_data(data, attitude, att_calc_p);

    att_calc_free(att_calc_p);
    gsl_vector_free(attitude);
  }

  return 0;
}
