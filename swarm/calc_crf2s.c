/*
 * calc_crf2s.c
 *
 * Compute Euler angles to go from CRF (star camera) frame to
 * s1,s2,s3 basis
 *
 * Output is a data file 'crf2s.dat' containing the time series of computed
 * Euler angles
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <sys/time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_eci.h>
#include <mainlib/ml_euler.h>
#include <mainlib/ml_oct.h>
#include <mainlib/ml_track.h>

#include "eulersat.c"
#include "filter.c"

/*
flag_data_times()
  Remove all current data flags with SATDATA_FLG_TIME, and
then flag all data points *not* inside [t0,t1] with
SATDATA_FLG_TIME. Used for processing 1 window of data at
a time.
*/

int
flag_data_times(const double t0, const double t1,
                satdata_mag *data)
{
  size_t i;

  for (i = 0; i < data->n; ++i)
    {
      if (data->flags[i] & SATDATA_FLG_TIME)
        data->flags[i] &= ~SATDATA_FLG_TIME;

      if (data->t[i] < t0 || data->t[i] > t1)
        data->flags[i] |= SATDATA_FLG_TIME;
    }

  return 0;
}

int
print_euler_angles(const char *filename, gsl_matrix *A)
{
  FILE *fp;
  const size_t n = A->size1;
  size_t i;
  size_t ndata = 0;
  double alpha_mean, beta_mean, gamma_mean;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_euler_angles: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  /* compute means of angles */
  alpha_mean = beta_mean = gamma_mean = 0.0;
  for (i = 0; i < n; ++i)
    {
      double alpha = gsl_matrix_get(A, i, 1);
      double beta = gsl_matrix_get(A, i, 2);
      double gamma = gsl_matrix_get(A, i, 3);

      /* check for missing data */
      if (alpha == 0.0 || beta == 0.0 || gamma == 0.0)
        continue;

      alpha_mean += alpha;
      beta_mean += beta;
      gamma_mean += gamma;
      ++ndata;
    }

  alpha_mean /= (double) ndata;
  beta_mean /= (double) ndata;
  gamma_mean /= (double) ndata;

  i = 1;
  fprintf(fp, "# Field %zu: time (years)\n", i++);
  fprintf(fp, "# Field %zu: alpha (deg)\n", i++);
  fprintf(fp, "# Field %zu: beta (deg)\n", i++);
  fprintf(fp, "# Field %zu: gamma (deg)\n", i++);
  fprintf(fp, "# Field %zu: alpha - mean (arcseconds)\n", i++);
  fprintf(fp, "# Field %zu: beta - mean (arcseconds)\n", i++);
  fprintf(fp, "# Field %zu: gamma - mean (arcseconds)\n", i++);

  for (i = 0; i < n; ++i)
    {
      double alpha = gsl_matrix_get(A, i, 1);
      double beta = gsl_matrix_get(A, i, 2);
      double gamma = gsl_matrix_get(A, i, 3);

      /* check for missing data */
      if (alpha == 0.0 || beta == 0.0 || gamma == 0.0)
        continue;

      /* matrix A is already in degrees */
      fprintf(fp, "%f %f %f %f %f %f %f\n",
              gsl_matrix_get(A, i, 0),
              alpha,
              beta,
              gamma,
              (alpha - alpha_mean) * 3600.0,
              (beta - beta_mean) * 3600.0,
              (gamma - gamma_mean) * 3600.0);
    }

  fclose(fp);

  return 0;
}

static int
crf2s_f(const gsl_vector *m, void *params, gsl_vector *f)
{
  int s = GSL_SUCCESS;
  eulersat_params *p = ((eulersat_params *) params);
  const satdata_mag *data = p->data;
  size_t i, j;
  size_t idx = 0;
  double Rs_data[9], Rq_data[9], R3_data[9], A_data[9];
  gsl_matrix_view Rs = gsl_matrix_view_array(Rs_data, 3, 3);
  gsl_matrix_view Rq = gsl_matrix_view_array(Rq_data, 3, 3);
  gsl_matrix_view R3 = gsl_matrix_view_array(R3_data, 3, 3);
  gsl_matrix_view A = gsl_matrix_view_array(A_data, 3, 3);
  const double alpha = gsl_vector_get(m, 0);
  const double beta = gsl_vector_get(m, 1);
  const double gamma = gsl_vector_get(m, 2);

  gsl_vector_set_zero(f);

  for (i = 0; i < data->n; ++i)
    {
      double *q = &(data->q[4*i]);

      /* skip flagged data */
      if (data->flags[i])
        continue;

      /* compute R_q and R_s for this position */
      euler_Rq(q, &Rq.matrix);
      eulersat_Rs(i, &Rs.matrix, data);

      /* compute Euler rotation matrix R_3 */
      euler_R3(alpha, beta, gamma, &R3.matrix);

      /* compute A = R_s^T R_q */
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &Rs.matrix,
                     &Rq.matrix, 0.0, &A.matrix);

      for (j = 0; j < 9; ++j)
        gsl_vector_set(f, idx++, R3.matrix.data[j] - A.matrix.data[j]);
    }

  assert(idx == p->n);

  return s;
} /* crf2s_f() */

/*
crf2s_proc()
  Compute Euler angles alpha,beta,gamma such that

R_s R_3(alpha,beta,gamma) B_VFM - B_main

is minimized in a least-squares sense.

Inputs: flags - EULER_FLG_xxx
        m     - (output) vector of 3 Euler angles
        data  - satellite data

Return: success/error
*/

static int
crf2s_proc(const size_t flags, gsl_vector *m, const satdata_mag *data)
{
  int s = 0;
  size_t n = data->n - satdata_nflagged(data);
  int info;
  gsl_multifit_fdfsolver *fdf_s;
  gsl_multifit_function_fdf f;
  const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
  eulersat_params params;

  if (n < 500)
    {
      fprintf(stderr, "crf2s_proc: insufficient number of data to process: %zu\n",
              n);
      return -1;
    }

  /* initial values */
  gsl_vector_set(m, EULER_IDX_ALPHA, 0.0 * M_PI / 180.0);
  gsl_vector_set(m, EULER_IDX_BETA, 0.0 * M_PI / 180.0);
  gsl_vector_set(m, EULER_IDX_GAMMA, 0.0 * M_PI / 180.0);

  f.f = &crf2s_f;
#if 1
  f.df = NULL;
#else
  f.df = &crf2s_df;
#endif
  f.n = 9 * n; /* 9 matrix elements */
  f.p = 3;
  f.params = &params;

  params.n = f.n;
  params.data = data;
  params.flags = flags;

  fprintf(stderr, "crf2s_proc: number of data = %zu\n", f.n);

  fdf_s = gsl_multifit_fdfsolver_alloc(T, f.n, f.p);

  fprintf(stderr, "crf2s_proc: initializing fdfsolver...");
  gsl_multifit_fdfsolver_set (fdf_s, &f, m);
  fprintf(stderr, "done\n");

  fprintf(stderr, "crf2s_proc: computing Euler angles...");
  s = eulersat_nonlinear_driver (fdf_s, 100, 1.0e-8, 1.0e-8, 0.0, &info);
  if (s != GSL_SUCCESS)
    fprintf(stderr, "crf2s_proc: error in fdfsolver: %d\n", s);
  fprintf(stderr, "done\n");

  eulersat_print_state(gsl_multifit_fdfsolver_niter(fdf_s), fdf_s);

  /* save angles */
  gsl_vector_memcpy(m, fdf_s->x);

  gsl_multifit_fdfsolver_free(fdf_s);

  return s;
} /* crf2s_proc() */

int
main_proc(const double period, satdata_mag *data)
{
  int s = 0;
  size_t i;
  gsl_vector *euler_m; /* Euler angles for VFM->SAT transformation */
  char *outfile = "crf2s.dat";
  size_t flags = 0;
  double t0, t1, t, dt; /* CDF_EPOCH */
  size_t nt;            /* number of time steps */

#if 0
  {
    double Rs_data[9], Rq_data[9], R_data[9];
    gsl_matrix_view Rs = gsl_matrix_view_array(Rs_data, 3, 3);
    gsl_matrix_view Rq = gsl_matrix_view_array(Rq_data, 3, 3);
    gsl_matrix_view R = gsl_matrix_view_array(R_data, 3, 3);

    for (i = 0; i < data->n; ++i)
      {
        double *q = &(data->q[4*i]);

        if (data->flags[i])
          continue;
        if (satdata_epoch2year(data->t[i]) < 2014.4)
          continue;

        euler_Rq(q, &Rq.matrix);
        eulersat_Rs(i, &Rs.matrix, data);

        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &Rs.matrix,
                       &Rq.matrix, 0.0, &R.matrix);
        print_octave(&R.matrix, NULL);
        getchar();
      }

    exit(1);
  }
#endif

  /*
   * nt-by-4 matrix, columns are:
   * time alpha beta gamma
   */
  gsl_matrix *A;

  euler_m = gsl_vector_alloc(3);

  /* search for first unflagged data time */
  for (i = 0; i < data->n; ++i)
    {
      if (data->flags[i])
        continue;

      t0 = data->t[i];
      break;
    }

  /* search for last unflagged data time */
  for (i = data->n; i > 0 && i--; )
    {
      if (data->flags[i])
        continue;

      t1 = data->t[i];
      break;
    }

  /* compute time step in ms */
  dt = period * 86400000;
  nt = (t1 - t0) / dt + 1;

  A = gsl_matrix_calloc(nt, 4);

  i = 0;
  for (t = t0; t <= t1; t += dt, i++)
    {
      flag_data_times(t, t + dt, data);

      /* compute Euler angles for times in [t,t+dt] */
      s = crf2s_proc(flags, euler_m, data);
      if (s)
        continue; /* insufficient number of data */

      /* save Euler angles in matrix A */
      gsl_matrix_set(A, i, 0, satdata_epoch2year(t + 0.5*dt));
      gsl_matrix_set(A, i, 1, wrap180(gsl_vector_get(euler_m, 0) * 180.0 / M_PI));
      gsl_matrix_set(A, i, 2, wrap180(gsl_vector_get(euler_m, 1) * 180.0 / M_PI));
      gsl_matrix_set(A, i, 3, wrap180(gsl_vector_get(euler_m, 2) * 180.0 / M_PI));
    }

  fprintf(stderr, "main_proc: writing euler angles to %s...", outfile);
  print_euler_angles(outfile, A);
  fprintf(stderr, "done\n");

  gsl_vector_free(euler_m);
  gsl_matrix_free(A);

  return s;
}

int
main(int argc, char *argv[])
{
  int c;
  satdata_mag *data = NULL;
  struct timeval tv0, tv1;
  double time_period = 3.0; /* time window in days */

  while ((c = getopt(argc, argv, "i:p:")) != (-1))
    {
      switch (c)
        {
          case 'i':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);
            data = satdata_swarm_read_idx(optarg, 0);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu points read, %g seconds)\n",
                    data->n, time_diff(tv0, tv1));
            break;

          case 'p':
            time_period = atof(optarg);
            break;

          default:
            fprintf(stderr, "Usage: %s <-i swarm_idx_file>\n", argv[0]);
            exit(1);
            break;
        }
    }

  if (!data)
    {
      fprintf(stderr, "Usage: %s <-i swarm_idx_file>\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "main: time window = %g [days]\n", time_period);

  /* filter by track rms */
  {
    size_t nrms;
    double thresh[] = { 30.0, 25.0, 17.0, 30.0 };
    track_workspace *track_p = track_alloc();

    track_init(data, NULL, track_p);

    nrms = track_filter_rms("swarm_rms.dat", thresh, data, track_p);
    fprintf(stderr, "main: flagged (%zu/%zu) (%.1f%%) points due to high rms\n",
            nrms, data->n, (double) nrms / (double) data->n * 100.0);
    track_free(track_p);
  }

  {
    const double lt_min = 6.0;
    const double lt_max = 18.0;
    size_t nlt;

    fprintf(stderr, "main: discarding local times in [%g,%g]...",
            lt_min, lt_max);
    gettimeofday(&tv0, NULL);
    nlt = swarm_filter_lt(lt_min, lt_max, data);
    gettimeofday(&tv1, NULL);
    fprintf(stderr, "done (%zu/%zu data flagged, %g seconds)\n",
            nlt, data->n, time_diff(tv0, tv1));
  }

  {
    const double qdlat_min = 0.0;
    const double qdlat_max = 60.0;
    size_t nqd;

    fprintf(stderr, "main: flagging data outside qdlat range [%g,%g]...",
            qdlat_min, qdlat_max);
    nqd = swarm_filter_qdlat(qdlat_min, qdlat_max, data);
    fprintf(stderr, "done (%zu/%zu data flagged)\n", nqd, data->n);
  }

  {
    size_t i;
    const int downsample = 20;

    fprintf(stderr, "main: downsampling data by factor %d...", downsample);

    for (i = 0; i < data->n; ++i)
      {
        if (i % downsample != 0)
          data->flags[i] |= SATDATA_FLG_OUTLIER;
      }

    fprintf(stderr, "done\n");
  }

#if 0
  /* discard data before 2014.4 (maneuvers etc), only for alt basis */
  if (use_altbasis)
    {
      size_t i;
      const double thresh = 2014.4;

      fprintf(stderr, "main: discarding data prior to %g...", thresh);

      for (i = 0; i < data->n; ++i)
        {
          double year;

          if (data->flags[i])
            continue;

          year = satdata_epoch2year(data->t[i]);
          if (year < thresh)
            data->flags[i] |= SATDATA_FLG_OUTLIER;
        }

      fprintf(stderr, "done\n");
    }
#endif

  main_proc(time_period, data);

  return 0;
}
