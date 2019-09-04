/*
 * calc_euler.c
 *
 * Compute Euler angle corrections to Swarm data every 3 days
 *
 * Usage:
 * ./calc_euler <-i swarm_idx_file> [-q]
 *   q: use star camera quaternions for VFM->NEC instead of s1,s2,s3 basis
 *
 * Output is a data file 'euler.dat' containing the time series of computed
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

#define ALPHA_SWA       (11.8127)
#define BETA_SWA        (-76.1945)
#define GAMMA_SWA       (-12.5843)

#define ALPHA_SWB       (-8.8763)
#define BETA_SWB        (-76.4296)
#define GAMMA_SWB       (-9.1003)

#define ALPHA_SWC       (1.8232)
#define BETA_SWC        (-76.8283)
#define GAMMA_SWC       (-1.9911)

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

int
main_proc(const int use_altbasis, const double period, satdata_mag *data)
{
  int s = 0;
  size_t i;
  gsl_vector *euler_sat; /* Euler angles for VFM->SAT transformation */
  char *eulersat_file = "eulersat.dat";
  char *outfile = "euler.dat";
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

  euler_sat = gsl_vector_alloc(3);

  /* use s1,s2,s3 basis instead of star camera */
  if (use_altbasis)
    flags = EULER_FLG_ALTBASIS;

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
      s = eulersat_proc(flags, euler_sat, data);
      if (s)
        continue; /* insufficient number of data */

      /* save Euler angles in matrix A */
      gsl_matrix_set(A, i, 0, satdata_epoch2year(t + 0.5*dt));
      gsl_matrix_set(A, i, 1, wrap180(gsl_vector_get(euler_sat, 0) * 180.0 / M_PI));
      gsl_matrix_set(A, i, 2, wrap180(gsl_vector_get(euler_sat, 1) * 180.0 / M_PI));
      gsl_matrix_set(A, i, 3, wrap180(gsl_vector_get(euler_sat, 2) * 180.0 / M_PI));
    }

  fprintf(stderr, "main_proc: writing euler angles to %s...", outfile);
  print_euler_angles(outfile, A);
  fprintf(stderr, "done\n");

#if 0
  fprintf(stderr, "main_proc: printing eulersat residuals to %s...",
          eulersat_file);
  eulersat_print_residuals(eulersat_file, flags, euler_sat, data);
  fprintf(stderr, "done\n");
#endif

  gsl_vector_free(euler_sat);
  gsl_matrix_free(A);

  return s;
}

int
main(int argc, char *argv[])
{
  int c;
  satdata_mag *data = NULL;
  struct timeval tv0, tv1;
  int use_altbasis = 1;
  double time_period = 3.0; /* time window in days */

  while ((c = getopt(argc, argv, "i:qp:")) != (-1))
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

          case 'q':
            use_altbasis = 0;
            break;

          case 'p':
            time_period = atof(optarg);
            break;

          default:
            fprintf(stderr, "Usage: %s <-i swarm_idx_file> [-q]\n", argv[0]);
            exit(1);
            break;
        }
    }

  if (!data)
    {
      fprintf(stderr, "Usage: %s <-i swarm_idx_file> [-q]\n", argv[0]);
      exit(1);
    }

  if (use_altbasis)
    fprintf(stderr, "main: using s1,s2,s3 basis\n");
  else
    fprintf(stderr, "main: using star camera quaternions\n");

  fprintf(stderr, "main: time window = %g [days]\n", time_period);

#if 0
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
#endif

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

  main_proc(use_altbasis, time_period, data);

  return 0;
}
