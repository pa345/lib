/*
 * crust.c
 *
 * Compare Swarm data with MF7 model
 *
 * ./crust <-i swarm_index_file>
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
#include <gsl/gsl_interp.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_indices.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_euler.h>
#include <mainlib/ml_msynth.h>
#include <mainlib/ml_pomme.h>
#include <mainlib/ml_track.h>

#if 0
/* Bangui ranges */
#define LON_MIN              (16.0)
#define LON_MAX              (22.0)
#define LAT_MIN              (-10.0)
#define LAT_MAX              (20.0)
#else
/* Kursk ranges */
#define LON_MIN              (33.0)
#define LON_MAX              (39.0)
#define LAT_MIN              (44.0)
#define LAT_MAX              (60.0)
#endif

#define USE_GRAD_DIFF        0

typedef struct
{
  double lt_min;
  double lt_max;
  double lon_min;
  double lon_max;
  double lat_min;
  double lat_max;
  double kp_min;
  double kp_max;
} print_parameters;

int
check_pomme_jump(const size_t idx,
                 const track_workspace *track_p,
                 const satdata_mag *data)
{
  int s = 0;
  const size_t start_idx = track_p->tracks[idx].start_idx;
  const size_t end_idx = track_p->tracks[idx].end_idx;
  const double thresh = 2.0;
  size_t i, j;

  for (i = start_idx; i < end_idx; ++i)
    {
      double B_ext[3], B_ext_next[3];

      B_ext[0] = SATDATA_VEC_X(data->B_ext, i);
      B_ext[1] = SATDATA_VEC_Y(data->B_ext, i);
      B_ext[2] = SATDATA_VEC_Z(data->B_ext, i);

      B_ext_next[0] = SATDATA_VEC_X(data->B_ext, i + 1);
      B_ext_next[1] = SATDATA_VEC_Y(data->B_ext, i + 1);
      B_ext_next[2] = SATDATA_VEC_Z(data->B_ext, i + 1);

      for (j = 0; j < 3; ++j)
        {
          if (fabs(B_ext_next[j] - B_ext[j]) > thresh)
            return 1;
        }
    }

  return s;
}

/* res = B_obs - B_main */
int
calc_residual(const size_t idx, double B[3], const satdata_mag *data)
{
  double B_obs[3], B_main[3], B_ext[3];
  size_t k;

  B_obs[0] = SATDATA_VEC_X(data->B, idx);
  B_obs[1] = SATDATA_VEC_Y(data->B, idx);
  B_obs[2] = SATDATA_VEC_Z(data->B, idx);

  B_main[0] = SATDATA_VEC_X(data->B_main, idx);
  B_main[1] = SATDATA_VEC_Y(data->B_main, idx);
  B_main[2] = SATDATA_VEC_Z(data->B_main, idx);

#if USE_GRAD_DIFF
  B_ext[0] = B_ext[1] = B_ext[2] = 0.0;
#else
  B_ext[0] = SATDATA_VEC_X(data->B_ext, idx);
  B_ext[1] = SATDATA_VEC_Y(data->B_ext, idx);
  B_ext[2] = SATDATA_VEC_Z(data->B_ext, idx);
#endif

  for (k = 0; k < 3; ++k)
    B[k] = B_obs[k] - B_main[k] - B_ext[k];

  return 0;
}

int
compute_offsets(const size_t idx, double offset[3],
                const track_workspace *track_p,
                const satdata_mag *data)
{
  int s = 0;
  const size_t start_idx = track_p->tracks[idx].start_idx;
  const size_t end_idx = track_p->tracks[idx].end_idx;
  gsl_matrix *A = gsl_matrix_alloc(10000, 3);
  gsl_vector_view v;
  size_t i;
  size_t row = 0;

  for (i = start_idx; i <= end_idx; ++i)
    {
      double B_res[3];

      if (data->latitude[i] < LAT_MIN || data->latitude[i] > LAT_MAX)
        continue;

      /* compute B_obs - B_main - B_ext */
      calc_residual(i, B_res, data);

      /* subtract B_crust */
      B_res[0] -= SATDATA_VEC_X(data->B_crust, i);
      B_res[1] -= SATDATA_VEC_Y(data->B_crust, i);
      B_res[2] -= SATDATA_VEC_Z(data->B_crust, i);

      gsl_matrix_set(A, row, 0, B_res[0]);
      gsl_matrix_set(A, row, 1, B_res[1]);
      gsl_matrix_set(A, row, 2, B_res[2]);
      ++row;
    }

  v = gsl_matrix_subcolumn(A, 0, 0, row);
  offset[0] = gsl_stats_mean(v.vector.data, v.vector.stride, row);

  v = gsl_matrix_subcolumn(A, 1, 0, row);
  offset[1] = gsl_stats_mean(v.vector.data, v.vector.stride, row);

  v = gsl_matrix_subcolumn(A, 2, 0, row);
  offset[2] = gsl_stats_mean(v.vector.data, v.vector.stride, row);

  gsl_matrix_free(A);

  return s;
}

size_t
print_track_single(FILE *fp, const size_t idx,
                   const double lat_min, const double lat_max,
                   double ssq[3],
                   double ssq_grad[3],
                   track_workspace *track_p,
                   const satdata_mag *data)
{
  const size_t start_idx = track_p->tracks[idx].start_idx;
  const size_t end_idx = track_p->tracks[idx].end_idx;
  size_t i, j, k;
  f107_workspace *f107_p = f107_alloc(F107_IDX_FILE);
  size_t grad_idx = 15;
  const double dt_min = (double) grad_idx - 2.0; /* time interval in s */
  const double dt_max = (double) grad_idx + 2.0;
  size_t ndata = 0;
  int sgn; /* sign for gradient/difference calculation */
  double offset[3] = { 0.0, 0.0, 0.0 };

  /* compute external field model correction for this track */
  track_dstcorr(idx, data, track_p);

  /* calculate component offsets between this track and MF7 */
  compute_offsets(idx, offset, track_p, data);

  for (i = start_idx; i <= end_idx; ++i)
    {
      time_t unix_time = satdata_epoch2timet(data->t[i]);
      double phi = data->longitude[i] * M_PI / 180.0;
      double lt = get_localtime(unix_time, phi);
      double qdlat = data->qdlat[i];
      double euvac, dt;
      double B_crust[3], B_res[3], B_res_corr[3];
      double B_crust2[3], B_res2[3], dB[3], dB_crust[3];

      if (data->latitude[i] < lat_min || data->latitude[i] > lat_max)
        continue;

      /* look for along-track difference */
      j = GSL_MIN(i + grad_idx, data->n - 1);
      dt = (data->t[j] - data->t[i]) / 1000.0; /* in s */
      if (dt < dt_min || dt > dt_max)
        continue; /* no gradient available */

      sgn = (data->latitude[j] > data->latitude[i]) ? 1 : -1;

      f107_get_euvac(unix_time, &euvac, f107_p);

      calc_residual(i, B_res, data);
      calc_residual(j, B_res2, data);

      /* compute residual including external field correction */
      track_residual(idx, i, B_res_corr, data, track_p);

      B_crust[0] = SATDATA_VEC_X(data->B_crust, i);
      B_crust[1] = SATDATA_VEC_Y(data->B_crust, i);
      B_crust[2] = SATDATA_VEC_Z(data->B_crust, i);

      B_crust2[0] = SATDATA_VEC_X(data->B_crust, j);
      B_crust2[1] = SATDATA_VEC_Y(data->B_crust, j);
      B_crust2[2] = SATDATA_VEC_Z(data->B_crust, j);

      for (k = 0; k < 3; ++k)
        {
          dB[k] = sgn * (B_res2[k] - B_res[k]);
          dB_crust[k] = sgn * (B_crust2[k] - B_crust[k]);

          /* update sum of squares difference between dB and dB_crust */
          ssq_grad[k] += pow(dB_crust[k] - dB[k], 2.0);

          ssq[k] += pow(B_res[k] - offset[k] - B_crust[k], 2.0);
        }

      ++ndata;

      printf("%ld %.8f %.2f %.2f %6.2f %5.1f %10.4f %8.4f %10.4f %10.4f %d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",
             unix_time,
             satdata_epoch2year(data->t[i]),
             get_ut(unix_time),
             lt,
             get_season(unix_time),
             euvac,
             data->longitude[i],
             data->latitude[i],
             data->altitude[i],
             qdlat,
             satdata_mag_satdir(i, data),
             B_res[0],
             B_res[1],
             B_res[2],
#if 1
             B_res[0] - offset[0],
             B_res[1] - offset[1],
             B_res[2] - offset[2],
#else
             B_res_corr[0] + B_crust[0],
             B_res_corr[1] + B_crust[1],
             B_res_corr[2] + B_crust[2],
#endif
             B_crust[0],
             B_crust[1],
             B_crust[2]);
    }

  fprintf(fp, "\n\n");

  f107_free(f107_p);

  return ndata;
} /* print_track_single() */

int
print_tracks(const print_parameters *params,
             track_workspace *track_p, const satdata_mag *data)
{
  int s = 0;
  size_t i = 1;
  size_t nrej_lt = 0;
  size_t nrej_rms = 0;
  size_t nrej_lon = 0;
  size_t nrej_jump = 0;
  size_t nrej_kp = 0;
  size_t nrej_tot = 0;
  FILE *fp = stdout;
  double rms[3] = { 0.0, 0.0, 0.0 };
  double rms_grad[3] = { 0.0, 0.0, 0.0 };
  size_t ndata = 0; /* number of data used for rms */
  kp_workspace *kp_p = kp_alloc(KP_IDX_FILE);

  fprintf(fp, "# Field %zu: time (UT)\n", i++);
  fprintf(fp, "# Field %zu: time (decimal year)\n", i++);
  fprintf(fp, "# Field %zu: UT (hours)\n", i++);
  fprintf(fp, "# Field %zu: local time (hours)\n", i++);
  fprintf(fp, "# Field %zu: season (day of year)\n", i++);
  fprintf(fp, "# Field %zu: EUVAC\n", i++);
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: altitude (km)\n", i++);
  fprintf(fp, "# Field %zu: QD latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: satellite direction\n", i++);
  fprintf(fp, "# Field %zu: X residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: Y residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: Z residual (nT)\n", i++);
  fprintf(fp, "# Field %zu: X residual (Dst corrected) (nT)\n", i++);
  fprintf(fp, "# Field %zu: Y residual (Dst corrected) (nT)\n", i++);
  fprintf(fp, "# Field %zu: Z residual (Dst corrected) (nT)\n", i++);
  fprintf(fp, "# Field %zu: X crust (nT)\n", i++);
  fprintf(fp, "# Field %zu: Y crust (nT)\n", i++);
  fprintf(fp, "# Field %zu: Z crust (nT)\n", i++);

  for (i = 0; i < track_p->n; ++i)
    {
      track_data *tptr = &(track_p->tracks[i]);
      time_t unix_time = satdata_epoch2timet(tptr->t_eq);
      double kp;
      double a, b;

      /* check for rms flag */
      if (tptr->flags)
        {
          ++nrej_rms;
          continue;
        }

      /* check local time */
      a = fmod(tptr->lt_eq - params->lt_min, 24.0);
      b = fmod(params->lt_max - params->lt_min, 24.0);

      if (a < 0.0)
        a += 24.0;
      if (b < 0.0)
        b += 24.0;

      if (a > b)
        {
          ++nrej_lt;
          continue;
        }
      else if (tptr->lon_eq < params->lon_min || tptr->lon_eq > params->lon_max)
        {
          ++nrej_lon;
          continue;
        }

      kp_get(unix_time, &kp, kp_p);
      if (kp < params->kp_min || kp > params->kp_max)
        {
          ++nrej_kp;
          continue;
        }

      /* only need to check POMME for vector measurements (not differences) */
#if !USE_GRAD_DIFF
      s = check_pomme_jump(i, track_p, data);
      if (s)
        {
          ++nrej_jump;
          continue;
        }
#endif

      ndata += print_track_single(fp, i, params->lat_min, params->lat_max, rms, rms_grad, track_p, data);
    }

  if (ndata > 0)
    {
      for (i = 0; i < 3; ++i)
        {
          rms[i] = sqrt(rms[i] / (double)ndata);
          rms_grad[i] = sqrt(rms_grad[i] / (double)ndata);
        }
    }

  nrej_tot = nrej_lt + nrej_lon + nrej_jump + nrej_kp + nrej_rms;

  fprintf(stderr, "print_tracks: rejected %zu/%zu tracks due to rms\n",
          nrej_rms, track_p->n);
  fprintf(stderr, "print_tracks: rejected %zu/%zu tracks due to LT\n",
          nrej_lt, track_p->n);
  fprintf(stderr, "print_tracks: rejected %zu/%zu tracks due to longitude\n",
          nrej_lon, track_p->n);
  fprintf(stderr, "print_tracks: rejected %zu/%zu tracks due to kp\n",
          nrej_kp, track_p->n);
  fprintf(stderr, "print_tracks: rejected %zu/%zu tracks due to POMME jumps\n",
          nrej_jump, track_p->n);
  fprintf(stderr, "print_tracks: rejected %zu/%zu total tracks (%zu remaining)\n",
          nrej_tot, track_p->n, track_p->n - nrej_tot);

  fprintf(stderr, "print_tracks: data available for rms: %zu\n", ndata);
  fprintf(stderr, "print_tracks: X rms = %.3f [nT]\n", rms[0]);
  fprintf(stderr, "print_tracks: Y rms = %.3f [nT]\n", rms[1]);
  fprintf(stderr, "print_tracks: Z rms = %.3f [nT]\n", rms[2]);
  fprintf(stderr, "print_tracks: DX rms = %.3f [nT]\n", rms_grad[0]);
  fprintf(stderr, "print_tracks: DY rms = %.3f [nT]\n", rms_grad[1]);
  fprintf(stderr, "print_tracks: DZ rms = %.3f [nT]\n", rms_grad[2]);

  kp_free(kp_p);

  return s;
} /* print_tracks() */

int
main(int argc, char *argv[])
{
  satdata_mag *data;
  struct timeval tv0, tv1;
  char *infile = NULL;
  print_parameters params;
  track_workspace *track_p;
  euler_workspace *euler_p = NULL;

  params.lt_min = 23.0;
  params.lt_max = 5.0;
  params.lon_min = LON_MIN;
  params.lon_max = LON_MAX;
  params.lat_min = LAT_MIN;
  params.lat_max = LAT_MAX;
  params.kp_min = 0.0;
  params.kp_max = 1.0;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "lt_min", required_argument, NULL, 'a' },
          { "lt_max", required_argument, NULL, 'b' },
          { "kp_min", required_argument, NULL, 'f' },
          { "kp_max", required_argument, NULL, 'g' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "i:e:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          case 'a':
            params.lt_min = atof(optarg);
            break;

          case 'b':
            params.lt_max = atof(optarg);
            break;

          case 'e':
            fprintf(stderr, "main: reading Euler angles from %s...", optarg);
            euler_p = euler_read(optarg);
            if (!euler_p)
              exit(1);
            fprintf(stderr, "done (%zu sets of angles read)\n", euler_p->n);
            break;

          case 'f':
            params.kp_min = atof(optarg);
            break;

          case 'g':
            params.kp_max = atof(optarg);
            break;

          default:
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i swarm_index_file> [-e euler_file] [--lt_min lt_min] [--lt_max lt_max] [--kp_min kp_min] [--kp_max kp_max]\n",
              argv[0]);
      exit(1);
    }

  fprintf(stderr, "input file = %s\n", infile);
  fprintf(stderr, "LT min  = %.2f\n", params.lt_min);
  fprintf(stderr, "LT max  = %.2f\n", params.lt_max);
  fprintf(stderr, "kp min  = %.2f\n", params.kp_min);
  fprintf(stderr, "kp max  = %.2f\n", params.kp_max);

  fprintf(stderr, "Reading %s...", infile);
  gettimeofday(&tv0, NULL);

  data = satdata_swarm_read_idx(infile, 0);
  if (!data)
    {
      fprintf(stderr, "main: error reading %s\n", infile);
      exit(1);
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu records read, %g seconds)\n", data->n,
          time_diff(tv0, tv1));

  if (euler_p)
    {
      fprintf(stderr, "main: rotating VFM measurements with new Euler angles...");
      euler_apply(data, euler_p);
      fprintf(stderr, "done\n");
    }

  track_p = track_alloc();

  {
    size_t nflag, nrms;
    double thresh[] = { 20.0, 25.0, 15.0, 15.0 };

    fprintf(stderr, "main: separating tracks...");
    nflag = track_init(data, NULL, track_p);
    fprintf(stderr, "done (%zu/%zu (%.1f%%) points flagged, %zu tracks found)\n",
            nflag, data->n, (double)nflag / (double)data->n * 100.0,
            track_p->n);

    nrms = track_flag_rms("swarm_rms.dat", thresh, data, track_p);
    fprintf(stderr, "main: flagged (%zu/%zu) (%.1f%%) points due to high rms\n",
            nrms, data->n, (double) nrms / (double) data->n * 100.0);
  }

  print_tracks(&params, track_p, data);

  satdata_mag_free(data);
  track_free(track_p);

  return 0;
} /* main() */
