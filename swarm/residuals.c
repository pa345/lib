/*
 * residuals.c
 *
 * Compare main field model with SWARM data
 *
 * usage: residuals -c ascii_coef_file -i swarm_index_file
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <assert.h>
#include <time.h>
#include <omp.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_msynth.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_solarpos.h>
#include <mainlib/ml_track.h>

#define SUBTRACT_MF7              1
#define DOWNSAMPLE                20

int
compute_rms(const char *outfile, msynth_workspace *w,
            const satdata_mag *data, double rms[4])
{
  size_t i;
  size_t n = 0;
  FILE *fp = fopen(outfile, "w");
  msynth_workspace *msynth_p;
  double *resid_f, *resid_x, *resid_y, *resid_z;
  struct timeval tv0, tv1;

  i = 1;
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: residual F (nT)\n", i++);
  fprintf(fp, "# Field %zu: residual X (nT)\n", i++);
  fprintf(fp, "# Field %zu: residual Y (nT)\n", i++);
  fprintf(fp, "# Field %zu: residual Z (nT)\n", i++);

  msynth_p = msynth_copy(w);

#if SUBTRACT_MF7
  msynth_set(1, 15, msynth_p);
#endif

  rms[0] = rms[1] = rms[2] = rms[3] = 0.0;

  /* initialize residuals array */
  resid_f = malloc(data->n * sizeof(double));
  resid_x = malloc(data->n * sizeof(double));
  resid_y = malloc(data->n * sizeof(double));
  resid_z = malloc(data->n * sizeof(double));
  for (i = 0; i < data->n; ++i)
    {
      resid_f[i] = -1.0e9;
      resid_x[i] = -1.0e9;
      resid_y[i] = -1.0e9;
      resid_z[i] = -1.0e9;
    }

  fprintf(stderr, "\n\t calculating residuals...");
  gettimeofday(&tv0, NULL);

  /* downsample to improve speed */
  for (i = 0; i < data->n; ++i)
    {
      size_t j;
      double t;
      double r = data->altitude[i] + data->R;
      double theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
      double phi = data->longitude[i] * M_PI / 180.0;
      double B_int[4], B_ext[4], B_crust[3], B_total[4];

      /* ignore flagged data */
      if (data->flags[i])
        continue;

      t = satdata_epoch2year(data->t[i]);

      msynth_eval(t, r, theta, phi, B_int, msynth_p);

      B_ext[0] = SATDATA_VEC_X(data->B_ext, i);
      B_ext[1] = SATDATA_VEC_Y(data->B_ext, i);
      B_ext[2] = SATDATA_VEC_Z(data->B_ext, i);

#if SUBTRACT_MF7
      B_crust[0] = SATDATA_VEC_X(data->B_crust, i);
      B_crust[1] = SATDATA_VEC_Y(data->B_crust, i);
      B_crust[2] = SATDATA_VEC_Z(data->B_crust, i);
#else
      B_crust[0] = B_crust[1] = B_crust[2] = 0.0;
#endif

      for (j = 0; j < 3; ++j)
        B_total[j] = B_int[j] + B_crust[j] + B_ext[j];

      B_total[3] = gsl_hypot3(B_total[0], B_total[1], B_total[2]);

      resid_x[i] = SATDATA_VEC_X(data->B, i) - B_total[0];
      resid_y[i] = SATDATA_VEC_Y(data->B, i) - B_total[1];
      resid_z[i] = SATDATA_VEC_Z(data->B, i) - B_total[2];
      resid_f[i] = data->F[i] - B_total[3];

      /* don't include high latitudes in rms */
      if (fabs(data->latitude[i]) <= 60.0)
        {
          rms[0] += pow(resid_x[i], 2.0);
          rms[1] += pow(resid_y[i], 2.0);
          rms[2] += pow(resid_z[i], 2.0);
          rms[3] += pow(resid_f[i], 2.0);
          ++n;
        }
    }

  for (i = 0; i < 4; ++i)
    rms[i] = sqrt(rms[i] / n);

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "\t done (%g seconds)", time_diff(tv0, tv1));
  fprintf(stderr, "\t ndata = %zu\n", n);
  fprintf(stderr, "\t rms X = %f [nT]\n", rms[0]);
  fprintf(stderr, "\t rms Y = %f [nT]\n", rms[1]);
  fprintf(stderr, "\t rms Z = %f [nT]\n", rms[2]);
  fprintf(stderr, "\t rms F = %f [nT]\n", rms[3]);

  /* write residual file */
  for (i = 0; i < data->n; ++i)
    {
      if (resid_f[i] < -1.0e8)
        continue; /* not calculated for this index */

      fprintf(fp, "%f %f %.12e %.12e %.12e %.12e\n",
              data->longitude[i],
              data->latitude[i],
              resid_f[i],
              resid_x[i],
              resid_y[i],
              resid_z[i]);
    }

  fprintf(stderr, "\t residual file = %s\n", outfile);

  fclose(fp);
  msynth_free(msynth_p);

  return GSL_SUCCESS;
}

/* flag any data within [lt_min,lt_max] */
size_t
flag_local_time(double lt_min, double lt_max, satdata_mag *data)
{
  size_t i;
  size_t cnt = 0;
  solarpos_workspace *work_sp = solarpos_alloc();

  for (i = 0; i < data->n; ++i)
    {
      time_t unix_time;
      double lt;

      /* ignore already flagged data for speed improvement */
      if (data->flags[i])
        continue;

      unix_time = satdata_epoch2timet(data->t[i]);
      lt = get_localtime(unix_time, data->longitude[i]*M_PI/180.0);

      /*
       * at high latitudes, calculate solar zenith angle for a better
       * measure of darkness/sunlight than local time
       */
      if (fabs(data->latitude[i]) > 60.0)
        {
          double lat = data->latitude[i] * M_PI / 180.0;
          double lon = wrappi(data->longitude[i] * M_PI / 180.0);
          double zenith;

          solarpos_calc_zenith(unix_time, lat, lon, &zenith, work_sp);
          zenith *= 180.0 / M_PI;
          assert(zenith >= 0.0);

          /* large zenith angle means darkness so allow it */
          if (zenith > 95.0)
            continue;
        }

      if (lt >= lt_min && lt <= lt_max)
        {
          ++cnt;
          data->flags[i] |= SATDATA_FLG_LT;
        }
    }

  solarpos_free(work_sp);

  return cnt;
}

int
main(int argc, char *argv[])
{
  int c;
  char *infile = NULL;
  satdata_mag *data;
  msynth_workspace *msynth_p = NULL;

  while ((c = getopt(argc, argv, "i:c:wgm:v:")) != (-1))
    {
      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          case 'c':
            fprintf(stderr, "main: reading %s...", optarg);
            msynth_p = msynth_read(optarg);
            fprintf(stderr, "done\n");
            break;

          case 'w':
            fprintf(stderr, "main: reading WMM coefficients...");
            msynth_p = msynth_wmm_read(MSYNTH_WMM_FILE);
            fprintf(stderr, "done\n");
            break;

          case 'g':
            fprintf(stderr, "main: reading IGRF coefficients...");
            msynth_p = msynth_igrf_read(MSYNTH_IGRF_FILE);
            fprintf(stderr, "done\n");
            break;

          case 'm':
            fprintf(stderr, "main: reading IGRF MF coefficients...");
            msynth_p = msynth_igrf_read_mf(optarg);
            fprintf(stderr, "done\n");
            break;

          case 'v':
            if (msynth_p == NULL)
              {
                fprintf(stderr, "main: error: must provide IGRF main field file\n");
                exit(1);
              }
            fprintf(stderr, "main: reading IGRF SV coefficients...");
            msynth_igrf_read_sv(optarg, msynth_p);
            fprintf(stderr, "done\n");
            break;

          default:
            break;
        }
    }

  if (!infile || !msynth_p)
    {
      fprintf(stderr, "usage: %s <-i swarm_index_file> <-c ascii_coef_file> [-w] [-g] [-m igrf_mf_candidate] [-v igrf_sv_candidate]\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "main: loading Swarm data from %s...", infile);
  data = satdata_swarm_read_idx(infile, 0);
  fprintf(stderr, "done (%zu data read)\n", data->n);

  {
    size_t i;

    fprintf(stderr, "main: downsampling data by factor %d...", DOWNSAMPLE);

    for (i = 0; i < data->n; ++i)
      {
        if (i % DOWNSAMPLE != 0)
          data->flags[i] |= SATDATA_FLG_OUTLIER;
      }

    fprintf(stderr, "done\n");
  }

  /* flag local time */
  {
    size_t nlt;
    const double lt_min = 4.0;
    const double lt_max = 21.0;

    fprintf(stderr, "main: flagging points inside LT window [%g,%g]...",
            lt_min, lt_max);

    nlt = flag_local_time(lt_min, lt_max, data);

    fprintf(stderr, "done (%zu/%zu data flagged)\n", nlt, data->n);
  }

  {
    size_t nflagged = satdata_nflagged(data);
    fprintf(stderr, "main: total flagged points: %zu/%zu\n",
            nflagged, data->n);
  }

  {
    double rms[4];
    struct timeval tv0, tv1;

    fprintf(stderr, "main: computing rms of Model/Swarm...");
    gettimeofday(&tv0, NULL);
    compute_rms("res.xyz", msynth_p, data, rms);
    gettimeofday(&tv1, NULL);
    fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));
  }

  msynth_free(msynth_p);

  return 0;
}
