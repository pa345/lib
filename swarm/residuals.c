/*
 * residuals.c
 *
 * Compare main field model with SWARM data
 *
 * usage: residuals -c ascii_coef_file -i magdata_file
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
#include <gsl/gsl_rstat.h>

#include <mainlib/ml_msynth.h>
#include <mainlib/ml_magdata.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_solarpos.h>
#include <mainlib/ml_track.h>

#define IGRF_EPOCH                (2015.0)

int
compute_rms(const char *outfile, const double tmin, const double tmax, msynth_workspace *w,
            const magdata *data)
{
  size_t i;
  FILE *fp = fopen(outfile, "w");
  msynth_workspace *msynth_p;
  gsl_rstat_workspace * rstat_x = gsl_rstat_alloc();
  gsl_rstat_workspace * rstat_y = gsl_rstat_alloc();
  gsl_rstat_workspace * rstat_z = gsl_rstat_alloc();
  gsl_rstat_workspace * rstat_f = gsl_rstat_alloc();
  struct timeval tv0, tv1;

  i = 1;
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: QD latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: residual X (nT)\n", i++);
  fprintf(fp, "# Field %zu: residual Y (nT)\n", i++);
  fprintf(fp, "# Field %zu: residual Z (nT)\n", i++);
  fprintf(fp, "# Field %zu: residual F (nT)\n", i++);

#if 1
  msynth_p = msynth_copy(w);
#else
  msynth_p = msynth_shc_read(MSYNTH_CHAOS_FILE);
#endif

  msynth_set(1, 15, msynth_p);

  fprintf(stderr, "\n\t calculating residuals...");
  gettimeofday(&tv0, NULL);

  for (i = 0; i < data->n; ++i)
    {
      double t;
      double r = data->r[i];
      double theta = data->theta[i];
      double phi = data->phi[i];
      double B_NEC[3], B_core[4], B_model[4], F;

      /* ignore flagged data */
      if (!MAGDATA_ExistVector(data->flags[i]) || !MAGDATA_ExistScalar(data->flags[i]))
        continue;

      t = satdata_epoch2year(data->t[i]);
      if (tmin > 0.0 && t < tmin)
        continue;
      if (tmax > 0.0 && t > tmax)
        continue;

      msynth_eval(t, r, theta, phi, B_core, msynth_p);

      B_model[0] = B_core[0] + data->Bx_model[i];
      B_model[1] = B_core[1] + data->By_model[i];
      B_model[2] = B_core[2] + data->Bz_model[i];
      B_model[3] = gsl_hypot3(B_model[0], B_model[1], B_model[2]);

      B_NEC[0] = data->Bx_nec[i];
      B_NEC[1] = data->By_nec[i];
      B_NEC[2] = data->Bz_nec[i];
      F = data->F[i];

      if (fabs(data->qdlat[i]) < 55.0)
        {
          gsl_rstat_add(B_NEC[0] - B_model[0], rstat_x);
          gsl_rstat_add(B_NEC[1] - B_model[1], rstat_y);
          gsl_rstat_add(B_NEC[2] - B_model[2], rstat_z);
          gsl_rstat_add(F - B_model[3], rstat_f);
        }

      fprintf(fp, "%f %f %f %.12e %.12e %.12e %.12e\n",
              data->phi[i] * 180.0 / M_PI,
              90.0 - data->theta[i] * 180.0 / M_PI,
              data->qdlat[i],
              B_NEC[0] - B_model[0],
              B_NEC[1] - B_model[1],
              B_NEC[2] - B_model[2],
              F - B_model[3]);
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

  fprintf(stderr, "%20s %20s %20s %20s %20s\n", "Component", "N", "mean (nT)", "sigma (nT)", "rms (nT)");

  fprintf(stderr, "%20s %20zu %20.2f %20.2f %20.2f\n",
          "X",
          gsl_rstat_n(rstat_x),
          gsl_rstat_mean(rstat_x),
          gsl_rstat_sd(rstat_x),
          gsl_rstat_rms(rstat_x));

  fprintf(stderr, "%20s %20zu %20.2f %20.2f %20.2f\n",
          "Y",
          gsl_rstat_n(rstat_y),
          gsl_rstat_mean(rstat_y),
          gsl_rstat_sd(rstat_y),
          gsl_rstat_rms(rstat_y));

  fprintf(stderr, "%20s %20zu %20.2f %20.2f %20.2f\n",
          "Z",
          gsl_rstat_n(rstat_z),
          gsl_rstat_mean(rstat_z),
          gsl_rstat_sd(rstat_z),
          gsl_rstat_rms(rstat_z));

  fprintf(stderr, "%20s %20zu %20.2f %20.2f %20.2f\n",
          "F",
          gsl_rstat_n(rstat_f),
          gsl_rstat_mean(rstat_f),
          gsl_rstat_sd(rstat_f),
          gsl_rstat_rms(rstat_f));

  fprintf(stderr, "\t residual file = %s\n", outfile);

  fclose(fp);
  msynth_free(msynth_p);
  gsl_rstat_free(rstat_x);
  gsl_rstat_free(rstat_y);
  gsl_rstat_free(rstat_z);

  return GSL_SUCCESS;
}

int
main(int argc, char *argv[])
{
  int c;
  magdata * data = NULL;
  char * outfile = "res.xyz";
  msynth_workspace *msynth_p = NULL;
#if 0
  double tmin = -1.0;
  double tmax = -1.0;
#elif 1
  double tmin = get_year(1418601600); /* Dec 15 2014 00:00:00 UTC */
  double tmax = get_year(1421280000); /* Jan 15 2015 00:00:00 UTC */
#elif 0
  double tmin = get_year(1416009600); /* Nov 15 2014 00:00:00 UTC */
  double tmax = get_year(1423958400); /* Feb 15 2015 00:00:00 UTC */
#endif
  struct timeval tv0, tv1;

  while ((c = getopt(argc, argv, "a:b:i:c:wgm:o:")) != (-1))
    {
      switch (c)
        {
          case 'a':
            tmin = atof(optarg);
            break;

          case 'b':
            tmax = atof(optarg);
            break;

          case 'i':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);

            data = magdata_read(optarg, NULL);
            if (!data)
              {
                fprintf(stderr, "main: error reading %s\n", optarg);
                exit(1);
              }

            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%zu records read, %g seconds)\n", data->n,
                    time_diff(tv0, tv1));
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
            fprintf(stderr, "main: reading IGRF MF coefficients from %s...", optarg);
            msynth_p = msynth_igrf_read_mf(optarg, IGRF_EPOCH);
            fprintf(stderr, "done\n");
            break;

          case 'o':
            outfile = optarg;
            break;

          default:
            break;
        }
    }

  if (!data || !msynth_p)
    {
      fprintf(stderr, "usage: %s <-i magdata_file> <-c ascii_coef_file> [-a tmin_years] [-b tmax_years] [-w] [-g] [-m igrf_mf_candidate] [-o output_file]\n", argv[0]);
      exit(1);
    }

  if (tmin > 0.0)
    fprintf(stderr, "main: tmin = %g\n", tmin);
  if (tmax > 0.0)
    fprintf(stderr, "main: tmax = %g\n", tmax);

  {
    struct timeval tv0, tv1;

    fprintf(stderr, "main: computing rms of Model/Swarm...");
    gettimeofday(&tv0, NULL);
    compute_rms(outfile, tmin, tmax, msynth_p, data);
    gettimeofday(&tv1, NULL);
    fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));
  }

  msynth_free(msynth_p);
  magdata_free(data);

  return 0;
}
