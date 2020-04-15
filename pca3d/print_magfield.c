/*
 * print_magfield.c
 *
 * Compute magnetic field from 3D TIEGCM current and print results
 *
 * ./print_magfield <-i tiegcm_nc_file>
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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rstat.h>

#include <mainlib/ml_common.h>
#include <mainlib/ml_bsearch.h>
#include <mainlib/ml_msynth.h>

#include <magfield/magfield.h>
#include <magfield/magfield_eval.h>

#include "magfit.h"
#include "tiegcm3d.h"

#define MAX_PTS            1000

/*
print_data()
  Print J grid for a fixed time and altitude

Inputs: data - tiegcm data
*/

int
print_data(const char *filename, const tiegcm3d_data *data, const int time_idx, const int ir, magfield_eval_workspace *w)
{
  int s = 0;
  const double r_B = R_EARTH_KM + 80.0;
  size_t i;
  size_t ilat, ilon;
  FILE *fp;

#if 0

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_data: unable to open %s: %s\n",
              filename, strerror(errno));
    }

  i = 1;
  fprintf(fp, "# Time: %ld (%.6f DOY)\n", data->t[time_idx], data->doy[time_idx] + data->ut[time_idx] / 24.0);
  fprintf(fp, "# Radius J: %.2f (km) [%.2f km altitude]\n", data->r[ir], data->r[ir] - R_EARTH_KM);
  fprintf(fp, "# Radius B: %.2f (km) [%.2f km altitude]\n", r_B, r_B - R_EARTH_KM);
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: data J_r (A/m^2)\n", i++);
  fprintf(fp, "# Field %zu: data J_t (A/m^2)\n", i++);
  fprintf(fp, "# Field %zu: data J_p (A/m^2)\n", i++);
  fprintf(fp, "# Field %zu: model J_r (A/m^2)\n", i++);
  fprintf(fp, "# Field %zu: model J_t (A/m^2)\n", i++);
  fprintf(fp, "# Field %zu: model J_p (A/m^2)\n", i++);
  fprintf(fp, "# Field %zu: B_r (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_t (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_p (nT)\n", i++);
  fprintf(fp, "# Field %zu: |B| (nT)\n", i++);

  for (ilon = 0; ilon < data->nlon; ++ilon)
    {
      double phi = data->glon[ilon] * M_PI / 180.0;

      for (ilat = 0; ilat < data->nlat; ++ilat)
        {
          size_t idx = TIEGCM3D_IDX(time_idx, ir, ilat, ilon, data);
          double theta = M_PI / 2.0 - data->glat[ilat] * M_PI / 180.0;
          double B[4], J[4];

          magfield_eval_B(r_B * 1.0e3, theta, phi, B, w);
          magfield_eval_J(data->r[ir] * 1.0e3, theta, phi, J, w);

          fprintf(fp, "%8.4f %8.4f %16.4e %16.4e %16.4e %16.4e %16.4e %16.4e %16.4e %16.4e %16.4e %16.4e\n",
                  data->glon[ilon],
                  data->glat[ilat],
                  data->Jr[idx],
                  data->Jt[idx],
                  data->Jp[idx],
                  J[0],
                  J[1],
                  J[2],
                  B[0] * 1.0e9,
                  B[1] * 1.0e9,
                  B[2] * 1.0e9,
                  B[3] * 1.0e9);
        }

      fprintf(fp, "\n");
    }

  fclose(fp);

#endif

#if 1

  {
    /* print out latitude profile at 80km and 280 deg longitude */
    const double epoch = 2020.0;
    const double r1 = R_EARTH_KM + 80.0;
    const double r2 = R_EARTH_KM + 450.0;
    const double r_J = R_EARTH_KM + 110.0;
    const double phi = 280.0 * M_PI / 180.0;
    const double dlat = 0.05;
    msynth_workspace * core_p = msynth_shc_read(MSYNTH_CHAOS_FILE);
    double lat;

    msynth_set(1, 15, core_p);
    fp = fopen("profile.txt", "w");

    i = 1;
    fprintf(fp, "# Time: %ld (%.6f DOY)\n", data->t[time_idx], data->doy[time_idx] + data->ut[time_idx] / 24.0);
    fprintf(fp, "# Longitude: %.2f (degrees)\n", phi * 180.0 / M_PI);
    fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
    fprintf(fp, "# Field %zu: J (A/m2) at %.2f (km)\n", i++, r_J - R_EARTH_KM);
    fprintf(fp, "# Field %zu: dF (nT) at %.2f (km)\n", i++, r1 - R_EARTH_KM);
    fprintf(fp, "# Field %zu: dF (nT) at %.2f (km)\n", i++, r2 - R_EARTH_KM);

    for (lat = -80.0; lat <= 80.0; lat += dlat)
      {
        double theta = M_PI / 2.0 - lat * M_PI / 180.0;
        double B1[4], B2[4], B_core1[4], B_core2[4], J[4];
        double tmp1[3], tmp2[3];
        double dF1, dF2;

        msynth_eval(epoch, r1, theta, phi, B_core1, core_p);
        msynth_eval(epoch, r2, theta, phi, B_core2, core_p);

        magfield_eval_B(r1 * 1.0e3, theta, phi, B1, w);
        magfield_eval_B(r2 * 1.0e3, theta, phi, B2, w);
        magfield_eval_J(r_J * 1.0e3, theta, phi, J, w);

        /* tmp1 = B_core1 + B1 */
        tmp1[0] = B_core1[0] - B1[1] * 1.0e9;
        tmp1[1] = B_core1[1] + B1[2] * 1.0e9;
        tmp1[2] = B_core1[2] - B1[0] * 1.0e9;

        /* tmp2 = B_core2 + B2 */
        tmp2[0] = B_core2[0] - B2[1] * 1.0e9;
        tmp2[1] = B_core2[1] + B2[2] * 1.0e9;
        tmp2[2] = B_core2[2] - B2[0] * 1.0e9;

        /* dF = | B_core + B_EEJ | - | B_core | */
        dF1 = gsl_hypot3(tmp1[0], tmp1[1], tmp1[2]) - B_core1[3];
        dF2 = gsl_hypot3(tmp2[0], tmp2[1], tmp2[2]) - B_core2[3];

        fprintf(fp, "%f %16.4e %16.4e %16.4e\n",
                lat,
                J[2],
                dF1,
                dF2);
      }

    fclose(fp);
  }

#endif

  return s;
}

/*
print_alt()
  Print altitude/latitude map of Jr/Jt/Jp for fixed time and longitude

Inputs: data - tiegcm data
*/

int
print_alt(const char *filename, const tiegcm3d_data *data, const int it, const int ilon, magfield_eval_workspace *w)
{
  int s = 0;
  const double phi = data->glon[ilon] * M_PI / 180.0;
  FILE *fp;
  size_t i;
  size_t ir, ilat;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_alt: unable to open %s: %s\n",
              filename, strerror(errno));
    }

  i = 1;
  fprintf(fp, "# Time: %ld (%.6f DOY)\n", data->t[it], data->doy[it] + data->ut[it] / 24.0);
  fprintf(fp, "# Longitude: %.2f (deg)\n", data->glon[ilon]);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: altitude (km)\n", i++);
  fprintf(fp, "# Field %zu: J_r (A/m^2)\n", i++);
  fprintf(fp, "# Field %zu: J_t (A/m^2)\n", i++);
  fprintf(fp, "# Field %zu: J_p (A/m^2)\n", i++);
  fprintf(fp, "# Field %zu: B_r (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_t (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_p (nT)\n", i++);
  fprintf(fp, "# Field %zu: |B| (nT)\n", i++);

  for (ilat = 0; ilat < data->nlat; ++ilat)
    {
      double theta = M_PI / 2.0 - data->glat[ilat] * M_PI / 180.0;

      for (ir = 0; ir < data->nr - 1; ++ir)
        {
          size_t idx = TIEGCM3D_IDX(it, ir, ilat, ilon, data);
          double B[4];

          magfield_eval_B(data->r[ir] * 1.0e3, theta, phi, B, w);

          fprintf(fp, "%8.4f %8.4f %16.4e %16.4e %16.4e %16.4e %16.4e %16.4e %16.4e\n",
                  data->glat[ilat],
                  data->r[ir] - R_EARTH_KM,
                  data->Jr[idx],
                  data->Jt[idx],
                  data->Jp[idx],
                  B[0] * 1.0e9,
                  B[1] * 1.0e9,
                  B[2] * 1.0e9,
                  B[3] * 1.0e9);
        }

      fprintf(fp, "\n");
    }

  fclose(fp);

  return s;
}

int
docalc(magfield_workspace * w)
{
  FILE *fp;
  const size_t nsat = 3;                  /* number of satellites in pearl-on-string */
  const double r_km = R_EARTH_KM + 85.0;
  const double r_m = r_km * 1.0e3;
  const double fac = 10.0;                /* fudge factor to make TIEGCM B field match better with ground variations */
  /*const double lon[4] = { 136.5, 141.0, 142.5, 145.5 };*/
  const double lon[4] = { 143.5, 148.5, 150.0, 152.5 };
  const double dlat = 8.0 / 110.0 * 10.0; /* 8 km/s / 110 km/deg =~ 0.072 degrees latitude per second */
  const double lat_min = -10.0;
  const double lat_max = 20.0;
  double phi[4];
  magfit_parameters params = magfit_default_parameters();
  magfit_workspace *magfit_p[4];
  double lat;
  size_t i, j, k, lonidx;
  double rnorm, snorm;
  gsl_rng *rng_p = gsl_rng_alloc(gsl_rng_default);
  const double sigma_r = 140.0;
  const double sigma_t = 64.0;
  size_t n = 0;                  /* number of points in track */
  double theta_array[MAX_PTS];   /* colatitudes of along-track positions */

  for (i = 0; i < 4; ++i)
    phi[i] = lon[i] * M_PI / 180.0;

  params.lat_min = -15.0;
  params.lat_max = 25.0;
  params.lat_spacing1d = 1.0;
  params.flags = MAGFIT_FLG_FIT_X | MAGFIT_FLG_FIT_Z;
  params.flags |= MAGFIT_FLG_SECS_FIT_DF;

  for (i = 0; i < 4; ++i)
    magfit_p[i] = magfit_alloc(magfit_secs1d, &params);

  /* store latitude points of each position along track */
  for (lat = lat_min; lat <= lat_max; lat += dlat)
    theta_array[n++] = M_PI / 2.0 - lat * M_PI / 180.0;

  fp = fopen("signal.dat", "w");

  for (k = 0; k < nsat; ++k)
    {
      for (lonidx = 0; lonidx < 4; ++lonidx)
        {
          if (k == 0)
            {
              i = 1;
              fprintf(fp, "# Field %zu: geocentric latitude (degrees)\n", i++);
              fprintf(fp, "# Field %zu: B_r original data for longitude %.2f degrees (nT)\n", i++, lon[lonidx]);
              fprintf(fp, "# Field %zu: B_r noise-added data for longitude %.2f degrees (nT)\n", i++, lon[lonidx]);
              fprintf(fp, "# Field %zu: B_t original data for longitude %.2f degrees (nT)\n", i++, lon[lonidx]);
              fprintf(fp, "# Field %zu: B_t noise-added data for longitude %.2f degrees (nT)\n", i++, lon[lonidx]);
            }

          for (i = 0; i < n; ++i)
            {
              double theta = theta_array[i];
              double B[4][4], B_noise[4][3], B_noise_nec[4][3];

              magfield_eval_B(r_m, theta, phi[lonidx], B[lonidx], w);

              /* convert to nT */
              for (j = 0; j < 3; ++j)
                B[lonidx][j] *= fac * 1.0e9;

              /* add noise */
              B_noise[lonidx][0] = B[lonidx][0] + gsl_ran_gaussian(rng_p, sigma_r);
              B_noise[lonidx][1] = B[lonidx][1] + gsl_ran_gaussian(rng_p, sigma_t);
              B_noise[lonidx][2] = B[lonidx][2];

              /* convert to NEC */
#if 1
              /* signal+noise */
              B_noise_nec[lonidx][0] = -B_noise[lonidx][1];
              B_noise_nec[lonidx][1] = B_noise[lonidx][2];
              B_noise_nec[lonidx][2] = -B_noise[lonidx][0];
#else
              /* signal, no noise */
              B_noise_nec[lonidx][0] = -B[lonidx][1];
              B_noise_nec[lonidx][1] = B[lonidx][2];
              B_noise_nec[lonidx][2] = -B[lonidx][0];
#endif

              magfit_add_datum(0.0, r_km, theta, phi[lonidx], 0.0, B_noise_nec[lonidx], magfit_p[lonidx]);

              if (k == 0)
                {
                  fprintf(fp, "%8.4f %8.4f %8.4f %8.4f %8.4f\n",
                          90.0 - theta * 180.0 / M_PI,
                          B[lonidx][0],
                          B_noise[lonidx][0],
                          B[lonidx][1],
                          B_noise[lonidx][1]);
                }
            }

          fprintf(fp, "\n\n");
        }
    }

  fclose(fp);

  /* fit EEJ model */
  for (i = 0; i < 4; ++i)
    magfit_fit(&rnorm, &snorm, magfit_p[i]);

  fp = fopen("fitted.dat", "w");

  for (lonidx = 0; lonidx < 4; ++lonidx)
    {
      i = 1;
      fprintf(fp, "# Field %zu: geocentric latitude (degrees)\n", i++);
      fprintf(fp, "# Field %zu: B_r fitted for longitude %.2f degrees (nT)\n", i++, lon[lonidx]);
      fprintf(fp, "# Field %zu: B_t fitted for longitude %.2f degrees (nT)\n", i++, lon[lonidx]);
      fprintf(fp, "# Field %zu: J_phi fitted for longitude %.2f degrees (A/m)\n", i++, lon[lonidx]);

      for (i = 0; i < n; ++i)
        {
          double theta = theta_array[i];
          double B1_fit_nec[3];
          double J1_nec[3];

          magfit_eval_B(0.0, r_km, theta, phi[lonidx], B1_fit_nec, magfit_p[lonidx]);
          magfit_eval_J(r_km, theta, phi[lonidx], J1_nec, magfit_p[lonidx]);

          fprintf(fp, "%8.4f %8.4f %8.4f %8.4f\n",
                  90.0 - theta * 180.0 / M_PI,
                  -B1_fit_nec[2],
                  -B1_fit_nec[0],
                  J1_nec[1]);
        }

      fprintf(fp, "\n\n");
    }

  fclose(fp);

  exit(1);

#if 0
  i = 1;
  fprintf(stdout, "# Field %zu: geocentric latitude (degrees)\n", i++);
  fprintf(stdout, "# Field %zu: B_r original data for longitude 145.69 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: B_t original data for longitude 145.69 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: B_r original data for longitude 149.32 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: B_t original data for longitude 149.32 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: B_r original data for longitude 150.68 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: B_t original data for longitude 150.68 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: B_r original data for longitude 152.95 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: B_t original data for longitude 152.95 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: B_r noisy data for longitude 145.69 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: B_t noisy data for longitude 145.69 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: B_r noisy data for longitude 149.32 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: B_t noisy data for longitude 149.32 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: B_r noisy data for longitude 150.68 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: B_t noisy data for longitude 150.68 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: B_r noisy data for longitude 152.95 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: B_t noisy data for longitude 152.95 degrees (nT)\n", i++);

  /* outer loop over number of satellite passes; inner loop over latitude range */
  for (k = 0; k < 3; ++k) {
  for (lat = -10.0; lat <= 20.0; lat += dlat)
    {
      double theta = M_PI / 2.0 - lat * M_PI / 180.0;
      double B1[4], B2[4], B3[4], B4[4];
      double B1_noise[4], B2_noise[4], B3_noise[4], B4_noise[4];
      double B1_nec[3], B2_nec[3], B3_nec[3], B4_nec[3];

      magfield_eval_B(r_m, theta, phi[0], B1, w);
      magfield_eval_B(r_m, theta, phi[1], B2, w);
      magfield_eval_B(r_m, theta, phi[2], B3, w);
      magfield_eval_B(r_m, theta, phi[3], B4, w);

      for (i = 0; i < 3; ++i)
        {
          B1[i] *= fac * 1.0e9;
          B2[i] *= fac * 1.0e9;
          B3[i] *= fac * 1.0e9;
          B4[i] *= fac * 1.0e9;

          B1_noise[i] = B1[i];
          B2_noise[i] = B2[i];
          B3_noise[i] = B3[i];
          B4_noise[i] = B4[i];
        }

      /* add noise */
#if 1
      B1_noise[0] += gsl_ran_gaussian(rng_p, sigma_Z);
      B2_noise[0] += gsl_ran_gaussian(rng_p, sigma_Z);
      B3_noise[0] += gsl_ran_gaussian(rng_p, sigma_Z);
      B4_noise[0] += gsl_ran_gaussian(rng_p, sigma_Z);

      B1_noise[1] += gsl_ran_gaussian(rng_p, sigma_X);
      B2_noise[1] += gsl_ran_gaussian(rng_p, sigma_X);
      B3_noise[1] += gsl_ran_gaussian(rng_p, sigma_X);
      B4_noise[1] += gsl_ran_gaussian(rng_p, sigma_X);
#endif

      /* convert spherical to NEC */

      B1_nec[0] = -B1_noise[1];
      B1_nec[1] = B1_noise[2];
      B1_nec[2] = -B1_noise[0];

      B2_nec[0] = -B2_noise[1];
      B2_nec[1] = B2_noise[2];
      B2_nec[2] = -B2_noise[0];

      B3_nec[0] = -B3_noise[1];
      B3_nec[1] = B3_noise[2];
      B3_nec[2] = -B3_noise[0];

      B4_nec[0] = -B4_noise[1];
      B4_nec[1] = B4_noise[2];
      B4_nec[2] = -B4_noise[0];

      magfit_add_datum(0.0, r_km, theta, phi[0], 0.0, B1_nec, magfit_p[0]);
      magfit_add_datum(0.0, r_km, theta, phi[1], 0.0, B2_nec, magfit_p[1]);
      magfit_add_datum(0.0, r_km, theta, phi[2], 0.0, B3_nec, magfit_p[2]);
      magfit_add_datum(0.0, r_km, theta, phi[3], 0.0, B4_nec, magfit_p[3]);

      if (k == 0)
        {
          printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
                 lat,
                 B1[0],
                 B1[1],
                 B2[0],
                 B2[1],
                 B3[0],
                 B3[1],
                 B4[0],
                 B4[1],
                 B1_noise[0],
                 B1_noise[1],
                 B2_noise[0],
                 B2_noise[1],
                 B3_noise[0],
                 B3_noise[1],
                 B4_noise[0],
                 B4_noise[1]);
        }
    }
    }

  /* fit EEJ model */
  for (i = 0; i < 4; ++i)
    magfit_fit(&rnorm, &snorm, magfit_p[i]);

  i = 1;
  fprintf(stdout, "\n\n");
  fprintf(stdout, "# Field %zu: geocentric latitude (degrees)\n", i++);
  fprintf(stdout, "# Field %zu: B_r fitted for longitude 145.69 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: B_t fitted for longitude 145.69 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: B_r fitted for longitude 149.32 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: B_t fitted for longitude 149.32 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: B_r fitted for longitude 150.68 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: B_t fitted for longitude 150.68 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: B_r fitted for longitude 152.95 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: B_t fitted for longitude 152.95 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: J_phi fitted for longitude 145.69 degrees (A/m)\n", i++);
  fprintf(stdout, "# Field %zu: J_phi fitted for longitude 149.32 degrees (A/m)\n", i++);
  fprintf(stdout, "# Field %zu: J_phi fitted for longitude 150.68 degrees (A/m)\n", i++);
  fprintf(stdout, "# Field %zu: J_phi fitted for longitude 152.95 degrees (A/m)\n", i++);

  for (lat = -10.0; lat <= 20.0; lat += dlat)
    {
      double theta = M_PI / 2.0 - lat * M_PI / 180.0;
      double B1_fit_nec[3], B2_fit_nec[3], B3_fit_nec[3], B4_fit_nec[3];
      double J1_nec[3], J2_nec[3], J3_nec[3], J4_nec[3];

      magfit_eval_B(0.0, r_km, theta, phi[0], B1_fit_nec, magfit_p[0]);
      magfit_eval_B(0.0, r_km, theta, phi[1], B2_fit_nec, magfit_p[1]);
      magfit_eval_B(0.0, r_km, theta, phi[2], B3_fit_nec, magfit_p[2]);
      magfit_eval_B(0.0, r_km, theta, phi[3], B4_fit_nec, magfit_p[3]);

      magfit_eval_J(r_km, theta, phi[0], J1_nec, magfit_p[0]);
      magfit_eval_J(r_km, theta, phi[1], J2_nec, magfit_p[1]);
      magfit_eval_J(r_km, theta, phi[2], J3_nec, magfit_p[2]);
      magfit_eval_J(r_km, theta, phi[3], J4_nec, magfit_p[3]);

      printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
             lat,
             -B1_fit_nec[2],
             -B1_fit_nec[0],
             -B2_fit_nec[2],
             -B2_fit_nec[0],
             -B3_fit_nec[2],
             -B3_fit_nec[0],
             -B4_fit_nec[2],
             -B4_fit_nec[0],
             J1_nec[1],
             J2_nec[1],
             J3_nec[1],
             J4_nec[1]);
    }
#endif

  for (i = 0; i < 4; ++i)
    magfit_free(magfit_p[i]);

  gsl_rng_free(rng_p);

  return 0;
}

int
main(int argc, char *argv[])
{
  tiegcm3d_data *data;
  struct timeval tv0, tv1;
  char *infile = NULL;
  char *outfile_map = "data_map.txt";
  char *outfile_alt = "data_alt.txt";
  int time_idx = 0;
  size_t lmax = 60;
  size_t mmax = 30;
  double lon = 150.0; /* desired longitude */
  int r_idx = 0;
  int lon_idx;
  magfield_workspace *magfield_p;
  magfield_eval_workspace *magfield_eval_p;
  magfield_params params;
  size_t i, j, k;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "lmax", required_argument, NULL, 'l' },
          { "mmax", required_argument, NULL, 'm' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "i:l:m:o:r:t:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'i':
            infile = optarg;
            break;

          case 'l':
            lmax = (size_t) atoi(optarg);
            break;

          case 'm':
            mmax = (size_t) atoi(optarg);
            break;

          case 'r':
            r_idx = atol(optarg);
            break;

          case 't':
            time_idx = atol(optarg);
            break;

          case 'o':
            outfile_map = optarg;
            break;

          default:
            break;
        }
    }

  if (!infile)
    {
      fprintf(stderr, "Usage: %s <-i tiegcm3d_nc_file> [-t time_idx] [-r r_idx] [-l lmax] [-m mmax] [-o output_file]\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "input file = %s\n", infile);
  fprintf(stderr, "lmax       = %zu\n", lmax);
  fprintf(stderr, "mmax       = %zu\n", mmax);

  fprintf(stderr, "main: reading %s...", infile);
  gettimeofday(&tv0, NULL);

  data = tiegcm3d_read(infile, NULL);
  if (!data)
    {
      fprintf(stderr, "main: error reading %s\n", infile);
      exit(1);
    }

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu records read, %g seconds)\n", data->nt,
          time_diff(tv0, tv1));

  params.lmin = 1;
  params.lmax = lmax;
  params.mmax = mmax;
  params.nr = data->nr;
  params.ntheta = data->nlat;
  params.nphi = data->nlon;
  params.rmin = data->r[0] * 1.0e3;
  params.rmax = data->r[data->nr - 1] * 1.0e3;
  params.R = R_EARTH_M;
  params.grid_type = MAGFIELD_GAUSS;

  fprintf(stderr, "rmin = %.1f [m]\n", params.rmin);
  fprintf(stderr, "rmax = %.1f [m]\n", params.rmax);

  fprintf(stderr, "\t allocating magfield workspace...");
  magfield_p = magfield_alloc(&params);
  magfield_eval_p = magfield_eval_alloc(&params);
  fprintf(stderr, "done\n");

  magfield_set_r(data->r, magfield_p);
  for (i = 0; i < data->nr; ++i)
    magfield_p->r[i] *= 1.0e3; /* convert to m */

  /* fill current grid */
  fprintf(stderr, "main: filling current grid...");

  for (i = 0; i < data->nr; ++i)
    {
      for (j = 0; j < data->nlat; ++j)
        {
          for (k = 0; k < data->nlon; ++k)
            {
              size_t idx = TIEGCM3D_IDX(0, i, j, k, data);

              magfield_grid_set(MAG_IDX_R, i, j, k, data->Jr[idx], magfield_p);
              magfield_grid_set(MAG_IDX_THETA, i, j, k, data->Jt[idx], magfield_p);
              magfield_grid_set(MAG_IDX_PHI, i, j, k, data->Jp[idx], magfield_p);
            }
        }
    }

  fprintf(stderr, "done\n");

  /* perform SH decomposition */
  fprintf(stderr, "main: performing SH decomposition...");
  magfield_decomp(magfield_p);
  fprintf(stderr, "done\n");

  /* locate index of desired longitude */
  lon_idx = bsearch_double(data->glon, lon, 0, data->nlon - 1);

#if 1
  fprintf(stderr, "main: lon_idx = %d (%.2f [deg])\n", lon_idx, data->glon[lon_idx]);

  magfield_eval_init(magfield_p->qtcoeff, magfield_p->qcoeff, magfield_p->pcoeff, magfield_eval_p);

  fprintf(stderr, "main: writing grid data to %s (time idx = %d, r idx = %d)...", outfile_map, time_idx, r_idx);
  print_data(outfile_map, data, time_idx, r_idx, magfield_eval_p);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: writing altitude data to %s...", outfile_alt);
  print_alt(outfile_alt, data, time_idx, lon_idx, magfield_eval_p);
  fprintf(stderr, "done\n");
#endif

  docalc(magfield_p);

  tiegcm3d_free(data);
  magfield_free(magfield_p);
  magfield_eval_free(magfield_eval_p);

  return 0;
}
