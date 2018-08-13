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

#include <common/common.h>
#include <common/bsearch.h>

#include <magfield/magfield.h>

#include "magfit.h"
#include "tiegcm3d.h"

/*
print_data()
  Print J grid for a fixed time and altitude

Inputs: data - tiegcm data
*/

int
print_data(const char *filename, const tiegcm3d_data *data, const int time_idx, const int ir, magfield_workspace *w)
{
  int s = 0;
  size_t i;
  size_t ilat, ilon;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_data: unable to open %s: %s\n",
              filename, strerror(errno));
    }

  i = 1;
  fprintf(fp, "# Time: %ld (%.6f DOY)\n", data->t[time_idx], data->doy[time_idx] + data->ut[time_idx] / 24.0);
  fprintf(fp, "# Radius: %.2f (km) [%.2f km altitude]\n", data->r[ir], data->r[ir] - 6371.2);
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: J_r (A/m^2)\n", i++);
  fprintf(fp, "# Field %zu: J_t (A/m^2)\n", i++);
  fprintf(fp, "# Field %zu: J_p (A/m^2)\n", i++);
  fprintf(fp, "# Field %zu: B_r (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_t (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_p (nT)\n", i++);
  fprintf(fp, "# Field %zu: |B| (nT)\n", i++);

  for (ilon = 0; ilon < data->nlon; ++ilon)
    {
      double phi = w->phi[ilon];

      for (ilat = 0; ilat < data->nlat; ++ilat)
        {
          size_t idx = TIEGCM3D_IDX(time_idx, ir, ilat, ilon, data);
          double theta = w->theta[ilat];
          double B[4];

          magfield_eval_B(data->r[ir] * 1.0e3, theta, phi, B, w);

          fprintf(fp, "%8.4f %8.4f %16.4e %16.4e %16.4e %16.4e %16.4e %16.4e %16.4e\n",
                  data->glon[ilon],
                  data->glat[ilat],
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

/*
print_alt()
  Print altitude/latitude map of Jr/Jt/Jp for fixed time and longitude

Inputs: data - tiegcm data
*/

int
print_alt(const char *filename, const tiegcm3d_data *data, const int it, const int ilon, magfield_workspace *w)
{
  int s = 0;
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
      double theta = w->theta[ilat];

      for (ir = 0; ir < data->nr; ++ir)
        {
          size_t idx = TIEGCM3D_IDX(it, ir, ilat, ilon, data);
          double B[4];

          magfield_eval_B(w->r[ir], theta, w->phi[ilon], B, w);

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
  const double r_km = R_EARTH_KM + 85.0;
  const double r_m = r_km * 1.0e3;
  const double lon[4] = { 145.69, 149.32, 150.68, 152.95 };
  const double dlat = 8.0 / 110.0 * 10.0; /* 8 km/s / 110 km/deg =~ 0.072 degrees latitude per second */
#if 0
  const double fac = 1.4;          /* factor to account for missing TIEGCM current below 109km */
#endif
  const double fac = 2.7;          /* factor to simulate solar max conditions */
  double phi[4];
  magfit_parameters params = magfit_default_parameters();
  magfit_workspace *magfit_p[4];
  double lat;
  size_t i, k;
  double rnorm, snorm;
  gsl_rng *rng_p = gsl_rng_alloc(gsl_rng_default);
  const double sigma = 60.0;

  for (i = 0; i < 4; ++i)
    phi[i] = lon[i] * M_PI / 180.0;

  params.lat_min = -15.0;
  params.lat_max = 25.0;
  params.lat_spacing1d = 1.0;
  params.flags = MAGFIT_FLG_FIT_X | MAGFIT_FLG_FIT_Z;
  params.flags |= MAGFIT_FLG_SECS_FIT_DF;

  for (i = 0; i < 4; ++i)
    magfit_p[i] = magfit_alloc(magfit_secs1d, &params);

  i = 1;
  fprintf(stdout, "# Field %zu: geocentric latitude (degrees)\n", i++);
  fprintf(stdout, "# Field %zu: |B| data for longitude 145.69 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: |B| data for longitude 149.32 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: |B| data for longitude 150.68 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: |B| data for longitude 152.95 degrees (nT)\n", i++);

  for (k = 0; k < 3; ++k) {
  for (lat = -10.0; lat <= 20.0; lat += dlat)
    {
      double theta = M_PI / 2.0 - lat * M_PI / 180.0;
      double B1[4], B2[4], B3[4], B4[4];
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

#if 1
          B1[i] += gsl_ran_gaussian(rng_p, sigma);
          B2[i] += gsl_ran_gaussian(rng_p, sigma);
          B3[i] += gsl_ran_gaussian(rng_p, sigma);
          B4[i] += gsl_ran_gaussian(rng_p, sigma);
#endif
        }

      B1[3] = gsl_hypot3(B1[0], B1[1], B1[2]);
      B2[3] = gsl_hypot3(B2[0], B2[1], B2[2]);
      B3[3] = gsl_hypot3(B3[0], B3[1], B3[2]);
      B4[3] = gsl_hypot3(B4[0], B4[1], B4[2]);

      /* convert spherical to NEC */

      B1_nec[0] = -B1[1];
      B1_nec[1] = B1[2];
      B1_nec[2] = -B1[0];

      B2_nec[0] = -B2[1];
      B2_nec[1] = B2[2];
      B2_nec[2] = -B2[0];

      B3_nec[0] = -B3[1];
      B3_nec[1] = B3[2];
      B3_nec[2] = -B3[0];

      B4_nec[0] = -B4[1];
      B4_nec[1] = B4[2];
      B4_nec[2] = -B4[0];

      magfit_add_datum(0.0, r_km, theta, phi[0], 0.0, B1_nec, magfit_p[0]);
      magfit_add_datum(0.0, r_km, theta, phi[1], 0.0, B2_nec, magfit_p[1]);
      magfit_add_datum(0.0, r_km, theta, phi[2], 0.0, B3_nec, magfit_p[2]);
      magfit_add_datum(0.0, r_km, theta, phi[3], 0.0, B4_nec, magfit_p[3]);

      if (k == 0)
        {
          printf("%8.4f %16.4f %16.4f %16.4f %16.4f\n",
                 lat,
                 -B1[1],
                 -B2[1],
                 -B3[1],
                 -B4[1]);
        }
    }
    }

  /* fit EEJ model */
  for (i = 0; i < 4; ++i)
    magfit_fit(&rnorm, &snorm, magfit_p[i]);

  i = 1;
  fprintf(stdout, "\n\n");
  fprintf(stdout, "# Field %zu: geocentric latitude (degrees)\n", i++);
  fprintf(stdout, "# Field %zu: |B| fitted for longitude 145.69 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: |B| fitted for longitude 149.32 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: |B| fitted for longitude 150.68 degrees (nT)\n", i++);
  fprintf(stdout, "# Field %zu: |B| fitted for longitude 152.95 degrees (nT)\n", i++);
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

      printf("%8.4f %16.4f %16.4f %16.4f %16.4f %16.4f %16.4f %16.4f %16.4f\n",
             lat,
#if 1
             B1_fit_nec[0],
             B2_fit_nec[0],
             B3_fit_nec[0],
             B4_fit_nec[0],
#else
             gsl_hypot3(B1_fit_nec[0], B1_fit_nec[1], B1_fit_nec[2]),
             gsl_hypot3(B2_fit_nec[0], B2_fit_nec[1], B2_fit_nec[2]),
             gsl_hypot3(B3_fit_nec[0], B3_fit_nec[1], B3_fit_nec[2]),
             gsl_hypot3(B4_fit_nec[0], B4_fit_nec[1], B4_fit_nec[2]),
#endif
             J1_nec[1],
             J2_nec[1],
             J3_nec[1],
             J4_nec[1]);
    }

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
  double lon = 150.0; /* desired longitude */
  int r_idx = 0;
  int lon_idx;
  magfield_workspace *magfield_p;
  magfield_params params;
  size_t i, j, k;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "i:o:r:t:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'i':
            infile = optarg;
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
      fprintf(stderr, "Usage: %s <-i tiegcm3d_nc_file> [-t time_idx] [-o output_file]\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "input file = %s\n", infile);

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
  params.lmax = 60;
  params.mmax = 30;
  params.nr = data->nr;
  params.ntheta = data->nlat;
  params.nphi = data->nlon;
  params.rmin = data->r[0] * 1.0e3;
  params.rmax = data->r[data->nr - 1] * 1.0e3;
  params.R = R_EARTH_M;
  params.grid_type = MAGFIELD_GAUSS;

  fprintf(stderr, "\t allocating magfield workspace...");
  magfield_p = magfield_alloc(&params);
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

              magfield_current_set(MAG_IDX_R, i, j, k, data->Jr[idx], magfield_p);
              magfield_current_set(MAG_IDX_THETA, i, j, k, data->Jt[idx], magfield_p);
              magfield_current_set(MAG_IDX_PHI, i, j, k, data->Jp[idx], magfield_p);
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

#if 0
  fprintf(stderr, "main: lon_idx = %d (%.2f [deg])\n", lon_idx, data->glon[lon_idx]);

  fprintf(stderr, "main: writing grid data to %s (time idx = %d, r idx = %d)...", outfile_map, time_idx, r_idx);
  print_data(outfile_map, data, time_idx, r_idx, magfield_p);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: writing altitude data to %s...", outfile_alt);
  print_alt(outfile_alt, data, time_idx, lon_idx, magfield_p);
  fprintf(stderr, "done\n");
#endif

  docalc(magfield_p);

  tiegcm3d_free(data);
  magfield_free(magfield_p);

  return 0;
}
