/*
 * print_fillion.c
 *
 * 1) Read Swarm file
 * 2) Read 3D TIEGCM grid file
 * 3) Compute magnetic field from 3D currents
 * 4) Output magnetic field values along Swarm orbit
 *
 * ./print_fillion <-i tiegcm_nc_file> <-s swarm_file>
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

#include <mainlib/ml_satdata.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_bsearch.h>
#include <mainlib/ml_interp.h>
#include <mainlib/ml_track.h>

#include <magfield/magfield.h>

#include "tiegcm3d.h"

#define MAX_PTS            1000

satdata_mag *
read_swarm(const char *filename)
{
  size_t nflag;
  satdata_mag *data;
  struct timeval tv0, tv1;

  fprintf(stderr, "read_swarm: reading %s...", filename);
  gettimeofday(&tv0, NULL);

  data = satdata_swarm_read_idx(filename, 0);

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%zu data read, %g seconds)\n",
          data->n, time_diff(tv0, tv1));

  /* check for instrument flags since we use Stage1 data */

  fprintf(stderr, "read_swarm: filtering for instrument flags...");
  nflag = satdata_swarm_filter_instrument(1, data);
  fprintf(stderr, "done (%zu/%zu (%.1f%%) data flagged)\n",
          nflag, data->n, (double)nflag / (double)data->n * 100.0);

  return data;
}

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

          /*magfield_eval_B(data->r[ir] * 1.0e3, theta, phi, B, w);*/
          magfield_eval_B((R_EARTH_KM + 450.0) * 1.0e3, theta, phi, B, w);

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
fill_grid(const size_t time_idx, tiegcm3d_data *data, magfield_workspace *magfield_p)
{
  size_t i, j, k;

  for (i = 0; i < data->nr; ++i)
    {
      for (j = 0; j < data->nlat; ++j)
        {
          for (k = 0; k < data->nlon; ++k)
            {
              size_t idx = TIEGCM3D_IDX(time_idx, i, j, k, data);

              magfield_grid_set(MAG_IDX_R, i, j, k, data->Jr[idx], magfield_p);
              magfield_grid_set(MAG_IDX_THETA, i, j, k, data->Jt[idx], magfield_p);
              magfield_grid_set(MAG_IDX_PHI, i, j, k, data->Jp[idx], magfield_p);
            }
        }
    }

  return GSL_SUCCESS;
}

int
print_track(FILE *fp, const size_t tidx, tiegcm3d_data * data, track_data * tptr,
            satdata_mag * data_sat, magfield_workspace * magfield_p,
            magfield_workspace * magfield_p2)
{
  size_t i, j;

  for (i = 0; i < tptr->n; ++i)
    {
      size_t didx = tptr->start_idx + i;
      time_t t = epoch2timet(data_sat->t[didx]);
      double fday = time2fday(t);
      double r = data_sat->r[didx];
      double theta = M_PI / 2.0 - data_sat->latitude[didx] * M_PI / 180.0;
      double phi = data_sat->longitude[didx] * M_PI / 180.0;
      struct tm *tmp = gmtime(&t);
      double B1[4], B2[4], B[4];

      tmp->tm_year = 2009 - 1900;
      t = mktime(tmp);

      if (t < data->t[tidx] || t > data->t[tidx+1])
        continue;

      magfield_eval_B(r * 1.0e3, theta, phi, B1, magfield_p);
      magfield_eval_B(r * 1.0e3, theta, phi, B2, magfield_p2);

      for (j = 0; j < 3; ++j)
        B[j] = interp1d((double) data->t[tidx], (double) data->t[tidx+1], B1[j], B2[j], (double) t);

      fprintf(fp, "%f %.4f %.4f %.4f %.4f %.4f %.4f\n",
              fday,
              r,
              data_sat->latitude[didx],
              data_sat->longitude[didx],
              B[0] * 1.0e9,
              B[1] * 1.0e9,
              B[2] * 1.0e9);
    }

  fprintf(fp, "\n\n");

  return GSL_SUCCESS;
}

int
main(int argc, char *argv[])
{
  tiegcm3d_data *data;
  satdata_mag *data_sat = NULL;
  track_workspace *track_p;
  struct timeval tv0, tv1;
  char *infile = NULL;
  char *outfile_map = "data_map.txt";
  char *outfile_alt = "data_alt.txt";
  char *outfileB = "satelliteB.txt";
  double lon = 150.0; /* desired longitude */
  time_t t;
  struct tm *tmp;
  size_t tidx = 0;
  int r_idx = 0;
  int lon_idx;
  magfield_workspace *magfield_p, *magfield_p2;
  magfield_params params;
  size_t i, j, k;
  FILE *fp;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "i:o:r:s:", long_options, &option_index);
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

          case 's':
            data_sat = read_swarm(optarg);
            break;

          case 'o':
            outfile_map = optarg;
            break;

          default:
            break;
        }
    }

  if (!infile || !data_sat)
    {
      fprintf(stderr, "Usage: %s <-i tiegcm3d_nc_file> <-s swarm_idx_file> [-o output_file]\n", argv[0]);
      exit(1);
    }

  fprintf(stderr, "input file = %s\n", infile);

  track_p = track_alloc();

  fprintf(stderr, "main: initializing tracks...");
  track_init(data_sat, NULL, track_p);
  fprintf(stderr, "done\n");

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
  magfield_p2 = magfield_alloc(&params);
  fprintf(stderr, "done\n");

  magfield_set_r(data->r, magfield_p);
  magfield_set_r(data->r, magfield_p2);
  for (i = 0; i < data->nr; ++i)
    {
      magfield_p->r[i] *= 1.0e3; /* convert to m */
      magfield_p2->r[i] *= 1.0e3; /* convert to m */
    }

  fp = fopen(outfileB, "w");

  i = 1;
  fprintf(fp, "# Field %zu: timestamp (MJD2000)\n", i++);
  fprintf(fp, "# Field %zu: geocentric radius (km)\n", i++);
  fprintf(fp, "# Field %zu: geocentric latitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: geocentric longitude (deg)\n", i++);
  fprintf(fp, "# Field %zu: B_r (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_t (nT)\n", i++);
  fprintf(fp, "# Field %zu: B_p (nT)\n", i++);

  t = epoch2timet(data_sat->t[0]);
  tmp = gmtime(&t);

  /* TIEGCM run is for 2009, change year and search for month/day in tiegcm3d data */
  tmp->tm_year = 2009 - 1900;
  t = mktime(tmp);

  tidx = bsearch_timet(data->t, t, 0, data->nt - 1);

  fprintf(stderr, "main: filling current grid for timestamp %zu...", tidx);
  fill_grid(tidx, data, magfield_p);
  fill_grid(tidx+1, data, magfield_p2);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: performing SH decomposition for timestamp %zu...", tidx);
  magfield_decomp(magfield_p);
  magfield_decomp(magfield_p2);
  fprintf(stderr, "done\n");

  for (i = 0; i < track_p->n; ++i)
    {
      track_data * tptr = &(track_p->tracks[i]);
      char buf[1024];

      sprintf(buf, "%s", ctime(&data->t[tidx]));
      buf[strlen(buf) - 1] = '\0';
      fprintf(stderr, "=== TRACK %zu: %s, LT = %.2f ===\n", i, buf, tptr->lt_eq);

      for (j = 0; j < tptr->n; ++j)
        {
          size_t didx = tptr->start_idx + j;
          double fday = time2fday(epoch2timet(data_sat->t[didx]));
          double r = data_sat->r[didx];
          double theta = M_PI / 2.0 - data_sat->latitude[didx] * M_PI / 180.0;
          double phi = data_sat->longitude[didx] * M_PI / 180.0;
          double B1[4], B2[4], B[4];

          t = epoch2timet(data_sat->t[didx]);
          tmp = gmtime(&t);
          tmp->tm_year = 2009 - 1900;
          t = mktime(tmp);

          if (t > data->t[tidx + 1])
            {
              tidx = bsearch_timet(data->t, t, 0, data->nt - 1);

              fprintf(stderr, "main: filling current grid for timestamp %zu...", tidx);
              fill_grid(tidx, data, magfield_p);
              fill_grid(tidx+1, data, magfield_p2);
              fprintf(stderr, "done\n");

              fprintf(stderr, "main: performing SH decomposition for timestamp %zu...", tidx);
              magfield_decomp(magfield_p);
              magfield_decomp(magfield_p2);
              fprintf(stderr, "done\n");
            }

          magfield_eval_B(r * 1.0e3, theta, phi, B1, magfield_p);
          magfield_eval_B(r * 1.0e3, theta, phi, B2, magfield_p2);

          for (k = 0; k < 3; ++k)
            B[k] = interp1d((double) data->t[tidx], (double) data->t[tidx+1], B1[k], B2[k], (double) t);

          fprintf(fp, "%f %.4f %.4f %.4f %.4f %.4f %.4f\n",
                  fday,
                  r,
                  data_sat->latitude[didx],
                  data_sat->longitude[didx],
                  B[0] * 1.0e9,
                  B[1] * 1.0e9,
                  B[2] * 1.0e9);
        }

      fprintf(fp, "\n\n");
    }

  fclose(fp);

#if 0
  /* locate index of desired longitude */
  lon_idx = bsearch_double(data->glon, lon, 0, data->nlon - 1);

  fprintf(stderr, "main: lon_idx = %d (%.2f [deg])\n", lon_idx, data->glon[lon_idx]);

  fprintf(stderr, "main: writing grid data to %s (time idx = %d, r idx = %d)...", outfile_map, tidx, r_idx);
  print_data(outfile_map, data, tidx, r_idx, magfield_p);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: writing altitude data to %s...", outfile_alt);
  print_alt(outfile_alt, data, tidx, lon_idx, magfield_p);
  fprintf(stderr, "done\n");
#endif

  tiegcm3d_free(data);
  magfield_free(magfield_p);
  track_free(track_p);

  return 0;
}
