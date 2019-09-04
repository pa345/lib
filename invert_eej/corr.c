/*
 * corr.c
 *
 * Correlate Swarm A/B data as a function of longitude
 *
 * Usage: ./corr [options]
 * -a profile_file   - ascii profile.dat file for 1st satellite
 * -b profile_file   - ascii profile.dat file for 2nd satellite
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics.h>

#include <mainlib/ml_indices.h>
#include <mainlib/ml_eef.h>
#include <mainlib/ml_satdata.h>
#include <mainlib/ml_bin2d.h>
#include <mainlib/ml_common.h>

static satdata_eef *
profile_read(const char *filename)
{
  FILE *fp;
  char buffer[2048];
  satdata_eef *data;
  size_t n = 0;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "profile_read: unable to open %s: %s\n",
              filename, strerror(errno));
      return NULL;
    }

  data = satdata_eef_alloc(100000);

  while (fgets(buffer, sizeof(buffer), fp) != 0)
    {
      size_t ntrack;
      time_t t;
      double lon, F2, lt, seas, kp, J;
      int dir;
      int c;

      if (*buffer == '#')
        continue;

      c = sscanf(buffer, "%zu %ld %lf %lf %lf %lf %lf %lf %d",
                 &ntrack,
                 &t,
                 &lon,
                 &lt,
                 &seas,
                 &F2,
                 &J,
                 &kp,
                 &dir);
      if (c < 9)
        continue;

      data->t[n] = satdata_timet2epoch(t);
      data->longitude[n] = lon;
      data->latitude[n] = 0.0;
      data->EEF[n] = 0.0;
      data->RelErr[n] = 0.0;
      data->J[n] = J;

      ++n;
    }

  fclose(fp);

  data->n = n;

  return data;
}

size_t
eef_search(const double t, const satdata_eef *eef)
{
  size_t i;
  size_t idx = 0;
  double dt_min = 1.0e12;

  for (i = 0; i < eef->n; ++i)
    {
      double dt = t - eef->t[i];

      if (fabs(dt) < dt_min)
        {
          dt_min = fabs(dt);
          idx = i;
        }
    }

  return idx;
}

size_t
docorr(const char *filename, const satdata_eef *eef1, const satdata_eef *eef2)
{
  size_t i, j;
  FILE *fp;
  size_t cnt = 0;
  bin2d_workspace *bin_p;
  const double tmin = -50.0;
  const double tmax = 50.0;
  const double lonmin = -50.0;
  const double lonmax = 20.0;
  const size_t nt = 10;
  const size_t nlon = 5;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "docorr: unable to open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  bin_p = bin2d_alloc(lonmin, lonmax, nlon, tmin, tmax, nt);

  i = 1;
  fprintf(fp, "# Field %zu: timestamp of satellite 1 (UT)\n", i++);
  fprintf(fp, "# Field %zu: delta t of equator crossings (min)\n", i++);
  fprintf(fp, "# Field %zu: delta longitude of equator crossings (degrees)\n", i++);
  fprintf(fp, "# Field %zu: peak J1 (A/m)\n", i++);
  fprintf(fp, "# Field %zu: peak J2 (A/m)\n", i++);

  for (i = 0; i < eef1->n; ++i)
    {
      double t1 = eef1->t[i];
      size_t idx = eef_search(t1, eef2);
      double t2 = eef2->t[idx];
      double dt_s = (t1 - t2) / 1000.0; /* convert to s */
      double dt_m = dt_s / 60.0;        /* convert to m */
      double dlon = wrap180(eef1->longitude[i] - eef2->longitude[idx]);

      if (fabs(dt_s) > 3600.0)
        continue;

      if (dlon < lonmin || dlon > lonmax)
        continue;

      if (dt_m < tmin || dt_m > tmax)
        continue;

      bin2d_add_element_corr(dlon, dt_m, eef1->J[i], eef2->J[idx], bin_p);

      fprintf(fp, "%ld %f %f %f %f\n",
              satdata_epoch2timet(t1),
              dt_m,
              dlon,
              eef1->J[i],
              eef2->J[idx]);
      ++cnt;
    }

  fclose(fp);

  i = 1;
  printf("# Field %zu: delta longitude (degrees)\n", i++);
  printf("# Field %zu: delta t (minutes)\n", i++);
  printf("# Field %zu: correlation\n", i++);
  printf("# Field %zu: number of points in bin\n", i++);

  for (i = 0; i < nlon; ++i)
    {
      for (j = 0; j < nt; ++j)
        {
          double dt, dlon;
          size_t n;
          double r;

          bin2d_xyval(i, j, &dlon, &dt, bin_p);
          n = bin2d_n(dlon, dt, bin_p);
          r = bin2d_correlation(dlon, dt, bin_p);

          printf("%f %f %f %zu\n",
                 dlon,
                 dt,
                 r,
                 n);
        }

      printf("\n");
    }

  bin2d_free(bin_p);

  return cnt;
}

int
main(int argc, char *argv[])
{
  char *outfile = "corr.dat";
  satdata_eef *eef1 = NULL;
  satdata_eef *eef2 = NULL;
  int c;

  while ((c = getopt(argc, argv, "a:b:")) != (-1))
    {
      switch (c)
        {
          case 'a':
            fprintf(stderr, "main: reading profile data from %s...", optarg);
            eef1 = profile_read(optarg);
            fprintf(stderr, "done (%zu data read)\n", eef1->n);
            break;

          case 'b':
            fprintf(stderr, "main: reading profile data from %s...", optarg);
            eef2 = profile_read(optarg);
            fprintf(stderr, "done (%zu data read)\n", eef2->n);
            break;

          default:
            exit(1);
        }
    }

  if (!eef1 || !eef2)
    {
      fprintf(stderr, "Usage: %s <-a profile_ascii_file> <-b profile_ascii_file>\n", argv[0]);
      exit(1);
    }

  {
    size_t ncorr;

    fprintf(stderr, "main: writing correlation data to %s...", outfile);
    ncorr = docorr(outfile, eef1, eef2);
    fprintf(stderr, "done (%zu equator crossing pairs found)\n", ncorr);
  }

  satdata_eef_free(eef1);
  satdata_eef_free(eef2);

  return 0;
} /* main() */
