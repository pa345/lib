/*
 * compare.c
 *
 * Compare Swarm EEF data with other observations (ie: JULIA)
 *
 * Usage: ./compare [options]
 * -e eef_file       - input file from invert_eej module
 * -s swarm_eef_file - input file in official Swarm EEF format
 * -j                - compare input with JULIA data
 *  -w wamnet_file   - compare with WAMNET dH values
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
#include <mainlib/ml_common.h>
#include <mainlib/ml_interp.h>

#include "grobs.h"

/* allowed difference in measurement times (minutes) */
#define TIME_WINDOW       (5.0)

/* allowed difference in longitude */
#define LON_WINDOW        (5.0)

#define MAX_KP            (2.0)

#define MAX_PAIRS         1000

#define MAX_DATA          8192

typedef struct
{
  time_t t[MAX_DATA];
  double lt[MAX_DATA];
  double longitude[MAX_DATA];
  double J[MAX_DATA];
  size_t n;
} current_data;

static current_data *
current_read(const char *filename)
{
  FILE *fp;
  char buffer[2048];
  current_data *data;
  size_t n = 0;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "current_read: unable to open %s: %s\n",
              filename, strerror(errno));
      return NULL;
    }

  data = calloc(1, sizeof(current_data));

  while (fgets(buffer, sizeof(buffer), fp) != 0)
    {
      time_t t;
      double lon, lt, J;
      int c;

      if (*buffer == '#')
        continue;

      c = sscanf(buffer, "%ld %lf %lf %lf",
                 &t,
                 &lt,
                 &lon,
                 &J);
      if (c < 4)
        continue;

      data->t[n] = t;
      data->lt[n] = lt;
      data->longitude[n] = lon;
      data->J[n] = J;

      ++n;
    }

  fclose(fp);

  data->n = n;

  return data;
}

#if 0
static int
compare_mag(const void *a, const void *b)
{
  const time_t t1 = *(const time_t *) a;
  jicmag_datum *db = (jicmag_datum *) b;
  const time_t t2 = db->t;
  double diff = (t2 - t1) / 60.0;

  if (fabs(diff) < TIME_WINDOW)
    return 0; /* found */
  else if (t1 < t2)
    return -1;
  else
    return 1;
}
#endif

/* longitude is lon of observatory pair in degrees */
static size_t
find_pairs_mag(const double longitude, const current_data *cdata, const grobs_data *mag,
               double *dH_mag, double *J)
{
  size_t i, j;
  size_t npairs = 0;
  kp_workspace *kp_p = kp_alloc(KP_IDX_FILE);

  i = 1;
  printf("# Field %zu: timestamp\n", i++);
  printf("# Field %zu: dH (nT)\n", i++);
  printf("# Field %zu: peak current density (A/km)\n", i++);

  for (i = 0; i < cdata->n; ++i)
    {
      time_t t = cdata->t[i];
      double dt, dH;
      double londiff = wrap180(cdata->longitude[i] - longitude);
      double kp;

      if (fabs(londiff) > LON_WINDOW)
        continue;

      /* reject bad current estimates */
      if (fabs(cdata->J[i]) > 5000.0)
        continue;

      kp_get(t, &kp, kp_p);
      if (kp > MAX_KP)
        continue;

      /* search for observatory timestamp corresponding to SECS measurement */
      j = bsearch_timet(mag->t, t, 0, mag->n - 1);
      dt = (t - mag->t[j]) / 60.0; /* convert to minutes */

      if (fabs(dt) > TIME_WINDOW)
        continue;

      /* interpolate observatory dH to time of SECS measurement */
      dH = interp1d((double) mag->t[j], (double) mag->t[j + 1],
                    mag->H[j], mag->H[j + 1], (double) t);

      printf("%ld %f %f\n",
             t,
             dH,
             cdata->J[i]);

      dH_mag[npairs] = dH;
      J[npairs] = cdata->J[i];

      ++npairs;
    }

  kp_free(kp_p);

  return npairs;
}

int
main(int argc, char *argv[])
{
  int c;
  grobs_data *data_wamnet = NULL;
  current_data *data_current = NULL;

  while ((c = getopt(argc, argv, "i:w:")) != (-1))
    {
      switch (c)
        {
          case 'i':
            fprintf(stderr, "main: reading SECS file %s...", optarg);
            data_current = current_read(optarg);
            fprintf(stderr, "done (%zu data read)\n", data_current->n);
            break;

          case 'w':
            fprintf(stderr, "main: reading WAMNET file %s...", optarg);
            data_wamnet = grobs_wamnet_read(optarg, NULL);
            fprintf(stderr, "done (%zu data read)\n", data_wamnet->n);

            break;

          default:
            exit(1);
        }
    }

  if (!data_wamnet)
    {
      fprintf(stderr, "Usage: %s <-i current_file> [-o output_file] [-w wamnet_file]\n", argv[0]);
      exit(1);
    }

  if (data_wamnet)
    {
      double dH[MAX_PAIRS], J[MAX_PAIRS];
      size_t npairs = find_pairs_mag(data_wamnet->glon, data_current, data_wamnet,
                                     dH, J);

      fprintf(stderr, "found %zu pairs\n", npairs);

      fprintf(stderr, "r(dH,J_peak) = %f\n",
              gsl_stats_correlation(dH, 1, J, 1, npairs));

      grobs_free(data_wamnet);
    }

  return 0;
}
