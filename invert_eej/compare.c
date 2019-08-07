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

#include <indices/indices.h>
#include <satdata/eef.h>
#include <satdata/satdata.h>
#include <common/common.h>
#include <common/interp.h>
#include <msynth/msynth.h>

#include "jicmag.h"
#include "julia.h"

#include "curvefit.h"
#include "efi.h"

/* allowed difference in measurement times (minutes) */
#define TIME_WINDOW       (5.0)

/* allowed difference in longitude */
#define LON_WINDOW        (10.0)

#define MIN_KP            (8.0)
#define MAX_KP            (9.0)

#define MAX_PAIRS         1000

static int
print_data(const char *filename, const satdata_eef *eef)
{
  size_t i;
  FILE *fp;

  fp = fopen(filename, "w");
  if (!fp)
    {
      fprintf(stderr, "print_data: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  i = 1;
  fprintf(fp, "# Field %zu: time (decimal year)\n", i++);
  fprintf(fp, "# Field %zu: local time (hours)\n", i++);
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: EEF (mV/m)\n", i++);
  fprintf(fp, "# Field %zu: RelErr\n", i++);

  for (i = 0; i < eef->n; ++i)
    {
      time_t t = satdata_epoch2timet(eef->t[i]);
      double phi = eef->longitude[i] * M_PI / 180.0;
      double lt = get_localtime(t, phi);

      fprintf(fp, "%f %f %f %f %f %f\n",
              satdata_epoch2year(eef->t[i]),
              lt,
              eef->longitude[i],
              eef->latitude[i],
              eef->EEF[i] * 1.0e3,
              eef->RelErr[i]);
    }

  fclose(fp);

  return 0;
}

static satdata_eef *
eef_read(const char *filename)
{
  FILE *fp;
  char buffer[2048];
  satdata_eef *data;
  size_t n = 0;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "eef_read: unable to open %s: %s\n",
              filename, strerror(errno));
      return NULL;
    }

  data = satdata_eef_alloc(100000);

  while (fgets(buffer, sizeof(buffer), fp) != 0)
    {
      time_t t;
      double lon, dummy, E, RelErr, J;
      int c;

      if (*buffer == '#')
        continue;

      c = sscanf(buffer, "%ld %lf %lf %lf %lf %lf %lf %lf",
                 &t,
                 &dummy,
                 &lon,
                 &dummy,
                 &dummy,
                 &E,
                 &RelErr,
                 &J);
      if (c < 7)
        continue;

      data->t[n] = satdata_timet2epoch(t);
      data->longitude[n] = lon;
      data->latitude[n] = 0.0;
      data->EEF[n] = E * 1.0e-3;
      data->RelErr[n] = RelErr;
      data->J[n] = J;

      ++n;
    }

  fclose(fp);

  data->n = n;

  return data;
}

jicmag_data *
wamnet_read(const char *filename)
{
  FILE *fp;
  char buf[JICMAG_MAX_BUFFER];
  size_t n;
  jicmag_data *data;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "wamnet_read: fopen: cannot open %s: %s\n",
              filename, strerror(errno));
      return 0;
    }

  data = jicmag_alloc(2000000);

  n = data->n;

  while (fgets(buf, sizeof(buf), fp) != NULL)
    {
      int c;
      double fday, dH;

      c = sscanf(buf, "%lf %lf", &fday, &dH);
      if (c < 2)
        continue;

      if (!isfinite(dH))
        continue;

      data->array[n].t = fday2timet(fday);
      data->array[n].H = dH;

      if (++n >= data->ntot)
        {
          fprintf(stderr, "wamnet_read: error: ntot too small\n");
          exit(1);
        }
    }

  fclose(fp);

  data->n = n;

  return data;
} /* wamnet_read() */

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
} /* compare_mag() */

/* longitude is lon of observatory pair in degrees */
static size_t
find_pairs_mag(const double longitude, const satdata_eef *eef, const jicmag_data *mag,
               double *E_eef, double *dH_mag, double *J)
{
  const double R = 6371.2;
  const double julia_theta = M_PI / 2.0 - JULIA_GEOCENTRIC_LAT * M_PI / 180.0;
  const double julia_phi = JULIA_LON * M_PI / 180.0;
  size_t i;
  size_t npairs = 0;
  msynth_workspace *msynth_p = msynth_igrf_read(MSYNTH_IGRF_FILE);
  f107_workspace *f107_p = f107_alloc(F107_IDX_FILE);
  kp_workspace *kp_p = kp_alloc(KP_IDX_FILE);

  i = 1;
  printf("# Field %zu: timestamp\n", i++);
  printf("# Field %zu: dH (nT)\n", i++);
  printf("# Field %zu: EEF (mV/m)\n", i++);
  printf("# Field %zu: peak current density (mA/m)\n", i++);
  printf("# Field %zu: EEF (Anderson formula) (mV/m)\n", i++);

  for (i = 0; i < eef->n; ++i)
    {
      time_t t = satdata_epoch2timet(eef->t[i]);
      double londiff = wrap180(eef->longitude[i] - longitude);
      jicmag_datum *ptr;
      double kp;

      if (fabs(londiff) > LON_WINDOW)
        continue;

      kp_get(t, &kp, kp_p);
      if ((kp < MIN_KP) || (kp > MAX_KP))
        continue;

      /* discard large outliers */
      if (fabs(eef->J[i]) > 150.0)
        continue;

      ptr = bsearch(&t, mag->array, mag->n, sizeof(jicmag_datum), compare_mag);
      if (ptr)
        {
          size_t idx = ptr - mag->array;
          double dH;
          double v, E, f107;
          double B[4];
          double year = get_year(ptr->t);

          assert(ptr->t == mag->array[idx].t);

#if 0
          dH = mag->array[idx].H;
#else
          dH = interp1d((double) mag->array[idx].t, (double) mag->array[idx + 1].t,
                        mag->array[idx].H, mag->array[idx + 1].H, t);
#endif

          if (!isfinite(dH))
            fprintf(stderr, "here\n");

          if (fabs(dH) > 170.0)
            continue;

          f107_get(ptr->t, &f107, f107_p);
          msynth_eval(year, R + 150.0, julia_theta, julia_phi, B, msynth_p);

          /* compute ExB drift with Anderson's formula */
          v = 12.26 - 0.0454*f107 + 0.1892*dH +
              0.00028*dH*dH - 0.0000022*dH*dH*dH;

          /* compute electric field in mV/m */
          E = v * B[3] * 1.0e-6;

          printf("%ld %f %f %f %f\n",
                 t,
                 dH,
                 eef->EEF[i] * 1.0e3,
                 eef->J[i],
                 E);

          E_eef[npairs] = eef->EEF[i] * 1.0e3;
          dH_mag[npairs] = dH;
          J[npairs] = eef->J[i];

          ++npairs;
        }
    }

  f107_free(f107_p);
  kp_free(kp_p);
  msynth_free(msynth_p);

  return npairs;
}

/*
find_pairs()
  Find data pairs
*/

static size_t
find_pairs(const satdata_eef *eef, const julia_data *julia,
           const satdata_efi *efi,
           double *E_eef, double *E_julia, double *E_efi)
{
  const double R = 6371.2;
  const double julia_theta = M_PI / 2.0 -
                             JULIA_GEOCENTRIC_LAT * M_PI / 180.0;
  size_t i;
  size_t n = 0;
  msynth_workspace *msynth_p = msynth_igrf_read(MSYNTH_IGRF_FILE);
  /*msynth_workspace *msynth_p = msynth_swarm_read(MSYNTH_CHAOS_FILE);*/

  i = 1;
  printf("# Field %zu: timestamp\n", i++);
  printf("# Field %zu: local time (hours)\n", i++);
  printf("# Field %zu: E_julia (mV/m)\n", i++);
  printf("# Field %zu: E_eef (mV/m)\n", i++);
  printf("# Field %zu: E_efi (mV/m)\n", i++);

  for (i = 0; i < eef->n; ++i)
    {
      time_t t = satdata_epoch2timet(eef->t[i]);
      double lt = get_localtime(t, JULIA_LON);
      double londiff = wrap180(eef->longitude[i] - JULIA_LON);
      int s;
      size_t julia_idx, efi_idx;
      double year;
      double B[4];

      if (fabs(londiff) > LON_WINDOW)
        continue;

      s = julia_search(t, TIME_WINDOW, julia, &julia_idx);
      if (s)
        continue; /* no JULIA data found */

#if 0
      if (efi)
        {
          s = efi_search(t, 0.1, efi, &efi_idx);
          if (s)
            continue; /* no EFI data found */
        }
#endif

      year = get_year(julia->array[julia_idx].t);
      msynth_eval(year, R + 150.0, julia_theta, JULIA_LON_RAD, B, msynth_p);

      /* compute E_julia in mV/m */
      E_julia[n] = julia->array[julia_idx].v_vert * B[3] * 1.0e-6;
      E_eef[n] = eef->EEF[i];

#if 0
      if (efi)
        E_efi[n] = efi->array[efi_idx].E_phi;
      else
#endif
        E_efi[n] = 0.0;

      /* JULIA seems to have some anomolous values, try to filter out */
      if ((E_julia[n] < 0.0 && E_eef[n] > 0.15) ||
          E_julia[n] > 1.2)
        continue;

      printf("%ld %f %f %f %f\n",
             t,
             lt,
             E_julia[n],
             E_eef[n],
             E_efi[n]);
      ++n;
    }

  msynth_free(msynth_p);

  return n;
}

static double
compute_rms(const double *x, const double *y, const size_t n)
{
  size_t i;
  double rms = 0.0;

  for (i = 0; i < n; ++i)
    rms += pow(x[i] - y[i], 2.0);

  rms = sqrt(rms / (double)n);
  return rms;
}

int
main(int argc, char *argv[])
{
  char *eef_file = NULL;
  char *swarm_eef_file = NULL;
  char *mag_file = NULL;
  char *out_file = NULL;
  satdata_eef *eef_data;
  satdata_efi *efi_data = NULL;
  struct timeval tv0, tv1;
  int c;
  julia_data *data_julia = NULL;
  jicmag_data *data_wamnet = NULL;

  while ((c = getopt(argc, argv, "e:f:m:jo:s:w:")) != (-1))
    {
      switch (c)
        {
          case 'e':
            eef_file = optarg;
            break;

          case 's':
            swarm_eef_file = optarg;
            break;

          case 'm':
            mag_file = optarg;
            break;

          case 'w':
            fprintf(stderr, "reading WAMNET file %s...", optarg);
            data_wamnet = wamnet_read(optarg);
            fprintf(stderr, "done (%zu data read)\n", data_wamnet->n);

            fprintf(stderr, "sorting WAMNET data...");
            jicmag_sort(data_wamnet);
            fprintf(stderr, "done\n");

            break;

          case 'j':
#if 1
            fprintf(stderr, "main: reading JULIA 5-minute average data...");
            data_julia = julia_read_avg_idx(JULIA_AVG_IDX_FILE, NULL);
#else
            fprintf(stderr, "main: reading JULIA data...");
            data_julia = julia_read_idx(JULIA_IDX_FILE, NULL);
#endif
            fprintf(stderr, "done (%zu data read)\n", data_julia->n);

            fprintf(stderr, "main: sorting JULIA data...");
            julia_sort(data_julia);
            fprintf(stderr, "done\n");

            break;

          case 'f':
            fprintf(stderr, "main: reading EFI data file %s...", optarg);
            efi_data = satdata_swarm_efi_read_idx(optarg, 0);
            fprintf(stderr, "done (%zu data read)\n", efi_data->n);
            break;

          case 'o':
            out_file = optarg;
            break;

          default:
            exit(1);
        }
    }

  if (!eef_file && !swarm_eef_file)
    {
      fprintf(stderr, "Usage: %s <-e eef_ascii_file> <-s swarm_eef_data_file> [-m jicmag_index_file] [-j] [-f] [-o output_file] [-w wamnet_file]\n", argv[0]);
      exit(1);
    }

  if (eef_file)
    {
      fprintf(stderr, "main: reading EEF file %s...", eef_file);
      gettimeofday(&tv0, NULL);
      eef_data = eef_read(eef_file);
      gettimeofday(&tv1, NULL);
      if (eef_data == NULL)
        exit(1);
      fprintf(stderr, "done (%zu records read, %g seconds)\n",
              eef_data->n, time_diff(tv0, tv1));
    }
  else if (swarm_eef_file)
    {
      fprintf(stderr, "main: reading Swarm EEF file %s...", swarm_eef_file);
      gettimeofday(&tv0, NULL);
      eef_data = satdata_eef_read_idx(swarm_eef_file, 0);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%zu records read, %g seconds)\n",
              eef_data->n, time_diff(tv0, tv1));
    }

  if (data_wamnet)
    {
      const double longitude_SAM = -5.77;
      const double longitude_JIC = -76.87;
      double E_eef[MAX_PAIRS], dH[MAX_PAIRS], J[MAX_PAIRS];
#if 1
      const double longitude = longitude_JIC;
#else
      const double longitude = longitude_SAM;
#endif
      size_t npairs = find_pairs_mag(longitude, eef_data, data_wamnet,
                                     E_eef, dH, J);

      fprintf(stderr, "found %zu pairs (longitude = %.2f)\n", npairs, longitude);

      fprintf(stderr, "r(E_eef,dH) = %f\n",
              gsl_stats_correlation(E_eef, 1, dH, 1, npairs));
      fprintf(stderr, "r(dH,J_peak) = %f\n",
              gsl_stats_correlation(dH, 1, J, 1, npairs));

      jicmag_free(data_wamnet);
    }
  else
    {
      double E_julia[MAX_PAIRS], E_eef[MAX_PAIRS], E_efi[MAX_PAIRS];
      size_t npairs;

      npairs = find_pairs(eef_data,
                          data_julia,
                          efi_data,
                          E_eef,
                          E_julia,
                          E_efi);

      fprintf(stderr, "main: found %zu pairs\n", npairs);

      if (npairs > 0)
        {
          const gsl_multifit_robust_type *T = gsl_multifit_robust_bisquare;
          curvefit_workspace *curve_p =
            curvefit_alloc(T, curvefit_poly, npairs, 2);

          if (data_julia)
            {
              double rms;

              fprintf(stderr, "r(E_eef,E_julia) = %f\n",
                      gsl_stats_correlation(E_eef, 1, E_julia, 1, npairs));

              rms = compute_rms(E_eef, E_julia, npairs);
              fprintf(stderr, "rms(E_eef,E_julia) = %f [mV/m]\n", rms);

              curvefit(0, E_julia, E_eef, curve_p);

              fprintf(stderr, "E_eef = %f E_julia + %f\n",
                      gsl_vector_get(curve_p->c, 1),
                      gsl_vector_get(curve_p->c, 0));
            }

          curvefit_free(curve_p);
        }
    }

  if (out_file)
    {
      fprintf(stderr, "main: printing EEF data to %s...", out_file);
      print_data(out_file, eef_data);
      fprintf(stderr, "done\n");
    }

  satdata_eef_free(eef_data);

  return 0;
} /* main() */
