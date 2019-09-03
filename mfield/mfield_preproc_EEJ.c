/*
 * mfield_preproc_EEJ.c
 *
 * Pre-process EEJ magnetic equator data and save in magdata format.
 *
 * Pre-processing steps are:
 *
 * 1. Apply impulse rejection filter to lat/lon EEJ magnetic equator positions to remove outliers
 * 2. Divide lat/lon data into time and longitude bins, compute the median of each bin to represent
 *    best knowledge of magnetic equator for that time period and longitude bin
 *
 * The result is an output file in magdata format containing all data points to
 * be used in the modeling. All data points will have a MAGDATA_FLG_FIT_xxx flag
 * set, and other flags will vary.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>
#include <complex.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <omp.h>
#include <libconfig.h>

#include <apex/apex.h>
#include <indices/indices.h>
#include <satdata/satdata.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_filter.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rstat.h>

#include <common/common.h>
#include <common/bsearch.h>

#include "magdata.h"

#define MAX_N 80000

#define EEJ_FLG_KP          (1 << 0)
#define EEJ_FLG_RC          (1 << 1)
#define EEJ_FLG_OUTLIER     (1 << 2)

typedef struct
{
  time_t t[MAX_N];         /* timestamp of magnetic equator measurement */
  double longitude[MAX_N]; /* longitude of magnetic equator (degrees) */
  double latitude[MAX_N];  /* geocentric latitude of magnetic equator (degrees) */
  double qdlat[MAX_N];     /* IGRF QD latitude of magnetic equator (degrees) */
  double lt[MAX_N];        /* local time (hours) */
  double kp[MAX_N];        /* kp index of magnetic equator measurement */
  double dRC[MAX_N];       /* dRC/dt of magnetic equator measurement (nT/hour) */
  size_t flags[MAX_N];     /* flags EEJ_FLG_xxx */
  double r;                /* radius of magnetic equator (km, same for all data) */
  size_t n;
} eej_data;

typedef struct
{
  double max_kp;          /* maximum kp */
  double max_dRC;         /* maximum dRC/dt (nT/hour) */
} preprocess_parameters;

/* sort EEJ dataset into ascending time order */
static int
eej_sort(eej_data * data)
{
  int s = 0;
  const size_t n = data->n;
  gsl_vector_long_view t = gsl_vector_long_view_array(data->t, n);
  gsl_vector_view lon = gsl_vector_view_array(data->longitude, n);
  gsl_vector_view lat = gsl_vector_view_array(data->latitude, n);
  gsl_vector_view qdlat = gsl_vector_view_array(data->qdlat, n);
  gsl_vector_view lt = gsl_vector_view_array(data->lt, n);
  gsl_vector_view kp = gsl_vector_view_array(data->kp, n);
  gsl_vector_view dRC = gsl_vector_view_array(data->dRC, n);
  gsl_permutation * p = gsl_permutation_alloc(n);

  gsl_sort_vector_long_index(p, &t.vector);

  gsl_permute_vector_long(p, &t.vector);
  gsl_permute_vector(p, &lon.vector);
  gsl_permute_vector(p, &lat.vector);
  gsl_permute_vector(p, &qdlat.vector);
  gsl_permute_vector(p, &lt.vector);
  gsl_permute_vector(p, &kp.vector);
  gsl_permute_vector(p, &dRC.vector);

  gsl_permutation_free(p);

  return s;
}

static int
eej_read(const char * filename, eej_data * data)
{
  int s = 0;
  FILE *fp;
  char buf[2048];
  size_t n = 0;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "eej_read: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  while (fgets(buf, 2048, fp) != NULL)
    {
      int c;
      size_t pnum;
      time_t t;
      double lt, r, lon, lat, lat_model, qdlat, F2, dlat, kp, dRC;
      int icode;

      if (*buf == '#')
        continue;

      c = sscanf(buf, "%zu %ld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",
                 &pnum,
                 &t,
                 &lt,
                 &r,
                 &lon,
                 &lat,
                 &lat_model,
                 &qdlat,
                 &F2,
                 &dlat,
                 &kp,
                 &dRC,
                 &icode);
      if (c < 13)
        continue;

      if (icode != 0)
        continue;

      data->t[n] = t;
      data->longitude[n] = lon;
      data->latitude[n] = lat;
      data->qdlat[n] = qdlat;
      data->lt[n] = lt;
      data->kp[n] = kp;
      data->dRC[n] = dRC;
      data->flags[n] = 0;

      if (++n >= MAX_N)
        {
          fprintf(stderr, "eej_read: MAX_N too small\n");
          break;
        }
    }

  data->r = R_EARTH_KM + 110.0;
  data->n = n;

  fclose(fp);

  /* make sure timestamps are sorted ascending (in case of combined A/B datasets) */
  s = eej_sort(data);

  return s;
}

/* read a EEJ peak data file which was smoothed with a B-spline (fiteq) */
static int
eej_read_smooth(const char * filename, eej_data * data)
{
  int s = 0;
  FILE *fp;
  char buf[2048];
  size_t n = 0;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "eej_read_smooth: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  while (fgets(buf, 2048, fp) != NULL)
    {
      int c;
      time_t t;
      double r, lon, lat, lat_model;

      if (*buf == '#')
        continue;

      c = sscanf(buf, "%ld %lf %lf %lf %lf",
                 &t,
                 &r,
                 &lon,
                 &lat,
                 &lat_model);
      if (c < 5)
        continue;

      data->t[n] = t;
      data->longitude[n] = lon;
      data->latitude[n] = lat;
      data->qdlat[n] = 0.0;
      data->lt[n] = 12.0;
      data->kp[n] = 0.0;
      data->dRC[n] = 0.0;
      data->flags[n] = 0;

      if (++n >= MAX_N)
        {
          fprintf(stderr, "eej_read_smooth: MAX_N too small\n");
          break;
        }
    }

  data->r = R_EARTH_KM + 110.0;
  data->n = n;

  fclose(fp);

  /* make sure timestamps are sorted ascending (in case of combined A/B datasets) */
  s = eej_sort(data);

  return s;
}

static int
eej_read_wavelet(const char * filename, eej_data * data)
{
  int s = 0;
  FILE *fp;
  char buf[2048];
  size_t n = 0;

  fp = fopen(filename, "r");
  if (!fp)
    {
      fprintf(stderr, "eej_read_wavelet: unable to open %s: %s\n",
              filename, strerror(errno));
      return -1;
    }

  while (fgets(buf, 2048, fp) != NULL)
    {
      int c;
      size_t pnum;
      double t_epoch, lon, lat, lat_model;

      if (*buf == '#')
        continue;

      c = sscanf(buf, "%zu %lf %lf %lf %lf",
                 &pnum,
                 &t_epoch,
                 &lon,
                 &lat,
                 &lat_model);
      if (c < 5)
        continue;

      data->t[n] = epoch2timet(t_epoch);
      data->longitude[n] = lon;
      data->latitude[n] = lat;
      data->qdlat[n] = 0.0;
      data->lt[n] = 12.0;
      data->kp[n] = 0.3;
      data->dRC[n] = 1.0;
      data->flags[n] = 0;

      if (++n >= MAX_N)
        {
          fprintf(stderr, "eej_read_wavelet: MAX_N too small\n");
          break;
        }
    }

  data->r = R_EARTH_KM + 110.0;
  data->n = n;

  fclose(fp);

  /* make sure timestamps are sorted ascending (in case of combined A/B datasets) */
  s = eej_sort(data);

  return s;
}

static int
check_parameters(preprocess_parameters * params)
{
  int s = 0;

  if (params->max_kp <= 0.0)
    {
      fprintf(stderr, "check_parameters: max_kp must be > 0\n");
      ++s;
    }

  if (params->max_dRC <= 0.0)
    {
      fprintf(stderr, "check_parameters: max_dRC must be > 0\n");
      ++s;
    }

  return s;
}

/* convert EEJ magnetic equator data to magdata format */
magdata *
copy_data(const size_t magdata_flags, const eej_data * data)
{
  magdata *mdata;
  magdata_datum datum;
  size_t ndata = 0;
  size_t i;

  for (i = 0; i < data->n; ++i)
    {
      if (data->flags[i])
        continue;

      ++ndata;
    }

  if (ndata == 0)
    return NULL;

  mdata = magdata_alloc(ndata, R_EARTH_KM);
  if (!mdata)
    return 0;

  mdata->global_flags = magdata_flags;

  magdata_datum_init(&datum);

  for (i = 0; i < data->n; ++i)
    {
      int s;
      size_t flags = 0;

      if (data->flags[i])
        continue;

      flags |= MAGDATA_FLG_Z;

      datum.t = satdata_timet2epoch(data->t[i]);
      datum.r = data->r;
      datum.theta = M_PI / 2.0 - data->latitude[i] * M_PI / 180.0;
      datum.phi = data->longitude[i] * M_PI / 180.0;
      datum.flags = flags;
      datum.lt = data->lt[i];
      datum.qdlat = data->qdlat[i];
      datum.B_nec[2] = 0.0; /* B_z = 0 */

      s = magdata_add(&datum, mdata);
      if (s)
        break;
    }

  mdata->global_flags |= MAGDATA_GLOBFLG_SATELLITE | MAGDATA_GLOBFLG_EEJ_MAGEQ;

  /* flag all EEJ data for main field modeling */
  for (i = 0; i < mdata->n; ++i)
    mdata->flags[i] |= MAGDATA_FLG_FIT_MF;

  return mdata;
}

int
copy_eej_data(const eej_data * data_in, eej_data * data_out)
{
  int s = 0;
  size_t i;

  for (i = 0; i < data_in->n; ++i)
    {
      data_out->t[i] = data_in->t[i];
      data_out->longitude[i] = data_in->longitude[i];
      data_out->latitude[i] = data_in->latitude[i];
      data_out->qdlat[i] = data_in->qdlat[i];
      data_out->lt[i] = data_in->lt[i];
      data_out->kp[i] = data_in->kp[i];
      data_out->dRC[i] = data_in->dRC[i];
      data_out->flags[i] = data_in->flags[i];
    }

  data_out->r = data_in->r;
  data_out->n = data_in->n;

  return s;
}

/*
preprocess_filter()
  Remove outliers from EEJ magnetic equator data with
impulse rejection filter

Inputs: data - EEJ data

Return: number of outliers detected
*/

size_t
preprocess_filter(const char * filename, eej_data * data)
{
  const size_t n = data->n;
  const size_t K = 51;      /* window size */
  const double t = 3.5;
  FILE *fp = fopen(filename, "w");
  gsl_vector * x = gsl_vector_alloc(n);
  gsl_vector * y = gsl_vector_alloc(n);
  gsl_vector * xmedian = gsl_vector_alloc(n);
  gsl_vector * xsigma = gsl_vector_alloc(n);
  gsl_vector_int * ioutlier = gsl_vector_int_alloc(n);
  gsl_filter_impulse_workspace * w = gsl_filter_impulse_alloc(K);
  gsl_permutation * p = gsl_permutation_alloc(n);
  gsl_vector_const_view lon = gsl_vector_const_view_array(data->longitude, n);
  gsl_vector_const_view lat = gsl_vector_const_view_array(data->latitude, n);
  size_t noutlier = 0;
  size_t i;

  /* compute permutation which sorts longitude into ascending order */
  gsl_sort_vector_index(p, &lon.vector);

  gsl_vector_memcpy(x, &lat.vector);
  gsl_permute_vector(p, x);

  /* apply impulse rejection filter */
  gsl_filter_impulse(GSL_FILTER_END_TRUNCATE, GSL_FILTER_SCALE_QN, t, x, y, xmedian, xsigma, &noutlier,
                     ioutlier, w);

  i = 1;
  fprintf(fp, "# Field %zu: longitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: geocentric latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: median geocentric latitude (degrees)\n", i++);
  fprintf(fp, "# Field %zu: lower error curve (degrees)\n", i++);
  fprintf(fp, "# Field %zu: upper error curve (degrees)\n", i++);
  fprintf(fp, "# Field %zu: outlier detected (0 or 1)\n", i++);

  for (i = 0; i < n; ++i)
    {
      size_t pi = gsl_permutation_get(p, i);
      double lon = data->longitude[pi];
      double lat = data->latitude[pi];
      double xmedi = gsl_vector_get(xmedian, i);
      double xsigmai = gsl_vector_get(xsigma, i);
      int outlier = gsl_vector_int_get(ioutlier, i);

      if (outlier)
        data->flags[pi] |= EEJ_FLG_OUTLIER;

      fprintf(fp, "%10.4f %10.4f %10.4f %10.4f %10.4f %d\n",
              lon,
              lat,
              xmedi,
              xmedi - t * xsigmai,
              xmedi + t * xsigmai,
              outlier);
    }

  fclose(fp);
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(xmedian);
  gsl_vector_free(xsigma);
  gsl_vector_int_free(ioutlier);
  gsl_permutation_free(p);
  gsl_filter_impulse_free(w);

  return noutlier;
}

/*
preprocess_vobs()
*/

int
preprocess_vobs(eej_data * data_in, eej_data * data_out)
{
  int s = 0;
  const size_t n = data_in->n;
  const double dlon = 2.0;
  const size_t nlon = (size_t) round(360.0 / dlon);
  const time_t t0 = data_in->t[0];
  const time_t t1 = data_in->t[n - 1];
  const time_t dt = 120 * 86400; /* time bin size in s */
  gsl_rstat_workspace ** rstat_p = malloc(nlon * sizeof(gsl_rstat_workspace));
  size_t nout = 0;
  time_t t;
  size_t i;

  for (i = 0; i < nlon; ++i)
    rstat_p[i] = gsl_rstat_alloc();

  for (t = t0; t <= t1; t += dt)
    {
      size_t idx1, idx2;

      idx1 = bsearch_timet(data_in->t, t, 0, n - 1);
      idx2 = bsearch_timet(data_in->t, t + dt, 0, n - 1);
      assert(data_in->t[idx1] <= t && t < data_in->t[idx2]);

      for (i = 0; i < nlon; ++i)
        gsl_rstat_reset(rstat_p[i]);

      for (i = idx1; i <= idx2; ++i)
        {
          size_t lon_bin;

          if (data_in->flags[i])
            continue;
          
          lon_bin = (size_t) (0.5 * (nlon - 1.0) * (data_in->longitude[i] / 180.0 + 1.0));
          assert(lon_bin < nlon);

          gsl_rstat_add(data_in->latitude[i], rstat_p[lon_bin]);

#if 0
          printf("%ld %f %f\n",
                 t,
                 data_in->longitude[i],
                 data_in->latitude[i]);
#endif
        }

      for (i = 0; i < nlon; ++i)
        {
          double lon = -180.0 + (i + 0.5) * dlon;

          if (gsl_rstat_n(rstat_p[i]) < 8)
            continue;

#if 0
          printf("%ld %f %f %zu\n",
                 t,
                 lon,
                 gsl_rstat_median(rstat_p[i]),
                 gsl_rstat_n(rstat_p[i]));
#endif

          data_out->t[nout] = t + dt / 2;
          data_out->latitude[nout] = gsl_rstat_median(rstat_p[i]);
          data_out->longitude[nout] = lon;
          data_out->lt[nout] = 0.0;
          data_out->qdlat[nout] = 0.0;
          data_out->kp[nout] = 0.0;
          data_out->dRC[nout] = 0.0;
          data_out->flags[nout] = 0;

          ++nout;
        }

#if 0
      printf("\n\n");
#endif
    }

  data_out->n = nout;
  data_out->r = data_in->r;

  for (i = 0; i < nlon; ++i)
    gsl_rstat_free(rstat_p[i]);

  free(rstat_p);

  return s;
}

/*
preprocess_data()
  Filter EEJ data for kp, |dRC/dt| and outliers

Inputs: params   - preprocess parameters
        data_in  - input EEJ dataset
        data_out - output EEJ dataset

Return: success/error
*/

int
preprocess_data(const preprocess_parameters *params, eej_data * data_in, eej_data * data_out)
{
  const char * outlier_file = "outlier.txt";
  size_t nkp = 0;
  size_t nRC = 0;
  size_t noutlier;
  size_t i;

  fprintf(stderr, "preprocess_data: filtering outliers...");
  noutlier = preprocess_filter(outlier_file, data_in);
  fprintf(stderr, "done (%zu outliers detected, wrote %s)\n", noutlier, outlier_file);

  fprintf(stderr, "preprocess_data: max_kp       = %g\n", params->max_kp);
  fprintf(stderr, "preprocess_data: max |dRC/dt| = %g\n", params->max_dRC);

#if 0
  for (i = 0; i < data_in->n; ++i)
    {
      if (data_in->kp[i] > params->max_kp)
        {
          ++nkp;
          data_in->flags[i] |= EEJ_FLG_KP;
        }

      if (fabs(data_in->dRC[i]) > params->max_dRC)
        {
          ++nRC;
          data_in->flags[i] |= EEJ_FLG_RC;
        }
    }
#endif

  fprintf(stderr, "\t\t flagged data due to outlier:   %zu (%.1f%%)\n", noutlier, (double) noutlier / (double) data_in->n * 100.0);
  fprintf(stderr, "\t\t flagged data due to kp:        %zu (%.1f%%)\n", nkp, (double) nkp / (double) data_in->n * 100.0);
  fprintf(stderr, "\t\t flagged data due to dRC/dt:    %zu (%.1f%%)\n", nRC, (double) nRC / (double) data_in->n * 100.0);

#if 0
  fprintf(stderr, "preprocess_data: calculating virtual observatory magnetic equator positions...");
  preprocess_vobs(data_in, data_out);
  fprintf(stderr, "done (%zu data computed)\n", data_out->n);
#else
  /* data_out = data_in */
  copy_eej_data(data_in, data_out);
#endif

  return 0;
}

static int
parse_config_file(const char *filename, preprocess_parameters *params)
{
  int s;
  config_t cfg;
  double fval;

  config_init(&cfg);

  s = config_read_file(&cfg, filename);
  if (s != CONFIG_TRUE)
    {
      fprintf(stderr, "parse_config_file: %s:%d - %s\n",
              config_error_file(&cfg),
              config_error_line(&cfg),
              config_error_text(&cfg));
      config_destroy(&cfg);
      return -1;
    }

  if (config_lookup_float(&cfg, "max_kp", &fval))
    params->max_kp = fval;
  if (config_lookup_float(&cfg, "max_dRC", &fval))
    params->max_dRC = fval;

  config_destroy(&cfg);

  return 0;
}

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options]\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --eej_file file         | -i file          - EEJ magnetic equator data file from gradient search analysis\n");
  fprintf(stderr, "\t --eej_smooth_file file  | -s file          - EEJ magnetic equator data file smoothed with periodic B-spline\n");
  fprintf(stderr, "\t --eej_wavelet_file file | -w file          - EEJ magnetic equator data file from wavelet analysis\n");
  fprintf(stderr, "\t --config_file           | -C config_file   - configuration file\n");
  fprintf(stderr, "\t --output_file           | -o output_file   - output file in magdata format\n");
}

int
main(int argc, char *argv[])
{
  int status;
  char *config_file = "EEJ_preproc.cfg";
  eej_data * data_in = malloc(sizeof(eej_data));
  eej_data * data_out = malloc(sizeof(eej_data));
  struct timeval tv0, tv1;
  preprocess_parameters params;
  size_t magdata_flags = 0;           /* MAGDATA_GLOBFLG_xxx */
  magdata * mdata;
  char * output_file = NULL;

  /* initialize parameters */
  params.max_kp = -1.0;
  params.max_dRC = -1.0;

  data_in->n = 0;
  data_out->n = 0;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "eej_file", required_argument, NULL, 'i' },
          { "eej_smooth_file", required_argument, NULL, 's' },
          { "config_file", required_argument, NULL, 'C' },
          { "output_file", required_argument, NULL, 'o' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "C:i:o:s:w:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'C':
            config_file = optarg;
            break;

          case 'i':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);
            eej_read(optarg, data_in);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%g seconds, %zu data read)\n",
                    time_diff(tv0, tv1), data_in->n);
            break;

          case 's':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);
            eej_read_smooth(optarg, data_in);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%g seconds, %zu data read)\n",
                    time_diff(tv0, tv1), data_in->n);
            break;

          case 'w':
            fprintf(stderr, "main: reading %s...", optarg);
            gettimeofday(&tv0, NULL);
            eej_read_wavelet(optarg, data_in);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%g seconds, %zu data read)\n",
                    time_diff(tv0, tv1), data_in->n);
            break;

          case 'o':
            output_file = optarg;
            break;

          default:
            break;
        }
    }

  /* parse configuration file */
  parse_config_file(config_file, &params);

  /* check parameters */
  status = check_parameters(&params);
  if (status)
    exit(1);

  if (data_in->n == 0)
    {
      print_help(argv);
      exit(1);
    }

  fprintf(stderr, "main: === PREPROCESSING EEJ DATA ===\n");
  preprocess_data(&params, data_in, data_out);

  if (data_out->n == 0)
    {
      fprintf(stderr, "main: error: no data available\n");
      exit(1);
    }

  mdata = copy_data(magdata_flags, data_out);

  if (output_file)
    {
      fprintf(stderr, "main: writing data to %s...", output_file);
      magdata_write(output_file, mdata);
      fprintf(stderr, "done (%zu data written)\n", mdata->n);
    }

  magdata_free(mdata);
  free(data_in);
  free(data_out);

  return 0;
}
