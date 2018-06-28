/*
 * mfield_preproc_obs.c
 *
 * Pre-process observatory data and save in magdata format.
 *
 * Pre-processing steps are:
 * 1. Instrument flags (recommended CHAMP flags except 1 star camera allowed)
 * 2. Track rms test
 * 3. filter for WMM criteria
 * 4. Downsample by factor 20
 * 5. Compute and store along-track gradients
 *
 * The result is an output file in magdata format containing all data point to
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
#include <flow/flow.h>
#include <indices/indices.h>
#include <obsdata/obsdata.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>

#include <common/common.h>
#include <msynth/msynth.h>

#include "magdata.h"
#include "magfit.h"
#include "track.h"

#include "mfield.h"
#include "mfield_data.h"

typedef struct
{
  size_t downsample;      /* downsampling factor */
  double min_LT;          /* minimum local time for field modeling */
  double max_LT;          /* maximum local time for field modeling */
  double qdlat_preproc_cutoff; /* QD latitude cutoff for high-latitudes */
  double min_zenith;      /* minimum zenith angle for high-latitude data selection */

  int subtract_B_main;    /* subtract a-priori main field from data */
  int subtract_B_crust;   /* subtract a-priori crustal field from data */
  int subtract_B_ext;     /* subtract a-priori external field from data */

  double max_kp;          /* maximum kp */
  double max_dRC;         /* maximum dRC/dt (nT/hour) */
} preprocess_parameters;

static int mfield_check_LT(const double lt, const double lt_min, const double lt_max);
static size_t model_flags(const size_t magdata_flags, const double t,
                          const double theta, const double phi, const double qdlat,
                          const preprocess_parameters * params);

#define MFIELD_IDX_X              0
#define MFIELD_IDX_Y              1
#define MFIELD_IDX_Z              2
#define MFIELD_IDX_F              3
#define MFIELD_IDX_DX_NS          4
#define MFIELD_IDX_DY_NS          5
#define MFIELD_IDX_DZ_NS          6
#define MFIELD_IDX_DF_NS          7
#define MFIELD_IDX_DX_EW          8
#define MFIELD_IDX_DY_EW          9
#define MFIELD_IDX_DZ_EW          10
#define MFIELD_IDX_DF_EW          11
#define MFIELD_IDX_B_EULER        12
#define MFIELD_IDX_END            13

static int
check_parameters(preprocess_parameters * params)
{
  int s = 0;

  if (params->downsample == 0)
    {
      fprintf(stderr, "check_parameters: downsample must be > 0\n");
      ++s;
    }

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

  if (params->min_LT < 0.0)
    {
      fprintf(stderr, "check_parameters: min_LT must be > 0\n");
      ++s;
    }

  if (params->max_LT < 0.0)
    {
      fprintf(stderr, "check_parameters: max_LT must be > 0\n");
      ++s;
    }

  if (params->qdlat_preproc_cutoff < 0.0)
    {
      fprintf(stderr, "check_parameters: qdlat_preproc_cutoff must be > 0\n");
      ++s;
    }

  if (params->min_zenith < 0.0)
    {
      fprintf(stderr, "check_parameters: min_zenith must be > 0\n");
      ++s;
    }

  if (params->subtract_B_main < 0)
    {
      fprintf(stderr, "check_parameters: subtract_B_main must be 0 or 1\n");
      ++s;
    }

  if (params->subtract_B_crust < 0)
    {
      fprintf(stderr, "check_parameters: subtract_B_crust must be 0 or 1\n");
      ++s;
    }

  if (params->subtract_B_ext < 0)
    {
      fprintf(stderr, "check_parameters: subtract_B_ext must be 0 or 1\n");
      ++s;
    }

  return s;
}

/* convert observatory SV data to magdata format */
magdata *
copy_data_SV(const size_t magdata_flags, const obsdata_station *station, preprocess_parameters * preproc_params)
{
  size_t ndata = station->n_sv_tot;
  magdata *mdata;
  magdata_params params;
  size_t npts[6] = { 0, 0, 0, 0, 0, 0 };
  size_t i;

  if (ndata == 0)
    return NULL;

  params.model_main = 0;
  params.model_crust = 0;
  params.model_ext = 0;

  mdata = magdata_alloc(ndata, R_EARTH_KM);
  if (!mdata)
    return 0;

  mdata->global_flags = magdata_flags;

  magdata_copy_station_SV(&params, station, mdata, npts);

  /* now determine which points in mdata will be used for MF modeling */
  for (i = 0; i < mdata->n; ++i)
    {
      size_t fitting_flags = 0;

      if (fabs(mdata->qdlat[i]) <= preproc_params->qdlat_preproc_cutoff)
        {
          /* mid-latitude point: check local time to determine whether to fit MF model */

          double LT = mdata->lt[i];
          int fit_MF = mfield_check_LT(LT, preproc_params->min_LT, preproc_params->max_LT);

          if (fit_MF)
            fitting_flags |= MAGDATA_FLG_FIT_MF;
        }
      else
        {
          /* high-latitude point - fit field model */
          fitting_flags |= MAGDATA_FLG_FIT_MF;
        }

      mdata->flags[i] |= fitting_flags;
    }

  return mdata;
}

/* convert observatory daily mean data to magdata format */
magdata *
copy_data(const size_t magdata_flags, const obsdata_station *station, preprocess_parameters * preproc_params)
{
  size_t ndata = station->n_mean_tot;
  magdata *mdata;
  magdata_params params;
  size_t npts[6] = { 0, 0, 0, 0, 0, 0 };
  size_t i;

  if (ndata == 0)
    return NULL;

  params.model_main = 0;
  params.model_ext = 0;

  /* subtract crustal field biases from data prior to modeling */
  if (preproc_params->subtract_B_crust)
    params.model_crust = 1;
  else
    params.model_crust = 0;

  mdata = magdata_alloc(ndata, R_EARTH_KM);
  if (!mdata)
    return 0;

  mdata->global_flags = magdata_flags;

  magdata_copy_station(&params, station, mdata, npts);

  /* now determine which points in mdata will be used for MF modeling */
  for (i = 0; i < mdata->n; ++i)
    {
      size_t fitting_flags = 0;

      if (fabs(mdata->qdlat[i]) <= preproc_params->qdlat_preproc_cutoff)
        {
          /* mid-latitude point: check local time to determine whether to fit MF model */

          double LT = mdata->lt[i];
          int fit_MF = mfield_check_LT(LT, preproc_params->min_LT, preproc_params->max_LT);

          if (fit_MF)
            fitting_flags |= MAGDATA_FLG_FIT_MF;
        }
      else
        {
          /* high-latitude point - fit field model */
          fitting_flags |= MAGDATA_FLG_FIT_MF;
        }

      mdata->flags[i] |= fitting_flags;
    }

  return mdata;
}

/*
mfield_check_LT()
  Check if a given LT is within [lt_min,lt_max] accounting for mod 24. So it
is possible to have input lt_min < lt_max in order to select data across midnight.

Example: [lt_min,lt_max] = [6,18] will select daytime data between 6am and 6pm
         [lt_min,lt_max] = [18,6] will select nighttime data between 6pm and 6am
         [lt_min,lt_max] = [22,5] will select nighttime data between 10pm and 5am
         [lt_min,lt_max] = [0,5] will select data between midnight and 5am

Return: 1 if LT \in [lt_min,lt_max] (mod 24); 0 otherwise
*/

static int
mfield_check_LT(const double lt, const double lt_min, const double lt_max)
{
  double a, b;

  b = fmod(lt_max - lt_min, 24.0);
  if (b < 0.0)
    b += 24.0;

  a = fmod(lt - lt_min, 24.0);
  if (a < 0.0)
    a += 24.0;

  if (a > b)
    return 0; /* invalid local time */

  /* valid local time */
  return 1;
}

/*
preprocess_data()
  Compute daily means and SV data for observatories

Inputs: params - preprocess parameters
          downsample - downsampling factor
        data   - observatory data

Return: pointer to sorted track workspace (should be freed by caller)
*/

int
preprocess_data(const preprocess_parameters *params, const size_t magdata_flags,
                obsdata *data)
{
  const double interval_SV = 0.5 * 365.25; /* interval for SV differences (days) */
  const size_t min_samples = 1;            /* minimum samples needed for daily means */
  obsdata_select_params obs_params = obsdata_select_default_params();
  struct timeval tv0, tv1;
  size_t i;

  fprintf(stderr, "preprocess_data: selecting quiet time data...");

  for (i = 0; i < data->nstation; ++i)
    {
      obsdata_station * station = data->stations[i];
      size_t nflagged_array[OBSDATA_IDX_END];
      size_t nflagged;

      fprintf(stderr, "\t selecting %s observatory data for quiet conditions...", station->name);
      gettimeofday(&tv0, NULL);
      nflagged = obsdata_station_select(&obs_params, station, data, nflagged_array);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds, %zu/%zu (%.1f%%) data flagged)\n",
              time_diff(tv0, tv1), nflagged, station->n, (double) nflagged / (double) station->n * 100.0);

      fprintf(stderr, "\t\t flagged data due to LT:        %zu (%.1f%%)\n", nflagged_array[OBSDATA_IDX_LT], (double) nflagged_array[OBSDATA_IDX_LT] / (double) station->n * 100.0);
      fprintf(stderr, "\t\t flagged data due to kp:        %zu (%.1f%%)\n", nflagged_array[OBSDATA_IDX_KP], (double) nflagged_array[OBSDATA_IDX_KP] / (double) station->n * 100.0);
      fprintf(stderr, "\t\t flagged data due to dRC/dt:    %zu (%.1f%%)\n", nflagged_array[OBSDATA_IDX_DRC], (double) nflagged_array[OBSDATA_IDX_DRC] / (double) station->n * 100.0);
      fprintf(stderr, "\t\t flagged data due to SMDL:      %zu (%.1f%%)\n", nflagged_array[OBSDATA_IDX_SMDL], (double) nflagged_array[OBSDATA_IDX_SMDL] / (double) station->n * 100.0);
      fprintf(stderr, "\t\t flagged data due to d/dt SMDL: %zu (%.1f%%)\n", nflagged_array[OBSDATA_IDX_DSMDL], (double) nflagged_array[OBSDATA_IDX_DSMDL] / (double) station->n * 100.0);

      fprintf(stderr, "\t computing %s daily means...", station->name);
      gettimeofday(&tv0, NULL);
      obsdata_station_daily(station);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds, %zu means computed)\n",
              time_diff(tv0, tv1), station->n_mean_tot);

      fprintf(stderr, "\t computing %s SV values...", station->name);
      gettimeofday(&tv0, NULL);
      obsdata_station_SV(min_samples, interval_SV, station);
      gettimeofday(&tv1, NULL);
      fprintf(stderr, "done (%g seconds, %zu SV values computed)\n",
              time_diff(tv0, tv1), station->n_sv_tot);
    }

  fprintf(stderr, "done\n");

#if 0
  /* select geomagnetically quiet data */
  fprintf(stderr, "preprocess_data: selecting geomagnetically quiet data...");
  mfield_preprocess_filter(magdata_flags, params, track_p, data);
  fprintf(stderr, "done\n");
#endif

#if 0
  /* downsample data */
  {
    size_t i;

    fprintf(stderr, "preprocess_data: downsampling data by factor %zu...", params->downsample);

    for (i = 0; i < data->n; ++i)
      {
        if (i % params->downsample != 0)
          data->flags[i] |= SATDATA_FLG_DOWNSAMPLE;
      }

    fprintf(stderr, "done\n");
  }
#endif

  return 0;
}

static int
parse_config_file(const char *filename, preprocess_parameters *params)
{
  int s;
  config_t cfg;
  double fval;
  int ival;

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
  if (config_lookup_float(&cfg, "min_LT", &fval))
    params->min_LT = fval;
  if (config_lookup_float(&cfg, "max_LT", &fval))
    params->max_LT = fval;
  if (config_lookup_float(&cfg, "qdlat_preproc_cutoff", &fval))
    params->qdlat_preproc_cutoff = fval;
  if (config_lookup_float(&cfg, "min_zenith", &fval))
    params->min_zenith = fval;

  if (config_lookup_int(&cfg, "downsample", &ival))
    params->downsample = (size_t) ival;

  if (config_lookup_int(&cfg, "subtract_B_main", &ival))
    params->subtract_B_main = (size_t) ival;
  if (config_lookup_int(&cfg, "subtract_B_crust", &ival))
    params->subtract_B_crust = (size_t) ival;
  if (config_lookup_int(&cfg, "subtract_B_ext", &ival))
    params->subtract_B_ext = (size_t) ival;

  config_destroy(&cfg);

  return 0;
}

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options]\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t --obs_data_file   | -O                        - Observatory data file (%s)\n", OBSDATA_BINARY_FILE);
  fprintf(stderr, "\t --downsample      | -d downsample             - downsampling factor\n");
  fprintf(stderr, "\t --config_file     | -C config_file            - configuration file\n");
}

int
main(int argc, char *argv[])
{
  int status;
  const char *path_dir = "/data/palken/lib/mfield/data/obs";
  char *datamap_file = "datamap.dat";
  char *data_file = "data.dat";
  char *config_file = "OBS.cfg";
  obsdata *data = NULL;
  struct timeval tv0, tv1;
  preprocess_parameters params;
  size_t downsample = 0;     /* downsample factor */
  size_t magdata_flags = 0;  /* MAGDATA_GLOBFLG_xxx */
  int flag_vec_rms = 1;
  char output_file[2048];
  size_t i;

  /* initialize parameters */
  params.downsample = 0;
  params.max_kp = -1.0;
  params.max_dRC = -1.0;
  params.min_LT = -1.0;
  params.max_LT = -1.0;
  params.qdlat_preproc_cutoff = -1.0;
  params.min_zenith = -1.0;
  params.subtract_B_main = -1;
  params.subtract_B_crust = -1;
  params.subtract_B_ext = -1;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "obs_file", no_argument, NULL, 'O' },
          { "downsample", required_argument, NULL, 'd' },
          { "config_file", required_argument, NULL, 'C' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "C:d:O", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'C':
            config_file = optarg;
            break;

          case 'd':
            downsample = (size_t) atoi(optarg);
            break;

          case 'O':
            fprintf(stderr, "main: reading %s...", OBSDATA_BINARY_FILE);
            gettimeofday(&tv0, NULL);
            data = obsdata_read(OBSDATA_BINARY_FILE);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%g seconds, %zu stations read, %zu total measurements)\n",
                    time_diff(tv0, tv1), data->nstation, obsdata_n(data));
            break;

          default:
            break;
        }
    }

  /* parse configuration file */
  parse_config_file(config_file, &params);

  /* replace config values with command-line arguments */
  if (downsample > 0)
    params.downsample = downsample;

  /* check parameters */
  status = check_parameters(&params);
  if (status)
    exit(1);

  if (!data)
    {
      print_help(argv);
      exit(1);
    }

  fprintf(stderr, "main: downsample       = %zu\n", params.downsample);
  fprintf(stderr, "main: LT minimum       = %.1f\n", params.min_LT);
  fprintf(stderr, "main: LT maximum       = %.1f\n", params.max_LT);

  fprintf(stderr, "main: === PREPROCESSING OBSERVATORY DATA ===\n");
  preprocess_data(&params, magdata_flags, data);

  for (i = 0; i < data->nstation; ++i)
    {
      obsdata_station *station = data->stations[i];
      magdata *mdata;

      /* copy SV data */

      mdata = copy_data_SV(magdata_flags, station, &params);

      if (mdata == NULL)
        {
          fprintf(stderr, "main: WARNING: station %s has no usable SV data\n", station->name);
          continue;
        }

      sprintf(output_file, "%s/%s_SV.dat", path_dir, station->name);
      fprintf(stderr, "main: writing data to %s...", output_file);
      magdata_write(output_file, mdata);
      fprintf(stderr, "done\n");

      magdata_free(mdata);

      /* copy daily mean data */
      mdata = copy_data(magdata_flags, station, &params);

      if (mdata == NULL)
        {
          fprintf(stderr, "main: WARNING: station %s has no usable daily mean data\n", station->name);
          continue;
        }

      sprintf(output_file, "%s/%s.dat", path_dir, station->name);
      fprintf(stderr, "main: writing data to %s...", output_file);
      magdata_write(output_file, mdata);
      fprintf(stderr, "done\n");

      magdata_free(mdata);
    }

  /* free data after copying arrays to free up memory */
  obsdata_free(data);

  exit(1);

#if 0
  magdata_init(mdata);

  fprintf(stderr, "main: computing spatial weighting of data...");
  gettimeofday(&tv0, NULL);
  magdata_calc(mdata);
  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));

#if 0
  fprintf(stderr, "main: writing data to %s...", data_file);
  magdata_print(data_file, mdata);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: writing data map to %s...", datamap_file);
  magdata_map(datamap_file, mdata);
  fprintf(stderr, "done\n");
#endif

  fprintf(stderr, "main: observatory rmin = %.1f (%.1f) [km]\n",
          mdata->rmin, mdata->rmin - mdata->R);
  fprintf(stderr, "main: observatory rmax = %.1f (%.1f) [km]\n",
          mdata->rmax, mdata->rmax - mdata->R);

  if (output_file)
    {
      fprintf(stderr, "main: writing data to %s...", output_file);
      magdata_write(output_file, mdata);
      fprintf(stderr, "done\n");
    }

  magdata_free(mdata);
#endif

  return 0;
}
