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

#include <flow/flow.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>

#include <mainlib/ml_apex.h>
#include <mainlib/ml_indices.h>
#include <mainlib/ml_obsdata.h>
#include <mainlib/ml_common.h>
#include <mainlib/ml_msynth.h>
#include <mainlib/ml_magdata.h>

#include "mfield.h"
#include "mfield_data.h"

typedef struct
{
  double qdlat_preproc_cutoff; /* QD latitude cutoff for high-latitudes */

  int subtract_B_main;    /* subtract a-priori main field from data */
  int subtract_B_crust;   /* subtract a-priori crustal field from data */
  int subtract_B_ext;     /* subtract a-priori external field from data */

  double max_kp;          /* maximum kp */
  double max_dRC;         /* maximum dRC/dt (nT/hour) */
} preprocess_parameters;

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

  if (params->qdlat_preproc_cutoff < 0.0)
    {
      fprintf(stderr, "check_parameters: qdlat_preproc_cutoff must be > 0\n");
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
copy_data_SV(const size_t magdata_flags, const obsdata_station *station)
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

  /*
   * the calculation of daily means and SV values already performed extensive data
   * selection (LT, MLT, etc). So just flag all mean/SV values for main field modeling
   */
  for (i = 0; i < mdata->n; ++i)
    mdata->flags[i] |= MAGDATA_FLG_FIT_MF;

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

  /* copy daily means into magdata struct */
  magdata_copy_station_means(&params, station, mdata, npts);

  /*
   * the calculation of daily means and SV values already performed extensive data
   * selection (LT, MLT, etc). So just flag all mean/SV values for main field modeling
   */
  for (i = 0; i < mdata->n; ++i)
    mdata->flags[i] |= MAGDATA_FLG_FIT_MF;

  return mdata;
}

/*
preprocess_data()
  Compute daily means and SV data for observatories

Inputs: params - preprocess parameters
        data   - observatory data

Return: pointer to sorted track workspace (should be freed by caller)
*/

int
preprocess_data(const preprocess_parameters *params, obsdata *data)
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
      fprintf(stderr, "\t\t flagged data due to MLT:       %zu (%.1f%%)\n", nflagged_array[OBSDATA_IDX_MLT], (double) nflagged_array[OBSDATA_IDX_MLT] / (double) station->n * 100.0);
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
  if (config_lookup_float(&cfg, "qdlat_preproc_cutoff", &fval))
    params->qdlat_preproc_cutoff = fval;

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
  fprintf(stderr, "\t --obs_data_file   | -O obs_binary_file        - Observatory data file (%s)\n", OBSDATA_BINARY_FILE);
  fprintf(stderr, "\t --config_file     | -C config_file            - configuration file\n");
}

int
main(int argc, char *argv[])
{
  int status;
  const char *path_dir = "/data/palken/lib/mfield/data/obs";
  char *datamap_file = "datamap.dat";
  char *config_file = "MF_preproc_OBS.cfg";
  char *obs_file = OBSDATA_BINARY_FILE;
  obsdata *data = NULL;
  struct timeval tv0, tv1;
  preprocess_parameters params;
  size_t magdata_flags = 0;  /* MAGDATA_GLOBFLG_xxx */
  char output_file[2048];
  size_t i;

  /* initialize parameters */
  params.max_kp = -1.0;
  params.max_dRC = -1.0;
  params.qdlat_preproc_cutoff = -1.0;
  params.subtract_B_main = -1;
  params.subtract_B_crust = -1;
  params.subtract_B_ext = -1;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "obs_file", required_argument, NULL, 'O' },
          { "config_file", required_argument, NULL, 'C' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "C:O:", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'C':
            config_file = optarg;
            break;

          case 'O':
            obs_file = optarg;
            break;

          default:
            print_help(argv);
            exit(1);
            break;
        }
    }

  /* parse configuration file */
  fprintf(stderr, "main: parsing configuration file %s...", config_file);
  parse_config_file(config_file, &params);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: reading %s...", obs_file);
  gettimeofday(&tv0, NULL);
  data = obsdata_read(obs_file);
  gettimeofday(&tv1, NULL);

  if (data == NULL)
    exit(1);

  fprintf(stderr, "done (%g seconds, %zu stations read, %zu total measurements)\n",
          time_diff(tv0, tv1), data->nstation, obsdata_n(data));

  /* check parameters */
  status = check_parameters(&params);
  if (status)
    exit(1);

  if (!data)
    {
      print_help(argv);
      exit(1);
    }

  fprintf(stderr, "main: === PREPROCESSING OBSERVATORY DATA ===\n");
  preprocess_data(&params, data);

  for (i = 0; i < data->nstation; ++i)
    {
      obsdata_station *station = data->stations[i];
      magdata *mdata;

      /* copy SV data */

      mdata = copy_data_SV(magdata_flags, station);

      if (mdata != NULL)
        {
          sprintf(output_file, "%s/%s_SV.dat", path_dir, station->name);
          fprintf(stderr, "main: writing data to %s...", output_file);
          magdata_write(output_file, mdata);
          fprintf(stderr, "done\n");

          magdata_free(mdata);
        }
      else
        {
          fprintf(stderr, "main: WARNING: station %s has no usable SV data\n", station->name);
        }

      /* copy daily mean data */
      mdata = copy_data(magdata_flags, station, &params);

      if (mdata != NULL)
        {
          sprintf(output_file, "%s/%s.dat", path_dir, station->name);
          fprintf(stderr, "main: writing data to %s...", output_file);
          magdata_write(output_file, mdata);
          fprintf(stderr, "done\n");

          magdata_free(mdata);
        }
      else
        {
          fprintf(stderr, "main: WARNING: station %s has no usable daily mean data\n", station->name);
        }
    }

  /* free data after copying arrays to free up memory */
  obsdata_free(data);

  return 0;
}
