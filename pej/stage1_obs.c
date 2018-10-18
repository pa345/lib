/*
 * stage1_obs.c
 *
 * Select observatory data according to IMF critera
 *
 * Steps are:
 * 1. Read observatory dataset
 * 2. Select data for IMF criteria
 * 3. convert to magdata format and store to disk
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

#include <obsdata/obsdata.h>
#include <indices/indices.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_statistics.h>

#include <common/bin3d.h>
#include <common/common.h>
#include <msynth/msynth.h>
#include <track/track.h>

#include "magdata.h"

typedef struct
{
  int all;            /* print all tracks */
  double kp_min;      /* minimum kp */
  double kp_max;      /* maximum kp */
  double IMF_Bz_min;  /* minimum IMF B_z (nT) */
  double IMF_Bz_max;  /* maximum IMF B_z (nT) */
  size_t downsample;  /* downsampling factor */
} preprocess_parameters;

/* subtract internal and external field models from data */
int
subtract_model(obsdata * data)
{
  int s = 0;
  size_t max_threads = (size_t) omp_get_max_threads();
  msynth_workspace **msynth_p = malloc(max_threads * sizeof(msynth_workspace *));
  size_t i;

  for (i = 0; i < max_threads; ++i)
    {
      msynth_p[i] = msynth_swarm_read(MSYNTH_CHAOS_FILE);
      msynth_set(1, 15, msynth_p[i]);
    }

#pragma omp parallel for private(i)
  for (i = 0; i < data->nstation; ++i)
    {
      int thread_id = omp_get_thread_num();
      obsdata_station *station = data->stations[i];
      double r = station->radius;
      double theta = M_PI / 2.0 - station->latitude * M_PI / 180.0;
      double phi = station->longitude * M_PI / 180.0;
      size_t j;

      for (j = 0; j < station->n; ++j)
        {
          double tyr = epoch2year(station->t[j]);
          double B_int[4];

          if (!OBSDATA_AvailableData(station->flags[j]))
            continue;

          msynth_eval(tyr, r, theta, phi, B_int, msynth_p[thread_id]);

          if (OBSDATA_ExistX(station->flags[j]))
            station->X[j] -= B_int[0] + station->X_ext[j];

          if (OBSDATA_ExistY(station->flags[j]))
            station->Y[j] -= B_int[1] + station->Y_ext[j];

          if (OBSDATA_ExistZ(station->flags[j]))
            station->Z[j] -= B_int[2] + station->Z_ext[j];
        }
    }

  for (i = 0; i < max_threads; ++i)
    msynth_free(msynth_p[i]);

  free(msynth_p);

  return s;
}

/*
preprocess_data()

Inputs: params - preprocess parameters
          lt_min     - minimum local time
          lt_max     - maximum local time
          downsample - downsampling factor
        data   - observatory data

Return: pointer to sorted track workspace (should be freed by caller)
*/

int
preprocess_data(const preprocess_parameters *params, obsdata *data)
{
  int s = 0;
  struct timeval tv0, tv1;

#if 0
  /* subtract internal and external field model from measurements */
  fprintf(stderr, "preprocess_data: subtracting internal/external models from data...");
  gettimeofday(&tv0, NULL);

  s = subtract_model(data);
  if (s)
    return s;

  gettimeofday(&tv1, NULL);
  fprintf(stderr, "done (%g seconds)\n", time_diff(tv0, tv1));
#endif

  if (params->all)
    return s;

  return s;
}

magdata *
copy_data(const obsdata *data)
{
  size_t ndata = obsdata_n(data);
  size_t npts[6] = { 0, 0, 0, 0, 0, 0 };
  magdata_params params;
  magdata *mdata = magdata_alloc(ndata, R_EARTH_KM);
  size_t i;

  if (!mdata)
    return 0;

  params.model_main = 0;
  params.model_crust = 0;
  params.model_ext = 0;

  fprintf(stderr, "\n");
  fprintf(stderr, "\t copy_data: copying stations to magdata format...");

  for (i = 0; i < data->nstation; ++i)
    {
      obsdata_station *station = data->stations[i];

      magdata_copy_station(&params, station, mdata, npts);
    }

  fprintf(stderr, "done (ndata = %zu mdata_n = %zu, mdata_ntot = %zu)\n", ndata, mdata->n, mdata->ntot);

  return mdata;
}

void
print_help(char *argv[])
{
  fprintf(stderr, "Usage: %s [options]\n", argv[0]);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t -b                                          - read BGS binary file for ionospheric field\n");
  fprintf(stderr, "\t --output_file | -o output_file              - output file\n");
  fprintf(stderr, "\t --kp_min | -v kp_min                        - kp minimum\n");
  fprintf(stderr, "\t --kp_max | -w kp_max                        - kp maximum\n");
  fprintf(stderr, "\t --imf_bz_min | -A imf_bz_min                - minimum IMF B_z (nT)\n");
  fprintf(stderr, "\t --imf_bz_max | -B imf_bz_max                - maximum IMF B_z (nT)\n");
}

int
main(int argc, char *argv[])
{
  char *data_file = "data/data.dat";
  obsdata *data = NULL;
  struct timeval tv0, tv1;
  preprocess_parameters params;
  magdata *mdata;

  /* defaults */
  params.kp_min = 0.0;
  params.kp_max = 20.0;
  params.IMF_Bz_min = -1.0e6;
  params.IMF_Bz_max = 1.0e6;

  while (1)
    {
      int c;
      int option_index = 0;
      static struct option long_options[] =
        {
          { "kp_min", required_argument, NULL, 'v' },
          { "kp_max", required_argument, NULL, 'w' },
          { "output_file", required_argument, NULL, 'o' },
          { "imf_bz_min", required_argument, NULL, 'A' },
          { "imf_bz_max", required_argument, NULL, 'B' },
          { 0, 0, 0, 0 }
        };

      c = getopt_long(argc, argv, "o:A:B:b", long_options, &option_index);
      if (c == -1)
        break;

      switch (c)
        {
          case 'b':
            fprintf(stderr, "main: reading %s...", OBSDATA_BINARY_IONO_FILE);
            gettimeofday(&tv0, NULL);
            data = obsdata_read(OBSDATA_BINARY_IONO_FILE);
            gettimeofday(&tv1, NULL);
            fprintf(stderr, "done (%g seconds, %zu stations read, %zu total measurements)\n",
                    time_diff(tv0, tv1), data->nstation, obsdata_n(data));
            break;

          case 'o':
            data_file = optarg;
            break;

          case 'v':
            params.kp_min = atof(optarg);
            break;

          case 'w':
            params.kp_max = atof(optarg);
            break;

          case 'A':
            params.IMF_Bz_min = atof(optarg);
            break;

          case 'B':
            params.IMF_Bz_max = atof(optarg);
            break;

          default:
            break;
        }
    }

  if (!data)
    {
      print_help(argv);
      exit(1);
    }

  fprintf(stderr, "main: kp minimum       = %.1f\n", params.kp_min);
  fprintf(stderr, "main: kp maximum       = %.1f\n", params.kp_max);
  fprintf(stderr, "main: IMF B_z minimum  = %.1f\n", params.IMF_Bz_min);
  fprintf(stderr, "main: IMF B_z maximum  = %.1f\n", params.IMF_Bz_max);

  preprocess_data(&params, data);

  {
#if 1
    obsdata_station *station = obsdata_station_find("KOU0", data);
    obsdata_station_print("kou.txt", station, data);
#else
    obsdata_station *station = obsdata_station_find("RES0", data);
    obsdata_station_print("res.txt", station);
#endif
  }

  fprintf(stderr, "main: converting to magdata format...");
  mdata = copy_data(data);
  fprintf(stderr, "done\n");

  fprintf(stderr, "main: writing data to %s...", data_file);
  magdata_write(data_file, mdata);
  fprintf(stderr, "done\n");

  magdata_free(mdata);
  obsdata_free(data);

  return 0;
}
